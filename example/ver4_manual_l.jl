using Printf, LinearAlgebra
using CompositeGrids
using ElectronGas
using Parameters, Random
using MCIntegration
using Lehmann

using ExpressionTree
using StaticArrays

Parquet = GWKT.Parquet
Manual = GWKT.Manual
DiagTree = GWKT.DiagTree

include("parameter.jl")
include("interaction.jl")


const steps = 1e7
const isF = true
const isProper = true #one interaction irreduble diagrams or not
const hasBubble = false #allow the bubble diagram or not
# const θgrid = collect(LinRange(0.1, π, Nk)) # external angle grid
const lgrid = [1, 2]
const Nk = length(lgrid)
# const ExtK = [@SVector [kF * cos(θ), kF * sin(θ), 0.0] for θ in θgrid]
const RefK = @SVector [kF, 0.0, 0.0]

struct Para{Q,T}
    dW0::Matrix{Float64}
    qgrid::Q
    τgrid::T # dedicated τgrid for dynamic interaction
    diag::Vector{DiagTree.Diagrams{Float64}}
    dir::Vector{Int}
    ex::Vector{Int}
    function Para(loopOrder::Int)
        qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 6 * kF], [0.0, 2kF], 16, 0.01 * kF, 8)
        τgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 16, β * 1e-4, 8)

        vqinv = [(q^2 + mass2) / (4π * e0^2) for q in qgrid.grid]
        dW0 = TwoPoint.dWRPA(vqinv, qgrid.grid, τgrid.grid, dim, EF, kF, β, spin, me) # dynamic part of the effective interaction

        # chan = [1, 2, 3]
        KinL = KoutL = [1, 0, 0]
        KinR = KoutR = [0, 1, 0]
        legK = [KinL, KoutL, KinR, KoutR]
        Gsym = [:mirror]
        Wsym = [:mirror, :timereversal]
        diag1, dir1, ex1 = Manual.build([1, 2, 3], legK, 3, spin, isProper, hasBubble, Gsym, Wsym)
        diag2, dir2, ex2 = Manual.build([2,], legK, 3, spin, isProper, hasBubble, Gsym, Wsym)
        diag3, dir3, ex3 = Manual.build([3,], legK, 3, spin, isProper, hasBubble, Gsym, Wsym)
        diag = [diag1, diag2, diag3]
        dir = [dir1, dir2, dir3]
        ex = [ex1, ex2, ex3]
        # DiagTree.showTree(diag1, dir1)
        # DiagTree.showTree(diag3, dir3)

        return new{typeof(qgrid),typeof(τgrid)}(dW0, qgrid, τgrid, diag, dir, ex)
    end
end

@inline function phase(varK, varT, extK, extT)
    tInL, tOutL, tInR, tOutR = varT[extT[INL]], varT[extT[OUTL]], varT[extT[INR]],
    varT[extT[OUTR]]
    if (isF)
        return cos(π / β * ((tInL + tOutL) - (tInR + tOutR)))
    else
        return cos(π / β * ((tInL - tOutL) + (tInR - tOutR)))
    end
end

@inline function evalPropagator(type, K, Tidx, varT, factor, para)
    τin, τout = varT[Tidx[1]], varT[Tidx[2]]
    if type == 1 #G
        ϵ = (dot(K, K) - kF^2) / (2me)
        return Spectral.kernelFermiT(τout - τin, ϵ, β) * factor
    else
        v, w = interactionDynamic(para, K, τin, τout)
        # v, w = interactionStatic(para, K, τin, τout)
        if type == 2 #v
            return -v * factor
        elseif type == 3 #W
            return -w * factor
        else
            error("not implemented!")
        end
    end
end

function integrand(config)
    if config.curr == 1
        x = config.var[3][1]
        chan = config.var[4][1]
        K, varT = config.var[1], config.var[2]
        ExtK = @SVector [kF * x, kF * sqrt(1 - x^2), 0.0]
        varK = [RefK, ExtK, K[1]]
        para = config.para
        wd, we = DiagTree.evalNaive(para.diag[chan], evalPropagator, varK, varT, [para.dir[chan], para.ex[chan]], phase, para)
        factor = 1 / (2π)^dim / β
        wd *= factor
        we *= factor
        # weight = eval_TU(config, RefK, RefK, ExtK[extKidx], ExtK[extKidx], false)
        # @assert abs(weight.d + wd) < 1e-10 "wd: $(weight.d) != $wd"
        # @assert abs(weight.e + we) < 1e-10 "we: $(weight.e) != $we"
        return Weight(wd, we)
    else
        error("impossible!")
    end
end

function phaseT(tInL, tOutL, tInR, tOutR, isT)
    if isT == false
        tOutL, tOutR = tOutR, tOutL
    end

    if (isF)
        return cos(π * ((tInL + tOutL) - (tInR + tOutR)) / β)
    else
        return cos(π * ((tInL - tOutL) + (tInR - tOutR)) / β)
    end
end

# function eval_TUC(config, KInL, KOutL, KInR, KOutR, isT)
#     para = config.para
#     K, T = config.var[1], config.var[2]
#     Qd = isT ? KInL - KOutL : KInL - KOutR
#     k1, k2 = K[1], K[1] - Qd
#     t1, t2 = [T[1], T[2]], [T[3], T[4]] # t1, t2 both have two tau variables

#     vld, wld, vle, wle = vertexDynamic(para, Qd, Vector(KInL) - k1, t1[1], t1[2])
#     vrd, wrd, vre, wre = vertexDynamic(para, Qd, Vector(KInR) - k2, t2[1], t2[2])

#     # possible green's functions on the top
#     ϵ1, ϵ2 = (dot(k1, k1) - kF^2) / (2me), (dot(k2, k2) - kF^2) / (2me)

#     gt1 = Spectral.kernelFermiT(t2[1] - t1[1], ϵ1, β)

#     gd1 = Spectral.kernelFermiT(t1[1] - t2[1], ϵ2, β)
#     G = gt1 * gd1 / (2π)^3 * phaseT(t1[1], t1[1], t2[1], t2[1], isT)
#     we += G * (vle * vre)
# end

function eval_TU(config, KInL, KOutL, KInR, KOutR, isT)
    para = config.para
    K, T, Ang = config.var[1], config.var[2], config.var[3]
    Qd = isT ? KInL - KOutL : KInL - KOutR
    k1, k2 = K[1], K[1] - Qd
    t1, t2 = [T[1], T[2]], [T[3], T[4]] # t1, t2 both have two tau variables
    # θ = para.extAngle[Ang[1]] # angle of the external momentum on the right
    # KInR = [kF * cos(θ), kF * sin(θ), 0.0]
    # Qd = zero(k1)

    vld, wld, vle, wle = vertexDynamic(para, Qd, Vector(KInL) - k1, t1[1], t1[2])
    vrd, wrd, vre, wre = vertexDynamic(para, Qd, Vector(KInR) - k2, t2[1], t2[2])

    # vld, wld, vle, wle = vertexStatic(para, Qd, KInL - k1, t1[1], t1[2])
    # vrd, wrd, vre, wre = vertexStatic(para, Qd, KInR - k2, t2[1], t2[2])

    wd, we = 0.0, 0.0

    # possible green's functions on the top
    ϵ1, ϵ2 = (dot(k1, k1) - kF^2) / (2me), (dot(k2, k2) - kF^2) / (2me)

    gt1 = Spectral.kernelFermiT(t2[1] - t1[1], ϵ1, β)
    gt2 = Spectral.kernelFermiT(t1[1] - t2[1], ϵ2, β)
    # wd += 1.0 / β * 1.0 / β * gt1 * gt2 / (2π)^3 * phase(t1[1], t1[1], t2[1], t2[1])
    # println(k1)
    # println(k2)
    # println(t1)
    # println(t2)

    # wd += 1.0 / β * 1.0 / β * gt1 * gt2 / (2π)^3

    # gt3 = Spectral.kernelFermiT(t1[1] - t2[2], ϵ2, β)
    # G = gt1 * gt3 / (2π)^3 * phase(t1[1], t1[1], t2[2], t2[1])
    # wd += G * (vld * wre)

    # wd += spin * (vld + wld) * (vrd + wrd) * gt1 * gt2 / (2π)^3 * phase(t1[1], t1[1], t2[1], t2[1])
    # println(vld, ", ", wld, "; ", vrd, ", ", wld, ", ", wd)

    ############## Diagram v x v ######################
    """
      KInL                      KInR
       |                         | 
  t1.L ↑     t1.L       t2.L     ↑ t2.L
       |-------------->----------|
       |       |    k1    |      |
       |   ve  |          |  ve  |
       |       |    k2    |      |
       |--------------<----------|
  t1.L ↑    t1.L        t2.L     ↑ t2.L
       |                         | 
      KInL                      KInR
"""
    gd1 = Spectral.kernelFermiT(t1[1] - t2[1], ϵ2, β)
    G = gt1 * gd1 / (2π)^3 * phaseT(t1[1], t1[1], t2[1], t2[1], isT)
    we += G * (vle * vre)
    # println(G * (vle * vre) * (2π)^3)
    ##################################################

    ############## Diagram w x v ######################
    """
      KInL                      KInR
       |                         | 
  t1.R ↑     t1.L       t2.L     ↑ t2.L
       |-------------->----------|
       |       |    k1    |      |
       |   we  |          |  ve  |
       |       |    k2    |      |
       |--------------<----------|
  t1.L ↑    t1.R        t2.L     ↑ t2.L
       |                         | 
      KInL                      KInR
    """
    gd2 = Spectral.kernelFermiT(t1[2] - t2[1], ϵ2, β)
    G = gt1 * gd2 / (2π)^3 * phaseT(t1[1], t1[2], t2[1], t2[1], isT)
    we += G * (wle * vre)
    # println(G * (wle * vre) * (2π)^3)
    ##################################################

    ############## Diagram v x w ######################
    """
      KInL                      KInR
       |                         | 
  t1.L ↑     t1.L       t2.L     ↑ t2.L
       |-------------->----------|
       |       |    k1    |      |
       |   ve  |          |  we  |
       |       |    k2    |      |
       |--------------<----------|
  t1.L ↑    t1.L        t2.R     ↑ t2.R
       |                         | 
      KInL                      KInR
    """
    gd3 = Spectral.kernelFermiT(t1[1] - t2[2], ϵ2, β)
    G = gt1 * gd3 / (2π)^3 * phaseT(t1[1], t1[1], t2[2], t2[1], isT)
    we += G * (vle * wre)
    # println(G * (vle * wre) * (2π)^3)
    ##################################################

    ############## Diagram w x w ######################
    """
      KInL                      KInR
       |                         | 
  t1.R ↑     t1.L       t2.L     ↑ t2.L
       |-------------->----------|
       |       |    k1    |      |
       |   we  |          |  we  |
       |       |    k2    |      |
       |--------------<----------|
  t1.L ↑    t1.R        t2.R     ↑ t2.R
       |                         | 
      KInL                      KInR
"""
    gd4 = Spectral.kernelFermiT(t1[2] - t2[2], ϵ2, β)
    G = gt1 * gd4 / (2π)^3 * phaseT(t1[1], t1[2], t2[2], t2[1], isT)
    we += G * (wle * wre)
    # println(G * (wle * wre) * (2π)^3)
    ##################################################

    # println(weight)
    # return Weight(wd, we)
    # println("weight ", wd)
    # exit(0)
    if isT
        return Weight(wd, we)
    else
        return Weight(-we, -wd)
    end
end

function measure(config)
    factor = 1.0 / config.reweight[config.curr]
    x = config.var[3][1]
    chan = config.var[4][1]
    # println(config.observable[1][1])
    if config.curr == 1
        weight = integrand(config)
        config.observable[1, chan, 1] += weight.d / 2 / abs(weight) * factor
        config.observable[1, chan, 2] += weight.e / 2 / abs(weight) * factor
        config.observable[2, chan, 1] += weight.d * x / 2 / abs(weight) * factor
        config.observable[2, chan, 2] += weight.e * x / 2 / abs(weight) * factor
    else
        return
    end
end

function MC()
    para = Para(1)

    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF)
    T = MCIntegration.Tau(β, β / 2.0)
    X = MCIntegration.Continuous([-1.0, 1.0], 0.2) #x=cos(θ)
    # ExtKidx = MCIntegration.Discrete(1, Nk)
    Chan = MCIntegration.Discrete(1, 3)

    # for (ti, t) in enumerate(T.data)
    #     t[1] = β * rand()
    #     t[2] = β * rand()
    # end

    dof = [[1, 4, 1, 1],] # K, T, ExtKidx
    obs = zeros(Nk, 3, 2) # observable for the Fock diagram 

    config = MCIntegration.Configuration(steps, (K, T, X, Chan), dof, obs; para = para)
    avg, std = MCIntegration.sample(config, integrand, measure; print = 0, Nblock = 16, reweight = 1e5)

    function info(idx, chan, di)
        if chan == 4 #print the sum of all channels
            return @sprintf("   %8.4f ±%8.4f", sum(avg[idx, :, di]), sum(std[idx, :, di]))
        else
            return @sprintf("   %8.4f ±%8.4f", avg[idx, chan, di], std[idx, chan, di])
        end
    end

    if isnothing(avg) == false
        avg *= NF
        std *= NF
        println("Direct ver4: ")
        for (ki, l) in enumerate(lgrid)
            println(@sprintf("%8.4f", l) * info(ki, 1, 1) * info(ki, 2, 1) * info(ki, 3, 1) * info(ki, 4, 1))
        end
        println("Exchange ver4: ")
        for (ki, l) in enumerate(lgrid)
            println(@sprintf("%8.4f", l) * info(ki, 1, 2) * info(ki, 2, 2) * info(ki, 3, 2) * info(ki, 4, 2))
        end

        Fp = Fm = cp = cm = 0.0
        θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 16)
        qs = [2 * kF * sin(θ / 2) for θ in θgrid.grid]
        Wp, Wm = KOstatic(Fp, Fm, cp, cm, 1.0, qs)

        Wp .*= NF
        Wm .*= NF

        Ws = -(Wp + 3 * Wm) / 2
        Wa = -(Wp - Wm) / 2
        # Ws .+= Fp - (cp - 3 * cm) / 2
        # Wa .+= Fm + (cp - 3 * cm) / 2
        # println(Ws)
        Ws0 = Interp.integrate1D(Ws .* sin.(θgrid.grid), θgrid) / 2
        Wa0 = Interp.integrate1D(Wa .* sin.(θgrid.grid), θgrid) / 2
        println("l=0:")
        # println("before flip: ", Ws0, Wa0)
        println("F0+=", -Ws0)
        println("F0-=", -Wa0)

        println("l=0 correction: ")
        dfd, dfd_err = avg[1, 1, 1], std[1, 1, 1]
        dfe, dfe_err = avg[1, 1, 2], std[1, 1, 2]
        println("dF0+=", dfd + dfe / 2, "  ± ", dfd_err + dfe_err / 2)
        println("dF0-=", dfe / 2, "  ± ", dfe_err / 2)

    end

end

MC()