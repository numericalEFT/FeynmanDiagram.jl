using Printf, LinearAlgebra
using CompositeGrids
using ElectronGas
using Parameters, Random
using MCIntegration
using Lehmann

using ExpressionTree
using StaticArrays

Parquet = GWKT.Parquet

include("parameter.jl")
include("interaction.jl")


const steps = 1e5
const isF = true
const Nk = 16
const θgrid = collect(LinRange(0.1, π, Nk)) # external angle grid
const ExtK = [@SVector [kF * cos(θ), kF * sin(θ), 0.0] for θ in θgrid]
const RefK = @SVector [kF, 0.0, 0.0]

struct Para{Q,T}
    dW0::Matrix{Float64}
    qgrid::Q
    τgrid::T # dedicated τgrid for dynamic interaction
    ver4::Parquet.Ver4
    function Para(loopOrder::Int)
        qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 6 * kF], [0.0, 2kF], 16, 0.01 * kF, 8)
        τgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 16, β * 1e-4, 8)

        vqinv = [(q^2 + mass2) / (4π * e0^2) for q in qgrid.grid]
        dW0 = TwoPoint.dWRPA(vqinv, qgrid.grid, τgrid.grid, dim, EF, kF, β, spin, me) # dynamic part of the effective interaction


        chan = [Parquet.T, Parquet.U]
        # chan = [Parquet.S]
        para = Parquet.Para(chan, [1, 2])
        ver4 = Parquet.Ver4{Weight}(loopOrder, 1, para)

        return new{typeof(qgrid),typeof(τgrid)}(dW0, qgrid, τgrid, ver4)
    end
end

function integrand(config)
    if config.curr == 1
        extKidx = config.var[3][1]
        KinL, KoutL, KinR, KoutR = RefK, RefK, ExtK[extKidx], ExtK[extKidx]
        eval(config, config.para.ver4, KinL, KoutL, KinR, KoutR, 1, true)
        ver4 = config.para.ver4
        w = ver4.weight
        wd, we = 0.0, 0.0
        for c in [Parquet.T, Parquet.U, Parquet.S]
            #     #     # println(c, ", ", Parquet.SymFactor[c])
            wd += w[c].d * Parquet.SymFactor[c]
            we += w[c].e * Parquet.SymFactor[c]
        end
        # wd = w[Parquet.S].d * Parquet.SymFactor[Parquet.S]
        # we = w[Parquet.S].e * Parquet.SymFactor[Parquet.S]
        # wd = 0.0
        # @assert Parquet.SymFactor[Parquet.S] ≈ -0.5
        # wd = w[Parquet.S].d * (0.5)
        # we = w[Parquet.S].e * (0.5)
        # @assert w[Parquet.S].d * Parquet.SymFactor[Parquet.S] ≈ wd
        # @assert w[Parquet.S].e * Parquet.SymFactor[Parquet.S] ≈ we
        # @assert Parquet.S
        # exit()
        return Weight(wd, we)
        # return Weight(0.0, config.para.ver4.weight[4].e * Parquet.SymFactor[4])
        # return Weight(0.0, config.para.ver4.weight[4].e * (-0.5))
    else
        error("impossible!")
    end
end

function measure(config)
    factor = 1.0 / config.reweight[config.curr]
    extKidx = config.var[3][1]
    # println(config.observable[1][1])
    if config.curr == 1
        weight = integrand(config)
        # println(weight.d, ", ", weight.e)
        # println(extKidx)
        config.observable[extKidx, 1] += weight.d / abs(weight) * factor
        config.observable[extKidx, 2] += weight.e / abs(weight) * factor
        # println(config.observable[:, 1])

    else
        return
    end
end

function evalG(G, K, varT)
    ϵ = (dot(K, K) - kF^2) / (2me)
    for (i, tpair) in enumerate(G.Tpair)
        G.weight[i] = Spectral.kernelFermiT(varT[tpair[2]] - varT[tpair[1]], ϵ, β)
    end
end

function phase(varT, Tpair)

    tInL, tOutL, tInR, tOutR = varT[Tpair[INL]], varT[Tpair[OUTL]], varT[Tpair[INR]],
    varT[Tpair[OUTR]]
    if (isF)
        return cos(π / β * ((tInL + tOutL) - (tInR + tOutR)))
    else
        return cos(π / β * ((tInL - tOutL) + (tInR - tOutR)))
    end
end

function eval(config, ver4, KinL, KoutL, KinR, KoutR, Kidx::Int, fast = false)
    para = config.para
    varK, varT = config.var[1], config.var[2]

    if ver4.loopNum == 0
        qd = KinL - KoutL
        qe = KinL - KoutR
        τIn, τOut = varT[ver4.Tidx], varT[ver4.Tidx+1]
        vd, wd, ve, we = vertexDynamic(para, qd, qe, τIn, τOut)

        # ver4.weight[1].d = 0.0
        # ver4.weight[1].e = 0.0
        # ver4.weight[2].d = 0.0
        # ver4.weight[2].e = 0.0
        # ver4.weight[3].d = 0.0
        # ver4.weight[3].e = 0.0

        ver4.weight[1].d = vd
        ver4.weight[1].e = ve
        ver4.weight[2].d = wd
        ver4.weight[2].e = 0.0
        ver4.weight[3].d = 0.0
        ver4.weight[3].e = we
        return
    end

    # LoopNum>=1
    for w in ver4.weight
        w.d, w.e = 0.0, 0.0 # initialize all weights
    end
    G = ver4.G
    K = varK[Kidx]
    evalG(G[1], K, varT)
    PhaseFactor = 1.0 / (2π)^dim
    Kt, Ku, Ks = similar(K), similar(K), similar(K) #Kt, Ku and Ks will be re-created later, slow in performance

    for c in ver4.chan
        if c == Parquet.T
            @. Kt = KoutL + K - KinL
            evalG(G[c], Kt, varT)
        elseif c == Parquet.U
            # can not be in box!
            @. Ku = KoutR + K - KinL
            evalG(G[c], Ku, varT)
        elseif c == Parquet.S
            # S channel, and cann't be in box!
            @. Ks = KinL + KinR - K
            evalG(G[c], Ks, varT)
        else
            error("not impossible!")
        end
    end
    for b in ver4.bubble
        c = b.chan
        Factor = Parquet.SymFactor[c] * PhaseFactor
        Llopidx = Kidx + 1
        Rlopidx = Kidx + 1 + b.Lver.loopNum

        if c == Parquet.T
            eval(config, b.Lver, KinL, KoutL, Kt, K, Llopidx)
            eval(config, b.Rver, K, Kt, KinR, KoutR, Rlopidx)
        elseif c == Parquet.U
            eval(config, b.Lver, KinL, KoutR, Ku, K, Llopidx)
            eval(config, b.Rver, K, Ku, KinR, KoutL, Rlopidx)
        elseif c == Parquet.S
            # S channel
            eval(config, b.Lver, KinL, Ks, KinR, K, Llopidx)
            eval(config, b.Rver, K, KoutL, Ks, KoutR, Rlopidx)
        else
            error("not implemented")
        end

        rN = length(b.Rver.weight)
        gWeight = 0.0
        for (l, Lw) in enumerate(b.Lver.weight)
            for (r, Rw) in enumerate(b.Rver.weight)
                map = b.map[(l-1)*rN+r]

                gWeight = G[1].weight[map.G0] * G[c].weight[map.Gx] * Factor

                if fast && ver4.level == 1
                    gWeight *= phase(varT, ver4.Tpair[map.ver])
                    w = ver4.weight[c]
                else
                    w = ver4.weight[map.ver]
                end

                if c == Parquet.T
                    # w.d += gWeight * (Lw.d * Rw.d * spin + Lw.d * Rw.e + Lw.e * Rw.d)
                    # w.d += gWeight * (Lw.d * Rw.e + Lw.e * Rw.d)
                    w.e += gWeight * Lw.e * Rw.e
                elseif c == Parquet.U
                    w.d += gWeight * Lw.e * Rw.e
                    # w.e += gWeight * (Lw.d * Rw.d * spin + Lw.d * Rw.e + Lw.e * Rw.d)
                    w.e += gWeight * (+Lw.d * Rw.e + Lw.e * Rw.d)
                elseif c == Parquet.S
                    # S channel,  see the note "code convention"
                    w.d += gWeight * (Lw.d * Rw.e + Lw.e * Rw.d)
                    w.e += gWeight * (Lw.d * Rw.d + Lw.e * Rw.e)
                else
                    error("not implemented")
                end

            end
        end

    end
end

function MC()
    para = Para(1)

    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF)
    T = MCIntegration.Tau(β, β / 2.0)
    ExtKidx = MCIntegration.Discrete(1, Nk)

    # for (ti, t) in enumerate(T.data)
    #     t[1] = β * rand()
    #     t[2] = β * rand()
    # end

    dof = [[1, 4, 1],] # K, T, ExtKidx
    obs = zeros(Nk, 2) # observable for the Fock diagram 

    config = MCIntegration.Configuration(steps, (K, T, ExtKidx), dof, obs; para = para)
    avg, std = MCIntegration.sample(config, integrand, measure; print = 10, Nblock = 16)
    if isnothing(avg) == false
        avg *= NF
        std *= NF
        println("Direct ver4: ")
        for (ki, theta) in enumerate(θgrid)
            @printf("%10.6f   %10.6f ± %10.6f\n", theta, avg[ki, 1], std[ki, 1])
        end
        println("Exchange ver4: ")
        for (ki, theta) in enumerate(θgrid)
            @printf("%10.6f   %10.6f ± %10.6f\n", theta, avg[ki, 2], std[ki, 2])
        end
    end

end

MC()