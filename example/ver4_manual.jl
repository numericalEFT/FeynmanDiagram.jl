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


const steps = 1e5
const isF = false
const irreducible = true
const Nk = 16
const θgrid = collect(LinRange(0.1, π, Nk)) # external angle grid
const ExtK = [@SVector [kF * cos(θ), kF * sin(θ), 0.0] for θ in θgrid]
const RefK = @SVector [kF, 0.0, 0.0]

struct Para{Q,T}
    dW0::Matrix{Float64}
    qgrid::Q
    τgrid::T # dedicated τgrid for dynamic interaction
    diag::DiagTree.Diagrams{Float64}
    dir::Int
    ex::Int
    function Para(loopOrder::Int)
        qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 6 * kF], [0.0, 2kF], 16, 0.01 * kF, 8)
        τgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 16, β * 1e-4, 8)

        vqinv = [(q^2 + mass2) / (4π * e0^2) for q in qgrid.grid]
        dW0 = TwoPoint.dWRPA(vqinv, qgrid.grid, τgrid.grid, dim, EF, kF, β, spin, me) # dynamic part of the effective interaction

        chan = [1,]
        KinL = KoutL = [1, 0, 0]
        KinR = KoutR = [0, 1, 0]
        legK = [KinL, KoutL, KinR, KoutR]
        Gsym = [:mirror]
        Wsym = [:mirror, :timereversal]
        diag, dir, ex = Manual.build(chan, legK, 3, spin, irreducible, Gsym, Wsym)
        DiagTree.showTree(diag, ex)

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
            return v * factor
        elseif type == 3 #W
            return w * factor
        else
            error("not implemented!")
        end
    end
end

function integrand(config)
    if config.curr == 1
        extKidx = config.var[3][1]
        K, varT = config.var[1], config.var[2]
        varK = [RefK, ExtK[extKidx], K[1]]
        para = config.para
        wd, we = DiagTree.evalNaive(para.diag, evalPropagator, varK, varT, [para.dir, para.ex], phase, para)
        factor = 1 / (2π)^dim
        return Weight(wd * factor, we * factor)
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
        config.observable[extKidx, 1] += weight.d / abs(weight) * factor
        config.observable[extKidx, 2] += weight.e / abs(weight) * factor
    else
        return
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