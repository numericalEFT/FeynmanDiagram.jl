using Printf, LinearAlgebra
using CompositeGrids
using ElectronGas
using Parameters, Random
using MCIntegration
using Lehmann

using ExpressionTree
using AbstractTrees
# using NewickTree
using StaticArrays

include("parameter.jl")
include("interaction.jl")


const steps = 1e6
const isF = false
const Nk = 16
const θgrid = [k * 2π / Nk for k = 0:Nk-1]
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


        chan = [Parquet.T, Parquet.U, Parquet.S]
        para = Parquet.Para(chan, [1, 2])
        ver4 = Parquet.Ver4{Weight}(loopOrder, 1, para)

        return new{typeof(qgrid),typeof(τgrid)}(dW0, qgrid, τgrid, ver4)
    end
end

function integrand(config)
    if config.curr == 1
        extKidx = config.var[3][1]
        KinL, KoutL, KinR, KoutR = RefK, RefK, ExtK[extKidx], ExtK[extKidx]
        eval(config, KinL, KoutL, KinR, KoutR, 1, true)
        w = config.para.ver4.weight
        return w[1] .+ w[2] .+ w[3] .+ w[4]
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
        config.observable[extKidx] += weight / abs(weight) * factor
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

function eval(config, KinL, KoutL, KinR, KoutR, Kidx::Int, fast = false)
    para = config.para
    ver4 = para.ver4
    varK, varT = config.var[1], config.var[2]

    if ver4.loopNum == 0
        qd = KinL - KoutL
        qe = KinL - KoutR
        τIn, τOut = varT[ver4.Tidx], varT[ver4.Tidx+1]
        vd, wd, ve, we = vertexDynamic(para, qd, qe, τIn, τOut)
        ver4.weight[1][DI] = vd
        ver4.weight[1][EX] = ve
        ver4.weight[2][DI] = wd
        ver4.weight[3][EX] = we
        return
    end

    # LoopNum>=1
    for w in ver4.weight
        w *= 0.0 # initialize all weights
    end
    G = ver4.G
    K = varK[Kidx]
    evalG(G[1], K, varT)
    PhaseFactor = 1.0 / (2π)^dim

    for c in ver4.chan
        if c == T
            Kt .= KoutL .+ K .- KinL
            evalG(G[T], Kt, varT)
        elseif c == U
            # can not be in box!
            Ku .= KoutR .+ K .- KinL
            evalG(G[U], Ku, varT)
        else
            # S channel, and cann't be in box!
            Ks .= KinL .+ KinR .- K
            evalG(G[S], Ks, varT)
        end
    end
    for b in ver4.bubble
        c = b.chan
        Factor = SymFactor[c] * PhaseFactor
        Llopidx = Kidx + 1
        Rlopidx = Kidx + 1 + b.Lver.loopNum

        if c == T
            eval(config, b.Lver, KinL, KoutL, Kt, K, Llopidx)
            eval(config, b.Rver, K, Kt, KinR, KoutR, Rlopidx)
        elseif c == U
            eval(config, b.Lver, KinL, KoutR, Ku, K, Llopidx)
            eval(config, b.Rver, K, Ku, KinR, KoutL, Rlopidx)
        else
            # S channel
            eval(config, b.Lver, KinL, Ks, KinR, K, Llopidx)
            eval(config, b.Rver, K, KoutL, Ks, KoutR, Rlopidx)
        end

        rN = length(b.Rver.weight)
        gWeight = 0.0
        for (l, Lw) in enumerate(b.Lver.weight)
            for (r, Rw) in enumerate(b.Rver.weight)
                map = b.map[(l-1)*rN+r]

                gWeight = G[1].weight[map.G] * G[c].weight[map.Gx] * Factor

                if fast && ver4.level == 1
                    pair = ver4.Tpair[map.ver]
                    if (isF)
                        dT =
                            varT[pair[INL]] + varT[pair[OUTL]] - varT[pair[INR]] -
                            varT[pair[OUTR]]
                    else
                        dT =
                            varT[pair[INL]] - varT[pair[OUTL]] + varT[pair[INR]] -
                            varT[pair[OUTR]]
                    end
                    gWeight *= cos(2.0 * pi / Beta * dT)
                    w = ver4.weight[c]
                else
                    w = ver4.weight[map.ver]
                end

                if c == T
                    w[DI] +=
                        gWeight *
                        (Lw[DI] * Rw[DI] * SPIN + Lw[DI] * Rw[EX] + Lw[EX] * Rw[DI])
                    w[EX] += gWeight * Lw[EX] * Rw[EX]
                elseif c == U
                    w[DI] += gWeight * Lw[EX] * Rw[EX]
                    w[EX] +=
                        gWeight *
                        (Lw[DI] * Rw[DI] * SPIN + Lw[DI] * Rw[EX] + Lw[EX] * Rw[DI])
                else
                    # S channel,  see the note "code convention"
                    w[DI] += gWeight * (Lw[DI] * Rw[EX] + Lw[EX] * Rw[DI])
                    w[EX] += gWeight * (Lw[DI] * Rw[DI] + Lw[EX] * Rw[EX])
                end

            end
        end

    end
end

function MC()
    para = Para(1)

    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF)
    T = MCIntegration.Tau(β, β / 2.0)
    ExtKidx = MCIntegration.Discrete(0, Nk - 1)

    dof = [[1, 4, 1],] # K, T, ExtKidx
    obs = zeros(2) # observable for the Fock diagram 

    config = MCIntegration.Configuration(steps, (K, T, ExtKidx), dof, obs; para = para)
    avg, std = MCIntegration.sample(config, integrand, measure; print = 0, Nblock = 16)
    if isnothing(avg) == false
        for (ki, theta) in θgrid
            @printf("%10.6f   %10.6f ± %10.6f\n", theta, avg[ki] / NF, std[ki] / NF)
        end
    end

end

MC()