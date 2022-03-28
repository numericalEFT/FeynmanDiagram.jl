using Printf, LinearAlgebra
using CompositeGrids
using ElectronGas
using Parameters, Random, DataFrames
using MCIntegration
using Lehmann

using FeynmanDiagram
using StaticArrays

const steps = 1e6
const isF = true

include("parameter.jl")
include("interaction.jl")
include("diagram.jl")


# println(dW0)
# exit(0)
const Nk = 16
const θgrid = collect(LinRange(0.1, π, Nk)) # external angle grid
const ExtK = [[kF * cos(θ), kF * sin(θ), 0.0] for θ in θgrid]

function measure(config)
    factor = 1.0 / config.reweight[config.curr]
    x = config.var[3][1]
    # println(config.observable[1][1])
    if config.curr == 1
        weight = integrand(config)
        config.observable[x, 1] += weight.d / abs(weight) * factor
        config.observable[x, 2] += weight.e / abs(weight) * factor
    else
        return
    end
end

function MC()
    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=2)
    K.data[:, 1] .= RefK
    T = MCIntegration.Tau(β, β / 2.0)
    # X = MCIntegration.Continuous([-1.0, 1.0], 0.2) #x=cos(θ)
    X = MCIntegration.Discrete(1, Nk)

    # for (ti, t) in enumerate(T.data)
    #     t[1] = β * rand()
    #     t[2] = β * rand()
    # end

    dof = [[para[o].innerLoopNum, para[o].totalTauNum, 1] for o in 1:Order] # K, T, ExtKidx
    obs = zeros(Nk, 2) # observable for the Fock diagram 

    config = MCIntegration.Configuration(steps, (K, T, X), dof, obs)
    avg, std = MCIntegration.sample(config, integrand, measure; print=0, Nblock=16)

    function info(idx, di)
        return @sprintf("   %8.4f ±%8.4f", avg[idx, di], std[idx, di])
    end

    if isnothing(avg) == false
        avg *= NF
        std *= NF
        N = size(avg)[1]
        grid = θgrid
        println("UpUp ver4: ")
        for li in 1:N
            @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], avg[li, 1], std[li, 1])
        end
        println("UpDown ver4: ")
        for li in 1:N
            @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], avg[li, 2], std[li, 2])
        end

        println("S ver4: ")
        for li in 1:N
            @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], (avg[li, 1] + avg[li, 2]) / 2, (std[li, 1] + std[li, 2]) / 2)
        end
        println("A ver4: ")
        for li in 1:N
            @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], (avg[li, 1] - avg[li, 2]) / 2, (std[li, 1] - std[li, 2]) / 2)
        end

    end

end

MC()