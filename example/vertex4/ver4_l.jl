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

include("interaction.jl")
include("ver4_diag.jl")


# println(dW0)
# exit(0)
const lgrid = [1, 2]
const Nl = length(lgrid)

function integrand(config)
    order = config.curr
    x = config.var[3][1]
    varK, varT = config.var[1], config.var[2]

    varK.data[:, 2] = [kF * x, kF * sqrt(1 - x^2), 0.0]

    ExprTree.evalNaive!(diag[order], varK.data, varT.data, eval)
    if !isempty(rootuu[order])
        wuu = sum(diag[order].node.current[root] * phase(varT, extTuu[order][ri]) for (ri, root) in enumerate(rootuu[order]))
    else
        wuu = 0.0
    end
    if !isempty(rootud[order])
        wud = sum(diag[order].node.current[root] * phase(varT, extTud[order][ri]) for (ri, root) in enumerate(rootud[order]))
    else
        wud = 0.0
    end
    # println(wuu, ",  ", wud)
    return Weight(wuu / β, wud / β)
end

function measure(config)
    factor = 1.0 / config.reweight[config.curr]
    x = config.var[3][1]
    # println(config.observable[1][1])
    if config.curr == 1
        weight = integrand(config)
        config.observable[1, 1] += weight.d / 2 / abs(weight) * factor
        config.observable[1, 2] += weight.e / 2 / abs(weight) * factor
        config.observable[2, 1] += weight.d * x / 2 / abs(weight) * factor
        config.observable[2, 2] += weight.e * x / 2 / abs(weight) * factor
    else
        return
    end
end

function MC()
    K = MCIntegration.FermiK(para.dim, kF, 0.2 * kF, 10.0 * kF, offset=2)
    K.data[:, 1] .= [kF, 0.0, 0.0]
    T = MCIntegration.Tau(β, β / 2.0)
    X = MCIntegration.Continuous([-1.0, 1.0], 0.2) #x=cos(θ)

    dof = [[diagpara[o].innerLoopNum, diagpara[o].totalTauNum, 1] for o in 1:Order] # K, T, ExtKidx
    obs = zeros(Nl, 2) # observable for the Fock diagram 

    config = MCIntegration.Configuration(steps, (K, T, X), dof, obs)
    avg, std = MCIntegration.sample(config, integrand, measure; print=0, Nblock=16)

    function info(idx, di)
        return @sprintf("   %8.4f ±%8.4f", avg[idx, di], std[idx, di])
    end

    if isnothing(avg) == false
        avg *= NF
        std *= NF
        N = size(avg)[1]
        grid = lgrid

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
            @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], (avg[li, 1] - avg[li, 2]) / 2, (std[li, 1] + std[li, 2]) / 2)
        end


    end

end

MC()