using Printf, LinearAlgebra
using CompositeGrids
using ElectronGas
using Parameters, Random
using MCIntegration
using Lehmann

const steps = 1e8

include("parameter.jl")
include("interaction.jl")

kgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 3kF], [0.0, kF], 8, 0.01 * kF, 8) # external K grid for sigma
dlr = DLR.DLRGrid(:fermi, 10EF, β, 1e-10)

struct Para{Q,T}
    dW0::Matrix{Float64}
    qgrid::Q
    τgrid::T # dedicated τgrid for dynamic interaction
    function Para()
        qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 6 * kF], [0.0, 2kF], 16, 0.01 * kF, 8)
        τgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 16, β * 1e-4, 8)

        vqinv = [(q^2 + mass2) / (4π * e0^2) for q in qgrid.grid]
        dW0 = TwoPoint.dWRPA(vqinv, qgrid.grid, τgrid.grid, dim, EF, kF, β, spin, me) # dynamic part of the effective interaction
        return new{typeof(qgrid),typeof(τgrid)}(dW0, qgrid, τgrid)
    end
end

function integrand(config)
    if config.curr == 1
        return eval2(config)
    else
        error("impossible!")
    end
end

function eval2(config)
    para = config.para

    K, T = config.var[1], config.var[2]
    k, τ = K[1], T[1]
    k0 = [0.0, 0.0, kF] # external momentum
    kq = k + k0
    ω = (dot(kq, kq) - kF^2) / (2me)
    g = Spectral.kernelFermiT(τ, ω, β)
    v, dW = interactionDynamic(para, k, 0.0, τ)
    phase = 1.0 / (2π)^3
    return g * dW * phase
end

function measure(config)
    factor = 1.0 / config.reweight[config.curr]
    τ = config.var[2][1]
    # println(config.observable[1][1])
    if config.curr == 1
        weight = integrand(config)
        config.observable[1] += weight * sin(-π / β * τ) / abs(weight) * factor
        config.observable[2] += weight * sin(π / β * τ) / abs(weight) * factor
    else
        return
    end
end

function fock(extn)
    para = Para()
    Ksize = length(kgrid.grid)

    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF)
    T = MCIntegration.Tau(β, β / 2.0)
    # Ext =
    dof = [[1, 1],] # degrees of freedom of the Fock diagram
    obs = zeros(2) # observable for the Fock diagram 

    config = MCIntegration.Configuration(steps, (K, T), dof, obs; para = para)
    avg, std = MCIntegration.sample(config, integrand, measure; print = 0, Nblock = 16)

    if isnothing(avg) == false
        @printf("%10.6f   %10.6f ± %10.6f\n", -1.0, avg[1], std[1])
        @printf("%10.6f   %10.6f ± %10.6f\n", 0.0, avg[2], std[2])

        dS_dw = (avg[1] - avg[2]) / (2π / β)
        error = (std[1] + std[2]) / (2π / β)
        println("dΣ/diω= $dS_dw ± $error")
        Z = (1 / (1 + dS_dw))
        Zerror = error / Z^2
        println("Z=  $Z ± $Zerror")
        # TODO: add errorbar estimation
        # println
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    # using Gaston
    ngrid = [-1, 0, 1]
    fock(ngrid)
end