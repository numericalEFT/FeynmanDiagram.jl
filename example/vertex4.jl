# This example demonstrated how to calculate one-loop diagram of free electrons using the Monte Carlo module
# Observable is normalized: Γ₄*N_F where N_F is the free electron density of states

using LinearAlgebra, Random, Printf, BenchmarkTools, InteractiveUtils, Parameters
using ElectronGas
using CompositeGrids
using MCIntegration
using Lehmann
include("parameter.jl")
include("interaction.jl")

const steps = 1e6 # MC steps of each worker

# claim all globals to be constant, otherwise, global variables could impact the efficiency
########################### parameters ##################################
const IsF = false # calculate quasiparticle interaction F or not
const AngSize = 16

########################## variables for MC integration ##################
const KInL = [kF, 0.0, 0.0] # incoming momentum of the left particle
const Qd = [0.0, 0.0, 0.0] # transfer momentum is zero in the forward scattering channel

struct Para{Q,T}
    extAngle::Vector{Float64}
    dW0::Matrix{Float64}
    qgrid::Q
    τgrid::T
    function Para(AngSize)
        extAngle = collect(LinRange(0.0, π, AngSize)) # external angle grid
        qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 6 * kF], [0.0, 2kF], 16, 0.01 * kF, 8)
        τgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 16, β * 1e-4, 8)

        vqinv = [(q^2 + mass2) / (4π * e0^2) for q in qgrid.grid]
        dW0 = TwoPoint.dWRPA(vqinv, qgrid.grid, τgrid.grid, dim, EF, kF, β, spin, me) # dynamic part of the effective interaction
        return new{typeof(qgrid),typeof(τgrid)}(extAngle, dW0, qgrid, τgrid)
    end
end

function phase(tInL, tOutL, tInR, tOutR)
    if (IsF)
        return cos(π * ((tInL + tOutL) - (tInR + tOutR)) / β)
    else
        return cos(π * ((tInL - tOutL) + (tInR - tOutR)) / β)
    end
end

function integrand(config)
    if config.curr == 1
        return eval_T(config)
    else
        error("Not implemented!")
    end
end

function eval_T(config)
    para = config.para
    T, K, Ang = config.var[1], config.var[2], config.var[3]
    k1, k2 = K[1], K[1] - Qd
    t1, t2 = T[1], T[2] # t1, t2 both have two tau variables
    θ = para.extAngle[Ang[1]] # angle of the external momentum on the right
    KInR = [kF * cos(θ), kF * sin(θ), 0.0]

    vld, wld, vle, wle = vertexDynamic(para, Qd, KInL - k1, t1[1], t1[2])
    vrd, wrd, vre, wre = vertexDynamic(para, Qd, KInR - k2, t2[1], t2[2])
    # wld, wle, wrd, wre = 0.0, 0.0, 0.0, 0.0

    wd, we = 0.0, 0.0

    # possible green's functions on the top
    ϵ1, ϵ2 = (dot(k1, k1) - kF^2) / (2me), (dot(k2, k2) - kF^2) / (2me)
    gt1 = Spectral.kernelFermiT(t2[1] - t1[1], ϵ1, β)


    gt2 = Spectral.kernelFermiT(t1[1] - t2[1], ϵ2, β)
    # wd += 1.0 / β * 1.0 / β * gt1 * gt2 / (2π)^3 * phase(t1[1], t1[1], t2[1], t2[1])

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
    G = gt1 * gd1 / (2π)^3 * phase(t1[1], t1[1], t2[1], t2[1])
    we += G * (vle * vre)
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
    G = gt1 * gd2 / (2π)^3 * phase(t1[1], t1[2], t2[1], t2[1])
    we += G * (wle * vre)
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
    G = gt1 * gd3 / (2π)^3 * phase(t1[1], t1[1], t2[2], t2[1])
    we += G * (vle * wre)
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
    G = gt1 * gd4 / (2π)^3 * phase(t1[1], t1[2], t2[2], t2[1])
    we += G * (wle * wre)
    ##################################################

    # println(weight)
    return Weight(wd, we)
end

function measure(config)
    angidx = config.var[3][1]
    factor = 1.0 / config.reweight[config.curr]
    if config.curr == 1
        weight = integrand(config)
        # println(weight)
        config.observable[angidx, 1] += weight.d / abs(weight) * factor
        config.observable[angidx, 2] += weight.e / abs(weight) * factor
    else
        error("Not implemented!")
    end
end

function run()
    T = MCIntegration.TauPair(β, β / 2.0)
    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF)
    Ext = MCIntegration.Discrete(1, AngSize) # external variable is specified

    dof = [[2, 1, 1],]
    obs = zeros(Float64, (AngSize, 2))

    para = Para(AngSize)

    config = MCIntegration.Configuration(steps, (T, K, Ext), dof, obs; para = para)
    avg, std = MCIntegration.sample(config, integrand, measure; print = 0, Nblock = 16)

    if isnothing(avg) == false
        println("NF = $NF")
        avg *= NF
        std *= NF
        println("Direct ver4: ")
        for (ki, theta) in enumerate(para.extAngle)
            @printf("%10.6f   %10.6f ± %10.6f\n", theta, avg[ki, 1], std[ki, 1])
        end
        println("Exchange ver4: ")
        for (ki, theta) in enumerate(para.extAngle)
            @printf("%10.6f   %10.6f ± %10.6f\n", theta, avg[ki, 2], std[ki, 2])
        end
    end
end

# @btime run(1, 10)
run()
# @time run(Repeat, totalStep)