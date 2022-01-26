# This example demonstrated how to calculate the bubble diagram of free electrons using the Monte Carlo module

using LinearAlgebra, Random, Printf, BenchmarkTools, InteractiveUtils, Parameters, StaticArrays
using MCIntegration
using FeynmanDiagram
using ElectronGas
using Lehmann

const steps = 1e7 # MC steps of each worker
const Order = 3
const rs = 1.0
const λ = 1.0
const beta = 25.0
const basic = reconstruct(Parameter.rydbergUnit(1 / beta, rs, 3), Λs = λ)
const β = basic.β
const kF = basic.kF
const me = basic.me
const spin = basic.spin

###################  build polarization diagram #######################
function getDiagPara(order)
    return GenericPara(diagType = PolarDiag,
        innerLoopNum = order,
        hasTau = true,
        loopDim = basic.dim,
        spin = basic.spin,
        interaction = [FeynmanDiagram.Interaction(ChargeCharge, Instant),],
        filter = [NoFock,]
    )
end

const diag_para = [getDiagPara(o) for o in 1:Order]
const polar = [mergeby(Parquet.polarization(diag_para[i])).diagram[1] for i in 1:Order]
# println(polar)
# plot_tree(polar)
const diag = [ExprTree.build(polar[o], 1)[1] for o in 1:Order] #different order has different set of K, T variables, thus must have different exprtrees
# println(diag)
# exit(0)

##################### parameters for MC ##############################
@with_kw struct ParaMC
    n::Int = 0 # external Matsubara frequency
    Qsize::Int = 8
    varK::Matrix = zeros(basic.dim, diag_para[end].totalLoopNum) # make sure that varK is large enough for all orders
    extQ::Vector{SVector{basic.dim,Float64}} = [@SVector [q, 0.0, 0.0] for q in LinRange(0.0, 3.0 * basic.kF, Qsize)]
end

##################### propagator and interaction evaluation ##############
function eval(id::GreenId, K, varT)
    τin, τout = varT[id.extT[1]], varT[id.extT[2]]
    # println(τBasis, ", ", varT)
    ϵ = dot(K, K) / (2me) - basic.μ
    τ = τout - τin
    if τ ≈ 0.0
        return Spectral.kernelFermiT(-1e-8, ϵ, basic.β)
    else
        return Spectral.kernelFermiT(τ, ϵ, basic.β)
    end
end

eval(id::InteractionId, K, varT) = (basic.e0)^2 / basic.ϵ0 / (dot(K, K) + basic.Λs)

################### interface to MC #########################################
function integrand(config)
    # if config.curr == 0
    #     error("impossible")
    # end
    @assert config.curr >= 0
    order = config.curr
    para = config.para

    T, K, Ext = config.var[1], config.var[2], config.var[3]
    extidx = Ext[1]

    para.varK[:, 1] .= para.extQ[extidx] # external momentum
    for i in 2:(diag_para[order].innerLoopNum+1)
        para.varK[:, i] .= K[i-1]
    end

    weight = ExprTree.evalNaive!(diag[order], para.varK, T, eval)
    w1 = weight[1] * cos(2π * para.n * (T[2] - T[1]) / β) / β
    # @assert w0 ≈ w1 "$w0 vesus $w1"

    # return weight[1] * cos(2π * para.n * (T[2] - T[1]) / β) / β
    return w1
end

function measure(config)
    obs = config.observable
    factor = 1.0 / config.reweight[config.curr]
    extidx = config.var[3][1]
    weight = integrand(config)
    obs[config.curr, extidx] += weight / abs(weight) * factor
end

function run(steps)

    para = ParaMC()
    @unpack extQ, Qsize = para
    @unpack kF, β = basic

    T = Tau(β, β / 2.0)
    K = FermiK(basic.dim, kF, 0.2 * kF, 10.0 * kF)

    Ext = MCIntegration.Discrete(1, length(extQ)) # external variable is specified

    dof = [[diag_para[o].totalTauNum, diag_para[o].innerLoopNum, 1] for o in 1:Order] # degrees of freedom of the normalization diagram and the bubble
    obs = zeros(Float64, (Order, Qsize)) # observable for the normalization diagram and the bubble

    config = Configuration(steps, (T, K, Ext), dof, obs; para = para)
    avg, std = MCIntegration.sample(config, integrand, measure; print = 0, Nblock = 16, reweight = 10000)
    # @profview MonteCarlo.sample(config, integrand, measure; print=0, Nblock=1)
    # sleep(100)

    if isnothing(avg) == false
        @unpack n, extQ = ParaMC()

        println("Order 1")
        for (idx, q) in enumerate(extQ)
            q = q[1]
            p = Polarization.Polarization0_ZeroTemp(q, para.n, basic)
            @printf("%10.6f  %10.6f ± %10.6f  %10.6f\n", q / basic.kF, avg[1, idx], std[1, idx], p)
        end

        for o = 2:Order
            println("Order $o")
            for (idx, q) in enumerate(extQ)
                q = q[1]
                @printf("%10.6f  %10.6f ± %10.6f\n", q / basic.kF, avg[o, idx], std[o, idx])
            end
        end
    end
end

run(steps)
# @time run(Steps)