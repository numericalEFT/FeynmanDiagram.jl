"""
This example demonstrated use NumericalEFT packages to calculate the multi-loop polarization diagrams
of electrons with Yukawa interaction.

If you want to run it with MPI, simply use:
mpiexec -n julia polarization.jl

Note that you may need to install julia MPI.jl following the link:
https://juliaparallel.github.io/MPI.jl/stable/configuration/#Julia-wrapper-for-mpiexec

With the parameters:
const steps = 1e7 # MC steps of each block
const Order = 3  #diagram order
const dim = 3
const rs = 1.0
const beta = 25.0  # β*E_F
const basic = Parameter.rydbergUnit(1 / beta, rs, dim, Λs = 1.0) # calculate all relevant parameters 

We expect the outputs:
Order 1
                q/kF             polarMC                  error               exact
            0.000000            0.097305 ±             0.000200           -0.048613
            0.428571            0.095600 ±             0.000220           -0.047862
            0.857143            0.091110 ±             0.000231           -0.045518
            1.285714            0.082607 ±             0.000290           -0.041234
            1.714286            0.067406 ±             0.000402           -0.033955
            2.142857            0.037449 ±             0.000305           -0.018655
            2.571429            0.022721 ±             0.000260           -0.011470
            3.000000            0.015731 ±             0.000219           -0.008007
Order 2
                q/kF             polarMC                  error               exact
            0.000000            0.022179 ±             0.000088
            0.428571            0.021614 ±             0.000136
            0.857143            0.019902 ±             0.000163
            1.285714            0.017330 ±             0.000189
            1.714286            0.012746 ±             0.000232
            2.142857            0.003328 ±             0.000081
            2.571429            0.001054 ±             0.000043
            3.000000            0.000446 ±             0.000023
Order 3
                q/kF             polarMC                  error               exact
            0.000000            0.007674 ±             0.000405
            0.428571            0.007842 ±             0.000414
            0.857143            0.007090 ±             0.000354
            1.285714            0.005093 ±             0.000366
            1.714286            0.003780 ±             0.000281
            2.142857            0.003119 ±             0.000291
            2.571429            0.001502 ±             0.000214
            3.000000            0.000768 ±             0.000151
"""

using Printf, LinearAlgebra
using MCIntegration, FeynmanDiagram, ElectronGas, Lehmann #NumericalEFT packages

##################### parameters for 3D UEG ##############################
const steps = 1e6 # MC steps of each block
const Order = 3  #diagram order
const dim = 3
const rs = 1.0
const beta = 25.0  # β*E_F
const basic = Parameter.rydbergUnit(1 / beta, rs, dim, Λs = 1.0) # calculate all relevant parameters 
const β, kF, μ, me, spin = basic.β, basic.kF, basic.μ, basic.me, basic.spin

##################### parameters for polarization ##############################
const n = 0        # external Matsubara frequency
const Qsize = 8    # numebr of external momentum points
const extQ = [[q, 0.0, 0.0] for q in LinRange(0.0, 3kF, Qsize)] #samples of external momentums

###################  parameter for polarization diagram #######################
diagPara(order) = GenericPara(diagType = PolarDiag, innerLoopNum = order, hasTau = true, loopDim = dim, spin = spin,
    interaction = [FeynmanDiagram.Interaction(ChargeCharge, Instant),],  #instant charge-charge interaction
    filter = [NoFock,])

println("Build the diagrams into an experssion tree ...")
const para = [diagPara(o) for o in 1:Order]
const polar = [mergeby(Parquet.polarization(para[i])).diagram[1] for i in 1:Order]   #diagram of different orders
#different order has different set of K, T variables, thus must have different exprtrees
const diag = [ExprTree.build(polar[o]) for o in 1:Order]    #experssion tree representation of diagrams 

##################### propagator and interaction evaluation ##############
function eval(id::BareGreenId, K, extT, varT)
    τin, τout = varT[extT[1]], varT[extT[2]]
    ϵ = dot(K, K) / (2me) - μ
    τ = τout - τin
    if τ ≈ 0.0
        return Spectral.kernelFermiT(-1e-8, ϵ, β)
    else
        return Spectral.kernelFermiT(τ, ϵ, β)
    end
end

eval(id::BareInteractionId, K, extT, varT) = (basic.e0)^2 / basic.ϵ0 / (dot(K, K) + basic.Λs)

################### interface to MC #########################################
function integrand(config)
    @assert config.curr >= 0
    order = config.curr
    #T: all internal Tau variables, K[>=2]: internal loop variables, Ext[1]: the index of the external momentum
    T, K, Ext = config.var[1], config.var[2], config.var[3]
    # copy the external momentum into the K.data matrix, which will be used to evaluate diagrams in the next line
    K.data[:, 1] .= extQ[Ext[1]]
    # K.data[:, 1]: external K, K.daata[:, >=2]: internal K, so that K.data contains all momentum
    ExprTree.evalNaive!(diag[order], K.data, T, eval) #evaluate the expression tree
    # there is an additional factor 1/β because we are integrating over both the incoming and the outing Tau variables of the poalrization
    return diag[order][1] * cos(2π * n * (T[2] - T[1]) / β) / β * spin
end

function measure(config)
    obs = config.observable
    factor = 1.0 / config.reweight[config.curr]
    weight = integrand(config)
    extidx = config.var[3][1]
    obs[config.curr, extidx] += weight / abs(weight) * factor
end

function run(steps)
    println("Initializing MC configuration ...")
    T = Tau(β, β / 2.0)  # Tau varaibles
    K = FermiK(dim, kF, 0.2 * kF, 10kF, offset = 1)  #momentum variables
    Ext = MCIntegration.Discrete(1, length(extQ)) # external variable is specified

    # degrees of freedom of the diagrams of different orders
    dof = [[para[o].totalTauNum, para[o].innerLoopNum, 1] for o in 1:Order]
    # observable for the diagrams of different orders
    obs = zeros(Float64, (Order, Qsize))

    config = Configuration(steps, (T, K, Ext), dof, obs)
    println("Start MC sampling ...")
    avg, std = MCIntegration.sample(config, integrand, measure; print = 0, Nblock = 16, reweight = 10000)

    if isnothing(avg) == false #if run with MPI, then only the master node has meaningful avg
        for o = 1:Order
            println("Order $o")
            @printf("%20s%20s   %20s%20s\n", "q/kF", "polarMC", "error", "exact")
            for (idx, q) in enumerate(extQ)
                if o == 1
                    p = Polarization.Polarization0_ZeroTemp(q[1], n, basic)
                    @printf("%20.6f%20.6f ± %20.6f%20.6f\n", q[1] / kF, avg[o, idx], std[o, idx], p)
                else
                    @printf("%20.6f%20.6f ± %20.6f\n", q[1] / kF, avg[o, idx], std[o, idx])
                end
            end
        end
    end
end

run(steps)
# @time run(Steps)