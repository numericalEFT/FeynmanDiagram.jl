"""
This example demonstrated use NumericalEFT packages to calculate the multi-loop sigma diagrams
of electrons with Yukawa interaction.

If you want to run it with MPI, simply use:
mpiexec -n julia polarization.jl

Note that you may need to install julia MPI.jl following the link:
https://juliaparallel.github.io/MPI.jl/stable/configuration/#Julia-wrapper-for-mpiexec
"""

using Printf, LinearAlgebra
using MCIntegration, FeynmanDiagram, ElectronGas, Lehmann #NumericalEFT packages

##################### parameters for 3D UEG ##############################
const steps = 1e6 # MC steps of each block
const Order = 3  #diagram order
const dim = 3
const rs = 1.0
const beta = 25.0  # β*E_F
const diagType = PolarDiag                 #build polarization diagram with Parquet algorithm
# const diagType = SigmaDiag                #build sigma diagram with Parquet algorithm
const basic = Parameter.rydbergUnit(1 / beta, rs, dim, Λs = 1.0) # calculate all relevant parameters 
const β, kF, μ, me, spin = basic.β, basic.kF, basic.μ, basic.me, basic.spin

##################### parameters for polarization ##############################
const n = 0        # external Matsubara frequency
const Qsize = 8    # numebr of external momentum points
const extQ = [[q, 0.0, 0.0] for q in LinRange(0.0, 3kF, Qsize)] #samples of external momentums

###################  parameter for polarization diagram #######################
diagPara(order) = GenericPara(diagType = diagType, innerLoopNum = order, hasTau = true, loopDim = dim, spin = spin,
    interaction = [FeynmanDiagram.Interaction(ChargeCharge, Instant),],  #instant charge-charge interaction
    filter = [NoFock,])

println("Build the diagrams into an experssion tree ...")
const para = [diagPara(o) for o in 1:Order]

#diagram of different orders
if diagType == SigmaDiag
    diags = [Parquet.sigma(para[i]) for i in 1:Order]
elseif diagType == PolarDiag
    diags = [Parquet.polarization(para[i]) for i in 1:Order]
end
#different order has different set of K, T variables, thus must have different exprtrees
const extT = [diags[o].extT for o in 1:Order]                        #external tau of each diagram
const tree = [ExprTree.build(diags[o].diagram) for o in 1:Order]     #experssion tree representation of diagrams 

##################### propagator and interaction evaluation ##############
function eval(id::GreenId, K, varT)
    τin, τout = varT[id.extT[1]], varT[id.extT[2]]
    ϵ = dot(K, K) / (2me) - μ
    τ = τout - τin
    if τ ≈ 0.0
        return Spectral.kernelFermiT(-1e-8, ϵ, β)
    else
        return Spectral.kernelFermiT(τ, ϵ, β)
    end
end

eval(id::InteractionId, K, varT) = (basic.e0)^2 / basic.ϵ0 / (dot(K, K) + basic.Λs)

# there is an additional factor 1/β because we are integrating over both the incoming and the outing Tau variables of the poalrization
if diagType == SigmaDiag
    phasefactor(extT) = sin((2n + 1) * π * (extT[2] - extT[1]) / β) / β
elseif diagType == PolarDiag
    phasefactor(extT) = cos(2n * π * (extT[2] - extT[1]) / β) * spin / β
end

################### interface to MC #########################################
function integrand(config)
    @assert config.curr >= 0
    order = config.curr
    #T: all internal Tau variables, K[>=2]: internal loop variables, Ext[1]: the index of the external momentum
    T, K, Ext = config.var[1], config.var[2], config.var[3]
    # copy the external momentum into the K.data matrix, which will be used to evaluate diagrams in the next line
    K.data[:, 1] .= extQ[Ext[1]]
    # K.data[:, 1]: external K, K.daata[:, >=2]: internal K, so that K.data contains all momentum
    weights = ExprTree.evalNaive!(tree[order], K.data, T, eval) #evaluate the expression tree
    return sum(w * phasefactor(extT[order][i]) for (i, w) in enumerate(weights))

    # return weight[1] * cos(2π * n * (T[2] - T[1]) / β) / β * spin
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
                    if diagType == PolarDiag
                        p = Polarization.Polarization0_ZeroTemp(q[1], n, basic)
                    else
                        p = 0.0
                    end
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