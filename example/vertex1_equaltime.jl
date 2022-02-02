"""
This example demonstrated use NumericalEFT packages to calculate the multi-loop equal-time sigma diagrams
of electrons with Yukawa interaction.

If you want to run it with MPI, simply use:
mpiexec -n julia polarization.jl

Note that you may need to install julia MPI.jl following the link:
https://juliaparallel.github.io/MPI.jl/stable/configuration/#Julia-wrapper-for-mpiexec
"""

using Printf, LinearAlgebra
using MCIntegration, FeynmanDiagram, ElectronGas, Lehmann, CompositeGrids #NumericalEFT packages
using JLD2

##################### parameters for 3D UEG ##############################
const steps = 1e8       # MC steps of each block
const Order = 3         #diagram order
const dim = 3
const rs = 1.0
const beta = 25.0       # β*E_F
# const diagType = PolarDiag                                         #build polarization diagram with Parquet algorithm
const diagType = SigmaDiag                                           #build sigma diagram with Parquet algorithm
@assert diagType == SigmaDiag || diagType == PolarDiag               #only support sigma or polarization
const isFermi = diagType == SigmaDiag ? true : false
const basic = Parameter.rydbergUnit(1 / beta, rs, dim, Λs = 1.0)     # calculate all relevant parameters 
const β, kF, μ, me, spin = basic.β, basic.kF, basic.μ, basic.me, basic.spin

##################### parameters for polarization ##############################
# const n = 0        # external Matsubara frequency
# numebr of external momentum points
const qgrid = CompositeGrid.CompositeLogGrid(:uniform, [0.0, 10kF], 4, 1.0, true, 8)
# qgrid from 0 to 10kF, 4 panels of log grid from dense to sparse and mininterval = 1, each panel contains 8 uniform points
println("qgrid: ", qgrid)
const Qsize = length(qgrid)
#samples of external momentums
const extQ = [[q, 0.0, 0.0] for q in qgrid]

###################  parameter for polarization diagram #######################
diagPara(order) = GenericPara(diagType = diagType, innerLoopNum = order, hasTau = true, loopDim = dim, spin = spin,
    interaction = [FeynmanDiagram.Interaction(ChargeCharge, Instant),],  #instant charge-charge interaction
    filter = [])
# filter = [Girreducible,])

println("Build the diagrams into an experssion tree ...")
const para = [diagPara(o) for o in 1:Order]

#diagram of different orders
if diagType == SigmaDiag
    diags = [Parquet.sigma(para[i]) for i in 1:Order]
elseif diagType == PolarDiag
    diags = [Parquet.polarization(para[i]) for i in 1:Order]
end
println(diags)
"""sh
output: 
DataFrames.DataFrame[1×4 DataFrame
 Row │ type       extT    diagram                            hash  
     │ Analytic…  Tuple…  Diagram…                           Int64 
─────┼─────────────────────────────────────────────────────────────
   1 │ Instant    (1, 1)  Σ Ins, t(1, 1)=0.0=4.031e-03⨁ (6…      9, 2×4 DataFrame
 Row │ type       extT    diagram                            hash  
     │ Analytic…  Tuple…  Diagram…                           Int64 
─────┼─────────────────────────────────────────────────────────────
   1 │ Instant    (1, 1)  Σ Ins, t(1, 1)=0.0=4.031e-03⨁ (2…     74
   2 │ Dynamic    (1, 2)  Σ Dyn, t(1, 2)=0.0=4.031e-03⨁ (7…     75, 3×4 DataFrame
 Row │ type       extT    diagram                            hash  
     │ Analytic…  Tuple…  Diagram…                           Int64 
─────┼─────────────────────────────────────────────────────────────
   1 │ Instant    (1, 1)  Σ Ins, t(1, 1)=0.0=4.031e-03⨁ (1…    630
   2 │ Dynamic    (1, 2)  Σ Dyn, t(1, 2)=0.0=4.031e-03⨁ (5…    631
   3 │ Dynamic    (1, 3)  Σ Dyn, t(1, 3)=0.0=4.031e-03⨁ (5…    632]
"""
#different order has different set of K, T variables, thus must have different exprtrees
# const t = [ExprTree.build(diags[o]).diagram for o in 1:Order]
t = [[ExprTree.build(diags[o].diagram[t]) for t in 1:o] for o in 1:Order]
# const extTidx = [[diags[o].extT[t][2] for t in 1:o] for o in 1:Order]
const tree = [t[1][1], t[2][1], t[2][2], t[3][1], t[3][2], t[3][3]]     #experssion tree representation of diagrams 
const extTidx = [1, 1, 2, 1, 2, 3]                        #external tau of each diagram
# println(tree)
# println(extTidx)
# exit(0)

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

################### interface to MC #########################################
function integrand(config)
    @assert config.curr >= 0
    index = config.curr
    #T: all internal Tau variables, K[>=2]: internal loop variables, Ext[1]: the index of the external momentum
    T, K, ExtQ = config.var[1], config.var[2], config.var[3]
    # copy the external momentum into the K.data matrix, which will be used to evaluate diagrams in the next line
    K.data[:, 1] .= extQ[ExtQ[1]]
    # K.data[:, 1]: external K, K.daata[:, >=2]: internal K, so that K.data contains all momentum
    if extTidx[index] == 1
        weights = ExprTree.evalNaive!(tree[index], K.data, T, eval) #evaluate the expression tree
        return weights[1] / β
    else
        Tc = deepcopy(T)
        Tc[extTidx[index]] = Tc[1] - 1e-8  #reset the external tout to be the same as tin, so that we measure the equal time property
        weights = ExprTree.evalNaive!(tree[index], K.data, Tc, eval) #evaluate the expression tree
        return weights[1] / β / β
    end
    # w = sum(w * phasefactor(T[tidx[i][1]], T[tidx[i][2]], n) for (i, w) in enumerate(weights))
    # return weights[1] * cos(2π * ngrid[nidx] * (T[2] - T[1]) / β) / β * spin
end

function measure(config)
    obs = config.observable
    factor = 1.0 / config.reweight[config.curr]
    weight = integrand(config)
    qidx = config.var[3][1]
    obs[qidx, config.curr] += weight / abs(weight) * factor
end

function run(steps)
    println("Initializing MC configuration ...")
    T = Tau(β, β / 2.0)  # Tau varaibles
    K = FermiK(dim, kF, 0.2 * kF, 10kF, offset = 1)  # momentum variables
    ExtQ = MCIntegration.Discrete(1, Qsize)    # external q and external wn

    # degrees of freedom of the diagrams of different orders
    # dof = [[para[o].totalTauNum - 1, para[o].innerLoopNum, 1] for o in 1:Order]
    dof = [[1, 1, 1], [2, 2, 1], [2, 2, 1], [3, 3, 1], [3, 3, 1], [3, 3, 1]]
    # observable for the diagrams of different orders
    obs = zeros((Qsize, length(dof)))

    # config = Configuration(steps, (T, K, ExtQ, ExtN), dof, obs; reweight = [0.01, 0.02])
    config = Configuration(steps, (T, K, ExtQ), dof, obs)
    println("Start MC sampling ...")
    avg, std = MCIntegration.sample(config, integrand, measure; print = 5, Nblock = 16)

    if isnothing(avg) == false #if run with MPI, then only the master node has meaningful avg

        q = [_q[1] for _q in extQ]
        jldsave("sigma.jld2"; basic, qgrid, avg, std)

        printstyled("zero frequency results: \n", color = :green)
        for o = 1:length(dof)
            printstyled("Group $o\n", color = :yellow)
            # @printf("%20s%20s%20s%20s%20s%20s\n", "q/kF", "real", "error", "imag", "error", "exact")
            @printf("%10s%14s%12s%14s\n", "q/kF", "sigma", "error", "exact")
            for (idx, q) in enumerate(qgrid)
                if o == 1
                    if diagType == PolarDiag
                        p = Polarization.Polarization0_ZeroTemp(q, n, basic)
                    else
                        p = SelfEnergy.Fock0_ZeroTemp(q, basic)
                    end
                    @printf("%10.6f%14.8f%12.8f%14.8f\n", q / kF, real(avg[idx, o]), real(std[idx, o]), p)
                else
                    @printf("%10.6f%14.8f%12.8f\n", q / kF, real(avg[idx, o]), real(std[idx, o]))
                    # @printf("%20.6f%20.6f%20.6f%20.6f%20.6f\n", q[1] / kF, real(avg[o, idx]), real(std[o, idx]), imag(avg[o, idx]), imag(std[o, idx]))
                end
            end

        end
    end

end

run(steps)
# @time run(Steps)