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
using JLD2

##################### parameters for 3D UEG ##############################
const steps = 1e6       # MC steps of each block
const Order = 2         #diagram order
const dim = 3
const rs = 1.0
const beta = 25.0       # β*E_F
const diagType = PolarDiag                                         #build polarization diagram with Parquet algorithm
# const diagType = SigmaDiag                                           #build sigma diagram with Parquet algorithm
@assert diagType == SigmaDiag || diagType == PolarDiag               #only support sigma or polarization
const isFermi = diagType == SigmaDiag ? true : false
const basic = Parameter.rydbergUnit(1 / beta, rs, dim, Λs = 1.0)     # calculate all relevant parameters 
const β, kF, μ, me, spin = basic.β, basic.kF, basic.μ, basic.me, basic.spin

##################### parameters for polarization ##############################
# const n = 0        # external Matsubara frequency
# numebr of external momentum points
const Qsize = 8
#samples of external momentums
const extQ = [[q, 0.0, 0.0] for q in LinRange(0.0, 10kF, Qsize)]
#construct the optimized basis using discrete Lehmann representation
const dlr = DLRGrid(Euv = 10 * basic.EF, β = β, rtol = 1e-8, isFermi = true)
# const Nsize = dlr.size
#get the optimal Matsubara frequencies
# const ngrid = dlr.n
# println(dlr)
const ngrid = collect(0:10)
const Nsize = length(ngrid)

###################  parameter for polarization diagram #######################
diagPara(order) = GenericPara(diagType = diagType, innerLoopNum = order, hasTau = true, loopDim = dim, spin = spin,
    interaction = [FeynmanDiagram.Interaction(ChargeCharge, Instant),],  #instant charge-charge interaction
    # filter = [])
    filter = [Girreducible,])

println("Build the diagrams into an experssion tree ...")
const para = [diagPara(o) for o in 1:Order]

#diagram of different orders
if diagType == SigmaDiag
    diags = [Parquet.sigma(para[i]) for i in 1:Order]
elseif diagType == PolarDiag
    diags = [Parquet.polarization(para[i]) for i in 1:Order]
end
# println("order 1")
# println(diags[1])
# plot_tree(diags[1].diagram, maxdepth = 9)
# println("order 2")
# println(diags[2])
# println("order 3")
# println(diags[3])
# plot_tree(diags[2].diagram, maxdepth = 9)
#different order has different set of K, T variables, thus must have different exprtrees
const extTidx = [diags[o].extT for o in 1:Order]                        #external tau of each diagram
const tree = [ExprTree.build(diags[o].diagram) for o in 1:Order]     #experssion tree representation of diagrams 
# println(extT)
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

# there is an additional factor 1/β because we are integrating over both the incoming and the outing Tau variables of the poalrization
if diagType == SigmaDiag
    phasefactor(tin, tout, n) = exp(-1im * (2n + 1) * π * (tout - tin) / β) / β
elseif diagType == PolarDiag
    phasefactor(tin, tout, n) = exp(-1im * 2n * π * (tout - tin) / β) * spin / β
end

################### interface to MC #########################################
function integrand(config)
    @assert config.curr >= 0
    order = config.curr
    #T: all internal Tau variables, K[>=2]: internal loop variables, Ext[1]: the index of the external momentum
    T, K, ExtQ, ExtN = config.var[1], config.var[2], config.var[3], config.var[4]
    # copy the external momentum into the K.data matrix, which will be used to evaluate diagrams in the next line
    K.data[:, 1] .= extQ[ExtQ[1]]
    # K.data[:, 1]: external K, K.daata[:, >=2]: internal K, so that K.data contains all momentum
    weights = ExprTree.evalNaive!(tree[order], K.data, T, eval) #evaluate the expression tree
    n = ngrid[ExtN[1]]
    # println(nidx)
    tidx = extTidx[order]
    # println(tidx)
    w = sum(w * phasefactor(T[tidx[i][1]], T[tidx[i][2]], n) for (i, w) in enumerate(weights))
    return w
    # return weights[1] * cos(2π * ngrid[nidx] * (T[2] - T[1]) / β) / β * spin
end

function measure(config)
    obs = config.observable
    factor = 1.0 / config.reweight[config.curr]
    weight = integrand(config)
    qidx = config.var[3][1]
    nidx = config.var[4][1]
    obs[nidx, qidx, config.curr] += weight / abs(weight) * factor
end

function run(steps)
    println("Initializing MC configuration ...")
    T = Tau(β, β / 2.0)  # Tau varaibles
    K = FermiK(dim, kF, 0.2 * kF, 10kF, offset = 1)  # momentum variables
    ExtQ = MCIntegration.Discrete(1, Qsize)    # external q and external wn
    ExtN = MCIntegration.Discrete(1, Nsize)    # external q and external wn

    # degrees of freedom of the diagrams of different orders
    dof = [[para[o].totalTauNum, para[o].innerLoopNum, 1, 1] for o in 1:Order]
    # observable for the diagrams of different orders
    obs = zeros(ComplexF64, (Nsize, Qsize, Order))
    # obs = zeros(Float64, (Nsize, Qsize, Order))

    # config = Configuration(steps, (T, K, ExtQ, ExtN), dof, obs; reweight = [0.01, 0.02])
    config = Configuration(steps, (T, K, ExtQ, ExtN), dof, obs)
    println("Start MC sampling ...")
    avg, std = MCIntegration.sample(config, integrand, measure; print = 5, Nblock = 16)

    if isnothing(avg) == false #if run with MPI, then only the master node has meaningful avg

        jldsave("sigma.jld2"; dlr, avg, std)
        nidx = 1
        printstyled("zero frequency results: \n", color = :green)
        for o = 1:Order
            printstyled("Order $o\n", color = :yellow)
            # @printf("%20s%20s%20s%20s%20s%20s\n", "q/kF", "real", "error", "imag", "error", "exact")
            @printf("%10s%14s%12s%14s%12s%14s\n", "q/kF", "real", "error", "imag", "error", "exact")
            # for (idx, q) in enumerate(extQ)
            #     if o == 1
            #         if diagType == PolarDiag
            #             p = Polarization.Polarization0_ZeroTemp(q[1], n, basic)
            #         else
            #             p = SelfEnergy.Fock0_ZeroTemp(q[1], basic)
            #         end
            #         @printf("%10.6f%14.8f%12.8f%14.8f%14.6f%14.8f\n", q[1] / kF, real(avg[nidx, idx, o]), real(std[nidx, idx, o]),
            #             imag(avg[nidx, idx, o]), imag(std[nidx, idx, o]), p)
            #     else
            #         @printf("%10.6f%14.8f%12.8f%14.8f%14.6f\n", q[1] / kF, real(avg[nidx, idx, o]), real(std[nidx, idx, o]),
            #             imag(avg[nidx, idx, o]), imag(std[nidx, idx, o]))
            #         # @printf("%20.6f%20.6f%20.6f%20.6f%20.6f\n", q[1] / kF, real(avg[o, idx]), real(std[o, idx]), imag(avg[o, idx]), imag(std[o, idx]))
            #     end
            # end

            idx = 1
            q = extQ[1][1]
            for (nidx, n) in enumerate(ngrid)
                if o == 1
                    if diagType == PolarDiag
                        p = Polarization.Polarization0_ZeroTemp(q, n, basic)
                    else
                        p = SelfEnergy.Fock0_ZeroTemp(q, basic)
                    end
                    @printf("%10.6f%14.8f%12.8f%14.8f%14.6f%14.8f\n", n, real(avg[nidx, idx, o]), real(std[nidx, idx, o]),
                        imag(avg[nidx, idx, o]), imag(std[nidx, idx, o]), p)
                else
                    @printf("%10.6f%14.8f%12.8f%14.8f%14.6f\n", n, real(avg[nidx, idx, o]), real(std[nidx, idx, o]),
                        imag(avg[nidx, idx, o]), imag(std[nidx, idx, o]))
                    # @printf("%20.6f%20.6f%20.6f%20.6f%20.6f\n", q[1] / kF, real(avg[o, idx]), real(std[o, idx]), imag(avg[o, idx]), imag(std[o, idx]))
                end
            end
        end

        # println()
        # avg = matfreq2tau(dlr, avg, [0.0,], axis = 1)
        # printstyled("equal time results: \n", color = :green)
        # for o = 1:Order
        #     printstyled("Order $o\n", color = :yellow)
        #     @printf("%10s%14s%14s\n", "q/kF", "real", "imag")
        #     for (idx, q) in enumerate(extQ)
        #         @printf("%10.6f%14.8f%14.8f\n", q[1] / kF, real(avg[1, idx, o]), imag(avg[1, idx, o]))
        #     end
        # end
    end

end

run(steps)
# @time run(Steps)