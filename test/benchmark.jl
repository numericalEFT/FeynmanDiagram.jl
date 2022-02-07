using FeynmanDiagram
using Lehmann
using LinearAlgebra

const diagType = Ver4Diag
const Order = 3
const Circle = 100000

const kF = 1.919
const β = 25.0 / kF^2
const Λs = 1.0

function benchmark(tree, N, varK, varT)
    for i in 1:N
        ExprTree.evalNaive!(tree, varK, varT, eval) #evaluate the expression tree
    end
end

function eval(id::GreenId, K, varT)
    τin, τout = varT[id.extT[1]], varT[id.extT[2]]
    ϵ = dot(K, K) - kF^2
    τ = τout - τin
    if τ ≈ 0.0
        return Spectral.kernelFermiT(-1e-8, ϵ, β)
    else
        return Spectral.kernelFermiT(τ, ϵ, β)
    end
end

eval(id::InteractionId, K, varT) = 8π / (dot(K, K) + Λs)

diagPara(order) = GenericPara(diagType = diagType, innerLoopNum = order, hasTau = true,
    interaction = [FeynmanDiagram.Interaction(ChargeCharge, Instant),],  #instant charge-charge interaction
    # filter = [NoFock,])
    filter = [Girreducible,])

println("Build the diagrams into an experssion tree ...")
const para = [diagPara(o) for o in 1:Order]

#diagram of different orders
if diagType == SigmaDiag
    diags = [Parquet.sigma(para[i]) for i in 1:Order]
elseif diagType == PolarDiag
    diags = [Parquet.polarization(para[i]) for i in 1:Order]
elseif diagType == Ver4Diag
    diags = [Parquet.vertex4(para[i]) for i in 1:Order]
else
    error("not implemented!")
end
const extT = [diags[o].extT for o in 1:Order]                        #external tau of each diagram
const tree = [ExprTree.build(mergeby(diags[o]).diagram[1]) for o in 1:Order]     #experssion tree representation of diagrams 
# const tree = [ExprTree.build(DiagTree.optimize!(mergeby(diags[o]).diagram[1])) for o in 1:Order]     #experssion tree representation of diagrams 

varK = rand(3, 16)
varT = [rand() for i in 1:8]

benchmark(tree[1], 100, varK, varT)
println("Run evaluation for $Circle times.")
for i in 1:Order
    @time benchmark(tree[i], Circle, varK, varT)
end