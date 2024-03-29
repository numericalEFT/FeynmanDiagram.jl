using FeynmanDiagram
using Lehmann
using LinearAlgebra

const type = Ver4Diag
const Order = 4
const Circle = 100000

const kF = 1.919
const β = 25.0 / kF^2
const Λs = 1.0
const dim = 3

function benchmark(tree, N, varK, varT)
    for i in 1:N
        ExprTree.evalKT!(tree, varK, varT) #evaluate the expression tree
        # ExprTree.evalTree!(tree; loopVar=varK, siteVar=varT) #evaluate the expression tree
    end
end

# normalize the momentum by assuming mirror symmetry
function normalize!(id::DiagramId)
    if id isa BareGreenId || id isa BareInteractionId
        idx = findfirst(x -> (x ≈ 0.0) == false, id.extK)
        sig = sign(id.extK[idx])
        id.extK ./= sig
    else
        error("$id is not expected!")
    end
end

@inbounds function DiagTree.eval(id::BareGreenId, K, extT, varT)
    # function green(id::BareGreenId, K, extT, varT)
    τin, τout = varT[extT[1]], varT[extT[2]]
    # ϵ = dot(K, K) - kF^2
    ϵ = K[1] * K[1] + K[2] * K[2] + K[3] * K[3] - kF^2
    τ = τout - τin
    if τ ≈ 0.0
        return Spectral.kernelFermiT(-1e-8, ϵ, β)
    else
        return Spectral.kernelFermiT(τ, ϵ, β)
    end
end

# interaction(id::BareInteractionId, K, extT, varT) = 8π / (K[1] * K[1] + K[2] * K[2] + K[3] * K[3] + Λs)
@inbounds DiagTree.eval(id::BareInteractionId, K, extT, varT) = 8π / (K[1] * K[1] + K[2] * K[2] + K[3] * K[3] + Λs)
# DiagTree.eval(id::BareInteractionId, K, extT, varT) = 8π / (dot(K, K) + Λs)
# DiagTree.eval(id, K, extT, varT) = 1.0

diagPara(order) = DiagParaF64(type=type, innerLoopNum=order, hasTau=true,
    interaction=[FeynmanDiagram.Interaction(ChargeCharge, Instant),],  #instant charge-charge interaction
    # filter = [NoFock,])
    filter=[NoHartree, Girreducible,])

println("Build the diagrams into an experssion tree ...")
const para = [diagPara(o) for o in 1:Order]

if type == SigmaDiag
    diags = [Parquet.sigma(para[i]) for i in 1:1]
elseif type == PolarDiag
    diags = [Parquet.polarization(para[i]) for i in 1:1]
elseif type == Ver4Diag
    diags = [Parquet.vertex4(para[i]) for i in 1:1]
else
    error("not implemented!")
end

#diagram of different orders
start = time()
if type == SigmaDiag
    diags = [Parquet.sigma(para[i]) for i in 1:Order]
elseif type == PolarDiag
    diags = [Parquet.polarization(para[i]) for i in 1:Order]
elseif type == Ver4Diag
    diags = [Parquet.vertex4(para[i]) for i in 1:Order]
else
    error("not implemented!")
end
println("Build DiagTree takes ", time() - start, " seconds")

const extT = [diags[o].extT for o in 1:Order]                        #external tau of each diagram

println("Building tree")
start = time()
# const tree = [ExprTree.build(mergeby(diags[o].diagram); verbose=1, normalize=normalize!) for o in 1:Order]     #experssion tree representation of diagrams 
const tree = [ExprTree.build(mergeby(diags[o].diagram, dim); verbose=1) for o in 1:Order]     #experssion tree representation of diagrams 
println("Done.")
# const tree = [ExprTree.build(DiagTree.optimize(mergeby(diags[o]).diagram[1])) for o in 1:Order]     #experssion tree representation of diagrams 
println("Build ExprTree takes ", time() - start, " seconds")

varK = rand(3, 16)
varT = [rand() for i in 1:8]

benchmark(tree[1], 100, varK, varT)
println("Run evaluation for $Circle times.")
for i in 1:Order
    @time benchmark(tree[i], Circle, varK, varT)
end