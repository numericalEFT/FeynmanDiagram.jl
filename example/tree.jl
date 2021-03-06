using FeynmanDiagram
using InteractiveUtils, LinearAlgebra
using Lehmann
# using Profile

function evalG(K, τin, τout)
    # println(τBasis, ", ", varT)
    kF, β = 1.0, 1.0
    ϵ = dot(K, K) / 2 - kF^2
    τ = τout - τin
    if τ ≈ 0.0
        return Spectral.kernelFermiT(-1e-8, ϵ, β)
    else
        return Spectral.kernelFermiT(τ, ϵ, β)
    end
end

evalV(K) = 8π / (dot(K, K) + 1)

eval(id::GreenId, K, varT) = evalG(K, varT[id.extT[1]], varT[id.extT[2]])
eval(id::InteractionId, K, varT) = evalV(K)


para = GenericPara(diagType = Ver4Diag, innerLoopNum = 3, filter = [Girreducible,], hasTau = true)
ver4 = Parquet.vertex4(para)
ver4 = mergeby(ver4, :response)
println(ver4)
tree, root = ExprTree.build(ver4.diagram)

varK = rand(3, para.totalLoopNum)
varT = [rand() for i in 1:para.totalTauNum]

ExprTree.evalNaive!(tree, varK, varT, eval)
@time ExprTree.evalNaive!(tree, varK, varT, eval)

# ExprTree.warn_type(tree, varK, varT, eval)

# open("typeinfo.txt", "w") do f
#     @code_warntype(f, ExprTree.evalNaive!(tree, varK, varT, evalPropagator))
# end
# @profile ExprTree.evalNaive!(tree, varK, varT, evalPropagator)

# evalDiagTree!(ver4, varK, varT, evalPropagator)
# @time evalDiagTree!(ver4, varK, varT, evalPropagator)
