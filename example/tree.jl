using FeynmanDiagram
using InteractiveUtils, LinearAlgebra
using Lehmann

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

evalPropagator(id::GreenId, K, Tbasis, varT) = evalG(K, varT[Tbasis[1]], varT[Tbasis[2]])
evalPropagator(id::InteractionId, K, Tbasis, varT) = evalV(K)

para = GenericPara(diagType = Ver4Diag, innerLoopNum = 3, filter = [Girreducible,], hasTau = true)
ver4 = Parquet.vertex4(para)
ver4 = mergeby(ver4, :response)
println(ver4)
tree, root = ExprTree.compile(ver4.diagram)

varK = rand(3, para.totalLoopNum)
varT = [rand() for i in 1:para.totalTauNum]

ExprTree.evalNaive!(tree, varK, varT, evalPropagator)
@time ExprTree.evalNaive!(tree, varK, varT, evalPropagator)

evalDiagTree!(ver4, varK, varT, evalPropagator)
@time evalDiagTree!(ver4, varK, varT, evalPropagator)
