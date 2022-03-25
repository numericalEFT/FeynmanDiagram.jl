@inline apply(o::Sum, diags::Vector{Diagram{W}}) where {W<:Number} = sum(d.weight for d in diags)
@inline apply(o::Prod, diags::Vector{Diagram{W}}) where {W<:Number} = prod(d.weight for d in diags)
@inline apply(o::Sum, diag::Diagram{W}) where {W<:Number} = diag.weight
@inline apply(o::Prod, diag::Diagram{W}) where {W<:Number} = diag.weight

@inline eval(d::DiagramId) = error("eval for $d has not yet implemented!")

function evalDiagNode!(diag::Diagram, varK, varT, evalBare::Function)
    if length(diag.subdiagram) == 0
        if hasproperty(diag.id, :extK)
            K = varK * diag.id.extK
            diag.weight = evalBare(diag.id, K, diag.id.extT, varT) * diag.factor
        else
            diag.weight = evalBare(diag.id, nothing, diag.id.extT, varT) * diag.factor
        end
    else
        diag.weight = apply(diag.operator, diag.subdiagram) * diag.factor
    end
    return diag.weight
end

function evalDiagTree!(diag::Diagram, varK, varT, evalBare::Function)
    for d in PostOrderDFS(diag)
        evalDiagNode!(d, varK, varT, evalBare)
    end
    return diag.weight
end

function evalDiagTree!(diags::Vector{Diagram{W}}, varK, varT, evalBare::Function) where {W}
    for d in diags
        evalDiagTree!(d, varK, varT, evalBare)
    end
    # return W[d.weight for d in diags]
end

function evalDiagTree!(df::DataFrame, varK, varT, evalBare::Function) where {W}
    for d in df[!, :diagram]
        evalDiagTree!(d, varK, varT, evalBare)
    end
    # return W[d.weight for d in df[!, :Diagram]]
end

function evalDiagNode!(diag::Diagram, evalBare::Function, vargs...)
    if length(diag.subdiagram) == 0
        diag.weight = evalBare(diag.id, vargs...) * diag.factor
    else
        diag.weight = apply(diag.operator, diag.subdiagram) * diag.factor
    end
    return diag.weight
end

function evalDiagTree!(diag::Diagram, evalBare::Function, vargs...)
    for d in PostOrderDFS(diag)
        evalDiagNode!(d, evalBare, vargs...)
    end
    return diag.weight
end

function evalDiagTree!(diags::Vector{Diagram{W}}, evalBare::Function, vargs...) where {W}
    for d in diags
        evalDiagTree!(d, evalBare, vargs...)
    end
    return W[d.weight for d in diags]
end

function evalDiagTree!(df::DataFrame, evalBare::Function, vargs...) where {W}
    for d in df[!, :diagram]
        evalDiagTree!(d, evalBare, vargs...)
    end
    # return W[d.weight for d in df[!, :Diagram]]
end