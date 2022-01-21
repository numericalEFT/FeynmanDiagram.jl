@inline apply(o::Sum, diags::Vector{Diagram{W}}) where {W<:Number} = sum(d.weight for d in diags)
@inline apply(o::Prod, diags::Vector{Diagram{W}}) where {W<:Number} = prod(d.weight for d in diags)
@inline apply(o::Sum, diag::Diagram{W}) where {W<:Number} = diag.weight
@inline apply(o::Prod, diag::Diagram{W}) where {W<:Number} = diag.weight

function evalDiagNode!(diag::Diagram, evalBare::Function, varK, varT)
    if length(diag.subdiagram) == 0
        diag.weight = evalBare(diag.id, varK, varT) * diag.factor
    else
        diag.weight = apply(diag.operator, diag.subdiagram) * diag.factor
    end
    return diag.weight
end

function evalDiagTree!(diag::Diagram, evalBare::Function, varK, varT)
    for d in PostOrderDFS(diag)
        evalDiagNode!(d, evalBare, varK, varT)
    end
    return diag.weight
end

function evalDiagTree!(diags::Vector{Diagram{W}}, evalBare::Function, varK, varT) where {W}
    for d in diags
        evalDiagTree!(d, evalBare, varK, varT)
    end
    # return W[d.weight for d in diags]
end

function evalDiagTree!(df::DataFrame, evalBare::Function, varK, varT) where {W}
    for d in df[!, :diagram]
        evalDiagTree!(d, evalBare, varK, varT)
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