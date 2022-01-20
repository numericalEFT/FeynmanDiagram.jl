@inline apply(o::Sum, diags::Vector{Diagram{W}}) where {W<:Number} = sum(d.weight for d in diags)
@inline apply(o::Prod, diags::Vector{Diagram{W}}) where {W<:Number} = prod(d.weight for d in diags)
@inline apply(o::Sum, diag::Diagram{W}) where {W<:Number} = diag.weight
@inline apply(o::Prod, diag::Diagram{W}) where {W<:Number} = diag.weight

function evalDiagNode!(diag::Diagram, evalBare::Function, vargs...; kwargs...)
    if isbare(diag)
        diag.weight = evalBare(diag.id, vargs...; kwargs...) * diag.factor
    else
        diag.weight = apply(diag.operator, diag.subdiagram) * diag.factor
    end
    return diag.weight
end

function evalDiagTree!(diag::Diagram, evalBare::Function, vargs...; kwargs...)
    for d in PostOrderDFS(diag)
        evalDiagNode!(d, evalBare, vargs...; kwargs...)
    end
    return diag.weight
end

function evalDiagTree!(diags::Vector{Diagram{W}}, evalBare::Function, vargs...; kwargs...) where {W}
    for d in diags
        evalDiagTree!(d, evalBare, vargs...; kwargs...)
    end
    return W[d.weight for d in diags]
end