@inline apply(o::Sum, diags::Vector{Diagram{W}}) where {W<:Number} = sum(d.weight for d in diags)
@inline apply(o::Prod, diags::Vector{Diagram{W}}) where {W<:Number} = prod(d.weight for d in diags)
@inline apply(o::Sum, diag::Diagram{W}) where {W<:Number} = diag.weight
@inline apply(o::Prod, diag::Diagram{W}) where {W<:Number} = diag.weight

@inline eval(d::DiagramId) = error("eval for $d has not yet implemented!")

######################### evaluator for KT representation ######################### 
function evalDiagNodeKT!(diag::Diagram, varK, varT, additional=nothing; eval=DiagTree.eval)
    if length(diag.subdiagram) == 0
        # if hasproperty(diag.id, :extK)
        if (isnothing(varK) == false) && (isnothing(varT) == false)
            K = varK * diag.id.extK
            if isnothing(additional)
                diag.weight = eval(diag.id, K, diag.id.extT, varT) * diag.factor
            else
                diag.weight = eval(diag.id, K, diag.id.extT, varT, additional) * diag.factor
            end
        elseif isnothing(varK)
            if isnothing(additional)
                diag.weight = eval(diag.id, diag.id.extT, varT) * diag.factor
            else
                diag.weight = eval(diag.id, diag.id.extT, varT, additional) * diag.factor
            end
        elseif isnothing(varT)
            K = varK * diag.id.extK
            if isnothing(additional)
                diag.weight = eval(diag.id, K) * diag.factor
            else
                diag.weight = eval(diag.id, K, additional) * diag.factor
            end
        else
            if isnothing(additional)
                diag.weight = eval(diag.id) * diag.factor
            else
                diag.weight = eval(diag.id, additional) * diag.factor
            end
        end
    else
        diag.weight = apply(diag.operator, diag.subdiagram) * diag.factor
    end
    return diag.weight
end

function evalKT!(diag::Diagram, varK, varT; eval=DiagTree.eval)
    for d in PostOrderDFS(diag)
        evalDiagNodeKT!(d, varK, varT; eval=eval)
    end
    return diag.weight
end

function evalKT!(diags::Vector{Diagram{W}}, varK, varT; eval=DiagTree.eval) where {W}
    for d in diags
        evalKT!(d, varK, varT; eval=eval)
    end
    # return W[d.weight for d in diags]
end

function evalKT!(df::DataFrame, varK, varT; eval=DiagTree.eval) where {W}
    for d in df[!, :diagram]
        evalKT!(d, varK, varT; eval=eval)
    end
    # return W[d.weight for d in df[!, :Diagram]]
end

######################### generic evaluator ######################### 
function evalDiagNode!(diag::Diagram, vargs...; eval=DiagTree.eval)
    if length(diag.subdiagram) == 0
        diag.weight = eval(diag.id, vargs...) * diag.factor
    else
        diag.weight = apply(diag.operator, diag.subdiagram) * diag.factor
    end
    return diag.weight
end

function eval!(diag::Diagram, vargs...; eval=DiagTree.eval)
    for d in PostOrderDFS(diag)
        evalDiagNode!(d, vargs...; eval=eval)
    end
    return diag.weight
end

function eval!(diags::Vector{Diagram{W}}, vargs...; eval=DiagTree.eval) where {W}
    for d in diags
        eval!(d, vargs...; eval=eval)
    end
    # return W[d.weight for d in diags]
end

function eval!(df::DataFrame, vargs...; eval=DiagTree.eval) where {W}
    for d in df[!, :diagram]
        eval!(d, vargs...; eval=eval)
    end
    # return W[d.weight for d in df[!, :Diagram]]
end