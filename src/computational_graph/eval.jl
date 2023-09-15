#@inline add( diags::Vector{Graph{F,W}}) where {F<:Number,W<:Number} = sum(d.weight for d in diags)
@inline apply(::Type{Sum}, diags::Vector{Graph{F,W}}, factors::Vector{F}) where {F<:Number,W<:Number} = sum(d.weight*d.factor*f for (d,f) in zip(diags,factors) )
@inline apply(::Type{Prod}, diags::Vector{Graph{F,W}}, factors::Vector{F}) where {F<:Number,W<:Number} = prod(d.weight*d.factor*f for (d,f) in zip(diags,factors) )
@inline apply(o::Sum, diag::Graph{F,W}) where {F<:Number,W<:Number} = diag.weight
@inline apply(o::Prod, diag::Graph{F,W}) where {F<:Number,W<:Number} = diag.weight

# #@inline eval(d::DiagramId) = error("eval for $d has not yet implemented!")
# #
function evalGraph!(g::Graph)
    result = nothing
    for node in PostOrderDFS(g)
        if isleaf(node)
            node.weight = 1.0 ##Currently set to 1 just for test. In the future, probably the whole function is not needed since we have the compiler.
        else
            node.weight = apply(node.operator, node.subgraphs, node.subgraph_factors) 
            #node.weight = add(node.subgraphs) * node.factor
        end
        result = node.weight * node.factor
    end
    return result
end

function evalGraph!(g::Number)
    return g
end
