@inline apply(::Type{Sum}, diags::Vector{Graph{F,W}}, factors::Vector{F}) where {F<:Number,W<:Number} = sum(d.weight * f for (d, f) in zip(diags, factors))
@inline apply(::Type{Prod}, diags::Vector{Graph{F,W}}, factors::Vector{F}) where {F<:Number,W<:Number} = prod(d.weight * f for (d, f) in zip(diags, factors))
@inline apply(::Type{Power{N}}, diags::Vector{Graph{F,W}}, factors::Vector{F}) where {N,F<:Number,W<:Number} = (diags[1].weight)^N * factors[1]
@inline apply(o::Sum, diag::Graph{F,W}) where {F<:Number,W<:Number} = diag.weight
@inline apply(o::Prod, diag::Graph{F,W}) where {F<:Number,W<:Number} = diag.weight
@inline apply(o::Power{N}, diag::Graph{F,W}) where {N,F<:Number,W<:Number} = diag.weight

@inline apply(::Type{Sum}, diags::Vector{FeynmanGraph{F,W}}, factors::Vector{F}) where {F<:Number,W<:Number} = sum(d.weight * f for (d, f) in zip(diags, factors))
@inline apply(::Type{Prod}, diags::Vector{FeynmanGraph{F,W}}, factors::Vector{F}) where {F<:Number,W<:Number} = prod(d.weight * f for (d, f) in zip(diags, factors))
@inline apply(::Type{Power{N}}, diags::Vector{FeynmanGraph{F,W}}, factors::Vector{F}) where {N,F<:Number,W<:Number} = (diags[1].weight)^N * factors[1]
@inline apply(o::Sum, diag::FeynmanGraph{F,W}) where {F<:Number,W<:Number} = diag.weight
@inline apply(o::Prod, diag::FeynmanGraph{F,W}) where {F<:Number,W<:Number} = diag.weight
@inline apply(o::Power{N}, diag::FeynmanGraph{F,W}) where {N,F<:Number,W<:Number} = diag.weight

function eval!(g::Graph{F,W}, leafmap::Dict{Int,Int}=Dict{Int,Int}(), leaf::Vector{W}=Vector{W}(); inherit=false, randseed::Int=-1) where {F,W}
    result = nothing
    if randseed > 0
        Random.seed!(randseed)
    end
    for node in PostOrderDFS(g)
        if isleaf(node)
            if !inherit
                if isempty(leafmap)
                    if randseed < 0
                        # node.weight = 1.0
                        node.weight = -3.0
                    else
                        node.weight = rand()
                    end
                else
                    node.weight = leaf[leafmap[node.id]]
                end
            end
        else
            node.weight = apply(node.operator, node.subgraphs, node.subgraph_factors)
        end
        result = node.weight
    end
    return result
end


function eval!(g::FeynmanGraph{F,W}, leafmap::Dict{Int,Int}=Dict{Int,Int}(), leaf::Vector{W}=Vector{W}(); inherit=false, randseed::Int=-1) where {F,W}
    result = nothing
    if randseed > 0
        Random.seed!(randseed)
    end
    for node in PostOrderDFS(g)
        if isleaf(node)
            if !inherit
                if isempty(leafmap)
                    if randseed < 0
                        node.weight = 1.0
                    else
                        node.weight = rand()
                    end
                else
                    node.weight = leaf[leafmap[node.id]]
                end
            end
        else
            node.weight = apply(node.operator, node.subgraphs, node.subgraph_factors)
        end
        result = node.weight
    end
    return result
end

function eval!(g::Number)
    return g
end

function eval!(nothing)
    return nothing
end
