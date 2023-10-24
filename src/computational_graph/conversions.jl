
"""
    function Base.convert(Graph, g::FeynmanGraph)
    
    Converts a FeynmanGraph `g` into a Graph, discarding its Feynman properties.
    After conversion, graph `g` is no longer guaranteed to be a valid (group of) Feynman diagram(s).

    # Arguments:
    - `g`  computational graph
"""
function Base.convert(::Type{Graph{F,W}}, g::FeynmanGraph{F,W}) where {F,W}
    return Graph{F,W}(g.subgraphs; subgraph_factors=g.subgraph_factors, name=g.name, operator=g.operator, orders=g.orders, factor=g.factor, weight=g.weight)
end

function Base.convert(::Type{FeynmanGraph{F,W}}, g::Graph{F,W}) where {F,W}
    error(
        "A set of Feynman properties (operator vertices, topology, etc.) must be specified to convert an object of type Graph to FeynmanGraph. " *
        "Please use constructor `FeynmanGraph(g::Graph, properties::FeynmanProperties)` instead."
    )
end

# Automatically promote FeynmanGraph to Graph for arithmetic operations
Base.promote_rule(::Type{Graph{F,W}}, ::Type{FeynmanGraph{F,W}}) where {F,W} = Graph{F,W}

# Arithmetic operations for mixed Graph/FeynmanGraph types
linear_combination(g1::Graph{F,W}, g2::FeynmanGraph{F,W}, c1, c2) where {F,W} = linear_combination(Base.promote(g1, g2)..., c1, c2)
linear_combination(g1::FeynmanGraph{F,W}, g2::Graph{F,W}, c1, c2) where {F,W} = linear_combination(Base.promote(g1, g2)..., c1, c2)
linear_combination(graphs::Vector{Union{Graph{F,W},FeynmanGraph{F,W}}}, constants::AbstractVector) where {F,W} = linear_combination(Base.promote(graphs)..., constants)
Base.:*(g1::Graph, g2::FeynmanGraph) = error("Multiplication of Feynman graphs is not well defined!")
Base.:*(g1::FeynmanGraph, g2::Graph) = error("Multiplication of Feynman graphs is not well defined!")
Base.:+(g1::Graph{F,W}, g2::FeynmanGraph{F,W}) where {F,W} = Base.:+(Base.promote(g1, g2)...)
Base.:+(g1::FeynmanGraph{F,W}, g2::Graph{F,W}) where {F,W} = Base.:+(Base.promote(g1, g2)...)
Base.:-(g1::Graph{F,W}, g2::FeynmanGraph{F,W}) where {F,W} = Base.:+(Base.promote(g1, g2)...)
Base.:-(g1::FeynmanGraph{F,W}, g2::Graph{F,W}) where {F,W} = Base.:+(Base.promote(g1, g2)...)
