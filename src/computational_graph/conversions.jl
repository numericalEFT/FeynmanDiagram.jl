
"""
    function Base.convert(::Type{G}, g::FeynmanGraph{F,W}) where {F,W,G<:Graph}
    
    Converts a FeynmanGraph `g` into a Graph, discarding its Feynman properties.
    After conversion, graph `g` is no longer guaranteed to be a valid (group of) Feynman diagram(s).

    # Arguments:
    - `g`  computational graph
"""
function Base.convert(::Type{G}, g::FeynmanGraph{F,W}) where {F,W,G<:Graph}
    return Graph(g.subgraphs; subgraph_factors=g.subgraph_factors, name=g.name, operator=g.operator(), orders=g.orders, ftype=F, wtype=W, factor=g.factor, weight=g.weight)
end

function Base.convert(::Type{FeynmanGraph}, g::Graph{F,W}) where {F,W}
    error(
        "A set of Feynman properties (operator vertices, topology, etc.) must be specified to convert an object of type Graph to FeynmanGraph. " *
        "Please use constructor `FeynmanGraph(g::Graph, properties::FeynmanProperties)` instead."
    )
end

# Automatically promote FeynmanGraph to Graph for vectors v::Vector{Union{Graph,FeynmanGraph}} and arithmetic operations
Base.promote_rule(::Type{Graph{F,W}}, ::Type{FeynmanGraph{F,W}}) where {F,W} = Graph{F,W}

# Arithmetic operations for mixed Graph/FeynmanGraph types
linear_combination(g1::Graph{F,W}, g2::FeynmanGraph{F,W}, c1=F(1), c2=F(1)) where {F,W} = linear_combination(Base.promote(g1, g2)..., c1, c2)
linear_combination(g1::FeynmanGraph{F,W}, g2::Graph{F,W}, c1=F(1), c2=F(1)) where {F,W} = linear_combination(Base.promote(g1, g2)..., c1, c2)
Base.:+(g1::Graph{F,W}, g2::FeynmanGraph{F,W}) where {F,W} = Base.:+(Base.promote(g1, g2)...)
Base.:+(g1::FeynmanGraph{F,W}, g2::Graph{F,W}) where {F,W} = Base.:+(Base.promote(g1, g2)...)
Base.:-(g1::Graph{F,W}, g2::FeynmanGraph{F,W}) where {F,W} = Base.:-(Base.promote(g1, g2)...)
Base.:-(g1::FeynmanGraph{F,W}, g2::Graph{F,W}) where {F,W} = Base.:-(Base.promote(g1, g2)...)

# No auto-promotion for graph products
multi_product(g1::Graph{F,W}, g2::FeynmanGraph{F,W}, c1=F(1), c2=F(1)) where {F,W} = error("Multiplication of Feynman graphs is not well defined!")
multi_product(g1::FeynmanGraph{F,W}, g2::Graph{F,W}, c1=F(1), c2=F(1)) where {F,W} = error("Multiplication of Feynman graphs is not well defined!")
Base.:*(g1::Graph{F,W}, g2::FeynmanGraph{F,W}) where {F,W} = error("Multiplication of Feynman graphs is not well defined!")
Base.:*(g1::FeynmanGraph{F,W}, g2::Graph{F,W}) where {F,W} = error("Multiplication of Feynman graphs is not well defined!")
