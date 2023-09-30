"""
    mutable struct Graph{F,W}
    
    A representation of a computational graph, e.g., an expression tree, with type stable node data.

# Members:
- `id::Int`  the unique hash id to identify the diagram
- `name::Symbol`  name of the diagram
- `subgraphs::Vector{Graph{F,W}}`  vector of sub-diagrams 
- `subgraph_factors::Vector{F}`  scalar multiplicative factors associated with each subgraph. Note that the subgraph factors may be manipulated algebraically. To associate a fixed multiplicative factor with this graph which carries some semantic meaning, use the `factor` argument instead.
- `operator::DataType`  node operation. Addition and multiplication are natively supported via operators Sum and Prod, respectively. Should be a concrete subtype of `AbstractOperator`.
- `factor::F`  total scalar multiplicative factor for the diagram
- `weight::W`  the weight of this node

# Example:
```julia-repl
julia> g1 = Graph([])
1=0.0

julia> g2 = Graph([]; factor=2)
2⋅2.0=0.0

julia> g = Graph([g1, g2]; operator=ComputationalGraphs.Sum())
3=0.0=⨁ (1,2)
```
"""
mutable struct Graph{F,W} <: AbstractGraph # Graph
    id::Int
    name::String # "" by default

    subgraphs::Vector{Graph{F,W}}
    subgraph_factors::Vector{F}

    operator::DataType
    factor::F
    weight::W

    """
        function Graph(subgraphs::AbstractVector; name="", operator::AbstractOperator=Sum(),
            ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype))
        
        Create a Graph struct from a set of subgraphs with the specified node data and operation.

    # Arguments:
    - `subgraphs`  vector of sub-diagrams 
    - `subgraph_factors`  scalar multiplicative factors associated with each subgraph. Note that the subgraph factors may be manipulated algebraically. To associate a fixed multiplicative factor with this graph which carries some semantic meaning, use the `factor` argument instead.
    - `name`  name of the diagram
    - `operator`  node operation, i.e., Sum, Prod, or a user-defined operator `Op <: AbstractOperator`
    - `ftype`  typeof(factor)
    - `wtype`  typeof(weight)
    - `factor`  fixed scalar multiplicative factor for this diagram (e.g., a permutation sign)
    - `weight`  the weight of this node
    """
    function Graph(subgraphs::AbstractVector; subgraph_factors=one.(eachindex(subgraphs)), name="", operator::AbstractOperator=Sum(),
        ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype)
    )
        return new{ftype,wtype}(uid(), name, subgraphs, subgraph_factors, typeof(operator), factor, weight)
    end
end

"""
    function constant_graph(factor=one(_dtype.factor))

    Returns a graph that represents a constant equal to f, where f is the factor with default value 1.

# Arguments:
- `f`:  constant factor
"""
function constant_graph(factor=one(_dtype.factor))
    return Graph([]; operator=Constant(), factor=factor, ftype=_dtype.factor, wtype=_dtype.weight, weight=one(_dtype.weight))
end


"""
    function Base.:*(g1::Graph{F,W}, c2::C) where {F,W,C}

    Returns a graph representing the scalar multiplication `g1*c2`.

# Arguments:
- `g1`  computational graph
- `c2`  scalar multiple
"""
function Base.:*(g1::Graph{F,W}, c2::C) where {F,W,C}
    g = Graph([g1,]; subgraph_factors=[F(c2),], operator=Prod(), ftype=F, wtype=W)
    # Merge multiplicative link
    if unary_istrivial(g1.operator) && onechild(g1)
        g.subgraph_factors[1] *= g1.subgraph_factors[1]
        g.subgraphs = g1.subgraphs
    end
    return g
end

"""
    function Base.:*(c1::C, g2::Graph{F,W}) where {F,W,C}

    Returns a graph representing the scalar multiplication `c1*g2`.

# Arguments:
- `c1`  scalar multiple
- `g2`  computational graph
"""
function Base.:*(c1::C, g2::Graph{F,W}) where {F,W,C}
    g = Graph([g2,]; subgraph_factors=[F(c1),], operator=Prod(), ftype=F, wtype=W)
    # Merge multiplicative link
    if unary_istrivial(g2.operator) && onechild(g2)
        g.subgraph_factors[1] *= g2.subgraph_factors[1]
        g.subgraphs = g2.subgraphs
    end
    return g
end

"""
    function linear_combination(g1::Graph{F,W}, g2::Graph{F,W}, c1::C, c2::C) where {F,W,C}

    Returns a graph representing the linear combination `c1*g1 + c2*g2`.

# Arguments:
- `g1`  first computational graph
- `g2`  second computational graph
- `c1`  first scalar multiple
- `c2`  second scalar multiple
"""
function linear_combination(g1::Graph{F,W}, g2::Graph{F,W}, c1::C, c2::C) where {F,W,C}
    g = Graph([g1, g2]; subgraph_factors=[F(c1), F(c2)], operator=Sum(), ftype=F, wtype=W)
    # Convert multiplicative links to in-place form
    if unary_istrivial(g1.operator) && onechild(g1)
        g.subgraph_factors[1] *= g1.subgraph_factors[1]
        g.subgraphs[1] = g1.subgraphs[1]
    end
    if unary_istrivial(g2.operator) && onechild(g2)
        g.subgraph_factors[2] *= g2.subgraph_factors[1]
        g.subgraphs[2] = g2.subgraphs[1]
    end
    return g
end

"""
    function linear_combination(graphs::Vector{Graph{F,W}}, constants::Vector{C}) where {F,W,C}

    Given a vector 𝐠 of graphs and an equally-sized vector 𝐜 of constants, returns a new
    graph representing the linear combination (𝐜 ⋅ 𝐠).

# Arguments:
- `graphs`  vector of computational graphs
- `constants`  vector of scalar multiples
"""
function linear_combination(graphs::Vector{Graph{F,W}}, constants::Vector{C}) where {F,W,C}
    g = Graph(graphs; subgraph_factors=constants, operator=Sum(), ftype=F, wtype=W)
    # Convert multiplicative links to in-place form
    for (i, sub_g) in enumerate(g.subgraphs)
        if unary_istrivial(sub_g.operator) && onechild(sub_g)
            g.subgraph_factors[i] *= sub_g.subgraph_factors[1]
            g.subgraphs[i] = sub_g.subgraphs[1]
        end
    end
    return g
end

# function Base.:+(c::C, g1::Graph{F,W}) where {F,W,C}
#     return linear_combination(g1, Unity, F(1), F(c))
# end
# function Base.:+(g1::Graph{F,W},c::C) where {F,W,C}
#     return linear_combination(g1, Unity, F(1), F(c))
# end

# function Base.:-(c::C, g1::Graph{F,W}) where {F,W,C}
#     return linear_combination(Unity, g1, F(c), F(-1))
# end
# function Base.:-(g1::Graph{F,W},c::C) where {F,W,C}
#     return linear_combination(g1, Unity, F(1), F(-c))
# end


"""
    function Base.:+(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}

    Returns a graph `g1 + g2` representing the addition of `g2` with `g1`.

# Arguments:
- `g1`  first computational graph
- `g2`  second computational graph
"""
function Base.:+(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
    return linear_combination(g1, g2, F(1), F(1))
end

"""
    function Base.:-(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}

    Returns a graph `g1 - g2` representing the subtraction of `g2` from `g1`.

# Arguments:
- `g1`  first computational graph
- `g2`  second computational graph
"""
function Base.:-(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
    return linear_combination(g1, g2, F(1), F(-1))
end


"""
    function multi_product(g1::Graph{F,W}, g2::Graph{F,W}, c1::C, c2::C) where {F,W,C}

    Returns a graph representing the multi product `c1*g1 * c2*g2`.

# Arguments:
- `g1`  first computational graph
- `g2`  second computational graph
- `c1`  first scalar multiple
- `c2`  second scalar multiple
"""
function multi_product(g1::Graph{F,W}, g2::Graph{F,W}, c1::C=1, c2::C=1) where {F,W,C}
    g = Graph([g1, g2]; subgraph_factors=[F(c1), F(c2)], operator=Prod(), ftype=F, wtype=W)
    # Convert multiplicative links to in-place form
    if unary_istrivial(g1.operator) && onechild(g1)
        g.subgraph_factors[1] *= g1.subgraph_factors[1]
        g.subgraphs[1] = g1.subgraphs[1]
    end
    if unary_istrivial(g2.operator) && onechild(g2)
        g.subgraph_factors[2] *= g2.subgraph_factors[1]
        g.subgraphs[2] = g2.subgraphs[1]
    end
    return g
end

"""
    function multi_product(graphs::Vector{Graph{F,W}}, constants::Vector{C}) where {F,W,C}

    Given a vector 𝐠 of graphs and an equally-sized vector 𝐜 of constants, returns a new
    graph representing the linear combination (𝐜 ⋅ 𝐠).

# Arguments:
- `graphs`  vector of computational graphs
- `constants`  vector of scalar multiples
"""
function multi_product(graphs::Vector{Graph{F,W}}, constants::Vector{C}=ones(C, length(graphs.subgraphs))) where {F,W,C}
    g = Graph(graphs; subgraph_factors=constants, operator=Prod(), ftype=F, wtype=W)
    # Convert multiplicative links to in-place form
    for (i, sub_g) in enumerate(g.subgraphs)
        if unary_istrivial(sub_g.operator) && onechild(sub_g)
            g.subgraph_factors[i] *= sub_g.subgraph_factors[1]
            g.subgraphs[i] = sub_g.subgraphs[1]
        end
    end
    return g
end


"""
    function Base.:*(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}

    Returns a graph `g1 * g2` representing the graph product between `g1` and `g2`.

# Arguments:
- `g1`  first computational graph
- `g2`  second computational graph
"""
function Base.:*(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
    return multi_product(g1, g2)
end
