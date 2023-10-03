"""
    mutable struct TaylorSeries{F,W,N}
    
    A representation of a computational graph, e.g., an expression tree, with type stable node data.

# Members:
- `id::Int`  the unique hash id to identify the diagram
- `name::Symbol`  name of the diagram
- `subgraphs::Vector{TaylorSeries{F,W,N}}`  vector of sub-diagrams 
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
mutable struct TaylorSeries{N,T}
    variable_number::Int
    expansion::Dict{SVector{N,Int},T}
    """
        function Graph(subgraphs::AbstractVector; name="", operator::AbstractOperator=Sum(),
            ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype)

    """
    function TaylorSeries(N::Int, T::DataType=Float64, expansion=Dict{SVector{N,Int},T}())
        return new{N,T}(N, expansion)
    end
end

# """
#     function constant_graph(factor=one(_dtype.factor))

#     Returns a graph that represents a constant equal to f, where f is the factor with default value 1.

# # Arguments:
# - `f`:  constant factor
# """
# function constant_graph(factor=one(_dtype.factor))
#     return Graph([]; operator=Constant(), factor=factor, ftype=_dtype.factor, wtype=_dtype.weight, weight=one(_dtype.weight))
# end


"""
    function Base.:*(g1::TaylorSeries{F,W,N}, c2::C) where {F,W,N,C}

    Returns a graph representing the scalar multiplication `g1*c2`.

# Arguments:
- `g1`  computational graph
- `c2`  scalar multiple
"""
function Base.:*(g1::TaylorSeries{N,T}, c2::Number) where {N,T}
    g = TaylorSeries(N, T)
    for (key, value) in g1.expansion
        g.expansion[key] = c2 * value
    end
    return g
end

"""
    function Base.:*(c1::C, g2::TaylorSeries{F,W,N}) where {F,W,N,C}

    Returns a graph representing the scalar multiplication `c1*g2`.

# Arguments:
- `c1`  scalar multiple
- `g2`  computational graph
"""
function Base.:*(c1::Number, g2::TaylorSeries{N,T}) where {N,T}
    g = TaylorSeries(N, T)
    for (key, value) in g2.expansion
        g.expansion[key] = c1 * value
    end
    return g
end


"""
    function Base.:+(g1::TaylorSeries{F,W,N}, g2::TaylorSeries{F,W,N}) where {F,W,N}

    Returns a graph `g1 + g2` representing the addition of `g2` with `g1`.

# Arguments:
- `g1`  first computational graph
- `g2`  second computational graph
"""
function Base.:+(g1::TaylorSeries{N,T}, g2::TaylorSeries{N,T}) where {N,T}
    g = TaylorSeries(N, T)
    for (key, value) in g1.expansion
        g.expansion[key] = value
    end
    for (key, value) in g2.expansion
        if haskey(g.expansion, key)
            g.expansion[key] += value
        else
            g.expansion[key] = value
        end
    end
    return g
end

"""
    function Base.:-(g1::TaylorSeries{F,W,N}, g2::TaylorSeries{F,W,N}) where {F,W,N}

    Returns a graph `g1 - g2` representing the subtraction of `g2` from `g1`.

# Arguments:
- `g1`  first computational graph
- `g2`  second computational graph
"""
function Base.:-(g1::TaylorSeries{N,T}, g2::TaylorSeries{N,T}) where {N,T}
    return g1 + (-1 * g2)
end


"""
    function Base.:*(g1::TaylorSeries{F,W,N}, g2::TaylorSeries{F,W,N}) where {F,W,N}

    Returns a graph `g1 * g2` representing the graph product between `g1` and `g2`.

# Arguments:
- `g1`  first computational graph
- `g2`  second computational graph
"""
function Base.:*(g1::TaylorSeries{N,T}, g2::TaylorSeries{N,T}) where {N,T}
    g = TaylorSeries(N, T)
    for (key1, value1) in g1.expansion
        for (key2, value2) in g2.expansion
            key = key1 + key2
            if haskey(g.expansion, key)
                g.expansion[key] += value1 * value2
            else
                g.expansion[key] = value1 * value2
            end
        end
    end
    return g
end
