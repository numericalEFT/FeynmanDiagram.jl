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
mutable struct TaylorSeries{T}
    id::Int
    expansion::Dict{Dict{Int,Int},T}
    variables::Set{TaylorSeries}
    """
        function Graph(subgraphs::AbstractVector; name="", operator::AbstractOperator=Sum(),
            ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype)

    """
    function TaylorSeries(T::DataType=Float64, expansion=Dict{Dict{Int,Int},T}(), variables=Set{TaylorSeries}())
        return new{T}(uid(), expansion, variables)
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

function identityseries(g::TaylorSeries{T}, value::T) where {T}
    gnew = TaylorSeries(T)
    push!(gnew.variables, g)
    gnew.expansion[Dict(g.id => 0)] = value
    gnew.expansion[Dict(g.id => 1)] = one(T)
    return gnew
end
"""
    function Base.:*(g1::TaylorSeries{F,W,N}, c2::C) where {F,W,N,C}

    Returns a graph representing the scalar multiplication `g1*c2`.

# Arguments:
- `g1`  computational graph
- `c2`  scalar multiple
"""
function Base.:*(g1::TaylorSeries{T}, c2::Number) where {T}
    g = TaylorSeries(T)
    for (key, value) in g1.expansion
        g.expansion[key] = c2 * value
    end
    g.variables = g1.variables
    return g
end

"""
    function Base.:*(c1::C, g2::TaylorSeries{F,W,N}) where {F,W,N,C}

    Returns a graph representing the scalar multiplication `c1*g2`.

# Arguments:
- `c1`  scalar multiple
- `g2`  computational graph
"""
function Base.:*(c1::Number, g2::TaylorSeries{T}) where {T}
    g = TaylorSeries(T)
    for (key, value) in g2.expansion
        g.expansion[key] = c1 * value
    end
    g.variables = g2.variables
    return g
end

"""
    function Base.:+(g1::TaylorSeries{F,W,N}, g2::TaylorSeries{F,W,N}) where {F,W,N}

    Returns a graph `g1 + g2` representing the addition of `g2` with `g1`.

# Arguments:
- `g1`  first computational graph
- `g2`  second computational graph
"""
function Base.:+(g1::TaylorSeries{T}, g2::TaylorSeries{T}) where {T}
    g = TaylorSeries(T)
    # for (key, value) in g1.expansion
    #     g.expansion[key] = value
    # end
    zero_order1 = Dict(var.id => 0 for var in g1.variables)
    zero_order2 = Dict(var.id => 0 for var in g2.variables)
    g.variables = union(g1.variables, g2.variables)
    for (key, value) in g1.expansion
        newkey = merge(key, zero_order2)
        g.expansion[newkey] = value
    end
    for (key, value) in g2.expansion
        newkey = merge(key, zero_order1)
        if haskey(g.expansion, newkey)
            g.expansion[newkey] += value
        else
            g.expansion[newkey] = value
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
function Base.:-(g1::TaylorSeries{T}, g2::TaylorSeries{T}) where {T}
    return g1 + (-1 * g2)
end


function merge_order(o1::Dict{Int,Int}, o2::Dict{Int,Int})
    o = copy(o1)
    for (id, order) in o2
        if haskey(o, id)
            o[id] += order
        else
            o[id] = order
        end
    end
    return o
end

"""
    function Base.:*(g1::TaylorSeries{F,W,N}, g2::TaylorSeries{F,W,N}) where {F,W,N}

    Returns a graph `g1 * g2` representing the graph product between `g1` and `g2`.

# Arguments:
- `g1`  first computational graph
- `g2`  second computational graph
"""
function Base.:*(g1::TaylorSeries{T}, g2::TaylorSeries{T}) where {T}
    g = TaylorSeries(T)
    g.variables = union(g1.variables, g2.variables)
    for (key1, value1) in g1.expansion
        for (key2, value2) in g2.expansion
            key = merge_order(key1, key2)
            if haskey(g.expansion, key)
                g.expansion[key] += value1 * value2
            else
                g.expansion[key] = value1 * value2
            end
        end
    end
    return g
end
