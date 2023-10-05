"""
    mutable struct TaylorSeries{F,W,N}
    
    A representation of a computational graph, e.g., an expression tree, with type stable node data.

# Members:
- `id::Int`  the unique hash id to identify the diagram
- `name::Symbol`  name of the diagram
- `expansion::Dict{Dict{Int,Int},T}`  The taylor expansion coefficients. The key Dict{Int,Int} labels the order with respect to each variables. 
- `variables::Set{V}`  Variables of the taylor series. Each variable must have an unique id. 

"""
mutable struct TaylorSeries{T,V}
    id::Int
    name::String # "" by default
    expansion::Dict{Dict{Int,Int},T}
    variables::Set{V}

    """
        function TaylorSeries(T::DataType=Float64, name="", expansion=Dict{Dict{Int,Int},T}(), variables=Set{V}())
            Create a TaylorSeries based on given expansion and variables.
    """
    function TaylorSeries(T::DataType=Float64, V::DataType=Float64, name="", expansion=Dict{Dict{Int,Int},T}(), variables=Set{V}())
        return new{T,V}(uid(), name, expansion, variables)
    end
end

"""
    function  identityseries(g::TaylorSeries{T,V}, value::T) where {T,V}
        For a given variable g,  create a TaylorSeries equal to f (g) = g. 
        Assign zero order taylor coefficient with given value, and first order coefficient as one.

# Arguments:
- 'g'  Variable
- ' value'  Zero order taylor coefficient of g.    
"""

function identityseries(g::V, value::T) where {T,V}
    gnew = TaylorSeries(T, V)
    push!(gnew.variables, g)
    gnew.expansion[Dict(g.id => 0)] = value
    gnew.expansion[Dict(g.id => 1)] = one(T)
    return gnew
end

"""
    function Base.:*(g1::TaylorSeries{T}, c2::Number) where {T}

    Returns a TaylorSeries representing the scalar multiplication `g1*c2`.

# Arguments:
- `g1`  TaylorSeries
- `c2`  scalar multiple
"""
function Base.:*(g1::TaylorSeries{T,V}, c2::Number) where {T,V}
    g = TaylorSeries(T, V)
    for (key, value) in g1.expansion
        g.expansion[key] = c2 * value
    end
    g.variables = g1.variables
    return g
end

"""
    function Base.:*(c1::Number, g2::TaylorSeries{T}) where {T}

    Returns a TaylorSeries representing the scalar multiplication `g2*c1`.

# Arguments:
- `g2`  TaylorSeries
- `c1`  scalar multiple
"""
function Base.:*(c1::Number, g2::TaylorSeries{T,V}) where {T,V}
    g = TaylorSeries(T, V)
    for (key, value) in g2.expansion
        g.expansion[key] = c1 * value
    end
    g.variables = g2.variables
    return g
end

"""
    function Base.:+(g1::TaylorSeries{T,V}, g2::TaylorSeries{T,V}) where {T,V}

    Returns a taylor series `g1 + g2` representing the addition of `g2` with `g1`.

# Arguments:
- `g1`  First taylor series
- `g2`  Second taylor series
"""
function Base.:+(g1::TaylorSeries{T,V}, g2::TaylorSeries{T,V}) where {T,V}
    g = TaylorSeries(T, V)
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
    function Base.:-(g1::TaylorSeries{T,V}, g2::TaylorSeries{T,V}) where {T,V}

    Returns a taylor series `g1 - g2` representing the difference of `g2` with `g1`.

# Arguments:
- `g1`  First taylor series
- `g2`  Second taylor series
"""
function Base.:-(g1::TaylorSeries{T}, g2::TaylorSeries{T}) where {T}
    return g1 + (-1 * g2)
end

"""
    function merge_order(o1::Dict{Int,Int}, o2::Dict{Int,Int})

    For two  dictionary representing  order of  two taylor coefficients c1 and c2, generate the order of  c1*c2.

    # Arguments:
    - `o1`  First order label
    - `o2`  Second order label
"""
function merge_order(o1::Dict{Int,Int}, o2::Dict{Int,Int})
    o = copy(o1)
    combination = 1.0
    for (id, order) in o2
        if haskey(o, id)
            o[id] += order
            combination *= binomial(o[id], order)
        else
            o[id] = order
        end
    end
    return o, combination
end

function taylor_combinatorial(o::Dict{Int,Int})
    coeff = 1.0
    for (id, order) in o
        coeff *= factorial(order)
    end
    return 1.0 / coeff
end
"""
    function Base.:*(g1::TaylorSeries{T,V}, g2::TaylorSeries{T,V}) where {T,V}

    Returns a taylor series `g1 * g2` representing the product of `g2` with `g1`.

# Arguments:
- `g1`  First taylor series
- `g2`  Second taylor series
"""
function Base.:*(g1::TaylorSeries{T,V}, g2::TaylorSeries{T,V}) where {T,V}
    g = TaylorSeries(T, V)
    g.variables = union(g1.variables, g2.variables)
    for (key1, value1) in g1.expansion
        for (key2, value2) in g2.expansion
            key, combination_coeff = merge_order(key1, key2)
            if haskey(g.expansion, key)
                g.expansion[key] += combination_coeff * value1 * value2
            else
                g.expansion[key] = combination_coeff * value1 * value2
            end
        end
    end
    return g
end
