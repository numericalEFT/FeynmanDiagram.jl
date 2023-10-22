"""
    function Base.:*(g1::TaylorSeries{T}, c2::Number) where {T}

    Returns a TaylorSeries representing the scalar multiplication `g1*c2`.

# Arguments:
- `g1`  TaylorSeries
- `c2`  scalar multiple
"""
function Base.:*(g1::TaylorSeries{T}, c2::Number) where {T}
    g = TaylorSeries{T}()
    for (order, coeff) in g1.coeffs
        g.coeffs[order] = c2 * coeff
    end
    return g
end

"""
    function Base.:*(c1::Number, g2::TaylorSeries{T}) where {T}

    Returns a TaylorSeries representing the scalar multiplication `g2*c1`.

# Arguments:
- `g2`  TaylorSeries
- `c1`  scalar multiple
"""
function Base.:*(c1::Number, g2::TaylorSeries{T}) where {T}
    g = TaylorSeries{T}()
    for (order, coeff) in g2.coeffs
        g.coeffs[order] = c1 * coeff
    end
    return g
end

"""
    function Base.:+(g1::TaylorSeries{T,V}, g2::TaylorSeries{T,V}) where {T,V}

    Returns a taylor series `g1 + g2` representing the addition of `g2` with `g1`.

# Arguments:
- `g1`  First taylor series
- `g2`  Second taylor series
"""
function Base.:+(g1::TaylorSeries{T}, g2::TaylorSeries{T}) where {T}
    g = TaylorSeries{T}()
    g.coeffs = copy(g1.coeffs)

    for (order, coeff) in g2.coeffs
        if haskey(g.coeffs, order)
            g.coeffs[order] += coeff
        else
            g.coeffs[order] = coeff
        end
    end
    return g
end


"""
    function Base.:+(g1::TaylorSeries{T,V}, g2::TaylorSeries{T,V}) where {T,V}

    Returns a taylor series `g1 + g2` representing the addition of `g2` with `g1`.

# Arguments:
- `g1`  First taylor series
- `g2`  Second taylor series
"""
function Base.:+(g1::TaylorSeries{T}, c::S) where {T,S<:Number}
    g = TaylorSeries{T}()
    g.coeffs = copy(g1.coeffs)

    zero_order = zeros(Int, get_numvars())
    if haskey(g.coeffs, zero_order)
        g.coeffs[order] += T(c)
    else
        g.coeffs[order] = T(c)
    end

    return g
end

function Base.:+(c::S, g1::TaylorSeries{T}) where {S<:Number,T}
    g = TaylorSeries{T}()
    g.coeffs = copy(g1.coeffs)

    zero_order = zeros(Int, get_numvars())
    if haskey(g.coeffs, zero_order)
        g.coeffs[zero_order] += T(c)
    else
        g.coeffs[zero_order] = T(c)
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
function Base.:-(c::S, g2::TaylorSeries{T}) where {T,S<:Number}
    return c + (-1 * g2)
end
function Base.:-(g1::TaylorSeries{T}, c::S) where {T,S<:Number}
    return g1 + (-1 * c)
end


function taylor_binomial(o1::Array{Int,1}, o2::Array{Int,1})
    @assert length(o1) == length(o2)
    result = 1
    for i in eachindex(o1)
        order = o1[i] + o2[i]
        if !iszero(order)
            result *= binomial(order, o1[i])
        end
    end
    return result
end

function taylor_factorial(o::Array{Int,1})
    result = 1
    for i in eachindex(o)
        result *= factorial(o[i])
    end
    return result
end
"""
    function Base.:*(g1::TaylorSeries{T,V}, g2::TaylorSeries{T,V}) where {T,V}

    Returns a taylor series `g1 * g2` representing the product of `g2` with `g1`.

# Arguments:
- `g1`  First taylor series
- `g2`  Second taylor series
"""
function Base.:*(g1::TaylorSeries{T}, g2::TaylorSeries{T}) where {T}
    g = TaylorSeries{T}()
    for (order1, coeff1) in g1.coeffs
        for (order2, coeff2) in g2.coeffs
            order = order1 + order2
            if sum(order) <= _params_Taylor_.order
                #combination_coeff = taylor_binomial(order1, order2)
                if haskey(g.coeffs, order)
                    #g.coeffs[order] += combination_coeff * coeff1 * coeff2
                    g.coeffs[order] += coeff1 * coeff2
                else
                    #g.coeffs[order] = combination_coeff * coeff1 * coeff2
                    g.coeffs[order] = coeff1 * coeff2
                end
            end
        end
    end
    return g
end

# function findidx(a::TaylorSeries, o::Array{Int,1})
#     @assert length(o) == get_numvars()
#     return findfirst(isequal(o), a.order)
# end

function getcoeff(g::TaylorSeries, order::Array{Int,1})
    if haskey(g.coeffs, order)
        return g.coeffs[order]
        #return 1 / taylor_factorial(order) * g.coeffs[order]
    else
        return nothing
    end
end

function getderivative(g::TaylorSeries, order::Array{Int,1})
    if haskey(g.coeffs, order)
        #return g.coeffs[order]
        return taylor_factorial(order) * g.coeffs[order]
    else
        return nothing
    end
end

function Base.one(x::TaylorSeries{T}) where {T}
    g = TaylorSeries{T}()
    g.coeffs[zeros(Int, get_numvars)] = one(T)
    return g
end
function Base.:^(x::TaylorSeries, p::Integer)
    p == 1 && return copy(x)
    p == 0 && return one(x)
    p == 2 && return square(x)
    p < 0 && throw(DomainError())
    return power_by_squaring(x, p)
end

function square(x::TaylorSeries)
    return x * x
end
function power_by_squaring(x::TaylorSeries, p::Integer)
    p == 1 && return copy(x)
    p == 0 && return one(x)
    p == 2 && return square(x)
    t = trailing_zeros(p) + 1
    p >>= t

    while (t -= 1) > 0
        x = square(x)
    end

    y = x
    while p > 0
        t = trailing_zeros(p) + 1
        p >>= t
        while (t -= 1) ≥ 0
            x = square(x)
        end
        y *= x
    end

    return y
end