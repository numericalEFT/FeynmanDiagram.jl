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
    function Base.:+(g1::TaylorSeries{T}, g2::TaylorSeries{T}) where {T}

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
    function Base.:+(g1::TaylorSeries{T}, c::S) where {T,S<:Number}

    Returns a taylor series `g1 + c` representing the addition of constant `c` with `g1`.
    
# Arguments:
- `g1`  Taylor series
- `c`  Constant
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


"""
    function Base.:+(c::S, g1::TaylorSeries{T}) where {S<:Number,T}

    Returns a taylor series `g1 + c` representing the addition of constant `c` with `g1`.
    
# Arguments:
- `g1`  Taylor series
- `c`  Constant
"""
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

"""
    function taylor_binomial(o1::Vector{Int}, o2::Vector{Int})

    Return the taylor binomial prefactor when product two high-order derivatives with order o1 and o2.

    # Arguments:
    - `o1`  Order of first derivative
    - `o2`  Order of second derivative
"""
function taylor_binomial(o1::Vector{Int}, o2::Vector{Int})
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


"""
    function taylor_factorial(o::Vector{Int})

    Return the taylor factorial prefactor with order o.

    # Arguments:
    - `o`  Order of the taylor coefficient
"""
function taylor_factorial(o::Vector{Int})
    result = 1
    for i in eachindex(o)
        result *= factorial(o[i])
    end
    return result
end

"""
    function Base.:*(g1::TaylorSeries{T}, g2::TaylorSeries{T}) where {T}

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
            if all(order .<= get_orders())
                #combination_coeff = taylor_binomial(order1, order2)
                if haskey(g.coeffs, order)
                    #g.coeffs[order] += combination_coeff * coeff1 * coeff2
                    # println(coeff1, coeff1.orders)
                    # println(coeff2, coeff2.orders)
                    # println(g.coeffs[order], g.coeffs[order].orders)
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

# """
#     function Base.:*(g1::TaylorSeries{T}, g2::TaylorSeries{T}) where {T}

#     Returns a taylor series `g1 * g2` representing the product of `g2` with `g1`.

# # Arguments:
# - `g1`  First taylor series
# - `g2`  Second taylor series
# """
# function multi_product(glist::Vector{TaylorSeries{T}}) where {T}
#     result = TaylorSeries{T}()
#     idxtuple = (keys(glist.coeffs) for g in glist)
#     for idx in collect(Iterators.product(idxtuple...))
#         orderidx = collect(orderidx)
#         totalorder = sum(glist[i].coeffs[orderidx[i]] for i in eachindex(glist))
#         if all(totalorder .<= get_orders())
#             #combination_coeff = taylor_binomial(order1, order2)
#             if haskey(g.coeffs, order)
#                 #g.coeffs[order] += combination_coeff * coeff1 * coeff2
#                 g.coeffs[order] += coeff1 * coeff2
#             else
#                 #g.coeffs[order] = combination_coeff * coeff1 * coeff2
#                 g.coeffs[order] = coeff1 * coeff2
#             end

#         end
#     end
#     return g
# end

"""
    function getcoeff(g::TaylorSeries, order::Vector{Int})

    Return the taylor coefficients with given order in taylor series g.

# Arguments:
- `g`  Taylor series
- `order`  Order of target coefficients
"""
function getcoeff(g::TaylorSeries, order::Vector{Int})
    if haskey(g.coeffs, order)
        return g.coeffs[order]
    else
        return nothing
    end
end

"""
    function getderivative(g::TaylorSeries, order::Vector{Int})

    Return the derivative with given order in taylor series g.

# Arguments:
- `g`  Taylor series
- `order`  Order of derivative
"""
function getderivative(g::TaylorSeries, order::Vector{Int})
    if haskey(g.coeffs, order)
        return taylor_factorial(order) * g.coeffs[order]
    else
        return nothing
    end
end


"""
    function  Base.one(g::TaylorSeries{T}) where {T}

    Return a constant one for a given taylor series.

# Arguments:
- `g`  Taylor series
"""
function Base.one(g::TaylorSeries{T}) where {T}
    unity = TaylorSeries{T}()
    unity.coeffs[zeros(Int, get_numvars)] = one(T)
    return unity
end


"""
    function Base.:^(x::TaylorSeries, p::Integer)

    Return the power of taylor series x^p, where p is an integer. 

# Arguments:
- `x`  Taylor series
- 'p' Power index
"""
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

# power_by_squaring; slightly modified from base/intfuncs.jl
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