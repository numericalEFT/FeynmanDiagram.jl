"""
    function Base.:*(g1::TaylorSeries{T}, c2::Number) where {T}

    Returns a TaylorSeries representing the scalar multiplication `g1*c2`.

# Arguments:
- `g1`  TaylorSeries
- `c2`  scalar multiple
"""
function Base.:*(g1::TaylorSeries{T}, c2::Number) where {T}
    g = TaylorSeries{T}()
    g.coeffs = c2 * g1.coeffs
    g.order = copy(g1.order)
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
    g.coeffs = c1 * g2.coeffs
    g.order = copy(g2.order)

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
    g.order = copy(g1.order)

    for (i, coeff) in enumerate(g2.coeffs)
        idx = findfirst(isequal(g2.order[i]), g.order)
        if isnothing(idx)
            push!(g.order, g2.order[i])
            push!(g.coeffs, coeff)
        else
            g.coeffs[idx] += coeff
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
function Base.:+(g1::TaylorSeries{T}, c::S) where {T,S}
    g = TaylorSeries{T}()
    g.coeffs = copy(g1.coeffs)
    g.order = copy(g1.order)
    zero_order = zeros(Int, get_numvars())
    idx = findfirst(isequal(zero_order), g.order)
    if isnothing(idx)
        push!(g.order, zero_order)
        push!(g.coeffs, T(c))
    else
        g.coeffs[idx] += T(c)
    end

    return g
end

function Base.:+(c::S, g1::TaylorSeries{T}) where {S,T}
    g = TaylorSeries{T}()
    g.coeffs = copy(g1.coeffs)
    g.order = copy(g1.order)
    zero_order = zeros(Int, get_numvars())
    idx = findfirst(isequal(zero_order), g.order)
    if isnothing(idx)
        push!(g.order, zero_order)
        push!(g.coeffs, T(c))
    else
        g.coeffs[idx] += T(c)
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
function Base.:-(c::T, g2::TaylorSeries{T}) where {T}
    return c + (-1 * g2)
end
function Base.:-(g1::TaylorSeries{T}, c::T) where {T}
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
    for (i1, coeff1) in enumerate(g1.coeffs)
        for (i2, coeff2) in enumerate(g2.coeffs)
            order = g1.order[i1] + g2.order[i2]
            if sum(order) <= _params_Taylor_.order
                idx = findfirst(isequal(order), g.order)
                combination_coeff = taylor_binomial(g1.order[i1], g2.order[i2])
                if isnothing(idx)
                    push!(g.order, order)
                    push!(g.coeffs, combination_coeff * coeff1 * coeff2)
                else
                    g.coeffs[idx] += combination_coeff * coeff1 * coeff2
                end
            end
        end
    end
    return g
end

function findidx(a::TaylorSeries, o::Array{Int,1})
    @assert length(o) == get_numvars()
    return findfirst(isequal(o), a.order)
end

function getcoeff(a::TaylorSeries, o::Array{Int,1})
    idx = findidx(a, o)
    if isnothing(idx)
        return nothing
    else
        return a.coeffs[idx] / taylor_factorial(o)
    end
end