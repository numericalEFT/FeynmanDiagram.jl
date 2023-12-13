"""
    ParamsTaylor

    DataType holding the current parameters for `TaylorSeries`. 
    This part of code is adopted from TaylorSeries.jl (https://github.com/JuliaDiff/TaylorSeries.jl)

**Fields:**

- `orders            :: Int`  Orders (degree) of the polynomials
- `num_vars         :: Int`  Number of variables
- `variable_names   :: Vector{String}`  Names of the variables
- `variable_symbols :: Vector{Symbol}`  Symbols of the variables

These parameters can be changed using [`set_variables`](@ref)
"""
mutable struct ParamsTaylor
    orders::Vector{Int}
    num_vars::Int
    variable_names::Vector{String}
    variable_symbols::Vector{Symbol}
end


ParamsTaylor(orders, num_vars, variable_names) = ParamsTaylor(orders, num_vars, variable_names, Symbol.(variable_names))

const _params_Taylor_ = ParamsTaylor([2, 2], 2, ["x₁", "x₂"])

## Utilities to get the maximum order, number of variables, their names and symbols
get_orders() = _params_Taylor_.orders
get_orders(idx::Int) = _params_Taylor_.orders[idx]
get_numvars() = _params_Taylor_.num_vars
get_variable_names() = _params_Taylor_.variable_names
get_variable_symbols() = _params_Taylor_.variable_symbols
function lookupvar(s::Symbol)
    ind = findfirst(x -> x == s, _params_Taylor_.variable_symbols)
    isa(ind, Nothing) && return 0
    return ind
end


"""
    set_variables([T::Type], names::String; [orders=get_orders(), numvars=-1])

Return a `TaylorSeries{T}` vector with each entry representing an
independent variable. `names` defines the output for each variable
(separated by a space). The default type `T` is `Float64`,
and the default for `orders` is the one defined globally.

If `numvars` is not specified, it is inferred from `names`. If only
one variable name is defined and `numvars>1`, it uses this name with
subscripts for the different variables.

```julia
julia> set_variables(Int, "x y z", orders=[4,4,4])
3-element Array{TaylorSeries.Taylor{Int},1}:
  1 x + 𝒪(x⁵y⁵z⁵)
  1 y + 𝒪(x⁵y⁵z⁵)
  1 z + 𝒪(x⁵y⁵z⁵)
```
"""
function set_variables(::Type{R}, names::Vector{T}; orders=get_orders()) where
{R,T<:AbstractString}
    # for o in orders
    #     o ≥ 1 || error("Order must be at least 1")
    # end
    num_vars = length(names)
    num_vars ≥ 1 || error("Number of variables must be at least 1")
    @assert length(orders) == num_vars "Input orders should have same length as number of variables."
    _params_Taylor_.variable_names = names
    _params_Taylor_.variable_symbols = Symbol.(names)
    if !(orders == get_orders() && num_vars == get_numvars())
        # if these are unchanged, no need to regenerate tables

        _params_Taylor_.orders = orders
        _params_Taylor_.num_vars = num_vars
    end
    # return a list of the new variables
    return TaylorSeries{R}[TaylorSeries(R, i) for i in 1:get_numvars()]
end


set_variables(::Type{R}, symbs::Vector{T}; orders=get_orders()) where
{R,T<:Symbol} = set_variables(R, string.(symbs), orders=orders)

set_variables(names::Vector{T}; orders=get_orders()) where {T<:AbstractString} =
    set_variables(Float64, names, orders=orders)
set_variables(symbs::Vector{T}; orders=get_orders()) where {T<:Symbol} =
    set_variables(Float64, symbs, orders=orders)

function set_variables(::Type{R}, names::T; orders=get_orders(), numvars=-1) where
{R,T<:AbstractString}

    variable_names = split(names)

    if length(variable_names) == 1 && numvars ≥ 1
        name = variable_names[1]
        variable_names = [string(name, subscriptify(i)) for i in 1:numvars]
    end

    set_variables(R, variable_names, orders=orders)
end
set_variables(::Type{R}, symbs::Symbol; orders=get_orders(), numvars=-1) where {R} =
    set_variables(R, string(symbs), orders=orders, numvars=numvars)

set_variables(names::T; orders=get_orders(), numvars=-1) where {T<:AbstractString} =
    set_variables(Float64, names, orders=orders, numvars=numvars)
set_variables(symbs::Symbol; orders=get_orders(), numvars=-1) =
    set_variables(Float64, string(symbs), orders=orders, numvars=numvars)