"""
    ParamsTaylor

DataType holding the current parameters for `TaylorSeries`. This part of code 

**Fields:**

- `order            :: Int`  Order (degree) of the polynomials
- `num_vars         :: Int`  Number of variables
- `variable_names   :: Vector{String}`  Names of the variables
- `variable_symbols :: Vector{Symbol}`  Symbols of the variables

These parameters can be changed using [`set_variables`](@ref)
"""
mutable struct ParamsTaylor
    order::Int
    num_vars::Int
    variable_names::Vector{String}
    variable_symbols::Vector{Symbol}
end


ParamsTaylor(order, num_vars, variable_names) = ParamsTaylor(order, num_vars, variable_names, Symbol.(variable_names))

const _params_Taylor_ = ParamsTaylor(6, 2, ["x₁", "x₂"])

## Utilities to get the maximum order, number of variables, their names and symbols
get_order() = _params_Taylor_.order
get_numvars() = _params_Taylor_.num_vars
get_variable_names() = _params_Taylor_.variable_names
get_variable_symbols() = _params_Taylor_.variable_symbols
function lookupvar(s::Symbol)
    ind = findfirst(x -> x == s, _params_Taylor_.variable_symbols)
    isa(ind, Nothing) && return 0
    return ind
end


"""
    set_variables([T::Type], names::String; [order=get_order(), numvars=-1])

Return a `Taylor{T}` vector with each entry representing an
independent variable. `names` defines the output for each variable
(separated by a space). The default type `T` is `Float64`,
and the default for `order` is the one defined globally.
Changing the `order` or `numvars` resets the hash_tables.

If `numvars` is not specified, it is inferred from `names`. If only
one variable name is defined and `numvars>1`, it uses this name with
subscripts for the different variables.

```julia
julia> set_variables(Int, "x y z", order=4)
3-element Array{TaylorSeries.Taylor{Int},1}:
  1 x + 𝒪(‖x‖⁵)
  1 y + 𝒪(‖x‖⁵)
  1 z + 𝒪(‖x‖⁵)

julia> set_variables("α", numvars=2)
2-element Array{TaylorSeries.Taylor{Float64},1}:
  1.0 α₁ + 𝒪(‖x‖⁵)
  1.0 α₂ + 𝒪(‖x‖⁵)

julia> set_variables("x", order=6, numvars=2)
2-element Array{TaylorSeries.Taylor{Float64},1}:
  1.0 x₁ + 𝒪(‖x‖⁷)
  1.0 x₂ + 𝒪(‖x‖⁷)
```
"""
function set_variables(::Type{R}, names::Vector{T}; order=get_order()) where
{R,T<:AbstractString}

    order ≥ 1 || error("Order must be at least 1")

    num_vars = length(names)
    num_vars ≥ 1 || error("Number of variables must be at least 1")

    _params_Taylor_.variable_names = names
    _params_Taylor_.variable_symbols = Symbol.(names)
    if !(order == get_order() && num_vars == get_numvars())
        # if these are unchanged, no need to regenerate tables

        _params_Taylor_.order = order
        _params_Taylor_.num_vars = num_vars
    end
    # return a list of the new variables
    return TaylorSeries{R}[TaylorSeries(R, i) for i in 1:get_numvars()]
end


set_variables(::Type{R}, symbs::Vector{T}; order=get_order()) where
{R,T<:Symbol} = set_variables(R, string.(symbs), order=order)

set_variables(names::Vector{T}; order=get_order()) where {T<:AbstractString} =
    set_variables(Float64, names, order=order)
set_variables(symbs::Vector{T}; order=get_order()) where {T<:Symbol} =
    set_variables(Float64, symbs, order=order)

function set_variables(::Type{R}, names::T; order=get_order(), numvars=-1) where
{R,T<:AbstractString}

    variable_names = split(names)

    if length(variable_names) == 1 && numvars ≥ 1
        name = variable_names[1]
        variable_names = [string(name, subscriptify(i)) for i in 1:numvars]
    end

    set_variables(R, variable_names, order=order)
end
set_variables(::Type{R}, symbs::Symbol; order=get_order(), numvars=-1) where {R} =
    set_variables(R, string(symbs), order=order, numvars=numvars)

set_variables(names::T; order=get_order(), numvars=-1) where {T<:AbstractString} =
    set_variables(Float64, names, order=order, numvars=numvars)
set_variables(symbs::Symbol; order=get_order(), numvars=-1) =
    set_variables(Float64, string(symbs), order=order, numvars=numvars)