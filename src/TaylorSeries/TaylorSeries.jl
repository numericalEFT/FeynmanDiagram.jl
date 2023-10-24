"""
    TaylorSeries

A Julia package for Taylor expansions in one or more independent variables.

The basic constructors is [`TaylorSeries`](@ref).

"""
module Taylor

using ..ComputationalGraphs
#using Markdown


#export show_params_TaylorN, show_monomials, displayBigO, use_show_default,    

include("parameter.jl")
include("constructors.jl")
include("print.jl")
include("arithmetic.jl")
export TaylorSeries

export get_order, get_numvars,
    set_variables, get_variables,
    get_variable_names, get_variable_symbols,
    displayBigO, use_show_default,
    getcoeff, taylor_factorial

end # module