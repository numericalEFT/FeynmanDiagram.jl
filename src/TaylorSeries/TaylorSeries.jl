"""
    TaylorSeries

A Julia package for Taylor expansions in one or more independent variables.

The basic constructors are [`Taylor1`](@ref) and [`TaylorN`](@ref);
see also [`HomogeneousPolynomial`](@ref).

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
    # jacobian, hessian, jacobian!, hessian!,
    displayBigO, use_show_default,
    getcoeff, taylor_factorial
# function __init__()
#     @static if !isdefined(Base, :get_extension)
#         @require IntervalArithmetic = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253" begin
#             include("../ext/TaylorSeriesIAExt.jl")
#         end
#     end
# end
end # module