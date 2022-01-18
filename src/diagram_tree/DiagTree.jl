module DiagTreeNew
using AbstractTrees
using Printf, PyCall, DataFrames

# @enum Operator begin
#     Add
#     Mul
# end

# Base.length(r::Operator) = 1
# Base.iterate(r::Operator) = (r, nothing)
# function Base.iterate(r::Operator, ::Nothing) end

include("traits.jl")
include("diagram.jl")
export DiagramId
export Diagram
export add_subdiagram!
export toDataFrame

# struct Cache and struct Pool
# include("pool.jl")

# struct Propagator, Node and Diagram
# include("object.jl")

# # IO operations
# include("io.jl")

# # diagram evaluation
# include("eval.jl")

end