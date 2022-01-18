module DiagTreeNew
using AbstractTrees
using Printf, PyCall, DataFrames

const ADD, MUL = 1, 2
export ADD, MUL

# include("interface.jl")
include("diagram.jl")

# struct Cache and struct Pool
include("pool.jl")

# struct Propagator, Node and Diagram
# include("object.jl")

# # IO operations
# include("io.jl")

# # diagram evaluation
# include("eval.jl")

end