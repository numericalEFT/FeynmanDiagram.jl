module ExprTreeNew
using AbstractTrees, LinearAlgebra
using ..DiagTree
# using Unrolled
# using InteractiveUtils

using Printf, PyCall

const ADD, MUL = 1, 2
export ADD, MUL

include("common.jl")

# struct Cache and struct Pool
include("pool.jl")

include("tree.jl")

# IO operations
include("io.jl")

# diagram evaluation
include("eval.jl")

end