module ExprTree
using AbstractTrees, LinearAlgebra, StaticArrays
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
export ExpressionTree

# IO operations
include("io.jl")

# diagram evaluation
include("eval.jl")

include("build.jl")

end