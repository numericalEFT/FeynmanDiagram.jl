module ExprTree
using AbstractTrees, LinearAlgebra
using ..DiagTree

using Printf, PyCall

const ADD, MUL = 1, 2
export ADD, MUL

include("common.jl")

# struct Cache and struct Pool
include("pool.jl")

# struct Propagator, Node and Diagram
include("object.jl")
export Component

# IO operations
include("io.jl")

# diagram evaluation
include("eval.jl")

include("compile.jl")

end