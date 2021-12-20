module ExpressionTree
using Random, LinearAlgebra

include("variable.jl")
export Var

include("GWKT/GWKT.jl")
export GWKT

include("diagram_tree.jl")
export DiagTree

# Write your package code here.

end
