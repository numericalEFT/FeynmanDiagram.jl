module FeynmanDiagram
using Random, LinearAlgebra

# include("variable.jl")
# export Var

# include("GWKT/GWKT.jl")
# export GWKT

include("diagramTree/DiagTree.jl")
export DiagTree

include("diagramBuilder/builder.jl")
export Builder
# Write your package code here.

end
