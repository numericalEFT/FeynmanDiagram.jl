module Compilers
using ..ComputationalGraphs
import ..ComputationalGraphs: FeynmanGraph

using ..AbstractTrees
using ..RuntimeGeneratedFunctions

RuntimeGeneratedFunctions.init(Compilers)

include("static.jl")

end