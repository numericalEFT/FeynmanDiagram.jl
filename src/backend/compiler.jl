module Compilers
using ..ComputationalGraphs
import ..ComputationalGraphs: Graph

using ..AbstractTrees
using ..RuntimeGeneratedFunctions

RuntimeGeneratedFunctions.init(Compilers)

include("static.jl")

end