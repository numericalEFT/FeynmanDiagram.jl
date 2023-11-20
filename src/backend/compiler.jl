module Compilers
using ..ComputationalGraphs
import ..ComputationalGraphs: id, name, set_name!, operator, subgraphs, subgraph_factors, factor

using ..AbstractTrees
using ..RuntimeGeneratedFunctions

RuntimeGeneratedFunctions.init(Compilers)

include("static.jl")
include("mindspore.jl")

end