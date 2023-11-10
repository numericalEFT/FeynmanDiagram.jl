module FrontEnds

import ..ComputationalGraphs
using LinearAlgebra
import ..QuantumOperators as Op
import ..ComputationalGraphs as IR
import ..ComputationalGraphs: FeynmanGraph
# import ..ComputationalGraphs: Graph
import ..ComputationalGraphs: _dtype

using ..DiagTree


include("pool.jl")
export LoopPool

include("LabelProduct.jl")
export LabelProduct

include("parquet.jl")
using .Parquet
# export Parquet

include("diagtree.jl")

end