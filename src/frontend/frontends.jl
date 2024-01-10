module FrontEnds

import ..ComputationalGraphs
using LinearAlgebra
import ..QuantumOperators as Op
import ..ComputationalGraphs as IR
import ..ComputationalGraphs: FeynmanGraph
import ..ComputationalGraphs: Graph
import ..ComputationalGraphs: _dtype

include("pool.jl")
export LoopPool

include("LabelProduct.jl")
export LabelProduct

# include("parquet/parquet.jl")
# using .Parquet
# export Parquet

end