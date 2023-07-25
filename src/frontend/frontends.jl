module FrontEnds

using LinearAlgebra
import ..QuantumOperators as Op
import ..ComputationalGraphs as IR
import ..ComputationalGraphs: Graph
import ..ComputationalGraphs: _dtype

include("pool.jl")
export LoopPool, CachePool

include("LabelProduct.jl")
export LabelProduct

include("parquet.jl")
using .Parquet
# export Parquet

end