module FrontEnds

import ..QuantumOperators as Op
import ..ComputationalGraphs as IR
import ..ComputationalGraphs: Graph
import ..ComputationalGraphs: _dtype

include("LabelProduct.jl")
export LabelProduct

include("parquet.jl")
using .Parquet
# export Parquet

include("readfile.jl")
export read_onediagram, read_diagrams


end