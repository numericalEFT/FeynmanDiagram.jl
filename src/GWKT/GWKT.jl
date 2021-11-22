module GWKT
# export PropagatorKT, Weight, addChild
include("parquet/parquet.jl")
export Parquet

include("diagtree.jl")
export DiagTree

end