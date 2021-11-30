module GWKT
# export PropagatorKT, Weight, addChild
include("diagtree.jl")
export DiagTree

include("parquet/parquet.jl")
export Parquet

include("manual/manual.jl")
export Manual

end