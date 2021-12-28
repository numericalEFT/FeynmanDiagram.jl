module Builder
using Parameters

include("parameter.jl")

using ..DiagTree

include("parquetAlg/parquet.jl")
export Parquet

end