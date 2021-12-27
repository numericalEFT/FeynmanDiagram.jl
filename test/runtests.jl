using ExpressionTree
using Lehmann
using Test, LinearAlgebra, Random, StaticArrays, Printf

include("diagtree.jl")

include("parquet_eval.jl")
using ..ParquetEval
include("parquet.jl")

@testset "ExpressionTree.jl" begin
    # Write your tests here.
end
