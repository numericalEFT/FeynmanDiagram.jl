using FeynmanDiagram
using Lehmann
using Test, LinearAlgebra, Random, StaticArrays, Printf, Parameters, Documenter
using AbstractTrees

# @testset "doctest" begin
#     DocMeta.setdocmeta!(FeynmanDiagram, :DocTestSetup, :(using FeynmanDiagram); recursive = true)
#     doctest(FeynmanDiagram)
# end
if isempty(ARGS)
    include("common.jl")
    include("diagram_tree.jl")
    include("expression_tree.jl")
    include("parquet_builder.jl")
else
    include(ARGS[1])
end

