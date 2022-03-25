using FeynmanDiagram
using Lehmann
using Test, LinearAlgebra, Random, StaticArrays, Printf, Parameters, Documenter
using AbstractTrees

# @testset "doctest" begin
#     DocMeta.setdocmeta!(FeynmanDiagram, :DocTestSetup, :(using FeynmanDiagram); recursive = true)
#     doctest(FeynmanDiagram)
# end

include("common.jl")
include("diagram_tree.jl")
include("expression_tree.jl")
include("parquet_builder.jl")

