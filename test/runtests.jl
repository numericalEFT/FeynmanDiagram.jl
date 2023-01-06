using FeynmanDiagram
using Lehmann
using Test, LinearAlgebra, Random, StaticArrays, Printf, Parameters, Documenter
using AbstractTrees

"""
Skip a testset

Use `@testset_skip` to replace `@testset` for some tests which should be skipped.

Usage
-----
Replace `@testset` with `@testset "reason"` where `"reason"` is a string saying why the
test should be skipped (which should come before the description string, if that is
present).
"""
macro testset_skip(args...)
    isempty(args) && error("No arguments to @testset_skip")
    length(args) < 2 && error("First argument to @testset_skip giving reason for "
                              *
                              "skipping is required")
    skip_reason = args[1]
    desc, testsettype, options = Test.parse_testset_args(args[2:end-1])
    ex = quote
        # record the reason for the skip in the description
        # and mark the tests as broken, but don't run tests
        local ts = Test.DefaultTestSet(string($desc, " - ", $skip_reason))
        push!(ts.results, Test.Broken(:skipped, "skipped tests"))
        local ret = Test.finish(ts)
        ret
    end
    return ex
end

# @testset "doctest" begin
#     DocMeta.setdocmeta!(FeynmanDiagram, :DocTestSetup, :(using FeynmanDiagram); recursive = true)
#     doctest(FeynmanDiagram)
# end
if isempty(ARGS)
    include("common.jl")
    include("quantum_operator.jl")
    include("computational_graph.jl")
    include("diagram_tree.jl")
    include("expression_tree.jl")
    include("parquet_builder.jl")
else
    include(ARGS[1])
end

