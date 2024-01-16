using FeynmanDiagram
using FeynmanDiagram: ComputationalGraphs as Graphs
using TensorOperations

@testset verbose = true "Tensorgraph" begin
    using FeynmanDiagram.ComputationalGraphs:
        TensorGraph, einsum
    g1 = TensorGraph([], dims=[5, 10, 15])
    g2 = TensorGraph([], dims=[15, 5, 20, 25])
    g3 = TensorGraph([], dims=[5, 10])

    g4 = einsum(g1, [1, 2, 3], g2, [3, 1, 4, 5], [2, 4, 5])
    @test g4.dims == [10, 20, 25]
    g5 = einsum(g1, [1, 2, 3], 2 * g2, [3, 1, 4, 5], [2, 3])
    @test g5.dims == [10, 15]
    g6 = einsum(g1, [1, 2, 3], g3, [1, 2], [2, 3])
    @test g6.dims == [10, 15]
    g7 = einsum((4 * g5 + 3 * g6), [2, 3], g2, [3, 1, 4, 5]) #if output axes is not specified, it is automatically derived following tensordot rule.
    @test g7.dims == [10, 5, 20, 25]

    #test from numpy.einsum
    a = collect(reshape(0:24, 5, 5))
    b = collect(0:4)
    ga = TensorGraph([], dims=collect(size(a)))
    gb = TensorGraph([], dims=collect(size(b)))

    @tensor result[i] := a[i, j] * b[j]
    @test collect(size(result)) == einsum(ga, [0, 1], gb, [1]).dims

    a = collect(reshape(0:59, 3, 4, 5))
    b = collect(reshape(0:23, 4, 3, 2))
    ga = TensorGraph([], dims=collect(size(a)))
    gb = TensorGraph([], dims=collect(size(b)))
    @tensor result[k, l] := a[i, j, k] * b[j, i, l]
    @test collect(size(result)) == einsum(ga, [0, 1, 2], gb, [1, 0, 3], [2, 3]).dims
end