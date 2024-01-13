using FeynmanDiagram
using FeynmanDiagram: ComputationalGraphs as Graphs


@testset verbose = true "Tensorgraph" begin
    using FeynmanDiagram.ComputationalGraphs:
        TensorGraph, einsum
    g1 = TensorGraph([], shape=[5, 10, 15])
    g2 = TensorGraph([], shape=[15, 5, 20, 25])
    g3 = TensorGraph([], shape=[5, 10])

    g4 = einsum(g1, [1, 2, 3], g2, [3, 1, 4, 5], [2, 4, 5])
    @test g4.shape == [10, 20, 25]
    g5 = einsum(g1, [1, 2, 3], g2, [3, 1, 4, 5], [2, 3], 1, -1)
    @test g5.shape == [10, 15]
    g6 = einsum(g1, [1, 2, 3], g3, [1, 2], [2, 3], 2, -2)
    @test g6.shape == [10, 15]
    g7 = einsum((4 * g5 + 3 * g6), [2, 3], g2, [3, 1, 4, 5], 1, -1) #if output axes is not specified, it is automatically derived following tensordot rule.
    @test g7.shape == [10, 5, 20, 25]
end