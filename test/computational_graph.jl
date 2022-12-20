@testset verbose = true "Graph" begin
    V = [𝑓⁺(1)𝑓⁻(2), 𝑓⁺(5)𝑓⁺(6)𝑓⁻(7)𝑓⁻(8), 𝑓⁺(3)𝑓⁻(4)]
    g1 = Graph(V, external=[1, 3])
    g2 = g1 * 2
    @testset "Graph equivalence" begin
        g1p = Graph(V, external=[1, 3])
        g2p = Graph(V, external=[1, 3], factor=2)
        # Test equivalence modulo fields id/factor
        @test isequiv(g1, g1p) == false
        @test isequiv(g1, g2p, :id) == false
        @test isequiv(g1, g2p, :factor) == false
        @test isequiv(g1, g1p, :id)
        @test isequiv(g1, g2p, :id, :factor)
        # Test inequivalence when subgraph lengths are different
        t = g1 + g1
        @test isequiv(t, g1, :id) == false
    end
    @testset "Scalar multiplication" begin
        @test vertices(g2) == vertices(g1)
        println(external(g2))
        println(external(g1))
        @test external(g2) == external(g1)
        @test g2.subgraph_factors == [2]
        @test g2.operator == ComputationalGraphs.Prod
        g2 = 2g1
        @test vertices(g2) == vertices(g1)
        @test external(g2) == external(g1)
        @test g2.subgraph_factors == [2]
        @test g2.operator == ComputationalGraphs.Prod
    end
    @testset "Graph addition" begin
        g3 = g1 + g2
        @test vertices(g3) == vertices(g1)
        @test external(g3) == external(g1)
        @test g3.factor == 1
        @test g3.subgraphs == [g1, g2]
        @test g3.subgraph_factors == [1, 1]
        @test isempty(g3.subgraphs[1].subgraph_factors)
        @test g3.subgraphs[2].subgraph_factors == [2]
        @test g3.operator == ComputationalGraphs.Sum
    end
    @testset "Graph subtraction" begin
        g4 = g1 - g2
        @test vertices(g4) == vertices(g1)
        @test external(g4) == external(g1)
        @test g4.factor == 1
        @test g4.subgraphs == [g1, g2]
        @test g4.subgraph_factors == [1, -1]
        @test isempty(g4.subgraphs[1].subgraph_factors)
        @test g4.subgraphs[2].subgraph_factors == [2]
        @test g4.subgraphs[2].subgraphs[1].factor == 1
        @test g4.operator == ComputationalGraphs.Sum
    end
    @testset "Linear combinations" begin
        # Binary form
        g5 = 3g1 + 5g2
        g5lc = ComputationalGraphs.linear_combination(g1, g2, 3, 5)
        @test g5.subgraph_factors == [1, 1]
        @test [g.subgraph_factors[1] for g in g5.subgraphs] == [3, 10]
        @test g5lc.subgraphs == [g1, g2]
        @test g5lc.subgraph_factors == [3, 5]
        # Requires optimization merge_prefactors on g5
        @test_broken isequiv(g5, g5lc, :id)
        # Vector form
        g6lc = ComputationalGraphs.linear_combination([g1, g2, g5, g2, g1], [3, 5, 7, 9, 11])
        @test g6lc.subgraphs == [g1, g2, g5, g2, g1]
        @test g6lc.subgraph_factors == [3, 5, 7, 9, 11]
    end
    @testset "Multiplicative chains" begin
        g6 = 7 * (5 * (3 * (2 * g1)))
        @test g6.subgraph_factors == [210]
        @test isempty(g6.subgraphs[1].subgraphs)
        @test isempty(g6.subgraphs[1].subgraph_factors)
        g7 = (((g1 * 2) * 3) * 5) * 7
        @test g7.subgraph_factors == [210]
        @test isempty(g7.subgraphs[1].subgraphs)
        @test isempty(g7.subgraphs[1].subgraph_factors)
    end
end

@testset "propagator" begin
    # g1 = propagator(𝑓⁺(1)𝑓⁻(2))
    g1 = propagator([𝑓⁺(1), 𝑓⁻(2)])
    @test g1.factor == 1
    @test g1.external == [1, 2]
    @test vertices(g1) == [𝑓⁺(1), 𝑓⁻(2)]
    @test external_with_ghost(g1) == external(g1) == [𝑓⁺(1), 𝑓⁻(2)]
    standardize_order!(g1)
    @test g1.factor == -1
    @test g1.external == [1, 2]
    @test external_with_ghost(g1) == external(g1) == [𝑓⁻(2), 𝑓⁺(1)]
    # @test vertices(g1) == [𝑓⁻(2)𝑓⁺(1)]

    g2 = propagator([𝑓⁺(1), 𝑓⁻(2), 𝑏⁺(1), 𝜙(1), 𝑓⁺(3), 𝑓⁻(1), 𝑓(1), 𝑏⁻(1), 𝜙(1)])
    @test vertices(g2) == external_with_ghost(g2) == external(g2) == [𝑓⁺(1), 𝑓⁻(2), 𝑏⁺(1), 𝜙(1), 𝑓⁺(3), 𝑓⁻(1), 𝑓(1), 𝑏⁻(1), 𝜙(1)]
    standardize_order!(g2)
    @test g2.factor == -1
    @test vertices(g2) == external_with_ghost(g2) == external(g2) == [𝑓⁻(1), 𝑏⁻(1), 𝜙(1), 𝑓⁻(2), 𝑓(1), 𝑓⁺(3), 𝜙(1), 𝑏⁺(1), 𝑓⁺(1)]
end

@testset verbose = true "feynman_diagram" begin
    @testset "Phi4" begin
        # phi theory 
        V1 = [𝜙(1)𝜙(1)𝜙(2)𝜙(2),]
        g1 = feynman_diagram(V1, [[1, 2], [3, 4]])    #vacuum diagram
        # g1 = feynman_diagram(V1, [1, 1, 2, 2])
        @test vertices(g1) == V1
        @test isempty(external(g1))
        @test g1.subgraph_factors == [1, 1]
    end
    @testset "Complex scalar field" begin
        #complex scalar field
        V2 = [𝑏⁺(1), 𝑏⁺(2)𝑏⁺(3)𝑏⁻(4)𝑏⁻(5), 𝑏⁺(6)𝑏⁺(7)𝑏⁻(8)𝑏⁻(9), 𝑏⁻(10)]
        # g2 = feynman_diagram(V2, [1, 2, 3, 4, 1, 4, 5, 2, 3, 5]; external=[1, 10])
        g2 = feynman_diagram(V2, [[1, 5], [2, 8], [3, 9], [4, 6], [7, 10]]; external=[1, 10])    # Green2
        @test vertices(g2) == V2
        @test external(g2) == [𝑏⁺(1), 𝑏⁻(10)]
        @test g2.subgraph_factors == [1, 1, 1, 1, 1]
    end
    @testset "Yukawa interaction" begin
        # Yukawa 
        V3 = [𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6)]
        # g3 = feynman_diagram(V3, [1, 2, 3, 2, 1, 3])
        g3 = feynman_diagram(V3, [[1, 5], [2, 4], [3, 6]])  #vacuum diagram
        @test vertices(g3) == V3
        @test isempty(external(g3))
        @test g3.factor == 1
        @test g3.subgraph_factors == [1, 1, 1]
        # @test internal_vertices(g3) == V3
        @test g3.subgraphs[1].factor == 1
        @test g3.subgraphs[1].vertices == external(g3.subgraphs[1]) == [𝑓⁺(1), 𝑓⁻(5)]
        standardize_order!(g3)
        @test g3.subgraphs[1].factor == -1
        @test external(g3.subgraphs[1]) == [𝑓⁻(5), 𝑓⁺(1)]

        V4 = [𝜙(13, true), 𝜙(14, true), 𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6), 𝑓⁺(7)𝑓⁻(8)𝜙(9), 𝑓⁺(10)𝑓⁻(11)𝜙(12),]
        g4 = feynman_diagram(V4, [[3, 10], [4, 12], [6, 13], [7, 9], [11, 14], [5, 1], [8, 2]], external=[5, 8]) # polarization diagram
        @test g4.factor == -1
        @test g4.subgraph_factors == [1, 1, 1, 1, 1]
        @test vertices(g4) == V4
        @test external(g4) == [𝜙(3), 𝜙(6)]
        @test external_with_ghost(g4) == [𝜙(13, true), 𝜙(14, true)]
        standardize_order!(g4)
        @test external(g4) == [𝜙(3), 𝜙(6)]
        # @test internal_vertices(g4) == V4[3:4]

        V5 = [𝑓⁻(1, true), 𝑓⁺(11, true), 𝜙(12, true), 𝑓⁺(2)𝑓⁻(3)𝜙(4), 𝑓⁺(5)𝑓⁻(6)𝜙(7), 𝑓⁺(8)𝑓⁻(9)𝜙(10)]
        g5 = feynman_diagram(V5, [[5, 7], [6, 12], [8, 10], [1, 4], [3, 9], [11, 2]], external=[4, 9, 11])  # vertex function
        @test g5.factor == 1
        @test g5.subgraph_factors == [1, 1, 1]
        @test vertices(g5) == V5
        @test external(g5) == [𝑓⁺(2), 𝜙(7), 𝑓⁻(9)]
        @test external_with_ghost(g5) == [𝑓⁻(1, true), 𝑓⁺(11, true), 𝜙(12, true)]
        # @test isempty(internal_vertices(g5))
        g5s = deepcopy(g5)
        standardize_order!(g5)
        @test g5.factor == -1
        @test external(g5) == [𝑓⁺(2), 𝜙(7), 𝑓⁻(9)]
        @test external_with_ghost(g5) == [𝑓⁺(11, true), 𝜙(12, true), 𝑓⁻(1, true)]
        # @test g5s == g5

        V5p = [𝑓⁺(11, true), 𝑓⁻(1, true), 𝜙(12, true), 𝑓⁺(2)𝑓⁻(3)𝜙(4), 𝑓⁺(5)𝑓⁻(6)𝜙(7), 𝑓⁺(8)𝑓⁻(9)𝜙(10)]
        g5p = feynman_diagram(V5p, [[5, 7], [6, 12], [8, 10], [2, 4], [3, 9], [11, 1]], external=[11, 9, 4])
        @test g5.factor ≈ g5p.factor    # reorder of external fake legs will not change the sign.
        @test g5p.subgraph_factors == [1, 1, 1]
        standardize_order!(g5p)
        @test external(g5p) == [𝑓⁻(9), 𝜙(7), 𝑓⁺(2)]
        @test g5p.factor ≈ g5.factor

        V6 = [𝑓⁻(8), 𝑓⁺(1), 𝑓⁺(2)𝑓⁻(3)𝜙(4), 𝑓⁺(5)𝑓⁻(6)𝜙(7)]
        g6 = feynman_diagram(V6, [[2, 4], [3, 7], [5, 8], [6, 1]], external=[1, 2])    # fermionic Green2
        @test g6.factor == -1
        @test g6.subgraph_factors == [1, 1, 1, 1]
        @test external(g6) == [𝑓⁻(8), 𝑓⁺(1)]
        standardize_order!(g6)
        @test g6.factor == 1
        @test vertices(g6)[1:2] == [𝑓⁺(1), 𝑓⁻(8)]
        @test Set(external(g6)) == Set([𝑓⁺(1), 𝑓⁻(8)])

        V7 = [𝑓⁻(7), 𝑓⁺ₑ(8), 𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6)]
        g7 = feynman_diagram(V7, [[3, 7], [5, 8], [6, 1], [4, 2]], external=[1, 4])     # sigma*G
        @test g7.factor == 1
        @test external(g7) == [𝑓⁻(7), 𝑓⁻(2)]
        @test external_with_ghost(g7) == [𝑓⁻(7), 𝑓⁺ₑ(8)]

        V8 = [𝑓⁻ₑ(1), 𝑓⁺(2), 𝑓⁻(12), 𝑓⁺ₑ(16), 𝑓⁺(3)𝑓⁻(4)𝜙(5), 𝑓⁺(6)𝑓⁻(7)𝜙(8), 𝑓⁺(9)𝑓⁻(10)𝜙(11), 𝑓⁺(13)𝑓⁻(14)𝜙(15)]
        g8 = feynman_diagram(V8, [[2, 6], [5, 9], [7, 16], [8, 15], [10, 13], [11, 3], [12, 4], [14, 1]], external=[2, 12, 3, 14])
        @test g8.factor == -1
        @test external(g8) == [𝑓⁺(2), 𝑓⁻(10), 𝑓⁻(12), 𝑓⁺(13)]
        @test external_with_ghost(g8) == [𝑓⁻ₑ(1), 𝑓⁺(2), 𝑓⁻(12), 𝑓⁺ₑ(16)]
        standardize_order!(g8)
        @test external_with_ghost(g8) == [𝑓⁺(2), 𝑓⁺ₑ(16), 𝑓⁻(12), 𝑓⁻ₑ(1)]
        @test Set(external(g8)) == Set([𝑓⁺(2), 𝑓⁻(10), 𝑓⁻(12), 𝑓⁺(13)])

        V8p = [𝑓⁺(2), 𝑓⁻ₑ(1), 𝑓⁻(12), 𝑓⁺ₑ(16), 𝑓⁺(3)𝑓⁻(4)𝜙(5), 𝑓⁺(6)𝑓⁻(7)𝜙(8), 𝑓⁺(9)𝑓⁻(10)𝜙(11), 𝑓⁺(13)𝑓⁻(14)𝜙(15)]
        g8p = feynman_diagram(V8p, [[1, 6], [5, 9], [7, 16], [8, 15], [10, 13], [11, 3], [12, 4], [14, 2]], external=[12, 1, 3, 14])
        @test g8p.factor == 1
    end
    @testset "f+f+f-f- interaction" begin
        V1 = [𝑓⁻ₑ(1), 𝑓⁻ₑ(2), 𝑓⁺(3), 𝑓⁺(4), 𝑓⁺(5)𝑓⁺(6)𝑓⁻(7)𝑓⁻(8), 𝑓⁺(9)𝑓⁺(10)𝑓⁻(11)𝑓⁻(12)]
        g1 = feynman_diagram(V1, [[1, 5], [2, 10], [3, 8], [4, 11], [6, 12], [7, 9]], external=[3, 4, 5, 10])
        g1p = feynman_diagram(V1, [[1, 10], [2, 5], [3, 8], [4, 11], [6, 12], [7, 9]], external=[3, 4, 5, 10])
        @test g1p.factor ≈ -g1.factor
        @test external(g1) == external(g1p)
        @test external_with_ghost(g1) == external_with_ghost(g1p)

        V2 = [𝑓⁻ₑ(1), 𝑓⁺ₑ(12), 𝑓⁺(2), 𝑓⁻(3), 𝑓⁺(4)𝑓⁺(5)𝑓⁻(6)𝑓⁻(7), 𝑓⁺(8)𝑓⁺(9)𝑓⁻(10)𝑓⁻(11)]
        g2 = feynman_diagram(V2, [[1, 9], [3, 8], [4, 5], [6, 12], [7, 10], [11, 2]], external=[3, 4, 9, 11])
        @test g2.factor == -1
        @test external(g2) == [𝑓⁺(2), 𝑓⁻(3), 𝑓⁺(8), 𝑓⁻(10)]
        @test external_with_ghost(g2) == [𝑓⁻ₑ(1), 𝑓⁺ₑ(12), 𝑓⁺(2), 𝑓⁻(3)]
        @test external_labels(g2) == [2, 3, 8, 10] # labels of external vertices
        @test external_with_ghost_labels(g2) == [1, 12, 2, 3] # labels of external vertices with ghost
        standardize_order!(g2)
        @test Set(external(g2)) == Set([𝑓⁺(2), 𝑓⁻(3), 𝑓⁺(8), 𝑓⁻(10)])
        @test external_with_ghost(g2) == [𝑓⁺ₑ(12), 𝑓⁺(2), 𝑓⁻(3), 𝑓⁻ₑ(1)]
        @test external_with_ghost_labels(g2) == [12, 2, 3, 1]
    end
    @testset "Multi-operator contractions" begin
        # multi-operator (>2) contractions
        Vm = [𝑓ₑ(1), 𝑓⁺(2)𝑓⁻(3)𝑏⁺(4), 𝜙(5)𝑓⁺(6)𝑓⁻(7), 𝑓(8)𝑏⁻(9)𝜙(10)]
        gm = feynman_diagram(Vm, [[2, 3, 4, 9], [5, 6, 7, 10], [8, 1]], external=[8])
        @test vertices(gm) == Vm
        @test gm.subgraph_factors == [1, 1]
        @test gm.subgraphs[1].vertices == external(gm.subgraphs[1]) == [𝑓⁺(2), 𝑓⁻(3), 𝑏⁺(4), 𝑏⁻(9)]
        @test gm.subgraphs[2].vertices == external(gm.subgraphs[2]) == [𝜙(5), 𝑓⁺(6), 𝑓⁻(7), 𝜙(10)]
        @test external_with_ghost(gm) == [𝑓ₑ(1)]
        @test external(gm) == [𝑓(8)]
        standardize_order!(gm)
        @test gm.subgraphs[1].factor == -1
        @test external(gm.subgraphs[1]) == [𝑓⁻(3), 𝑏⁻(9), 𝑏⁺(4), 𝑓⁺(2)]
        @test gm.subgraphs[2].factor == -1
        @test external(gm.subgraphs[2]) == [𝜙(5), 𝑓⁻(7), 𝜙(10), 𝑓⁺(6)]

        ggm = deepcopy(gm)
        ggm.id = 1000
        @test isequiv(gm, ggm, :id)
    end
    @testset "Construct feynman diagram from sub-diagrams" begin
        V1 = [𝜙(5), 𝜙(6), 𝜙(7), 𝜙(8)]
        g1 = feynman_diagram(V1, [[1, 2, 3, 4],], external=[1, 2, 3, 4])    #vacuum diagram
        V2 = [𝜙(9), 𝜙(10), 𝜙(11), 𝜙(12)]
        g2 = feynman_diagram(V2, [[1, 2, 3, 4],], external=[1, 2, 3, 4])    #vacuum diagram

        g = feynman_diagram([𝜙(1), 𝜙(2), 𝜙(3), 𝜙(4), g1, g2], [[1, 5], [2, 6], [7, 9], [8, 10], [3, 11], [4, 12]]; external=[1, 2, 3, 4])

        @test g.vertices[1:4] == [𝜙(1), 𝜙(2), 𝜙(3), 𝜙(4)]
        @test external(g) == [𝜙(1), 𝜙(2), 𝜙(3), 𝜙(4)]
        @test g.vertices[5] == 𝜙(5)𝜙(6)𝜙(7)𝜙(8)
        @test g.vertices[6] == 𝜙(9)𝜙(10)𝜙(11)𝜙(12)
        @test g.subgraphs[end-1] == g1
        @test g.subgraphs[end] == g2

        V3 = [𝑓⁻ₑ(1), 𝑓⁺ₑ(12), 𝑓⁺(2), 𝑓⁻(3), 𝑓⁺(4)𝑓⁺(5)𝑓⁻(6)𝑓⁻(7), 𝑓⁺(8)𝑓⁺(9)𝑓⁻(10)𝑓⁻(11)]
        g3 = feynman_diagram(V3, [[1, 9], [3, 8], [4, 5], [6, 12], [7, 10], [11, 2]], external=[3, 4, 9, 11])

        g4 = feynman_diagram([𝑓⁻ₑ(13), 𝑓⁺ₑ(14), 𝑓⁻(15), 𝑓⁺(16), g3],
            [[1, 5], [2, 6], [3, 7], [4, 8]],
            external=[3, 4, 5, 6]
        )
        @test g4.vertices[5] == 𝑓⁺(2)𝑓⁻(3)𝑓⁺(8)𝑓⁻(10)
        @test external(g4) == [𝑓⁻(15), 𝑓⁺(16), 𝑓⁺(2), 𝑓⁻(3)]
        @test external_with_ghost(g4) == [𝑓⁻ₑ(13), 𝑓⁺ₑ(14), 𝑓⁻(15), 𝑓⁺(16)]
    end

end

@testset "relabel and standardize_labels" begin
    using FeynmanDiagram.ComputationalGraphs

    @testset "relabel" begin
        # construct a graph
        V = [𝑓ₑ(1), 𝑓⁺(2)𝑓⁻(3)𝑏⁺(4), 𝜙(5)𝑓⁺(6)𝑓⁻(7), 𝑓(8)𝑏⁻(9)𝜙(10)]
        g1 = feynman_diagram(V, [[2, 3, 4, 9], [5, 6, 7, 10], [8, 1]], external=[8])

        map = Dict(4 => 1, 6 => 1, 8 => 1, 9 => 1, 10 => 1)
        g2 = relabel(g1, map)
        uniqlabels = ComputationalGraphs.collect_labels(g2)
        @test uniqlabels == [1, 2, 3, 5, 7]

        map = Dict([i => 1 for i in 2:10])
        g3 = relabel(g1, map)
        uniqlabels = ComputationalGraphs.collect_labels(g3)
        @test uniqlabels == [1,]
    end

    @testset "standardize_labels" begin
        V = [𝑓ₑ(1), 𝑓⁺(2)𝑓⁻(3)𝑏⁺(4), 𝜙(5)𝑓⁺(6)𝑓⁻(7), 𝑓(8)𝑏⁻(9)𝜙(10)]
        g1 = feynman_diagram(V, [[2, 3, 4, 9], [5, 6, 7, 10], [8, 1]], external=[8])

        map = Dict([i => (11 - i) for i in 1:5])
        g2 = relabel(g1, map)

        g3 = standardize_labels(g2)
        uniqlabels = ComputationalGraphs.collect_labels(g3)
        @test uniqlabels == [1, 2, 3, 4, 5]
    end
end

@testset "graph vector" begin
    import FeynmanDiagram.ComputationalGraphs as Graphs

    p1 = Graphs.propagator([𝑓⁺(1), 𝑓⁻(2)])
    p2 = Graphs.propagator([𝑓⁺(1), 𝑓⁻(3)])
    p3 = Graphs.propagator([𝑓⁺(2), 𝑓⁻(3)])

    gv = [p1, p2, p3]

    g1 = Graphs.group(gv, [1,])
    @test Set(g1[[𝑓⁺(1),]]) == Set([p1, p2])
    @test Set(g1[[𝑓⁺(2),]]) == Set([p3,])

    g2 = Graphs.group(gv, [2,])
    @test Set(g2[[𝑓⁻(2),]]) == Set([p1,])
    @test Set(g2[[𝑓⁻(3),]]) == Set([p2, p3])

    g3 = Graphs.group(gv, [1, 2])
    @test Set(g3[[𝑓⁺(1), 𝑓⁻(2)]]) == Set([p1,])
    @test Set(g3[[𝑓⁺(1), 𝑓⁻(3)]]) == Set([p2,])
    @test Set(g3[[𝑓⁺(2), 𝑓⁻(3)]]) == Set([p3,])
end
