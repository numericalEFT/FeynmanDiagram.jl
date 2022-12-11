@testset verbose = true "Graph" begin
    V = [𝑓⁺(1)𝑓⁻(2), 𝑓⁺(5)𝑓⁺(6)𝑓⁻(7)𝑓⁻(8), 𝑓⁺(3)𝑓⁻(4)]
    g1 = Graph(V, external=[1, 3])
    g2 = g1 * 2
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
    g1 = propagator(𝑓⁺(1)𝑓⁻(2))
    @test g1.factor == 1
    @test g1.external == [1, 2]
    @test vertices(g1) == [𝑓⁺(1)𝑓⁻(2)]
    standardize_order!(g1)
    @test g1.factor == -1
    @test g1.external == [2, 1]
    @test OperatorProduct(external(g1)) == 𝑓⁻(2)𝑓⁺(1)
    # @test vertices(g1) == [𝑓⁻(2)𝑓⁺(1)]

    g2 = propagator(𝑓⁺(1)𝑓⁻(2)𝑏⁺(1)𝜙(1)𝑓⁺(3)𝑓⁻(1)𝑓(1)𝑏⁻(1)𝜙(1))
    @test vertices(g2) == [𝑓⁺(1)𝑓⁻(2)𝑏⁺(1)𝜙(1)𝑓⁺(3)𝑓⁻(1)𝑓(1)𝑏⁻(1)𝜙(1)]
    standardize_order!(g2)
    @test g2.factor == -1
    @test OperatorProduct(external(g2)) == 𝑓⁻(1)𝑏⁻(1)𝜙(1)𝑓⁻(2)𝑓(1)𝑓⁺(3)𝜙(1)𝑏⁺(1)𝑓⁺(1)
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
        # @test internal_vertices(g1) == V1
        # @test isequiv(g1, gg1, :id)
    end
    @testset "Complex scalar field" begin
        #complex scalar field
        V2 = [𝑏⁺(1), 𝑏⁺(2)𝑏⁺(3)𝑏⁻(4)𝑏⁻(5), 𝑏⁺(6)𝑏⁺(7)𝑏⁻(8)𝑏⁻(9), 𝑏⁻(10)]
        # g2 = feynman_diagram(V2, [1, 2, 3, 4, 1, 4, 5, 2, 3, 5]; external=[1, 10])
        g2 = feynman_diagram(V2, [[1, 5], [2, 8], [3, 9], [4, 6], [7, 10]]; external=[1, 10])    # Green2
        @test vertices(g2) == V2
        @test external(g2) == OperatorProduct(V2)[[1, 10]]
        @test g2.subgraph_factors == [1, 1, 1, 1, 1]
        # @test internal_vertices(g2) == V2[2:3]
        # @test isequiv(g2, gg2, :id)
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
        @test g3.subgraphs[1].vertices == [𝑓⁺(1)𝑓⁻(5)]
        @test g3.subgraphs[1].factor == 1
        @test OperatorProduct(external(g3.subgraphs[1])) == 𝑓⁺(1)𝑓⁻(5)
        standardize_order!(g3)
        @test g3.subgraphs[1].factor == -1
        @test OperatorProduct(external(g3.subgraphs[1])) == 𝑓⁻(5)𝑓⁺(1)
        # @test gg3.subgraphs[1].factor == 1
        # @test !isequiv(g3, gg3, :id)

        V4 = [𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6), 𝑓⁺(7)𝑓⁻(8)𝜙(9), 𝑓⁺(10)𝑓⁻(11)𝜙(12)]
        gg4 = feynman_diagram(V4, [[1, 8], [2, 10], [4, 11], [5, 7], [9, 12]]) # polarization diagram
        @test gg4.factor == -1
        @test gg4.subgraph_factors == [1, 1, 1, 1, 1]
        @test vertices(gg4) == V4
        @test external(gg4) == OperatorProduct(V4)[[3, 6]]
        standardize_order!(gg4)
        @test external(gg4) == OperatorProduct(V4)[[3, 6]]
        # @test internal_vertices(g4) == V4[3:4]
        # g4 = feynman_diagram(V4, [1, 2, 0, 3, 4, 0, 4, 1, 5, 2, 3, 5], external=[3, 6])
        # @test isequiv(g4, gg4, :id)

        V5 = [𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6), 𝑓⁺(7)𝑓⁻(8)𝜙(9)]
        g5 = feynman_diagram(V5, [[2, 4], [3, 9], [5, 7]], external=[1, 6, 8])  # vertex function
        @test g5.factor == 1
        @test g5.subgraph_factors == [1, 1, 1]
        @test vertices(g5) == V5
        @test external(g5) == OperatorProduct(V5)[[1, 6, 8]]
        # @test isempty(internal_vertices(g5))
        g5s = deepcopy(g5)
        standardize_order!(g5)
        @test external(g5) == OperatorProduct(V5)[[1, 6, 8]]
        @test g5s == g5

        gg5 = feynman_diagram(V5, [[2, 4], [3, 9], [5, 7]], external=[8, 6, 1])
        @test g5.factor ≈ gg5.factor * (-1)
        @test gg5.subgraph_factors == [1, 1, 1]
        standardize_order!(gg5)
        @test external(gg5) == OperatorProduct(V5)[[6, 1, 8]]
        @test gg5.factor ≈ g5.factor
        # gg5 = feynman_diagram(V5, [1, 2, 1, 3, 3, 2], external=[1, 2, 3])
        ggg5 = feynman_diagram(V5, [[2, 4], [3, 9], [5, 7]])
        @test isequiv(g5, ggg5, :id, :type)

        V6 = [𝑓⁺(1), 𝑓⁺(2)𝑓⁻(3)𝜙(4), 𝑓⁺(5)𝑓⁻(6)𝜙(7), 𝑓⁻(8)]
        g6 = feynman_diagram(V6, [[1, 3], [2, 6], [4, 7], [5, 8]], external=[8, 1])    # fermionic Green2
        @test g6.factor == -1
        @test g6.subgraph_factors == [1, 1, 1, 1]
        @test external(g6) == OperatorProduct(V6)[[8, 1]]
        standardize_order!(g6)
        @test g6.factor == 1
        @test external(g6) == OperatorProduct(V6)[[1, 8]]

        V7 = [𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6), 𝑓⁻(7)]
        g7 = feynman_diagram(V7, [[1, 5], [3, 6], [4, 7]], external=[2, 7])     # sigma*G
        @test g7.factor == 1

        V8 = [𝑓⁺(1), 𝑓⁺(2)𝑓⁻(3)𝜙(4), 𝑓⁺(5)𝑓⁻(6)𝜙(7), 𝑓⁺(8)𝑓⁻(9)𝜙(10), 𝑓⁻(11), 𝑓⁺(12)𝑓⁻(13)𝜙(14)]
        g8 = feynman_diagram(V8, [[1, 3], [2, 6], [4, 14], [5, 13], [7, 10], [8, 11]], external=[1, 9, 11, 12])
        @test g8.factor == -1
        gg8 = feynman_diagram(V8, [[1, 3], [2, 6], [4, 14], [5, 13], [7, 10], [8, 11]], external=[9, 1, 11, 12])
        @test gg8.factor == -1
    end
    @testset "Multi-operator contractions" begin
        # multi-operator (>2) contractions
        Vm = [𝑓⁺(1)𝑓⁻(2)𝑏⁺(3), 𝜙(4)𝑓⁺(5)𝑓⁻(6), 𝑓(7)𝑏⁻(8)𝜙(9)]
        gm = feynman_diagram(Vm, [[1, 2, 3, 8], [4, 5, 6, 9]])
        @test vertices(gm) == Vm
        @test gm.subgraph_factors == [1, 1]
        @test gm.subgraphs[1].vertices == [𝑓⁺(1)𝑓⁻(2)𝑏⁺(3)𝑏⁻(8)]
        @test gm.subgraphs[2].vertices == [𝜙(4)𝑓⁺(5)𝑓⁻(6)𝜙(9)]
        @test OperatorProduct(external(gm.subgraphs[1])) == 𝑓⁺(1)𝑓⁻(2)𝑏⁺(3)𝑏⁻(8)
        @test OperatorProduct(external(gm.subgraphs[2])) == 𝜙(4)𝑓⁺(5)𝑓⁻(6)𝜙(9)
        standardize_order!(gm)
        @test gm.subgraphs[1].factor == -1
        @test OperatorProduct(external(gm.subgraphs[1])) == 𝑓⁻(2)𝑏⁻(8)𝑏⁺(3)𝑓⁺(1)
        @test gm.subgraphs[2].factor == -1
        @test OperatorProduct(external(gm.subgraphs[2])) == 𝜙(4)𝑓⁻(6)𝜙(9)𝑓⁺(5)

        ggm = deepcopy(gm)
        ggm.id = 1000
        @test isequiv(gm, ggm, :id)
    end
end
