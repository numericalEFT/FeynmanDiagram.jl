using FeynmanDiagram: ComputationalGraphs as Graphs

# 𝓞 represents a non-trivial unary operation
struct O <: Graphs.AbstractOperator end

# 𝓞1, 𝓞2, and 𝓞3 represent trivial unary operations
struct O1 <: Graphs.AbstractOperator end
struct O2 <: Graphs.AbstractOperator end
struct O3 <: Graphs.AbstractOperator end
Graphs.unary_istrivial(::Type{O}) where {O<:Union{O1,O2,O3}} = true

@testset verbose = true "Graph" begin
    @testset verbose = true "Operations" begin
        g1 = Graph([])
        g2 = 2 * g1
        g2p = Graph([]; factor=2)
        @testset "Equivalence" begin
            g1_new_instance = Graph([])
            # Test equivalence modulo fields id/factor
            @test isequiv(g1, g1_new_instance) == false
            @test isequiv(g1, g1_new_instance, :id)
            @test isequiv(g1, g2p, :id) == false
            @test isequiv(g1, g2p, :factor) == false
            @test isequiv(g1, g2p, :id, :factor)
            # Test inequivalence when subgraph lengths are different
            t = g1 + g1
            @test isequiv(t, g1, :id) == false
        end
        @testset "Scalar multiplication" begin
            @test g2.subgraph_factors == [2]
            @test g2.operator == Graphs.Prod
            g3 = g1 * 2
            @test g3.subgraph_factors == [2]
            @test g3.operator == Graphs.Prod
        end
        @testset "Addition" begin
            g3 = g1 + g2
            @test g3.factor == 1
            @test g3.subgraphs == [g1, g1]
            @test g3.subgraph_factors == [1, 2]
            @test g3.subgraphs[1].subgraph_factors == g1.subgraph_factors
            @test g3.operator == Graphs.Sum
        end
        @testset "Subtraction" begin
            g4 = g1 - g2
            @test g4.factor == 1
            @test g4.subgraphs == [g1, g1]
            @test g4.subgraph_factors == [1, -2]
            @test g4.subgraphs[1].subgraph_factors == g1.subgraph_factors
            @test g4.subgraphs[2].subgraph_factors == g1.subgraph_factors
            @test g4.operator == Graphs.Sum
        end
        @testset "Linear combinations" begin
            # Binary form
            # NOTE: since g2 = 2 * g1, 5g2 ↦ 10g1 in final expressions
            g5 = 3g1 + 5g2
            g5lc = linear_combination(g1, g2, 3, 5)
            @test g5lc.subgraphs == [g1, g1]
            @test g5lc.subgraph_factors == [3, 10]
            @test isequiv(g5, g5lc, :id)
            # Vector form
            g6lc = linear_combination([g1, g2, g5, g2, g1], [3, 5, 7, 9, 11])
            @test g6lc.subgraphs == [g1, g1, g5, g1, g1]
            @test g6lc.subgraph_factors == [3, 10, 7, 18, 11]
            # Test one-level merging of multiplicative chains
            g7lc = g1 + 2 * (3 * g1 + 5 * g2p)
            g7lc_expect = g1 + 2 * linear_combination([g1, g2p], [3, 5])
            @test isequiv(g7lc, g7lc_expect, :id)
        end
        @testset "Multiplicative chains" begin
            g6 = 7 * (5 * (3 * (2 * g1)))
            @test g6.subgraph_factors == [210]
            @test g6.subgraphs[1].subgraphs == g1.subgraphs
            @test g6.subgraphs[1].subgraph_factors == g1.subgraph_factors
            g7 = (((g1 * 2) * 3) * 5) * 7
            @test g7.subgraph_factors == [210]
            @test g7.subgraphs[1].subgraphs == g1.subgraphs
            @test g7.subgraphs[1].subgraph_factors == g1.subgraph_factors
        end
    end
    @testset verbose = true "Transformations" begin
        @testset "Replace subgraph" begin
            g1 = Graph([])
            g2 = Graph([]; factor=2)
            g3 = Graph([]; factor=3)
            gsum = g2 + g3
            groot = g1 + gsum
            replace_subgraph!(groot, g2, g3)
            @test isequiv(gsum.subgraphs[1], gsum.subgraphs[2])
            gnew = replace_subgraph(groot, g2, g3)
            @test isequiv(gnew, g1 + (g3 + g3), :id)
        end
        @testset "Prune trivial unary operations" begin
            g1 = Graph([])
            # +g1
            g2 = Graph([g1,]; operator=Graphs.Sum())
            # +(+g1)
            g3 = Graph([g2,]; operator=Graphs.Sum())
            # +2(+g1)
            g3p = Graph([g2,]; subgraph_factors=[2,], operator=Graphs.Sum())
            # +(+(+g1))
            g4 = Graph([g3,]; operator=Graphs.Sum())
            # +(+2(+g1))
            g4p = Graph([g3p,]; operator=Graphs.Sum())
            @test Graphs.unary_istrivial(Graphs.Prod)
            @test Graphs.unary_istrivial(Graphs.Sum)
            @test Graphs.merge_factorless_chain(g2) == g1
            @test Graphs.merge_factorless_chain(g3) == g1
            @test Graphs.merge_factorless_chain(g4) == g1
            @test Graphs.merge_factorless_chain(g3p) == g3p
            @test Graphs.merge_factorless_chain(g4p) == g3p
            g5 = Graph([g1,]; operator=O())
            @test Graphs.unary_istrivial(O) == false
            @test Graphs.merge_factorless_chain(g5) == g5
        end
        g1 = Graph([])
        g2 = Graph([g1,]; subgraph_factors=[5,], operator=Graphs.Prod())
        g3 = Graph([g2,]; subgraph_factors=[3,], operator=Graphs.Prod())
        # g = 2*(3*(5*g1))
        g = Graph([g3,]; subgraph_factors=[2,], operator=Graphs.Prod())
        # gp = 2*(3*(g1 + 5*g1))
        g2p = g1 + g2
        g3p = Graph([g2p,]; subgraph_factors=[3,], operator=Graphs.Prod())
        gp = Graph([g3p,]; subgraph_factors=[2,], operator=Graphs.Prod())
        @testset "Merge chains" begin
            # g ↦ 30*(*(*g1))
            g_merged = Graphs.merge_chain_prefactors(g)
            @test g_merged.subgraph_factors == [30,]
            @test all(isfactorless(node) for node in PreOrderDFS(eldest(g_merged)))
            # in-place form
            gc = deepcopy(g)
            Graphs.merge_chain_prefactors!(gc)
            @test isequiv(gc, g_merged, :id)
            # gp ↦ 6*(*(g1 + 5*g1))
            gp_merged = Graphs.merge_chain_prefactors(gp)
            @test gp_merged.subgraph_factors == [6,]
            @test isfactorless(eldest(gp)) == false
            @test isfactorless(eldest(gp_merged))
            @test eldest(eldest(gp_merged)) == g2p
            # g ↦ 30*g1
            g_merged = merge_chains(g)
            @test isequiv(g_merged, 30 * g1, :id)
            # in-place form
            merge_chains!(g)
            @test isequiv(g, 30 * g1, :id)
            # gp ↦ 6*(g1 + 5*g1)
            gp_merged = merge_chains(gp)
            @test isequiv(gp_merged, 6 * g2p, :id)
            # Test a generic trivial unary chain
            # *(O3(5 * O2(3 * O1(2 * h)))) ↦ 30 * h
            h = Graph([])
            h1 = Graph([h,]; subgraph_factors=[2,], operator=O1())
            h2 = Graph([h1,]; subgraph_factors=[3,], operator=O2())
            h3 = Graph([h2,]; subgraph_factors=[5,], operator=O3())
            h4 = Graph([h3,]; operator=Graphs.Prod())
            h4_merged = merge_chains(h4)
            @test isequiv(h4_merged, 30 * h, :id)
            # in-place form
            merge_chains!(h4)
            @test isequiv(h4, 30 * h, :id)
        end
        @testset "Merge prefactors" begin
            g1 = propagator(𝑓⁺(1)𝑓⁻(2))
            h1 = linear_combination(g1, g1, 1, 2)
            @test h1.subgraph_factors == [1, 2]
            h2 = merge_linear_combination(h1)
            @test h2.subgraph_factors == [3]
            @test length(h2.subgraphs) == 1
            @test h2.subgraphs[1] == g1
            g2 = propagator(𝑓⁺(1)𝑓⁻(2), factor=2)
            h3 = linear_combination(g1, g2, 1, 2)
            h4 = merge_linear_combination(h3)
            @test isequiv(h3, h4, :id)
            h5 = linear_combination([g1, g2, g2, g1], [3, 5, 7, 9])
            h6 = merge_linear_combination(h5)
            @test length(h6.subgraphs) == 2
            @test h6.subgraphs == [g1, g2]
            @test h6.subgraph_factors == [12, 12]
            g3 = 2 * g1
            h7 = linear_combination([g1, g3, g3, g1], [3, 5, 7, 9])
            h8 = merge_linear_combination(h7)
            @test length(h8.subgraphs) == 1
            @test h8.subgraphs == [g1]
            @test h8.subgraph_factors == [36]
        end
    end
    @testset verbose = true "Optimizations" begin
        @testset "Remove one-child parents" begin
            # h = O(7 * (5 * (3 * (2 * g)))) ↦ O(210 * g)
            g1 = Graph([])
            g2 = 2 * g1
            g3 = Graph([g2,]; subgraph_factors=[3,], operator=Graphs.Prod())
            g4 = Graph([g3,]; subgraph_factors=[5,], operator=Graphs.Prod())
            h = Graph([g4,]; subgraph_factors=[7,], operator=O())
            hvec = repeat([deepcopy(h)], 3)
            # Test on a single graph
            Graphs.merge_all_chains!(h)
            @test h.operator == O
            @test h.subgraph_factors == [210,]
            @test eldest(h) == g1
            # Test on a vector of graphs
            Graphs.merge_all_chains!(hvec)
            @test all(h.operator == O for h in hvec)
            @test all(h.subgraph_factors == [210,] for h in hvec)
            @test all(eldest(h) == g1 for h in hvec)

            g2 = 2 * g1
            g3 = Graph([g2,]; subgraph_factors=[3,], operator=Graphs.Prod())
            g4 = Graph([g3,]; subgraph_factors=[5,], operator=Graphs.Prod())
            h0 = Graph([g1, g4]; subgraph_factors=[2, 7], operator=O())
            Graphs.merge_all_chains!(h0)
            @test h0.subgraph_factors == [2, 210]
            @test h0.subgraphs[2] == g1

            h1 = Graph([h0]; subgraph_factors=[3,], operator=Graphs.Prod())
            h2 = Graph([h1]; subgraph_factors=[5,], operator=Graphs.Prod())
            h = Graph([h2]; subgraph_factors=[7,], operator=O())
            Graphs.merge_all_chains!(h)
            @test h.subgraph_factors == [105]
            @test eldest(h) == h0
        end
        @testset "merge all linear combinations" begin
            g1 = Graph([])
            g2 = 2 * g1
            g3 = Graph([], factor=3.0)
            h = Graph([g1, g1, g3], subgraph_factors=[-1, 3, 1])
            h0 = Graph([deepcopy(h), g2])
            _h = Graph([g1, g3], subgraph_factors=[2, 1])
            hvec = repeat([deepcopy(h)], 3)
            # Test on a single graph
            Graphs.merge_all_linear_combinations!(h)
            @test isequiv(h, _h, :id)
            # Test on a vector of graphs
            Graphs.merge_all_linear_combinations!(hvec)
            @test all(isequiv(h, _h, :id) for h in hvec)

            Graphs.merge_all_linear_combinations!(h0)
            @test isequiv(h0.subgraphs[1], _h, :id)
        end
        @testset "optimize" begin
            g1 = Graph([])
            g2 = 2 * g1
            g3 = Graph([g2,]; subgraph_factors=[3,], operator=Graphs.Prod())
            g4 = Graph([g3,]; subgraph_factors=[5,], operator=Graphs.Prod())
            g5 = Graph([], factor=3.0)
            h0 = Graph([g1, g4, g5], subgraph_factors=[2, -1, 1])
            h1 = Graph([h0], operator=Graphs.Prod(), subgraph_factors=[2])
            h = Graph([h1, g5])
            _h = Graph([Graph([g1, g5], subgraph_factors=[-28, 1]), g5], subgraph_factors=[2, 1])

            hvec_op, leafMap = Graphs.optimize(repeat([deepcopy(h)], 3))
            leaf = rand(2)
            @test all(isequiv(h, _h, :id) for h in hvec_op)
            @test Graphs.eval!(hvec_op[1], leafMap, leaf) ≈ Graphs.eval!(h, leafMap, leaf)

            leafMap1 = Graphs.optimize!([h])
            @test isequiv(h, _h, :id, :weight)
        end
    end
end

@testset verbose = true "FeynmanGraph" begin
    @testset verbose = true "Operations" begin
        V = [interaction(𝑓⁺(1)𝑓⁻(2)𝑓⁺(3)𝑓⁻(4)), interaction(𝑓⁺(5)𝑓⁺(6)𝑓⁻(7)𝑓⁻(8)),
            external_vertex(𝑓⁺(9)), external_vertex(𝑓⁺(10))]
        g1 = FeynmanGraph(V; topology=[[2, 6], [3, 7], [4, 9], [8, 10]],
            external_indices=[1, 5, 9, 10], external_legs=[false, false, true, true])
        g2 = g1 * 2
        g2p = FeynmanGraph(V; topology=[[2, 6], [3, 7], [4, 9], [8, 10]],
            external_indices=[1, 5, 9, 10], external_legs=[false, false, true, true], factor=2)
        @testset "Properties" begin
            @test diagram_type(g1) == Graphs.GenericDiag
            @test orders(g1) == zeros(Int, 16)
            @test vertices(g1) == OperatorProduct[𝑓⁺(1)𝑓⁻(2)𝑓⁺(3)𝑓⁻(4), 𝑓⁺(5)𝑓⁺(6)𝑓⁻(7)𝑓⁻(8), 𝑓⁺(9), 𝑓⁺(10)]
            @test topology(g1) == [[2, 6], [3, 7], [4, 9], [8, 10]]
            @test external_indices(g1) == [1, 5, 9, 10]
            @test external_operators(g1) == 𝑓⁺(1)𝑓⁺(5)𝑓⁺(9)𝑓⁺(10)
            @test external_legs(g1) == [false, false, true, true]
            parameters = FeynmanProperties(
                diagram_type(g1),
                vertices(g1),
                topology(g1),
                external_indices(g1),
                external_legs(g1),
            )
            parameters_no_topology = FeynmanProperties(
                diagram_type(g1),
                vertices(g1),
                [],
                external_indices(g1),
                external_legs(g1),
            )
            @test parameters == g1.properties
            @test parameters != parameters_no_topology
            @test parameters_no_topology == drop_topology(g1.properties)
        end
        @testset "Equivalence" begin
            g1_new_instance = FeynmanGraph(V; topology=[[2, 6], [3, 7], [4, 9], [8, 10]],
                external_indices=[1, 5, 9, 10], external_legs=[false, false, true, true])
            g1_from_parameters = FeynmanGraph(V, g1.properties)
            # Test equivalence modulo fields id/factor
            @test isequiv(g1, g1_new_instance) == false
            @test isequiv(g1, g1_from_parameters) == false
            @test isequiv(g1, g2p, :id) == false
            @test isequiv(g1, g2p, :factor) == false
            @test isequiv(g1, g1_new_instance, :id)
            @test isequiv(g1, g1_from_parameters, :id)
            @test isequiv(g1, g2p, :id, :factor)
            # Test inequivalence when subgraph lengths are different
            t = g1 + g1
            @test isequiv(t, g1, :id) == false
        end
        @testset "Scalar multiplication" begin
            @test vertices(g2) == vertices(g1)
            println(external_operators(g2))
            println(external_operators(g1))
            @test external_operators(g2) == external_operators(g1)
            @test g2.subgraph_factors == [2]
            @test g2.operator == Graphs.Prod
            g2 = 2g1
            @test vertices(g2) == vertices(g1)
            @test external_operators(g2) == external_operators(g1)
            @test g2.subgraph_factors == [2]
            @test g2.operator == Graphs.Prod
        end
        @testset "Addition" begin
            g3 = g1 + g2
            @test vertices(g3) == vertices(g1)
            @test external_operators(g3) == external_operators(g1)
            @test g3.factor == 1
            @test g3.subgraphs == [g1, g1]
            @test g3.subgraph_factors == [1, 2]
            @test g3.subgraphs[1].subgraph_factors == g1.subgraph_factors
            @test g3.operator == Graphs.Sum
        end
        @testset "Subtraction" begin
            g4 = g1 - g2
            @test vertices(g4) == vertices(g1)
            @test external_operators(g4) == external_operators(g1)
            @test g4.factor == 1
            @test g4.subgraphs == [g1, g1]
            @test g4.subgraph_factors == [1, -2]
            @test g4.subgraphs[1].subgraph_factors == g1.subgraph_factors
            @test g4.subgraphs[2].subgraph_factors == g1.subgraph_factors
            @test g4.operator == Graphs.Sum
        end
        @testset "Linear combinations" begin
            # Binary form
            # NOTE: since g2 = 2 * g1, 5g2 ↦ 10g1 in final expressions
            g5 = 3g1 + 5g2
            g5lc = linear_combination(g1, g2, 3, 5)
            @test g5lc.subgraphs == [g1, g1]
            @test g5lc.subgraph_factors == [3, 10]
            @test isequiv(g5, g5lc, :id)
            # Vector form
            g6lc = linear_combination([g1, g2, g5, g2, g1], [3, 5, 7, 9, 11])
            @test g6lc.subgraphs == [g1, g1, g5, g1, g1]
            @test g6lc.subgraph_factors == [3, 10, 7, 18, 11]
            # Test one-level merging of multiplicative chains
            g7lc = g1 + 2 * (3 * g1 + 5 * g2p)
            g7lc_expect = g1 + 2 * linear_combination([g1, g2p], [3, 5])
            @test isequiv(g7lc, g7lc_expect, :id)
        end
        @testset "Multiplicative chains" begin
            g6 = 7 * (5 * (3 * (2 * g1)))
            @test g6.subgraph_factors == [210]
            @test g6.subgraphs[1].subgraphs == g1.subgraphs
            @test g6.subgraphs[1].subgraph_factors == g1.subgraph_factors
            g7 = (((g1 * 2) * 3) * 5) * 7
            @test g7.subgraph_factors == [210]
            @test g7.subgraphs[1].subgraphs == g1.subgraphs
            @test g7.subgraphs[1].subgraph_factors == g1.subgraph_factors
        end
    end

    @testset verbose = true "Transformations" begin
        @testset "Relabel" begin
            # construct a graph
            V = [𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6), 𝑓⁺(7)𝑓⁻(8)𝜙(9)]
            g1 = feynman_diagram(interaction.(V), [[1, 5], [3, 9], [4, 8]])

            map = Dict(3 => 1, 4 => 1, 5 => 1, 9 => 1, 8 => 1)
            g2 = relabel(g1, map)
            uniqlabels = Graphs.collect_labels(g2)
            @test uniqlabels == [1, 2, 6, 7]

            map = Dict([i => 1 for i in 2:9])
            g3 = relabel(g1, map)
            uniqlabels = Graphs.collect_labels(g3)
            @test uniqlabels == [1,]
        end
        @testset "Standardize labels" begin
            V = [𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6), 𝑓⁺(7)𝑓⁻(8)𝜙(9), 𝑓⁺(10)]
            g1 = feynman_diagram([interaction.(V[1:3]); external_vertex(V[end])], [[1, 5], [3, 9], [4, 8], [2, 10]])

            map = Dict([i => (11 - i) for i in 1:5])
            g2 = relabel(g1, map)

            g3 = standardize_labels(g2)
            uniqlabels = Graphs.collect_labels(g3)
            @test uniqlabels == [1, 2, 3, 4, 5]
        end
        @testset "Replace subgraph" begin
            V2 = [external_vertex(𝜙(1)), interaction(𝜙(2)𝜙(3)), external_vertex(𝜙(4))]
            g1 = feynman_diagram(V2, [[1, 2], [3, 4]])
            g2 = feynman_diagram(V2, [[1, 3], [2, 4]])
            g3 = feynman_diagram(V2, [[1, 4], [2, 3]])
            gsum = g2 + g3
            groot = g1 + gsum
            replace_subgraph!(groot, g2, g3)
            @test isequiv(gsum.subgraphs[1], gsum.subgraphs[2])
            gnew = replace_subgraph(groot, g2, g3)
            @test isequiv(gnew, g1 + (g3 + g3), :id)
        end
        @testset "Prune trivial unary operations" begin
            g1 = propagator(𝑓⁺(1)𝑓⁻(2))
            # +g1
            g2 = FeynmanGraph([g1,], drop_topology(g1.properties); operator=Graphs.Sum())
            # +(+g1)
            g3 = FeynmanGraph([g2,], drop_topology(g2.properties); operator=Graphs.Sum())
            # +2(+g1)
            g3p = FeynmanGraph([g2,], drop_topology(g2.properties); subgraph_factors=[2,], operator=Graphs.Sum())
            # +(+(+g1))
            g4 = FeynmanGraph([g3,], drop_topology(g3.properties); operator=Graphs.Sum())
            # +(+2(+g1))
            g4p = FeynmanGraph([g3p,], drop_topology(g3p.properties); operator=Graphs.Sum())
            @test Graphs.unary_istrivial(Graphs.Prod)
            @test Graphs.unary_istrivial(Graphs.Sum)
            @test Graphs.merge_factorless_chain(g2) == g1
            @test Graphs.merge_factorless_chain(g3) == g1
            @test Graphs.merge_factorless_chain(g4) == g1
            @test Graphs.merge_factorless_chain(g3p) == g3p
            @test Graphs.merge_factorless_chain(g4p) == g3p
            g5 = FeynmanGraph([g1,], drop_topology(g1.properties); operator=O())
            @test Graphs.unary_istrivial(O) == false
            @test Graphs.merge_factorless_chain(g5) == g5
        end
        g1 = propagator(𝑓⁻(1)𝑓⁺(2))
        g2 = FeynmanGraph([g1,], g1.properties; subgraph_factors=[5,], operator=Graphs.Prod())
        g3 = FeynmanGraph([g2,], g2.properties; subgraph_factors=[3,], operator=Graphs.Prod())
        # g = 2*(3*(5*g1))
        g = FeynmanGraph([g3,], g3.properties; subgraph_factors=[2,], operator=Graphs.Prod())
        # gp = 2*(3*(g1 + 5*g1))
        g2p = g1 + g2
        g3p = FeynmanGraph([g2p,], g2p.properties; subgraph_factors=[3,], operator=Graphs.Prod())
        gp = FeynmanGraph([g3p,], g3p.properties; subgraph_factors=[2,], operator=Graphs.Prod())
        @testset "Merge chains" begin
            # g ↦ 30*(*(*g1))
            g_merged = Graphs.merge_chain_prefactors(g)
            @test g_merged.subgraph_factors == [30,]
            @test all(isfactorless(node) for node in PreOrderDFS(eldest(g_merged)))
            # in-place form
            gc = deepcopy(g)
            Graphs.merge_chain_prefactors!(gc)
            @test isequiv(gc, g_merged, :id)
            # gp ↦ 6*(*(g1 + 5*g1))
            gp_merged = Graphs.merge_chain_prefactors(gp)
            @test gp_merged.subgraph_factors == [6,]
            @test isfactorless(eldest(gp)) == false
            @test isfactorless(eldest(gp_merged))
            @test isequiv(eldest(eldest(gp_merged)), g2p, :id)
            # g ↦ 30*g1
            g_merged = merge_chains(g)
            @test isequiv(g_merged, 30 * g1, :id)
            # in-place form
            merge_chains!(g)
            @test isequiv(g, 30 * g1, :id)
            # gp ↦ 6*(g1 + 5*g1)
            gp_merged = merge_chains(gp)
            @test isequiv(gp_merged, 6 * g2p, :id)
            # Test a generic trivial unary chain
            # *(O3(5 * O2(3 * O1(2 * h)))) ↦ 30 * h
            h = propagator(𝑓⁻(1)𝑓⁺(2))
            h1 = FeynmanGraph([h,], h.properties; subgraph_factors=[2,], operator=O1())
            h2 = FeynmanGraph([h1,], h1.properties; subgraph_factors=[3,], operator=O2())
            h3 = FeynmanGraph([h2,], h2.properties; subgraph_factors=[5,], operator=O3())
            h4 = FeynmanGraph([h3,], h3.properties; operator=Graphs.Prod())
            h4_merged = merge_chains(h4)
            @test isequiv(h4_merged, 30 * h, :id)
            # in-place form
            merge_chains!(h4)
            @test isequiv(h4, 30 * h, :id)
        end
        @testset "Merge prefactors" begin
            g1 = propagator(𝑓⁺(1)𝑓⁻(2))
            h1 = linear_combination(g1, g1, 1, 2)
            @test h1.subgraph_factors == [1, 2]
            h2 = merge_linear_combination(h1)
            @test h2.subgraph_factors == [3]
            @test length(h2.subgraphs) == 1
            @test isequiv(h2.subgraphs[1], g1, :id)
            g2 = propagator(𝑓⁺(1)𝑓⁻(2), factor=2)
            h3 = linear_combination(g1, g2, 1, 2)
            h4 = merge_linear_combination(h3)
            @test isequiv(h3, h4, :id)
            h5 = linear_combination([g1, g2, g2, g1], [3, 5, 7, 9])
            h6 = merge_linear_combination(h5)
            @test length(h6.subgraphs) == 2
            @test h6.subgraphs == [g1, g2]
            @test h6.subgraph_factors == [12, 12]
            g3 = 2 * g1
            h7 = linear_combination([g1, g3, g3, g1], [3, 5, 7, 9])
            h8 = merge_linear_combination(h7)
            @test length(h8.subgraphs) == 1
            @test h8.subgraphs == [g1]
            @test h8.subgraph_factors == [36]
        end
    end

    @testset verbose = true "Optimizations" begin
        @testset "Remove one-child parents" begin
            g1 = propagator(𝑓⁻(1)𝑓⁺(2))
            g2 = 2 * g1
            # h = O(7 * (5 * (3 * (2 * g)))) ↦ O(210 * g)
            g3 = FeynmanGraph([g2,], g2.properties; subgraph_factors=[3,], operator=Graphs.Prod())
            g4 = FeynmanGraph([g3,], g3.properties; subgraph_factors=[5,], operator=Graphs.Prod())
            h = FeynmanGraph([g4,], drop_topology(g4.properties); subgraph_factors=[7,], operator=O())
            hvec = repeat([h], 3)
            # Test on a single graph
            Graphs.merge_all_chains!(h)
            @test h.operator == O
            @test h.subgraph_factors == [210,]
            @test isequiv(eldest(h), g1, :id)
            # Test on a vector of graphs
            Graphs.merge_all_chains!(hvec)
            @test all(h.operator == O for h in hvec)
            @test all(h.subgraph_factors == [210,] for h in hvec)
            @test all(isequiv(eldest(h), g1, :id) for h in hvec)

            g2 = 2 * g1
            g3 = FeynmanGraph([g2,], g2.properties; subgraph_factors=[3,], operator=Graphs.Prod())
            g4 = FeynmanGraph([g3,], g3.properties; subgraph_factors=[5,], operator=Graphs.Prod())
            h = FeynmanGraph([g1, g4], drop_topology(g4.properties); subgraph_factors=[2, 7], operator=O())
            Graphs.merge_all_chains!(h)
            @test h.subgraph_factors == [2, 210]
        end
        @testset "optimize" begin
            g1 = propagator(𝑓⁻(1)𝑓⁺(2))
            g2 = 2 * g1
            g3 = FeynmanGraph([g2,], g2.properties; subgraph_factors=[3,], operator=Graphs.Prod())
            g4 = FeynmanGraph([g3,], g3.properties; subgraph_factors=[5,], operator=Graphs.Prod())
            g5 = propagator(𝑓⁻(1)𝑓⁺(2), factor=3.0)
            h0 = FeynmanGraph([g1, g4, g5], subgraph_factors=[2, -1, 1])
            h1 = FeynmanGraph([h0], operator=Graphs.Prod(), subgraph_factors=[2])
            h = FeynmanGraph([h1, g5])
            _h = FeynmanGraph([FeynmanGraph([g1, g5], subgraph_factors=[-28, 1]), g5], subgraph_factors=[2, 1])

            hvec_op, leafMap = Graphs.optimize(repeat([deepcopy(h)], 3))
            leaf = rand(2)
            @test all(isequiv(h, _h, :id) for h in hvec_op)
            @test Graphs.eval!(hvec_op[1], leafMap, leaf) ≈ Graphs.eval!(h, leafMap, leaf)

            leafMap1 = Graphs.optimize!([h])
            @test isequiv(h, _h, :id, :weight)
        end
    end

    @testset "FeynmanGraphVector" begin
        p1 = propagator(𝑓⁺(1)𝑓⁻(2))
        p2 = propagator(𝑓⁺(1)𝑓⁻(3))
        p3 = propagator(𝑓⁺(2)𝑓⁻(3))

        gv = [p1, p2, p3]

        g1 = Graphs.group(gv, [2,])
        @test Set(g1[[𝑓⁺(1),]]) == Set([p1, p2])
        @test Set(g1[[𝑓⁺(2),]]) == Set([p3,])

        g2 = Graphs.group(gv, [1,])
        @test Set(g2[[𝑓⁻(2),]]) == Set([p1,])
        @test Set(g2[[𝑓⁻(3),]]) == Set([p2, p3])

        g3 = Graphs.group(gv, [2, 1])
        @test Set(g3[[𝑓⁺(1), 𝑓⁻(2)]]) == Set([p1,])
        @test Set(g3[[𝑓⁺(1), 𝑓⁻(3)]]) == Set([p2,])
        @test Set(g3[[𝑓⁺(2), 𝑓⁻(3)]]) == Set([p3,])
    end

    @testset "Propagator" begin
        g1 = propagator(𝑓⁺(1)𝑓⁻(2))
        @test g1.factor == -1
        @test external_indices(g1) == [2, 1]
        @test vertices(g1) == [𝑓⁺(1), 𝑓⁻(2)]
        @test external_operators(g1) == 𝑓⁻(2)𝑓⁺(1)
        @test external_labels(g1) == [2, 1]
    end

    @testset "Interaction" begin
        ops = 𝑓⁺(1)𝑓⁻(2)𝑓⁻(3)𝑓⁺(4)𝜙(5)
        g1 = interaction(ops)
        @test g1.factor == 1
        @test external_indices(g1) == [1, 2, 3, 4, 5]
        @test vertices(g1) == [ops]
        @test external_operators(g1) == ops
        @test external_labels(g1) == [1, 2, 3, 4, 5]

        g2 = interaction(ops, reorder=normal_order)
        @test g2.factor == -1
        @test vertices(g2) == [ops]
        @test external_operators(g2) == 𝑓⁺(1)𝑓⁺(4)𝜙(5)𝑓⁻(3)𝑓⁻(2)
        @test external_labels(g2) == [1, 4, 5, 3, 2]
    end

    @testset verbose = true "Feynman diagram" begin
        @testset "Phi4" begin
            # phi theory 
            V1 = [interaction(𝜙(1)𝜙(2)𝜙(3)𝜙(4))]
            g1 = feynman_diagram(V1, [[1, 2], [3, 4]])    #vacuum diagram
            @test vertices(g1) == [𝜙(1)𝜙(2)𝜙(3)𝜙(4)]
            @test isempty(external_operators(g1))
            @test g1.subgraph_factors == [1, 1, 1]
        end
        @testset "Complex scalar field" begin
            #complex scalar field
            V2 = [𝑏⁺(1), 𝑏⁺(2)𝑏⁺(3)𝑏⁻(4)𝑏⁻(5), 𝑏⁺(6)𝑏⁺(7)𝑏⁻(8)𝑏⁻(9), 𝑏⁻(10)]
            g2V = [external_vertex(V2[1]), interaction(V2[2]), interaction(V2[3]), external_vertex(V2[4])]
            g2 = feynman_diagram(g2V, [[1, 5], [2, 8], [3, 9], [4, 6], [7, 10]])    # Green2
            @test vertices(g2) == V2
            @test external_operators(g2) == 𝑏⁺(1)𝑏⁻(10)
            @test g2.subgraph_factors == ones(Int, 9)
        end
        @testset "Yukawa interaction" begin
            # Yukawa 
            V3 = [𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6)]
            g3 = feynman_diagram(interaction.(V3), [[1, 5], [2, 4], [3, 6]])  #vacuum diagram
            @test vertices(g3) == V3
            @test isempty(external_operators(g3))
            @test g3.factor == 1
            @test g3.subgraph_factors == ones(Int, 5)
            @test g3.subgraphs[3].factor == -1
            @test vertices(g3.subgraphs[3]) == [𝑓⁺(1), 𝑓⁻(5)]
            @test external_operators(g3.subgraphs[3]) == 𝑓⁻(5)𝑓⁺(1)

            V4 = [𝑓⁺(1)𝑓⁻(2), 𝑓⁺(3)𝑓⁻(4)𝜙(5), 𝑓⁺(6)𝑓⁻(7)𝜙(8), 𝑓⁺(9)𝑓⁻(10)]
            g4 = feynman_diagram([external_vertex(V4[1]), interaction.(V4[2:3])..., external_vertex(V4[4])],
                [[1, 4], [2, 6], [3, 10], [5, 8], [7, 9]]) # polarization diagram
            @test g4.factor == -1
            @test g4.subgraph_factors == ones(Int, 9)
            @test vertices(g4) == V4
            @test external_operators(g4) == 𝑓⁺(1)𝑓⁻(2)𝑓⁺(9)𝑓⁻(10)

            V5 = [𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6), 𝑓⁺(7)𝑓⁻(8)𝜙(9)]
            g5 = feynman_diagram(interaction.(V5), [[1, 5], [3, 9], [4, 8]])  # vertex function
            @test g5.factor == -1
            @test g5.subgraph_factors == ones(Int, 6)
            @test vertices(g5) == V5
            @test external_operators(g5) == 𝑓⁻(2)𝜙(6)𝑓⁺(7)
            g5p = feynman_diagram(interaction.(V5), [[1, 5], [3, 9], [4, 8]], [3, 1, 2])
            @test g5.factor ≈ -g5p.factor    # reorder of external fake legs will not change the sign.
            @test g5p.subgraph_factors == ones(Int, 6)
            @test external_operators(g5p) == 𝑓⁺(7)𝑓⁻(2)𝜙(6)

            V6 = [𝑓⁻(8), 𝑓⁺(1), 𝑓⁺(2)𝑓⁻(3)𝜙(4), 𝑓⁺(5)𝑓⁻(6)𝜙(7)]
            g6 = feynman_diagram([external_vertex.(V6[1:2]); interaction.(V6[3:4])], [[2, 4], [3, 7], [5, 8], [6, 1]])    # fermionic Green2
            @test g6.factor == -1
            @test g6.subgraph_factors == ones(Int, 8)
            @test external_operators(g6) == 𝑓⁻(8)𝑓⁺(1)

            V7 = [𝑓⁻(7), 𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6)]
            g7 = feynman_diagram([external_vertex(V7[1]), interaction.(V7[2:3])...], [[2, 6], [4, 7], [5, 1]])     # sigma*G
            @test g7.factor == 1
            @test external_operators(g7) == 𝑓⁻(7)𝑓⁻(2)

            V8 = [𝑓⁺(2), 𝑓⁻(12), 𝑓⁺(3)𝑓⁻(4)𝜙(5), 𝑓⁺(6)𝑓⁻(7)𝜙(8), 𝑓⁺(9)𝑓⁻(10)𝜙(11), 𝑓⁺(13)𝑓⁻(14)𝜙(15)]
            g8 = feynman_diagram([external_vertex.(V8[1:2]); interaction.(V8[3:end])], [[1, 4], [3, 7], [5, 14], [6, 13], [8, 11], [9, 2]])
            @test g8.factor == -1
            @test vertices(g8) == V8
            @test external_operators(g8) == 𝑓⁺(2)𝑓⁻(12)𝑓⁻(10)𝑓⁺(13)

            g8p = feynman_diagram([external_vertex.(V8[1:2]); interaction.(V8[3:end])],
                [[1, 4], [3, 7], [5, 14], [6, 13], [8, 11], [9, 2]], [2, 1])
            @test g8p.factor == 1
            @test external_operators(g8p) == 𝑓⁺(2)𝑓⁻(12)𝑓⁺(13)𝑓⁻(10)
        end
        @testset "f+f+f-f- interaction" begin
            V1 = [𝑓⁺(3), 𝑓⁺(4), 𝑓⁺(5)𝑓⁺(6)𝑓⁻(7)𝑓⁻(8), 𝑓⁺(9)𝑓⁺(10)𝑓⁻(11)𝑓⁻(12)]
            g1 = feynman_diagram([external_vertex.(V1[1:2]); interaction.(V1[3:4])], [[1, 6], [2, 9], [4, 10], [5, 7]])
            g1p = feynman_diagram([external_vertex.(V1[2:-1:1]); interaction.(V1[3:4])],
                [[2, 6], [1, 9], [4, 10], [5, 7]], [2, 1])
            @test g1p.factor ≈ g1.factor
            @test external_operators(g1) == 𝑓⁺(3)𝑓⁺(4)𝑓⁺(5)𝑓⁺(10)
            @test vertices(g1p) == [𝑓⁺(4), 𝑓⁺(3), 𝑓⁺(5)𝑓⁺(6)𝑓⁻(7)𝑓⁻(8), 𝑓⁺(9)𝑓⁺(10)𝑓⁻(11)𝑓⁻(12)]
            @test external_operators(g1p) == 𝑓⁺(4)𝑓⁺(3)𝑓⁺(10)𝑓⁺(5)

            V2 = [𝑓⁺(2), 𝑓⁻(3), 𝑓⁺(4)𝑓⁺(5)𝑓⁻(6)𝑓⁻(7), 𝑓⁺(8)𝑓⁺(9)𝑓⁻(10)𝑓⁻(11)]
            g2 = feynman_diagram([external_vertex.(V2[1:2]); interaction.(V2[3:4])], [[1, 6], [2, 3], [4, 10], [5, 8]])
            @test g2.factor == -1
            @test external_operators(g2) == 𝑓⁺(2)𝑓⁻(3)𝑓⁺(8)𝑓⁻(10)
            @test external_labels(g2) == [2, 3, 8, 10] # labels of external vertices    
        end
        @testset "Construct feynman diagram from sub-diagrams" begin
            V1 = [𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6)]
            g1 = feynman_diagram(interaction.(V1), [[3, 6]])
            V2 = [𝑓⁺(7)𝑓⁻(8)𝜙(9), 𝑓⁺(10)𝑓⁻(11)𝜙(12)]
            g2 = feynman_diagram(interaction.(V2), [[3, 6]])

            V3 = [𝑓⁻(13), 𝑓⁻(14), 𝑓⁺(15), 𝑓⁺(16)]
            g = feynman_diagram([g1, g2, external_vertex.(V3)...], [[1, 6], [2, 12], [3, 9], [4, 5], [7, 10], [8, 11]])

            @test vertices(g) == [𝑓⁺(1)𝑓⁻(2)𝑓⁺(4)𝑓⁻(5), 𝑓⁺(7)𝑓⁻(8)𝑓⁺(10)𝑓⁻(11), V3...]
            @test external_operators(g) == reduce(*, V3)
        end
    end
end

@testset verbose = true "Tree properties" begin
    using FeynmanDiagram.ComputationalGraphs:
        haschildren, onechild, isleaf, isbranch, ischain, isfactorless, eldest
    # Leaves: gᵢ
    g1 = propagator(𝑓⁻(1)𝑓⁺(2))
    g2 = propagator(𝑓⁻(1)𝑓⁺(2), factor=2)
    # Branches: Ⓧ --- gᵢ
    g3 = 1 * g1
    g4 = 1 * g2
    g5 = 2 * g1
    # Chains: Ⓧ --- Ⓧ --- gᵢ (simplified by default)
    g6 = FeynmanGraph([g5,], g5.properties; subgraph_factors=[1,], operator=Graphs.Prod())
    g7 = FeynmanGraph([g3,], g3.properties; subgraph_factors=[2,], operator=Graphs.Prod())
    # General trees
    g8 = 2 * (3 * g1 + 5 * g2)
    g9 = g1 + 2 * (3 * g1 + 5 * g2)  # ↦ g1 + 2 * linear_combination([g1, g2], [3, 5])
    @testset "Leaves" begin
        @test haschildren(g1) == false
        @test onechild(g1) == false
        @test isleaf(g1)
        @test isbranch(g1) == false
        @test ischain(g1)
        @test isfactorless(g1)
        @test isfactorless(g2) == false
        @test_throws AssertionError eldest(g1)
    end
    @testset "Branches" begin
        @test haschildren(g3)
        @test onechild(g3)
        @test isleaf(g3) == false
        @test isbranch(g3)
        @test ischain(g3)
        @test isfactorless(g3)
        @test isfactorless(g4)
        @test isfactorless(g5) == false
        @test isleaf(eldest(g3))
    end
    @testset "Chains" begin
        @test haschildren(g6)
        @test onechild(g6)
        @test isleaf(g6) == false
        @test isbranch(g6) == false
        @test ischain(g6)
        @test isfactorless(g6)
        @test isfactorless(g7) == false
        @test isbranch(eldest(g6))
    end
    @testset "General" begin
        @test haschildren(g8)
        @test onechild(g8)
        @test isleaf(g8) == false
        @test isbranch(g8) == false
        @test ischain(g8) == false
        @test isfactorless(g8) == false
        @test onechild(eldest(g8)) == false
    end
    @testset "Iteration" begin
        count_pre = sum(1 for node in PreOrderDFS(g9))
        count_post = sum(1 for node in PostOrderDFS(g9))
        @test count_pre == 5
        @test count_post == 5
    end
end
