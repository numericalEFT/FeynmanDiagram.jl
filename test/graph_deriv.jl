using FeynmanDiagram: ComputationalGraphs as Graphs

# 𝓞 represents a non-trivial unary operation
struct O <: Graphs.AbstractOperator end

# 𝓞1, 𝓞2, and 𝓞3 represent trivial unary operations
struct O1 <: Graphs.AbstractOperator end
struct O2 <: Graphs.AbstractOperator end
struct O3 <: Graphs.AbstractOperator end
Graphs.unary_istrivial(::Type{O}) where {O<:Union{O1,O2,O3}} = true

@testset verbose = true "Auto Differentiation" begin
    using FeynmanDiagram.ComputationalGraphs:
        eval!, forwardAD, node_derivative, backAD, forwardAD_root!, build_all_leaf_derivative, build_derivative_graph
    g1 = Graph([])
    g2 = Graph([])
    g3 = Graph([], factor=2.0)
    G3 = g1
    G4 = 4 * g1 * g1
    G5 = 4 * (2 * G3 + 3 * G4)
    G6 = (2 * g1 + 3 * g2) * (4 * g1 + g3) * g1
    #G6 = (g1 + g2) * (g1 + g2) * g1
    G7 = (3 * g1 + 4 * g2 + 5 * g3) * 3 * g1

    @testset "node_derivative" begin
        F1 = g1 * g1
        F2 = (3 * g1) * (4 * g1)
        F3 = (2 * g1 * g2) * (3 * g1)
        F4 = (2 * g1 + 3 * g2) + g1
        @test eval!(node_derivative(F1, g1)) == 2
        @test eval!(node_derivative(F2, g1)) == 24
        @test eval!(node_derivative(F1, g2)) == nothing
        @test eval!(node_derivative(F3, g1)) == 6 #The derivative is local, and only considers the children at root 
        print(node_derivative(F4, g1), "\n")
        @test eval!(node_derivative(F4, g1)) == 1
    end
    @testset "Eval" begin
        # Current test assign all green's function equal to 1 for simplicity.
        # print(eval!(forwardAD(G5, g1.id)),"\n")
        # print(eval!(forwardAD(G3, g1.id)),"\n")
        # print(eval!(forwardAD(G3, g2.id)),"\n")
        # print(eval!(forwardAD(G6, g1.id)),"\n")
        # print(eval!(forwardAD(forwardAD(G6, g1.id), g2.id)),"\n")
        # print(eval!(forwardAD(forwardAD(G6, g1.id), g3.id)),"\n")
        # gs = Compilers.to_julia_str([forwardAD(G5, g1.id),], name="eval_graph!")
        # println(gs,"\n")
        @test eval!(forwardAD(G3, g1.id)) == 1
        @test eval!(forwardAD(G4, g1.id)) == 8
        @test eval!(forwardAD(G5, g1.id)) == 104
        @test eval!(forwardAD(G6, g1.id)) == 62
        @test eval!(forwardAD(G6, g3.id)) == 5
        @test eval!(forwardAD(forwardAD(G6, g1.id), g2.id)) == 30
        #backAD(G5, true)
        for (i, G) in enumerate([G3, G4, G5, G6, G7])
            back_deriv = backAD(G)
            for (id_pair, value_back) in back_deriv
                # gs = Compilers.to_julia_str([value,], name="eval_graph!")
                # println("id:$(key)", gs, "\n")
                value_forward = forwardAD(G, id_pair[2])
                @test eval!(value_back) == eval!(value_forward)
                # print("value:$(i+2) $(eval!(value_forward))\n")
            end
        end
        # gs = Compilers.to_julia_str([G6,], name="eval_graph!")
        # println("G6  ", gs, "\n")
        # for (id, G) in backAD(G6)
        #     gs = Compilers.to_julia_str([G,], name="eval_graph!")
        #     println("first order derive id:$(id)", gs, "\n")
        #     back_deriv = backAD(G)
        #     for (id_pair, value_back) in back_deriv
        #         gs = Compilers.to_julia_str([value_back,], name="eval_graph!")
        #         println("second order derive id:$(id_pair)", gs, "\n")
        #         value_forward = forwardAD(G, id_pair[2])
        #         @test eval!(value_back) == eval!(value_forward)
        #         print("value:$(id_pair) $(eval!(value_forward))\n")
        #     end
        # end

        # for (order_vec, graph) in build_all_leaf_derivative(G6, 3)
        #     print("$(order_vec), $(eval!(graph)) \n")
        # end
    end
    @testset "forwardAD_root!" begin
        F3 = g1 + g2
        F2 = linear_combination([g1, g3, F3], [2, 1, 3])
        F1 = Graph([g1, F2, F3], operator=Graphs.Prod(), subgraph_factors=[3.0, 1.0, 1.0])

        kg1, kg2, kg3 = (g1.id, (1,)), (g2.id, (1,)), (g3.id, (1,))
        kF1, kF2, kF3 = (F1.id, (1,)), (F2.id, (1,)), (F3.id, (1,))

        dual = forwardAD_root!(F1)  # auto-differentation!
        @test dual[kF3].subgraphs == [dual[kg1], dual[kg2]]
        @test dual[kF2].subgraphs == [dual[kg1], dual[kg3], dual[kF3]]

        leafmap = Dict{Int,Int}()
        leafmap[g1.id], leafmap[g2.id], leafmap[g3.id] = 1, 2, 3
        leafmap[dual[kg1].id] = 4
        leafmap[dual[kg2].id] = 5
        leafmap[dual[kg3].id] = 6
        leaf = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0]   # d F1 / d g1
        @test eval!(dual[kF1], leafmap, leaf) == 120.0
        @test eval!(dual[kF2], leafmap, leaf) == 5.0
        @test eval!(dual[kF3], leafmap, leaf) == 1.0

        leaf = [5.0, -1.0, 2.0, 0.0, 1.0, 0.0]  # d F1 / d g2
        @test eval!(dual[kF1], leafmap, leaf) == 570.0
        @test eval!(dual[kF2], leafmap, leaf) == 3.0
        @test eval!(dual[kF3], leafmap, leaf) == 1.0

        leaf = [5.0, -1.0, 2.0, 0.0, 0.0, 1.0]  # d F1 / d g3
        @test eval!(dual[kF1], leafmap, leaf) == 60.0
        @test eval!(dual[kF2], leafmap, leaf) == 1.0
        @test eval!(dual[kF3], leafmap, leaf) == 0.0

        F0 = F1 * F3
        kF0 = (F0.id, (1,))
        dual1 = forwardAD_root!(F0)
        leafmap[dual1[kg1].id] = 4
        leafmap[dual1[kg2].id] = 5
        leafmap[dual1[kg3].id] = 6

        leaf = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0]
        @test eval!(dual1[kF0], leafmap, leaf) == 300.0
        leaf = [5.0, -1.0, 2.0, 0.0, 1.0, 0.0]
        @test eval!(dual1[kF0], leafmap, leaf) == 3840.0
        leaf = [5.0, -1.0, 2.0, 0.0, 0.0, 1.0]
        @test eval!(dual1[kF0], leafmap, leaf) == 240.0
        @test isequiv(dual[kF1], dual1[kF1], :id, :weight, :vertices)

        F0_r1 = F1 + F3
        kF0_r1 = (F0_r1.id, (1,))
        dual = forwardAD_root!([F0, F0_r1])
        leafmap[dual[kg1].id] = 4
        leafmap[dual[kg2].id] = 5
        leafmap[dual[kg3].id] = 6
        @test eval!(dual[kF0], leafmap, leaf) == 240.0
        @test eval!(dual[kF0_r1], leafmap, leaf) == 60.0
        @test isequiv(dual[kF0], dual1[kF0], :id, :weight)
        @test isequiv(dual[kF1], dual1[kF1], :id, :weight)
    end
    @testset "build_derivative_graph" begin
        F3 = g1 + g2
        F2 = linear_combination([g1, g3, F3], [2, 1, 3])
        F1 = Graph([g1, F2, F3], operator=Graphs.Prod(), subgraph_factors=[3.0, 1.0, 1.0])

        leafmap = Dict{Int,Int}()
        leafmap[g1.id], leafmap[g2.id], leafmap[g3.id] = 1, 2, 3
        orders = (3, 2, 2)
        dual = build_derivative_graph(F1, orders)

        leafmap[dual[(g1.id, (1, 0, 0))].id], leafmap[dual[(g2.id, (0, 1, 0))].id], leafmap[dual[(g3.id, (0, 0, 1))].id] = 4, 5, 6

        for order in Iterators.product((0:x for x in orders)...)
            order == (0, 0, 0) && continue
            for g in [g1, g2, g3]
                if !haskey(leafmap, dual[(g.id, order)].id)
                    leafmap[dual[(g.id, order)].id] = 7
                end
            end
        end
        leaf = [5.0, -1.0, 2.0, 1.0, 1.0, 1.0, 0.0]
        @test eval!(dual[(F1.id, (1, 0, 0))], leafmap, leaf) == 1002
        @test eval!(dual[(F1.id, (2, 0, 0))], leafmap, leaf) == 426
        @test eval!(dual[(F1.id, (3, 0, 0))], leafmap, leaf) == 90
        @test eval!(dual[(F1.id, (3, 1, 0))], leafmap, leaf) == 0
    end
end
