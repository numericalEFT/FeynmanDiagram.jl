@testset "Graph" begin
    V = [𝑓⁺(1)𝑓⁻(2), 𝑓⁺(5)𝑓⁺(6)𝑓⁻(7)𝑓⁻(8), 𝑓⁺(3)𝑓⁻(4)]
    g1 = Graph(V, external=[1, 3])
    g2 = g1 * 2
    @test vertices(g2) == vertices(g1)
    println(external(g2))
    println(external(g1))
    @test external(g2) == external(g1)
    @test g2.factor == 2
    @test g2.subgraph_factors == [1]
    @test g2.operator == FeynmanDiagram.ComputationalGraphs.Prod
    g2 = 2g1
    @test vertices(g2) == vertices(g1)
    @test external(g2) == external(g1)
    @test g2.factor == 2
    @test g2.subgraph_factors == [1]
    @test g2.operator == FeynmanDiagram.ComputationalGraphs.Prod
    g3 = g1 + g2
    @test vertices(g3) == vertices(g1)
    @test external(g3) == external(g1)
    @test g3.factor == 1
    @test g3.subgraph_factors == [1, 2]
    @test g3.subgraphs == [g1, g2]
    @test g3.operator == FeynmanDiagram.ComputationalGraphs.Sum
    g4 = g1 - g2
    @test vertices(g4) == vertices(g1)
    @test external(g4) == external(g1)
    @test g4.factor == 1
    @test g4.subgraph_factors == [1, -2]
    @test g4.subgraphs[2].factor == -2
    @test g4.subgraphs[2].subgraphs[1].factor == 2
    @test g4.operator == FeynmanDiagram.ComputationalGraphs.Sum
end

# @testset "Contractions" begin
#     # Test 1: Scalar fields with Wick crossings, parity = +1
#     vertices1 = [𝜙(1)𝜙(2), 𝜙(3)𝜙(4)𝜙(5)𝜙(6), 𝜙(7)𝜙(8)]
#     parity1, ind1, edges1 = contractions_to_edges(vertices1, [1, 2, 3, 4, 1, 3, 4, 2])
#     ops = reduce(*, vertices1)
#     @test ind1 == [[1, 5], [2, 8], [3, 6], [4, 7]]
#     @test Set(edges1) == Set([(ops[1], ops[5]), (ops[2], ops[8]), (ops[3], ops[6]), (ops[4], ops[7])])
#     @test parity1 == 1

#     # Test 2: Bosons with Wick crossings, parity = +1
#     vertices2 = [𝑏⁺(1)𝑏⁺(2)𝑏⁻(3), 𝑏⁻(4)𝑏⁺(5), 𝑏⁻(6)𝑏⁺(7)𝑏⁻(8)]
#     parity2, ind2, edges2 = contractions_to_edges(vertices2, [1, 2, 3, 4, 3, 1, 4, 2])
#     ops = reduce(*, vertices2)
#     @test ind2 == [[1, 6], [2, 8], [3, 5], [4, 7]]
#     @test Set(edges2) == Set([(ops[1], ops[6]), (ops[2], ops[8]), (ops[3], ops[5]), (ops[4], ops[7])])
#     @test parity2 == 1

#     # Test 3: Indistinguishable Majoranas with no Wick crossings, parity = +1
#     vertices3 = [𝑓(1)𝑓(1)𝑓(1)𝑓(1)𝑓(1)𝑓(1)𝑓(1)𝑓(1),]
#     parity3, ind3, edges3 = contractions_to_edges(vertices3, [1, 2, 3, 4, 4, 3, 2, 1])
#     ops = reduce(*, vertices3)
#     @test ind3 == [[1, 8], [2, 7], [3, 6], [4, 5]]
#     @test Set(edges3) == Set([(ops[1], ops[8]), (ops[2], ops[7]), (ops[3], ops[6]), (ops[4], ops[5])])
#     # P = (1 8 2 7 3 6 4 5) = (1)(5 3 2 8)(4 7)(6) => parity = +1
#     @test parity3 == 1

#     # Test 4: Fermions with Wick crossings, parity = -1
#     vertices4 = [𝑓⁺(1)𝑓⁻(2), 𝑓⁺(5)𝑓⁺(6)𝑓⁻(7)𝑓⁻(8), 𝑓⁺(3)𝑓⁻(4),]
#     parity4, ind4, edges4 = contractions_to_edges(vertices4, [1, 2, 2, 3, 1, 4, 4, 3])
#     ops = reduce(*, vertices4)
#     @test ind4 == [[1, 5], [2, 3], [4, 8], [6, 7]]
#     @test Set(edges4) == Set([(ops[1], ops[5]), (ops[2], ops[3]), (ops[4], ops[8]), (ops[6], ops[7])])
#     # P = (1 5 2 3 4 8 6 7) = (1)(2 5 4)(3)(6 8)(7) => parity = -1
#     @test parity4 == -1

#     # Test 5: Mixed bosonic/classical/fermionic operators, parity = -1
#     vertices5 = [𝑏⁺(1)𝑓⁺(2)𝜙(3), 𝑓⁻(4)𝑓⁻(5), 𝑏⁻(6)𝑓⁺(7)𝜙(8)]
#     ops = reduce(*, vertices5)
#     parity5, ind5, edges5 = contractions_to_edges(vertices5, [1, 2, 3, 2, 4, 1, 4, 3])
#     @test ind5 == [[1, 6], [2, 4], [3, 8], [5, 7]]
#     @test Set(edges5) == Set([(ops[1], ops[6]), (ops[2], ops[4]), (ops[3], ops[8]), (ops[5], ops[7])])
#     # Flattened fermionic edges: [2, 4, 5, 7]
#     # => P = (1 2 3 4) = (1)(2)(3 4) => parity = 1
#     @test parity5 == 1
# end

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

@testset "feynman_diagram" begin
    # phi theory 
    V1 = [𝜙(1)𝜙(1)𝜙(2)𝜙(2),]
    g1 = feynman_diagram(V1, [[1, 2], [3, 4]])    #vacuum diagram
    # g1 = feynman_diagram(V1, [1, 1, 2, 2])
    @test vertices(g1) == V1
    @test isempty(external(g1))
    # @test internal_vertices(g1) == V1
    # @test isequiv(g1, gg1, :id)


    #complex scalar field
    V2 = [𝑏⁺(1), 𝑏⁺(2)𝑏⁺(3)𝑏⁻(4)𝑏⁻(5), 𝑏⁺(6)𝑏⁺(7)𝑏⁻(8)𝑏⁻(9), 𝑏⁻(10)]
    # g2 = feynman_diagram(V2, [1, 2, 3, 4, 1, 4, 5, 2, 3, 5]; external=[1, 10])
    g2 = feynman_diagram(V2, [[1, 5], [2, 8], [3, 9], [4, 6], [7, 10]]; external=[1, 10], type=:green)    # Green2
    @test vertices(g2) == V2
    @test external(g2) == OperatorProduct(V2)[[1, 10]]
    standardize_order!(g2)
    @test external(g2) == OperatorProduct(V2)[[10, 1]]
    # @test internal_vertices(g2) == V2[2:3]
    # @test isequiv(g2, gg2, :id)

    # Yukawa 
    V3 = [𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6)]
    # g3 = feynman_diagram(V3, [1, 2, 3, 2, 1, 3])
    g3 = feynman_diagram(V3, [[1, 5], [2, 4], [3, 6]])  #vacuum diagram
    @test vertices(g3) == V3
    @test isempty(external(g3))
    @test g3.factor == 1
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
    gg4 = feynman_diagram(V4, [[1, 8], [2, 10], [4, 11], [5, 7], [9, 12]], type=:sigma) # polarization diagram
    @test gg4.factor == -1
    @test vertices(gg4) == V4
    @test external(gg4) == OperatorProduct(V4)[[3, 6]]
    standardize_order!(gg4)
    @test external(gg4) == OperatorProduct(V4)[[3, 6]]
    # @test internal_vertices(g4) == V4[3:4]
    # g4 = feynman_diagram(V4, [1, 2, 0, 3, 4, 0, 4, 1, 5, 2, 3, 5], external=[3, 6])
    # @test isequiv(g4, gg4, :id)

    V5 = [𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6), 𝑓⁺(7)𝑓⁻(8)𝜙(9)]
    g5 = feynman_diagram(V5, [[2, 4], [3, 9], [5, 7]], external=[1, 6, 8], type=:interaction)  # vertex function
    @test g5.factor == 1
    @test vertices(g5) == V5
    @test external(g5) == OperatorProduct(V5)[[1, 6, 8]]
    # @test isempty(internal_vertices(g5))
    g5s = deepcopy(g5)
    standardize_order!(g5)
    @test external(g5) == OperatorProduct(V5)[[1, 6, 8]]
    @test g5s == g5

    gg5 = feynman_diagram(V5, [[2, 4], [3, 9], [5, 7]], external=[8, 6, 1], type=:interaction)
    @test g5.factor ≈ gg5.factor * (-1)
    standardize_order!(gg5)
    @test external(gg5) == OperatorProduct(V5)[[6, 1, 8]]
    @test gg5.factor ≈ g5.factor
    # gg5 = feynman_diagram(V5, [1, 2, 1, 3, 3, 2], external=[1, 2, 3])
    ggg5 = feynman_diagram(V5, [[2, 4], [3, 9], [5, 7]])
    @test isequiv(g5, ggg5, :id, :type)


    V6 = [𝑓⁺(1), 𝑓⁺(2)𝑓⁻(3)𝜙(4), 𝑓⁺(5)𝑓⁻(6)𝜙(7), 𝑓⁻(8)]
    g6 = feynman_diagram(V6, [[1, 3], [2, 6], [4, 7], [5, 8]], external=[1, 8], type=:green)    # fermionic Green2
    @test g6.factor == 1
    @test external(g6) == OperatorProduct(V6)[[1, 8]]
    standardize_order!(g6)
    @test g6.factor == -1
    gg6 = feynman_diagram(V6, [[1, 3], [2, 6], [4, 7], [5, 8]], external=[8, 1], type=:green)
    @test gg6.factor ≈ g6.factor

    V7 = [𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6), 𝑓⁻(7)]
    g7 = feynman_diagram(V7, [[1, 5], [3, 6], [4, 7]], external=[2, 7])     # sigma*G
    @test g7.factor == 1

    V8 = [𝑓⁺(1), 𝑓⁺(2)𝑓⁻(3)𝜙(4), 𝑓⁺(5)𝑓⁻(6)𝜙(7), 𝑓⁺(8)𝑓⁻(9)𝜙(10), 𝑓⁻(11), 𝑓⁺(12)𝑓⁻(13)𝜙(14)]
    g8 = feynman_diagram(V8, [[1, 3], [2, 6], [4, 14], [5, 13], [7, 10], [8, 11]], external=[1, 9, 11, 12])
    @test g8.factor == -1
    gg8 = feynman_diagram(V8, [[1, 3], [2, 6], [4, 14], [5, 13], [7, 10], [8, 11]], external=[9, 1, 11, 12])
    @test gg8.factor == -1


    # multi-oeprators (>2) contractions
    Vm = [𝑓⁺(1)𝑓⁻(2)𝑏⁺(3), 𝜙(4)𝑓⁺(5)𝑓⁻(6), 𝑓(7)𝑏⁻(8)𝜙(9)]
    gm = feynman_diagram(Vm, [[1, 2, 3, 8], [4, 5, 6, 9]])
    @test vertices(gm) == Vm
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


    # construct Feynman diagram from Graphs
    # g1 = ComputationalGraphs.propagator(𝑓⁺(1)𝑓⁻(2),)
    # g2 = ComputationalGraphs.propagator(𝑓⁺(2)𝑓⁻(1),)
    # g = feynman_diagram([g1, g2], [1, 2, 2, 1]; external=[1, 2]) #build Feynman diagram from Graphs with Wick's contractions
    # @test external(g) == [external(g1); external(g2)]
    # @test isempty(internal_vertices(g))

    # g = feynman_diagram([g1, g2], [1, 2, 2, 1]; external=[1, 2]) #build Feynman diagram from Graphs with topology
    # @test external(g) == [external(g1); external(g2)]
    # @test isempty(internal_vertices(g))
end