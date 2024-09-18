@testset "Contractions" begin
    # Test 1: Scalar fields with Wick crossings, parity = +1
    vertices1 = [𝜙(1)𝜙(2), 𝜙(3)𝜙(4)𝜙(5)𝜙(6), 𝜙(7)𝜙(8)]
    parity1, ind1, edges1 = contractions_to_edges(vertices1, [1, 2, 3, 4, 1, 3, 4, 2])
    ops = reduce(*, vertices1)
    @test ind1 == [[1, 5], [2, 8], [3, 6], [4, 7]]
    @test Set(edges1) == Set([(ops[1], ops[5]), (ops[2], ops[8]), (ops[3], ops[6]), (ops[4], ops[7])])
    @test parity1 == 1

    # Test 2: Bosons with Wick crossings, parity = +1
    vertices2 = [𝑏⁺(1)𝑏⁺(2)𝑏⁻(3), 𝑏⁻(4)𝑏⁺(5), 𝑏⁻(6)𝑏⁺(7)𝑏⁻(8)]
    parity2, ind2, edges2 = contractions_to_edges(vertices2, [1, 2, 3, 4, 3, 1, 4, 2])
    ops = reduce(*, vertices2)
    @test ind2 == [[1, 6], [2, 8], [3, 5], [4, 7]]
    @test Set(edges2) == Set([(ops[1], ops[6]), (ops[2], ops[8]), (ops[3], ops[5]), (ops[4], ops[7])])
    @test parity2 == 1

    # Test 3: Indistinguishable Majoranas with no Wick crossings, parity = +1
    vertices3 = [𝑓(1)𝑓(1)𝑓(1)𝑓(1)𝑓(1)𝑓(1)𝑓(1)𝑓(1),]
    parity3, ind3, edges3 = contractions_to_edges(vertices3, [1, 2, 3, 4, 4, 3, 2, 1])
    ops = reduce(*, vertices3)
    @test ind3 == [[1, 8], [2, 7], [3, 6], [4, 5]]
    @test Set(edges3) == Set([(ops[1], ops[8]), (ops[2], ops[7]), (ops[3], ops[6]), (ops[4], ops[5])])
    # P = (1 8 2 7 3 6 4 5) = (1)(5 3 2 8)(4 7)(6) => parity = +1
    @test parity3 == 1

    # Test 4: Fermions with Wick crossings, parity = -1
    vertices4 = [𝑓⁺(1)𝑓⁻(2), 𝑓⁺(5)𝑓⁺(6)𝑓⁻(7)𝑓⁻(8), 𝑓⁺(3)𝑓⁻(4),]
    parity4, ind4, edges4 = contractions_to_edges(vertices4, [1, 2, 2, 3, 1, 4, 4, 3])
    ops = reduce(*, vertices4)
    @test ind4 == [[1, 5], [2, 3], [4, 8], [6, 7]]
    @test Set(edges4) == Set([(ops[1], ops[5]), (ops[2], ops[3]), (ops[4], ops[8]), (ops[6], ops[7])])
    # P = (1 5 2 3 4 8 6 7) = (1)(2 5 4)(3)(6 8)(7) => parity = -1
    @test parity4 == -1

    # Test 5: Mixed bosonic/classical/fermionic operators, parity = -1
    vertices5 = [𝑏⁺(1)𝑓⁺(2)𝜙(3), 𝑓⁻(4)𝑓⁻(5), 𝑏⁻(6)𝑓⁺(7)𝜙(8)]
    ops = reduce(*, vertices5)
    parity5, ind5, edges5 = contractions_to_edges(vertices5, [1, 2, 3, 2, 4, 1, 4, 3])
    @test ind5 == [[1, 6], [2, 4], [3, 8], [5, 7]]
    @test Set(edges5) == Set([(ops[1], ops[6]), (ops[2], ops[4]), (ops[3], ops[8]), (ops[5], ops[7])])
    # Flattened fermionic edges: [2, 4, 5, 7]
    # => P = (1 2 3 4) = (1)(2)(3 4) => parity = 1
    @test parity5 == 1
end

@testset "feynman_diagram from Wick" begin
    # construct Feynman diagram from FeynmanGraphs
    g1 = ComputationalGraphs.propagator(𝑓⁺(1)𝑓⁻(2),)
    g2 = ComputationalGraphs.propagator(𝑓⁺(2)𝑓⁻(1),)
    g = feynman_diagram([g1, g2], [1, 2, 2, 1]; external=[1, 2]) #build Feynman diagram from FeynmanGraphs with Wick's contractions
    @test external(g) == [external(g1); external(g2)]
    @test isempty(internal_vertices(g))

    g = feynman_diagram([g1, g2], [1, 2, 2, 1]; external=[1, 2]) #build Feynman diagram from FeynmanGraphs with topology
    @test external(g) == [external(g1); external(g2)]
    @test isempty(internal_vertices(g))
end
