@testset "Diagram" begin
    # electron gas
    g1 = 𝐺ᶠ(1, 2) * 𝑊(1, 2) * 𝐺ᶠ(2, 1)          # vaccum diagram
    @test g1.couplings[1] == Coupling_yukawa
    @test isempty(g1.external_vertices)
    @test g1.internal_vertices[1] == InternalVertex(1, 0, Coupling_yukawa)
    @test g1.internal_vertices[2] == InternalVertex(2, 0, Coupling_yukawa)
    @test checkVertices(g1)

    g2 = 𝐺ᶠ(1, 1) * 𝑊(1, 2) * 𝐺ᶠ(2, 2)
    @test g2.couplings[1] == Coupling_yukawa
    @test isempty(g2.external_vertices)
    @test g2.internal_vertices[1] == InternalVertex(1, 0, Coupling_yukawa)
    @test g2.internal_vertices[2] == InternalVertex(2, 0, Coupling_yukawa)
    @test checkVertices(g2)

    g3 = 𝐺ᶠ(1, 2) * 𝐺ᶠ(2, 3) * 𝑊(2, 3) * 𝐺ᶠ(3, 4)
    @test g3.couplings[1] == Coupling_yukawa
    @test g3.external_vertices[1] == ExternalVertex(1, 0, CompositeOperator([𝑓dag]))
    @test g3.external_vertices[2] == ExternalVertex(4, 0, CompositeOperator([𝑓]))
    @test g3.internal_vertices[1] == InternalVertex(2, 0, Coupling_yukawa)
    @test g3.internal_vertices[2] == InternalVertex(3, 0, Coupling_yukawa)
    @test checkVertices(g3)

    # phi4 theory
    g4 = 𝐺ᵠ(1, 1) * 𝐺ᵠ(1, 1) * 𝑊(1, Coupling_phi4)
    @test g4.couplings[1] == Coupling_phi4
    @test isempty(g4.external_vertices)
    @test g4.internal_vertices[1] == InternalVertex(1, 0, Coupling_phi4)
    @test checkVertices(g4)
    # @test !checkVertices(𝐺ᵠ(1, 1) * 𝐺ᵠ(1, 1) * 𝐺ᵠ(1, 1) * 𝑊(1, Coupling_phi4))

    g5 = 𝐺ᵠ(1, 2) * 𝑊(2, Coupling_phi4) * 𝐺ᵠ(2, 3) * 𝐺ᵠ(2, 3) * 𝐺ᵠ(2, 3) * 𝑊(3, Coupling_phi4) * 𝐺ᵠ(3, 4)
    @test g5.couplings[1] == Coupling_phi4
    @test g5.external_vertices[1] == ExternalVertex(1, 0, CompositeOperator([ϕ]))
    @test g5.external_vertices[2] == ExternalVertex(4, 0, CompositeOperator([ϕ]))
    @test g5.internal_vertices[1] == InternalVertex(2, 0, Coupling_phi4)
    @test g5.internal_vertices[2] == InternalVertex(3, 0, Coupling_phi4)
    @test checkVertices(g5)
end