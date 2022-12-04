@testset "QuantumExpr" begin
    @test 𝑓(1) == QuantumExpr(QuantumOperator(:f, 1))
    @test QuantumOperators.isfermionic(𝑓(1)[1])
    @test QuantumOperators.isfermionic(𝑓⁺(1)[1])
    @test QuantumOperators.isfermionic(𝑓⁻(1)[1])
    @test QuantumOperators.iscreation(𝑓⁺(1)[1])
    @test QuantumOperators.iscreation(𝑏⁺(1)[1])

    qe1 = QuantumExpr([QuantumOperator(:f⁺, 1), QuantumOperator(:f⁻, 2), QuantumOperator(:ϕ, 3)])
    qe2 = QuantumExpr([QuantumOperator(:f⁺, 1), QuantumOperator(:f⁻, 2),
        QuantumOperator(:ϕ, 3), QuantumOperator(:b⁻, 4)])
    qe3 = QuantumExpr([QuantumOperator(:b⁻, 4), QuantumOperator(:f⁺, 1), QuantumOperator(:f⁻, 2),
        QuantumOperator(:ϕ, 3)])
    @test QuantumOperator(:f⁺, 1) * QuantumOperator(:f⁻, 2) * QuantumOperator(:ϕ, 3) == qe1
    @test 𝑓⁺(1)𝑓⁻(2)𝜙(3) == qe1
    @test qe1 * 𝑏⁻(4) == qe2
    @test qe1 * QuantumOperator(:b⁻, 4) == qe2
    @test QuantumOperator(:b⁻, 4) * qe1 == qe3
    @test QuantumExpr(qe1) == qe1.operators

end

@testset "Parity" begin
    # P = (1) => sgn(P) = 1
    p1 = [1]
    @test QuantumOperators.parity(p1) == 1
    @test QuantumOperators.parity_old(p1) == 1

    # P = (2 3 1 5 6 4) = (1 2 3) (4 5 6) => sgn(P) = 1
    p2 = [2, 3, 1, 5, 6, 4]
    @test QuantumOperators.parity(p2) == 1
    @test QuantumOperators.parity_old(p2) == 1

    # P = (3 4 1 2) = (1 3) (2 4) => sgn(P) = 1
    p3 = [3, 4, 1, 2]
    @test QuantumOperators.parity(p3) == 1
    @test QuantumOperators.parity_old(p3) == 1

    # P = (3 5 1 2 4 6 7) = (1 3) (2 5 4) (6) (7) => sgn(P) = -1
    p4 = [3, 5, 1, 2, 4, 6, 7]
    @test QuantumOperators.parity(p4) == -1
    @test QuantumOperators.parity_old(p4) == -1
end