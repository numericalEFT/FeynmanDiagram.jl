@testset "OperatorProduct" begin
    @test 𝑓(1) == OperatorProduct(QuantumOperator(:f, 1))
    @test isfermionic(𝑓(1)[1])
    @test isfermionic(𝑓⁺(1)[1])
    @test isfermionic(𝑓⁻(1)[1])
    @test QuantumOperators.iscreation(𝑓⁺(1)[1])
    @test QuantumOperators.iscreation(𝑏⁺(1)[1])
    @test (𝑓⁻(1)[1])' == 𝑓⁺(1)[1]

    qe1 = OperatorProduct([QuantumOperator(:f⁺, 1), QuantumOperator(:f⁻, 2), QuantumOperator(:ϕ, 3)])
    qe2 = OperatorProduct([QuantumOperator(:f⁺, 1), QuantumOperator(:f⁻, 2),
        QuantumOperator(:ϕ, 3), QuantumOperator(:f⁻, 4)])
    qe3 = OperatorProduct([QuantumOperator(:b⁻, 4), QuantumOperator(:f⁺, 1), QuantumOperator(:f⁻, 2),
        QuantumOperator(:ϕ, 3)])
    @test QuantumOperator(:f⁺, 1) * QuantumOperator(:f⁻, 2) * QuantumOperator(:ϕ, 3) == qe1
    @test 𝑓⁺(1)𝑓⁻(2)𝜙(3) == qe1
    @test qe1 * 𝑓⁻(4) == qe2
    @test qe1 * QuantumOperator(:f⁻, 4) == qe2
    @test QuantumOperator(:b⁻, 4) * qe1 == qe3
    @test OperatorProduct(qe1) == qe1.operators
    @test !isfermionic(qe1)
    @test isfermionic(qe2)
    @test !isfermionic(qe3)
    @test qe1' == 𝜙(3)𝑓⁺(2)𝑓⁻(1)
    @test qe3' == 𝜙(3)𝑓⁺(2)𝑓⁻(1)𝑏⁺(4)
end

@testset "correlator order" begin
    o1 = 𝑓⁺(1)𝑓⁻(2)𝑓⁺(5)𝑓⁺(6)𝑓⁻(1)𝑓⁻(5)
    sign, o1 = correlator_order(o1)
    @test sign == 1
    @test o1 == 𝑓⁻(1)𝑓⁻(5)𝑓⁻(2)𝑓⁺(6)𝑓⁺(5)𝑓⁺(1)

    o2 = 𝑓⁺(1)𝑓⁻(2)𝑏⁺(1)𝜙(1)𝑓⁺(6)𝑓⁺(5)𝑓⁻(1)𝑓⁻(5)𝑏⁻(1)
    sign, o2 = correlator_order(o2)
    @test sign == -1
    @test o2 == 𝑓⁻(1)𝑏⁻(1)𝑓⁻(5)𝑓⁻(2)𝜙(1)𝑓⁺(6)𝑓⁺(5)𝑏⁺(1)𝑓⁺(1)

    o3 = 𝑓⁺(1)𝑓⁻(2)𝑏⁺(1)𝜙(1)𝑓⁺(3)𝑓⁻(1)𝑓(1)𝑏⁻(1)𝜙(1)
    sign, o3 = correlator_order(o3)
    @test sign == -1
    @test o3 == 𝑓⁻(1)𝑏⁻(1)𝜙(1)𝑓⁻(2)𝑓(1)𝑓⁺(3)𝜙(1)𝑏⁺(1)𝑓⁺(1)
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