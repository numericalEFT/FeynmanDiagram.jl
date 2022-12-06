"""
    struct QuantumOperator

    struct of a quantum operator.

# Members:
- `operator::Symbol`  symbol of quantum operator, supports :f⁺, :f⁻, :f, :b⁺, :b⁻, :ϕ
- `label::Int`  label of the operator indices. It could represent spacetime, spin, momentum, flavor, etc.
"""
struct QuantumOperator
    operator::Symbol
    label::Int
    function QuantumOperator(operator::Symbol, label::Int)
        @assert label > 0
        return new(operator, label)
    end
end
Base.isequal(a::QuantumOperator, b::QuantumOperator) = ((a.operator == b.operator) && (a.label == b.label))
Base.:(==)(a::QuantumOperator, b::QuantumOperator) = Base.isequal(a, b)

function Base.show(io::IO, o::QuantumOperator)
    print(io, "$(String(o.operator))($(o.label))")
end

Base.show(io::IO, ::MIME"text/plain", o::QuantumOperator) = Base.show(io, o)

"""
    Base.adjoint(operator::QuantumOperator)

    Return the conjuated quantum operator of `operator`.
"""
function Base.adjoint(operator::QuantumOperator)
    if operator.operator in [:f, :ϕ]
        return operator
    elseif operator.operator == :f⁺
        return QuantumOperator(:f⁻, operator.label)
    elseif operator.operator == :f⁻
        return QuantumOperator(:f⁺, operator.label)
    elseif operator.operator == :b⁺
        return QuantumOperator(:b⁻, operator.label)
    elseif operator.operator == :b⁻
        return QuantumOperator(:b⁺, operator.label)
    end
end

"""
    function isfermionic(operator::QuantumOperator)

    Check if `operator` is a fermionic operator.
"""
function isfermionic(operator::QuantumOperator)
    operator.operator in [:f⁺, :f⁻, :f]
end

function iscreation(operator::QuantumOperator)
    operator.operator in [:f⁺, :b⁺]
end

function isannihilation(operator::QuantumOperator)
    operator.operator in [:f⁻, :b⁻]
end