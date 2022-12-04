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

function isfermionic(operator::QuantumOperator)
    operator.operator in [:f⁺, :f⁻, :f]
end

function iscreation(operator::QuantumOperator)
    operator.operator in [:f⁺, :b⁺]
end
