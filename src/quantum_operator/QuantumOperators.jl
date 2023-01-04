module QuantumOperators

include("operator.jl")

export QuantumOperator

include("expression.jl")

export OperatorProduct, isfermionic, iscreation, isannihilation
export 𝑓⁻, 𝑓⁺, 𝑓, 𝑏⁻, 𝑏⁺, 𝜙
# export 𝑓⁻ₑ, 𝑓⁺ₑ, 𝑓ₑ, 𝑏⁻ₑ, 𝑏⁺ₑ, 𝜙ₑ
export fermionic_annihilation, fermionic_creation, majorana
export bosonic_annihilation, bosonic_creation, real_classic
export normal_order, correlator_order

end