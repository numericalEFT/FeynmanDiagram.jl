module QuantumOperators

include("operator.jl")

export QuantumOperator

include("expression.jl")

export OperatorProduct, isfermionic, iscreation, isannihilation, isghost
export 𝑓⁻, 𝑓⁺, 𝑓, 𝑏⁻, 𝑏⁺, 𝜙
export fermionic_annihilation, fermionic_creation, majorana
export bosonic_annihilation, bosonic_creation, real_classic
export normal_order, correlator_order

end