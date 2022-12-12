module ComputationalGraphs

using AbstractTrees
using Printf, PyCall, DataFrames

import ..QuantumOperators: QuantumOperator, OperatorProduct, 𝑓⁻, 𝑓⁺, 𝑓, 𝑏⁻, 𝑏⁺, 𝜙, iscreation, isannihilation, isfermionic, isghost, parity, normal_order, correlator_order

include("common.jl")
export labelreset

include("graph.jl")
export Graph, isequiv
export feynman_diagram, contractions_to_edges, propagator, standardize_order!
export is_external, is_internal, external, vertices, real_legs, fake_legs
# export 𝐺ᶠ, 𝐺ᵇ, 𝐺ᵠ, 𝑊, Green2, Interaction

# include("tree.jl")
# include("operation.jl")

include("io.jl")
# plot_tree

# include("eval.jl")
# include("optimize.jl")

end