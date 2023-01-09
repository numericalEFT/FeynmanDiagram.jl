module ComputationalGraphs

using AbstractTrees
using Printf, PyCall, DataFrames

macro todo()
    return :(error("Not yet implemented!"))
end

import ..QuantumOperators: QuantumOperator, OperatorProduct, 𝑓⁻, 𝑓⁺, 𝑓, 𝑏⁻, 𝑏⁺, 𝜙, iscreation, isannihilation, isfermionic, parity, normal_order, correlator_order
# import ..QuantumOperators: 𝑓⁻ₑ, 𝑓⁺ₑ, 𝑓ₑ, 𝑏⁻ₑ, 𝑏⁺ₑ, 𝜙ₑ

include("common.jl")
export labelreset
export _dtype
export set_datatype

include("graph.jl")
# export AbstractOperator, Prod, Sum
export Graph, isequiv, linear_combination
# export GraphType, Interaction, ExternalVertex, Propagator, SelfEnergy, VertexDiag, GreenDiag, GenericDiag
export feynman_diagram, propagator, interaction, external_vertex
# export standardize_order!
export is_external, is_internal, vertices, external
export external_labels
# export 𝐺ᶠ, 𝐺ᵇ, 𝐺ᵠ, 𝑊, Green2, Interaction

include("tree_properties.jl")
export haschildren, onechild, isleaf, isbranch, ischain, isfactorless, eldest

# include("operation.jl")
include("graphvector.jl")
include("io.jl")
# plot_tree

# include("eval.jl")
# include("optimize.jl")

include("transform.jl")
export relabel!, standardize_labels!, replace_subgraph!, merge_prodchain_subfactors!, inplace_prod!
export relabel, standardize_labels, replace_subgraph, merge_prodchain_subfactors, inplace_prod
export prune_trivial_unary, merge_prefactors

end
