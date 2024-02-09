module FeynmanDiagram
using Random, LinearAlgebra, Parameters, AbstractTrees, RuntimeGeneratedFunctions

macro todo()
    return :(error("Not yet implemented!"))
end

include("quantum_operator/QuantumOperators.jl")

using .QuantumOperators
export QuantumOperators
export QuantumOperator, OperatorProduct, isfermionic
export 𝑓⁻, 𝑓⁺, 𝑓, 𝑏⁻, 𝑏⁺, 𝜙
# export 𝑓⁻ₑ, 𝑓⁺ₑ, 𝑓ₑ, 𝑏⁻ₑ, 𝑏⁺ₑ, 𝜙ₑ
export fermionic_annihilation, fermionic_creation, majorana
export bosonic_annihilation, bosonic_creation, real_classic
export correlator_order, normal_order

include("computational_graph/ComputationalGraph.jl")
using .ComputationalGraphs
export ComputationalGraphs
export labelreset, parity
# export AbstractOperator, Prod, Sum

export AbstractGraph, AbstractOperator
export Graph, FeynmanGraph, FeynmanProperties
export plot_tree

export isequiv, drop_topology, is_external, is_internal, diagram_type, orders, vertices, topology
export external_legs, external_indices, external_operators, external_labels
export multi_product, linear_combination, feynman_diagram, propagator, interaction, external_vertex

export haschildren, onechild, isleaf, isbranch, ischain, has_zero_subfactors, eldest, count_operation
export relabel!, standardize_labels!, replace_subgraph!, merge_linear_combination!, merge_multi_product!, remove_zero_valued_subgraphs!
export relabel, standardize_labels, replace_subgraph, merge_linear_combination, merge_multi_product, remove_zero_valued_subgraphs
export open_parenthesis, open_parenthesis!, flatten_prod!, flatten_prod, flatten_sum!, flatten_sum, flatten_chains!, flatten_chains
export optimize!, optimize, flatten_all_chains!, merge_all_linear_combinations!, merge_all_multi_products!, remove_all_zero_valued_subgraphs!, remove_duplicated_leaves!, remove_duplicated_nodes!

include("TaylorSeries/TaylorSeries.jl")
using .Taylor
export Taylor

include("utility.jl")
using .Utility
export Utility
export taylorexpansion!

include("frontend/frontends.jl")
using .FrontEnds
export FrontEnds
export LabelProduct
export Filter, Wirreducible, Girreducible, NoBubble, NoHartree, NoFock, Proper, DirectOnly
using .GV
export GV, diagdictGV, diagdictGV_ver4, leafstates
using .Parquet
export Parquet, diagdict_parquet, diagdict_parquet_extraAD

include("backend/compiler.jl")
using .Compilers
export Compilers

end
