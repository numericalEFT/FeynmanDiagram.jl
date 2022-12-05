module ComputationalGraphs

using AbstractTrees
using Printf, PyCall, DataFrames, Parameters

@enum TwoBodyChannel Alli = 1 PHr PHEr PPr AnyChan
@enum Permutation Di = 1 Ex DiEx

# export TwoBodyChannel, Alli, PHr, PHEr, PPr, AnyChan
# export Permutation, Di, Ex, DiEx

Base.length(r::TwoBodyChannel) = 1
Base.iterate(r::TwoBodyChannel) = (r, nothing)
function Base.iterate(r::TwoBodyChannel, ::Nothing) end

Base.length(r::Permutation) = 1
Base.iterate(r::Permutation) = (r, nothing)
function Base.iterate(r::Permutation, ::Permutation) end

const INL, OUTL, INR, OUTR = 1, 2, 3, 4

import ..QuantumOperators: QuantumOperator, QuantumExpr, 𝑓⁻, 𝑓⁺, 𝑓, 𝑏⁻, 𝑏⁺, 𝜙, iscreation, isfermionic, parity, correlator_order, correlator_order!

include("common.jl")
export labelreset, parity, parity_old

include("graph.jl")
export Graph
export feynman_diagram, contractions_to_edges, propagator, standardize_order!
export is_external, is_internal, external_vertices, internal_vertices, vertices
# export 𝐺ᶠ, 𝐺ᵇ, 𝐺ᵠ, 𝑊, Green2, Interaction
# export Coupling_yukawa, Coupling_phi3, Coupling_phi4, Coupling_phi6

# include("tree.jl")
# include("operation.jl")

include("io.jl")
# plot_tree

# include("eval.jl")
# include("optimize.jl")

# export addSubDiagram!
# export evalDiagTree!
# export evalDiagTreeKT!
# export Operator, Sum, Prod
# export uidreset
# export toDataFrame, mergeby, plot_tree

end