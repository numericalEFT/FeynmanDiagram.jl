module DiagTreeNew
using AbstractTrees
using Printf, PyCall, DataFrames

@enum Channel I = 1 T U S ITUS
@enum Permutation Di = 1 Ex DiEx

export Channel, I, T, U, S, ITUS
export Permutation, Di, Ex, DiEx

Base.length(r::Channel) = 1
Base.iterate(r::Channel) = (r, nothing)
function Base.iterate(r::Channel, ::Nothing) end

Base.length(r::Permutation) = 1
Base.iterate(r::Permutation) = (r, nothing)
function Base.iterate(r::Permutation, ::Permutation) end

include("common.jl")
include("traits.jl")
include("diagram.jl")
include("diagram_io.jl")
export Diagram
export addSubDiagram!
export toDataFrame
export evalDiagNode!
export evalDiagTree!
export Operator, Sum, Prod
export DiagramId, Ver4Id, Ver3Id, GreenId, SigmaId, PolarId, InteractionId

end