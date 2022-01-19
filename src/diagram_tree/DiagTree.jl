module DiagTreeNew
using AbstractTrees
using Printf, PyCall, DataFrames

@enum Channel I = 1 T U S Ts Us Ic Tc Uc Sc Tsc Usc
@enum DiEx Di = 1 Ex

export Channel, I, T, U, S
export DiEx, Di, Ex

Base.length(r::Channel) = 1
Base.iterate(r::Channel) = (r, nothing)
function Base.iterate(r::Channel, ::Nothing) end

Base.length(r::DiEx) = 1
Base.iterate(r::DiEx) = (r, nothing)
function Base.iterate(r::DiEx, ::DiEx) end

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