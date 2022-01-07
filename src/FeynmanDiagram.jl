module FeynmanDiagram
using Random, LinearAlgebra, Parameters

@enum DiagramType begin
    SigmaDiag          #self-energy
    GreenDiag          #green's function
    PolarDiag          #polarization
    Ver3Diag           #3-point vertex function
    Ver4Diag           #4-point vertex function
end

@enum Filter begin
    Wirreducible  #remove all polarization subdiagrams
    Girreducible  #remove all self-energy inseration
    NoHatree
    NoFock
    NoBubble  # true to remove all bubble subdiagram
    Proper  #ver4, ver3, and polarization diagrams may require to be irreducible along the transfer momentum/frequency
end

@enum InteractionName begin
    Composite
    ChargeCharge
    SpinSpin
    ProperChargeCharge
    ProperSpinSpin
end

@enum InteractionType begin
    Instant
    Dynamic
    D_Instant #derivative of instant interaction
    D_Dynamic #derivative of the dynamic interaction
end

export SigmaDiag, PolarDiag, Ver3Diag, Ver4Diag
export Wirreducible, Girreducible, NoBubble, NoHatree, Proper
export InteractionName, ChargeCharge, SpinSpin
export InteractionType, Instant, Dynamic, D_Instant, D_Dynamic

include("parameter.jl")
export GenericPara, Interaction

include("diagramTree/DiagTree.jl")
export DiagTree

include("diagramBuilder/builder.jl")
export Builder
# Write your package code here.

end
