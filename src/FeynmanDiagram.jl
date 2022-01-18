module FeynmanDiagram
using Random, LinearAlgebra, Parameters

@enum DiagramType begin
    SigmaDiag          #self-energy
    GreenDiag          #green's function
    PolarDiag          #polarization
    Ver3Diag           #3-point vertex function
    Ver4Diag           #4-point vertex function
end
Base.length(r::DiagramType) = 1
Base.iterate(r::DiagramType) = (r, nothing)
function Base.iterate(r::DiagramType, ::Nothing) end

@enum Filter begin
    Wirreducible  #remove all polarization subdiagrams
    Girreducible  #remove all self-energy inseration
    NoHatree
    NoFock
    NoBubble  # true to remove all bubble subdiagram
    Proper  #ver4, ver3, and polarization diagrams may require to be irreducible along the transfer momentum/frequency
end

Base.length(r::Filter) = 1
Base.iterate(r::Filter) = (r, nothing)
function Base.iterate(r::Filter, ::Nothing) end

@enum Response begin
    Composite
    ChargeCharge
    SpinSpin
    ProperChargeCharge
    ProperSpinSpin
    UpUp
    UpDown
end

Base.length(r::Response) = 1
Base.iterate(r::Response) = (r, nothing)
function Base.iterate(r::Response, ::Nothing) end

@enum AnalyticProperty begin
    Instant
    Dynamic
    D_Instant #derivative of instant interaction
    D_Dynamic #derivative of the dynamic interaction
end

Base.length(r::AnalyticProperty) = 1
Base.iterate(r::AnalyticProperty) = (r, nothing)
function Base.iterate(r::AnalyticProperty, ::Nothing) end

export SigmaDiag, PolarDiag, Ver3Diag, Ver4Diag
export Wirreducible, Girreducible, NoBubble, NoHatree, Proper
export InteractionName, ChargeCharge, SpinSpin, UpUp, UpDown
export InteractionType, Instant, Dynamic, D_Instant, D_Dynamic

include("common.jl")
export GenericPara, Interaction

include("diagram_tree/DiagTree.jl")
export DiagTreeNew
export DiagramId, Diagram, add_subdiagram!, toDataFrame

include("diagramTree/DiagTree.jl")
export DiagTree
export Component, Diagrams
export addpropagator!, addnode!
export setroot!, addroot!
export evalNaive, showTree

include("parquetBuilder/parquet.jl")
export ParquetNew

include("diagramBuilder/builder.jl")
export Builder
# Write your package code here.

end
