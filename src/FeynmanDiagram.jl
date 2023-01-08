module FeynmanDiagram
using Random, LinearAlgebra, Parameters

macro todo()
    return :(error("Not yet implemented!"))
end

@enum DiagramType begin
    VacuumDiag         #vaccum diagram for the free energy
    SigmaDiag          #self-energy
    GreenDiag          #green's function
    PolarDiag          #polarization
    Ver3Diag           #3-point vertex function
    Ver4Diag           #4-point vertex function
    GnDiag             #n-point Green's function
    GcDiag             #n-point connected Green's function
end
Base.length(r::DiagramType) = 1
Base.iterate(r::DiagramType) = (r, nothing)
function Base.iterate(r::DiagramType, ::Nothing) end

abstract type DiagType end
abstract type Vacuum <: DiagType end
abstract type Tadpole <: DiagType end
abstract type FermiPropagator <: DiagType end
abstract type BosePropagator <: DiagType end
abstract type FermiSelfEnergy <: DiagType end
abstract type BoseSelfEnergy <: DiagType end
abstract type VertexDiag <: DiagType end
abstract type GncDiag <: DiagType end
abstract type GndDiag <: DiagType end

@enum Filter begin
    Wirreducible  #remove all polarization subdiagrams
    Girreducible  #remove all self-energy inseration
    NoHartree
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

export SigmaDiag, PolarDiag, Ver3Diag, Ver4Diag, GreenDiag
export VacuumDiag, GnDiag, GcDiag
export Wirreducible, Girreducible, NoBubble, NoHartree, NoFock, Proper
export Response, ChargeCharge, SpinSpin, UpUp, UpDown
export AnalyticProperty, Instant, Dynamic, D_Instant, D_Dynamic

export DiagType
export FermiPropagator, BosePropagator, FermiSelfEnergy, BoseSelfEnergy, VertexDiag
export Vacuum, Tadpole, GncDiag, GndDiag

include("common.jl")
export DiagPara, DiagParaF64
export Interaction, interactionTauNum, innerTauNum

include("common_new.jl")
export DiagramPara, DiagramParaF64
# export Interaction, interactionTauNum, innerTauNum

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
export labelreset, parity, parity_old
# export AbstractOperator, Prod, Sum
export Graph, isequiv, linear_combination
# export GraphType, Interaction, ExternalVertex, Propagator, SelfEnergy, VertexDiag, GreenDiag, GenericDiag
export feynman_diagram, propagator, interaction, external_vertex
# export standardize_order!
export is_external, is_internal, vertices, external
export external_labels
# export reducibility, connectivity
# export 𝐺ᶠ, 𝐺ᵇ, 𝐺ᵠ, 𝑊, Green2, Interaction
# export Coupling_yukawa, Coupling_phi3, Coupling_phi4, Coupling_phi6

include("backend/compiler.jl")
using .Compilers
export Compilers

include("frontend/frontends.jl")
using .FrontEnds
export FrontEnds

include("diagram_tree/DiagTree.jl")
using .DiagTree
export DiagTree
export TwoBodyChannel, Alli, PHr, PHEr, PPr, AnyChan
export Permutation, Di, Ex, DiEx
export Diagram, addSubDiagram!, toDataFrame
export evalDiagNode!, evalDiagTree!, evalDiagTreeKT!
export Operator, Sum, Prod
export DiagramId, GenericId, Ver4Id, Ver3Id, GreenId, SigmaId, PolarId
export PropagatorId, BareGreenId, BareInteractionId
export BareGreenNId, BareHoppingId, GreenNId, ConnectedGreenNId
export uidreset, toDataFrame, mergeby, plot_tree

include("parquet_builder/parquet.jl")
using .Parquet
export Parquet
export ParquetBlocks

include("strong_coupling_expansion_builder/strong_coupling_expansion")
using .SCE
export SCE
export Gn

include("expression_tree/ExpressionTree.jl")
using .ExprTree
export ExprTree
export Component, ExpressionTree
# export Propagator
export addpropagator!, addnode!
export setroot!, addroot!
export evalNaive, showTree

##################### precompile #######################
# precompile as the final step of the module definition:
if ccall(:jl_generating_output, Cint, ()) == 1   # if we're precompiling the package
    let
        para = DiagParaF64(type=Ver4Diag, innerLoopNum=2, hasTau=true)
        # ver4 = Parquet.vertex4(para)  # this will force precompilation
        ver4 = Parquet.build(para)  # this will force precompilation

        mergeby(ver4, [:response])
        mergeby(ver4.diagram)
        mergeby(ver4.diagram, [:response]; idkey=[:extT, :response])

        para = DiagParaF64(type=SigmaDiag, innerLoopNum=2, hasTau=true)
        Parquet.build(para)  # this will force precompilation
        para = DiagParaF64(type=GreenDiag, innerLoopNum=2, hasTau=true)
        Parquet.green(para)  # this will force precompilation
        para = DiagParaF64(type=PolarDiag, innerLoopNum=2, hasTau=true)
        # Parquet.polarization(para)  # this will force precompilation
        Parquet.build(para)  # this will force precompilation
        para = DiagParaF64(type=Ver3Diag, innerLoopNum=2, hasTau=true)
        # Parquet.vertex3(para)  # this will force precompilation
        Parquet.build(para)  # this will force precompilation

        DiagTree.removeHartreeFock!(ver4.diagram)
        DiagTree.derivative(ver4.diagram, BareGreenId)
        DiagTree.derivative(ver4.diagram, BareInteractionId)
        # DiagTree.removeHartreeFock!(ver4.diagram)
        ExprTree.build(ver4.diagram)
    end
end

end