module GV

import ..QuantumOperators as Op
import ..ComputationalGraphs as IR
import ..ComputationalGraphs: FeynmanGraph
import ..ComputationalGraphs: Graph
import ..ComputationalGraphs: _dtype
import ..Taylor: set_variables
import ..FrontEnds
import ..FrontEnds: LabelProduct
import ..FrontEnds: Filter, NoHartree, NoFock, DirectOnly
import ..FrontEnds: Wirreducible  #remove all polarization subdiagrams
import ..FrontEnds: Girreducible  #remove all self-energy inseration
import ..FrontEnds: NoBubble  # true to remove all bubble subdiagram
import ..FrontEnds: Proper  #ver4, ver3, and polarization diagrams may require to be irreducible along the transfer momentum/frequency
import ..FrontEnds: Response, Composite, ChargeCharge, SpinSpin, UpUp, UpDown
import ..FrontEnds: AnalyticProperty, Instant, Dynamic
import ..FrontEnds: TwoBodyChannel, Alli, PHr, PHEr, PPr, AnyChan
import ..FrontEnds: DiagramId, Ver4Id, Ver3Id, GreenId, SigmaId, PolarId, GenericId, BareGreenId, BareInteractionId

using AbstractTrees

export diagsGV, diagsGV_ver4

include("GV_diagrams/readfile.jl")


"""
    function diagsGV(type::Symbol, order::Int, GOrder::Int=0, VerOrder::Int=0;
        labelProd::Union{Nothing,LabelProduct}=nothing, spinPolarPara::Float64=0.0,
        tau_labels::Union{Nothing,Vector{Int}}=nothing, filter::Vector{Filter}=[NoHartree])

    Generates a `Vector{FeynmanGraph}`: the given-`type` diagrams with static interactions of a given order, where the actual order of diagrams equals to `order + VerOrder + 2 * GOrder`.
    Generates a `LabelProduct`: `labelProd` with inputs `tau_labels` and all the possible momenta-loop basis. 
    Generates external tau labels Vector{Vector{Int}}. The i-th labels (Vector{Int}) corresponds to the i-th `FeynmanGraph` in `Vector{FeynmanGraph}`.

# Arguments:
- `type` (Symbol): The type of the diagrams, including `:spinPolar`, `:chargePolar`, `:sigma`, `:green`, or `:freeEnergy`.
- `order` (Int): The order of the diagrams without counterterms.
- `GOrder` (Int): The order of self-energy counterterms.
- `VerOrder` (Int): The order of interaction counterterms.
- `labelProd` (Union{Nothing,LabelProduct}=nothing, optional): The initial cartesian QuantumOperator.label product (defaults to `nothing`).
- `spinPolarPara` (Float64, optional): The spin-polarization parameter (n_up - n_down) / (n_up + n_down) (defaults to `0.0`).
- `tau_labels`(Union{Nothing, Vector{Int}}, optional): The labels for the discrete time of each vertex. (defaults to `nothing`).

# Returns
A tuple `(diagrams, fermi_labelProd, bose_labelProd, extT_labels)` where 
- `diagrams` is a `Vector{FeynmanGraph}` object representing the diagrams, 
- `labelProd` is a `LabelProduct` object containing the labels for the leaves of graphs, 
- `extT_labels` is a `Vector{Union{Tuple,Vector{Int}}}` object containing the external tau labels for each `FeynmanGraph` in `diagrams`.
"""
function diagsGV(type::Symbol, order::Int, GOrder::Int, VerOrder::Int;
    labelProd::Union{Nothing,LabelProduct}=nothing, spinPolarPara::Float64=0.0,
    tau_labels::Union{Nothing,Vector{Int}}=nothing)
    if type == :spinPolar
        filename = string(@__DIR__, "/GV_diagrams/groups_spin/Polar$(order)_$(VerOrder)_$(GOrder).diag")
    elseif type == :chargePolar
        filename = string(@__DIR__, "/GV_diagrams/groups_charge/Polar$(order)_$(VerOrder)_$(GOrder).diag")
    elseif type == :sigma
        filename = string(@__DIR__, "/GV_diagrams/groups_sigma/Sigma$(order)_$(VerOrder)_$(GOrder).diag")
    elseif type == :green
        filename = string(@__DIR__, "/GV_diagrams/groups_green/Green$(order)_$(VerOrder)_$(GOrder).diag")
    elseif type == :freeEnergy
        filename = string(@__DIR__, "/GV_diagrams/groups_free_energy/FreeEnergy$(order)_$(VerOrder)_$(GOrder).diag")
    else
        error("no support for $type diagram")
    end
    # println("Reading ", filename)

    if isnothing(labelProd)
        return read_diagrams(filename; tau_labels=tau_labels, diagType=type, spinPolarPara=spinPolarPara)
    else
        return read_diagrams(filename; labelProd=labelProd, diagType=type, spinPolarPara=spinPolarPara)
    end
end

function diagsGV(type::Symbol, order::Int; spinPolarPara::Float64=0.0, filter::Vector{Filter}=[NoHartree])
    if type == :spinPolar
        filename = string(@__DIR__, "/GV_diagrams/groups_spin/Polar$(order)_0_0.diag")
    elseif type == :chargePolar
        filename = string(@__DIR__, "/GV_diagrams/groups_charge/Polar$(order)_0_0.diag")
    elseif type == :sigma
        filename = string(@__DIR__, "/GV_diagrams/groups_sigma/Sigma$(order)_0_0.diag")
    elseif type == :green
        filename = string(@__DIR__, "/GV_diagrams/groups_green/Green$(order)_0_0.diag")
    elseif type == :freeEnergy
        filename = string(@__DIR__, "/GV_diagrams/groups_free_energy/FreeEnergy$(order)_0_0.diag")
    else
        error("no support for $type diagram")
    end

    return read_diagrams(filename, type, filter=filter, spinPolarPara=spinPolarPara)
end

"""
    function diagsGV_ver4(order::Int; spinPolarPara::Float64=0.0, filter::Vector{Filter}=[NoHartree], channels=[PHr, PHEr, PPr, Alli])

    Generates a `Vector{Graph}`: the 4-point vertex diagrams with static interactions of a given order.

# Arguments:
- `order` (Int): The order of the diagrams without counterterms.
- `spinPolarPara` (Float64, optional): The spin-polarization parameter (n_up - n_down) / (n_up + n_down) (defaults to `0.0`).
- `channels` (optional): The channels for the diagrams, defaults to `[PHr, PHEr, PPr, Alli]`.
- `filter` (optional): Filter criteria for the diagrams, defaults to `[NoHartree]`.
"""
function diagsGV_ver4(order::Int; spinPolarPara::Float64=0.0, channels=[PHr, PHEr, PPr, Alli], filter::Vector{Filter}=[NoHartree])
    if channels == [Alli]
        filename = string(@__DIR__, "/GV_diagrams/groups_vertex4/Vertex4I$(order)_0_0.diag")
        return read_vertex4diagrams(filename, spinPolarPara=spinPolarPara, channels=channels, filter=filter)
    else
        filename = string(@__DIR__, "/GV_diagrams/groups_vertex4/Vertex4$(order)_0_0.diag")
        return read_vertex4diagrams(filename, spinPolarPara=spinPolarPara, channels=channels, filter=filter)
    end
end

end