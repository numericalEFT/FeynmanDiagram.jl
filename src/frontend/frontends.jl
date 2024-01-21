module FrontEnds

import ..ComputationalGraphs
import ..ComputationalGraphs: Graph, _dtype
import ..QuantumOperators
import ..Taylor
import ..Utility
using LinearAlgebra

@enum TwoBodyChannel Alli = 1 PHr PHEr PPr AnyChan

@enum Filter begin
    Wirreducible  #remove all polarization subdiagrams
    Girreducible  #remove all self-energy inseration
    NoHartree
    NoFock
    NoBubble  # true to remove all bubble subdiagram
    Proper  #ver4, ver3, and polarization diagrams may require to be irreducible along the transfer momentum/frequency
    DirectOnly # only direct interaction, this can be useful for debug purpose
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
end

Base.length(r::AnalyticProperty) = 1
Base.iterate(r::AnalyticProperty) = (r, nothing)
function Base.iterate(r::AnalyticProperty, ::Nothing) end

function short(name::Response)
    if name == ChargeCharge
        return "cc"
    elseif name == SpinSpin
        return "σσ"
    elseif name == UpUp
        return "↑↑"
    elseif name == UpDown
        return "↑↓"
    else
        @error("$name is not implemented!")
    end
end

function short(type::AnalyticProperty)
    if type == Instant
        return "Ins"
    elseif type == Dynamic
        return "Dyn"
    else
        @error("$type is not implemented!")
    end
end

function symbol(name::Response, type::AnalyticProperty, addition=nothing)
    if isnothing(addition)
        return Symbol("$(short(name))$(short(type))")
    else
        return Symbol("$(short(name))$(short(type))$(addition)")
    end

end

include("diagram_id.jl")
export DiagramId, GenericId, Ver4Id, Ver3Id, GreenId, SigmaId, PolarId, GreenNId, ConnectedGreenNId
export PropagatorId, BareGreenId, BareInteractionId, BareHoppingId, BareGreenNId

include("pool.jl")
export LoopPool

include("LabelProduct.jl")
export LabelProduct

include("parquet/parquet.jl")
export Parquet

include("GV.jl")
export GV

export get_ver4I

const vertex4I_diags = Dict{Int,Vector{Graph}}()
function initialize_vertex4I_diags(; filter=[NoHartree], spinPolarPara::Float64=0.0)
    dict_graphs = GV.diagdictGV_AD(:vertex4I, [(3, 0, 0), (4, 0, 0)], filter=filter, spinPolarPara=spinPolarPara)
    vertex4I_diags[3] = dict_graphs[(3, 0, 0)][1]
    vertex4I_diags[4] = dict_graphs[(4, 0, 0)][1]
end

get_ver4I() = vertex4I_diags

initialize_vertex4I_diags()

end