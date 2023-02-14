module GW

import ..QuantumOperators as Op
import ..ComputationalGraphs as IR
import ..ComputationalGraphs: Graph
import ..ComputationalGraphs: _dtype
using ..FrontEnds

export PolarEachOrder, PolarDiagrams

include("GW_diagrams/readfile.jl")

"""
    function PolarEachOrder(type::Symbol, order::Int, VerOrder::Int=0, SigmaOrder::Int=0; loopPool::Union{LoopPool,Nothing}=nothing,
        tau_labels::Union{Nothing,Vector{Int}}=nothing, GTypes::Union{Nothing,Vector{Int}}=nothing, WTypes::Union{Nothing,Vector{Int}}=nothing)
 
    Generates a `Graph`: the polarization diagrams of a given order, where the actual order of diagrams equals to `order + VerOrder + 2 * SigmaOrder`.
    Generates fermionic/bosonic `LabelProduct`: `fermi_labelProd`/`bose_labelProd` with inputs `tau_labels`, `GTypes`/`WTypes`, and updated `loopPool`. 

# Arguments:
- `type` (Symbol): The type of the diagrams, either `:spin` or `:charge`.
- `order` (Int): The order of the diagrams without counterterms.
- `VerOrder` (Int, optional): The order of interaction counterterms (defaults to 0).
- `SigmaOrder` (Int, optional): The order of self-energy counterterms (defaults to 0).
- `dim` (Int, optional): The dimension of the system (defaults to 3).
- `loopPool` (Union{LoopPool,Nothing}=nothing, optional): The initial pool of loop momenta (defaults to `nothing`).
- `tau_labels`(Union{Nothing, Vector{Int}}, optional): The labels for the discrete time of each vertex. (defaults to `nothing`).
- `GTypes`: The types of fermion propagators `G` in the diagrams (defaults to `collect(0:SigmaOrder)`).
- `WTypes`: The types of boson interaction `W` in the diagrams (defaults to `collect(0:VerOrder)`).

# Returns
A tuple `(diagrams, fermi_labelProd, bose_labelProd)` where 
- `diagrams` is a `Graph` objects representing the diagrams, 
- `fermi_labelProd` is a `LabelProduct` object containing the labels for the fermionic `G` objects in the diagrams, 
- `bose_labelProd` is a `LabelProduct` object containing the labels for the bosonic `W` objects in the diagrams.
"""
function PolarEachOrder(type::Symbol, order::Int, VerOrder::Int=0, SigmaOrder::Int=0; dim::Int=3, loopPool::Union{LoopPool,Nothing}=nothing,
    tau_labels::Union{Nothing,Vector{Int}}=nothing, GTypes::Union{Nothing,Vector{Int}}=nothing, WTypes::Union{Nothing,Vector{Int}}=nothing)
    if type == :spin
        filename = string(@__DIR__, "/GW_diagrams/groups_spin/Polar$(order)_$(VerOrder)_$(SigmaOrder).diag")
    elseif type == :charge
        filename = string(@__DIR__, "/GW_diagrams/groups_charge/Polar$(order)_$(VerOrder)_$(SigmaOrder).diag")
    end

    if isnothing(GTypes)
        GTypes = collect(0:SigmaOrder)
    end
    if isnothing(WTypes)
        WTypes = collect(0:VerOrder)
    end
    return read_diagrams(filename; dim=dim, loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, WTypes=WTypes)
end

"""
    function PolarDiagrams(type::Symbol, MaxOrder::Int, has_counterterm::Bool=false, dim::Int=3)

    Generates a `Graph`: the `dim`-dimensional polarization diagrams in a given `type`, to a given maximum order `MaxOrder`, with switchable couterterms. 
    Generates fermionic/bosonic `LabelProduct`: `fermi_labelProd`/`bose_labelProd` for this `Graph`.

# Arguments:
- `type` (Symbol): The type of the diagrams, either `:spin` or `:charge`.
- `Maxorder` (Int): The maximum actual order of the diagrams.
- `has_counterterm` (Bool, optional): `false` for G0W0, `true` for GW with interaction and self-energy counterterms.
- `dim` (Int, optional): The dimension of the system (defaults to 3).

# Returns
A tuple `(diagrams, fermi_labelProd, bose_labelProd)` where 
- `diagrams` is a `Graph` objects representing the diagrams, 
- `fermi_labelProd` is a `LabelProduct` object containing the labels for the fermionic `G` objects in the diagrams, 
- `bose_labelProd` is a `LabelProduct` object containing the labels for the bosonic `W` objects in the diagrams.
"""
function PolarDiagrams(type::Symbol, MaxOrder::Int, has_counterterm::Bool=false, dim::Int=3)
    graphs = Graph{_dtype.factor,_dtype.weight}[]
    tau_labels = collect(1:MaxOrder+1)
    loopPool = LoopPool(:K, dim, MaxOrder + 1, Float64)
    if has_counterterm
        GTypes = collect(0:(MaxOrder-1)÷2)
        WTypes = collect(0:MaxOrder-2)
        for order in 1:MaxOrder
            for VerOrder in WTypes
                order == 1 && VerOrder > 0 && continue
                for SigmaOrder in GTypes
                    order + VerOrder + 2 * SigmaOrder > MaxOrder && continue
                    g, fermi_labelProd, bose_labelProd = PolarEachOrder(type, order, VerOrder, SigmaOrder;
                        dim=dim, loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, WTypes=WTypes)
                    push!(graphs, g)
                    loopPool = fermi_labelProd.labels[3]
                end
            end
        end
    else
        GTypes, WTypes = [0], [0]
        for order in 1:MaxOrder
            g, fermi_labelProd, bose_labelProd = PolarEachOrder(type, order;
                loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, WTypes=WTypes)
            push!(graphs, g)
            loopPool = fermi_labelProd.labels[3]
        end
    end
    fermi_labelProd = LabelProduct(tau_labels, GTypes, loopPool)
    bose_labelProd = LabelProduct(tau_labels, WTypes, loopPool)

    return IR.linear_combination(graphs, ones(_dtype.factor, length(graphs))), fermi_labelProd, bose_labelProd
end

end