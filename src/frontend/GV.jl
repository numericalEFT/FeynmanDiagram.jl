module GV

import ..QuantumOperators as Op
import ..ComputationalGraphs as IR
import ..ComputationalGraphs: Graph
import ..ComputationalGraphs: _dtype
import ..ComputationalGraphs: group
using ..FrontEnds

export PolarEachOrder, PolarDiagrams, SigmaDiagrams

include("GV_diagrams/readfile.jl")

"""
    function PolarEachOrder(type::Symbol, order::Int, VerOrder::Int=0, SigmaOrder::Int=0; loopPool::Union{LoopPool,Nothing}=nothing,
        tau_labels::Union{Nothing,Vector{Int}}=nothing, GTypes::Union{Nothing,Vector{Int}}=nothing, VTypes::Union{Nothing,Vector{Int}}=nothing)
 
    Generates a `Graph`: the polarization diagrams with static interactions of a given order, where the actual order of diagrams equals to `order + VerOrder + 2 * SigmaOrder`.
    Generates fermionic/bosonic `LabelProduct`: `fermi_labelProd`/`bose_labelProd` with inputs `tau_labels`, `GTypes`/`VTypes`, and updated `loopPool`. 

# Arguments:
- `type` (Symbol): The type of the diagrams, either `:spin` or `:charge`.
- `order` (Int): The order of the diagrams without counterterms.
- `VerOrder` (Int, optional): The order of interaction counterterms (defaults to 0).
- `SigmaOrder` (Int, optional): The order of self-energy counterterms (defaults to 0).
- `dim` (Int, optional): The dimension of the system (defaults to 3).
- `loopPool` (Union{LoopPool,Nothing}=nothing, optional): The initial pool of loop momenta (defaults to `nothing`).
- `tau_labels`(Union{Nothing, Vector{Int}}, optional): The labels for the discrete time of each vertex. (defaults to `nothing`).
- `GTypes`: The types of fermion propagators `G` in the diagrams (defaults to `collect(0:SigmaOrder)`).
- `VTypes`: The types of boson static interaction `V` in the diagrams (defaults to `collect(0:VerOrder)`).

# Returns
A tuple `(diagrams, fermi_labelProd, bose_labelProd)` where 
- `diagrams` is a `Graph` objects representing the diagrams, 
- `fermi_labelProd` is a `LabelProduct` object containing the labels for the fermionic `G` objects in the diagrams, 
- `bose_labelProd` is a `LabelProduct` object containing the labels for the bosonic `W` objects in the diagrams.
"""
function PolarEachOrder(type::Symbol, order::Int, VerOrder::Int=0, SigmaOrder::Int=0; dim::Int=3, loopPool::Union{LoopPool,Nothing}=nothing,
    tau_labels::Union{Nothing,Vector{Int}}=nothing, GTypes::Union{Nothing,Vector{Int}}=nothing, VTypes::Union{Nothing,Vector{Int}}=nothing)
    diagtype = :polar
    if type == :spin
        filename = string(@__DIR__, "/GV_diagrams/groups_spin/Polar$(order)_$(VerOrder)_$(SigmaOrder).diag")
    elseif type == :charge
        filename = string(@__DIR__, "/GV_diagrams/groups_charge/Polar$(order)_$(VerOrder)_$(SigmaOrder).diag")
    elseif type == :sigma
        diagtype = type
        filename = string(@__DIR__, "/GV_diagrams/groups_sigma/Sigma$(order)_$(VerOrder)_$(SigmaOrder).diag")
    end

    # println("Reading ", filename)

    if isnothing(GTypes)
        GTypes = collect(0:SigmaOrder)
        type == :sigma && append!(GTypes, [-2, -3])
    end
    isnothing(VTypes) && (VTypes = collect(0:VerOrder))
    return read_diagrams(filename; dim=dim, loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, VTypes=VTypes, diagType=diagtype)
end

"""
    function PolarDiagrams(type::Symbol, MaxOrder::Int, has_counterterm::Bool=false, dim::Int=3)

    Generates a `Graph`: the `dim`-dimensional polarization diagrams with static interactions in a given `type`, to a given maximum order `MaxOrder`, with switchable couterterms. 
    Generates fermionic/bosonic `LabelProduct`: `fermi_labelProd`/`bose_labelProd` for this `Graph`.

# Arguments:
- `type` (Symbol): The type of the diagrams, either `:spin`, `:charge`, or `:sigma`.
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
    if type == :sigma
        MaxLoopNum = MaxOrder + 2
        tau_labels = collect(1:MaxOrder)
    else
        tau_labels = collect(1:MaxOrder+1)
    end
    # MaxLoopNum = MaxOrder + 1
    loopPool = LoopPool(:K, dim, MaxLoopNum, Float64)
    if has_counterterm
        GTypes = collect(0:MaxOrder-1)
        type == :sigma && append!(GTypes, [-2, -3])
        VTypes = collect(0:MaxOrder-1)
        for order in 1:MaxOrder
            # order != 1 && continue
            for VerOrder in VTypes
                type != :sigma && order == 1 && VerOrder > 0 && continue
                for SigmaOrder in 0:MaxOrder-1
                    order + VerOrder + SigmaOrder > MaxOrder && continue
                    g, fermi_labelProd, bose_labelProd = PolarEachOrder(type, order, VerOrder, SigmaOrder;
                        dim=dim, loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, VTypes=VTypes)
                    g.name = "$(order)$(VerOrder)$(SigmaOrder)"
                    push!(graphs, g)
                    loopPool = fermi_labelProd.labels[3]
                end
            end
        end
    else
        GTypes, VTypes = [0], [0]
        type == :sigma && append!(GTypes, [-2, -3])
        for order in 1:MaxOrder
            g, fermi_labelProd, bose_labelProd = PolarEachOrder(type, order;
                loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, VTypes=VTypes)
            push!(graphs, g)
            loopPool = fermi_labelProd.labels[3]
        end
    end
    fermi_labelProd = LabelProduct(tau_labels, GTypes, loopPool)
    bose_labelProd = LabelProduct(tau_labels, VTypes, loopPool)

    return IR.linear_combination(graphs, ones(_dtype.factor, length(graphs))), fermi_labelProd, bose_labelProd
end

function SigmaDiagrams(type::Symbol, MaxOrder::Int, has_counterterm::Bool=false, dim::Int=3)
    dict_graphs = Dict{Tuple{Int,Int,Int},Vector{Graph{_dtype.factor,_dtype.weight}}}()
    MaxLoopNum = MaxOrder + 2
    tau_labels = collect(1:MaxOrder)
    # MaxLoopNum = MaxOrder + 1
    loopPool = LoopPool(:K, dim, MaxLoopNum, Float64)
    if has_counterterm
        GTypes = collect(0:MaxOrder-1)
        append!(GTypes, [-2, -3])
        VTypes = collect(0:MaxOrder-1)
        for order in 1:MaxOrder
            for VerOrder in VTypes
                for SigmaOrder in 0:MaxOrder-1
                    order + VerOrder + SigmaOrder > MaxOrder && continue
                    gvec, fermi_labelProd, bose_labelProd = PolarEachOrder(:sigma, order, VerOrder, SigmaOrder;
                        dim=dim, loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, VTypes=VTypes)
                    # g.name = "$(order)$(VerOrder)$(SigmaOrder)"
                    key = (order, VerOrder, SigmaOrder)
                    dict_graphs[key] = gvec
                    # push!(graphs, gvec)
                    loopPool = fermi_labelProd.labels[3]
                end
            end
        end
    else
        GTypes, VTypes = [0], [0]
        append!(GTypes, [-2, -3])
        for order in 1:MaxOrder
            gvec, fermi_labelProd, bose_labelProd = PolarEachOrder(:sigma, order;
                loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, VTypes=VTypes)
            key = (order, VerOrder, SigmaOrder)
            dict_graphs[key] = gvec
            # push!(graphs, gvec)
            loopPool = fermi_labelProd.labels[3]
        end
    end
    fermi_labelProd = LabelProduct(tau_labels, GTypes, loopPool)
    bose_labelProd = LabelProduct(tau_labels, VTypes, loopPool)

    # return IR.linear_combination(graphs_eqT, ones(_dtype.factor, length(graphs_eqT))), fermi_labelProd, bose_labelProd
    return dict_graphs, fermi_labelProd, bose_labelProd
end

end