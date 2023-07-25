module GV

import ..QuantumOperators as Op
import ..ComputationalGraphs as IR
import ..ComputationalGraphs: Graph
import ..ComputationalGraphs: _dtype
using ..FrontEnds
using AbstractTrees

export PolarEachOrder, PolarDiagrams, SigmaDiagrams, LeavesState

include("GV_diagrams/readfile.jl")

"""
    function PolarEachOrder(type::Symbol, order::Int, VerOrder::Int=0, GOrder::Int=0; loopPool::Union{LoopPool,Nothing}=nothing,
        tau_labels::Union{Nothing,Vector{Int}}=nothing, GTypes::Union{Nothing,Vector{Int}}=nothing, VTypes::Union{Nothing,Vector{Int}}=nothing)
 
    Generates a `Graph`: the polarization diagrams with static interactions of a given order, where the actual order of diagrams equals to `order + VerOrder + 2 * GOrder`.
    Generates fermionic/bosonic `LabelProduct`: `fermi_labelProd`/`bose_labelProd` with inputs `tau_labels`, `GTypes`/`VTypes`, and updated `loopPool`. 

# Arguments:
- `type` (Symbol): The type of the diagrams, either `:spin` or `:charge`.
- `order` (Int): The order of the diagrams without counterterms.
- `VerOrder` (Int, optional): The order of interaction counterterms (defaults to 0).
- `GOrder` (Int, optional): The order of self-energy counterterms (defaults to 0).
- `dim` (Int, optional): The dimension of the system (defaults to 3).
- `loopPool` (Union{LoopPool,Nothing}=nothing, optional): The initial pool of loop momenta (defaults to `nothing`).
- `tau_labels`(Union{Nothing, Vector{Int}}, optional): The labels for the discrete time of each vertex. (defaults to `nothing`).
- `GTypes`: The types of fermion propagators `G` in the diagrams (defaults to `collect(0:GOrder)`).
- `VTypes`: The types of boson static interaction `V` in the diagrams (defaults to `collect(0:VerOrder)`).

# Returns
A tuple `(diagrams, fermi_labelProd, bose_labelProd)` where 
- `diagrams` is a `Graph` objects representing the diagrams, 
- `fermi_labelProd` is a `LabelProduct` object containing the labels for the fermionic `G` objects in the diagrams, 
- `bose_labelProd` is a `LabelProduct` object containing the labels for the bosonic `W` objects in the diagrams.
"""
function PolarEachOrder(type::Symbol, order::Int, VerOrder::Int=0, GOrder::Int=0; dim::Int=3, loopPool::Union{LoopPool,Nothing}=nothing,
    tau_labels::Union{Nothing,Vector{Int}}=nothing, GTypes::Union{Nothing,Vector{Int}}=nothing, VTypes::Union{Nothing,Vector{Int}}=nothing)
    diagtype = :polar
    if type == :spin
        filename = string(@__DIR__, "/GV_diagrams/groups_spin/Polar$(order)_$(VerOrder)_$(GOrder).diag")
    elseif type == :charge
        filename = string(@__DIR__, "/GV_diagrams/groups_charge/Polar$(order)_$(VerOrder)_$(GOrder).diag")
    elseif type == :sigma
        diagtype = type
        filename = string(@__DIR__, "/GV_diagrams/groups_sigma/Sigma$(order)_$(VerOrder)_$(GOrder).diag")
    end

    # println("Reading ", filename)

    if isnothing(GTypes)
        GTypes = collect(0:GOrder)
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
- `type` (Symbol): The type of the polarization diagrams, either `:spin` or `:charge`.
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
    MaxLoopNum = MaxOrder + 1
    tau_labels = collect(1:MaxOrder+1)
    loopPool = LoopPool(:K, dim, MaxLoopNum, Float64)
    if has_counterterm
        GTypes = collect(0:MaxOrder-1)
        VTypes = collect(0:MaxOrder-1)
        for order in 1:MaxOrder
            for VerOrder in VTypes
                order == 1 && VerOrder > 0 && continue
                for GOrder in 0:MaxOrder-1
                    order + VerOrder + GOrder > MaxOrder && continue
                    g, fermi_labelProd, bose_labelProd = PolarEachOrder(type, order, VerOrder, GOrder;
                        dim=dim, loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, VTypes=VTypes)
                    g.name = "$(order)$(VerOrder)$(GOrder)"
                    push!(graphs, g)
                    loopPool = fermi_labelProd.labels[3]
                end
            end
        end
    else
        GTypes, VTypes = [0], [0]
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

function SigmaDiagrams(MaxOrder::Int, has_counterterm::Bool=false, dim::Int=3)
    # dict_graphs = Dict{Tuple{Int,Int,Int},Vector{Graph{_dtype.factor,_dtype.weight}}}()
    dict_graphs = Dict{Tuple{Int,Int,Int},Tuple{Vector{Graph{_dtype.factor,_dtype.weight}},Vector{Vector{Int}}}}()
    MaxLoopNum = MaxOrder + 2
    # tau_labels = collect(1:MaxOrder)
    tau_labels = collect(1:MaxOrder+2)
    loopPool = LoopPool(:K, dim, MaxLoopNum, Float64)
    # graphvector = Vector{Graph{_dtype.factor,_dtype.weight}}()

    if has_counterterm
        GTypes = collect(0:MaxOrder-1)
        append!(GTypes, [-2, -3])
        VTypes = collect(0:MaxOrder-1)
        for order in 1:MaxOrder
            for VerOrder in VTypes
                for GOrder in 0:MaxOrder-1
                    order + VerOrder + GOrder > MaxOrder && continue
                    gvec, fermi_labelProd, bose_labelProd, extT_labels = PolarEachOrder(:sigma, order, VerOrder, GOrder;
                        dim=dim, loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, VTypes=VTypes)
                    # g.name = "$(order)$(VerOrder)$(GOrder)"
                    key = (order, VerOrder, GOrder)
                    dict_graphs[key] = (gvec, extT_labels)
                    # push!(graphs, gvec)
                    loopPool = fermi_labelProd.labels[3]
                end
            end
        end
    else
        GTypes, VTypes = [0], [0]
        append!(GTypes, [-2, -3])
        for order in 1:MaxOrder
            gvec, fermi_labelProd, bose_labelProd, extT_labels = PolarEachOrder(:sigma, order;
                loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, VTypes=VTypes)
            key = (order, 0, 0)
            dict_graphs[key] = (gvec, extT_labels)
            loopPool = fermi_labelProd.labels[3]
        end
    end
    # IR.optimize!(graphvector)
    fermi_labelProd = LabelProduct(tau_labels, GTypes, loopPool)
    bose_labelProd = LabelProduct(tau_labels, VTypes, loopPool)

    # return IR.linear_combination(graphs_eqT, ones(_dtype.factor, length(graphs_eqT))), fermi_labelProd, bose_labelProd
    return dict_graphs, fermi_labelProd, bose_labelProd
end

function SigmaDiagrams(gkeys::Vector{Tuple{Int,Int,Int}}, dim::Int=3)
    dict_graphs = Dict{Tuple{Int,Int,Int},Tuple{Vector{Graph{_dtype.factor,_dtype.weight}},Vector{Vector{Int}}}}()
    MaxLoopNum = maximum([key[1] for key in gkeys]) + 2
    MaxVerOrder = maximum([key[2] for key in gkeys])
    MaxGOrder = maximum([key[3] for key in gkeys])

    tau_labels = collect(1:MaxLoopNum)
    loopPool = LoopPool(:K, dim, MaxLoopNum, Float64)
    GTypes = collect(0:MaxGOrder)
    append!(GTypes, [-2, -3])
    VTypes = collect(0:MaxVerOrder)

    # graphvector = Vector{_dtype.factor,_dtype.weight}()
    propagatorMap, interactionMap = Dict{eltype(gkeys),Dict{Int,Int}}(), Dict{eltype(gkeys),Dict{Int,Int}}()
    for key in gkeys
        gvec, fermi_labelProd, bose_labelProd, extT_labels = PolarEachOrder(:sigma, key...;
            dim=dim, loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, VTypes=VTypes)
        dict_graphs[key] = (gvec, extT_labels)
        loopPool = fermi_labelProd.labels[3]
        propagatorMap[key], interactionMap[key] = IR.optimize!(gvec)
        # append!(graphvector, gvec)
    end
    # IR.optimize!(graphvector)

    fermi_labelProd = LabelProduct(tau_labels, GTypes, loopPool)
    bose_labelProd = LabelProduct(tau_labels, VTypes, loopPool)
    return dict_graphs, fermi_labelProd, bose_labelProd, (propagatorMap, interactionMap)
end

function LeavesState(FeynGraphs::Dict{T,Tuple{Vector{G},Vector{Vector{Int}}}},
    FermiLabel::LabelProduct, BoseLabel::LabelProduct, graph_keys::Vector{T}) where {T,G<:Graph}
    #read information of each leaf from the generated graph and its LabelProduct, the information include type, loop momentum, imaginary time.
    num_g = length(graph_keys)
    ExtT_index = [Vector{Vector{Int}}() for _ in 1:num_g]

    PropagatorType = [Vector{Int}() for _ in 1:num_g]
    PropagatorInTau = [Vector{Int}() for _ in 1:num_g]
    PropagatorOutTau = [Vector{Int}() for _ in 1:num_g]
    PropagatorLoopIndex = [Vector{Int}() for _ in 1:num_g]
    PropagatorValue = [Vector{Float64}() for _ in 1:num_g]
    InteractionType = [Vector{Int}() for _ in 1:num_g]
    InteractionInTau = [Vector{Int}() for _ in 1:num_g]
    InteractionOutTau = [Vector{Int}() for _ in 1:num_g]
    InteractionLoopIndex = [Vector{Int}() for _ in 1:num_g]
    InteractionValue = [Vector{Float64}() for _ in 1:num_g]

    for (ikey, key) in enumerate(graph_keys)
        ExtT_index[ikey] = FeynGraphs[key][2]  # external tau variables
        for j in eachindex(ExtT_index[ikey])
            for g in Leaves(FeynGraphs[key][1][j])
                g.name == "visited" && continue
                if g.type == IR.Interaction
                    push!(InteractionType[ikey], 0)
                    In = Out = g.vertices[1][1].label
                    push!(InteractionLoopIndex[ikey], 1)
                    push!(InteractionInTau[ikey], FermiLabel[In][1])
                    push!(InteractionOutTau[ikey], FermiLabel[Out][1])
                    push!(InteractionValue[ikey], 1.0)
                elseif g.type == IR.Propagator
                    if (Op.isfermionic(g.vertices[1]))
                        In, Out = g.vertices[2][1].label, g.vertices[1][1].label
                        if FermiLabel[In][2] in [-2, -3]
                            push!(PropagatorType[ikey], 0)
                            push!(PropagatorLoopIndex[ikey], 1)
                        else
                            push!(PropagatorType[ikey], FermiLabel[In][2] * 2 + 1)
                            push!(PropagatorLoopIndex[ikey], FrontEnds.linear_to_index(FermiLabel, In)[end]) #the label of LoopPool for each fermionic leaf
                        end
                        push!(PropagatorInTau[ikey], FermiLabel[In][1])
                        push!(PropagatorOutTau[ikey], FermiLabel[Out][1])
                    else
                        In, Out = g.vertices[2][1].label, g.vertices[1][1].label
                        push!(PropagatorType[ikey], BoseLabel[In][2] * 2 + 2)
                        push!(PropagatorLoopIndex[ikey], FrontEnds.linear_to_index(BoseLabel, In)[end]) #the label of LoopPool for each bosonic leaf
                        push!(PropagatorInTau[ikey], BoseLabel[In][1])
                        push!(PropagatorOutTau[ikey], BoseLabel[Out][1])
                    end
                    push!(PropagatorValue[ikey], 1.0)
                end
                g.name = "visited"
            end
        end
    end
    return (PropagatorValue, PropagatorType, PropagatorInTau, PropagatorOutTau, PropagatorLoopIndex),
    (InteractionValue, InteractionType, InteractionInTau, InteractionOutTau, InteractionLoopIndex), ExtT_index
end

# function LeavesState(FeynGraphs::Dict{T,Tuple{Vector{G},Vector{Vector{Int}}}},
#     FermiLabel::LabelProduct, BoseLabel::LabelProduct, graph_keys::Vector{T}) where {T,G<:Graph}
#     #read information of each leaf from the generated graph and its LabelProduct, the information include type, loop momentum, imaginary time.
#     ExtT_index = [Vector{Vector{Int}}() for _ in eachindex(graph_keys)]
#     dict_pools = Dict{T,CachePool{Tuple{Int,Int,Int,Int},Float64}}()

#     for (ikey, key) in enumerate(graph_keys)
#         ExtT_index[ikey] = FeynGraphs[key][2]  # external tau variables
#         pool = CachePool(Tuple{Int,Int,Int,Int}, Float64)
#         for j in eachindex(ExtT_index[ikey])
#             for g in Leaves(FeynGraphs[key][1][j])
#                 if g.type == IR.Interaction
#                     LeafType, LeafLoopIndex = 0, 1
#                     In = Out = g.vertices[1][1].label
#                 elseif (Op.isfermionic(g.vertices[1]))
#                     In, Out = g.vertices[2][1].label, g.vertices[1][1].label
#                     if FermiLabel[In][2] in [-2, -3]
#                         LeafType, LeafLoopIndex = 0, 1
#                     else
#                         LeafType = FermiLabel[In][2] * 2 + 1  # 1,3,5...
#                         LeafLoopIndex = FrontEnds.linear_to_index(FermiLabel, In)[end] #the label of LoopPool for each fermionic leaf
#                     end
#                 else
#                     In, Out = g.vertices[2][1].label, g.vertices[1][1].label
#                     LeafType = BoseLabel[In][2] * 2 + 2  # 2,4,6...
#                     LeafLoopIndex = FrontEnds.linear_to_index(BoseLabel, In)[end] #the label of LoopPool for each bosonic leaf
#                 end
#                 id = FrontEnds.FastHashInt(g.id)
#                 if haskey(pool.object, id)
#                     @assert pool.object[id] == (LeafType, LeafLoopIndex, FermiLabel[In][1], FermiLabel[Out][1]) "the leaf with the same id should be the same!"
#                 else
#                     pool.object[id] = (LeafType, LeafLoopIndex, FermiLabel[In][1], FermiLabel[Out][1])
#                     pool[g.id] = 1.0
#                 end
#             end
#         end
#         dict_pools[key] = pool
#     end
#     return dict_pools, ExtT_index
# end

end