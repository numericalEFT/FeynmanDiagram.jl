module FrontEnds

import ..ComputationalGraphs as IR
import ..ComputationalGraphs: Graph, FeynmanGraph, _dtype
import ..QuantumOperators
import ..Taylor
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

include("GV.jl")
export GV

include("parquet/parquet.jl")
export Parquet

# include("strong_coupling_expansion_builder/strong_coupling_expansion.jl")
# export SCE

"""
    function leafstates(leaf_maps::Vector{Dict{Int,G}}, labelProd::LabelProduct) where {G<:Union{Graph,FeynmanGraph}}

    Extracts leaf information from the leaf mapping from the leaf value's index to the leaf node for all graph partitions
    and their associated LabelProduct data (`labelProd`). 
    The information includes their initial value, type, orders, in/out time, and loop momenta.
    
# Arguments:
- `leaf_maps`: A vector of the dictionary mapping the leaf value's index to the FeynmanGraph/Graph of this leaf. 
               Each dict corresponds to a graph partition, such as (order, Gorder, Vorder).
- `labelProd`: A LabelProduct used to label the leaves of graphs.

# Returns
- A tuple of vectors containing information about the leaves of graphs, including their initial values, types, orders, input and output time indexes, and loop-momenta indexes.
"""
function leafstates(leaf_maps::Vector{Dict{Int,G}}, labelProd::LabelProduct) where {G<:Union{Graph,FeynmanGraph}}
    #read information of each leaf from the generated graph and its LabelProduct, the information include type, loop momentum, imaginary time.
    num_g = length(leaf_maps)
    leafType = [Vector{Int}() for _ in 1:num_g]
    leafOrders = [Vector{Vector{Int}}() for _ in 1:num_g]
    leafInTau = [Vector{Int}() for _ in 1:num_g]
    leafOutTau = [Vector{Int}() for _ in 1:num_g]
    leafLoopIndex = [Vector{Int}() for _ in 1:num_g]
    leafValue = [Vector{Float64}() for _ in 1:num_g]

    for (ikey, leafmap) in enumerate(leaf_maps)
        len_leaves = length(keys(leafmap))
        sizehint!(leafType[ikey], len_leaves)
        sizehint!(leafOrders[ikey], len_leaves)
        sizehint!(leafInTau[ikey], len_leaves)
        sizehint!(leafOutTau[ikey], len_leaves)
        sizehint!(leafLoopIndex[ikey], len_leaves)
        leafValue[ikey] = ones(Float64, len_leaves)

        for idx in 1:len_leaves
            g = leafmap[idx]
            vertices = g.properties.vertices
            if IR.diagram_type(g) == IR.Interaction
                In = Out = vertices[1][1].label
                push!(leafType[ikey], 0)
                push!(leafLoopIndex[ikey], 1)
            elseif IR.diagram_type(g) == IR.Propagator
                if (Op.isfermionic(vertices[1]))
                    In, Out = vertices[2][1].label, vertices[1][1].label
                    # push!(leafType[ikey], g.orders[1] * 2 + 1)
                    push!(leafType[ikey], 1)
                    push!(leafLoopIndex[ikey], FrontEnds.linear_to_index(labelProd, In)[end]) #the label of LoopPool for each fermionic leaf
                else
                    In, Out = vertices[2][1].label, vertices[1][1].label
                    # push!(leafType[ikey], g.orders[2] * 2 + 2)
                    push!(leafType[ikey], 2)
                    push!(leafLoopIndex[ikey], FrontEnds.linear_to_index(labelProd, In)[end]) #the label of LoopPool for each bosonic leaf
                end
            end
            push!(leafOrders[ikey], g.orders)
            push!(leafInTau[ikey], labelProd[In][1])
            push!(leafOutTau[ikey], labelProd[Out][1])
        end
    end
    return (leafValue, leafType, leafOrders, leafInTau, leafOutTau, leafLoopIndex)
end

"""
    function leafstates(leaf_maps::Vector{Dict{Int,G}}, maxloopNum::Int)

    Extracts leaf information from the leaf mapping from the leaf value's index to the leaf node for all graph partitions. 
    The information includes their initial value, type, orders, in/out time, and loop momentum index.
    The loop basis is also obtained for all the graphs.
    
# Arguments:
- `leaf_maps`: A vector of the dictionary mapping the leaf value's index to the Graph of this leaf. 
               Each dict corresponds to a graph partition, such as (order, Gorder, Vorder).
- `maxloopNum`: The maximum loop-momentum number.

# Returns
- A tuple of vectors containing information about the leaves of graphs, including their initial values, types, orders, input and output time indexes, and loop-momenta indexes.
- Loop-momentum basis (`::Vector{Vector{Float64}}`) for all the graphs.
"""
function leafstates(leaf_maps::Vector{Dict{Int,G}}, maxloopNum::Int) where {G<:Graph}

    num_g = length(leaf_maps)
    leafType = [Vector{Int}() for _ in 1:num_g]
    leafOrders = [Vector{Vector{Int}}() for _ in 1:num_g]
    leafInTau = [Vector{Int}() for _ in 1:num_g]
    leafOutTau = [Vector{Int}() for _ in 1:num_g]
    leafLoopIndex = [Vector{Int}() for _ in 1:num_g]
    leafValue = [Vector{Float64}() for _ in 1:num_g]

    loopbasis = Vector{Float64}[]
    for (ikey, leafmap) in enumerate(leaf_maps)
        len_leaves = length(keys(leafmap))
        sizehint!(leafType[ikey], len_leaves)
        sizehint!(leafOrders[ikey], len_leaves)
        sizehint!(leafInTau[ikey], len_leaves)
        sizehint!(leafOutTau[ikey], len_leaves)
        sizehint!(leafLoopIndex[ikey], len_leaves)
        leafValue[ikey] = ones(Float64, len_leaves)

        for idx in 1:len_leaves
            leaf = leafmap[idx]
            @assert IR.isleaf(leaf)
            diagId, leaf_orders = leaf.properties, leaf.orders
            loopmom = copy(diagId.extK)
            len = length(loopmom)
            @assert maxloopNum >= len

            if maxloopNum > length(loopmom)
                Base.append!(loopmom, zeros(Float64, maxloopNum - len))
            end
            flag = true
            for bi in eachindex(loopbasis)
                if loopbasis[bi] ≈ loopmom
                    push!(leafLoopIndex[ikey], bi)
                    flag = false
                    break
                end
            end
            if flag
                push!(loopbasis, loopmom)
                push!(leafLoopIndex[ikey], length(loopbasis))
            end

            push!(leafInTau[ikey], diagId.extT[1])
            push!(leafOutTau[ikey], diagId.extT[2])

            push!(leafOrders[ikey], leaf_orders)
            push!(leafType[ikey], FrontEnds.index(typeof(diagId)))

        end
    end

    return (leafValue, leafType, leafOrders, leafInTau, leafOutTau, leafLoopIndex), loopbasis
end

export leafstates

end