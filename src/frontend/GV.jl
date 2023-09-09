module GV

import ..QuantumOperators as Op
import ..ComputationalGraphs as IR
import ..ComputationalGraphs: Graph
import ..ComputationalGraphs: _dtype
using ..FrontEnds
using AbstractTrees

export eachorder_diag, diagdictGV, leafstates

include("GV_diagrams/readfile.jl")

"""
    function eachorder_diag(type::Symbol, order::Int, VerOrder::Int=0, GOrder::Int=0; loopPool::Union{LoopPool,Nothing}=nothing,
        tau_labels::Union{Nothing,Vector{Int}}=nothing, GTypes::Union{Nothing,Vector{Int}}=nothing, VTypes::Union{Nothing,Vector{Int}}=nothing)
 
    Generates a `Vector{Graph}`: the polarization diagrams with static interactions of a given order, where the actual order of diagrams equals to `order + VerOrder + 2 * GOrder`.
    Generates fermionic/bosonic `LabelProduct`: `fermi_labelProd`/`bose_labelProd` with inputs `tau_labels`, `GTypes`/`VTypes`, and updated `loopPool`. 
    Generates external tau labels Vector{Vector{Int}}. The i-th labels (Vector{Int}) corresponds to the i-th `Graph` in `Vector{Graph}`.

# Arguments:
- `type` (Symbol): The type of the diagrams, including `:spinPolar`, `:chargePolar`, `:sigma`, `:green`, or `:freeEnergy`.
- `order` (Int): The order of the diagrams without counterterms.
- `VerOrder` (Int, optional): The order of interaction counterterms (defaults to 0).
- `GOrder` (Int, optional): The order of self-energy counterterms (defaults to 0).
- `dim` (Int, optional): The dimension of the system (defaults to 3).
- `spinPolarPara` (Float64, optional): The spin-polarization parameter (n_up - n_down) / (n_up + n_down) (defaults to `0.0`).
- `loopPool` (Union{LoopPool,Nothing}=nothing, optional): The initial pool of loop momenta (defaults to `nothing`).
- `tau_labels`(Union{Nothing, Vector{Int}}, optional): The labels for the discrete time of each vertex. (defaults to `nothing`).
- `GTypes`: The types of fermion propagators `G` in the diagrams (defaults to `collect(0:GOrder)`).
- `VTypes`: The types of boson static interaction `V` in the diagrams (defaults to `collect(0:VerOrder)`).

# Returns
A tuple `(diagrams, fermi_labelProd, bose_labelProd, extT_labels)` where 
- `diagrams` is a `Vector{Graph}` object representing the diagrams, 
- `fermi_labelProd` is a `LabelProduct` object containing the labels for the fermionic `G` objects in the diagrams, 
- `bose_labelProd` is a `LabelProduct` object containing the labels for the bosonic `W` objects in the diagrams.
- `extT_labels` is a `Vector{Vector{Int}}` object containing the external tau labels for each `Graph` in `diagrams`.
"""
function eachorder_diag(type::Symbol, order::Int, GOrder::Int=0, VerOrder::Int=0; dim::Int=3, spinPolarPara::Float64=0.0,
    loopPool::Union{LoopPool,Nothing}=nothing, tau_labels::Union{Nothing,Vector{Int}}=nothing,
    GTypes::Union{Nothing,Vector{Int}}=nothing, VTypes::Union{Nothing,Vector{Int}}=nothing)
    diagtype = :polar
    if type == :spinPolar
        filename = string(@__DIR__, "/GV_diagrams/groups_spin/Polar$(order)_$(VerOrder)_$(GOrder).diag")
    elseif type == :chargePolar
        filename = string(@__DIR__, "/GV_diagrams/groups_charge/Polar$(order)_$(VerOrder)_$(GOrder).diag")
    elseif type == :sigma
        diagtype = type
        filename = string(@__DIR__, "/GV_diagrams/groups_sigma/Sigma$(order)_$(VerOrder)_$(GOrder).diag")
    elseif type == :green
        diagtype = type
        filename = string(@__DIR__, "/GV_diagrams/groups_green/Green$(order)_$(VerOrder)_$(GOrder).diag")
    elseif type == :freeEnergy
        diagtype = type
        filename = string(@__DIR__, "/GV_diagrams/groups_free_energy/FreeEnergy$(order)_$(VerOrder)_$(GOrder).diag")
    end

    # println("Reading ", filename)

    if isnothing(GTypes)
        GTypes = collect(0:GOrder)
        type == :sigma && append!(GTypes, [-2, -3])
        type == :green && push!(GTypes, -2)
    end
    isnothing(VTypes) && (VTypes = collect(0:VerOrder))
    return read_diagrams(filename; dim=dim, loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, VTypes=VTypes,
        diagType=diagtype, spinPolarPara=spinPolarPara)
end

"""
    function diagdictGV(type::Symbol, MaxOrder::Int, has_counterterm::Bool=false, dim::Int=3)

    Generates a Graph Dict: the `dim`-dimensional spin/charge polarization or self-energy diagrams with static interactions in a given `type`, to a given maximum order `MaxOrder`, with switchable couterterms. 
    Generates fermionic/bosonic `LabelProduct`: `fermi_labelProd`/`bose_labelProd` for these Graphs.
    Generates a Tuple (propagatorMap, interactionMap) for mapping `g.id` to the index of unique proapgators and interactions, respectively. 

# Arguments:
- `type` (Symbol): The type of the Feynman diagrams, including `:spinPolar`, `:chargePolar`, `:sigma`, `:green`, or `:freeEnergy`.
- `Maxorder` (Int): The maximum actual order of the diagrams.
- `has_counterterm` (Bool): `false` for G0W0, `true` for GW with self-energy and interaction counterterms (defaults to `false`).
- `dim` (Int): The dimension of the system (defaults to 3).
- `spinPolarPara` (Float64, optional): The spin-polarization parameter (n_up - n_down) / (n_up + n_down) (defaults to `0.0`).

# Returns
A tuple `(diagrams, fermi_labelProd, bose_labelProd)` where 
- `diagrams` is a `Dict{Tuple{Int,Int,Int},Tuple{Vector{Graph},Vector{Vector{Int}}}}` object representing the diagrams. 
   The key is (order, Gorder, Vorder). The element is a Tuple (diagrams, extT_labels).
- `fermi_labelProd` is a `LabelProduct` object containing the labels for the fermionic `G` objects in the diagrams, 
- `bose_labelProd` is a `LabelProduct` object containing the labels for the bosonic `W` objects in the diagrams.
- `(propagatorMap, interactionMap)` maps `g.id` to the index of unique proapgators and interactions, respectively. 
"""
function diagdictGV(type::Symbol, MaxOrder::Int, has_counterterm::Bool=false, dim::Int=3; spinPolarPara::Float64=0.0)
    dict_graphs = Dict{Tuple{Int,Int,Int},Tuple{Vector{Graph{_dtype.factor,_dtype.weight}},Vector{Vector{Int}}}}()
    if type == :sigma
        MaxLoopNum = MaxOrder + 2
        tau_labels = collect(1:MaxLoopNum)
    elseif type in [:chargePolar, :spinPolar, :green]
        MaxLoopNum = MaxOrder + 1
        tau_labels = collect(1:MaxLoopNum)
    elseif type == :freeEnergy
        MaxLoopNum = MaxOrder + 1
        tau_labels = collect(1:MaxLoopNum-1)
    else
        error("no support for $type diagram")
    end
    loopPool = LoopPool(:K, dim, MaxLoopNum, Float64)

    propagatorMap, interactionMap = Dict{Tuple{Int,Int,Int},Dict{Int,Int}}(), Dict{Tuple{Int,Int,Int},Dict{Int,Int}}()
    if has_counterterm
        GTypes = collect(0:MaxOrder-1)
        type == :sigma && append!(GTypes, [-2, -3])
        type == :green && push!(GTypes, -2)
        VTypes = collect(0:MaxOrder-1)
        for order in 1:MaxOrder
            for VerOrder in VTypes
                type in [:chargePolar, :spinPolar] && order == 1 && VerOrder > 0 && continue
                for GOrder in 0:MaxOrder-1
                    order + VerOrder + GOrder > MaxOrder && continue
                    gvec, fermi_labelProd, bose_labelProd, extT_labels = eachorder_diag(type, order, GOrder, VerOrder;
                        dim=dim, loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, VTypes=VTypes, spinPolarPara=spinPolarPara)
                    # push!(graphs, g)
                    key = (order, GOrder, VerOrder)
                    dict_graphs[key] = (gvec, extT_labels)
                    loopPool = fermi_labelProd.labels[3]
                    propagatorMap[key], interactionMap[key] = IR.optimize!(gvec)
                end
            end
        end
    else
        GTypes, VTypes = [0], [0]
        type == :sigma && append!(GTypes, [-2, -3])
        type == :green && push!(GTypes, -2)
        for order in 1:MaxOrder
            gvec, fermi_labelProd, bose_labelProd, extT_labels = eachorder_diag(type, order;
                loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, VTypes=VTypes, spinPolarPara=spinPolarPara)
            # push!(graphs, g)
            key = (order, 0, 0)
            dict_graphs[key] = (gvec, extT_labels)
            loopPool = fermi_labelProd.labels[3]
            propagatorMap[key], interactionMap[key] = IR.optimize!(gvec)
        end
    end
    fermi_labelProd = LabelProduct(tau_labels, GTypes, loopPool)
    bose_labelProd = LabelProduct(tau_labels, VTypes, loopPool)

    # return IR.linear_combination(graphs, ones(_dtype.factor, length(graphs))), fermi_labelProd, bose_labelProd
    return dict_graphs, fermi_labelProd, bose_labelProd, (propagatorMap, interactionMap)
end

"""
    function diagdictGV(type::Symbol, MaxOrder::Int, has_counterterm::Bool=false, dim::Int=3)

    Generates a Graph Dict: the `dim`-dimensional spin/charge polarization or self-energy diagrams with static interactions in a given `type`, to a given maximum order `MaxOrder`, with switchable couterterms. 
    Generates fermionic/bosonic `LabelProduct`: `fermi_labelProd`/`bose_labelProd` for these Graphs.
    Generates a Tuple (propagatorMap, interactionMap) for mapping `g.id` to the index of unique proapgators and interactions, respectively. 

# Arguments:
- `type` (Symbol): The type of the Feynman diagrams, including `:spinPolar`, `:chargePolar`, `:sigma`, `:green`, or `:freeEnergy`.
- `gkeys` (Vector{Tuple{Int,Int,Int}}): The (order, Gorder, Vorder) of the diagrams. Gorder is the order of self-energy counterterms, and Vorder is the order of interaction counterterms. 
- `dim` (Int): The dimension of the system (defaults to 3).
- `spinPolarPara` (Float64, optional): The spin-polarization parameter (n_up - n_down) / (n_up + n_down) (defaults to `0.0`).

# Returns
A tuple `(diagrams, fermi_labelProd, bose_labelProd)` where 
- `diagrams` is a `Dict{Tuple{Int,Int,Int},Tuple{Vector{Graph},Vector{Vector{Int}}}}` object representing the diagrams. 
   The key is (order, Gorder, Vorder). The element is a Tuple (diagrams, extT_labels).
- `fermi_labelProd` is a `LabelProduct` object containing the labels for the fermionic `G` objects in the diagrams, 
- `bose_labelProd` is a `LabelProduct` object containing the labels for the bosonic `W` objects in the diagrams.
- `(propagatorMap, interactionMap)` maps `g.id` to the index of unique proapgators and interactions, respectively. 
"""
function diagdictGV(type::Symbol, gkeys::Vector{Tuple{Int,Int,Int}}, dim::Int=3; spinPolarPara::Float64=0.0)
    dict_graphs = Dict{Tuple{Int,Int,Int},Tuple{Vector{Graph{_dtype.factor,_dtype.weight}},Vector{Vector{Int}}}}()
    if type == :sigma
        MaxLoopNum = maximum([key[1] for key in gkeys]) + 2
        tau_labels = collect(1:MaxLoopNum)
    elseif type in [:chargePolar, :spinPolar, :green]
        MaxLoopNum = maximum([key[1] for key in gkeys]) + 1
        tau_labels = collect(1:MaxLoopNum)
    elseif type == :freeEnergy
        MaxLoopNum = maximum([key[1] for key in gkeys]) + 1
        tau_labels = collect(1:MaxLoopNum-1)
    else
        error("no support for $type diagram")
    end
    MaxGOrder = maximum([key[2] for key in gkeys])
    MaxVerOrder = maximum([key[3] for key in gkeys])

    loopPool = LoopPool(:K, dim, MaxLoopNum, Float64)
    GTypes = collect(0:MaxGOrder)
    type == :sigma && append!(GTypes, [-2, -3])
    type == :green && push!(GTypes, -2)
    VTypes = collect(0:MaxVerOrder)

    # graphvector = Vector{_dtype.factor,_dtype.weight}()
    propagatorMap, interactionMap = Dict{eltype(gkeys),Dict{Int,Int}}(), Dict{eltype(gkeys),Dict{Int,Int}}()
    for key in gkeys
        gvec, fermi_labelProd, bose_labelProd, extT_labels = eachorder_diag(type, key...;
            dim=dim, loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, VTypes=VTypes, spinPolarPara=spinPolarPara)
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

"""
    function leafstates(
        FeynGraphs::Dict{T, Tuple{Vector{G}, Vector{Vector{Int}}}},
        FermiLabel::LabelProduct, BoseLabel::LabelProduct,
        graph_keys::Vector{T}
    ) where {T, G <: Graph}

    Extracts leaf information from a Dict collection of Feynman graphs (`FeynGraphs` with its keys `graph_keys`)
    and their associated LabelProduct data (`FermiLabel` and `BoseLabel`). 
    The information includes their initial value, type, in/out time, and loop momenta.
    
# Arguments:
- `FeynGraphs`: A dictionary mapping keys of type T to tuples containing a vector of `Graph` objects and a vector of external time labels.
- `FermiLabel`: A LabelProduct used to label the fermionic `G` objects in the graphs.
- `BoseLabel`: A LabelProduct used to label bosonic `W` objects in the graphs.
- `graph_keys`: A vector containing keys of type `T`, specifying which graphs to analyze.

# Returns
- A tuple of vectors containing information about the propagators in the graphs, including their initial values, types, input and output time indexes, and loop-momenta indexes.
- A tuple of vectors containing information about the interactions in the graphs, including their initial values, types, input and output time indexes, and loop-momenta indexes.
- A Vector{Vector{Int}} representing the external tau variables of each vector of graph corresponding to each key of type `T`.
"""
function leafstates(FeynGraphs::Dict{T,Tuple{Vector{G},Vector{Vector{Int}}}},
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

        leaves = Vector{G}()
        for graph in FeynGraphs[key][1]
            append!(leaves, collect(Leaves(graph)))
        end
        sort!(leaves, by=x -> x.id) #sort the id of the leaves in an asscend order
        unique!(x -> x.id, leaves) #filter out the leaves with the same id number  

        for g in leaves
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
    return (PropagatorValue, PropagatorType, PropagatorInTau, PropagatorOutTau, PropagatorLoopIndex),
    (InteractionValue, InteractionType, InteractionInTau, InteractionOutTau, InteractionLoopIndex), ExtT_index
end

end