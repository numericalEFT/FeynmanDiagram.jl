module GV

import ..QuantumOperators as Op
import ..ComputationalGraphs as IR
import ..ComputationalGraphs: FeynmanGraph
import ..ComputationalGraphs: _dtype
using ..FrontEnds
using AbstractTrees

export eachorder_diag, diagdictGV, leafstates

include("GV_diagrams/readfile.jl")

"""
    function eachorder_diag(type::Symbol, order::Int, VerOrder::Int=0, GOrder::Int=0; loopPool::Union{LoopPool,Nothing}=nothing,
        tau_labels::Union{Nothing,Vector{Int}}=nothing, GTypes::Union{Nothing,Vector{Int}}=nothing, VTypes::Union{Nothing,Vector{Int}}=nothing)
 
    Generates a `Vector{FeynmanGraph}`: the polarization diagrams with static interactions of a given order, where the actual order of diagrams equals to `order + VerOrder + 2 * GOrder`.
    Generates fermionic/bosonic `LabelProduct`: `fermi_labelProd`/`bose_labelProd` with inputs `tau_labels`, `GTypes`/`VTypes`, and updated `loopPool`. 
    Generates external tau labels Vector{Vector{Int}}. The i-th labels (Vector{Int}) corresponds to the i-th `FeynmanGraph` in `Vector{FeynmanGraph}`.

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
- `diagrams` is a `Vector{FeynmanGraph}` object representing the diagrams, 
- `fermi_labelProd` is a `LabelProduct` object containing the labels for the fermionic `G` objects in the diagrams, 
- `bose_labelProd` is a `LabelProduct` object containing the labels for the bosonic `W` objects in the diagrams.
- `extT_labels` is a `Vector{Vector{Int}}` object containing the external tau labels for each `FeynmanGraph` in `diagrams`.
"""
function eachorder_diag(type::Symbol, labelProd::LabelProduct, order::Int,
    GOrder::Int=0, VerOrder::Int=0; spinPolarPara::Float64=0.0, tau_labels::Union{Nothing,Vector{Int}}=nothing)
    # GTypes::Union{Nothing,Vector{Int}}=nothing, VTypes::Union{Nothing,Vector{Int}}=nothing)
    diagtype = :polar
    if type == :spinPolar
        filename = string(@__DIR__, "/GV_diagrams/groups_spin/Polar$(order)_$(VerOrder)_$(GOrder).diag")
    elseif type == :chargePolar
        filename = string(@__DIR__, "/GV_diagrams/groups_charge/Polar$(order)_$(VerOrder)_$(GOrder).diag")
    elseif type == :sigma_old
        diagtype = type
        filename = string(@__DIR__, "/GV_diagrams/groups_sigma_old/Sigma$(order)_$(VerOrder)_$(GOrder).diag")
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

    Gorders = collect(0:GOrder)
    # type == :sigma_old && append!(Gorders, [-2, -3])
    type in [:green, :sigma] && pushfirst!(Gorders, -2)
    Vorders = collect(0:VerOrder)
    # isnothing(VTypes) && (VTypes = collect(0:VerOrder))
    GVorders = [Gorders, Vorders]

    return read_diagrams(filename, GVorders; labelProd=labelProd, tau_labels=tau_labels, diagType=diagtype, spinPolarPara=spinPolarPara)
end

"""
    function diagdictGV(type::Symbol, MaxOrder::Int, has_counterterm::Bool=false, dim::Int=3)

    Generates a FeynmanGraph Dict: the `dim`-dimensional spin/charge polarization or self-energy diagrams with static interactions in a given `type`, to a given maximum order `MaxOrder`, with switchable couterterms. 
    Generates fermionic/bosonic `LabelProduct`: `fermi_labelProd`/`bose_labelProd` for these FeynmanGraphs.
    Generates a leafMap for mapping `g.id` to the index of unique leaf.

# Arguments:
- `type` (Symbol): The type of the Feynman diagrams, including `:spinPolar`, `:chargePolar`, `:sigma_old`, `:green`, or `:freeEnergy`.
- `Maxorder` (Int): The maximum actual order of the diagrams.
- `has_counterterm` (Bool): `false` for G0W0, `true` for GW with self-energy and interaction counterterms (defaults to `false`).
- `dim` (Int): The dimension of the system (defaults to 3).
- `spinPolarPara` (Float64, optional): The spin-polarization parameter (n_up - n_down) / (n_up + n_down) (defaults to `0.0`).

# Returns
A tuple `(dict_graphs, fermi_labelProd, bose_labelProd, leafMap)` where 
- `dict_graphs` is a `Dict{Tuple{Int,Int,Int},Tuple{Vector{FeynmanGraph},Vector{Vector{Int}}}}` object representing the diagrams. 
   The key is (order, Gorder, Vorder). The element is a Tuple (graphVector, extT_labels).
- `fermi_labelProd` is a `LabelProduct` object containing the labels for the fermionic `G` objects in the diagrams, 
- `bose_labelProd` is a `LabelProduct` object containing the labels for the bosonic `W` objects in the diagrams.
- `leafMap` maps `g.id` to the index of unique leaf. 
"""
function diagdictGV(type::Symbol, MaxOrder::Int, has_counterterm::Bool=false, dim::Int=3;
    MinOrder::Int=1, spinPolarPara::Float64=0.0)
    dict_graphs = Dict{Tuple{Int,Int,Int},Tuple{Vector{FeynmanGraph{_dtype.factor,_dtype.weight}},Vector{Vector{Int}}}}()
    if type == :sigma_old
        MaxLoopNum = MaxOrder + 2
        tau_labels = collect(1:MaxLoopNum)
    elseif type == :sigma
        MaxLoopNum = MaxOrder + 1
        tau_labels = collect(1:MaxLoopNum-1)
    elseif type in [:chargePolar, :spinPolar, :green]
        MaxLoopNum = MaxOrder + 1
        tau_labels = collect(1:MaxLoopNum)
        if type == :spinPolar
            @assert iszero(spinPolarPara) "no support for the spin polarization in the spin-polarized systems"
        end
    elseif type == :freeEnergy
        MaxLoopNum = MaxOrder + 1
        tau_labels = collect(1:MaxLoopNum-1)
        MaxLoopNum == 1 && (tau_labels = [1])  # must set a tau label
    else
        error("no support for $type diagram")
    end
    # loopPool = LoopPool(:K, dim, MaxLoopNum, Float64)
    loopbasis = [vcat([1.0], [0.0 for _ in 2:MaxLoopNum])]
    # Create label product
    fermi_labelProd = LabelProduct(tau_labels, GTypes, loopbasis)
    bose_labelProd = LabelProduct(tau_labels, VTypes, loopbasis)


    leafMap = Dict{Tuple{Int,Int,Int},Dict{Int,Int}}()
    if has_counterterm
        GTypes = collect(0:MaxOrder-MinOrder)
        type == :sigma_old && append!(GTypes, [-2, -3])
        type in [:green, :sigma] && push!(GTypes, -2)
        type == :freeEnergy && push!(GTypes, -1)
        VTypes = collect(0:MaxOrder-1)
        for order in MinOrder:MaxOrder
            for VerOrder in VTypes
                type in [:chargePolar, :spinPolar] && order == 1 && VerOrder > 0 && continue
                order == 0 && VerOrder > 0 && continue
                for GOrder in collect(0:MaxOrder-MinOrder)
                    order + VerOrder + GOrder > MaxOrder && continue
                    gvec, fermi_labelProd, bose_labelProd, extT_labels = eachorder_diag(type, order, GOrder, VerOrder;
                        tau_labels=tau_labels, GTypes=GTypes, VTypes=VTypes, spinPolarPara=spinPolarPara)
                    key = (order, GOrder, VerOrder)
                    dict_graphs[key] = (gvec, extT_labels)
                    # loopPool = fermi_labelProd.labels[3]
                    leafMap[key] = IR.optimize!(gvec)
                end
            end
        end
    else
        GTypes, VTypes = [0], [0]
        type == :sigma_old && append!(GTypes, [-2, -3])
        type in [:green, :sigma] && push!(GTypes, -2)
        for order in 1:MaxOrder
            gvec, fermi_labelProd, bose_labelProd, extT_labels = eachorder_diag(type, order;
                tau_labels=tau_labels, GTypes=GTypes, VTypes=VTypes, spinPolarPara=spinPolarPara)
            key = (order, 0, 0)
            dict_graphs[key] = (gvec, extT_labels)
            # loopPool = fermi_labelProd.labels[3]
            leafMap[key] = IR.optimize!(gvec)
        end
    end
    # fermi_labelProd = LabelProduct(tau_labels, GTypes, loopPool)
    # bose_labelProd = LabelProduct(tau_labels, VTypes, loopPool)

    return dict_graphs, fermi_labelProd, bose_labelProd, leafMap
end

"""
    function diagdictGV(type::Symbol, MaxOrder::Int, has_counterterm::Bool=false, dim::Int=3)

    Generates a FeynmanGraph Dict: the `dim`-dimensional spin/charge polarization or self-energy diagrams with static interactions in a given `type`, to a given maximum order `MaxOrder`, with switchable couterterms. 
    Generates fermionic/bosonic `LabelProduct`: `fermi_labelProd`/`bose_labelProd` for these FeynmanGraphs.
    Generates a leafMap for mapping `g.id` to the index of unique leaf.

# Arguments:
- `type` (Symbol): The type of the Feynman diagrams, including `:spinPolar`, `:chargePolar`, `:sigma_old`, `:green`, or `:freeEnergy`.
- `gkeys` (Vector{Tuple{Int,Int,Int}}): The (order, Gorder, Vorder) of the diagrams. Gorder is the order of self-energy counterterms, and Vorder is the order of interaction counterterms. 
- `dim` (Int): The dimension of the system (defaults to 3).
- `spinPolarPara` (Float64, optional): The spin-polarization parameter (n_up - n_down) / (n_up + n_down) (defaults to `0.0`).

# Returns
A tuple `(dict_graphs, fermi_labelProd, bose_labelProd, leafMap)` where 
- `dict_graphs` is a `Dict{Tuple{Int,Int,Int},Tuple{Vector{FeynmanGraph},Vector{Vector{Int}}}}` object representing the diagrams. 
   The key is (order, Gorder, Vorder). The element is a Tuple (graphVector, extT_labels).
- `fermi_labelProd` is a `LabelProduct` object containing the labels for the fermionic `G` objects in the diagrams, 
- `bose_labelProd` is a `LabelProduct` object containing the labels for the bosonic `W` objects in the diagrams.
- `leafMap` maps `g.id` to the index of unique leaf. 
"""
function diagdictGV(type::Symbol, gkeys::Vector{Tuple{Int,Int,Int}}, dim::Int=3; spinPolarPara::Float64=0.0)
    dict_graphs = Dict{Tuple{Int,Int,Int},Tuple{Vector{FeynmanGraph{_dtype.factor,_dtype.weight}},Vector{Vector{Int}}}}()
    if type == :sigma_old
        MaxLoopNum = maximum([key[1] for key in gkeys]) + 2
        tau_labels = collect(1:MaxLoopNum)
    elseif type == :sigma
        MaxLoopNum = maximum([key[1] for key in gkeys]) + 1
        tau_labels = collect(1:MaxLoopNum-1)
    elseif type in [:chargePolar, :spinPolar, :green]
        MaxLoopNum = maximum([key[1] for key in gkeys]) + 1
        tau_labels = collect(1:MaxLoopNum)
        if type == :spinPolar
            @assert iszero(spinPolarPara) "no support for the spin polarization in the spin-polarized systems"
        end
    elseif type == :freeEnergy
        MaxLoopNum = maximum([key[1] for key in gkeys]) + 1
        tau_labels = collect(1:MaxLoopNum-1)
        MaxLoopNum == 1 && (tau_labels = [1])  # must set a tau label
    else
        error("no support for $type diagram")
    end
    MaxGOrder = maximum([key[2] for key in gkeys])
    MaxVerOrder = maximum([key[3] for key in gkeys])

    loopPool = LoopPool(:K, dim, MaxLoopNum, Float64)
    GTypes = collect(0:MaxGOrder)
    type == :sigma_old && append!(GTypes, [-2, -3])
    type in [:green, :sigma] && push!(GTypes, -2)
    type == :freeEnergy && push!(GTypes, -1)
    VTypes = collect(0:MaxVerOrder)

    # graphvector = Vector{_dtype.factor,_dtype.weight}()
    leafMap = Dict{eltype(gkeys),Dict{Int,Int}}()
    for key in gkeys
        gvec, fermi_labelProd, bose_labelProd, extT_labels = eachorder_diag(type, key...;
            dim=dim, loopPool=loopPool, tau_labels=tau_labels, GTypes=GTypes, VTypes=VTypes, spinPolarPara=spinPolarPara)
        dict_graphs[key] = (gvec, extT_labels)
        loopPool = fermi_labelProd.labels[3]
        leafMap[key] = IR.optimize!(gvec)
        # append!(graphvector, gvec)
    end
    # IR.optimize!(graphvector)

    fermi_labelProd = LabelProduct(tau_labels, GTypes, loopPool)
    bose_labelProd = LabelProduct(tau_labels, VTypes, loopPool)
    return dict_graphs, fermi_labelProd, bose_labelProd, leafMap
end

"""
    function leafstates(
        FeynGraphs::Dict{T, Tuple{Vector{G}, Vector{Vector{Int}}}},
        FermiLabel::LabelProduct, BoseLabel::LabelProduct,
        graph_keys::Vector{T}
    ) where {T, G<:FeynmanGraph}

    Extracts leaf information from a Dict collection of Feynman graphs (`FeynGraphs` with its keys `graph_keys`)
    and their associated LabelProduct data (`FermiLabel` and `BoseLabel`). 
    The information includes their initial value, type, in/out time, and loop momenta.
    
# Arguments:
- `FeynGraphs`: A dictionary mapping keys of type T to tuples containing a vector of `FeynmanGraph` objects and a vector of external time labels.
- `FermiLabel`: A LabelProduct used to label the fermionic `G` objects in the graphs.
- `BoseLabel`: A LabelProduct used to label bosonic `W` objects in the graphs.
- `graph_keys`: A vector containing keys of type `T`, specifying which graphs to analyze.

# Returns
- A tuple of vectors containing information about the leaves in the graphs, including their initial values, types, input and output time indexes, and loop-momenta indexes.
- A Vector{Vector{Int}} representing the external tau variables of each vector of graph corresponding to each key of type `T`.
"""
function leafstates(FeynGraphs::Dict{T,Tuple{Vector{G},Vector{Vector{Int}}}},
    FermiLabel::LabelProduct, BoseLabel::LabelProduct, graph_keys::Vector{T}) where {T,G<:FeynmanGraph}
    #read information of each leaf from the generated graph and its LabelProduct, the information include type, loop momentum, imaginary time.
    num_g = length(graph_keys)
    ExtT_index = [Vector{Vector{Int}}() for _ in 1:num_g]

    leafType = [Vector{Int}() for _ in 1:num_g]
    leafInTau = [Vector{Int}() for _ in 1:num_g]
    leafOutTau = [Vector{Int}() for _ in 1:num_g]
    leafLoopIndex = [Vector{Int}() for _ in 1:num_g]
    leafValue = [Vector{Float64}() for _ in 1:num_g]

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
            vertices = IR.vertices(g)
            if IR.diagram_type(g) == IR.Interaction
                push!(leafType[ikey], 0)
                In = Out = vertices[1][1].label
                push!(leafLoopIndex[ikey], 1)
                push!(leafInTau[ikey], FermiLabel[In][1])
                push!(leafOutTau[ikey], FermiLabel[Out][1])
                push!(leafValue[ikey], 1.0)
            elseif IR.diagram_type(g) == IR.Propagator
                if (Op.isfermionic(vertices[1]))
                    In, Out = vertices[2][1].label, vertices[1][1].label
                    if FermiLabel[In][2] in [-2, -3]
                        push!(leafType[ikey], 0)
                        push!(leafLoopIndex[ikey], 1)
                    else
                        push!(leafType[ikey], FermiLabel[In][2] * 2 + 1)
                        push!(leafLoopIndex[ikey], FrontEnds.linear_to_index(FermiLabel, In)[end]) #the label of LoopPool for each fermionic leaf
                    end
                    push!(leafInTau[ikey], FermiLabel[In][1])
                    push!(leafOutTau[ikey], FermiLabel[Out][1])
                else
                    In, Out = vertices[2][1].label, vertices[1][1].label
                    push!(leafType[ikey], BoseLabel[In][2] * 2 + 2)
                    push!(leafLoopIndex[ikey], FrontEnds.linear_to_index(BoseLabel, In)[end]) #the label of LoopPool for each bosonic leaf
                    push!(leafInTau[ikey], BoseLabel[In][1])
                    push!(leafOutTau[ikey], BoseLabel[Out][1])
                end
                push!(leafValue[ikey], 1.0)
            end
            g.name = "visited"
        end
    end
    return (leafValue, leafType, leafInTau, leafOutTau, leafLoopIndex), ExtT_index
end

end