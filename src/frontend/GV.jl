module GV

import ..QuantumOperators as Op
import ..ComputationalGraphs as IR
import ..ComputationalGraphs: FeynmanGraph
import ..ComputationalGraphs: Graph
import ..ComputationalGraphs: _dtype
import ..Parquet
import ..Taylor
using ..DiagTree
using ..FrontEnds
using AbstractTrees

import ..Utility: taylorexpansion!

import ..Filter
import ..Wirreducible  #remove all polarization subdiagrams
import ..Girreducible  #remove all self-energy inseration
import ..NoHartree
import ..NoFock
import ..NoBubble  # true to remove all bubble subdiagram
import ..Proper  #ver4, ver3, and polarization diagrams may require to be irreducible along the transfer momentum/frequency
import ..DirectOnly

import ..DiagramType
import ..VacuumDiag
import ..GreenDiag
import ..SigmaDiag
import ..PolarDiag
import ..Ver3Diag
import ..Ver4Diag

import ..Composite
import ..ChargeCharge
import ..SpinSpin
import ..UpUp
import ..UpDown
import ..Response

import ..Instant
import ..Dynamic
import ..D_Instant
import ..D_Dynamic
import ..AnalyticProperty

import ..Interaction
import ..DiagPara
import ..DiagParaF64

import ..DiagramId
import ..Ver4Id
import ..Ver3Id
import ..GreenId
import ..SigmaId
import ..PolarId
import ..BareInteractionId
import ..BareGreenId

export eachorder_diag, diagdictGV, diagdict_parquet, leafstates, leafstates_diagtree

include("GV_diagrams/readfile.jl")


"""
    function diagdictGV(type::Symbol, MaxOrder::Int, has_counterterm::Bool=false;
        MinOrder::Int=1, spinPolarPara::Float64=0.0)

    Generates a FeynmanGraph Dict: the Feynman diagrams with static interactions in a given `type`, and 
    spin-polarizaition parameter `spinPolarPara`, to given minmimum/maximum orders `MinOrder/MaxOrder`, with switchable couterterms. 
    Generates a `LabelProduct`: `labelProd` for these FeynmanGraphs.

# Arguments:
- `type` (Symbol): The type of the Feynman diagrams, including `:spinPolar`, `:chargePolar`, `:sigma_old`, `:green`, or `:freeEnergy`.
- `Maxorder` (Int): The maximum actual order of the diagrams.
- `has_counterterm` (Bool): `false` for G0W0, `true` for GW with self-energy and interaction counterterms (defaults to `false`).
- `MinOrder` (Int, optional): The minmimum actual order of the diagrams (defaults to `1`).
- `spinPolarPara` (Float64, optional): The spin-polarization parameter (n_up - n_down) / (n_up + n_down) (defaults to `0.0`).

# Returns
A tuple `(dict_graphs, fermi_labelProd, bose_labelProd, leafMap)` where 
- `dict_graphs` is a `Dict{Tuple{Int,Int,Int},Tuple{Vector{FeynmanGraph},Vector{Vector{Int}}}}` object representing the diagrams. 
   The key is (order, Gorder, Vorder). The element is a Tuple (graphVector, extT_labels).
- `labelProd` is a `LabelProduct` object containing the labels for the leaves of graphs, 
"""
function diagdictGV(type::Symbol, MaxOrder::Int, has_counterterm::Bool=false;
    MinOrder::Int=1, spinPolarPara::Float64=0.0)
    dict_graphs = Dict{Tuple{Int,Int,Int},Tuple{Vector{FeynmanGraph{_dtype.factor,_dtype.weight}},Vector{Vector{Int}}}}()
    if type == :sigma
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
    loopbasis = [vcat([1.0], [0.0 for _ in 2:MaxLoopNum])]
    # Create label product
    labelProd = LabelProduct(tau_labels, loopbasis)

    if has_counterterm
        Gorders = 0:MaxOrder-MinOrder
        Vorders = 0:MaxOrder-1
        for order in MinOrder:MaxOrder
            for VerOrder in Vorders
                type in [:chargePolar, :spinPolar] && order == 1 && VerOrder > 0 && continue
                order == 0 && VerOrder > 0 && continue
                for GOrder in Gorders
                    order + VerOrder + GOrder > MaxOrder && continue
                    gvec, labelProd, extT_labels = eachorder_diag(type, order, GOrder, VerOrder;
                        labelProd=labelProd, spinPolarPara=spinPolarPara)
                    key = (order, GOrder, VerOrder)
                    dict_graphs[key] = (gvec, extT_labels)
                    IR.optimize!(gvec)
                    IR.optimize!(gvec)
                end
            end
        end
    else
        for order in 1:MaxOrder
            gvec, labelProd, extT_labels = eachorder_diag(type, order;
                labelProd=labelProd, spinPolarPara=spinPolarPara)
            key = (order, 0, 0)
            dict_graphs[key] = (gvec, extT_labels)
            IR.optimize!(gvec)
        end
    end

    return dict_graphs, labelProd
end

"""
    function diagdictGV(type::Symbol, gkeys::Vector{Tuple{Int,Int,Int}}; spinPolarPara::Float64=0.0)

    Generates a FeynmanGraph Dict: the Feynman diagrams with static interactions in the given `type` and 
    spin-polarizaition parameter `spinPolarPara`, with given couterterm-orders (from `gkeys`). 
    Generates a `LabelProduct`: `labelProd` for these FeynmanGraphs.

# Arguments:
- `type` (Symbol): The type of the Feynman diagrams, including `:spinPolar`, `:chargePolar`, `:sigma_old`, `:green`, or `:freeEnergy`.
- `gkeys` (Vector{Tuple{Int,Int,Int}}): The (order, Gorder, Vorder) of the diagrams. Gorder is the order of self-energy counterterms, and Vorder is the order of interaction counterterms. 
- `spinPolarPara` (Float64, optional): The spin-polarization parameter (n_up - n_down) / (n_up + n_down) (defaults to `0.0`).

# Returns
A tuple `(dict_graphs, fermi_labelProd, bose_labelProd, leafMap)` where 
- `dict_graphs` is a `Dict{Tuple{Int,Int,Int},Tuple{Vector{FeynmanGraph},Vector{Vector{Int}}}}` object representing the diagrams. 
   The key is (order, Gorder, Vorder). The element is a Tuple (graphVector, extT_labels).
- `labelProd` is a `LabelProduct` object containing the labels for the leaves of graphs, 
"""
function diagdictGV(type::Symbol, gkeys::Vector{Tuple{Int,Int,Int}}; spinPolarPara::Float64=0.0)
    dict_graphs = Dict{Tuple{Int,Int,Int},Tuple{Vector{FeynmanGraph{_dtype.factor,_dtype.weight}},Vector{Vector{Int}}}}()
    if type == :sigma
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

    loopbasis = [vcat([1.0], [0.0 for _ in 2:MaxLoopNum])]
    # Create label product
    labelProd = LabelProduct(tau_labels, loopbasis)

    # graphvector = Vector{_dtype.factor,_dtype.weight}()
    for key in gkeys
        gvec, labelProd, extT_labels = eachorder_diag(type, key...;
            labelProd=labelProd, spinPolarPara=spinPolarPara)
        dict_graphs[key] = (gvec, extT_labels)
        IR.optimize!(gvec)
        IR.optimize!(gvec)
        # append!(graphvector, gvec)
    end
    # IR.optimize!(graphvector)

    return dict_graphs, labelProd
end

"""
    function eachorder_diag(type::Symbol, order::Int, VerOrder::Int=0, GOrder::Int=0; loopPool::Union{LoopPool,Nothing}=nothing,
        tau_labels::Union{Nothing,Vector{Int}}=nothing, GTypes::Union{Nothing,Vector{Int}}=nothing, VTypes::Union{Nothing,Vector{Int}}=nothing)
 
    Generates a `Vector{FeynmanGraph}`: the given-`type` diagrams with static interactions of a given order, where the actual order of diagrams equals to `order + VerOrder + 2 * GOrder`.
    Generates a `LabelProduct`: `labelProd` with inputs `tau_labels` and all the possible momenta-loop basis. 
    Generates external tau labels Vector{Vector{Int}}. The i-th labels (Vector{Int}) corresponds to the i-th `FeynmanGraph` in `Vector{FeynmanGraph}`.

# Arguments:
- `type` (Symbol): The type of the diagrams, including `:spinPolar`, `:chargePolar`, `:sigma`, `:green`, or `:freeEnergy`.
- `order` (Int): The order of the diagrams without counterterms.
- `GOrder` (Int, optional): The order of self-energy counterterms (defaults to 0).
- `VerOrder` (Int, optional): The order of interaction counterterms (defaults to 0).
- `labelProd` (Union{Nothing,LabelProduct}=nothing, optional): The initial cartesian QuantumOperator.label product (defaults to `nothing`).
- `spinPolarPara` (Float64, optional): The spin-polarization parameter (n_up - n_down) / (n_up + n_down) (defaults to `0.0`).
- `tau_labels`(Union{Nothing, Vector{Int}}, optional): The labels for the discrete time of each vertex. (defaults to `nothing`).

# Returns
A tuple `(diagrams, fermi_labelProd, bose_labelProd, extT_labels)` where 
- `diagrams` is a `Vector{FeynmanGraph}` object representing the diagrams, 
- `labelProd` is a `LabelProduct` object containing the labels for the leaves of graphs, 
- `extT_labels` is a `Vector{Vector{Int}}` object containing the external tau labels for each `FeynmanGraph` in `diagrams`.
"""
function eachorder_diag(type::Symbol, order::Int, GOrder::Int=0, VerOrder::Int=0;
    labelProd::Union{Nothing,LabelProduct}=nothing, spinPolarPara::Float64=0.0, tau_labels::Union{Nothing,Vector{Int}}=nothing)
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

    if isnothing(labelProd)
        return read_diagrams(filename; tau_labels=tau_labels, diagType=diagtype, spinPolarPara=spinPolarPara)
    else
        return read_diagrams(filename; labelProd=labelProd, diagType=diagtype, spinPolarPara=spinPolarPara)
    end
end


function diagdict_parquet(type::Symbol, MaxOrder::Int, has_counterterm::Bool=true; MinOrder::Int=1,
    spinPolarPara::Float64=0.0, isDynamic=false, filter=[NoHartree])
    # spinPolarPara::Float64=0.0, isDynamic=false, channel=[PHr, PHEr, PPr], filter=[NoHartree])

    diagtype = _diagtype(type)
    spin = 2.0 / (spinPolarPara + 1)
    dict_graphs = Dict{Tuple{Int,Int,Int},Tuple{Vector{Graph},Vector{Vector{Int}}}}()
    # dict_graphs = Dict{Tuple{Int,Int,Int},Tuple{Vector{Graph},Vector{Tuple{Vararg{Int}}}}}()

    KinL, KoutL, KinR = zeros(16), zeros(16), zeros(16)
    KinL[1], KoutL[2], KinR[3] = 1.0, 1.0, 1.0

    if has_counterterm
        for order in MinOrder:MaxOrder
            Taylor.set_variables("x y"; orders=[MaxOrder - order, MaxOrder - order])
            para = diagPara(diagtype, isDynamic, spin, order, filter, KinL - KoutL)
            # legK = [DiagTree.getK(para.totalLoopNum + 3, 1), DiagTree.getK(para.totalLoopNum + 3, 2), DiagTree.getK(para.totalLoopNum + 3, 3)]
            # d::Vector{Diagram{Float64}} = Parquet.vertex4(para, legK, channel).diagram
            # diags::Vector{Diagram{Float64}} = Parquet.build(para).diagram
            parquet_builder = Parquet.build(para)
            diags, extT = parquet_builder.diagram, parquet_builder.extT

            propagator_var = Dict(BareGreenId => [true, false], BareInteractionId => [false, true]) # Specify variable dependence of fermi (first element) and bose (second element) particles.
            taylor_vec, taylormap = taylorexpansion!(diags, propagator_var)

            for t in taylor_vec
                for (o, graph) in t.coeffs
                    key = (order, o...)
                    sum(key) > MaxOrder && continue
                    if haskey(dict_graphs, key)
                        push!(dict_graphs[key][1], graph)
                    else
                        dict_graphs[key] = ([graph,], collect.(extT))
                    end
                end
            end
        end
    else
        Taylor.set_variables("x y"; orders=[0, 0])
        for order in MinOrder:MaxOrder
            para = diagPara(diagtype, isDynamic, spin, order, filter, KinL - KoutL)
            parquet_builder = Parquet.build(para)
            diags, extT = parquet_builder.diagram, parquet_builder.extT

            propagator_var = Dict(BareGreenId => [true, false], BareInteractionId => [false, true]) # Specify variable dependence of fermi (first element) and bose (second element) particles.
            taylor_vec, taylormap = taylorexpansion!(diags, propagator_var)
            for t in taylor_vec
                graph = t.coeffs[[0, 0]]
                key = (order, 0, 0)
                if haskey(dict_graphs, key)
                    push!(dict_graphs[key][1], graph)
                else
                    dict_graphs[key] = ([graph,], extT)
                end
            end
        end
    end

    for gvec in values(dict_graphs)
        IR.optimize!(gvec[1])
    end
    return dict_graphs
end

function diagdict_parquet(type::Symbol, gkeys::Vector{Tuple{Int,Int,Int}};
    spinPolarPara::Float64=0.0, isDynamic=false, filter=[NoHartree])
    # spinPolarPara::Float64=0.0, isDynamic=false, channel=[PHr, PHEr, PPr], filter=[NoHartree])

    diagtype = _diagtype(type)
    spin = 2.0 / (spinPolarPara + 1)
    dict_graphs = Dict{Tuple{Int,Int,Int},Tuple{Vector{Graph},Vector{Vector{Int}}}}()

    KinL, KoutL, KinR = zeros(16), zeros(16), zeros(16)
    KinL[1], KoutL[2], KinR[3] = 1.0, 1.0, 1.0

    MinOrder = minimum([p[1] for p in gkeys])
    MaxOrder = maximum([p[1] for p in gkeys])
    for order in MinOrder:MaxOrder
        Taylor.set_variables("x y"; orders=[MaxOrder - order, MaxOrder - order])
        para = diagPara(diagtype, isDynamic, spin, order, filter, KinL - KoutL)
        parquet_builder = Parquet.build(para)
        diags, extT = parquet_builder.diagram, parquet_builder.extT

        propagator_var = Dict(BareGreenId => [true, false], BareInteractionId => [false, true]) # Specify variable dependence of fermi (first element) and bose (second element) particles.
        taylor_vec, taylormap = taylorexpansion!(diags, propagator_var)

        for t in taylor_vec
            for (o, graph) in t.coeffs
                key = (order, o...)
                key ∉ gkeys && continue
                if haskey(dict_graphs, key)
                    push!(dict_graphs[key][1], graph)
                else
                    dict_graphs[key] = ([graph,], collect.(extT))
                end
            end
        end
    end

    for gvec in values(dict_graphs)
        IR.optimize!(gvec[1])
    end
    return dict_graphs
end

function diagPara(type, isDynamic::Bool, spin, order, filter, transferLoop)
    # function diagPara(type, isDynamic, spin, order, filter, transferLoop)
    inter = [Interaction(ChargeCharge, isDynamic ? [Instant, Dynamic] : [Instant,]),]  #instant charge-charge interaction
    if type == VacuumDiag
        innerLoopNum = order + 1
    elseif type == Ver4Diag
        innerLoopNum = order - 1
    else
        innerLoopNum = order
    end

    return DiagParaF64(
        type=type,
        innerLoopNum=innerLoopNum,
        hasTau=true,
        spin=spin,
        # firstLoopIdx=4,
        interaction=inter,
        filter=filter,
        transferLoop=transferLoop
    )
end

function _diagtype(type::Symbol)
    if type == :freeEnergy
        diagtype = VacuumDiag
    elseif type == :sigma
        diagtype = SigmaDiag
    elseif type == :green
        diagtype = GreenDiag
    elseif type == :chargePolar
        diagtype = PolarDiag
    elseif type == :vertex3
        diagtype = Ver3Diag
    elseif type == :vertex4
        diagtype = Ver4Diag
    end
end

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

function leafstates_diagtree(leaf_maps::Vector{Dict{Int,Graph}}, maxloopNum::Int)

    num_g = length(leaf_maps)
    leafType = [Vector{Int}() for _ in 1:num_g]
    leafOrders = [Vector{Vector{Int}}() for _ in 1:num_g]
    leafInTau = [Vector{Int}() for _ in 1:num_g]
    leafOutTau = [Vector{Int}() for _ in 1:num_g]
    leafLoopIndex = [Vector{Int}() for _ in 1:num_g]
    leafValue = [Vector{Float64}() for _ in 1:num_g]

    loopbasis = Vector{Float64}[]
    # tau_labels = Vector{Int}[]
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

            # push!(tau_labels, collect(diagId.extT))
            push!(leafInTau[ikey], diagId.extT[1])
            push!(leafOutTau[ikey], diagId.extT[2])

            push!(leafOrders[ikey], leaf_orders)
            push!(leafType[ikey], DiagTree.index(typeof(diagId)))

        end
    end

    return (leafValue, leafType, leafOrders, leafInTau, leafOutTau, leafLoopIndex), loopbasis
end

end