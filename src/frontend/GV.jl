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
        # append!(graphvector, gvec)
    end
    # IR.optimize!(graphvector)

    return dict_graphs, labelProd
end

"""
    function leafstates(leaf_maps::Vector{Dict{Int,G}}, labelProd::LabelProduct) where {T,G<:FeynmanGraph}

    Extracts leaf information from the leaf mapping from the leaf value's index to the leaf node for all graph partitions
    and their associated LabelProduct data (`labelProd`). 
    The information includes their initial value, type, in/out time, and loop momenta.
    
# Arguments:
- `leaf_maps`: A vector of the dictionary mapping the leaf value's index to the FeynmanGraph of this leaf. 
               Each dict corresponds to a graph partition, such as (order, Gorder, Vorder).
- `labelProd`: A LabelProduct used to label the leaves of graphs.

# Returns
- A tuple of vectors containing information about the leaves of graphs, including their initial values, types, input and output time indexes, and loop-momenta indexes.
"""
function leafstates(leaf_maps::Vector{Dict{Int,G}}, labelProd::LabelProduct) where {G<:FeynmanGraph}
    #read information of each leaf from the generated graph and its LabelProduct, the information include type, loop momentum, imaginary time.
    num_g = length(leaf_maps)
    leafType = [Vector{Int}() for _ in 1:num_g]
    leafInTau = [Vector{Int}() for _ in 1:num_g]
    leafOutTau = [Vector{Int}() for _ in 1:num_g]
    leafLoopIndex = [Vector{Int}() for _ in 1:num_g]
    leafValue = [Vector{Float64}() for _ in 1:num_g]

    for (ikey, leafmap) in enumerate(leaf_maps)
        len_leaves = length(keys(leafmap))
        sizehint!(leafType[ikey], len_leaves)
        sizehint!(leafInTau[ikey], len_leaves)
        sizehint!(leafOutTau[ikey], len_leaves)
        sizehint!(leafLoopIndex[ikey], len_leaves)
        leafValue[ikey] = ones(Float64, len_leaves)

        for idx in 1:len_leaves
            g = leafmap[idx]
            vertices = IR.vertices(g)
            if IR.diagram_type(g) == IR.Interaction
                push!(leafType[ikey], 0)
                In = Out = vertices[1][1].label
                push!(leafLoopIndex[ikey], 1)
                push!(leafInTau[ikey], labelProd[In][1])
                push!(leafOutTau[ikey], labelProd[Out][1])
            elseif IR.diagram_type(g) == IR.Propagator
                if (Op.isfermionic(vertices[1]))
                    In, Out = vertices[2][1].label, vertices[1][1].label
                    push!(leafType[ikey], g.orders[1] * 2 + 1)
                    push!(leafLoopIndex[ikey], FrontEnds.linear_to_index(labelProd, In)[end]) #the label of LoopPool for each fermionic leaf
                    push!(leafInTau[ikey], labelProd[In][1])
                    push!(leafOutTau[ikey], labelProd[Out][1])
                else
                    In, Out = vertices[2][1].label, vertices[1][1].label
                    push!(leafType[ikey], g.orders[2] * 2 + 2)
                    push!(leafLoopIndex[ikey], FrontEnds.linear_to_index(labelProd, In)[end]) #the label of LoopPool for each bosonic leaf
                    push!(leafInTau[ikey], labelProd[In][1])
                    push!(leafOutTau[ikey], labelProd[Out][1])
                end
            end
        end
    end
    return (leafValue, leafType, leafInTau, leafOutTau, leafLoopIndex)
end


function diagPara(type, isDynamic, spin, order, filter, transferLoop)
    inter = [FeynmanDiagram.Interaction(ChargeCharge, isDynamic ? [Instant, Dynamic] : [Instant,]),]  #instant charge-charge interaction
    return DiagParaF64(
        type=type,
        innerLoopNum=order - 1,
        hasTau=true,
        spin=spin,
        # firstLoopIdx=4,
        interaction=inter,
        filter=filter,
        transferLoop=transferLoop
    )
end

# function parquet(type::Symbol, gkeys::Vector{Tuple{Int,Int,Int}}; spinPolarPara::Float64=0.0,
function parquet(type::Symbol, MaxOrder::Int, has_counterterm::Bool=true; MinOrder::Int=1,
    spinPolarPara::Float64=0.0, isDynamic=false, filter=[NoHartree])
    # spinPolarPara::Float64=0.0, isDynamic=false, channel=[PHr, PHEr, PPr], filter=[NoHartree])

    if type == :freeEnergy
        diagtype = VaccumDiag
    elseif type == :sigma
        diagtype = SigmaDiag
    elseif type == :green
        diagtype = :GreenDiag
    elseif type == :chargePolar
        diagtype = PolarDiag
    end

    spin = 2.0 / (spinPolarPara + 1)
    diagpara = Vector{DiagParaF64}()
    dict_graphs = Dict{Tuple{Int,Int,Int},Tuple{Vector{Graph{_dtype.factor,_dtype.weight}},Vector{Vector{Int}}}}()

    KinL, KoutL, KinR = zeros(16), zeros(16), zeros(16)
    KinL[1], KoutL[2], KinR[3] = 1.0, 1.0, 1.0

    if has_counterterm
        set_variables("x y"; orders=[0, 0])
        for order in MinOrder:MaxOrder
            para = diagpara(isDynamic, spin, order, filter, KinL - KoutL)
            # legK = [DiagTree.getK(para.totalLoopNum + 3, 1), DiagTree.getK(para.totalLoopNum + 3, 2), DiagTree.getK(para.totalLoopNum + 3, 3)]
            # d::Vector{Diagram{Float64}} = Parquet.vertex4(para, legK, channel).diagram
            d::Vector{Diagram{Float64}} = Parquet.build(para).diagram
            propagator_var = Dict(DiagTree.BareGreenId => [true, false], DiagTree.BareInteractionId => [false, true]) # Specify variable dependence of fermi (first element) and bose (second element) particles.
            t, taylormap, from_coeff_map = taylorexpansion!(root, propagator_var)
            # t, taylormap = taylorexpansion!(root, propagator_var)
        end
    else
        for order in MinOrder:MaxOrder
            set_variables("x y"; orders=[MaxOrder - order, MaxOrder - order])
            para = diagpara(isDynamic, spin, order, filter, KinL - KoutL)
            # diagpara = Vector{DiagParaF64}()
            # legK = [DiagTree.getK(para.totalLoopNum + 3, 1), DiagTree.getK(para.totalLoopNum + 3, 2), DiagTree.getK(para.totalLoopNum + 3, 3)]
            # d::Vector{Diagram{Float64}} = Parquet.vertex4(para, legK, channel).diagram
            d::Vector{Diagram{Float64}} = Parquet.vertex4(para).diagram

            propagator_var = Dict(DiagTree.BareGreenId => [true, false], DiagTree.BareInteractionId => [false, true]) # Specify variable dependence of fermi (first element) and bose (second element) particles.
            t, taylormap, from_coeff_map = taylorexpansion!(root, propagator_var)
        end
    end

    return t, taylormap, from_coeff_map
end

end