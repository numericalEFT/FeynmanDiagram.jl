"""
    function diagdict_parquet(diagtype::Union{DiagramType,Symbol}, MaxOrder::Int, MinOrder::Int=1;
        has_counterterm::Bool=true, spinPolarPara::Float64=0.0, optimize_level=0, channels=[PHr, PHEr, PPr, Alli],
        isDynamic=false, filter=[NoHartree], transferLoop=nothing, extK=nothing
    )

    Generates a Graph Dict: the Feynman diagrams with dynamic/instant interactions in a given `diagtype`, and 
    spin-polarizaition parameter `spinPolarPara`, to given minmimum/maximum orders `MinOrder/MaxOrder`, with switchable couterterms. 

# Arguments:
- `diagtype` (Union{DiagramType,Symbol}): The type of the Feynman diagrams, including  `:freeEnergy`/`VacuumDiag`, `:sigma`/`SigmaDiag`, 
    `:green`/`GreenDiag`, `:polar`/`PolarDiag`, `:vertex3`/`Ver3Diag`, and `:vertex4`/`Ver4Diag`.
- `Maxorder` (Int): The maximum actual order of the diagrams.
- `MinOrder` (Int): The minmimum actual order of the diagrams (defaults to `1`).
- `has_counterterm` (Bool, optional): `false` for G0W0, `true` for GW with self-energy and interaction counterterms (defaults to `false`).
- `spinPolarPara` (Float64, optional): The spin-polarization parameter (n_up - n_down) / (n_up + n_down) (defaults to `0.0`).
- `optimize_level` (Int, optional): The optimization level for the diagrams, defaults to `0`.
- `channels` (optional): The channels for the diagrams, defaults to `[PHr, PHEr, PPr, Alli]`.
- `isDynamic` (Bool, optional): Flag to specify if the interactions are dynamic, defaults to false.
- `filter` (optional): Filter criteria for the diagrams, defaults to `[NoHartree]`.
- `transferLoop` (optional): Transfer loop for the diagrams, defaults to `nothing`.
- `extK` (optional): External k-vector for the diagrams, defaults to `nothing`.

# Returns
- `dict_graphs` is a `Dict` object representing the diagrams. The key is (order, Gorder, Vorder). 
    The element is a Tuple (graphVector, extT_labels) or (graphVector, extT_labels, responseVector).
    `graphVector::Vector{Graph}` is a vector of `Graph` objects; `extT_labels::Vector{Vector{Int}}` is a vector of external tau labels;
    `responseVector::Vector{Response}` is a vector of `Response` objects.
"""
function diagdict_parquet(diagtype::Union{DiagramType,Symbol}, MaxOrder::Int, MinOrder::Int=1;
    has_counterterm::Bool=true, spinPolarPara::Float64=0.0, optimize_level=0, channels=[PHr, PHEr, PPr, Alli],
    isDynamic=false, filter=[NoHartree], transferLoop=nothing, extK=nothing
)
    if diagtype isa Symbol
        diagtype = _diagtype(diagtype)
    end
    if diagtype in [VacuumDiag, SigmaDiag, GreenDiag]
        return diagdict_parquet_noresponse(diagtype, MaxOrder, MinOrder, has_counterterm=has_counterterm,
            spinPolarPara=spinPolarPara, optimize_level=optimize_level, isDynamic=isDynamic, filter=filter,
            transferLoop=transferLoop, extK=extK)
    else
        return diagdict_parquet_response(diagtype, MaxOrder, MinOrder, has_counterterm=has_counterterm,
            spinPolarPara=spinPolarPara, optimize_level=optimize_level, channels=channels, isDynamic=isDynamic,
            filter=filter, transferLoop=transferLoop, extK=extK)
    end
end

"""
    function diagdict_parquet(diagtype::Union{DiagramType,Symbol}, MaxOrder::Int, MinOrder::Int=1;
        has_counterterm::Bool=true, spinPolarPara::Float64=0.0, optimize_level=0, channels=[PHr, PHEr, PPr, Alli],
        isDynamic=false, filter=[NoHartree], transferLoop=nothing, extK=nothing
    )

    Generates a Graph Dict: the Feynman diagrams with dynamic/instant interactions in a given `diagtype`, and 
    spin-polarizaition parameter `spinPolarPara`, with given couterterm-orders (from `gkeys`).  

# Arguments:
- `diagtype` (Union{DiagramType,Symbol}): The type of the Feynman diagrams, including  `:freeEnergy`/`VacuumDiag`, `:sigma`/`SigmaDiag`, 
    `:green`/`GreenDiag`, `:polar`/`PolarDiag`, `:vertex3`/`Ver3Diag`, and `:vertex4`/`Ver4Diag`.
- `gkeys` (Vector{Tuple{Int,Int,Int}}): The (order, Gorder, Vorder) of the diagrams. Gorder is the order of self-energy counterterms, and Vorder is the order of interaction counterterms. 
- `spinPolarPara` (Float64, optional): The spin-polarization parameter (n_up - n_down) / (n_up + n_down) (defaults to `0.0`).
- `optimize_level` (Int, optional): The optimization level for the diagrams, defaults to `0`.
- `channels` (optional): The channels for the diagrams, defaults to `[PHr, PHEr, PPr, Alli]`.
- `isDynamic` (Bool, optional): Flag to specify if the interactions are dynamic, defaults to false.
- `filter` (optional): Filter criteria for the diagrams, defaults to `[NoHartree]`.
- `transferLoop` (optional): Transfer loop for the diagrams, defaults to `nothing`.
- `extK` (optional): External k-vector for the diagrams, defaults to `nothing`.

# Returns
- `dict_graphs` is a `Dict` object representing the diagrams. The key is (order, Gorder, Vorder). 
    The element is a Tuple (graphVector, extT_labels) or (graphVector, extT_labels, responseVector).
    `graphVector::Vector{Graph}` is a vector of `Graph` objects; `extT_labels::Vector{Vector{Int}}` is a vector of external tau labels;
    `responseVector::Vector{Response}` is a vector of `Response` objects.
"""
function diagdict_parquet(diagtype::Union{DiagramType,Symbol}, gkeys::Vector{Tuple{Int,Int,Int}};
    spinPolarPara::Float64=0.0, optimize_level=0, channels=[PHr, PHEr, PPr, Alli],
    isDynamic=false, filter=[NoHartree], transferLoop=nothing, extK=nothing
)
    if diagtype isa Symbol
        diagtype = _diagtype(diagtype)
    end
    if diagtype in [VacuumDiag, SigmaDiag, GreenDiag]
        return diagdict_parquet_noresponse(diagtype, gkeys,
            spinPolarPara=spinPolarPara, optimize_level=optimize_level, isDynamic=isDynamic, filter=filter,
            transferLoop=transferLoop, extK=extK)
    else
        return diagdict_parquet_response(diagtype, gkeys,
            spinPolarPara=spinPolarPara, optimize_level=optimize_level, channels=channels, isDynamic=isDynamic,
            filter=filter, transferLoop=transferLoop, extK=extK)
    end
end

"""
    function diagdict_parquet_extraAD(diagtype::Union{DiagramType,Symbol}, gkeys::Vector{Tuple{Int,Int,Int}}, extra_variables::Dict{String,Int};
        spinPolarPara::Float64=0.0, optimize_level=0, channels=[PHr, PHEr, PPr, Alli],
        isDynamic=false, filter=[NoHartree], transferLoop=nothing, extK=nothing
    )

    Generates a Graph Dict: the Feynman diagrams with dynamic/instant interactions in a given `diagtype`, and 
    spin-polarizaition parameter `spinPolarPara`, with given couterterm-orders (from `gkeys`) and extra derivative variables (from `extra_variables`).  

# Arguments:
- `diagtype` (Union{DiagramType,Symbol}): The type of the Feynman diagrams, including  `:freeEnergy`/`VacuumDiag`, `:sigma`/`SigmaDiag`, 
    `:green`/`GreenDiag`, `:polar`/`PolarDiag`, `:vertex3`/`Ver3Diag`, and `:vertex4`/`Ver4Diag`.
- `gkeys` (Vector{Tuple{Int,Int,Int}}): The (order, Gorder, Vorder) of the diagrams. Gorder is the order of self-energy counterterms, and Vorder is the order of interaction counterterms. 
- `extra_variables` (Dict{String,Int}): The extra derivative variables, with their orders. 
- `spinPolarPara` (Float64, optional): The spin-polarization parameter (n_up - n_down) / (n_up + n_down) (defaults to `0.0`).
- `optimize_level` (Int, optional): The optimization level for the diagrams, defaults to `0`.
- `channels` (optional): The channels for the diagrams, defaults to `[PHr, PHEr, PPr, Alli]`.
- `isDynamic` (Bool, optional): Flag to specify if the interactions are dynamic, defaults to false.
- `filter` (optional): Filter criteria for the diagrams, defaults to `[NoHartree]`.
- `transferLoop` (optional): Transfer loop for the diagrams, defaults to `nothing`.
- `extK` (optional): External k-vector for the diagrams, defaults to `nothing`.

# Returns
- `dict_graphs` is a `Dict` object representing the diagrams. The key is (order, Gorder, Vorder). 
    The element is a Tuple (graphVector, extT_labels) or (graphVector, extT_labels, responseVector).
    `graphVector::Vector{Graph}` is a vector of `Graph` objects, `extT_labels::Vector{Vector{Int}}` is a vector of external tau labels.
    `responseVector::Vector{Response}` is a vector of `Response` objects.
"""
function diagdict_parquet_extraAD(diagtype::Union{DiagramType,Symbol}, gkeys::Vector{Tuple{Int,Int,Int}}, extra_variables::Dict{String,Int};
    spinPolarPara::Float64=0.0, optimize_level=0, channels=[PHr, PHEr, PPr, Alli],
    isDynamic=false, filter=[NoHartree], transferLoop=nothing, extK=nothing
)
    if diagtype isa Symbol
        diagtype = _diagtype(diagtype)
    end
    if diagtype in [VacuumDiag, SigmaDiag, GreenDiag]
        return diagdict_parquet_noresponse_extraAD(diagtype, gkeys, extra_variables,
            spinPolarPara=spinPolarPara, optimize_level=optimize_level, isDynamic=isDynamic, filter=filter,
            transferLoop=transferLoop, extK=extK)
    else
        return diagdict_parquet_response_extraAD(diagtype, gkeys, extra_variables,
            spinPolarPara=spinPolarPara, optimize_level=optimize_level, channels=channels, isDynamic=isDynamic,
            filter=filter, transferLoop=transferLoop, extK=extK)
    end
end

function diagdict_parquet_noresponse(diagtype::DiagramType, MaxOrder::Int, MinOrder::Int=1;
    has_counterterm::Bool=true, spinPolarPara::Float64=0.0, optimize_level=0,
    isDynamic=false, filter=[NoHartree], transferLoop=nothing, extK=nothing)

    # diagtype = _diagtype(type)
    spin = 2.0 / (spinPolarPara + 1)
    dict_graphs = Dict{NTuple{3,Int},Tuple{Vector{Graph},Vector{Vector{Int}}}}()
    # dict_graphs = Dict{Tuple{Int,Int,Int},Tuple{Vector{Graph},Vector{Tuple{Vararg{Int}}}}}()

    if has_counterterm
        for order in MinOrder:MaxOrder
            set_variables("x y"; orders=[MaxOrder - order, MaxOrder - order])
            para = diagPara(diagtype, isDynamic, spin, order, filter, transferLoop)
            parquet_builder = Parquet.build(para, extK)
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
        set_variables("x y"; orders=[0, 0])
        for order in MinOrder:MaxOrder
            para = diagPara(diagtype, isDynamic, spin, order, filter, transferLoop)
            parquet_builder = Parquet.build(para, extK)
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
            dict_graphs[(order, 0, 0)] = (diags, collect.(extT))
        end
    end

    for gvec in values(dict_graphs)
        IR.optimize!(gvec[1], level=optimize_level)
    end
    return dict_graphs
end

function diagdict_parquet_response(diagtype::DiagramType, MaxOrder::Int, MinOrder::Int=1;
    has_counterterm::Bool=true, spinPolarPara::Float64=0.0, optimize_level=0, channels=[PHr, PHEr, PPr, Alli],
    isDynamic=false, filter=[NoHartree], transferLoop=nothing, extK=nothing)

    # diagtype = _diagtype(type)
    spin = 2.0 / (spinPolarPara + 1)
    dict_graphs = Dict{NTuple{3,Int},Tuple{Vector{Graph},Vector{Vector{Int}},Vector{Response}}}()
    # dict_graphs = Dict{Tuple{Int,Int,Int},Tuple{Vector{Graph},Vector{Tuple{Vararg{Int}}}}}()

    if has_counterterm
        for order in MinOrder:MaxOrder
            set_variables("x y"; orders=[MaxOrder - order, MaxOrder - order])
            para = diagPara(diagtype, isDynamic, spin, order, filter, transferLoop)
            parquet_builder = Parquet.build(para, extK, channels=channels)
            diags, extT, responses = parquet_builder.diagram, parquet_builder.extT, parquet_builder.response

            propagator_var = Dict(BareGreenId => [true, false], BareInteractionId => [false, true]) # Specify variable dependence of fermi (first element) and bose (second element) particles.
            taylor_vec, taylormap = taylorexpansion!(diags, propagator_var)

            for t in taylor_vec
                for (o, graph) in t.coeffs
                    key = (order, o...)
                    sum(key) > MaxOrder && continue
                    if haskey(dict_graphs, key)
                        push!(dict_graphs[key][1], graph)
                    else
                        dict_graphs[key] = ([graph,], collect.(extT), responses)
                    end
                end
            end
        end
    else
        set_variables("x y"; orders=[0, 0])
        for order in MinOrder:MaxOrder
            para = diagPara(diagtype, isDynamic, spin, order, filter, transferLoop)
            parquet_builder = Parquet.build(para, extK, channels=channels)
            diags, extT, responses = parquet_builder.diagram, parquet_builder.extT, parquet_builder.response

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
            dict_graphs[(order, 0, 0)] = (diags, collect.(extT), responses)
        end
    end

    for gvec in values(dict_graphs)
        IR.optimize!(gvec[1], level=optimize_level)
    end
    return dict_graphs
end

function diagdict_parquet_noresponse(diagtype::DiagramType, gkeys::Vector{Tuple{Int,Int,Int}};
    spinPolarPara::Float64=0.0, optimize_level=0,
    isDynamic=false, filter=[NoHartree], transferLoop=nothing, extK=nothing)

    spin = 2.0 / (spinPolarPara + 1)
    dict_graphs = Dict{NTuple{3,Int},Tuple{Vector{Graph},Vector{Vector{Int}}}}()
    diag_orders = unique([p[1] for p in gkeys])
    maxOrder = maximum(diag_orders)

    for order in diag_orders
        GVorder = maxOrder - order
        set_variables("x y"; orders=[GVorder, GVorder])
        para = diagPara(diagtype, isDynamic, spin, order, filter, transferLoop)
        parquet_builder = Parquet.build(para, extK)

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
        IR.optimize!(gvec[1], level=optimize_level)
    end
    return dict_graphs
end

function diagdict_parquet_response(diagtype::DiagramType, gkeys::Vector{Tuple{Int,Int,Int}};
    spinPolarPara::Float64=0.0, optimize_level=0, channels=[PHr, PHEr, PPr, Alli],
    isDynamic=false, filter=[NoHartree], transferLoop=nothing, extK=nothing)

    spin = 2.0 / (spinPolarPara + 1)
    dict_graphs = Dict{NTuple{3,Int},Tuple{Vector{Graph},Vector{Vector{Int}},Vector{Response}}}()
    diag_orders = unique([p[1] for p in gkeys])
    maxOrder = maximum(diag_orders)

    for order in diag_orders
        GVorder = maxOrder - order
        set_variables("x y"; orders=[GVorder, GVorder])
        para = diagPara(diagtype, isDynamic, spin, order, filter, transferLoop)
        parquet_builder = Parquet.build(para, extK, channels=channels)

        diags, extT, responses = parquet_builder.diagram, parquet_builder.extT, parquet_builder.response

        propagator_var = Dict(BareGreenId => [true, false], BareInteractionId => [false, true]) # Specify variable dependence of fermi (first element) and bose (second element) particles.
        taylor_vec, taylormap = taylorexpansion!(diags, propagator_var)

        for t in taylor_vec
            for (o, graph) in t.coeffs
                key = (order, o...)
                key ∉ gkeys && continue
                if haskey(dict_graphs, key)
                    push!(dict_graphs[key][1], graph)
                else
                    dict_graphs[key] = ([graph,], collect.(extT), responses)
                end
            end
        end
    end

    for gvec in values(dict_graphs)
        IR.optimize!(gvec[1], level=optimize_level)
    end
    return dict_graphs
end

function diagdict_parquet_noresponse_extraAD(diagtype::DiagramType, gkeys::Vector{Tuple{Int,Int,Int}}, extra_variables::Dict{String,Int};
    spinPolarPara::Float64=0.0, optimize_level=0, isDynamic=false, filter=[NoHartree], transferLoop=nothing, extK=nothing)

    spin = 2.0 / (spinPolarPara + 1)
    # num_vars = 3 + length(keys(extra_variables))
    dict_graphs = Dict{NTuple{3,Int},Tuple{Vector{Graph},Vector{Vector{Int}}}}()

    extra_varnames = ""
    extra_orders = Int[]
    for (var_name, order) in extra_variables
        extra_varnames *= " $var_name"
        push!(extra_orders, order)
    end
    diag_orders = unique([p[1] for p in gkeys])
    maxOrder = maximum(diag_orders)

    for order in diag_orders
        GVorder = maxOrder - order
        # set_variables("x y k"; orders=[MaxOrder - order, MaxOrder - order, 1])
        set_variables("x y" * extra_varnames; orders=[GVorder, GVorder, extra_orders...])
        para = diagPara(diagtype, isDynamic, spin, order, filter, transferLoop)
        parquet_builder = Parquet.build(para, extK)
        diags, extT = parquet_builder.diagram, parquet_builder.extT

        var_dependence = Dict{Int,Vector{Bool}}()
        for diag in diags
            for leaf in Leaves(diag)
                if leaf.id isa BareGreenId
                    if leaf.id.extK[1] != 0
                        var_dependence[leaf.hash] = [true, false, true]
                    else
                        var_dependence[leaf.hash] = [true, false, false]
                    end
                elseif leaf.id isa BareInteractionId
                    if leaf.id.extK[1] != 0
                        var_dependence[leaf.hash] = [false, true, true]
                    else
                        var_dependence[leaf.hash] = [false, true, false]
                    end
                end
            end
        end

        taylor_vec, taylormap = taylorexpansion!(diags, var_dependence)

        for t in taylor_vec
            for (o, graph) in t.coeffs
                o[3:end] != extra_orders && continue
                key = (order, o[1], o[2])
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
        IR.optimize!(gvec[1], level=optimize_level)
    end
    return dict_graphs
end

function diagdict_parquet_response_extraAD(diagtype::DiagramType, gkeys::Vector{Tuple{Int,Int,Int}}, extra_variables::Dict{String,Int};
    spinPolarPara::Float64=0.0, optimize_level=0, channels=[PHr, PHEr, PPr, Alli],
    isDynamic=false, filter=[NoHartree], transferLoop=nothing, extK=nothing)

    spin = 2.0 / (spinPolarPara + 1)
    # num_vars = 3 + length(keys(extra_variables))
    dict_graphs = Dict{NTuple{3,Int},Tuple{Vector{Graph},Vector{Vector{Int}},Vector{Response}}}()

    extra_varnames = ""
    extra_orders = Int[]
    for (var_name, order) in extra_variables
        extra_varnames *= " $var_name"
        push!(extra_orders, order)
    end
    diag_orders = unique([p[1] for p in gkeys])
    maxOrder = maximum(diag_orders)

    for order in diag_orders
        GVorder = maxOrder - order
        # set_variables("x y k"; orders=[MaxOrder - order, MaxOrder - order, 1])
        set_variables("x y" * extra_varnames; orders=[GVorder, GVorder, extra_orders...])
        para = diagPara(diagtype, isDynamic, spin, order, filter, transferLoop)
        parquet_builder = Parquet.build(para, extK, channels=channels)
        diags, extT, responses = parquet_builder.diagram, parquet_builder.extT, parquet_builder.response

        var_dependence = Dict{Int,Vector{Bool}}()
        for diag in diags
            for leaf in Leaves(diag)
                if leaf.id isa BareGreenId
                    if leaf.id.extK[1] != 0
                        var_dependence[leaf.hash] = [true, false, true]
                    else
                        var_dependence[leaf.hash] = [true, false, false]
                    end
                elseif leaf.id isa BareInteractionId
                    if leaf.id.extK[1] != 0
                        var_dependence[leaf.hash] = [false, true, true]
                    else
                        var_dependence[leaf.hash] = [false, true, false]
                    end
                end
            end
        end

        taylor_vec, taylormap = taylorexpansion!(diags, var_dependence)

        for t in taylor_vec
            for (o, graph) in t.coeffs
                o[3:end] != extra_orders && continue
                key = (order, o[1], o[2])
                key ∉ gkeys && continue
                if haskey(dict_graphs, key)
                    push!(dict_graphs[key][1], graph)
                else
                    dict_graphs[key] = ([graph,], collect.(extT), responses)
                end
            end
        end
    end

    for gvec in values(dict_graphs)
        IR.optimize!(gvec[1], level=optimize_level)
    end
    return dict_graphs
end

function _diagtype(type::Symbol)
    if type == :freeEnergy
        return VacuumDiag
    elseif type == :sigma
        return SigmaDiag
    elseif type == :green
        return GreenDiag
    elseif type == :chargePolar
        return PolarDiag
    elseif type == :vertex3
        return Ver3Diag
    elseif type in [:vertex4, :vertex4I]
        return Ver4Diag
    else
        error("$type is not implemented")
    end
end

function diagPara(type, isDynamic::Bool, spin, order, filter, transferLoop=nothing)
    inter = [Interaction(ChargeCharge, isDynamic ? [Instant, Dynamic] : [Instant,]),]  #instant charge-charge interaction
    if type == VacuumDiag
        innerLoopNum = order + 1
    else
        innerLoopNum = order
    end

    if Proper in filter
        @assert !isnothing(transferLoop) "transferLoop must be provided if Proper is in filter"
    end

    if isnothing(transferLoop)
        return DiagPara(
            type=type,
            innerLoopNum=innerLoopNum,
            hasTau=true,
            spin=spin,
            interaction=inter,
            filter=filter,
        )
    else
        return DiagPara(
            type=type,
            innerLoopNum=innerLoopNum,
            hasTau=true,
            spin=spin,
            interaction=inter,
            filter=filter,
            transferLoop=transferLoop
        )
    end
end