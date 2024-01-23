"""
    function diagdict_parquet(type::Symbol, MaxOrder::Int, has_counterterm::Bool=true; MinOrder::Int=1,
        spinPolarPara::Float64=0.0, isDynamic=false, filter=[NoHartree])

    Generates a Graph Dict: the Feynman diagrams with dynamic/instant interactions in a given `type`, and 
    spin-polarizaition parameter `spinPolarPara`, to given minmimum/maximum orders `MinOrder/MaxOrder`, with switchable couterterms. 

# Arguments:
- `type` (Symbol): The type of the Feynman diagrams, including  `:sigma`, `:chargePolar`, `:green`, `vertex3`, `vertex4`, or `:freeEnergy`.
- `Maxorder` (Int): The maximum actual order of the diagrams.
- `MinOrder` (Int): The minmimum actual order of the diagrams (defaults to `1`).
- `has_counterterm` (Bool, optional): `false` for G0W0, `true` for GW with self-energy and interaction counterterms (defaults to `false`).
- `spinPolarPara` (Float64, optional): The spin-polarization parameter (n_up - n_down) / (n_up + n_down) (defaults to `0.0`).
- `isDynamic` (Bool, optional): Flag to specify if the interactions are dynamic, defaults to false.
- `filter` (optional): Filter criteria for the diagrams, defaults to `[NoHartree]`.

# Returns
- `dict_graphs` is a `Dict{Tuple{Int,Int,Int},Tuple{Vector{Graph},Vector{Vector{Int}}}}` object representing the diagrams. 
   The key is (order, Gorder, Vorder). The element is a Tuple (graphVector, extT_labels).
"""
function diagdict_parquet(type::Symbol, MaxOrder::Int, MinOrder::Int=1;
    has_counterterm::Bool=true, spinPolarPara::Float64=0.0, optimize_level=0, channels=[PHr, PHEr, PPr, Alli],
    isDynamic=false, filter=[NoHartree], transferLoop=nothing, extK=nothing)

    diagtype = _diagtype(type)
    spin = 2.0 / (spinPolarPara + 1)
    dict_graphs = Dict{NTuple{3,Int},Tuple{Vector{Graph},Vector{Vector{Int}}}}()
    # dict_graphs = Dict{Tuple{Int,Int,Int},Tuple{Vector{Graph},Vector{Tuple{Vararg{Int}}}}}()

    if has_counterterm
        for order in MinOrder:MaxOrder
            set_variables("x y"; orders=[MaxOrder - order, MaxOrder - order])
            para = diagPara(diagtype, isDynamic, spin, order, filter, transferLoop)
            parquet_builder = Parquet.build(para, extK, channels=channels)
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
            parquet_builder = Parquet.build(para, extK, channels=channels)
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

function diagdict_parquet(type::Symbol, gkeys::Vector{Tuple{Int,Int,Int}};
    spinPolarPara::Float64=0.0, optimize_level=0, channels=[PHr, PHEr, PPr, Alli],
    isDynamic=false, filter=[NoHartree], transferLoop=nothing, extK=nothing)

    diagtype = _diagtype(type)
    spin = 2.0 / (spinPolarPara + 1)
    dict_graphs = Dict{NTuple{3,Int},Tuple{Vector{Graph},Vector{Vector{Int}}}}()
    diag_orders = unique([p[1] for p in gkeys])
    Gorder = maximum([p[2] for p in gkeys])
    Vorder = maximum([p[3] for p in gkeys])

    for order in diag_orders
        set_variables("x y"; orders=[Gorder, Vorder])
        para = diagPara(diagtype, isDynamic, spin, order, filter, transferLoop)
        parquet_builder = Parquet.build(para, extK, channels=channels)

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

function diagdict_parquet(type::Symbol, gkeys::Vector{Tuple{Int,Int,Int}}, extra_variables::Dict{String,Int};
    spinPolarPara::Float64=0.0, optimize_level=0, channels=[PHr, PHEr, PPr, Alli],
    isDynamic=false, filter=[NoHartree], transferLoop=nothing, extK=nothing)

    diagtype = _diagtype(type)
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
    Gorder = maximum([p[2] for p in gkeys])
    Vorder = maximum([p[3] for p in gkeys])

    for order in diag_orders
        # set_variables("x y k"; orders=[MaxOrder - order, MaxOrder - order, 1])
        set_variables("x y" * extra_varnames; orders=[Gorder, Vorder, extra_orders...])
        para = diagPara(diagtype, isDynamic, spin, order, filter, transferLoop)
        parquet_builder = Parquet.build(para, extK, channels=channels)
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