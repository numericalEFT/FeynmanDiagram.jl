function _combinegroups(groups, getid, operator, name)
    # combine diagrams in a group into one composite diagram
    gdf = combine(groups) do group # for each group in groups
        # check the documentation of ``combine" for details https://dataframes.juliadata.org/stable/man/split_apply_combine/
        # id = isnothing(getid) ? GenericId(group.diagram[1].id.para, Tuple(group[1, fields])) : getid(group)
        id = getid(group)

        if nrow(group) == 1
            # if there is only one diagram in df, and the new id is either GenericId or the id of the existing diagram, 
            # then simply return the current df without creating a new diagram
            # ! the new factor will be multiplied to the factor of the exisiting diagram!
            if id isa GenericId || typeof(id) == typeof(group.diagram[1].properties)
                # diag = deepcopy(group[1, :diagram])
                diag = group.diagram[1]
                return (diagram=diag, hash=diag.id)
            end
        end
        diag = Graph(group.diagram; properties=id, operator=operator, name=name)
        return (diagram=diag, hash=diag.id)
    end
    return gdf
end

function _mergediag(group, id, operator, name)
    if nrow(group) == 1
        # if there is only one diagram in df, and the new id is either GenericId or the id of the existing diagram, 
        # then simply return the current df without creating a new diagram
        # ! the new factor will be multiplied to the factor of the exisiting diagram!
        if id isa GenericId || typeof(id) == typeof(group.diagram[1].properties)
            # diag = deepcopy(group[1, :diagram])
            diag = group.diagram[1]
            return diag
        end
    end
    return Graph(group.diagram; properties=id, operator=operator, name=name)
end

function _combine(groups, getid, operator, name)
    """
    # if fields = [:response, :extT], then

    # 1. groups.cols is like: Vector{Symbol}[:response, :extT]

    # 2. groups.keymap is like: 

    #     Dict{Any, Int64} with 2 entries:
    #     (UpDown, (1, 1, 1, 1)) => 2
    #     (UpUp, (1, 1, 1, 1))   => 1
    # """
    d = Dict{Symbol,Any}()
    _keys = keys(groups)
    for col in groupcols(groups)
        d[col] = [key[col] for key in _keys]
    end
    d[:diagram] = [_mergediag(groups[key], getid(groups[key]), operator, name) for key in _keys]
    d[:hash] = [diag.id for diag in d[:diagram]]
    return DataFrame(d, copycols=false)
end

# function mergeby(df::DataFrame, fields=Vector{Symbol}();
#     operator=Sum(), name::Symbol=:none,
#     getid::Function=g -> GenericId(g[1, :diagram].properties.para, Tuple(g[1, fields]))
# )
#     if isempty(df)
#         return df
#     else
#         return mergeby(df, fields; operator=operator, name=name, getid=getid)
#     end
# end

function mergeby(df::DataFrame, fields=Vector{Symbol}();
    operator=Sum(), name::Symbol=:none,
    getid::Function=g -> GenericId(g[1, :diagram].properties.para, Tuple(g[1, fields]))
)
    if isempty(df)
        return df
    else
        if all(x -> typeof(x.properties) == typeof(df.diagram[1].properties), df[!, :diagram]) == false
            @warn "Not all DiagramIds in $df are the same!"
        end
        groups = DataFrames.groupby(df, fields, sort=true)
        ########  less memory usage but can not pass the test right now ##############
        d = _combine(groups, getid, operator, name)
        ######## alternative approach (more memory)  ##################
        # d = _combinegroups(groups, getid, factor, operator, name)
        # println("old\n$d \n new\n$cd")
        return d
    end
end

function mergeby(diags::Union{Graph,Tuple,AbstractVector}, fields=nothing; idkey=nothing, kwargs...)
    if diags isa Graph
        return diags
    else
        if isempty(diags)
            return diags
        else
            W = typeof(diags[1].weight)
            @assert all(x -> (x.weight isa W), diags) "all diagrams should be of the same type. \n$diags"
            diags = collect(diags)
            if isnothing(fields) && isnothing(idkey)
                return mergeby(diags; kwargs...)
            else
                return mergeby(diags, fields; idkey=idkey, kwargs...)
            end
        end
    end
end

# function mergeby(diags::AbstractVector, fields=[]; idkey::Vector{Symbol}=[], kwargs...)
function mergeby(diags::Vector{Graph{F,W}}, fields; idkey=Vector{Symbol}(), kwargs...) where {F,W}
    if isempty(diags)
        return diags
    else
        df = toDataFrame(diags, idkey)
        mergedf = mergeby(df, fields; kwargs...)
        return Vector{Graph{F,W}}(mergedf.diagram)
    end
end

function mergeby(diags::Vector{Graph{F,W}};
    operator=Sum(), name::Symbol=:none,
    getid::Function=d -> GenericId(d[1].properties.para)
) where {F,W}
    if isempty(diags)
        return diags
    else
        id = getid(diags)
        if length(diags) == 1 && (id isa GenericId || typeof(id) == typeof(diags[1].properties))
            # if there is only one diagram, and the new id is either GenericId or the id of the existing diagram, 
            # then simply return the current diagram without creating a new diagram
            # ! the new factor will be multiplied to the factor of the exisiting diagram!
            return diags
        end
        diag = Graph(diags; properties=id, operator=operator, name=name)
        return [diag,]
    end
end
# mergeby(df::DataFrame; kwargs...) = mergeby(df, []; kwargs...)
# mergeby(diags::Vector{Graph}; kwargs...) = mergeby(diags, []; kwargs...)

"""
    update_extKT!(diags::Vector{Graph}, para::DiagPara, legK::Vector{Vector{Float64}})

    Update the external momenta (`extK`) and external times (`extT`) of all the nodes in a vector of graphs in-place.

# Arguments
- `diags::Vector{Graph}`: A vector of `Graph` objects.
- `para::DiagPara`: parameters reconstructed in the graphs. Its `firstTauIdx` will update the `extT` of graphs.
- `legK::Vector{Vector{Float64}}`: basus of the external momenta for the legs of the diagram as [left in, left out, right in, right out]. 
"""
function update_extKT!(diags::Vector{Graph}, para::DiagPara, legK::Vector{Vector{Float64}})
    visited = Set{Int}()
    tauIdx = para.firstTauIdx
    # len_extK = para.totalLoopNum
    # num_extK = len_extK - para.innerLoopNum
    # extK = [k[1:len_extK] for k in legK[1:num_extK]]
    len_extK = length(legK[1])
    num_extK = length(legK) - 1
    extK = legK[1:end-1]

    sumK = zeros(len_extK)
    _K = zeros(len_extK)

    for graph in diags
        tau_shift = tauIdx - graph.properties.extT[1]
        for node in PreOrderDFS(graph)
            node.id in visited && continue
            node.id = IR.uid()
            push!(visited, node.id)
            prop = IR.properties(node)
            if !hasproperty(prop, :extK) || !hasproperty(prop, :extT)
                continue
            end
            K = prop.extK
            T = prop.extT
            if prop isa Ver4Id || prop isa Ver3Id
                for i in eachindex(K)
                    resize!(K[i], len_extK)
                    K[i] .= legK[i][1:len_extK]
                end
                if tau_shift != 0
                    node.properties = FrontEnds.reconstruct(prop, :para => para, :extT => Tuple(t + tau_shift for t in T))
                else
                    node.properties = FrontEnds.reconstruct(prop, :para => para)
                end
            elseif prop isa PropagatorId || prop isa GreenId || prop isa SigmaId || prop isa PolarId
                original_len_K = length(K)
                if length(K) < len_extK
                    resize!(K, len_extK)
                    K[original_len_K+1:end] .= 0.0
                else
                    resize!(K, len_extK)
                end

                _K[num_extK+1:end] .= K[num_extK+1:end]
                for i in eachindex(extK)
                    sumK .+= K[i] * extK[i]
                end
                K .= sumK .+ _K
                fill!(sumK, 0.0)
                if tau_shift != 0
                    node.properties = FrontEnds.reconstruct(prop, :extT => Tuple(t + tau_shift for t in T))
                end
            end
        end
    end
end

"""
    update_extKT(diags::Vector{Graph}, para::DiagPara, legK::Vector{Vector{Float64}}) -> Vector{Graph}

    Returns a new vector of graphs with updated external momenta (`extK`) and external times (`extT`), 
    based on the provided graphs, parameters, and external legs' momenta.

# Arguments
- `diags::Vector{Graph}`: A vector of `Graph` objects.
- `para::DiagPara`: parameters reconstructed in the graphs. Its `firstTauIdx` will update the `extT` of graphs.
- `legK::Vector{Vector{Float64}}`: basus of the external momenta for the legs of the diagram as [left in, left out, right in, right out]. 

# Returns
- `Vector{Graph}`: A new vector of `Graph` objects with updated `extK`, `extT`, and `para` (if existed) properties for each node.
"""
function update_extKT(diags::Vector{Graph}, para::DiagPara, legK::Vector{Vector{Float64}})
    graphs = deepcopy(diags)
    update_extKT!(graphs, para, legK)
    return graphs
end