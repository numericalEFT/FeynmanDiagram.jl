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


function update_extK!(diag::Graph, extK::Vector{Vector{Float64}})
    visited = Set{Int}()
    num_extK = length(extK)
    len_extK = length(extK[1])

    sumK = zeros(len_extK)
    _K = zeros(len_extK)
    for leaf in Leaves(diag)
        if !(leaf.id in visited)
            push!(visited, leaf.id)
            prop = IR.properties(leaf)
            K = prop.extK

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
            # K[1:end] = sum([K[i] * extK[i] for i in eachindex(extK)]) + _K

            # if prop isa BareGreenId
            #     new_properties = BareGreenId(prop.type, k=K, t=prop.extT)
            # elseif prop isa BareInteractionId
            #     new_properties = BareInteractionId(prop.response, prop.type, k=K, t=prop.extT)
            # else
            #     error("unexpected property type $prop")
            # end
            # IR.set_properties!(leaf, new_properties)
        end
    end
end

function update_extK!(diags::Vector{Graph}, extK::Vector{Vector{Float64}})
    visited = Set{Int}()
    num_extK = length(extK)
    len_extK = length(extK[1])

    sumK = zeros(len_extK)
    _K = zeros(len_extK)
    for diag in diags
        for leaf in Leaves(diag)
            if !(leaf.id in visited)
                push!(visited, leaf.id)
                prop = IR.properties(leaf)
                K = prop.extK

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
            end
        end
    end
end