# Base.hash(d::DiagramId) = hash(d) % 1000000
"""
    mutable struct Diagram{W}
    
    struct of a diagram. A diagram of a sum or produce of various subdiagrams.

# Members
- hash::Int          : the unique hash number to identify the diagram
- name::Symbol       : name of the diagram
- id::DiagramId      : diagram id 
- operator::Operator : operation, support Sum() and Prod()
- factor::W          : additional factor of the diagram
- subdiagram::Vector{Diagram{W}}   : vector of sub-diagrams
- weight::W          : weight of the diagram
"""
mutable struct Diagram{W}
    hash::Int
    name::Symbol
    id::DiagramId
    operator::Operator
    factor::W
    subdiagram::Vector{Diagram{W}}

    weight::W
    # parent::Diagram

    """
        function Diagram{W}(id::DiagramId, operator::Operator = Sum(), subdiagram = []; name = :none, factor = W(1), weight = W(0)) where {W}

        construct Diagram struct.

    # Arguments
    - id            : DiagramId of the diagram
    - operator      : Sum() or Prod()
    - subdiagram    : subdiagrams of the diagram, it should be a vector of Diagram struct
    - name          : name of the diaram
    - factor        : the additional factor of the diagram
    - weight        : the initial weight
    """
    function Diagram{W}(id::DiagramId, operator::Operator=Sum(), subdiagram=[];
        name=:none, factor=W(1), weight=W(0)) where {W}
        # return new{W}(uid(), name, id, operator, factor, deepcopy(subdiagram), weight)
        return new{W}(uid(), name, id, operator, factor, subdiagram, weight)
    end
    # function Diagram(id::DiagramId, operator::Operator=Sum(), subdiagram=[]; type::DataType=id.para.weightType,
    #     name=:none, factor=one(type), weight=zero(type))
    #     # return new{type}(uid(), name, id, operator, factor, deepcopy(subdiagram), weight)
    #     return new{type}(uid(), name, id, operator, factor, subdiagram, weight)
    # end
end

isbare(diag::Diagram) = isempty(diag.subdiagram)

# function addSubDiagram!(parent::Diagram, child::Diagram)
#     for c in parent.subdiagram
#         if c.id == child.id
#             return false
#         end
#     end
#     push!(parent.subdiagram, deepcopy(child))
# end

# function addSubDiagram!(parent::Diagram, child::Vector{Diagram{W}}) where {W}
#     for d in child
#         addSubDiagram!(parent, d)
#     end
# end

# _diagram(df, index) = df[index, :Diagram]

function _combinegroups(groups, getid, factor, operator, name)
    # combine diagrams in a group into one composite diagram
    gdf = combine(groups) do group # for each group in groups
        # check the documentation of ``combine" for details https://dataframes.juliadata.org/stable/man/split_apply_combine/
        # id = isnothing(getid) ? GenericId(group.diagram[1].id.para, Tuple(group[1, fields])) : getid(group)
        id = getid(group)

        if nrow(group) == 1
            # if there is only one diagram in df, and the new id is either GenericId or the id of the existing diagram, 
            # then simply return the current df without creating a new diagram
            # ! the new factor will be multiplied to the factor of the exisiting diagram!
            if id isa GenericId || typeof(id) == typeof(group.diagram[1].id)
                # diag = deepcopy(group[1, :diagram])
                diag = group.diagram[1]
                diag.factor *= factor
                return (diagram=diag, hash=diag.hash)
            end
        end
        W = typeof(group.diagram[1].weight)
        diag = Diagram{W}(id, operator, group.diagram, name=name, factor=factor)
        return (diagram=diag, hash=diag.hash)
    end
    return gdf
end

function _mergediag(::Type{W}, group, factor, id, operator, name) where {W}
    if nrow(group) == 1
        # if there is only one diagram in df, and the new id is either GenericId or the id of the existing diagram, 
        # then simply return the current df without creating a new diagram
        # ! the new factor will be multiplied to the factor of the exisiting diagram!
        if id isa GenericId || typeof(id) == typeof(group.diagram[1].id)
            # diag = deepcopy(group[1, :diagram])
            diag::Diagram{W} = group.diagram[1]
            diag.factor *= factor
            return diag
        end
    end
    return Diagram{W}(id, operator, group.diagram, name=name, factor=factor)
end

function _combine(::Type{W}, groups, factor, getid, operator, name) where {W}
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
    d[:diagram] = [_mergediag(W, groups[key], factor, getid(groups[key]), operator, name) for key in _keys]
    d[:hash] = [diag.hash for diag in d[:diagram]]
    return DataFrame(d, copycols=false)
end

function mergeby(df::DataFrame, fields=Vector{Symbol}();
    operator=Sum(), name::Symbol=:none, factor=1.0,
    getid::Function=g -> GenericId(g[1, :diagram].id.para, Tuple(g[1, fields]))
)
    if isempty(df)
        return df
    else
        W = typeof(df.diagram[1].weight)
        return mergeby(W, df, fields; operator=operator, name=name, factor=factor, getid=getid)
    end
end

function mergeby(::Type{W}, df::DataFrame, fields=Vector{Symbol}();
    operator=Sum(), name::Symbol=:none, factor=1.0,
    getid::Function=g -> GenericId(g[1, :diagram].id.para, Tuple(g[1, fields]))
) where {W}
    if isempty(df)
        return df
    else
        if all(x -> typeof(x.id) == typeof(df.diagram[1].id), df[!, :diagram]) == false
            @warn "Not all DiagramIds in $df are the same!"
        end
        groups = DataFrames.groupby(df, fields, sort=true)
        ########  less memory usage but can not pass the test right now ##############
        d = _combine(W, groups, factor, getid, operator, name)
        ######## alternative approach (more memory)  ##################
        # d = _combinegroups(groups, getid, factor, operator, name)
        # println("old\n$d \n new\n$cd")
        return d
    end
end

function mergeby(diags::Union{Diagram,Tuple,AbstractVector}, fields=nothing; idkey=nothing, kwargs...)
    if diags isa Diagram
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
function mergeby(diags::Vector{Diagram{W}}, fields; idkey=Vector{Symbol}(), kwargs...) where {W}
    if isempty(diags)
        return diags
    else
        df = toDataFrame(diags, idkey)
        mergedf = mergeby(df, fields; kwargs...)
        return Vector{Diagram{W}}(mergedf.diagram)
    end
end

function mergeby(diags::Vector{Diagram{W}};
    operator=Sum(), name::Symbol=:none, factor=1.0,
    getid::Function=d -> GenericId(d[1].id.para::DiagPara{W})
) where {W}
    if isempty(diags)
        return diags
    else
        id = getid(diags)
        if length(diags) == 1 && (id isa GenericId || typeof(id) == typeof(diags[1].id))
            # if there is only one diagram, and the new id is either GenericId or the id of the existing diagram, 
            # then simply return the current diagram without creating a new diagram
            # ! the new factor will be multiplied to the factor of the exisiting diagram!
            diags[1].factor *= factor
            return diags
        end
        diag = Diagram{W}(id, operator, diags, name=name, factor=factor)
        return [diag,]
    end
end
# mergeby(df::DataFrame; kwargs...) = mergeby(df, []; kwargs...)
# mergeby(diags::Vector{Diagram{W}}; kwargs...) where {W} = mergeby(diags, []; kwargs...)



#####################  interface to AbstractTrees ########################### 
function AbstractTrees.children(diag::Diagram)
    return diag.subdiagram
end

## Things that make printing prettier
AbstractTrees.printnode(io::IO, diag::Diagram) = print(io, "\u001b[32m$(diag.hash)\u001b[0m : $diag")
# AbstractTrees.printnode(io::IO, diag::Diagram) = print(io, "$(diag)")

######### define the following for the type stability #########################
# AbstractTrees.childrentype(diag::Diagram{W}) where {W} = Vector{Diagram{W}}

# AbstractTrees.NodeType(::Diagram{W}) where {W} = HasNodeType()
AbstractTrees.nodetype(::Diagram{W}) where {W} = Diagram{W}

## Optional enhancements
# These next two definitions allow inference of the item type in iteration.
# (They are not sufficient to solve all internal inference issues, however.)
Base.eltype(::Type{<:TreeIterator{Diagram{W}}}) where {W} = Diagram{W}
Base.IteratorEltype(::Type{<:TreeIterator{Diagram{W}}}) where {W} = Base.HasEltype()
