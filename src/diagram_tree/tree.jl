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
    function Diagram{W}(id::DiagramId, operator::Operator = Sum(), subdiagram = [];
        name = :none, factor = W(1), weight = W(0)) where {W}
        return new{W}(uid(), name, id, operator, factor, deepcopy(subdiagram), weight)
    end
    function Diagram(id::DiagramId, operator::Operator = Sum(), subdiagram = []; type::DataType = id.para.weightType,
        name = :none, factor = one(type), weight = zero(type))
        return new{type}(uid(), name, id, operator, factor, deepcopy(subdiagram), weight)
    end
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

function mergeby(df::DataFrame, fields = [];
    operator = Sum(), name::Symbol = :none,
    factor = isempty(df) ? nothing : one(df[1, :diagram].factor),
    getid::Function = g -> isempty(g) ? nothing : GenericId(g[1, :diagram].id.para, Tuple(g[1, fields]))
)
    if isempty(df)
        return df
    end

    # df = toDataFrame(diags, verbose = verbose)
    if all(x -> typeof(x.id) == typeof(df[1, :diagram].id), df[!, :diagram]) == false
        @warn "Not all DiagramIds in $diags are the same!"
    end

    groups = DataFrames.groupby(df, fields, sort = true)

    # combine diagrams in a group into one composite diagram
    gdf = combine(groups) do group # for each group in groups
        # check the documentation of ``combine" for details https://dataframes.juliadata.org/stable/man/split_apply_combine/

        if nrow(group) == 1
            # if there is only one diagram in df, and the new id is either GenericId or the id of the existing diagram, 
            # then simply return the current df without creating a new diagram
            # ! the new factor will be multiplied to the factor of the exisiting diagram!
            id = getid(group)
            if id isa GenericId || typeof(id) == typeof(group[1, :diagram].id)
                diag = deepcopy(group[1, :diagram])
                diag.factor *= factor
                return (diagram = diag, hash = diag.hash)
            end
        end
        diag = Diagram(getid(group), operator, group[:, :diagram], name = name, factor = factor)
        (diagram = diag, hash = diag.hash)
    end
    return gdf
end

function mergeby(diags::AbstractVector, fields = []; expand::Bool = false, kwargs...)
    df = toDataFrame(diags, expand = expand)
    return mergeby(df, fields; kwargs...)
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

## Optional enhancements
# These next two definitions allow inference of the item type in iteration.
# (They are not sufficient to solve all internal inference issues, however.)
# Base.eltype(::Type{<:TreeIterator{BinaryNode{T}}}) where {T} = BinaryNode{T}
# Base.IteratorEltype(::Type{<:TreeIterator{BinaryNode{T}}}) where {T} = Base.HasEltype()
