# Base.hash(d::DiagramId) = hash(d) % 1000000

mutable struct Diagram{W}
    hash::Int
    name::Symbol
    id::DiagramId
    operator::Operator
    factor::W
    subdiagram::Vector{Diagram{W}}

    weight::W
    # parent::Diagram

    function Diagram(id::DiagramId, operator::Operator = Sum(), subdiagram = []; type::DataType = id.para.weightType,
        name = :none, factor = one(type), weight = zero(type))
        return new{type}(uid(), name, id, operator, factor, subdiagram, weight)
    end

    #constructor for DiagramId without a field of GenericPara
    function Diagram{W}(id::DiagramId, operator::Operator = Sum(), subdiagram = [];
        name = :none, factor = W(1), weight = W(0)) where {W}
        return new{W}(uid(), name, id, operator, factor, subdiagram, weight)
    end
end

isbare(diag::Diagram) = isempty(diag.subdiagram)

function addSubDiagram!(parent::Diagram, child::Diagram)
    for c in parent.subdiagram
        if c.id == child.id
            return false
        end
    end
    push!(parent.subdiagram, child)
end

function addSubDiagram!(parent::Diagram, child::Vector{Diagram{W}}) where {W}
    for d in child
        addSubDiagram!(parent, d)
    end
end

# _diagram(df, index) = df[index, :Diagram]

function mergeby(df::DataFrame, fields = [];
    operator = Sum(), name::Symbol = :none, factor = one(df[1, :diagram].factor),
    getid::Function = g -> GenericId(g[1, :diagram].id.para, g[1, fields])
)

    # df = toDataFrame(diags, verbose = verbose)
    if all(x -> typeof(x.id) == typeof(df[1, :diagram].id), df[!, :diagram]) == false
        @warn "Not all DiagramIds in $diags are the same!"
    end

    groups = DataFrames.groupby(df, fields, sort = true)

    # combine diagrams in a group into one composite diagram
    gdf = combine(groups) do group # for each group in groups
        # check the documentation of ``combine" for details https://dataframes.juliadata.org/stable/man/split_apply_combine/
        (diagram = Diagram(getid(group), operator, group[:, :diagram], name = name, factor = factor),)
    end
    return gdf
end

function mergeby(diags::Vector{Diagram{W}}, fields = []; expand::Bool = false, kwargs...) where {W}
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
