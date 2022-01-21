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

function mergeby(df::DataFrame, fields;
    verbose::Int = 0, operator = Sum(), factor = one(df.Diagram[1].factor), name::Symbol = :none,
    getid::Function = g -> GenericId(g.Diagram[1].id.para, fields),
    kwargs...)

    # df = toDataFrame(diags, verbose = verbose)
    if all(x -> typeof(x.id) == typeof(df.Diagram[1].id), df.Diagram) == false
        @warn "Not all DiagramIds in $diags are the same!"
    end

    groups = DataFrames.groupby(df, fields)

    # return combine(groups, [:id, :Diagram] =>
    #     ((id, diagrams) -> Diagram(getid(id[1]), operator, diagrams, name = name, factor = factor)) => :Diagram)
    gdf = combine(groups) do g
        (Diagram = Diagram(getid(g), operator, g[:, :Diagram], name = name, factor = factor),)
    end
    return gdf
end
function mergeby(diags::Vector{Diagram{W}}, fields; kwargs...) where {W}
    df = toDataFrame(diags, verbose = kwargs[:verbose])
    return mergeby(df, fields; kwargs...)
end



#####################  interface to AbstractTrees ########################### 
function AbstractTrees.children(diag::Diagram)
    return diag.subdiagram
end

## Things that make printing prettier
# AbstractTrees.printnode(io::IO, diag::Diagram) = print(io, "\u001b[32m$(diag.hash)\u001b[0m : $(diag.id)")
AbstractTrees.printnode(io::IO, diag::Diagram) = print(io, "$(diag)")

## Optional enhancements
# These next two definitions allow inference of the item type in iteration.
# (They are not sufficient to solve all internal inference issues, however.)
# Base.eltype(::Type{<:TreeIterator{BinaryNode{T}}}) where {T} = BinaryNode{T}
# Base.IteratorEltype(::Type{<:TreeIterator{BinaryNode{T}}}) where {T} = Base.HasEltype()
