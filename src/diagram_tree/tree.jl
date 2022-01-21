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

_diagram(df, index) = df[index, :Diagram]

function mergeby(diags::Vector{Diagram{W}}, fields; verbose::Int = 0, operator = Sum(), factor = one(W),
    getid::Function = g -> GenericId(_diagram(g, 1).id.para, g[1, fields]),
    getname::Function = g -> :none
) where {W}

    df = toDataFrame(diags, verbose = verbose)
    if all(x -> x == df[1, :id], df[!, :id]) == false
        @warn "Not all DiagramIds in $diags are the same!"
    end

    group = DataFrames.groupby(df, fields)

    d = Dict{Any,Diagram{W}}()

    for g in group
        key = g[1, fields]
        if length(key) > 1
            key = Tuple(key)
        end
        d[key] = Diagram(getid(g), operator, g[:, :Diagram], name = getname(g), factor = factor)
    end
    return d
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