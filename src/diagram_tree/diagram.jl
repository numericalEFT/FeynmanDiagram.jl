abstract type DiagramId end

# toDict(d::DiagramId) = error("toDict not implemented!")
Base.show(io::IO, d::DiagramId) = error("Base.show not implemented!")
Base.isequal(a::DiagramId, b::DiagramId) = error("Base.isequal not implemented!")
Base.:(==)(a::DiagramId, b::DiagramId) = Base.isequal(a, b)

# Base.hash(d::DiagramId) = hash(d) % 1000000

struct Diagram{T<:DiagramId}
    hash::Int
    id::T
    subdiagram::Vector{Diagram}
    # parent::Diagram

    Diagram(id::T) where {T<:DiagramId} = new{T}(hash(id) % 1000000, id, [])
end

function add_subdiagram!(parent::Diagram, child::Diagram)
    for c in parent.subdiagram
        if c.id == child.id
            return false
        end
    end
    push!(parent.subdiagram, child)
end

function toDataFrame(diagVec::AbstractVector,
    id2Dict::Function = x -> Dict{Symbol,Any}([fn => getfield(x, fn) for fn ∈ fieldnames(typeof(x))]))
    diags = []
    for diag in diagVec
        d = id2Dict(diag.id)
        d[:hash] = diag.hash
        d[:DiagramId] = diag.id
        d[:subdiagram] = [d.id for d in diag.subdiagram]
        push!(diags, d)
    end
    return DataFrame(diags)
end

## Things we need to define
function AbstractTrees.children(diag::Diagram)
    return diag.subdiagram
end

## Things that make printing prettier
AbstractTrees.printnode(io::IO, diag::Diagram) = print(io, "\u001b[32m$(diag.hash)\u001b[0m : $(diag.id)")

## Optional enhancements
# These next two definitions allow inference of the item type in iteration.
# (They are not sufficient to solve all internal inference issues, however.)
# Base.eltype(::Type{<:TreeIterator{BinaryNode{T}}}) where {T} = BinaryNode{T}
# Base.IteratorEltype(::Type{<:TreeIterator{BinaryNode{T}}}) where {T} = Base.HasEltype()

## Let's test it. First build a tree.