abstract type DiagramId end

# toDict(d::DiagramId) = error("toDict not implemented!")
Base.Dict(x::DiagramId) = Dict{Symbol,Any}([fn => getfield(x, fn) for fn ∈ fieldnames(typeof(x))])
Base.show(io::IO, d::DiagramId) = error("Base.show not implemented!")
Base.isequal(a::DiagramId, b::DiagramId) = error("Base.isequal not implemented!")
Base.:(==)(a::DiagramId, b::DiagramId) = Base.isequal(a, b)
eval(d::DiagramId) = error("eval for $d has not yet implemented!")

# Base.hash(d::DiagramId) = hash(d) % 1000000

struct Diagram{T<:DiagramId}
    hash::Int # Two diagram MAY have the same hash number, only provide a inttuiative way to distinish two diagrams. 
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

function toDataFrame(diagVec::AbstractVector)
    diags = []
    for diag in diagVec
        d = Dict{Symbol,Any}(Dict(diag.id))
        # d = id2Dict(diag.id)
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