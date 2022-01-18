# Base.hash(d::DiagramId) = hash(d) % 1000000

struct Diagram{T<:DiagramId,W}
    hash::Int # Two diagram MAY have the same hash number, only provide a inttuiative way to distinish two diagrams. 
    id::T
    operator::Operator
    factor::W
    subdiagram::Vector{Any}

    weight::W
    # parent::Diagram

    function Diagram(id::T, operator::Operator = Add(), factor::W = 1.0) where {T<:DiagramId,W}
        return new{T,W}(hash(id) % 1000000, id, operator, factor, [], zero(W))
    end
end

function add_subdiagram!(parent::Diagram, child::Diagram)
    for c in parent.subdiagram
        if c.id == child.id
            return false
        end
    end
    push!(parent.subdiagram, child)
end

function toDataFrame(diag::Diagram, maxdepth::Int = 1)
    @assert maxdepth == 1 "deep convert has not yet been implemented!"
    d = Dict{Symbol,Any}(Dict(diag.id))
    # d = id2Dict(diag.id)
    d[:hash] = diag.hash
    d[:operator] = diag.operator
    d[:factor] = diag.factor
    d[:weight] = diag.weight
    d[:DiagramId] = diag.id
    d[:subdiagram] = Tuple(d.id for d in diag.subdiagram)
    return DataFrame(d)
end

function toDataFrame(diagVec::AbstractVector, maxdepth::Int = 1)
    diags = []
    for diag in diagVec
        push!(diags, toDataFrame(diag, maxdepth))
    end
    return vcat(diags)
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