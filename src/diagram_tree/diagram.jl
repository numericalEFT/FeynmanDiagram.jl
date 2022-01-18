# Base.hash(d::DiagramId) = hash(d) % 1000000

struct Diagram{W}
    hash::Int # Two diagram MAY have the same hash number, only provide a inttuiative way to distinish two diagrams. 
    id::DiagramId
    operator::Operator
    factor::W
    subdiagram::Vector{Diagram{W}}

    weight::W
    # parent::Diagram

    function Diagram(id::DiagramId, operator::Operator = Sum(), factor::W = 1.0) where {W}
        return new{W}(hash(id) % 1000000, id, operator, factor, [], zero(W))
    end
end

function Base.show(io::IO, diag::Diagram)

    if length(diag.subdiagram) == 0
        typestr = "bare"
    elseif length(diag.subdiagram) == 1
        typestr = "#$(diag.hash)"
    else
        subdiag = prod(["#$(d.hash), " for d in diag.subdiagram[1:end-1]])
        subdiag *= "#$(diag.subdiagram[end].hash)"
        typestr = "$(diag.operator)($subdiag)"
    end
    print(io, "\u001b[32m#$(diag.hash)\u001b[0m : $(diag.id) = $(diag.factor)x$typestr = $(diag.weight)")
end

isbare(diag::Diagram) = isempty(diag.subdiagram)

function add_subdiagram!(parent::Diagram, child::Diagram)
    for c in parent.subdiagram
        if c.id == child.id
            return false
        end
    end
    push!(parent.subdiagram, child)
end


@inline apply(o::Sum, diags::Vector{Diagram{W}}) where {W<:Number} = sum(d.weight for d in diags)
@inline apply(o::Prod, diags::Vector{Diagram{W}}) where {W<:Number} = prod(d.weight for d in diags)
@inline apply(o::Sum, diag::Diagram{W}) where {W<:Number} = diag.weight
@inline apply(o::Prod, diag::Diagram{W}) where {W<:Number} = diag.weight

function eval!(diag::Diagram)
    if isbare(diag)
        diag.weight = apply(diag.operator, diag.subdiagram) * factor
    else
        diag.weight = eval(diag.id) * factor
    end
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

function toDataFrame(diagVec::Vector{Diagram{W}}, maxdepth::Int = 1) where {W}
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
# AbstractTrees.printnode(io::IO, diag::Diagram) = print(io, "\u001b[32m$(diag.hash)\u001b[0m : $(diag.id)")
AbstractTrees.printnode(io::IO, diag::Diagram) = print(io, "$(diag)")

## Optional enhancements
# These next two definitions allow inference of the item type in iteration.
# (They are not sufficient to solve all internal inference issues, however.)
# Base.eltype(::Type{<:TreeIterator{BinaryNode{T}}}) where {T} = BinaryNode{T}
# Base.IteratorEltype(::Type{<:TreeIterator{BinaryNode{T}}}) where {T} = Base.HasEltype()

## Let's test it. First build a tree.