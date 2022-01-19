# Base.hash(d::DiagramId) = hash(d) % 1000000

mutable struct Diagram{W}
    hash::Int # Two diagram MAY have the same hash number, only provide a inttuiative way to distinish two diagrams. 
    id::DiagramId
    operator::Operator
    factor::W
    subdiagram::Vector{Diagram{W}}

    weight::W
    # parent::Diagram

    function Diagram{W}(operator::Operator = Sum(), subdiagram = []; id::DiagramId = GenericId(), factor = W(1), weight = W(0)) where {W}
        return new{W}(hash(id) % 1000000, id, operator, factor, subdiagram, weight)
    end
    function Diagram(operator::Operator = Sum(), subdiagram = []; id::DiagramId = GenericId(), factor = Float64(1), weight = Float64(0))
        return new{Float64}(hash(id) % 1000000, id, operator, factor, subdiagram, weight)
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


@inline apply(o::Sum, diags::Vector{Diagram{W}}) where {W<:Number} = sum(d.weight for d in diags)
@inline apply(o::Prod, diags::Vector{Diagram{W}}) where {W<:Number} = prod(d.weight for d in diags)
@inline apply(o::Sum, diag::Diagram{W}) where {W<:Number} = diag.weight
@inline apply(o::Prod, diag::Diagram{W}) where {W<:Number} = diag.weight

function evalDiagNode!(diag::Diagram, evalBare::Function, vargs...; kwargs...)
    if isbare(diag)
        diag.weight = evalBare(diag.id, vargs...; kwargs...) * diag.factor
    else
        diag.weight = apply(diag.operator, diag.subdiagram) * diag.factor
    end
    return diag.weight
end

function evalDiagTree!(diag::Diagram, evalBare::Function, vargs...; kwargs...)
    for d in PostOrderDFS(diag)
        evalDiagNode!(d, evalBare, vargs...; kwargs...)
    end
    return diag.weight
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