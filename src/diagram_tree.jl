module DiagTree
include("pool.jl")
using ..Var


struct Propagator{PARA,F}
    para::PARA
    order::Int
    factor::F
    variable::Vector{Int}
    function Propagator(order, variable = [], factor::F = 1.0, para::P = 0) where {F,P}
        return new{P,F}(para, order, factor, variable)
    end
    function Propagator{P,F}(order, variable = [], factor = 1.0, para = 0) where {P,F}
        return new{P,F}(para, order, factor, variable)
    end
end

function Base.isequal(a::Propagator{P,F}, b::Propagator{P,F}) where {P,F}
    if (isequal(a.para, b.para) == false) || (a.order != b.order) || (a.variable != b.variable)
        return false
    else
        return true
    end
end
Base.:(==)(a::Propagator{P,F}, b::Propagator{P,F}) where {P,F} = Base.isequal(a, b)

struct Node{PARA,F}
    para::PARA
    operation::Int #1: multiply, 2: add, ...
    factor::F
    components::Vector{Vector{Int}}
    child::Vector{Int}
    parent::Int # parent id
    # child::SubArray{CachedObject{Node{PARA, P},W},1,Vector{CachedObject{NODE,W}},Tuple{Vector{Int64}},false}

    function Node(operation::Int, components = [[]], child = [], factor::F = 1.0, parent = 0, para::P = 0) where {F,P}
        return new{P,F}(para, operation, factor, components, child, parent)
    end
    function Node{P,F}(operation::Int, components = [[]], child = [], factor = 1.0, parent = 0, para = 0) where {F,P}
        return new{P,F}(para, operation, factor, components, child, parent)
    end
end

function Base.isequal(a::Node{P}, b::Node{P}) where {P}
    # only parent is allowed to be different
    if (isequal(a.para, b.para) == false) || (a.operation != b.operation) || (a.components != b.components) || (a.child != b.child) || (a.factor != b.factor)
        return false
    else
        return true
    end
end
Base.:(==)(a::Node{P}, b::Node{P}) where {P} = Base.isequal(a, b)

mutable struct Diagrams{V,P,PARA,F,W}
    basisPool::V
    propagatorPool::P
    nodePool::Pool{Node{PARA,F},W}
    root::Vector{Int}
    # root::SubArray{CachedObject{NODE,W},1,Vector{CachedObject{NODE,W}},Tuple{Vector{Int64}},false}
    #SubArray has 5 type parameters. The first two are the standard element type and dimensionality. 
    #The next is the type of the parent AbstractArray. The most heavily-used is the fourth parameter, 
    #a Tuple of the types of the indices for each dimension. The final one, L, is only provided as a convenience for dispatch;
    #it's a boolean that represents whether the index types support fast linear indexing.
    # loopNum::Int
    # tauNum::Int
    function Diagrams{PARA,F,W}(basis::V, propagator::P) where {V,P,PARA,F,W}
        return new{V,P,PARA,F,W}(basis, propagator, Pool{Node{PARA,F},W}(), [])
    end
    function Diagrams{F,W}(basis::V, propagator::P) where {V,P,F,W}
        PARA = Int
        return new{V,P,PARA,F,W}(basis, propagator, Pool{Node{PARA,F},W}(), [])
    end
end

function addPropagator(diag::Diagrams, index::Int, order::Int, basis::AbstractVector, factor = 1, para = 0, currWeight = 0)
    basisPool = diag.basisPool
    propagatorPool = diag.propagatorPool
    # @assert length(basis) == length(variablePool) == length(currVar) "$(length(basis)) == $(length(variablePool)) == $(length(currVar)) breaks"

    PROPAGATOR_POOL = typeof(propagatorPool[index])
    CACHEDPROPAGATOR = eltype(fieldtype(PROPAGATOR_POOL, :pool))
    PROPAGATOR = fieldtype(CACHEDPROPAGATOR, :object)
    PARA = fieldtype(PROPAGATOR, :para)
    F = fieldtype(PROPAGATOR, :factor)

    vidx = zeros(length(basisPool))
    for (bi, b) in enumerate(basis)
        # b[1]: basis, b[2]: initialize variable (curr)
        vidx[bi] = append(basisPool[bi], b[1], b[2])
    end
    prop = Propagator{PARA,F}(order, vidx, factor, para)
    return append(diag.propagatorPool[index], prop, currWeight)
end

# function Node{F, P}(operation::Int, components = [[]], child = [], parent = 0, factor = 1.0, para = 0) where {F,P}
# function addNode(diag::Diagrams, operator, components::Vector{AbstractVector}, childNodes::AbstractVector; factor = 1.0, parent = 0, para = 0, currWeight = 0)
function addNode(diag::Diagrams, operator, components, childNodes; factor = 1.0, parent = 0, para = 0, currWeight = 0.0)
    nodePool = diag.nodePool
    @assert length(components) == length(diag.propagatorPool) "each element of the components is an index vector of the corresponding propagator"

    _NodePool = typeof(nodePool)
    _CachedNode = eltype(fieldtype(_NodePool, :pool))
    _Node = fieldtype(_CachedNode, :object)
    PARA = fieldtype(_Node, :para)
    F = fieldtype(_Node, :factor)

    node = Node{PARA,F}(operator, components, childNodes, factor, parent, para)

    # println("para: ", PARA)
    # println("F: ", F)
    # println("node: ", node)
    nidx = append(nodePool, node, currWeight)
    return nidx
end

end