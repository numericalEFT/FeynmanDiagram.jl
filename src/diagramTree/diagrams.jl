struct Propagator{PARA,F}
    para::PARA
    order::Int
    factor::F
    basis::Vector{Int}
    function Propagator(order, basis = [], factor::F = 1.0, para::P = 0) where {F,P}
        return new{P,F}(para, order, factor, basis)
    end
    function Propagator{P,F}(order, basis = [], factor = 1.0, para = 0) where {P,F}
        return new{P,F}(para, order, factor, basis)
    end
end

function Base.isequal(a::Propagator{P,F}, b::Propagator{P,F}) where {P,F}
    if (isequal(a.para, b.para) == false) || (a.order != b.order) || (a.basis != b.basis)
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
    childNodes::Vector{Int}
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
    function Diagrams{PARA,F,W}(basisPool::V, propagatorPool::P) where {V,P,PARA,F,W}
        return new{V,P,PARA,F,W}(basisPool, propagatorPool, Pool{Node{PARA,F},W}(), [])
    end
    function Diagrams{F,W}(basisPool::V, propagatorPool::P) where {V,P,F,W}
        PARA = Int
        return new{V,P,PARA,F,W}(basisPool, propagatorPool, Pool{Node{PARA,F},W}(), [])
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

    vidx = zeros(length(basis))
    for (bi, b) in enumerate(basis)
        # b[1]: basis, b[2]: initialize variable (curr)
        vidx[bi] = append(basisPool[bi], b[1], b[2])
    end
    prop = Propagator{PARA,F}(order, vidx, factor, para)
    return append(diag.propagatorPool[index], prop, currWeight)
end

function addNode(diag::Diagrams, operator, components, childNodes; factor = 1.0, parent = 0, para = 0, currWeight = 0.0)
    nodePool = diag.nodePool
    @assert length(components) == length(diag.propagatorPool) "each element of the components is an index vector of the corresponding propagator"

    _NodePool = typeof(nodePool)
    _CachedNode = eltype(fieldtype(_NodePool, :pool))
    _Node = fieldtype(_CachedNode, :object)
    PARA = fieldtype(_Node, :para)
    F = fieldtype(_Node, :factor)

    node = Node{PARA,F}(operator, components, childNodes, factor, parent, para)

    nidx = append(nodePool, node, currWeight)
    return nidx
end