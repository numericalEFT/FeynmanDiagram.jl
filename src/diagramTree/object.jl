"""
    mutable struct Propagator{PARA,F}

     Object of all kinds of propagators (many-body Green's functions, interactions, vertex functions, ...)

# Members
- para::PARA : User-defined parameters, which will be used to evaluate the factor and the weight of the propagator
- order::Int : The order of the propagators (for example, in a Feynman diagrammatic expansion, the one-body Green's function usually has order zero, and an interaction usually has order one)
- factor::F : Additional factor of the propagator
- basis::Vector{Int} : Index to the cached basis stored in certain pool. They are essentail for the weight evaluation of the propagator (the first index for the first basis type, the second index for the second basis type, etc.). 
"""
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

function propagatorPool(weightType::DataType, factorType::DataType, paraType::DataType = Int)
    propagatorType = Propagator{paraType,factorType}
    return Pool{Cache{propagatorType,weightType}}()
end

"""
    mutable struct Node{PARA,F}

    Node Object, which is the building block of the diagram tree. Each node is a collection of CACHED proapgator objects and other child CACHED node objects

# Members
- para::PARA     : user-defined parameters, which will be used to evaluate the factor and the weight of the node (e.g., if the node represents a vertex function, then the parameter may be the momentum basis of the external legs)
- operation::Int : #1: multiply, 2: add, ...
- factor::F      : additional factor of the node
- components::Vector{Vector{Int}}  : Index to the cached propagators stored in certain pools. Each Vector{Int} is for one kind of propagator.
- childNodes::Vector{Int}  : Indices to the cached nodes stored in certain pool. They are the child of the current node in the diagram tree.
- parent::Int : Index to the cached nodes which is the parent of the current node.
"""
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

"""
    mutable struct Diagrams{V,P,PARA,F,W}

    Diagram Object represents a set of Feynman diagrams in a tree (forest) structure

# Members
- basisPool::V      : Tuple of pools of cached basis  in a format of (BasisPool1, BasisPool2, ...)
- propagatorPool::P : Tuple of pools of cached propagators in a format of (PropagatorPool1, PropagatorPool2, ...)
- nodePool::Pool{Node{PARA,F},W} : Pool of the nodes in the diagram tree
- root::Vector{Int} : indices of the cached nodes that are the root(s) of the diagram tree. Each element corresponds to one root.
"""
mutable struct Diagrams{V,P,PARA,F,W}
    basisPool::V
    propagatorPool::P
    nodePool::Pool{Cache{Node{PARA,F},W}}
    root::Vector{Int}
    # function Diagrams{PARA,F,W}(basisPool::V, propagatorPool::P) where {V,P,PARA,F,W}
    #     return new{V,P,PARA,F,W}(basisPool, propagatorPool, Pool{Cache{Node{PARA,F},W}}(), [])
    # end
    # function Diagrams{F,W}(basisPool::V, propagatorPool::P) where {V,P,F,W}
    #     PARA = Int
    #     return new{V,P,PARA,F,W}(basisPool, propagatorPool, Pool{Cache{Node{PARA,F},W}}(), [])
    # end
    function Diagrams(basisPool::V, propagatorPool::P, nodeWeightType::DataType, nodeFactorType::DataType, nodeParaType::DataType = Int) where {V,P}
        @assert V <: Tuple "Tuple is required for efficiency!"
        @assert P <: Tuple "Tuple is required for efficiency!"
        # println(basisPool)
        # println(propagatorPool)
        nodePool = Pool{Cache{Node{nodeParaType,nodeFactorType},nodeWeightType}}()
        return new{V,P,nodeParaType,nodeFactorType,nodeWeightType}(basisPool, propagatorPool, nodePool, [])
    end
end

"""
    function addPropagator(diag::Diagrams, index::Int, order::Int, basis::AbstractVector, factor = 1, para = 0, currWeight = 0)

    Add a propagator into the diagram tree.

# Arguments
- diag::Diagrams : Diagram tree.
- index::Int     : Index of the propagator in the diagram tree.
- order::Int     : Order of the propagator
- basis::AbstractVector : Vector of the pair [variable basis, initial variable]. For example, if a propagator involves both a momentum and a time variable, then the basis should be [[Kbasis, K0], [Tbasis, T0]], where K0, T0 are the initial values of the momentum and the time.
- factor = 1     : Factor of the propagator
- para = 0       : Additional paramenter required to evaluate the propagator. If not needed, simply leave it as an integer.
- currWeight = 0 : Initial weight of the propagator
"""
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
        if isCached(basisPool[bi])
            vidx[bi] = append(basisPool[bi], b[1], b[2], true)
        else
            vidx[bi] = append(basisPool[bi], b[1], b[2], false)
        end
    end
    prop = Propagator{PARA,F}(order, vidx, factor, para)
    return append(diag.propagatorPool[index], prop, currWeight, true)
end

"""
    function addNode(diag::Diagrams, operator, components, childNodes; factor = 1.0, parent = 0, para = 0, currWeight = 0.0)

    Add a node into the diagram tree.

# Arguments
- diag::Diagrams : Diagram tree.
- operator::Int  : #1: multiply, 2: add, ...
- components     : Index to the cached propagators stored in certain pools. It should be in the format of Vector{Vector{Int}}, where each Vector{Int} is for one kind of propagator.
- childNodes     : Indices to the cached nodes stored in certain pool. They are the child of the current node in the diagram tree. It should be in the format of Vector{Int}.
- factor = 1     : Factor of the node
- para = 0       : Additional paramenter required to evaluate the node. If not needed, simply leave it as an integer.
- currWeight = 0 : Initial weight of the node.
"""
function addNode(diag::Diagrams, operator, components, childNodes; factor = 1.0, parent = 0, para = 0, currWeight = 0.0)
    nodePool = diag.nodePool
    @assert length(components) == length(diag.propagatorPool) "each element of the components is an index vector of the corresponding propagator"

    _NodePool = typeof(nodePool)
    _CachedNode = eltype(fieldtype(_NodePool, :pool))
    _Node = fieldtype(_CachedNode, :object)
    PARA = fieldtype(_Node, :para)
    F = fieldtype(_Node, :factor)

    node = Node{PARA,F}(operator, components, childNodes, factor, parent, para)

    nidx = append(nodePool, node, currWeight, true)
    return nidx
end