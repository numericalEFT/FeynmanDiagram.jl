"""
    mutable struct Propagator{PARA, F, N}

     Object of all kinds of propagators (many-body Green's functions, interactions, vertex functions, ...)

# Members
- para::PARA : User-defined parameters, which will be used to evaluate the factor and the weight of the propagator
- order::Int : The order of the propagators (for example, in a Feynman diagrammatic expansion, the one-body Green's function usually has order zero, and an interaction usually has order one)
- factor::F : Additional factor of the propagator
- basis::Vector{Int} : Index to the cached basis stored in certain pool. They are essentail for the weight evaluation of the propagator (the first index for the first basis type, the second index for the second basis type, etc.). 
"""
struct Propagator{PARA,F}
    name::Symbol
    para::PARA
    order::Int
    factor::F
    loopIdx::Int
    siteBasis::Vector{Int}
    # function Propagator(order, basis = [], factor::F = 1.0, para::P = 0) where {F,P}
    #     return new{P,F}(para, order, factor, basis)
    # end
    function Propagator{P,F}(name::Symbol, order::Int, para, factor, loopidx::Int, sitebasis) where {P,F}
        return new{P,F}(name, para, order, F(factor), loopidx, sitebasis)
    end
end

function Base.isequal(a::Propagator{P,F}, b::Propagator{P,F}) where {P,F}
    if (isequal(a.para, b.para) == false) || (a.order != b.order) || (a.siteBasis != b.siteBasis) || (a.loopIdx != b.loopIdx) || (a.factor ≈ b.factor) == false
        return false
    else
        return true
    end
end
Base.:(==)(a::Propagator{P,F}, b::Propagator{P,F}) where {P,F} = Base.isequal(a, b)

function propagatorPool(poolName::Symbol, weightType::DataType; factorType::DataType = weightType, paraType::DataType = Nothing)
    propagatorType = Propagator{paraType,factorType}
    return CachedPool(poolName, propagatorType, weightType)
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
    name::Symbol
    para::PARA
    operation::Int #1: multiply, 2: add, ...
    factor::F
    # order::Int
    propagators::Vector{Int}
    childNodes::Vector{Int}
    parent::Int # parent id
    # child::SubArray{CachedObject{Node{PARA, P},W},1,Vector{CachedObject{NODE,W}},Tuple{Vector{Int64}},false}

    # function Node(operation::Int, components = [[]], child = [], factor::F = 1.0, parent = 0, para::P = :N) where {F,P}
    #     return new{P,F}(para, operation, factor, components, child, parent)
    # end
    function Node{F}(name::Symbol, operation::Int, para::P, propagators = [], child = [], factor = 1.0, parent = 0) where {F,P}
        # @assert typeof(para) == P
        return new{P,F}(name, para, operation, F(factor), propagators, child, parent)
    end
    function Node{P,F}(name::Symbol, operation::Int, para, propagators = [], child = [], factor = 1.0, parent = 0) where {F,P}
        # @assert typeof(para) == P
        return new{P,F}(name, para, operation, F(factor), propagators, child, parent)
    end
end

function Base.isequal(a::Node{P}, b::Node{P}) where {P}
    # only parent is allowed to be different
    if length(a.propagators) != length(b.propagators)
        return false
    end
    for i in 1:length(a.propagators)
        if Set(a.propagators[i]) != Set(b.propagators[i])
            return false
        end
    end

    if (isequal(a.para, b.para) == false) || (a.operation != b.operation) || (Set(a.childNodes) != Set(b.childNodes)) || (a.factor ≈ b.factor) == false
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
    name::Symbol
    basisPool::V
    propagatorPool::P
    nodePool::CachedPool{Node{PARA,F},W}
    root::Vector{Int}
    # function Diagrams{PARA,F,W}(basisPool::V, propagatorPool::P) where {V,P,PARA,F,W}
    #     return new{V,P,PARA,F,W}(basisPool, propagatorPool, Pool{Cache{Node{PARA,F},W}}(), [])
    # end
    # function Diagrams{F,W}(basisPool::V, propagatorPool::P) where {V,P,F,W}
    #     PARA = Int
    #     return new{V,P,PARA,F,W}(basisPool, propagatorPool, Pool{Cache{Node{PARA,F},W}}(), [])
    # end
    function Diagrams(basisPool::V, propagatorPool::P, nodeWeightType::DataType; nodeFactorType = nodeWeightType, nodeParaType::DataType = Nothing, name = :none) where {V,P}
        nodePool = CachedPool(:node, Node{nodeParaType,nodeFactorType}, nodeWeightType)
        return new{V,P,nodeParaType,nodeFactorType,nodeWeightType}(name, basisPool, propagatorPool, nodePool, [])
    end
end

Base.getindex(diag::Diagrams, i) = diag.nodePool.current[diag.root[i]]
Base.firstindex(diag::Diagrams) = 1
Base.lastindex(diag::Diagrams) = length(diag.root)

struct Component
    index::Int
    isNode::Bool
    poolName::Symbol
    object::Any
end

Base.show(io::IO, c::Component) = print(io, "#$(c.index) in $(c.poolName) Pool")

Base.zero(::Type{Component}) = Component(0, false, :none, nothing)

function setroot!(diag, rootVec::Vector{Component})
    for r in rootVec
        @assert r.isNode "root of Diagrams must be a node!"
    end
    diag.root = [r.index for r in rootVec]
end

function addroot!(diag, root::Component)
    @assert root.isNode "root of Diagrams must be a node!"
    push!(diag.root, r.index)
end

"""
    function addPropagator!(diag::Diagrams, index::Int, order::Int, name, factor = 1; site = [], loop = nothing, para = nothing)

    Add a propagator into the diagram tree.

# Arguments
- diag::Diagrams : Diagram tree.
- index::Int     : Index of the propagator in the diagram tree.
- order::Int     : Order of the propagator.
- name = :none   : name of the propagator.
- factor = 1     : Factor of the propagator.
- site = []      : site basis (e.g, time and space coordinate) of the propagator.
- loop = nothing : loop basis (e.g, momentum and frequency) of the propagator.
- para = nothing : Additional paramenter required to evaluate the propagator.
"""
function addPropagator!(diag::Diagrams, name, factor = 1.0; site = [], loop = nothing, para = nothing, order::Int = 0)
    loopPool = diag.basisPool
    propagatorPool = diag.propagatorPool
    # @assert length(basis) == length(variablePool) == length(currVar) "$(length(basis)) == $(length(variablePool)) == $(length(currVar)) breaks"

    # @assert 0 < index <= length(propagatorPool) "proapgator Pool index $index is illegal!"

    PROPAGATOR_POOL = typeof(propagatorPool)
    PROPAGATOR = eltype(fieldtype(PROPAGATOR_POOL, :object))
    PARA = fieldtype(PROPAGATOR, :para)
    F = fieldtype(PROPAGATOR, :factor)

    @assert typeof(para) <: PARA "Type of $para is $(typeof(para)), while we expect $PARA"

    # function Propagator{P,F}(order, para, factor, loopbasis, localbasis) where {P,F}
    loopidx = 0
    if isnothing(loop) == false
        @assert typeof(loop) <: AbstractVector "LoopBasis should be a Vector!"
        loopidx = append(loopPool, loop)
    end
    prop = Propagator{PARA,F}(name, order, para, factor, loopidx, collect(site))
    pidx = append(diag.propagatorPool, prop)
    # return component(pidx, false, propagatorPool[index].name)
    return pidx
end

"""
    function addNode!(diag::Diagrams, operator, name, factor = 1.0; propagator = nothing, child = [], parent = 0, para = nothing)

    Add a node into the diagram tree.

# Arguments
- diag::Diagrams : Diagram tree.
- operator::Int  : #1: multiply, 2: add, ...
- name           : name of the node
- factor = 1.0   : Factor of the node
- propagator = nothing  : Index to the cached propagators stored in certain pools. It should be in the format of Vector{Vector{Int}}, where each Vector{Int} is for one kind of propagator.
- child = []        : Indices to the cached nodes stored in certain pool. They are the child of the current node in the diagram tree. It should be in the format of Vector{Int}.
- para = nothing    : Additional paramenter required to evaluate the node. Set to nothing by default.
"""
function addNode!(diag::Diagrams, operator, name, factor = 1.0; propagator = nothing, child = [], para = nothing)
    nodePool = diag.nodePool

    if isnothing(propagator)
        propagator = []
    end

    # @assert length(propagator) == length(diag.propagatorPool) "each element of the propagator is an index vector of the corresponding propagator"

    filterZero(list) = [l for l in list if l != 0]

    #filter the propagators and nodes with index 0, they are the empty object
    propagator = filterZero(propagator)
    childNodes = filterZero(child)

    empty = true
    for c in propagator
        if isempty(c) == false
            empty = false
        end
    end
    # if all components are empty and the childnodes are empty, then no need to create new node, simply return 0
    if (empty == true) && isempty(childNodes)
        return 0
    end
    # if all components are empty && there is only one child node && factor=1, then no need to create new node, simply return the child node
    if (empty == true) && length(childNodes) == 1 && (factor ≈ 1)
        return childNodes[1]
    end

    _NodePool = typeof(nodePool)
    _Node = eltype(fieldtype(_NodePool, :object))
    PARA = fieldtype(_Node, :para)
    F = fieldtype(_Node, :factor)
    # println("node PARA: ", PARA)
    # println("node F: ", F)
    # @assert PARA == typeof(para) "Type of $para is not $PARA"

    for pidx in propagator
        @assert pidx <= length(diag.propagatorPool) "Failed to add node with propagator = $propagator, and child =$childNodes. $pidx is not in GW pool (length = $(length(diag.propagatorPool)))."
    end
    for nidx in childNodes
        @assert nidx <= length(diag.nodePool) "Failed to add node with propagator = $propagator, and child =$childNodes. $nidx is not in nodePool."
    end


    node = Node{PARA,F}(name, operator, para, propagator, childNodes, factor, 0)

    nidx = append(nodePool, node)
    return nidx
    # return component(nidx, true, :none)
end

"""
    function getNode(diag::Diagrams, nidx::Int)
    
    get Node in the diag with the index nidx.
"""
function getNode(diag::Diagrams, nidx::Int)
    return diag.nodePool.object[nidx]
end
function getNode(diag::Diagrams, n::Component)
    @assert n.isNode == true
    return diag.nodePool.object[n.index]
end

# function getPropagator(diag::Diagrams, pidx::Int)
#     for p in diag.propagatorPool
#         if p.name == poolName
#             return p.object[pidx]
#         end
#     end
#     error("$poolName propagator pool doesn't exist!")
# end

"""
    function getNodeWeight(diag::Diagrams, nidx::Int)
    
    get Node weight in the diag with the index nidx.
"""
function getNodeWeight(diag::Diagrams, nidx::Int)
    return diag.nodePool.current[nidx]
end


function addpropagator!(diag::Diagrams, name, factor = 1.0; site = [], loop = nothing, para = nothing, order::Int = 0)
    pidx = addPropagator!(diag, name, factor; site = site, loop = loop, para = para, order = order)
    @assert pidx > 0
    return Component(pidx, false, :propagator, diag.propagatorPool.object[pidx])
end

function addnode!(diag::Diagrams, operator, name, components, factor = 1.0; para = nothing)
    _components = [c for c in components if c.index > 0]
    if operator == MUL
        if length(_components) < length(components)
            return zero(Component) #if some of the components doesn't exist, then product of the components doens't exist
        end
    elseif operator == ADD
        if length(_components) == 0
            return zero(Component) #if all of the components doesn't exist, then sum of the components doens't exist
        end
    end

    child = []
    propagator = []
    for c in collect(_components)
        @assert c isa Component "$c is not a DiagTree.Component"
        if c.isNode
            push!(child, c.index)
        else
            push!(propagator, c.index)
        end
    end
    nidx = addNode!(diag, operator, name, factor; propagator = propagator, child = child, para = para)
    return Component(nidx, true, diag.nodePool.name, getNode(diag, nidx))
end