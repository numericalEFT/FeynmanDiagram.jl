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
    name::Symbol
    para::PARA
    order::Int
    factor::F
    loopIdx::Int
    siteBasis::Vector{Int}
    # function Propagator(order, basis = [], factor::F = 1.0, para::P = 0) where {F,P}
    #     return new{P,F}(para, order, factor, basis)
    # end
    function Propagator{F}(name::Symbol, order::Int, para::P, factor, loopidx::Int, sitebasis) where {F,P}
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
    propagators::Vector{Vector{Int}}
    childNodes::Vector{Int}
    parent::Int # parent id
    # child::SubArray{CachedObject{Node{PARA, P},W},1,Vector{CachedObject{NODE,W}},Tuple{Vector{Int64}},false}

    # function Node(operation::Int, components = [[]], child = [], factor::F = 1.0, parent = 0, para::P = :N) where {F,P}
    #     return new{P,F}(para, operation, factor, components, child, parent)
    # end
    function Node{F}(name::Symbol, operation::Int, para::P, propagators = [[]], child = [], factor = 1.0, parent = 0) where {F,P}
        @assert typeof(para) == P
        return new{P,F}(name, para, operation, F(factor), propagators, child, parent)
    end
end

function Base.isequal(a::Node{P}, b::Node{P}) where {P}
    # only parent is allowed to be different
    if (isequal(a.para, b.para) == false) || (a.operation != b.operation) || (a.propagators != b.propagators) || (a.childNodes != b.childNodes) || (a.factor != b.factor)
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
        # typeof(nothing) = Nothing
        # @assert V <: Tuple "Tuple is required for efficiency!"
        @assert P <: Tuple "Tuple is required for efficiency!"
        # println(basisPool)
        # println(propagatorPool)
        poolname = []
        for p in propagatorPool
            push!(poolname, p.name)
        end
        @assert length(poolname) == length(Set(poolname)) "All propagatorPool should have different names, yet got $poolname!"
        nodePool = CachedPool(:node, Node{nodeParaType,nodeFactorType}, nodeWeightType)
        return new{V,P,nodeParaType,nodeFactorType,nodeWeightType}(name, basisPool, propagatorPool, nodePool, [])
    end
end

struct component
    index::Int
    isNode::Bool
    propagatorPoolName::Symbol
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
function addPropagator!(diag::Diagrams, index::Int, order::Int, name, factor = 1.0; site = [], loop = nothing, para = nothing)
    loopPool = diag.basisPool
    propagatorPool = diag.propagatorPool
    # @assert length(basis) == length(variablePool) == length(currVar) "$(length(basis)) == $(length(variablePool)) == $(length(currVar)) breaks"

    PROPAGATOR_POOL = typeof(propagatorPool[index])
    PROPAGATOR = eltype(fieldtype(PROPAGATOR_POOL, :object))
    PARA = fieldtype(PROPAGATOR, :para)
    F = fieldtype(PROPAGATOR, :factor)

    @assert PARA == typeof(para) "Type of $para is not $PARA"

    # function Propagator{P,F}(order, para, factor, loopbasis, localbasis) where {P,F}
    loopidx = 0
    if isnothing(loop) == false
        @assert typeof(loop) <: AbstractVector "LoopBasis should be a Vector!"
        loopidx = append(loopPool, loop)
    end
    prop = Propagator{F}(name, order, para, factor, loopidx, collect(site))
    pidx = append(diag.propagatorPool[index], prop)
    # return component(pidx, false, propagatorPool[index].name)
    return pidx
end
"""
    function addPropagator!(diag::Diagrams, poolName::Symbol, order::Int, name, factor = 1; site = [], loop = nothing, para = nothing)

    Add a propagator into the diagram tree referenced by the pool name.

# Arguments
- diag::Diagrams : Diagram tree.
- poolName::Symbol : name of the propagator pool
- order::Int     : Order of the propagator.
- name = :none   : name of the propagator.
- factor = 1     : Factor of the propagator.
- site = []      : site basis (e.g, time and space coordinate) of the propagator.
- loop = nothing : loop basis (e.g, momentum and frequency) of the propagator.
- para = nothing : Additional paramenter required to evaluate the propagator.
"""
function addPropagator!(diag::Diagrams, poolName::Symbol, order::Int, name, factor = 1.0; site = [], loop = nothing, para = nothing)
    for (idx, pool) in enumerate(diag.propagatorPool)
        if pool.name == poolName
            return addPropagator!(diag, idx, order, name, factor; site = site, loop = loop, para = para)
        end
    end
end

function addpropagator!(diag::Diagrams, poolName::Symbol, order::Int, name, factor = 1.0; site = [], loop = nothing, para = nothing)
    pidx = addPropagator!(diag, poolName, order, name, factor; site = site, loop = loop, para = para)
    return component(pidx, false, poolName)
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
        propagator = [[] for p in diag.propagatorPool]
    end

    @assert length(propagator) == length(diag.propagatorPool) "each element of the propagator is an index vector of the corresponding propagator"

    filterZero(list) = [l for l in list if l != 0]

    #filter the propagators and nodes with index 0, they are the empty object
    propagator = [filterZero(c) for c in propagator]
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
    @assert PARA == typeof(para) "Type of $para is not $PARA"

    for (i, p) in enumerate(propagator)
        for pidx in p
            @assert pidx <= length(diag.propagatorPool[i]) "Failed to add node with propagator = $propagator, and child =$childNodes. $pidx is not in pool $i (length = $(length(diag.propagatorPool[i])))."
        end
    end
    for nidx in childNodes
        @assert nidx <= length(diag.nodePool) "Failed to add node with propagator = $propagator, and child =$childNodes. $nidx is not in nodePool."
    end


    node = Node{F}(name, operator, para, propagator, childNodes, factor, 0)

    nidx = append(nodePool, node)
    return nidx
    # return component(nidx, true, :none)
end
"""
    function addNodeByName(diag::Diagrams, operator, name = :none; child = [], factor = 1.0, parent = 0, para = nothing, kwargs...)

    An alternative way to add node. Instead of providing the components vector, one specify the propagators in kwargs. 

# Arguments
- diag::Diagrams : Diagram tree.
- operator::Int  : #1: multiply, 2: add, ...
- name           : name of the node
- factor = 1.0   : Factor of the node
- child = []        : Indices to the cached nodes stored in certain pool. They are the child of the current node in the diagram tree. It should be in the format of Vector{Int}.
- para = nothing    : Additional paramenter required to evaluate the node. Set to nothing by default.
- kwargs...      : Accept arbitrary pair of parameters as: PropagatorPoolName = [index of the propagator1, index of the propagator2, ...]
"""
function addNodeByName!(diag::Diagrams, operator, name, factor = 1.0; child = [], para = nothing, kwargs...)
    components = [[] for p in diag.propagatorPool]
    for key in keys(kwargs)
        for (idx, p) in enumerate(diag.propagatorPool)
            if p.name == key
                append!(components[idx], kwargs[key])
            end
        end
    end
    # println(kwargs, " got ", components)
    return addNode!(diag, operator, name, factor; propagator = components, child = child, para = para)
end

function addnode!(diag::Diagrams, operator, name, factor = 1.0; components = Vector{component}([]), parent = 0, para = nothing)
    child = []
    propagator = [[] for p in diag.propagatorPool]
    if isempty(components) == false
        for c in components
            if c.isNode
                push!(child, c.index)
            else
                for (idx, p) in enumerate(diag.propagatorPool)
                    if p.name == c.propagatorPoolName
                        append!(propagator[idx], c.index)
                    end
                end
            end
        end
    end
    nidx = addNode!(diag, operator, name, factor; propagator = propagator, child = child, para = para)
    return component(nidx, true, :none)
end

"""
    function getNode(diag::Diagrams, nidx::Int)
    
    get Node in the diag with the index nidx.
"""
function getNode(diag::Diagrams, nidx::Int)
    return diag.nodePool.object[nidx]
end

"""
    function getNodeWeight(diag::Diagrams, nidx::Int)
    
    get Node weight in the diag with the index nidx.
"""
function getNodeWeight(diag::Diagrams, nidx::Int)
    return diag.nodePool.current[nidx]
end