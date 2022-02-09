"""
    mutable struct Diagrams{V,P,PARA,F,W}

    Diagram Object represents a set of Feynman diagrams in a tree (forest) structure

# Members
- basisPool::V      : Tuple of pools of cached basis  in a format of (BasisPool1, BasisPool2, ...)
- propagatorPool::P : Tuple of pools of cached propagators in a format of (PropagatorPool1, PropagatorPool2, ...)
- nodePool::Pool{Node{PARA,F},W} : Pool of the nodes in the diagram tree
- root::Vector{Int} : indices of the cached nodes that are the root(s) of the diagram tree. Each element corresponds to one root.
"""
mutable struct ExpressionTree{V,pPARA,nPARA,F,W}
    name::Symbol
    loopBasis::V
    propagator::CachedPool{Propagator{pPARA,F},W}
    # node::CachedPool{Node{nPARA,F},W}
    node::CachedPool{Any,W}
    root::Vector{Int}
    function ExpressionTree(; loopBasis::V, weight::DataType, factor::DataType = weight, nodePara::DataType = Nothing, propagatorPara::DataType = Nothing, name = :none) where {V,P}
        # nodePool = CachedPool(:node, Node{nodePara,factor}, weight)
        nodePool = CachedPool(:node, Any, weight)
        propagatorPool = CachedPool(:propagator, Propagator{propagatorPara,factor}, weight)
        return new{V,propagatorPara,nodePara,factor,weight}(name, loopBasis, propagatorPool, nodePool, [])
    end
end

weight(tree::ExpressionTree) = tree.node.current

Base.getindex(diag::ExpressionTree, i) = diag.node.current[diag.root[i]]
Base.firstindex(diag::ExpressionTree) = 1
Base.lastindex(diag::ExpressionTree) = length(diag.root)

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
    function addPropagator!(diag::ExpressionTree, name, factor = 1.0; site = [], loop = nothing, para = nothing, order::Int = 0)

    Add a propagator into the diagram tree.

# Arguments
- diag           : diagrammatic experssion tree.
- order::Int = 0 : Order of the propagator.
- name = :none   : name of the propagator.
- factor = 1     : Factor of the propagator.
- site = []      : site basis (e.g, time and space coordinate) of the propagator.
- loop = nothing : loop basis (e.g, momentum and frequency) of the propagator.
- para = nothing : Additional paramenter required to evaluate the propagator.
"""
function addPropagator!(diag::ExpressionTree{V,pPARA,nPARA,F,W}, name, factor = 1.0; site = [], loop = nothing, para = nothing, order::Int = 0) where {V,pPARA,nPARA,F,W}
    loopPool = diag.loopBasis
    propagatorPool = diag.propagator
    # @assert length(basis) == length(variablePool) == length(currVar) "$(length(basis)) == $(length(variablePool)) == $(length(currVar)) breaks"

    # @assert 0 < index <= length(propagatorPool) "proapgator Pool index $index is illegal!"

    # PROPAGATOR_POOL = typeof(propagatorPool)
    # PROPAGATOR = eltype(fieldtype(PROPAGATOR_POOL, :object))
    # PARA = fieldtype(PROPAGATOR, :para)
    # F = fieldtype(PROPAGATOR, :factor)

    # @assert typeof(para) <: PARA "Type of $para is $(typeof(para)), while we expect $PARA"

    # function Propagator{P,F}(order, para, factor, loopbasis, localbasis) where {P,F}
    loopidx = 0
    if isnothing(loop) == false
        @assert typeof(loop) <: AbstractVector "LoopBasis should be a Vector!"
        loopidx = append(loopPool, loop)
    end
    prop = Propagator{pPARA,F}(name, order, para, factor, loopidx, collect(site))
    # pidx = append(diag.propagator, prop)
    pidx = append(diag.node, prop)
    # return component(pidx, false, propagatorPool[index].name)
    return pidx
end

"""
    function addNode!(diag::ExpressionTree, operator, name, factor = 1.0; propagator = nothing, child = [], para = nothing)

    Add a node into the diagram tree.

# Arguments
- diag::ExpressionTree  : diagrammatic experssion tree.
- operator::Int         : #1: multiply, 2: add, ...
- name                  : name of the node
- factor = 1.0          : Factor of the node
- propagator = nothing  : Index to the cached propagators stored in certain pools. It should be in the format of Vector{Vector{Int}}, where each Vector{Int} is for one kind of propagator.
- child = []            : Indices to the cached nodes stored in certain pool. They are the child of the current node in the diagram tree. It should be in the format of Vector{Int}.
- para = nothing        : Additional paramenter required to evaluate the node. Set to nothing by default.
"""
function addNode!(diag::ExpressionTree{V,pPARA,nPARA,F,W}, operator, name, factor = 1.0; propagator = nothing, child = [], para = nothing) where {V,pPARA,nPARA,F,W}
    nodePool = diag.node

    if isnothing(propagator)
        propagator = []
    end

    # @assert length(propagator) == length(diag.propagatorPool) "each element of the propagator is an index vector of the corresponding propagator"

    # filterZero(list) = [l for l in list if l != 0]

    # #filter the propagators and nodes with index 0, they are the empty object
    # propagator = filterZero(propagator)
    # childNodes = filterZero(child)

    # empty = true
    # for c in propagator
    #     if isempty(c) == false
    #         empty = false
    #     end
    # end
    # # if all components are empty and the childnodes are empty, then no need to create new node, simply return 0
    # if (empty == true) && isempty(childNodes)
    #     return 0
    # end
    # # if all components are empty && there is only one child node && factor=1, then no need to create new node, simply return the child node
    # if (empty == true) && length(childNodes) == 1 && (factor ≈ 1)
    #     return childNodes[1]
    # end

    # _NodePool = typeof(nodePool)
    # println(_NodePool)
    # _Node = eltype(fieldtype(_NodePool, :object))
    # println(_Node)
    # PARA = fieldtype(_Node, :para)
    # F = fieldtype(_Node, :factor)
    # # println("node PARA: ", PARA)
    # # println("node F: ", F)
    # # @assert PARA == typeof(para) "Type of $para is not $PARA"

    # for pidx in propagator
    #     @assert pidx <= length(diag.propagator) "Failed to add node with propagator = $propagator, and child =$childNodes. $pidx is not in GW pool (length = $(length(diag.propagator)))."
    # end

    children = deepcopy(propagator)
    append!(children, child)
    for nidx in children
        @assert nidx <= length(diag.node) "Failed to add node with propagator = $propagator, and child =$children. $nidx is not in nodePool."
    end

    node = Node{nPARA,F}(name, operator, para, [], children, factor, 0)

    nidx = append(nodePool, node)
    return nidx
    # return component(nidx, true, :none)
end

"""
    function getNode(diag::Diagrams, nidx::Int)
    
    get Node in the diag with the index nidx.
"""
function getNode(diag, nidx::Int)
    return diag.node.object[nidx]
end
function getNode(diag, n::Component)
    @assert n.isNode == true
    return diag.node.object[n.index]
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
    function getNodeWeight(tree, nidx::Int)
    
    get Node weight in the diagram experssion tree with the index nidx.
"""
function getNodeWeight(tree, nidx::Int)
    return tree.node.current[nidx]
end


function addpropagator!(diag, name, factor = 1.0; site = [], loop = nothing, para = nothing, order::Int = 0)
    pidx = addPropagator!(diag, name, factor; site = site, loop = loop, para = para, order = order)
    @assert pidx > 0
    return Component(pidx, false, :propagator, diag.node.object[pidx])
end

function addnode!(diag, operator, name, components, factor = 1.0; para = nothing)
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
    return Component(nidx, true, diag.node.name, getNode(diag, nidx))
end