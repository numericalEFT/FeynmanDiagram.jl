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
    loopidx::Int
    siteidx::Vector{Int}
    children::Vector{Int}
    function Node{P,F}(name::Symbol, para; loopidx = 0, siteidx = [], operator = ADD, children = [], factor = 1.0) where {F,P}
        # @assert typeof(para) == P
        return new{P,F}(name, para, operator, F(factor), loopidx, siteidx, children)
    end
end

function Base.isequal(a::Node{P}, b::Node{P}) where {P}
    # only parent is allowed to be different
    if (isequal(a.para, b.para) == false) || (a.operation != b.operation) || (Set(a.children) != Set(b.children)) || (a.factor ≈ b.factor) == false || a.loopidx != b.loopidx || a.siteidx != b.siteidx
        return false
    else
        return true
    end
end
Base.:(==)(a::Node{P}, b::Node{P}) where {P} = Base.isequal(a, b)

"""
    mutable struct ExpressionTree{V,PARA,F,W}

    Diagram Object represents a set of Feynman diagrams in an experssion tree (forest) structure

# Members
- name::Symbol                     : Name of the tree
- loopBasis::V                     : Tuple of pools of cached basis  in a format of (BasisPool1, BasisPool2, ...)
- node::CachedPool{Node{PARA,F},W} : Pool of the nodes in the diagram tree
- root::Vector{Int}                : indices of the cached nodes that are the root(s) of the diagram tree. Each element corresponds to one root.
"""
mutable struct ExpressionTree{V,PARA,F,W}
    name::Symbol
    loopBasis::V
    node::CachedPool{Node{PARA,F},W}
    root::Vector{Int}
    function ExpressionTree(; loopBasis::V, weight::DataType, factor::DataType = weight, nodePara::DataType = Nothing, name = :none) where {V}
        nodePool = CachedPool(:node, Node{nodePara,factor}, weight)
        return new{V,nodePara,factor,weight}(name, loopBasis, nodePool, [])
    end
end

weight(tree::ExpressionTree) = tree.node.current

Base.getindex(diag::ExpressionTree, i) = diag.node.current[diag.root[i]]
Base.firstindex(diag::ExpressionTree) = 1
Base.lastindex(diag::ExpressionTree) = length(diag.root)

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
function addpropagator!(diag::ExpressionTree{V,PARA,F,W}, name, factor = 1.0; site = [], loop = nothing, para = nothing, order::Int = 0) where {V,PARA,F,W}
    loopPool = diag.loopBasis
    loopidx = 0
    if isnothing(loop) == false
        @assert typeof(loop) <: AbstractVector "LoopBasis should be a Vector!"
        loopidx = append(loopPool, loop)
    end
    # prop = Propagator{pPARA,F}(name, order, para, factor, loopidx, collect(site))
    prop = Node{PARA,F}(name, para; factor = factor, loopidx = loopidx, siteidx = collect(site))
    # pidx = append(diag.propagator, prop)
    pidx = append(diag.node, prop)
    # return component(pidx, false, propagatorPool[index].name)
    return pidx
end

"""
    function addnode!(diag::ExpressionTree{V,PARA,F,W}, operator, name, children::Union{Tuple, AbstractVector}, factor = 1.0; para = nothing) where {V,PARA,F,W}

    Add a node into the expression tree.

# Arguments
- diag::ExpressionTree  : diagrammatic experssion tree.
- operator::Int         : #1: multiply, 2: add, ...
- name                  : name of the node
- children              : Indices to the cached nodes stored in certain pool. They are the child of the current node in the diagram tree. It should be in the format of Vector{Int}.
- factor = 1.0          : Factor of the node
- para = nothing        : Additional paramenter required to evaluate the node. Set to nothing by default.
"""
function addnode!(diag::ExpressionTree{V,PARA,F,W}, operator, name, children::Union{Tuple,AbstractVector}, factor = 1.0; para = nothing) where {V,PARA,F,W}
    nodePool = diag.node

    for nidx in children
        @assert nidx <= length(diag.node) "Failed to add node with propagator = $propagator, and child =$children. $nidx is not in nodePool."
    end

    node = Node{PARA,F}(name, para; operator = operator, children = children, factor = factor)

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


# function addpropagator!(diag, name, factor = 1.0; site = [], loop = nothing, para = nothing, order::Int = 0)
#     pidx = addPropagator!(diag, name, factor; site = site, loop = loop, para = para, order = order)
#     @assert pidx > 0
#     return Component(pidx, false, :propagator, diag.node.object[pidx])
# end

# function addnode!(diag, operator, name, components, factor = 1.0; para = nothing)
#     _components = [c for c in components if c.index > 0]
#     if operator == MUL
#         if length(_components) < length(components)
#             return zero(Component) #if some of the components doesn't exist, then product of the components doens't exist
#         end
#     elseif operator == ADD
#         if length(_components) == 0
#             return zero(Component) #if all of the components doesn't exist, then sum of the components doens't exist
#         end
#     end

#     child = []
#     propagator = []
#     for c in collect(_components)
#         @assert c isa Component "$c is not a DiagTree.Component"
#         if c.isNode
#             push!(child, c.index)
#         else
#             push!(propagator, c.index)
#         end
#     end
#     nidx = addNode!(diag, operator, name, factor; propagator = propagator, child = child, para = para)
#     return Component(nidx, true, diag.node.name, getNode(diag, nidx))
# end