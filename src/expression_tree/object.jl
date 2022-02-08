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

struct Component
    index::Int
    isNode::Bool
    poolName::Symbol
    object::Any
end

Base.show(io::IO, c::Component) = print(io, "#$(c.index) in $(c.poolName) Pool")

Base.zero(::Type{Component}) = Component(0, false, :none, nothing)
