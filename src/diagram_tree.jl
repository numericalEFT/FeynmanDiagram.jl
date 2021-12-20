module DiagTree
include("pool.jl")
using ..Var


struct Propagator{PARA}
    para::PARA
    order::Int
    variable::Vector{Int}
    function Propagator(order, variable = [], para::P = 0) where {P}
        return new{P}(para, variable, order)
    end
end

function Base.isequal(a::Propagator{P}, b::Propagator{P}) where {P}
    if (isequal(a.para, b.para) == false) || (a.order != b.order) || (a.variable != b.variable)
        return false
    else
        return true
    end
end
Base.:(==)(a::Propagator{P}, b::Propagator{P}) where {P} = Base.isequal(a, b)


struct Node{PARA}
    para::PARA
    operation::Int #1: multiply, 2: add, ...
    components::Vector{Vector{Int}}
    child::Vector{Int}
    parent::Int # parent id
    # child::SubArray{CachedObject{Node{PARA, P},W},1,Vector{CachedObject{NODE,W}},Tuple{Vector{Int64}},false}

    function Node(operation::Int; components = [[]], child = [], parent = 0, para::P = 0) where {P}
        return new{P}(para, operation, components, child, parent)
    end
end

function Base.isequal(a::Node{P}, b::Node{P}) where {P}
    # only parent is allowed to be different
    if (isequal(a.para, b.para) == false) || (a.operation != b.operation) || (a.components != b.components) || (a.child != b.child)
        return false
    else
        return true
    end
end
Base.:(==)(a::Node{P}, b::Node{P}) where {P} = Base.isequal(a, b)

mutable struct Diagrams{V,P,NODE,W}
    variable::V
    propagator::P
    tree::Pool{NODE,W}
    root::Vector{Int}
    # root::SubArray{CachedObject{NODE,W},1,Vector{CachedObject{NODE,W}},Tuple{Vector{Int64}},false}
    #SubArray has 5 type parameters. The first two are the standard element type and dimensionality. 
    #The next is the type of the parent AbstractArray. The most heavily-used is the fourth parameter, 
    #a Tuple of the types of the indices for each dimension. The final one, L, is only provided as a convenience for dispatch;
    #it's a boolean that represents whether the index types support fast linear indexing.
    # loopNum::Int
    # tauNum::Int
    function Diagrams{NODE,W}(var::V, propagator::P) where {V,P,NODE,W}
        # momenta = Vector{Momentum}(undef, 0)
        # propagators = Vector{Propagator}(undef, 0)
        # tree = Vector{Node{W}}(undef, 0)
        return new{V,P,NODE,W}(var, propagator, Pool{NODE,W}(), [])
    end
end

end