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
end

function Base.isequal(a::Propagator{P,F}, b::Propagator{P,F}) where {P,F}
    if (isequal(a.para, b.para) == false) || (a.order != b.order) || (a.variable != b.variable)
        return false
    else
        return true
    end
end
Base.:(==)(a::Propagator{P,F}, b::Propagator{P,F}) where {P,F} = Base.isequal(a, b)

struct Node{PARA,W}
    para::PARA
    operation::Int #1: multiply, 2: add, ...
    factor::W
    components::Vector{Vector{Int}}
    child::Vector{Int}
    parent::Int # parent id
    # child::SubArray{CachedObject{Node{PARA, P},W},1,Vector{CachedObject{NODE,W}},Tuple{Vector{Int64}},false}

    function Node(operation::Int; components = [[]], child = [], parent = 0, factor::W = 1.0, para::P = 0) where {W,P}
        return new{W,P}(para, operation, factor, components, child, parent)
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

mutable struct Diagrams{V,P,PARA,W}
    variable::V
    propagator::P
    tree::Pool{Node{PARA,W},W}
    root::Vector{Int}
    # root::SubArray{CachedObject{NODE,W},1,Vector{CachedObject{NODE,W}},Tuple{Vector{Int64}},false}
    #SubArray has 5 type parameters. The first two are the standard element type and dimensionality. 
    #The next is the type of the parent AbstractArray. The most heavily-used is the fourth parameter, 
    #a Tuple of the types of the indices for each dimension. The final one, L, is only provided as a convenience for dispatch;
    #it's a boolean that represents whether the index types support fast linear indexing.
    # loopNum::Int
    # tauNum::Int
    function Diagrams{PARA,W}(var::V, propagator::P) where {V,P,PARA,W}
        return new{V,P,NODE,W}(var, propagator, Pool{Node{PARA,W},W}(), [])
    end
    function Diagrams{W}(var::V, propagator::P) where {V,P,W}
        PARA = Int
        return new{V,P,PARA,W}(var, propagator, Pool{Node{PARA,W},W}(), [])
    end
end

function addPropagator(diag::Diagrams, index, order, basis, currVar, factor = 1, para = 0, curr = 0)
    variablePool = diag.variable
    propagator = diag.propagator
    # @assert length(basis) == length(variablePool) == length(currVar) "$(length(basis)) == $(length(variablePool)) == $(length(currVar)) breaks"

    PARA = fieldtype(fieldtype(eltype(fieldtype(typeof(propagator[index]), :pool)), :object), :para)
    FACTOR = fieldtype(eltype(fieldtype(typeof(propagator[index]), :pool)), :curr)

    # PARA = fieldtype(typeof(propagator[index]), 1)
    # FACTOR = fieldtype(typeof(propagator[index]), 2)
    factor = FACTOR(factor)
    para = PARA(para)
    curr = FACTOR(curr)

    vidx = zeros(length(basis))

    # if length(basis) == 1
    #     vidx[1], _ = append(variablePool, basis, currVar)
    # else
    vidx = zeros(length(basis))
    for (bi, b) in enumerate(basis)
        # TYPE = fieldtype(eltype(fieldtype(typeof(variablePool[bi]), :pool)), :curr)
        vidx[bi], _ = append(variablePool[bi], b, currVar[bi])
    end
    # end
    return append(diag.propagator[index], Propagator(order, vidx, factor, para), curr)
end

end