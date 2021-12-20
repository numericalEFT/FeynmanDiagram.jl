module DiagTree
include("pool.jl")
using ..Var


struct Propagator{PARA,F}
    para::PARA
    order::Int
    factor::F
    variable::Vector{Int}
    function Propagator(order, variable = [], factor::F = 1.0, para::P = 0) where {F,P}
        return new{P,F}(para, variable, order)
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

# function addPropagator(diag::Diagrams, index, basis, curr, order, factor::F = 1.0, para::P = 0) where {F,P}
#     @assert length(basis) == length(curr) == length(variable)

#     PARA = fieldtype(propagator[i], 1)
#     FACTOR = fieldtype(propagator[i], 2)
#     factor = FACTOR(factor)
#     para = PARA(para)

#     if length(basis) == 1
#         vidx = append(variable, basis, curr)
#         return append(diag.propagator[index], Propagator{PARA,FACTOR}(order, [vidx,], factor, para), curr)
#     else
#         for (b, bi) in enumerate(basis)
#             # @assert typeof(b) ==  
#             append(variable[bi], basis[bi], curr[bi])
#         end
#     end
# end



end