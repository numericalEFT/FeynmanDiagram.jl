module DiagTree
include("pool.jl")

# struct Node{PARA,P}
#     para::PARA
#     operation::Int #1: multiply, 2: add, ...
#     propagator::P # tuple of CachedObject that host possible components (Nodes, propagators, ...)

#     parent::Int # parent id
#     # child::SubArray{CachedObject{Node{PARA, P},W},1,Vector{CachedObject{NODE,W}},Tuple{Vector{Int64}},false}

#     function Node(para::P, operation::Int, branch::B, parent = 0) where {P,B}
#         return new{P,B}(para, operation, branch, parent)
#     end
# end

# mutable struct Diagrams{V,P,NODE,W}
#     variable::V
#     propagator::P
#     tree::Pool{NODE,W}
#     # root::SubArray{CachedObject{NODE,W},1,Vector{CachedObject{NODE,W}},Tuple{Vector{Int64}},false}
#     #SubArray has 5 type parameters. The first two are the standard element type and dimensionality. 
#     #The next is the type of the parent AbstractArray. The most heavily-used is the fourth parameter, 
#     #a Tuple of the types of the indices for each dimension. The final one, L, is only provided as a convenience for dispatch;
#     #it's a boolean that represents whether the index types support fast linear indexing.
#     # loopNum::Int
#     # tauNum::Int
#     function Diagrams(var::V, propagator::V) where {W}
#         momenta = Vector{Momentum}(undef, 0)
#         propagators = Vector{Propagator}(undef, 0)
#         tree = Vector{Node{W}}(undef, 0)
#         return new(momenta, propagators, tree, [])
#     end
# end

end