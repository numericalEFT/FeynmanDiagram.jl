module DiagTree
include("pool.jl")

# mutable struct Node{W,PARA,COM}
#     para::PARA
#     factor::Float64 #symmetry factor, Fermi factor, spin factor
#     operation::Int #1: multiply, 2: add

#     #### link to the other nodes and components ##########
#     parent::Int # parent id
#     nodes::Vector{Int} # store the indices of the child nodes
#     components::COM # tuple of Vector{Int} that host other components

#     function Node{W}(_parent = 0) where {W}
#         new{W}(1, 0, 1.0, 0, false, [], [], W(0), W(0), _parent, [], [])
#     end
#     function Node{W}(id, operation::Int, factor, propagators = [], nodes = []; extK = [], extT = [], parent = 0) where {W}
#         new{W}(id, operation, factor, 0, false, collect(extK), collect(extT), W(0), W(0), parent, propagators, nodes)
#     end
# end

end