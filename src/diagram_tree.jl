module DiagTree
include("pool.jl")

mutable struct Node{W,PARA,COM}
    para::PARA
    factor::Float64 #symmetry factor, Fermi factor, spin factor
    operation::Int #1: multiply, 2: add

    #### link to the other nodes and components ##########
    parent::Int # parent id
    nodes::Vector{Int} # store the indices of the child nodes
    components::COM # tuple of Vector{Int} that host other components

    function Node{W}(_parent = 0) where {W}
        new{W}(1, 0, 1.0, 0, false, [], [], W(0), W(0), _parent, [], [])
    end
    function Node{W}(id, operation::Int, factor, propagators = [], nodes = []; extK = [], extT = [], parent = 0) where {W}
        new{W}(id, operation, factor, 0, false, collect(extK), collect(extT), W(0), W(0), parent, propagators, nodes)
    end
end

isleaf(node) = (length(node.nodes) == 0) #if a node doesn't have child nodes, it must be a leaf (it must be ponting to propagators)

mutable struct Diagrams{W}
    momenta::Vector{Momentum}
    propagators::Vector{Propagator}
    tree::Vector{Node{W}}
    root::Vector{Int} #index of the root of the diagram tree
    # loopNum::Int
    # tauNum::Int
    function Diagrams{W}() where {W}
        momenta = Vector{Momentum}(undef, 0)
        propagators = Vector{Propagator}(undef, 0)
        tree = Vector{Node{W}}(undef, 0)
        return new(momenta, propagators, tree, [])
    end
end


end