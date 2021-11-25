module DiagTree
# export PropagatorKT, Weight, addChild
export Diagrams, Momentum, Propagator, addMomentum!, addPropagator!, Node, addChild

struct Momentum
    basis::Vector{Int} #loop basis of the momentum
    curr::Vector{Float64}
    new::Vector{Float64}

    version::Int128
    excited::Bool #if set to excited, then the current weight needs to be replaced with the new weight

    function Momentum(loopbasis)
        _current = loopbasis .* 0
        _new = _current
        return new(loopbasis, _current, _new, 1, false)
    end
end

struct Propagator{W}
    type::Int #1: Green's function, 2: interaction
    order::Int #the propagator may have an internal order (say, a Green's function diagram with multiple self-energy sub-diagrams)
    Kidx::Int #loop basis of the momentum
    Tidx::Tuple{Int,Int}

    version::Int128
    excited::Bool #if set to excited, then the current weight needs to be replaced with the new weight
    curr::W
    new::W

    function Propagator{W}(_type::Int, _order, _Kidx, _Tidx)
        return new(_type, _order, _Kidx, Tuple(_Tidx), 0, false, W(0), W(0))
    end
end

mutable struct Node{W}
    operation::Int #0: multiply, 1: add, 2: subtract
    factor::W #symmetry factor, Fermi factor, spin factor
    version::Int128
    excited::Bool #if set to excited, then the current weight needs to be replaced with the new weight
    curr::W
    new::W

    #### link to the other nodes ##########
    parent::Int
    propagators::Vector{Int}
    nodes::Vector{Int} #if the Node is a leaf, then child stores the index of propagator, otherwise, it stores the indices of the child nodes

    function Node{W}(_parent) where {W}
        new{W}(0, 1.0, 0, false, W(0), W(0), _parent, [], [])
    end
end

struct Diagrams{W}
    momenta::Vector{Momentum}
    propagators::Vector{Propagator}
    tree::Vector{Node{W}}
    rootIdx::Int #index of the root of the diagram tree
    # loopNum::Int
    # tauNum::Int
    function Diagrams{W}() where {W}
        momenta = Vector{Momentum}(undef, 0)
        propagators = Vector{Propagator}(undef, 0)
        tree = Vector{Node{W}}(undef, 0)
        return new(momenta, propagators, tree, 0)
    end
end



Base.show(io::IO, w::Node{W}) where {W} = print(io, "Type: $(w.type), Prop: $(w.propagatorIdx), curr: $(w.current)")

function addChild(tree::Vector{Node{W}}, _parent) where {W}
    push!(tree, Node{W}(_parent))
    idx = length(tree) #the index of the last element, which is the new weight node
    push!(tree[_parent].child, idx)
    return idx
end

function addNode!(diagrams::Diagrams{W}, propagatorIdx) where {W}
    tree = diagrams.tree
    node = Node{W}(-1)
    node.propagatorIdx = propagatorIdx
    push!(tree, node)
    diagrams.rootIdx = length(tree)
    return length(tree) #index of the new node
end

compareTidx(Tidx1, Tidx2, hasTimeReversal) = hasTimeReversal ? ((Tidx1 == Tidx2) || (Tidx1 == (Tidx2[2], Tidx2[1]))) : Tidx1 == Tidx2
compareKidx(Kidx1, Kidx2, hasMirrorSymmetry) = hasMirrorSymmetry ? ((Kidx1 ≈ Kidx2) || (Kidx1 ≈ -Kidx2)) : Kidx1 ≈ Kidx2

function addMomentum!(diagrams::Diagrams, _Kbasis, _symmetry)
    momenta = diagrams.momenta
    for (i, K) in enumerate(momenta)
        if compareKidx(K.basis, _Kbasis, :mirror in _symmetry)
            return i, false #existing momentum
        end
    end
    push!(momenta, Momentum(_Kbasis))
    return length(momenta), true #new momentum
end

# add new propagators to the propagator list
function addPropagator!(diagrams::Diagrams, type::Int, order, _Kbasis, _Tidx, _symmetry = [])
    propagators = diagrams.propagators

    _Kidx, isNewK = addMomentum!(diagrams, _Kbasis, _symmetry)

    # for (i, n) in enumerate(diagrams.tree)
    #     p = n.
    # end

    for (i, p) in enumerate(propagators)
        if p.type == type && p.order == order
            Tflag = compareTidx(p.Tidx, _Tidx, (:particlehole in _symmetry || :timereversal in _symmetry))
            if Tflag && p.Kidx == _Kidx
                return i, false #existing propagator
            end
        end
    end
    push!(propagators, Propagator(type, order, _Kidx, _Tidx))
    return length(propagators), true #new propagator
end


end