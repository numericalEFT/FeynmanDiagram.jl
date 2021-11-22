module DiagTree
# export PropagatorKT, Weight, addChild
export Diagrams, Momentum, Propagator, addMomentum!, addPropagator!, Node, addChild

struct Momentum
    loop::Vector{Int} #loop basis of the momentum
    current::Vector{Float64}
    new::Vector{Float64}

    version::Int128
    excited::Bool #if set to excited, then the current weight needs to be replaced with the new weight

    function Momentum(loopbasis)
        _current = loopbasis .* 0
        _new = current
        return new(loopbasis, _current, _new, 0, false)
    end
end

struct Propagator
    type::Symbol #1: Green's function, 2: interaction
    order::Int #the propagator may have an internal order (say, a Green's function diagram with multiple self-energy sub-diagrams)
    Kidx::Int #loop basis of the momentum
    Tidx::Tuple{Int,Int}

    function Propagator(_type, _order, _Kidx, _Tidx)
        return new(_type, _order, _Kidx, _Tidx)
    end
end

struct Diagrams{W}
    Gsymmetry::Vector{Symbol}
    Wsymmetry::Vector{Symbol}
    momenta::Vector{Momentum}
    propagators::Vector{Propagator}
    tree::Vector{Node{W}}
    rootIdx::Int #index of the root of the diagram tree
    # loopNum::Int
    # tauNum::Int
    function Diagrams{W}(Gsymmetry, Wsymmetry) where {W}
        momenta = Vector{Momentum}(undef, 0)
        propagators = Vector{propagators}(undef, 0)
        tree::Vector{Node{W}}(undef, 0)
        return new(Gsymmetry, Wsymmetry, momenta, propagators, tree, 0)
    end
end


mutable struct Node{W<:Number}
    type::Int #type of the weight, Green's function, interaction, node of some intermediate step
    propagatorIdx::Int #if the weight is for a propagator, then this is the index of the propagator in the propagator table
    version::Int128
    excited::Bool #if set to excited, then the current weight needs to be replaced with the new weight
    current::W
    new::W

    #### link to the other nodes ##########
    parent::Int
    child::Vector{Int} #if the Node is a leaf, then child stores the index of propagator, otherwise, it stores the indices of the child nodes
    operation::Int #0: multiply, 1: add, 2: subtract
    factor::W #symmetry factor, Fermi factor, spin factor

    function Node{W}(_parent) where {W}
        new{W}(0, 0, 0, false, W(0), W(0), _parent, [])
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
    return length(tree) #index of the new node
end



compareTidx(Tidx1, Tidx2, hasTimeReversal) = hasTimeReversal ? ((Tidx1 == Tidx2) || (Tidx1 == (Tidx2[2], Tidx2[1]))) : Tidx1 == Tidx2
compareKidx(Kidx1, Kidx2, hasMirrorSymmetry) = hasMirrorSymmetry ? ((Kidx1 ≈ Kidx2) || (Kidx1 ≈ -Kidx2)) : Kidx1 ≈ Kidx2

function addMomentum!(diagrams::Diagrams, _Kbasis)
    momenta = diagrams.momenta
    for (i, K) in enumerate(momenta)
        Kbasis = momenta[K.loop]
        if compareKidx(Kbasis, _Kbasis, :mirror in _symmetry)
            return i
        end
    end
    push!(momenta, Momentum(_Kbasis))
    return length(momenta)
end

# add new propagators to the propagator list
function addPropagator!(diagrams::Diagrams, type::Symbol, order, _Kidx, _Tidx)
    propagators = diagrams.propagators
    momenta = diagrams.momenta
    _Kbasis = momenta[_Kidx]
    if type == :G
        _symmetry = diagrams.Gsymmetry
    elseif type == :V || type == :W
        _symmetry = diagrams.Wsymmetry
    else
        error("not implemented!")
    end
    for (i, p) in enumerate(propagators)
        Kbasis = momenta[p.Kidx]
        if p.type == type && p.order == order
            Tflag = compareTidx(p.Tidx, _Tidx, (:particlehole in _symmetry || :timereversal in _symmetry))
            Kflag = compareKidx(Kbasis, _Kbasis, :mirror in _symmetry)
            if Tflag && Kflag
                return i
            end
        end
    end
    push!(propagators, Propagator(type, order, _Kidx, _Tidx))
    return length(propagators)
end


end