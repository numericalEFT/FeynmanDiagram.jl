module DiagTree
using PyCall, Printf
# export PropagatorKT, Weight, addChild
export Diagrams, Momentum, Propagator, addMomentum!, addPropagator!, Node, addChild

mutable struct Momentum
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

mutable struct Propagator{W}
    type::Int #1: Green's function, 2: interaction
    order::Int #the propagator may have an internal order (say, a Green's function diagram with multiple self-energy sub-diagrams)
    Kidx::Int #loop basis of the momentum
    Tidx::Tuple{Int,Int}
    factor::Float64 #additional factor, for example, if the propagator is an exchange interaction, then factor = -1.0

    version::Int128
    excited::Bool #if set to excited, then the current weight needs to be replaced with the new weight
    curr::W
    new::W

    function Propagator{W}(_type::Int, _order::Int, _Kidx, _Tidx, _factor = 1.0) where {W}
        return new{W}(_type, _order, _Kidx, Tuple(_Tidx), _factor, 0, false, W(0), W(0))
    end
end

mutable struct Node{W}
    id::Int
    operation::Int #1: multiply, 2: add
    factor::Float64 #symmetry factor, Fermi factor, spin factor
    version::Int128
    excited::Bool #if set to excited, then the current weight needs to be replaced with the new weight
    extK::Vector{Int}
    extT::Vector{Int}
    curr::W
    new::W

    #### link to the other nodes ##########
    parent::Int
    propagators::Vector{Int}
    nodes::Vector{Int} #if the Node is a leaf, then child stores the index of propagator, otherwise, it stores the indices of the child nodes

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



Base.show(io::IO, w::Node{W}) where {W} = print(io, "Type: $(w.type), Prop: $(w.propagatorIdx), curr: $(w.current)")

function addChild(tree::Vector{Node{W}}, _parent) where {W}
    push!(tree, Node{W}(_parent))
    idx = length(tree) #the index of the last element, which is the new weight node
    push!(tree[_parent].child, idx)
    return idx
end

function addNode!(diagrams::Diagrams{W}, operation::Int, factor, propagators = [], nodes = []; isRoot = false, extT = [], extK = []) where {W}
    if isempty(propagators) && isempty(nodes)
        return -1
    end
    tree = diagrams.tree
    idx = length(tree) + 1
    # make sure that the propagators already exist
    for p in propagators
        @assert (1 <= p <= length(diagrams.propagators)) || p == -1
    end
    # make sure that the nodes already exist
    for n in nodes
        @assert (1 <= n <= length(diagrams.tree)) || n == -1
    end
    node = Node{W}(idx, operation, factor, propagators, nodes; extT = extT, extK = extK)
    push!(tree, node)
    if isRoot
        push!(diagrams.root, idx)
    end
    return idx #index of the new node
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
function addPropagator!(diagrams::Diagrams{W}, type::Int, order::Int, _Kbasis, _Tidx, symmetry = [], factor = 1.0) where {W}
    propagators = diagrams.propagators
    _Kidx, isNewK = addMomentum!(diagrams, _Kbasis, symmetry)

    for (i, p) in enumerate(propagators)
        if p.type == type && p.order == order
            Tflag = compareTidx(p.Tidx, _Tidx, (:particlehole in symmetry || :timereversal in symmetry))
            if Tflag && p.Kidx == _Kidx && p.factor ≈ factor
                return i, false #existing propagator
            end
        end
    end
    push!(propagators, Propagator{W}(type, order, _Kidx, _Tidx, factor))
    return length(propagators), true #new propagator
end

# function addPropagator!(diagrams::Diagrams{W}, type::Int, order::Int, _Kbasis, _Tidx::Vector{Tuple{Int,Int}}, _symmetry = [], _factor = 1.0) where {W}
#     @assert length(_Kbasis) == length(_Tidx)
#     pidx = []
#     for (ti, tidx) in enumerate(_Tidx)
#         kbasis = _Kbasis[ti]
#     end
# end

"""
    showTree(ver4, para::Para; verbose=0, depth=999)

    Visualize the diagram tree using ete3 python package

#Arguments
- `ver4`: the 4-vertex diagram tree to visualize
- `para`: parameters
- `verbose=0`: the amount of information to show
- `depth=999`: deepest level of the diagram tree to show
"""
function showTree(diag::Diagrams, _root = diag.root[end]; verbose = 0, depth = 999)
    tree = diag.tree
    K = diag.momenta
    root = tree[_root]

    # pushfirst!(PyVector(pyimport("sys")."path"), @__DIR__) #comment this line if no need to load local python module
    ete = PyCall.pyimport("ete3")

    function info(node)
        s = "N$(node.id):"
        if isempty(node.extT) == false
            s *= "T$(Tuple(node.extT)), "
        end

        if node.operation == 1
            if (node.factor ≈ 1.0) == false
                s *= @sprintf("%3.1f", node.factor)
            end
            s *= "× "
        elseif node.operation == 2
            @assert node.factor ≈ 1.0
            s *= "+ "
        else
            error("not implemented!")
        end
        return s
    end


    function treeview(node, level, t = nothing)
        if isnothing(t)
            t = ete.Tree(name = " ")
        end

        nt = t.add_child(name = info(node))
        name_face = ete.TextFace(nt.name, fgcolor = "black", fsize = 10)
        nt.add_face(name_face, column = 0, position = "branch-top")

        for child in node.nodes
            if child != -1
                treeview(tree[child], level + 1, nt)
            else
                nnt = nt.add_child(name = "0")
            end
        end
        for pidx in node.propagators
            if pidx != -1
                p = diag.propagators[pidx]
                factor = @sprintf("%3.1f", p.factor)
                nnt = nt.add_child(name = "P$pidx: typ $(p.type), K$(K[p.Kidx].basis), T$(p.Tidx), $factor")
            else
                nnt = nt.add_child(name = "0")
            end

            # name_face = ete.TextFace(nnt.name, fgcolor = "black", fsize = 10)
            # nnt.add_face(name_face, column = 0, position = "branch-top")
        end

        return t
    end

    t = treeview(root, 1)
    # style = ete.NodeStyle()
    # style["bgcolor"] = "Khaki"
    # t.set_style(style)


    ts = ete.TreeStyle()
    ts.show_leaf_name = true
    # ts.show_leaf_name = True
    # ts.layout_fn = my_layout
    ####### show tree vertically ############
    # ts.rotation = 90 #show tree vertically

    ####### show tree in an arc  #############
    # ts.mode = "c"
    # ts.arc_start = -180
    # ts.arc_span = 180
    # t.write(outfile="/home/kun/test.txt", format=8)
    t.show(tree_style = ts)
end

include("treeEval.jl")
end