module DiagTree
using Printf, PyCall
include("pool.jl")
using ..Var


struct Propagator{PARA,F}
    para::PARA
    order::Int
    factor::F
    basis::Vector{Int}
    function Propagator(order, basis = [], factor::F = 1.0, para::P = 0) where {F,P}
        return new{P,F}(para, order, factor, basis)
    end
    function Propagator{P,F}(order, basis = [], factor = 1.0, para = 0) where {P,F}
        return new{P,F}(para, order, factor, basis)
    end
end

function Base.isequal(a::Propagator{P,F}, b::Propagator{P,F}) where {P,F}
    if (isequal(a.para, b.para) == false) || (a.order != b.order) || (a.basis != b.basis)
        return false
    else
        return true
    end
end
Base.:(==)(a::Propagator{P,F}, b::Propagator{P,F}) where {P,F} = Base.isequal(a, b)

struct Node{PARA,F}
    para::PARA
    operation::Int #1: multiply, 2: add, ...
    factor::F
    components::Vector{Vector{Int}}
    childNodes::Vector{Int}
    parent::Int # parent id
    # child::SubArray{CachedObject{Node{PARA, P},W},1,Vector{CachedObject{NODE,W}},Tuple{Vector{Int64}},false}

    function Node(operation::Int, components = [[]], child = [], factor::F = 1.0, parent = 0, para::P = 0) where {F,P}
        return new{P,F}(para, operation, factor, components, child, parent)
    end
    function Node{P,F}(operation::Int, components = [[]], child = [], factor = 1.0, parent = 0, para = 0) where {F,P}
        return new{P,F}(para, operation, factor, components, child, parent)
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

mutable struct Diagrams{V,P,PARA,F,W}
    basisPool::V
    propagatorPool::P
    nodePool::Pool{Node{PARA,F},W}
    root::Vector{Int}
    function Diagrams{PARA,F,W}(basisPool::V, propagatorPool::P) where {V,P,PARA,F,W}
        return new{V,P,PARA,F,W}(basisPool, propagatorPool, Pool{Node{PARA,F},W}(), [])
    end
    function Diagrams{F,W}(basisPool::V, propagatorPool::P) where {V,P,F,W}
        PARA = Int
        return new{V,P,PARA,F,W}(basisPool, propagatorPool, Pool{Node{PARA,F},W}(), [])
    end
end

function addPropagator(diag::Diagrams, index::Int, order::Int, basis::AbstractVector, factor = 1, para = 0, currWeight = 0)
    basisPool = diag.basisPool
    propagatorPool = diag.propagatorPool
    # @assert length(basis) == length(variablePool) == length(currVar) "$(length(basis)) == $(length(variablePool)) == $(length(currVar)) breaks"

    PROPAGATOR_POOL = typeof(propagatorPool[index])
    CACHEDPROPAGATOR = eltype(fieldtype(PROPAGATOR_POOL, :pool))
    PROPAGATOR = fieldtype(CACHEDPROPAGATOR, :object)
    PARA = fieldtype(PROPAGATOR, :para)
    F = fieldtype(PROPAGATOR, :factor)

    vidx = zeros(length(basis))
    for (bi, b) in enumerate(basis)
        # b[1]: basis, b[2]: initialize variable (curr)
        vidx[bi] = append(basisPool[bi], b[1], b[2])
    end
    prop = Propagator{PARA,F}(order, vidx, factor, para)
    return append(diag.propagatorPool[index], prop, currWeight)
end

function addNode(diag::Diagrams, operator, components, childNodes; factor = 1.0, parent = 0, para = 0, currWeight = 0.0)
    nodePool = diag.nodePool
    @assert length(components) == length(diag.propagatorPool) "each element of the components is an index vector of the corresponding propagator"

    _NodePool = typeof(nodePool)
    _CachedNode = eltype(fieldtype(_NodePool, :pool))
    _Node = fieldtype(_CachedNode, :object)
    PARA = fieldtype(_Node, :para)
    F = fieldtype(_Node, :factor)

    node = Node{PARA,F}(operator, components, childNodes, factor, parent, para)

    nidx = append(nodePool, node, currWeight)
    return nidx
end

function printBasisPool(diag::Diagrams)
    Nmax = maximum([length(b) for b in diag.basisPool])
    print("index  ")
    for (bi, b) in enumerate(diag.basisPool)
        print("basis#$bi  ")
    end
    print("\n")
    for i = 1:Nmax
        print("$i   ")
        for b in diag.basisPool
            if length(b) >= i
                cachedobj = b[i]
                print("$(cachedobj.object)  ")
            end
        end
        print("\n")
    end
end

"""
    showTree(diag::Diagrams, _root = diag.root[end]; verbose = 0, depth = 999)

    Visualize the diagram tree using ete3 python package

#Arguments
- `diag`: the Diagrams struct to visualize
- `_root`: the index of the root node to visualize
- `verbose=0`: the amount of information to show
- `depth=999`: deepest level of the diagram tree to show
"""
function showTree(diag::Diagrams, _root = diag.root[end]; verbose = 0, depth = 999)
    tree = diag.nodePool
    basisPool = diag.basisPool
    root = tree[_root]

    printBasisPool(diag)

    # pushfirst!(PyVector(pyimport("sys")."path"), @__DIR__) #comment this line if no need to load local python module
    ete = PyCall.pyimport("ete3")

    function info(cachedNode)
        s = "N$(cachedNode.id):"
        node = cachedNode.object
        s *= sprint(show, node.para)
        s *= ", "

        if node.operation == 1
            if (node.factor ≈ 1.0) == false
                s *= @sprintf("%3.1f", node.factor)
            end
            s *= "x "
        elseif node.operation == 2
            @assert node.factor ≈ 1.0
            s *= "+ "
        else
            error("not implemented!")
        end
        return s
    end


    function treeview(cachedNode, level, t = nothing)
        if isnothing(t)
            t = ete.Tree(name = " ")
        end

        nt = t.add_child(name = info(cachedNode))
        name_face = ete.TextFace(nt.name, fgcolor = "black", fsize = 10)
        nt.add_face(name_face, column = 0, position = "branch-top")

        for child in cachedNode.object.childNodes
            if child != -1
                treeview(tree[child], level + 1, nt)
            else
                nnt = nt.add_child(name = "0")
            end
        end

        for (ci, component) in enumerate(cachedNode.object.components)
            # println(component, ", ", ci)
            propagatorPool = diag.propagatorPool[ci]
            for pidx in component
                p = propagatorPool[pidx].object #Propagator
                factor = @sprintf("%3.1f", p.factor)
                # nnt = nt.add_child(name = "P$pidx: typ $ci, K$(K[p.Kidx].basis), T$(p.Tidx), $factor")
                nnt = nt.add_child(name = "P$pidx: typ $ci, basis $(p.basis), factor: $factor")
            end
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

end