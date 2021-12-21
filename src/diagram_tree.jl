module DiagTree
using Printf, PyCall
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
    function Propagator{P,F}(order, variable = [], factor = 1.0, para = 0) where {P,F}
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
    # root::SubArray{CachedObject{NODE,W},1,Vector{CachedObject{NODE,W}},Tuple{Vector{Int64}},false}
    #SubArray has 5 type parameters. The first two are the standard element type and dimensionality. 
    #The next is the type of the parent AbstractArray. The most heavily-used is the fourth parameter, 
    #a Tuple of the types of the indices for each dimension. The final one, L, is only provided as a convenience for dispatch;
    #it's a boolean that represents whether the index types support fast linear indexing.
    # loopNum::Int
    # tauNum::Int
    function Diagrams{PARA,F,W}(basis::V, propagator::P) where {V,P,PARA,F,W}
        return new{V,P,PARA,F,W}(basis, propagator, Pool{Node{PARA,F},W}(), [])
    end
    function Diagrams{F,W}(basis::V, propagator::P) where {V,P,F,W}
        PARA = Int
        return new{V,P,PARA,F,W}(basis, propagator, Pool{Node{PARA,F},W}(), [])
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

# function Node{F, P}(operation::Int, components = [[]], child = [], parent = 0, factor = 1.0, para = 0) where {F,P}
# function addNode(diag::Diagrams, operator, components::Vector{AbstractVector}, childNodes::AbstractVector; factor = 1.0, parent = 0, para = 0, currWeight = 0)
function addNode(diag::Diagrams, operator, components, childNodes; factor = 1.0, parent = 0, para = 0, currWeight = 0.0)
    # println("components: ", components)
    nodePool = diag.nodePool
    @assert length(components) == length(diag.propagatorPool) "each element of the components is an index vector of the corresponding propagator"

    _NodePool = typeof(nodePool)
    _CachedNode = eltype(fieldtype(_NodePool, :pool))
    _Node = fieldtype(_CachedNode, :object)
    PARA = fieldtype(_Node, :para)
    F = fieldtype(_Node, :factor)

    node = Node{PARA,F}(operator, components, childNodes, factor, parent, para)

    # println("para: ", PARA)
    # println("F: ", F)
    # println("node: ", node)
    nidx = append(nodePool, node, currWeight)
    return nidx
end

"""
    showTree(diag::Diagrams, _root = diag.root[end]; verbose = 0, depth = 999)

    Visualize the diagram tree using ete3 python package

#Arguments
- `ver4`: the 4-vertex diagram tree to visualize
- `para`: parameters
- `verbose=0`: the amount of information to show
- `depth=999`: deepest level of the diagram tree to show
"""
function showTree(diag::Diagrams, _root = diag.root[end]; verbose = 0, depth = 999)
    # println(diag)
    # println(diag.nodePool)
    tree = diag.nodePool
    basisPool = diag.basisPool
    root = tree[_root]

    # pushfirst!(PyVector(pyimport("sys")."path"), @__DIR__) #comment this line if no need to load local python module
    ete = PyCall.pyimport("ete3")

    function info(cachedNode)
        s = "N$(cachedNode.id):"
        node = cachedNode.object
        s *= sprint(show, node)
        # if isempty(node.extT) == false
        #     s *= "T$(Tuple(node.extT)), "
        # end

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

        # println("node: ", cachedNode.object.components)

        for (ci, component) in enumerate(cachedNode.object.components)
            # println(component, ", ", ci)
            propagatorPool = diag.propagatorPool[ci]
            for pidx in component
                p = propagatorPool[pidx].object #Propagator
                factor = @sprintf("%3.1f", p.factor)
                # nnt = nt.add_child(name = "P$pidx: typ $ci, K$(K[p.Kidx].basis), T$(p.Tidx), $factor")
                nnt = nt.add_child(name = "P$pidx: typ $ci, var $(p.variable) $factor")
            end
        end

        # for pidx in cacheNode.propagators
        #     if pidx != -1
        #         p = diag.propagators[pidx]
        #         factor = @sprintf("%3.1f", p.factor)
        #         nnt = nt.add_child(name = "P$pidx: typ $(p.type), K$(K[p.Kidx].basis), T$(p.Tidx), $factor")
        #     else
        #         nnt = nt.add_child(name = "0")
        #     end

        #     # name_face = ete.TextFace(nnt.name, fgcolor = "black", fsize = 10)
        #     # nnt.add_face(name_face, column = 0, position = "branch-top")
        # end

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