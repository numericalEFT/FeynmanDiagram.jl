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
                # cachedobj = b[i]
                # print("$(cachedobj.object)  ")
                print(sprint(show, b[i]) * "   ")
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