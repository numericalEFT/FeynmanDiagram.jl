function printBasisPool(diag::Diagrams, io = Base.stdout)
    printstyled(io, "Loop Basis ($(length(diag.basisPool)) in total)\n", color = :blue)
    title = @sprintf("%5s%40s\n", "index", "loopBasis")
    printstyled(io, title, color = :green)
    for i = 1:length(diag.basisPool)
        b = diag.basisPool.basis[:, i]
        @printf(io, "%5i%40s\n", i, "$b")
    end
    println(io)
end

function printPropagator(diag::Diagrams, io = Base.stdout)
    for pool in diag.propagatorPool
        printstyled(io, "Propagator $(pool.name) ($(length(pool)) in total)\n", color = :blue)
        title = @sprintf("%5s%40s%40s%40s\n", "index", "para", "loop", "site")
        printstyled(io, title, color = :green)
        for (idx, p) in enumerate(pool.object)
            site = isempty(p.siteBasis) ? "" : "$(p.siteBasis)"
            loop = p.loopIdx <= 0 ? "" : "$(diag.basisPool[p.loopIdx])"
            # loop = p.loopIdx <= 0 ? "" : "$(p.loopIdx)"
            @printf(io, "%5i%40s%40s%40s\n", idx, "$(p.para)", loop, site)
        end
        println(io)
    end
    println(io)
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
function showTree(diag::Diagrams, _root::Int; verbose = 0, depth = 999)
    tree = diag.nodePool
    basisPool = diag.basisPool
    root = tree[_root]
    id = _root

    printBasisPool(diag)
    printPropagator(diag)

    # pushfirst!(PyVector(pyimport("sys")."path"), @__DIR__) #comment this line if no need to load local python module
    ete = PyCall.pyimport("ete3")

    function factor(f)
        if f ≈ 1
            return ""
        end
        s = "$f"
        if length(s) <= 4
            return s
        else
            return @sprintf("%6.3e", f)
        end
    end

    function info(node, id)
        s = "N$(id) $(node.para):"
        # s *= sprint(show, node.para)
        # s *= ": "

        s *= factor(node.factor)

        if node.operation == MUL
            s *= ", x"
        elseif node.operation == ADD
            s *= ", +"
        else
            error("not implemented!")
        end
        return s
    end


    function treeview(node, id, level, t = nothing)
        if isnothing(t)
            t = ete.Tree(name = " ")
        end

        nt = t.add_child(name = info(node, id))
        name_face = ete.TextFace(nt.name, fgcolor = "black", fsize = 10)
        nt.add_face(name_face, column = 0, position = "branch-top")

        for child in node.childNodes
            if child != -1
                treeview(tree[child], child, level + 1, nt)
            else
                nnt = nt.add_child(name = "0")
            end
        end

        for (ci, component) in enumerate(node.components)
            # println(component, ", ", ci)
            propagatorPool = diag.propagatorPool[ci]
            for pidx in component
                p = propagatorPool.object[pidx] #Propagator
                site = isempty(p.siteBasis) ? "" : " site $(p.siteBasis),"
                # loop = p.loopIdx <= 0 ? "" : "loop $(diag.basisPool[p.loopIdx])"
                loop = p.loopIdx <= 0 ? "" : "$(p.loopIdx)"
                nnt = nt.add_child(name = "P$pidx $(p.para): $loop,$site $(factor(p.factor))")
            end
        end

        return t
    end

    t = treeview(root, id, 1)
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