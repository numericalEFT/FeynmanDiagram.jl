function printBasisPool(diag, io=Base.stdout)
    printstyled(io, "Loop Basis ($(length(diag.loopBasis)) in total)\n", color=:blue)
    title = @sprintf("%5s%40s\n", "index", "loopBasis")
    printstyled(io, title, color=:green)
    for i = 1:length(diag.loopBasis)
        b = diag.loopBasis.basis[:, i]
        @printf(io, "%5i%40s\n", i, "$b")
    end
    println(io)
end

function printNodes(diag, io=Base.stdout)
    printstyled(io, "Node ($(length(diag.node)) in total)\n", color=:blue)
    title = @sprintf("%5s%5s%40s%40s\n", "index", "name", "para", "child")
    printstyled(io, title, color=:green)
    for (idx, n) in enumerate(diag.node.object)
        # site = isempty(p.siteBasis) ? "" : "$(p.siteBasis)"
        # loop = p.loopIdx <= 0 ? "" : "$(diag.basisPool[p.loopIdx])"
        # loop = p.loopIdx <= 0 ? "" : "$(p.loopIdx)"
        if n isa Propagator
            site = isempty(n.siteBasis) ? "" : "$(n.siteBasis)"
            loop = n.loopIdx <= 0 ? "" : "$(diag.loopBasis[n.loopIdx])"
            @printf(io, "%5i%5s%40s%40s%40s\n", idx, "$(n.name)", "$(n.para)", loop, site)
        else
            @printf(io, "%5i%5s%40s%40s\n", idx, "$(n.name)", "$(n.para)", "$(n.childNodes)")
        end
    end
    println(io)
end

function Base.show(io::IO, tree::ExpressionTree)
    print(io, "ExprTree: $(tree.name) with root $(tree.root)")
    # print(io, "$(short(v.response))#$(v.order), k$(v.extK), t$(v.extT)")
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
function showTree(tree::ExpressionTree, _root::Int; verbose=0, depth=999)

    # pushfirst!(PyVector(pyimport("sys")."path"), @__DIR__) #comment this line if no need to load local python module
    ete = PyCall.pyimport("ete3")

    function name_para(p)
        name = (p.name == :none) ? "" : " $(p.name)"
        para = isnothing(p.para) ? "" : " $(p.para)"
        return name * para
    end

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

    function info(node, idx)
        s = "N$(idx)$(name_para(node)) = $(tree.node.current[idx])"
        # s *= sprint(show, node.para)
        # s *= ": "

        if node.operation == MUL
            s *= " = $(factor(node.factor))x"
        elseif node.operation == ADD
            s *= " = $(factor(node.factor))+"
        else
            error("not implemented!")
        end
        return s
    end


    function treeview(idx::Int, level, t=nothing)
        if isnothing(t)
            t = ete.Tree(name=" ")
        end

        if tree.node.object[idx] isa PropagatorId
            p = tree.node.object[idx]    #Propagator
            site = isempty(p.siteBasis) ? "" : " t$(p.siteBasis),"
            loop = p.loopIdx <= 0 ? "" : "k$(tree.loopBasis[p.loopIdx])"
            # loop = p.loopIdx <= 0 ? "" : "$(p.loopIdx)"
            nnt = t.add_child(name="P$(idx)$(name_para(p)): $loop,$site $(factor(p.factor))")
        else # composite node
            nt = t.add_child(name=info(tree.node.object[idx], idx))
            name_face = ete.TextFace(nt.name, fgcolor="black", fsize=10)
            nt.add_face(name_face, column=0, position="branch-top")

            for child in tree.node.object[idx].children
                if child != -1
                    treeview(child, level + 1, nt)
                else
                    nnt = nt.add_child(name="0")
                end
            end
        end

        return t
    end

    t = treeview(_root, 1)
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
    t.show(tree_style=ts)
end
# function showTree(diag::ExpressionTree, _root::Component; kwargs...)
#     @assert _root.isNode "Can not visualize $_root, because it is not a Node!"
#     return showTree(diag, _root.index; kwargs...)
# end