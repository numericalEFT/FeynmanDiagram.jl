function _summary(diag::Diagram, color = true)
    function factor()
        f = diag.factor
        if f isa Number && f ≈ 1
            return ""
        end
        s = "$(f)"
        if f isa Float64
            return length(s) <= 4 ? s : @sprintf("%6.3e", f)
        else
            return s
        end
    end
    function hash()
        return color ? "\u001b[32m#$(diag.hash)\u001b[0m" : "#$(diag.hash)"
    end
    if length(diag.subdiagram) == 0
        s = "$(hash()): $(diag.id)"
        f = factor()
        return isempty(f) ? "$s " : "$s, $f x "
    elseif length(diag.subdiagram) == 1
        return "$(hash()): $(diag.id) = $(factor())"
    else
        return "$(hash()): $(diag.id) = $(factor())$(diag.operator) "
    end
end

function Base.show(io::IO, diag::Diagram)
    if length(diag.subdiagram) == 0
        typestr = "bare"
    elseif length(diag.subdiagram) == 1
        typestr = "#$(diag.hash)"
    else
        subdiag = prod(["#$(d.hash), " for d in diag.subdiagram[1:end-1]])
        subdiag *= "#$(diag.subdiagram[end].hash)"
        typestr = "($subdiag)"
    end
    print(io, "$(_summary(diag, true))$typestr = $(diag.weight)")
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
function plot_tree(diag::Diagram; verbose = 0, maxdepth = 999)

    # pushfirst!(PyVector(pyimport("sys")."path"), @__DIR__) #comment this line if no need to load local python module
    ete = PyCall.pyimport("ete3")

    function treeview(node, level, t = ete.Tree(name = " "))
        nt = t.add_child(name = _summary(node, false))
        name_face = ete.TextFace(nt.name, fgcolor = "black", fsize = 10)
        nt.add_face(name_face, column = 0, position = "branch-top")
        for child in node.subdiagram
            treeview(child, level + 1, nt)
        end
        return t
    end

    t = treeview(diag, 1)
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