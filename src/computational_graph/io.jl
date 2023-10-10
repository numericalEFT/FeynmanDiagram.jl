function _ops_to_str(ops::Vector{OperatorProduct})
    strs = ["$(o)" for o in ops]
    return join(strs, "|")
end

function short(factor, ignore=nothing)
    if isnothing(ignore) == false && applicable(isapprox, factor, ignore) && factor ≈ ignore
        return ""
    end
    str = "$(factor)"
    if factor isa Float64
        return length(str) <= 4 ? str : @sprintf("%6.3e", factor)
    elseif factor isa Vector{Float64}
        return length(str) <= 4 ? str : reduce(*, [@sprintf("%6.3e", f) for f in factor])
    else
        return str
    end
end

function short_orders(orders)
    orders_no_trailing_zeros = ""
    idx_last_set = findlast(x -> x != 0, orders)
    if isnothing(idx_last_set) == false
        orders_no_trailing_zeros *= string(orders[1:idx_last_set])
    end
    return orders_no_trailing_zeros
end

function _stringrep(graph::AbstractGraph, color=true)
    namestr = isempty(graph.name) ? "" : "-$(graph.name)"
    idstr = "$(graph.id)$namestr"
    if graph isa FeynmanGraph
        idstr *= ":$(_ops_to_str(vertices(graph)))"
    end
    fstr = short(graph.factor, one(graph.factor))
    wstr = short(graph.weight)
    ostr = short_orders(orders(graph))
    # =$(node.weight*(2π)^(3*node.id.para.innerLoopNum))

    if length(graph.subgraphs) == 0
        return isempty(fstr) ? "$(idstr)$(ostr)=$wstr" : "$(idstr)⋅$(fstr)=$wstr"
    else
        return "$(idstr)$(ostr)=$wstr=$(fstr)$(graph.operator) "
    end
end

"""
    show(io::IO, graph::AbstractGraph; kwargs...)

    Write a text representation of an AbstractGraph `graph` to the output stream `io`.
"""
function Base.show(io::IO, graph::AbstractGraph; kwargs...)
    if length(graph.subgraphs) == 0
        typestr = ""
    else
        typestr = join(["$(g.id)" for g in graph.subgraphs], ",")
        typestr = "($typestr)"
    end
    print(io, "$(_stringrep(graph, true))$typestr")
end
Base.show(io::IO, ::MIME"text/plain", graph::AbstractGraph; kwargs...) = Base.show(io, graph; kwargs...)

"""
    function plot_tree(graph::AbstractGraph; verbose = 0, maxdepth = 6)

    Visualize the computational graph as a tree using ete3 python package

#Arguments
- `graph::AbstractGraph`        : the computational graph struct to visualize
- `verbose=0`   : the amount of information to show
- `maxdepth=6`  : deepest level of the computational graph to show
"""
function plot_tree(graph::AbstractGraph; verbose=0, maxdepth=6)

    # pushfirst!(PyVector(pyimport("sys")."path"), @__DIR__) #comment this line if no need to load local python module
    ete = PyCall.pyimport("ete3")

    function treeview(node, level, t=ete.Tree(name=" "))
        if level > maxdepth
            return
        end
        name = "$(_stringrep(node, false))"
        nt = t.add_child(name=name)

        if length(node.subgraphs) > 0
            name_face = ete.TextFace(nt.name, fgcolor="black", fsize=10)
            nt.add_face(name_face, column=0, position="branch-top")
            for child in node.subgraphs
                treeview(child, level + 1, nt)
            end
        end

        return t
    end

    t = treeview(graph, 1)

    # NOTE: t.set_style does not update the original PyObject as expected, i.e.,
    #       `t.set_style(ete.NodeStyle(bgcolor="Khaki"))` does not modify t.
    #
    # The low-level approach circumvents this by directly updating the original PyObject `t."img_style"`
    PyCall.set!(t."img_style", "bgcolor", "Khaki")

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
function plot_tree(graphs::Vector{<:AbstractGraph}; kwargs...)
    for graph in graphs
        plot_tree(graph; kwargs...)
    end
end
