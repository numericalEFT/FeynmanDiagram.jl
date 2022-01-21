function toDict(diag::Diagram; maxdepth::Int)
    @assert maxdepth == 1 "deep convert has not yet been implemented!"

    d = Dict{Symbol,Any}()
    d[:hash] = diag.hash
    d[:id] = diag.id
    d[:name] = diag.name
    d[:diagram] = diag
    d[:subdiagram] = Tuple(d.hash for d in diag.subdiagram)
    d[:operator] = diag.operator
    d[:factor] = diag.factor
    d[:weight] = diag.weight

    return d
end

function toDict(v::DiagramId)
    d = Dict{Symbol,Any}()
    for field in fieldnames(typeof(v))
        data = getproperty(v, field)
        #DataFrame will expand a vector into multiple rows. To prevent it, we transform all vectors into tuples
        d[field] = data isa AbstractVector ? Tuple(data) : data
    end
    return d
end

# function toDataFrame(diagVec::AbstractVector; expand::Bool = false)

# end

function toDataFrame(diagVec::AbstractVector; expand::Bool = false, maxdepth::Int = 1)
    vec_of_diag_dict = [toDict(d, maxdepth = maxdepth) for d in diagVec]
    vec_of_id_dict = [toDict(d.id) for d in diagVec]
    names = Set(reduce(union, keys(d) for d in vec_of_diag_dict))
    if expand
        idnames = Set(reduce(union, keys(d) for d in vec_of_id_dict))
        @assert isempty(intersect(names, idnames)) "collision of diagram names $names and id names $idnames"
        names = union(names, idnames)

        for (di, d) in enumerate(vec_of_diag_dict)
            merge!(d, vec_of_id_dict[di]) #add id dict into the diagram dict
        end
    end
    # println(names)
    df = DataFrame([name => [] for name in names])

    for dict in vec_of_diag_dict
        append!(df, dict, cols = :union)
    end
    return df
end

function _summary(diag::Diagram{W}, color = true) where {W}
    function factor()
        factor = diag.factor
        if factor isa Number && factor ≈ one(W)
            return ""
        end
        fstr = "$(factor)"
        if factor isa Float64
            return length(fstr) <= 4 ? fstr : @sprintf("%6.3e", factor)
        else
            return fstr
        end
    end

    namestr = diag.name == :none ? "" : "$(diag.name) "
    idstr = "$namestr$(diag.id)"

    if length(diag.subdiagram) == 0
        f = factor()
        return isempty(f) ? "$idstr " : "$idstr, $f⋅"
    else
        return "$idstr=$(factor())$(diag.operator) "
    end
end

function Base.show(io::IO, diag::Diagram)
    if length(diag.subdiagram) == 0
        typestr = "bare"
    else
        subdiag = prod(["$(d.hash), " for d in diag.subdiagram[1:end-1]])
        subdiag *= "$(diag.subdiagram[end].hash)"
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
        nt = t.add_child(name = "$(node.hash): " * _summary(node, false))

        if length(node.subdiagram) > 0
            name_face = ete.TextFace(nt.name, fgcolor = "black", fsize = 10)
            nt.add_face(name_face, column = 0, position = "branch-top")
            for child in node.subdiagram
                treeview(child, level + 1, nt)
            end
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