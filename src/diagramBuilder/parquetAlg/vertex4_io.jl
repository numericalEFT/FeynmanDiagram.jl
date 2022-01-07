################## implement AbstractTrees interface #######################
# refer to https://github.com/JuliaCollections/AbstractTrees.jl for more details
function AbstractTrees.children(ver4::Ver4)
    return ver4.bubble
end

function AbstractTrees.children(bubble::Bubble)
    return (bubble.Lver, bubble.Rver)
end

function iterate(ver4::Ver4{W}) where {W}
    if length(ver4.bubble) == 0
        return nothing
    else
        return (ver4.bubble[1], 1)
    end
end

function iterate(bub::Bubble)
    return (bub.Lver, false)
end

function iterate(ver4::Ver4{W}, state) where {W}
    if state >= length(ver4.bubble) || length(ver4.bubble) == 0
        return nothing
    else
        return (ver4.bubble[state+1], state + 1)
    end
end

function iterate(bub::Bubble, state::Bool)
    state && return nothing
    return (bub.Rver, true)
end

Base.IteratorSize(::Type{Ver4{W}}) where {W} = Base.SizeUnknown()
Base.eltype(::Type{Ver4{W}}) where {W} = Ver4{W}

Base.IteratorSize(::Type{Bubble{Ver4{W}}}) where {W} = Base.SizeUnknown()
Base.eltype(::Type{Bubble{Ver4{W}}}) where {W} = Bubble{Ver4{W}}

AbstractTrees.printnode(io::IO, ver4::Ver4) = print(io, tpair(ver4))
AbstractTrees.printnode(io::IO, bub::Bubble) = print(io,
    "\u001b[32m$(bub.id): $(bub.chan) $(bub.Lver.para.innerLoopNum)Ⓧ $(bub.Rver.para.innerLoopNum)\u001b[0m")

function tpair(ver4, MaxT = 18)
    s = "\u001b[31m$(ver4.id):\u001b[0m"
    if ver4.para.innerLoopNum > 0
        s *= "$(ver4.para.innerLoopNum)lp, T$(length(ver4.Tpair))⨁ "
    else
        s *= "⨁ "
    end
    # if ver4.loopNum <= 1
    for (ti, T) in enumerate(ver4.Tpair)
        if ti <= MaxT
            s *= "($(T[1]),$(T[2]),$(T[3]),$(T[4]))"
        else
            s *= "..."
            break
        end
    end
    # end
    return s
end

##### pretty print of Bubble and Ver4  ##########################
Base.show(io::IO, bub::Bubble) = AbstractTrees.printnode(io::IO, bub)
Base.show(io::IO, ver4::Ver4) = AbstractTrees.printnode(io::IO, ver4)


"""
convert Ver4 tree struct to a string in the newick format
"""
function newick(ver4::Ver4)
    return "$(newickVer4(ver4));"
end

"""
convert Bubble tree struct to a string in the newick format
"""
function newick(bub::Bubble)
    return "$(newickBuble(bub));"
end

function newickBubble(bub::Bubble)
    # Recursive version
    # Practically a postorder tree traversal
    left = newickVer4(bub.Lver)
    right = newickVer4(bub.Rver)
    return "($left,$right)$(bub.id)_$(ChanName[bub.chan])_$(bub.Lver.para.innerLoopNum)Ⓧ$(bub.Rver.para.innerLoopNum)"
end


function newickVer4(ver4::Ver4)
    # Recursive version
    # Practically a postorder tree traversal

    function tpairNewick(ver4)
        if ver4.para.innerLoopNum > 0
            s = "$(ver4.id):lp$(ver4.para.innerLoopNum)_T$(length(ver4.Tpair))⨁"
        else
            s = "$(Ver4.id):⨁"
        end
        # if ver4.loopNum <= 1
        for (ti, T) in enumerate(ver4.Tpair)
            if ti <= 5
                s *= "⟨$(T[1])-$(T[2])-$(T[3])-$(T[4])⟩"
            else
                s *= "…"
                break
            end
        end
        # end
        return s
    end

    if ver4.para.innerLoopNum == 0
        return tpairNewick(ver4)
    else
        s = "("
        for (bi, bub) in enumerate(ver4.bubble)
            s *= newickBubble(bub)
            if bi != length(ver4.bubble)
                s *= ","
            else
                s *= ")"
            end
        end
        return s * tpairNewick(ver4)
    end
end

"""
    showTree(ver4, para::Para; verbose=0, depth=999)

    Visualize the diagram tree using ete3 python package

#Arguments
- `ver4`: the 4-vertex diagram tree to visualize
- `para`: parameters
- `verbose=0`: the amount of information to show
- `depth=999`: deepest level of the diagram tree to show
"""
function showTree(ver4; verbose = 0, depth = 999)

    # pushfirst!(PyVector(pyimport("sys")."path"), @__DIR__) #comment this line if no need to load local python module
    ete = pyimport("ete3")

    function tpairETE(ver4, depth)
        s = "$(ver4.id):"
        if ver4.para.innerLoopNum > 0
            s *= "$(ver4.para.innerLoopNum)lp, T$(length(ver4.Tpair))⨁ "
        else
            s *= "⨁ "
        end
        if ver4.para.innerLoopNum == 0 || ver4.level > depth
            MaxT = Inf
        else
            MaxT = 1
        end

        # if ver4.loopNum <= 1
        for (ti, T) in enumerate(ver4.Tpair)
            if ti <= MaxT
                s *= "($(T[1]),$(T[2]),$(T[3]),$(T[4]))"
            else
                s *= "..."
                break
            end
        end
        # end
        return s
    end


    function treeview(ver4, t = nothing)
        if isnothing(t)
            t = ete.Tree(name = " ")
        end

        if ver4.para.innerLoopNum == 0 || ver4.level > depth
            nt = t.add_child(name = tpairETE(ver4, depth))
            return t
        else
            # prefix = "$(ver4.id): $(ver4.loopNum) lp, $(length(ver4.Tpair)) elem"
            # nt = t.add_child(name=prefix * ", ⨁")
            nt = t.add_child(name = tpairETE(ver4, depth))
            name_face = ete.TextFace(nt.name, fgcolor = "black", fsize = 10)
            nt.add_face(name_face, column = 0, position = "branch-top")
        end

        for bub in ver4.bubble
            nnt = nt.add_child(name = "$(bub.id): $(bub.chan) $(bub.Lver.para.innerLoopNum)Ⓧ$(bub.Rver.para.innerLoopNum)")

            name_face = ete.TextFace(nnt.name, fgcolor = "black", fsize = 10)
            nnt.add_face(name_face, column = 0, position = "branch-top")

            treeview(bub.Lver, nnt)
            treeview(bub.Rver, nnt)
        end

        return t
    end

    t = treeview(ver4)
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