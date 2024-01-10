function to_dotstatic(operator::Type, id::Int, subgraphs::AbstractVector{<:AbstractGraph}, subgraph_factors::AbstractVector)
    error(
        "Static representation for computational graph nodes with operator $(operator) not yet implemented! "
    )
end

function to_dotstatic(::Type{ComputationalGraphs.Sum}, id::Int, subgraphs::Vector{Graph{F,W}}, subgraph_factors::Vector{F}) where {F,W}
    node_temp = ""
    arrow_temp = ""
    # opr_node = "g$(id)[shape=box, label = <($factor)*&oplus;>, style=filled, fillcolor=cyan,fontsize=18]"
    opr_node = "g$(id)[shape=box, label = <&oplus;>, style=filled, fillcolor=cyan,fontsize=18]"
    node_temp *= opr_node
    for (g, gfactor) in zip(subgraphs, subgraph_factors)
        if gfactor != 1
            arrow_temp *= "g$(g.id)->g$(id)[arrowhead=vee,label=$gfactor,fontsize=16]\n"
        else
            arrow_temp *= "g$(g.id)->g$(id)[arrowhead=vee,]\n"
        end
    end
    return node_temp, arrow_temp
end

function to_dotstatic(::Type{ComputationalGraphs.Prod}, id::Int, subgraphs::Vector{Graph{F,W}}, subgraph_factors::Vector{F}) where {F,W}
    node_temp = ""
    arrow_temp = ""
    # opr_node = "g$id[shape=box, label = <($factor)&otimes;>, style=filled, fillcolor=cornsilk,fontsize=18]\n"
    opr_node = "g$id[shape=box, label = <&otimes;>, style=filled, fillcolor=cornsilk,fontsize=18]\n"
    node_temp *= opr_node
    for (g, gfactor) in zip(subgraphs, subgraph_factors)
        if gfactor != 1
            arrow_temp *= "g$(g.id)->g$(id)[arrowhead=vee,label=$gfactor,fontsize=16]\n"
        else
            arrow_temp *= "g$(g.id)->g$(id)[arrowhead=vee,]\n"
        end
    end
    # end
    return node_temp, arrow_temp
end

function to_dotstatic(::Type{ComputationalGraphs.Power{N}}, id::Int, subgraphs::Vector{Graph{F,W}}, subgraph_factors::Vector{F}) where {N,F,W}
    node_temp = ""
    arrow_temp = ""
    # opr_node = "g$id[shape=box, label = <($factor)*Pow($N)>, style=filled, fillcolor=darkolivegreen,fontsize=18]\n"
    opr_node = "g$id[shape=box, label = <Pow($N)>, style=filled, fillcolor=darkolivegreen,fontsize=18]\n"
    node_temp *= opr_node
    if subgraph_factors[1] != 1
        arrow_temp *= "g$(subgraphs[1].id)->$opr_name[arrowhead=vee,label=$gfactor,fontsize=16]\n"
    else
        arrow_temp *= "g$(subgraphs[1].id)->$opr_name[arrowhead=vee,]\n"
    end
    return node_temp, arrow_temp
end

function to_dotstatic(::Type{ComputationalGraphs.Sum}, id::Int, subgraphs::Vector{FeynmanGraph{F,W}}, subgraph_factors::Vector{F}) where {F,W}
    node_temp = ""
    arrow_temp = ""
    # opr_node = "g$(id)[shape=box, label = <($factor)*&oplus;>, style=filled, fillcolor=cyan,fontsize=18]"
    opr_node = "g$(id)[shape=box, label = <&oplus;>, style=filled, fillcolor=cyan,fontsize=18]"
    opr_name = "g$id"
    node_temp *= opr_node
    for (g, gfactor) in zip(subgraphs, subgraph_factors)
        if gfactor != 1
            arrow_temp *= "g$(g.id)->$opr_name[arrowhead=vee,label=$gfactor,fontsize=16]\n"
        else
            arrow_temp *= "g$(g.id)->$opr_name[arrowhead=vee,]\n"
        end
    end
    return node_temp, arrow_temp
end

function to_dotstatic(::Type{ComputationalGraphs.Prod}, id::Int, subgraphs::Vector{FeynmanGraph{F,W}}, subgraph_factors::Vector{F}) where {F,W}
    node_temp = ""
    arrow_temp = ""
    # opr_node = "g$id[shape=box, label = <($factor)&otimes;>, style=filled, fillcolor=cornsilk,fontsize=18]\n"
    opr_node = "g$id[shape=box, label = <&otimes;>, style=filled, fillcolor=cornsilk,fontsize=18]\n"
    node_temp *= opr_node
    for (g, gfactor) in zip(subgraphs, subgraph_factors)
        if gfactor != 1
            arrow_temp *= "g$(g.id)->g$(id)[arrowhead=vee,label=$gfactor,fontsize=16]\n"
        else
            arrow_temp *= "g$(g.id)->g$(id)[arrowhead=vee,]\n"
        end
    end
    # end
    return node_temp, arrow_temp
end

function to_dotstatic(::Type{ComputationalGraphs.Power{N}}, id::Int, subgraphs::Vector{FeynmanGraph{F,W}}, subgraph_factors::Vector{F}) where {N,F,W}
    node_temp = ""
    arrow_temp = ""
    # opr_node = "g$id[shape=box, label = <($factor)*Pow($N)>, style=filled, fillcolor=darkolivegreen,fontsize=18]\n"
    opr_node = "g$id[shape=box, label = <Pow($N)>, style=filled, fillcolor=darkolivegreen,fontsize=18]\n"
    node_temp *= opr_node
    if subgraph_factors[1] != 1
        arrow_temp *= "g$(subgraphs[1].id)->$opr_name[arrowhead=vee,label=$(subgraph_factors[1]),fontsize=16]\n"
    else
        arrow_temp *= "g$(subgraphs[1].id)->$opr_name[arrowhead=vee,]\n"
    end
    return node_temp, arrow_temp
end

"""
    function to_dot_str(graphs::AbstractVector{<:AbstractGraph}, name::String="")

    Compile a list of graphs into a string for dot language.

    # Arguments:
    - `graphs`  vector of computational graphs
    - `title`   The name of the compiled function (defaults to nothing)
"""
function to_dot_str(graphs::AbstractVector{<:AbstractGraph}, name::String="")
    head = "digraph ComputationalGraph { \nlabel=\"$name\"\n"
    head *= "ReturnNode[shape=box, label = \"Return\", style=filled, fillcolor=darkorange,fontsize=18]\n"
    body_node = ""
    body_arrow = ""
    leafidx = 1
    root = [id(g) for g in graphs]
    inds_visitedleaf = Int[]
    inds_visitednode = Int[]
    rootidx = 1
    for graph in graphs
        for g in PostOrderDFS(graph) #leaf first search
            g_id = id(g)
            isroot = false
            if g_id in root
                isroot = true
            end
            if isempty(subgraphs(g)) #leaf
                g_id in inds_visitedleaf && continue
                leafname = get_leafname(g, leafidx)
                gnode_str = "g$g_id[label=<$leafname>, style=filled, fillcolor=paleturquoise,fontsize=18]\n"
                body_node *= gnode_str
                leafidx += 1
                push!(inds_visitedleaf, g_id)
            else
                g_id in inds_visitednode && continue
                temp_node, temp_arrow = to_dotstatic(operator(g), g_id, subgraphs(g), subgraph_factors(g))
                body_node *= temp_node
                body_arrow *= temp_arrow
                push!(inds_visitednode, g_id)
            end
            if isroot
                body_arrow *= "g$(g_id)->ReturnNode[arrowhead=vee,]\n"
                rootidx += 1
            end
        end
    end
    tail = "   }\n"
    expr = head * body_node * body_arrow * tail
    # println(expr)
    return expr
end

function compile_dot(graphs::AbstractVector{<:AbstractGraph}, filename::String; graph_name="")
    dot_string = to_dot_str(graphs, graph_name)
    open(filename, "w") do f
        write(f, dot_string)
    end
end

function get_leafname(g, leafidx)
    leaftype = Nothing
    if g isa FeynmanGraph
        leaftype = g.properties.diagtype
    elseif g isa Graph
        leaftype = typeof(g.properties)
    else
        error("Unknown graph type: $(typeof(g))")
    end

    if leaftype == BareGreenId
        leafname = "<i>G</i><sub>$leafidx</sub>"
        println(leaftype, ": ", leafidx, " ", g.properties.extK, " ", g.properties.extT)
    elseif leaftype == BareInteractionId
        println(leaftype, ": ", leafidx, " ", g.properties.response, " ", g.properties.type, " ", g.properties.extK, " ", g.properties.extT)
        leafname = "<i>V</i><sub>$leafidx</sub>"
    elseif leaftype == PolarId
        leafname = "&Pi;<sub>$leafidx</sub>"
    elseif leaftype == Ver3Id
        leafname = "&Gamma;<sup>(3)</sup><sub>$leafidx</sub>"
    elseif leaftype == Ver4Id
        leafname = "&Gamma;<sup>(4)</sup><sub>$leafidx</sub>"
    elseif leaftype == ComputationalGraphs.Propagator
        if isfermionic(g.properties.vertices[1])
            leafname = "<i>G</i><sub>$leafidx</sub>"
        else
            leafname = "<i>V</i><sub>$leafidx</sub>"
        end
    elseif leaftype == ComputationalGraphs.Interaction
        leafname = "<i>Ver</i><sub>$leafidx</sub>"
    else
        println("Unknown leaf type: $leaftype")
        leafname = "L<sub>$leafidx</sub>"
    end
    return leafname
end
