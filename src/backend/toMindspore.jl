# ms = pyimport("mindspore")

"""
    function to_pystatic(operator::Type, subgraphs::AbstractVector{<:AbstractGraph}, subgraph_factors::AbstractVector)

Returns the static representation of a computational graph node `g` with operator `operator`, subgraphs `subgraphs`, and subgraph factors `subgraph_factors` in python.
"""
function to_pystatic(operator::Type, subgraphs::AbstractVector{<:AbstractGraph}, subgraph_factors::AbstractVector)
    error(
        "Static representation for computational graph nodes with operator $(operator) not yet implemented! " *
        "Please define a method `to_static(::Type{$(operator)}, subgraphs::$(typeof(subgraphs)), subgraph_factors::$(typeof(subgraph_factors)))`."
    )
end

function to_pystatic(::Type{ComputationalGraphs.Sum}, subgraphs::Vector{Graph{F,W}}, subgraph_factors::Vector{F}) where {F,W}
    if length(subgraphs) == 1
        factor_str = subgraph_factors[1] == 1 ? "" : " * $(subgraph_factors[1])"
        return "(g$(subgraphs[1].id)$factor_str)"
    else
        terms = ["g$(g.id)" * (gfactor == 1 ? "" : " * $gfactor") for (g, gfactor) in zip(subgraphs, subgraph_factors)]
        return "(" * join(terms, " + ") * ")"
    end
end

function to_pystatic(::Type{ComputationalGraphs.Prod}, subgraphs::Vector{Graph{F,W}}, subgraph_factors::Vector{F}) where {F,W}
    if length(subgraphs) == 1
        factor_str = subgraph_factors[1] == 1 ? "" : " * $(subgraph_factors[1])"
        return "(g$(subgraphs[1].id)$factor_str)"
    else
        terms = ["g$(g.id)" * (gfactor == 1 ? "" : " * $gfactor") for (g, gfactor) in zip(subgraphs, subgraph_factors)]
        return "(" * join(terms, " * ") * ")"
        # return "(" * join(["g$(g.id)" for g in subgraphs], " * ") * ")"
    end
end

function to_pystatic(::Type{ComputationalGraphs.Power{N}}, subgraphs::Vector{Graph{F,W}}, subgraph_factors::Vector{F}) where {N,F,W}
    factor_str = subgraph_factors[1] == 1 ? "" : " * $(subgraph_factors[1])"
    return "((g$(subgraphs[1].id))**$N$factor_str)"
end

function to_pystatic(::Type{ComputationalGraphs.Sum}, subgraphs::Vector{FeynmanGraph{F,W}}, subgraph_factors::Vector{F}) where {F,W}
    if length(subgraphs) == 1
        factor_str = subgraph_factors[1] == 1 ? "" : " * $(subgraph_factors[1])"
        return "(g$(subgraphs[1].id)$factor_str)"
    else
        terms = ["g$(g.id)" * (gfactor == 1 ? "" : " * $gfactor") for (g, gfactor) in zip(subgraphs, subgraph_factors)]
        return "(" * join(terms, " + ") * ")"
    end
end

function to_pystatic(::Type{ComputationalGraphs.Prod}, subgraphs::Vector{FeynmanGraph{F,W}}, subgraph_factors::Vector{F}) where {F,W}
    if length(subgraphs) == 1
        factor_str = subgraph_factors[1] == 1 ? "" : " * $(subgraph_factors[1])"
        return "(g$(subgraphs[1].id)$factor_str)"
    else
        terms = ["g$(g.id)" * (gfactor == 1 ? "" : " * $gfactor") for (g, gfactor) in zip(subgraphs, subgraph_factors)]
        return "(" * join(terms, " * ") * ")"
    end
end

function to_pystatic(::Type{ComputationalGraphs.Power{N}}, subgraphs::Vector{FeynmanGraph{F,W}}, subgraph_factors::Vector{F}) where {N,F,W}
    factor_str = subgraph_factors[1] == 1 ? "" : " * $(subgraph_factors[1])"
    return "((g$(subgraphs[1].id))**$N$factor_str)"
end

"""
    function to_julia_str(graphs::AbstractVector{<:AbstractGraph})
    
Compile a list of graphs into a string for a python static function and output a python script which support the static graph representation in mindspore framework.
"""

function to_python_str_ms(graphs::AbstractVector{<:AbstractGraph})
    head = "import mindspore as ms\n@ms.jit\n"
    head *= "def graphfunc(leaf):\n"
    body = "    graph_list = []\n"
    leafidx = 1
    root = [id(g) for g in graphs]
    inds_visitedleaf = Int[]
    inds_visitednode = Int[]
    for graph in graphs
        for g in PostOrderDFS(graph) #leaf first search
            g_id = id(g)
            target = "g$(g_id)"
            isroot = false
            if g_id in root
                isroot = true
            end
            if isempty(subgraphs(g)) #leaf
                g_id in inds_visitedleaf && continue
                factor_str = factor(g) == 1 ? "" : " * $(factor(g))"
                body *= "    $target = ms.Tensor(leaf[$(leafidx-1)])$factor_str\n"
                leafidx += 1
                push!(inds_visitedleaf, g_id)
            else
                g_id in inds_visitednode && continue
                factor_str = factor(g) == 1 ? "" : " * $(factor(g))"
                body *= "    $target = $(to_pystatic(operator(g), subgraphs(g), subgraph_factors(g)))$factor_str\n"
                push!(inds_visitednode, g_id)
            end
            if isroot
                body *= "    graph_list.append($target)\n"
            end
        end
    end
    tail = "    return graph_list\n"
    tail*= "def to_StaticGraph(leaf)\n"
    tail*= "    output = graphfunc(leaf)\n"
    tail*= "    return output"
    expr = head * body * tail
    println(expr)
    # return head * body * tail
    f = open("GraphFunc.py", "w")
    write(f, expr)
    return expr
end

# function to_mindspore_graph(graphs::AbstractVector{<:AbstractGraph})
#     pyexpr = to_python_str_ms(graphs)
#     py"""
#         import mindspore as ms
#         exec($pyexpr)
#         ms_graph = jit(fn=graphfunc)
#         out = ms_graph()
#     """
#     return py"out"
# end
