"""
    function to_python_str(graphs::AbstractVector{<:AbstractGraph})

    Compile a list of graphs into a string for a python static function and output a python script.

# Arguments:
- `graphs`  vector of computational graphs
"""
function to_python_str(graphs::AbstractVector{<:AbstractGraph};
    root::AbstractVector{Int}=[id(g) for g in graphs], name::String="eval_graph", in_place::Bool=false)
    head = ""
    body = ""
    leafidx = 0
    inds_visitedleaf = Int[]
    inds_visitednode = Int[]
    map_validx_leaf = Dict{Int,eltype(graphs)}()  # mapping from the index of the leafVal to the leaf graph 
    for graph in graphs
        for g in PostOrderDFS(graph) #leaf first search
            g_id = id(g)
            target = "g$(g_id)"
            isroot = false
            if g_id in root
                target_root = "root[:, $(findfirst(x -> x == g_id, root)-1)]"
                isroot = true
            end
            if isempty(subgraphs(g)) #leaf
                g_id in inds_visitedleaf && continue
                body *= "    $target = leafVal[:, $(leafidx)]\n"
                leafidx += 1
                map_validx_leaf[leafidx] = g
                push!(inds_visitedleaf, g_id)
            else
                g_id in inds_visitednode && continue
                body *= "    $target = $(to_static(operator(g), subgraphs(g), subgraph_factors(g), lang=:python))\n"
                push!(inds_visitednode, g_id)
            end
            if isroot
                body *= "    $target_root = $target\n"
            end
        end
    end
    if in_place
        head *= "def $name(root, leafVal):\n"
    else
        head *= "import torch\n"
        head *= "def $name(leafVal):\n"
        head *= "    root = torch.empty(leafVal.shape[0], $(length(graphs)), dtype=leafVal.dtype, device=leafVal.device)\n"
    end
    tail = "    return root\n\n"

    return head * body * tail, map_validx_leaf
end
function compile_Python(graphs::AbstractVector{<:AbstractGraph}, filename::String;
    root::AbstractVector{Int}=[id(g) for g in graphs], func_name="eval_graph")
    py_string, leafmap = to_python_str(graphs, root=root, name=func_name)
    open(filename, "a") do f
        write(f, py_string)
    end
    return leafmap
end
