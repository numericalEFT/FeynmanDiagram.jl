"""
    function to_python_str(graphs::AbstractVector{<:AbstractGraph})

    Compile a list of graphs into a string for a python static function and output a python script.

# Arguments:
- `graphs`  vector of computational graphs
"""
function to_python_str(graphs::AbstractVector{<:AbstractGraph})
    head = ""
    body = ""
    leafidx = 0
    root = [id(g) for g in graphs]
    inds_visitedleaf = Int[]
    inds_visitednode = Int[]
    gid_to_leafid = Dict{String,Int64}()
    rootidx = 0
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
                body *= "    $target = leaf[$(leafidx)]\n"
                gid_to_leafid[target] = leafidx
                leafidx += 1
                push!(inds_visitedleaf, g_id)
            else
                g_id in inds_visitednode && continue
                body *= "    $target = $(to_static(operator(g), subgraphs(g), subgraph_factors(g), lang=:python))\n"
                push!(inds_visitednode, g_id)
            end
            if isroot
                body *= "    root$(rootidx) = $target\n"
                rootidx += 1
            end
        end
    end
    head *= "def graphfunc(leaf):\n"
    output = ["root$(i)" for i in 0:rootidx-1]
    output = join(output, ",")
    tail = "    return $output\n\n"

    expr = head * body * tail

    return expr, gid_to_leafid
end
function compile_Python(graphs::AbstractVector{<:AbstractGraph}, filename::String="GraphFunc.py")
    py_string, leafmap = to_python_str(graphs)
    open(filename, "w") do f
        write(f, py_string)
    end
    return leafmap
end
