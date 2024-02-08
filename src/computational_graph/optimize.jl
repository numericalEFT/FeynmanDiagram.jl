"""
    function optimize!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; level=0, verbose=0, normalize=nothing)

    In-place optimization of given `graphs`. Removes duplicated leaves, flattens chains, 
    merges linear combinations, and removes zero-valued subgraphs. When `level > 0`, also removes duplicated intermediate nodes.

# Arguments:
- `graphs`: A tuple or vector of graphs.
- `level`: Optimization level (default: 0). A value greater than 0 triggers more extensive but slower optimization processes, such as removing duplicated intermediate nodes.
- `verbose`: Level of verbosity (default: 0).
- `normalize`: Optional function to normalize the graphs (default: nothing).

# Returns
- Returns the optimized graphs. If the input graphs is empty, it returns nothing.
"""
function optimize!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; level=0, verbose=0, normalize=nothing)
    if isempty(graphs)
        return nothing
    else
        if level > 0
            if graphs isa Tuple
                root = Graph(collect(graphs))
            else
                root = Graph(graphs)
            end
            remove_duplicated_nodes!(root, verbose=verbose)
        else
            remove_duplicated_leaves!(graphs, verbose=verbose, normalize=normalize)
        end

        flatten_all_chains!(graphs, verbose=verbose)
        merge_all_linear_combinations!(graphs, verbose=verbose)
        remove_all_zero_valued_subgraphs!(graphs, verbose=verbose)
        return graphs
    end
end

"""
    function optimize(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; level=0, verbose=0, normalize=nothing)

    Optimizes a copy of given `graphs`. Removes duplicated nodes (when `level > 0`) or leaves, flattens chains, 
    merges linear combinations, and removing zero-valued subgraphs.

# Arguments:
- `graphs`: A tuple or vector of graphs.
- `level`: Optimization level (default: 0). A value greater than 0 triggers more extensive but slower optimization processes, such as removing duplicated nodes.
- `verbose`: Level of verbosity (default: 0).
- `normalize`: Optional function to normalize the graphs (default: nothing).

# Returns:
- A tuple/vector of optimized graphs.
"""
function optimize(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; level=0, verbose=0, normalize=nothing)
    graphs_new = deepcopy(graphs)
    optimize!(graphs_new, level=level, verbose=verbose, normalize=normalize)
    return graphs_new
end

"""
    function flatten_all_chains!(g::AbstractGraph; verbose=0)
F
    Flattens all nodes representing trivial unary chains in-place in the given graph `g`. 

# Arguments:
- `graphs`: The graph to be processed.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- The mutated graph `g` with all chains flattened.
"""
function flatten_all_chains!(g::AbstractGraph; verbose=0)
    verbose > 0 && println("flatten all nodes representing trivial unary chains.")
    for sub_g in g.subgraphs
        flatten_all_chains!(sub_g)
        flatten_chains!(sub_g)
    end
    flatten_chains!(g)
    return g
end

"""
    function flatten_all_chains!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0)

    Flattens all nodes representing trivial unary chains in-place in the given graphs.

# Arguments:
- `graphs`: A collection of graphs to be processed.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- The mutated collection `graphs` with all chains in each graph flattened.
"""
function flatten_all_chains!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0)
    verbose > 0 && println("flatten all nodes representing trivial unary chains.")
    # Post-order DFS
    for g in graphs
        flatten_all_chains!(g.subgraphs)
        flatten_chains!(g)
    end
    return graphs
end

"""
    function remove_all_zero_valued_subgraphs!(g::AbstractGraph; verbose=0)

    Recursively removes all zero-valued subgraph(s) in-place in the given graph `g`.

# Arguments:
- `g`: An AbstractGraph.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graph.
# 
"""
function remove_all_zero_valued_subgraphs!(g::AbstractGraph; verbose=0)
    verbose > 0 && println("merge nodes representing a linear combination of a non-unique list of graphs.")
    # Post-order DFS
    for sub_g in subgraphs(g)
        remove_all_zero_valued_subgraphs!(sub_g)
        remove_zero_valued_subgraphs!(sub_g)
    end
    remove_zero_valued_subgraphs!(g)
    return g
end

"""
    function remove_all_zero_valued_subgraphs!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0)

    Recursively removes all zero-valued subgraph(s) in-place in the given graphs.

# Arguments:
- `graphs`: A collection of graphs to be processed.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graphs.
# 
"""
function remove_all_zero_valued_subgraphs!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0)
    verbose > 0 && println("merge nodes representing a linear combination of a non-unique list of graphs.")
    # Post-order DFS
    for g in graphs
        remove_all_zero_valued_subgraphs!(subgraphs(g))
        remove_zero_valued_subgraphs!(g)
    end
    return graphs
end

"""
    function merge_all_linear_combinations!(g::AbstractGraph; verbose=0)

    Merges all nodes representing a linear combination of a non-unique list of subgraphs in-place in the given graph `g`.

# Arguments:
- `g`: An AbstractGraph.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graph.
# 
"""
function merge_all_linear_combinations!(g::AbstractGraph; verbose=0)
    verbose > 0 && println("merge nodes representing a linear combination of a non-unique list of graphs.")
    # Post-order DFS
    for sub_g in subgraphs(g)
        merge_all_linear_combinations!(sub_g)
        merge_linear_combination!(sub_g)
    end
    merge_linear_combination!(g)
    return g
end

"""
    function merge_all_linear_combinations!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0)

    Merges all nodes representing a linear combination of a non-unique list of subgraphs in-place in the given graphs. 

# Arguments:
- `graphs`: A collection of graphs to be processed.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graphs.
# 
"""
function merge_all_linear_combinations!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0)
    verbose > 0 && println("merge nodes representing a linear combination of a non-unique list of graphs.")
    # Post-order DFS
    for g in graphs
        merge_all_linear_combinations!(subgraphs(g))
        merge_linear_combination!(g)
    end
    return graphs
end

"""
    function merge_all_multi_products!(g::Graph; verbose=0)

    Merges all nodes representing a multi product of a non-unique list of subgraphs in-place in the given graph `g`.

# Arguments:
- `g::Graph`: A Graph.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graph.
# 
"""
function merge_all_multi_products!(g::Graph; verbose=0)
    verbose > 0 && println("merge nodes representing a multi product of a non-unique list of graphs.")
    # Post-order DFS
    for sub_g in g.subgraphs
        merge_all_multi_products!(sub_g)
        merge_multi_product!(sub_g)
    end
    merge_multi_product!(g)
    return g
end

"""
    function merge_all_multi_products!(graphs::Union{Tuple,AbstractVector{<:Graph}}; verbose=0)

    Merges all nodes representing a multi product of a non-unique list of subgraphs in-place in the given graphs. 

# Arguments:
- `graphs`: A collection of graphs to be processed.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graphs.
# 
"""
function merge_all_multi_products!(graphs::Union{Tuple,AbstractVector{<:Graph}}; verbose=0)
    verbose > 0 && println("merge nodes representing a multi product of a non-unique list of graphs.")
    # Post-order DFS
    for g in graphs
        merge_all_multi_products!(g.subgraphs)
        merge_multi_product!(g)
    end
    return graphs
end

"""
    function unique_nodes!(graphs::AbstractVector{<:AbstractGraph})

    Identifies and retrieves unique nodes from a set of graphs.

# Arguments:
- `graphs`: A collection of graphs to be processed.

# Returns:
- A mapping dictionary from the id of each leaf to the unique leaf node.
"""
function unique_nodes!(graphs::AbstractVector{<:AbstractGraph}, mapping::Dict{Int,<:AbstractGraph}=Dict{Int,eltype(graphs)}())
    # function unique_nodes!(graphs::AbstractVector{<:AbstractGraph})
    ############### find the unique Leaves #####################
    # unique_graphs = []
    # mapping = Dict{Int,eltype(graphs)}()
    unique_graphs = collect(values(mapping))

    for g in graphs
        flag = true
        for e in unique_graphs
            if isequiv(e, g, :id, :name, :weight)
                mapping[id(g)] = e
                flag = false
                break
            end
        end
        if flag
            push!(unique_graphs, g)
            mapping[id(g)] = g
        end
    end
    return mapping
end

"""
    function remove_duplicated_leaves!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0, normalize=nothing, kwargs...)

    Removes duplicated leaf nodes in-place from a collection of graphs. It also provides optional normalization for these leaves.

# Arguments:
- `graphs`: A collection of graphs to be processed.
- `verbose`: Level of verbosity (default: 0).
- `normalize`: Optional function to normalize the graphs (default: nothing).
"""
function remove_duplicated_leaves!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0, normalize=nothing, kwargs...)
    verbose > 0 && println("remove duplicated leaves.")
    leaves = Vector{eltype(graphs)}()
    for g in graphs
        append!(leaves, collect(Leaves(g)))
    end
    if isnothing(normalize) == false
        @assert normalize isa Function "a function call is expected for normalize"
        for leaf in leaves
            normalize(id(leaf))
        end
    end
    sort!(leaves, by=x -> id(x)) #sort the id of the leaves in an asscend order
    unique!(x -> id(x), leaves) #filter out the leaves with the same id number

    mapping = unique_nodes!(leaves)

    for g in graphs
        for n in PreOrderDFS(g)
            for (si, sub_g) in enumerate(subgraphs(n))
                if isleaf(sub_g)
                    set_subgraph!(n, mapping[id(sub_g)], si)
                end
            end
        end
    end

    return graphs
end

function remove_duplicated_nodes!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0, kwargs...)
    verbose > 0 && println("remove duplicated nodes.")

    nodes_all = Vector{eltype(graphs)}()
    for g in graphs
        for node in PostOrderDFS(g)
            push!(nodes_all, node)
        end
    end

    sort!(nodes_all, by=x -> id(x)) #sort the id of the leaves in an asscend order
    unique!(x -> id(x), nodes_all) #filter out the leaves with the same id number

    mapping = unique_nodes!(nodes_all)

    for g in graphs
        for n in PreOrderDFS(g)
            for (si, sub_g) in enumerate(subgraphs(n))
                set_subgraph!(n, mapping[id(sub_g)], si)
            end
        end
    end

    return graphs
end

function remove_duplicated_nodes!(root::G; verbose=0) where {G<:AbstractGraph}
    verbose > 0 && println("remove duplicated nodes.")
    # A dictionary to keep track of unique nodes based on a key (like id, or a hash of properties)

    # remove_duplicated_leaves!([root])

    unique_nodes = Dict{Int,G}()
    # for l in Leaves(root)
    #     if !haskey(unique_nodes, id(l))
    #         unique_nodes[id(l)] = l
    #     end
    # end

    # Helper function to process a node
    function process_node(node)
        # Compute a key for the node (here, I'm using `id` for simplicity)
        node_key = id(node)

        # Check if a node with the same key already exists
        if haskey(unique_nodes, node_key)
            return unique_nodes[node_key]
        else
            # Check if the node is equivalent to any existing unique node
            for g in values(unique_nodes)
                if isequiv(node, g, :id, :name, :weight)
                    return g
                end
            end

            # Process child nodes if the node is unique
            for (i, child) in enumerate(subgraphs(node))
                unique_child = process_node(child)
                set_subgraph!(node, unique_child, i)
            end

            # Add the (now potentially updated) node to the unique_nodes dictionary
            unique_nodes[node_key] = node
            return node
        end
    end

    # Start processing from the root
    process_node(root)

    return root
end

"""
    function burn_from_targetleaves!(graphs::AbstractVector{G}, targetleaves_id::AbstractVector{Int}; verbose=0) where {G <: AbstractGraph}

    Removes all nodes connected to the target leaves in-place via "Prod" operators.

# Arguments:
- `graphs`: A vector of graphs.
- `targetleaves_id::AbstractVector{Int}`: Vector of target leafs' id.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- The id of a constant graph with a zero factor if any graph in `graphs` was completely burnt; otherwise, `nothing`.
"""
function burn_from_targetleaves!(graphs::AbstractVector{G}, targetleaves_id::AbstractVector{Int}; verbose=0) where {G<:AbstractGraph}
    verbose > 0 && println("remove all nodes connected to the target leaves via Prod operators.")

    graphs_sum = linear_combination(graphs, one.(eachindex(graphs)))
    ftype = eltype(subgraph_factors(graphs[1]))

    for leaf in Leaves(graphs_sum)
        if !isdisjoint(id(leaf), targetleaves_id)
            set_name!(leaf, "BURNING")
        end
    end

    for node in PostOrderDFS(graphs_sum)
        if any(x -> name(x) == "BURNING", subgraphs(node))
            if operator(node) == Prod || operator(node) <: Power
                set_subgraphs!(node, G[])
                set_subgraph_factors!(node, ftype[])
                set_name!(node, "BURNING")
            else
                _subgraphs = G[]
                _subgraph_factors = ftype[]
                for (i, subg) in enumerate(subgraphs(node))
                    if name(subg) != "BURNING"
                        push!(_subgraphs, subg)
                        push!(_subgraph_factors, subgraph_factor(node, i))
                    end
                end
                set_subgraphs!(node, _subgraphs)
                set_subgraph_factors!(node, _subgraph_factors)
                if isempty(_subgraph_factors)
                    set_name!(node, "BURNING")
                end
            end
        end
    end

    # g_c0 = constant_graph(ftype(0))
    g_c1 = constant_graph(ftype(1))
    has_c0 = false
    for g in graphs
        if name(g) == "BURNING"
            has_c0 = true
            set_id!(g, id(g_c1))
            set_operator!(g, Unitary)
            # set_subgraphs!(g, subgraphs(g_c0))
            # set_subgraph_factors!(g, subgraph_factors(g_c0))
            set_weight!(g, 0.0)
        end
    end

    has_c0 ? (return id(g_c1)) : (return nothing)
end