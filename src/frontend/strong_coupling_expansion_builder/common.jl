
function parity(p)
    """
    calculate the parity of a given permutation of the array [1, 2, 3, ...]
    """
    n = length(p)
    not_seen = Set{Int}(1:n)
    seen = Set{Int}()
    cycles = Array{Int,1}[]
    while !isempty(not_seen)
        cycle = Int[]
        x = pop!(not_seen)
        while !in(x, seen)
            push!(cycle, x)
            push!(seen, x)
            x = p[x]
            pop!(not_seen, x, 0)
        end
        push!(cycles, cycle)
    end
    cycle_lengths = map(length, cycles)
    even_cycles = filter(i -> i % 2 == 0, cycle_lengths)
    length(even_cycles) % 2 == 0 ? 1 : -1
end

struct Partition
    """
    2-partition of 2N-legs of Green's functions
    """
    l::Vector{Vector{Int}}
    r::Vector{Vector{Int}}
    sign::Vector{Int}
    function Partition(order)
        N = 2^order - 2
        function addincoming(l, r)
            # map 1, 2, 3, ... to order+1, order+2, ...
            l = l .+ order
            r = r .+ order
            nl, nr = length(l), length(r)
            inl = [i for i in 1:nl]
            inr = [i for i in nl+1:order]
            append!(inl, l)
            append!(inr, r)
            return inl, inr
        end

        left = Vector{Vector{Int}}(undef, N)
        right = similar(left)
        sign = Vector{Int}(undef, N)
        for (si, s) in enumerate(partitions(1:order, 2))
            # generate 2-partition of the outgoing legs
            # all incoming legs are kept the same
            s1, s2 = addincoming(s[1], s[2])
            left[2si-1], right[2si-1] = s1, s2
            sign[2si-1] = parity(vcat(s1, s2))

            q1, q2 = addincoming(s[2], s[1])
            left[2si], right[2si] = q1, q2
            sign[2si] = parity(vcat(q1, q2))
            # println("($s1, $s2) -> $p1, ($q1, $q2) -> $p2")
        end
        return new(left, right, sign)
    end
end

"""
    build_graph_data(edge_list)

Processes a list of directed edges to build necessary graph structures.

# Arguments
- `edge_list`: Vector of tuples `(u, v)` representing edges `u -> v`.

# Returns
- `nodes`: Set of all unique nodes.
- `adj_undirected`: Adjacency list treating edges as undirected (for connectivity checks).
- `in_degrees`: Dictionary of in-degrees (for closed check).
- `out_degrees`: Dictionary of out-degrees (for closed check).
"""
function build_graph_data(edge_list::Vector{Tuple{Int,Int}})
    nodes = Set{Int}()
    adj_undirected = Dict{Int,Vector{Int}}()
    in_degrees = Dict{Int,Int}()
    out_degrees = Dict{Int,Int}()

    for (u, v) in edge_list
        push!(nodes, u)
        push!(nodes, v)

        # Build undirected adjacency (for traversal/connectivity)
        push!(get!(adj_undirected, u, Int[]), v)
        push!(get!(adj_undirected, v, Int[]), u)

        # Build directed degrees (for is_closed check)
        out_degrees[u] = get(out_degrees, u, 0) + 1
        in_degrees[v] = get(in_degrees, v, 0) + 1
    end

    return nodes, adj_undirected, in_degrees, out_degrees
end

"""
    is_connected(edge_list)

Checks if the graph is weakly connected (all nodes are reachable from each other 
if we ignore edge direction).
"""
function is_connected(edge_list::Vector{Tuple{Int,Int}})
    nodes, adj_undirected, _, _ = build_graph_data(edge_list)

    if length(nodes) <= 1
        return true
    end

    # BFS Traversal (treating graph as undirected)
    visited = Set{Int}()
    queue = [first(nodes)] # Start from an arbitrary node

    while !isempty(queue)
        u = popfirst!(queue)
        if u ∉ visited
            push!(visited, u)
            # Visit neighbors
            for v in get(adj_undirected, u, Int[])
                if v ∉ visited
                    push!(queue, v)
                end
            end
        end
    end

    return length(visited) == length(nodes)
end

"""
    is_closed(edge_list)

Checks if the graph is "closed" according to the directed definition:
For every node, its in-degree must equal its out-degree.
"""
function is_closed(edge_list::Vector{Tuple{Int,Int}})
    nodes, _, in_degrees, out_degrees = build_graph_data(edge_list)

    if isempty(nodes)
        return true
    end

    for node in nodes
        n_in = get(in_degrees, node, 0)
        n_out = get(out_degrees, node, 0)

        if n_in != n_out
            return false
        end
    end

    return true
end

"""
    get_connected_component_indices(edge_list)

Splits the graph into weakly connected components.

# Returns
- `Vector{Vector{Int}}`: A list of components. Each component is a list 
  of indices (1-based) referring to the edges in the original `edge_list`.
"""
function get_connected_component_indices(edge_list::Vector{Tuple{Int,Int}})::Vector{Vector{Int}}
    nodes, adj_undirected, _, _ = build_graph_data(edge_list)

    visited_nodes = Set{Int}()
    components_indices = Vector{Vector{Int}}()

    for start_node in nodes
        if start_node ∈ visited_nodes
            continue
        end

        # --- New Component Found ---
        current_component_nodes = Set{Int}()
        queue = [start_node]

        # BFS to find all nodes in this component
        while !isempty(queue)
            u = popfirst!(queue)
            if u ∉ current_component_nodes
                push!(current_component_nodes, u)
                push!(visited_nodes, u)

                for v in get(adj_undirected, u, Int[])
                    if v ∉ current_component_nodes
                        push!(queue, v)
                    end
                end
            end
        end

        # --- Map Edges to this Component ---
        # Find which edge INDICES belong to the nodes we just found
        current_indices = Int[]
        for (i, (u, v)) in enumerate(edge_list)
            # If a node is in this component, the edge belongs to it.
            # (We only need to check 'u' because if 'u' is in, 'v' must be too).
            if u ∈ current_component_nodes
                push!(current_indices, i)
            end
        end

        push!(components_indices, current_indices)
    end

    return components_indices
end

# """
#     build_graph_info(edge_list)

# Builds an adjacency list, a set of nodes, and a degree count for each node
# from a list of edges.

# # Arguments
# - `edge_list::Vector{Tuple{Int, Int}}`: A Vector where each element is a 
#   `Tuple{Int, Int}` representing an edge.

# # Returns
# - `adj::Dict{Int, Vector{Int}}`: Adjacency list
# - `nodes::Set{Int}`: Set of all unique nodes
# - `degrees::Dict{Int, Int}`: Dictionary mapping each node to its degree
# """
# function build_graph_info(edge_list::Vector{Tuple{Int,Int}})
#     adj = Dict{Int,Vector{Int}}()
#     nodes = Set{Int}()
#     degrees = Dict{Int,Int}()

#     for (u, v) in edge_list
#         # 1. Add both nodes to the node set
#         push!(nodes, u)
#         push!(nodes, v)

#         # 2. Build the adjacency list (for an undirected graph)
#         # get!(dict, key, default) is useful:
#         # If key exists, return its value; otherwise, create it with 
#         # the default value and return that default.
#         push!(get!(adj, u, Int[]), v)
#         push!(get!(adj, v, Int[]), u)

#         # 3. Calculate degrees
#         # get(dict, key, default) returns the value for key, or default 
#         # if key is not present.
#         degrees[u] = get(degrees, u, 0) + 1
#         degrees[v] = get(degrees, v, 0) + 1
#     end

#     return adj, nodes, degrees
# end

# """
#     is_connected(edge_list)

# Checks if the graph defined by the edge list is connected.
# """
# function is_connected(edge_list::Vector{Tuple{Int,Int}})
#     # We only need the adjacency list and the set of nodes
#     adj, nodes, _ = build_graph_info(edge_list)

#     # Edge case: A graph with 0 or 1 nodes is considered connected.
#     if length(nodes) <= 1
#         return true
#     end

#     visited = Set{Int}()
#     # Use a Vector as a queue (push! to add to the end, popfirst! to remove from the front)
#     queue = [first(nodes)] # Start the search from an arbitrary node

#     while !isempty(queue)
#         u = popfirst!(queue)

#         if u ∉ visited
#             push!(visited, u)

#             # Add all unvisited neighbors to the queue
#             # Use get(adj, u, []) to safely handle nodes with no neighbors
#             for v in get(adj, u, Int[])
#                 if v ∉ visited
#                     push!(queue, v)
#                 end
#             end
#         end
#     end

#     # The graph is connected if the number of visited nodes
#     # equals the total number of nodes.
#     return length(visited) == length(nodes)
# end

# """
#     is_closed(edge_list)

# Checks if the graph is "closed," defined as every node having a degree of at least 2.
# """
# function is_closed(edge_list::Vector{Tuple{Int,Int}}; is_even=true)
#     # We only need the set of nodes and the degree counts
#     _, nodes, degrees = build_graph_info(edge_list)

#     # Edge case: An empty graph (no nodes) vacuously satisfies the 
#     # condition "all nodes...".
#     if isempty(nodes)
#         return true
#     end

#     for node in nodes
#         # Get the degree for the current node.
#         # Our helper function guarantees all nodes in the `nodes` set are also in the `degrees` dictionary.
#         node_degree = get(degrees, node, 0) # Use get for safety

#         if node_degree < 2
#             # Found a node with degree < 2, so the graph is not "closed"
#             return false
#         end
#         if is_even && isodd(node_degree)
#             return false
#         end
#     end

#     # If the loop finishes, all nodes had a degree of 2 or more.
#     return true
# end

# """
#     get_connected_component_indices(edge_list)
# Takes a graph's edge list and divides the graph into its connected components.

# # Arguments
# - `edge_list::Vector{Tuple{Int, Int}}`: The edge list of the full graph.

# # Returns
# - `Vector{Vector{Int}}`: A Vector where each element is a `Vector{Int}` containing the 1-based indices
#   of the edges (from the original `edge_list`) that belong to that component.
# """
# function get_connected_component_indices(edge_list::Vector{Tuple{Int,Int}})::Vector{Vector{Int}}

#     # 1. Build the graph structure for traversal
#     adj, all_nodes, _ = build_graph_info(edge_list)

#     # Tracks nodes already assigned to a component
#     visited_nodes = Set{Int}()

#     # Final list of [component_indices_1, component_indices_2, ...]
#     all_components_indices = Vector{Vector{Int}}()

#     # 2. Iterate through every node in the graph
#     for start_node in all_nodes

#         # If we've already processed this node, skip it
#         if start_node ∈ visited_nodes
#             continue
#         end

#         # --- Found a new component ---

#         # 3. Find all nodes in this new component (using BFS)
#         current_component_nodes = Set{Int}()
#         queue = [start_node]
#         push!(visited_nodes, start_node)

#         while !isempty(queue)
#             u = popfirst!(queue)
#             push!(current_component_nodes, u) # Add node to this component's set

#             for v in get(adj, u, Int[])
#                 if v ∉ visited_nodes
#                     push!(visited_nodes, v) # Mark as globally visited
#                     push!(queue, v)
#                 end
#             end
#         end

#         # 4. Now, find all *indices* of edges that belong to this component.
#         #    We iterate through the original edge_list with their indices.
#         current_component_indices = Int[]

#         # `enumerate` in Julia is 1-based, perfect for this.
#         for (idx, (u, v)) in enumerate(edge_list)

#             # If either node of the edge is in our component set,
#             # the edge belongs to this component.
#             if u ∈ current_component_nodes
#                 push!(current_component_indices, idx)
#             end
#         end

#         # 5. Add this component's index list to our main list
#         push!(all_components_indices, current_component_indices)
#     end

#     return all_components_indices
# end

function _get_canonical_details(v, colored_edges)
    in_degrees = zeros(Int, v)
    out_degrees = zeros(Int, v)
    for (u, v_node) in colored_edges
        out_degrees[u] += 1
        in_degrees[v_node] += 1
    end
    partitions = Dict()
    for i in 1:v
        key = (in_degrees[i], out_degrees[i])
        if !haskey(partitions, key)
            partitions[key] = []
        end
        push!(partitions[key], i)
    end

    sorted_partition_keys = sort(collect(keys(partitions)))
    vertex_partitions = [partitions[key] for key in sorted_partition_keys]

    min_repr = ""
    automorphism_count = 0
    canonical_mapping = Dict()
    partition_perms = [permutations(p) for p in vertex_partitions]
    for combo_of_perms in Iterators.product(partition_perms...)
        relabel_map = Dict()
        for (original_partition, permuted_partition) in zip(vertex_partitions, combo_of_perms)
            for (original_vertex, permuted_vertex) in zip(original_partition, permuted_partition)
                relabel_map[original_vertex] = permuted_vertex
            end
        end
        edge_strings = sort!(["$(relabel_map[u])-$((relabel_map[v_node]))" for (u, v_node) in colored_edges])
        current_repr = join(edge_strings, ";")
        if min_repr == "" || current_repr < min_repr
            min_repr = current_repr
            automorphism_count = 1
            canonical_mapping = relabel_map
        elseif current_repr == min_repr
            automorphism_count += 1
        end
    end
    return min_repr, automorphism_count, canonical_mapping
end

function get_symmetry_factor(v, skeleton_edges)
    unique_colored_graphs = Dict()
    current_colored_edges = [(edge[1], edge[2]) for (i, edge) in enumerate(skeleton_edges)]
    canonical_repr, sym_factor_v, _ = _get_canonical_details(v, current_colored_edges)
    if !haskey(unique_colored_graphs, canonical_repr)
        edge_counts = Dict()
        for edge in current_colored_edges
            edge_counts[edge] = get(edge_counts, edge, 0) + 1
        end
        sym_factor_e = 1
        for count in values(edge_counts)
            sym_factor_e *= factorial(count)
        end
        unique_colored_graphs[canonical_repr] = (current_colored_edges, (sym_factor_v, sym_factor_e), v)
    end

    # return values(unique_colored_graphs)
    return sym_factor_e * sym_factor_v
end