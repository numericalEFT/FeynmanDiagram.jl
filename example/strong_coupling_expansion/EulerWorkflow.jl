module EulerWorkflow

export run_workflow, print_results

using IterTools
using Combinatorics

#=============================================================================#
# SECTION 1: GRAPH GENERATION LOGIC (Internal Implementation)
#=============================================================================#

# --- Custom Graph Representation and Functions ---

function _is_weakly_connected(v_count, edges)
    if v_count == 0
        return false
    end
    if isempty(edges)
        return v_count == 1
    end
    adj = [[] for _ in 1:v_count]
    for (u, v) in edges
        push!(adj[u], v)
        push!(adj[v], u)
    end
    q = [1]
    visited = falses(v_count)
    visited[1] = true
    count = 1
    while !isempty(q)
        u = popfirst!(q)
        for neighbor in adj[u]
            if !visited[neighbor]
                visited[neighbor] = true
                count += 1
                push!(q, neighbor)
            end
        end
    end
    return count == v_count
end

function _is_eulerian(v_count, edges)
    if v_count == 0
        return true
    end
    in_degrees = zeros(Int, v_count)
    out_degrees = zeros(Int, v_count)
    for (u, v) in edges
        if u > v_count || v > v_count
            return false
        end
        out_degrees[u] += 1
        in_degrees[v] += 1
    end
    return all(in_degrees .== out_degrees)
end

function _compositions(n, k)
    if k < 1
        return n == 0 ? [[]] : []
    end
    if k == 1
        return [[n]]
    end
    result = []
    for i in 0:n
        for p in _compositions(n - i, k - 1)
            push!(result, vcat([i], p))
        end
    end
    return result
end

# --- Recursive Skeleton Generation ---

function _get_skeleton_canonical_form(v_count, edges)
    v_labels = 1:v_count
    min_repr = ""
    for p in permutations(v_labels)
        relabel_map = Dict(zip(v_labels, p))
        edge_strings = sort!(["$(relabel_map[u])-$((relabel_map[v]))" for (u, v) in edges])
        current_repr = join(edge_strings, ";")
        if min_repr == "" || current_repr < min_repr
            min_repr = current_repr
        end
    end
    return min_repr
end

function _generate_recursive(m_target, current_edges, found_skeletons, found_canonical_forms)
    if length(current_edges) == m_target
        v_count = isempty(current_edges) ? 0 : maximum(e -> max(e[1], e[2]), current_edges)
        if _is_weakly_connected(v_count, current_edges) && _is_eulerian(v_count, current_edges)
            canonical_form = _get_skeleton_canonical_form(v_count, current_edges)
            if !(canonical_form in found_canonical_forms)
                push!(found_canonical_forms, canonical_form)
                push!(found_skeletons, (v_count, current_edges))
            end
        end
        return
    end
    max_v = isempty(current_edges) ? 0 : maximum(e -> max(e[1], e[2]), current_edges)
    for u in 1:max_v
        for v in 1:max_v
            _generate_recursive(m_target, [current_edges; (u, v)], found_skeletons, found_canonical_forms)
        end
    end
    if max_v > 0
        for u in 1:max_v
            _generate_recursive(m_target, [current_edges; (u, max_v + 1)], found_skeletons, found_canonical_forms)
            _generate_recursive(m_target, [current_edges; (max_v + 1, u)], found_skeletons, found_canonical_forms)
        end
    end
    if isempty(current_edges)
        _generate_recursive(m_target, [(1, 1)], found_skeletons, found_canonical_forms)
        _generate_recursive(m_target, [(1, 2)], found_skeletons, found_canonical_forms)
    end
end

function _generate_uncolored_skeletons(m)
    found_skeletons = []
    found_canonical_forms = Set()
    _generate_recursive(m, [], found_skeletons, found_canonical_forms)
    return found_skeletons
end

# --- Graph Coloring and Symmetry Factor Calculation ---

function _get_canonical_details(v, colored_edges)
    in_degrees = zeros(Int, v)
    out_degrees = zeros(Int, v)
    for (u, v_node, _) in colored_edges
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
    for combo_of_perms in product(partition_perms...)
        relabel_map = Dict()
        for (original_partition, permuted_partition) in zip(vertex_partitions, combo_of_perms)
            for (original_vertex, permuted_vertex) in zip(original_partition, permuted_partition)
                relabel_map[original_vertex] = permuted_vertex
            end
        end
        edge_strings = sort!(["$(relabel_map[u])-$((relabel_map[v_node]))-$(order)" for (u, v_node, order) in colored_edges])
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

function _color_graph(skeleton, orders)
    v, skeleton_edges = skeleton
    unique_colored_graphs = Dict()
    for p_orders in unique(permutations(orders))
        current_colored_edges = [(edge[1], edge[2], p_orders[i]) for (i, edge) in enumerate(skeleton_edges)]
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
    end
    return values(unique_colored_graphs)
end

#=============================================================================#
# SECTION 2: GRAPH DECOMPOSITION LOGIC (Internal Implementation)
#=============================================================================#

function _find_decomposition_paths(graph_data, graph_library)
    decomposition_paths = []

    primitive_splits = _find_primitive_and_remainder(graph_data, graph_library)

    for (id_A, mapping_A, data_B) in primitive_splits

        if isempty(data_B[1])
            continue
        else
            paths_for_B = _find_decomposition_paths(data_B, graph_library)

            if isempty(paths_for_B)
                id_B, map_B = _get_id_and_mapping(data_B, graph_library)
                if id_B != -1
                    new_path = sort([(id_A, mapping_A), (id_B, map_B)], by=x -> x[1])
                    push!(decomposition_paths, new_path)
                end
            else
                for path_B in paths_for_B
                    new_path = sort([(id_A, mapping_A); path_B], by=x -> x[1])
                    push!(decomposition_paths, new_path)
                end
            end
        end
    end

    return unique(decomposition_paths)
end

function _find_primitive_and_remainder(graph_data, graph_library)
    colored_edges, _, _ = graph_data
    num_edges = length(colored_edges)
    found_splits = []
    if num_edges == 0
        return found_splits
    end

    sorted_edges = sort(colored_edges, by=e -> (e[1], e[2], e[3]))
    first_edge = sorted_edges[1]
    other_edges = sorted_edges[2:end]

    for k in 0:length(other_edges)
        for subset_of_others_indices in Combinatorics.combinations(1:length(other_edges), k)

            subset_A_edges = [first_edge; other_edges[subset_of_others_indices]]

            if length(subset_A_edges) == num_edges
                continue
            end

            is_A_valid, id_A, mapping_A = _is_valid_subgraph(subset_A_edges, graph_library)

            if is_A_valid
                remainder_indices = setdiff(1:length(other_edges), subset_of_others_indices)
                subset_B_edges = other_edges[remainder_indices]

                v_count_B = isempty(subset_B_edges) ? 0 : maximum(e -> max(e[1], e[2]), subset_B_edges)
                data_B = (subset_B_edges, (0, 0), v_count_B)

                push!(found_splits, (id_A, mapping_A, data_B))
            end
        end
    end
    return unique(found_splits)
end

function _is_valid_subgraph(edges, graph_library)
    if isempty(edges)
        return (false, -1, [])
    end
    parent_verts = sort(unique(vcat([[e[1], e[2]] for e in edges]...)))
    v_count_dense = length(parent_verts)
    parent_to_dense_map = Dict(v => i for (i, v) in enumerate(parent_verts))
    edges_dense = [(parent_to_dense_map[e[1]], parent_to_dense_map[e[2]], e[3]) for e in edges]
    if !_is_weakly_connected(v_count_dense, [(e[1], e[2]) for e in edges_dense]) || !_is_eulerian(v_count_dense, [(e[1], e[2]) for e in edges_dense])
        return (false, -1, [])
    end
    canonical_form, _, dense_to_canonical_map = _get_canonical_details(v_count_dense, edges_dense)
    if haskey(graph_library, canonical_form)
        canonical_to_parent_map = Dict()
        for p_v in parent_verts
            dense_v = parent_to_dense_map[p_v]
            canon_v = dense_to_canonical_map[dense_v]
            canonical_to_parent_map[canon_v] = p_v
        end
        ordered_mapping = [canonical_to_parent_map[i] for i in 1:v_count_dense]
        return (true, graph_library[canonical_form].id, ordered_mapping)
    else
        return (false, -1, [])
    end
end

function _get_id_and_mapping(graph_data, graph_library)
    colored_edges, _, v_count = graph_data
    if isempty(colored_edges)
        return -1, []
    end

    canonical_form, _, dense_to_canonical_map = _get_canonical_details(v_count, colored_edges)

    if haskey(graph_library, canonical_form)
        canonical_to_original_map = Dict()
        for original_v in 1:v_count
            canon_v = dense_to_canonical_map[original_v]
            canonical_to_original_map[canon_v] = original_v
        end

        ordered_mapping = [canonical_to_original_map[i] for i in 1:v_count]
        return (graph_library[canonical_form].id, ordered_mapping)
    end

    return -1, []
end

#=============================================================================#
# SECTION 3: PUBLIC API FUNCTIONS
#=============================================================================#

"""
    run_workflow(M, N)

Performs all calculations and returns the graph catalog and decomposition results.
"""
function run_workflow(M, N)
    # --- Generation ---
    println("--- Generating Graph Catalog ---")
    all_graphs = Dict{Tuple{Int,Int},Vector}()
    for m in 1:M
        skeletons = _generate_uncolored_skeletons(m)
        for n_val in 0:N
            # Initialize the entry for this (m,n) if not exists
            if !haskey(all_graphs, (m, n_val))
                all_graphs[(m, n_val)] = []
            end

            order_partitions = unique([sort(p) for p in _compositions(n_val, m)])
            for orders in order_partitions
                println("orders: ", orders)
                for skeleton in skeletons
                    println(skeleton)
                    colored_graphs = _color_graph(skeleton, orders)
                    for graph_data in colored_graphs
                        push!(all_graphs[(m, n_val)], graph_data)
                    end
                end
            end
        end
    end

    # Calculate total graphs count
    total_graphs = sum(length(v) for v in values(all_graphs))
    println("Finished Generation. Found $total_graphs total unique graphs.")

    # --- Decomposition ---
    println("\n--- Analyzing Graph Decompositions ---")
    graph_library = Dict()

    # Build graph library with all graphs
    graph_id = 1
    for ((m, n), graphs_list) in all_graphs
        for graph_data in graphs_list
            colored_edges, _, v = graph_data
            canonical_form, _, _ = _get_canonical_details(v, colored_edges)
            graph_library[canonical_form] = (id=graph_id, data=graph_data)
            graph_id += 1
        end
    end
    println("Graph library created with $(length(graph_library)) entries.")

    all_decomposition_results = Dict{Tuple{Int,Int},Vector}()
    graph_id = 1
    for ((m, n), graphs_list) in all_graphs
        for graph_data in graphs_list
            all_paths = _find_decomposition_paths(graph_data, graph_library)
            if !isempty(all_paths)
                if !haskey(all_decomposition_results, (m, n))
                    all_decomposition_results[(m, n)] = []
                end
                push!(all_decomposition_results[(m, n)], (target_id=graph_id, paths=all_paths))
            end
            graph_id += 1
        end
    end

    # Calculate total decomposable graphs count
    total_decomposable = sum(length(v) for v in values(all_decomposition_results))
    println("Finished Decomposition Analysis. Found $total_decomposable decomposable graphs.")

    return all_graphs, all_decomposition_results
end

"""
    print_results(graph_catalog, decomposition_results)

Prints the graph catalog and decomposition relationships in a human-readable format.
"""
function print_results(graph_catalog, decomposition_results)
    println("\n========================================")
    println("          FINAL RESULTS")
    println("========================================")

    println("\n--- GRAPH CATALOG ---")
    if isempty(graph_catalog)
        println("No graphs were generated.")
    else
        # Sort the catalog by (m,n) keys for consistent printing
        sorted_keys = sort(collect(keys(graph_catalog)), by=x -> (x[1], x[2]))
        for (m, n) in sorted_keys
            graphs_list = graph_catalog[(m, n)]
            if !isempty(graphs_list)
                println("\n[m=$m, n=$n]:")
                for (i, graph_data) in enumerate(graphs_list)
                    # Find the global index of this graph
                    global_index = 0
                    current_count = 0
                    for (key, list) in graph_catalog
                        if key == (m, n)
                            global_index = current_count + 1
                            break
                        else
                            current_count += length(list)
                        end
                    end
                    global_index += i - 1

                    colored_edges, (sym_factor_v, sym_factor_e), v = graph_data
                    total_sym_factor = sym_factor_v * sym_factor_e
                    actual_m = length(colored_edges)
                    actual_n = sum(e[3] for e in colored_edges)

                    println("\n  - Graph #$global_index (m=$actual_m, n=$actual_n, v=$v)")
                    println("    Symmetry Factor: $total_sym_factor (Vertex: $sym_factor_v, Edge: $sym_factor_e)")
                    println("    Edges:")
                    sorted_edges = sort(colored_edges, by=e -> (e[1], e[2], e[3]))
                    for (u, v_node, order) in sorted_edges
                        println("      ($u, $v_node) with order $order")
                    end
                end
            end
        end
    end

    println("\n--- DECOMPOSITION RELATIONSHIPS ---")
    if isempty(decomposition_results)
        println("No decomposable graphs were found.")
    else
        # Sort the results by (m,n) keys for consistent printing
        sorted_keys = sort(collect(keys(decomposition_results)), by=x -> (x[1], x[2]))
        for (m, n) in sorted_keys
            results_list = decomposition_results[(m, n)]
            println("\n[m=$m, n=$n]:")
            for result in results_list
                # Find the global index of this graph
                global_index = 0
                current_count = 0
                for (key, list) in decomposition_results
                    if key == (m, n)
                        global_index = current_count + 1
                        break
                    else
                        current_count += length(list)
                    end
                end
                global_index += result.target_id - 1

                println("\n  ----------------------------------------")
                println("  Decomposition(s) for Graph #$global_index")
                for path in result.paths
                    if length(path) > 1
                        path_str = join(["Graph #$(id)($(join(mapping, ","))) " for (id, mapping) in path], " + ")
                        println("    - Path: " * path_str)
                    end
                end
                println("  ----------------------------------------")
            end
        end
    end
end

end # end module