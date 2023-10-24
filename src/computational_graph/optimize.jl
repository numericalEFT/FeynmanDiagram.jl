"""
    function optimize!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0, normalize=nothing)

    In-place optimization of given `graphs`. Removes duplicated leaves, merges chains, and merges linear combinations.

# Arguments:
- `graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}`: A tuple or vector of graphs.
- `verbose`: Level of verbosity (default: 0).
- `normalize`: Optional function to normalize the graphs (default: nothing).

# Returns:
- A mapping dictionary from the id of each unique leaf node to its index in collect(1:length(leafs)).
"""
function optimize!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0, normalize=nothing)
    if isempty(graphs)
        return nothing
    else
        graphs = collect(graphs)
        leaf_mapping = remove_duplicated_leaves!(graphs, verbose=verbose, normalize=normalize)
        merge_all_chains!(graphs, verbose=verbose)
        merge_all_linear_combinations!(graphs, verbose=verbose)
        return leaf_mapping
    end
end

"""
    function optimize(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0, normalize=nothing)

    Optimize a copy of given `graphs`. Removes duplicated leaves, merges chains, and merges linear combinations.

# Arguments:
- `graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}`: A tuple or vector of graphs.
- `verbose`: Level of verbosity (default: 0).
- `normalize`: Optional function to normalize the graphs (default: nothing).

# Returns:
- A tuple/vector of optimized graphs.
- A mapping dictionary from the id of each unique leaf node to its index in collect(1:length(leafs)).
"""
function optimize(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0, normalize=nothing)
    graphs_new = deepcopy(graphs)
    leaf_mapping = optimize!(graphs_new, verbose=verbose, normalize=normalize)
    return graphs_new, leaf_mapping
end

"""
    function merge_all_chain_prefactors!(g::AbstractGraph; verbose=0)

    In-place merge prefactors of all nodes representing trivial unary chains towards the root level for a single graph.

# Arguments:
- `g::AbstractGraph`: An AbstractGraph.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graph.
# 
"""
function merge_all_chain_prefactors!(g::AbstractGraph; verbose=0)
    verbose > 0 && println("merge prefactors of all nodes representing trivial unary chains toward root level.")
    # Post-order DFS
    for sub_g in subgraphs(g)
        merge_all_chain_prefactors!(sub_g)
        merge_chain_prefactors!(sub_g)
    end
    merge_chain_prefactors!(g)
    return g
end

"""
    function merge_all_chain_prefactors!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0)

    In-place merge prefactors of all nodes representing trivial unary chains towards the root level for given graphs.

# Arguments:
- `graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}`: A tuple or vector of graphs.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graphs.
# 
"""
function merge_all_chain_prefactors!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0)
    verbose > 0 && println("merge prefactors of all nodes representing trivial unary chains toward root level.")
    # Post-order DFS
    for g in graphs
        merge_all_chain_prefactors!(subgraphs(g))
        merge_chain_prefactors!(g)
    end
    return graphs
end

"""
    function merge_all_factorless_chains!(g::AbstractGraph; verbose=0)

    In-place merge all nodes representing factorless trivial unary chains within a single graph.

# Arguments:
- `g::AbstractGraph`: An AbstractGraph.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graph.
# 
"""
function merge_all_factorless_chains!(g::AbstractGraph; verbose=0)
    verbose > 0 && println("merge all nodes representing factorless trivial unary chains.")
    # Post-order DFS
    for sub_g in subgraphs(g)
        merge_all_factorless_chains!(sub_g)
        merge_factorless_chain!(sub_g)
    end
    merge_factorless_chain!(g)
    return g
end

"""
    function merge_all_factorless_chains!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0)

    In-place merge all nodes representing factorless trivial unary chains within given graphs.

# Arguments:
- `graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}`: A tuple or vector of graphs.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graphs.
# 
"""
function merge_all_factorless_chains!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0)
    verbose > 0 && println("merge all nodes representing factorless trivial unary chains.")
    # Post-order DFS
    for g in graphs
        merge_all_factorless_chains!(subgraphs(g))
        merge_factorless_chain!(g)
    end
    return graphs
end

"""
    function merge_all_chains!(g::AbstractGraph; verbose=0)

    In-place merge all nodes representing trivial unary chains within a single graph.
    This function consolidates both chain prefactors and factorless chains.

# Arguments:
- `g::AbstractGraph`: An AbstractGraph.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graph.
# 
"""
function merge_all_chains!(g::AbstractGraph; verbose=0)
    verbose > 0 && println("merge all nodes representing trivial unary chains.")
    merge_all_chain_prefactors!(g, verbose=verbose)
    merge_all_factorless_chains!(g, verbose=verbose)
    return g
end

"""
    function merge_all_chains!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0)

    In-place merge all nodes representing trivial unary chains in given graphs. 
    This function consolidates both chain prefactors and factorless chains.

# Arguments:
- `graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}`: A tuple or vector of graphs.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graphs.
# 
"""
function merge_all_chains!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0)
    verbose > 0 && println("merge all nodes representing trivial unary chains.")
    merge_all_chain_prefactors!(graphs, verbose=verbose)
    merge_all_factorless_chains!(graphs, verbose=verbose)
    return graphs
end

"""
    function merge_all_linear_combinations!(g::AbstractGraph; verbose=0)

    In-place merge all nodes representing a linear combination of a non-unique list of subgraphs within a single graph.

# Arguments:
- `g::AbstractGraph`: An AbstractGraph.
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

    In-place merge all nodes representing a linear combination of a non-unique list of subgraphs in given graphs. 

# Arguments:
- `graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}`: A tuple or vector of graphs.
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

    In-place merge all nodes representing a multi product of a non-unique list of subgraphs within a single graph.

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
    function merge_all_multi_products!(graphs::Union{Tuple,AbstractVector{Graph}}; verbose=0)

    In-place merge all nodes representing a multi product of a non-unique list of subgraphs in given graphs. 

# Arguments:
- `graphs::Union{Tuple,AbstractVector{Graph}}`: A tuple or vector of graphs.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graphs.
# 
"""
function merge_all_multi_products!(graphs::Union{Tuple,AbstractVector{Graph}}; verbose=0)
    verbose > 0 && println("merge nodes representing a multi product of a non-unique list of graphs.")
    # Post-order DFS
    for g in graphs
        merge_all_multi_products!(g.subgraphs)
        merge_multi_product!(g)
    end
    return graphs
end

"""
    function unique_leaves(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}})

    Identify and retrieve unique leaf nodes from a set of graphs.

# Arguments:
- `graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}`: A tuple or vector of graphs.

# Returns:
- The vector of unique leaf nodes.
- A mapping dictionary from the id of each unique leaf node to its index in collect(1:length(leafs)).
"""
function unique_leaves(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}})
    ############### find the unique Leaves #####################
    unique_graphs = []
    mapping = Dict{Int,Int}()

    idx = 1
    for g in graphs
        flag = true
        for (ie, e) in enumerate(unique_graphs)
            if isequiv(e, g, :id)
                mapping[id(g)] = ie
                flag = false
                break
            end
        end
        if flag
            push!(unique_graphs, g)
            mapping[id(g)] = idx
            idx += 1
        end
    end
    return unique_graphs, mapping
end

"""
    function remove_duplicated_leaves!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0, normalize=nothing, kwargs...)

    In-place remove duplicated leaf nodes from a collection of graphs. It also provides optional normalization for these leaves.

# Arguments:
- `graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}`: A tuple or vector of graphs.
- `verbose`: Level of verbosity (default: 0).
- `normalize`: Optional function to normalize the graphs (default: nothing).

# Returns:
- A mapping dictionary from the id of each unique leaf node to its index in collect(1:length(leafs)).
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

    _unique_leaves, leafmap = unique_leaves(leaves)
    verbose > 0 && length(leaves) > 0 && println("Number of independent Leaves $(length(leaves)) → $(length(_unique_leaves))")

    for g in graphs
        for n in PreOrderDFS(g)
            for (si, sub_g) in enumerate(subgraphs(n))
                if isleaf(sub_g)
                    set_subgraph!(n, _unique_leaves[leafmap[id(sub_g)]], si)
                end
            end
        end
    end

    return leafmap
end

"""
    function burn_from_targetleaves!(graphs::AbstractVector{<:AbstractGraph}, targetleaves_id::AbstractVector{Int}; verbose=0)

    In-place remove all nodes connected to the target leaves via "Prod" operators.

# Arguments:
- `graphs::AbstractVector{<:AbstractGraph}`: A vector of graphs.
- `targetleaves_id::AbstractVector{Int}`: Vector of target leafs' id.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- The id of a constant graph with a zero factor if any graph in `graphs` was completely burnt; otherwise, `nothing`.
"""
function burn_from_targetleaves!(graphs::AbstractVector{<:AbstractGraph}, targetleaves_id::AbstractVector{Int}; verbose=0)
    verbose > 0 && println("remove all nodes connected to the target leaves via Prod operators.")

    graphs_sum = linear_combination(graphs, one.(eachindex(graphs)))
    ftype = typeof(factor(graphs[1]))

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

    g_c0 = constant_graph(ftype(0))
    has_c0 = false
    for g in graphs
        if name(g) == "BURNING"
            has_c0 = true
            set_id!(g, id(g_c0))
            set_operator!(g, Constant)
            set_factor!(g, ftype(0))
        end
    end

    has_c0 ? (return id(g_c0)) : (return nothing)
end