function optimize!(graphs::Union{Tuple,AbstractVector{G}}; verbose=0, normalize=nothing) where {G<:AbstractGraph}
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

function optimize(graphs::Union{Tuple,AbstractVector{G}}; verbose=0, normalize=nothing) where {G<:AbstractGraph}
    graphs_new = deepcopy(graphs)
    leaf_mapping = optimize!(graphs_new)
    return graphs_new, leaf_mapping
end

function merge_all_chain_prefactors!(g::G; verbose=0) where {G<:AbstractGraph}
    verbose > 0 && println("merge prefactors of all nodes representing trivial unary chains toward root level.")
    # Post-order DFS
    for sub_g in g.subgraphs
        merge_all_chain_prefactors!(sub_g)
        merge_chain_prefactors!(sub_g)
    end
    merge_chain_prefactors!(g)
    return g
end

function merge_all_chain_prefactors!(graphs::AbstractVector{G}; verbose=0) where {G<:AbstractGraph}
    verbose > 0 && println("merge prefactors of all nodes representing trivial unary chains toward root level.")
    # Post-order DFS
    for g in graphs
        merge_all_chain_prefactors!(g.subgraphs)
        merge_chain_prefactors!(g)
    end
    return graphs
end

function merge_all_factorless_chains!(g::G; verbose=0) where {G<:AbstractGraph}
    verbose > 0 && println("merge all nodes representing factorless trivial unary chains.")
    # Post-order DFS
    for sub_g in g.subgraphs
        merge_all_factorless_chains!(sub_g)
        merge_factorless_chain!(sub_g)
    end
    merge_factorless_chain!(g)
    return g
end

function merge_all_factorless_chains!(graphs::AbstractVector{G}; verbose=0) where {G<:AbstractGraph}
    verbose > 0 && println("merge all nodes representing factorless trivial unary chains.")
    # Post-order DFS
    for g in graphs
        merge_all_factorless_chains!(g.subgraphs)
        merge_factorless_chain!(g)
    end
    return graphs
end

function merge_all_chains!(g::G; verbose=0) where {G<:AbstractGraph}
    verbose > 0 && println("merge all nodes representing trivial unary chains.")
    merge_all_chain_prefactors!(g, verbose=verbose)
    merge_all_factorless_chains!(g, verbose=verbose)
    return g
end

function merge_all_chains!(graphs::AbstractVector{G}; verbose=0) where {G<:AbstractGraph}
    verbose > 0 && println("merge all nodes representing trivial unary chains.")
    merge_all_chain_prefactors!(graphs, verbose=verbose)
    merge_all_factorless_chains!(graphs, verbose=verbose)
    return graphs
end

function merge_all_linear_combinations!(g::G; verbose=0) where {G<:AbstractGraph}
    verbose > 0 && println("merge nodes representing a linear combination of a non-unique list of graphs.")
    # Post-order DFS
    for sub_g in g.subgraphs
        merge_all_linear_combinations!(sub_g)
        merge_linear_combination!(sub_g)
    end
    merge_linear_combination!(g)
    return g
end

function merge_all_linear_combinations!(graphs::AbstractVector{G}; verbose=0) where {G<:AbstractGraph}
    verbose > 0 && println("merge nodes representing a linear combination of a non-unique list of graphs.")
    # Post-order DFS
    for g in graphs
        merge_all_linear_combinations!(g.subgraphs)
        merge_linear_combination!(g)
    end
    return graphs
end

function unique_leaves(_graphs::AbstractVector{G}) where {G<:AbstractGraph}
    ############### find the unique Leaves #####################
    uniqueGraph = []
    mapping = Dict{Int,Int}()

    idx = 1
    for g in _graphs
        flag = true
        for (ie, e) in enumerate(uniqueGraph)
            if isequiv(e, g, :id)
                mapping[g.id] = ie
                flag = false
                break
            end
        end
        if flag
            push!(uniqueGraph, g)
            mapping[g.id] = idx
            idx += 1
        end
    end
    return uniqueGraph, mapping
end

function remove_duplicated_leaves!(graphs::AbstractVector{G}; verbose=0, normalize=nothing, kwargs...) where {G<:AbstractGraph}
    verbose > 0 && println("remove duplicated leaves.")
    leaves = Vector{G}()
    for g in graphs
        append!(leaves, collect(Leaves(g)))
    end
    if isnothing(normalize) == false
        @assert normalize isa Function "a function call is expected for normalize"
        for leaf in leaves
            normalize(leaf.id)
        end
    end
    sort!(leaves, by=x -> x.id) #sort the id of the leaves in an asscend order
    unique!(x -> x.id, leaves) #filter out the leaves with the same id number

    uniqueLeaf, leafMap = unique_leaves(leaves)
    verbose > 0 && length(leaves) > 0 && println("Number of independent Leaves $(length(leaves)) → $(length(uniqueLeaf))")

    for g in graphs
        for n in PreOrderDFS(g)
            for (si, sub_g) in enumerate(n.subgraphs)
                if isleaf(sub_g)
                    n.subgraphs[si] = uniqueLeaf[leafMap[sub_g.id]]
                end
            end
        end
    end

    return leafMap
end
