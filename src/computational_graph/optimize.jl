function optimize!(graphs::Union{Tuple,AbstractVector{G}}; verbose=0, normalize=nothing) where {G<:AbstractGraph}
    if isempty(graphs)
        return nothing
    else
        graphs = collect(graphs)
        leaf_mapping = remove_duplicated_leaves!(graphs, verbose=verbose, normalize=normalize)
        remove_onechild_parents!(graphs, verbose=verbose)
        merge_nodes_prefactors!(graphs, verbose=verbose)
        return leaf_mapping
    end
end

function optimize(graphs::Union{Tuple,AbstractVector{G}}; verbose=0, normalize=nothing) where {G<:AbstractGraph}
    graphs_new = deepcopy(graphs)
    leaf_mapping = optimize!(graphs_new)
    return graphs_new, leaf_mapping
end

function remove_onechild_parents!(g::G; verbose=0) where {G<:AbstractGraph}
    verbose > 0 && println("remove nodes with only one child.")
    for sub_g in g.subgraphs
        remove_onechild_parents!(sub_g)
        inplace_prod!(sub_g)
    end
    inplace_prod!(g)
    return g
end

function remove_onechild_parents!(graphs::AbstractVector{G}; verbose=0) where {G<:AbstractGraph}
    verbose > 0 && println("remove nodes with only one child.")
    for g in graphs
        remove_onechild_parents!(g.subgraphs)
        for sub_g in g.subgraphs
            inplace_prod!(sub_g)
        end
        inplace_prod!(g)
    end
    return graphs
end

function merge_nodes_prefactors!(g::G; verbose=0) where {G<:AbstractGraph}
    verbose > 0 && println("merge nodes' subgraphs with multiplicative prefactors.")
    for sub_g in g.subgraphs
        merge_nodes_prefactors!(sub_g)
        merge_prefactors!(sub_g)
    end
    merge_prefactors!(g)
    return g
end

function merge_nodes_prefactors!(graphs::AbstractVector{G}; verbose=0) where {G<:AbstractGraph}
    verbose > 0 && println("merge nodes' subgraphs with multiplicative prefactors.")
    for g in graphs
        merge_nodes_prefactors!(g.subgraphs)
        for sub_g in g.subgraphs
            merge_prefactors!(sub_g)
        end
        merge_prefactors!(g)
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
