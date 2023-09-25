function optimize!(graphs::Union{Tuple,AbstractVector}; verbose=0, normalize=nothing)
    if isempty(graphs)
        return nothing
    else
        graphs = collect(graphs)
        removeOneChildParent!(graphs, verbose=verbose)
        mappings = removeDuplicatedLeaves!(graphs, verbose=verbose, normalize=normalize)
        return mappings
    end
end

function removeOneChildParent!(g::G; verbose=0) where {G<:AbstractGraph}
    verbose > 0 && println("remove nodes with only one child.")
    for sub_g in g.subgraphs
        removeOneChildParent!(sub_g)
        inplace_prod!(sub_g)
    end
    inplace_prod!(g)
    return g
end

function removeOneChildParent!(graphs::AbstractVector{G}; verbose=0) where {G<:AbstractGraph}
    verbose > 0 && println("remove nodes with only one child.")
    for g in graphs
        removeOneChildParent!(g.subgraphs)
        for sub_g in g.subgraphs
            inplace_prod!(sub_g)
        end
        inplace_prod!(g)
    end
    return graphs
end

function uniqueLeaves(_graphs::AbstractVector{G}) where {G<:AbstractGraph}
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

function removeDuplicatedLeaves!(graphs::AbstractVector{G}; verbose=0, normalize=nothing, kwargs...) where {G<:AbstractGraph}
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

    uniqueLeaf, leafMap = uniqueLeaves(leaves)
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
