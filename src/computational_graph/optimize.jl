function optimize!(graphs::Union{Tuple,AbstractVector}; verbose=0, normalize=nothing)
    if isempty(graphs)
        # return graphs
        return nothing
    else
        graphs = collect(graphs)
        removeOneChildParent!(graphs, verbose=verbose)
        mappings = removeDuplicatedLeaves!(graphs, verbose=verbose, normalize=normalize)
        return mappings
        # return graphs, mappings
    end
end

function removeOneChildParent!(graphs::AbstractVector{G}; verbose=0) where {G<:Graph}
    verbose > 0 && println("remove nodes with only one child.")
    for g in graphs
        removeOneChildParent!(g.subgraphs)
        for sub_g in g.subgraphs
            merge_prodchain_subfactors!(sub_g)
        end
    end
    return graphs
end

function removeDuplicatedLeaves!(graphs::AbstractVector{G}; verbose=0, normalize=nothing, kwargs...) where {G<:Graph}
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

    for l in leaves
        #make sure all leaves are either propagators or interactions
        @assert l.type in [Interaction, Propagator]
    end

    function uniqueLeaves(_graphs::Vector{G}) where {G}
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
                # push!(mapping, length(uniqueGraph))
                mapping[g.id] = idx
                idx += 1
            end
        end
        return uniqueGraph, mapping
    end

    green = [l for l in leaves if l.type == Propagator]
    interaction = [l for l in leaves if l.type == Interaction]

    uniqueGreen, greenMap = uniqueLeaves(green)
    uniqueInteraction, interactionMap = uniqueLeaves(interaction)

    verbose > 0 && length(green) > 0 && println("Number of independent Propagators $(length(green)) → $(length(uniqueGreen))")
    verbose > 0 && length(interaction) > 0 && println("Number of independent Interactions $(length(interaction)) → $(length(uniqueInteraction))")

    for g in graphs
        for n in PreOrderDFS(g)
            for (si, sub_g) in enumerate(n.subgraphs)
                if sub_g.type == Propagator
                    n.subgraphs[si] = uniqueGreen[greenMap[sub_g.id]]
                elseif sub_g.type == Interaction
                    n.subgraphs[si] = uniqueInteraction[interactionMap[sub_g.id]]
                end
            end
        end
    end

    # return uniqueGreen, uniqueInteraction
    return greenMap, interactionMap
end

# function removeDuplicatedLeaves!(graphs::AbstractVector{G}; verbose=0, normalize=nothing, kwargs...) where {G<:Graph}
#     verbose > 0 && println("remove duplicated leaves.")
#     uniqueGreen, uniqueInteraction = Vector{G}(), Vector{G}()
#     for g in graphs
#         for l in Leaves(g)
#             if l.type == Interaction
#                 loc = findfirst(x -> isequiv(x, l, :id), uniqueInteraction)
#                 if !isnothing(loc)
#                     l.id = uniqueInteraction[loc].id
#                 else
#                     push!(uniqueInteraction, l)
#                 end
#             elseif l.type == Propagator
#                 loc = findfirst(x -> isequiv(x, l, :id), uniqueGreen)
#                 if !isnothing(loc)
#                     l.id = uniqueGreen[loc].id
#                 else
#                     push!(uniqueGreen, l)
#                 end
#             else
#                 error("the leaf's type cannot be $(l.type)!")
#             end
#         end
#     end
#     return uniqueGreen, uniqueInteraction
# end