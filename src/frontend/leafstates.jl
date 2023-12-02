function leafstates(leaf_maps::Vector{Dict{Int,G}}, coeff_map::Vector{Dict{Int,Tuple{Diagram{W},Vector{Int}}}},
    maxloopNum::Int) where {G<:Graph,W}

    num_g = length(leaf_maps)
    leafType = [Vector{Int}() for _ in 1:num_g]
    leafOrders = [Vector{Vector{Int}}() for _ in 1:num_g]
    leafInTau = [Vector{Int}() for _ in 1:num_g]
    leafOutTau = [Vector{Int}() for _ in 1:num_g]
    leafLoopIndex = [Vector{Int}() for _ in 1:num_g]
    leafValue = [Vector{Float64}() for _ in 1:num_g]

    loopbasis = Vector{Float64}[]
    # tau_labels = Vector{Int}[]
    for (ikey, leafmap) in enumerate(leaf_maps)
        len_leaves = length(keys(leafmap))
        sizehint!(leafType[ikey], len_leaves)
        sizehint!(leafOrders[ikey], len_leaves)
        sizehint!(leafInTau[ikey], len_leaves)
        sizehint!(leafOutTau[ikey], len_leaves)
        sizehint!(leafLoopIndex[ikey], len_leaves)
        leafValue[ikey] = ones(W, len_leaves)

        for idx in 1:len_leaves
            leaf = leafmap[idx]
            diag, orders = coeff_map[leaf.id]
            diagId = diag.id
            loopmom = copy(diagId.extK)
            len = length(loopmom)
            @assert DiagTree.isbare(diag)
            @assert maxloopNum >= len

            if maxloopNum > length(loopmom)
                append!(loopmom, zeros(Float64, maxloopNum - len))
            end
            flag = true
            for bi in eachindex(loopbasis)
                if loopbasis[bi] ≈ loopmom
                    push!(leafLoopIndex[ikey], bi)
                    flag = false
                    break
                end
            end
            if flag
                push!(loopbasis, loopmom)
                push!(leafLoopIndex[ikey], length(loopbasis))
            end

            # push!(tau_labels, collect(diagId.extT))
            push!(leafInTau[ikey], diagId.extT[1])
            push!(leafOutTau[ikey], diagId.extT[2])

            push!(leafOrders[ikey], orders)
            push!(leafType, DiagTree.index(diagId))

        end
    end

    return (leafValue, leafType, leafOrders, leafInTau, leafOutTau, leafLoopIndex), loopbasis
end