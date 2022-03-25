function optimize(diag::Union{Diagram,Tuple,AbstractVector}, optlevel = 1; verbose = 0, kwargs...)
    single = false
    diag = deepcopy(diag)
    if diag isa Diagram
        single = true
        diag = [diag,]
    end
    removeOneChildParent!(diag, verbose = verbose)
    removeDuplicatedLeaves!(diag, verbose = verbose)
    return single ? diag[1] : diag
end

"""
    removeOneChildParent!(diags::AbstractVector; verbose = 0)

    remove duplicated nodes such as:  ---> ver4 ---> InteractionId. Leaf will not be touched!
"""
function removeOneChildParent!(diags::AbstractVector; verbose = 0)
    verbose > 0 && println("remove nodes with only one child.")
    for diag in diags
        #deep first search, remove one-child parent from the leaf level first
        removeOneChildParent!(diag.subdiagram)
        #then remove the one-child subdiagram of the current diagram
        for (si, subdiag) in enumerate(diag.subdiagram)
            if length(subdiag.subdiagram) == 1
                subdiag.subdiagram[1].factor *= subdiag.factor
                diag.subdiagram[si] = subdiag.subdiagram[1]
            end
        end
    end
    return diags
end

"""
    removeDuplicatedLeaves!(diags::AbstractVector; verbose = 0)

    remove duplicated nodes such as:  ---> ver4 ---> InteractionId. Leaf will not be touched!
"""
function removeDuplicatedLeaves!(diags::AbstractVector; verbose = 0)
    verbose > 0 && println("remove duplicated leaves.")
    leaves = []
    for diag in diags
        #leaves must be the propagators
        append!(leaves, collect(Leaves(diag)))
    end
    # println([d.hash for d in leaves])
    sort!(leaves, by = x -> x.hash) #sort the hash of the leaves in an asscend order
    unique!(x -> x.hash, leaves) #filter out the leaves with the same hash number

    for l in leaves
        #make sure all leaves are either Green's functions or interactions
        @assert l.id isa PropagatorId
    end

    function uniqueLeaves(_diags::AbstractVector)
        ############### find the unique Leaves #####################
        uniqueDiag = []
        mapping = Dict{Int,Any}()
        for diag in _diags
            flag = true
            for (ei, e) in enumerate(uniqueDiag)
                if e.factor ≈ diag.factor && e.id == diag.id
                    mapping[diag.hash] = e
                    flag = false
                    break
                end
            end
            if flag
                push!(uniqueDiag, diag)
                # push!(mapping, length(uniqueDiag))
                mapping[diag.hash] = diag
            end
        end
        return uniqueDiag, mapping
    end

    # println(leaves)
    green = [l for l in leaves if l.id isa BareGreenId]
    interaction = [l for l in leaves if l.id isa BareInteractionId]
    greenN = [l for l in leaves if l.id isa BareGreenNId]
    hopping = [l for l in leaves if l.id isa BareHoppingId]
    # println(green)
    # println(interaction)

    uniqueGreen, greenMap = uniqueLeaves(green)
    uniqueGreenN, greenNMap = uniqueLeaves(greenN)
    uniqueInteraction, interactionMap = uniqueLeaves(interaction)
    uniqueHopping, hoppingMap = uniqueLeaves(hopping)
    # println(uniqueInteraction)
    # display(greenMap)

    verbose > 0 && length(green) > 0 && println("Number of independent Greens $(length(green)) → $(length(uniqueGreen))")
    verbose > 0 && length(greenN) > 0 && println("Number of independent GreenNs $(length(greenN)) → $(length(uniqueGreenN))")
    verbose > 0 && length(interaction) > 0 && println("Number of independent Interactions $(length(interaction)) → $(length(uniqueInteraction))")
    verbose > 0 && length(hopping) > 0 && println("Number of independent Hopping $(length(hopping)) → $(length(uniqueHopping))")

    for diag in diags
        for n in PreOrderDFS(diag)
            for (si, subdiag) in enumerate(n.subdiagram)
                @assert (n.id isa PropagatorId) == false "the diagram $n with subdiagrams cannot be a proapgator!"

                if subdiag.id isa PropagatorId
                    if subdiag.id isa BareGreenId
                        n.subdiagram[si] = greenMap[subdiag.hash]
                    elseif subdiag.id isa BareInteractionId
                        n.subdiagram[si] = interactionMap[subdiag.hash]
                    elseif subdiag.id isa BareGreenNId
                        n.subdiagram[si] = greenNMap[subdiag.hash]
                    elseif subdiag.id isa BareHoppingId
                        n.subdiagram[si] = hoppingMap[subdiag.hash]
                    else
                        error("not implemented!")
                    end
                end
            end
        end
    end

    return uniqueGreen, uniqueInteraction
    # return diags
end