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
    for diag in diags
        for (si, subdiag) in enumerate(diag.subdiagram)
            if length(subdiag.subdiagram) == 1
                subdiag.subdiagram[1].factor *= subdiag.factor
                diag.subdiagram[si] = subdiag.subdiagram[1]
            end
        end
        removeOneChildParent!(diag.subdiagram)
    end
    return diags
end

"""
    removeDuplicatedLeaves!(diags::AbstractVector; verbose = 0)

    remove duplicated nodes such as:  ---> ver4 ---> InteractionId. Leaf will not be touched!
"""
function removeDuplicatedLeaves!(diags::AbstractVector; verbose = 0)
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
        @assert (l.id isa BareGreenId) || (l.id isa BareInteractionId)
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
    # println(green)
    # println(interaction)

    uniqueGreen, greenMap = uniqueLeaves(green)
    uniqueInteraction, interactionMap = uniqueLeaves(interaction)
    # println(uniqueInteraction)
    # display(greenMap)

    verbose > 0 && println("Number of independent Greens $(length(green)) → $(length(uniqueGreen))")
    verbose > 0 && println("Number of independent Interactions $(length(interaction)) → $(length(uniqueInteraction))")

    for diag in diags
        for n in PreOrderDFS(diag)
            for (si, subdiag) in enumerate(n.subdiagram)
                @assert (n.id isa PropagatorId) == false "a diagram with subdiagrams cannot be a proapgator!"

                if subdiag.id isa PropagatorId
                    if subdiag.id isa BareGreenId
                        n.subdiagram[si] = greenMap[subdiag.hash]
                    elseif subdiag.id isa BareInteractionId
                        n.subdiagram[si] = interactionMap[subdiag.hash]
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