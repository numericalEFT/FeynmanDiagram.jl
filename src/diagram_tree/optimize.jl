function optimize(diag::Union{Tuple,AbstractVector}, optlevel = 1; kwargs...) end

"""
    removeduplication(diag::Union{Tuple,AbstractVector})

    remove duplicated nodes such as:  ---> ver4 ---> InteractionId. Leaf will not be touched!

"""
function removeOneChildParent(diag::Union{Diagram,Tuple,AbstractVector})
    if diag isa Diagram
        diag = [diag,]
    end
    for d in diag
        for (si, subdiag) in enumerate(d.subdiagram)
            if length(subdiag.subdiagram) == 1
                d.subdiagram[si] = subdiag.subdiagram[1]
            end
        end
    end
end