function oneOrderHigher(diag::Diagram{W}, ::Type{Id}, subdiagram = []) where {W,Id}
    if diag.id isa PropagatorId && (diag.id isa Id) == false
        #for bare propagator, a derivative of different id vanishes
        return nothing
    end
    order = deepcopy(diag.order)
    if Id == BareGreenId
        order[1] += 1
    elseif Id == BareInteractionId
        order[2] += 1
    else
        error("not implemented!")
    end

    d = Diagram{W}(diag.id, diag.operator, subdiagram; name = diag.name, factor = diag.factor, weight = diag.weight, order = order)
    return d
end

function derivative(diags::Union{Diagram,Tuple,AbstractVector}, ::Type{ID}) where {ID<:PropagatorId}
    # use a dictionary to host the dual diagram of a diagram for a given hash number
    # a dual diagram is defined as the derivative of the original diagram

    single = false
    if diags isa Diagram
        diags = [diags,]
        single = true
    end
    dual = Dict{Int,Any}()
    for diag in diags
        for d in PostOrderDFS(diag)
            if haskey(dual, d.hash)
                continue
            end
            id = d.id
            if id isa PropagatorId
                # for propagators like bare Green's function and interaction, derivative simply means increase an order by one
                dual[d.hash] = oneOrderHigher(d, ID)
            else # composite diagram
                if d.operator isa Sum
                    # for a diagram which is a sum of subdiagrams, derivative means a sub of derivative subdiagrams
                    children = [dual[sub.hash] for sub in d.subdiagram if isnothing(dual[sub.hash]) == false]
                    if isempty(children)
                        dual[d.hash] = nothing
                    else
                        dual[d.hash] = oneOrderHigher(d, ID, children)
                    end
                elseif d.operator isa Prod
                    # d = s1xs2x... = s1'xs2x... + s1xs2'x... + ...
                    terms = []
                    for (si, sub) in enumerate(d.subdiagram)
                        if isnothing(dual[sub.hash])
                            continue
                        end
                        children = deepcopy(d.subdiagram)
                        children[si] = dual[sub.hash]
                        dd = oneOrderHigher(d, ID, children)
                        if isnothing(dd) == false
                            push!(terms, dd)
                        end
                    end
                    if isempty(terms)
                        dual[d.hash] = nothing
                    else
                        dual[d.hash] = Diagram{id.para.weightType}(id, Sum(), terms, name = Symbol("$(diag.name)'"), order = terms[1].order)
                    end
                else
                    error("not implemented!")
                end
            end
        end
    end
    if single
        return dual[diags[1].hash]
    else
        return [dual[diag.hash] for diag in diags]
    end
end