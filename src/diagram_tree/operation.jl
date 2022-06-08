function oneOrderHigher(diag::Diagram{W}, ::Type{Id}, subdiagram=[]) where {W,Id}
    if diag.id isa PropagatorId && (diag.id isa Id) == false
        #for bare propagator, a derivative of different id vanishes
        return nothing
    end
    id = deepcopy(diag.id)
    if Id == BareGreenId
        id.order[1] += 1
    elseif Id == BareInteractionId
        id.order[2] += 1
    else
        error("not implemented!")
    end

    d = Diagram{W}(id, diag.operator, subdiagram; name=diag.name, factor=diag.factor, weight=diag.weight)
    return d
end

"""
    function derivative(diags::Union{Diagram,Tuple,AbstractVector}, ::Type{ID}) where {ID<:PropagatorId}
    
    Automatic differentiation derivative on the diagrams

# Arguments
- diags     : diagrams to take derivative
- ID        : DiagramId to apply the differentiation
"""
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
                        newid = deepcopy(id)
                        newid.order .= terms[1].id.order
                        dual[d.hash] = Diagram{id.para.weightType}(newid, Sum(), terms, name=Symbol("$(d.name)'"))
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

"""
    function derivative(diags::Union{Diagram,Tuple,AbstractVector}, ::Type{ID}, order::Int) where {ID<:PropagatorId}
    
    Automatic differentiation derivative on the diagrams

# Arguments
- diags      : diagrams to take derivative
- ID         : DiagramId to apply the differentiation
- order::Int : derivative order
"""
function derivative(diags::Union{Diagram,Tuple,AbstractVector}, ::Type{ID}, order::Int) where {ID<:PropagatorId}
    @assert order >= 0
    if order == 0
        return diags
    end
    result = diags
    for o in 1:order
        result = derivative(result, ID)
    end
    return result
end

"""
    function removeHatreeFock!(diags::Union{Diagram,Tuple,AbstractVector})
    
    Remove the Hatree-Fock insertions that without any derivatives on the propagator and the interaction. 

# Arguments
- diags      : diagrams to remove the Fock insertion

# Remarks
- The operations removeHatreeFock! and taking derivatives doesn't commute with each other! 
- If the input diagram is a Hatree-Fock diagram, then the overall weight will become zero! 
"""
function removeHatreeFock!(diags::Union{Diagram,Tuple,AbstractVector})
    single = false
    if diags isa Diagram
        diags = [diags,]
        single = true
    end
    for diag in diags
        for d in PreOrderDFS(diag)
            # for subdiag in d.subdiagram
            if d.id isa SigmaId
                # println(d, " with ", d.id.para.innerLoopNum)
                if isempty(d.id.order) || all(x -> x == 0, d.id.order) #eithr order is empty or a vector of zeros
                    if d.id.para.innerLoopNum == 1
                        d.factor = 0.0
                    end
                end
            end
            # end
        end
    end
    if single
        return diags[1]
    else
        return diags
    end
end