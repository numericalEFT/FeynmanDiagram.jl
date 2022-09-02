function oneOrderHigher(diag::Diagram{W}, ::Type{Id}, subdiagram=[];
    index::Int=index(Id), #which index to increase the order
    operator::Operator=diag.operator,
    name::Symbol=diag.name, factor::W=diag.factor) where {W,Id}
    # if diag.id isa PropagatorId && (diag.id isa Id) == false
    #     #for bare propagator, a derivative of different id vanishes
    #     return nothing
    # end
    id = deepcopy(diag.id)
    @assert index <= length(id.order) "$(id) only supports derivatives up to the index $(length(id.order))!"
    id.order[index] += 1
    d = Diagram{W}(id, operator, subdiagram; name=name, factor=factor)
    return d
end

function hasOrderHigher(diag::Diagram{W}, ::Type{ID}) where {W,ID<:PropagatorId}
    if diag.id isa ID
        #for bare propagator, a derivative of different id vanishes
        return true
    else
        return false
    end
end

"""
    function derivative(diags::Union{Tuple,AbstractVector}, ::Type{ID}; index::Int=index(ID)) where {W,ID<:PropagatorId}
    function derivative(diags::Vector{Diagram{W}}, ::Type{ID}; index::Int=index(ID)) where {W,ID<:PropagatorId}
    
    Automatic differentiation derivative on the diagrams

# Arguments
- diags     : diagrams to take derivative
- ID        : DiagramId to apply the differentiation
- index     : index of the id.order array element to increase the order
"""
function derivative(diags::Union{Tuple,AbstractVector}, ::Type{ID}; index::Int=index(ID)) where {W,ID<:PropagatorId}
    if isempty(diags)
        return diags
    else
        diags = collect(diags)
        diags = derivative(diags, ID; index=index)
        return diags
    end
end

function derivative(diags::Vector{Diagram{W}}, ::Type{ID}; index::Int=index(ID)) where {W,ID<:PropagatorId}
    # use a dictionary to host the dual diagram of a diagram for a given hash number
    # a dual diagram is defined as the derivative of the original diagram

    dual = Dict{Int,Diagram{W}}()
    for diag in diags
        for d::Diagram{W} in PostOrderDFS(diag) # postorder traversal will visit all subdiagrams of a diagram first
            if haskey(dual, d.hash)
                continue
            end
            if d.id isa PropagatorId
                # for propagators like bare Green's function and interaction, derivative simply means increase an order by one
                if hasOrderHigher(d, ID)
                    dual[d.hash] = oneOrderHigher(d, ID; index=index)
                end
            else # composite diagram
                if d.operator isa Sum
                    # for a diagram which is a sum of subdiagrams, derivative means a sub of derivative subdiagrams
                    children = [dual[sub.hash] for sub in d.subdiagram if haskey(dual, sub.hash)]
                    if isempty(children) == false
                        dual[d.hash] = oneOrderHigher(d, ID, children; index=index)
                    end
                elseif d.operator isa Prod
                    # d = s1xs2x... = s1'xs2x... + s1xs2'x... + ...
                    terms = Vector{Diagram{W}}()
                    for (si, sub) in enumerate(d.subdiagram)
                        if haskey(dual, sub.hash) == false
                            continue
                        end
                        children = [si == sj ? dual[sub.hash] : sub for (sj, sub) in enumerate(d.subdiagram)]
                        if isempty(children) == false
                            push!(terms, oneOrderHigher(d, ID, children; index=index))
                        end
                    end

                    if isempty(terms) == false
                        # !the summed dual diagram must have factor = 1.0
                        dual[d.hash] = oneOrderHigher(d, ID, terms; index=index, name=Symbol("$(d.name)'"), factor=W(1), operator=Sum())
                    end
                else
                    error("not implemented!")
                end
            end
        end
    end
    return [dual[diag.hash] for diag in diags if haskey(dual, diag.hash)]
end

"""
    function derivative(diags::Union{Diagram,Tuple,AbstractVector}, ::Type{ID}, order::Int) where {ID<:PropagatorId}
    
    Automatic differentiation derivative on the diagrams

# Arguments
- diags      : diagrams to take derivative
- ID         : DiagramId to apply the differentiation
- order::Int : derivative order
"""
function derivative(diags::Union{Tuple,AbstractVector}, ::Type{ID}, order::Int; index::Int=index(ID)) where {ID<:PropagatorId}
    @assert order >= 0
    if order == 0
        return diags
    end
    result = diags
    for o in 1:order
        result = derivative(result, ID; index=index)
    end
    return result
end

"""
    function removeHartreeFock!(diag::Diagram{W}) where {W}
    function removeHartreeFock!(diags::Union{Tuple,AbstractVector})
    
    Remove the Hartree-Fock insertions that without any derivatives on the propagator and the interaction. 

# Arguments
- diags      : diagrams to remove the Fock insertion

# Remarks
- The operations removeHartreeFock! and taking derivatives doesn't commute with each other! 
- If the input diagram is a Hartree-Fock diagram, then the overall weight will become zero! 
- The return value is always nothing
"""
function removeHartreeFock!(diag::Diagram{W}) where {W}
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
    end
end
function removeHartreeFock!(diags::Union{Tuple,AbstractVector})
    for diag in diags
        removeHartreeFock!(diag)
    end
end