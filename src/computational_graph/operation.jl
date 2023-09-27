# function oneOrderHigher(diag::Graph{F,W}, ::Type{Id}, subdiagram=[];
#     index::Int=index(Id), #which index to increase the order
#     operator::Operator=diag.operator,
#     name::Symbol=diag.name, factor::W=diag.factor) where {F, W,Id}
#     # if diag.id isa PropagatorId && (diag.id isa Id) == false
#     #     #for bare propagator, a derivative of different id vanishes
#     #     return nothing
#     # end
#     id = deepcopy(diag.id)
#     @assert index <= length(id.order) "$(id) only supports derivatives up to the index $(length(id.order))!"
#     id.order[index] += 1
#     d = Graph{F,W}(id, operator, subdiagram; name=name, factor=factor)
#     return d
# end

# function hasOrderHigher(diag::Graph{F,W}, ::Type{ID}) where {F, W,ID<:PropagatorId}
#     if diag.id isa ID
#         #for bare propagator, a derivative of different id vanishes
#         return true
#     else
#         return false
#     end
# end

# """
#     function derivative(diags::Union{Tuple,AbstractVector}, ::Type{ID}; index::Int=index(ID)) where {W,ID<:PropagatorId}
#     function derivative(diags::Vector{Graph{F,W}}, ::Type{ID}; index::Int=index(ID)) where {W,ID<:PropagatorId}

#     Automatic differentiation derivative on the diagrams

# # Arguments
# - diags     : diagrams to take derivative
# - ID        : DiagramId to apply the differentiation
# - index     : index of the id.order array element to increase the order
# """
# function derivative(diags::Union{Tuple,AbstractVector}, ::Type{ID}) where {F, W,ID<:PropagatorId}
#     if isempty(diags)
#         return diags
#     else
#         diags = collect(diags)
#         diags = derivative(diags, ID; index=index)
#         return diags
#     end
# end
function separate_number_graph(g::Vector{Union{F,Graph{F,W}}}) where {F,W}
    subgraphs = Vector{Graph{F,W}}()
    subnumber = nothing
    for child in g
        if typeof(child) <: Number
            if isnothing(subnumber)
                subnumber = child
            else
                subnumber += child
            end
        elseif typeof(child) <: Graph{F,W}
            push!(subgraphs, child)
        else
            error("The type of subgraphs in derivative is incorrect!")
        end
    end
    return subgraphs, subnumber
end

function separate_number_graph(g::Vector{Union{F,Graph{F,W}}}, coeff::Vector{F}) where {F,W}
    @assert length(g) == length(coeff)
    subgraphs = Vector{Graph{F,W}}()
    subnumber = nothing
    subcoeff = Vector{F}()
    for (i, child) in enumerate(g)
        if typeof(child) <: Number
            if isnothing(subnumber)
                subnumber = child * coeff[i]
            else
                subnumber += child * coeff[i]
            end
        elseif typeof(child) <: Graph{F,W}
            push!(subgraphs, child)
            push!(subcoeff, coeff[i])
        else
            error("The type of subgraphs in derivative is incorrect!")
        end
    end
    return subgraphs, subnumber, subcoeff
end

function frontAD(diag::Graph{F,W}, ID::Int) where {F,W}
    # use a dictionary to host the dual diagram of a diagram for a given hash number
    # a dual diagram is defined as the derivative of the original diagram
    rootid = -1
    dual = Dict{Int,Union{F,Graph{F,W}}}()
    for d::Graph{F,W} in PostOrderDFS(diag) # postorder traversal will visit all subdiagrams of a diagram first
        rootid = d.id
        if haskey(dual, d.id)
            continue
        end
        if isleaf(d)
            # For leaves, derivative with respect to same leaf is 1, otherwise is 0 (no dual graph in this case).
            if d.id == ID
                dual[d.id] = F(1) * d.factor
            end
        else # composite diagram
            if d.operator == Sum
                children = Vector{Union{F,Graph{F,W}}}()
                coeff = Vector{F}()
                for (i, sub) in enumerate(d.subgraphs)
                    if haskey(dual, sub.id)
                        push!(children, dual[sub.id])
                        push!(coeff, d.subgraph_factors[i])
                    end
                end

                subgraphs, subnumber, subcoeff = separate_number_graph(children, coeff)
                if isempty(subgraphs) == false
                    if !isnothing(subnumber)
                        push!(subgraphs, constant_graph(F(subnumber)))     #If both numbers and graphs appear in derivative, convert number to a unity graph, and asign the number to subgraph_factors of parent node.
                        push!(subcoeff, 1.0)
                    end
                    dual[d.id] = linear_combination(subgraphs, subcoeff)
                    dual[d.id].factor *= d.factor
                elseif !isnothing(subnumber)  #if only numbers appear in derivative, return a number
                    dual[d.id] = subnumber * d.factor
                end
            elseif d.operator == Prod
                # d = s1xs2x... = s1'xs2x... + s1xs2'x... + ...
                factor = 1.0
                children = Vector{Union{F,Graph{F,W}}}()
                for (si, sub) in enumerate(d.subgraphs) # First generate each addend s1'xs2x...
                    if haskey(dual, sub.id) == false
                        continue
                    end
                    factor *= d.subgraph_factors[si]
                    child = dual[sub.id]
                    for (sj, sub) in enumerate(d.subgraphs)
                        if !(si == sj)
                            child = child * d.subgraphs[sj]
                        end
                    end
                    push!(children, child)
                end

                subgraphs, subnumber = separate_number_graph(children)

                if isempty(subgraphs) == false
                    if !isnothing(subnumber)
                        push!(subgraphs, constant_graph(F(subnumber)))     #If both numbers and graphs appear in derivative, convert number to a unity graph, and asign the number to subgraph_factors of parent node.
                        push!(subcoeff, 1.0)
                    end
                    subcoeff = ones(F, length(subgraphs))
                    dual[d.id] = linear_combination(subgraphs, subcoeff)
                    dual[d.id].factor *= d.factor * factor
                elseif !isnothing(subnumber)  #if only numbers appear in derivative, return a number
                    dual[d.id] = subnumber * d.factor * factor
                end
            else
                error("not implemented!")
            end
        end
    end
    if isempty(dual)
        return 0.0
    end
    return dual[rootid]
end


# function backAD(diag::Graph{F,W}, ID::Int; unity = Graph([]; type = diag.type(), ftype = F, wtype = W,  weight = W(1.0)) )where {F, W}
#     # use a dictionary to host the dual diagram of a diagram for a given hash number
#     # a dual diagram is defined as the derivative of the original diagram
#     rootid = -1
#     AD_array = Dict{Int,Graph{F, W}}()
#     dual = Dict{Int,Union{F, Graph{F, W}}}()
#     for d::Graph{F,W} in PreOrderDFS(diag) # postorder traversal will visit all subdiagrams of a diagram first
#         if isleaf(d)
#             AD_array[d.id] = dual[d.id]
#         end
#         if haskey(dual, d.id)
#             continue
#         end
#         if isleaf(d) 
#             # For leaves, derivative with respect to same leaf is 1, otherwise is 0 (no dual graph in this case).
#             if d.id == ID
#                 dual[d.id] = F(1)*d.factor
#             end
#         else # composite diagram
#             if d.operator == Sum
#                 children = Vector{Union{F,Graph{F,W}}}()
#                 coeff = Vector{F}()
#                 for (i,sub) in enumerate(d.subgraphs)
#                     if haskey(dual, sub.id)
#                         push!( children, dual[sub.id])
#                         push!(coeff, d.subgraph_factors[i])
#                     end
#                 end

#                 subgraphs, subnumber, subcoeff = separate_number_graph(children, coeff)
#                 if isempty(subgraphs) == false
#                     if !isnothing(subnumber)
#                         push!(subgraphs, unity)     #If both numbers and graphs appear in derivative, convert number to a unity graph, and asign the number to subgraph_factors of parent node.
#                         push!(subcoeff, subnumber) 
#                     end
#                     dual[d.id] =linear_combination(subgraphs,subcoeff)
#                     dual[d.id].factor *= d.factor
#                 elseif !isnothing(subnumber)  #if only numbers appear in derivative, return a number
#                     dual[d.id] = subnumber*d.factor
#                 end
#             elseif d.operator == Prod
#                 # d = s1xs2x... = s1'xs2x... + s1xs2'x... + ...
#                 factor = 1.0
#                 children = Vector{Union{F,Graph{F,W}}}()
#                 for (si, sub) in enumerate(d.subgraphs) # First generate each addend s1'xs2x...
#                     if haskey(dual, sub.id) == false
#                         continue
#                     end
#                     factor *= d.subgraph_factors[si]
#                     child = dual[sub.id]
#                     for (sj, sub) in enumerate(d.subgraphs)
#                         if !(si == sj)
#                             child = child * d.subgraphs[sj]
#                         end
#                     end
#                     push!(children, child)
#                 end

#                 subgraphs, subnumber = separate_number_graph(children)

#                 if isempty(subgraphs) == false
#                     if !isnothing(subnumber)
#                         push!(subgraphs, unity)     #If both numbers and graphs appear in derivative, convert number to a unity graph, and asign the number to subgraph_factors of parent node.
#                         push!(subcoeff, subnumber) 
#                     end
#                     subcoeff = ones(F, length(subgraphs))
#                     dual[d.id] =linear_combination(subgraphs,subcoeff)
#                     dual[d.id].factor *= d.factor*factor
#                 elseif !isnothing(subnumber)  #if only numbers appear in derivative, return a number
#                     dual[d.id] = subnumber*d.factor*factor
#                 end
#             else
#                 error("not implemented!")
#             end
#         end
#     end
#     if isempty(dual)
#         return 0.0
#     end
#     return dual[rootid]
# end



# """
#     function derivative(diags::Union{Diagram,Tuple,AbstractVector}, ::Type{ID}, order::Int) where {ID<:PropagatorId}

#     Automatic differentiation derivative on the diagrams

# # Arguments
# - diags      : diagrams to take derivative
# - ID         : DiagramId to apply the differentiation
# - order::Int : derivative order
# """
# function derivative(diags::Union{Tuple,AbstractVector}, ::Type{ID}, order::Int; index::Int=index(ID)) where {ID<:PropagatorId}
#     @assert order >= 0
#     if order == 0
#         return diags
#     end
#     result = diags
#     for o in 1:order
#         result = derivative(result, ID; index=index)
#     end
#     return result
# end

# """
#     function removeHartreeFock!(diag::Graph{F,W}) where {F,W}
#     function removeHartreeFock!(diags::Union{Tuple,AbstractVector})

#     Remove the Hartree-Fock insertions that without any derivatives on the propagator and the interaction. 

# # Arguments
# - diags      : diagrams to remove the Fock insertion

# # Remarks
# - The operations removeHartreeFock! and taking derivatives doesn't commute with each other! 
# - If the input diagram is a Hartree-Fock diagram, then the overall weight will become zero! 
# - The return value is always nothing
# """
# function removeHartreeFock!(diag::Graph{F,W}) where {F,W}
#     for d in PreOrderDFS(diag)
#         # for subdiag in d.subdiagram
#         if d.id isa SigmaId
#             # println(d, " with ", d.id.para.innerLoopNum)
#             if isempty(d.id.order) || all(x -> x == 0, d.id.order) #eithr order is empty or a vector of zeros
#                 if d.id.para.innerLoopNum == 1
#                     d.factor = 0.0
#                 end
#             end
#         end
#     end
# end
# function removeHartreeFock!(diags::Union{Tuple,AbstractVector})
#     for diag in diags
#         removeHartreeFock!(diag)
#     end
# end

# """
#     function MatsubaraSum(diags::Vector{Graph{F,W}}) where {F,W}

#     Automatic MatsubaraSum on the diagrams

# # Arguments
# - diags     : diagrams to take MatsubaraSum
# """
# function MatsubaraSum(diags::Vector{Graph{F,W}}) where {F,W}
# end
# # function MatsubaraSum(diags::Union{Tuple,AbstractVector}) where {F,W}
# # end
