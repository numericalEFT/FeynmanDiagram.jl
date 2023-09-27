"""
    function linear_combination_number_with_graph(g::Vector{Union{F,Graph{F,W}}},  coeff::Vector{F}=ones(F, length(g))) where {F,W}

    Given a vector 𝐠 of mixed numbers and graphs,  and an equally-sized vector 𝐜 of constants (by default are all unity), and return a new
    graph representing the linear combination (𝐜 ⋅ 𝐠). All numbers are treated as graphs, unless all members of g are numbers, in which case simply return (𝐜 ⋅ 𝐠) as a number.
# Arguments:
- `g::Vector{Union{F,Graph{F,W}}}`:  Vector of numbers and graphs  𝐠 to be combined
- `coeff::Vector{F}`: coefficient vector  𝐜
"""

function linear_combination_number_with_graph(g::Vector{Union{F,Graph{F,W}}}, coeff::Vector{F}=ones(F, length(g))) where {F,W}
    @assert length(g) == length(coeff)
    subgraphs = Vector{Graph{F,W}}()
    subnumber = nothing
    subcoeff = Vector{F}()
    result = nothing
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

    if isempty(subgraphs) == false
        if !isnothing(subnumber)
            push!(subgraphs, constant_graph(F(subnumber)))     #If both numbers and graphs appear in derivative, convert number to a unity graph, and asign the number to subgraph_factors of parent node.
            push!(subcoeff, 1.0)
        end
        result = linear_combination(subgraphs, subcoeff)
        #result.factor *= d.factor
    elseif !isnothing(subnumber)  #if only numbers appear in derivative, return a number
        result = subnumber #* d.factor
    end
    return result
    #return subgraphs, subnumber, subcoeff
end
"""
    function forwardAD(diag::Graph{F,W}, ID::Int) where {F,W}
    
    Given a  graph G and an id of graph g, calculate the derivative d G / d g by forward propagation algorithm.
# Arguments:
- `diag::Graph{F,W}`:  Graph to be differentiated
- `ID::Int`: ID of the graph that is  
"""
function forwardAD(diag::Graph{F,W}, ID::Int) where {F,W}
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
                dual[d.id] = F(1)
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
                dum = linear_combination_number_with_graph(children, coeff)
                if !isnothing(dum)
                    dual[d.id] = dum
                end
                # subgraphs, subnumber, subcoeff = linear_combination_number_with_graph(children, coeff)
                # if isempty(subgraphs) == false
                #     if !isnothing(subnumber)
                #         push!(subgraphs, constant_graph(F(subnumber)))     #If both numbers and graphs appear in derivative, convert number to a unity graph, and asign the number to subgraph_factors of parent node.
                #         push!(subcoeff, 1.0)
                #     end
                #     dual[d.id] = linear_combination(subgraphs, subcoeff)
                #     dual[d.id].factor *= d.factor
                # elseif !isnothing(subnumber)  #if only numbers appear in derivative, return a number
                #     dual[d.id] = subnumber * d.factor
                # end
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
                dum = linear_combination_number_with_graph(children)
                if !isnothing(dum)
                    dual[d.id] = factor * dum
                end

                #subgraphs, subnumber = linear_combination_number_with_graph(children)

                # if isempty(subgraphs) == false
                #     if !isnothing(subnumber)
                #         push!(subgraphs, constant_graph(F(subnumber)))     #If both numbers and graphs appear in derivative, convert number to a constant graph, and asign the number to subgraph_factors of parent node.
                #     end
                #     subcoeff = ones(F, length(subgraphs))
                #     dual[d.id] = linear_combination(subgraphs, subcoeff)
                #     dual[d.id].factor *= d.factor * factor
                # elseif !isnothing(subnumber)  #if only numbers appear in derivative, return a number
                #     dual[d.id] = subnumber * d.factor * factor
                # end
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

function all_parent(diag::Graph{F,W}) where {F,W}
    result = Dict{Int,Vector{Graph{F,W}}}()
    for d in PostOrderDFS(diag)
        if !haskey(result, d.id)
            parents = Vector{Graph{F,W}}()
            for g in PostOrderDFS(diag)
                if d.id in [sub.id for sub in g.subgraphs]
                    push!(parents, g)
                end
            end
            result[d.id] = parents
        end
    end
    return result
end

function node_derivative(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W} #return d g1/ d g2 
    if isleaf(g1)
        return nothing
    elseif g1.operator == Sum
        sum_factor = 0.0
        exist = false #track if g2 exist in g1 subgraphs.
        for i in 1:length(g1.subgraphs)
            if g1.subgraphs[i].id == g2.id
                exist = true
                sum_factor += g1.subgraph_factors[i]
            end
        end
        if exist
            return F(sum_factor)
        else
            return nothing
        end
    elseif g1.operator == Prod
        count = 0 # count how many times g2 appears in g1, i.e. power of g2 in g1
        subgraphs = []
        subgraphfactors = []
        factor = nothing
        first_time = true #Track if its the first time we find g2 in g1 subgraph.
        for i in 1:length(g1.subgraphs)
            if g1.subgraphs[i].id == g2.id
                if first_time # We should remove the first g2 in g1
                    first_time = false
                    factor = g1.subgraph_factors[i] #save the factor of first g2.
                    count += 1
                else
                    count += 1
                    push!(subgraphs, g1.subgraphs[i]) #Keep other g2 in g1
                    push!(subgraphfactors, g1.subgraph_factors[i])
                end
            else
                push!(subgraphs, g1.subgraphs[i])
                push!(subgraphfactors, g1.subgraph_factors[i])
            end
        end
        if count == 0
            return nothing
        end
        if isempty(subgraphs)
            return factor
        else
            if !isnothing(factor)
                subgraphfactors[1] *= count * factor
            end
            g = deepcopy(g1)
            g.subgraphs = subgraphs
            g.subgraph_factors = subgraphfactors
            return g
        end
    else
        return nothing
    end
end

function recursive_backAD!(diag::Graph{F,W}, parents::Dict{Int,Vector{Graph{F,W}}}, dual::Dict{Int,Union{F,Graph{F,W}}}, result::Dict{Tuple{Int,Int},Graph{F,W}}, root_id::Int) where {F,W}
    if !haskey(dual, diag.id)
        derivative_list = Vector{Union{F,Graph{F,W}}}()
        if isempty(parents[diag.id]) # if d is the root node, the backward derivative is just 1.
            dual[diag.id] = F(1)
        else
            for parent in parents[diag.id]
                parent_AD = recursive_backAD!(parent, parents, dual, result, root_id)
                d_node = node_derivative(parent, diag)
                if !isnothing(d_node) && !isnothing(parent_AD)
                    push!(derivative_list, d_node * parent_AD)
                end
            end
            dum = linear_combination_number_with_graph(derivative_list)
            dual[diag.id] = dum
        end
    end
    if isleaf(diag) # The final result is derivative with respect to each leaf 
        if typeof(dual[diag.id]) <: Number
            result[(root_id, diag.id)] = constant_graph(F(dual[diag.id]))
        else
            result[(root_id, diag.id)] = dual[diag.id]
        end
    end
    return dual[diag.id]
end

function backAD(diag::Graph{F,W}, debug::Bool=false) where {F,W}
    # use a dictionary to host the dual diagram of a diagram for a given hash number
    # a dual diagram is defined as the derivative of the original diagram
    dual = Dict{Int,Union{F,Graph{F,W}}}()
    result = Dict{Tuple{Int,Int},Graph{F,W}}()
    parents = all_parent(diag)
    for d in Leaves(diag)#PreOrderDFS(diag) # preorder traversal will visit all parents first
        #print("type: $(d.type)\n")
        if d.type == Constant || haskey(dual, d.id)
            continue
        end
        recursive_backAD!(d, parents, dual, result, diag.id)
    end
    return result
end

# function derivative_variable!(diag::Graph{F,W}, derivative_mapping::Dict{Int, Tuple{Int,Vector{Int}}} , variable_idx::Int) where {F,W}
#     children = Vector{Graph{F,W}}()
#     for (id, AD) in backAD(diag)
#         if haskey((id, ))
#         new_node = Graph([], ftype=F, wtype=W) # Create the new green's function that represents d g/ d x, where x is the variable, and g is the current leaf.
#         if !isnothing(AD)
#             push!(children, new_node * AD)
#         end
#     end
#     if isempty(children)
#         return nothing
#     else
#         return linear_combination(children, ones(F, length(children)))
#     end
# end

function build_all_leaf_derivative(diag::Graph{F,W}, max_order::Int) where {F,W}
    result = Dict{Vector{Int},Graph{F,W}}()
    current_func = backAD(diag)
    leave_number = length(current_func)
    order_dict = Dict{Int,Vector{Int}}()

    leafmap = Dict{Int,Int}()
    count = 1
    for (id_pair, func) in current_func
        leafmap[id_pair[2]] = count
        order = zeros(Int, leave_number)
        order[count] += 1
        order_dict[func.id] = order
        result[order] = func
        count += 1
    end

    for i in 2:max_order
        new_func = Dict{Tuple{Int,Int},Graph{F,W}}()
        for (id_pair, func) in current_func
            AD = backAD(func)
            print("$(i) $(AD)\n")
            for (id_pair_AD, func_AD) in AD
                # print("$(func_AD)\n")
                order = copy(order_dict[func.id])
                order[leafmap[id_pair_AD[2]]] += 1

                if !(order in values(order_dict))
                    new_func[id_pair_AD] = func_AD
                    order_dict[func_AD.id] = order
                    result[order] = func_AD
                end
            end
        end
        current_func = new_func
    end
    return result
end

# function build_all_variable_derivative(diag::Graph{F,W}, max_order::Int, variable_number::Int) where {F,W}
#     leaf_derivative = build_all_leaf_derivative(diag, max_order)


# end

# function backAD(diag::Graph{F,W}, debug::Bool=false) where {F,W}
#     # use a dictionary to host the dual diagram of a diagram for a given hash number
#     # a dual diagram is defined as the derivative of the original diagram
#     dual = Dict{Int,Union{F,Graph{F,W}}}()
#     result = Dict{Int,Union{F,Graph{F,W}}}()
#     parents = all_parent(diag)
#     for d in reverse([node for node in PostOrderDFS(diag)]) # preorder traversal will visit all parents first
#         if debug
#             print("node: $(d.id) $(d.subgraphs) $(d.subgraph_factors)\n")
#         end
#         if haskey(dual, d.id)
#             continue
#         end

#         if isempty(parents[d.id]) # if d is the root node, the backward derivative is just 1.
#             dual[d.id] = F(1)
#         else
#             derivative_list = Vector{Union{F,Graph{F,W}}}()
#             for parent in parents[d.id]
#                 if debug
#                     print("\t parent: $(parent.id) $(parent.subgraphs) $(parent.subgraph_factors)\n")
#                 end
#                 if haskey(dual, parent.id)
#                     d_node = node_derivative(parent, d)
#                     if !isnothing(d_node)
#                         push!(derivative_list, d_node * dual[parent.id])
#                     end
#                     if debug
#                         if isnothing(d_node) || typeof(d_node) <: Number
#                             print("\tderivative: $(d_node)\n")
#                         else
#                             print("\tderivative: $(d_node.id) $(d_node.subgraphs) $(d_node.subgraph_factors)\n")
#                         end
#                         if isnothing(dual[parent.id]) || typeof(dual[parent.id]) <: Number
#                             print("\tparent derivative: $(dual[parent.id])\n")
#                         else
#                             print("\tparent derivative: $(dual[parent.id].id) $(dual[parent.id].subgraphs) $(dual[parent.id].subgraph_factors)\n")
#                         end
#                     end
#                 end
#             end
#             dum = linear_combination_number_with_graph(derivative_list)
#             if !isnothing(dum)
#                 dual[d.id] = dum
#             end
#         end
#         if isleaf(d) && haskey(dual, d.id) # The final result is derivative with respect to each leaf 
#             result[d.id] = dual[d.id]
#         end
#     end
#     return result
# end


function forwardAD_root(diags::AbstractVector{G}) where {G<:Graph}
    dual = Dict{Int,G}()
    # println("rootID: ", diag.id)
    for diag in diags
        for node in PreOrderDFS(diag)
            visited = false
            if haskey(dual, node.id)
                dual[node.id].name != "None" && continue
                visited = true
            end
            # println("Node: ", node.id)

            if node.operator == Sum
                nodes_deriv = G[]
                for sub_node in node.subgraphs
                    if haskey(dual, sub_node.id)
                        # println("subNode haskey: ", sub_node.id)
                        push!(nodes_deriv, dual[sub_node.id])
                    else
                        # println("subNode nokey: ", sub_node.id)
                        g_dual = Graph(G[]; name="None")
                        push!(nodes_deriv, g_dual)
                        dual[sub_node.id] = g_dual
                    end
                end
                if visited
                    dual[node.id].subgraphs = nodes_deriv
                    dual[node.id].subgraph_factors = node.subgraph_factors
                    dual[node.id].name = node.name
                else
                    dual[node.id] = Graph(nodes_deriv; subgraph_factors=node.subgraph_factors, factor=node.factor)
                end
            elseif node.operator == Prod
                nodes_deriv = G[]
                for (i, sub_node) in enumerate(node.subgraphs)
                    if haskey(dual, sub_node.id)
                        # println("subNode haskey: ", sub_node.id)
                        subgraphs = [j == i ? dual[subg.id] : subg for (j, subg) in enumerate(node.subgraphs)]
                        push!(nodes_deriv, Graph(subgraphs; operator=Prod(), subgraph_factors=node.subgraph_factors))
                    else
                        # println("subNode nokey: ", sub_node.id)
                        g_dual = Graph(G[]; name="None")
                        dual[sub_node.id] = g_dual
                        subgraphs = [j == i ? g_dual : subg for (j, subg) in enumerate(node.subgraphs)]
                        push!(nodes_deriv, Graph(subgraphs; operator=Prod(), subgraph_factors=node.subgraph_factors))

                    end
                end
                if visited
                    dual[node.id].subgraphs = nodes_deriv
                    dual[node.id].subgraph_factors = one.(eachindex(nodes_deriv))
                    dual[node.id].name = node.name
                else
                    dual[node.id] = Graph(nodes_deriv; factor=node.factor)
                end
            end
        end
    end
    return dual
end

forwardAD_root(diag::Graph) = forwardAD_root([diag])