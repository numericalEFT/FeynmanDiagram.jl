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
- `diag::Graph{F,W}`:  Graph G to be differentiated
- `ID::Int`: ID of the graph g  
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

"""
    function all_parent(diag::Graph{F,W}) where {F,W}
    
    Given a  graph, find all parents node of each node, and return them in a dictionary.
# Arguments:
- `diag::Graph{F,W}`:  Target graph
"""

function all_parent(diag::Graph{F,W}) where {F,W}
    result = Dict{Int,Vector{Graph{F,W}}}()
    for d in PostOrderDFS(diag)
        if !haskey(result, d.id)
            parents = Vector{Graph{F,W}}()
            parents_id = Vector{Int}()
            for g in PostOrderDFS(diag)
                if !(g.id in parents_id) && d.id in [sub.id for sub in g.subgraphs]
                    push!(parents, g)
                    push!(parents_id, g.id)
                end
            end
            result[d.id] = parents
        end
    end
    return result
end

"""
function node_derivative(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W} 
    
    Return the local derivative d g1/ dg2 at node g1. The local derivative only considers the subgraph of node g1, and ignores g2 that appears in deeper layers.
    Example: For g1 = G *g2, and G = g3*g2,  return d g1/ dg2 = G = g3*g2 instead of 2 g2*g3. 
# Arguments:
- `diag::Graph{F,W}`:  Target graph
"""

function node_derivative(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W} #return d g1/ d g2 
    if isleaf(g1)
        return nothing
    elseif g1.operator == Sum
        sum_factor = 0.0
        exist = false #track if g2 exist in g1 subgraphs.
        for i in eachindex(g1.subgraphs)
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
        for i in eachindex(g1.subgraphs)
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
        if d.operator == Constant || haskey(dual, d.id)
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
#     leaf_derivative, leafmap = build_all_leaf_derivative(diag, max_order)
#     for (id, idx) in leafmap
#         for order in max_order
#             for 
#         end
#     end
# end

function insert_dualDict!(dict_kv::Dict{Tk,Tv}, dict_vk::Dict{Tv,Tk}, key::Tk, value::Tv) where {Tk,Tv}
    dict_kv[key] = value
    if !haskey(dict_vk, value)
        dict_vk[value] = Set()
    end
    push!(dict_vk[value], key)
end


"""
    function forwardAD_root!(graphs::AbstractVector{G}, idx::Int=1,
        dual::Dict{Tuple{Int,NTuple{N,Bool}},G}=Dict{Tuple{Int,Tuple{Bool}},G}()) where {G<:Graph,N}
 
    Computes the forward automatic differentiation (AD) of the given graphs beginning from the roots.

# Arguments:
- `graphs`: A vector of graphs.
- `idx`: Index for differentiation (default: 1).
- `dual`: A dictionary that holds the result of differentiation.

# Returns:
- The dual dictionary populated with all differentiated graphs, including the intermediate AD.
"""
function forwardAD_root!(graphs::AbstractVector{G}, idx::Int=1,
    dual::Dict{Tuple{Int,NTuple{N,Bool}},G}=Dict{Tuple{Int,Tuple{Bool}},G}()) where {G<:Graph,N}
    # dualinv::Dict{G,Tuple{Int,NTuple{N,Int}}}=Dict{G,Tuple{Int,Tuple{Int}}}()) where {G<:Graph,N}

    @assert idx <= N "the differential variable's index must be described by the key of dual."
    # dual = Dict{Int,G}()
    # println("rootID: ", diag.id)
    key2 = Tuple(i == idx ? true : false for i in 1:N)
    for diag in graphs
        for node in PreOrderDFS(diag)
            visited = false
            # if haskey(dual, node.id)
            if any(key[1] == node.id && key[2] == key2 for key in keys(dual))
                dual[(node.id, key2)].name != "None" && continue
                visited = true
            end
            # println("Node: ", node.id)

            if node.operator == Sum
                nodes_deriv = G[]
                for sub_node in node.subgraphs
                    key = (sub_node.id, key2)
                    if haskey(dual, key)
                        # println("subNode haskey: ", sub_node.id)
                        push!(nodes_deriv, dual[key])
                    else
                        # println("subNode nokey: ", sub_node.id)
                        g_dual = Graph(G[]; name="None")
                        push!(nodes_deriv, g_dual)
                        dual[key] = g_dual
                    end
                end
                key_node = (node.id, key2)
                if visited
                    dual[key_node].subgraphs = nodes_deriv
                    dual[key_node].subgraph_factors = node.subgraph_factors
                    dual[key_node].name = node.name
                else
                    dual[key_node] = Graph(nodes_deriv; subgraph_factors=node.subgraph_factors, factor=node.factor)
                end
            elseif node.operator == Prod
                nodes_deriv = G[]
                for (i, sub_node) in enumerate(node.subgraphs)
                    key = (sub_node.id, key2)
                    if haskey(dual, key)
                        # println("subNode haskey: ", sub_node.id)
                        subgraphs = [j == i ? dual[key] : subg for (j, subg) in enumerate(node.subgraphs)]
                        push!(nodes_deriv, Graph(subgraphs; operator=Prod(), subgraph_factors=node.subgraph_factors))
                    else
                        # println("subNode nokey: ", sub_node.id)
                        g_dual = Graph(G[]; name="None")
                        dual[key] = g_dual
                        subgraphs = [j == i ? g_dual : subg for (j, subg) in enumerate(node.subgraphs)]
                        push!(nodes_deriv, Graph(subgraphs; operator=Prod(), subgraph_factors=node.subgraph_factors))

                    end
                end
                key_node = (node.id, key2)
                if visited
                    dual[key_node].subgraphs = nodes_deriv
                    dual[key_node].subgraph_factors = one.(eachindex(nodes_deriv))
                    dual[key_node].name = node.name
                else
                    dual[key_node] = Graph(nodes_deriv; factor=node.factor)
                end
            end
        end
    end
    return dual
end

function forwardAD_root!(diag::G, idx::Int=1,
    dual::Dict{Tuple{Int,NTuple{N,Bool}},G}=Dict{Tuple{Int,Tuple{Bool}},G}()) where {G<:Graph,N}
    return forwardAD_root!([diag], idx, dual)
end


@inline function find_last_neighbor(item)
    # loc = findfirst(val -> val > 0, item)
    loc = findlast(val -> val > 0, item)
    if isnothing(loc)
        return nothing
    else
        return ntuple(j -> j == loc ? item[j] - 1 : item[j], length(item))
    end
end

"""
    function build_derivative_graph(graphs::AbstractVector{G}, orders::NTuple{N,Int};
        nodes_id=nothing) where {G<:Graph,N}
 
    Constructs a derivative graph using forward automatic differentiation with given graphs and derivative orders.

# Arguments:
- `graphs`: A vector of graphs.
- `orders::NTuple{N,Int}`: A tuple indicating the orders of differentiation. `N` represents the number of independent variables to be differentiated.
- `nodes_id`: Optional node IDs to indicate saving their derivative graph.

# Returns:
- A dictionary containing the dual derivative graphs for all indicated nodes. 
If `isnothing(nodes_id)`, indicated nodes include all leaf and root nodes. Otherwise, indicated nodes include all root nodes and other nodes from `nodes_id`.
"""
function build_derivative_graph(graphs::AbstractVector{G}, orders::NTuple{N,Int};
    nodes_id=nothing) where {G<:Graph,N}

    roots_id = Set{Int}()
    for g in graphs
        push!(roots_id, g.id)
    end
    if isnothing(nodes_id)
        nodes_id = Set{Int}()
        for g in graphs
            for leaf in Leaves(g)
                push!(nodes_id, leaf.id)
            end
        end
    end

    dual_oneorder = Dict{Tuple{Int,NTuple{N,Bool}},G}()
    cumsum_orders = cumsum(orders)
    idx0 = findfirst(val -> 1 <= val, cumsum_orders)
    first_order = ntuple(j -> j == idx0 ? 1 : 0, N)

    # generate dual 1-order derivative graph 
    dual_oneorder = forwardAD_root!(graphs, idx0, dual_oneorder)
    dual_graphs = [dual_oneorder[(g.id, first_order)] for g in graphs]
    for x in 2:sum(orders)
        idx = findfirst(val -> x <= val, cumsum_orders)
        dual_oneorder = forwardAD_root!(dual_graphs, idx, dual_oneorder)
        dual_graphs = [dual_oneorder[(g.id, ntuple(j -> j == idx, N))] for g in dual_graphs]
    end

    dual = Dict{Tuple{Int,NTuple{N,Int}},G}()
    # generate dual derivative graph Dict for all nodes (except for root nodes)
    iter_orders = Tuple(0:x for x in orders)
    for node_id in nodes_id
        for order in Iterators.product(iter_orders...)
            order == ntuple(j -> 0, N) && continue
            prev_order = find_last_neighbor(order)
            if prev_order == ntuple(j -> 0, N)
                dual[(node_id, order)] = dual_oneorder[(node_id, prev_order .!= order)]
            else
                dual[(node_id, order)] = dual_oneorder[(dual[(node_id, prev_order)].id, prev_order .!= order)]
            end
        end
    end

    # generate dual derivative graph Dict for all root nodes
    _cum = pushfirst!(collect(cumsum_orders), 0)
    for root_id in roots_id
        dual[(root_id, first_order)] = dual_oneorder[(root_id, first_order)]
        prev_order = first_order
        for x in 2:sum(orders)
            idx = findfirst(val -> x <= val, cumsum_orders)
            order = ntuple(j -> j == idx ? x - _cum[idx] : (j < idx ? orders[j] : 0), N)
            dual[(root_id, order)] = dual_oneorder[(dual[(root_id, prev_order)].id, prev_order .!= order)]
            prev_order = order
        end
    end

    return dual
end

function build_derivative_graph(graph::G, orders::NTuple{N,Int};
    nodes_id=nothing) where {G<:Graph,N}
    return build_derivative_graph([graph], orders; nodes_id=nodes_id)
end