
# 1. 定义叶子的唯一特征 Key
struct LeafKey
    operator::DataType
    orders::Vector{Int}
    properties::Any
end

function Base.hash(v::LeafKey, h::UInt)
    h = hash(v.operator, h)
    h = hash(v.orders, h)
    h = hash(v.properties, h)
    # h = hash(string(v.properties), h)
    return h
end

function Base.isequal(a::LeafKey, b::LeafKey)
    if a.operator != b.operator || a.orders != b.orders
        return false
    end
    # return string(a.properties) == string(b.properties)
    return isequal(a.properties, b.properties)
end

"""
    optimize_randomized!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; 
                         digits=10, verbose=0, normalize=nothing)

基于随机化求值的高效去重优化。
1. 首先去除重复 Leaf。
2. 给所有 Unique Leaf 赋予随机复数权值。
3. 自底向上计算图的“随机指纹”。
4. 根据指纹合并重复的中间节点。
"""
function optimize!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}};
    digits=10, verbose=0, normalize=nothing, seed=1234)
    verbose > 0 && println("Starting randomized optimization...")

    # 1. 先用精确哈希去重 Leaf (复用之前的优化函数)
    #    这是必须的，保证同一个物理变量对应同一个随机值
    remove_duplicated_leaves!(graphs, verbose=verbose, normalize=normalize)

    # 2. 为每个 Unique Leaf 分配随机值
    #    使用 Dict 存储: Leaf ID -> Random Value
    rng = MersenneTwister(seed) # 固定种子以保证可复现性
    leaf_vals = Dict{Int,ComplexF64}()

    # 收集所有 graphs 中的 unique leaf 并赋值
    # 这里不需要显式收集，我们在遍历时动态赋值即可，或者先遍历一遍收集
    # 为了简单，我们在 recursive_eval 中动态查表

    # 3. 缓存表
    # val_map: 节点 ID -> 计算出的随机值 (防止 DAG 重复计算)
    val_map = Dict{Int,ComplexF64}()

    # unique_table: (Operator, RoundedValue) -> 规范化节点
    # Key 使用 Tuple{DataType, ComplexF64}，其中 Value 经过 round 处理
    unique_table = Dict{Any,eltype(graphs)}()

    # 访问记录，用于构建 DAG
    memo_node = Dict{Int,eltype(graphs)}()

    # --- 核心递归函数 ---
    function recursive_build(g::AbstractGraph)
        # 1. 记忆化：如果已经处理过该节点结构，直接返回规范化节点
        if haskey(memo_node, g.id)
            return memo_node[g.id]
        end

        # 2. 计算当前节点的随机值 (Value)
        current_val = zero(ComplexF64)

        if isleaf(g)
            # 如果是 Leaf，获取或分配随机值
            if !haskey(leaf_vals, g.id)
                # 生成一个随机复数，实部虚部都在 [0, 100) 之间避免太小
                leaf_vals[g.id] = complex(rand(rng) * 100, rand(rng) * 100)
            end
            current_val = leaf_vals[g.id]
        else
            # 如果是 Branch，递归处理子节点
            # 先递归优化子节点，确保存储在 g.subgraphs 中的已经是规范化节点
            for (i, sub_g) in enumerate(subgraphs(g))
                canonical_sub = recursive_build(sub_g)
                set_subgraph!(g, canonical_sub, i)

                sub_val = val_map[canonical_sub.id]
                factor = subgraph_factor(g, i)

                current_val += sub_val * factor
            end
        end

        # 记录计算出的值
        val_map[g.id] = current_val

        # 3. 查表去重 (Dedup)
        # 将复数 round 一下，作为 Key
        # round 的目的是消除 1.000000000000001 和 1.0 的区别
        rounded_val = round(current_val, digits=digits)

        # 构造指纹 Key: (Operator, 随机计算值, Orders)
        # 注意：Orders 也要相同才算相等
        dedup_key = (g.operator, rounded_val, g.orders)

        if haskey(unique_table, dedup_key)
            canonical_node = unique_table[dedup_key]
        else
            unique_table[dedup_key] = g
            canonical_node = g
        end

        memo_node[g.id] = canonical_node
        return canonical_node
    end

    if graphs isa AbstractVector
        for (i, g) in enumerate(graphs)
            graphs[i] = recursive_build(g)
        end
    else
        for g in graphs
            recursive_build(g)
        end
    end

    return graphs
end

"""
    function optimize!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; level=0, verbose=0, normalize=nothing)

    In-place optimization of given `graphs`. Removes duplicated leaves, flattens chains, 
    merges linear combinations, and removes zero-valued subgraphs. When `level > 0`, also removes duplicated intermediate nodes.

# Arguments:
- `graphs`: A tuple or vector of graphs.
- `level`: Optimization level (default: 0). A value greater than 0 triggers more extensive but slower optimization processes, such as removing duplicated intermediate nodes.
- `verbose`: Level of verbosity (default: 0).
- `normalize`: Optional function to normalize the graphs (default: nothing).

# Returns
- Returns the optimized graphs. If the input graphs is empty, it returns nothing.
"""
# function optimize!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; level=0, verbose=0, normalize=nothing)
#     if isempty(graphs)
#         return nothing
#     else
#         if level > 0
#             if graphs isa Tuple
#                 root = Graph(collect(graphs))
#             else
#                 root = Graph(graphs)
#             end
#             remove_duplicated_nodes!(root, verbose=verbose)
#         else
#             remove_duplicated_leaves!(graphs, verbose=verbose, normalize=normalize)
#         end

#         flatten_all_chains!(graphs, verbose=verbose)
#         merge_all_linear_combinations!(graphs, verbose=verbose)
#         remove_all_zero_valued_subgraphs!(graphs, verbose=verbose)
#         return graphs
#     end
# end

"""
    function optimize!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; level=0, verbose=0, normalize=nothing)

    In-place optimization of given `graphs`. 
    Uses a fused post-order traversal with memoization to perform flattening, merging, and zero-removal in a single pass.
"""
function optimize_v0!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; level=0, verbose=0, normalize=nothing)
    if isempty(graphs)
        return nothing
    end

    # 1. Structural Simplification (Global)
    if level > 0
        if graphs isa Tuple
            root = Graph(collect(graphs))
        else
            root = Graph(graphs)
        end
        remove_duplicated_nodes!(root, verbose=verbose)
    else
        remove_duplicated_leaves!(graphs, verbose=verbose, normalize=normalize)
    end

    println("structural simplification done.")

    # 2. Algebraic Simplification (Local, Fused)
    visited = Set{Int}()

    for g in graphs
        recursive_optimize!(g, visited)
    end

    return graphs
end

"""
    recursive_optimize!(g::AbstractGraph, visited::Set{Int})

    Helper function to perform fused optimization operations in post-order.
"""
function recursive_optimize!(g::AbstractGraph, visited::Set{Int})
    if g.id in visited
        return
    end

    # Post-order DFS
    for sub_g in subgraphs(g)
        recursive_optimize!(sub_g, visited)
    end

    # Apply operations bottom-up
    # 1. 先移除零值子图，可能会使当前节点变成零或简化结构。
    remove_zero_valued_subgraphs!(g)

    # 2. 扁平化链条。如果去零后导致出现了单链 (e.g., Sum(0, g1) -> Sum(g1) -> g1)，这里可以进一步简化。
    # flatten_chains!(g)
    _flatten_chains_shallow!(g)

    # 3. 合并线性组合。这通常是在同一层级上操作，应该在子图稳定后进行。
    merge_linear_combination!(g)

    # 标记当前节点为已处理
    push!(visited, g.id)
end

"""
    _flatten_chains_shallow!(g::AbstractGraph)

flatten_chains! 的非递归版本。仅检查直接子节点是否为 trivial unary 链。
"""
function _flatten_chains_shallow!(g::AbstractGraph)
    # 遍历当前子节点
    # 注意：直接修改 g.subgraphs 可能影响迭代，但在 enumerate 下，
    # 只要我们是一对一替换 (sub_g -> child)，索引是安全的。
    for (i, sub_g) in enumerate(subgraphs(g))
        # 检查子节点是否是 "单传" 节点 (Trivial Unary)
        # 且该子节点已经在 recursive_optimize! 中被处理过，是最简形式
        if unary_istrivial(sub_g) && onechild(sub_g)
            # 提升孙子节点 (Grandchild) 替换子节点
            child = sub_g.subgraphs[1]
            factor = sub_g.subgraph_factors[1]

            # 原地修改图连接
            set_subgraph!(g, child, i)

            # 更新因子: parent_factor * sub_factor
            new_factor = subgraph_factor(g, i) * factor
            set_subgraph_factor!(g, new_factor, i)
        end
    end
    return g
end

"""
    function optimize(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; level=0, verbose=0, normalize=nothing)

    Optimizes a copy of given `graphs`. Removes duplicated nodes (when `level > 0`) or leaves, flattens chains, 
    merges linear combinations, and removing zero-valued subgraphs.

# Arguments:
- `graphs`: A tuple or vector of graphs.
- `level`: Optimization level (default: 0). A value greater than 0 triggers more extensive but slower optimization processes, such as removing duplicated nodes.
- `verbose`: Level of verbosity (default: 0).
- `normalize`: Optional function to normalize the graphs (default: nothing).

# Returns:
- A tuple/vector of optimized graphs.
"""
function optimize(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0, normalize=nothing)
    graphs_new = deepcopy(graphs)
    # optimize!(graphs_new, level=level, verbose=verbose, normalize=normalize)
    optimize!(graphs_new, verbose=verbose, normalize=normalize)
    return graphs_new
end

"""
    function flatten_all_chains!(g::AbstractGraph; verbose=0)
F
    Flattens all nodes representing trivial unary chains in-place in the given graph `g`. 

# Arguments:
- `graphs`: The graph to be processed.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- The mutated graph `g` with all chains flattened.
"""
function flatten_all_chains!(g::AbstractGraph; verbose=0)
    verbose > 0 && println("flatten all nodes representing trivial unary chains.")
    for sub_g in g.subgraphs
        flatten_all_chains!(sub_g)
        flatten_chains!(sub_g)
    end
    flatten_chains!(g)
    return g
end

"""
    function flatten_all_chains!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0)

    Flattens all nodes representing trivial unary chains in-place in the given graphs.

# Arguments:
- `graphs`: A collection of graphs to be processed.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- The mutated collection `graphs` with all chains in each graph flattened.
"""
function flatten_all_chains!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0)
    verbose > 0 && println("flatten all nodes representing trivial unary chains.")
    # Post-order DFS
    for g in graphs
        flatten_all_chains!(g.subgraphs)
        flatten_chains!(g)
    end
    return graphs
end

"""
    function remove_all_zero_valued_subgraphs!(g::AbstractGraph; verbose=0)

    Recursively removes all zero-valued subgraph(s) in-place in the given graph `g`.

# Arguments:
- `g`: An AbstractGraph.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graph.
# 
"""
function remove_all_zero_valued_subgraphs!(g::AbstractGraph; verbose=0)
    verbose > 0 && println("merge nodes representing a linear combination of a non-unique list of graphs.")
    # Post-order DFS
    for sub_g in subgraphs(g)
        remove_all_zero_valued_subgraphs!(sub_g)
        remove_zero_valued_subgraphs!(sub_g)
    end
    remove_zero_valued_subgraphs!(g)
    return g
end

"""
    function remove_all_zero_valued_subgraphs!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0)

    Recursively removes all zero-valued subgraph(s) in-place in the given graphs.

# Arguments:
- `graphs`: A collection of graphs to be processed.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graphs.
# 
"""
function remove_all_zero_valued_subgraphs!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0)
    verbose > 0 && println("merge nodes representing a linear combination of a non-unique list of graphs.")
    # Post-order DFS
    for g in graphs
        remove_all_zero_valued_subgraphs!(subgraphs(g))
        remove_zero_valued_subgraphs!(g)
    end
    return graphs
end

"""
    function merge_all_linear_combinations!(g::AbstractGraph; verbose=0)

    Merges all nodes representing a linear combination of a non-unique list of subgraphs in-place in the given graph `g`.

# Arguments:
- `g`: An AbstractGraph.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graph.
# 
"""
function merge_all_linear_combinations!(g::AbstractGraph; verbose=0)
    verbose > 0 && println("merge nodes representing a linear combination of a non-unique list of graphs.")
    # Post-order DFS
    for sub_g in subgraphs(g)
        merge_all_linear_combinations!(sub_g)
        merge_linear_combination!(sub_g)
    end
    merge_linear_combination!(g)
    return g
end

"""
    function merge_all_linear_combinations!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0)

    Merges all nodes representing a linear combination of a non-unique list of subgraphs in-place in the given graphs. 

# Arguments:
- `graphs`: A collection of graphs to be processed.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graphs.
# 
"""
function merge_all_linear_combinations!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0)
    verbose > 0 && println("merge nodes representing a linear combination of a non-unique list of graphs.")
    # Post-order DFS
    for g in graphs
        merge_all_linear_combinations!(subgraphs(g))
        merge_linear_combination!(g)
    end
    return graphs
end

"""
    function merge_all_multi_products!(g::Graph; verbose=0)

    Merges all nodes representing a multi product of a non-unique list of subgraphs in-place in the given graph `g`.

# Arguments:
- `g::Graph`: A Graph.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graph.
# 
"""
function merge_all_multi_products!(g::Graph; verbose=0)
    verbose > 0 && println("merge nodes representing a multi product of a non-unique list of graphs.")
    # Post-order DFS
    for sub_g in g.subgraphs
        merge_all_multi_products!(sub_g)
        merge_multi_product!(sub_g)
    end
    merge_multi_product!(g)
    return g
end

"""
    function merge_all_multi_products!(graphs::Union{Tuple,AbstractVector{<:Graph}}; verbose=0)

    Merges all nodes representing a multi product of a non-unique list of subgraphs in-place in the given graphs. 

# Arguments:
- `graphs`: A collection of graphs to be processed.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- Optimized graphs.
# 
"""
function merge_all_multi_products!(graphs::Union{Tuple,AbstractVector{<:Graph}}; verbose=0)
    verbose > 0 && println("merge nodes representing a multi product of a non-unique list of graphs.")
    # Post-order DFS
    for g in graphs
        merge_all_multi_products!(g.subgraphs)
        merge_multi_product!(g)
    end
    return graphs
end

"""
    function unique_nodes!(graphs::AbstractVector{<:AbstractGraph})

    Identifies and retrieves unique nodes from a set of graphs.

# Arguments:
- `graphs`: A collection of graphs to be processed.

# Returns:
- A mapping dictionary from the id of each leaf to the unique leaf node.
"""
function unique_nodes!(graphs::AbstractVector{<:AbstractGraph}, mapping::Dict{Int,<:AbstractGraph}=Dict{Int,eltype(graphs)}())
    # function unique_nodes!(graphs::AbstractVector{<:AbstractGraph})
    ############### find the unique Leaves #####################
    # unique_graphs = []
    # mapping = Dict{Int,eltype(graphs)}()
    unique_graphs = collect(values(mapping))

    for g in graphs
        flag = true
        for e in unique_graphs
            if isequiv(e, g, :id, :name, :weight)
                mapping[id(g)] = e
                flag = false
                break
            end
        end
        if flag
            push!(unique_graphs, g)
            mapping[id(g)] = g
        end
    end
    return mapping
end

"""
    function remove_duplicated_leaves!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0, normalize=nothing, kwargs...)

    Removes duplicated leaf nodes in-place from a collection of graphs. It also provides optional normalization for these leaves.

# Arguments:
- `graphs`: A collection of graphs to be processed.
- `verbose`: Level of verbosity (default: 0).
- `normalize`: Optional function to normalize the graphs (default: nothing).
"""
# function remove_duplicated_leaves!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0, normalize=nothing, kwargs...)
#     verbose > 0 && println("remove duplicated leaves.")
#     leaves = Vector{eltype(graphs)}()
#     for g in graphs
#         append!(leaves, collect(Leaves(g)))
#     end
#     if isnothing(normalize) == false
#         @assert normalize isa Function "a function call is expected for normalize"
#         for leaf in leaves
#             normalize(id(leaf))
#         end
#     end
#     sort!(leaves, by=x -> id(x)) #sort the id of the leaves in an asscend order
#     unique!(x -> id(x), leaves) #filter out the leaves with the same id number

#     mapping = unique_nodes!(leaves)

#     for g in graphs
#         for n in PreOrderDFS(g)
#             for (si, sub_g) in enumerate(subgraphs(n))
#                 if isleaf(sub_g)
#                     set_subgraph!(n, mapping[id(sub_g)], si)
#                 end
#             end
#         end
#     end

#     return graphs
# end

"""
    remove_duplicated_leaves!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0, normalize=nothing)

使用哈希表原地去重叶子节点。复杂度从 O(N^2) 降低为 O(N)。
遵循 isequiv(a, b, :id, :name, :weight) 规则，即忽略 id, name, weight 进行比较。
"""
# function remove_duplicated_leaves!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0, normalize=nothing)
#     verbose > 0 && println("Optimizing leaves with Hash Map...")


#     # 2. 全局叶子缓存：LeafKey -> 唯一的 Leaf 节点对象
#     leaf_cache = Dict{LeafKey,eltype(graphs)}()

#     # 3. 遍历去重
#     visited = Set{Int}()

#     function _optimize_leaves_recursive!(g::AbstractGraph)
#         if g.id in visited
#             return
#         end

#         for (i, sub_g) in enumerate(subgraphs(g))
#             if isleaf(sub_g)
#                 key = LeafKey(sub_g.operator, sub_g.orders, sub_g.properties)

#                 if haskey(leaf_cache, key)
#                     canonical_leaf = leaf_cache[key]
#                     set_subgraph!(g, canonical_leaf, i)
#                 else
#                     # if !isnothing(normalize)
#                     #     normalize(sub_g.id)
#                     # end
#                     leaf_cache[key] = sub_g
#                 end
#             else
#                 _optimize_leaves_recursive!(sub_g)
#             end
#         end

#         push!(visited, g.id)
#     end

#     for g in graphs
#         _optimize_leaves_recursive!(g)
#     end

#     return graphs
# end

function remove_duplicated_leaves!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0, normalize=nothing)
    verbose > 0 && println("Optimizing leaves with Hash Map...")

    leaf_cache = Dict{LeafKey,eltype(graphs)}()
    visited = Set{Int}()

    function _process_node_children!(g::AbstractGraph)
        if g.id in visited
            return
        end
        push!(visited, g.id)

        for (i, sub_g) in enumerate(subgraphs(g))
            if isleaf(sub_g)
                key = LeafKey(sub_g.operator, sub_g.orders, sub_g.properties)
                if haskey(leaf_cache, key)
                    set_subgraph!(g, leaf_cache[key], i)
                else
                    # if !isnothing(normalize); normalize(sub_g.id); end
                    leaf_cache[key] = sub_g
                end
                # --------------------
            else
                _process_node_children!(sub_g)
            end
        end
    end

    if graphs isa AbstractVector
        for (i, g) in enumerate(graphs)
            if isleaf(g)
                key = LeafKey(g.operator, g.orders, g.properties)
                if haskey(leaf_cache, key)
                    graphs[i] = leaf_cache[key]
                else

                    leaf_cache[key] = g
                end
            else
                _process_node_children!(g)
            end
        end
    else
        for g in graphs
            if !isleaf(g)
                _process_node_children!(g)
            end
        end
    end

    return graphs
end

function remove_duplicated_nodes!(graphs::Union{Tuple,AbstractVector{<:AbstractGraph}}; verbose=0, kwargs...)
    verbose > 0 && println("remove duplicated nodes.")

    nodes_all = Vector{eltype(graphs)}()
    for g in graphs
        for node in PostOrderDFS(g)
            push!(nodes_all, node)
        end
    end

    sort!(nodes_all, by=x -> id(x)) #sort the id of the leaves in an asscend order
    unique!(x -> id(x), nodes_all) #filter out the leaves with the same id number

    mapping = unique_nodes!(nodes_all)

    for g in graphs
        for n in PreOrderDFS(g)
            for (si, sub_g) in enumerate(subgraphs(n))
                set_subgraph!(n, mapping[id(sub_g)], si)
            end
        end
    end

    return graphs
end

function remove_duplicated_nodes!(root::G; verbose=0) where {G<:AbstractGraph}
    verbose > 0 && println("remove duplicated nodes.")
    # A dictionary to keep track of unique nodes based on a key (like id, or a hash of properties)

    # remove_duplicated_leaves!([root])

    unique_nodes = Dict{Int,G}()
    # for l in Leaves(root)
    #     if !haskey(unique_nodes, id(l))
    #         unique_nodes[id(l)] = l
    #     end
    # end

    # Helper function to process a node
    function process_node(node)
        # Compute a key for the node (here, I'm using `id` for simplicity)
        node_key = id(node)

        # Check if a node with the same key already exists
        if haskey(unique_nodes, node_key)
            return unique_nodes[node_key]
        else
            # Check if the node is equivalent to any existing unique node
            for g in values(unique_nodes)
                if isequiv(node, g, :id, :name, :weight)
                    return g
                end
            end

            # Process child nodes if the node is unique
            for (i, child) in enumerate(subgraphs(node))
                unique_child = process_node(child)
                set_subgraph!(node, unique_child, i)
            end

            # Add the (now potentially updated) node to the unique_nodes dictionary
            unique_nodes[node_key] = node
            return node
        end
    end

    # Start processing from the root
    process_node(root)

    return root
end

"""
    function burn_from_targetleaves!(graphs::AbstractVector{G}, targetleaves_id::AbstractVector{Int}; verbose=0) where {G <: AbstractGraph}

    Removes all nodes connected to the target leaves in-place via "Prod" operators.

# Arguments:
- `graphs`: A vector of graphs.
- `targetleaves_id::AbstractVector{Int}`: Vector of target leafs' id.
- `verbose`: Level of verbosity (default: 0).

# Returns:
- The id of a constant graph with a zero factor if any graph in `graphs` was completely burnt; otherwise, `nothing`.
"""
function burn_from_targetleaves!(graphs::AbstractVector{G}, targetleaves_id::AbstractVector{Int}; verbose=0) where {G<:AbstractGraph}
    verbose > 0 && println("remove all nodes connected to the target leaves via Prod operators.")

    graphs_sum = linear_combination(graphs, one.(eachindex(graphs)))
    ftype = eltype(subgraph_factors(graphs[1]))

    for leaf in Leaves(graphs_sum)
        if !isdisjoint(id(leaf), targetleaves_id)
            set_name!(leaf, "BURNING")
        end
    end

    for node in PostOrderDFS(graphs_sum)
        if any(x -> name(x) == "BURNING", subgraphs(node))
            if operator(node) == Prod || operator(node) <: Power
                set_subgraphs!(node, G[])
                set_subgraph_factors!(node, ftype[])
                set_name!(node, "BURNING")
            else
                _subgraphs = G[]
                _subgraph_factors = ftype[]
                for (i, subg) in enumerate(subgraphs(node))
                    if name(subg) != "BURNING"
                        push!(_subgraphs, subg)
                        push!(_subgraph_factors, subgraph_factor(node, i))
                    end
                end
                set_subgraphs!(node, _subgraphs)
                set_subgraph_factors!(node, _subgraph_factors)
                if isempty(_subgraph_factors)
                    set_name!(node, "BURNING")
                end
            end
        end
    end

    # g_c0 = constant_graph(ftype(0))
    g_c1 = constant_graph(ftype(1))
    has_c0 = false
    for g in graphs
        if name(g) == "BURNING"
            has_c0 = true
            set_id!(g, id(g_c1))
            set_operator!(g, Unitary)
            # set_subgraphs!(g, subgraphs(g_c0))
            # set_subgraph_factors!(g, subgraph_factors(g_c0))
            set_weight!(g, 0.0)
        end
    end

    has_c0 ? (return id(g_c1)) : (return nothing)
end