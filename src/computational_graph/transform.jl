# this file is included in ComputationalGraphs.jl

"""
    function relabel!(g::FeynmanGraph, map::Dict{Int,Int})

    Relabels the quantum operators in g and its subgraphs according to `map`.
    For example, `map = {1=>2, 3=>2}`` will find all quantum operators with labels 1 and 3, and then map them to 2.

# Arguments:
- `g::FeynmanGraph`: graph to be modified
- `map`: mapping from old labels to the new ones
"""
function relabel!(g::FeynmanGraph, map::Dict{Int,Int})

    for i in eachindex(vertices(g))
        op = vertices(g)[i]
        for j in eachindex(op.operators)
            qo = op.operators[j]
            if haskey(map, qo.label)
                op.operators[j] = QuantumOperator(qo, map[qo.label])
            end
        end
    end

    for i in eachindex(g.subgraphs)
        relabel!(g.subgraphs[i], map)
    end

    return g
end

"""
    function relabel(g::FeynmanGraph, map::Dict{Int,Int})

    Returns a copy of g with quantum operators in g and its subgraphs relabeled according to `map`.
    For example, `map = {1=>2, 3=>2}` will find all quantum operators with labels 1 and 3, and then map them to 2.

# Arguments:
- `g::FeynmanGraph`: graph to be modified
- `map`: mapping from old labels to the new ones
"""
relabel(g::FeynmanGraph, map::Dict{Int,Int}) = relabel!(deepcopy(g), map)

"""
    function collect_labels(g::FeynmanGraph)

    Returns the list of sorted unique labels in graph g.

# Arguments:
- `g::FeynmanGraph`: graph to find labels for
"""
function collect_labels(g::FeynmanGraph)
    labels = Vector{Int}([])
    for i in eachindex(vertices(g))
        op = vertices(g)[i]
        for j in eachindex(op.operators)
            qo = op.operators[j]
            if !(qo.label in labels)
                push!(labels, qo.label)
            end
        end
    end

    uniqlables = sort(unique(labels))
    return uniqlables
end

"""
    function standardize_labels!(g::FeynmanGraph)

    Finds all labels involved in g and its subgraphs and 
    modifies g by relabeling in standardized order, e.g.,
    (1, 4, 5, 7, ...) ↦ (1, 2, 3, 4, ....)

# Arguments:
- `g::FeynmanGraph`: graph to be relabeled
"""
function standardize_labels!(g::FeynmanGraph)
    #TBD
    uniqlabels = collect_labels(g)
    map = Dict{Int,Int}()
    for i in eachindex(uniqlabels)
        push!(map, uniqlabels[i] => i)
    end
    return relabel!(g, map)
end

"""
    function standardize_labels!(g::FeynmanGraph)

    Finds all labels involved in g and its subgraphs and returns 
    a copy of g relabeled in a standardized order, e.g.,
    (1, 4, 5, 7, ...) ↦ (1, 2, 3, 4, ....)

# Arguments:
- `g::FeynmanGraph`: graph to be relabeled
"""
standardize_labels(g::FeynmanGraph) = standardize_labels!(deepcopy(g))

"""
    function replace_subgraph!(g::AbstractGraph, w::AbstractGraph, m::AbstractGraph)

    Modifies g by replacing the subgraph w with a new graph m.
    For Feynman diagrams, subgraphs w and m should have the same diagram type, orders, and external indices.

# Arguments:
- `g::AbstractGraph`: graph to be modified
- `w::AbstractGraph`: subgraph to replace
- `m::AbstractGraph`: new subgraph
"""
function replace_subgraph!(g::AbstractGraph, w::AbstractGraph, m::AbstractGraph)
    if g isa FeynmanGraph
        @assert w isa FeynmanGraph && m isa FeynmanGraph "Feynman diagrams should be replaced with Feynman diagrams"
        @assert isleaf(g) == false "Target parent graph cannot be a leaf"
        @assert diagram_type(w) == diagram_type(m) "Old and new subgraph should have the same diagram type"
        @assert orders(w) == orders(m) "Old and new subgraph should have the same orders"
        @assert external_indices(w) == external_indices(m) "Old and new subgraph should have the same external indices"
    end
    for node in PreOrderDFS(g)
        for (i, child) in enumerate(children(node))
            if isequiv(child, w, :id)
                node.subgraphs[i] = m
                return
            end
        end
    end
end

"""
    function replace_subgraph(g::AbstractGraph, w::AbstractGraph, m::AbstractGraph)

    Creates a modified copy of g by replacing the subgraph w with a new graph m.
    For Feynman diagrams, subgraphs w and m should have the same diagram type, orders, and external indices.

# Arguments:
- `g::AbstractGraph`: graph to be modified
- `w::AbstractGraph`: subgraph to replace
- `m::AbstractGraph`: new subgraph
"""
function replace_subgraph(g::AbstractGraph, w::AbstractGraph, m::AbstractGraph)
    if g isa FeynmanGraph
        @assert w isa FeynmanGraph && m isa FeynmanGraph "Feynman diagrams should be replaced with Feynman diagrams"
        @assert isleaf(g) == false "Target parent graph cannot be a leaf"
        @assert diagram_type(w) == diagram_type(m) "Old and new subgraph should have the same diagram type"
        @assert orders(w) == orders(m) "Old and new subgraph should have the same orders"
        @assert external_indices(w) == external_indices(m) "Old and new subgraph should have the same external indices"
    end
    g_new = deepcopy(g)
    for node in PreOrderDFS(g_new)
        for (i, child) in enumerate(children(node))
            if isequiv(child, w, :id)
                node.subgraphs[i] = m
                break
            end
        end
    end
    return g_new
end

"""
    function prune_trivial_unary(g::AbstractGraph)

    Returns a simplified copy of g if it represents a trivial unary chain.
    Otherwise, returns the original graph. For example, +(+(+g)) ↦ g.
    Does nothing unless g has the following structure: Ⓧ --- ⋯ --- Ⓧ ⋯ (!),
    where the stop-case (!) represents a leaf, an operator 𝓞' != Ⓧ, or a non-unary Ⓧ node.

# Arguments:
- `g::AbstractGraph`: graph to be modified
"""
function prune_trivial_unary(g::AbstractGraph)
    while unary_istrivial(g.operator) && onechild(g) && isfactorless(g)
        g = eldest(g)
    end
    return g
end

"""
    function merge_prodchain_subfactors!(g::AbstractGraph)

    Simplifies the subgraph factors of a graph g representing a unary Prod
    chain by merging them at root level, e.g., 2*(3*(5*g)) ↦ 30*(*(*g)). 
    Does nothing unless g has the following structure: 𝓞 --- Ⓧ --- ⋯ --- Ⓧ ⋯ (!),
    where the stop-case (!) represents a leaf, an operator 𝓞' != Ⓧ, or a non-unary Ⓧ node.

# Arguments:
- `g::AbstractGraph`: graph to be modified
"""
function merge_prodchain_subfactors!(g::AbstractGraph)
    if isleaf(g) || onechild(g) == false
        return g
    end
    child = eldest(g)
    children_factor = 1
    while onechild(child)
        # Break case: end of Prod chain, found 𝓞' != Ⓧ
        child.operator != Prod && break
        # Move this subfactor to running total
        children_factor *= child.subgraph_factors[1]
        child.subgraph_factors[1] = 1
        # Descend one level
        child = eldest(child)
    end
    # Update g subfactor with total factors from children
    g.subgraph_factors[1] *= children_factor
    return g
end

"""
    function merge_prodchain_subfactors(g::AbstractGraph)

    Returns a copy of a graph g representing a unary Prod chain with subgraph factors
    simplified by merging them at the root level, e.g., 2*(3*(5*g)) ↦ 30*(*(*g)). 
    Does nothing unless g has the following structure: 𝓞 --- Ⓧ --- ⋯ --- Ⓧ ⋯ (!),
    where the stop-case (!) represents a leaf, an operator 𝓞' != Ⓧ, or a non-unary Ⓧ node.

# Arguments:
- `g::AbstractGraph`: graph to be modified
"""
merge_prodchain_subfactors(g::AbstractGraph) = merge_prodchain_subfactors!(deepcopy(g))

"""
    function inplace_prod!(g::AbstractGraph)

    Converts a graph g representing a unary Prod chain to in-place form by merging its subgraph factors at
    root level and pruning the resultant unary product operation(s), e.g., 2*(3*(5*g)) ↦ 30*(*(*g)) ↦ 30*g.
    Does nothing unless g has the following structure: 𝓞 --- Ⓧ --- ⋯ --- Ⓧ ⋯ (!),
    where the stop-case (!) represents a leaf, an operator 𝓞' != Ⓧ, or a non-unary Ⓧ node.

# Arguments:
- `g::AbstractGraph`: graph to be modified
"""
function inplace_prod!(g::AbstractGraph)
    if isleaf(g) || onechild(g) == false
        return g
    end
    # First shift subfactors to root level, then prune left-over trivial unary operations.
    merge_prodchain_subfactors!(g)
    g.subgraphs[1] = prune_trivial_unary(eldest(g))
    return g
end

"""
    function inplace_prod(g::AbstractGraph)

    Returns a copy of a graph g representing a unary Prod chain converted to in-place form by merging its subgraph 
    factors at root level and pruning the resultant unary product operation(s), e.g., 2*(3*(5*g)) ↦ 30*(*(*g)) ↦ 30*g.
    Does nothing unless g has the following structure: 𝓞 --- Ⓧ --- ⋯ --- Ⓧ ⋯ (!),
    where the stop-case (!) represents a leaf, an operator 𝓞' != Ⓧ, or a non-unary Ⓧ node.

# Arguments:
- `g::AbstractGraph`: graph to be modified
"""
inplace_prod(g::AbstractGraph) = inplace_prod!(deepcopy(g))

"""
    function merge_prefactors(g::Graph)
   
    Returns a copy of graph g with multiplicative prefactors factorized,
    e.g., 3*g1 + 5*g2 + 7*g1 + 9*g2 ↦ 10*g1 + 14*g2. Does nothing if the
    graph g does not represent a Sum operation.

# Arguments:
- `g::Graph`: graph to be modified
"""
function merge_prefactors(g::Graph{F,W}) where {F,W}
    if g.operator == Sum
        added = falses(length(g.subgraphs))
        subg_fac = eltype(g.subgraph_factors)[]
        subg = eltype(g.subgraphs)[]
        k = 0
        for i in eachindex(added)
            added[i] && continue
            push!(subg, g.subgraphs[i])
            push!(subg_fac, g.subgraph_factors[i])
            added[i] = true
            k += 1
            for j in (i+1):length(g.subgraphs)
                if added[j] == false && isequiv(g.subgraphs[i], g.subgraphs[j], :id)
                    added[j] = true
                    subg_fac[k] += g.subgraph_factors[j]
                end
            end
        end
        nodedata = union(getproperty.(subg, :nodedata))
        g_merged = Graph(subg; subgraph_factors=subg_fac,
            nodedata=nodedata, operator=Sum(), ftype=F, wtype=W)
        return g_merged
    else
        return g
    end
end

"""
    function merge_prefactors(g::FeynmanGraph)
   
    Returns a copy of Feynman graph g with multiplicative prefactors factorized,
    e.g., 3*g1 + 5*g2 + 7*g1 + 9*g2 ↦ 10*g1 + 14*g2. Does nothing if the
    graph g does not represent a Sum operation.

# Arguments:
- `g::FeynmanGraph`: graph to be modified
"""
function merge_prefactors(g::FeynmanGraph{F,W}) where {F,W}
    if g.operator == Sum
        added = falses(length(g.subgraphs))
        subg_fac = eltype(g.subgraph_factors)[]
        subg = eltype(g.subgraphs)[]
        k = 0
        for i in eachindex(added)
            added[i] && continue
            push!(subg, g.subgraphs[i])
            push!(subg_fac, g.subgraph_factors[i])
            added[i] = true
            k += 1
            for j in (i+1):length(g.subgraphs)
                if added[j] == false && isequiv(g.subgraphs[i], g.subgraphs[j], :id)
                    added[j] = true
                    subg_fac[k] += g.subgraph_factors[j]
                end
            end
        end
        g_merged = FeynmanGraph(subg, g.properties; subgraph_factors=subg_fac, operator=Sum(), ftype=F, wtype=W)
        return g_merged
    else
        return g
    end
end
