# this file is included in ComputationalGraphs.jl

# relabel constructor for QuantumOperator
QuantumOperator(qo::QuantumOperator, label::Int) = QuantumOperator(qo.operator(), label)
# relabel constructor for OperatorProduct

"""
    function relabel!(g::Graph, map::Dict{Int,Int})

    This function maps the labels of the quantum operators in g and its subgraphs to other labels. 
    For example, map = {1=>2, 3=>2} will find all quantum operators with labels 1 and 3, and then map them to 2. The graph g is modified.

# Arguments:
- `g::Graph`: graph to be modified
- `map`: mapping from old labels to the new ones
"""
function relabel!(g::Graph, map::Dict{Int,Int})

    for i in eachindex(g.vertices)
        op = g.vertices[i]
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
    function relabel(g::Graph, map::Dict{Int,Int})

    This function maps the labels of the quantum operators in g and its subgraphs to other labels. 
    For example, map = {1=>2, 3=>2} will find all quantum operators with labels 1 and 3, and then map them to 2. Returns a modified copy of g.

# Arguments:
- `g::Graph`: graph to be modified
- `map`: mapping from old labels to the new ones
"""
relabel(g::Graph, map::Dict{Int,Int}) = relabel!(deepcopy(g), map)

"""
    function collect_labels(g::Graph)

Return sorted unique labels of graph g.
"""
function collect_labels(g::Graph)
    labels = Vector{Int}([])
    for i in eachindex(g.vertices)
        op = g.vertices[i]
        for j in eachindex(op.operators)
            qo = op.operators[j]
            if !(qo.label in labels)
                push!(labels, qo.label)
            end
        end
    end

    uniqlables = sort(unique(labels))
end

"""
    function standardize_labels!(g::Graph)

    This function first finds all labels involved in g and its subgraphs (for example, 1, 4, 5, 7, ...), 
    then relabel them in the order 1, 2, 3, 4, ....
    The graph g is modified.
"""
function standardize_labels!(g::Graph)
    #TBD
    uniqlabels = collect_labels(g)
    map = Dict{Int,Int}()
    for i in eachindex(uniqlabels)
        push!(map, uniqlabels[i] => i)
    end
    return relabel!(g, map)
end

"""
    function standardize_labels!(g::Graph)

    This function first finds all labels involved in g and its subgraphs (for example, 1, 4, 5, 7, ...), 
    then relabel them in the order 1, 2, 3, 4, ....
    Returns a standardized copy of g.
"""
standardize_labels(g::Graph) = standardize_labels!(deepcopy(g))

"""
    function replace_subgraph!(g::Graph, w::Graph, m::graph)

    In place function that replaces the children graph w in graph g with a new graph m.
    Graphs w and m should have the same internal and external vertices, and topology
"""
function replace_subgraph!(g::Graph, w::Graph, m::Graph)
    @assert !isleaf(g) "Target parent graph can not be a leaf"
    @assert w.vertices == m.vertices "Old and new subgraph should have the same vertices"
    @assert w.external == m.external "Old and new subgraph should have the same external vertices"
    print("isleaf $(isleaf(g))\n")
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
    function replace_subgraph(g::Graph, w::Graph, m::graph)

    Generate a copy of graph g, with the children graph w replaced by a new graph m.
    Graph w and m should have the same internal and external vertices, and topology
"""
function replace_subgraph(g::Graph, w::Graph, m::Graph)
    @assert w.vertices == m.vertices "Old and new subgraph should have the same vertices"
    @assert w.external == m.external "Old and new subgraph should have the same external vertices"
    g0 = deepcopy(g)
    for node in PreOrderDFS(g0)
        for (i, child) in enumerate(children(node))
            if isequiv(child, w, :id)
                node.subgraphs[i] = m
                break
            end
        end
    end
    return g0
end

"""
    function prune_trivial_unary(g::Graph)

    Returns a simplified copy of g if it represents a trivial unary operation.
    Otherwise, returns the original graph.
"""
function prune_trivial_unary(g::Graph)
    # No-op; g is not a branch (depth-1, one-child tree)
    if isbranch(g) == false
        return g
    end
    # Prune trivial unary operation
    if unary_istrivial(g.operator) && isfactorless(g)
        return eldest(g)
    else
        return g
    end
end

"""
    function inplace_prod(g::Graph)

    Converts a unary Prod link to in-place form by propagating subgraph_factors up a level.
"""
function inplace_prod(g::Graph)
    if onechild(g) == false
        return g
    end
    child = eldest(g)
    if onechild(child) && child.operator == Prod
        # Merge subgraph factors at parent tree level
        g.subgraph_factors[1] *= child.subgraph_factors[1]
        child.subgraph_factors[1] = 1
    end
    return g
end

"""
    function merge_prefactors(g::Graph)
   
    Factorize the prefactors of a multiplicative graph g.
"""
function merge_prefactors(g0::Graph{F,W}) where {F,W}
    if g0.operator == Sum
        added = falses(length(g0.subgraphs))
        subg_fac = eltype(g0.subgraph_factors)[]
        subg = eltype(g0.subgraphs)[]
        k = 0
        for i in eachindex(added)
            if added[i]
                continue
            end
            push!(subg, g0.subgraphs[i])
            push!(subg_fac, g0.subgraph_factors[i])
            added[i] = true
            k += 1
            for j in (i+1):length(g0.subgraphs)
                if added[j] == false && isequiv(g0.subgraphs[i], g0.subgraphs[j], :id)
                    added[j] = true
                    subg_fac[k] += g0.subgraph_factors[j]
                end
            end
        end
        g = Graph(subg; topology=g0.topology, vertices=g0.vertices, external=g0.external, hasLeg=g0.hasLeg,
            subgraph_factors=subg_fac, type=g0.type(), operator=g0.operator())
        return g
    else
        return g0
    end
end

############LEGACY BELOW################

# function readDiag(io::IO)

#     return Diagram
# end

# function printDiag(diag::Union{Tuple,AbstractVector}, kwargs...)

#     f = open("./DiagFiles/Diag.txt", "w")
#     writedlm(f, [])
# end
