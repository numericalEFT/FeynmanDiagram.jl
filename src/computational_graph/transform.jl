# this file is included in ComputationalGraphs.jl

# relabel constructor for QuantumOperator
QuantumOperator(qo::QuantumOperator, label::Int) = QuantumOperator(qo.operator(), label, qo.is_ghost)
# relabel constructor for OperatorProduct

"""
    function relabel!(g::Graph, map::Dict{Int,Int})

    This function maps the labels of the quantum operators in g and its subgraphs to other labels. 
For example, map = {1=>2, 3=>2} will find all quantum operators with labels 1 and 3, and then map them to 2.
The graph g is modified.

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
For example, map = {1=>2, 3=>2} will find all quantum operators with labels 1 and 3, and then map them to 2.
Return a new copy of modified g.

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
Return a copy of stantardized g.
"""
standardize_labels(g::Graph) = standardize_labels!(deepcopy(g))

"""
    function prune_trivial_unary(g::Graph)

Simplifies a graph g if it represents a trivial unary operation. Otherwise, returns the original graph.
"""
function prune_trivial_unary(g::Graph)
    # No-op; g is not a branch (depth-1, one-child tree)
    if isbranch(g) == false
        return g
    end
    # Prune trivial unary operations
    if unary_istrivial(g.operator) && isfactorless(g)
        return eldest(g)
    else
        return g
    end
end

"""
    function simplify_products(g::Graph)

Simplifies subgraph factors for a graph g by shifting them up to root level and merging the link.
"""
function simplify_products(g::Graph)
    # No-op for non-multiplicative node operations, branches/leaves, and non-chain graphs
    if g.operator != Prod || isleaf(g) || isbranch(g) || onechild(g) == false
        return g
    end
    # Shift multiplicative subfactors to root level
    gs = deepcopy(g)
    child_gs = eldest(g)
    while onechild(g)
        gs.subgraph_factors[1] *= child_gs.subgraph_factors[1]
        child_gs.subgraph_factors[1] = 1
        g = child_gs
    end
    # 
end

"""
    function inplace_prod(g::Graph{F,W}) where {F,W}

Converts a unary Prod chain to in-place form by propagating subgraph_factors up the chain.
"""
function inplace_prod(g::Graph{F,W}) where {F,W}
    # Find unary Prod chain subgraphs of g
    gt = deepcopy(g)
    for sg in gt

    end

    if g.operator == Prod && ischain(g)
        gs = g.subgraphs[1]
        return Graph(gs.vertices; external=gs.external, type=gs.type, topology=gs.topology, subgraphs=gs.subgraphs,
            factor=g.subgraph_factors[1] * g.factor * gs.factor, operator=gs.operator, ftype=F, wtype=W)
    else
        return g
    end
end

"""
    function merge_chain(g::Graph)

Simplifies a graph g with a unary operator chain as its root.
"""
function merge_chain(g::Graph)
    if g.depth == 0
        return g
    elseif g.depth == 1
        # A branch is mergeable iff the unary operation is trivial
        return prune_trivial_unary(g)
    end
    # First, simplify subgraph factors for multiplicative chains
    gs = g.operator == Prod ? simplify_subfactors(g) : deepcopy(g)
    # Then, merge operator chains of length > 1
    while onechild(gs)
        # Break case: last branch found
        if node.depth == 1
            # A branch is mergeable iff the unary operation is trivial
            return prune_trivial_unary(g)
        else
            # If the branch is not factorless, propagate the subgraph factor up the chain
        end

        isunary = length(gs) == 1
        isprunable = g.subgraph_factors[1] == 1 && g.factor == 1 && g.operator in [Prod, Sum]
        if isunary && isprunable
            return gs[1]
        else
            return g
        end
    end
    return error("Encountered an unexpected error.")
end

"""
    function merge_prefactors(g::Graph)
        Factorize the prefactors of a multiplicative graph g.
"""
function merge_prefactors(g0::Graph)
    if (g0.operator==Sum)
        added = falses(length(g0.subgraphs))
        subg_fac = (eltype(g0.subgraph_factors))[]
        subg = (eltype(g0.subgraphs))[]
        k = 0
        for i in eachindex(added)
            if added[i] 
                continue
            end
            push!(subg,g0.subgraphs[i])
            push!(subg_fac,g0.subgraph_factors[i])
            added[i] = true
            k += 1
            for j in i+1:length(g0.subgraphs)
                if(added[j] == false && isequiv(g0.subgraphs[i], g0.subgraphs[j], :id))
                    added[j] = true
                    subg_fac[k] += g0.subgraph_factors[j]
                end
            end
        end
        g = Graph(g0.vertices; external=g0.external, type=g0.type, topology=g0.topology,
        subgraphs=subg, subgraph_factors= subg_fac, operator= g0.operator())
        return g
    else
        return g0
    end
end

# function prune_unary(g::Graph{F, W}) where {F, W}
#     if (g.operator in [Prod, Sum] && length(g.subgraph) = 1 && g.subgraph_factors[1]≈ one(F) && g.factor ≈ one(F))
#         return g.subgraph[1]
#     else
#         return g
#     end
# end


# function inplace_prod(g1::Graph) 
#     if (length(g1.subgraphs)==1 && length(g1.subgraphs[1].subgraphs) ==1 && (g1.subgraphs[1].operator == Prod))
#         g1.subgraph_factors[1] *= g1.subgraphs[1].subgraph_factors[1] 
#         g1.subgraphs[1].subgraph_factors[1] = 1
#     end
#     return g1
# end

######

# """Converts a unary Prod node to in-place form by merging factors and subgraph_factors."""
# function inplace_prod(g::Graph{F,W}) where {F,W}
#     if g.operator == Prod && length(g.subgraphs) == 1
#         gs = g.subgraphs[1]
#         return Graph(gs.vertices; external=gs.external, type=gs.type, topology=gs.topology, subgraphs=gs.subgraphs,
#             factor=g.subgraph_factors[1] * g.factor * gs.factor, operator=gs.operator, ftype=F, wtype=W)
#     else
#         return g
#     end
# end


############LEGACY BELOW################

# function readDiag(io::IO)

#     return Diagram
# end

# function printDiag(diag::Union{Tuple,AbstractVector}, kwargs...)

#     f = open("./DiagFiles/Diag.txt", "w")
#     writedlm(f, [])
# end
