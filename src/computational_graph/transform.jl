function readDiag(io::IO)

    return Diagram
end

function printDiag(diag::Union{Tuple,AbstractVector}, kwargs...)

    f = open("./DiagFiles/Diag.txt", "w")
    writedlm(f, [])
end

"""
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

"""Converts a unary Prod chain to in-place form by propagating subgraph_factors up the chain."""
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

"""Simplifies a graph g with a unary operator chain as its root."""
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

# function merge_prefactors(g0::Graph{F,W}) where {F,W}
#     if (g1.operator==Sum && length(g1.subgraphs)==2 && isequiv(g1.subgraphs[1], g1.subgraphs[2], :factor, :id, :subgraph_factors))
#         g1 = g0.subgraphs[1]
#         g2 = g0.subgraphs[2]
#         g_subg = Graph(g1.vertices; external=g1.external, type=g1.type, topology=g1.topology,
#         subgraphs=g1.subgraphs, operator=g1.operator(), ftype=F, wtype=W)
#         g = Graph(g1.vertices; external=g1.external, type=g1.type, topology=g1.topology,
#         subgraphs=[g_subg,], operator=Prod(), ftype=F, wtype=W)
#         g.subgraph_factors[1] = (g1.subgraph_factors[1]*g1.factor+g1.subgraph_factors[2]*g1.subgraphs[2].factor) * g0.factor
#         return g
#     else
#         return g1
#     end
# end
