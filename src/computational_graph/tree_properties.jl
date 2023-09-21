##################### AbstractTrees interface for AbstracGraphs ########################### 

## Things that make printing prettier
AbstractTrees.printnode(io::IO, g::AbstractGraph) = print(io, "\u001b[32m$(g.id)\u001b[0m : $g")

## Guarantee type-stable tree iteration for Graphs, StableGraphs, and FeynmanGraphs
AbstractTrees.NodeType(::Graph) = HasNodeType()
AbstractTrees.NodeType(::StableGraph) = HasNodeType()
AbstractTrees.NodeType(::FeynmanGraph) = HasNodeType()
AbstractTrees.nodetype(::Graph{F,W}) where {F,W} = Graph{F,W}
AbstractTrees.nodetype(::StableGraph{F,W,NT}) where {F,W,NT} = StableGraph{F,W,NT}
AbstractTrees.nodetype(::FeynmanGraph{F,W}) where {F,W} = FeynmanGraph{F,W}

## Optional enhancements
# These next two definitions allow inference of the item type in iteration.
# (They are not sufficient to solve all internal inference issues, however.)
Base.IteratorEltype(::Type{<:TreeIterator{Graph{F,W}}}) where {F,W} = Base.HasEltype()
Base.eltype(::Type{<:TreeIterator{Graph{F,W}}}) where {F,W} = Graph{F,W}
Base.IteratorEltype(::Type{<:TreeIterator{StableGraph{F,W,NT}}}) where {F,W,NT} = Base.HasEltype()
Base.eltype(::Type{<:TreeIterator{StableGraph{F,W,NT}}}) where {F,W,NT} = StableGraph{F,W,NT}
Base.IteratorEltype(::Type{<:TreeIterator{FeynmanGraph{F,W}}}) where {F,W} = Base.HasEltype()
Base.eltype(::Type{<:TreeIterator{FeynmanGraph{F,W}}}) where {F,W} = FeynmanGraph{F,W}

function AbstractTrees.children(g::AbstractGraph)
    return g.subgraphs
end

##################### Tree properties ########################### 

"""
    function haschildren(g::AbstractGraph)

    Returns whether the graph has any children (subgraphs).

# Arguments:
- `g::AbstractGraph`: graph to be analyzed
"""
haschildren(g::AbstractGraph) = isempty(g.subgraphs) == false

"""
    function onechild(g::AbstractGraph)

    Returns whether the graph g has only one child (subgraph).

# Arguments:
- `g::AbstractGraph`: graph to be analyzed
"""
onechild(g::AbstractGraph) = length(children(g)) == 1

"""
    function isleaf(g::AbstractGraph)

    Returns whether the graph g is a leaf (terminating tree node).

# Arguments:
- `g::AbstractGraph`: graph to be analyzed
"""
isleaf(g::AbstractGraph) = isempty(g.subgraphs)

"""
    function isbranch(g::AbstractGraph)

    Returns whether the graph g is a branch-type (depth-1 and one-child) graph.

# Arguments:
- `g::AbstractGraph`: graph to be analyzed
"""
isbranch(g::AbstractGraph) = onechild(g) && isleaf(eldest(g))

"""
    function ischain(g::AbstractGraph)

    Returns whether the graph g is a chain-type graph.

# Arguments:
- `g::AbstractGraph`: graph to be analyzed
"""
function ischain(g::AbstractGraph)
    isleaf(g) && return true
    while onechild(g)
        g = eldest(g)
        isleaf(g) && return true
    end
    return false
end

"""
    function isfactorless(g)

    Returns whether the graph g is factorless, i.e., has unity factor and, if applicable,
    subgraph factor(s). Note that this function does not recurse through subgraphs of g, so 
    that one may have, e.g., `isfactorless(g) == true` but `isfactorless(eldest(g)) == false`.

# Arguments:
- `g::AbstractGraph`: graph to be analyzed
"""
function isfactorless(g::AbstractGraph)
    if isleaf(g)
        return isapprox_one(g.factor)
    else
        return all(isapprox_one.([g.factor; g.subgraph_factors]))
    end
end

"""
    function eldest(g::AbstractGraph)

    Returns the first child (subgraph) of a graph g.

# Arguments:
- `g::AbstractGraph`: graph for which to find the first child
"""
function eldest(g::AbstractGraph)
    @assert haschildren(g) "Graph has no children!"
    return children(g)[1]
end
