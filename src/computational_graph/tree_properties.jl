##################### AbstractTrees interface for Graphs ########################### 

## Things that make printing prettier
AbstractTrees.printnode(io::IO, g::Graph) = print(io, "\u001b[32m$(g.id)\u001b[0m : $g")

## Guarantee type-stable tree iteration for Graphs
AbstractTrees.NodeType(::Graph) = HasNodeType()
AbstractTrees.nodetype(::Graph{F,W}) where {F,W} = Graph{F,W}

## Optional enhancements
# These next two definitions allow inference of the item type in iteration.
# (They are not sufficient to solve all internal inference issues, however.)
Base.IteratorEltype(::Type{<:TreeIterator{Graph{F,W}}}) where {F,W} = Base.HasEltype()
Base.eltype(::Type{<:TreeIterator{Graph{F,W}}}) where {F,W} = Graph{F,W}

function AbstractTrees.children(g::Graph)
    return g.subgraphs
end

##################### Tree properties ########################### 

"""
    function haschildren(g::Graph)

    Returns whether the graph has any children (subgraphs).

# Arguments:
- `g::Graph`: graph to be analyzed
"""
haschildren(g::Graph) = isempty(g.subgraphs) == false

"""
    function onechild(g::Graph)

    Returns whether the graph g has only one child (subgraph).

# Arguments:
- `g::Graph`: graph to be analyzed
"""
onechild(g::Graph) = length(children(g)) == 1

"""
    function isleaf(g::Graph)

    Returns whether the graph g is a leaf (terminating tree node).

# Arguments:
- `g::Graph`: graph to be analyzed
"""
isleaf(g::Graph) = isempty(g.subgraphs)

"""
    function isbranch(g::Graph)

    Returns whether the graph g is a branch-type (depth-1 and one-child) graph.

# Arguments:
- `g::Graph`: graph to be analyzed
"""
isbranch(g::Graph) = onechild(g) && isleaf(eldest(g))

"""
    function ischain(g::Graph)

    Returns whether the graph g is a chain-type graph.

# Arguments:
- `g::Graph`: graph to be analyzed
"""
function ischain(g::Graph)
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
- `g::Graph`: graph to be analyzed
"""
function isfactorless(g)
    if isleaf(g)
        return isapprox_one(g.factor)
    else
        return all(isapprox_one.([g.factor; g.subgraph_factors]))
    end
end

"""
    function eldest(g::Graph)

    Returns the first child (subgraph) of a graph g.

# Arguments:
- `g::Graph`: graph for which to find the first child
"""
function eldest(g::Graph)
    @assert haschildren(g) "Graph has no children!"
    return children(g)[1]
end

"""
    function totaloperation(g::Graph)

    Returns the total number of  additions and multiplications in the graph.

# Arguments:
- `g::Graph`: graph for which to find the total number of  operations.
"""
function totaloperation(g::Graph)
    totalsum = 0
    totalprod = 0
    for node in PreOrderDFS(g)
        if length(node.subgraphs) > 0 
            if node.operator == Prod
                totalprod += length(node.subgraphs) - 1
            elseif node.operator == Sum
                totalsum += length(node.subgraphs) - 1 
            end
        end
    end
    return [totalsum, totalprod]
end
