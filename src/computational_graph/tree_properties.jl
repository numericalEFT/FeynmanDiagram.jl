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

# Does the graph have any children?
haschildren(g::Graph) = isempty(g.subgraphs) == false

# Is the graph a leaf?
isleaf(g::Graph) = isempty(g.subgraphs)

# Does the graph have only one child?
onechild(g::Graph) = length(children(g)) == 1

# Is the graph a branch (depth-1 and one-child)?
isbranch(g::Graph) = onechild(g) && isleaf(eldest(g))

# Get the first child of a graph
function eldest(g::Graph)
    @assert haschildren(g) "Graph has no children!"
    return children(g)[1]
end

# Is the graph a chain?
function ischain(g::Graph)
    isleaf(g) && return true
    while onechild(g)
        isleaf(g) && return true
        g = eldest(g)
    end
    return false
end

# Is the graph factorless?
function isfactorless(g)
    if isleaf(g)
        return isapprox_one(g.factor)
    else
        return all(isapprox_one.([g.factor; g.subgraph_factors]))
    end
end
