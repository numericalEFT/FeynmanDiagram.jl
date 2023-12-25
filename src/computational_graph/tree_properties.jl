##################### AbstractTrees interface for AbstractGraphs ########################### 

## Things that make printing prettier
AbstractTrees.printnode(io::IO, g::AbstractGraph) = print(io, "\u001b[32m$(id(g))\u001b[0m : $g")

## Guarantee type-stable tree iteration for Graphs and FeynmanGraphs
AbstractTrees.NodeType(::Graph) = HasNodeType()
AbstractTrees.NodeType(::FeynmanGraph) = HasNodeType()
AbstractTrees.nodetype(::Graph{F,W}) where {F,W} = Graph{F,W}
AbstractTrees.nodetype(::FeynmanGraph{F,W}) where {F,W} = FeynmanGraph{F,W}

## Optional enhancements
# These next two definitions allow inference of the item type in iteration.
# (They are not sufficient to solve all internal inference issues, however.)
Base.IteratorEltype(::Type{<:TreeIterator{Graph{F,W}}}) where {F,W} = Base.HasEltype()
Base.eltype(::Type{<:TreeIterator{Graph{F,W}}}) where {F,W} = Graph{F,W}
Base.IteratorEltype(::Type{<:TreeIterator{FeynmanGraph{F,W}}}) where {F,W} = Base.HasEltype()
Base.eltype(::Type{<:TreeIterator{FeynmanGraph{F,W}}}) where {F,W} = FeynmanGraph{F,W}

function AbstractTrees.children(g::AbstractGraph)
    return subgraphs(g)
end

##################### Tree properties ########################### 

"""
    function haschildren(g::AbstractGraph)

    Returns whether the graph has any children (subgraphs).

# Arguments:
- `g::AbstractGraph`: graph to be analyzed
"""
haschildren(g::AbstractGraph) = isempty(subgraphs(g)) == false

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
isleaf(g::AbstractGraph) = isempty(subgraphs(g))

"""
    function isbranch(g::AbstractGraph)

    Returns whether the graph g is a branch-type (depth-1 and one-child) graph.

# Arguments:
- `g::AbstractGraph`: graph to be analyzed
"""
isbranch(g::AbstractGraph) = onechild(g) && isleaf(eldest(g))

"""
    function ischain(g::AbstractGraph)

    Returns whether the graph g is a chain-type graph (i.e., a unary string).

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
        return isapprox_one(factor(g))
    else
        return all(isapprox_one.([factor(g); subgraph_factors(g)]))
    end
end

"""
    function has_zero_subfactors(g)

    Returns whether the graph g has only zero-valued subgraph factor(s). 
    Note that this function does not recurse through subgraphs of g, so that one may have, e.g.,
    `isfactorless(g) == true` but `isfactorless(eldest(g)) == false`.
    By convention, returns `false` if g is a leaf.

# Arguments:
- `g::AbstractGraph`: graph to be analyzed
"""
function has_zero_subfactors(g::AbstractGraph)
    if isleaf(g)
        return false  # convention: subgraph_factors = [] ⟹ subfactorless = false
    else
        return iszero(subgraph_factors(g))
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
    return subgraph(g)
end

"""
    function count_operation(g::Graph)

    Returns the total number of  additions and multiplications in the graph.

# Arguments:
- `g::Graph`: graph for which to find the total number of  operations.
"""
function count_operation(g::G) where {G<:AbstractGraph}
    totalsum = 0
    totalprod = 0
    # totalpower = 0
    for node in PreOrderDFS(g)
        if length(node.subgraphs) > 0
            if node.operator == Prod
                totalprod += length(node.subgraphs) - 1
            elseif node.operator == Sum
                totalsum += length(node.subgraphs) - 1
                # elseif node.operator <: Power
                #     totalpower += 1
            end
        end
    end
    return [totalsum, totalprod]
end

function count_operation(g::Vector{G}) where {G<:AbstractGraph}
    visited = Set{Int}()
    totalsum = 0
    totalprod = 0
    for graph in g
        for node in PreOrderDFS(graph)
            if !(node.id in visited)
                push!(visited, node.id)
                if length(node.subgraphs) > 0
                    if node.operator == Prod
                        totalprod += length(node.subgraphs) - 1
                    elseif node.operator == Sum
                        totalsum += length(node.subgraphs) - 1
                    end
                end
            end
        end
    end
    return [totalsum, totalprod]
end


function count_operation(g::Dict{Vector{Int},G}) where {G<:AbstractGraph}
    visited = Set{Int}()
    totalsum = 0
    totalprod = 0
    for (order, graph) in g
        for node in PreOrderDFS(graph)
            if !(node.id in visited)
                push!(visited, node.id)
                if length(node.subgraphs) > 0
                    if node.operator == Prod
                        totalprod += length(node.subgraphs) - 1
                    elseif node.operator == Sum
                        totalsum += length(node.subgraphs) - 1
                    end
                end
            end
        end
    end
    return [totalsum, totalprod]
end

function count_operation(g::Number)
    return [0, 0]
end


function count_operation(nothing)
    return [0, 0]
end

function count_expanded_operation(g::G) where {G<:AbstractGraph}
    totalsum = 0
    totalprod = 0

    len_subg = length(subgraphs(g))
    subgraphs_sum = zeros(Int, len_subg)
    subgraphs_prod = zeros(Int, len_subg)
    for (i, subg) in enumerate(subgraphs(g))
        subgraphs_sum[i], subgraphs_prod[i] = count_expanded_operation(subg)
    end

    if isleaf(g)
        return [0, 0]
    else
        if operator(g) == Sum
            totalsum = sum(subgraphs_sum) + len_subg - 1
            totalprod = sum(subgraphs_prod)
        elseif operator(g) == Prod
            totalsum = prod(subgraphs_sum .+ 1) - 1
            innerprod = 0
            for i in 1:len_subg
                innerprod += subgraphs_prod[i] * prod([subgraphs_sum[j] + 1 for j in 1:len_subg if j != i])
            end
            totalprod = innerprod + (totalsum + 1) * (len_subg - 1)
        end
    end
    return [totalsum, totalprod]
end