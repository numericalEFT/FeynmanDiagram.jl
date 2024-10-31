##################### AbstractTrees interface for AbstractGraphs ########################### 

## Things that make printing prettier
AbstractTrees.printnode(io::IO, g::AbstractGraph) = print(io, _stringrep(g; color=true))

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
    function has_zero_subfactors(g::AbstractGraph, operator_type::Type{<:AbstractOperator})

    Determines whether the graph `g` has only zero-valued subgraph factors based on the specified operator type.
    This function does not recurse through the subgraphs of `g`, so it only checks the immediate subgraph factors.
    If `g` is a leaf (i.e., has no subgraphs), the function returns `false` by convention.

    The behavior of the function depends on the operator type:
    - `Sum`: Checks if all subgraph factors are zero.
    - `Prod`: Checks if any subgraph factor is zero.
    - `Power{N}`: Checks if the first subgraph factor is zero.
    - Other `AbstractOperator`: Defaults to return `false`.
# Arguments:
- `g::AbstractGraph`: graph to be analyzed
- `operator`: the operator used in graph `g`
"""
function has_zero_subfactors(g::AbstractGraph, ::Type{Sum})
    @assert g.operator == Sum "Operator must be Sum"
    return iszero(subgraph_factors(g))
end

function has_zero_subfactors(g::AbstractGraph, ::Type{Prod})
    @assert g.operator == Prod "Operator must be Prod"
    return 0 in subgraph_factors(g)
end

function has_zero_subfactors(g::AbstractGraph, ::Type{Power{N}}) where {N}
    @assert g.operator <: Power "Operator must be a Power"
    return iszero(subgraph_factors(g)[1])
end

function has_zero_subfactors(g::AbstractGraph, ::Type{<:AbstractOperator})
    @info "has_zero_subfactors: Operator type $operator is not specifically defined. Defaults to return false."
    return false
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
    function count_leaves(g::G) where {G<:AbstractGraph}

    Returns the total number of leaves with unique id in the graph.

# Arguments:
- `g::Graph`: graph for which to find the total number of leaves.
"""
function count_leaves(g::G) where {G<:AbstractGraph}
    leaves = collect(Leaves(g))
    unique!(x -> x.id, leaves)

    return length(leaves)
end

function count_leaves(graphs::Vector{G}) where {G<:AbstractGraph}
    leaves = Vector{G}()
    for g in graphs
        append!(leaves, collect(Leaves(g)))
    end
    unique!(x -> x.id, leaves)

    return length(leaves)
end

"""
    function count_operation(g::G) where {G<:AbstractGraph}

    Returns the total number of  additions and multiplications in the graph.

# Arguments:
- `g::Graph`: graph for which to find the total number of  operations.
"""
function count_operation(g::G) where {G<:AbstractGraph}
    visited = Set{Int}()
    totalsum = 0
    totalprod = 0
    # totalpower = 0
    for node in PreOrderDFS(g)
        if !(node.id in visited)
            push!(visited, node.id)
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

"""
    function count_expanded_operation(g::G) where {G<:AbstractGraph}

    Returns the total number of operations in the totally expanded version (without any parentheses in the mathematical expression) of the graph.

# Arguments:
- `g::Graph`: graph for which to find the total number of operations in its expanded version.
"""
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