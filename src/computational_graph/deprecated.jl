module Deprecated

# """Type alias for a directed graph edge e = (a₁⁺, a₂⁻) from e[1] to e[2]."""
const EdgeType = Tuple{QuantumOperator,QuantumOperator}

function Base.:+(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
    @assert g1.type == g2.type "g1 and g2 are not of the same type."
    # TODO: more check
    @assert Set(vertices(g1)) == Set(vertices(g2)) "g1 and g2 have different vertices."
    @assert Set(external(g1)) == Set(external(g2)) "g1 and g2 have different external vertices."
    @assert g1.orders == g2.orders "g1 and g2 have different orders."

    return Graph(g1.vertices; external=g1.external, type=g1.type, subgraphs=[g1, g2], operator=Sum(), ftype=F, wtype=W)
end

function Base.:-(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
    return g1 + (-F(1)) * g2
end

"""
    function internal_vertices(g::Graph)

    Return all internal vertices (::Vector{OperatorProduct}) of Graph `g`.
"""
internal_vertices(g::Graph) = g.vertices[setdiff(eachindex(g.vertices), g.external)]

"""
    function feynman_diagram(vertices::Vector{OperatorProduct}, contractions::Vector{Int};
        external=[], factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", type=:generic)

    Create a Graph representing feynman diagram from all vertices and Wick contractions.

# Arguments:
- `vertices::Vector{OperatorProduct}`  vertices of the diagram
- `contractions::Vector{Int}`  contraction-index vector respresnting Wick contractions
- `external`  index of external vertices
- `factor`  scalar multiplicative factor for the diagram
- `weight`  weight of the diagram
- `name`  name of the diagram
- `type`  type of the diagram

# Example:
```julia-repl
julia> g = feynman_diagram([𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6)], [1, 2, 3, 2, 1, 3])
1: generic graph from f⁺(1)f⁻(2)ϕ(3)|f⁺(4)f⁻(5)ϕ(6)

julia> g.subgraphs
3-element Vector{Graph{Float64, Float64}}:
 2: propagtor graph from f⁺(1)f⁻(5)
 3: propagtor graph from f⁻(2)f⁺(4)
 4: propagtor graph from ϕ(3)ϕ(6)
```
"""
function feynman_diagram(vertices::Vector{OperatorProduct}, contractions::Vector{Int};
    external=[], factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", type=:generic)

    contraction_sign, topology, edges = contractions_to_edges(vertices, contractions)
    g = Graph(vertices; external=external, topology=topology, name=name, type=type, operator=Prod(),
        factor=factor * contraction_sign, weight=weight)
    for edge in edges
        push!(g.subgraphs, propagator(reduce(*, edge)))
    end
    return g
end
function feynman_diagram(graphs::Vector{Graph{F,W}}, contractions::Vector{Int};
    external=[], factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", type=:generic) where {F,W}

    vertices = [v for g in graphs for v in external_vertices(g)]
    return feynman_diagram(vertices, contractions; external=external, factor=factor, weight=weight, name=name, type=type)
end

function feynman_diagram(graphs::Vector{Graph{F,W}}, topology::Vector{Vector{Int}};
    external=[], factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", type=:generic) where {F,W}

    vertices = reduce(*, [v for g in graphs for v in external_vertices(g)])
    return feynman_diagram(vertices, topology; external=external, factor=factor, weight=weight, name=name, type=type)
end

"""
    function contractions_to_edges(vertices::Vector{OperatorProduct}, contractions::Vector{Int})

    Converts a list of Wick contractions associated with a list of vertices
    to a list of edges e = (i, j) directed from `operators[i]` to `operators[j]`, where
    `operators` is a flattened list of operators associated with the specified `vertices`.

# Example: 
```julia-repl
julia> vertices = [𝑏⁺(1)𝑓⁺(2)𝜙(3), 𝑓⁻(4)𝑓⁻(5), 𝑏⁻(6)𝑓⁺(7)𝜙(8)];

julia> edges, sign = contractions_to_edges(vertices, [1, 2, 3, 2, 4, 1, 4, 3])
(Tuple{QuantumOperator, QuantumOperator}[(b⁺(1), b⁻(6)), (f⁺(2), f⁻(4)), (ϕ(3), ϕ(8)), (f⁻(5), f⁺(7))], 1)
```
- flattened fermionic edges: [2, 4, 5, 7]
- sign = parity([1, 2, 3, 4]) = 1
"""
function contractions_to_edges(vertices::Vector{OperatorProduct}, contractions::Vector{Int})
    #TODO: only works for weak-coupling expansion with Wick's theorem for now.
    # Obtain the flattened list of non-composite operators
    operators = [o for v in vertices for o in v.operators]
    # Filter some illegal contractions
    is_invalid =
        isodd(length(contractions)) ||
        # isodd(length(unique(contractions))) ||
        length(operators) != length(contractions)
    if is_invalid
        throw(
            ArgumentError(
                "Input $contractions does not specify a legal set of Wick contractions",
            ),
        )
    end
    # Loop over operators and pair
    next_pairing = 1
    edges = Vector{EdgeType}()
    topology = Vector{Int}[]
    permutation = Int[]

    for (i, wick_i) in enumerate(contractions)
        if i < next_pairing
            continue  # Operator already paired
        end
        for (js, wick_j) in enumerate(contractions[(i+1):end])
            j = i + js  # Iterating for (j > i)
            if j < next_pairing
                continue  # Operator already paired
            end
            if wick_i == wick_j
                @debug "Found Wick contraction #$wick_i, adding edge between operators $i and $j"
                @assert operators[j]'.operator == operators[i].operator
                isfermionic(operators[j]) && append!(permutation, [i, j])
                push!(edges, (operators[i], operators[j]))
                push!(topology, [i, j])
                # Move on to next pair
                next_pairing += 1
                break
            end
        end
    end
    # Deduce the parity of the contraction
    # Get the permutation parity for the flattened list of fermionic edges, e.g.,
    # permutation = [1, 5, 4, 8, 7, 6] => P = (1 3 2 6 5 4) => sign = +1
    sign = isempty(permutation) ? 1 : parity(sortperm(permutation))

    return sign, topology, edges
end

end # module Deprecated