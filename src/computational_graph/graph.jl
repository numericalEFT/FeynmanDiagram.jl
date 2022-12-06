abstract type Operator end
struct Sum <: Operator end
struct Prod <: Operator end
Base.isequal(a::Operator, b::Operator) = (typeof(a) == typeof(b))
Base.:(==)(a::Operator, b::Operator) = Base.isequal(a, b)
apply(o::Operator, diags) = error("not implemented!")

Base.show(io::IO, o::Operator) = print(io, typeof(o))
Base.show(io::IO, ::Type{Sum}) = print(io, "⨁")
Base.show(io::IO, ::Type{Prod}) = print(io, "Ⓧ")

"""Type alias for a directed graph edge e = (a₁⁺, a₂⁻) from e[1] to e[2]."""
const EdgeType = Tuple{QuantumOperator,QuantumOperator}

"""
    mutable struct Graph{F,W}
    
    Computational Graph representation of a collection of Feynman diagrams. All Feynman diagrams should share the same set of external and internal vertices.

# Members:
- `id::Int`  the unique hash id to identify the diagram
- `name::Symbol`  name of the diagram
- `type::Symbol`  type of the diagram, support :propagator, :interaction, :sigma, :green, :generic
- `orders::Vector{Int}`  orders of the diagram, e.g. loop order, derivative order, etc.
- `external::Vector{Int}`  index of external vertices
- `vertices::Vector{OperatorProduct}`  vertices of the diagram. Each index is composited by the product of quantum operators.
- `topology::Vector{Vector{Int}}` topology of the diagram. Each Vector{Int} stores vertices' index connected with each other (as a propagator). 
- `subgraph::Vector{Graph{F,W}}`  vector of sub-diagrams 
- `operator::Operator`  node operation, support Sum() and Prod()
- `factor::F`  additional factor of the diagram
- `weight::W`  weight of the diagram

# Example:
```julia-repl
julia> g = Graph([𝑓⁺(1)𝑓⁻(2), 𝑓⁺(3)𝑓⁻(4)], external=[1, 2], subgraph=[Graph([𝑓⁺(1)𝑓⁻(4)], []), Graph([𝑓⁻(2)𝑓⁺(3)], [])])
3:f⁺(1)f⁻(2)|f⁺(3)f⁻(4)=0.0=⨁ (1,2)

julia> g.subgraph
2-element Vector{Graph{Float64, Float64}}:
 1:f⁺(1)f⁻(4)=0.0
 2:f⁻(2)f⁺(3)=0.0
```
"""
mutable struct Graph{F,W} # Graph
    id::Int
    name::String # "" by default
    type::Symbol # :propagator, :interaction, :sigma, :green, :generic
    orders::Vector{Int}

    external::Vector{Int} # index of external vertices
    vertices::Vector{OperatorProduct} # vertices of the diagram
    topology::Vector{Vector{Int}}

    subgraph::Vector{Graph{F,W}}

    operator::DataType
    factor::F
    weight::W

    """
        function Graph(vertices::Vector{OperatorProduct}; external=[], subgraph=[],
            name="", type=:generic, operator::Operator=Sum(), orders=zeros(Int, 16),
            ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype))
        
        Create a Graph struct from vertices and external indices.

    # Arguments:
    - `vertices::Vector{OperatorProduct}`  vertices of the diagram
    - `external`  index of external vertices, by default, all vertices are external
    - `topology` topology of the diagram
    - `subgraph`  vector of sub-diagrams 
    - `name`  name of the diagram
    - `type`  type of the diagram
    - `operator::Operator`  node operation
    - `orders`  orders of the diagram
    - `ftype`  typeof(factor)
    - `wtype`  typeof(weight)
    - `factor`  additional factor of the diagram
    - `weight`  weight of the diagram
    """
    function Graph(vertices::AbstractVector; external=collect(1:length(vertices)), subgraph=[], topology=[],
        name="", type=:generic, operator::Operator=Sum(), orders=zeros(Int, 16),
        ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype)
    )
        vertices = [OperatorProduct(v) for v in vertices]
        return new{ftype,wtype}(uid(), name, type, orders, external, vertices, topology, subgraph, typeof(operator), factor, weight)
    end

    """
        function Graph(extV::AbstractVector, intV::AbstractVector; subgraph=[],
            name="", type=:generic, operator::Operator=Sum(), orders=zeros(Int, 16),
            ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype))
        
        Create a Graph struct from external and internal vertices.

    # Arguments:
    - `extV::AbstractVector`  external vertices of the diagram
    - `intV::AbstractVector`  internal vertices of the diagram
    - `topology` topology of the diagram
    - `subgraph`  vector of sub-diagrams 
    - `name`  name of the diagram
    - `type`  type of the diagram
    - `operator::Datatype`  node operation, Sum, Prod, etc.
    - `orders`  orders of the diagram
    - `ftype`  typeof(factor)
    - `wtype`  typeof(weight)
    - `factor`  additional factor of the diagram
    - `weight`  weight of the diagram
    """
    function Graph(extV::AbstractVector, intV::AbstractVector; topology=[], subgraph=[],
        name="", type=:generic, operator::Operator=Sum(), orders=zeros(Int, 16),
        ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype)
    )
        vertices = [extV..., intV...]
        ext = collect(1:length(extV))
        return new{ftype,wtype}(uid(), name, type, orders, ext, vertices, topology, subgraph, typeof(operator), factor, weight)
    end
end

function Base.isequal(a::Graph, b::Graph)
    typeof(a) != typeof(b) && return false
    for field in fieldnames(typeof(a))
        if field == :weight
            (getproperty(a, :weight) ≈ getproperty(b, :weight)) == false && return false
        else
            getproperty(a, field) != getproperty(b, field) && return false
        end
    end
    return true
end
Base.:(==)(a::Graph, b::Graph) = Base.isequal(a, b)
# isbare(diag::Graph) = isempty(diag.subgraph)

"""
    function isequiv(a::Graph, b::Graph, args...)

    Determine whether `a` is equivalent to `b` without considering fields in `args`.
"""
function isequiv(a::Graph, b::Graph, args...)
    typeof(a) != typeof(b) && return false
    for field in fieldnames(typeof(a))
        field in [args...] && continue
        if field == :weight
            (getproperty(a, :weight) ≈ getproperty(b, :weight)) == false && return false
        elseif field == :subgraph
            !all(isequiv.(getproperty(a, field), getproperty(b, field), args...)) && return false
        else
            getproperty(a, field) != getproperty(b, field) && return false
        end
    end
    return true
end

"""
    function is_external(g::Graph, i::Int) 

    Check if `i::Int` in the external indices of Graph `g`.
"""
is_external(g::Graph, i::Int) = i in g.external

"""
    function is_internal(g::Graph, i::Int) 

    Check if `i::Int` in the internal indices of Graph `g`.
"""
is_internal(g::Graph, i::Int) = (i in g.external) == false

"""
    function external_vertices(g::Graph)

    Return all external vertices (::Vector{OperatorProduct}) of Graph `g`.
"""
external_vertices(g::Graph) = g.vertices[g.external]

"""
    function internal_vertices(g::Graph)

    Return all internal vertices (::Vector{OperatorProduct}) of Graph `g`.
"""
internal_vertices(g::Graph) = g.vertices[setdiff(1:length(g.vertices), g.external)]

"""
    function vertices(g::Graph)

    Return all vertices (::Vector{OperatorProduct}) of Graph `g`.
"""
vertices(g::Graph) = g.vertices

#TODO: add function return reducibility of Graph. 
function reducibility(g::Graph)
    return (OneFermiIrreducible,)
end

#TODO: add function for connected diagram check. 
function connectivity(g::Graph)
    isempty(g.subgraph) && return true
end

function Base.:*(g1::Graph{F,W}, c2::C) where {F,W,C}
    return Graph(g1.vertices; external=g1.external, type=g1.type, topology=g1.topology, subgraph=[g1,], operator=Prod(), ftype=F, wtype=W, factor=F(c2))
end

function Base.:*(c1::C, g2::Graph{F,W}) where {F,W,C}
    return Graph(g2.vertices; external=g2.external, type=g2.type, topology=g2.topology, subgraph=[g2,], operator=Prod(), ftype=F, wtype=W, factor=F(c1))
end

function Base.:+(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
    @assert g1.type == g2.type "g1 and g2 are not of the same type."
    # TODO: more check
    @assert Set(vertices(g1)) == Set(vertices(g2)) "g1 and g2 have different vertices."
    @assert Set(external_vertices(g1)) == Set(external_vertices(g2)) "g1 and g2 have different external vertices."
    @assert g1.orders == g2.orders "g1 and g2 have different orders."

    return Graph(g1.vertices; external=g1.external, type=g1.type, subgraph=[g1, g2], operator=Sum(), ftype=F, wtype=W)
end

function Base.:-(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
    return g1 + (-F(1)) * g2
end

"""
    function feynman_diagram(vertices::Vector{OperatorProduct}, contractions::Vector{Int};
        external=[], factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", type=:generic)
    
    Create a Graph representing feynman diagram from all vertices and Wick contractions.

# Arguments:
- `vertices::Vector{OperatorProduct}`  vertices of the diagram
- `contractions::Vector{Int}`  contraction-index vector respresnting Wick contractions
- `external`  index of external vertices
- `factor`  additional factor of the diagram
- `weight`  weight of the diagram
- `name`  name of the diagram
- `type`  type of the diagram

# Example:
```julia-repl
julia> g = feynman_diagram([𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6)], [1, 2, 3, 2, 1, 3])
1: generic graph from f⁺(1)f⁻(2)ϕ(3)|f⁺(4)f⁻(5)ϕ(6)

julia> g.subgraph
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
        push!(g.subgraph, propagator(reduce(*, edge)))
    end
    return g
end
# function feynman_diagram(graphs::Vector{Graph{F,W}}, contractions::Vector{Int};
#     external=[], factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", type=:generic) where {F,W}

#     vertices = [v for g in graphs for v in external_vertices(g)]
#     return feynman_diagram(vertices, contractions; external=external, factor=factor, weight=weight, name=name, type=type)
# end

"""
    function feynman_diagram(vertices::Vector{OperatorProduct}, topology::Vector{Vector{Int}};
        external=[], factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", type=:generic)
    
    Create a Graph representing feynman diagram from all vertices and topology (connections between vertices).

# Arguments:
- `vertices::Vector{OperatorProduct}`  vertices of the diagram
- `topology::Vector{Vector{Int}}` topology of the diagram. Each Vector{Int} stores vertices' index connected with each other (as a propagator). 
- `external`  index of external vertices
- `factor`  additional factor of the diagram
- `weight`  weight of the diagram
- `name`  name of the diagram
- `type`  type of the diagram

# Example:
```julia-repl
julia> g = feynman_diagram([𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6)], [[5, 1], [2, 4], [3, 6]])
1: generic graph from f⁺(1)f⁻(2)ϕ(3)|f⁺(4)f⁻(5)ϕ(6)

julia> g.subgraph
3-element Vector{Graph{Float64, Float64}}:
 2: propagtor graph from f⁻(5)f⁺(1)
 3: propagtor graph from f⁻(2)f⁺(4)
 4: propagtor graph from ϕ(3)ϕ(6)
```
"""
function feynman_diagram(vertices::Vector{OperatorProduct}, topology::Vector{Vector{Int}};
    external=[], factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", type=:generic)

    operators = [o for v in vertices for o in v.operators]
    permutation = collect(Iterators.flatten(topology))
    filter!(p -> p ∉ findall(x -> !x, isfermionic.(operators)), permutation)
    sign = isempty(permutation) ? 1 : parity(sortperm(permutation))

    g = Graph(vertices; external=external, topology=topology, name=name, type=type, operator=Prod(),
        factor=factor * sign, weight=weight)
    for connection in topology
        push!(g.subgraph, propagator(reduce(*, operators[connection])))
    end
    return g
end
# function feynman_diagram(graphs::Vector{Graph{F,W}}, topology::Vector{Vector{Int}};
#     external=[], factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", type=:generic) where {F,W}

#     vertices = reduce(*, [v for g in graphs for v in external_vertices(g)])
#     return feynman_diagram(vertices, topology; external=external, factor=factor, weight=weight, name=name, type=type)
# end

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

"""
    function propagator(ops::OperatorProduct;
        name="", diagtype=:propagtor, factor=one(_dtype.factor), weight=zero(_dtype.weight), operator=Sum())

    Create a propagator-type Graph from given OperatorProduct `ops`.
"""
function propagator(ops::OperatorProduct;
    name="", diagtype=:propagtor, factor=one(_dtype.factor), weight=zero(_dtype.weight), operator=Sum())
    return Graph([ops], []; type=diagtype, name=name, operator=operator, factor=factor, weight=weight)
end

"""
    function standardize_order!(g::Graph)

    Standardize the order of all leaves (propagators) of Graph by correlator ordering.

# Example: 
```julia-repl
julia> g = propagator(𝑓⁺(1)𝑏⁺(2)𝜙(3)𝑓⁻(1)𝑏⁻(2))
1: propagtor graph from f⁺(1)b⁺(2)ϕ(3)f⁻(1)b⁻(2)

julia> standardize_order!(g)

julia> g, g.factor
(1: propagtor graph from f⁻(1)b⁻(2)ϕ(3)b⁺(2)f⁺(1), -1.0)
```
"""
function standardize_order!(g::Graph)
    for leaf in Leaves(g)
        for (i, vertex) in enumerate(leaf.vertices)
            sign, newvertex = correlator_order(vertex)
            leaf.vertices[i] = OperatorProduct(newvertex)
            leaf.factor *= sign
        end
    end
end

#####################  interface to AbstractTrees ########################### 
function AbstractTrees.children(diag::Graph)
    return diag.subgraph
end

## Things that make printing prettier
AbstractTrees.printnode(io::IO, diag::Graph) = print(io, "\u001b[32m$(diag.id)\u001b[0m : $diag")
AbstractTrees.nodetype(::Graph{F,W}) where {F,W} = Graph{F,W}

## Optional enhancements
# These next two definitions allow inference of the item type in iteration.
# (They are not sufficient to solve all internal inference issues, however.)
Base.IteratorEltype(::Type{<:TreeIterator{Graph{F,W}}}) where {F,W} = Base.HasEltype()
Base.eltype(::Type{<:TreeIterator{Graph{F,W}}}) where {F,W} = Graph{F,W}
