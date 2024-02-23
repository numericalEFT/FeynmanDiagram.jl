abstract type DiagramType end
struct Interaction <: DiagramType end
struct ExternalVertex <: DiagramType end
struct Propagator <: DiagramType end
struct SelfEnergy <: DiagramType end
struct VertexDiag <: DiagramType end
struct GreenDiag <: DiagramType end
struct GenericDiag <: DiagramType end

"""
    mutable struct FeynmanProperties

    Diagrammatic properties associated with a given Feynman diagram.

# Members:
- `diagtype::DataType`  classification of the Feynman diagram. Should be one of the following supported DiagramTypes: Interaction, ExternalVertex, Propagator, SelfEnergy, VertexDiag, GreenDiag, or GenericDiag.
- `vertices::Vector{OperatorProduct}`  vertices of the diagram. Each index is composited by the product of quantum operators. 
- `topology::Vector{Vector{Int}}` topology of the diagram. Each Vector{Int} stores vertices' index connected with each other (as a propagator). 
- `external_indices::Vector{Int}`  indices of actual external vertices in terms of QuantumOperators
- `external_legs::Vector{Bool}` indicates which external vertices have real legs (true: real leg, false: fake leg)
"""
# TODO: add additional properties, e.g., isconnected::Bool and isreducible::Bool
mutable struct FeynmanProperties
    diagtype::DataType # :propagator, :interaction, :sigma, :green, :generic, ...
    vertices::Vector{OperatorProduct}
    topology::Vector{Vector{Int}}
    external_indices::Vector{Int}
    external_legs::Vector{Bool}
end

function Base.isequal(a::FeynmanProperties, b::FeynmanProperties)
    for field in fieldnames(FeynmanProperties)
        getproperty(a, field) != getproperty(b, field) && return false
    end
    return true
end
Base.:(==)(a::FeynmanProperties, b::FeynmanProperties) = Base.isequal(a, b)

"""
    function drop_topology(p::FeynmanProperties)

    Returns a copy of the given FeynmanProperties `p` modified to have no topology.
"""
drop_topology(p::FeynmanProperties) = FeynmanProperties(p.diagtype, p.vertices, [], p.external_indices, p.external_legs)
"""
    mutable struct FeynmanGraph{F<:Number,W}
    
    Computational graph representation of a (collection of) Feynman diagram(s). All Feynman diagrams should share the same set of external and internal vertices.

# Members:
- `id::Int`  the unique hash id to identify the diagram
- `name::Symbol`  name of the diagram
- `orders::Vector{Int}`  orders associated with the Feynman graph, e.g., loop/derivative orders
- `properties::FeynmanProperties`  diagrammatic properties, e.g., the operator vertices and topology
- `subgraphs::Vector{FeynmanGraph{F,W}}`  vector of sub-diagrams 
- `subgraph_factors::Vector{F}`  scalar multiplicative factors associated with each subdiagram
- `operator::DataType`  node operation (Sum, Prod, etc.)
- `weight::W`  weight of the diagram

# Example:
```julia-repl
julia> g1 = FeynmanGraph([]; vertices=[𝑓⁺(1),𝑓⁻(2)], external_indices=[1,2], external_legs=[true,true])
1:f⁺(1)|f⁻(2)=0.0

julia> g2 = FeynmanGraph([]; vertices=[𝑓⁺(3),𝑓⁻(4)], external_indices=[1,2], external_legs=[true,true])
2:f⁺(3)|f⁻(4)=0.0

julia> g = FeynmanGraph([g1,g2]; vertices=[𝑓⁺(1),𝑓⁻(2),𝑓⁺(3),𝑓⁻(4)], operator=ComputationalGraphs.Prod(), external_indices=[1,2,3,4], external_legs=[true,true,true,true])
3:f⁺(1)|f⁻(2)|f⁺(3)|f⁻(4)=0.0=Ⓧ (1,2)
```
"""
mutable struct FeynmanGraph{F<:Number,W} <: AbstractGraph # FeynmanGraph
    id::Int
    name::String # "" by default
    orders::Vector{Int}
    properties::FeynmanProperties

    subgraphs::Vector{FeynmanGraph{F,W}}
    subgraph_factors::Vector{F}

    operator::DataType
    weight::W

    """
        function FeynmanGraph(subgraphs::AbstractVector; topology=[], vertices::Union{Vector{OperatorProduct},Nothing}=nothing, external_indices=[], external_legs=[],
            subgraph_factors=one.(eachindex(subgraphs)), name="", diagtype::DiagramType=GenericDiag(), operator::AbstractOperator=Sum(),
            orders=zeros(Int, 16), ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype))
        
        Create a FeynmanGraph struct from a set of subgraphs, vertices and external indices.

    # Arguments:
    - `subgraphs`  vector of sub-diagrams 
    - `topology` topology of the diagram
    - `vertices::Union{Vector{OperatorProduct},Nothing}`  vertices of the diagram, nothing by default
    - `external_indices`  indices of actual external vertices in terms of QuantumOperators, empty by default
    - `external_legs` indicates which external indices correspond to real legs (true: real leg, false: fake leg)
    - `subgraph_factors`  scalar multiplicative factors associated with each subdiagram
    - `name`  name of the diagram
    - `diagtype::DiagramType`  type of the diagram
    - `operator::AbstractOperator`  node operation (Sum, Prod, etc.)
    - `orders`  orders associated with the Feynman graph, e.g., loop/derivative orders
    - `ftype`  typeof(factor)
    - `wtype`  typeof(weight)
    - `factor`  overall scalar multiplicative factor for this diagram (e.g., permutation sign)
    - `weight`  weight of the diagram
    """
    function FeynmanGraph(subgraphs::AbstractVector; topology=[], vertices::Union{Vector{OperatorProduct},Nothing}=nothing, external_indices=[], external_legs=[],
        subgraph_factors=one.(eachindex(subgraphs)), name="", diagtype::DiagramType=GenericDiag(), operator::AbstractOperator=Sum(),
        orders=zeros(Int, 16), ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype)
    )
        @assert length(external_indices) == length(external_legs)
        if typeof(operator) <: Power
            @assert length(subgraphs) == 1 "FeynmanGraph with Power operator must have one and only one subgraph."
        end
        # @assert allunique(subgraphs) "all subgraphs must be distinct."
        if isnothing(vertices)
            vertices = [external_operators(g) for g in subgraphs if diagram_type(g) != Propagator]
        end
        properties = FeynmanProperties(typeof(diagtype), vertices, topology, external_indices, external_legs)

        g = new{ftype,wtype}(uid(), String(name), orders, properties, subgraphs, subgraph_factors, typeof(operator), weight)
        if factor ≈ one(ftype)
            return g
        else
            return new{ftype,wtype}(uid(), String(name), orders, properties, [g,], [factor,], Prod, weight * factor)
        end
    end

    """
        function FeynmanGraph(subgraphs::AbstractVector, properties::FeynmanProperties;
            subgraph_factors=one.(eachindex(subgraphs)), name="", operator::AbstractOperator=Sum(),
            ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype))
        
        Create a FeynmanGraph struct from a given set of subgraphs and Feynman properties.

    # Arguments:
    - `subgraphs`  vector of sub-diagram
    - `properties::FeynmanProperties`  diagrammatic properties, e.g., the operator vertices and topology 
    - `subgraph_factors`  scalar multiplicative factors associated with each subdiagram
    - `name`  name of the diagram
    - `operator::AbstractOperator`  node operation (Sum, Prod, etc.)
    - `ftype`  typeof(factor)
    - `wtype`  typeof(weight)
    - `factor`  overall scalar multiplicative factor for this diagram (e.g., permutation sign)
    - `weight`  weight of the diagram
    """
    function FeynmanGraph(subgraphs::AbstractVector, properties::FeynmanProperties;
        subgraph_factors=one.(eachindex(subgraphs)), name="", operator::AbstractOperator=Sum(),
        orders=zeros(Int, 16), ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype)
    )
        @assert length(properties.external_indices) == length(properties.external_legs)
        if typeof(operator) <: Power
            @assert length(subgraphs) == 1 "FeynmanGraph with Power operator must have one and only one subgraph."
        end
        # @assert allunique(subgraphs) "all subgraphs must be distinct."
        g = new{ftype,wtype}(uid(), String(name), orders, properties, subgraphs, subgraph_factors, typeof(operator), weight)
        if factor ≈ one(ftype)
            return g
        else
            return new{ftype,wtype}(uid(), String(name), orders, properties, [g,], [factor,], Prod, weight * factor)
        end
    end

    """
        function FeynmanGraph(g::Graph, properties::FeynmanProperties)

        Create a Feynman graph given a graph `g` and the Feynman properties (external vertices, topology, etc.) to endow it with.

    # Arguments:
    - `g`  computational graph
    - `properties::FeynmanProperties`  diagrammatic properties, e.g., the operator vertices and topology 
    """
    function FeynmanGraph(g::Graph{F,W}, properties::FeynmanProperties) where {F,W}
        @assert length(properties.external_indices) == length(properties.external_legs)
        # @assert allunique(subgraphs) "all subgraphs must be distinct."
        # return new{F,W}(uid(), g.name, g.orders, properties, g.subgraphs, g.subgraph_factors, g.operator, g.weight)
        return new{F,W}(uid(), g.name, g.orders, properties, [FeynmanGraph(subg, subg.properties) for subg in g.subgraphs], g.subgraph_factors, g.operator, g.weight)
    end
end

### AbstractGraph interface for FeynmanGraph ###

# Getters
id(g::FeynmanGraph) = g.id
name(g::FeynmanGraph) = g.name
orders(g::FeynmanGraph) = g.orders
operator(g::FeynmanGraph) = g.operator
weight(g::FeynmanGraph) = g.weight
subgraph(g::FeynmanGraph, i=1) = g.subgraphs[i]
subgraphs(g::FeynmanGraph) = g.subgraphs
subgraphs(g::FeynmanGraph, indices::AbstractVector{Int}) = g.subgraphs[indices]
subgraph_factor(g::FeynmanGraph, i=1) = g.subgraph_factors[i]
subgraph_factors(g::FeynmanGraph) = g.subgraph_factors
subgraph_factors(g::FeynmanGraph, indices::AbstractVector{Int}) = g.subgraph_factors[indices]

# Setters
set_id!(g::FeynmanGraph, id::Int) = (g.id = id)
set_name!(g::FeynmanGraph, name::String) = (g.name = name)
set_orders!(g::FeynmanGraph, orders::Vector{Int}) = (g.orders = orders)
set_operator!(g::FeynmanGraph, operator::Type{<:AbstractOperator}) = (g.operator = operator)
set_operator!(g::FeynmanGraph, operator::AbstractOperator) = (g.operator = typeof(operator))
set_weight!(g::FeynmanGraph{F,W}, weight) where {F,W} = (g.weight = W(weight))
set_subgraph!(g::FeynmanGraph{F,W}, subgraph::FeynmanGraph{F,W}, i=1) where {F,W} = (g.subgraphs[i] = subgraph)
set_subgraphs!(g::FeynmanGraph{F,W}, subgraphs::Vector{FeynmanGraph{F,W}}) where {F,W} = (g.subgraphs = subgraphs)
set_subgraphs!(g::FeynmanGraph{F,W}, subgraphs::Vector{FeynmanGraph{F,W}}, indices::AbstractVector{Int}) where {F,W} = (g.subgraphs[indices] = subgraphs)
set_subgraph_factor!(g::FeynmanGraph{F,W}, subgraph_factor, i=1) where {F,W} = (g.subgraph_factors[i] = F(subgraph_factor))
set_subgraph_factors!(g::FeynmanGraph{F,W}, subgraph_factors::AbstractVector) where {F,W} = (g.subgraph_factors = Vector{F}(subgraph_factors))
set_subgraph_factors!(g::FeynmanGraph{F,W}, subgraph_factors::AbstractVector, indices::AbstractVector{Int}) where {F,W} = (g.subgraph_factors[indices] = Vector{F}(subgraph_factors))

###############################

"""
function is_external_operators(g::FeynmanGraph, i) 

Check if `i` in the external indices of FeynmanGraph `g`.
"""
is_external(g::FeynmanGraph, i) = i in g.properties.external_indices

"""
    function is_internal(g::FeynmanGraph, i) 

Check if `i` in the internal indices of FeynmanGraph `g`.
"""
is_internal(g::FeynmanGraph, i) = (i in g.properties.external_indices) == false

"""
    function diagram_type(g::FeynmanGraph)

Returns the diagram type (::DiagramType) of FeynmanGraph `g`.
"""
diagram_type(g::FeynmanGraph) = g.properties.diagtype

"""
    function vertices(g::FeynmanGraph)

Returns all vertices (::Vector{OperatorProduct}) of FeynmanGraph `g`.
"""
vertices(g::FeynmanGraph) = g.properties.vertices

"""
    function vertex(g::FeynmanGraph, i=1)

Returns the `i`th vertex (::OperatorProduct) of FeynmanGraph `g`.
Defaults to the first vertex if an index `i` is not supplied.
"""
vertex(g::FeynmanGraph, i=1) = g.properties.vertices[i]

"""
    function topology(g::FeynmanGraph)

Returns the topology (::Vector{Vector{Int}}) of FeynmanGraph `g`.
"""
topology(g::FeynmanGraph) = g.properties.topology

"""
    function external_legs(g::FeynmanGraph)

Returns a list of Boolean indices external_legs (::Vector{Bool}) indicating which external vertices of FeynmanGraph `g` have real legs (true: real leg, false: fake leg).
"""
external_legs(g::FeynmanGraph) = g.properties.external_legs

"""
    function external_indices(g::FeynmanGraph)

Returns a list of indices (::Vector{Int}}) to the external vertices of the FeynmanGraph `g`.
"""
external_indices(g::FeynmanGraph) = g.properties.external_indices

"""
    function external_operators(g::FeynmanGraph)

Returns all physical external operators (::OperatorProduct}) of FeynmanGraph `g`.
"""
external_operators(g::FeynmanGraph) = OperatorProduct(OperatorProduct(g.properties.vertices)[g.properties.external_indices])

"""
    function external_labels(g::FeynmanGraph)

Returns the labels of all physical external vertices of FeynmanGraph `g`.
"""
external_labels(g::FeynmanGraph) = [o.label for o in external_operators(g)]

function reducibility(g::FeynmanGraph)
    #TODO: add function return reducibility of FeynmanGraph. 
    @todo
    return (OneFermiIrreducible,)
end

function connectivity(g::FeynmanGraph)
    #TODO: add function for connected diagram check. 
    @todo
    isleaf(g) && return true
end

"""
    function Base.:*(g1::Graph{F,W}, c2) where {F,W}

Returns a graph representing the scalar multiplication `g1*c2`.

# Arguments:
- `g1`  Feynman graph
- `c2`  scalar multiple
"""
function Base.:*(g1::FeynmanGraph{F,W}, c2) where {F,W}
    g = FeynmanGraph([g1,], g1.properties; subgraph_factors=[F(c2),], operator=Prod(), orders=orders(g1), ftype=F, wtype=W)
    # Convert trivial unary link to in-place form
    if unary_istrivial(g1) && onechild(g1)
        g.subgraph_factors[1] *= g1.subgraph_factors[1]
        g.subgraphs = g1.subgraphs
    end
    return g
end

"""
    function Base.:*(c1, g2::Graph{F,W}) where {F,W}

Returns a graph representing the scalar multiplication `c1*g2`.

# Arguments:
- `c1`  scalar multiple
- `g2`  Feynman graph
"""
function Base.:*(c1, g2::FeynmanGraph{F,W}) where {F,W}
    g = FeynmanGraph([g2,], g2.properties; subgraph_factors=[F(c1),], operator=Prod(), orders=orders(g2), ftype=F, wtype=W)
    # Convert trivial unary link to in-place form
    if unary_istrivial(g2) && onechild(g2)
        g.subgraph_factors[1] *= g2.subgraph_factors[1]
        g.subgraphs = g2.subgraphs
    end
    return g
end

"""
    function linear_combination(g1::FeynmanGraph{F,W}, g2::FeynmanGraph{F,W}, c1, c2) where {F,W}

Returns a graph representing the linear combination `c1*g1 + c2*g2`. If `g1 == g2`, it will return a graph representing `(c1+c2)*g1`
Feynman Graphs `g1` and `g2` must have the same diagram type, orders, and external vertices.

# Arguments:
- `g1`  first Feynman graph
- `g2`  second Feynman graph
- `c1`:  first scalar multiple (defaults to 1).
- `c2`:  second scalar multiple (defaults to 1).
"""
function linear_combination(g1::FeynmanGraph{F,W}, g2::FeynmanGraph{F,W}, c1=F(1), c2=F(1)) where {F,W}
    @assert diagram_type(g1) == diagram_type(g2) "g1 and g2 are not of the same graph type."
    @assert orders(g1) == orders(g2) "g1 and g2 have different orders."
    @assert Set(external_operators(g1)) == Set(external_operators(g2)) "g1 and g2 have different external vertices."
    empty_topology = []  # No topology for Sum nodes
    total_vertices = union(vertices(g1), vertices(g2))
    properties = FeynmanProperties(diagram_type(g1), total_vertices, empty_topology, external_indices(g1), external_legs(g1))

    f1 = typeof(c1) == F ? c1 : F(c1)
    f2 = typeof(c2) == F ? c2 : F(c2)
    subgraphs = [g1, g2]
    subgraph_factors = [f1, f2]
    # Convert trivial unary links to in-place form
    if unary_istrivial(g1) && onechild(g1)
        subgraph_factors[1] *= g1.subgraph_factors[1]
        subgraphs[1] = g1.subgraphs[1]
    end
    if unary_istrivial(g2) && onechild(g2)
        subgraph_factors[2] *= g2.subgraph_factors[1]
        subgraphs[2] = g2.subgraphs[1]
    end

    if subgraphs[1].id == subgraphs[2].id
        # if isequiv(subgraphs[1], subgraphs[2], :id)
        g = FeynmanGraph([subgraphs[1]], properties; subgraph_factors=[sum(subgraph_factors)], operator=Sum(), orders=orders(g1), ftype=F, wtype=W)
    else
        g = FeynmanGraph(subgraphs, properties; subgraph_factors=subgraph_factors, operator=Sum(), orders=orders(g1), ftype=F, wtype=W)
    end
    return g
end

"""
    function linear_combination(graphs::Vector{FeynmanGraph{F,W}}, constants::AbstractVector=ones(F, length(graphs))) where {F,W}

Given a vector 𝐠 of graphs each with the same type and external/internal vertices and 
an equally-sized vector 𝐜 of constants, returns a new graph representing the linear combination (𝐜 ⋅ 𝐠). 
The function identifies unique graphs from the input `graphs` and sums their associated `constants`.
All input Graphs must have the same diagram type, orders, and external vertices.

# Arguments:
- `graphs`  vector of input FeymanGraphs
- `constants`  vector of scalar multiples (defaults to ones(F, length(graphs))).

# Returns:
- A new `FeynmanGraph{F,W}` object representing the linear combination of the unique input `graphs` weighted by the constants, 
where duplicate graphs in the input `graphs` are combined by summing their associated constants. 

# Example:
    Given graphs `g1`, `g2`, `g1` and constants `c1`, `c2`, `c3`, the function computes `(c1+c3)*g1 + c2*g2`.
"""
function linear_combination(graphs::Vector{FeynmanGraph{F,W}}, constants::AbstractVector=ones(F, length(graphs))) where {F,W}
    @assert alleq(diagram_type.(graphs)) "Graphs are not all of the same graph type."
    @assert alleq(orders.(graphs)) "Graphs do not all have the same order."
    @assert alleq(Set.(external_operators.(graphs))) "Graphs do not share the same set of external vertices."
    g1 = graphs[1]
    empty_topology = []  # No topology for Sum nodes
    total_vertices = union(Iterators.flatten(vertices.(graphs)))
    properties = FeynmanProperties(diagram_type(g1), total_vertices, empty_topology, external_indices(g1), external_legs(g1))

    subgraphs = graphs
    subgraph_factors = eltype(constants) == F ? constants : Vector{F}(constants)
    # Convert trivial unary links to in-place form
    for (i, sub_g) in enumerate(graphs)
        if unary_istrivial(sub_g) && onechild(sub_g)
            subgraph_factors[i] *= sub_g.subgraph_factors[1]
            subgraphs[i] = sub_g.subgraphs[1]
        end
    end
    unique_graphs = FeynmanGraph{F,W}[]
    unique_factors = F[]
    for (idx, g) in enumerate(subgraphs)
        i = findfirst(isequal(g.id), id.(unique_graphs))
        if isnothing(i)
            push!(unique_graphs, g)
            push!(unique_factors, subgraph_factors[idx])
        else
            unique_factors[i] += subgraph_factors[idx]
        end
    end

    g = FeynmanGraph(unique_graphs, properties; subgraph_factors=unique_factors, operator=Sum(), orders=orders(g1), ftype=F, wtype=W)
    return g
end

"""
    function Base.:+(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}

Returns a graph `g1 + g2` representing the addition of two Feynman diagrams `g2` with `g1`.
Diagrams `g1` and `g2` must have the same diagram type, orders, and external vertices.

# Arguments:
- `g1`  first Feynman graph
- `g2`  second Feynman graph
"""
function Base.:+(g1::FeynmanGraph{F,W}, g2::FeynmanGraph{F,W}) where {F,W}
    return linear_combination(g1, g2, F(1), F(1))
end

"""
    function Base.:-(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}

Returns a graph `g1 - g2` representing the subtraction of `g2` from `g1`.
Diagrams `g1` and `g2` must have the same diagram type, orders, and external vertices.

# Arguments:
- `g1`  first Feynman graph
- `g2`  second Feynman graph
"""
function Base.:-(g1::FeynmanGraph{F,W}, g2::FeynmanGraph{F,W}) where {F,W}
    return linear_combination(g1, g2, F(1), F(-1))
end

function Base.:*(g1::FeynmanGraph, g2::FeynmanGraph)
    error("Multiplication of Feynman graphs is not well defined!")
end

"""
    function feynman_diagram(subgraphs::Vector{FeynmanGraph{F,W}}, topology::Vector{Vector{Int}}, perm_noleg::Union{Vector{Int},Nothing}=nothing;
        factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", diagtype::DiagramType=GenericDiag()) where {F,W}
    
Create a FeynmanGraph representing feynman diagram from all subgraphs and topology (connections between vertices),
where each ExternalVertex is given in `vertices`, 
while internal vertices are constructed with external legs of graphs in `vertices`, or simply OperatorProduct in `vertices`.
    
# Arguments:
- `subgraphs::Vector{FeynmanGraph{F,W}}` all subgraphs of the diagram. All external operators of subgraphs constitute all operators of the new diagram.
- `topology::Vector{Vector{Int}}` topology of the diagram. Each Vector{Int} stores operators' index connected with each other (as a propagator). 
- `perm_noleg::Union{Vector{Int},Nothing}=nothing` permutation of all the nonleg external operators. By default, setting nothing means to use the default order from subgraphs.
- `factor::F`  overall scalar multiplicative factor for this diagram (e.g., permutation sign)
- `weight`  weight of the diagram
- `name`  name of the diagram
- `diagtype`  type of the diagram

# Example:
```julia-repl
julia> V = [𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6), 𝑓⁺(7)𝑓⁻(8)𝜙(9)];
julia> g = feynman_diagram(interaction.(V), [[1, 5], [3, 9], [4, 8]], [3, 1, 2])
7:f⁺(1)f⁻(2)ϕ(3)|f⁺(4)f⁻(5)ϕ(6)|f⁺(7)f⁻(8)ϕ(9)=0.0=Ⓧ (1,2,3,4,5,6)

julia> g.subgraphs
6-element Vector{FeynmanGraph{Float64, Float64}}:
 1:f⁺(1)f⁻(2)ϕ(3)=0.0
 2:f⁺(4)f⁻(5)ϕ(6)=0.0
 3:f⁺(7)f⁻(8)ϕ(9)=0.0
 4:f⁺(1)|f⁻(5)⋅-1.0=0.0
 5:ϕ(3)|ϕ(9)=0.0
 6:f⁺(4)|f⁻(8)⋅-1.0=0.0
```
"""
function feynman_diagram(subgraphs::Vector{FeynmanGraph{F,W}}, topology::Vector{Vector{Int}}, perm_noleg::Union{Vector{Int},Nothing}=nothing;
    contraction_orders::Union{Nothing,Vector{Vector{Int}}}=nothing, factor=one(F), weight=zero(W),
    name="", diagtype::DiagramType=GenericDiag(), is_signed::Bool=false) where {F,W}

    # external_ops = OperatorProduct(operators[external]) # the external operators for the building diagram after contractions
    contraction = collect(Iterators.flatten(topology))
    @assert length(unique(contraction)) == length(contraction)  # no repeated index

    vertices, all_external_legs = OperatorProduct[], Bool[]
    external_leg, external_noleg = Int[], Int[] # index all leg/nonleg external operators
    ind = 0

    subgraphs = deepcopy(subgraphs)
    orders_length = length(orders(subgraphs[1]))
    diag_orders = zeros(Int, orders_length)
    for g in subgraphs
        diag_orders += orders(g)
        diagram_type(g) == Propagator && continue  # exclude propagator subgraph to avoid double counting.
        push!(vertices, external_operators(g))
        append!(all_external_legs, external_legs(g))
        if diagram_type(g) == ExternalVertex
            append!(external_leg, external_indices(g) .+ ind) # ExternalVertex will be legged after contraction.
        else
            gext = setdiff(external_indices(g) .+ ind, contraction) # select all external operators
            gextLeg = external_legs(g)[gext.-ind]
            # the selected gext[i] with gextLeg[i]==true is the external vertice with a leg
            append!(external_leg, gext[gextLeg])
            append!(external_noleg, gext[gextLeg.==false])
        end
        ind += length(external_indices(g))
    end

    @assert !any(all_external_legs[setdiff(eachindex(all_external_legs), external_noleg)]) "all contracted operators should have no leg."
    @assert external_leg ⊆ contraction
    @assert isempty(intersect(contraction, external_noleg)) "all nonleg external operators should not be contracted"
    if !isnothing(perm_noleg)
        @assert length(unique(perm_noleg)) == length(perm_noleg) == length(external_noleg)
        external_noleg = external_noleg[perm_noleg]
    end

    operators = OperatorProduct(vertices) # all external operators from subgraphs
    permutation = union(contraction, external_noleg)
    @assert Set(permutation) == Set(eachindex(operators)) # permutation must exhaust all operators

    if !is_signed
        fermionic_operators = isfermionic.(operators)
        filter!(p -> fermionic_operators[p], permutation)
        sign = isempty(permutation) ? 1 : parity(sortperm(permutation))
    else
        sign = 1
    end

    # subgraphs_noVer = FeynmanGraph{F,W}[]
    if isnothing(contraction_orders)
        for (i, connection) in enumerate(topology)
            push!(subgraphs, propagator(operators[connection]; orders=zeros(Int, orders_length)))
            # push!(subgraphs_noVer, propagator(operators[connection]; orders=zeros(Int, orders_length)))
        end
    else
        for (i, connection) in enumerate(topology)
            propagator_orders = zeros(Int, orders_length)
            propagator_orders[eachindex(contraction_orders[i])] = contraction_orders[i]
            push!(subgraphs, propagator(operators[connection]; orders=propagator_orders))
            # push!(subgraphs_noVer, propagator(operators[connection]; orders=propagator_orders))
            diag_orders += propagator_orders
        end
    end
    _external_indices = union(external_leg, external_noleg)
    _external_legs = append!([true for i in eachindex(external_leg)], [false for i in eachindex(external_noleg)])
    # return FeynmanGraph(subgraphs_noVer; topology=topology, external_indices=_external_indices, external_legs=_external_legs, vertices=vertices,
    return FeynmanGraph(subgraphs; topology=topology, external_indices=_external_indices, external_legs=_external_legs, vertices=vertices,
        orders=diag_orders, name=name, diagtype=diagtype, operator=Prod(), factor=factor * sign, weight=weight)
end

# do nothing when already a OperatorProduct; 
_extract_vertex(::Type{<:OperatorProduct}, g) = g
# helper functions extracting external legs from g::FeynmanGraph to form a vertex 
_extract_vertex(::Type{<:FeynmanGraph}, g) = OperatorProduct(external_operators(g))

"""
    function propagator(ops::Union{OperatorProduct,Vector{QuantumOperator}};
        name="", factor=one(_dtype.factor), weight=zero(_dtype.weight), operator=Sum())

Create a Propagator-type FeynmanGraph from given OperatorProduct or Vector{QuantumOperator} `ops`, including two quantum operators.
"""
function propagator(ops::Union{OperatorProduct,Vector{QuantumOperator}}; orders::Union{Nothing,Vector{Int}}=nothing,
    name="", factor=one(_dtype.factor), weight=zero(_dtype.weight), operator=Sum())
    @assert length(ops) == 2
    @assert adjoint(ops[1].operator) == ops[2].operator
    sign, perm = correlator_order(OperatorProduct(ops))
    if isnothing(orders)
        return FeynmanGraph(FeynmanGraph[]; topology=[[1, 2]], external_indices=perm, external_legs=[true, true], vertices=OperatorProduct.(ops),
            diagtype=Propagator(), name=name, operator=operator, factor=factor * sign, weight=weight)
    else
        return FeynmanGraph(FeynmanGraph[]; topology=[[1, 2]], external_indices=perm, external_legs=[true, true], vertices=OperatorProduct.(ops),
            orders=orders, diagtype=Propagator(), name=name, operator=operator, factor=factor * sign, weight=weight)
    end
end

"""
    function interaction(ops::OperatorProduct; name="", reorder::Union{Function,Nothing}=nothing,
        factor=one(_dtype.factor), weight=zero(_dtype.weight), operator=Sum())

Create a Interaction-type FeynmanGraph from given OperatorProduct `ops`, including several quantum operators for a vertex.
One can call a reorder function for the operators ordering.  
"""
function interaction(ops::OperatorProduct; name="", reorder::Union{Function,Nothing}=nothing,
    factor=one(_dtype.factor), weight=zero(_dtype.weight), operator=Sum())
    @assert !isfermionic(ops) "interaction OperatorProduct must be bosonic."
    if !isnothing(reorder)
        sign, perm = reorder(ops)
        return FeynmanGraph(FeynmanGraph[]; external_indices=perm, external_legs=[false for i in eachindex(perm)],
            vertices=[OperatorProduct(ops)], diagtype=Interaction(), name=name, operator=operator, factor=factor * sign, weight=weight)
    end
    _external_indices = collect(eachindex(ops))
    return FeynmanGraph(FeynmanGraph[]; external_indices=_external_indices, external_legs=[false for i in eachindex(_external_indices)],
        vertices=[ops], diagtype=Interaction(), name=name, operator=operator, factor=factor, weight=weight)
end

"""
    function external_vertex(ops::OperatorProduct;
        name="", factor=one(_dtype.factor), weight=zero(_dtype.weight), operator=Sum())

Create a ExternalVertex-type FeynmanGraph from given OperatorProduct `ops`, including several quantum operators for an purely external vertex.
"""
function external_vertex(ops::OperatorProduct;
    name="", factor=one(_dtype.factor), weight=zero(_dtype.weight), operator=Sum())
    external_indices = collect(eachindex(ops))
    return FeynmanGraph(FeynmanGraph[]; external_indices=external_indices, external_legs=[false for i in external_indices],
        vertices=[ops], diagtype=ExternalVertex(), name=name, operator=operator, factor=factor, weight=weight)
end

"""
    function group(gv::AbstractVector{G}, indices::Vector{Int}) where {G<:FeynmanGraph}

Group the graphs in `gv` by the external operators at the indices `indices`. Return a dictionary of `Vector{OperatorProduct}` to `GraphVector`.

# Example

```julia-repl
julia> p1 = propagator(𝑓⁺(1)𝑓⁻(2));

julia> p2 = propagator(𝑓⁺(1)𝑓⁻(3));

julia> p3 = propagator(𝑓⁺(2)𝑓⁻(3));

julia> gv = [p1, p2, p3];

julia> ComputationalGraphs.group(gv, [1, 2])
Dict{Vector{OperatorProduct}, Vector{FeynmanGraph{Float64, Float64}}} with 3 entries:
  [f⁻(2), f⁺(1)] => [1:f⁺(1)|f⁻(2)⋅-1.0=0.0]
  [f⁻(3), f⁺(1)] => [2:f⁺(1)|f⁻(3)⋅-1.0=0.0]
  [f⁻(3), f⁺(2)] => [3:f⁺(2)|f⁻(3)⋅-1.0=0.0]

julia> ComputationalGraphs.group(gv, [1, ])
Dict{Vector{OperatorProduct}, Vector{FeynmanGraph{Float64, Float64}}} with 2 entries:
  [f⁻(3)] => [2:f⁺(1)|f⁻(3)⋅-1.0=0.0, 3:f⁺(2)|f⁻(3)⋅-1.0=0.0]
  [f⁻(2)] => [1:f⁺(1)|f⁻(2)⋅-1.0=0.0]

julia> ComputationalGraphs.group(gv, [2, ])
Dict{Vector{OperatorProduct}, Vector{FeynmanGraph{Float64, Float64}}} with 2 entries:
  [f⁺(2)] => [3:f⁺(2)|f⁻(3)⋅-1.0=0.0]
  [f⁺(1)] => [1:f⁺(1)|f⁻(2)⋅-1.0=0.0, 2:f⁺(1)|f⁻(3)⋅-1.0=0.0]
```
"""
function group(gv::AbstractVector{G}, indices::Vector{Int}) where {G<:FeynmanGraph}
    l = length(external_indices(gv[1]))
    @assert all(x -> length(external_indices(x)) == l, gv)
    groups = Dict{Vector{OperatorProduct},Vector{G}}()
    for t in gv
        ext = external_operators(t)
        key = [OperatorProduct(ext[i]) for i in indices]
        if haskey(groups, key)
            push!(groups[key], t)
        else
            groups[key] = [t,]
        end
    end
    return groups
end
