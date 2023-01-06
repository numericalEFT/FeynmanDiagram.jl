abstract type AbstractOperator end
struct Sum <: AbstractOperator end
struct Prod <: AbstractOperator end
Base.isequal(a::AbstractOperator, b::AbstractOperator) = (typeof(a) == typeof(b))
Base.:(==)(a::AbstractOperator, b::AbstractOperator) = Base.isequal(a, b)
apply(o::AbstractOperator, diags) = error("not implemented!")

Base.show(io::IO, o::AbstractOperator) = print(io, typeof(o))
Base.show(io::IO, ::Type{Sum}) = print(io, "⨁")
Base.show(io::IO, ::Type{Prod}) = print(io, "Ⓧ")

# Is the unary operation trivial (𝓞g = g)?
unary_istrivial(::Type{AbstractOperator}) = false
unary_istrivial(::Type{O}) where {O<:Union{Sum,Prod}} = true  # (+g) ≡ g and (*g) ≡ g

# Is the operation associative: a 𝓞 (b 𝓞 c) = (a 𝓞 b) 𝓞 c = a 𝓞 b 𝓞 c?
isassociative(::Type{AbstractOperator}) = false
isassociative(::Type{Sum}) = true
# NOTE: Associativity of Prod (graph composition)
#       requires Base.*(g1, g2) and Base./(g1, g2)
# isassociative(::Type{Prod}) = true

abstract type GraphType end
struct Interaction <: GraphType end
struct ExternalVertex <: GraphType end
struct Propagator <: GraphType end
struct SelfEnergy <: GraphType end
struct VertexDiag <: GraphType end
struct GreenDiag <: GraphType end
struct GenericDiag <: GraphType end

"""
    mutable struct Graph{F,W}
    
    Computational Graph representation of a collection of Feynman diagrams. All Feynman diagrams should share the same set of external and internal vertices.

# Members:
- `id::Int`  the unique hash id to identify the diagram
- `name::Symbol`  name of the diagram
- `type::DataType`  type of the diagram, support Interaction, ExternalVertex, Propagator, SelfEnergy, VertexDiag, GreenDiag, and GenericDiag.
- `orders::Vector{Int}`  orders of the diagram, e.g. loop order, derivative order, etc.
- `vertices::Vector{OperatorProduct}`  vertices of the diagram. Each index is composited by the product of quantum operators. 
- `topology::Vector{Vector{Int}}` topology of the diagram. Each Vector{Int} stores vertices' index connected with each other (as a propagator). 
- `external::Vector{Int}`  index of ACTUAL external vertices (as QuantumOperators)
- `hasLeg::Vector{Bool}` index of each external operator (true: legged, false: nonleg)
- `subgraphs::Vector{Graph{F,W}}`  vector of sub-diagrams 
- `subgraph_factors::Vector{F}`  scalar multiplicative factors associated with each subdiagram
- `operator::DataType`  node operation, support Sum and Prod
- `factor::F`  total scalar multiplicative factor for the diagram
- `weight::W`  weight of the diagram

# Example:
```julia-repl
julia> g1 = Graph([], vertices=[𝑓⁺(1),𝑓⁻(2)], external=[1,2], hasLeg=[true,true])
1:f⁺(1)|f⁻(2)=0.0

julia> g2 = Graph([], vertices=[𝑓⁺(3),𝑓⁻(4)], external=[1,2], hasLeg=[true,true])
2:f⁺(3)|f⁻(4)=0.0

julia> g = Graph([g1,g2], vertices=[𝑓⁺(1),𝑓⁻(2),𝑓⁺(3),𝑓⁻(4)], operator=ComputationalGraphs.Prod(), external=[1,2,3,4], hasLeg=[true,true,true,true])
3:f⁺(1)|f⁻(2)|f⁺(3)|f⁻(4)=0.0=Ⓧ (1,2)
```
"""
mutable struct Graph{F,W} # Graph
    id::Int
    name::String # "" by default
    # type::Symbol # :propagator, :interaction, :sigma, :green, :generic
    type::DataType
    orders::Vector{Int}

    vertices::Vector{OperatorProduct} # vertices of the diagram
    topology::Vector{Vector{Int}}
    external::Vector{Int} # index of external operators
    hasLeg::Vector{Bool} # Bool indexes for all external operators (true: real leg, false: fake leg)

    subgraphs::Vector{Graph{F,W}}
    subgraph_factors::Vector{F}

    operator::DataType
    factor::F
    weight::W

    """
        function Graph(vertices::Vector{OperatorProduct}; external=[], subgraphs=[],
            name="", type=:generic, operator::AbstractOperator=Sum(), orders=zeros(Int, 16),
            ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype))
        
        Create a Graph struct from vertices and external indices.

    # Arguments:
    - `subgraphs`  vector of sub-diagrams 
    - `topology` topology of the diagram
    - `vertices::Union{Vector{OperatorProduct},Nothing}`  vertices of the diagram, nothing by default
    - `external`  index of actual external vertices in terms of QuantumOperators, empty by default
    - `subgraph_factors::Vector{F}`  scalar multiplicative factors associated with each subdiagram
    - `name`  name of the diagram
    - `type::GraphType`  type of the diagram
    - `operator::DataType`  node operation, Sum, Prod, etc.
    - `orders`  orders of the diagram
    - `ftype`  typeof(factor)
    - `wtype`  typeof(weight)
    - `factor::F`  overall scalar multiplicative factor for this diagram (e.g., permutation sign)
    - `weight`  weight of the diagram
    """
    function Graph(subgraphs::AbstractVector; topology=[], vertices::Union{Vector{OperatorProduct},Nothing}=nothing, external=[], hasLeg=[],
        subgraph_factors=one.(eachindex(subgraphs)), name="", type::GraphType=GenericDiag(), operator::AbstractOperator=Sum(),
        orders=zeros(Int, 16), ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype)
    )
        @assert length(external) == length(hasLeg)
        # ftype = promote_type([typeof(g.factor) for g in subgraphs]...)
        # wtype = promote_type([typeof(g.weight) for g in subgraphs]...)
        if isnothing(vertices)
            vertices = [OperatorProduct(OperatorProduct(g.vertices)[g.external]) for g in subgraphs if g.type != Propagator]
        end
        return new{ftype,wtype}(uid(), name, typeof(type), orders, vertices, topology, external,
            hasLeg, subgraphs, subgraph_factors, typeof(operator), factor, weight)
    end
    # function Graph(vertices::AbstractVector; external=[], subgraphs=[], subgraph_factors=one.(eachindex(subgraphs)),
    #     topology=[], name="", type=:generic, operator::AbstractOperator=Sum(), orders=zeros(Int, 16),
    #     ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype)
    # )
    #     vertices = [OperatorProduct(v) for v in vertices]
    #     return new{ftype,wtype}(uid(), name, type, orders, external, vertices, topology,
    #         subgraphs, subgraph_factors, typeof(operator), factor, weight)
    # end
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
# isbare(g::Graph) = isempty(g.subgraphs)

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
        elseif field == :subgraphs
            length(a.subgraphs) != length(b.subgraphs) && return false
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

# """
#     function isghost(g::Graph, i::Int) 

#     Check if `i::Int` in the ghost operator's indices of Graph `g`.
# """
# isghost(g::Graph, i::Int) = isghost(OperatorProduct(g.vertices)[i])

"""
    function vertices(g::Graph)

    Return all vertices (::Vector{OperatorProduct}) of Graph `g`.
"""
vertices(g::Graph) = g.vertices

"""
    function external(g::Graph)

    Return all physical external vertices (::OperatorProduct}) of Graph `g`.
"""
external(g::Graph) = OperatorProduct(OperatorProduct(g.vertices)[g.external])
# external(g::Graph) = OperatorProduct.(OperatorProduct(g.vertices)[g.external])

"""
    function external_labels(g::Graph)

    Return the labels of all physical external vertices of Graph `g`.
"""
external_labels(g::Graph) = [o.label for o in external(g)]

# """
#     function external_with_ghost(g::Graph)

#     Return all the external vertices (::Vector{OperatorProduct}), including real legs and ghost legs.
# """
# external_with_ghost(g::Graph) = OperatorProduct.(OperatorProduct(g.vertices)[eachindex(g.external)])

# """
#     function external_with_ghost_labels(g::Graph)

#     Return the labels of all external vertices, including both real legs and ghost legs.
# """
# external_with_ghost_labels(g::Graph) = [o[1].label for o in external_with_ghost(g)]

#TODO: add function return reducibility of Graph. 
function reducibility(g::Graph)
    return (OneFermiIrreducible,)
end

#TODO: add function for connected diagram check. 
function connectivity(g::Graph)
    isempty(g.subgraphs) && return true
end

function Base.:*(g1::Graph{F,W}, c2::C) where {F,W,C}
    g = Graph([g1,]; topology=g1.topology, external=g1.external, hasLeg=g1.hasLeg, vertices=g1.vertices,
        type=g1.type(), subgraph_factors=[F(c2),], operator=Prod(), ftype=F, wtype=W)
    # Merge multiplicative link
    if g1.operator == Prod && length(g1.subgraph_factors) == 1
        g.subgraph_factors[1] *= g1.subgraph_factors[1]
        g.subgraphs = g1.subgraphs
    end
    return g
end

function Base.:*(c1::C, g2::Graph{F,W}) where {F,W,C}
    g = Graph([g2,]; topology=g2.topology, external=g2.external, hasLeg=g2.hasLeg, vertices=g2.vertices,
        type=g2.type(), subgraph_factors=[F(c1),], operator=Prod(), ftype=F, wtype=W)
    # Merge multiplicative link
    if g2.operator == Prod && length(g2.subgraph_factors) == 1
        g.subgraph_factors[1] *= g2.subgraph_factors[1]
        g.subgraphs = g2.subgraphs
    end
    return g
end

"""Returns a graph representing the linear combination `c1*g1 + c2*g2`."""
function linear_combination(g1::Graph{F,W}, g2::Graph{F,W}, c1::C, c2::C) where {F,W,C}
    # TODO: more check
    @assert g1.type == g2.type "g1 and g2 are not of the same type."
    @assert g1.orders == g2.orders "g1 and g2 have different orders."
    @assert Set(vertices(g1)) == Set(vertices(g2)) "g1 and g2 have different vertices."
    @assert Set(external(g1)) == Set(external(g2)) "g1 and g2 have different external vertices."
    return Graph([g1, g2]; external=g1.external, hasLeg=g1.hasLeg, vertices=g1.vertices, type=g1.type(),
        subgraph_factors=[F(c1), F(c2)], operator=Sum(), ftype=F, wtype=W)
end

"""
Given a vector `graphs` of graphs each with the same type and external/internal
vertices and an equally-sized vector `constants` of constants, returns a new
graph representing the linear combination ⟨`graphs`, `constants`⟩.
"""
function linear_combination(graphs::Vector{Graph{F,W}}, constants::Vector{C}) where {F,W,C}
    # TODO: more check
    @assert alleq(getproperty.(graphs, :type)) "Graphs are not all of the same type."
    @assert alleq(getproperty.(graphs, :orders)) "Graphs do not all have the same order."
    @assert alleq(Set.(vertices.(graphs))) "Graphs do not share the same set of vertices."
    @assert alleq(Set.(external.(graphs))) "Graphs do not share the same set of external vertices."
    g1 = graphs[1]
    return Graph(graphs; external=g1.external, hasLeg=g1.hasLeg, vertices=g1.vertices, type=g1.type(),
        subgraph_factors=constants, operator=Sum(), ftype=F, wtype=W)
end

function Base.:+(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
    return linear_combination(g1, g2, F(1), F(1))
end

function Base.:-(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
    return linear_combination(g1, g2, F(1), F(-1))
end

# function Base.:+(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
#     return linear_combination([g1, g2], [F(1), F(1)])
# end

# function Base.:-(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
#     return linear_combination([g1, g2], [F(1), F(-1)])
# end

"""
    function feynman_diagram(subgraphs::Vector{Graph{F,W}}, topology::Vector{Vector{Int}}, perm_noleg::Union{Vector{Int},Nothing}=nothing;
        factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", diagtype::GraphType=GenericDiag()) where {F,W}
    
    Create a Graph representing feynman diagram from all subgraphs and topology (connections between vertices),
    where each ExternalVertex is given in `vertices`, 
    while internal vertices are constructed with external legs of graphs in `vertices`, or simply OperatorProduct in `vertices`.
    
# Arguments:
- `subgraphs::Vector{Graph{F,W}}` all subgraphs of the diagram. All external operators of subgraphs constitute all operators of the new diagram.
- `topology::Vector{Vector{Int}}` topology of the diagram. Each Vector{Int} stores operators' index connected with each other (as a propagator). 
- `perm_noleg::Union{Vector{Int},Nothing}=nothing` permutation of all the nonleg external operators. By default, setting nothing means to use the default order from subgraphs.
- `factor::F`  overall scalar multiplicative factor for this diagram (e.g., permutation sign)
- `weight`  weight of the diagram
- `name`  name of the diagram
- `type`  type of the diagram

# Example:
```julia-repl
julia> V = [𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6), 𝑓⁺(7)𝑓⁻(8)𝜙(9)];
julia> g = feynman_diagram(interaction.(V), [[1, 5], [3, 9], [4, 8]], perm_noleg=[3, 1, 2])
7:f⁺(1)f⁻(2)ϕ(3)|f⁺(4)f⁻(5)ϕ(6)|f⁺(7)f⁻(8)ϕ(9)=0.0=Ⓧ (1,2,3,4,5,6)

julia> g.subgraphs
6-element Vector{Graph{Float64, Float64}}:
 1:f⁺(1)f⁻(2)ϕ(3)=0.0
 2:f⁺(4)f⁻(5)ϕ(6)=0.0
 3:f⁺(7)f⁻(8)ϕ(9)=0.0
 4:f⁺(1)|f⁻(5)⋅-1.0=0.0
 5:ϕ(3)|ϕ(9)=0.0
 6:f⁺(4)|f⁻(8)⋅-1.0=0.0
```
"""
function feynman_diagram(subgraphs::Vector{Graph{F,W}}, topology::Vector{Vector{Int}}, perm_noleg::Union{Vector{Int},Nothing}=nothing;
    factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", diagtype::GraphType=GenericDiag()) where {F,W}

    # external_ops = OperatorProduct(operators[external]) # the external operators for the building diagram after contractions
    contraction = collect(Iterators.flatten(topology))
    @assert length(unique(contraction)) == length(contraction)  # no repeated index

    vertices, all_hasLeg = OperatorProduct[], Bool[]
    external_leg, external_noleg = Int[], Int[] # index all leg/nonleg external operators
    ind = 0
    for g in subgraphs
        g.type == Propagator && continue  # exclude propagator subgraph to avoid double counting.
        push!(vertices, external(g))
        append!(all_hasLeg, g.hasLeg)
        if g.type == ExternalVertex
            append!(external_leg, g.external .+ ind) # ExternalVertex will be legged after contraction.
        else
            gext = setdiff(g.external .+ ind, contraction) # select all external operators
            gextLeg = g.hasLeg[gext.-ind]
            # the selected gext[i] with gextLeg[i]==true is the external vertice with a leg
            append!(external_leg, gext[gextLeg])
            append!(external_noleg, gext[gextLeg.==false])
        end
        ind += length(g.external)
    end

    @assert !any(all_hasLeg[setdiff(eachindex(all_hasLeg), external_noleg)]) "all contracted operators should have no leg."
    @assert external_leg ⊆ contraction
    @assert isempty(intersect(contraction, external_noleg)) "all nonleg external operators should not be contracted"
    if !isnothing(perm_noleg)
        @assert length(unique(perm_noleg)) == length(perm_noleg) == length(external_noleg)
        external_noleg = external_noleg[perm_noleg]
    end

    operators = OperatorProduct(vertices) # all external operators from subgraphs
    permutation = union(contraction, external_noleg)
    @assert Set(permutation) == Set(eachindex(operators)) # permutation must exhaust all operators

    fermionic_operators = isfermionic.(operators)
    filter!(p -> fermionic_operators[p], permutation)
    sign = isempty(permutation) ? 1 : parity(sortperm(permutation))

    for connection in topology
        push!(subgraphs, propagator(operators[connection]))
    end
    _external = union(external_leg, external_noleg)
    _hasLeg = append!([true for i in eachindex(external_leg)], [false for i in eachindex(external_noleg)])
    return Graph(subgraphs; topology=topology, external=_external, hasLeg=_hasLeg, vertices=vertices,
        name=name, type=diagtype, operator=Prod(), factor=factor * sign, weight=weight)
end

# """
#     function feynman_diagram(vertices::AbstractVector, topology::Vector{Vector{Int}};
#         external=[], factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", type=:generic)

#     Create a Graph representing feynman diagram from all vertices and topology (connections between vertices),
#     where external vertices are given at the beginning of `vertices`, 
#     while internal vertices are constructed with external legs of graphs in `vertices`, or simply OperatorProduct in `vertices`.

# # Arguments:
# - `vertices::AbstractVector` all vertices of the diagram. Formal external operators at the beginning, internal vertices given from Graph or OperatorProduct.
# - `topology::Vector{Vector{Int}}` topology of the diagram. Each Vector{Int} stores vertices' index connected with each other (as a propagator). 
# - `external`  index of external vertices. They are the actual external quantum operators, not the ghost operators.
# - `factor::F`  overall scalar multiplicative factor for this diagram (e.g., permutation sign)
# - `weight`  weight of the diagram
# - `name`  name of the diagram
# - `type`  type of the diagram

# # Example:
# ```julia-repl
# julia> g = feynman_diagram([𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6)], [[5, 1], [2, 4], [3, 6]])
# 4:f⁺(1)f⁻(2)ϕ(3)|f⁺(4)f⁻(5)ϕ(6)=0.0=-1.0Ⓧ (1,2,3)

# julia> g.subgraphs
# 3-element Vector{Graph{Float64, Float64}}:
# 1:f⁻(5)|f⁺(1)=0.0
# 2:f⁻(2)|f⁺(4)=0.0
# 3:ϕ(3)|ϕ(6)=0.0
# ```
# """
# function feynman_diagram(vertices::AbstractVector,
#     topology::Vector{Vector{Int}};
#     external=[], factor=one(_dtype.factor),
#     weight=zero(_dtype.weight),
#     name="", type=:generic)

#     # convert all Graphs to vertices by grouping all its external legs as a vertex
#     converted_vertices = [_extract_vertex(typeof(v), v) for v in vertices]
#     # collect graphs in vertices as subdiagrams
#     subdiagrams = [v for v in vertices if isa(v, Graph)]

#     operators = [o for v in converted_vertices for o in v.operators]
#     permutation = collect(Iterators.flatten(topology))
#     ind_ops = collect(eachindex(operators))


#     @assert length(unique(permutation)) == length(permutation) # no repeated index
#     @assert length(unique(external)) == length(external) # no repeated index
#     # @assert Set(permutation) == Set(ind_ops) # permutation must exhaust all operators
#     ind_ghost = filter(p -> isghost(operators[p]), ind_ops)
#     @assert all(ind_ghost .<= length(external)) # external real/fake legs must be placed at the beginning of vertices.

#     ind_fakeleg = Int[]
#     subgraphs = Graph[]
#     for connection in topology
#         if isempty(intersect(connection, ind_ghost))
#             push!(subgraphs, propagator(operators[connection]))
#         else
#             @assert length(connection) == 2 "Ghost external operator can only be connected to a single internal operator"
#             ind_fop = setdiff(connection, ind_ghost)
#             append!(ind_fakeleg, ind_fop)
#         end
#     end
#     @assert ind_fakeleg ⊆ external "external operators are not consistent with ghost operators in vertices. Ghost leg indices: $ind_fakeleg not in external: $external"

#     fermionic_operators = isfermionic.(operators)
#     filter!(p -> fermionic_operators[p], permutation)
#     sign = isempty(permutation) ? 1 : parity(sortperm(permutation))

#     append!(subgraphs, subdiagrams)

#     g = Graph(converted_vertices; external=external, subgraphs=subgraphs, topology=topology, name=name,
#         type=type, operator=Prod(), factor=factor * sign, weight=weight)
#     return g
# end

# do nothing when already a OperatorProduct; 
_extract_vertex(::Type{<:OperatorProduct}, g) = g
# helper functions extracting external legs from g::Graph to form a vertex 
_extract_vertex(::Type{<:Graph}, g) = OperatorProduct(external(g))

"""
    function propagator(ops::Union{OperatorProduct,Vector{QuantumOperator}};
        name="", factor=one(_dtype.factor), weight=zero(_dtype.weight), operator=Sum())

    Create a Propagator-type Graph from given OperatorProduct or Vector{QuantumOperator} `ops`, including two quantum operators.
"""
function propagator(ops::Union{OperatorProduct,Vector{QuantumOperator}};
    name="", factor=one(_dtype.factor), weight=zero(_dtype.weight), operator=Sum())
    @assert length(ops) == 2
    @assert adjoint(ops[1].operator) == ops[2].operator
    sign, perm = correlator_order(OperatorProduct(ops))
    return Graph(Graph[]; topology=[[1, 2]], external=perm, hasLeg=[true, true], vertices=OperatorProduct.(ops),
        type=Propagator(), name=name, operator=operator, factor=factor * sign, weight=weight)
end

"""
    function interaction(ops::OperatorProduct; name="", reorder::Union{Function,Nothing}=nothing,
        factor=one(_dtype.factor), weight=zero(_dtype.weight), operator=Sum())
    Create a Interaction-type Graph from given OperatorProduct `ops`, including several quantum operators for a vertex.
    One can call a reorder function for the operators ordering.  
"""
function interaction(ops::OperatorProduct; name="", reorder::Union{Function,Nothing}=nothing,
    factor=one(_dtype.factor), weight=zero(_dtype.weight), operator=Sum())
    @assert !isfermionic(ops) "interaction OperatorProduct must be bosonic."
    if !isnothing(reorder)
        sign, perm = reorder(ops)
        return Graph(Graph[]; external=perm, hasLeg=[false for i in eachindex(perm)],
            vertices=[OperatorProduct(ops)], type=Interaction(), name=name, operator=operator, factor=factor * sign, weight=weight)
    end
    _external = collect(eachindex(ops))
    return Graph(Graph[]; external=_external, hasLeg=[false for i in eachindex(_external)],
        vertices=[ops], type=Interaction(), name=name, operator=operator, factor=factor, weight=weight)
end

"""
    function external_vertex(ops::OperatorProduct;
        name="", factor=one(_dtype.factor), weight=zero(_dtype.weight), operator=Sum())
    Create a ExternalVertex-type Graph from given OperatorProduct `ops`, including several quantum operators for an purely external vertex.
"""
function external_vertex(ops::OperatorProduct;
    name="", factor=one(_dtype.factor), weight=zero(_dtype.weight), operator=Sum())
    external = collect(eachindex(ops))
    return Graph(Graph[]; external=external, hasLeg=[false for i in external],
        vertices=[ops], type=ExternalVertex(), name=name, operator=operator, factor=factor, weight=weight)
end

# function propagator(ops::OperatorProduct;
#     name="", diagtype=:propagator, factor=one(_dtype.factor), weight=zero(_dtype.weight), operator=Sum())
#     return Graph([ops,]; external=collect(eachindex(ops)), type=diagtype, name=name, operator=operator, factor=factor, weight=weight)
# end

"""
    function standardize_order!(g::Graph)

    Standardize the external operators' order of Graph. 
    Reorder all leaves (propagators) of Graph by correlator ordering. 
    Reorder all non-leaves of Graph by normal ordering.

# Example: 
```julia-repl
julia> g = propagator([𝑓⁺(1), 𝑏⁺(2), 𝜙(3), 𝑓⁻(1), 𝑏⁻(2)])
1:f⁺(1)|b⁺(2)|ϕ(3)|f⁻(1)|b⁻(2)=0.0

julia> standardize_order!(g)

julia> g
11:f⁻(1)|b⁻(2)|ϕ(3)|b⁺(2)|f⁺(1)⋅-1.0=0.0
```
"""
# function standardize_order!(g::Graph)
#     for node in PreOrderDFS(g)
#         extL = external_with_ghost(node)
#         if isempty(node.subgraphs)
#             sign, perm = correlator_order(OperatorProduct(extL))
#             # node.external = node.external[perm]
#         else
#             sign, perm = normal_order(OperatorProduct(extL))
#             inds_real = [i for (i, op) in enumerate(extL) if !isghost(op[1])]
#             node.external = union(sortperm(perm)[inds_real], setdiff(node.external, perm))
#             for connection in node.topology
#                 for (i, ind) in enumerate(connection)
#                     ind in perm && (connection[i] = perm[ind])
#                 end
#             end
#         end
#         node.vertices[eachindex(node.external)] = node.vertices[perm]
#         node.factor *= sign
#     end
# end

# function attach_legs(graph::Graph, ops::OperatorProduct, attachInds::Vector{Int};
#     factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", type=:generic)
#     @assert length(ops) == length(attachInds)
#     subgraphs = [graph]
#     for (i, ind) in enumerate(attachInds)
#         push!(subgraphs, propagator(OperatorProduct([external(graph)[ind], ops[i]])))
#     end

#     return Graph(subgraphs,)
# end
