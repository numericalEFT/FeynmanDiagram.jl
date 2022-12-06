abstract type Operator end
struct Sum <: Operator end
struct Prod <: Operator end
Base.isequal(a::Operator, b::Operator) = (typeof(a) == typeof(b))
Base.:(==)(a::Operator, b::Operator) = Base.isequal(a, b)
apply(o::Operator, diags) = error("not implemented!")

Base.show(io::IO, o::Operator) = print(io, typeof(o))
Base.show(io::IO, o::Sum) = print(io, "⨁")
Base.show(io::IO, o::Prod) = print(io, "Ⓧ")

# function vstr(r, c)
#     N = length(r)
#     # cstr(x) = x ? "⁺" : "⁻"
#     s = ""
#     for i = 1:N-1
#         s *= "$(r[i])$c"
#     end
#     s *= "$(r[end])$c"
#     return s
# end

# function vcstr(r, creation)
#     N = length(r)
#     # cstr(x) = x ? "⁺" : "⁻"
#     s = ""
#     for i = 1:N-1
#         if creation[i]
#             s *= "$(r[i])⁺"
#         else
#             s *= "$(r[i])⁻"
#         end
#     end
#     if creation[end]
#         s *= "$(r[end])⁺"
#     else
#         s *= "$(r[end])⁻"
#     end
#     return s
# end

# @enum Reducibility begin
#     OneFermiIrreducible
#     OneBoseIrreducible
#     ParticleHoleIrreducible
#     ParticleParticleIrreducible
# end

"""Type alias for a directed graph edge e = (a₁⁺, a₂⁻) from e[1] to e[2]."""
const EdgeType = Tuple{QuantumOperator,QuantumOperator}

"""
    mutable struct Graph{F,W}
    
    mutable struct of a Feynman diagram. A diagram of a sum or produce of various subgraphs.

# Members:
- `id::Int`  the unique hash id to identify the diagram
- `name::Symbol`  name of the diagram
- `type::Symbol`  type of the diagram, support :propagator, :interaction, :sigma, :green, :generic
- `orders::Vector{Int}`  orders of the diagram, e.g. loop order, derivative order, etc.
- `external::Vector{Int}`  index of external vertices
- `vertices::Vector{OperatorProduct}`  vertices of the diagram. Each index is composited by the product of quantum operators.
- `subgraph::Vector{Graph{F,W}}`  vector of sub-diagrams 
- `operator::Operator`  node operation, support Sum() and Prod()
- `factor::F`  additional factor of the diagram
- `weight::W`  weight of the diagram

# Example:
```julia-repl
julia> g = Graph([𝑓⁺(1)𝑓⁻(2), 𝑓⁺(3)𝑓⁻(4)], external=[1, 2], subgraph=[Graph([𝑓⁺(1)𝑓⁻(4)], []), Graph([𝑓⁻(2)𝑓⁺(3)], [])])
3: generic graph from f⁺(1)f⁻(2)|f⁺(3)f⁻(4)

julia> g.subgraph
2-element Vector{Graph{Float64, Float64}}:
 1: generic graph from f⁺(1)f⁻(4)
 2: generic graph from f⁻(2)f⁺(3)
```
"""
mutable struct Graph{F,W} # Graph
    id::Int
    name::String # "" by default
    type::Symbol # :propagator, :interaction, :sigma, :green, :generic
    orders::Vector{Int}

    external::Vector{Int} # index of external vertices
    vertices::Vector{OperatorProduct} # vertices of the diagram

    subgraph::Vector{Graph{F,W}}

    operator::Operator
    factor::F
    weight::W

    """
        function Graph(vertices::Vector{OperatorProduct}; external=[], subgraph=[],
            name="", type=:generic, operator::Operator=Sum(), orders=zeros(Int, 16),
            ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype))
        
        Create a Graph struct from vertices and external indices.

    # Arguments:
    - `vertices::Vector{OperatorProduct}`  vertices of the diagram
    - `external`  index of external vertices
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
    function Graph(vertices::Vector{OperatorProduct}; external=[], subgraph=[],
        name="", type=:generic, operator::Operator=Sum(), orders=zeros(Int, 16),
        ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype)
    )
        return new{ftype,wtype}(uid(), name, type, orders, external, vertices, subgraph, operator, factor, weight)
    end

    """
        function Graph(extV::AbstractVector, intV::AbstractVector; subgraph=[],
            name="", type=:generic, operator::Operator=Sum(), orders=zeros(Int, 16),
            ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype))
        
        Create a Graph struct from external and internal vertices.

    # Arguments:
    - `extV::AbstractVector`  external vertices of the diagram
    - `intV::AbstractVector`  internal vertices of the diagram
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
    function Graph(extV::AbstractVector, intV::AbstractVector; subgraph=[],
        name="", type=:generic, operator::Operator=Sum(), orders=zeros(Int, 16),
        ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype)
    )
        vertices = [extV; intV]
        ext = collect(1:length(extV))
        return new{ftype,wtype}(uid(), name, type, orders, ext, vertices, subgraph, operator, factor, weight)
    end
end

function _ops_to_str(ops::Vector{OperatorProduct})
    strs = ["$(o)" for o in ops]
    return join(strs, "|")
end

"""
    show(io::IO, g::Graph)

    Write a text representation of `Graph` to the output stream `io`.
"""
function Base.show(io::IO, g::Graph)
    #TODO: improve a text representation of Graph to the output stream.
    if isempty(g.name)
        print(io, "$(g.id): $(g.type) graph from $(_ops_to_str(g.vertices))")
    else
        print(io, "$(g.id), $(g.name): $(g.type) graph from $(_ops_to_str(g.vertices))")
    end
end
Base.show(io::IO, ::MIME"text/plain", g::Graph) = Base.show(io, g)

function Base.isequal(a::Graph, b::Graph)
    typeof(a) != typeof(b) && return false
    for field in fieldnames(typeof(a))
        if field == :weight
            (getproperty(a, :weight) ≈ getproperty(b, :weight)) == false && return false
        end
        getproperty(a, field) != getproperty(b, field) && return false
    end
    return true
end
Base.:(==)(a::Graph, b::Graph) = Base.isequal(a, b)
# isbare(diag::Graph) = isempty(diag.subgraph)

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

# function Base.:*(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
#     type = :generic
#     for v1 in g1.internal_vertices
#         for v2 in g2.internal_vertices
#             @assert v1.point != v2.point "g1 and g2 have the same internal vertex point."
#         end
#     end
#     couplings = union(g1.couplings, g2.couplings)
#     extV, intV = _getVertices(g1.external_vertices, g2.external_vertices, couplings)
#     return Graph{F,W}(extV, [g1.internal_vertices; g2.internal_vertices; intV]; type=type, couplings=couplings, subgraph=[g1, g2], operator=Prod())
# end

function Base.:*(g1::Graph{F,W}, c2::C) where {F,W,C}
    # factor = F(1) * c2
    # ftype = typeof(factor)
    return Graph(g1.vertices; external=g1.external, type=g1.type, subgraph=[g1,], operator=Prod(), ftype=F, wtype=W, factor=F(c2))
end

function Base.:*(c1::C, g2::Graph{F,W}) where {F,W,C}
    return Graph(g2.vertices; external=g2.external, type=g2.type, subgraph=[g2,], operator=Prod(), ftype=F, wtype=W, factor=F(c1))
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

    g = Graph(vertices; external=external, name=name, type=type, operator=Prod(), factor=factor, weight=weight)
    edges, contraction_sign = contractions_to_edges(vertices, contractions)
    g.factor *= contraction_sign
    for edge in edges
        push!(g.subgraph, propagator(edge[1] * edge[2]))
    end

    return g
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

    return edges, sign
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

# function checkVertices(g::Graph{F,W}) where {F,W}
#     @assert !isempty(g.couplings) "Graph.couplings must be defined to check the legality of vertices"
#     for v in g.internal_vertices
#         @assert v.operator in g.couplings "internal vertex point $(v.point) is illegal."
#     end
#     for v in g.external_vertices
#         is_legal = false
#         for ops in g.couplings
#             @assert v.operator != ops "external vertex point $(v.point) is illegal."
#             num_pops, num_ops = _countervector(v.operator.operators), _countervector(ops.operators)
#             all([num_ops[op] >= num_pops[op] for op in v.operator.operators]) && (is_legal = true)
#         end
#         @assert is_legal "external vertex point $(v.point) is illegal."
#     end
#     return true
# end

# """
#     function _getVertices(V1::Vector{ExternalVertex}, V2::Vector{ExternalVertex}, couplings::Couplings)

#     Give internal and external vertices when two diagrams combine.

#     # Arguments
#     - V1::Vector{ExternalVertex}            : External vertices of diagram I.
#     - V2::Vector{ExternalVertex}            : External vertices of diagram II.
#     - couplings::Vector{OperatorProduct}  : all the vertex couplings in the Feynman rule. 
# """
# function _getVertices(V1::Vector{ExternalVertex}, V2::Vector{ExternalVertex}, couplings::Vector{OperatorProduct})
#     V1_ind = [v.point for v in V1]
#     V2_ind = [v.point for v in V2]
#     common = intersect(V1_ind, V2_ind)
#     total = union(V1, V2)
#     intV = []
#     extV = [v for v in total if v.point ∉ common]
#     for point in common
#         # ifinternal, ifexternal = false, false
#         i1 = findfirst(isequal(point), V1_ind)
#         i2 = findfirst(isequal(point), V2_ind)
#         point_ops = OperatorProduct([V1[i1].operator.operators; V2[i2].operator.operators])
#         if isempty(couplings)
#             append!(extV, [ExternalVertex(point, V1[i1].current, point_ops)])
#         else
#             for ops in couplings
#                 if point_ops == ops
#                     append!(intV, [InternalVertex(point, V1[i1].current, point_ops)])
#                 else
#                     append!(extV, [ExternalVertex(point, V1[i1].current, point_ops)])
#                 end
#             end
#         end
#         # @assert ifinternal || ifexternal "point $point is illegal."
#     end
#     return extV, intV
# end

# function 𝐺ᶠ(point_in::Int, point_out::Int, current::Int=0; kwargs...)
#     return Green2(point_in, point_out, current; isFermi=true, kwargs...)
# end

# function 𝐺ᵇ(point_in::Int, point_out::Int, current::Int=0; kwargs...)
#     return Green2(point_in, point_out, current; isFermi=false, kwargs...)
# end

# function 𝐺ᵠ(point_in::Int, point_out::Int, current::Int=0; kwargs...)
#     return Green2(point_in, point_out, current; isFermi=false, isComplex=false, kwargs...)
# end

# function Green2(point_in::Int, point_out::Int, current::Int=0;
#     isFermi=true, isComplex=true, flavor::Int=1, couplings=[], dtype=Float64,
#     factor=zero(dtype), weight=zero(dtype), name="G2", subgraph=[], operator=Sum())
#     if isFermi && isComplex
#         opin, opout = fermionic_creation(flavor), fermionic_annihilation(flavor)
#     elseif isFermi && !isComplex
#         opin, opout = majorana(flavor), majorana(flavor)
#     elseif !isFermi && isComplex
#         opin, opout = bosonic_creation(flavor), bosonic_annihilation(flavor)
#     else
#         opin, opout = real_scalar(flavor), real_scalar(flavor)
#     end
#     if point_in == point_out
#         extV = [ExternalVertex(point_in, current, OperatorProduct([opin, opout]))]
#     else
#         extV = [ExternalVertex(point_in, current, OperatorProduct([opin])),
#             ExternalVertex(point_out, current, OperatorProduct([opout]))]
#     end
#     if isnothing(subgraph)
#         diagtype = :propagator
#     else
#         diagtype = :green
#     end
#     return Graph{dtype,dtype}(extV, type=diagtype, couplings=couplings, subgraph=subgraph,
#         name=name, operator=operator, factor=factor, weight=weight)
# end

# const 𝑊 = Interaction

# function Interaction(point_in::Int, point_out::Int, current::Int=0; flavor::Int=1,
#     couplings=[Coupling_yukawa,], dtype=Float64, factor=zero(dtype), weight=zero(dtype), name="W")
#     ext_in = ExternalVertex(point_in, current, OperatorProduct([real_scalar(flavor),]))
#     ext_out = ExternalVertex(point_out, current, OperatorProduct([real_scalar(flavor),]))
#     diagtype = :interaction2
#     return Graph{dtype,dtype}([ext_in, ext_out], type=diagtype, couplings=couplings,
#         name=name, factor=factor, weight=weight)
# end

# function Interaction(point::Int, coupling::OperatorProduct, current::Int=0;
#     dtype=Float64, factor=zero(dtype), weight=zero(dtype), name="W", operators=[])
#     extV = ExternalVertex(point, current, OperatorProduct(operators))
#     # diagtype = Symbol("bareVertex$(length(operator))"...)
#     diagtype = :interaction
#     return Graph{dtype,dtype}([extV,], type=diagtype, couplings=[coupling,], name=name, factor=factor, weight=weight)
# end