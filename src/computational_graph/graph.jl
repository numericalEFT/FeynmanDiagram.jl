abstract type Operator end
struct Sum <: Operator end
struct Prod <: Operator end
# struct Diff <: Operator end
# struct Integral <: Operator end
Base.isequal(a::Operator, b::Operator) = (typeof(a) == typeof(b))
Base.:(==)(a::Operator, b::Operator) = Base.isequal(a, b)
apply(o::Operator, diags) = error("not implemented!")

Base.show(io::IO, o::Operator) = print(io, typeof(o))
Base.show(io::IO, o::Sum) = print(io, "⨁")
Base.show(io::IO, o::Prod) = print(io, "Ⓧ")
# Base.show(io::IO, o::Diff) = print(io, "d")
# Base.show(io::IO, o::Integral) = print(io, "∫")

function vstr(r, c)
    N = length(r)
    # cstr(x) = x ? "⁺" : "⁻"
    s = ""
    for i = 1:N-1
        s *= "$(r[i])$c"
    end
    s *= "$(r[end])$c"
    return s
end

function vcstr(r, creation)
    N = length(r)
    # cstr(x) = x ? "⁺" : "⁻"
    s = ""
    for i = 1:N-1
        if creation[i]
            s *= "$(r[i])⁺"
        else
            s *= "$(r[i])⁻"
        end
    end
    if creation[end]
        s *= "$(r[end])⁺"
    else
        s *= "$(r[end])⁻"
    end
    return s
end

# @enum Reducibility begin
#     OneFermiIrreducible
#     OneBoseIrreducible
#     ParticleHoleIrreducible
#     ParticleParticleIrreducible
# end

struct QuantumOperator
    operator::Symbol
    label::Int
    function QuantumOperator(operator::Symbol, label::Int)
        @assert label > 0
        return new(operator, label)
    end
end
Base.isequal(a::QuantumOperator, b::QuantumOperator) = ((a.operator == b.operator) && (a.label == b.label))
Base.:(==)(a::QuantumOperator, b::QuantumOperator) = Base.isequal(a, b)

function Base.show(io::IO, o::QuantumOperator)
    print(io, "$(String(o.operator))($(o.label))")
end

Base.show(io::IO, ::MIME"text/plain", o::QuantumOperator) = Base.show(io, o)

fermionic_annihilation(i) = QuantumOperator(:f⁻, i)
fermionic_creation(i) = QuantumOperator(:f⁺, i)
majorana(i) = QuantumOperator(:f, i)
bosonic_annihilation(i) = QuantumOperator(:b⁻, i)
bosonic_creation(i) = QuantumOperator(:b⁺, i)
real_classic(i) = QuantumOperator(:phi, i)

const 𝑓⁻ = fermionic_annihilation
const 𝑓⁺ = fermionic_creation
const 𝑓 = majorana
const 𝑏⁻ = bosonic_annihilation
const 𝑏⁺ = bosonic_creation
const 𝜙 = real_classic

struct CompositeOperator <: AbstractVector{QuantumOperator}
    operators::Vector{QuantumOperator}
    function CompositeOperator(operators::Vector{QuantumOperator})
        return new(operators)
    end
    function CompositeOperator(operator::QuantumOperator)
        return new([operator,])
    end
    function CompositeOperator(operators::CompositeOperator)
        return new([operators...,])
    end
end

#TODO: make compositeoperator norm ordered when it is created

Base.eltype(::Type{CompositeOperator}) = QuantumOperator
Base.getindex(o::CompositeOperator, i::Int) = o.operators[i]
Base.setindex!(o::CompositeOperator, v::QuantumOperator, i::Int) = o.operators[i] = v
Base.length(o::CompositeOperator) = length(o.operators)
Base.size(o::CompositeOperator) = size(o.operators)

Base.show(io::IO, o::CompositeOperator) = print(io, reduce(*, ["$o" for o in o.operators]))
Base.show(io::IO, ::MIME"text/plain", o::CompositeOperator) = Base.show(io, o)

function Base.:*(o1::QuantumOperator, o2::QuantumOperator)
    return CompositeOperator([o1, o2])
end

function Base.:*(o1::CompositeOperator, o2::QuantumOperator)
    return CompositeOperator([o1.operators; o2])
end

function Base.:*(o1::QuantumOperator, o2::CompositeOperator)
    return CompositeOperator([o1; o2.operators])
end

function Base.:*(o1::CompositeOperator, o2::CompositeOperator)
    return CompositeOperator([o1.operators; o2.operators])
end

"""
    Converts a CompositeOperator to normal-ordered form in place and returns the associated statistical sign.
"""
function normal_order!(operator::CompositeOperator)
    sign = 1
    return sign
end

"""
    Computes the permutation required to convert a CompositeOperator to normal-ordered form. 
    Returns the associated statistical sign and permutation.
"""
function normal_order(operator::CompositeOperator)
    sign = 1
    permutation = collect(eachindex(operator.operators))
    return sign, permutation
end

"""Type alias for a directed graph edge e = (a₁⁺, a₂⁻) from e[1] to e[2]."""
const EdgeType = Tuple{Int,Int}

"""
    mutable struct Graph{F,W}
    
    struct of a Feynman diagram. A diagram of a sum or produce of various subgraphs.

# Members
- hash::Int            : the unique hash number to identify the diagram
- name::Symbol         : name of the diagram
- para::GraphPara    : internal parameters of the diagram
- orders::Vector{Int}  : orders of the diagram, loop order, derivative order, etc.
# - couplings::Couplings : all the vertex couplings in the Feynman rule. 
# - internal_points::Vector{Int} : internal points in the diagram
# - currents::Vector{Float64} : independent currents in the diagram
- external_vertices::Vector{ExternalVertex}    : external vertices of the diagram
- internal_vertices::Vector{InternalVertex}    : internal vertices of the diagram
# - isConnected::Bool    : connected or disconnected Green's function
# - isAmputated::Bool    : amputated Green's function or not
- subgraph::Vector{Graph{W}}   : vector of sub-diagrams 
- operator::Operator   : operation, support Sum() and Prod()
- factor::F            : additional factor of the diagram
- weight::W            : weight of the diagram
"""
mutable struct Graph{F,W} # Graph
    id::Int
    name::String # "" by default
    type::Symbol # :propagator, :interaction, :sigma, :green, :generic
    orders::Vector{Int}

    external::Vector{Int} # index of external vertices
    vertices::Vector{CompositeOperator} # vertices of the diagram

    subgraph::Vector{Graph{F,W}}

    operator::Operator
    factor::F
    weight::W

    function Graph(extV::Vector{CompositeOperator}, intV::Vector{CompositeOperator}; subgraph=[],
        name="", type=:generic, operator::Operator=Sum(), orders=zeros(Int, 16),
        ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype)
    )
        vertices = [extV..., intV...]
        ext = collect(1:length(extV))
        return new{ftype,wtype}(uid(), name, type, orders, ext, int, vertices, subgraph, operator, factor, weight)
    end
    function Graph(vertices::Vector{CompositeOperator}; external::Vector{Int}=[], subgraph=[],
        name="", type=:generic, operator::Operator=Sum(), orders=zeros(Int, 16),
        ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype)
    )
        return new{ftype,wtype}(uid(), name, type, orders, external, vertices, subgraph, operator, factor, weight)
    end
end

#TODO: improve a text representation of Graph to the output stream.
function Base.show(io::IO, g::Graph)
    if isempty(g.name)
        print(io, "$(g.id): $(g.type) graph from $(g.vertices)")
    else
        print(io, "$(g.id), $(g.name): $(g.type) graph from $(g.vertices)")
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

is_external(g::Graph, i::Int) = i in g.external
is_internal(g::Graph, i::Int) = (i in g.external) == false
external_vertices(g::Graph) = g.vertices[g.external]
internal_vertices(g::Graph) = g.vertices[setdiff(1:length(g.vertices), g.external)]
vertices(g::Graph) = g.vertices

#TODO: add function return reducibility of Graph. 
function reducibility(g::Graph)
    return (OneFermiIrreducible,)
end

#TODO: add function for connected diagram check. 
function connectivity(g::Graph)
    isempty(g.subgraph) && return true
end
function connectivity!(g::Graph)
    g.isConnected = connectivity(g)
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

function Base.:*(g1::Graph{F,W}, c2::C) where {F,W,C<:F}
    return Graph{F,W}(g1.external_vertices, g1.internal_vertices; type=g1.type, subgraph=[g1,], operator=Prod(), factor=c2)
end

function Base.:*(c1::C, g2::Graph{F,W}) where {F,W,C<:F}
    return Graph{F,W}(g2.external_vertices, g2.internal_vertices; type=g2.type, subgraph=[g2,], operator=Prod(), factor=c1)
end

function Base.:+(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
    @assert g1.type == g2.type "g1 and g2 are not of the same type."
    # TODO: more check
    type = g1.type
    @assert Set(vertices(g1)) == Set(vertices(g2)) "g1 and g2 have different external vertices."
    @assert Set(external_vertices(g1)) == Set(external_vertices(g2)) "g1 and g2 have different internal vertices."
    @assert g1.orders == g2.orders "g1 and g2 have different orders."

    return Graph{F,W}(g1.external_vertices, g1.internal_vertices; type=type, subgraph=[g1, g2], operator=Sum())
end

function Base.:-(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
    return g1 + (-1) * g2
end

function feynman_diagram(vertices::Vector{CompositeOperator}, contractions::Vector{Int};
    external::Vector{Int}=[],
    factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", type=:generic)

    operators = [vertices[i].operator for i in 1:length(extV)]
    comp_op = reduce(*, operators)
    sign_normal, permutation = normal_order(comp_op)
    edges, contraction_sign = contraction_to_edges(comp_op, contractions)

    return feynman_diagram(vertices, edges; external=external, factor=factor .* contraction_sign, weight=weight, name=name, type=type)
end

function feynman_diagram(vertices::Vector{CompositeOperator}, edges::Vector{Tuple{Int,Int}};
    external::Vector{Int}=[],
    factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", type=:generic)

    g = Graph(vertices; external=external, name=name, type=type, operator=Prod(), factor=factor, weight=weight)

    for edge in edges
        append!(g.subgraph, [])  # problem: how to determine vertices in the propagator? 
    end
end

"""
    Converts an ordered list of Wick contractions associated with an unspecified normal-ordered `CompositeOperator` to a list of edges.
    For a list of N Wick contractions, the first N entries of `contractions` should be: `1, 2, 3, ..., N`.

    Example: [1, 2, 3, 4, 1, 3, 4, 2] => [(1, 5), (2, 8), (3, 6), (4, 7)]

"""
function contractions_to_edges(operator::CompositeOperator, contractions::Vector{Int})
    # Filter illegal contractions
    is_invalid =
        isodd(length(contractions)) ||
        isodd(length(unique(contractions))) ||
        length(operator.operators) != length(contractions)
    if is_invalid
        throw(
            ArgumentError(
                "Input $contractions does not specify a legal set of Wick contractions",
            ),
        )
    end
    N = length(contractions) ÷ 2
    if contractions[1:N] != collect(1:N)
        throw(ArgumentError("Input $contractions is not normal-ordered"))
    end
    edges = Vector{EdgeType}()
    # Loop over annihilation operators (N < j ≤ 2N) and pair
    for (js, i) in enumerate(contractions[(N+1):end])
        j = js + N
        @debug "Wick contraction #$i, adding edge: ($i, $j)"
        push!(edges, (i, j))
    end
    contraction_sign = 1
    # TODO: Deduce the statistical sign associated with this normal-ordered contraction.
    #       We need to track the sign as we put the annihilation operators in anti-sorted
    return edges, contraction_sign
end

# function add_edge!(g::Graph{F,W}, contraction::Dict{Int,Vector{Int}}) where {F,W}
#     for vi in keys(contraction)
#         if vi <= length(g.external_vertices)
#             v = g.external_vertices[vi]
#         else
#             v = g.internal_vertices[vi-length(g.external_vertices)]
#         end
#         @assert length(v.operator) == length(contraction[vi])
#         for (loc_v, contr) in enumerate(contraction[vi])
#             contr == -1 && continue
#             for (ui, contrs) in contraction
#                 if contr in contrs
#                     loc_u = findlast(isequal(contr), contrs)
#                     #TODO: check the legality of the contraction
#                     if ui == vi
#                         if loc_v != loc_u
#                             v1 = ExternalVertex(v.current, [v.operator[loc_v], v.operator[loc_u]])
#                             append!(g.subgraph, [propagtor(v1, v1; ftype=F, wtype=W)])
#                         end
#                     else
#                         u = ui <= length(g.external_vertices) ? g.external_vertices[ui] : g.internal_vertices[ui-length(g.external_vertices)]
#                         v1 = ExternalVertex(v.current, [v.operator[loc_v]])
#                         v2 = ExternalVertex(u.current, [u.operator[loc_u]])
#                         append!(g.subgraph, [propagtor(v1, v2; ftype=F, wtype=W)])
#                     end
#                     contraction[vi][loc_v], contraction[ui][loc_u] = -1, -1
#                 end
#             end
#         end
#     end
# end

function propagator(v1::QuantumOperator, v2::QuantumOperator; name="propagator", diagtype=:propagtor,
    factor=one(_dtype.factor), weight=zero(_dtype.weight), operator=Sum())
    if v1 == v2
        extV = [v1] # Hatree bubble
    else
        extV = [v1, v2]
    end
    return Graph(extV, [], type=diagtype, name=name, operator=operator, factor=factor, weight=weight)
end

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

"""
    function _getVertices(V1::Vector{ExternalVertex}, V2::Vector{ExternalVertex}, couplings::Couplings)

    Give internal and external vertices when two diagrams combine.
  
    # Arguments
    - V1::Vector{ExternalVertex}            : External vertices of diagram I.
    - V2::Vector{ExternalVertex}            : External vertices of diagram II.
    - couplings::Vector{CompositeOperator}  : all the vertex couplings in the Feynman rule. 
"""
# function _getVertices(V1::Vector{ExternalVertex}, V2::Vector{ExternalVertex}, couplings::Vector{CompositeOperator})
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
#         point_ops = CompositeOperator([V1[i1].operator.operators; V2[i2].operator.operators])
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
#         extV = [ExternalVertex(point_in, current, CompositeOperator([opin, opout]))]
#     else
#         extV = [ExternalVertex(point_in, current, CompositeOperator([opin])),
#             ExternalVertex(point_out, current, CompositeOperator([opout]))]
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
#     ext_in = ExternalVertex(point_in, current, CompositeOperator([real_scalar(flavor),]))
#     ext_out = ExternalVertex(point_out, current, CompositeOperator([real_scalar(flavor),]))
#     diagtype = :interaction2
#     return Graph{dtype,dtype}([ext_in, ext_out], type=diagtype, couplings=couplings,
#         name=name, factor=factor, weight=weight)
# end

# function Interaction(point::Int, coupling::CompositeOperator, current::Int=0;
#     dtype=Float64, factor=zero(dtype), weight=zero(dtype), name="W", operators=[])
#     extV = ExternalVertex(point, current, CompositeOperator(operators))
#     # diagtype = Symbol("bareVertex$(length(operator))"...)
#     diagtype = :interaction
#     return Graph{dtype,dtype}([extV,], type=diagtype, couplings=[coupling,], name=name, factor=factor, weight=weight)
# end