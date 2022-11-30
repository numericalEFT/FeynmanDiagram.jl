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

struct QuantumOperator#{T} <: AbstractVector{T}
    operator::Symbol
    flavor::Int
end
Base.isequal(a::QuantumOperator, b::QuantumOperator) = ((a.operator == b.operator) && (a.flavor == b.flavor))
Base.:(==)(a::QuantumOperator, b::QuantumOperator) = Base.isequal(a, b)

# const CompositeOperator = Vector{QuantumOperator}
const 𝑓 = QuantumOperator(:f⁻, 1)    #fermionic annihilation
const 𝑓dag = QuantumOperator(:f⁺, 1) #fermionic creation
const γ = QuantumOperator(:f, 1)    #fermionic real field (Majorana fermion)
const 𝑏 = QuantumOperator(:b⁻, 1)    #bosonic annihilation
const 𝑏dag = QuantumOperator(:b⁺, 1) #bosonic creation
const ϕ = QuantumOperator(:phi, 1) #classical real field (it can be paired with f⁺f⁻, b⁺b⁻, or other classical fields):w

# f = [QuantumOperator(:f⁻, i) for i in 1:3]
# f = fermionic_ann() # f is a function
# f(1), f(2), f(3)
fermionic_annihilation(flavor) = QuantumOperator(:f⁻, flavor)
fermionic_annihilation() = flavor -> QuantumOperator(:f⁻, flavor)
fermionic_creation(flavor) = QuantumOperator(:f⁺, flavor)
fermionic_creation() = flavor -> QuantumOperator(:f⁺, flavor)
majorana(flavor) = QuantumOperator(:f⁺, flavor)
majorana() = flavor -> QuantumOperator(:f⁺, flavor)
bosonic_annihilation(flavor) = QuantumOperator(:b⁻, flavor)
bosonic_annihilation() = flavor -> QuantumOperator(:b⁻, flavor)
bosonic_creation(flavor) = QuantumOperator(:b⁺, flavor)
bosonic_creation() = flavor -> QuantumOperator(:b⁺, flavor)
real_scalar(flavor) = QuantumOperator(:phi, flavor)
real_scalar() = flavor -> QuantumOperator(:phi, flavor)

function _countervector(it)
    y = Dict{eltype(it),Int}()
    for i in it
        y[i] = get(y, i, 0) + 1
    end
    return y
end

struct CompositeOperator
    operators::Vector{QuantumOperator}
end
Base.isequal(a::CompositeOperator, b::CompositeOperator) = (_countervector(a.operators) == _countervector(b.operators))
Base.:(==)(a::CompositeOperator, b::CompositeOperator) = Base.isequal(a, b)
const Coupling_yukawa = CompositeOperator([𝑓, 𝑓dag, ϕ])
const Coupling_phi3 = CompositeOperator([ϕ, ϕ, ϕ])
const Coupling_phi4 = CompositeOperator([ϕ, ϕ, ϕ, ϕ])
const Coupling_phi6 = CompositeOperator([ϕ, ϕ, ϕ, ϕ, ϕ, ϕ])

# """
#     struct Couplings

#     struct of all the vertex couplings in the Feynman rule. 
# # Members
# - operators::Vector{CompositeOperator}     : each CompositeOperator element stores the quantum operators inclued in one kind of vertex coupling 
# """
# struct Couplings
#     operators::Vector{CompositeOperator}
#     # TODO: function to testify the legality of the operators
# end
# const simp_coups = Couplings([[f, fdag, phi],])

abstract type Vertex end

struct ExternalVertex <: Vertex
    point::Int
    current::Int
    operator::CompositeOperator
    # operators::Vector{Symbol} # list of the composite operators,
    # :f for fermionic real field (Majorana fermion)
    # :b for bosonic real field
    # :phi for classical real field (it can be paired with f⁺f⁻, b⁺b⁻, or other classical fields)
    # :f⁺, :f⁻ for complex fermionic field
    # :b⁺, :b⁻ for complex bosonic field
    # flavors::Vector{Int} # flavor of each operator, it allows the field to be scalar, vector or even tensor
end
Base.isequal(a::ExternalVertex, b::ExternalVertex) = ((a.point == b.point) && (a.current == b.current) && (a.operator == b.operator))
Base.:(==)(a::ExternalVertex, b::ExternalVertex) = Base.isequal(a, b)

struct InternalVertex <: Vertex
    point::Int
    current::Int
    operator::CompositeOperator
end
Base.isequal(a::InternalVertex, b::InternalVertex) = ((a.point == b.point) && (a.current == b.current) && (a.operator == b.operator))
Base.:(==)(a::InternalVertex, b::InternalVertex) = Base.isequal(a, b)


"""
    mutable struct Diagram{F,W}
    
    struct of a Feynman diagram. A diagram of a sum or produce of various subdiagrams.

# Members
- hash::Int            : the unique hash number to identify the diagram
- name::Symbol         : name of the diagram
- para::DiagramPara    : internal parameters of the diagram
- orders::Vector{Int}  : orders of the diagram, loop order, derivative order, etc.
# - couplings::Couplings : all the vertex couplings in the Feynman rule. 
# - internal_points::Vector{Int} : internal points in the diagram
# - currents::Vector{Float64} : independent currents in the diagram
- external_vertices::Vector{ExternalVertex}    : external vertices of the diagram
- internal_vertices::Vector{InternalVertex}    : internal vertices of the diagram
# - isConnected::Bool    : connected or disconnected Green's function
# - isAmputated::Bool    : amputated Green's function or not
- subdiagram::Vector{Diagram{W}}   : vector of sub-diagrams 
- operator::Operator   : operation, support Sum() and Prod()
- factor::F            : additional factor of the diagram
- weight::W            : weight of the diagram
"""
mutable struct Diagram{F,W} # Diagram
    id::Int
    name::String # "" by default
    type::Symbol # :propagator, :interaction, :sigma, :green, :generic
    orders::Vector{Int}
    couplings::Vector{CompositeOperator}

    external_vertices::Vector{ExternalVertex}
    internal_vertices::Vector{InternalVertex}
    # is_connected::Bool
    subdiagram::Vector{Diagram{W}}

    operator::Operator
    factor::F
    weight::W

    function Diagram{F,W}(extV=[], intV=[]; subdiagram=[], couplings=[],
        name="", type=:generic, operator::Operator=Sum(), factor=F(1), weight=W(0)) where {F,W}
        orders = zeros(Int, 16)
        return new{F,W}(uid(), name, type, orders, couplings, extV, intV, subdiagram, operator, factor, weight)
    end
end

#TODO: improve a text representation of Diagram to the output stream.
function Base.show(io::IO, g::Diagram)
    # strc = g.is_connected ? "connected " : "disconnected "
    # stra = g.isAmputated ? "amputated " : " "
    print(io, "id $(g.id): Green's function $(g.name)")
end

function Base.isequal(a::Diagram, b::Diagram)
    typeof(a) != typeof(b) && return false
    for field in fieldnames(typeof(a))
        if field == :weight
            (getproperty(a, :weight) ≈ getproperty(b, :weight)) == false && return false
        end
        getproperty(a, field) != getproperty(b, field) && return false
    end
    return true
end
Base.:(==)(a::Diagram, b::Diagram) = Base.isequal(a, b)
# isbare(diag::Diagram) = isempty(diag.subdiagram)

#TODO: add function return reducibility of Diagram. 
function reducibility(g::Diagram)
    return (OneFermiIrreducible,)
end

#TODO: add function for connected diagram check. 
function connectivity(g::Diagram)
    isempty(g.subdiagram) && return true
end
function connectivity!(g::Diagram)
    g.isConnected = connectivity(g)
end

function Base.:*(g1::Diagram{F,W}, g2::Diagram{F,W}) where {F,W}
    type = :generic
    for v1 in g1.internal_vertices
        for v2 in g2.internal_vertices
            @assert v1.point != v2.point "g1 and g2 have the same internal vertex point."
        end
    end
    # ext1_ind = [v.point for v in g1.external_vertices]
    # ext2_ind = [v.point for v in g2.external_vertices]
    # has_common = !isempty(intersect(ext1_ind, ext2_ind))
    # total = union(g1.external_vertices, g2.external_vertices)
    # ext = [v for v in total if v.point ∉ common]
    couplings = union(g1.couplings, g2.couplings)
    extV, intV = _getVertices(g1.external_vertices, g2.external_vertices, couplings)
    return Diagram{F,W}(extV, [g1.internal_vertices; g2.internal_vertices; intV]; type=type, couplings=couplings, subdiagram=[g1, g2], operator=Prod())
end

function Base.:*(g1::Diagram{F,W}, c2::Number) where {F,W}
    return Diagram{F,W}(g1.external_vertices, g1.internal_vertices; type=g1.type, couplings=g1.couplings, subdiagram=[g1,], operator=Prod(), factor=c2)
end

function Base.:*(c1::Number, g2::Diagram{F,W}) where {F,W}
    return Diagram{F,W}(g2.external_vertices, g2.internal_vertices; type=g2.type, couplings=g2.couplings, subdiagram=[g2,], operator=Prod(), factor=c1)
end

function Base.:+(g1::Diagram{F,W}, g2::Diagram{F,W}) where {F,W}
    @assert g1.type == g2.type "g1 and g2 are not of the same type."
    # @assert g1.isAmputated == g2.isAmputated "g1 and g2 are not of the same amputated status."
    # TODO: more check
    type = g1.type
    @assert Set(g1.external_vertices) == Set(g2.external_vertices) "g1 and g2 have different external vertices."
    @assert Set(g1.internal_vertices) == Set(g2.internal_vertices) "g1 and g2 have different internal vertices."
    #TODO: add external vertices creation/annihilation check
    return Diagram{F,W}(g1.external_vertices, g1.internal_vertices; type=type, couplings=g1.couplings, subdiagram=[g1, g2], operator=Sum())
end

function Base.:-(g1::Diagram{F,W}, g2::Diagram{F,W}) where {F,W}
    return g1 + (-1) * g2
end

function checkVertices(g::Diagram{F,W}) where {F,W}
    @assert !isempty(g.couplings) "Diagram.couplings must be defined to check the legality of vertices"
    for v in g.internal_vertices
        @assert v.operator in g.couplings "internal vertex point $(v.point) is illegal."
    end
    for v in g.external_vertices
        is_legal = false
        for ops in g.couplings
            @assert v.operator != ops "external vertex point $(v.point) is illegal."
            num_pops, num_ops = _countervector(v.operator.operators), _countervector(ops.operators)
            all([num_ops[op] >= num_pops[op] for op in v.operator.operators]) && (is_legal = true)
        end
        @assert is_legal "external vertex point $(v.point) is illegal."
    end
    return true
end

"""
    function _getVertices(V1::Vector{ExternalVertex}, V2::Vector{ExternalVertex}, couplings::Couplings)

    Give internal and external vertices when two diagrams combine.
  
    # Arguments
    - V1::Vector{ExternalVertex}            : External vertices of diagram I.
    - V2::Vector{ExternalVertex}            : External vertices of diagram II.
    - couplings::Vector{CompositeOperator}  : all the vertex couplings in the Feynman rule. 
"""
function _getVertices(V1::Vector{ExternalVertex}, V2::Vector{ExternalVertex}, couplings::Vector{CompositeOperator})
    V1_ind = [v.point for v in V1]
    V2_ind = [v.point for v in V2]
    common = intersect(V1_ind, V2_ind)
    total = union(V1, V2)
    intV = []
    extV = [v for v in total if v.point ∉ common]
    for point in common
        # ifinternal, ifexternal = false, false
        i1 = findfirst(isequal(point), V1_ind)
        i2 = findfirst(isequal(point), V2_ind)
        point_ops = CompositeOperator([V1[i1].operator.operators; V2[i2].operator.operators])
        if isempty(couplings)
            append!(extV, [ExternalVertex(point, V1[i1].current, point_ops)])
        else
            for ops in couplings
                if point_ops == ops
                    append!(intV, [InternalVertex(point, V1[i1].current, point_ops)])
                else
                    append!(extV, [ExternalVertex(point, V1[i1].current, point_ops)])
                end
                # num_pops, num_ops = _countervector(point_ops), _countervector(ops.operators)
                # if num_pops == num_ops
                #     ifinternal = true
                #     break
                # end
                # if point_ops ⊆ ops.operators
                #     if all([num_ops[op] >= num_pops[op] for op in point_ops])
                #         append!(extV, [ExternalVertex(point, V1[i1].current, CompositeOperator(point_ops))])
                #         # ifexternal = true
                #         break
                #     end
                # end
            end
        end
        # @assert ifinternal || ifexternal "point $point is illegal."
    end
    return extV, intV
end

function 𝐺ᶠ(point_in::Int, point_out::Int, current::Int=0; kwargs...)
    return Green2(point_in, point_out, current; isFermi=true, kwargs...)
end

function 𝐺ᵇ(point_in::Int, point_out::Int, current::Int=0; kwargs...)
    return Green2(point_in, point_out, current; isFermi=false, kwargs...)
end

function 𝐺ᵠ(point_in::Int, point_out::Int, current::Int=0; kwargs...)
    return Green2(point_in, point_out, current; isFermi=false, isComplex=false, kwargs...)
end

function Green2(point_in::Int, point_out::Int, current::Int=0;
    isFermi=true, isComplex=true, flavor::Int=1, couplings=[], dtype=Float64,
    factor=zero(dtype), weight=zero(dtype), name="G2", subdiagram=[], operator=Sum())
    if isFermi && isComplex
        opin, opout = fermionic_creation(flavor), fermionic_annihilation(flavor)
    elseif isFermi && !isComplex
        opin, opout = majorana(flavor), majorana(flavor)
    elseif !isFermi && isComplex
        opin, opout = bosonic_creation(flavor), bosonic_annihilation(flavor)
    else
        opin, opout = real_scalar(flavor), real_scalar(flavor)
    end
    if point_in == point_out
        extV = [ExternalVertex(point_in, current, CompositeOperator([opin, opout]))]
    else
        extV = [ExternalVertex(point_in, current, CompositeOperator([opin])),
            ExternalVertex(point_out, current, CompositeOperator([opout]))]
    end
    if isnothing(subdiagram)
        diagtype = :propagator
    else
        diagtype = :green2
    end
    return Diagram{dtype,dtype}(extV, type=diagtype, couplings=couplings, subdiagram=subdiagram,
        name=name, operator=operator, factor=factor, weight=weight)
end

const 𝑊 = Interaction

function Interaction(point_in::Int, point_out::Int, current::Int=0; flavor::Int=1,
    couplings=[Coupling_yukawa,], dtype=Float64, factor=zero(dtype), weight=zero(dtype), name="W")
    ext_in = ExternalVertex(point_in, current, CompositeOperator([real_scalar(flavor),]))
    ext_out = ExternalVertex(point_out, current, CompositeOperator([real_scalar(flavor),]))
    diagtype = :interaction2
    return Diagram{dtype,dtype}([ext_in, ext_out], type=diagtype, couplings=couplings,
        name=name, factor=factor, weight=weight)
end

function Interaction(point::Int, coupling::CompositeOperator, current::Int=0;
    dtype=Float64, factor=zero(dtype), weight=zero(dtype), name="W", operators=[])
    extV = ExternalVertex(point, current, CompositeOperator(operators))
    # diagtype = Symbol("bareVertex$(length(operator))"...)
    diagtype = :interaction
    return Diagram{dtype,dtype}([extV,], type=diagtype, couplings=[coupling,], name=name, factor=factor, weight=weight)
end

# function Diagram{W}(::Type{Vacuum}; isConnected=false, isAmputated=true, extV=[],
#     subdiagram=[], reducibility=[], name=:VacuumDiagram, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
#     @assert para.type == Vacuum "types from input and para are inconsistent."
#     @assert length(extV) == 0 "input parameters do not support Vacuum diagram."
#     return Diagram{W}(Vaccum, isConnected, isAmputated, [], subdiagram, reducibility=reducibility,
#         name=name, operator=operator, factor=factor, weight=weight)
# end

# function Diagram{W}(::Type{Tadpole}, para::DiagramPara; isConnected=true, isAmputated=false, extV=[],
#     subdiagram=[], reducibility=[], name=:TadpoleDiagram, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
#     @assert para.type == Tadpole "types from input and para are inconsistent."
#     @assert length(extV) == 1 "input parameters do not support Tadpole diagram."
#     return Diagram{W}(para, isConnected, isAmputated, extV, subdiagram, reducibility=reducibility,
#         name=name, operator=operator, factor=factor, weight=weight)
# end

# function Diagram{W}(::Type{FermiPropagator}, para::DiagramPara; isConnected=true, isAmputated=false, extV=[],
#     subdiagram=[], reducibility=[OneFermiIrreducible,], name=:FermiPropagator, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
#     @assert para.type == FermiPropagator "types from input and para are inconsistent."
#     @assert length(extV) == 2 && extV[1].isFermi && extV[2].isFermi &&
#             ((extV[1].isCreation && !extV[2].isCreation) || (!extV[1].isCreation && extV[2].isCreation)) "input parameters do not support FermiPropagator."
#     return Diagram{W}(para, isConnected, isAmputated, extV, subdiagram, reducibility=reducibility,
#         name=name, operator=operator, factor=factor, weight=weight)
# end

# function Diagram{W}(::Type{BosePropagator}, para::DiagramPara; isConnected=true, isAmputated=false, extV=[],
#     subdiagram=[], reducibility=[OneBoseIrreducible,], name=:BosePropagator, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
#     @assert para.type == BosePropagator "types from input and para are inconsistent."
#     @assert length(extV) == 2 && !extV[1].isFermi && !extV[2].isFermi &&
#             ((extV[1].isCreation && !extV[2].isCreation) || (!extV[1].isCreation && extV[2].isCreation)) "input parameters do not support BosePropagator."
#     return Diagram{W}(para, isConnected, isAmputated, extV, subdiagram, reducibility=reducibility,
#         name=name, operator=operator, factor=factor, weight=weight)
# end

# function Diagram{W}(::Type{FermiSelfEnergy}, para::DiagramPara; isConnected=true, isAmputated=true, extV=[],
#     subdiagram=[], reducibility=[OneFermiIrreducible,], name=:FermiSelfEnergy, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
#     @assert para.type == FermiSelfEnergy "types from input and para are inconsistent."
#     @assert length(extV) == 2 && extV[1].isFermi && extV[2].isFermi "input parameters do not support FermiSelfEnergy."
#     return Diagram{W}(para, isConnected, isAmputated, extV, subdiagram, reducibility=reducibility,
#         name=name, operator=operator, factor=factor, weight=weight)
# end

# function Diagram{W}(::Type{BoseSelfEnergy}, para::DiagramPara; isConnected=true, isAmputated=true, extV=[],
#     subdiagram=[], reducibility=[OneFermiIrreducible,], name=:BoseSelfEnergy, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
#     @assert para.type == BoseSelfEnergy "types from input and para are inconsistent."
#     @assert length(extV) == 2 && !extV[1].isFermi && !extV[2].isFermi "input parameters do not support BoseSelfEnergy."
#     return Diagram{W}(para, isConnected, isAmputated, extV, subdiagram, reducibility=reducibility,
#         name=name, operator=operator, factor=factor, weight=weight)
# end

# function Diagram{W}(::Type{VertexDiag}, para::DiagramPara; Nf=2, Nb=1, isConnected=true, isAmputated=true, extV=[],
#     subdiagram=[], reducibility=[OneFermiIrreducible, OneBoseirreducible], name=:VertexDiagram, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
#     @assert para.type == VertexDiag "types from input and para are inconsistent."
#     @assert length(extV) > 2 && length(extV) == Nf + Nb "input parameters do not support Vertex diagram."
#     return Diagram{W}(para, isConnected, isAmputated, extV, subdiagram, reducibility=reducibility,
#         name=name, operator=operator, factor=factor, weight=weight)
# end

# function Diagram{W}(::Type{GncDiag}, para::DiagramPara; N=4, isConnected=true, isAmputated=false, extV=[],
#     subdiagram=[], reducibility=[], name=:ConnectedNDiagram, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
#     @assert para.type == GncDiag "types from input and para are inconsistent."
#     @assert length(extV) > 2 && length(extV) == N "input parameters do not support Connected N-point diagram."
#     return Diagram{W}(para, isConnected, isAmputated, extV, subdiagram, reducibility=reducibility,
#         name=name, operator=operator, factor=factor, weight=weight)
# end

# function Diagram{W}(::Type{GndDiag}, para::DiagramPara; N=4, isConnected=false, isAmputated=false, extV=[],
#     subdiagram=[], reducibility=[], name=:DisconnectedNDiagram, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
#     @assert para.type == GndDiag "types from input and para are inconsistent."
#     @assert length(extV) > 2 && length(extV) == N "input parameters do not support Disconnected N-point diagram."
#     return Diagram{W}(para, isConnected, isAmputated, extV, subdiagram, reducibility=reducibility,
#         name=name, operator=operator, factor=factor, weight=weight)
# end
