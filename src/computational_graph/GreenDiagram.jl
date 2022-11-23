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

@enum Reducibility begin
    OneFermiIrreducible
    OneBoseIrreducible
    ParticleHoleIrreducible
    ParticleParticleIrreducible
end

struct ExternalVertice
    point::Int
    current::Int
    isCreation::Bool
    isFermi::Bool
end

"""
    mutable struct GreenDiagram{W}
    
    struct of a Feynman diagram. A diagram of a sum or produce of various subdiagrams.

# Members
- hash::Int           : the unique hash number to identify the diagram
- name::Symbol        : name of the diagram
- para::DiagramPara   : internal parameters of the diagram
- orders::Vector{Int} : orders of the diagram, loop order, derivative order, etc.
- internal_points::Vector{Int} : internal points in the diagram
- currents::Vector{Float64} : independent currents in the diagram
- extVertices::Vector{ExternalVertice}    : external vertices of the diagram
- isConnected::Bool   : connected or disconnected Green's function
- isAmputated::Bool   : amputated Green's function or not
- subdiagram::Vector{GreenDiagram{W}}   : vector of sub-diagrams 
- operator::Operator  : operation, support Sum() and Prod()
- factor::W           : additional factor of the diagram
- weight::W           : weight of the diagram
"""
mutable struct GreenDiagram{W} # GreenDiagram
    hash::Int
    name::Symbol
    type::DataType
    para::DiagramPara
    orders::Vector{Int}
    internal_points::Vector{Int}
    currents::Vector{Float64}

    extVertices::Vector{ExternalVertice}
    isConnected::Bool
    isAmputated::Bool
    reducibility::Vector{Reducibility}
    subdiagram::Vector{GreenDiagram{W}}

    operator::Operator
    factor::W
    weight::W

    function GreenDiagram{W}(para::DiagramPara, isConnected, isAmputated, extV=[], subdiagram=[];
        type=para.type, reducibility=[], name=:GreenDiagram, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
        @assert type <: DiagType "$type is not implemented in DiagType."
        orders = zeros(Int, 16)
        g = new{W}(uid(), name, type, para, orders, [], [], extV, isConnected, isAmputated, reducibility, subdiagram, operator, factor, weight)
        reducibility!(g)
        return g
    end
end

function Base.show(io::IO, g::GreenDiagram)
    strc = g.isConnected ? "connected " : "disconnected "
    stra = g.isAmputated ? "amputated " : " "
    print(io, "hash $(g.hash): " * strc * stra * "Green's function $(g.name)")
end

Base.:(==)(a::GreenDiagram, b::GreenDiagram) = Base.isequal(a, b)
function Base.isequal(a::GreenDiagram, b::GreenDiagram)
    typeof(a) != typeof(b) && return false
    for field in fieldnames(typeof(a))
        if field == :weight
            (getproperty(a, :weight) ≈ getproperty(b, :weight)) == false && return false
        end
        getproperty(a, field) != getproperty(b, field) && return false
    end
    return true
end
# isbare(diag::GreenDiagram) = isempty(diag.subdiagram)

function reducibility(g::GreenDiagram{W}) where {W}
    return ()
end

function reducibility!(g::GreenDiagram{W}) where {W}
    # @assert g.reducibility ⊆ reducibility(g) "wrong reducibility in g (hash: $g.hash)."
    g.reducibility = union(reducibility(g), g.reducibility)
end

function Base.:*(g1::GreenDiagram{W}, g2::GreenDiagram{W}) where {W}
end

function Base.:*(g1::GreenDiagram{W}, c2::Number) where {W}
end

function Base.:*(c1::Number, g2::GreenDiagram{W}) where {W}
end

function Base.:+(g1::GreenDiagram{W}, g2::GreenDiagram{W}) where {W}
end

function Base.:-(g1::GreenDiagram{W}, g2::GreenDiagram{W}) where {W}
end
# g = Sigma(...)
# w = W(...)
# ver4 = Vertex4(...)

# graph = g*w+ver4*0.5

function GreenDiagram{W}(::Type{Vacuum}, para::DiagramPara; isConnected=false, isAmputated=true, extV=[],
    subdiagram=[], reducibility=[], name=:VacuumDiagram, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
    @assert para.type == Vacuum "types from input and para are inconsistent."
    @assert length(extV) == 0 "input parameters do not support Vacuum diagram."
    return GreenDiagram{W}(para, isConnected, isAmputated, [], subdiagram, reducibility=reducibility,
        name=name, operator=operator, factor=factor, weight=weight)
end

function GreenDiagram{W}(::Type{Tadpole}, para::DiagramPara; isConnected=true, isAmputated=false, extV=[],
    subdiagram=[], reducibility=[], name=:TadpoleDiagram, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
    @assert para.type == Tadpole "types from input and para are inconsistent."
    @assert length(extV) == 1 "input parameters do not support Tadpole diagram."
    return GreenDiagram{W}(para, isConnected, isAmputated, extV, subdiagram, reducibility=reducibility,
        name=name, operator=operator, factor=factor, weight=weight)
end

function GreenDiagram{W}(::Type{FermiPropagator}, para::DiagramPara; isConnected=true, isAmputated=false, extV=[],
    subdiagram=[], reducibility=[OneFermiIrreducible,], name=:FermiPropagator, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
    @assert para.type == FermiPropagator "types from input and para are inconsistent."
    @assert length(extV) == 2 && extV[1].isFermi && extV[2].isFermi &&
            ((extV[1].isCreation && !extV[2].isCreation) || (!extV[1].isCreation && extV[2].isCreation)) "input parameters do not support FermiPropagator."
    return GreenDiagram{W}(para, isConnected, isAmputated, extV, subdiagram, reducibility=reducibility,
        name=name, operator=operator, factor=factor, weight=weight)
end

function GreenDiagram{W}(::Type{BosePropagator}, para::DiagramPara; isConnected=true, isAmputated=false, extV=[],
    subdiagram=[], reducibility=[OneBoseIrreducible,], name=:BosePropagator, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
    @assert para.type == BosePropagator "types from input and para are inconsistent."
    @assert length(extV) == 2 && !extV[1].isFermi && !extV[2].isFermi &&
            ((extV[1].isCreation && !extV[2].isCreation) || (!extV[1].isCreation && extV[2].isCreation)) "input parameters do not support BosePropagator."
    return GreenDiagram{W}(para, isConnected, isAmputated, extV, subdiagram, reducibility=reducibility,
        name=name, operator=operator, factor=factor, weight=weight)
end

function GreenDiagram{W}(::Type{FermiSelfEnergy}, para::DiagramPara; isConnected=true, isAmputated=true, extV=[],
    subdiagram=[], reducibility=[OneFermiIrreducible,], name=:FermiSelfEnergy, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
    @assert para.type == FermiSelfEnergy "types from input and para are inconsistent."
    @assert length(extV) == 2 && extV[1].isFermi && extV[2].isFermi "input parameters do not support FermiSelfEnergy."
    return GreenDiagram{W}(para, isConnected, isAmputated, extV, subdiagram, reducibility=reducibility,
        name=name, operator=operator, factor=factor, weight=weight)
end

function GreenDiagram{W}(::Type{BoseSelfEnergy}, para::DiagramPara; isConnected=true, isAmputated=true, extV=[],
    subdiagram=[], reducibility=[OneFermiIrreducible,], name=:BoseSelfEnergy, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
    @assert para.type == BoseSelfEnergy "types from input and para are inconsistent."
    @assert length(extV) == 2 && !extV[1].isFermi && !extV[2].isFermi "input parameters do not support BoseSelfEnergy."
    return GreenDiagram{W}(para, isConnected, isAmputated, extV, subdiagram, reducibility=reducibility,
        name=name, operator=operator, factor=factor, weight=weight)
end

function GreenDiagram{W}(::Type{VertexDiag}, para::DiagramPara; Nf=2, Nb=1, isConnected=true, isAmputated=true, extV=[],
    subdiagram=[], reducibility=[OneFermiIrreducible, OneBoseirreducible], name=:VertexDiagram, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
    @assert para.type == VertexDiag "types from input and para are inconsistent."
    @assert length(extV) > 2 && length(extV) == Nf + Nb "input parameters do not support Vertex diagram."
    return GreenDiagram{W}(para, isConnected, isAmputated, extV, subdiagram, reducibility=reducibility,
        name=name, operator=operator, factor=factor, weight=weight)
end

function GreenDiagram{W}(::Type{GncDiag}, para::DiagramPara; N=4, isConnected=true, isAmputated=false, extV=[],
    subdiagram=[], reducibility=[], name=:ConnectedNDiagram, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
    @assert para.type == GncDiag "types from input and para are inconsistent."
    @assert length(extV) > 2 && length(extV) == N "input parameters do not support Connected N-point diagram."
    return GreenDiagram{W}(para, isConnected, isAmputated, extV, subdiagram, reducibility=reducibility,
        name=name, operator=operator, factor=factor, weight=weight)
end

function GreenDiagram{W}(::Type{GndDiag}, para::DiagramPara; N=4, isConnected=false, isAmputated=false, extV=[],
    subdiagram=[], reducibility=[], name=:DisconnectedNDiagram, operator::Operator=Sum(), factor=W(1), weight=W(0)) where {W}
    @assert para.type == GndDiag "types from input and para are inconsistent."
    @assert length(extV) > 2 && length(extV) == N "input parameters do not support Disconnected N-point diagram."
    return GreenDiagram{W}(para, isConnected, isAmputated, extV, subdiagram, reducibility=reducibility,
        name=name, operator=operator, factor=factor, weight=weight)
end
