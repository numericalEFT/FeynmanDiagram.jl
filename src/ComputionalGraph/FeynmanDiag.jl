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


@enum DiagType begin
    Vacuum
    Tadpole
    Full
    Connected
    Amputated
    OneFermiIrreducible
    OneBoseIrreducible
    ParticleHoleIrreducible
    ParticleParticleIrreducible
end

@enum Filter begin
    NoHartree
    NoFock
    NoBubble  # true to remove all bubble subdiagram
    Proper  #ver4, ver3, and polarization diagrams may require to be irreducible along the transfer momentum/frequency
end

struct Leg
    point::Int
    current::Vector{Float64}
    isCreation::Bool
    isFermi::Bool
    type::DiagType
end

#abstract type DiagType end

#struct Full<:DiagType
#struct OneFermiIrreducible <: DiagType
#struct OneBoseIrreducible <: DiagType

struct GenericGreen{P,W}
    hash::Int
    name::Symbol
    para::P
    legs::Vector{Leg}
    type::Vector{DiagType}
    subdiagram::Vector{GenericGreen{P}}

    factor::W
    weight::W
    function GenericGreen{P}(para::P, legs=[], subdiagram=[]; name=:none, type=[],
        factor=W(1), weight=W(0)) where {P,W}
        return new{P}(uid(), name, para, legs, type, subdiagram, factor, weight)
    end
    #para::DiagPara
    #point::Vector{Int}
    #current::Vector{C}
    #isCreation::Vector{Bool}
    # orbital::Vector{Int}
    #isFermi::Vector{Bool}
    #isAmputate::Vector{Bool}
    #isRenorm::Vector{Bool}
    #filter::Any
    # function ConnectedGreenNId(para::DiagPara, N::Int; orbital=[], t=[], creation=[], r=[], isFermi=[], isAmputate=[], isRenorm=[], filter=nothing)
    #     @assert length(orbital) == length(t) == length(r) == length(creation) == length(isFermi) == length(isAmputate) == length(isRenorm) == N
    #     return new(para, N, r, t, creation, orbital, isFermi, isAmputate, isRenorm)
    # end
end

Base.show(io::IO, g::GenericGreen) = print(io, "... Diagram")
isbare(diag::GenericGreen) = isempty(diag.subdiagram)
Base.:(==)(a::GenericGreen, b::GenericGreen) = Base.isequal(a, b)

# Base.show(io::IO, v::DiagramId) = print(io, "($(vstr(v.site, "ᵣ"))|$(vstr(v.orbital, "ₒ"))|$(vcstr(v.extT, v.creation)))")

function Base.:*(g1::GenericGreen{P,W}, g2::GenericGreen{P}) where {P,W}
end

function Base.:*(g1::GenericGreen{P,W}, c2::Number) where {P,W}
end

function Base.:*(c1::Number, g2::GenericGreen{P,W}) where {P,W}
end

function Base.:+(g1::GenericGreen{P,W}, g2::GenericGreen{P,W}) where {P,W}
end

function Base.:-(g1::GenericGreen{P,W}, g2::GenericGreen{P,W}) where {P,W}
end
# g = Sigma(...)
# w = W(...)
# ver4 = Vertex4(...)

# graph = g*w+ver4*0.5

#struct PropagatorId <: DiagramId end
#struct SigmaId <: DiagramId end
#struct PolarId <: DiagramId end
#struct Ver3Id <: VertexId end
##struct Ver4Id <: VertexId end
#struct VerNId <: VertexId end

function FermiPropagator{P,W}(para::P, legs, subdiagram=[]; name=:none, extratype=[], factor::W, weight::W) where {P,W}
    @assert length(legs) == 2 && legs[1].isFermi && legs[2].isFermi && ((legs[1].isCreation && !legs[2].isCreation)
                                                                        ||
                                                                        (!legs[1].isCreation && legs[2].isCreation)) "Error: input parameters do not support FermiPropagator."
    type = union([Connected, OneFermiIrreducible], extratype)
    return GenericGreen{P,W}(para, legs, subdiagram; name, type, factor, weight)
end

function BosePropagator{P,W}(para::P, legs, subdiagram=[]; name=:none, extratype=[], factor::W, weight::W) where {P,W}
    @assert length(legs) == 2 && !legs[1].isFermi && !legs[2].isFermi && ((legs[1].isCreation && !legs[2].isCreation)
                                                                          ||
                                                                          (!legs[1].isCreation && legs[2].isCreation)) "Error: input parameters do not support BosePropagator."
    type = union([Connected, OneBoseIrreducible], extratype)
    return GenericGreen{P,W}(para, legs, subdiagram; name, type, factor, weight)
end

function FermiSelfEnergy{P,W}(para::P, legs, subdiagram=[]; name=:none, extratype=[], factor::W, weight::W) where {P,W}
    @assert length(legs) == 2 && legs[1].isFermi && legs[2].isFermi && Amputated in legs[1].type &&
            Amputated in legs[2].type "Error: input parameters do not support FermiSelfEnergy."
    type = union([Connected, Amputated, OneFermiIrreducible], extratype)
    return GenericGreen{P,W}(para, legs, subdiagram; name, type, factor, weight)
end

function BoseSelfEnergy{P,W}(para::P, legs, subdiagram=[]; name=:none, extratype=[], factor::W, weight::W) where {P,W}
    @assert length(legs) == 2 && !legs[1].isFermi && !legs[2].isFermi && Amputated in legs[1].type &&
            Amputated in legs[2].type "Error: input parameters do not support BoseSelfEnergy (Polarization)."
    type = union([Connected, Amputated, OneBoseIrreducible], extratype)
    return GenericGreen{P,W}(para, legs, subdiagram; name, type, factor, weight)
end


function Vertex{P,W}(para::P, N::Int, legs, subdiagram=[]; name=:none, extratype=[], factor::W, weight::W) where {P,W}
    @assert length(legs) == N && Amputated in legs[1].type "Error: input parameters do not support Vertex."
    # or (!legs[1].isCreation && legs[2].isCreation)) "Error: input parameters do not support FermiPropagator."
    type = union([Connected, Amputated, OneFermiIrreducible], extratype)
    return GenericGreen{P,W}(para, legs, subdiagram; name, type, factor, weight)
end


function Base.isequal(a::GenericGreen, b::GenericGreen)
    if typeof(a) != typeof(b)
        return false
    end
    for field in fieldnames(typeof(a))
        if field == :extK
            if (getproperty(a, :extK) ≈ getproperty(b, :extK)) == false
                return false
            end
        end
        if getproperty(a, field) != getproperty(b, field)
            return false
        end
    end
    return true
end
