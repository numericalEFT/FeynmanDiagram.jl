"""
    abstract type DiagramId end

    The abstract type of all diagrams/subdiagrams/bare propagators
"""
abstract type DiagramId end

"""
    abstract type PropagatorId <: DiagramId end

    The abstract type of all bare propagators
"""
abstract type PropagatorId <: DiagramId end

# Base.Dict(x::DiagramId) = Dict{Symbol,Any}([fn => getfield(x, fn) for fn ∈ fieldnames(typeof(x))])
# Base.show(io::IO, d::DiagramId) = error("Base.show not implemented!")
# Base.isequal(a::DiagramId, b::DiagramId) = error("Base.isequal not implemented!")
Base.:(==)(a::DiagramId, b::DiagramId) = Base.isequal(a, b)

struct BareGreenId <: PropagatorId
    para::DiagPara
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function BareGreenId(para::DiagPara, type::AnalyticProperty=Dynamic; k, t)
        return new(para, type, k, Tuple(t))
    end
end
Base.show(io::IO, v::BareGreenId) = print(io, "$(short(v.type)), k$(v.extK), t$(v.extT)")

struct BareInteractionId <: PropagatorId # bare W-type interaction, with only one extK
    para::DiagPara
    response::Response #UpUp, UpDown, ...
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function BareInteractionId(para::DiagPara, response::Response, type::AnalyticProperty=Instant; k, t=(0, 0))
        return new(para, response, type, k, Tuple(t))
    end
end
Base.show(io::IO, v::BareInteractionId) = print(io, "$(short(v.response))$(short(v.type)), k$(v.extK), t$(v.extT)")

struct GenericId <: DiagramId
    para::DiagPara
    extra::Any
    GenericId(para::DiagPara, extra=Nothing) = new(para, extra)
end
Base.show(io::IO, v::GenericId) = print(io, v.extra == Nothing ? "" : "$(v.extra)")

struct GreenId <: DiagramId
    para::DiagPara
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function GreenId(para::DiagPara, type::AnalyticProperty=Dynamic; k, t)
        return new(para, type, k, Tuple(t))
    end
end
Base.show(io::IO, v::GreenId) = print(io, "$(short(v.type)), k$(v.extK), t$(v.extT)")

struct SigmaId <: DiagramId
    para::DiagPara
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function SigmaId(para::DiagPara, type::AnalyticProperty; k, t=(0, 0))
        return new(para, type, k, t)
    end
end
Base.show(io::IO, v::SigmaId) = print(io, "$(short(v.type))#$(v.order), t$(v.extT)")

struct PolarId <: DiagramId
    para::DiagPara
    response::Response #UpUp, UpDown, ...
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    order::Vector{Int}
    function PolarId(para::DiagPara, response::Response; k, t=(0, 0))
        return new(para, response, k, t)
    end
end
Base.show(io::IO, v::PolarId) = print(io, "$(short(v.response)), k$(v.extK), t$(v.extT)")

struct Ver3Id <: DiagramId
    para::DiagPara
    response::Response #UpUp, UpDown, ...
    extK::Vector{Vector{Float64}}
    extT::Tuple{Int,Int,Int} #all possible extT from different interactionType
    function Ver3Id(para::DiagPara, response::Response; k, t=(0, 0, 0))
        return new(para, response, k, Tuple(t))
    end
end
Base.show(io::IO, v::Ver3Id) = print(io, "$(short(v.response)),t$(v.extT)")

struct Ver4Id <: DiagramId
    para::DiagPara
    response::Response #UpUp, UpDown, ...
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    channel::TwoBodyChannel # particle-hole, particle-hole exchange, particle-particle, irreducible
    extK::Vector{Vector{Float64}}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
    function Ver4Id(para::DiagPara, response::Response, type::AnalyticProperty=Dynamic; k, t=(0, 0, 0, 0), chan::TwoBodyChannel=AnyChan)
        return new(para, response, type, chan, k, Tuple(t))
    end
end
Base.show(io::IO, v::Ver4Id) = print(io, (v.channel == AnyChan ? "" : "$(v.channel) ") * "$(short(v.response))$(short(v.type)),t$(v.extT)")

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

function Base.isequal(a::DiagramId, b::DiagramId)
    if typeof(a) != typeof(b)
        return false
    end
    for field in fieldnames(typeof(a))
        # field in [:para, :permutation] && continue #both parameter and permutation needs to be compared
        # if field == :extK
        #     if !(getproperty(a, :extK) ≈ getproperty(b, :extK)) && !(getproperty(a, :extK) ≈ -getproperty(b, :extK))
        #         return false
        #     end
        #     continue
        # end
        if getproperty(a, field) != getproperty(b, field)
            return false
        end
    end
    return true
end

function index(type)
    if type == BareGreenId
        return 1
    elseif type == BareInteractionId
        return 2
    else
        error("Not Implemented!")
    end
end


