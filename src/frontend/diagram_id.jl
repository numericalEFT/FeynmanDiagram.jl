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
    type::AnalyticProperty #Instant, Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function BareGreenId(type::AnalyticProperty=Dynamic; k, t)
        idx = findfirst(!iszero, k)
        if isnothing(idx) || k[idx] > 0
            return new(type, k, Tuple(t))
        else
            return new(type, -k, Tuple(t))
        end
    end
end
Base.show(io::IO, v::BareGreenId) = print(io, "$(short(v.type)), k$(v.extK), t$(v.extT)")

struct BareInteractionId <: PropagatorId # bare W-type interaction, with only one extK
    response::Response #UpUp, UpDown, ...
    type::AnalyticProperty #Instant, Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function BareInteractionId(response::Response, type::AnalyticProperty=Instant; k, t=(0, 0))
        idx = findfirst(!iszero, k)
        if isnothing(idx) || k[idx] > 0
            return new(response, type, k, Tuple(t))
        else
            return new(response, type, -k, Tuple(t))
        end
    end
end
Base.show(io::IO, v::BareInteractionId) = print(io, "$(short(v.response))$(short(v.type)), k$(v.extK), t$(v.extT)")

function Base.isequal(a::BareInteractionId, b::BareInteractionId)
    # Check if response, type, and extK are not equal
    if (a.response != b.response) || (a.type != b.type) || ((a.extK ≈ b.extK) == false)
        return false
    end

    # Check the conditions for Instant and Dynamic types
    # both Instant or Dynamic can have extT = [1, 1] or [1, 2]
    # This is because that Instant interaction may need an auxiliary time index to increase the number of the internal time variables to two.

    # if extT[1] == extT[2], that means the interaction is not time-dependent, then the specific time is not important

    # For example, if a.extT = [1, 1] and b.extT = [2, 2], then return true
    # Or, if a.extT = [1, 2] and b.extT = [1, 2], then return true
    # otherwise, return false

    return ((a.extT[1] == a.extT[2]) && (b.extT[1] == b.extT[2])) || (a.extT == b.extT)

    # If none of the conditions are met, return false
    return false
end

struct GenericId{P} <: DiagramId
    para::P
    extra::Any
    GenericId(para::P, extra=Nothing) where {P} = new{P}(para, extra)
end
Base.show(io::IO, v::GenericId) = print(io, v.extra == Nothing ? "" : "$(v.extra)")

struct GreenId{P} <: DiagramId
    para::P
    type::AnalyticProperty #Instant, Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function GreenId(para::P, type::AnalyticProperty=Dynamic; k, t) where {P}
        idx = findfirst(!iszero, k)
        if isnothing(idx) || k[idx] > 0
            return new{P}(para, type, k, Tuple(t))
        else
            return new{P}(para, type, -k, Tuple(t))
        end
    end
end
Base.show(io::IO, v::GreenId) = print(io, "$(short(v.type)), k$(v.extK), t$(v.extT)")

struct SigmaId{P} <: DiagramId
    para::P
    type::AnalyticProperty #Instant, Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function SigmaId(para::P, type::AnalyticProperty; k, t=(0, 0)) where {P}
        idx = findfirst(!iszero, k)
        if isnothing(idx) || k[idx] > 0
            return new{P}(para, type, k, Tuple(t))
        else
            return new{P}(para, type, -k, Tuple(t))
        end
    end
end
Base.show(io::IO, v::SigmaId) = print(io, "$(short(v.type))#$(v.order), t$(v.extT)")

struct PolarId{P} <: DiagramId
    para::P
    response::Response #UpUp, UpDown, ...
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    order::Vector{Int}
    function PolarId(para::P, response::Response; k, t=(0, 0)) where {P}
        idx = findfirst(!iszero, k)
        if isnothing(idx) || k[idx] > 0
            return new{P}(para, response, k, Tuple(t))
        else
            return new{P}(para, response, -k, Tuple(t))
        end
    end
end
Base.show(io::IO, v::PolarId) = print(io, "$(short(v.response)), k$(v.extK), t$(v.extT)")

struct Ver3Id{P} <: DiagramId
    para::P
    response::Response #UpUp, UpDown, ...
    extK::Vector{Vector{Float64}}
    extT::Tuple{Int,Int,Int} #all possible extT from different interactionType
    function Ver3Id(para::P, response::Response; k, t=(0, 0, 0)) where {P}
        return new{P}(para, response, k, Tuple(t))
    end
end
Base.show(io::IO, v::Ver3Id) = print(io, "$(short(v.response)),t$(v.extT)")

struct Ver4Id{P} <: DiagramId
    para::P
    response::Response #UpUp, UpDown, ...
    type::AnalyticProperty #Instant, Dynamic
    channel::TwoBodyChannel # particle-hole, particle-hole exchange, particle-particle, irreducible
    extK::Vector{Vector{Float64}}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
    function Ver4Id(para::P, response::Response, type::AnalyticProperty=Dynamic;
        k, t=(0, 0, 0, 0), chan::TwoBodyChannel=AnyChan) where {P}
        return new{P}(para, response, type, chan, k, Tuple(t))
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


"""
hopping function c⁺c⁻
"""
struct BareHoppingId{P} <: PropagatorId
    para::P
    site::Tuple{Int,Int}
    orbital::Tuple{Int,Int}
    extT::Tuple{Int,Int}
    function BareHoppingId(para::P, orbital, t, r) where {P}
        return new{P}(para, r, orbital, t)
    end
end
Base.show(io::IO, v::BareHoppingId) = print(io, "($(vstr(v.site, "ᵣ"))|$(vstr(v.orbital, "ₒ"))|$(vcstr(v.extT, [true, false])))")

"""
time-ordered N-point Bare Green's function
"""
struct BareGreenNId{P} <: PropagatorId
    para::P
    site::Int
    creation::Vector{Bool}
    orbital::Vector{Int}
    extT::Vector{Int}
    N::Int
    function BareGreenNId(para::P; orbital=[], t=[], creation=[], r=0) where {P}
        @assert length(orbital) == length(t) == length(creation)
        return new{P}(para, r, creation, orbital, t, length(orbital))
    end
end
Base.show(io::IO, v::BareGreenNId) = print(io, "($(v.site)ᵣ|$(vstr(v.orbital, "ₒ"))|$(vcstr(v.extT, v.creation)))")

"""
time-ordered N-point Composite Green's function
"""
struct GreenNId{P} <: DiagramId
    para::P
    site::Vector{Int}
    creation::Vector{Bool}
    orbital::Vector{Int}
    extT::Vector{Int}
    N::Int
    function GreenNId(para::P; orbital=[], t=[], creation=[], r=[]) where {P}
        @assert length(orbital) == length(t) == length(r) == length(creation)
        return new{P}(para, r, creation, orbital, t, length(orbital))
    end
end
Base.show(io::IO, v::GreenNId) = print(io, "($(vstr(v.site, "ᵣ"))|$(vstr(v.orbital, "ₒ"))|$(vcstr(v.extT, v.creation)))")

"""
time-ordered N-point Composite Green's function
"""
struct ConnectedGreenNId{P} <: DiagramId
    para::P
    site::Vector{Int}
    creation::Vector{Bool}
    orbital::Vector{Int}
    extT::Vector{Int}
    N::Int
    function ConnectedGreenNId(para::P; orbital=[], t=[], creation=[], r=[]) where {P}
        @assert length(orbital) == length(t) == length(r) == length(creation)
        return new{P}(para, r, creation, orbital, t, length(orbital))
    end
end
Base.show(io::IO, v::ConnectedGreenNId) = print(io, "($(vstr(v.site, "ᵣ"))|$(vstr(v.orbital, "ₒ"))|$(vcstr(v.extT, v.creation)))")


function Base.isequal(a::DiagramId, b::DiagramId)
    if typeof(a) != typeof(b)
        return false
    end
    for field in fieldnames(typeof(a))
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
    elseif type == BareGreenNId
        return 3
    elseif type == BareHoppingId
        return 4
    else
        error("Not Implemented!")
    end
end


