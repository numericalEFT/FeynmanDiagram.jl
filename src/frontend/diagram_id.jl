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
Base.:(==)(a::DiagramId, b::DiagramId) = Base.isequal(a, b)

struct BareGreenId <: PropagatorId
    type::AnalyticProperty #Instant, Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function BareGreenId(type::AnalyticProperty, k::Vector{T}, t::Tuple{Int,Int}) where {T<:Real}
        return new(type, mirror_symmetrize(k), t)
    end
    function BareGreenId(type::AnalyticProperty=Dynamic; k, t)
        return new(type, mirror_symmetrize(k), Tuple(t))
    end
end
Base.show(io::IO, v::BareGreenId) = print(io, "$(short(v.type)), k$(v.extK), t$(v.extT)")
function Base.hash(v::BareGreenId, h::UInt)
    h = hash(BareGreenId, h)
    h = hash(v.type, h)
    h = hash(v.extK, h)
    h = hash(v.extT, h)
    return h
end
function Base.isequal(a::BareGreenId, b::BareGreenId)
    return a.type == b.type && isequal(a.extK, b.extK) && a.extT == b.extT
end

struct BareInteractionId <: PropagatorId # bare W-type interaction, with only one extK
    response::Response #UpUp, UpDown, ...
    type::AnalyticProperty #Instant, Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function BareInteractionId(response::Response, type::AnalyticProperty, k::Vector{T}, t::Tuple{Int,Int}) where {T<:Real}
        return new(response, type, mirror_symmetrize(k), t)
    end
    function BareInteractionId(response::Response, type::AnalyticProperty=Instant; k, t=(0, 0))
        return new(response, type, mirror_symmetrize(k), Tuple(t))
    end
end
Base.show(io::IO, v::BareInteractionId) = print(io, "$(short(v.response))$(short(v.type)), k$(v.extK), t$(v.extT)")
function Base.hash(v::BareInteractionId, h::UInt)
    h = hash(BareInteractionId, h)
    h = hash(v.response, h)
    h = hash(v.type, h)
    h = hash(v.extK, h)
    # h = hash(round.(v.extK, digits=8), h)

    if v.extT[1] == v.extT[2]
        # the interaction is not time-dependent, then the specific time is not important
        h = hash(:time_independent, h)
    else
        h = hash(v.extT, h)
    end

    return h
end
function Base.isequal(a::BareInteractionId, b::BareInteractionId)
    # Check if response, type, and extK are not equal
    if (a.response != b.response) || (a.type != b.type) || !isequal(a.extK, b.extK)
        return false
    end

    return ((a.extT[1] == a.extT[2]) && (b.extT[1] == b.extT[2])) || (a.extT == b.extT)
end

struct GenericId{P} <: DiagramId
    para::P
    extra::Any
    GenericId(para::P, extra=nothing) where {P} = new{P}(para, extra)
end
Base.show(io::IO, v::GenericId) = print(io, isnothing(v.extra) ? "" : "$(v.extra)")
function Base.hash(v::GenericId, h::UInt)
    h = hash(GenericId, h)
    h = hash(v.para, h)
    h = hash(v.extra, h)
    return h
end
function Base.isequal(a::GenericId, b::GenericId)
    return a.para == b.para && a.extra == b.extra
end

function mirror_symmetrize(k::Vector{T}) where {T<:Number}
    idx = findfirst(!iszero, k)
    if isnothing(idx) || k[idx] > 0
        return k
    else
        mk = -k
        if T <: Real
            for i in 1:length(mk)
                if mk[i] == -T(0)
                    mk[i] = T(0)
                end
            end
        end
        return mk
    end
end

struct GreenId{P} <: DiagramId
    para::P
    type::AnalyticProperty #Instant, Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function GreenId(para::P, type::AnalyticProperty, k::Vector{T}, t::Tuple{Int,Int}) where {P,T<:Real}
        return new{P}(para, type, mirror_symmetrize(k), t)
    end
    function GreenId(para::P, type::AnalyticProperty=Dynamic; k, t) where {P}
        return new{P}(para, type, mirror_symmetrize(k), Tuple(t))
    end
end
Base.show(io::IO, v::GreenId) = print(io, "$(short(v.type)), k$(v.extK), t$(v.extT)")
function Base.hash(v::GreenId, h::UInt)
    h = hash(GreenId, h)
    h = hash(v.para, h)
    h = hash(v.type, h)
    h = hash(v.extK, h)
    h = hash(v.extT, h)
    return h
end
function Base.isequal(a::GreenId, b::GreenId)
    return a.type == b.type && a.extT == b.extT && isequal(a.extK, b.extK) && a.para == b.para
end

struct VacuumId{P} <: DiagramId
    para::P
    function VacuumId(para::P) where {P}
        return new{P}(para)
    end
end
Base.show(io::IO, v::VacuumId) = print(io, "vacuum")
function Base.hash(v::VacuumId, h::UInt)
    h = hash(VacuumId, h)
    h = hash(v.para, h)
    return h
end
function Base.isequal(a::VacuumId, b::VacuumId)
    return a.para == b.para
end

struct SigmaId{P} <: DiagramId
    para::P
    type::AnalyticProperty #Instant, Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function SigmaId(para::P, type::AnalyticProperty, k::Vector{T}, t::Tuple{Int,Int}) where {P,T<:Real}
        return new{P}(para, type, mirror_symmetrize(k), t)
    end
    function SigmaId(para::P, type::AnalyticProperty; k, t=(0, 0)) where {P}
        return new{P}(para, type, mirror_symmetrize(k), Tuple(t))
    end
end
Base.show(io::IO, v::SigmaId) = print(io, "$(short(v.type)), k$(v.extK), t$(v.extT)")
function Base.hash(v::SigmaId, h::UInt)
    h = hash(SigmaId, h)
    h = hash(v.para, h)
    h = hash(v.type, h)
    h = hash(v.extK, h)
    h = hash(v.extT, h)
    return h
end
function Base.isequal(a::SigmaId, b::SigmaId)
    if typeof(a) != typeof(b)
        return false
    end
    return a.type == b.type && a.extT == b.extT && isequal(a.extK, b.extK) && a.para == b.para
end

struct PolarId{P} <: DiagramId
    para::P
    response::Response #UpUp, UpDown, ...
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function PolarId(para::P, response::Response, k::Vector{T}, t::Tuple{Int,Int}) where {P,T<:Real}
        return new{P}(para, response, mirror_symmetrize(k), t)
    end
    function PolarId(para::P, response::Response; k, t=(0, 0)) where {P}
        return new{P}(para, response, mirror_symmetrize(k), Tuple(t))
    end
end
Base.show(io::IO, v::PolarId) = print(io, "$(short(v.response)), k$(v.extK), t$(v.extT)")
function Base.hash(v::PolarId, h::UInt)
    h = hash(PolarId, h)
    h = hash(v.para, h)
    h = hash(v.response, h)
    h = hash(v.extK, h)
    h = hash(v.extT, h)
    return h
end
function Base.isequal(a::PolarId, b::PolarId)
    if typeof(a) != typeof(b)
        return false
    end
    return a.response == b.response && a.extT == b.extT && isequal(a.extK, b.extK) && a.para == b.para
end

struct Ver3Id{P} <: DiagramId
    para::P
    response::Response #UpUp, UpDown, ...
    extK::Vector{Vector{Float64}}
    extT::Tuple{Int,Int,Int} #all possible extT from different interactionType
    function Ver3Id(para::P, response::Response, k::Vector{Vector{T}}, t::Tuple{Int,Int,Int}) where {P,T<:Real}
        return new{P}(para, response, k, t)
    end
    function Ver3Id(para::P, response::Response; k, t=(0, 0, 0)) where {P}
        return new{P}(para, response, k, Tuple(t))
    end
end
Base.show(io::IO, v::Ver3Id) = print(io, "$(short(v.response)),t$(v.extT)")
function Base.hash(v::Ver3Id, h::UInt)
    h = hash(Ver3Id, h)
    h = hash(v.para, h)
    h = hash(v.response, h)
    h = hash(v.extK, h)
    h = hash(v.extT, h)
    return h
end
function Base.isequal(a::Ver3Id, b::Ver3Id)
    if typeof(a) != typeof(b)
        return false
    end
    return a.response == b.response && a.extT == b.extT && isequal(a.extK, b.extK) && a.para == b.para
end

struct Ver4Id{P} <: DiagramId
    para::P
    response::Response #UpUp, UpDown, ...
    type::AnalyticProperty #Instant, Dynamic
    channel::TwoBodyChannel # particle-hole, particle-hole exchange, particle-particle, irreducible
    extK::Vector{Vector{Float64}}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
    function Ver4Id(para::P, response::Response, type::AnalyticProperty, chan::TwoBodyChannel, k::Vector{Vector{T}}, t::NTuple{4,Int}) where {P,T<:Real}
        return new{P}(para, response, type, chan, k, t)
    end
    function Ver4Id(para::P, response::Response, type::AnalyticProperty=Dynamic;
        k, t=(0, 0, 0, 0), chan::TwoBodyChannel=AnyChan) where {P}
        return new{P}(para, response, type, chan, k, Tuple(t))
    end
end
Base.show(io::IO, v::Ver4Id) = print(io, (v.channel == AnyChan ? "" : "$(v.channel) ") * "$(short(v.response))$(short(v.type)),t$(v.extT)")
function Base.hash(v::Ver4Id, h::UInt)
    h = hash(Ver4Id, h)
    h = hash(v.para, h)
    h = hash(v.response, h)
    h = hash(v.type, h)
    h = hash(v.channel, h)
    h = hash(v.extK, h)
    h = hash(v.extT, h)
    return h
end
function Base.isequal(a::Ver4Id, b::Ver4Id)
    if typeof(a) != typeof(b)
        return false
    end
    return a.response == b.response && a.type == b.type && a.channel == b.channel && a.extT == b.extT && isequal(a.extK, b.extK) && a.para == b.para
end

function vstr(r, c)
    N = length(r)
    # cstr(x) = x ? "⁺" : "⁻"
    s = ""
    for i ∈ 1:(N-1)
        s *= "$(r[i])$c"
    end
    s *= "$(r[end])$c"
    return s
end

function vcstr(r, creation)
    N = length(r)
    # cstr(x) = x ? "⁺" : "⁻"
    s = ""
    for i ∈ 1:(N-1)
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


struct OperatorId{P} <: DiagramId
    para::P
    site::Int
    creation::Bool
    orbital::Int
    extT::Int
    function OperatorId(para::P, r::Int, c::Bool, orbital::Int, t::Int) where {P}
        return new{P}(para, r, c, orbital, t)
    end
end
function Base.hash(v::OperatorId, h::UInt)
    h = hash(OperatorId, h)
    h = hash(v.para, h)
    h = hash(v.site, h)
    h = hash(v.creation, h)
    h = hash(v.orbital, h)
    h = hash(v.extT, h)
    return h
end
function Base.isequal(a::OperatorId, b::OperatorId)
    if typeof(a) != typeof(b)
        return false
    end
    return a.site == b.site && a.creation == b.creation && a.orbital == b.orbital && a.extT == b.extT && a.para == b.para
end


"""
hopping function c⁺c⁻
"""
struct BareHoppingId{P} <: PropagatorId
    para::P
    site::Tuple{Int,Int}
    orbital::Tuple{Int,Int}
    extT::Tuple{Int,Int}
    function BareHoppingId(para::P, r::Tuple{Int,Int}, orbital::Tuple{Int,Int}, t::Tuple{Int,Int}) where {P}
        return new{P}(para, r, orbital, t)
    end
end
Base.show(io::IO, v::BareHoppingId) = print(io, "($(vstr(v.site, "ᵣ"))|$(vstr(v.orbital, "ₒ"))|$(vcstr(v.extT, [true, false])))")
function Base.hash(v::BareHoppingId, h::UInt)
    h = hash(BareHoppingId, h)
    h = hash(v.para, h)
    h = hash(v.site, h)
    h = hash(v.orbital, h)
    h = hash(v.extT, h)
    return h
end
function Base.isequal(a::BareHoppingId, b::BareHoppingId)
    if typeof(a) != typeof(b)
        return false
    end
    return a.site == b.site && a.orbital == b.orbital && a.extT == b.extT && a.para == b.para
end

struct DetHoppingId{P} <: DiagramId
    para::P
    site::Vector{Int}
    extT::Vector{Int}
    orbital::Int
    N::Int
    function DetHoppingId(para::P, r::Vector{Int}, orbital::Int, t::Vector{Int}, N=length(r)) where {P}
        @assert length(r) == length(t) == N
        return new{P}(para, r, t, orbital, N)
    end
    function DetHoppingId(para::P; orbital=1, t=[], r=[]) where {P}
        @assert length(t) == length(r)
        return new{P}(para, r, t, orbital, length(r))
    end
end
function Base.hash(v::DetHoppingId, h::UInt)
    h = hash(DetHoppingId, h)
    # h = hash(v.para, h)
    h = hash(v.site, h)
    # h = hash(v.extT, h)
    h = hash(v.orbital, h)
    h = hash(v.N, h)
    return h
end
function Base.isequal(a::DetHoppingId, b::DetHoppingId)
    if typeof(a) != typeof(b)
        return false
    end
    # return a.N == b.N && a.orbital == b.orbital && a.site == b.site && a.extT == b.extT && a.para == b.para
    return a.site == b.site && a.orbital == b.orbital
end

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
    function BareGreenNId(para::P, r::Int, creation::Vector{Bool}, orbital::Vector{Int}, t::Vector{Int}, N::Int=length(orbital)) where {P}
        @assert length(orbital) == length(t) == length(creation) == N
        return new{P}(para, r, creation, orbital, t, N)
    end
    function BareGreenNId(para::P; orbital=[], t=[], creation=[], r=0) where {P}
        @assert length(orbital) == length(t) == length(creation)
        return new{P}(para, r, creation, orbital, t, length(orbital))
    end
end
Base.show(io::IO, v::BareGreenNId) = print(io, "($(v.site)ᵣ|$(vstr(v.orbital, "ₒ"))|$(vcstr(v.extT, v.creation)))")
function Base.hash(v::BareGreenNId, h::UInt)
    h = hash(BareGreenNId, h)
    h = hash(v.para, h)
    h = hash(v.site, h)
    h = hash(v.creation, h)
    h = hash(v.orbital, h)
    h = hash(v.extT, h)
    h = hash(v.N, h)
    return h
end
function Base.isequal(a::BareGreenNId, b::BareGreenNId)
    if typeof(a) != typeof(b)
        return false
    end
    return a.N == b.N && a.site == b.site && a.creation == b.creation && a.orbital == b.orbital && a.extT == b.extT && a.para == b.para
end

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
    function GreenNId(para::P, r::Vector{Int}, creation::Vector{Bool}, orbital::Vector{Int}, t::Vector{Int}, N::Int=length(orbital)) where {P}
        @assert N > 0
        @assert length(orbital) == length(t) == length(r) == length(creation) == N
        return new{P}(para, r, creation, orbital, t, N)
    end
    function GreenNId(para::P; orbital=[], t=[], creation=[], r=[]) where {P}
        @assert length(orbital) > 0
        @assert length(orbital) == length(t) == length(r) == length(creation)
        return new{P}(para, r, creation, orbital, t, length(orbital))
    end
end
Base.show(io::IO, v::GreenNId) = print(io, "($(vstr(v.site, "ᵣ"))|$(vstr(v.orbital, "ₒ"))|$(vcstr(v.extT, v.creation)))")
function Base.hash(v::GreenNId, h::UInt)
    h = hash(GreenNId, h)
    h = hash(v.para, h)
    h = hash(v.site, h)
    h = hash(v.creation, h)
    h = hash(v.orbital, h)
    h = hash(v.extT, h)
    h = hash(v.N, h)
    return h
end
function Base.isequal(a::GreenNId, b::GreenNId)
    if typeof(a) != typeof(b)
        return false
    end
    return a.N == b.N && a.site == b.site && a.creation == b.creation && a.orbital == b.orbital && a.extT == b.extT && a.para == b.para
end

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
    function ConnectedGreenNId(para::P, r::Vector{Int}, creation::Vector{Bool}, orbital::Vector{Int}, t::Vector{Int}, N::Int=length(orbital)) where {P}
        @assert N > 0
        @assert length(orbital) == length(t) == length(r) == length(creation) == N
        return new{P}(para, r, creation, orbital, t, N)
    end
    function ConnectedGreenNId(para::P; orbital=[], t=[], creation=[], r=[]) where {P}
        @assert length(orbital) > 0
        @assert length(orbital) == length(t) == length(r) == length(creation)
        return new{P}(para, r, creation, orbital, t, length(orbital))
    end
end
Base.show(io::IO, v::ConnectedGreenNId) = print(io, "($(vstr(v.site, "ᵣ"))|$(vstr(v.orbital, "ₒ"))|$(vcstr(v.extT, v.creation)))")
function Base.hash(v::ConnectedGreenNId, h::UInt)
    h = hash(ConnectedGreenNId, h)
    h = hash(v.para, h)
    h = hash(v.site, h)
    h = hash(v.creation, h)
    h = hash(v.orbital, h)
    h = hash(v.extT, h)
    h = hash(v.N, h)
    return h
end
function Base.isequal(a::ConnectedGreenNId, b::ConnectedGreenNId)
    if typeof(a) != typeof(b)
        return false
    end
    return a.N == b.N && a.site == b.site && a.creation == b.creation && a.orbital == b.orbital && a.extT == b.extT && a.para == b.para
end

# Must define Base.hash and Base.isequal for new DiagramId types!
# function Base.isequal(a::DiagramId, b::DiagramId)
#     if typeof(a) != typeof(b)
#         return false
#     end
#     for field in fieldnames(typeof(a))
#         if getproperty(a, field) != getproperty(b, field)
#             return false
#         end
#     end
#     return true
# end

function index(type)
    if type == BareGreenId
        return 1
    elseif type == BareInteractionId
        return 2
    elseif type <: BareGreenNId
        return 3
    elseif type <: BareHoppingId
        return 4
    elseif type <: GreenNId
        return 5
    elseif type <: DetHoppingId
        return 6
    else
        # error("Not Implemented!")
        return 0
    end
end

"""
	reconstruct(instance::DiagramId, updates::Pair{Symbol}...)

Create a new instance of the same type as `instance`, with specified fields updated to new values.

# Usage
new_instance = reconstruct(old_instance, :field1 => new_value1, :field2 => new_value2)
"""
function reconstruct(instance::DiagramId, updates::Pair{Symbol}...)
    # Get the type of the instance
    T = typeof(instance)

    # Extract field names and values from the instance
    field_names = fieldnames(T)
    field_values = [getfield(instance, fn) for fn in field_names]

    # Update fields based on the updates provided
    for (field, new_value) in updates
        field_idx = findfirst(==(field), field_names)
        if field_idx !== nothing
            field_values[field_idx] = new_value
        else
            throw(ArgumentError("Field $field does not exist in type $T"))
        end
    end

    # Construct a new instance with the updated field values
    return Base.typename(T).wrapper(field_values...)
end
