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
    function BareGreenId(type::AnalyticProperty, k::Vector{T}, t::Tuple{Int,Int}) where {T<:Real}
        return new(type, mirror_symmetrize(k), t)
    end
    function BareGreenId(type::AnalyticProperty=Dynamic; k, t)
        return new(type, mirror_symmetrize(k), Tuple(t))
    end
end
Base.show(io::IO, v::BareGreenId) = print(io, "$(short(v.type)), k$(v.extK), t$(v.extT)")
function Base.isequal(a::BareGreenId, b::BareGreenId)
    return a.type == b.type && a.extT == b.extT && a.extK == b.extK
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
    GenericId(para::P, extra=nothing) where {P} = new{P}(para, extra)
end
Base.show(io::IO, v::GenericId) = print(io, v.extra == Nothing ? "" : "$(v.extra)")
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
function Base.isequal(a::GreenId, b::GreenId)
    return a.type == b.type && a.extT == b.extT && a.extK == b.extK && a.para == b.para
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
Base.show(io::IO, v::SigmaId) = print(io, "$(short(v.type))#$(v.order), t$(v.extT)")
function Base.isequal(a::SigmaId, b::SigmaId)
    if typeof(a) != typeof(b)
        return false
    end
    return a.type == b.type && a.extT == b.extT && a.extK == b.extK && a.para == b.para
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
function Base.isequal(a::PolarId, b::PolarId)
    if typeof(a) != typeof(b)
        return false
    end
    return a.response == b.response && a.extT == b.extT && a.extK == b.extK && a.para == b.para
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
function Base.isequal(a::Ver3Id, b::Ver3Id)
    if typeof(a) != typeof(b)
        return false
    end
    return a.response == b.response && a.extT == b.extT && a.extK == b.extK && a.para == b.para
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
function Base.isequal(a::Ver4Id, b::Ver4Id)
    if typeof(a) != typeof(b)
        return false
    end
    return a.response == b.response && a.type == b.type && a.channel == b.channel && a.extT == b.extT && a.extK == b.extK && a.para == b.para
end

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
    function BareHoppingId(para::P, r::Tuple{Int,Int}, orbital::Tuple{Int,Int}, t::Tuple{Int,Int}) where {P}
        return new{P}(para, r, orbital, t)
    end
end
Base.show(io::IO, v::BareHoppingId) = print(io, "($(vstr(v.site, "ᵣ"))|$(vstr(v.orbital, "ₒ"))|$(vcstr(v.extT, [true, false])))")
function Base.isequal(a::BareHoppingId, b::BareHoppingId)
    if typeof(a) != typeof(b)
        return false
    end
    return a.site == b.site && a.orbital == b.orbital && a.extT == b.extT && a.para == b.para
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
        @assert length(orbital) == length(t) == length(r) == length(creation) == N
        return new{P}(para, r, creation, orbital, t, N)
    end
    function GreenNId(para::P; orbital=[], t=[], creation=[], r=[]) where {P}
        @assert length(orbital) == length(t) == length(r) == length(creation)
        return new{P}(para, r, creation, orbital, t, length(orbital))
    end
end
Base.show(io::IO, v::GreenNId) = print(io, "($(vstr(v.site, "ᵣ"))|$(vstr(v.orbital, "ₒ"))|$(vcstr(v.extT, v.creation)))")
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
        @assert length(orbital) == length(t) == length(r) == length(creation) == N
        return new{P}(para, r, creation, orbital, t, N)
    end
    function ConnectedGreenNId(para::P; orbital=[], t=[], creation=[], r=[]) where {P}
        @assert length(orbital) == length(t) == length(r) == length(creation)
        return new{P}(para, r, creation, orbital, t, length(orbital))
    end
end
Base.show(io::IO, v::ConnectedGreenNId) = print(io, "($(vstr(v.site, "ᵣ"))|$(vstr(v.orbital, "ₒ"))|$(vcstr(v.extT, v.creation)))")
function Base.isequal(a::ConnectedGreenNId, b::ConnectedGreenNId)
    if typeof(a) != typeof(b)
        return false
    end
    return a.N == b.N && a.site == b.site && a.creation == b.creation && a.orbital == b.orbital && a.extT == b.extT && a.para == b.para
end

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