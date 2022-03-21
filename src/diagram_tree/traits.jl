abstract type Operator end
struct Sum <: Operator end
struct Prod <: Operator end
Base.isequal(a::Operator, b::Operator) = (typeof(a) == typeof(b))
Base.:(==)(a::Operator, b::Operator) = Base.isequal(a, b)
apply(o::Operator, diags) = error("not implemented!")

Base.show(io::IO, o::Operator) = print(io, typeof(o))
Base.show(io::IO, o::Sum) = print(io, "⨁")
Base.show(io::IO, o::Prod) = print(io, "Ⓧ")

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
    para::GenericPara
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    order::Vector{Int}
    function BareGreenId(para::GenericPara, type::AnalyticProperty = Dynamic, order = [0, 0]; k, t)
        return new(para, type, k, Tuple(t), order)
    end
end
Base.show(io::IO, v::BareGreenId) = print(io, "$(short(v.type))#$(v.order), k$(v.extK), t$(v.extT)")

struct BareInteractionId <: PropagatorId
    para::GenericPara
    response::Response #UpUp, UpDown, ...
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    permutation::Permutation
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    order::Vector{Int}
    function BareInteractionId(para::GenericPara, response::Response, type::AnalyticProperty = Instant, order = [0, 0]; k, t = (0, 0), permu::Permutation = DiEx)
        return new(para, response, type, permu, k, Tuple(t), order)
    end
end
Base.show(io::IO, v::BareInteractionId) = print(io, "$(short(v.response))$(short(v.type))$(v.permutation)#$(v.order), k$(v.extK), t$(v.extT)")

struct GenericId <: DiagramId
    para::GenericPara
    extra::Any
    order::Vector{Int}
    GenericId(para::GenericPara, extra = Nothing, order = [0, 0]) = new(para, extra, order)
end
Base.show(io::IO, v::GenericId) = print(io, v.extra == Nothing ? "#$(v.order)" : "$(v.extra)#$(v.order)")

struct GreenId <: DiagramId
    para::GenericPara
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    order::Vector{Int}
    function GreenId(para::GenericPara, type::AnalyticProperty = Dynamic, order = [0, 0]; k, t)
        return new(para, type, k, Tuple(t), order)
    end
end
Base.show(io::IO, v::GreenId) = print(io, "$(short(v.type))#$(v.order), k$(v.extK), t$(v.extT)")

struct SigmaId <: DiagramId
    para::GenericPara
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    order::Vector{Int}
    function SigmaId(para::GenericPara, type::AnalyticProperty, order = [0, 0]; k, t = (0, 0))
        return new(para, type, k, t, order)
    end
end
Base.show(io::IO, v::SigmaId) = print(io, "$(short(v.type))#$(v.order), t$(v.extT)")

struct PolarId <: DiagramId
    para::GenericPara
    response::Response #UpUp, UpDown, ...
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    order::Vector{Int}
    function PolarId(para::GenericPara, response::Response, order = [0, 0]; k, t = (0, 0))
        return new(para, response, k, t, order)
    end
end
Base.show(io::IO, v::PolarId) = print(io, "$(short(v.response))#$(v.order), k$(v.extK), t$(v.extT)")

struct Ver3Id <: DiagramId
    para::GenericPara
    response::Response #UpUp, UpDown, ...
    extK::Vector{Vector{Float64}}
    extT::Tuple{Int,Int,Int} #all possible extT from different interactionType
    order::Vector{Int}
    function Ver3Id(para::GenericPara, response::Response, order = [0, 0]; k, t = (0, 0, 0))
        return new(para, response, k, Tuple(t), order)
    end
end
Base.show(io::IO, v::Ver3Id) = print(io, "$(short(v.response))#$(v.order),t$(v.extT)")

struct Ver4Id <: DiagramId
    para::GenericPara
    response::Response #UpUp, UpDown, ...
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    channel::TwoBodyChannel # particle-hole, particle-hole exchange, particle-particle, irreducible
    extK::Vector{Vector{Float64}}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
    order::Vector{Int}
    function Ver4Id(para::GenericPara, response::Response, type::AnalyticProperty = Dynamic, order = [0, 0]; k, t = (0, 0, 0, 0), chan::TwoBodyChannel = AnyChan)
        return new(para, response, type, chan, k, Tuple(t), order)
    end
end
Base.show(io::IO, v::Ver4Id) = print(io, (v.channel == AnyChan ? "" : "$(v.channel) ") * "$(short(v.response))$(short(v.type))#$(v.order),t$(v.extT)")

"""
time-ordered N-point Bare Green's function
"""
struct BareGreenNId <: PropagatorId
    para::GenericPara
    site::Int
    creation::Vector{Bool}
    orbital::Vector{Int}
    extT::Vector{Int}
    N::Int
    function BareGreenNId(para::GenericPara, creation = [], orbital = [], t = [], r = 0)
        @assert length(orbital) == length(t) == length(creation)
        return new(para, r, creation, orbital, t, length(orbital))
    end
end
function fieldstr(creation, orbital, extT)
    cstr(x) = x ? "⁺" : "⁻"
    N = length(creation)
    return reduce(*, ["$(extT[i])_$(orbital[i])$(cstr(creation[i]))" for i in 1:N])
end
Base.show(io::IO, v::BareGreenNId) = print(io, "gn$(v.N)($(v.site)ᵣ|$(fieldstr(v.creation, v.orbital, v.extT)))")

"""
time-ordered N-point Composite Green's function
"""
struct GreenNId <: PropagatorId
    para::GenericPara
    site::Vector{Int}
    creation::Vector{Bool}
    orbital::Vector{Int}
    extT::Vector{Int}
    N::Int
    function GreenNId(para::GenericPara, creation = [], orbital = [], t = [], r = [])
        @assert length(orbital) == length(t) == length(creation)
        return new(para, r, creation, orbital, t, length(orbital))
    end
end
Base.show(io::IO, v::GreenNId) = print(io, "Gn$(v.N)($(v.site)ᵣ|$(fieldstr(v.creation, v.orbital, v.extT)))")

"""
time-ordered N-point Composite Green's function
"""
struct ConnectedGreenNId <: PropagatorId
    para::GenericPara
    site::Vector{Int}
    creation::Vector{Bool}
    orbital::Vector{Int}
    extT::Vector{Int}
    N::Int
    function GreenNId(para::GenericPara, creation = [], orbital = [], t = [], r = [])
        @assert length(orbital) == length(t) == length(creation)
        return new(para, r, creation, orbital, t, length(orbital))
    end
end
Base.show(io::IO, v::ConnectedGreenNId) = print(io, "Gc$(v.N)($(v.site)ᵣ|$(fieldstr(v.creation, v.orbital, v.extT)))")

function Base.isequal(a::DiagramId, b::DiagramId)
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


