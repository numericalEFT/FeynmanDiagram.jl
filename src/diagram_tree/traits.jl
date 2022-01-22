abstract type Operator end
struct Sum <: Operator end
struct Prod <: Operator end
Base.isequal(a::Operator, b::Operator) = (typeof(a) == typeof(b))
Base.:(==)(a::Operator, b::Operator) = Base.isequal(a, b)
apply(o::Operator, diags) = error("not implemented!")

Base.show(io::IO, o::Operator) = print(io, typeof(o))
Base.show(io::IO, o::Sum) = print(io, "⨁")
Base.show(io::IO, o::Prod) = print(io, "Ⓧ")


abstract type DiagramId end

# Base.Dict(x::DiagramId) = Dict{Symbol,Any}([fn => getfield(x, fn) for fn ∈ fieldnames(typeof(x))])
# Base.show(io::IO, d::DiagramId) = error("Base.show not implemented!")
# Base.isequal(a::DiagramId, b::DiagramId) = error("Base.isequal not implemented!")
# eval(d::DiagramId) = error("eval for $d has not yet implemented!")
Base.:(==)(a::DiagramId, b::DiagramId) = Base.isequal(a, b)

# function summary(d::DiagramId, verbose::Int = 0, color = false)
#     for field in fieldnames(typeof(d))
#     end
# end

struct GenericId <: DiagramId
    para::GenericPara
    extra::Any
    GenericId(para::GenericPara, extra = Nothing) = new(para, extra)
end
Base.show(io::IO, v::GenericId) = print(io, v.extra == Nothing ? "" : "$(v.extra)")

struct GreenId <: DiagramId
    para::GenericPara
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function GreenId(para::GenericPara, type::AnalyticProperty = Dynamic; k, t)
        return new(para, type, k, Tuple(t))
    end
end
Base.show(io::IO, v::GreenId) = print(io, "$(short(v.type)), k$(v.extK), t$(v.extT)")

struct InteractionId <: DiagramId
    para::GenericPara
    response::Response #UpUp, UpDown, ...
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    permutation::Permutation
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function InteractionId(para::GenericPara, response::Response, type::AnalyticProperty = Instant; k, t = (0, 0), permu::Permutation = DiEx)
        return new(para, response, type, permu, k, Tuple(t))
    end
end
Base.show(io::IO, v::InteractionId) = print(io, "$(short(v.response))$(short(v.type))$(v.permutation), k$(v.extK), t$(v.extT)")

struct SigmaId <: DiagramId
    para::GenericPara
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function SigmaId(para::GenericPara, type::AnalyticProperty; k, t = (0, 0))
        return new(para, type, k, t)
    end
end
Base.show(io::IO, v::SigmaId) = print(io, "$(short(v.type)), t$(v.extT)")

struct PolarId <: DiagramId
    para::GenericPara
    response::Response #ChargeCharge, SpinSpin, ...
    extK::Vector{Float64}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
    function PolarId(para::GenericPara, response::Response; k, t = (0, 0))
        return new(para, response, k, t)
    end
end
Base.show(io::IO, v::PolarId) = print(io, "$(short(v.response)), k$(v.extK), t$(v.extT)")

struct Ver3Id <: DiagramId
    para::GenericPara
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Vector{Float64}}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
end
Base.show(io::IO, v::Ver3Id) = print(io, "$(short(v.type)), t$(v.extT)")

struct Ver4Id <: DiagramId
    para::GenericPara
    response::Response #UpUp, UpDown, ...
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    channel::TwoBodyChannel # particle-hole (T), particle-hole exchange (U), particle-particle (S), irreducible (I)
    extK::Vector{Vector{Float64}}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
    function Ver4Id(para::GenericPara, response::Response, type::AnalyticProperty = Dynamic; k, t = (0, 0, 0, 0), chan::TwoBodyChannel = AnyChan)
        return new(para, response, type, chan, k, Tuple(t))
    end
end
Base.show(io::IO, v::Ver4Id) = print(io, (v.channel == AnyChan ? "" : "$(v.channel) ") * "$(short(v.response))$(short(v.type)),t$(v.extT)")


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


