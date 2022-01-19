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
    GenericId(para::GenericPara) = new(para)
end
Base.show(io::IO, v::GenericId) = print(io, "")

struct GreenId <: DiagramId
    para::GenericPara
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function GreenId(para::GenericPara, type::AnalyticProperty = Dynamic; k, t)
        return new(para, type, k, Tuple(t))
    end
end
Base.show(io::IO, v::GreenId) = print(io, "$(v.type), k$(v.extK), t$(v.extT)")

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
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
end
Base.show(io::IO, v::SigmaId) = print(io, "$(short(v.type)), k$(v.extK), t$(v.extT)")

struct PolarId <: DiagramId
    response::Response #ChargeCharge, SpinSpin, ...
    extK::Vector{Float64}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
    para::GenericPara
end
Base.show(io::IO, v::PolarId) = print(io, "$(short(v.response))$(short(v.type)), k$(v.extK), t$(v.extT)")

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

function toDict(v::DiagramId; verbose::Int)
    d = Dict{Symbol,Any}()
    for field in fieldnames(typeof(v))
        if verbose > 1 && field == :extT
            tidx = getproperty(v, :extT)
            if length(tidx) == 2 # for sigma, polar
                d[:TinL], d[:ToutL] = tidx[1], tidx[2]
            elseif length(tidx) == 3 # vertex3
                d[:TinL], d[:ToutL], d[:TinR] = tidx[INL], tidx[OUTL], tidx[INR]
            elseif length(tidx) == 4 # vertex4
                d[:TinL], d[:ToutL], d[:TinR], d[:ToutR] = tidx[INL], tidx[OUTL], tidx[INR], tidx[OUTR]
            else
                error("not implemented!")
            end
        else
            # elseif expand && field == :extK
            #     k = getproperty(v, :extK)
            #     if length(k) == 1 # for sigma, polar
            #         d[:KinL] = k[1]
            #     elseif length(k) == 2 # for sigma, polar
            #         d[:KinL], d[:KoutL] = k[1], k[2]
            #     elseif length(k) == 3 # vertex3
            #         d[:KinL], d[:KoutL], d[:KinR] = k[INL], k[OUTL], k[INR]
            #     elseif length(k) == 4 # vertex4
            #         d[:KinL], d[:KoutL], d[:KinR], d[:KoutR] = k[INL], k[OUTL], k[INR], k[OUTR]
            #     else
            #         error("not implemented!")
            #     end
            data = getproperty(v, field)
            #DataFrame will expand a vector into multiple rows. To prevent it, we transform all vectors into tuples
            d[field] = data isa AbstractVector ? Tuple(data) : data
        end
    end
    return d
end


