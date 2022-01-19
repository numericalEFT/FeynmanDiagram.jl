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
    uid::Int
    name::Symbol
    GenericId(name::Symbol = :none) = new(uid(), name)
end
Base.show(io::IO, v::GenericId) = print(io, "$(v.name == :none ? "" : v.name)")

struct GreenId <: DiagramId
    uid::Int
    name::Symbol
    para::GenericPara
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function GreenId(para::GenericPara, type::AnalyticProperty = Dynamic; k, t, name = :none)
        return new(uid(), name, para, type, k, Tuple(t))
    end
end
Base.show(io::IO, v::GreenId) = print(io, "G$(v.name == :none ? "" : v.name) $(v.type), k$(v.extK), t$(v.extT)")

struct InteractionId <: DiagramId
    uid::Int
    name::Symbol
    para::GenericPara
    response::Response #UpUp, UpDown, ...
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    function InteractionId(para::GenericPara, response::Response, type::AnalyticProperty = Instant; k, t = (0, 0), name = :none)
        return new(uid(), name, para, response, type, k, Tuple(t))
    end
end
Base.show(io::IO, v::InteractionId) = print(io, "W$(v.name == :none ? "" : v.name) $(short(v.response))$(short(v.type)), k$(v.extK), t$(v.extT)")

struct SigmaId <: DiagramId
    uid::Int
    name::Symbol
    para::GenericPara
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
end
Base.show(io::IO, v::SigmaId) = print(io, "Σ$(v.name == :none ? "" : v.name) $(short(v.type)), k$(v.extK), t$(v.extT)")

struct PolarId <: DiagramId
    uid::Int
    name::Symbol
    response::Response #ChargeCharge, SpinSpin, ...
    extK::Vector{Float64}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
    para::GenericPara
end
Base.show(io::IO, v::PolarId) = print(io, "Π$(v.name == :none ? "" : v.name) $(short(v.response))$(short(v.type)), k$(v.extK), t$(v.extT)")

struct Ver3Id <: DiagramId
    uid::Int
    name::Symbol
    para::GenericPara
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Vector{Float64}}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
end
Base.show(io::IO, v::Ver3Id) = print(io, "Γ3$(v.name == :none ? "" : v.name) $(short(v.type)), t$(v.extT)")

struct Ver4Id <: DiagramId
    uid::Int
    name::Symbol
    para::GenericPara
    response::Response #UpUp, UpDown, ...
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    DiEx::Int # 1 for direct, 2 for exchange
    extK::Vector{Vector{Float64}}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
end
Base.show(io::IO, v::Ver4Id) = print(io, "Γ4$(v.name == :none ? "" : v.name) $(short(v.response))$(short(v.type))$(v.DiEx == DI ? :D : (v.DiEx == EX ? :E : :DE)),t$(v.extT)")


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

function Base.Dict(v::DiagramId)
    d = Dict{Symbol,Any}()
    for field in fieldnames(typeof(v))
        if field == :extT
            tidx = getproperty(v, :extT)
            if length(tidx) == 2 # for sigma, polar
                d[:Tin], d[:Tout] = tidx[1], tidx[2]
            elseif length(tidx) == 3 # vertex3
                d[:TinL], t[:ToutL], t[:TinR] = tidx[INL], tidx[OUTL], tidx[INR]
            elseif length(tidx) == 4 # vertex4
                d[:TinL], t[:ToutL], t[:TinR], t[:ToutR] = tidx[INL], tidx[OUTL], tidx[INR], tidx[OUTR]
            else
                error("not implemented!")
            end
        elseif field == :extK
            k = getproperty(v, :extK)
            if length(k) == 2 # for sigma, polar
                d[:Kin], d[:Kout] = k[1], k[2]
            elseif length(k) == 3 # vertex3
                d[:KinL], d[:KoutL], d[:KinR] = k[INL], k[OUTL], k[INR]
            elseif length(k) == 4 # vertex4
                d[:KinL], d[:KoutL], d[:KinR], d[:KoutR] = k[INL], k[OUTL], k[INR], k[OUTR]
            else
                error("not implemented!")
            end
        else
            d[field] = getproperty(v, field)
        end
    end
    return d
end


