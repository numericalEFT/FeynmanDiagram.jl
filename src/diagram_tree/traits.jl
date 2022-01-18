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
eval(d::DiagramId) = error("eval for $d has not yet implemented!")
Base.:(==)(a::DiagramId, b::DiagramId) = Base.isequal(a, b)

struct Vertex4 <: DiagramId
    ######### properties that defines a unique ver4 ###################
    response::Response #UpUp, UpDown, ...
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    DiEx::Int # 1 for direct, 2 for exchange
    extK::Vector{Vector{Float64}}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
    para::GenericPara
end
Base.show(io::IO, v::Vertex4) = print(io, "Γ4 $(v.response)$(v.type)$(v.DiEx == DI ? :D : (v.DiEx == EX ? :E : :DE)),t$(v.extT)")

struct Green <: DiagramId
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    para::GenericPara
end
Base.show(io::IO, v::Green) = print(io, "G $(v.type), k$(v.extK), t$(v.extT)")

struct BareInteraction <: DiagramId
    response::Response #UpUp, UpDown, ...
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int} #all possible extT from different interactionType
    para::GenericPara
end
Base.show(io::IO, v::BareInteraction) = print(io, "W $(v.response)$(v.type), k$(v.extK), t$(v.extT)")

struct Sigma <: DiagramId
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
    para::GenericPara
end
Base.show(io::IO, v::Sigma) = print(io, "Σ $(v.type), k$(v.extK), t$(v.extT)")

struct Vertex3 <: DiagramId
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Vector{Float64}}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
    para::GenericPara
end
Base.show(io::IO, v::Vertex3) = print(io, "Γ3 $(v.type), t$(v.extT)")

struct Polar <: DiagramId
    response::Response #ChargeCharge, SpinSpin, ...
    extK::Vector{Float64}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
    para::GenericPara
end
Base.show(io::IO, v::Polar) = print(io, "Π $(v.response)$(v.type), k$(v.extK), t$(v.extT)")

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


