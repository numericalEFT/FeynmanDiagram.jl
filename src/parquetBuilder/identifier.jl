abstract type Identifier end

struct Vertex4 <: Identifier
    ######### properties that defines a unique ver4 ###################
    name::ResponseName #ChargeCharge, SpinSpin, ...
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    DiEx::Int # 1 for direct, 2 for exchange
    extK::Vector{Vector{Float64}}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
end
function Base.isequal(a::Vertex4, b::Vertex4)
    return (a.name == b.name) && (a.type == b.type) && (a.extT == b.extT) && (a.DiEx == b.DiEx) && (a.extK ≈ b.extK)
end

struct Sigma <: Identifier
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Float64}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
end

struct Vertex3 <: Identifier
    type::AnalyticProperty #Instant, Dynamic, D_Instant, D_Dynamic
    extK::Vector{Vector{Float64}}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
end

struct Polarization <: Identifier
    name::ResponseName #ChargeCharge, SpinSpin, ...
    extK::Vector{Float64}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
end

struct Nodes{I<:Identifier}
    id::Identifier
    nodes::Vector{Component}
    function Nodes(id::I, nodes = []) where {I}
        return new{I}(id, collect(nodes))
    end
    function Nodes(id::I, nodes::Component) where {I}
        return new{I}(id, [nodes,])
    end
end

Base.:(==)(a::Identifier, b::Identifier) = Base.isequal(a, b)

function add!(nodesVec::Vector{Nodes{I}}, newId::I, nodes, compare::Function = Base.isequal) where {I<:Identifier}
    for (ni, n) in enumerate(nodesVec)
        if compare(n.id, newId)
            append!(n.nodes, nodes)
            return ni
        end
    end
    push!(nodesVec, Nodes(newId, nodes))
    return length(nodesVec)
end

function merge(nodesVec::Vector{Nodes{I}}, compare::Function = Base.isequal) where {I<:Identifier}
    merged = Vector{Nodes{I}}([])
    for n in nodesVec
        add!(merged, n.id, n.nodes, compare)
    end
    return merged
end

function merge(nodesVec::Vector{Nodes{I}}, syms::Symbol...) where {I<:Identifier}
    function compare(id1, id2)
        for s in syms
            if s != :extK
                if getproperty(id1, s) != getproperty(id2, s)
                    return false
                end
            else
                if (getproperty(id1, s) ≈ getproperty(id2, s)) == false
                    return false
                end
            end
        end
        return true
    end
    return merge(nodesVec, compare)
end