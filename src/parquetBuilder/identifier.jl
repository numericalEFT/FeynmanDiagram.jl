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
Base.show(io::IO, v::Vertex4) = print(io, "$(symbol(v.name, v.type))_$(v.DiEx ==1 ? :di : :ex), T$(v.extT)")

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

mutable struct Node{I<:Identifier}
    id::I
    node::Component
    children::Vector{Component}
    # function Node(id::I, node, children) where {I}
    #     return new{I}(id, node, collect(children))
    # end
    function Node(id::I; node = zero(Component), children = []) where {I}
        return new{I}(id, node, collect(children))
    end
    # function Node(id::I; node = zero(Component), children::Component) where {I}
    #     return new{I}(id, node, [children,])
    # end
end

Base.:(==)(a::Identifier, b::Identifier) = Base.isequal(a, b)

function add!(nodesVec::Vector{Node{I}}, newId::I, children, compare::Function = Base.isequal) where {I<:Identifier}
    for (ni, n) in enumerate(nodesVec)
        if compare(n.id, newId)
            append!(n.children, children)
            return ni
        end
    end
    push!(nodesVec, Node(newId, children = children))
    return length(nodesVec)
end

# function merge(nodesVec::Vector{Node{I}}, compare::Function = Base.isequal) where {I<:Identifier}
#     merged = Vector{Node{I}}([])
#     for n in nodesVec
#         add!(merged, n.id, n.children, compare)
#     end
#     return merged
# end

function classify(nodesVec::Vector{Node{I}}, comparedSyms::Symbol...) where {I<:Identifier}
    function compare(id1, id2)
        for s in comparedSyms
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

    mergedNodes = Vector{Node{I}}([])
    for n in nodesVec
        add!(mergedNodes, n.id, n.children, compare)
    end


    group = []
    for mn in mergedNodes
        key = Tuple(Base.getproperty(mn.id, sym) for sym in comparedSyms)
        push!(group, (key, mn.children))
    end

    # group[]
    # group = [mn.children for mn in mergedNodes]
    return group
end

function classify!(diag::DiagTree.Diagrams, nodesVec::Vector{Node{I}}, comparedSyms::Symbol...) where {I<:Identifier}
    group = classify(nodesVec, comparedSyms...)
    # println(group)
    componentgroup = []
    name = reduce(*, sym for sym in comparedSyms)
    for g in group
        component = DiagTree.addnode!(diag, ADD, name, g[2], para = g[1])
        push!(componentgroup, (g[1], component))
    end
    # println("merge", componentgroup)
    return componentgroup
end

# function merge(diag, nodesVec::Vector{Node{I}}, comparedSyms::Symbol...) where {I<:Identifier}
#     # if one of the comparedSyms is different, two objects are different 

#     group = classify(nodesVec, comparedSyms...)

#     merged = []
#     for g in group
#     end
#     return merged
#     # return merge(nodesVec, compare)
# end