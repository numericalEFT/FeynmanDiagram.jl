mutable struct Element{T}
    extT::Tuple{Int,Int,Int,Int}
    node::T
end
Base.zero(::Type{Element{T}}) where {T} = Element{T}((0, 0, 0, 0), T(0))

mutable struct Ver4Identity{name<:InteractionName}
    ######### properties that defines a unique ver4 ###################
    # name::InteractionName #composite, chargecharge, spinspin, ...
    legK::Vector{Vector{Float64}}
    extT::Set{Tuple{Int,Int,Int,Int}} #all possible extT from different interactionType
end

function Base.isequal(a::Ver4Identity, b::Ver4Identity)
    # only parent is allowed to be different
    if (typeof(a) != typeof(b)) || (a.legK != b.legK) || (a.extT != b.extT)
        return false
    else
        return true
    end
end
Base.:(==)(a::Ver4Identity, b::Ver4Identity) = Base.isequal(a, b)

mutable struct Ver4Node{T}
    identity::Any

    ######## possible interactionType of the ver4 ############################# 
    instant::Tuple{Element{T},Element{T}}
    dynamic::Tuple{Element{T},Element{T}}
    d_instant::Tuple{Element{T},Element{T}} #derivative
    d_dynamic::Tuple{Element{T},Element{T}} #derivative
    function Ver4Info{T}(name::InteractionName, legK = [[],], extT = []) where {T}
        identity = Ver4Identity{name}(legK, extT)
        info = new{T}(identity, zero(Element{T}), zero(Element{T}), zero(Element{T}), zero(Element{T}))
        updateIdentity(info)
        return info
    end
end

function Base.isequal(a::Ver4Node{T}, b::Ver4Node{T}) where {T}
    # only parent is allowed to be different
    if a.identity != b.identity
        return false
    else
        return true
    end
end
Base.:(==)(a::Ver4Node{T}, b::Ver4Node{T}) where {T} = Base.isequal(a, b)

function updateIdentity(node::Ver4Node)
    node.extT = Set([])

    function push(element)
        if element[1].extT != (0, 0, 0, 0)
            push!(node.extT, element[1].extT)
        end
        if element[2].extT != (0, 0, 0, 0)
            push!(node.extT, element[2].extT)
        end
    end
    push(node.instant)
    push(node.dynamic)
    push(node.d_instant)
    push(node.d_dynamic)
end

