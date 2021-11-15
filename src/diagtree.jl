module DiagTree
using ..Var

abstract type AbstractWeight end

struct PropagatorKT{Weight} <: AbstractWeight
    id::Int
    type::Int
    order::Int
    version::Int
    Kidx::Int
    Tidx::Tuple{Int,Int}
    weight::Weight
end

# ##### pretty print of Bubble and Ver4  ##########################
# Base.show(io::IO, bub::Bubble) = AbstractTrees.printnode(io::IO, bub)
# Base.show(io::IO, ver4::Ver4) = AbstractTrees.printnode(io::IO, ver4)

# ################## implement AbstractTrees interface #######################
# # refer to https://github.com/JuliaCollections/AbstractTrees.jl for more details
# function AbstractTrees.children(ver4::Ver4)
#     return ver4.bubble
# end

# function AbstractTrees.children(bubble::Bubble)
#     return (bubble.Lver, bubble.Rver)
# end

# function iterate(ver4::Ver4{W}) where {W}
#     if length(ver4.bubble) == 0
#         return nothing
#     else
#         return (ver4.bubble[1], 1)
#     end
# end

# function iterate(bub::Bubble{Ver4{W},W}) where {W}
#     return (bub.Lver, false)
# end

# function iterate(ver4::Ver4{W}, state) where {W}
#     if state >= length(ver4.bubble) || length(ver4.bubble) == 0
#         return nothing
#     else
#         return (ver4.bubble[state+1], state + 1)
#     end
# end

# function iterate(bub::Bubble{Ver4{W},W}, state::Bool) where {W}
#     state && return nothing
#     return (bub.Rver, true)
# end

# Base.IteratorSize(::Type{Ver4{W}}) where {W} = Base.SizeUnknown()
# Base.eltype(::Type{Ver4{W}}) where {W} = Ver4{W}

# Base.IteratorSize(::Type{Bubble{Ver4{W},W}}) where {W} = Base.SizeUnknown()
# Base.eltype(::Type{Bubble{Ver4{W},W}}) where {W} = Bubble{Ver4{W},W}

# AbstractTrees.printnode(io::IO, ver4::Ver4) = print(io, tpair(ver4))
# AbstractTrees.printnode(io::IO, bub::Bubble) = print(io, "\u001b[32m$(bub.id): $(ChanName[bub.chan]) $(bub.Lver.loopNum)Ⓧ $(bub.Rver.loopNum)\u001b[0m")
end