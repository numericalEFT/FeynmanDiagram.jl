module Var
using StaticArrays

abstract type DerivedVariable end

# mutable struct Momentum{D<:Int,T<:Number} <: DerivedVariable
#     id::Int
#     version::Int
#     basis::Vector{Int}
#     K::SVector{D,T}
#     function Momentum{D,T}(_id, _basis, _K = @SVector rand(T, D), _version = 0) where {D<:Int,T<:Number}
#         return new(_id, _version, _basis, _K)
#     end
# end

end