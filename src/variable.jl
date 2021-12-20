module Var
using StaticArrays

abstract type CachedVariable end

mutable struct Momentum{D,T} <: CachedVariable
    symmetry::Symbol
    basis::Vector{T}
    function Momentum{D,T}(_basis, _symmetry = :none) where {D,T}
        return new(_symmetry, _basis)
    end
end

# mutable struct TauPair{N<:Int}

# const CachedMomentum = CachedObject{Momentum,MVector}

# macro 
end