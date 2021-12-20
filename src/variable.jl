module Var
using StaticArrays
using LinearAlgebra

# const SymOperator1d = Tuple{Float64,Float64}
# const SymOperatorNd = Tuple{Float64,Float64}

# mutable struct Momentum{D,T}
#     symmetry::Symbol
#     basis::Vector{T}
#     function Momentum{D,T}(_basis, _symmetry = :none) where {D,T}
#         return new(_symmetry, _basis)
#     end
# end

# mutable struct TauPair{N}
#     symmetry::Vector{Tuple{Float64,Float64}}
#     basis::Vector{}
# end

function refection(Type, D::Int)
    return (-diagm(ones(D)), zeros(Type, (D, D)))
end

# function particleHole(::type, N::Int, β)
#     return Tuple{Matrix{type},Matrix{type}}(-Diagonal(ones(D)), zeros(type, (D, D)))
# end

struct VectorVariable{T}
    symmetry::Vector{Tuple{Matrix{T},Matrix{T}}}
    basis::Vector{T}
    function VectorVariable(_basis, _symmetry = [])
        for sym in _symmetry
            @assert size(sym[1])[1] == length(_basis) "rotation operator size $(size(sym[1])) doesn't match with the basis size $(length(_basis))"
            @assert size(sym[1])[2] == length(_basis) "rotation operator size $(size(sym[1])) doesn't match with the basis size $(length(_basis))"
            @assert size(sym[2])[1] == length(_basis) "translation operator size $(size(sym[2])) doesn't match with the basis size $(length(_basis))"
            @assert size(sym[2])[2] == length(_basis) "translation operator size $(size(sym[2])) doesn't match with the basis size $(length(_basis))"
        end
        return new{eltype(_basis)}(_symmetry, _basis)
    end
end

function Base.isequal(a::VectorVariable{T}, b::VectorVariable{T}) where {T}
    if a ≈ b
        return true
    end

    for sym in symmetry
        rotation, translation = sym[1], sym[2]
        if rotation * a + translation ≈ b
            return true
        end
    end
    return false
end

end