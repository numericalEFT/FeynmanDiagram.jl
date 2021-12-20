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
    return (-diagm(ones(D)), zeros(Type, D))
end

# function particleHole(::type, N::Int, β)
#     return Tuple{Matrix{type},Matrix{type}}(-Diagonal(ones(D)), zeros(type, (D, D)))
# end

struct VectorVariable{T}
    symmetry::Vector{Tuple{Matrix{T},Vector{T}}}
    basis::Vector{T}
    function VectorVariable(_basis, _symmetry = [])
        for sym in _symmetry
            @assert size(sym[1])[1] == length(_basis) "rotation operator size $(size(sym[1])) doesn't match with the basis size $(length(_basis))"
            @assert size(sym[1])[2] == length(_basis) "rotation operator size $(size(sym[1])) doesn't match with the basis size $(length(_basis))"
            @assert length(sym[2]) == length(_basis) "translation operator size $(length(sym[2])) doesn't match with the basis size $(length(_basis))"
        end
        return new{eltype(_basis)}(_symmetry, _basis)
    end
    function VectorVariable{T}(_basis, _symmetry = []) where {T}
        for sym in _symmetry
            @assert size(sym[1])[1] == length(_basis) "rotation operator size $(size(sym[1])) doesn't match with the basis size $(length(_basis))"
            @assert size(sym[1])[2] == length(_basis) "rotation operator size $(size(sym[1])) doesn't match with the basis size $(length(_basis))"
            @assert length(sym[2]) == length(_basis) "translation operator size $(length(sym[2])) doesn't match with the basis size $(length(_basis))"
        end
        return new{T}(T.(_symmetry), T.(_basis))
    end
end
"""
Base.:(==)(a::VectorVariable, b::VectorVariable) = Base.isequal(a, b)
function Base.isequal(a::VectorVariable{T}, b::VectorVariable{T}) where {T}

Compare two vectors with the function isequal or the operator ==
"""
Base.:(==)(a::VectorVariable, b::VectorVariable) = Base.isequal(a, b)
function Base.isequal(a::VectorVariable{T}, b::VectorVariable{T}) where {T}
    if length(a.symmetry) != length(b.symmetry)
        return false
    end

    for (si, sym) in enumerate(a.symmetry)
        if !(sym[1] ≈ b.symmetry[si][1]) || !(sym[2] ≈ b.symmetry[si][2])
            return false
        end
    end

    if a.basis ≈ b.basis
        return true
    end

    for sym in a.symmetry

        rotation, translation = sym[1], sym[2]
        if rotation * a.basis .+ translation ≈ b.basis
            return true
        end
    end
    return false
end

struct ScalarVariable{T}
    symmetry::Vector{Tuple{T,T}}
    basis::T
    function ScalarVariable(_basis, _symmetry = [])
        return new{typeof(_basis)}(_symmetry, _basis)
    end
end
"""
Base.:(==)(a::VectorVariable, b::VectorVariable) = Base.isequal(a, b)
function Base.isequal(a::VectorVariable{T}, b::VectorVariable{T}) where {T}

Compare two vectors with the function isequal or the operator ==
"""
Base.:(==)(a::ScalarVariable, b::ScalarVariable) = Base.isequal(a, b)
function Base.isequal(a::ScalarVariable{T}, b::ScalarVariable{T}) where {T}
    if length(a.symmetry) != length(b.symmetry)
        return false
    end

    for (si, sym) in enumerate(a.symmetry)
        if !(sym[1] ≈ b.symmetry[si][1]) || !(sym[2] ≈ b.symmetry[si][2])
            return false
        end
    end

    if a.basis ≈ b.basis
        return true
    end

    for sym in a.symmetry

        rotation, translation = sym[1], sym[2]
        if rotation * a.basis .+ translation ≈ b.basis
            return true
        end
    end
    return false
end

end