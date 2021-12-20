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

# function refection(Type, D::Int)
#     return (-diagm(ones(D)), zeros(Type, D))
# end

# function particleHole(::type, N::Int, β)
#     return Tuple{Matrix{type},Matrix{type}}(-Diagonal(ones(D)), zeros(type, (D, D)))
# end

struct VectorVariable{PARA,T}
    para::PARA
    basis::Vector{T}
    function VectorVariable(_basis, para::P = 0) where {P}
        return new{P,eltype(_basis)}(para, _basis)
    end
    function VectorVariable{P,T}(_basis, para = 0) where {P,T}
        return new{P,T}(P(para), T.(_basis))
    end
end
"""
Base.:(==)(a::VectorVariable, b::VectorVariable) = Base.isequal(a, b)
function Base.isequal(a::VectorVariable{T}, b::VectorVariable{T}) where {T}

Compare two vectors with the function isequal or the operator ==
"""
Base.:(==)(a::VectorVariable, b::VectorVariable) = Base.isequal(a, b)
function Base.isequal(a::VectorVariable{P,T}, b::VectorVariable{P,T}) where {P,T}
    if applicable(isequal, a, b)
        return isequal(a.basis, b.basis)
    elseif applicable(isapprox, a, b)
        return a.basis ≈ b.basis
    else
        error("Comparison between ScalarVariable $a and $b has not yet been implemented!")
    end
end

struct ScalarVariable{PARA,T}
    para::PARA
    basis::T
    function ScalarVariable(_basis, para::P = 0) where {P}
        return new{P,typeof(_basis)}(para, _basis)
    end
    function ScalarVariable{P,T}(_basis, para = 0) where {P,T}
        return new{P,T}(P(para), T(_basis))
    end
end
"""
Base.:(==)(a::VectorVariable, b::VectorVariable) = Base.isequal(a, b)
function Base.isequal(a::VectorVariable{T}, b::VectorVariable{T}) where {T}

Compare two vectors with the function isequal or the operator ==
"""
Base.:(==)(a::ScalarVariable, b::ScalarVariable) = Base.isequal(a, b)
function Base.isequal(a::ScalarVariable{P,T}, b::ScalarVariable{P,T}) where {P,T}
    if applicable(isequal, a, b)
        return isequal(a.basis, b.basis)
    elseif applicable(isapprox, a, b)
        return a.basis ≈ b.basis
    else
        error("Comparison between ScalarVariable $a and $b has not yet been implemented!")
    end
end

end