"""
    mutable struct TaylorSeries{T}
    
    A representation of a taylor series. 

# Members:
- `name::Symbol`  name of the diagram
- `coeffs::Dict{Vector{Int},T}`  The taylor expansion coefficients. The integer array define the order of corresponding coefficient. 
"""
mutable struct TaylorSeries{T}
    name::String # "" by default
    coeffs::Dict{Vector{Int},T}

    """
        function TaylorSeries{T}(coeffs::Dict{Vector{Int},T}=Dict{Vector{Int},T}(), name::String="") where {T}
            Create a TaylorSeries based on given coefficients.
    """
    function TaylorSeries{T}(coeffs::Dict{Vector{Int},T}=Dict{Vector{Int},T}(), name::String="") where {T}
        return new{T}(name, coeffs)
    end
end

"""
    function TaylorSeries(::Type{T}, nv::Int) where {T}
    
    Create a taylor series equal to variable with index nv. For example, if global variables are "x y", in put nv=2 generate series t=y.

# Arguments:
- `::Type{T}`  DataType of coefficients in taylor series.
- `nv::Int`  Index of variable. 
"""
function TaylorSeries(::Type{T}, nv::Int) where {T}
    @assert 0 < nv ≤ get_numvars()
    v = zeros(Int, get_numvars())
    @inbounds v[nv] = 1
    return TaylorSeries{T}(Dict{Vector{Int},T}(v => one(T)))
end
TaylorSeries(nv::Int) = TaylorSeries(Float64, nv)


