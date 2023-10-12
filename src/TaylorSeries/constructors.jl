"""
    mutable struct TaylorSeries{F,W,N}
    
    A representation of a computational graph, e.g., an expression tree, with type stable node data.

# Members:
- `id::Int`  the unique hash id to identify the diagram
- `name::Symbol`  name of the diagram
- `expansion::Dict{Dict{Int,Int},T}`  The taylor expansion coefficients. The key Dict{Int,Int} labels the order with respect to each variables. 
- `variables::Set{V}`  Variables of the taylor series. Each variable must have an unique id. 
- `truncate::Dict{Int, Int}` For each variable, the taylor series is truncated to certain order. If empty, the function must be a polinomial, with all none-zero partial derivatives saved. 
"""
mutable struct TaylorSeries{T}
    name::String # "" by default
    coeffs::Array{T,1}
    order::Array{Array{Int,1},1}
    """
        function TaylorSeries(T::DataType=Float64, name="", expansion=Dict{Dict{Int,Int},T}(), variables=Set{V}())
            Create a TaylorSeries based on given expansion and variables.
    """
    function TaylorSeries{T}(coeffs::Array{T,1}=Array{T,1}(), order::Array{Array{Int,1},1}=Array{Array{Int,1},1}(), name::String="") where {T}
        return new{T}(name, coeffs, order)
    end
end


function TaylorSeries(::Type{T}, nv::Int) where {T}
    @assert 0 < nv ≤ get_numvars()
    v = zeros(Int, get_numvars())
    @inbounds v[nv] = 1
    return TaylorSeries{T}([one(T)], [v])
end
TaylorSeries(nv::Int) = TaylorSeries(Float64, nv)


