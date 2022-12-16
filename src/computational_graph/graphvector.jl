"""
struct GraphVector{F,W} <: AbstractVector
"""
struct GraphVector{F,W} <: AbstractVector
    graphs::Vector{Graph{F,W}}
    GraphVector(graphs::AbstractVector) = new{_dtype.factor,_dtype.weight}(graphs)
    GraphVector{F,W}(graphs::AbstractVector) = new{F,W}(graphs)
end



