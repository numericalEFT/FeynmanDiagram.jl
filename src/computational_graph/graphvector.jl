# """
#     struct GraphVector{F,W} <: AbstractVector{Graph{F,W}}

# A vector of graphs. All graphs must have the same number of external legs, but can have different external labels.

# Support two constructors:
#     GraphVector(graphs::AbstractVector) # construct from a vector of Graph,  use the default factor and weight types
#     GraphVector{F,W}(graphs::AbstractVector) # construct from a vector of Graph,  use the specified factor and weight types `F` and `W`
# """
# struct GraphVector{F,W} <: AbstractVector{Graph{F,W}}
#     graphs::Vector{Graph{F,W}}
#     function GraphVector{F,W}(graphs::AbstractVector) where {F,W}
#         l = length(graphs[1].external)
#         @assert all(x -> length(x.external) == l, graphs)
#         return new{F,W}(graphs)
#     end
#     GraphVector(graphs::AbstractVector) = GraphVector{_dtype.factor,_dtype.weight}(graphs)
# end

# ########## iterator interface  ##########
# Base.eltype(::Type{GraphVector{F,W}}) where {F,W} = Graph{F,W}
# Base.length(v::GraphVector) = length(v.graphs)
# Base.size(v::GraphVector) = size(v.graphs)
# Base.getindex(v::GraphVector, i::Int) = v.graphs[i]
# Base.setindex!(v::GraphVector, g::Graph, i::Int) = v.graphs[i] = g
# Base.firstindex(v::GraphVector) = firstindex(v.graphs)
# Base.lastindex(v::GraphVector) = lastindex(v.graphs)
# Base.resize!(v::GraphVector, n::Int) = resize!(v.graphs, n)


# function merge(gv::GraphVector, indices::Vector{Int})
#     return [GraphVector([gv[i] for i in range(1, length(gv)) if i ∉ indices]) for i in indices]
# end

"""
    function group(gv::Vector{Graph}, indices::Vector{Int})

Group the graphs in `gv` by the external operators at the indices `indices`. Return a dictionary of `Vector{OperatorProduct}` to `GraphVector`.

# Example

```julia-repl
julia> p1 = propagator([𝑓⁺(1), 𝑓⁻(2)]);

julia> p2 = propagator([𝑓⁺(1), 𝑓⁻(3)]);

julia> p3 = propagator([𝑓⁺(2), 𝑓⁻(3)]);

julia> gv = [p1, p2, p3];

julia> ComputationalGraphs.group(gv, [1, 2])
Dict{Vector{OperatorProduct}, Vector{Graph{Float64, Float64}}} with 3 entries:
  [f⁺(1), f⁻(3)] => [2:f⁺(1)|f⁻(3)=0.0]
  [f⁺(1), f⁻(2)] => [1:f⁺(1)|f⁻(2)=0.0]
  [f⁺(2), f⁻(3)] => [3:f⁺(2)|f⁻(3)=0.0]

julia> ComputationalGraphs.group(gv, [1, ])
Dict{Vector{OperatorProduct}, Vector{Graph{Float64, Float64}}} with 2 entries:
  [f⁺(2)] => [3:f⁺(2)|f⁻(3)=0.0]
  [f⁺(1)] => [1:f⁺(1)|f⁻(2)=0.0, 2:f⁺(1)|f⁻(3)=0.0]

julia> ComputationalGraphs.group(gv, [2, ])
Dict{Vector{OperatorProduct}, Vector{Graph{Float64, Float64}}} with 2 entries:
  [f⁻(3)] => [2:f⁺(1)|f⁻(3)=0.0, 3:f⁺(2)|f⁻(3)=0.0]
  [f⁻(2)] => [1:f⁺(1)|f⁻(2)=0.0]
```

"""
function group(gv::AbstractVector{G}, indices::Vector{Int}) where {G<:Graph}
  l = length(gv[1].external)
  @assert all(x -> length(x.external) == l, gv)
  groups = Dict{Vector{OperatorProduct},Vector{G}}()
  for t in gv
    ext = external(t)
    key = [OperatorProduct(ext[i]) for i in indices]
    if haskey(groups, key)
      push!(groups[key], t)
    else
      groups[key] = [t,]
    end
  end
  return groups
end





