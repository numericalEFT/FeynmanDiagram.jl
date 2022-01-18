abstract type DiagramId end

Base.Dict(x::DiagramId) = Dict{Symbol,Any}([fn => getfield(x, fn) for fn ∈ fieldnames(typeof(x))])
Base.show(io::IO, d::DiagramId) = error("Base.show not implemented!")
Base.isequal(a::DiagramId, b::DiagramId) = error("Base.isequal not implemented!")
Base.:(==)(a::DiagramId, b::DiagramId) = Base.isequal(a, b)
eval(d::DiagramId) = error("eval for $d has not yet implemented!")

abstract type Operator end
struct Add <: Operator end
struct Mutiply <: Operator end
Base.show(io::IO, o::Operator) = print(io, typeof(o))
Base.isequal(a::Operator, b::Operator) = (typeof(a) == typeof(b))
Base.:(==)(a::Operator, b::Operator) = Base.isequal(a, b)

apply(o::Operator, weights) = error("not implemented!")
apply(o::Add, weights::Vector{W}) where {W<:Number} = sum(weights)
apply(o::Add, weight::Number) = weight
apply(o::Mutiply, weights::Number) = weight


