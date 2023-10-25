abstract type AbstractGraph end

abstract type AbstractOperator end
struct Sum <: AbstractOperator end
struct Prod <: AbstractOperator end
struct Constant <: AbstractOperator end
struct Power{N} <: AbstractOperator
    function Power(N::Real)
        @assert N ∉ [0, 1] "Power{$N} makes no sense."
        new{N}()
    end
end
Base.eltype(::Type{<:Power{N}}) where {N} = N
decrement_power(::Type{<:Power{N}}) where {N} = N == 2 ? Sum() : Power(N - 1)
# struct Power <: AbstractOperator
#     exponent::Real
# end
Base.isequal(a::AbstractOperator, b::AbstractOperator) = (typeof(a) == typeof(b))
Base.:(==)(a::AbstractOperator, b::AbstractOperator) = Base.isequal(a, b)
apply(o::AbstractOperator, diags) = error("not implemented!")

Base.show(io::IO, o::AbstractOperator) = print(io, typeof(o))
Base.show(io::IO, ::Type{Sum}) = print(io, "⨁")
Base.show(io::IO, ::Type{Prod}) = print(io, "Ⓧ")
Base.show(io::IO, ::Type{Constant}) = print(io, "C")
Base.show(io::IO, ::Type{Power{N}}) where {N} = print(io, "^$N")

# Is the unary form of operator 𝓞 trivial: 𝓞(G) ≡ G?
# NOTE: this property implies that 𝓞(c * G) = c * G = c * 𝓞(G), so
#       we may propagate the subgraph factor c up to the parent graph.
unary_istrivial(::Type{<:AbstractOperator}) = false
unary_istrivial(::Type{<:Union{Sum,Prod}}) = true  # (+g) ≡ g and (*g) ≡ g

# Is the operation associative: a 𝓞 (b 𝓞 c) = (a 𝓞 b) 𝓞 c = a 𝓞 b 𝓞 c?
isassociative(::Type{<:AbstractOperator}) = false
isassociative(::Type{Sum}) = true
# NOTE: Associativity of Prod (graph composition)
#       requires Base.*(g1, g2) and Base./(g1, g2)
# isassociative(::Type{Prod}) = true

function Base.isequal(a::AbstractGraph, b::AbstractGraph)
    typeof(a) != typeof(b) && return false
    for field in fieldnames(typeof(a))
        if field == :weight
            (getproperty(a, :weight) ≈ getproperty(b, :weight)) == false && return false
        else
            getproperty(a, field) != getproperty(b, field) && return false
        end
    end
    return true
end
Base.:(==)(a::AbstractGraph, b::AbstractGraph) = Base.isequal(a, b)

"""
    function isequiv(a::AbstractGraph, b::AbstractGraph, args...)

    Determine whether `a` is equivalent to `b` without considering fields in `args`.
"""
function isequiv(a::AbstractGraph, b::AbstractGraph, args...)
    typeof(a) != typeof(b) && return false
    for field in fieldnames(typeof(a))
        field in args && continue
        if field == :weight
            (getproperty(a, :weight) ≈ getproperty(b, :weight)) == false && return false
        elseif field == :subgraphs
            length(a.subgraphs) != length(b.subgraphs) && return false
            !all(isequiv.(getproperty(a, field), getproperty(b, field), args...)) && return false
        else
            getproperty(a, field) != getproperty(b, field) && return false
        end
    end
    return true
end