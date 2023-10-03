abstract type AbstractGraph end

abstract type AbstractOperator end
struct Sum <: AbstractOperator end
struct Prod <: AbstractOperator end
struct Constant <: AbstractOperator end
Base.isequal(a::AbstractOperator, b::AbstractOperator) = (typeof(a) == typeof(b))
Base.:(==)(a::AbstractOperator, b::AbstractOperator) = Base.isequal(a, b)
apply(o::AbstractOperator, diags) = error("not implemented!")

Base.show(io::IO, o::AbstractOperator) = print(io, typeof(o))
Base.show(io::IO, ::Type{Sum}) = print(io, "⨁")
Base.show(io::IO, ::Type{Prod}) = print(io, "Ⓧ")
Base.show(io::IO, ::Type{Constant}) = print(io, "C")

# Is the unary form of operator 𝓞 trivial: 𝓞(G) ≡ G?
# NOTE: this property implies that 𝓞(c * G) = c * G = c * 𝓞(G), so
#       we may propagate the subgraphs factor c up to the parent graph.
unary_istrivial(::Type{O}) where {O<:AbstractOperator} = false
unary_istrivial(::Type{O}) where {O<:Union{Sum,Prod}} = true  # (+g) ≡ g and (*g) ≡ g

# Is the operation associative: a 𝓞 (b 𝓞 c) = (a 𝓞 b) 𝓞 c = a 𝓞 b 𝓞 c?
isassociative(::Type{O}) where {O<:AbstractOperator} = false
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

"""
    function Base.:*(g1::AbstractGraph, c2::C) where {C}

    Returns a graph representing the scalar multiplication `g1*c2`.

# Arguments:
- `g1`  computational graph
- `c2`  scalar multiple
"""
# function Base.:*(g1::AbstractGraph, c2::C) where {C}
#     g = similar(g1)
#     g.subgraphs = [g1,]
#     g.subgraph_factors = [c2,]
#     g.operator = Prod

#     # Merge multiplicative link
#     if unary_istrivial(g1.operator) && onechild(g1)
#         g.subgraph_factors[1] *= g1.subgraph_factors[1]
#         g.subgraphs = g1.subgraphs
#     end
#     return g
# end

# """
#     function Base.:*(c1::C, g2::AbstractGraph) where {C}

#     Returns a graph representing the scalar multiplication `c1*g2`.

# # Arguments:
# - `c1`  scalar multiple
# - `g2`  computational graph
# """
# function Base.:*(c1::C, g2::AbstractGraph) where {C}
#     g = similar(g2)
#     g.subgraphs = [g2,]
#     g.subgraph_factors = [c1,]
#     g.operator = Prod
#     # Merge multiplicative link
#     if unary_istrivial(g2.operator) && onechild(g2)
#         g.subgraph_factors[1] *= g2.subgraph_factors[1]
#         g.subgraphs = g2.subgraphs
#     end
#     return g
# end

# """
#     function linear_combination(g1::AbstractGraph, g2::AbstractGraph, c1::C, c2::C) where {C}

#     Returns a graph representing the linear combination `c1*g1 + c2*g2`.

# # Arguments:
# - `g1`  first computational graph
# - `g2`  second computational graph
# - `c1`  first scalar multiple
# - `c2`  second scalar multiple
# """
# function linear_combination(g1::AbstractGraph, g2::AbstractGraph, c1::C, c2::C) where {C}
#     g = similar(g1)
#     g.subgraphs = [g1, g2]
#     g.subgraph_factors = [c1, c2]
#     g.operator = Sum
#     # Convert multiplicative links to in-place form
#     if unary_istrivial(g1.operator) && onechild(g1)
#         g.subgraph_factors[1] *= g1.subgraph_factors[1]
#         g.subgraphs[1] = g1.subgraphs[1]
#     end
#     if unary_istrivial(g2.operator) && onechild(g2)
#         g.subgraph_factors[2] *= g2.subgraph_factors[1]
#         g.subgraphs[2] = g2.subgraphs[1]
#     end
#     return g
# end

# """
#     function linear_combination(graphs::Vector{AbstractGraph}, constants::Vector{C}) where {C}

#     Given a vector 𝐠 of graphs and an equally-sized vector 𝐜 of constants, returns a new
#     graph representing the linear combination (𝐜 ⋅ 𝐠).

# # Arguments:
# - `graphs`  vector of computational graphs
# - `constants`  vector of scalar multiples
# """
# function linear_combination(graphs::Vector{AbstractGraph}, constants::Vector{C}) where {C}
#     if isempty(graphs)
#         return 
#     g = similar(g1)
#     g.subgraphs = graphs
#     g.subgraph_factors = constants
#     g.operator = Sum

#     # Convert multiplicative links to in-place form
#     for (i, sub_g) in enumerate(g.subgraphs)
#         if unary_istrivial(sub_g.operator) && onechild(sub_g)
#             g.subgraph_factors[i] *= sub_g.subgraph_factors[1]
#             g.subgraphs[i] = sub_g.subgraphs[1]
#         end
#     end
#     return g
# end

# # function Base.:+(c::C, g1::AbstractGraph) where {C}
# #     return linear_combination(g1, Unity, F(1), F(c))
# # end
# # function Base.:+(g1::AbstractGraph,c::C) where {C}
# #     return linear_combination(g1, Unity, F(1), F(c))
# # end

# # function Base.:-(c::C, g1::AbstractGraph) where {C}
# #     return linear_combination(Unity, g1, F(c), F(-1))
# # end
# # function Base.:-(g1::AbstractGraph,c::C) where {C}
# #     return linear_combination(g1, Unity, F(1), F(-c))
# # end


# """
#     function Base.:+(g1::AbstractGraph, g2::AbstractGraph) 

#     Returns a graph `g1 + g2` representing the addition of `g2` with `g1`.

# # Arguments:
# - `g1`  first computational graph
# - `g2`  second computational graph
# """
# function Base.:+(g1::AbstractGraph, g2::AbstractGraph)
#     return linear_combination(g1, g2, 1, 1)
# end

# """
#     function Base.:-(g1::AbstractGraph, g2::AbstractGraph) 

#     Returns a graph `g1 - g2` representing the subtraction of `g2` from `g1`.

# # Arguments:
# - `g1`  first computational graph
# - `g2`  second computational graph
# """
# function Base.:-(g1::AbstractGraph, g2::AbstractGraph)
#     return linear_combination(g1, g2, 1, -1)
# end


# """
#     function multi_product(g1::AbstractGraph, g2::AbstractGraph, c1::C, c2::C) where {C}

#     Returns a graph representing the multi product `c1*g1 * c2*g2`.

# # Arguments:
# - `g1`  first computational graph
# - `g2`  second computational graph
# - `c1`  first scalar multiple
# - `c2`  second scalar multiple
# """
# function multi_product(g1::AbstractGraph, g2::AbstractGraph, c1::C=1, c2::C=1) where {C}
#     g = deepcopy(g1)
#     g.subgraphs = [g1, g2]
#     g.subgraph_factors = [c1, c2]
#     g.operator = Prod

#     # Convert multiplicative links to in-place form
#     if unary_istrivial(g1.operator) && onechild(g1)
#         g.subgraph_factors[1] *= g1.subgraph_factors[1]
#         g.subgraphs[1] = g1.subgraphs[1]
#     end
#     if unary_istrivial(g2.operator) && onechild(g2)
#         g.subgraph_factors[2] *= g2.subgraph_factors[1]
#         g.subgraphs[2] = g2.subgraphs[1]
#     end
#     return g
# end

# """
#     function multi_product(graphs::Vector{AbstractGraph}, constants::Vector{C}) where {C}

#     Given a vector 𝐠 of graphs and an equally-sized vector 𝐜 of constants, returns a new
#     graph representing the linear combination (𝐜 ⋅ 𝐠).

# # Arguments:
# - `graphs`  vector of computational graphs
# - `constants`  vector of scalar multiples
# """
# function multi_product(graphs::Vector{AbstractGraph}, constants::Vector{C}=ones(C, length(graphs.subgraphs))) where {C}
#     g = deepcopy(g1)
#     g.subgraphs = graphs
#     g.subgraph_factors = constants
#     g.operator = Prod

#     # Convert multiplicative links to in-place form
#     for (i, sub_g) in enumerate(g.subgraphs)
#         if unary_istrivial(sub_g.operator) && onechild(sub_g)
#             g.subgraph_factors[i] *= sub_g.subgraph_factors[1]
#             g.subgraphs[i] = sub_g.subgraphs[1]
#         end
#     end
#     return g
# end


# """
#     function Base.:*(g1::AbstractGraph, g2::AbstractGraph) 

#     Returns a graph `g1 * g2` representing the graph product between `g1` and `g2`.

# # Arguments:
# - `g1`  first computational graph
# - `g2`  second computational graph
# """
# function Base.:*(g1::AbstractGraph, g2::AbstractGraph)
#     return multi_product(g1, g2)
# end
