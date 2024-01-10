abstract type AbstractGraph end

abstract type AbstractOperator end
struct Sum <: AbstractOperator end
struct Prod <: AbstractOperator end
struct Unitary <: AbstractOperator end
struct Power{N} <: AbstractOperator
    function Power(N::Int)
        @assert N ∉ [0, 1] "Power{$N} makes no sense."
        new{N}()
    end
end
Base.eltype(::Type{<:Power{N}}) where {N} = N
decrement_power(::Type{<:Power{N}}) where {N} = N == 2 ? Sum() : Power(N - 1)
Base.isequal(a::AbstractOperator, b::AbstractOperator) = (typeof(a) == typeof(b))
Base.:(==)(a::AbstractOperator, b::AbstractOperator) = Base.isequal(a, b)
apply(o::AbstractOperator, diags) = error("not implemented!")

Base.show(io::IO, o::AbstractOperator) = print(io, typeof(o))
Base.show(io::IO, ::Type{Sum}) = print(io, "⨁")
Base.show(io::IO, ::Type{Prod}) = print(io, "Ⓧ")
Base.show(io::IO, ::Type{Unitary}) = print(io, "𝟙")
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

"""
    function unary_istrivial(g::AbstractGraph)

    Returns true if the unary form of the graph operation of node `g` is trivial.
    Otherwise, returns false.
"""
unary_istrivial(g::AbstractGraph) = unary_istrivial(operator(g))

"""
    function isassociative(g::AbstractGraph)

    Returns true if the graph operation of node `g` is associative.
    Otherwise, returns false.
"""
isassociative(g::AbstractGraph) = isassociative(operator(g))


### Getters/Setters ###

"""
    function id(g::AbstractGraph)

    Returns the unique hash id of computational graph `g`.
"""
function id(g::AbstractGraph) end

"""
    function name(g::AbstractGraph)

    Returns the name of computational graph `g`.
"""
function name(g::AbstractGraph) end

"""
    function orders(g::AbstractGraph)

    Returns the derivative orders (::Vector{Int}) of computational graph `g`.
"""
function orders(g::AbstractGraph) end

"""
function operator(g::AbstractGraph)
    
    Returns the operation associated with computational graph node `g`.
"""
function operator(g::AbstractGraph) end

"""
    function weight(g::AbstractGraph)

    Returns the weight of the computational graph `g`.
"""
function weight(g::AbstractGraph) end

"""
function subgraph(g::AbstractGraph, i=1)
    
    Returns a copy of the `i`th subgraph of computational graph `g`.
    Defaults to the first subgraph if an index `i` is not supplied.
"""
function subgraph(g::AbstractGraph, i=1) end

"""
function subgraphs(g::AbstractGraph)
    
    Returns the subgraphs of computational graph `g`.
"""
function subgraphs(g::AbstractGraph) end

"""
function subgraphs(g::AbstractGraph, indices::AbstractVector{Int})

    Returns the subgraphs of computational graph `g` at indices `indices`.
    By default, calls `subgraph(g, i)` for each `i` in `indices`. 
"""
function subgraphs(g::AbstractGraph, indices::AbstractVector{Int})
    return [subgraph(g, i) for i in indices]
end

"""
function subgraph_factor(g::AbstractGraph, i=1)
    
    Returns a copy of the `i`th subgraph factor of computational graph `g`.
    Defaults to the first subgraph factor if an index `i` is not supplied.
"""
function subgraph_factor(g::AbstractGraph, i=1) end

"""
function subgraph_factors(g::AbstractGraph)
    
    Returns the subgraph factors of computational graph `g`.
"""
function subgraph_factors(g::AbstractGraph) end

"""
function subgraph_factors(g::AbstractGraph, indices::AbstractVector{Int})

    Returns the subgraph factors of computational graph `g` at indices `indices`.
    By default, calls `subgraph_factor(g, i)` for each `i` in `indices`. 
"""
function subgraph_factors(g::AbstractGraph, indices::AbstractVector{Int})
    return [subgraph_factor(g, i) for i in indices]
end


### Setters ###

"""
function set_id!(g::AbstractGraph, id)

    Update the id of graph `g` to `id`.
"""
function set_id!(g::AbstractGraph, id) end

"""
function set_name!(g::AbstractGraph, name::AbstractString)

    Update the name of graph `g` to `name`.
"""
function set_name!(g::AbstractGraph, name::AbstractString) end

"""
function set_orders!(g::AbstractGraph, orders::AbstractVector)

    Update the orders of graph `g` to `orders`.
"""
function set_orders!(g::AbstractGraph, orders::AbstractVector) end

"""
function set_operator!(g::AbstractGraph, operator::AbstractOperator)

    Update the operator of graph `g` to `typeof(operator)`.
"""
function set_operator!(g::AbstractGraph, operator::AbstractOperator) end

"""
function set_operator!(g::AbstractGraph, operator::Type{<:AbstractOperator})

    Update the operator of graph `g` to `operator`.
"""
function set_operator!(g::AbstractGraph, operator::Type{<:AbstractOperator}) end

"""
function set_weight!(g::AbstractGraph, weight)

    Update the weight of graph `g` to `weight`.
"""
function set_weight!(g::AbstractGraph, weight) end

"""
function set_subgraph!(g::AbstractGraph, subgraph::AbstractGraph, i=1)

    Update the `i`th subgraph of graph `g` to `subgraph`.
    Defaults to the first subgraph factor if an index `i` is not supplied.
"""
function set_subgraph!(g::AbstractGraph, subgraph::AbstractGraph, i=1) end

"""
function set_subgraphs!(g::AbstractGraph, subgraphs::AbstractVector{<:AbstractGraph})

    Update the full list of subgraphs of graph `g` to `subgraphs`.
"""
function set_subgraphs!(g::AbstractGraph, subgraphs::AbstractVector{<:AbstractGraph}) end

"""
function set_subgraphs!(g::AbstractGraph, subgraphs::AbstractVector{<:AbstractGraph}, indices::AbstractVector{Int})

    Update the specified subgraphs of graph `g` at indices `indices` to `subgraphs`.
    By default, calls `set_subgraph!(g, subgraphs[i], i)` for each `i` in `indices`. 
"""
function set_subgraphs!(g::AbstractGraph, subgraphs::AbstractVector{<:AbstractGraph}, indices::AbstractVector{Int})
    @assert length(subgraphs) == length(indices)
    for (i, subg) in zip(indices, subgraphs)
        set_subgraph!(g, subg, i)
    end
end

"""
function set_subgraph_factor!(g::AbstractGraph, subgraph_factor, i=1)

    Update the `i`th subgraph factor of graph `g` to `subgraph_factor`.
    Defaults to the first subgraph factor if an index `i` is not supplied.
"""
function set_subgraph_factor!(g::AbstractGraph, subgraph_factor, i=1) end

"""
function set_subgraph_factors!(g::AbstractGraph, subgraph_factors::AbstractVector)

    Update the subgraph factors of graph `g` to `subgraphs`.
"""
function set_subgraph_factors!(g::AbstractGraph, subgraph_factors::AbstractVector) end

"""
function set_subgraph_factors!(g::AbstractGraph, subgraph_factors::AbstractVector, indices::AbstractVector{Int})

    Update the specified subgraph factors of graph `g` at indices `indices` to `subgraph_factors`.
    By default, calls `set_subgraph_factor!(g, subgraph_factors[i], i)` for each `i` in `indices`. 
"""
function set_subgraph_factors!(g::AbstractGraph, subgraph_factors::AbstractVector, indices::AbstractVector{Int})
    @assert length(subgraph_factors) == length(indices)
    for (i, subg_fac) in zip(indices, subgraph_factors)
        set_subgraph_factor!(g, subg_fac, i)
    end
end

"""
function disconnect_subgraphs!(g::G) where {G<:AbstractGraph}

    Empty the subgraphs and subgraph_factors of graph `g`. Any child nodes of g
    not referenced elsewhere in the full computational graph are effectively deleted.
"""
function disconnect_subgraphs!(g::AbstractGraph)
    empty!(subgraphs(g))
    empty!(subgraph_factors(g))
end

### Methods ###

# Tests for exact equality between two abstract graphs
function Base.isequal(a::AbstractGraph, b::AbstractGraph)
    typeof(a) != typeof(b) && return false
    (weight(a) ≈ weight(b)) == false && return false  # check graph weights for approximate equality
    length(subgraphs(a)) != length(subgraphs(b)) && return false

    pa = sortperm(subgraphs(a), by=x -> id(x))
    pb = sortperm(subgraphs(b), by=x -> id(x))
    subgraph_factors(a)[pa] != subgraph_factors(b)[pb] && return false
    subgraphs(a)[pa] != subgraphs(b)[pb] && return false

    for field in fieldnames(typeof(a))
        if field in [:weight, :subgraphs, :subgraph_factors]
            continue
            # if field == :weight && getproperty(a, :weight) == weight(a) && getproperty(b, :weight) == weight(b)
            #     continue  # skip graph weights if already accounted for
            # elseif field == :subgraphs && getproperty(a, :subgraphs) == subgraphs(a) && getproperty(b, :subgraphs) == subgraphs(b)
            #     continue  # skip subgraphs if already accounted for
            # elseif field == :subgraph_factors && getproperty(a, :subgraph_factors) == subgraph_factors(a) && getproperty(b, :subgraph_factors) == subgraph_factors(b)
            #     continue  # skip subgraph_factors if already accounted for
        else
            getproperty(a, field) != getproperty(b, field) && return false
        end
    end
    return true
end
Base.:(==)(a::AbstractGraph, b::AbstractGraph) = Base.isequal(a, b)

"""
    function isequiv(a::AbstractGraph, b::AbstractGraph, args...)

    Determine whether graph `a` is equivalent to graph `b` without considering fields in `args`.
"""
function isequiv(a::AbstractGraph, b::AbstractGraph, args...)
    typeof(a) != typeof(b) && return false
    # Check graph weights for approximate equality where applicable
    if :weight ∉ args
        (weight(a) ≈ weight(b)) == false && return false
    end
    # Check that all subgraphs are equivalent modulo `args`
    length(subgraphs(a)) != length(subgraphs(b)) && return false

    # if :id ∉ args
    #     pa = sortperm(subgraphs(a), by=x -> id(x))
    #     pb = sortperm(subgraphs(b), by=x -> id(x))
    #     subgraph_factors(a)[pa] != subgraph_factors(b)[pb] && return false
    #     !all(isequiv.(subgraphs(a)[pa], subgraphs(b)[pb], args...)) && return false
    # else
    a_pairs = collect(zip(subgraphs(a), subgraph_factors(a)))
    b_pairs = collect(zip(subgraphs(b), subgraph_factors(b)))
    for (suba, suba_factor) in a_pairs
        match_found = false
        for (idx, (subb, subb_factor)) in enumerate(b_pairs)
            if suba_factor == subb_factor && isequiv(suba, subb, args...)
                deleteat!(b_pairs, idx)
                match_found = true
                break
            end
        end
        !match_found && return false
    end
    # end

    for field in fieldnames(typeof(a))
        if field in [:weight, :subgraphs, :subgraph_factors, args...]
            continue
            # if field == :weight && getproperty(a, :weight) == weight(a) && getproperty(b, :weight) == weight(b)
            #     continue  # skip graph weights if already accounted for
            # elseif field == :subgraphs && getproperty(a, :subgraphs) == subgraphs(a) && getproperty(b, :subgraphs) == subgraphs(b)
            #     continue  # skip subgraphs if already accounted for
            # elseif field == :subgraph_factors && getproperty(a, :subgraph_factors) == subgraph_factors(a) && getproperty(b, :subgraph_factors) == subgraph_factors(b)
            #     continue  # skip subgraph_factors if already accounted for
        else
            # getproperty(a, field) != getproperty(b, field) && return false
            if getproperty(a, field) != getproperty(b, field)
                # println(field)
                return false
            end
        end
    end
    return true
end


### Arithmetic operations ###

errmsg(G::Type) = "Method not yet implemented for user-defined graph type $G."
linear_combination(g1::G, g2::G, c1, c2) where {G<:AbstractGraph} = error(errmsg(G))
linear_combination(graphs::AbstractVector{G}, constants::AbstractVector) where {G<:AbstractGraph} = error(errmsg(G))
multi_product(g1::G, g2::G, c1, c2) where {G<:AbstractGraph} = error(errmsg(G))
multi_product(graphs::AbstractVector{G}, constants::AbstractVector) where {G<:AbstractGraph} = error(errmsg(G))
Base.:*(c1, g2::G) where {G<:AbstractGraph} = error(errmsg(G))
Base.:*(g1::G, c2) where {G<:AbstractGraph} = error(errmsg(G))
Base.:*(g1::G, g2::G) where {G<:AbstractGraph} = error(errmsg(G))
Base.:+(g1::G, g2::G) where {G<:AbstractGraph} = error(errmsg(G))
Base.:-(g1::G, g2::G) where {G<:AbstractGraph} = error(errmsg(G))
