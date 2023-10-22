abstract type AbstractGraph end

abstract type AbstractOperator end
struct Sum <: AbstractOperator end
struct Prod <: AbstractOperator end
Base.isequal(a::AbstractOperator, b::AbstractOperator) = (typeof(a) == typeof(b))
Base.:(==)(a::AbstractOperator, b::AbstractOperator) = Base.isequal(a, b)
apply(o::AbstractOperator, diags) = error("not implemented!")

Base.show(io::IO, o::AbstractOperator) = print(io, typeof(o))
Base.show(io::IO, ::Type{Sum}) = print(io, "⨁")
Base.show(io::IO, ::Type{Prod}) = print(io, "Ⓧ")

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
function factor(g::AbstractGraph)

    Returns the fixed scalar-multiplicative factor of the computational graph `g`.
"""
function factor(g::AbstractGraph) end

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
function subgraphs(g::AbstractGraph, indices::AbstractVector{Int})

    Returns the subgraph factors of computational graph `g` at indices `indices`.
    By default, calls `subgraph_factor(g, i)` for each `i` in `indices`. 
"""
function subgraph_factors(g::AbstractGraph, indices::AbstractVector{Int})
    return [subgraph_factor(g, i) for i in indices]
end


### Setters ###

"""
function set_name!(g::AbstractGraph, name::AbstractString)

    Update the name of graph `g` to `name`.
"""
function set_name!(g::AbstractGraph, name::AbstractString) end

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

### Methods ###

# Tests for exact equality between two abstract graphs
function Base.isequal(a::AbstractGraph, b::AbstractGraph)
    typeof(a) != typeof(b) && return false
    (weight(a) ≈ weight(b)) == false && return false  # graph weights approximately equal
    for field in fieldnames(typeof(a))
        if field == :weight && getproperty(a, :weight) == weight(a) && getproperty(b, :weight) == weight(b)
            continue  # skip graph weights if already accounted for
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
    (weight(a) ≈ weight(b)) == false && return false  # graph weights approximately equal
    # Check that all subgraphs are equivalent modulo `args`
    length(subgraphs(a)) != length(subgraphs(b)) && return false
    !all(isequiv.(subgraphs(a), subgraphs(b), args...)) && return false
    for field in fieldnames(typeof(a))
        field in args && continue
        if field == :weight && getproperty(a, :weight) == weight(a) && getproperty(b, :weight) == weight(b)
            continue  # skip graph weights if already accounted for
        elseif field == :subgraphs && getproperty(a, :subgraphs) == subgraphs(a) && getproperty(b, :subgraphs) == subgraphs(b)
            continue  # skip subgraphs if already accounted for
        else
            getproperty(a, field) != getproperty(b, field) && return false
        end
    end
    return true
end
