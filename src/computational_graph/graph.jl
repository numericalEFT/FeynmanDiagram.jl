"""
    mutable struct Graph{F,W}
    
    A representation of a computational graph, e.g., an expression tree, with type stable node data.

# Members:
- `id::Int`  the unique hash id to identify the diagram
- `name::Symbol`  name of the diagram
- `orders::Vector{Int}`  orders associated with the graph, e.g., derivative orders
- `subgraphs::Vector{Graph{F,W}}`  vector of sub-diagrams 
- `subgraph_factors::Vector{F}`  scalar multiplicative factors associated with each subgraph. Note that the subgraph factors may be manipulated algebraically. To associate a fixed multiplicative factor with this graph which carries some semantic meaning, use the `factor` argument instead.
- `operator::DataType`  node operation. Addition and multiplication are natively supported via operators Sum and Prod, respectively. Should be a concrete subtype of `AbstractOperator`.
- `factor::F`  total scalar multiplicative factor for the diagram
- `weight::W`  the weight of this node

# Example:
```julia-repl
julia> g1 = Graph([])
1=0.0

julia> g2 = Graph([]; factor=2)
2⋅2.0=0.0

julia> g = Graph([g1, g2]; operator=ComputationalGraphs.Sum())
3=0.0=⨁ (1,2)
```
"""
mutable struct Graph{F,W} <: AbstractGraph # Graph
    id::Int
    name::String # "" by default
    orders::Vector{Int}

    subgraphs::Vector{Graph{F,W}}
    subgraph_factors::Vector{F}
    parent_graphs::Vector{Graph{F,W}}

    operator::DataType
    factor::F
    weight::W

    """
        function Graph(subgraphs::AbstractVector; name="", operator::AbstractOperator=Sum(),
            ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype))
        
        Create a Graph struct from a set of subgraphs with the specified node data and operation.

    # Arguments:
    - `subgraphs`  vector of sub-diagrams 
    - `subgraph_factors`  scalar multiplicative factors associated with each subgraph. Note that the subgraph factors may be manipulated algebraically. To associate a fixed multiplicative factor with this graph which carries some semantic meaning, use the `factor` argument instead.
    - `name`  name of the diagram
    - `orders`  orders associated with the graph, e.g., derivative orders
    - `operator`  node operation, i.e., Sum, Prod, or a user-defined operator `Op <: AbstractOperator`
    - `ftype`  typeof(factor)
    - `wtype`  typeof(weight)
    - `factor`  fixed scalar multiplicative factor for this diagram (e.g., a permutation sign)
    - `weight`  the weight of this node
    """
    function Graph(subgraphs::AbstractVector; subgraph_factors=one.(eachindex(subgraphs)), parent_graphs::AbstractVector=eltype(subgraphs)[],
        name="", operator::AbstractOperator=Sum(), orders=zeros(Int, 16),
        ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype)
    )
        if typeof(operator) <: Power
            @assert length(subgraphs) == 1 "Graph with Power operator must have one and only one subgraph."
        end
        # @assert allunique(subgraphs) "all subgraphs must be distinct."
        # return new{ftype,wtype}(uid(), name, orders, subgraphs, subgraph_factors, parent_graphs, typeof(operator), factor, weight)
        g = new{ftype,wtype}(uid(), name, orders, subgraphs, subgraph_factors, parent_graphs, typeof(operator), factor, weight)
        for sub_g in subgraphs
            g ∉ sub_g.parent_graphs && push!(sub_g.parent_graphs, g)
        end
        return g
    end
end

"""
    function constant_graph(factor=one(_dtype.factor))

    Returns a graph that represents a constant equal to f, where f is the factor with default value 1.

# Arguments:
- `f`:  constant factor
"""
function constant_graph(factor=one(_dtype.factor))
    return Graph([]; operator=Constant(), factor=factor, ftype=_dtype.factor, wtype=_dtype.weight, weight=one(_dtype.weight))
end

"""
    function orders(g::Graph)

    Returns the derivative orders (::Vector{Int}) of Graph `g`.
"""
orders(g::Graph) = g.orders

"""
    function Base.:*(g1::Graph{F,W}, c2::C) where {F,W,C}

    Returns a graph representing the scalar multiplication `g1*c2`.

# Arguments:
- `g1`  computational graph
- `c2`  scalar multiple
"""
function Base.:*(g1::Graph{F,W}, c2::C) where {F,W,C}
    g = Graph([g1,]; subgraph_factors=[F(c2),], operator=Prod(), orders=orders(g1), ftype=F, wtype=W)
    # Merge multiplicative link
    if unary_istrivial(g1.operator) && onechild(g1)
        g.subgraph_factors[1] *= g1.subgraph_factors[1]
        g.subgraphs = g1.subgraphs
        pop!(g1.parent_graphs)
        push!(g1.subgraphs[1].parent_graphs, g)
    end
    return g
end

"""
    function Base.:*(c1::C, g2::Graph{F,W}) where {F,W,C}

    Returns a graph representing the scalar multiplication `c1*g2`.

# Arguments:
- `c1`  scalar multiple
- `g2`  computational graph
"""
function Base.:*(c1::C, g2::Graph{F,W}) where {F,W,C}
    g = Graph([g2,]; subgraph_factors=[F(c1),], operator=Prod(), orders=orders(g2), ftype=F, wtype=W)
    # Merge multiplicative link
    if unary_istrivial(g2.operator) && onechild(g2)
        g.subgraph_factors[1] *= g2.subgraph_factors[1]
        g.subgraphs = g2.subgraphs
        pop!(g2.parent_graphs)
        push!(g2.subgraphs[1].parent_graphs, g)
    end
    return g
end

"""
    function linear_combination(g1::Graph{F,W}, g2::Graph{F,W}, c1::C, c2::C) where {F,W,C}

    Returns a graph representing the linear combination `c1*g1 + c2*g2`.
    If `g1 == g2`, it will return a graph representing `(c1+c2)*g1`.
    Graphs `g1` and `g2` must have the same orders.

# Arguments:
- `g1`  first computational graph
- `g2`  second computational graph
- `c1`  first scalar multiple
- `c2`  second scalar multiple
"""
function linear_combination(g1::Graph{F,W}, g2::Graph{F,W}, c1::C=1, c2::C=1) where {F,W,C}
    @assert orders(g1) == orders(g2) "g1 and g2 have different orders."
    subgraphs = [g1, g2]
    subgraph_factors = [F(c1), F(c2)]
    # Convert multiplicative links to in-place form
    if unary_istrivial(g1.operator) && onechild(g1)
        subgraph_factors[1] *= g1.subgraph_factors[1]
        subgraphs[1] = g1.subgraphs[1]
    end
    if unary_istrivial(g2.operator) && onechild(g2)
        subgraph_factors[2] *= g2.subgraph_factors[1]
        subgraphs[2] = g2.subgraphs[1]
    end

    if subgraphs[1] == subgraphs[2]
        g = Graph([subgraphs[1]]; subgraph_factors=[sum(subgraph_factors)], operator=Sum(), orders=orders(g1), ftype=F, wtype=W)
    else
        g = Graph(subgraphs; subgraph_factors=subgraph_factors, operator=Sum(), orders=orders(g1), ftype=F, wtype=W)
    end
    return g
end

"""
    function linear_combination(graphs::Vector{Graph{F,W}}, constants::Vector{C}) where {F,W,C}

    Given a vector 𝐠 of graphs and an equally-sized vector 𝐜 of constants, returns a new
    graph representing the linear combination (𝐜 ⋅ 𝐠). 
    The function identifies unique graphs from the input `graphs` and sums their associated `constants`.
    All input graphs must have the same orders.

# Arguments:
- `graphs`  vector of computational graphs
- `constants`  vector of scalar multiples (defaults to ones(C, length(graphs))).

# Returns:
- A new `Graph{F,W}` object representing the linear combination of the unique input `graphs` weighted by the constants, 
where duplicate graphs in the input `graphs` are combined by summing their associated constants. 

# Example:
    Given graphs `g1`, `g2`, `g1` and constants `c1`, `c2`, `c3`, the function computes `(c1+c3)*g1 + c2*g2`.
"""
function linear_combination(graphs::Vector{Graph{F,W}}, constants::Vector{C}=ones(C, length(graphs))) where {F,W,C}
    @assert alleq(orders.(graphs)) "Graphs do not all have the same order."
    subgraphs, subgraph_factors = graphs, constants
    # parameters = union(getproperty.(graphs, :parameters))
    # Convert multiplicative links to in-place form
    for (i, sub_g) in enumerate(graphs)
        if unary_istrivial(sub_g.operator) && onechild(sub_g)
            subgraph_factors[i] *= sub_g.subgraph_factors[1]
            subgraphs[i] = sub_g.subgraphs[1]
        end
    end

    unique_graphs = Graph{F,W}[]
    unique_factors = F[]
    for (idx, g) in enumerate(subgraphs)
        i = findfirst(isequal(g), unique_graphs)
        if isnothing(i)
            push!(unique_graphs, g)
            push!(unique_factors, subgraph_factors[idx])
        else
            unique_factors[i] += subgraph_factors[idx]
        end
    end
    g = Graph(unique_graphs; subgraph_factors=unique_factors, operator=Sum(), orders=orders(graphs[1]), ftype=F, wtype=W)

    return g
end

# function Base.:+(c::C, g1::Graph{F,W}) where {F,W,C}
#     return linear_combination(g1, Unity, F(1), F(c))
# end
# function Base.:+(g1::Graph{F,W},c::C) where {F,W,C}
#     return linear_combination(g1, Unity, F(1), F(c))
# end

# function Base.:-(c::C, g1::Graph{F,W}) where {F,W,C}
#     return linear_combination(Unity, g1, F(c), F(-1))
# end
# function Base.:-(g1::Graph{F,W},c::C) where {F,W,C}
#     return linear_combination(g1, Unity, F(1), F(-c))
# end


"""
    function Base.:+(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}

    Returns a graph `g1 + g2` representing the addition of `g2` with `g1`.
    Graphs `g1` and `g2` must have the same orders.

# Arguments:
- `g1`  first computational graph
- `g2`  second computational graph
"""
function Base.:+(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
    return linear_combination(g1, g2, F(1), F(1))
end

"""
    function Base.:-(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}

    Returns a graph `g1 - g2` representing the subtraction of `g2` from `g1`.
    Graphs `g1` and `g2` must have the same orders.

# Arguments:
- `g1`  first computational graph
- `g2`  second computational graph
"""
function Base.:-(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
    return linear_combination(g1, g2, F(1), F(-1))
end


"""
    function multi_product(g1::Graph{F,W}, g2::Graph{F,W}, c1::C=1, c2::C=1) where {F,W,C}

    Returns a graph representing the multi product `c1*g1 * c2*g2`.
    If `g1 == g2`, it will return a graph representing `c1*c2 * (g1)^2` with `Power(2)` operator.

# Arguments:
- `g1`:  first computational graph
- `g2`:  second computational graph
- `c1`:  first scalar multiple (defaults to 1).
- `c2`:  second scalar multiple (defaults to 1).
"""
function multi_product(g1::Graph{F,W}, g2::Graph{F,W}, c1::C=1, c2::C=1) where {F,W,C}
    # g = Graph([g1, g2]; subgraph_factors=[F(c1), F(c2)], operator=Prod(), ftype=F, wtype=W)
    subgraphs = [g1, g2]
    subgraph_factors = [F(c1), F(c2)]
    # Convert multiplicative links to in-place form
    if unary_istrivial(g1.operator) && onechild(g1)
        subgraph_factors[1] *= g1.subgraph_factors[1]
        subgraphs[1] = g1.subgraphs[1]
    end
    if unary_istrivial(g2.operator) && onechild(g2)
        subgraph_factors[2] *= g2.subgraph_factors[1]
        subgraphs[2] = g2.subgraphs[1]
    end

    if subgraphs[1] == subgraphs[2]
        g = Graph([subgraphs[1]]; subgraph_factors=[prod(subgraph_factors)], operator=Power(2), ftype=F, wtype=W)
    else
        g = Graph(subgraphs; subgraph_factors=subgraph_factors, operator=Prod(), ftype=F, wtype=W)
    end
    return g
end

"""
    multi_product(graphs::Vector{Graph{F,W}}, constants::Vector{C}=ones(C, length(graphs.subgraphs))) where {F,W,C}

    Construct a product graph from multiple input graphs, where each graph can be weighted by a constant. 
    For graphs that are repeated more than once, it adds a power operator to the subgraph to represent the repetition.
    Moreover, it optimizes any trivial unary operators in the resulting product graph.

# Arguments:
- `graphs::Vector{Graph{F,W}}`: A vector of input graphs to be multiplied.
- `constants::Vector{C}`: A vector of scalar multiples. If not provided, it defaults to a vector of ones of the same length as `graphs`.

Returns:
- A new product graph with the unique subgraphs (or powered versions thereof) and the associated constants as subgraph factors.

# Example:
    Given graphs `g1`, `g2`, `g1` and constants `c1`, `c2`, `c3`, the function computes `(c1*c3)*(g1)^2 * c2*g2`.
"""
function multi_product(graphs::Vector{Graph{F,W}}, constants::Vector{C}=ones(C, length(graphs))) where {F,W,C}
    subgraphs, subgraph_factors = graphs, constants
    # Convert multiplicative links to in-place form
    for (i, sub_g) in enumerate(graphs)
        if unary_istrivial(sub_g.operator) && onechild(sub_g)
            subgraph_factors[i] *= sub_g.subgraph_factors[1]
            subgraphs[i] = sub_g.subgraphs[1]
        end
    end

    unique_graphs = Vector{Graph{F,W}}()
    unique_factors = F[]
    repeated_counts = Int[]
    for (idx, g) in enumerate(subgraphs)
        loc = findfirst(isequal(g), unique_graphs)
        if isnothing(loc)
            push!(unique_graphs, g)
            push!(unique_factors, subgraph_factors[idx])
            push!(repeated_counts, 1)
        else
            unique_factors[loc] *= subgraph_factors[idx]
            repeated_counts[loc] += 1
        end
    end

    if length(unique_factors) == 1
        g = Graph(unique_graphs; subgraph_factors=unique_factors, operator=Power(repeated_counts[1]), ftype=F, wtype=W)
    else
        subgraphs = Vector{Graph{F,W}}()
        for (idx, g) in enumerate(unique_graphs)
            if repeated_counts[idx] == 1
                push!(subgraphs, g)
            else
                push!(subgraphs, Graph([g], operator=Power(repeated_counts[idx]), ftype=F, wtype=W))
            end
        end
        g = Graph(subgraphs; subgraph_factors=unique_factors, operator=Prod(), ftype=F, wtype=W)
    end
    return g
end


"""
    function Base.:*(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}

    Returns a graph `g1 * g2` representing the graph product between `g1` and `g2`.

# Arguments:
- `g1`  first computational graph
- `g2`  second computational graph
"""
function Base.:*(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
    return multi_product(g1, g2)
end
