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
    function Graph(subgraphs::AbstractVector; subgraph_factors=one.(eachindex(subgraphs)), name="", operator::AbstractOperator=Sum(),
        orders=zeros(Int, 16), ftype=_dtype.factor, wtype=_dtype.weight, factor=one(ftype), weight=zero(wtype)
    )
        if typeof(operator) <: Power
            @assert length(subgraphs) == 1 "Graph with Power operator must have one and only one subgraph."
        end
        # @assert allunique(subgraphs) "all subgraphs must be distinct."
        return new{ftype,wtype}(uid(), name, orders, subgraphs, subgraph_factors, typeof(operator), factor, weight)
    end
end

### AbstractGraph interface for Graph ###

# Getters
id(g::Graph) = g.id
name(g::Graph) = g.name
orders(g::Graph) = g.orders
operator(g::Graph) = g.operator
factor(g::Graph) = g.factor
weight(g::Graph) = g.weight
subgraph(g::Graph, i=1) = g.subgraphs[i]
subgraphs(g::Graph) = g.subgraphs
subgraphs(g::Graph, indices::AbstractVector{Int}) = g.subgraphs[indices]
subgraph_factor(g::Graph, i=1) = g.subgraph_factors[i]
subgraph_factors(g::Graph) = g.subgraph_factors
subgraph_factors(g::Graph, indices::AbstractVector{Int}) = g.subgraph_factors[indices]

# Setters
set_id!(g::Graph, id::Int) = (g.id = id)
set_name!(g::Graph, name::String) = (g.name = name)
set_orders!(g::Graph, orders::Vector{Int}) = (g.orders = orders)
set_operator!(g::Graph, operator::Type{<:AbstractOperator}) = (g.operator = operator)
set_operator!(g::Graph, operator::AbstractOperator) = (g.operator = typeof(operator))
set_factor!(g::Graph{F,W}, factor) where {F,W} = (g.factor = F(factor))
set_weight!(g::Graph{F,W}, weight) where {F,W} = (g.weight = W(weight))
set_subgraph!(g::Graph{F,W}, subgraph::Graph{F,W}, i=1) where {F,W} = (g.subgraphs[i] = subgraph)
set_subgraphs!(g::Graph{F,W}, subgraphs::Vector{Graph{F,W}}) where {F,W} = (g.subgraphs = subgraphs)
set_subgraphs!(g::Graph{F,W}, subgraphs::Vector{Graph{F,W}}, indices::AbstractVector{Int}) where {F,W} = (g.subgraphs[indices] = subgraphs)
set_subgraph_factor!(g::Graph{F,W}, subgraph_factor, i=1) where {F,W} = (g.subgraph_factors[i] = F(subgraph_factor))
set_subgraph_factors!(g::Graph{F,W}, subgraph_factors::AbstractVector) where {F,W} = (g.subgraph_factors = Vector{F}(subgraph_factors))
set_subgraph_factors!(g::Graph{F,W}, subgraph_factors::AbstractVector, indices::AbstractVector{Int}) where {F,W} = (g.subgraph_factors[indices] = Vector{F}(subgraph_factors))

###############################

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
    function Base.:*(g1::Graph{F,W}, c2) where {F,W}

    Returns a graph representing the scalar multiplication `g1*c2`.

# Arguments:
- `g1`  computational graph
- `c2`  scalar multiple
"""
function Base.:*(g1::Graph{F,W}, c2) where {F,W}
    g = Graph([g1,]; subgraph_factors=[F(c2),], operator=Prod(), orders=orders(g1), ftype=F, wtype=W)
    # Convert trivial unary link to in-place form
    if unary_istrivial(g1) && onechild(g1)
        g.subgraph_factors[1] *= g1.subgraph_factors[1]
        # g.subgraph_factors[1] *= g1.subgraph_factors[1] * g1.factor
        g.subgraphs = g1.subgraphs
    end
    return g
end

"""
    function Base.:*(c1, g2::Graph{F,W}) where {F,W}

    Returns a graph representing the scalar multiplication `c1*g2`.

# Arguments:
- `c1`  scalar multiple
- `g2`  computational graph
"""
function Base.:*(c1, g2::Graph{F,W}) where {F,W}
    g = Graph([g2,]; subgraph_factors=[F(c1),], operator=Prod(), orders=orders(g2), ftype=F, wtype=W)
    # Convert trivial unary link to in-place form
    if unary_istrivial(g2) && onechild(g2)
        g.subgraph_factors[1] *= g2.subgraph_factors[1]
        # g.subgraph_factors[1] *= g2.subgraph_factors[1] * g2.factor
        g.subgraphs = g2.subgraphs
    end
    return g
end

"""
    function linear_combination(g1::Graph{F,W}, g2::Graph{F,W}, c1, c2) where {F,W}

    Returns a graph representing the linear combination `c1*g1 + c2*g2`.
    If `g1 == g2`, it will return a graph representing `(c1+c2)*g1`.

# Arguments:
- `g1`  first computational graph
- `g2`  second computational graph
- `c1`  first scalar multiple
- `c2`  second scalar multiple
"""
function linear_combination(g1::Graph{F,W}, g2::Graph{F,W}, c1=F(1), c2=F(1)) where {F,W}
    @assert orders(g1) == orders(g2) "g1 and g2 have different orders."
    f1 = typeof(c1) == F ? c1 : F(c1)
    f2 = typeof(c2) == F ? c2 : F(c2)
    subgraphs = [g1, g2]
    subgraph_factors = [f1, f2]
    # Convert trivial unary links to in-place form
    if unary_istrivial(g1) && onechild(g1)
        subgraph_factors[1] *= g1.subgraph_factors[1]
        # subgraph_factors[1] *= g1.subgraph_factors[1] * g1.factor
        subgraphs[1] = g1.subgraphs[1]
    end
    if unary_istrivial(g2) && onechild(g2)
        subgraph_factors[2] *= g2.subgraph_factors[1]
        # subgraph_factors[2] *= g2.subgraph_factors[1] * g2.factor
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
    function linear_combination(graphs::Vector{Graph{F,W}}, constants::AbstractVector=ones(F, length(graphs))) where {F,W}

    Given a vector 𝐠 of graphs and an equally-sized vector 𝐜 of constants, returns a new
    graph representing the linear combination (𝐜 ⋅ 𝐠). 
    The function identifies unique graphs from the input `graphs` and sums their associated `constants`.

# Arguments:
- `graphs`  vector of computational graphs
- `constants`  vector of scalar multiples (defaults to ones(F, length(graphs))).

# Returns:
- A new `Graph{F,W}` object representing the linear combination of the unique input `graphs` weighted by the constants, 
where duplicate graphs in the input `graphs` are combined by summing their associated constants. 

# Example:
    Given graphs `g1`, `g2`, `g1` and constants `c1`, `c2`, `c3`, the function computes `(c1+c3)*g1 + c2*g2`.
"""
function linear_combination(graphs::Vector{Graph{F,W}}, constants::AbstractVector=ones(F, length(graphs))) where {F,W}
    @assert alleq(orders.(graphs)) "Graphs do not all have the same order."
    subgraphs = graphs
    subgraph_factors = eltype(constants) == F ? constants : Vector{F}(constants)
    # Convert trivial unary links to in-place form
    for (i, sub_g) in enumerate(graphs)
        if unary_istrivial(sub_g) && onechild(sub_g)
            subgraph_factors[i] *= sub_g.subgraph_factors[1]
            # subgraph_factors[i] *= sub_g.subgraph_factors[1] * sub_g.factor
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
    if isempty(unique_graphs)
        return nothing
    end
    g = Graph(unique_graphs; subgraph_factors=unique_factors, operator=Sum(), orders=orders(graphs[1]), ftype=F, wtype=W)

    return g
end

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
    function multi_product(g1::Graph{F,W}, g2::Graph{F,W}, c1=F(1), c2=F(1)) where {F,W,C}

    Returns a graph representing the multi product `c1*g1 * c2*g2`.
    If `g1 == g2`, it will return a graph representing `c1*c2 * (g1)^2` with `Power(2)` operator.

# Arguments:
- `g1`:  first computational graph
- `g2`:  second computational graph
- `c1`:  first scalar multiple (defaults to 1).
- `c2`:  second scalar multiple (defaults to 1).
"""
function multi_product(g1::Graph{F,W}, g2::Graph{F,W}, c1=F(1), c2=F(1)) where {F,W}
    @assert orders(g1) == orders(g2) "g1 and g2 have different orders."
    f1 = typeof(c1) == F ? c1 : F(c1)
    f2 = typeof(c2) == F ? c2 : F(c2)
    subgraphs = [g1, g2]
    subgraph_factors = [f1, f2]
    # Convert trivial unary links to in-place form
    if unary_istrivial(g1) && onechild(g1)
        subgraph_factors[1] *= g1.subgraph_factors[1]
        # subgraph_factors[1] *= g1.subgraph_factors[1] * g1.factor
        subgraphs[1] = g1.subgraphs[1]
    end
    if unary_istrivial(g2) && onechild(g2)
        subgraph_factors[2] *= g2.subgraph_factors[1]
        # subgraph_factors[2] *= g2.subgraph_factors[1] * g2.factor
        subgraphs[2] = g2.subgraphs[1]
    end

    if subgraphs[1] == subgraphs[2]
        g = Graph([subgraphs[1]]; subgraph_factors=[prod(subgraph_factors)], operator=Power(2), ftype=F, wtype=W)
    else
        g = Graph(subgraphs; subgraph_factors=subgraph_factors, operator=Prod(), orders=orders(g1), ftype=F, wtype=W)
    end
    return g
end

"""
    multi_product(graphs::Vector{Graph{F,W}}, constants::AbstractVector=ones(F, length(graphs))) where {F,W,C}

    Construct a product graph from multiple input graphs, where each graph can be weighted by a constant. 
    For graphs that are repeated more than once, it adds a power operator to the subgraph to represent the repetition.
    Moreover, it optimizes any trivial unary operators in the resulting product graph.

# Arguments:
- `graphs::Vector{Graph{F,W}}`: A vector of input graphs to be multiplied.
- `constants::AbstractVector`: A vector of scalar multiples. If not provided, it defaults to a vector of ones of the same length as `graphs`.

Returns:
- A new product graph with the unique subgraphs (or powered versions thereof) and the associated constants as subgraph factors.

# Example:
    Given graphs `g1`, `g2`, `g1` and constants `c1`, `c2`, `c3`, the function computes `(c1*c3)*(g1)^2 * c2*g2`.
"""
function multi_product(graphs::Vector{Graph{F,W}}, constants::AbstractVector=ones(F, length(graphs))) where {F,W}
    @assert alleq(orders.(graphs)) "Graphs do not all have the same order."
    g1 = graphs[1]
    subgraphs = graphs
    subgraph_factors = eltype(constants) == F ? constants : Vector{F}(constants)
    # Convert trivial unary links to in-place form
    for (i, sub_g) in enumerate(graphs)
        if unary_istrivial(sub_g) && onechild(sub_g)
            subgraph_factors[i] *= sub_g.subgraph_factors[1]
            # subgraph_factors[i] *= sub_g.subgraph_factors[1] * sub_g.factor
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

    if isempty(unique_graphs)
        return nothing
    end

    if length(unique_factors) == 1
        g = Graph(unique_graphs; subgraph_factors=unique_factors, operator=Power(repeated_counts[1]), orders=orders(g1), ftype=F, wtype=W)
    else
        subgraphs = Vector{Graph{F,W}}()
        for (idx, g) in enumerate(unique_graphs)
            if repeated_counts[idx] == 1
                push!(subgraphs, g)
            else
                push!(subgraphs, Graph([g], operator=Power(repeated_counts[idx]), orders=orders(g1), ftype=F, wtype=W))
            end
        end
        g = Graph(subgraphs; subgraph_factors=unique_factors, operator=Prod(), orders=orders(g1), ftype=F, wtype=W)
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
