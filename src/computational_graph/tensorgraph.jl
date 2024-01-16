struct EinSum <: AbstractOperator
    input_axes1::Vector{Int}
    input_axes2::Vector{Int}
    output_axes::Vector{Int}
    function EinSum(input_axes1::Vector{Int}, input_axes2::Vector{Int}, output_axes::Vector{Int})
        new(input_axes1, input_axes2, output_axes)
    end
end



"""
    mutable struct TensorGraph{F<:Number,W}
    
    A representation of a TensorGraph, e.g., an expression tree of tensors, with type stable node data.

# Members:
- `id::Int`  the unique hash id to identify the diagram
- `name::Symbol`  name of the diagram
- `orders::Vector{Int}`  orders associated with the TensorGraph, e.g., derivative orders
- `subgraphs::Vector{TensorGraph{F,W}}`  vector of sub-diagrams 
- `subgraph_factors::Vector{F}`  scalar multiplicative factors associated with each subgraph. Note that the subgraph factors may be manipulated algebraically. To associate a fixed multiplicative factor with this graph which carries some semantic meaning, use the `factor` argument instead.
- `operator::DataType`  node operation. Addition and multiplication are natively supported via operators Sum and Prod, respectively. Should be a concrete subtype of `AbstractOperator`.
- `weight::Vecctor{W}`  the weight of this node.  Should be a tensor that has the 
- `dims::Vecctor{Int}` the shape of this tensor node. Default value is [1].
- `properties::Any` extra information of Green's functions. Default value is nothing.

# Example:
```julia-repl
julia> g1 = TensorGraph([])
1=0.0

julia> g = TensorGraph([g1, g2]; operator=ComputationalGraphs.Sum())
3=0.0=⨁ (1,2)
```
"""
mutable struct TensorGraph{F<:Number,W} <: AbstractGraph # Graph
    id::Int
    name::String # "" by default
    orders::Vector{Int}

    subgraphs::Vector{TensorGraph{F,W}}
    subgraph_factors::Vector{F}

    operator::DataType
    dims::Vector{Int}
    weight::Array{W}
    properties::Any
    """
        function TensorGraph(subgraphs::AbstractVector; name="", operator::AbstractOperator=Sum(),
            ftype=_dtype.factor, wtype=_dtype.weight, weight=zero(wtype))
        
        Create a TensorGraph struct from a set of subgraphs with the specified node data and operation.

    # Arguments:
    - `subgraphs`  vector of sub-diagrams 
    - `subgraph_factors`  scalar multiplicative factors associated with each subgraph. Note that the subgraph factors may be manipulated algebraically. To associate a fixed multiplicative factor with this graph which carries some semantic meaning, use the `factor` argument instead.
    - `name`  name of the diagram
    - `orders`  orders associated with the graph, e.g., derivative orders
    - `operator`  node operation, i.e., Sum, Prod, or a user-defined operator `Op <: AbstractOperator`
    - `ftype`  typeof(factor)
    - `wtype`  typeof(weight)
    - `weight`  the weight of this node
    - `dims::Vecctor{Int}` the shape of this tensor node. Default value is [1].
    - `properties::Any` extra information of Green's functions. Default value is nothing.
    """
    function TensorGraph(subgraphs::AbstractVector; factor=one(_dtype.factor), subgraph_factors=one.(eachindex(subgraphs)), name="", operator::AbstractOperator=Sum(),
        orders=zeros(Int, 16), ftype=_dtype.factor, wtype=_dtype.weight, dims=[1], weight=zeros(wtype, dims...), properties=nothing
    )

        g = new{ftype,wtype}(uid(), name, orders, subgraphs, subgraph_factors, typeof(operator), dims, weight, properties)
        if factor ≈ one(ftype)
            return g
        else
            return new{ftype,wtype}(uid(), name, orders, [g,], [factor,], Prod, dims, weight * factor, properties)
        end
    end
end

### AbstractGraph interface for Graph ###

# Getters
id(g::TensorGraph) = g.id
name(g::TensorGraph) = g.name
orders(g::TensorGraph) = g.orders
operator(g::TensorGraph) = g.operator
weight(g::TensorGraph) = g.weight
properties(g::TensorGraph) = g.properties
subgraph(g::TensorGraph, i=1) = g.subgraphs[i]
subgraphs(g::TensorGraph) = g.subgraphs
subgraphs(g::TensorGraph, indices::AbstractVector{Int}) = g.subgraphs[indices]
subgraph_factor(g::TensorGraph, i=1) = g.subgraph_factors[i]
subgraph_factors(g::TensorGraph) = g.subgraph_factors
subgraph_factors(g::TensorGraph, indices::AbstractVector{Int}) = g.subgraph_factors[indices]
dims(g::TensorGraph) = g.dims
# Setters
set_id!(g::TensorGraph, id::Int) = (g.id = id)
set_name!(g::TensorGraph, name::String) = (g.name = name)
set_orders!(g::TensorGraph, orders::Vector{Int}) = (g.orders = orders)
set_operator!(g::TensorGraph, operator::Type{<:AbstractOperator}) = (g.operator = operator)
set_operator!(g::TensorGraph, operator::AbstractOperator) = (g.operator = typeof(operator))
set_weight!(g::TensorGraph{F,W}, weight) where {F,W} = (g.weight = W(weight))
set_properties!(g::TensorGraph, properties) = (g.properties = properties)
set_subgraph!(g::TensorGraph{F,W}, subgraph::TensorGraph{F,W}, i=1) where {F,W} = (g.subgraphs[i] = subgraph)
set_subgraphs!(g::TensorGraph{F,W}, subgraphs::Vector{TensorGraph{F,W}}) where {F,W} = (g.subgraphs = subgraphs)
set_subgraphs!(g::TensorGraph{F,W}, subgraphs::Vector{TensorGraph{F,W}}, indices::AbstractVector{Int}) where {F,W} = (g.subgraphs[indices] = subgraphs)
set_subgraph_factor!(g::TensorGraph{F,W}, subgraph_factor, i=1) where {F,W} = (g.subgraph_factors[i] = F(subgraph_factor))
set_subgraph_factors!(g::TensorGraph{F,W}, subgraph_factors::AbstractVector) where {F,W} = (g.subgraph_factors = Vector{F}(subgraph_factors))
set_subgraph_factors!(g::TensorGraph{F,W}, subgraph_factors::AbstractVector, indices::AbstractVector{Int}) where {F,W} = (g.subgraph_factors[indices] = Vector{F}(subgraph_factors))



#Internal function that derives the output axes label based on the input, following the tensordot contraction rule.
function derive_output_axes(input_axes1::Vector{Int}, input_axes2::Vector{Int})
    result = Vector{Int}()
    copy2 = copy(input_axes2)
    for (idx1, label) in enumerate(input_axes1)
        idx2 = findfirst(x -> x == label, copy2)
        if !isnothing(idx2)
            deleteat!(copy2, idx2)
        else
            push!(result, label)
        end
    end

    return vcat(result, copy2)
end

#Internal function that generates the dimension of output tenosr in an Einstein summation.
function dims_contraction(dims1::Vector{Int}, input_axes1::Vector{Int}, dims2::Vector{Int}, input_axes2::Vector{Int}, output_axes::Vector{Int})
    @assert length(dims1) == length(input_axes1) "Einsum subscript does not have  the right number of axes for operand 1. "
    @assert length(dims2) == length(input_axes2) "Einsum subscript does not have  the right number of axes for operand 2. "
    visited = [false for _ in input_axes2]
    result = zeros(Int, length(output_axes))
    count = 0
    for (idx1, label) in enumerate(input_axes1)
        idx2 = findfirst(x -> x == label, input_axes2)
        if !isnothing(idx2)
            @assert dims1[idx1] == dims2[idx2] "Axis[$(idx1)] in first operand do not have same length as axis[$(idx2)] in second operand for contraction."
            visited[idx2] = true
        end

        idx = findfirst(x -> x == label, output_axes)
        if !isnothing(idx)
            result[idx] = dims1[idx1]
            count += 1
        end
    end

    for (idx2, label) in enumerate(input_axes2)
        if visited[idx2]
            continue
        else
            idx = findfirst(x -> x == label, output_axes)
            if !isnothing(idx)
                result[idx] = dims2[idx2]
                count += 1
            end
        end
    end
    @assert count == length(result) "Output has wrong number of axes."
    return result
end

"""
    function Base.:*(g1::TensorGraph{F,W}, c2) where {F,W}

    Returns a graph representing the scalar multiplication `g1*c2`.

# Arguments:
- `g1`  tensor graph
- `c2`  scalar multiple
"""
function Base.:*(g1::TensorGraph{F,W}, c2::Number) where {F,W}
    g = TensorGraph([g1,]; subgraph_factors=[F(c2),], operator=Prod(), orders=orders(g1), ftype=F, wtype=W, dims=g1.dims)
    # Convert trivial unary link to in-place form
    if unary_istrivial(g1) && onechild(g1)
        g.subgraph_factors[1] *= g1.subgraph_factors[1]
        # g.subgraph_factors[1] *= g1.subgraph_factors[1] * g1.factor
        g.subgraphs = g1.subgraphs
    end

    return g
end

"""
    function Base.:*(c1, g2::TensorGraph{F,W}) where {F,W}

    Returns a graph representing the scalar multiplication `c1*g2`.

# Arguments:
- `c1`  scalar multiple
- `g2`  tensor graph
"""
function Base.:*(c1::Number, g2::TensorGraph{F,W}) where {F,W}
    g = TensorGraph([g2,]; subgraph_factors=[F(c1),], operator=Prod(), orders=orders(g2), ftype=F, wtype=W, dims=g2.dims)
    # Convert trivial unary link to in-place form
    if unary_istrivial(g2) && onechild(g2)
        g.subgraph_factors[1] *= g2.subgraph_factors[1]
        # g.subgraph_factors[1] *= g2.subgraph_factors[1] * g2.factor
        g.subgraphs = g2.subgraphs
    end
    return g
end

"""
    function linear_combination(g1::TensorGraph{F,W}, g2::TensorGraph{F,W}, c1, c2) where {F,W}

    Returns a graph representing the linear combination `c1*g1 + c2*g2`.
    If `g1 == g2`, it will return a graph representing `(c1+c2)*g1`.
    TensorGraphs `g1` and `g2` must have the same orders.

# Arguments:
- `g1`  first tensor graph
- `g2`  second tensor graph
- `c1`  first scalar multiple
- `c2`  second scalar multiple
"""
function linear_combination(g1::TensorGraph{F,W}, g2::TensorGraph{F,W}, c1::Number=F(1), c2::Number=F(1)) where {F,W}
    @assert dims(g1) == dims(g2) "g1 and g2 have different dimensions."
    if length(g1.orders) > length(g2.orders)
        g2.orders = [orders(g2); zeros(Int, length(g1.orders) - length(g2.orders))]
    else
        g1.orders = [orders(g1); zeros(Int, length(g2.orders) - length(g1.orders))]
    end
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
        g = TensorGraph([subgraphs[1]]; subgraph_factors=[sum(subgraph_factors)], operator=Sum(), orders=orders(g1), ftype=F, wtype=W, dims=subgraphs[1].dims)
    else
        g = TensorGraph(subgraphs; subgraph_factors=subgraph_factors, operator=Sum(), orders=orders(g1), ftype=F, wtype=W, dims=subgraphs[1].dims)
    end

    return g
end

"""
    function linear_combination(graphs::Vector{TensorGraph{F,W}}, constants::AbstractVector=ones(F, length(graphs))) where {F,W}

    Given a vector 𝐠 of graphs and an equally-sized vector 𝐜 of constants, returns a new
    graph representing the linear combination (𝐜 ⋅ 𝐠). 
    The function identifies unique graphs from the input `graphs` and sums their associated `constants`.
    All input graphs must have the same orders.

# Arguments:
- `graphs`  vector of tensor graphs
- `constants`  vector of scalar multiples (defaults to ones(F, length(graphs))).

# Returns:
- A new `TensorGraph{F,W}` object representing the linear combination of the unique input `graphs` weighted by the constants, 
where duplicate graphs in the input `graphs` are combined by summing their associated constants. 

# Example:
    Given graphs `g1`, `g2`, `g1` and constants `c1`, `c2`, `c3`, the function computes `(c1+c3)*g1 + c2*g2`.
"""
function linear_combination(graphs::Vector{TensorGraph{F,W}}, constants::AbstractVector=ones(F, length(graphs))) where {F,W}
    @assert alleq(dims.(graphs)) "TensorGraphs do not all have the same dimesnions."
    maxlen_orders = maximum(length.(orders.(graphs)))
    for g in graphs
        g.orders = [orders(g); zeros(Int, maxlen_orders - length(orders(g)))]
    end
    @assert alleq(orders.(graphs)) "TensorGraphs do not all have the same order."

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

    unique_graphs = TensorGraph{F,W}[]
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
    g = TensorGraph(unique_graphs; subgraph_factors=unique_factors, operator=Sum(), orders=orders(graphs[1]), ftype=F, wtype=W, dims=unique_graphs.dims)
    return g
end

"""
    function Base.:+(g1::TensorGraph{F,W}, g2::TensorGraph{F,W}) where {F,W}

    Returns a graph `g1 + g2` representing the addition of `g2` with `g1`.
    TensorGraphs `g1` and `g2` must have the same orders.

# Arguments:
- `g1`  first tensor graph
- `g2`  second tensor graph
"""
function Base.:+(g1::TensorGraph{F,W}, g2::TensorGraph{F,W}) where {F,W}
    return linear_combination(g1, g2, F(1), F(1))
end

"""
    function Base.:-(g1::TensorGraph{F,W}, g2::TensorGraph{F,W}) where {F,W}

    Returns a graph `g1 - g2` representing the subtraction of `g2` from `g1`.
    TensorGraphs `g1` and `g2` must have the same orders.

# Arguments:
- `g1`  first tensor graph
- `g2`  second tensor graph
"""
function Base.:-(g1::TensorGraph{F,W}, g2::TensorGraph{F,W}) where {F,W}
    return linear_combination(g1, g2, F(1), F(-1))
end

"""
function einsum(g1::TensorGraph{F,W}, input_axes1::Vector{Int}, g2::TensorGraph{F,W}, input_axes2::Vector{Int}, output_axes::Vector{Int}) where {F,W}

    Returns the Einstein summation of two tensors. If the output axes are not specified, the function falls back to the tensor dot contraction rule.

# Arguments:
- `g1`:  first tensor graph
- `g2`:  second tensor graph
- `input_axes1`:  axes of first tensor. Should be an integer array, with each axis labeled by a unique integer.
- `input_axes2`:  axes of second tensor. 
- `input_axes2`:  axes of output tensor.
"""
function einsum(g1::TensorGraph{F,W}, input_axes1::Vector{Int}, g2::TensorGraph{F,W}, input_axes2::Vector{Int}, output_axes::Vector{Int}) where {F,W}
    # @assert orders(g1) == orders(g2) "g1 and g2 have different orders."
    subgraphs = [g1, g2]
    subgraph_factors = [one(F), one(F)]
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


    if length(g1.orders) > length(g2.orders)
        g2.orders = [orders(g2); zeros(Int, length(g1.orders) - length(g2.orders))]
    else
        g1.orders = [orders(g1); zeros(Int, length(g2.orders) - length(g1.orders))]
    end
    dims = dims_contraction(subgraphs[1].dims, input_axes1, subgraphs[2].dims, input_axes2, output_axes)
    g = TensorGraph(subgraphs; subgraph_factors=subgraph_factors, operator=EinSum(input_axes1, input_axes2, output_axes), orders=orders(g1) + orders(g2), ftype=F, wtype=W, dims=dims)

    return g
end

function einsum(g1::TensorGraph{F,W}, input_axes1::Vector{Int}, g2::TensorGraph{F,W}, input_axes2::Vector{Int}) where {F,W}
    output_axes = derive_output_axes(input_axes1, input_axes2)
    return einsum(g1, input_axes1, g2, input_axes2, output_axes)
end
# """
#     multi_product(graphs::Vector{TensorGraph{F,W}}, constants::AbstractVector=ones(F, length(graphs))) where {F,W,C}

#     Construct a product graph from multiple input graphs, where each graph can be weighted by a constant. 
#     For graphs that are repeated more than once, it adds a power operator to the subgraph to represent the repetition.
#     Moreover, it optimizes any trivial unary operators in the resulting product graph.

# # Arguments:
# - `graphs::Vector{TensorGraph{F,W}}`: A vector of input graphs to be multiplied.
# - `constants::AbstractVector`: A vector of scalar multiples. If not provided, it defaults to a vector of ones of the same length as `graphs`.

# Returns:
# - A new product graph with the unique subgraphs (or powered versions thereof) and the associated constants as subgraph factors.

# # Example:
#     Given graphs `g1`, `g2`, `g1` and constants `c1`, `c2`, `c3`, the function computes `(c1*c3)*(g1)^2 * c2*g2`.
# """
# function multi_product(graphs::Vector{TensorGraph{F,W}}, constants::AbstractVector=ones(F, length(graphs))) where {F,W}
#     # @assert alleq(orders.(graphs)) "Graphs do not all have the same order."
#     g1 = graphs[1]
#     subgraphs = graphs
#     subgraph_factors = eltype(constants) == F ? constants : Vector{F}(constants)

#     maxlen_orders = maximum(length.(orders.(graphs)))
#     g_orders = zeros(Int, maxlen_orders)
#     # Convert trivial unary links to in-place form
#     for (i, sub_g) in enumerate(graphs)
#         if unary_istrivial(sub_g) && onechild(sub_g)
#             subgraph_factors[i] *= sub_g.subgraph_factors[1]
#             # subgraph_factors[i] *= sub_g.subgraph_factors[1] * sub_g.factor
#             subgraphs[i] = sub_g.subgraphs[1]
#         end
#         sub_g.orders = [orders(sub_g); zeros(Int, maxlen_orders - length(orders(sub_g)))]
#         g_orders += orders(sub_g)
#     end

#     unique_graphs = Vector{TensorGraph{F,W}}()
#     unique_factors = F[]
#     repeated_counts = Int[]
#     for (idx, g) in enumerate(subgraphs)
#         loc = findfirst(isequal(g), unique_graphs)
#         if isnothing(loc)
#             push!(unique_graphs, g)
#             push!(unique_factors, subgraph_factors[idx])
#             push!(repeated_counts, 1)
#         else
#             unique_factors[loc] *= subgraph_factors[idx]
#             repeated_counts[loc] += 1
#         end
#     end

#     if isempty(unique_graphs)
#         return nothing
#     end

#     if length(unique_factors) == 1
#         g = TensorGraph(unique_graphs; subgraph_factors=unique_factors, operator=Power(repeated_counts[1]), orders=g_orders, ftype=F, wtype=W, dims=dims_contraction([uni_graph.dims for uni_graph in unique_graphs]))
#     else
#         subgraphs = Vector{TensorGraph{F,W}}()
#         for (idx, g) in enumerate(unique_graphs)
#             if repeated_counts[idx] == 1
#                 push!(subgraphs, g)
#             else
#                 push!(subgraphs, TensorGraph([g], operator=Power(repeated_counts[idx]), orders=orders(g1) * repeated_counts[idx], ftype=F, wtype=W))
#             end
#         end
#         g = TensorGraph(subgraphs; subgraph_factors=unique_factors, operator=Prod(), orders=g_orders, ftype=F, wtype=W, dims=dims_contraction([subgraph.dims for subgraph in subgraphs]))
#     end
#     return g
# end

# """
#     function Base.:*(g1::TensorGraph{F,W}, g2::TensorGraph{F,W}) where {F,W}

#     Returns a graph `g1 * g2` representing the graph product between `g1` and `g2`.

# # Arguments:
# - `g1`  first tensor graph
# - `g2`  second tensor graph
# """
# function Base.:*(g1::TensorGraph{F,W}, g2::TensorGraph{F,W}) where {F,W}
#     return multi_product(g1, g2)
# end
