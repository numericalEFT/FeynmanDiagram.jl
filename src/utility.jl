module Utility
using ..ComputationalGraphs
#using ..ComputationalGraphs: Sum, Prod, Power, decrement_power
using ..ComputationalGraphs: decrement_power
using ..ComputationalGraphs: build_all_leaf_derivative, eval!, isfermionic
import ..ComputationalGraphs: count_operation
using ..ComputationalGraphs.AbstractTrees
using ..DiagTree
using ..DiagTree: Diagram, PropagatorId, BareGreenId, BareInteractionId
using ..FrontEnds: LabelProduct
using ..Taylor

@inline apply(::Type{ComputationalGraphs.Sum}, diags::Vector{T}, factors::Vector{F}) where {T<:TaylorSeries,F<:Number} = sum(d * f for (d, f) in zip(diags, factors))
@inline apply(::Type{ComputationalGraphs.Prod}, diags::Vector{T}, factors::Vector{F}) where {T<:TaylorSeries,F<:Number} = prod(d * f for (d, f) in zip(diags, factors))
@inline apply(::Type{ComputationalGraphs.Power{N}}, diags::Vector{T}, factors::Vector{F}) where {N,T<:TaylorSeries,F<:Number} = (diags[1])^N * factors[1]

@inline apply(::Type{DiagTree.Sum}, diags::Vector{T}, factors::Vector{F}) where {T<:TaylorSeries,F<:Number} = sum(d * f for (d, f) in zip(diags, factors))
@inline apply(::Type{DiagTree.Prod}, diags::Vector{T}, factors::Vector{F}) where {T<:TaylorSeries,F<:Number} = prod(d * f for (d, f) in zip(diags, factors))


"""
    function taylorexpansion!(graph::G, var_dependence::Dict{Int,Vector{Bool}}=Dict{Int,Vector{Bool}}(); to_coeff_map::Dict{Int,TaylorSeries{G}}=Dict{Int,TaylorSeries{G}}()) where {G<:Graph}
    
    Return a taylor series of graph g, together with a map of between nodes of g and correponding taylor series.
# Arguments:
- `graph`  Target graph 
- `var_dependence::Dict{Int,Vector{Bool}}` A dictionary that specifies the variable dependence of target graph leaves. Should map the id of each leaf to a Bool vector. 
    The length of the vector should be the same as number of variables.
- `to_coeff_map::Dict{Int,TaylorSeries}` A dicitonary that maps id of each node of target graph to its correponding taylor series.
"""
function taylorexpansion!(graph::G, var_dependence::Dict{Int,Vector{Bool}}=Dict{Int,Vector{Bool}}(); to_coeff_map::Dict{Int,TaylorSeries{G}}=Dict{Int,TaylorSeries{G}}()) where {G<:Graph}
    if haskey(to_coeff_map, graph.id) #If already exist, use taylor series in to_coeff_map.
        return to_coeff_map[graph.id], to_coeff_map

    elseif isleaf(graph)
        if haskey(var_dependence, graph.id)
            var = var_dependence[graph.id]
        else
            var = fill(false, get_numvars()) #if dependence not provided, assume the graph depends on no variables
        end
        ordtuple = ((var[idx]) ? (0:get_orders(idx)) : (0:0) for idx in 1:get_numvars())
        result = TaylorSeries{G}()
        for order in collect(Iterators.product(ordtuple...)) #varidx specifies the variables graph depends on. Iterate over all taylor coefficients of those variables.
            o = collect(order)
            if sum(o) == 0      # For a graph the zero order taylor coefficient is just itself.
                result.coeffs[o] = graph
            else
                coeff = Graph([]; operator=ComputationalGraphs.Sum(), factor=graph.factor)
                result.coeffs[o] = coeff
            end
        end
        to_coeff_map[graph.id] = result
        return result, to_coeff_map
    else
        to_coeff_map[graph.id] = graph.factor * apply(graph.operator, [taylorexpansion!(sub, var_dependence; to_coeff_map=to_coeff_map)[1] for sub in graph.subgraphs], graph.subgraph_factors)
        return to_coeff_map[graph.id], to_coeff_map
    end
end

"""
    function taylorexpansion!(graph::FeynmanGraph{F,W}, var_dependence::Dict{Int,Vector{Bool}}=Dict{Int,Vector{Bool}}(); to_coeff_map::Dict{Int,TaylorSeries{G}}=Dict{Int,TaylorSeries{G}}()) where {F,W}
    
    Return a taylor series of FeynmanGraph g, together with a map of between nodes of g and correponding taylor series.
# Arguments:
- `graph`  Target FeynmanGraph 
- `var_dependence::Dict{Int,Vector{Bool}}` A dictionary that specifies the variable dependence of target graph leaves. Should map the id of each leaf to a Bool vector. 
    The length of the vector should be the same as number of variables.
- `to_coeff_map::Dict{Int,TaylorSeries}` A dicitonary that maps id of each node of target graph to its correponding taylor series.
- `from_coeff_map::Dict{Int,Tuple{Int,Vector{Bool}}}` A dicitonary that maps a taylor coefficient to its owner FeynmanGraph. The key should be the id of coefficient graph, and value should be a tuple of (feynmangraph.id, order).
"""
function taylorexpansion!(graph::FeynmanGraph{F,W}, var_dependence::Dict{Int,Vector{Bool}}=Dict{Int,Vector{Bool}}(); to_coeff_map::Dict{Int,TaylorSeries{Graph{F,W}}}=Dict{Int,TaylorSeries{Graph{F,W}}}(), from_coeff_map::Dict{Int,Tuple{Int,Vector{Int}}}=Dict{Int,Tuple{Int,Vector{Int}}}()) where {F,W}
    if haskey(to_coeff_map, graph.id) #If already exist, use taylor series in to_coeff_map.
        if isleaf(graph)
            for (order, coeff) in to_coeff_map[graph.id].coeffs
                if haskey(from_coeff_map, coeff.id)
                    @assert from_coeff_map[coeff.id] == (graph.id, order) "The graph g$(graph.id) is mapped to two different leaf taylor series!"
                else
                    from_coeff_map[coeff.id] = (graph.id, order)
                end
            end
        end
        return to_coeff_map[graph.id], to_coeff_map, from_coeff_map

    elseif isleaf(graph)
        if haskey(var_dependence, graph.id)
            var = var_dependence[graph.id]
        else
            var = fill(false, get_numvars()) #if dependence not provided, assume the graph depends on no variables
        end
        ordtuple = ((var[idx]) ? (0:get_orders(idx)) : (0:0) for idx in 1:get_numvars())
        result = TaylorSeries{Graph{F,W}}()
        for order in collect(Iterators.product(ordtuple...)) #varidx specifies the variables graph depends on. Iterate over all taylor coefficients of those variables.
            o = collect(order)
            coeff = Graph([]; operator=ComputationalGraphs.Sum(), factor=graph.factor)
            result.coeffs[o] = coeff
            from_coeff_map[coeff.id] = (graph.id, o)
        end
        to_coeff_map[graph.id] = result
        return result, to_coeff_map, from_coeff_map
    else
        to_coeff_map[graph.id] = graph.factor * apply(graph.operator, [taylorexpansion!(sub, var_dependence; to_coeff_map=to_coeff_map, from_coeff_map=from_coeff_map)[1] for sub in graph.subgraphs], graph.subgraph_factors)
        return to_coeff_map[graph.id], to_coeff_map, from_coeff_map
    end
end

"""
    function taylorexpansion!(graph::Diagram{W}, var_dependence::Dict{Int,Vector{Bool}}=Dict{Int,Vector{Bool}}(); to_coeff_map::Dict{Int,TaylorSeries{G}}=Dict{Int,TaylorSeries{G}}()) where {W}
    
    Return a taylor series of Diagram g, together with a map of between nodes of g and correponding taylor series.
# Arguments:
- `graph`  Target diagram 
- `var_dependence::Dict{Int,Vector{Bool}}` A dictionary that specifies the variable dependence of target diagram leaves. Should map the id of each leaf to a Bool vector. 
    The length of the vector should be the same as number of variables.
- `to_coeff_map::Dict{Int,TaylorSeries}` A dicitonary that maps id of each node of target diagram to its correponding taylor series.
"""
function taylorexpansion!(graph::Diagram{W}, var_dependence::Dict{Int,Vector{Bool}}=Dict{Int,Vector{Bool}}(); to_coeff_map::Dict{Int,TaylorSeries{Graph{W,W}}}=Dict{Int,TaylorSeries{Graph{W,W}}}()) where {W}
    if haskey(to_coeff_map, graph.hash) #If already exist, use taylor series in to_coeff_map.
        return to_coeff_map[graph.hash], to_coeff_map

    elseif isempty(graph.subdiagram)
        if haskey(var_dependence, graph.hash)
            var = var_dependence[graph.hash]
        else
            var = fill(false, get_numvars()) #if dependence not provhashed, assume the graph depends on no variables
        end
        ordtuple = ((var[idx]) ? (0:get_orders(idx)) : (0:0) for idx in 1:get_numvars())
        result = TaylorSeries{Graph{W,W}}()
        for order in collect(Iterators.product(ordtuple...)) #varidx specifies the variables graph depends on. Iterate over all taylor coefficients of those variables.
            o = collect(order)
            coeff = Graph([]; operator=ComputationalGraphs.Sum(), factor=graph.factor)
            result.coeffs[o] = coeff
        end
        to_coeff_map[graph.hash] = result
        return result, to_coeff_map
    else
        to_coeff_map[graph.hash] = graph.factor * apply(typeof(graph.operator), [taylorexpansion!(sub, var_dependence; to_coeff_map=to_coeff_map)[1] for sub in graph.subdiagram], ones(W, length(graph.subdiagram)))
        return to_coeff_map[graph.hash], to_coeff_map
    end
end

"""
    function taylorexpansion!(graph::FeynmanGraph{F,W}, propagator_var::Tuple{Vector{Bool},Vector{Bool}}, label::Tuple{LabelProduct,LabelProduct}; to_coeff_map::Dict{Int,TaylorSeries{Graph{F,W}}}=Dict{Int,TaylorSeries{Graph{F,W}}}()) where {F,W}
    
    Return a taylor series of FeynmanGraph g, together with a map of between nodes of g and correponding taylor series. In this set up, the leaves that are the same type of propagators (fermi or bose) depend on the same set of variables, 
    whereas other types of Feynman diagrams (such as vertices) depends on no variables that need to be differentiated (for AD purpose, they are just constants).
# Arguments:
- `graph`  Target FeynmanGraph
- `propagator_var::Tuple{Vector{Bool},Vector{Bool}}` A Tuple that specifies the variable dependence of fermi (first element) and bose propagator (second element).
    The dependence is given by a vector of the length same as the number of variables.
- `label::Tuple{LabelProduct,LabelProduct}` A Tuple fermi (first element) and bose LabelProduct (second element).
- `to_coeff_map::Dict{Int,TaylorSeries}` A dicitonary that maps id of each node of target diagram to its correponding taylor series.
"""
function taylorexpansion!(graph::FeynmanGraph{F,W}, propagator_var::Tuple{Vector{Bool},Vector{Bool}}; to_coeff_map::Dict{Int,TaylorSeries{Graph{F,W}}}=Dict{Int,TaylorSeries{Graph{F,W}}}(), from_coeff_map::Dict{Int,Tuple{Int,Vector{Int}}}=Dict{Int,Tuple{Int,Vector{Int}}}()) where {F,W}
    var_dependence = Dict{Int,Vector{Bool}}()
    for leaf in Leaves(graph)
        if ComputationalGraphs.diagram_type(leaf) == ComputationalGraphs.Propagator
            # In = leaf.properties.vertices[2][1].label
            if isfermionic(leaf.properties.vertices[1])
                # if label[1][In][2] >= 0 #For fake propagator, this label is smaller than zero, and those propagators should not be differentiated.
                var_dependence[leaf.id] = [propagator_var[1][idx] ? true : false for idx in 1:get_numvars()]
                # end
            else
                var_dependence[leaf.id] = [propagator_var[2][idx] ? true : false for idx in 1:get_numvars()]
            end
        end
    end
    return taylorexpansion!(graph, var_dependence; to_coeff_map=to_coeff_map, from_coeff_map=from_coeff_map)
end

"""
    function taylorexpansion!(graph::Diagram{W}, propagator_var::Dict{DataType,Vector{Bool}}; to_coeff_map::Dict{Int,TaylorSeries{Graph{W,W}}}=Dict{Int,TaylorSeries{Graph{W,W}}}()) where {W}
    
    Return a taylor series of Diagram g, together with a map of between nodes of g and correponding taylor series. In this set up, the leaves that are the same type of diagrams (such as Green functions) depend on the same set of variables.
    
# Arguments:
- `graph`  Target Diagram
- `propagator_var::Dict{DataType,Vector{Bool}}` A dictionary that specifies the variable dependence of different types of diagrams. Should be a map between DataTypes in DiagramID and Bool vectors.
    The dependence is given by a vector of the length same as the number of variables.
- `to_coeff_map::Dict{Int,TaylorSeries}` A dicitonary that maps id of each node of target diagram to its correponding taylor series.
"""
function taylorexpansion!(graph::Diagram{W}, propagator_var::Dict{DataType,Vector{Bool}}; to_coeff_map::Dict{Int,TaylorSeries{Graph{W,W}}}=Dict{Int,TaylorSeries{Graph{W,W}}}()) where {W}
    var_dependence = Dict{Int,Vector{Bool}}()
    for leaf in Leaves(graph)
        if haskey(propagator_var, typeof(leaf.id))
            var_dependence[leaf.hash] = [propagator_var[typeof(leaf.id)][idx] ? true : false for idx in 1:get_numvars()]
        end
    end
    return taylorexpansion!(graph, var_dependence; to_coeff_map=to_coeff_map)
end

function taylorexpansion!(graphs::Vector{G}, var_dependence::Dict{Int,Vector{Bool}}=Dict{Int,Vector{Bool}}(); to_coeff_map::Dict{Int,TaylorSeries{G}}=Dict{Int,TaylorSeries{G}}()) where {G<:Graph}
    result = Vector{TaylorSeries{G}}()
    for graph in graphs
        taylor, _ = taylorexpansion!(graph, var_dependence; to_coeff_map=to_coeff_map)
        push!(result, taylor)
    end
    return result, to_coeff_map
end

function taylorexpansion!(graphs::Vector{Diagram{W}}, propagator_var::Dict{DataType,Vector{Bool}}; to_coeff_map::Dict{Int,TaylorSeries{Graph{W,W}}}=Dict{Int,TaylorSeries{Graph{W,W}}}()) where {W}
    result = Vector{TaylorSeries{Graph{W,W}}}()
    for graph in graphs
        taylor, _ = taylorexpansion!(graph, propagator_var; to_coeff_map=to_coeff_map)
        push!(result, taylor)
    end
    return result, to_coeff_map
end
"""
    taylorexpansion_withmap(g::G; coeffmode=true, var::Vector{Int}=collect(1:get_numvars())) where {G<:Graph}
    
    Return a taylor series of graph g, together with a map of chain relationships between generated derivatives.
    This function is only internally used for constructing high order derivatives by naive nested forward AD.
    It is only for banch mark purpose and not exported.
# Arguments:
- `g`  Target graph 
- `coeffmode` If true, the generated taylor series saves taylor coefficients with the factorial prefactor. If false, the taylor series saves derivatives instead
- `var` The index of variables graph depends on
"""
function taylorexpansion_withmap(g::G; coeffmode=true, var::Vector{Bool}=fill(true, get_numvars())) where {G<:Graph}
    @assert isleaf(g)
    chainrule_map_leaf = Dict{Int,Dict{Int,G}}()
    maxorder = get_orders()
    current_func = Dict(zeros(Int, get_numvars()) => g)
    result = TaylorSeries{G}()
    result.coeffs[zeros(Int, get_numvars())] = g

    for i in 1:sum(get_orders())
        new_func = Dict{Vector{Int},G}()
        for (order, func) in current_func
            if !haskey(chainrule_map_leaf, func.id)
                chainrule_map_leaf[func.id] = Dict{Int,G}()
            end
            for idx in eachindex(var)
                if var[idx]
                    ordernew = copy(order)
                    ordernew[idx] += 1
                    if ordernew[idx] <= get_orders(idx)
                        if !haskey(result.coeffs, ordernew)
                            if coeffmode
                                funcAD = Graph([]; operator=ComputationalGraphs.Sum(), factor=g.factor)
                            else
                                #funcAD = taylor_factorial(ordernew) * Graph([]; operator=ComputationalGraphs.Sum(), factor=g.factor)
                                funcAD = Graph([]; operator=ComputationalGraphs.Sum(), factor=taylor_factorial(ordernew) * g.factor)
                            end
                            new_func[ordernew] = funcAD
                            result.coeffs[ordernew] = funcAD
                            chainrule_map_leaf[func.id][idx] = funcAD
                        else
                            chainrule_map_leaf[func.id][idx] = result.coeffs[ordernew]
                        end
                    end
                end
            end
        end
        current_func = new_func
    end

    return result, chainrule_map_leaf
end



#Functions below generate high order derivatives with naive nested forward AD. This part would be significantly refactored later with 
# Taylor Series API.

function build_derivative_backAD!(g::G, leaftaylor::Dict{Int,TaylorSeries{G}}=Dict{Int,TaylorSeries{G}}()) where {G<:Graph}
    chainrule_map_leaf = Dict{Int,Dict{Int,G}}()
    for leaf in Leaves(g)
        if !haskey(leaftaylor, leaf.id)
            leaftaylor[leaf.id], map = taylorexpansion_withmap(leaf; coeffmode=false)
            chainrule_map_leaf = merge(chainrule_map_leaf, map)
        end
    end

    leafAD, chainrule_map = build_all_leaf_derivative(g)
    current_func = Dict(zeros(Int, get_numvars()) => g)

    result = TaylorSeries{G}()
    result.coeffs[zeros(Int, get_numvars())] = g
    for i in 1:sum(get_orders())
        new_func = Dict{Vector{Int},G}()
        for (order, func) in current_func
            for idx in 1:get_numvars()
                ordernew = copy(order)
                ordernew[idx] += 1
                if ordernew[idx] <= get_orders(idx)
                    if !haskey(result.coeffs, ordernew)
                        funcAD = forwardAD_taylor(func, idx, chainrule_map, chainrule_map_leaf, leaftaylor)
                        if !isnothing(funcAD)
                            new_func[ordernew] = funcAD
                            result.coeffs[ordernew] = funcAD
                        end
                    end
                end
            end
        end
        current_func = new_func
    end
    return result, leaftaylor
end


function forwardAD_taylor(g::G, varidx::Int, chainrule_map::Dict{Int,Array{G,1}}, chainrule_map_leaf::Dict{Int,Dict{Int,G}}, leaftaylor::Dict{Int,TaylorSeries{G}}) where {G<:Graph}
    # if haskey(chainrule_map, g.id)
    #     return chainrule!(varidx, chainrule_map[g.id], leaftaylor)
    # elseif haskey(chainrule_map_leaf, g.id)
    if haskey(chainrule_map_leaf, g.id)
        map = chainrule_map_leaf[g.id]
        if haskey(map, varidx)
            return map[varidx]
        else
            return nothing
        end
    elseif g.operator == ComputationalGraphs.Sum
        children = Array{G,1}()
        for graph in g.subgraphs
            dgraph = forwardAD_taylor(graph, varidx, chainrule_map, chainrule_map_leaf, leaftaylor)
            if !isnothing(dgraph)
                push!(children, dgraph)
            end
        end
        if isempty(children)
            return nothing
        else
            return linear_combination(children, g.subgraph_factors)
        end
    elseif g.operator == ComputationalGraphs.Prod
        children = Array{G,1}()
        for (i, graph) in enumerate(g.subgraphs)
            dgraph = forwardAD_taylor(graph, varidx, chainrule_map, chainrule_map_leaf, leaftaylor)
            if !isnothing(dgraph)
                subgraphs = [j == i ? dgraph : subg for (j, subg) in enumerate(g.subgraphs)]
                push!(children, Graph(subgraphs; operator=ComputationalGraphs.Prod(), subgraph_factors=g.subgraph_factors))
            end
        end
        if isempty(children)
            return nothing
        else
            return linear_combination(children)
        end
    elseif g.operator <: ComputationalGraphs.Power

        dgraph = forwardAD_taylor(g.subgraphs[1], varidx, chainrule_map, chainrule_map_leaf, leaftaylor)
        if isnothing(dgraph)
            return nothing
        else
            power = eltype(g.operator)
            if power == 1
                return dgraph
            else
                return dgraph * Graph(g.subgraphs; subgraph_factors=power * g.subgraph_factors, operator=decrement_power(g.operator))
            end
        end
    end
end

function chainrule!(varidx::Int, dg::Array{G,1}, leaftaylor::Dict{Int,TaylorSeries{G}}) where {G<:Graph}
    children = Array{G,1}()
    order = zeros(Int, get_numvars())
    order[varidx] += 1
    for i in 1:length(dg)÷2
        taylor = leaftaylor[dg[2*i-1].id]
        if haskey(taylor.coeffs, order)
            coeff = taylor.coeffs[order]
            push!(children, coeff * dg[2*i])
        end
    end
    if isempty(children)
        return nothing
    else
        return linear_combination(children)
    end
end

function count_operation(g::TaylorSeries{G}) where {G<:Graph}
    return count_operation(g.coeffs)
end

function count_operation(graphs::Vector{TaylorSeries{G}}) where {G<:Graph}
    if length(graphs) == 0
        return [0, 0]
    else
        allcoeffs = Vector{G}()
        for g in graphs
            for (order, coeffs) in g.coeffs
                push!(allcoeffs, coeffs)
            end
        end
        return count_operation(allcoeffs)
    end
end

function count_operation(graphs::Vector{TaylorSeries{G}}, order::Vector{Int}) where {G<:Graph}
    if length(graphs) == 0
        return [0, 0]
    else
        allcoeffs = Vector{G}()
        for g in graphs
            push!(allcoeffs, g.coeffs[order])
        end
        return count_operation(allcoeffs)
    end
end

end