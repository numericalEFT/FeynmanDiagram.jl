module Utility
using ..ComputationalGraphs
using ..ComputationalGraphs: Sum, Prod, Power
using ..ComputationalGraphs: build_all_leaf_derivative, eval!
using ..ComputationalGraphs.AbstractTrees
using ..Taylor


"""
    function graph_to_taylor(graph::G, varidx::Array{Int,1}=collect(1:get_numvars())) where {G<:Graph}
"""
function graph_to_taylor(graph::G, varidx::Array{Int,1}=collect(1:get_numvars())) where {G<:Graph}
    @assert isleaf(graph)
    maxorder = get_order()
    ordtuple = ((idx in varidx) ? (0:maxorder) : (0:0) for idx in 1:get_numvars())
    result = TaylorSeries{G}()
    for order in collect(Iterators.product(ordtuple...)) #varidx specifies the variables graph depends on. Iterate over all taylor coefficients of those variables.
        o = collect(order)
        if sum(o) <= get_order()
            coeff = Graph([]; operator=Sum(), factor=graph.factor)
            result.coeffs[o] = coeff
        end
    end
    return result
end


function graph_to_taylor_withmap!(g::G; coeffmode=true, varidx::Array{Int,1}=collect(1:get_numvars()), chainrule_map_leaf::Dict{Int,Dict{Int,G}}=Dict{Int,Dict{Int,G}}()) where {G<:Graph}
    @assert isleaf(g)
    maxorder = get_order()
    current_func = Dict(zeros(Int, get_numvars()) => g)
    result = TaylorSeries{G}()
    result.coeffs[zeros(Int, get_numvars())] = g

    for i in 1:get_order()
        new_func = Dict{Array{Int,1},G}()
        for (order, func) in current_func
            if !haskey(chainrule_map_leaf, func.id)
                chainrule_map_leaf[func.id] = Dict{Int,G}()
            end
            for idx in varidx
                ordernew = copy(order)
                ordernew[idx] += 1
                if !haskey(result.coeffs, ordernew)
                    if coeffmode
                        funcAD = Graph([]; operator=Sum(), factor=g.factor)
                    else
                        #funcAD = taylor_factorial(ordernew) * Graph([]; operator=Sum(), factor=g.factor)
                        funcAD = Graph([]; operator=Sum(), factor=taylor_factorial(ordernew) * g.factor)
                    end
                    new_func[ordernew] = funcAD
                    result.coeffs[ordernew] = funcAD
                    chainrule_map_leaf[func.id][idx] = funcAD
                else
                    chainrule_map_leaf[func.id][idx] = result.coeffs[ordernew]
                end
            end
        end
        current_func = new_func
    end

    return result
end

@inline apply(::Type{Sum}, diags::Vector{T}, factors::Vector{F}) where {T<:TaylorSeries,F<:Number} = sum(d * f for (d, f) in zip(diags, factors))
@inline apply(::Type{Prod}, diags::Vector{T}, factors::Vector{F}) where {T<:TaylorSeries,F<:Number} = prod(d * f for (d, f) in zip(diags, factors))
@inline apply(::Type{Power{N}}, diags::Vector{T}, factors::Vector{F}) where {N,T<:TaylorSeries,F<:Number} = (diags[1])^N * factors[1]

function taylorexpansion!(graph::G, taylormap::Dict{Int,T}=Dict{Int,TaylorSeries{G}}()) where {G<:Graph,T<:TaylorSeries}
    if isempty(taylormap)
        for g in Leaves(graph)
            if !haskey(taylormap, g.id)
                #taylormap[g.id] = graph_to_taylor(g)
                taylormap[g.id] = graph_to_taylor_withmap!(g)
            end
        end
    end
    rootid = -1
    for g in PostOrderDFS(graph) # postorder traversal will visit all subdiagrams of a diagram first
        rootid = g.id
        if isleaf(g) || haskey(taylormap, g.id)
            continue
        end
        taylormap[g.id] = apply(g.operator, [taylormap[sub.id] for sub in g.subgraphs], g.subgraph_factors)
    end
    return taylormap[rootid]
end

function build_derivative_backAD!(g::G, leaftaylor::Dict{Int,TaylorSeries{G}}=Dict{Int,TaylorSeries{G}}()) where {G<:Graph}
    chainrule_map_leaf = Dict{Int,Dict{Int,G}}()
    for leaf in Leaves(g)
        if !haskey(leaftaylor, leaf.id)
            leaftaylor[leaf.id] = graph_to_taylor_withmap!(leaf; coeffmode=false, chainrule_map_leaf=chainrule_map_leaf)
        end
    end
    for (i, leaf) in chainrule_map_leaf
        for (j, coeff) in leaf
            print("test $(i) $(j) $(coeff) $(eval!(coeff))\n")
        end
        print("\n")
    end
    leafAD, chainrule_map = build_all_leaf_derivative(g)
    current_func = Dict(zeros(Int, get_numvars()) => g)

    result = TaylorSeries{G}()
    result.coeffs[zeros(Int, get_numvars())] = g
    for i in 1:get_order()
        new_func = Dict{Array{Int,1},G}()
        for (order, func) in current_func
            for idx in 1:get_numvars()
                ordernew = copy(order)
                ordernew[idx] += 1
                if !haskey(result.coeffs, ordernew)
                    funcAD = forwardAD_taylor(func, idx, chainrule_map, chainrule_map_leaf, leaftaylor)
                    if !isnothing(funcAD)
                        new_func[ordernew] = funcAD
                        result.coeffs[ordernew] = funcAD
                    end
                end
            end
        end
        current_func = new_func
    end
    return result
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
    elseif g.operator == Sum
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
    elseif g.operator == Prod
        children = Array{G,1}()
        for (i, graph) in enumerate(g.subgraphs)
            dgraph = forwardAD_taylor(graph, varidx, chainrule_map, chainrule_map_leaf, leaftaylor)
            if !isnothing(dgraph)
                subgraphs = [j == i ? dgraph : subg for (j, subg) in enumerate(g.subgraphs)]
                push!(children, Graph(subgraphs; operator=Prod(), subgraph_factors=g.subgraph_factors))
            end
        end
        if isempty(children)
            return nothing
        else
            return linear_combination(children)
        end
    elseif g.operator <: Power

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

end