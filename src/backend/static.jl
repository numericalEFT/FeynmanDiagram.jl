function _to_static(::Type{ComputationalGraphs.Sum}, subgraphs::Vector{Graph{F,W}}, subgraph_factors::Vector{F}) where {F,W}
    if length(subgraphs) == 1
        factor_str = subgraph_factors[1] == 1 ? "" : " * $(subgraph_factors[1])"
        return "(g$(subgraphs[1].id)$factor_str)"
    else
        terms = ["g$(g.id)" * (gfactor == 1 ? "" : " * $gfactor") for (g, gfactor) in zip(subgraphs, subgraph_factors)]
        return "(" * join(terms, " + ") * ")"
    end
end

function _to_static(::Type{ComputationalGraphs.Prod}, subgraphs::Vector{Graph{F,W}}, subgraph_factors::Vector{F}) where {F,W}
    if length(subgraphs) == 1
        factor_str = subgraph_factors[1] == 1 ? "" : " * $(subgraph_factors[1])"
        return "(g$(subgraphs[1].id)$factor_str)"
    else
        terms = ["g$(g.id)" * (gfactor == 1 ? "" : " * $gfactor") for (g, gfactor) in zip(subgraphs, subgraph_factors)]
        return "(" * join(terms, " * ") * ")"
        # return "(" * join(["g$(g.id)" for g in subgraphs], " * ") * ")"
    end
end

function _to_static(::Type{ComputationalGraphs.Power{N}}, subgraphs::Vector{Graph{F,W}}, subgraph_factors::Vector{F}) where {N,F,W}
    factor_str = subgraph_factors[1] == 1 ? "" : " * $(subgraph_factors[1])"
    return "((g$(subgraphs[1].id))^$N$factor_str)"
end

function _to_static(::Type{ComputationalGraphs.Sum}, subgraphs::Vector{FeynmanGraph{F,W}}, subgraph_factors::Vector{F}) where {F,W}
    if length(subgraphs) == 1
        factor_str = subgraph_factors[1] == 1 ? "" : " * $(subgraph_factors[1])"
        return "(g$(subgraphs[1].id)$factor_str)"
    else
        terms = ["g$(g.id)" * (gfactor == 1 ? "" : " * $gfactor") for (g, gfactor) in zip(subgraphs, subgraph_factors)]
        return "(" * join(terms, " + ") * ")"
    end
end

function _to_static(::Type{ComputationalGraphs.Prod}, subgraphs::Vector{FeynmanGraph{F,W}}, subgraph_factors::Vector{F}) where {F,W}
    if length(subgraphs) == 1
        factor_str = subgraph_factors[1] == 1 ? "" : " * $(subgraph_factors[1])"
        return "(g$(subgraphs[1].id)$factor_str)"
    else
        terms = ["g$(g.id)" * (gfactor == 1 ? "" : " * $gfactor") for (g, gfactor) in zip(subgraphs, subgraph_factors)]
        return "(" * join(terms, " * ") * ")"
    end
end

function _to_static(::Type{ComputationalGraphs.Power{N}}, subgraphs::Vector{FeynmanGraph{F,W}}, subgraph_factors::Vector{F}) where {N,F,W}
    factor_str = subgraph_factors[1] == 1 ? "" : " * $(subgraph_factors[1])"
    return "((g$(subgraphs[1].id))^$N$factor_str)"
end

"""
    function to_julia_str(graphs::AbstractVector{G}; root::AbstractVector{Int}=[g.id for g in graphs], name::String="eval_graph!") where {G<:AbstractGraph}
    
Compile a list of graphs into a string for a julia static function. The function takes two arguments: `root` and `leaf`. 
`root` is a vector of the root node ids of the graphs, and `leaf` is a vector of the leaf nodes' weights of the graphs. 
"""
function to_julia_str(graphs::AbstractVector{G}; root::AbstractVector{Int}=[g.id for g in graphs], name::String="eval_graph!") where {G<:AbstractGraph}
    head = "function $name(root::AbstractVector, leaf::AbstractVector)\n "
    body = ""
    leafidx = 1
    inds_visitedleaf = Int[]
    inds_visitednode = Int[]
    for graph in graphs
        for g in PostOrderDFS(graph) #leaf first search
            if g.id in root
                target = "root[$(findfirst(x -> x == g.id, root))]"
            else
                target = "g$(g.id)"
            end
            if isempty(g.subgraphs) #leaf
                g.id in inds_visitedleaf && continue
                factor_str = g.factor == 1 ? "" : " * $(g.factor)"
                body *= "    $target = leaf[$leafidx]$factor_str\n "
                leafidx += 1
                push!(inds_visitedleaf, g.id)
            else
                g.id in inds_visitednode && continue
                factor_str = g.factor == 1 ? "" : " * $(g.factor)"
                body *= "    $target = $(_to_static(g.operator, g.subgraphs, g.subgraph_factors))$factor_str\n "
                push!(inds_visitednode, g.id)
            end
        end
    end
    tail = "end"
    return head * body * tail
end

"""
    function to_julia_str(graphs::AbstractVector{G}, leafMap::Dict{Int,Int}; root::AbstractVector{Int}=[g.id for g in graphs],
        name::String="eval_graph!") where {G<:AbstractGraph}
    
Compile a list of Feynman graphs into a string for a julia static function. The complied function takes two arguments: `root` and `leafVal`. 
`root` is a vector of the root node ids of the graphs, and `leafVal` is a vector of the leaf nodes' weights of the graphs. 

# Arguments:
- `graphs` (AbstractVector{G}): The vector object representing the Feynman graphs,
- `leafMap (Dict{Int,Int})`: The mapping dictionary from the id of each leaf to the index of the leaf weight's table `leafVal`.
- `root` (AbstractVector{Int}, optional): The vector of the root node ids of the graphs (defaults to `[g.id for g in graphs]`).
- `name` (String,optional): The name of the complied function (defaults to `"eval_graph!"`).  
"""
function to_julia_str(graphs::AbstractVector{G}, leafMap::Dict{Int,Int}; root::AbstractVector{Int}=[g.id for g in graphs],
    name::String="eval_graph!") where {G<:AbstractGraph}
    head = "function $name(root::AbstractVector, leafVal::AbstractVector)\n "
    body = ""
    inds_visitedleaf = Int[]
    inds_visitednode = Int[]
    for graph in graphs
        for g in PostOrderDFS(graph) #leaf first search
            if g.id in root
                target = "root[$(findfirst(x -> x == g.id, root))]"
            else
                target = "g$(g.id)"
            end
            if isempty(g.subgraphs) #leaf
                g.id in inds_visitedleaf && continue
                factor_str = g.factor == 1 ? "" : " * $(g.factor)"
                body *= "    $target = leafVal[$(leafMap[g.id])]$factor_str\n "
                push!(inds_visitedleaf, g.id)
            else
                g.id in inds_visitednode && continue
                factor_str = g.factor == 1 ? "" : " * $(g.factor)"
                body *= "    $target = $(_to_static(g.operator, g.subgraphs, g.subgraph_factors))$factor_str\n "
                push!(inds_visitednode, g.id)
            end
        end
    end
    tail = "end"
    return head * body * tail
end

"""
    function compile(graphs::AbstractVector{G}; root::AbstractVector{Int}=[g.id for g in graphs]) where {G<:AbstractGraph}
    
Compile a list of graphs into a julia static function. 
The function takes two arguments: `root` and `leaf`. `root` is a vector of the root node ids of the graphs, and `leaf` is a vector of the leaf node ids of the graphs. 
This function calls to_julia_str and generate a defined function using RuntimeGeneratedFunctions.
Comparing to eval(Meta.parse(to_julia_str(...))), this function does not leak out the function name into global scope.

# Example:
```julia
factor = 1.5
V1 = [𝑓⁺(1)𝑓⁻(2), 𝑓⁺(3)𝑓⁻(4)]
subgraphs = [external_vertex(V1[1]), external_vertex(V1[2])]
g = FeynmanGraph(subgraphs; factor=factor)
# println(g)
eval_graph! = Compilers.compile([g,])
root = [0.0,]
leaf = [1.0, 2.0]

@assert eval_graph!(root, leaf) ≈ (leaf[1] + leaf[2]) * factor
```
"""
function compile(graphs::AbstractVector{G};
    root::AbstractVector{Int}=[g.id for g in graphs]) where {G<:AbstractGraph}
    # this function return a runtime generated function defined by compile()
    func_string = to_julia_str(graphs; root=root, name="func_name!")
    func_expr = Meta.parse(func_string)
    return @RuntimeGeneratedFunction(func_expr)
end

function compile(graphs::AbstractVector{G}, leafMap::Dict{Int,Int};
    root::AbstractVector{Int}=[g.id for g in graphs]) where {G<:AbstractGraph}
    # this function return a runtime generated function defined by compile()
    func_string = to_julia_str(graphs, leafMap; root=root, name="func_name!")
    func_expr = Meta.parse(func_string)
    return @RuntimeGeneratedFunction(func_expr)
end