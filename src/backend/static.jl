function _to_static(::Type{ComputationalGraphs.Sum}, subgraphs::Vector{Graph{F,W}}, subgraph_factors::Vector{F}) where {F,W}
    if length(subgraphs) == 1
        return "(g$(subgraphs[1].id) * $(subgraph_factors[1]))"
    else
        return "(" * join(["g$(g.id) * $gfactor" for (g, gfactor) in zip(subgraphs, subgraph_factors)], " + ") * ")"
    end
end

function _to_static(::Type{ComputationalGraphs.Prod}, subgraphs::Vector{Graph{F,W}}, subgraph_factors::Vector{F}) where {F,W}
    if length(subgraphs) == 1
        return "(g$(subgraphs[1].id))"
    else
        return "(" * join(["g$(g.id)" for g in subgraphs], " * ") * ")"
    end
end

"""
    function to_julia_str(graphs::AbstractVector; root::AbstractVector{Int}=[g.id for g in graphs], name::String="eval_graph!")
    
Compile a list of graphs into a string for a julia static function. The function takes two arguments: `root` and `leaf`. `root` is a vector of the root node ids of the graphs, and `leaf` is a vector of the leaf node ids of the graphs. 
"""
function to_julia_str(graphs::AbstractVector; root::AbstractVector{Int}=[g.id for g in graphs], name::String="eval_graph!")
    head = "function $name(root::AbstractVector, leaf::AbstractVector)\n "
    body = ""
    leafidx = 1
    for graph in graphs
        for g in PostOrderDFS(graph) #leaf first search
            if g.id in root
                target = "root[$(findfirst(x -> x == g.id, root))]"
            else
                target = "g$(g.id)"
            end
            if isempty(g.subgraphs) #leaf
                body *= "    $target = leaf[$leafidx]\n "
                leafidx += 1
            else
                body *= "    $target = $(_to_static(g.operator, g.subgraphs, g.subgraph_factors))*$(g.factor)\n "
            end
        end
    end
    tail = "end"
    return head * body * tail
end

function to_julia_str(graphs::AbstractVector, propagatorMap::Dict{Int,Int}, interactionMap::Dict{Int,Int}; root::AbstractVector{Int}=[g.id for g in graphs], name::String="eval_graph!")
    head = "function $name(root::AbstractVector, propagatorVal::AbstractVector, interactionVal::AbstractVector)\n "
    body = ""
    pIdx, iIdx = 1, 1
    for graph in graphs
        for g in PostOrderDFS(graph) #leaf first search
            if g.id in root
                target = "root[$(findfirst(x -> x == g.id, root))]"
            else
                target = "g$(g.id)"
            end
            if isempty(g.subgraphs) #leaf
                if g.type == ComputationalGraphs.Propagator
                    body *= "    $target = propagatorVal[$(propagatorMap[g.id])]\n "
                    pIdx += 1
                elseif g.type == ComputationalGraphs.Interaction
                    body *= "    $target = interactionVal[$(interactionMap[g.id])]\n "
                    iIdx += 1
                end
            else
                body *= "    $target = $(_to_static(g.operator, g.subgraphs, g.subgraph_factors))*$(g.factor)\n "
            end
        end
    end
    tail = "end"
    return head * body * tail
end

"""
    function compile(graphs::AbstractVector; root::AbstractVector{Int}=[g.id for g in graphs])
    
Compile a list of graphs into a julia static function. 
The function takes two arguments: `root` and `leaf`. `root` is a vector of the root node ids of the graphs, and `leaf` is a vector of the leaf node ids of the graphs. 
This function calls to_julia_str and generate a defined function using RuntimeGeneratedFunctions.
Comparing to eval(Meta.parse(to_julia_str(...))), this function does not leak out the function name into global scope.

# Example:
```julia
factor = 1.5
V1 = [𝑓⁺(1)𝑓⁻(2), 𝑓⁺(3)𝑓⁻(4)]
subgraphs = [external_vertex(V1[1]), external_vertex(V1[2])]
g = Graph(subgraphs; factor=factor)
# println(g)
eval_graph! = Compilers.compile([g,])
root = [0.0,]
leaf = [1.0, 2.0]

@assert eval_graph!(root, leaf) ≈ (leaf[1] + leaf[2]) * factor
```
"""
function compile(graphs::AbstractVector;
    root::AbstractVector{Int}=[g.id for g in graphs])
    # this function return a runtime generated function defined by compile()
    func_string = to_julia_str(graphs; root=root, name="func_name!")
    func_expr = Meta.parse(func_string)
    return @RuntimeGeneratedFunction(func_expr)
end

function compile(graphs::AbstractVector, propagatorMap::Dict{Int,Int}, interactionMap::Dict{Int,Int};
    root::AbstractVector{Int}=[g.id for g in graphs])
    # this function return a runtime generated function defined by compile()
    func_string = to_julia_str(graphs, propagatorMap, interactionMap; root=root, name="func_name!")
    func_expr = Meta.parse(func_string)
    return @RuntimeGeneratedFunction(func_expr)
end