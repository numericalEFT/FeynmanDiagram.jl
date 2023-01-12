function _to_static(::Type{ComputationalGraphs.Sum}, subgraphs::Vector{Graph{F,W}}) where {F,W}
    if length(subgraphs) == 1
        return "(g$(subgraphs[1].id))"
    else
        return "(" * join(["g$(g.id)" for g in subgraphs], " + ") * ")"
    end
end

function _to_static(::Type{ComputationalGraphs.Prod}, subgraphs::Vector{Graph{F,W}}) where {F,W}
    if length(subgraphs) == 1
        return "(g$(subgraphs[1].id))"
    else
        return "(" * join(["g$(g.id)" for g in subgraphs], " * ") * ")"
    end
end

"""
    function static_graph(graphs::AbstractVector; root::AbstractVector{Int}=[g.id for g in graphs], name::String="eval_graph!")
    
Compile a list of graphs into a string for a julia static function. The function takes two arguments: `root` and `leaf`. `root` is a vector of the root node ids of the graphs, and `leaf` is a vector of the leaf node ids of the graphs. 

# Example:
```julia-repl
julia> g = Graph([𝑓⁺(1)𝑓⁻(2), 𝑓⁺(3)𝑓⁻(4)], external=[1, 2], subgraphs=[Graph([𝑓⁺(1)𝑓⁻(4)]), Graph([𝑓⁻(2)𝑓⁺(3)])])
3:f⁺(1)f⁻(2)|f⁺(3)f⁻(4)=0.0=⨁ (1,2)

julia> gs = Compilers.static_graph([g, ])
"function eval_graph!(root::AbstractVector, leaf::AbstractVector)\n     g1 = leaf[1]\n     g2 = leaf[2]\n     root[1] = (g1 + g2)*1.0\n end"

julia> eval(Meta.parse(gs)) #compile the string into a callable function `eval_graph!(root, leaf)`
eval_graph! (generic function with 1 method)

julia> leaf = [1.0, 2.0]; root = [0.0,];

julia> eval_graph!(root, leaf)
3.0
"""
function static_graph(graphs::AbstractVector; root::AbstractVector{Int}=[g.id for g in graphs], name::String="eval_graph!")
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
                body *= "    $target = $(_to_static(g.operator, g.subgraphs))*$(g.factor)\n "
            end
        end
    end
    tail = "end"
    return head * body * tail
end
"""
    function static_graph_rgf(graphs::AbstractVector; root::AbstractVector{Int}=[g.id for g in graphs])
    
Compile a list of graphs into a string for a julia static function. The function takes two arguments: `root` and `leaf`. `root` is a vector of the root node ids of the graphs, and `leaf` is a vector of the leaf node ids of the graphs. 
This function calls static_graph and generate a defined function using RuntimeGeneratedFunctions.
Comparing to eval(Meta.parse(static_graph(...))), this function does not leak out the function name into global scope.
"""
function static_graph_rgf(graphs::AbstractVector;
    root::AbstractVector{Int}=[g.id for g in graphs])
    # this function return a runtime generated function defined by static_graph()
    func_string = static_graph(graphs; root=root, name="func_name!")
    func_expr = Meta.parse(func_string)
    return @RuntimeGeneratedFunction(func_expr)
end