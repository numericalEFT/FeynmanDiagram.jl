"""
    function to_static(operator::Type, subgraphs::AbstractVector{<:AbstractGraph}, subgraph_factors::AbstractVector)

Returns the static representation of a computational graph node `g` with operator `operator`, subgraphs `subgraphs`, and subgraph factors `subgraph_factors`.
"""
function to_static(operator::Type, subgraphs::AbstractVector{<:AbstractGraph}, subgraph_factors::AbstractVector)
    error(
        "Static representation for computational graph nodes with operator $(operator) not yet implemented! " *
        "Please define a method `to_static(::Type{$(operator)}, subgraphs::$(typeof(subgraphs)), subgraph_factors::$(typeof(subgraph_factors)))`."
    )
end

function to_static(::Type{ComputationalGraphs.Sum}, subgraphs::Vector{Graph{F,W}}, subgraph_factors::Vector{F}) where {F,W}
    if length(subgraphs) == 1
        return "(g$(subgraphs[1].id) * $(subgraph_factors[1]))"
    else
        return "(" * join(["g$(g.id) * $gfactor" for (g, gfactor) in zip(subgraphs, subgraph_factors)], " + ") * ")"
    end
end

function to_static(::Type{ComputationalGraphs.Prod}, subgraphs::Vector{Graph{F,W}}, subgraph_factors::Vector{F}) where {F,W}
    if length(subgraphs) == 1
        return "(g$(subgraphs[1].id))"
    else
        return "(" * join(["g$(g.id)" for g in subgraphs], " * ") * ")"
    end
end

function to_static(::Type{ComputationalGraphs.Sum}, subgraphs::Vector{FeynmanGraph{F,W}}, subgraph_factors::Vector{F}) where {F,W}
    if length(subgraphs) == 1
        return "(g$(subgraphs[1].id) * $(subgraph_factors[1]))"
    else
        return "(" * join(["g$(g.id) * $gfactor" for (g, gfactor) in zip(subgraphs, subgraph_factors)], " + ") * ")"
    end
end

function to_static(::Type{ComputationalGraphs.Prod}, subgraphs::Vector{FeynmanGraph{F,W}}, subgraph_factors::Vector{F}) where {F,W}
    if length(subgraphs) == 1
        return "(g$(subgraphs[1].id))"
    else
        return "(" * join(["g$(g.id)" for g in subgraphs], " * ") * ")"
    end
end


"""
    function to_julia_str(graphs::AbstractVector{<:AbstractGraph}; root::AbstractVector{Int}=[id(g) for g in graphs], name::String="eval_graph!")
    
Compile a list of graphs into a string for a julia static function. The function takes two arguments: `root` and `leaf`. 
`root` is a vector of the root node ids of the graphs, and `leaf` is a vector of the leaf nodes' weights of the graphs. 
"""
function to_julia_str(graphs::AbstractVector{<:AbstractGraph}; root::AbstractVector{Int}=[id(g) for g in graphs], name::String="eval_graph!")
    head = "function $name(root::AbstractVector, leaf::AbstractVector)\n "
    body = ""
    leafidx = 1
    for graph in graphs
        for g in PostOrderDFS(graph) #leaf first search
            if id(g) in root
                target = "root[$(findfirst(x -> x == id(g), root))]"
            else
                target = "g$(id(g))"
            end
            if isempty(subgraphs(g)) #leaf
                body *= "    $target = leaf[$leafidx]\n "
                leafidx += 1
            else
                body *= "    $target = $(to_static(operator(g), subgraphs(g), subgraph_factors(g)))*$(factor(g))\n "
            end
        end
    end
    tail = "end"
    return head * body * tail
end

"""
    function to_julia_str(graphs::AbstractVector{<:AbstractGraph}, leafMap::Dict{Int,Int}; root::AbstractVector{Int}=[id(g) for g in graphs],
        name::String="eval_graph!")
    
Compile a list of Feynman graphs into a string for a julia static function. The complied function takes two arguments: `root` and `leafVal`. 
`root` is a vector of the root node ids of the graphs, and `leafVal` is a vector of the leaf nodes' weights of the graphs. 

# Arguments:
- `graphs` (AbstractVector{G}): The vector object representing the Feynman graphs,
- `leafMap (Dict{Int,Int})`: The mapping dictionary from the id of each leaf to the index of the leaf weight's table `leafVal`.
- `root` (AbstractVector{Int}, optional): The vector of the root node ids of the graphs (defaults to `[id(g) for g in graphs]`).
- `name` (String,optional): The name of the complied function (defaults to `"eval_graph!"`).  
"""
function to_julia_str(graphs::AbstractVector{<:AbstractGraph}, leafMap::Dict{Int,Int}; root::AbstractVector{Int}=[id(g) for g in graphs],
    name::String="eval_graph!")
    head = "function $name(root::AbstractVector, leafVal::AbstractVector)\n "
    body = ""
    for graph in graphs
        for g in PostOrderDFS(graph) #leaf first search
            if id(g) in root
                target = "root[$(findfirst(x -> x == id(g), root))]"
            else
                target = "g$(id(g))"
            end
            if isempty(subgraphs(g)) #leaf
                name(g) == "compiled" && continue
                body *= "    $target = leafVal[$(leafMap[id(g)])]\n "
                set_name!(g, "compiled")
            else
                body *= "    $target = $(to_static(operator(g), subgraphs(g), subgraph_factors(g)))*$(factor(g))\n "
            end
        end
    end
    tail = "end"
    return head * body * tail
end

"""
    function compile(graphs::AbstractVector{<:AbstractGraph}; root::AbstractVector{Int}=[id(g) for g in graphs])
    
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
function compile(graphs::AbstractVector{<:AbstractGraph};
    root::AbstractVector{Int}=[id(g) for g in graphs])
    # this function return a runtime generated function defined by compile()
    func_string = to_julia_str(graphs; root=root, name="func_name!")
    func_expr = Meta.parse(func_string)
    return @RuntimeGeneratedFunction(func_expr)
end

function compile(graphs::AbstractVector{<:AbstractGraph}, leafMap::Dict{Int,Int};
    root::AbstractVector{Int}=[id(g) for g in graphs])
    # this function return a runtime generated function defined by compile()
    func_string = to_julia_str(graphs, leafMap; root=root, name="func_name!")
    func_expr = Meta.parse(func_string)
    return @RuntimeGeneratedFunction(func_expr)
end