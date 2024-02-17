"""
    function to_static(operator::Type, subgraphs::AbstractVector{<:AbstractGraph}, subgraph_factors::AbstractVector)

    Returns the static representation of a computational graph node `g` with operator `operator`, subgraphs `subgraphs`, and subgraph factors `subgraph_factors`.
"""
function to_static(operator::Type, subgraphs::AbstractVector{<:AbstractGraph}, subgraph_factors::AbstractVector; lang::Symbol=:julia)
    error(
        "Static representation for computational graph nodes with operator $(operator) not yet implemented! " *
        "Please define a method `to_static(::Type{$(operator)}, subgraphs::$(typeof(subgraphs)), subgraph_factors::$(typeof(subgraph_factors)))`."
    )
end

function to_static(::Type{ComputationalGraphs.Sum}, subgraphs::Vector{Graph{F,W}}, subgraph_factors::Vector{F}; lang::Symbol=:julia) where {F,W}
    if length(subgraphs) == 1
        factor_str = subgraph_factors[1] == 1 ? "" : " * $(subgraph_factors[1])"
        return "(g$(subgraphs[1].id)$factor_str)"
    else
        terms = ["g$(g.id)" * (gfactor == 1 ? "" : " * $gfactor") for (g, gfactor) in zip(subgraphs, subgraph_factors)]
        return "(" * join(terms, " + ") * ")"
    end
end

function to_static(::Type{ComputationalGraphs.Prod}, subgraphs::Vector{Graph{F,W}}, subgraph_factors::Vector{F}; lang::Symbol=:julia) where {F,W}
    if length(subgraphs) == 1
        factor_str = subgraph_factors[1] == 1 ? "" : " * $(subgraph_factors[1])"
        return "(g$(subgraphs[1].id)$factor_str)"
    else
        terms = ["g$(g.id)" * (gfactor == 1 ? "" : " * $gfactor") for (g, gfactor) in zip(subgraphs, subgraph_factors)]
        return "(" * join(terms, " * ") * ")"
        # return "(" * join(["g$(g.id)" for g in subgraphs], " * ") * ")"
    end
end

function to_static(::Type{ComputationalGraphs.Power{N}}, subgraphs::Vector{Graph{F,W}}, subgraph_factors::Vector{F}; lang::Symbol=:julia) where {N,F,W}
    factor_str = subgraph_factors[1] == 1 ? "" : " * $(subgraph_factors[1])"
    if lang == :julia
        op_str = "^"
    elseif lang == :c || lang == :python
        op_str = "**"
    else
        error("Unsupported language")
    end
    return "((g$(subgraphs[1].id))$(op_str)$N$factor_str)"
end

function to_static(::Type{ComputationalGraphs.Sum}, subgraphs::Vector{FeynmanGraph{F,W}}, subgraph_factors::Vector{F}; lang::Symbol=:julia) where {F,W}
    if length(subgraphs) == 1
        factor_str = subgraph_factors[1] == 1 ? "" : " * $(subgraph_factors[1])"
        return "(g$(subgraphs[1].id)$factor_str)"
    else
        terms = ["g$(g.id)" * (gfactor == 1 ? "" : " * $gfactor") for (g, gfactor) in zip(subgraphs, subgraph_factors)]
        return "(" * join(terms, " + ") * ")"
    end
end

function to_static(::Type{ComputationalGraphs.Prod}, subgraphs::Vector{FeynmanGraph{F,W}}, subgraph_factors::Vector{F}; lang::Symbol=:julia) where {F,W}
    if length(subgraphs) == 1
        factor_str = subgraph_factors[1] == 1 ? "" : " * $(subgraph_factors[1])"
        return "(g$(subgraphs[1].id)$factor_str)"
    else
        terms = ["g$(g.id)" * (gfactor == 1 ? "" : " * $gfactor") for (g, gfactor) in zip(subgraphs, subgraph_factors)]
        return "(" * join(terms, " * ") * ")"
    end
end

function to_static(::Type{ComputationalGraphs.Power{N}}, subgraphs::Vector{FeynmanGraph{F,W}}, subgraph_factors::Vector{F}; lang::Symbol=:julia) where {N,F,W}
    factor_str = subgraph_factors[1] == 1 ? "" : " * $(subgraph_factors[1])"
    if lang == :julia
        op_str = "^"
    elseif lang == :c || lang == :python
        op_str = "**"
    else
        error("Unsupported language")
    end
    return "((g$(subgraphs[1].id))$(op_str)$N$factor_str)"
end

"""
    function to_julia_str(graphs::AbstractVector{<:AbstractGraph}, leafMap::Dict{Int,Int}; root::AbstractVector{Int}=[id(g) for g in graphs],
        name::String="eval_graph!")
    
Compile a list of Feynman graphs into a string for a julia static function. The complied function takes two arguments: `root` and `leafVal`. 
`root` is a vector of the root node ids of the graphs, and `leafVal` is a vector of the leaf nodes' weights of the graphs. 

# Arguments:
- `graphs` (AbstractVector{G}): The vector object representing the Feynman graphs,
- `root` (AbstractVector{Int}, optional): The vector of the root node ids of the graphs (defaults to `[id(g) for g in graphs]`).
- `name` (String,optional): The name of the complied function (defaults to `"eval_graph!"`). 

# Returns:
- A String representing the compiled Julia function.
- `leafMap (Dict{Int,G})`: A dictionary that maps the index of the leaf weight's table `leafVal` to the leaf graph.
"""
function to_julia_str(graphs::AbstractVector{<:AbstractGraph}; root::AbstractVector{Int}=[id(g) for g in graphs],
    name::String="eval_graph!")
    head = "\nfunction $name(root::AbstractVector, leafVal::AbstractVector)\n"
    body = ""
    inds_visitedleaf = Int[]
    inds_visitednode = Int[]
    idx_leafVal = 1
    map_validx_leaf = Dict{Int,eltype(graphs)}()  # mapping from the index of the leafVal to the leaf graph 
    for graph in graphs
        for g in PostOrderDFS(graph) #leaf first search
            g_id = id(g)
            target = "g$(g_id)"
            isroot = false
            if g_id in root
                target_root = "root[$(findfirst(x -> x == g_id, root))]"
                isroot = true
            end
            if isempty(subgraphs(g)) #leaf
                g_id in inds_visitedleaf && continue
                body *= "    $target = leafVal[$idx_leafVal]\n"
                map_validx_leaf[idx_leafVal] = g
                idx_leafVal += 1
                push!(inds_visitedleaf, g_id)
            else
                g_id in inds_visitednode && continue
                body *= "    $target = $(to_static(operator(g), subgraphs(g), subgraph_factors(g)))\n"
                push!(inds_visitednode, g_id)
            end
            if isroot
                body *= "    $target_root = $target\n"
            end
        end
    end
    tail = "end"
    return head * body * tail, map_validx_leaf
end

function julia_to_C_typestr(type::DataType)
    if type == Float64
        return "double "
    elseif type == Float32
        return "float "
    elseif type == Int64
        return "long long "
    elseif type == Int32
        return "int "
    elseif type == ComplexF32
        return "complex float "
    elseif type == ComplexF64
        return "complex double "
    elseif type <: Array
        return julia_to_C_typestr(eltype(type)) * "*"
    else
        error("Unsupported type")
    end
end

function to_Cstr(graphs::AbstractVector{<:AbstractGraph}; root::AbstractVector{Int}=[id(g) for g in graphs],
    datatype::DataType=_dtype.weight, name::String="eval_graph")
    # head = "#include <stdio.h>"
    ctype_str = julia_to_C_typestr(datatype)
    head = "\nvoid $name($ctype_str*root, $ctype_str*leafVal)\n{\n"

    declare = "    $ctype_str"
    body = ""
    inds_visitedleaf = Int[]
    inds_visitednode = Int[]
    idx_leafVal = 0
    map_validx_leaf = Dict{Int,eltype(graphs)}()  # mapping from the index of the leafVal to the leaf graph 
    for graph in graphs
        for g in PostOrderDFS(graph) #leaf first search
            g_id = id(g)
            target = "g$(g_id)"
            isroot = false
            if g_id in root
                target_root = "root[$(findfirst(x -> x == g_id, root)-1)]"
                isroot = true
            end
            if isempty(subgraphs(g)) #leaf
                g_id in inds_visitedleaf && continue
                declare *= " g$g_id,"
                body *= "    $target = leafVal[$idx_leafVal];\n"
                idx_leafVal += 1
                map_validx_leaf[idx_leafVal] = g
                push!(inds_visitedleaf, g_id)
            else
                g_id in inds_visitednode && continue
                declare *= " g$g_id,"
                body *= "    $target = $(to_static(operator(g), subgraphs(g), subgraph_factors(g), lang=:c));\n"
                push!(inds_visitednode, g_id)
            end
            if isroot
                body *= "    $target_root = $target;\n"
            end
        end
    end
    declare = chop(declare) * ";\n"
    tail = "}"
    return head * declare * body * tail, map_validx_leaf
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
    func_string, leafmap = to_julia_str(graphs; root=root, name="eval_graph!")
    func_expr = Meta.parse(func_string)
    return @RuntimeGeneratedFunction(func_expr), leafmap
end

"""
    function compile_Julia(graphs::AbstractVector{<:AbstractGraph}, filename::String; 
        root::AbstractVector{Int}=[id(g) for g in graphs], func_name="eval_graph!")

    Compiles a set of graphs into Julia code and append the generated code to a specified file. 

# Arguments
- `graphs::AbstractVector{<:AbstractGraph}`: An array of graph objects. These graphs are processed to generate Julia code.
- `filename::String`: The name of the file to which the generated code will be appended. The file is created if it does not exist.
- `root::AbstractVector{Int}` (keyword): An array of integers representing root nodes for each graph in `graphs`. By default, it is an array of IDs obtained by calling `id(g)` for each graph `g` in `graphs`.
- `func_name::String` (keyword): The base name for the function(s) to be generated. Defaults to `"eval_graph!"`.

# Returns
- A dictionary (`leafmap`) that maps the index of the leaf weight's table `leafVal` to the leaf graph.
"""
function compile_Julia(graphs::AbstractVector{<:AbstractGraph}, filename::String;
    root::AbstractVector{Int}=[id(g) for g in graphs], func_name="eval_graph!")
    func_string, leafmap = to_julia_str(graphs; root=root, name=func_name)
    open(filename, "a") do f
        write(f, func_string)
    end
    return leafmap
end

"""
    function compile_C(graphs::AbstractVector{<:AbstractGraph}, filename::String; 
        datatype::DataType=_dtype.weight, root::AbstractVector{Int}=[id(g) for g in graphs], func_name="eval_graph")

    Compiles a set of graphs into C language code and append the generated code to a specified file. 

# Arguments
- `datatype::DataType`: This type is used for variables types in the generated C code.
- `graphs::AbstractVector{<:AbstractGraph}`: An array of graph objects. These graphs are processed to generate Julia code.
- `filename::String`: The name of the file to which the generated code will be appended. The file is created if it does not exist.
- `root::AbstractVector{Int}` (keyword): An array of integers representing root nodes for each graph in `graphs`. By default, it is an array of IDs obtained by calling `id(g)` for each graph `g` in `graphs`.
- `func_name::String` (keyword): The base name for the function(s) to be generated. Defaults to `"eval_graph"`.

# Returns
- A dictionary (`leafmap`) that maps the index of the leaf weight's table `leafVal` to the leaf graph.
"""
function compile_C(graphs::AbstractVector{<:AbstractGraph}, filename::String;
    datatype::DataType=_dtype.weight, root::AbstractVector{Int}=[id(g) for g in graphs], func_name="eval_graph")
    func_string, leafmap = to_Cstr(graphs; datatype=datatype, root=root, name=func_name)
    open(filename, "a") do f
        write(f, func_string)
    end
    return leafmap
end