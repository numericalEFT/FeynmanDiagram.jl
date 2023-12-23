"""
    function Graph!(d::DiagTree.Diagram{W}; map=Dict{Int,DiagTree.DiagramId}()) where {W}
    
Converts a DiagTree `d` into a Graph, storing the diagram information (which is a DiagramId object) in a Graph.id to DiagramId dictionary ``map".

# Arguments:
- `d`:  DiagTree.Diagram object.
- `map`:  map between the Graph.id and DiagramId. It will be updated with the new nodes and leaves contained in the DiagTree.Diagram `d`.


# Example:
```julia-repl
julia> para = DiagParaF64(type = Ver4Diag, innerLoopNum=2);

julia> ver4=Parquet.build(para)
2×5 DataFrame
 Row │ diagram                            extT          hash   response  type      
     │ Diagram…                           Tuple…        Int64  Response  Analytic… 
─────┼─────────────────────────────────────────────────────────────────────────────
   1 │ 5978:↑↑Dyn#[0, 0, 0, 0],t(1, 1, …  (1, 1, 1, 1)   5978  UpUp      Dynamic
   2 │ 5979:↑↓Dyn#[0, 0, 0, 0],t(1, 1, …  (1, 1, 1, 1)   5979  UpDown    Dynamic

julia> d = ver4.diagram[1]  # take the first diagram
5978:↑↑Dyn#[0, 0, 0, 0],t(1, 1, 1, 1)=0.0=⨁ (5026, 5071, 5137, 5146, 5175, 5220, 5312, 5321, 5350, 5396, 5463, 5473, 5503, 5549, 5642, 5652, 5682, 5793, 5831, 5968)

julia> root = FrontEnds.Graph!(d)
```

"""
# function Graph!(d::DiagTree.Diagram{W}; map=Dict{Int,DiagTree.DiagramId}()) where {W}
function Graph!(d::DiagTree.Diagram{W}) where {W}

    function op(o)
        if o isa DiagTree.Sum
            return ComputationalGraphs.Sum()
        elseif o isa DiagTree.Prod
            return ComputationalGraphs.Prod()
        else
            error("Unknown operator: $o")
        end
    end

    subgraphs = ComputationalGraphs.Graph{W,W}[]
    for g in d.subdiagram
        # res, map = Graph!(g; map=map)
        res = Graph!(g)
        push!(subgraphs, res)
    end

    if isempty(subgraphs)
        root = ComputationalGraphs.Graph(subgraphs; subgraph_factors=ones(W, length(subgraphs)), factor=d.factor, name=String(d.name),
            operator=op(d.operator), orders=d.id.order, ftype=W, wtype=W, weight=d.weight, properties=d.id)
    else
        tree = ComputationalGraphs.Graph(subgraphs; subgraph_factors=ones(W, length(subgraphs)),
            operator=op(d.operator), orders=d.id.order, ftype=W, wtype=W, weight=d.weight)
        root = ComputationalGraphs.Graph([tree,]; subgraph_factors=[d.factor,], orders=tree.orders,
            ftype=W, wtype=W, weight=d.weight * d.factor)
    end

    return root
    # @assert haskey(map, root.id) == false "DiagramId already exists in map: $(root.id)"
    # @assert haskey(map, tree.id) == false "DiagramId already exists in map: $(tree.id)"
    # map[root.id] = d.id
    # map[tree.id] = d.id

    # return root, map
end

"""
    function extract_var_dependence(map::Dict{Int,DiagTree.DiagramId}, ::Type{ID}, numvars::Int)
    
    Given a map between graph id and DiagramId, extract the variable dependence of all graphs.

# Arguments:
- `map::Dict{Int,DiagTree.DiagramId}`:  A dictionary mapping graph ids to DiagramIds. DiagramId stores the diagram information of the corresponding graph. 
- `ID`:  The particular type of ID that has the given variable dependence. 
- `numvars`: The number of variables which the diagram depends on.
"""
function extract_var_dependence(map::Dict{Int,DiagTree.DiagramId}, ::Type{ID}; numvars::Int=1) where {ID<:PropagatorId}
    var_dependence = Dict{Int,Vector{Bool}}()
    for (id, diagID) in map
        if diagID isa ID
            var_dependence[id] = [true for _ in 1:numvars]
        else
            var_dependence[id] = [false for _ in 1:numvars]
        end
    end
    return var_dependence
end