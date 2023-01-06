module Deprecated

# """Type alias for a directed graph edge e = (a₁⁺, a₂⁻) from e[1] to e[2]."""
const EdgeType = Tuple{QuantumOperator,QuantumOperator}

function Base.:+(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
    @assert g1.type == g2.type "g1 and g2 are not of the same type."
    # TODO: more check
    @assert Set(vertices(g1)) == Set(vertices(g2)) "g1 and g2 have different vertices."
    @assert Set(external(g1)) == Set(external(g2)) "g1 and g2 have different external vertices."
    @assert g1.orders == g2.orders "g1 and g2 have different orders."

    return Graph(g1.vertices; external=g1.external, type=g1.type, subgraphs=[g1, g2], operator=Sum(), ftype=F, wtype=W)
end

function Base.:-(g1::Graph{F,W}, g2::Graph{F,W}) where {F,W}
    return g1 + (-F(1)) * g2
end

"""
    function internal_vertices(g::Graph)

    Return all internal vertices (::Vector{OperatorProduct}) of Graph `g`.
"""
internal_vertices(g::Graph) = g.vertices[setdiff(eachindex(g.vertices), g.external)]

"""
    function feynman_diagram(vertices::Vector{OperatorProduct}, contractions::Vector{Int};
        external=[], factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", type=:generic)

    Create a Graph representing feynman diagram from all vertices and Wick contractions.

# Arguments:
- `vertices::Vector{OperatorProduct}`  vertices of the diagram
- `contractions::Vector{Int}`  contraction-index vector respresnting Wick contractions
- `external`  index of external vertices
- `factor`  scalar multiplicative factor for the diagram
- `weight`  weight of the diagram
- `name`  name of the diagram
- `type`  type of the diagram

# Example:
```julia-repl
julia> g = feynman_diagram([𝑓⁺(1)𝑓⁻(2)𝜙(3), 𝑓⁺(4)𝑓⁻(5)𝜙(6)], [1, 2, 3, 2, 1, 3])
1: generic graph from f⁺(1)f⁻(2)ϕ(3)|f⁺(4)f⁻(5)ϕ(6)

julia> g.subgraphs
3-element Vector{Graph{Float64, Float64}}:
 2: propagtor graph from f⁺(1)f⁻(5)
 3: propagtor graph from f⁻(2)f⁺(4)
 4: propagtor graph from ϕ(3)ϕ(6)
```
"""
function feynman_diagram(vertices::Vector{OperatorProduct}, contractions::Vector{Int};
    external=[], factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", type=:generic)

    contraction_sign, topology, edges = contractions_to_edges(vertices, contractions)
    g = Graph(vertices; external=external, topology=topology, name=name, type=type, operator=Prod(),
        factor=factor * contraction_sign, weight=weight)
    for edge in edges
        push!(g.subgraphs, propagator(reduce(*, edge)))
    end
    return g
end
function feynman_diagram(graphs::Vector{Graph{F,W}}, contractions::Vector{Int};
    external=[], factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", type=:generic) where {F,W}

    vertices = [v for g in graphs for v in external_vertices(g)]
    return feynman_diagram(vertices, contractions; external=external, factor=factor, weight=weight, name=name, type=type)
end

function feynman_diagram(graphs::Vector{Graph{F,W}}, topology::Vector{Vector{Int}};
    external=[], factor=one(_dtype.factor), weight=zero(_dtype.weight), name="", type=:generic) where {F,W}

    vertices = reduce(*, [v for g in graphs for v in external_vertices(g)])
    return feynman_diagram(vertices, topology; external=external, factor=factor, weight=weight, name=name, type=type)
end

"""
    function contractions_to_edges(vertices::Vector{OperatorProduct}, contractions::Vector{Int})

    Converts a list of Wick contractions associated with a list of vertices
    to a list of edges e = (i, j) directed from `operators[i]` to `operators[j]`, where
    `operators` is a flattened list of operators associated with the specified `vertices`.

# Example: 
```julia-repl
julia> vertices = [𝑏⁺(1)𝑓⁺(2)𝜙(3), 𝑓⁻(4)𝑓⁻(5), 𝑏⁻(6)𝑓⁺(7)𝜙(8)];

julia> edges, sign = contractions_to_edges(vertices, [1, 2, 3, 2, 4, 1, 4, 3])
(Tuple{QuantumOperator, QuantumOperator}[(b⁺(1), b⁻(6)), (f⁺(2), f⁻(4)), (ϕ(3), ϕ(8)), (f⁻(5), f⁺(7))], 1)
```
- flattened fermionic edges: [2, 4, 5, 7]
- sign = parity([1, 2, 3, 4]) = 1
"""
function contractions_to_edges(vertices::Vector{OperatorProduct}, contractions::Vector{Int})
    #TODO: only works for weak-coupling expansion with Wick's theorem for now.
    # Obtain the flattened list of non-composite operators
    operators = [o for v in vertices for o in v.operators]
    # Filter some illegal contractions
    is_invalid =
        isodd(length(contractions)) ||
        # isodd(length(unique(contractions))) ||
        length(operators) != length(contractions)
    if is_invalid
        throw(
            ArgumentError(
                "Input $contractions does not specify a legal set of Wick contractions",
            ),
        )
    end
    # Loop over operators and pair
    next_pairing = 1
    edges = Vector{EdgeType}()
    topology = Vector{Int}[]
    permutation = Int[]

    for (i, wick_i) in enumerate(contractions)
        if i < next_pairing
            continue  # Operator already paired
        end
        for (js, wick_j) in enumerate(contractions[(i+1):end])
            j = i + js  # Iterating for (j > i)
            if j < next_pairing
                continue  # Operator already paired
            end
            if wick_i == wick_j
                @debug "Found Wick contraction #$wick_i, adding edge between operators $i and $j"
                @assert operators[j]'.operator == operators[i].operator
                isfermionic(operators[j]) && append!(permutation, [i, j])
                push!(edges, (operators[i], operators[j]))
                push!(topology, [i, j])
                # Move on to next pair
                next_pairing += 1
                break
            end
        end
    end
    # Deduce the parity of the contraction
    # Get the permutation parity for the flattened list of fermionic edges, e.g.,
    # permutation = [1, 5, 4, 8, 7, 6] => P = (1 3 2 6 5 4) => sign = +1
    sign = isempty(permutation) ? 1 : parity(sortperm(permutation))

    return sign, topology, edges
end

# function addSubDiagram!(parent::Diagram, child::Diagram)
#     for c in parent.subdiagram
#         if c.id == child.id
#             return false
#         end
#     end
#     push!(parent.subdiagram, deepcopy(child))
# end

# function addSubDiagram!(parent::Diagram, child::Vector{Diagram{W}}) where {W}
#     for d in child
#         addSubDiagram!(parent, d)
#     end
# end

# _diagram(df, index) = df[index, :Diagram]

function _combinegroups(groups, factor, operator, name)
    # combine diagrams in a group into one composite diagram
    gdf = combine(groups) do group # for each group in groups
        # check the documentation of ``combine" for details https://dataframes.juliadata.org/stable/man/split_apply_combine/
        # id = isnothing(getid) ? GenericId(group.diagram[1].id.para, Tuple(group[1, fields])) : getid(group)

        if nrow(group) == 1
            # if there is only one diagram in df, and the new id is either GenericId or the id of the existing diagram, 
            # then simply return the current df without creating a new diagram
            # ! the new factor will be multiplied to the factor of the exisiting diagram!
            if group isa Diagram
                # diag = deepcopy(group[1, :diagram])
                diag = group[1]
                diag.factor *= factor
                return (diagram=diag, hash=diag.hash)
            end
        end
        W = typeof(group[1].para), typeof(group[1].weight)
        # generate new para, legs and type.
        diag = Diagram{W}(newpara, newlegs, groups; operator=operator, name=name, factor=factor)
        return (diagram=diag, hash=diag.hash)
    end
    return gdf
end

function _mergediag(::Type{W}, group, factor, operator, name) where {W}
    if nrow(group) == 1
        # if there is only one diagram in df, and the new id is either GenericId or the id of the existing diagram, 
        # then simply return the current df without creating a new diagram
        # ! the new factor will be multiplied to the factor of the exisiting diagram!
        # if id isa GenericId || typeof(id) == typeof(group.diagram[1].id)
        # diag = deepcopy(group[1, :diagram])
        diag = group[1]
        diag.factor *= factor
        return diag
        # end
    end
    # return Diagram(group[1].para, group[1].legs, group.subdiagram; operator=operator, name=name, factor=factor)
end

function _combine(::Type{W}, groups, factor, getid, operator, name) where {W}
    """
    # if fields = [:response, :extT], then

    # 1. groups.cols is like: Vector{Symbol}[:response, :extT]

    # 2. groups.keymap is like: 

    #     Dict{Any, Int64} with 2 entries:
    #     (UpDown, (1, 1, 1, 1)) => 2
    #     (UpUp, (1, 1, 1, 1))   => 1
    # """
    d = Dict{Symbol,Any}()
    _keys = keys(groups)
    for col in groupcols(groups)
        d[col] = [key[col] for key in _keys]
    end
    d[:diagram] = [_mergediag(W, groups[key], factor, operator, name) for key in _keys]
    d[:hash] = [diag.hash for diag in d[:diagram]]
    return DataFrame(d, copycols=false)
end

function mergeby(df::DataFrame, fields=Vector{Symbol}();
    operator=Sum(), name::Symbol=:none, factor=1.0,
    getid::Function=g -> GenericId(g[1, :diagram].id.para, Tuple(g[1, fields]))
)
    if isempty(df)
        return df
    else
        W = typeof(df.diagram[1].weight)
        return mergeby(W, df, fields; operator=operator, name=name, factor=factor, getid=getid)
    end
end

function mergeby(::Type{W}, df::DataFrame, fields=Vector{Symbol}();
    operator=Sum(), name::Symbol=:none, factor=1.0,
    getid::Function=g -> GenericId(g[1, :diagram].id.para, Tuple(g[1, fields]))
) where {W}
    if isempty(df)
        return df
    else
        if all(x -> typeof(x.id) == typeof(df.diagram[1].id), df[!, :diagram]) == false
            @warn "Not all DiagramIds in $df are the same!"
        end
        groups = DataFrames.groupby(df, fields, sort=true)
        ########  less memory usage but can not pass the test right now ##############
        d = _combine(W, groups, factor, getid, operator, name)
        ######## alternative approach (more memory)  ##################
        # d = _combinegroups(groups, getid, factor, operator, name)
        # println("old\n$d \n new\n$cd")
        return d
    end
end

function mergeby(diags::Union{Diagram,Tuple,AbstractVector}, fields=nothing; idkey=nothing, kwargs...)
    if diags isa Diagram
        return diags
    else
        if isempty(diags)
            return diags
        else
            W = typeof(diags[1].weight)
            @assert all(x -> (x.weight isa W), diags) "all diagrams should be of the same type. \n$diags"
            diags = collect(diags)
            if isnothing(fields) && isnothing(idkey)
                return mergeby(diags; kwargs...)
            else
                return mergeby(diags, fields; idkey=idkey, kwargs...)
            end
        end
    end
end

# function mergeby(diags::AbstractVector, fields=[]; idkey::Vector{Symbol}=[], kwargs...)
function mergeby(diags::Vector{Diagram{W}}, fields; idkey=Vector{Symbol}(), kwargs...) where {W}
    if isempty(diags)
        return diags
    else
        df = toDataFrame(diags, idkey)
        mergedf = mergeby(df, fields; kwargs...)
        return Vector{Diagram{W}}(mergedf.diagram)
    end
end

function mergeby(diags::Vector{Diagram{W}};
    operator=Sum(), name::Symbol=:none, factor=1.0,
    getid::Function=d -> GenericId(d[1].id.para::DiagPara{W})
) where {W}
    if isempty(diags)
        return diags
    else
        id = getid(diags)
        if length(diags) == 1 && (id isa GenericId || typeof(id) == typeof(diags[1].id))
            # if there is only one diagram, and the new id is either GenericId or the id of the existing diagram, 
            # then simply return the current diagram without creating a new diagram
            # ! the new factor will be multiplied to the factor of the exisiting diagram!
            diags[1].factor *= factor
            return diags
        end
        diag = Diagram{W}(id, operator, diags, name=name, factor=factor)
        return [diag,]
    end
end
# mergeby(df::DataFrame; kwargs...) = mergeby(df, []; kwargs...)
# mergeby(diags::Vector{Diagram{W}}; kwargs...) where {W} = mergeby(diags, []; kwargs...)

end # module Deprecated