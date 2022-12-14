# this file is included in ComputationalGraphs.jl

# relabel constructor for QuantumOperator
QuantumOperator(qo::QuantumOperator, label::Int) = QuantumOperator(qo.operator(), label, qo.is_ghost)
# relabel constructor for OperatorProduct

function relabel!(g::Graph, map::Dict{Int,Int})

    for i in 1:length(g.external)
        if haskey(map, g.external[i])
            g.external[i] = map[g.external[i]]
        end
    end

    for i in 1:length(g.vertices)
        op = g.vertices[i]
        for j in 1:length(op.operators)
            qo = op.operators[j]
            if haskey(map, qo.label)
                op.operators[j] = QuantumOperator(qo, map[qo.label])
            end
        end
    end

    for i in 1:length(g.subgraphs)
        relabel!(g.subgraphs[i], map)
    end

    return g
end

relabel(g::Graph, map::Dict{Int,Int}) = relabel!(deepcopy(g), map)

function collect_labels(g::Graph)
    labels = Vector{Int}([])
    for i in 1:length(g.vertices)
        op = g.vertices[i]
        for j in 1:length(op.operators)
            qo = op.operators[j]
            if !(qo.label in labels)
                push!(labels, qo.label)
            end
        end
    end

    uniqlables = sort(unique(labels))
end

function standardize_labels!(g::Graph)
    #TBD
    uniqlabels = collect_labels(g)
    map = Dict{Int,Int}()
    for i in 1:length(uniqlabels)
        push!(map, uniqlabels[i] => i)
    end
    return relabel!(g, map)
end

standardize_labels(g::Graph) = standardize_labels!(deepcopy(g))

############LEGACY BELOW################

# function readDiag(io::IO)

#     return Diagram
# end

# function printDiag(diag::Union{Tuple,AbstractVector}, kwargs...)

#     f = open("./DiagFiles/Diag.txt", "w")
#     writedlm(f, [])
# end
