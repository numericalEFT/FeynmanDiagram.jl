
# this is a standalone file implementing:
# Basic IR transform required by the Parquet Algorithm

using Test

module IRTransform

using FeynmanDiagram
using FeynmanDiagram.ComputationalGraphs
using FeynmanDiagram.QuantumOperators
using FeynmanDiagram.QuantumOperators: QuantumOperator

# relabel constructor for QuantumOperator
QuantumOperator(qo::QuantumOperator, label::Int) = QuantumOperator(qo.operator(), label, qo.is_ghost)
# relabel constructor for OperatorProduct

export relabel!, standardize_labels!
export relabel, standardize_labels

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

    for i in 1:length(g.topology)
        for j in 1:length(g.topology[i])
            if haskey(map, g.topology[i][j])
                g.topology[i][j] = map[g.topology[i][j]]
            end
        end
        g.topology[i] = unique(g.topology[i])
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

end



@testset "IRTransform" begin
    using .IRTransform
    using .IRTransform.FeynmanDiagram
    using .IRTransform.FeynmanDiagram.ComputationalGraphs

    @testset "relabel" begin
        # construct a graph
        V = [𝑓ₑ(1), 𝑓⁺(2)𝑓⁻(3)𝑏⁺(4), 𝜙(5)𝑓⁺(6)𝑓⁻(7), 𝑓(8)𝑏⁻(9)𝜙(10)]
        g1 = feynman_diagram(V, [[2, 3, 4, 9], [5, 6, 7, 10], [8, 1]], external=[8])

        map = Dict(4 => 1, 6 => 1, 8 => 1, 9 => 1, 10 => 1)
        g2 = relabel(g1, map)
        uniqlabels = IRTransform.collect_labels(g2)
        @test uniqlabels == [1, 2, 3, 5, 7]

        map = Dict([i => 1 for i in 2:10])
        g3 = relabel(g1, map)
        uniqlabels = IRTransform.collect_labels(g3)
        @test uniqlabels == [1,]
    end

    @testset "standardize_labels" begin
        V = [𝑓ₑ(1), 𝑓⁺(2)𝑓⁻(3)𝑏⁺(4), 𝜙(5)𝑓⁺(6)𝑓⁻(7), 𝑓(8)𝑏⁻(9)𝜙(10)]
        g1 = feynman_diagram(V, [[2, 3, 4, 9], [5, 6, 7, 10], [8, 1]], external=[8])

        map = Dict([i => (11 - i) for i in 1:5])
        g2 = relabel(g1, map)

        g3 = standardize_labels(g2)
        uniqlabels = IRTransform.collect_labels(g3)
        @test uniqlabels == [1, 2, 3, 4, 5]
    end


end
