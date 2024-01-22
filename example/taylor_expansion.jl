using FeynmanDiagram
using FeynmanDiagram.Taylor
using FeynmanDiagram.ComputationalGraphs:
    eval!, forwardAD, node_derivative, backAD, build_all_leaf_derivative
using FeynmanDiagram.Utility:
    taylorexpansion!, build_derivative_backAD!, count_operation

function benchmark_AD(glist::Vector{T}) where {T<:Graph}
    #taylormap = Dict{Int,TaylorSeries{T}}()
    totaloperation = [0, 0]
    taylorlist = Vector{TaylorSeries{T}}()
    for g in glist
        var_dependence = Dict{Int,Vector{Bool}}()
        for leaf in FeynmanDiagram.Leaves(g)
            var_dependence[leaf.id] = [true for _ in 1:get_numvars()]
        end
        @time t, taylormap, from_coeff_map = taylorexpansion!(g, var_dependence)


        operation = count_operation(t)
        totaloperation = totaloperation + operation
        push!(taylorlist, t)
        print("operation number: $(operation)\n")
        t_compare, leaftaylor = build_derivative_backAD!(g)
        for (order, coeff) in (t_compare.coeffs)
            @assert (eval!(coeff)) == (eval!(Taylor.taylor_factorial(order) * t.coeffs[order]))
        end
    end

    total_uniqueoperation = count_operation(taylorlist)
    print(" total operation number: $(length(taylorlist)) $(totaloperation) $(total_uniqueoperation)\n")
    return total_uniqueoperation
end
g1 = Graph([])
g2 = Graph([])
g3 = Graph([]) #, factor=2.0)
g4 = Graph([])
g5 = Graph([])
g6 = Graph([])
G3 = g1
G4 = 1.0 * g1 * g2
G5 = 1.0 * (3.0 * G3 + 0.5 * G4)
G6 = (1.0 * g1 + 2.0 * g2) * (g1 + g3)

using FeynmanDiagram.Taylor:
    TaylorSeries, getcoeff, set_variables

set_variables("x y", orders=[3, 2])

benchmark_AD([G3, G4, G5, G6])






