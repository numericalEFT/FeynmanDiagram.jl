using FeynmanDiagram
using FeynmanDiagram.ComputationalGraphs:
    eval!, forwardAD, node_derivative, backAD, build_all_leaf_derivative, Leaves, taylorexpansion!, count_operation, build_derivative_backAD!
g1 = Graph([])
g2 = Graph([])
g3 = Graph([], factor=2.0)
G3 = g1
G4 = 1.0 * g1 * g2
G5 = 1.0 * (3.0 * G3 + 0.5 * G4)
G6 = (1.0 * g1 + 0.5 * g2) * (1.5 * g1 + g3) * g1
using FeynmanDiagram.Taylor:
    TaylorSeries, getcoeff, findidx, set_variables

set_variables("x y", order=4)
T5 = taylorexpansion!(G5)
print("$(count_operation(T5.coeffs))\n")
for i in eachindex(T5.coeffs)
    # gs = Compilers.to_julia_str([T5.coeffs[i],], name="eval_graph!")
    # println("T5.coeffs[$(i)]  ", gs, "\n")
    print("$(T5.order[i]) $(eval!(T5.coeffs[i])) $(count_operation(T5.coeffs[i]))\n")
end

T5_compare = build_derivative_backAD!(G5)
print("$(count_operation(T5_compare.coeffs))\n")
for i in eachindex(T5_compare.coeffs)
    # gs = Compilers.to_julia_str([T5_compare.coeffs[i],], name="eval_graph!")
    # println("T5_compare.coeffs[$(i)]  ", gs, "\n")
    print("$(T5_compare.order[i]) $(eval!(T5_compare.coeffs[i])) $(count_operation(T5_compare.coeffs[i]))\n")
end

# derive = build_all_leaf_derivative(G6, 3)
# T6 = TaylorSeries(typeof(G6), typeof(G6), "", derive, var)
# for (order_vec, graph) in derive
#     print("$(order_vec), $(eval!(graph)) \n")
# end

# derive = build_all_leaf_derivative(G5, 3)
# #var = Set(leaf for leaf in Leaves(G5))
# var = Dict(leaf.id => leaf for leaf in Leaves(G5))
# T5 = TaylorSeries(typeof(G5), typeof(G5), "", derive, var)
# for (order_vec, graph) in derive
#     print("$(order_vec), $(eval!(graph)) \n")
# end

# T = T5 * T6
# for (order_vec, graph) in T.expansion
#     print("T5*T6 $(order_vec), $(eval!(graph)) \n")
# end

# G = G5 * G6
# derive = build_all_leaf_derivative(G, 5)
# #var = Set(leaf for leaf in Leaves(G))
# var = Dict(leaf.id => leaf for leaf in Leaves(G))
# T_compare = TaylorSeries(typeof(G), typeof(G), "", derive, var)

# for (order_vec, graph) in T_compare.expansion
#     value = eval!(graph)
#     @assert value == eval!(T.expansion[order_vec])
#     print("compare  T5*T6 $(order_vec), $(value) \n")
# end


