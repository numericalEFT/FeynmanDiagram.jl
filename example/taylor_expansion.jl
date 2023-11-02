using FeynmanDiagram
using FeynmanDiagram.Taylor
using FeynmanDiagram.ComputationalGraphs:
    eval!, forwardAD, node_derivative, backAD, build_all_leaf_derivative
using FeynmanDiagram.Utility:
    taylorexpansion!, taylorexpansion, build_derivative_backAD!, count_operation

function banchmark_AD(glist::Vector{T}) where {T<:Graph}
    taylormap = Dict{Int,TaylorSeries{T}}()
    totaloperation = [0, 0]
    taylorlist = Vector{TaylorSeries{T}}()
    for g in glist
        @time t, taylormap = taylorexpansion!(g; taylormap=taylormap)


        operation = count_operation(t)
        totaloperation = totaloperation + operation
        push!(taylorlist, t)
        print("operation number: $(operation)\n")
        t_compare = build_derivative_backAD!(g)
        for (order, coeff) in (t_compare.coeffs)
            @assert (eval!(coeff)) == (eval!(Taylor.taylor_factorial(order) * t.coeffs[order]))
            # gs = Compilers.to_julia_str([coeff,], name="eval_graph!")
            # println("$(order)  ", gs, "\n")
            # print("$(order) $(eval!(coeff)) $(eval!(getderivative(t5,order))) $(count_operation(coeff))\n")
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
#G5 = g5 * g6
#G6 = (0.5 * g1 * g2 + 0.5 * g2 * g3)
#G6 = (g1 + g2) * (0.5 * g1 + g3) * g1 #   (0.5 * g1 + g3)
#G6 = g1 * g2 * g3 * g4 * g5 * g6
G6 = (1.0 * g1 + 2.0 * g2) * (g1 + g3)
#G6 = 1.5 * g1*g1 + 0.5 * g2 * 1.5 * g1 +    0.5*g2*g3
using FeynmanDiagram.Taylor:
    TaylorSeries, getcoeff, set_variables


# set_variables("x y", order=3)
# @time T5 = taylorexpansion(G5)
# print(T5)
set_variables("x y", orders=[3, 2])
#set_variables("x y z a", order=[1, 2, 3, 2])


#@time taylormap = taylorexpansion(G6)


banchmark_AD([G3, G4, G5, G6])


# T5 = taylormap[G6.id]
# #order = [0, 0, 0, 0, 0, 0]
# #@time print(T5.coeffs[order])
# print("$(count_operation(T5.coeffs))\n")
# # for (order, coeff) in (T5.coeffs)
# #     #gs = Compilers.to_julia_str([coeff,], name="eval_graph!")
# #     #println("$(order)  ", gs, "\n")
# #     print("$(order) $(eval!(coeff)) $(eval!(getcoeff(T5,order))) $(coeff.id) $(count_operation(coeff))\n")
# # end

# print("TaylorSeries $(T5)\n")

# @time T5_compare = build_derivative_backAD!(G6)
# print("$(count_operation(T5_compare.coeffs))\n")
# for (order, coeff) in (T5_compare.coeffs)
#     @assert (eval!(coeff)) == (eval!(Taylor.taylor_factorial(order) * T5.coeffs[order]))
#     # gs = Compilers.to_julia_str([coeff,], name="eval_graph!")
#     # println("$(order)  ", gs, "\n")
#     # print("$(order) $(eval!(coeff)) $(eval!(getderivative(T5,order))) $(count_operation(coeff))\n")
# end





