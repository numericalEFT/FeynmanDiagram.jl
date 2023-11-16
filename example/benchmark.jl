using FeynmanDiagram
using FeynmanDiagram.Taylor
using FeynmanDiagram.ComputationalGraphs:
    eval!, Leaves
using FeynmanDiagram.Utility:
    taylorexpansion!, count_operation


function assign_leaves(g::FeynmanGraph, taylormap)
    leafmap = Dict{Int,Int}()
    leafvec = Vector{Float64}()
    idx = 0
    for leaf in Leaves(g)
        taylor = taylormap[leaf.id]
        for (order, coeff) in taylor.coeffs
            idx += 1
            push!(leafvec, 1.0 / taylor_factorial(order))
            leafmap[coeff.id] = idx
            print("assign $(order) $(coeff.id)  $(taylor_factorial(order)) $(leafvec[idx])\n")
        end
    end
    return leafmap, leafvec
end

#dict_g, fl, bl, leafmap = diagdictGV(:sigma, [(2, 0, 0), (2, 0, 1), (2, 0, 2), (2, 1, 0), (2, 1, 1), (2, 2, 0), (2, 1, 2), (2, 2, 2)], 3)
dict_g, fl, bl, leafmap = diagdictGV(:sigma, [(3, 0, 0), (3, 0, 3), (3, 0, 2), (3, 0, 1)], 3)
g = dict_g[(3, 0, 0)]

set_variables("x y", orders=[1, 3])
propagator_var = ([true, false], [false, true]) # Specify variable dependence of fermi (first element) and bose (second element) particles.
t, taylormap = taylorexpansion!(g[1][1], propagator_var, (fl, bl))

for (order, graph) in dict_g
    if graph[2][1] == g[2][1]
        idx = 1
    else
        idx = 2
    end
    print("$(count_operation(t.coeffs[[order[2],order[3]]]))\n")
    print("$(count_operation(graph[1][idx]))\n")
    print("$(order) $(eval!(graph[1][idx])) $(eval!(t.coeffs[[order[2],order[3]]]))\n")
end

