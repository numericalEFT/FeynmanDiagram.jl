using FeynmanDiagram
import FeynmanDiagram.FrontEnds: NoHartree
import FeynmanDiagram.Compilers
using FeynmanDiagram.ComputationalGraphs:
        eval!
using FeynmanDiagram.ComputationalGraphs.AbstractTrees

function visualize_graph(g::Graph, fname::String)
    Compilers.compile_dot([g], "$fname.dot")
    run(`dot -Tpdf -o $fname.pdf $fname.dot`)
    println("Generated $fname.pdf")
end

para = Parquet.DiagPara(type = Parquet.SigmaDiag, innerLoopNum = 2, hasTau = true, filter=[NoHartree,]);
sigmadf = Parquet.build(para) 

optimize!(sigmadf.diagram)
print(typeof(sigmadf.diagram)) 

renormalization_orders = [0, 3];

leaf_dep_funcs = [pr -> pr isa FrontEnds.BareGreenId, pr -> pr isa FrontEnds.BareInteractionId];
dict_sigma = taylorAD([sigmadf.diagram[2]], renormalization_orders, leaf_dep_funcs)
dict_sigma_nest = taylorAD_nest([sigmadf.diagram[2]], renormalization_orders, leaf_dep_funcs)




for (order, graph) in dict_sigma
    graph_nest = dict_sigma_nest[order]
    print("order $(order)\n")
    print("taylor $(count_operation(graph[1])), $(eval!(graph[1]))\n")
    print("nest $(count_operation(graph_nest[1])), $(eval!(graph_nest[1]))\n")
    # Visualize each graph for this order
    # order_str = join(order, "_")
    # for (i, (g_taylor, g_nest)) in enumerate(zip(graph, graph_nest))
    #     visualize_graph(g_taylor, "sigma_taylor_order$(order_str)_graph$(i)")
    #     visualize_graph(g_nest,   "sigma_nest_order$(order_str)_graph$(i)")
    # end
end
