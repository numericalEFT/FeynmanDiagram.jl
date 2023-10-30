using FeynmanDiagram: Taylor as Taylor

@testset verbose = true "TaylorSeries" begin
    using FeynmanDiagram.Taylor:
        getcoeff, set_variables, taylor_factorial
    a, b, c, d, e = set_variables("a b c d e", order=3)
    F1 = (a + b) * (a + b) * (a + b)
    print("$(F1)")
    @test getcoeff(F1, [2, 1, 0, 0, 0]) == 3.0
    @test getcoeff(F1, [1, 2, 0, 0, 0]) == 3.0
    @test getcoeff(F1, [3, 0, 0, 0, 0]) == 1.0
    @test getcoeff(F1, [0, 3, 0, 0, 0]) == 1.0
    F2 = (1 + a) * (3 + 2c)
    @test getcoeff(F2, [0, 0, 0, 0, 0]) == 3.0
    @test getcoeff(F2, [1, 0, 0, 0, 0]) == 3.0
    @test getcoeff(F2, [0, 0, 1, 0, 0]) == 2.0
    @test getcoeff(F2, [1, 0, 1, 0, 0]) == 2.0
    F3 = (a + b)^3
    @test getcoeff(F1, [2, 1, 0, 0, 0]) == 3.0
    @test getcoeff(F1, [1, 2, 0, 0, 0]) == 3.0
    @test getcoeff(F1, [3, 0, 0, 0, 0]) == 1.0
    @test getcoeff(F1, [0, 3, 0, 0, 0]) == 1.0
    using FeynmanDiagram.ComputationalGraphs:
        eval!, forwardAD, node_derivative, backAD, build_all_leaf_derivative, count_operation
    using FeynmanDiagram.Utility:
        taylorexpansion!, build_derivative_backAD!
    g1 = Graph([])
    g2 = Graph([])
    g3 = Graph([], factor=2.0)
    G3 = g1
    G4 = 1.0 * g1 * g1
    G5 = 1.0 * (3.0 * G3 + 0.5 * G4)
    G6 = (1.0 * g1 + 2.0 * g2) * (g1 + g3)
    using FeynmanDiagram.Taylor:
        TaylorSeries, getcoeff, set_variables

    set_variables("x y z", order=5)
    for G in [G3, G4, G5, G6]
        T = taylorexpansion!(G)
        T_compare = build_derivative_backAD!(G)
        for (order, coeff) in T_compare.coeffs
            @test eval!(coeff) == eval!(taylor_factorial(order) * T.coeffs[order])
        end
    end
end