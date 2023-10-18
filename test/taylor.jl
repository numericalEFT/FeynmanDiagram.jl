using FeynmanDiagram: Taylor as Taylor

@testset verbose = true "TaylorSeries" begin
    using FeynmanDiagram.Taylor:
        getcoeff
    a, b, c, d, e = set_variables("a b c d e", order=3)
    F1 = (a + b) * (a + b) * (a + b)
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
end