@testset "DiagTree" begin
    DiagTree = GWKT.DiagTree
    # Write your tests here.
    Gsymmetry = [:mirror, :particlehole]
    Wsymmetry = [:mirror, :timereversal]
    kbasis = [1, 0, 1, 1]
    Tidx = [(1, 2), (1, 2), (2, 1), (2, 1)]
    Kidx = [kbasis, -kbasis, kbasis, -kbasis]
    N = length(Tidx)
    Gtype, Wtype = 1, 2

    ############# Test G with symmetry ############################
    diag = DiagTree.Diagrams{Float64}()
    idx = [DiagTree.addPropagator!(diag, Gtype, 1, Kidx[i], Tidx[i], Gsymmetry)[1] for i = 1:N]
    @test idx == [1, 1, 1, 1]
    ############# Test G without symmetry ############################
    diag = DiagTree.Diagrams{Float64}()
    idx = [DiagTree.addPropagator!(diag, Gtype, 1, Kidx[i], Tidx[i], [])[1] for i = 1:N]
    @test idx == [1, 2, 3, 4]
    ############# Test W with symmetry ############################
    diag = DiagTree.Diagrams{Float64}()
    idx = [DiagTree.addPropagator!(diag, Wtype, 1, Kidx[i], Tidx[i], Wsymmetry)[1] for i = 1:N]
    @test idx == [1, 1, 1, 1]
    ############# Test W without symmetry ############################
    diag = DiagTree.Diagrams{Float64}()
    idx = [DiagTree.addPropagator!(diag, Wtype, 1, Kidx[i], Tidx[i], [])[1] for i = 1:N]
    @test idx == [1, 2, 3, 4]
end

@testset "Diagrams" begin
    """
        k1-k3                     k2+k3 
        |                         | 
    t1.L ↑     t1.L       t2.L     ↑ t2.L
        |-------------->----------|
        |       |  k3+k4   |      |
        |   v   |          |  v   |
        |       |    k4    |      |
        |--------------<----------|
    t1.L ↑    t1.L        t2.L     ↑ t2.L
        |                         | 
        k1                        k2
    """
    DiagTree = GWKT.DiagTree
    Gsym = [:mirror]
    Wsym = [:mirror, :timereversal]
    Gtype, Wtype = 1, 2
    diag = DiagTree.Diagrams{Float64}()
    g1, _ = DiagTree.addPropagator!(diag, Gtype, 0, [0, 0, 1, 1], [1, 2], Gsym)
    g2, _ = DiagTree.addPropagator!(diag, Gtype, 0, [0, 0, 0, 1], [2, 1], Gsym)
    v1d, _ = DiagTree.addPropagator!(diag, Wtype, 0, [0, 0, 1, 0], [1, 1], Wsym)
    v1e, _ = DiagTree.addPropagator!(diag, Wtype, 0, [1, 0, -1, -1], [1, 1], Wsym)
    v2d, _ = DiagTree.addPropagator!(diag, Wtype, 0, [0, 0, 1, 0], [2, 2], Wsym)
    v2e, _ = DiagTree.addPropagator!(diag, Wtype, 0, [0, 1, 0, -1], [2, 2], Wsym)

    MUL, ADD = 0, 1
    factor = 1.0 / (2π)^3

    # root = addNode!(diag, 0)
    # addChild!(diag, root)
    gg_n = addNode!(diag, MUL, 1.0, [g1, g2], [])
    vdd = addNode!(diag, MUL, 1.0, [v1d, v2d], [])
    vde = addNode!(diag, MUL, 1.0, [v1d, v2e], [])
    ved = addNode!(diag, MUL, 1.0, [v1e, v2d], [])
    vee = addNode!(diag, MUL, 1.0, [v1e, v2e], [])
    vsum = addNode!(diag, ADD, 1.0, [], [vdd, vde, ved, vee])
    root = addNode!(diag, MUL, factor, [], [gg_n, vsum])

end