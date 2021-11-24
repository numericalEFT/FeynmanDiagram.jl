@testset "DiagTree" begin
    DiagTree = GWKT.DiagTree
    # Write your tests here.
    Gsymmetry = [:mirror, :particlehole]
    Wsymmetry = [:mirror, :timereversal]
    kbasis = [1, 0, 1, 1]
    Tidx = [(1, 2), (1, 2), (2, 1), (2, 1)]
    Kidx = [kbasis, -kbasis, kbasis, -kbasis]
    N = length(Tidx)

    ############# Test G with symmetry ############################
    diag = DiagTree.Diagrams{Float64}(Gsymmetry, Wsymmetry)
    idx = [DiagTree.addPropagator!(diag, :G, 1, Kidx[i], Tidx[i])[1] for i = 1:N]
    @test idx == [1, 1, 1, 1]
    ############# Test G without symmetry ############################
    diag = DiagTree.Diagrams{Float64}([], Wsymmetry)
    idx = [DiagTree.addPropagator!(diag, :G, 1, Kidx[i], Tidx[i])[1] for i = 1:N]
    @test idx == [1, 2, 3, 4]
    ############# Test W with symmetry ############################
    diag = DiagTree.Diagrams{Float64}(Gsymmetry, Wsymmetry)
    idx = [DiagTree.addPropagator!(diag, :W, 1, Kidx[i], Tidx[i])[1] for i = 1:N]
    @test idx == [1, 1, 1, 1]
    ############# Test W without symmetry ############################
    diag = DiagTree.Diagrams{Float64}(Gsymmetry, [])
    idx = [DiagTree.addPropagator!(diag, :W, 1, Kidx[i], Tidx[i])[1] for i = 1:N]
    @test idx == [1, 2, 3, 4]
end