@testset "DiagTree" begin
    # Write your tests here.
    propagators = Vector{DiagTree.PropagatorKT}(undef, 0)
    Gsymmetry = [:mirror, :particlehole]
    Wsymmetry = [:mirror, :timereversal]
    kbasis = [1, 0, 1, 1]
    Tidx = [(1, 2), (1, 2), (2, 1), (2, 1)]
    Kidx = [kbasis, -kbasis, kbasis, -kbasis]
    N = length(Tidx)

    ############# Test G with symmetry ############################
    idx = [DiagTree.addPropagator(propagators, :G, 1, Tidx[i], Kidx[i], Gsymmetry) for i = 1:N]
    @test idx == [1, 1, 1, 1]
    ############# Test G without symmetry ############################
    idx = [DiagTree.addPropagator(propagators, :G, 1, Tidx[i], Kidx[i], []) for i = 1:N]
    @test idx == [1, 2, 3, 4]
    ############# Test W with symmetry ############################
    idx = [DiagTree.addPropagator(propagators, :W, 1, Tidx[i], Kidx[i], Wsymmetry) for i = 1:N]
    @test idx == [5, 5, 5, 5]
    ############# Test W without symmetry ############################
    idx = [DiagTree.addPropagator(propagators, :W, 1, Tidx[i], Kidx[i], []) for i = 1:N]
    @test idx == [5, 6, 7, 8]
end