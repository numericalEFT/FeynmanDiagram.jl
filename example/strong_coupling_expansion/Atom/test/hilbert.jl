@testset "Spin-1/2 Fock" begin
    # test six states with N=2
    # (UP and DOWN)
    n = [([1, 1], [0, 0]), ([0, 1], [0, 1]), ([1, 0], [0, 1]), ([0, 1], [1, 0]), ([1, 0], [1, 0]), ([0, 0], [1, 1])]
    idx = [3, 5, 6, 9, 10, 12]

    for i in 1:length(n)
    @test Hilbert.state2idx(n[i]...) == idx[i]
    @test Hilbert.idx2state(idx[i], 2) == n[i]
    @test Hilbert.occupation(idx[i], 2, 1, UP) == n[i][1][1]
    @test Hilbert.occupation(idx[i], 2, 1, DOWN) == n[i][2][1]
    @test Hilbert.occupation(idx[i], 2, 2, UP) == n[i][1][2]
    @test Hilbert.occupation(idx[i], 2, 2, DOWN) == n[i][2][2]
    @test Hilbert.occupationTotal(idx[i], 2) == (sum(n[i][1]), sum(n[i][2]))

end
    basis = Hilbert.BinaryFock(2, 2, :all)
    @test basis.idx == idx
    # basis = BinaryFock(1, :all, :all)
    # println(creation(basis, 1, UP))
    # println(creation(basis, 1, UP))
end