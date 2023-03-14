@testset "LoopPool" begin
    loopbasis = [[1.0, 1.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0], [1.0, 0.0, -1.0, 0.0]]
    loopPool = FrontEnds.LoopPool(:K, 3, loopbasis)
    @test size(loopPool) == (length(loopbasis),)
    @test loopPool[1] == [1.0, 1.0, 0.0, 0.0]
    @test loopPool[end] == [1.0, 0.0, -1.0, 0.0]
    loopPool[2] = [1.0, 0.0, -1.0, 0.0]
    @test loopPool[2] == [1.0, 0.0, -1.0, 0.0]

    dim, N = 3, 4
    loopPool = FrontEnds.LoopPool(:K, dim, N, Float64)
    basis1 = [1.0, 0.0, 0.0, 1.0]
    basis2 = [1.0, 1.0, 0.0, 0.0]
    basis3 = [1.0, 0.0, -1.0, 1.0]
    idx1 = FrontEnds.append(loopPool, basis1)
    idx2 = FrontEnds.append(loopPool, basis2)
    idx3 = FrontEnds.append(loopPool, basis2)
    idx4 = FrontEnds.append(loopPool, basis1)
    idx5 = FrontEnds.append(loopPool, basis3)
    @test length(loopPool) == 3
    @test idx1 == idx4
    @test idx2 == idx3

    varK = rand(dim, N)
    FrontEnds.update(loopPool, varK)
    @test FrontEnds.loop(loopPool, 1) ≈ varK * basis1
    @test FrontEnds.loop(loopPool, 2) ≈ varK * basis2
    @test FrontEnds.loop(loopPool, 3) ≈ varK * basis3
end

@testset "LabelProduct" begin
    flavors = [1, 2, 3]
    tau_labels = collect(1:5)
    loopbasis = [[1.0, 1.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0], [1.0, 0.0, -1.0, 0.0]]
    loopPool = FrontEnds.LoopPool(:K, 3, loopbasis)
    labelProd = LabelProduct(flavors, tau_labels, loopPool)

    @test length(labelProd) == 3 * 5 * 5
    @test size(labelProd) == (3, 5, 5)
    @test FrontEnds.index_to_linear(labelProd, 2, 4, 3) == 2 + 3 * 3 + 3 * 5 * 2
    @test FrontEnds.linear_to_index(labelProd, 41) == (2, 4, 3)
    @test FrontEnds.linear_to_index(labelProd.dims, 41) == (2, 4, 3)
    @test labelProd[38] == labelProd[2, 3, 3] == (2, 3, [0.0, 0.0, 1.0, 0.0])

    @test eltype(typeof(labelProd)) == (eltype(typeof(flavors)), eltype(typeof(tau_labels)), eltype(typeof(loopPool)))
    @test FrontEnds._find_label(typeof(labelProd.labels), FrontEnds.LoopPool) == 3
end