@testset "Parameter" begin
    p = GenericPara(diagType = Ver4Diag, innerLoopNum = 1)
    q = GenericPara(diagType = Ver4Diag, innerLoopNum = 2)
    a = GenericPara(diagType = Ver4Diag, innerLoopNum = 2)

    @test p != q
    @test q == a

    aa = reconstruct(a, transferLoop = [0.0, 0.0, 0.0])
    @test a != aa

    aaa = reconstruct(a, interaction = [])
    @test a != aaa
end