@testset "Parameter" begin
    p = GenericPara(diagType = Ver4Diag, innerLoopNum = 1)
    q = GenericPara(diagType = Ver4Diag, innerLoopNum = 2)
    a = GenericPara(diagType = Ver4Diag, innerLoopNum = 2)

    @test p != q
    @test q == a

    aa = reconstruct(a, transferLoop = [0.0, 0.0, 0.0])
    @test a != aa

    # reconstruct with the same diagType but different interaction
    aaa = reconstruct(a, interaction = [])
    @test a != aaa

    # reconstruct with different diagram type leads to different parameter
    #reconstructed GenericPara uses the old parameters such as firstLoopIdx, totalLoopNum etc., which are different between diagTypes
    s = reconstruct(a, diagType = SigmaDiag)
    ss = GenericPara(diagType = SigmaDiag, innerLoopNum = 2)
    @test s != ss

end