@testset "Parameter" begin
    p = DiagParaF64(type=Ver4Diag, innerLoopNum=1)
    q = DiagParaF64(type=Ver4Diag, innerLoopNum=2)
    a = DiagParaF64(type=Ver4Diag, innerLoopNum=2)

    @test p != q
    @test q == a

    aa = reconstruct(a, transferLoop=[0.0, 0.0, 0.0])
    @test a != aa

    # reconstruct with the same type but different interaction
    aaa = reconstruct(a, interaction=[])
    @test a != aaa

    # reconstruct with different diagram type leads to different parameter
    #reconstructed DiagPara uses the old parameters such as firstLoopIdx, totalLoopNum etc., which are different between types
    s = reconstruct(a, type=SigmaDiag)
    ss = DiagParaF64(type=SigmaDiag, innerLoopNum=2)
    @test s != ss

end