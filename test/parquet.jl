@testset "Parquet" begin
    function testDiagWeigt(loopNum, chan, Kdim = 3, spin = 2, interactionTauNum = 1)

        Parquet = Builder.Parquet
        K0 = zeros(2 + loopNum)
        KinL, KoutL, KinR, KoutR = deepcopy(K0), deepcopy(K0), deepcopy(K0), deepcopy(K0)
        KinL[1] = KoutL[1] = 1
        KinR[2] = KoutR[2] = 1
        legK = [KinL, KoutL, KinR, KoutR]

        loopBasisNum = loopNum + 2
        varK = rand(Kdim, loopBasisNum)
        varT = [rand() for i in 1:2*(loopNum+1)]

        #################### DiagTree ####################################
        diag, ver4, dir, ex = Parquet.build(Float64, chan, loopNum, legK, Kdim, 3, interactionTauNum, spin)
        # the weighttype of the returned ver4 is Float64
        rootDir = DiagTree.addNode!(diag, DiagTree.ADD, :dir; child = dir, para = (0, 0, 0, 0))
        rootEx = DiagTree.addNode!(diag, DiagTree.ADD, :ex; child = ex, para = (0, 0, 0, 0))
        diag.root = [rootDir, rootEx]

        # Parquet.print_tree(ver4)

        w1 = DiagTree.evalNaive(diag, varK, varT, ParquetEval.evalPropagator)

        printstyled("naive DiagTree evaluator cost:", color = :green)
        @time DiagTree.evalNaive(diag, varK, varT, ParquetEval.evalPropagator)

        ##################### lower level subroutines  #######################################
        para = Parquet.Para(Float64, Kdim, loopNum, loopBasisNum, chan, interactionTauNum, spin)
        ver4 = Parquet.Ver4{ParquetEval.Weight}(para, loopNum, 1)

        KinL, KoutL, KinR, KoutR = varK[:, 1], varK[:, 1], varK[:, 2], varK[:, 2]
        ParquetEval.eval(ver4, varK, varT, KinL, KoutL, KinR, KoutR, 3, spin, true)

        printstyled("parquet evaluator cost:", color = :green)
        @time ParquetEval.eval(ver4, varK, varT, KinL, KoutL, KinR, KoutR, 3, spin, true)
        w2 = ver4.weight[1]

        # Parquet.print_tree(ver4)
        # DiagTree.showTree(diag, diag.root[1])
        # Parquet.showTree(ver4)
        # DiagTree.printBasisPool(diag)
        # DiagTree.printPropagator(diag)
        # println(diag.propagatorPool[1].object[2])

        @test w1[1] ≈ w2[1]
        @test w1[2] ≈ w2[2]
    end

    Parquet = Builder.Parquet
    for l = 1:3
        testDiagWeigt(l, [Parquet.T,])
        testDiagWeigt(l, [Parquet.U,])
        testDiagWeigt(l, [Parquet.S,])
        testDiagWeigt(l, [Parquet.T, Parquet.U, Parquet.S])
    end
end