@testset "Parquet" begin

    function evalG(K, τin, τout)
        # println(τBasis, ", ", varT)
        kF, β = 1.0, 1.0
        ϵ = dot(K, K) / 2 - kF^2
        τ = τout - τin
        return Spectral.kernelFermiT(τ, ϵ, β)
    end

    evalV(K) = 8π / (dot(K, K) + 1)

    function evalPropagator(idx, object, K, varT, diag)
        if idx == 1 #GPool
            # println(object.siteBasis)
            return evalG(K, varT[object.siteBasis[1]], varT[object.siteBasis[2]])
        elseif idx == 2 #VPool
            return evalV(K)
        else
            error("object with name = $(object.name) is not implemented")
        end
    end

    function testDiagWeigt(loopNum, chan, Kdim = 3, spin = 2, interactionTauNum = 1; filter = [], timing = false, eval = true)
        println("$(Int.(chan)) Channel Test")

        K0 = zeros(loopNum + 2)
        KinL, KoutL, KinR, KoutR = deepcopy(K0), deepcopy(K0), deepcopy(K0), deepcopy(K0)
        KinL[1] = KoutL[1] = 1
        KinR[2] = KoutR[2] = 1
        legK = [KinL, KoutL, KinR, KoutR]

        para = Builder.GenericPara(
            loopDim = Kdim,
            interactionTauNum = interactionTauNum,
            innerLoopNum = loopNum,
            totalLoopNum = length(KinL),
            totalTauNum = (loopNum + 1) * interactionTauNum,
            spin = spin,
            weightType = Float64,
            firstLoopIdx = 3,
            firstTauIdx = 1,
            filter = filter,
            transferLoop = KinL - KoutL
        )

        Parquet = Builder.Parquet

        F = [Parquet.U, Parquet.S]
        V = [Parquet.T, Parquet.U]

        varK = rand(Kdim, para.totalLoopNum)
        varT = [rand() for i in 1:para.totalTauNum]

        #################### DiagTree ####################################
        diag, ver4, dir, ex = Parquet.buildVer4(para, legK, chan, F, V)
        # the weighttype of the returned ver4 is Float64
        rootDir = DiagTree.addNode!(diag, DiagTree.ADD, :dir; child = dir, para = (0, 0, 0, 0))
        rootEx = DiagTree.addNode!(diag, DiagTree.ADD, :ex; child = ex, para = (0, 0, 0, 0))
        diag.root = [rootDir, rootEx]

        ver4 = Parquet.Ver4{Parquet.Weight}(para, chan, F, V)
        # Parquet.print_tree(ver4)

        if eval
            w1 = DiagTree.evalNaive(diag, varK, varT, evalPropagator)

            if timing
                printstyled("naive DiagTree evaluator cost:", color = :green)
                @time DiagTree.evalNaive(diag, varK, varT, evalPropagator)
            end

            ##################### lower level subroutines  #######################################

            KinL, KoutL, KinR, KoutR = varK[:, 1], varK[:, 1], varK[:, 2], varK[:, 2]
            Parquet.eval(para, ver4, varK, varT, [KinL, KoutL, KinR, KoutR], evalG, evalV, true)

            if timing
                printstyled("parquet evaluator cost:", color = :green)
                @time Parquet.eval(para, ver4, varK, varT, [KinL, KoutL, KinR, KoutR], evalG, evalV, true)
            end

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

        return para, diag, ver4
    end

    Parquet = Builder.Parquet
    for l = 1:3
        testDiagWeigt(l, [Parquet.T,])
        testDiagWeigt(l, [Parquet.U,])
        testDiagWeigt(l, [Parquet.S,])
        testDiagWeigt(l, [Parquet.T, Parquet.U, Parquet.S]; timing = true)
    end

    para, diag, ver4 = testDiagWeigt(2, [Parquet.T, Parquet.U, Parquet.S]; filter = [Builder.Proper], eval = false)
    for i in 1:length(diag.basisPool)
        @test (diag.basisPool.basis[:, i] ≈ para.transferLoop) == false
    end
end