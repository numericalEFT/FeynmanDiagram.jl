@testset "Partition" begin
    p = Builder.Parquet.orderedPartition(5, 2)
    expect = [[4, 1], [1, 4], [2, 3], [3, 2]]
    @test Set(p) == Set(expect)

    p = Builder.Parquet.orderedPartition(3, 2, 0)
    expect = [[3, 0], [0, 3], [1, 2], [2, 1]]
    @test Set(p) == Set(expect)
end

@testset "FindFirstIdx" begin
    # p = Builder.Parquet.orderedPartition(5, 4, 0)

    function testLoopIdx(partition, firstidx, expected)
        firstLoopIdx, total = Builder.Parquet.findFirstLoopIdx(partition, firstidx)
        @test firstLoopIdx == expected
        totalExp = sum(partition) + firstidx - 1
        @test total == totalExp
    end

    testLoopIdx([1, 1, 2, 1], 1, [1, 2, 3, 5])
    testLoopIdx([1, 1, 2, 1], 0, [0, 1, 2, 4])
    testLoopIdx([1, 0, 2, 0], 1, [1, 2, 2, 4])
    testLoopIdx([1,], 1, [1,])

    function testTauIdx(partition, isG, firstidx, tauNum, expected)
        firstIdx, total = Builder.Parquet.findFirstTauIdx(partition, isG, firstidx, tauNum)
        @test firstIdx == expected
    end
    tauNum = 1
    # isG = [false, true, false, true]
    isG = [Builder.Ver4Diag, Builder.GreenDiag, Builder.Ver4Diag, Builder.GreenDiag]
    testTauIdx([1, 1, 2, 1], isG, 1, tauNum, [1, 3, 4, 7])
    testTauIdx([1, 1, 2, 1], isG, 0, tauNum, [0, 2, 3, 6])
    testTauIdx([1, 0, 2, 0], isG, 1, tauNum, [1, 3, 3, 6])

end

function evalG(K, τin, τout)
    # println(τBasis, ", ", varT)
    kF, β = 1.0, 1.0
    ϵ = dot(K, K) / 2 - kF^2
    τ = τout - τin
    if τ ≈ 0.0
        return Spectral.kernelFermiT(-1e-8, ϵ, β)
    else
        return Spectral.kernelFermiT(τ, ϵ, β)
    end
end

evalV(K) = 8π / (dot(K, K) + 1)

function evalPropagator(idx, object, K, varT, diag)
    if idx == 1 #GPool
        return evalG(K, varT[object.siteBasis[1]], varT[object.siteBasis[2]])
        # return evalFakeG(K, varT[object.siteBasis[1]], varT[object.siteBasis[2]])
    elseif idx == 2 #VPool
        return evalV(K)
    else
        error("object with name = $(object.name) is not implemented")
    end
end

function evalGfixK(K, τin, τout)
    # println(τBasis, ", ", varT)
    K = zero(K)
    kF, β = 1.0, 1.0
    ϵ = dot(K, K) / 2 - kF^2
    τ = τout - τin
    if τ ≈ 0.0
        return Spectral.kernelFermiT(-1e-8, ϵ, β)
    else
        return Spectral.kernelFermiT(τ, ϵ, β)
    end
end

evalVfixK(K) = 1.0

function evalPropagatorfixK(idx, object, K, varT, diag)
    if idx == 1 #GPool
        return evalGfixK(K, varT[object.siteBasis[1]], varT[object.siteBasis[2]])
        # return evalFakeG(K, varT[object.siteBasis[1]], varT[object.siteBasis[2]])
    elseif idx == 2 #VPool
        return evalVfixK(K)
    else
        error("object with name = $(object.name) is not implemented")
    end
end

function evalFakePropagator(idx, object, K, varT, diag)
    return 1.0
end

function evalFakeG(K, τin, τout)
    return 1.0
end

function evalFakeV(K)
    return 1.0
end

evalFake(id::DiagramId, varK, varT) = 1.0

function evalfixK(id::GreenId, varK, varT)
    K = varK * id.extK
    tin, tout = id.extT
    return evalGfixK(K, varT[tin], varT[tout])
end

function evalfixK(id::InteractionId, varK, varT)
    K = varK * id.extK
    return evalVfixK(K)
end

function eval(id::GreenId, varK, varT)
    K = varK * id.extK
    tin, tout = id.extT
    return evalG(K, varT[tin], varT[tout])
end

function eval(id::InteractionId, varK, varT)
    K = varK * id.extK
    return evalV(K)
end


@testset "Parquet Ver4" begin
    Benchmark = ParquetNew.Benchmark
    Parquet = ParquetNew

    function getfunction(type)
        if type == :physical
            return evalPropagator, evalG, evalV
        elseif type == :fixK
            return evalPropagatorfixK, evalGfixK, evalVfixK
        elseif type == :fake
            return evalFakePropagator, evalFakeG, evalFakeV
        else
            error("not implemented")
        end
    end

    function testVertex4(loopNum, chan, type::Symbol; filter = [], timing = false, eval = true)
        println("$(Int.(chan)) Channel Test")
        Kdim, spin = 3, 2
        interactionTauNum = 1
        isFermi = true

        K0 = zeros(loopNum + 2)
        KinL, KoutL, KinR, KoutR = deepcopy(K0), deepcopy(K0), deepcopy(K0), deepcopy(K0)
        KinL[1] = KoutL[1] = 1
        KinR[2] = KoutR[2] = 1
        legK = [KinL, KoutL, KinR, KoutR]

        para = GenericPara(
            diagType = Ver4Diag,
            loopDim = Kdim,
            isFermi = isFermi,
            hasTau = true,
            innerLoopNum = loopNum,
            totalLoopNum = length(KinL),
            totalTauNum = (loopNum + 1) * interactionTauNum,
            spin = spin,
            weightType = Float64,
            firstLoopIdx = 3,
            firstTauIdx = 1,
            filter = union(filter, [Girreducible,]), #ver4 evaluation only support one-particle-irreducible diagram
            transferLoop = KinL - KoutL,
            interaction = [Interaction(ChargeCharge, Instant),]
        )

        F = [Parquet.U, Parquet.S]
        V = [Parquet.T, Parquet.U]
        All = [Parquet.T, Parquet.U, Parquet.S]
        # F = [Parquet.U]
        # V = []
        # All = [Parquet.T, Parquet.U]

        varK = rand(Kdim, para.totalLoopNum)
        varT = [rand() for i in 1:para.totalTauNum]

        #################### DiagTree ####################################
        diag, nodes = Parquet.buildVer4(para, legK, chan, F = F, V = V, All = All)
        d = Parquet.groupby!(diag, nodes, :response)
        DiagTree.setroot!(diag, [d[UpUp], d[UpDown]])
        # DiagTree.showTree(diag, d[UpUp])
        # DiagTree.showTree(diag, d[UpDown])

        ver4 = Benchmark.Ver4{Benchmark.Weight}(para, Int.(chan), Int.(F), Int.(V))
        # Parquet.print_tree(ver4)

        if eval
            # w1 = DiagTree.evalNaive(diag, varK, varT, evalPropagator)
            w1 = DiagTree.evalNaive(diag, varK, varT, evalPropagator)
            # println(w1)

            if timing
                printstyled("naive DiagTree evaluator cost:", color = :green)
                @time DiagTree.evalNaive(diag, varK, varT, evalPropagator)
            end

            ##################### lower level subroutines  #######################################

            KinL, KoutL, KinR, KoutR = varK[:, 1], varK[:, 1], varK[:, 2], varK[:, 2]
            # Benchmark.eval(para, ver4, varK, varT, [KinL, KoutL, KinR, KoutR], evalG, evalV, true)
            Benchmark.eval(para, ver4, varK, varT, [KinL, KoutL, KinR, KoutR], evalG, evalV, true)

            if timing
                printstyled("parquet evaluator cost:", color = :green)
                @time Benchmark.eval(para, ver4, varK, varT, [KinL, KoutL, KinR, KoutR], evalG, evalV, true)
            end

            w2 = ver4.weight[1]

            # Parquet.print_tree(ver4)
            # DiagTree.showTree(diag, diag.root[1])
            # Parquet.showTree(ver4)
            # DiagTree.printBasisPool(diag)
            # DiagTree.printPropagator(diag)
            # println(diag.propagatorPool[1].object[2])
            # println(w1, " vs ", w2)

            # The upup channel of charge-charge vertex4 == Direct + exchange 
            @test w1[1] ≈ w2[1] + w2[2]
            # The updown channel of charge-charge vertex4 == Direct
            @test w1[2] ≈ w2[1]
        end

        return para, diag, ver4
    end

    function testEval(type)
        for l = 1:3
            testVertex4(l, [Parquet.T,], type)
            testVertex4(l, [Parquet.U,], type)
            testVertex4(l, [Parquet.S,], type)
            testVertex4(l, [Parquet.T, Parquet.U, Parquet.S], type; timing = true)
        end
    end

    # testEval(:fake)
    # testEval(:fixK)
    testEval(:physical)

    #test only proper diagrams are generated if the switch is turned on
    para, diag, ver4 = testVertex4(3, [Parquet.T, Parquet.U, Parquet.S], :physical; filter = [Builder.Proper], eval = false)
    for i in 1:length(diag.basisPool)
        @test (diag.basisPool.basis[:, i] ≈ para.transferLoop) == false
    end
end

@testset "ParquetNew Ver4" begin
    Benchmark = ParquetNew.Benchmark
    Parquet = ParquetNew

    function getfunction(type)
        if type == :physical
            return eval, evalG, evalV
        elseif type == :fixK
            return evalfixK, evalGfixK, evalVfixK
        elseif type == :fake
            return evalFake, evalFakeG, evalFakeV
        else
            error("not implemented")
        end
    end

    function testVertex4(loopNum, chan, type::Symbol; filter = [], timing = false, toeval = true)
        println("$(Int.(chan)) Channel Test")
        Kdim, spin = 3, 2
        interactionTauNum = 1
        isFermi = true

        K0 = zeros(loopNum + 2)
        KinL, KoutL, KinR, KoutR = deepcopy(K0), deepcopy(K0), deepcopy(K0), deepcopy(K0)
        KinL[1] = KoutL[1] = 1
        KinR[2] = KoutR[2] = 1
        legK = [KinL, KoutL, KinR, KoutR]

        blocks = ParquetBlocks(phi = [PHEr, PPr], ppi = [PHr, PHEr])

        para = GenericPara(
            diagType = Ver4Diag,
            loopDim = Kdim,
            isFermi = isFermi,
            hasTau = true,
            innerLoopNum = loopNum,
            totalLoopNum = length(KinL),
            totalTauNum = (loopNum + 1) * interactionTauNum,
            spin = spin,
            weightType = Float64,
            firstLoopIdx = 3,
            firstTauIdx = 1,
            filter = union(filter, [Girreducible,]), #ver4 evaluation only support one-particle-irreducible diagram
            transferLoop = KinL - KoutL,
            interaction = [Interaction(ChargeCharge, Instant),],
            extra = blocks
        )

        varK = rand(Kdim, para.totalLoopNum)
        varT = [rand() for i in 1:para.totalTauNum]

        #################### DiagTree ####################################
        diags = Parquet.vertex4(para, legK, chan)
        diags = mergeby(diags, :response)
        # DiagTreeNew.plot_tree(diags[1])
        # DiagTreeNew.plot_tree(diags[2])

        ver4 = Benchmark.Ver4{Benchmark.Weight}(para, Int.(chan), Int.(blocks.phi), Int.(blocks.ppi))

        if toeval

            eval, evalG, evalV = getfunction(type)

            # w1 = DiagTree.evalNaive(diag, varK, varT, evalPropagator)
            evalDiagTree!(diags, eval, varK, varT)
            w1 = [diags.diagram[1].weight, diags.diagram[2].weight]
            # println(w1)

            if timing
                printstyled("naive DiagTree evaluator cost:", color = :green)
                @time evalDiagTree!(diags, eval, varK, varT)
            end

            ##################### lower level subroutines  #######################################

            KinL, KoutL, KinR, KoutR = varK[:, 1], varK[:, 1], varK[:, 2], varK[:, 2]
            # Benchmark.eval(para, ver4, varK, varT, [KinL, KoutL, KinR, KoutR], evalG, evalV, true)
            Benchmark.eval(para, ver4, varK, varT, [KinL, KoutL, KinR, KoutR], evalG, evalV, true)

            if timing
                printstyled("parquet evaluator cost:", color = :green)
                @time Benchmark.eval(para, ver4, varK, varT, [KinL, KoutL, KinR, KoutR], evalG, evalV, true)
            end

            w2 = ver4.weight[1]

            # println(w1, " vs ", w2)

            # The upup channel of charge-charge vertex4 == Direct + exchange 
            @test w1[1] ≈ w2[1] + w2[2]
            # The updown channel of charge-charge vertex4 == Direct
            @test w1[2] ≈ w2[1]
        end

        return para, diags, ver4
    end

    function testEval(type)
        for l = 1:3
            testVertex4(l, [PHr,], type)
            testVertex4(l, [PHEr,], type)
            testVertex4(l, [PPr,], type)
            testVertex4(l, [PHr, PHEr, PPr], type; timing = true)
        end
    end

    # testEval(:fake)
    # testEval(:fixK)
    testEval(:physical)

    #test only proper diagrams are generated if the switch is turned on
    # para, diag, ver4 = testVertex4(3, [Parquet.T, Parquet.U, Parquet.S], :physical; filter = [Builder.Proper], eval = false)
    # for i in 1:length(diag.basisPool)
    #     @test (diag.basisPool.basis[:, i] ≈ para.transferLoop) == false
    # end
end

@testset "Parquet Sigma" begin
    Parquet = ParquetNew
    function getSigma(loopNum; Kdim = 3, spin = 2, interactionTauNum = 1, filter = [], isFermi = true, subdiagram = false)
        println("LoopNum =$loopNum Sigma Test")

        para = GenericPara(
            diagType = SigmaDiag,
            loopDim = Kdim,
            hasTau = true,
            innerLoopNum = loopNum,
            totalLoopNum = loopNum + 1,
            totalTauNum = loopNum * interactionTauNum,
            isFermi = isFermi,
            spin = spin,
            weightType = Float64,
            firstLoopIdx = 2,
            firstTauIdx = 1,
            filter = filter,
            interaction = [Interaction(ChargeCharge, Instant),],
            extra = ParquetBlocks(phi = [PHEr, PPr], ppi = [PHr, PHEr])
        )

        extK = zeros(para.totalLoopNum)
        extK[1] = 1.0

        varK = rand(Kdim, para.totalLoopNum)
        varT = [rand() for i in 1:para.totalTauNum]

        #################### DiagTree ####################################
        diag = Parquet.sigma(para, extK, subdiagram)
        diag = mergeby(diag)
        # print_tree(diag.diagram[1])

        return para, diag.diagram[1], varK, varT
    end


    function testDiagramNumber(para, diag, varK, varT)
        # w = DiagTree.evalNaive(diag, varK, varT, evalFakePropagator)
        w = evalDiagTree!(diag, evalFake, varK, varT)
        # plot_tree(diag, maxdepth = 7)
        factor = (1 / (2π)^para.loopDim)^para.innerLoopNum
        num = w / factor
        @test num ≈ sigma_G2v(para.innerLoopNum, para.spin)
    end


    ##################  G^2*v expansion #########################################
    for l = 1:4
        # ret = getSigma(l, spin = 1, isFermi = false, filter = [Builder.Girreducible,])
        # testDiagramNumber(ret...)
        ret = getSigma(l, spin = 2, isFermi = false, filter = [Builder.Girreducible,])
        testDiagramNumber(ret...)
    end

    # para, diag, varK, varT = getSigma(1, spin = 2, isFermi = false, filter = [Builder.NoFock,], subdiagram = true)
    # @test isempty(diag.root)

end

# @testset "Green" begin
#     Parquet = Builder.Parquet

#     function buildG(loopNum, extT, firstTauIdx; Kdim = 3, spin = 2, interactionTauNum = 1, filter = [], isFermi = true)
#         para = Builder.GenericPara(
#             loopDim = Kdim,
#             # interactionTauNum = interactionTauNum,
#             hasTau = true,
#             innerLoopNum = loopNum,
#             totalLoopNum = loopNum + 1,
#             totalTauNum = loopNum * interactionTauNum + 2,
#             isFermi = isFermi,
#             spin = spin,
#             weightType = Float64,
#             firstLoopIdx = 2,
#             firstTauIdx = firstTauIdx,
#             filter = filter,
#             interaction = [Builder.Interaction(Builder.ChargeCharge, Builder.Instant),]
#         )
#         extK = zeros(para.totalLoopNum)
#         extK[1] = 1.0
#         diag, Gidx = Parquet.buildG(para, extK, extT)
#         return diag, Gidx
#     end
#     # diag, Gidx = buildG(2, [1, 2], 3; filter = [])
#     # DiagTree.showTree(diag, Gidx)

#     # If G is irreducible, then only loop-0 G exist
#     diag, G = buildG(1, [1, 2], 3; filter = [Builder.Girreducible,])
#     @test G.index == 0

#     # If Fock diagram is not allowed, then one-loop G diagram should not be exist
#     diag, G = buildG(1, [1, 2], 3; filter = [Builder.NoFock,])
#     @test G.index == 0
#     # Even if Fock diagram is not allowed, then loopNum>=1 G diagram can exist
#     diag, G = buildG(2, [1, 2], 3; filter = [Builder.NoFock,])
#     @test G.index > 0

# end


@testset "Parquet Vertex3" begin
    Parquet = ParquetNew
    function getGamma3(loopNum; Kdim = 3, spin = 2, interactionTauNum = 1, filter = [Girreducible, Proper,], isFermi = true, subdiagram = false)
        println("LoopNum =$loopNum Vertex3 Test")

        para = GenericPara(
            diagType = Ver3Diag,
            loopDim = Kdim,
            innerLoopNum = loopNum,
            isFermi = isFermi,
            hasTau = true,
            filter = filter,
            interaction = [Interaction(ChargeCharge, Instant),]
        )

        K0 = zeros(para.totalLoopNum)
        KinL, Q = deepcopy(K0), deepcopy(K0)
        Q[1] = 1
        KinL[2] = 1
        legK = [Q, KinL]

        varK = rand(Kdim, para.totalLoopNum)
        varT = [rand() for i in 1:para.totalTauNum]

        #################### DiagTree ####################################
        vertex3 = Parquet.vertex3(para, legK)
        diag = mergeby(vertex3)
        # print_tree(diag.diagram[1])

        return para, diag.diagram[1], varK, varT
    end


    function testDiagramNumber(para, diag, varK, varT)
        # w = DiagTree.evalNaive(diag, varK, varT, evalFakePropagator)
        w = evalDiagTree!(diag, evalFake, varK, varT)
        # plot_tree(diag, maxdepth = 9)
        factor = (1 / (2π)^para.loopDim)^para.innerLoopNum
        num = w / factor
        @test num ≈ gamma3_G2v(para.innerLoopNum, para.spin)
    end


    ##################  G^2*v expansion #########################################
    for l = 1:3
        # ret = getSigma(l, spin = 1, isFermi = false, filter = [Builder.Girreducible,])
        # testDiagramNumber(ret...)
        ret = getGamma3(l, isFermi = false, filter = [Girreducible, Proper])
        testDiagramNumber(ret...)
    end

    # para, diag, varK, varT = getSigma(1, spin = 2, isFermi = false, filter = [Builder.NoFock,], subdiagram = true)
    # @test isempty(diag.root)

end