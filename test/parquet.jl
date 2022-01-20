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
    # return evalG(zero(K), τin, τout)
    # return evalG(K, τin, τout)
    return 1.0
end

function evalFakeV(K)
    return 1.0
    # return evalV(K)
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
        # DiagTree.showTree(diag, uu[2].index)
        # DiagTree.showTree(diag, ud[2].index)

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

    testEval(:fake)
    testEval(:fixK)
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
        diags = Parquet.buildVer4(para, legK, chan)
        d = Parquet.groupby!(diag, nodes, :response)
        DiagTree.setroot!(diag, [d[UpUp], d[UpDown]])
        # DiagTree.showTree(diag, uu[2].index)
        # DiagTree.showTree(diag, ud[2].index)

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

    testEval(:fake)
    testEval(:fixK)
    testEval(:physical)

    #test only proper diagrams are generated if the switch is turned on
    para, diag, ver4 = testVertex4(3, [Parquet.T, Parquet.U, Parquet.S], :physical; filter = [Builder.Proper], eval = false)
    for i in 1:length(diag.basisPool)
        @test (diag.basisPool.basis[:, i] ≈ para.transferLoop) == false
    end
end

# @testset "Parquet Sigma" begin
#     function getSigma(loopNum; Kdim = 3, spin = 2, interactionTauNum = 1, filter = [], isFermi = true, subdiagram = false)
#         println("LoopNum =$loopNum Sigma Test")

#         para = Builder.GenericPara(
#             loopDim = Kdim,
#             # interactionTauNum = interactionTauNum,
#             hasTau = true,
#             innerLoopNum = loopNum,
#             totalLoopNum = loopNum + 1,
#             totalTauNum = loopNum * interactionTauNum,
#             isFermi = isFermi,
#             spin = spin,
#             weightType = Float64,
#             firstLoopIdx = 2,
#             firstTauIdx = 1,
#             filter = filter,
#             interaction = [Builder.Interaction(Builder.ChargeCharge, Builder.Instant),]
#         )

#         extK = zeros(para.totalLoopNum)
#         extK[1] = 1.0

#         Parquet = Builder.Parquet

#         varK = rand(Kdim, para.totalLoopNum)
#         varT = [rand() for i in 1:para.totalTauNum]

#         #################### DiagTree ####################################
#         diag, instant, dynamic = Parquet.buildSigma(para, extK, subdiagram)
#         # the weighttype of the returned ver4 is Float64
#         sumRoot = DiagTree.addnode!(diag, DiagTree.ADD, :sum, vcat(instant, dynamic); para = [0, 0])
#         if sumRoot.index != 0
#             push!(diag.root, sumRoot.index)
#         end

#         return para, diag, varK, varT
#     end


#     function testDiagramNumber(para, diag, varK, varT)
#         w = DiagTree.evalNaive(diag, varK, varT, evalFakePropagator)
#         factor = (1 / (2π)^para.loopDim)^para.innerLoopNum
#         num = w / factor
#         @test num[1] ≈ sigma_G2v(para.innerLoopNum, para.spin)
#     end


#     Parquet = Builder.Parquet

#     ##################  G^2*v expansion #########################################
#     for l = 1:4
#         ret = getSigma(l, spin = 1, isFermi = false, filter = [Builder.Girreducible,])
#         testDiagramNumber(ret...)
#         ret = getSigma(l, spin = 2, isFermi = false, filter = [Builder.Girreducible,])
#         testDiagramNumber(ret...)
#     end
#     # para, diag, _, _ = getSigma(3, spin = 2, isFermi = false, filter = [])
#     # for r in diag.root
#     #     DiagTree.showTree(diag, r)
#     # end

#     para, diag, varK, varT = getSigma(1, spin = 2, isFermi = false, filter = [Builder.NoFock,], subdiagram = true)
#     @test isempty(diag.root)

# end

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