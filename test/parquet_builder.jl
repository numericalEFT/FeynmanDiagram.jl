@testset "Partition" begin
    p = Parquet.orderedPartition(5, 2)
    expect = [[4, 1], [1, 4], [2, 3], [3, 2]]
    @test Set(p) == Set(expect)

    p = Parquet.orderedPartition(3, 2, 0)
    expect = [[3, 0], [0, 3], [1, 2], [2, 1]]
    @test Set(p) == Set(expect)
end

@testset "FindFirstIdx" begin
    function testLoopIdx(partition, firstidx, expected)
        firstLoopIdx, total = Parquet.findFirstLoopIdx(partition, firstidx)
        @test firstLoopIdx == expected
        totalExp = sum(partition) + firstidx - 1
        @test total == totalExp
    end

    testLoopIdx([1, 1, 2, 1], 1, [1, 2, 3, 5])
    testLoopIdx([1, 1, 2, 1], 0, [0, 1, 2, 4])
    testLoopIdx([1, 0, 2, 0], 1, [1, 2, 2, 4])
    testLoopIdx([1,], 1, [1,])

    function testTauIdx(partition, isG, firstidx, tauNum, expected)
        firstIdx, total = Parquet.findFirstTauIdx(partition, isG, firstidx, tauNum)
        @test firstIdx == expected
    end
    tauNum = 1
    # isG = [false, true, false, true]
    isG = [Ver4Diag, GreenDiag, Ver4Diag, GreenDiag]
    testTauIdx([1, 1, 2, 1], isG, 1, tauNum, [1, 3, 4, 7])
    testTauIdx([1, 1, 2, 1], isG, 0, tauNum, [0, 2, 3, 6])
    testTauIdx([1, 0, 2, 0], isG, 1, tauNum, [1, 3, 3, 6])

end

@testset "Filter" begin

    # for G irreducible diagrams, only 0-loop G is allowed
    @test Parquet.isValidG([Girreducible,], 0) == true
    @test Parquet.isValidG([Girreducible,], 1) == false
    @test Parquet.isValidG([Girreducible,], 2) == false

    # for Fock irreducible diagrams, only 0-loop or 2, 3, 4...-loop G is allowed
    @test Parquet.isValidG([NoFock,], 0) == true
    @test Parquet.isValidG([NoFock,], 1) == true
    #one-loop G diagram becomes invalid only if both Hartree and Fock are filtered
    @test Parquet.isValidG([NoFock, NoHartree], 1) == false
    @test Parquet.isValidG([NoFock, NoHartree], 2) == true

    # for G irreducible diagrams, no sigma subdiagram is allowed
    @test Parquet.isValidSigma([Girreducible,], 0, true) == false
    @test Parquet.isValidSigma([Girreducible,], 1, true) == false
    @test Parquet.isValidSigma([Girreducible,], 2, true) == false

    @test Parquet.isValidSigma([Girreducible,], 0, false) == false
    @test Parquet.isValidSigma([Girreducible,], 1, false) == true
    @test Parquet.isValidSigma([Girreducible,], 2, false) == true

    # for Fock irreducible diagrams, no Fock sigma subdiagram is allowed
    @test Parquet.isValidSigma([NoFock,], 0, true) == false
    #one-loop sigma diagram can be either Hartree or Fock diagram
    #one-loop sigma sub-diagram becomes invalid only if both Hartree and Fock are filtered
    @test Parquet.isValidSigma([NoFock,], 1, true) == true
    @test Parquet.isValidSigma([NoFock, NoHartree], 1, true) == false
    @test Parquet.isValidSigma([NoFock, NoHartree], 2, true) == true

    @test Parquet.isValidSigma([NoFock,], 0, false) == false
    @test Parquet.isValidSigma([NoFock,], 1, false) == true
    @test Parquet.isValidSigma([NoFock, NoHartree], 1, false) == true
    @test Parquet.isValidSigma([NoFock,], 2, false) == true
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

evalGfixK(K, τin, τout) = evalG(zero(K), τin, τout)
evalVfixK(K) = 1.0

evalFakeG(K, τin, τout) = 1.0
evalFakeV(K) = 1.0

################## api for expression tree ##############################
evalPropagator(id::BareGreenId, K, extT, varT) = evalG(K, varT[extT[1]], varT[extT[2]])
evalPropagator(id::BareInteractionId, K, extT, varT) = evalV(K)
evalPropagatorfixK(id::BareGreenId, K, extT, varT) = evalGfixK(K, varT[extT[1]], varT[extT[2]])
evalPropagatorfixK(id::BareInteractionId, K, extT, varT) = evalVfixK(K)
evalFakePropagator(id::PropagatorId, K, extT, varT) = 1.0


@testset "ParquetNew Ver4" begin
    Benchmark = Parquet.Benchmark

    function getfunction(type)
        if type == :physical
            return evalG, evalV, evalPropagator
        elseif type == :fixK
            return evalGfixK, evalVfixK, evalPropagatorfixK
        elseif type == :fake
            return evalFakeG, evalFakeV, evalFakePropagator
        else
            error("not implemented")
        end
    end

    function testVertex4(loopNum, chan, type::Symbol; filter=[NoHartree,], timing=false, toeval=true)
        println("$(Int.(chan)) Channel Test")
        Kdim, spin = 3, 2
        interactionTauNum = 1
        isFermi = true

        K0 = zeros(loopNum + 2)
        KinL, KoutL, KinR, KoutR = deepcopy(K0), deepcopy(K0), deepcopy(K0), deepcopy(K0)
        KinL[1] = KoutL[1] = 1
        KinR[2] = KoutR[2] = 1
        legK = [KinL, KoutL, KinR, KoutR]

        blocks = ParquetBlocks(phi=[PHEr, PPr], ppi=[PHr, PHEr])

        para = DiagParaF64(
            type=Ver4Diag,
            loopDim=Kdim,
            isFermi=isFermi,
            hasTau=true,
            innerLoopNum=loopNum,
            totalLoopNum=length(KinL),
            totalTauNum=(loopNum + 1) * interactionTauNum,
            spin=spin,
            firstLoopIdx=3,
            firstTauIdx=1,
            filter=union(filter, [Girreducible,]), #ver4 evaluation only support one-particle-irreducible diagram
            transferLoop=KinL - KoutL,
            interaction=[Interaction(ChargeCharge, Instant),],
            extra=blocks
        )

        varK = rand(Kdim, para.totalLoopNum)
        varT = [rand() for i in 1:para.totalTauNum]

        #################### DiagTree ####################################
        diags = Parquet.vertex4(para, legK, chan)
        diags = mergeby(diags, :response)
        # DiagTreeNew.plot_tree(diags[1])
        # DiagTreeNew.plot_tree(diags[2])

        ################### ExprTree ###################################
        tree = ExprTree.build(diags.diagram)
        # println("root", root)

        ################### original Parquet builder ###################################
        ver4 = Benchmark.Ver4{Benchmark.Weight}(para, Int.(chan), Int.(blocks.phi), Int.(blocks.ppi))


        if toeval

            evalG, evalV, evalPropagator = getfunction(type)

            # w1 = DiagTree.evalNaive(diag, varK, varT, evalPropagator)
            DiagTree.evalKT!(diags, varK, varT; eval=evalPropagator)
            w1 = [diags.diagram[1].weight, diags.diagram[2].weight]
            if timing
                printstyled("naive DiagTree evaluator cost:", color=:green)
                @time DiagTree.evalKT!(diags, varK, varT; eval=evalPropagator)
            end

            ExprTree.evalNaive!(tree, varK, varT; eval=evalPropagator)
            w1e = [tree[1], tree[2]]
            if timing
                printstyled("naive ExprTree cost:", color=:green)
                @time ExprTree.evalKT!(tree, varK, varT; eval=evalPropagator)
            end


            # optdiags = DiagTree.optimize!(diags.diagram)
            opttree = ExprTree.build(diags.diagram)
            ExprTree.evalKT!(opttree, varK, varT; eval=evalPropagator)
            w1eopt = [opttree[1], opttree[2]]

            if timing
                printstyled("naive optimized ExprTree cost:", color=:green)
                @time ExprTree.evalKT!(opttree, varK, varT; eval=evalPropagator)
            end

            ##################### lower level subroutines  #######################################

            KinL, KoutL, KinR, KoutR = varK[:, 1], varK[:, 1], varK[:, 2], varK[:, 2]
            legK = [KinL, KoutL, KinR, KoutR]
            # Benchmark.eval(para, ver4, varK, varT, [KinL, KoutL, KinR, KoutR], evalG, evalV, true)
            Benchmark.eval(para, ver4, varK, varT, legK, evalG, evalV, true)

            if timing
                printstyled("parquet evaluator cost:", color=:green)
                # @btime sin(p, ver4, var) setup = (x = rand())
                # @time Benchmark.eval(para, ver4, varK, varT, [KinL, KoutL, KinR, KoutR], evalG, evalV, true)
                @time Benchmark.eval(para, ver4, varK, varT, legK, evalG, evalV, true)
                # @btime Benchmark.eval(p, v4, vK, vT, lK, eG, eV, flag) setup = (p = para, v4 = ver4, vK = varK, vT = varT, l = legK, eG = evalG, eV = evalV, flag = true)
            end

            w2 = ver4.weight[1]

            # println(w1, " vs ", w1e, " vs ", w2)

            @assert w1 ≈ w1e
            @assert w1 ≈ w1eopt

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
            testVertex4(l, [PHr, PHEr, PPr], type; timing=true)
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
    function getSigma(loopNum; Kdim=3, spin=2, interactionTauNum=1, filter=[NoHartree,], isFermi=true, subdiagram=false)
        println("LoopNum =$loopNum Sigma Test")

        para = DiagParaF64(
            type=SigmaDiag,
            loopDim=Kdim,
            hasTau=true,
            innerLoopNum=loopNum,
            totalLoopNum=loopNum + 1,
            totalTauNum=loopNum * interactionTauNum,
            isFermi=isFermi,
            spin=spin,
            firstLoopIdx=2,
            firstTauIdx=1,
            filter=filter,
            interaction=[Interaction(ChargeCharge, Instant),],
            extra=ParquetBlocks(phi=[PHEr, PPr], ppi=[PHr, PHEr])
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
        w = DiagTree.evalKT!(diag, varK, varT; eval=evalFakePropagator)
        # plot_tree(diag, maxdepth = 7)
        # factor = (1 / (2π)^para.loopDim)^para.innerLoopNum
        factor = 1.0
        num = w / factor
        @test num * (-1)^(para.innerLoopNum) ≈ Parquet.Benchmark.count_sigma_G2v(para.innerLoopNum, para.spin)
    end


    ##################  G^2*v expansion #########################################
    for l = 1:4
        # ret = getSigma(l, spin = 1, isFermi = false, filter = [Builder.Girreducible,])
        # testDiagramNumber(ret...)
        ret = getSigma(l, spin=2, isFermi=false, filter=[NoHartree, Girreducible,])
        testDiagramNumber(ret...)
    end

    # para, diag, varK, varT = getSigma(1, spin = 2, isFermi = false, filter = [Builder.NoFock,], subdiagram = true)
    # @test isempty(diag.root)

end

@testset "Green" begin
    function buildG(loopNum, extT; Kdim=3, spin=2, interactionTauNum=1, filter=[NoHartree,], isFermi=true)
        para = DiagParaF64(
            type=GreenDiag,
            loopDim=Kdim,
            hasTau=true,
            innerLoopNum=loopNum,
            isFermi=isFermi,
            spin=spin,
            filter=filter,
            interaction=[Interaction(ChargeCharge, Instant),]
        )
        extK = zeros(para.totalLoopNum)
        extK[1] = 1.0
        if Parquet.isValidG(para)
            G = Parquet.green(para, extK, extT)
            return G
        else
            return nothing
        end
    end
    # diag, Gidx = buildG(2, [1, 2], 3; filter = [])
    # DiagTree.showTree(diag, Gidx)

    # If G is irreducible, then only loop-0 G exist for main diagram, and no G exist for subdiagram
    G = buildG(0, [1, 2]; filter=[NoHartree, Girreducible,])
    @test G isa Diagram
    G = buildG(1, [1, 2]; filter=[NoHartree, Girreducible,])
    @test isnothing(G)
    G = buildG(2, [1, 2]; filter=[NoHartree, Girreducible,])
    @test isnothing(G)

    # If Fock diagram is not allowed, then one-loop G diagram should not be exist for subdiagram
    G = buildG(0, [1, 2]; filter=[NoHartree, NoFock,])
    @test G isa Diagram
    G = buildG(1, [1, 2]; filter=[NoHartree, NoFock,])
    @test isnothing(G)
    G = buildG(2, [1, 2]; filter=[NoHartree, NoFock,]) #high order subdiagram is allowed
    @test G isa Diagram

end


@testset "Parquet Vertex3" begin
    function getGamma3(loopNum; Kdim=3, spin=2, interactionTauNum=1, filter=[NoHartree, Girreducible, Proper,], isFermi=true, subdiagram=false)
        println("LoopNum =$loopNum Vertex3 Test")

        para = DiagParaF64(
            type=Ver3Diag,
            loopDim=Kdim,
            innerLoopNum=loopNum,
            isFermi=isFermi,
            hasTau=true,
            filter=filter,
            interaction=[Interaction(ChargeCharge, Instant),]
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
        w = DiagTree.evalKT!(diag, varK, varT; eval=evalFakePropagator)
        # plot_tree(diag, maxdepth = 9)
        # factor = (1 / (2π)^para.loopDim)^para.innerLoopNum
        factor = 1.0
        num = w / factor
        @test num * (-1)^(para.innerLoopNum) ≈ Parquet.Benchmark.count_ver3_G2v(para.innerLoopNum, para.spin)
    end


    ##################  G^2*v expansion #########################################
    for l = 1:3
        # ret = getSigma(l, spin = 1, isFermi = false, filter = [Builder.Girreducible,])
        # testDiagramNumber(ret...)
        ret = getGamma3(l, isFermi=false, filter=[NoHartree, Girreducible, Proper])
        testDiagramNumber(ret...)
    end

    # para, diag, varK, varT = getSigma(1, spin = 2, isFermi = false, filter = [Builder.NoFock,], subdiagram = true)
    # @test isempty(diag.root)

end


@testset "Parquet Polarization" begin
    function getPolar(loopNum; Kdim=3, spin=2, interactionTauNum=1, filter=[NoHartree, Girreducible,], isFermi=true, subdiagram=false)
        println("LoopNum =$loopNum Polarization Test")

        para = DiagParaF64(
            type=PolarDiag,
            loopDim=Kdim,
            innerLoopNum=loopNum,
            isFermi=isFermi,
            hasTau=true,
            filter=filter,
            interaction=[Interaction(ChargeCharge, Instant),]
        )

        Q = zeros(para.totalLoopNum)
        Q[1] = 1

        varK = rand(Kdim, para.totalLoopNum)
        varT = [rand() for i in 1:para.totalTauNum]

        #################### DiagTree ####################################
        diag = Parquet.polarization(para, Q)
        # print_tree(diag.diagram[1])
        return para, diag, varK, varT
    end

    # Test polarization Parquet builder when filter 'Proper' is specified explicitly
    getPolar(1, filter=[Proper, NoHartree, NoFock,])

    ##################  G^2*v expansion #########################################
    for l = 1:4
        para, diag, varK, varT = getPolar(l, isFermi=false, filter=[NoHartree, Girreducible,])
        diag = mergeby(diag).diagram[1]
        w = DiagTree.evalKT!(diag, varK, varT; eval=evalFakePropagator)
        # factor = (1 / (2π)^para.loopDim)^para.innerLoopNum
        factor = 1.0
        num = w / factor
        # println(num * para.spin)
        @test num * para.spin * (-1)^(para.innerLoopNum - 1) ≈ Parquet.Benchmark.count_polar_G2v(para.innerLoopNum, para.spin)
    end

    ##################  g^2*v expansion #########################################
    for l = 1:4
        para, diag, varK, varT = getPolar(l, isFermi=false, filter=[NoHartree, NoFock,])
        diag = mergeby(diag).diagram[1]
        w = DiagTree.evalKT!(diag, varK, varT, eval=evalFakePropagator)
        # factor = (1 / (2π)^para.loopDim)^para.innerLoopNum
        factor = 1.0
        num = w / factor
        # println(num * para.spin)
        @test num * para.spin * (-1)^(para.innerLoopNum - 1) ≈ Parquet.Benchmark.count_polar_g2v_noFock(para.innerLoopNum, para.spin)
    end

    ##################  g^2*v expansion for the upup polarization #########################################
    for l = 1:4
        para, diag, varK, varT = getPolar(l, isFermi=false, filter=[NoHartree, NoFock,])
        w = DiagTree.evalKT!(diag.diagram[1], varK, varT, eval=evalFakePropagator)
        # factor = (1 / (2π)^para.loopDim)^para.innerLoopNum
        factor = 1.0
        num = w / factor
        # println(num * para.spin)
        # println("$diag")
        @test num * para.spin * (-1)^(para.innerLoopNum - 1) ≈ Parquet.Benchmark.count_polar_g2v_noFock_upup(para.innerLoopNum, para.spin)
    end
end