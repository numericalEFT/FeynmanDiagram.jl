@testset "More on Parquet" begin
    function testVer4(loopNum, chan, interaction, Kdim = 3, spin = 2; hasTau = true, filter = [], timing = false, eval = true)
        println("$(Int.(chan)) Channel Test")

        K0 = zeros(loopNum + 2)
        KinL, KoutL, KinR, KoutR = deepcopy(K0), deepcopy(K0), deepcopy(K0), deepcopy(K0)
        KinL[1] = KoutL[1] = 1
        KinR[2] = KoutR[2] = 1
        legK = [KinL, KoutL, KinR, KoutR]

        interactionTauNum = Builder.interactionTauNum(hasTau, interaction)

        para = Builder.GenericPara(
            loopDim = Kdim,
            innerLoopNum = loopNum,
            totalLoopNum = length(KinL),
            spin = spin,
            weightType = Float64,
            firstLoopIdx = 3,
            filter = union(filter, [Builder.Girreducible,]), #ver4 evaluation only support one-particle-irreducible diagram
            transferLoop = KinL - KoutL,
            interaction = interaction,
            hasTau = hasTau,
            firstTauIdx = 1,
            totalTauNum = (loopNum + 1) * interactionTauNum
        )

        Parquet = Builder.ParquetNew

        F = [Parquet.U, Parquet.S]
        V = [Parquet.T, Parquet.U]

        varK = rand(Kdim, para.totalLoopNum)
        varT = [rand() for i in 1:para.totalTauNum]

        #################### DiagTree ####################################
        diag, nodes = Parquet.Ver4(para, legK, chan, F, V)
        # the weighttype of the returned ver4 is Float64
        dir = [n.node.index for n in nodes if n.DiEx == 1]
        ex = [n.node.index for n in nodes if n.DiEx == 2]
        println(nodes)

        rootDir = DiagTree.addnode!(diag, DiagTree.ADD, :dir, dir; para = [0, 0, 0, 0])
        rootEx = DiagTree.addnode!(diag, DiagTree.ADD, :ex, ex; para = [0, 0, 0, 0])
        diag.root = [rootDir.index, rootEx.index]

        # println(diag.root)
        DiagTree.showTree(diag, diag.root[1])
    end

    Parquet = Builder.ParquetNew
    testVer4(0, [Parquet.T, Parquet.U, Parquet.S], [Builder.Interaction(Builder.ChargeCharge, [Builder.Instant,]),]; hasTau = true)

end