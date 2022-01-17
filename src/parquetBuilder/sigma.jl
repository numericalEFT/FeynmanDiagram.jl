"""
    function buildSigma(para, externLoop; F = [I, U, S], V = [I, T, U], All = union(F, V), diag = newDiagTree(para, :Sigma), subdiagram = false)
    
    build sigma diagram. 
    When sigma is created as a subdiagram, then no Fock diagram is generated if para.filter contains NoFock, and no sigma diagram is generated if para.filter contains Girreducible

"""
function buildSigma(para, externLoop, subdiagram = false;
    F = [I, U, S], V = [I, T, U], All = union(F, V),
    diag = newDiagTree(para, :Sigma))
    @assert para.innerLoopNum >= 1
    @assert length(externLoop) == para.totalLoopNum
    tright = para.firstTauIdx - 1 + para.innerLoopNum * para.interactionTauNum
    @assert para.totalTauNum >= tright "totalTauNum = $(para.totalTauNum) is not enough, sigma requires $tright\npara=$para"
    @assert para.totalLoopNum >= para.firstLoopIdx -1 + para.innerLoopNum

    if isValidSigma(para.filter, para.innerLoopNum, subdiagram) == false
        return diag, Vector{Component}([]), Vector{Component}([])
    end

    K = zero(externLoop)
    K[para.firstLoopIdx] = 1.0
    t0 = para.firstTauIdx

    function collapse!(dict, diag, nodes, factor)
        for n in nodes
            node = diag.nodePool.object[n.index]
            extT = node.para
            sigmaT = (extT[INL], extT[OUTR])
            if haskey(dict, sigmaT) == false
                dict[sigmaT] = []
            end
            push!(dict[sigmaT], (n, extT, factor))
        end
    end

    instant, dynamic = [], []

    ############### Fock-type diagram #######################################
    #  /=== W ===\
    # /           \
    # ----- G -----
    ##########################################################################
    qe = K - externLoop
    factor = 1 / (2π)^para.loopDim

    # if para.interactionTauNum == 1 || para.interactionTauNum == 0
    #     paraG = reconstruct(para, innerLoopNum = para.innerLoopNum - 1, firstLoopIdx = para.firstLoopIdx + 1, firstTauIdx = t0 + 1)
    #     if isValidG(paraG)
    #         diag, g = buildG(paraG, K, [t0, t0]; diag = diag)
    #         @assert g.index != 0
    #         # v = DiagTree.addpropagator!(diag, :Vpool, 1, :V; loop = qe)
    #         push!(instant, DiagTree.addnode!(diag, MUL, :fockΣ, [g, v], factor; para = [t0, t0]))
    #     end
    # elseif para.interactionTauNum == 2
    #     paraG = reconstruct(para, innerLoopNum = para.innerLoopNum - 1, firstLoopIdx = para.firstLoopIdx + 1, firstTauIdx = t0 + 2)
    #     if isValidG(paraG)
    #         diag, gv = buildG(paraG, K, [t0, t0]; diag = diag)
    #         @assert gv.index != 0
    #         v = DiagTree.addpropagator!(diag, :Vpool, 1, :V; loop = qe, site = (t0, t0 + 1))
    #         push!(instant, DiagTree.addnode!(diag, MUL, :fockΣ, [gv, v], factor; para = (t0, t0)))

    #         diag, gw = buildG(paraG, K, [t0, t0 + 1]; diag = diag)
    #         @assert gw.index != 0
    #         w = DiagTree.addpropagator!(diag, :Wpool, 1, :W; loop = qe, site = (t0, t0 + 1))
    #         push!(dynamic, DiagTree.addnode!(diag, MUL, :fockΣ, [gw, w], factor; para = [t0, t0 + 1]))
    #     end
    # else
    #     error("not implemented!")
    # end

    #if interaction is dynamic, then first two tau variables are reversed for the in and out vertices
    paraG = reconstruct(para, diagType = GreenDiag, innerLoopNum = para.innerLoopNum - 1,
        firstLoopIdx = para.firstLoopIdx + 1, firstTauIdx = t0 + para.interactionTauNum)
    paraW = reconstruct(para, diagType = Ver4Diag, innerLoopNum = 0, firstTauIdx = t0)
    #TODO: add validation for paraW
    if isValidG(paraG)
        legK = [externLoop, K, K, externLoop]
        ver4 = toDataFrame(bareVer4!(diag, paraW, legK, [DI, EX]))
        println(ver4)
        println(groupby(ver4, :response))
    end

    # if para.innerLoopNum >= 2
    #     ################# Sigma beyond the Fock diagram #################################
    #     #all self-energy diagram will be countered twice, thus a factor 1/2 is needed.
    #     factor = 1 / (2π)^para.loopDim * (1 / 2)
    #     KinL, KoutR = externLoop, externLoop
    #     KinR, KoutL = K, K
    #     for ver4LoopNum in 1:para.innerLoopNum-1
    #         gLoopNum = para.innerLoopNum - 1 - ver4LoopNum
    #         # if G doesn't exist, continue without creating ver4 node in the diagram
    #         if isValidG(para.filter, gLoopNum) == false
    #             continue
    #         end

    #         ver4Para = reconstruct(para, firstLoopIdx = para.firstLoopIdx + 1, innerLoopNum = ver4LoopNum)
    #         diag, ver4, dir, ex = buildVer4(ver4Para, [KinL, KoutL, KinR, KoutR],
    #             [T,], F, V, All; Fouter = [], Allouter = All, diag = diag)

    #         dict = Dict{Tuple{Int,Int},Vector{Any}}()
    #         collapse!(dict, diag, dir, 1.0)
    #         collapse!(dict, diag, ex, para.spin)
    #         #exchange ver4 has an additional Fermi loop compared to the direct counterpart when plugged into sigma

    #         for key in keys(dict)
    #             nodes = []
    #             for (n, extT, spinFactor) in dict[key]

    #                 tpair = [extT[INL], extT[OUTR]]
    #                 paraG = reconstruct(para,
    #                     innerLoopNum = gLoopNum,
    #                     firstLoopIdx = para.firstLoopIdx + 1 + ver4LoopNum,
    #                     firstTauIdx = maxTauIdx(ver4) + 1)
    #                 diag, g = buildG(paraG, K, tpair; diag = diag)
    #                 @assert g.index != 0
    #                 push!(nodes, DiagTree.addnode!(diag, MUL, :Σ, [n, g], factor * spinFactor; para = tpair))
    #             end

    #             push!(dynamic, DiagTree.addnode!(diag, DiagTree.ADD, :Σ, nodes; para = collect(key))) #key = (extT[INL], extT[OUTR])
    #             @assert dynamic[end].index != 0
    #         end
    #     end
    # end

    # make sure all incoming Tau idx is equal
    # for nidx in dynamic
    #     node = DiagTree.getNode(diag, nidx)
    #     @assert node.para[1] == para.firstTauIdx
    #     @assert node.para[2] <= tright
    # end
    # @assert isempty(dynamic) == false || isempty(instant) == false "Sigma diagram doesn't exist for \n$para"
    return diag, instant, dynamic

end