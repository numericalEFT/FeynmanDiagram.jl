function buildVer3(para, legK, subdiagram = false; F = [I, U, S], V = [I, T, U], All = union(F, V), chan = All, diag = newDiagTree(para, :Ver3))
    KinR, KoutR = legK[1], legK[2]
    q = KinR - KoutR
    K = zero(KinR)
    loopIdx = para.firstLoopIdx
    tleft = para.firstTauIdx
    K[loopIdx] = 1.0
    tauNum = para.interactionTauNum
    maxLoopIdx = para.firstLoopIdx - 1 + para.innerLoopNum
    maxTauIdx = para.firstTauIdx - 1 + (1 + para.innerLoopNum * para.interactionTauNum)

    partition = orderedPartition(para.innerLoopNum - 1, 3, 0)
    for p in partition
        g1L, verL, g2L = p

        # idx, maxidx = findFirstLoopIdx(p, [true, false, true], loopIdx + 1)
        # g1firstLoopIdx, verfirstLoopIdx, g2firstLoopIdx = idx
        # @assert maxidx == maxLoopIdx

        # idx, maxidx = findFirstTauIdx(p, [true, false, true], tleft + 1, tauNum)
        # g1firstTauIdx, verfirstTauIdx, g2firstTauIdx = idx
        # @assert maxidx == maxTauIdx

        # g1Para = reconstruct(para, innerLoopNum = g1L, firstLoopIdx = g1firstLoopIdx, firstTauIdx = g1firstTauIdx)
        # g2Para = reconstruct(para, innerLoopNum = g2L, firstLoopIdx = g2firstLoopIdx, firstTauIdx = g2firstTauIdx)
        # verPara = reconstruct(para, innerLoopNum = verL, findFirstLoopIdx = verfirstLoopIdx, firstTauIdx = verfirstTauIdx)

        # if isValidG(g1Para) == false || isValidG(g2Para) == false
        #     continue
        # end

        # Kx = K + KinR - KoutR

        # g1, isnode1 = buildG(g1Para, K, [tleft, tleft + 1]; diag = diag)
        # ver4 = buildVer4(verPara, [K, Kx, KinR, KoutR], chan, F, V, All; diag = diag)

        # g1, isnode1 = buildG(g1Para, K, [tleft, tleft]; diag = diag)
    end
end