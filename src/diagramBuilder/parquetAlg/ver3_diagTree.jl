function buildVer3(para, legK, subdiagram = false; F = [I, U, S], V = [I, T, U], All = union(F, V), diag = newDiagTree(para, :Ver3))
    KinR, KoutR = legK[1], legK[2]
    q = k1 - k2
    K = zero(k1)
    K[para.firstLoopIdx] = 1.0
    tleft = para.firstTauIdx
    tauNum = para.interactionTauNum
    maxLoopIdx = para.firstLoopIdx - 1 + para.innerLoopNum
    maxTauIdx = para.firstTauIdx - 1 + (1 + para.innerLoopNum * para.interactionTauNum)

    partition = orderedPartition(para.innerLoopNum - 1, 3, 0)
    for (g1L, verL, g2L) in partition
        # # example 1: (1, 2, 1), 
        # # example 2: (0, 2, 0)
        # G1firstLoopIdx = para.firstLoopIdx + 1 #first inner loop index of G1 (may not exist if G1 loop number is zero) 
        # verfirstLoopIdx = para.firstLoopIdx + g1L
        # G2firstLoopIdx = verfirstLoopIdx + verL #first inner loop index of G2 (may not exist if G2 loop number is zero)
        # loopRight = verfirstLoopIdx + verL + g2L - 1
        # # eg1: 4+2+1-1 == 6, eg2: 3+2+0-1 == 4
        # @assert loopRight == maxLoopIdx

        # G1firstTauIdx = tleft + 1
        # verfirstTauIdx = tleft + g1L * tauNum
        # G2firstTauIdx = verfirstTauIdx + (verL + 1) * TauNum
        # tauRight = G2firstTauIdx + oGx * TauNum - 1
        # # eg1: 4+3+1-1==7, eg2: 3+3+0-1==5
        # @assert tauRight == maxTauIdx(ver4)

    end
end