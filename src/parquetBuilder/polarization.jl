
function buildPolarization(para, extK, subdiagram = false; name = :Π)

    (subdiagram == false) && uidreset()
    @assert para.diagType == PolarDiag
    @assert para.innerLoopNum >= 1
    @assert length(extK) == para.totalLoopNum
    tright = para.firstTauIdx - 1 + para.innerLoopNum * para.interactionTauNum
    @assert para.totalTauNum >= tright "totalTauNum = $(para.totalTauNum) is not enough, sigma requires $tright\npara=$para"
    @assert para.totalLoopNum >= para.firstLoopIdx -1 + para.innerLoopNum

    if (para.extra isa ParquetBlocks) == false
        parquetblocks = ParquetBlocks(phi = [PPr, PHEr], ppi = [PHr, PHEr], Γ4 = [PPr, PHr, PHEr])
        para = reconstruct(para, extra = parquetblocks)
    end

    K = zero(extK)
    LoopIdx = para.firstLoopIdx
    K[LoopIdx] = 1.0
    @assert (K ≈ extK) == false
    t0 = para.firstTauIdx
    extT = (t0, t0 + 1)

    polar = DataFrame()

    ################# Π0 = GG #######################
    for (oG1, oG2) in orderedPartition(para.innerLoopNum - 1, 2, 0)

        idx, maxLoop = findFirstLoopIdx([oG1, oG2], LoopIdx + 1)
        @assert maxLoop <= para.totalLoopNum "maxLoop = $maxLoop > $(para.totalLoopNum)"
        G1firstLoopIdx, G2firstLoopIdx = idx

        idx, maxTau = findFirstTauIdx([oG1, oG2], [GreenDiag, GreenDiag], para.firstTauIdx, para.interactionTauNum)
        @assert maxTau <= para.totalTauNum
        G1firstTauIdx, G2firstTauIdx = idx

        paraG1 = reconstruct(para, diagType = GreenDiag, innerLoopNum = oG1,
            firstLoopIdx = G1firstLoopIdx, firstTauIdx = G1firstTauIdx)
        paraG2 = reconstruct(para, diagType = GreenDiag, innerLoopNum = oG2,
            firstLoopIdx = G2firstLoopIdx, firstTauIdx = G2firstTauIdx)

        if isValidG(paraG1) && isValidG(paraG2)
            g1 = buildG(paraG1, extK .+ K, (t0, t0 + 1), true, name = :Gp)
            g2 = buildG(paraG2, K, (t0 + 1, t0), true, name = :Gh)
            Π0uu = Diagram(PolarId(para, UpUp, k = extK, t = extT), Prod(), [g1, g2]; name = :Π0)
            push!(polar, (response = UpUp, diagram = Π0uu))
        end
    end

    # legK = [extK, K, K, extK]
    polar = mergeby(polar, [:response,]; name = name,
        getid = g -> PolarId(para, g[1, :response], k = extK, t = extT)
    )
    return polar
end