"""
    function buildSigma(para, extK; F = [I, U, S], V = [I, T, U], All = union(F, V), diag = newDiagTree(para, :Sigma), subdiagram = false)
    
    build sigma diagram. 
    When sigma is created as a subdiagram, then no Fock diagram is generated if para.filter contains NoFock, and no sigma diagram is generated if para.filter contains Girreducible

"""
function buildSigma(para, extK, subdiagram = false; name = :none)
    subdiagram == false && uidreset()
    @assert para.innerLoopNum >= 1
    @assert length(extK) == para.totalLoopNum
    tright = para.firstTauIdx - 1 + para.innerLoopNum * para.interactionTauNum
    @assert para.totalTauNum >= tright "totalTauNum = $(para.totalTauNum) is not enough, sigma requires $tright\npara=$para"
    @assert para.totalLoopNum >= para.firstLoopIdx -1 + para.innerLoopNum

    diags = Diagram{para.weightType}[]

    if isValidSigma(para.filter, para.innerLoopNum, subdiagram) == false
        return diags
    end

    K = zero(extK)
    LoopIdx = para.firstLoopIdx
    K[LoopIdx] = 1.0
    t0 = para.firstTauIdx
    factor = 1 / (2π)^para.loopDim
    qe = K - extK
    legK = [extK, K, K, extK]


    function toSigma!(ver4)
        df = toDataFrame(ver4, verbose = 0)
        @assert all(d -> typeof(d.id) == Ver4Id, df.diagram) == false
        @assert all(x -> x == UpUp || x == UpDown, df[:, :response])



        for df in groupby(df, [:response, :type, :TinL, :ToutR])
            response, type = df[:response], df[:type]
            #type: Instant or Dynamic
            sid = SigmaId(para, type, k = extK, t = (df[:TinL], df[:ToutR]))
            for g in buildG(paraG, K, (df[:ToutL], df[:TinR]); name = :Gfock)
                # Sigma = G*(2 W↑↑ - W↑↓)
                spinfactor = (response == UpUp) ? 2 : -1
                push!(diags, Diagram(sid, Prod(), [g, df[:Diagram]], factor = spinfactor, name = name))
            end
        end
    end

    for (oG, oW) in orderedPartition(para.innerLoopNum, 2, 0)

        idx, maxLoop = findFirstLoopIdx([oG, oW], LoopIdx + 1)
        @assert maxLoop <= para.totalLoopNum
        GfirstLoopIdx, WfirstLoopIdx = idx

        idx, maxTau = findFirstTauIdx([oG, oW], [GreenDiag, Ver4Diag], para.firstTauIdx, para.interactionTauNum)
        @assert maxTau <= para.totalTauNum
        GfirstTauIdx, WfirstTauIdx = idx

        paraG = reconstruct(para, diagType = GreenDiag, innerLoopNum = oG,
            firstLoopIdx = GfirstLoopIdx, firstTauIdx = GfirstTauIdx)
        paraW = reconstruct(para, diagType = Ver4Diag, innerLoopNum = oW,
            firstLoopIdx = WfirstLoopIdx, firstTauIdx = WfirstTauIdx)

        #TODO: add validation for paraW
        if isValidG(paraG)
            if oW == 0 # Fock-type Σ
                ver4 = bareVer4(paraW, legK, [Ex,])
            else # composite Σ
                ver4 = buildVer4(paraW, [PHr,], true, phi_toplevel = [], Γ4_toplevel = paraW.extra.Γ4)
            end
            toSigma!(ver4)
        end
    end

    return diags

end