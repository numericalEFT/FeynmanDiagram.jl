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

    if isValidSigma(para.filter, para.innerLoopNum, subdiagram) == false
        return nothing
    end

    K = zero(extK)
    LoopIdx = para.firstLoopIdx
    K[LoopIdx] = 1.0
    t0 = para.firstTauIdx
    factor = 1 / (2π)^para.loopDim
    qe = K - extK
    legK = [extK, K, K, extK]

    function GWwithGivenExTtoΣ(group)
        allsame(group, [:response, :type, :GT])
        #type: Instant or Dynamic
        response, type = group[1, :response], group[1, :type]
        sid = SigmaId(para, type, k = extK, t = group[1, :extT])
        g = buildG(paraG, K, group[1, :GT]; name = :Gfock) #there is only one G diagram for a extT
        @assert g isa Diagram
        # Sigma = G*(2 W↑↑ - W↑↓)
        spinfactor = (response == UpUp) ? 2 : -1
        return Diagram(sid, Prod(), [g, group[1, :diagram]], factor = spinfactor, name = name)
    end

    compositeSigma = DataFrame()

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

            df = toDataFrame(ver4, expand = true)
            allsame(group, :id)
            @assert all(x -> x == UpUp || x == UpDown, df[:, :response])
            #transform extT coloum intwo extT for Σ and extT for G
            df = transform(df, :extT => ByRow(x -> [(x[INL], x[OUTR]), (x[OUTL], x[INR])]) => [:extT, :GT])
            println(df)
            for group in groupby(df, [:response, :type, :GT])
                newsigma = (type = group[1, :type], extT = group[1, :extT], diagram = GWwithGivenExTtoΣ(group))
                push!(compositeSigma, newsigma)
            end
        end
    end

    return mergeby(compositeSigma, [:type, :extT], name = name,
        getid = g -> SigmaId(para, g[1, :type], k = extK, t = g[1, :extT]))

end