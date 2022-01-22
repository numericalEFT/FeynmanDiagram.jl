"""
    function buildSigma(para, extK; F = [I, U, S], V = [I, T, U], All = union(F, V), diag = newDiagTree(para, :Sigma), subdiagram = false)
    
    build sigma diagram. 
    When sigma is created as a subdiagram, then no Fock diagram is generated if para.filter contains NoFock, and no sigma diagram is generated if para.filter contains Girreducible

"""
function buildSigma(para, extK, subdiagram = false; name = :Σ)
    (subdiagram == false) && uidreset()
    @assert para.diagType == SigmaDiag
    @assert para.innerLoopNum >= 1
    @assert length(extK) == para.totalLoopNum
    tright = para.firstTauIdx - 1 + para.innerLoopNum * para.interactionTauNum
    @assert para.totalTauNum >= tright "totalTauNum = $(para.totalTauNum) is not enough, sigma requires $tright\npara=$para"
    @assert para.totalLoopNum >= para.firstLoopIdx -1 + para.innerLoopNum

    if isValidSigma(para.filter, para.innerLoopNum, subdiagram) == false
        return nothing
    end

    if (para.extra isa ParquetBlocks) == false
        parquetblocks = ParquetBlocks(phi = [PPr, PHEr], ppi = [PHr, PHEr], Γ4 = [PPr, PHr, PHEr])
        para = reconstruct(para, extra = parquetblocks)
    end

    K = zero(extK)
    LoopIdx = para.firstLoopIdx
    K[LoopIdx] = 1.0
    @assert (K ≈ extK) == false
    legK = [extK, K, K, extK]

    function GWwithGivenExTtoΣ(group, oW, paraG)
        # println(group)
        # @assert length(group[:, :diagram]) == 1
        # allsame(group, [:response, :type, :GT])
        @assert group[:response] == UpUp || group[:response] == UpDown
        #type: Instant or Dynamic
        response, type = group[:response], group[:type]
        sid = SigmaId(para, type, k = extK, t = group[:extT])
        g = buildG(paraG, K, group[:GT]; name = oW == 0 ? :Gfock : :G_Σ) #there is only one G diagram for a extT
        @assert g isa Diagram
        # Sigma = G*(2 W↑↑ - W↑↓)
        # ! The sign of ↑↓ is from the spin symmetry, not from the fermionic statistics!
        spinfactor = (response == UpUp) ? 2 : -1
        # spinfactor = (response == UpUp) ? 0 : 1
        if oW > 0 # oW are composte Sigma, there is a symmetry factor 1/2
            spinfactor *= 0.5
        end
        # plot_tree(mergeby(DataFrame(group)), maxdepth = 7)
        sigmadiag = Diagram(sid, Prod(), [g, group[:diagram]], factor = spinfactor, name = name)
        # plot_tree(sigmadiag, maxdepth = 7)
        return (type = type, extT = group[:extT], diagram = sigmadiag)
    end

    compositeSigma = DataFrame()

    for (oG, oW) in orderedPartition(para.innerLoopNum - 1, 2, 0)

        idx, maxLoop = findFirstLoopIdx([oG, oW], LoopIdx + 1)
        @assert maxLoop <= para.totalLoopNum "maxLoop = $maxLoop > $(para.totalLoopNum)"
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
                paraW0 = reconstruct(paraW, filter = union(paraW.filter, Proper), transferLoop = zero(K))
                ver4 = buildVer4(paraW0, legK, [], true)
            else # composite Σ
                ver4 = buildVer4(paraW, legK, [PHr,], true, phi_toplevel = [], Γ4_toplevel = [PHr, PHEr, PPr,])
                # plot_tree(mergeby(ver4).diagram[1])
            end
            #transform extT coloum intwo extT for Σ and extT for G
            df = transform(ver4, :extT => ByRow(x -> [(x[INL], x[OUTR]), (x[OUTL], x[INR])]) => [:extT, :GT])

            groups = mergeby(df, [:response, :type, :GT, :extT], operator = Sum())
            for mergedVer4 in eachrow(groups)
                push!(compositeSigma, GWwithGivenExTtoΣ(mergedVer4, oW, paraG))
            end
        end
    end

    Factor = 1 / (2π)^para.loopDim
    sigmadf = mergeby(compositeSigma, [:type, :extT], name = name, factor = Factor,
        getid = g -> SigmaId(para, g[1, :type], k = extK, t = g[1, :extT]))
    return sigmadf
end