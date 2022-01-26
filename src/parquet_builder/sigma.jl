"""
    function sigma(para, extK = DiagTree.getK(para.totalLoopNum, 1), subdiagram = false; name = :Σ)
    
    Build sigma diagram. 
    When sigma is created as a subdiagram, then no Fock diagram is generated if para.filter contains NoFock, and no sigma diagram is generated if para.filter contains Girreducible

# Arguments
- `para`            : parameters. It should provide internalLoopNum, interactionTauNum, firstTauIdx
- `extK`            : basis of external loop. 
- `subdiagram`      : a sub-vertex or not
- `name`            : name of the diagram

# Output
- A DataFrame with fields `:type`, `:extT`, `:diagram`, `:hash`
- All sigma share the same incoming Tau index, but not the outgoing one
"""
function sigma(para, extK = DiagTree.getK(para.totalLoopNum, 1), subdiagram = false; name = :Σ)
    (subdiagram == false) && uidreset()
    @assert para.diagType == SigmaDiag
    @assert para.innerLoopNum >= 1
    @assert length(extK) == para.totalLoopNum

    if isValidSigma(para.filter, para.innerLoopNum, subdiagram) == false
        return DataFrame(type = [], extT = [], diagram = [])
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
        g = green(paraG, K, group[:GT], true; name = oW == 0 ? :Gfock : :G_Σ) #there is only one G diagram for a extT
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

        idx, maxLoop = findFirstLoopIdx([oW, oG], LoopIdx + 1)
        @assert maxLoop <= para.totalLoopNum "maxLoop = $maxLoop > $(para.totalLoopNum)"
        WfirstLoopIdx, GfirstLoopIdx = idx

        # it is important to do W first, because the left in of W is also the incoming leg of sigma, they have the same Tidx
        idx, maxTau = findFirstTauIdx([oW, oG], [Ver4Diag, GreenDiag], para.firstTauIdx, para.interactionTauNum)
        @assert maxTau <= para.totalTauNum
        WfirstTauIdx, GfirstTauIdx = idx

        paraG = reconstruct(para, diagType = GreenDiag, innerLoopNum = oG,
            firstLoopIdx = GfirstLoopIdx, firstTauIdx = GfirstTauIdx)
        paraW = reconstruct(para, diagType = Ver4Diag, innerLoopNum = oW,
            firstLoopIdx = WfirstLoopIdx, firstTauIdx = WfirstTauIdx)

        #TODO: add validation for paraW
        if isValidG(paraG)
            if oW == 0 # Fock-type Σ
                paraW0 = reconstruct(paraW, filter = union(paraW.filter, Proper), transferLoop = zero(K))
                ver4 = vertex4(paraW0, legK, [], true)
            else # composite Σ
                ver4 = vertex4(paraW, legK, [PHr,], true, phi_toplevel = [], Γ4_toplevel = [PHr, PHEr, PPr,])
                # plot_tree(mergeby(ver4).diagram[1])
            end
            #transform extT coloum intwo extT for Σ and extT for G
            # plot_tree(ver4)
            df = transform(ver4, :extT => ByRow(x -> [(x[INL], x[OUTR]), (x[OUTL], x[INR])]) => [:extT, :GT])

            groups = mergeby(df, [:response, :type, :GT, :extT], operator = Sum())
            for mergedVer4 in eachrow(groups)
                push!(compositeSigma, GWwithGivenExTtoΣ(mergedVer4, oW, paraG))
            end
        end
    end

    if isempty(compositeSigma)
        return DataFrame(type = [], extT = [], diagram = [])
    end

    Factor = 1 / (2π)^para.loopDim
    sigmadf = mergeby(compositeSigma, [:type, :extT], name = name, factor = Factor,
        getid = g -> SigmaId(para, g[1, :type], k = extK, t = g[1, :extT]))

    @assert all(x -> x[1] == para.firstTauIdx, sigmadf.extT) "all sigma should share the same in Tidx\n$sigmadf"
    # println(sigmadf)
    # plot_tree(sigmadf)
    return sigmadf
end