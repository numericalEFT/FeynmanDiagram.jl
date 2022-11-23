"""
    function sigma(para, extK = DiagTree.getK(para.totalLoopNum, 1), subdiagram = false; name = :Σ, resetuid = false, blocks::ParquetBlocks=ParquetBlocks())
    
    Build sigma diagram. 
    When sigma is created as a subdiagram, then no Fock diagram is generated if para.filter contains NoFock, and no sigma diagram is generated if para.filter contains Girreducible

# Arguments
- `para`            : parameters. It should provide internalLoopNum, interactionTauNum, firstTauIdx
- `extK`            : basis of external loop. 
- `subdiagram`      : a sub-vertex or not
- `name`            : name of the diagram
- `resetuid`        : restart uid count from 1
- `blocks`          : building blocks of the Parquet equation. See the struct ParquetBlocks for more details.

# Output
- A DataFrame with fields `:type`, `:extT`, `:diagram`, `:hash`
- All sigma share the same incoming Tau index, but not the outgoing one
"""
function sigma(para::DiagPara{W}, extK=DiagTree.getK(para.totalLoopNum, 1), subdiagram=false; name=:Σ, resetuid=false, blocks::ParquetBlocks=ParquetBlocks()) where {W}
    resetuid && uidreset()
    (para.type == SigmaDiag) || error("$para is not for a sigma diagram")
    (para.innerLoopNum >= 1) || error("sigma must has more than one inner loop")
    # @assert length(extK) == para.totalLoopNum
    # @assert (para.innerLoopNum <= 1) || ((NoBubble in para.filter) == false) "The many-body correction sigma only accounts for half of the bubble counterterm right now."
    if (para.innerLoopNum > 1) && (NoBubble in para.filter)
        @warn "Sigma with two or more loop orders still contain bubble subdiagram even if NoBubble is turned on in para.filter!"
    end

    (length(extK) >= para.totalLoopNum) || error("expect dim of extK>=$(para.totalLoopNum), got $(length(extK))")
    extK = extK[1:para.totalLoopNum]

    compositeSigma = DataFrame(type=AnalyticProperty[], extT=Tuple{Int,Int}[], diagram=Diagram{W}[])

    if isValidSigma(para.filter, para.innerLoopNum, subdiagram) == false
        # return DataFrame(type=[], extT=[], diagram=[])
        return compositeSigma
    end

    # if (para.extra isa ParquetBlocks) == false
    #     parquetblocks = ParquetBlocks(phi=[PPr, PHEr], ppi=[PHr, PHEr], Γ4=[PPr, PHr, PHEr])
    #     para::DiagPara = reconstruct(para, extra=parquetblocks)
    # end

    K = zero(extK)
    LoopIdx = para.firstLoopIdx
    K[LoopIdx] = 1.0
    (isapprox(K, extK) == false) || error("K and extK can not be the same")
    legK = [extK, K, K, extK]

    function GWwithGivenExTtoΣ(group, oW, paraG)
        # println(group)
        # @assert length(group[:, :diagram]) == 1
        # allsame(group, [:response, :type, :GT])
        (group[:response] == UpUp || group[:response] == UpDown) || error("GW with given ExT to Σ only works for UpUp or UpDown")
        #type: Instant or Dynamic
        response, type = group[:response], group[:type]
        sid = SigmaId(para, type, k=extK, t=group[:extT])
        g = green(paraG, K, group[:GT], true; name=(oW == 0 ? :Gfock : :G_Σ), blocks=blocks) #there is only one G diagram for a extT
        (g isa Diagram) || error("green function must return a Diagram")
        # Sigma = G*(2 W↑↑ - W↑↓)
        # ! The sign of ↑↓ is from the spin symmetry, not from the fermionic statistics!
        spinfactor = (response == UpUp) ? 2 : -1
        # spinfactor = (response == UpUp) ? 0 : 1
        if oW > 0 # oW are composte Sigma, there is a symmetry factor 1/2
            spinfactor *= 0.5
            # elseif oW == 0 # the Fock diagram requires an additional minus sign, because the interaction is currently a direct one, but is expected to be exchange.
            #     spinfactor *= paraG.isFermi ? -1.0 : 1.0
        end
        # plot_tree(mergeby(DataFrame(group)), maxdepth = 7)
        sigmadiag = Diagram{W}(sid, Prod(), [g, group[:diagram]], factor=spinfactor, name=name)
        # plot_tree(sigmadiag, maxdepth = 7)
        return (type=type, extT=group[:extT], diagram=sigmadiag)
    end

    for (oG, oW) in orderedPartition(para.innerLoopNum - 1, 2, 0)

        idx, maxLoop = findFirstLoopIdx([oW, oG], LoopIdx + 1)
        (maxLoop <= para.totalLoopNum) || error("maxLoop = $maxLoop > $(para.totalLoopNum)")
        WfirstLoopIdx, GfirstLoopIdx = idx

        # it is important to do W first, because the left in of W is also the incoming leg of sigma, they have the same Tidx
        idx, maxTau = findFirstTauIdx([oW, oG], [Ver4Diag, GreenDiag], para.firstTauIdx, interactionTauNum(para))
        (maxTau <= para.totalTauNum) || error("maxTau = $maxTau > $(para.totalTauNum)")
        WfirstTauIdx, GfirstTauIdx = idx

        paraG = reconstruct(para, type=GreenDiag, innerLoopNum=oG,
            firstLoopIdx=GfirstLoopIdx, firstTauIdx=GfirstTauIdx)
        paraW = reconstruct(para, type=Ver4Diag, innerLoopNum=oW,
            firstLoopIdx=WfirstLoopIdx, firstTauIdx=WfirstTauIdx)

        #TODO: add validation for paraW
        # println("oG: $oG, oW: $oW  -> ", isValidG(paraG))
        if isValidG(paraG)
            # println(paraW)
            if oW == 0 # Fock-type Σ
                if NoHartree in paraW.filter
                    paraW0 = reconstruct(paraW, filter=union(paraW.filter, Proper), transferLoop=zero(K))
                    ver4 = vertex4(paraW0, legK, [], true)
                else
                    ver4 = vertex4(paraW, legK, [], true)
                end
                # println(ver4)
            else # composite Σ
                # paraW0 = reconstruct(paraW, filter=union(paraW.filter, Proper), transferLoop=extK-K)
                ver4 = vertex4(paraW, legK, [PHr,], true; blocks=blocks, blockstoplevel=ParquetBlocks(phi=[], Γ4=[PHr, PHEr, PPr]))
                # plot_tree(mergeby(ver4).diagram[1])
            end
            #transform extT coloum intwo extT for Σ and extT for G
            # plot_tree(ver4)
            df = transform(ver4, :extT => ByRow(x -> [(x[INL], x[OUTR]), (x[OUTL], x[INR])]) => [:extT, :GT])

            groups = mergeby(W, df, [:response, :type, :GT, :extT], operator=Sum())
            for mergedVer4 in eachrow(groups)
                push!(compositeSigma, GWwithGivenExTtoΣ(mergedVer4, oW, paraG))
            end
        end
    end


    if isempty(compositeSigma)
        return compositeSigma
    else
        # Factor = 1 / (2π)^para.loopDim
        Factor = 1.0
        # sigmaid(g) = SigmaId(para, g[1, :type], k=extK, t=g[1, :extT])
        # sigmadf = mergeby(compositeSigma, [:type, :extT], name=name, factor=Factor, getid=sigmaid)
        sigmadf = mergeby(W, compositeSigma, [:type, :extT], name=name, factor=Factor,
            getid=g -> SigmaId(para, g[1, :type], k=extK, t=g[1, :extT]))

        all(x -> x[1] == para.firstTauIdx, sigmadf.extT) || error("all sigma should share the same in Tidx\n$sigmadf")
        # println(sigmadf)
        # plot_tree(sigmadf)
        return sigmadf
    end
end