"""
    green(para::DiagPara, extK = DiagTree.getK(para.totalLoopNum, 1), extT = para.hasTau ? (1, 2) : (0, 0), subdiagram = false;
        name = :G, resetuid = false, blocks::ParquetBlocks=ParquetBlocks())
    
    Build composite Green's function.
    By definition, para.firstTauIdx is the first Tau index of the left most self-energy subdiagram.

# Arguments
- `para`            : parameters. It should provide internalLoopNum, interactionTauNum, firstTauIdx
- `extK`            : basis of external loop. 
- `extT`: [Tau index of the left leg, Tau index of the right leg]
- `subdiagram`      : a sub-vertex or not
- `name`            : name of the diagram
- `resetuid`        : restart uid count from 1
- `blocks`          : building blocks of the Parquet equation. See the struct ParquetBlocks for more details.


# Output
- A Diagram object or nothing if the Green's function is illegal. 
"""
function green(para::DiagPara{W}, extK=DiagTree.getK(para.totalLoopNum, 1), extT=para.hasTau ? (1, 2) : (0, 0), subdiagram=false;
    name=:G, resetuid=false, blocks::ParquetBlocks=ParquetBlocks()) where {W}

    @assert isValidG(para) "$para doesn't gives a valid Green's function"
    @assert para.type == GreenDiag
    @assert para.innerLoopNum >= 0
    # @assert length(extK) == para.totalLoopNum
    @assert length(extT) == 2

    @assert length(extK) >= para.totalLoopNum "expect dim of extK>=$(para.totalLoopNum), got $(length(extK))"
    extK = extK[1:para.totalLoopNum]

    resetuid && uidreset()

    tin, tout = extT[1], extT[2]
    t0 = para.firstTauIdx

    # if isValidG(para) == false
    #     return nothing
    # end

    if para.innerLoopNum == 0
        return Diagram{W}(BareGreenId(para, k=extK, t=extT), name=name)
    end

    # ################# after this step, the Green's function must be nontrivial! ##################

    # if (para.extra isa ParquetBlocks) == false
    #     parquetblocks = ParquetBlocks(phi=[PPr, PHEr], ppi=[PHr, PHEr], Γ4=[PPr, PHr, PHEr])
    #     para::DiagPara = reconstruct(para, extra=parquetblocks)
    # end

    function ΣG(group, oG, Tidx, Kidx, ΣTidx)
        #type: Instant or Dynamic
        paraG = reconstruct(para, type=GreenDiag, firstTauIdx=Tidx, firstLoopIdx=Kidx, innerLoopNum=oG)
        G = Parquet.green(paraG, extK, group[:GT], true; blocks=blocks)
        # G = Diagram(GreenId(paraG, k = extK, t = group[:GT]), name = Symbol("g#$li")) #there is only one G diagram for a extT
        @assert G isa Diagram
        # println(group)
        pairT = (t=(ΣTidx, (group[:GT][2])),)
        return Diagram{W}(GenericId(para, pairT), Prod(), [group[:diagram], G], name=:ΣG)
    end

    para0 = reconstruct(para, innerLoopNum=0) #parameter for g0
    g0 = Diagram{W}(BareGreenId(para0, k=extK, t=(tin, t0)), name=:g0)
    ΣGpairs = Vector{Diagram{W}}()
    for p in orderedPartition(para.innerLoopNum, 2, 0)
        oΣ, oG = p

        if (isValidSigma(para.filter, oΣ, true) == false) || (isValidG(para.filter, oG) == false)
            continue
        end

        idx, maxTau = findFirstTauIdx(p, [SigmaDiag, GreenDiag], t0, interactionTauNum(para))
        @assert maxTau <= para.totalTauNum "maxTau = $maxTau > $(para.totalTauNum)"
        if para.hasTau
            @assert tin < t0 || tin > maxTau "external T index cann't be with in [$t0, $maxTau]"
            @assert tout < t0 || tout > maxTau "external T index cann't be with in [$t0, $maxTau]"
        end
        ΣfirstTidx, GfirstTidx = idx

        idx, maxLoop = findFirstLoopIdx(p, para.firstLoopIdx)
        @assert maxLoop <= para.totalLoopNum "maxLoop = $maxLoop > $(para.totalLoopNum)"
        ΣfirstKidx, GfirstKidx = idx

        sigmaPara = reconstruct(para, type=SigmaDiag, firstTauIdx=ΣfirstTidx, firstLoopIdx=ΣfirstKidx, innerLoopNum=oΣ)
        # println(ΣfirstTidx)
        sigma = Parquet.sigma(sigmaPara, extK, true, name=:Σ, blocks=blocks)
        @assert all(x -> x[1] == ΣfirstTidx, sigma.extT) "all sigma should share the same in Tidx\n$sigma"

        #combine sigmas with the same out Tidx
        df = transform(sigma, :extT => ByRow(x -> [x[1], (x[2], extT[2]),]) => [:Tin, :GT,])
        allsame(df, :Tin)
        groups = mergeby(df, :GT, operator=Sum())

        #for each sigma group, attach a G, all pair ΣG share the same in and out Tidx
        append!(ΣGpairs, [ΣG(g, oG, GfirstTidx, GfirstKidx, ΣfirstTidx) for g in eachrow(groups)])
    end
    # println(ΣGpairs)
    # println(operator)
    ΣGmerged = mergeby(ΣGpairs; operator=Sum(), name=:gΣG)[1]
    compositeG = Diagram{W}(GreenId(para, k=extK, t=extT), Prod(), [g0, ΣGmerged], name=:G)
    return compositeG
end