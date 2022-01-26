"""
    function buildG(para, extK, extT, subdiagram = false; F = [I, U, S], V = [I, T, U], All = union(F, V), diag = newDiagTree(para, :G))
    
    Build composite Green's function.
    By definition, para.firstTauIdx is the first Tau index of the left most self-energy subdiagram.

- `extT`: [Tau index of the left leg, Tau index of the right leg]

"""
function green(para, extK = DiagTree.getK(para.totalLoopNum, 1), extT = para.hasTau ? (1, 2) : (0, 0), subdiagram = false; name = :G)
    @assert para.diagType == GreenDiag
    @assert para.innerLoopNum >= 0
    @assert length(extK) == para.totalLoopNum
    @assert length(extT) == 2

    (subdiagram == false) && uidreset()

    tin, tout = extT[1], extT[2]
    t0 = para.firstTauIdx

    if isValidG(para) == false
        return nothing
    end

    if para.innerLoopNum == 0
        return Diagram(GreenId(para, k = extK, t = extT), name = name)
    end

    # ################# after this step, the Green's function must be nontrivial! ##################

    if (para.extra isa ParquetBlocks) == false
        parquetblocks = ParquetBlocks(phi = [PPr, PHEr], ppi = [PHr, PHEr], Γ4 = [PPr, PHr, PHEr])
        para = reconstruct(para, extra = parquetblocks)
    end

    function ΣG(group, oG, Tidx, Kidx, ΣTidx)
        #type: Instant or Dynamic
        paraG = reconstruct(para, diagType = GreenDiag, firstTauIdx = Tidx, firstLoopIdx = Kidx, innerLoopNum = oG)
        G = Parquet.green(paraG, extK, group[:GT], true)
        # G = Diagram(GreenId(paraG, k = extK, t = group[:GT]), name = Symbol("g#$li")) #there is only one G diagram for a extT
        @assert isnothing(G) == false
        # println(group)
        pairT = (t = (ΣTidx, (group[:GT][2])),)
        return Diagram(GenericId(para, pairT), Prod(), [group[:diagram], G], name = :ΣG)
    end

    para0 = reconstruct(para, innerLoopNum = 0) #parameter for g0
    g0 = Diagram(GreenId(para0, k = extK, t = (tin, t0)), name = :g0)
    ΣGpairs = []
    for p in orderedPartition(para.innerLoopNum, 2, 0)
        oΣ, oG = p

        if (isValidSigma(para.filter, oΣ, true) == false) || (isValidG(para.filter, oG) == false)
            continue
        end

        idx, maxTau = findFirstTauIdx(p, [SigmaDiag, GreenDiag], t0, para.interactionTauNum)
        @assert maxTau <= para.totalTauNum "maxTau = $maxTau > $(para.totalTauNum)"
        if para.hasTau
            @assert tin < t0 || tin > maxTau "external T index cann't be with in [$t0, $maxTau]"
            @assert tout < t0 || tout > maxTau "external T index cann't be with in [$t0, $maxTau]"
        end
        ΣfirstTidx, GfirstTidx = idx

        idx, maxLoop = findFirstLoopIdx(p, para.firstLoopIdx)
        @assert maxLoop <= para.totalLoopNum "maxLoop = $maxLoop > $(para.totalLoopNum)"
        ΣfirstKidx, GfirstKidx = idx

        sigmaPara = reconstruct(para, diagType = SigmaDiag, firstTauIdx = ΣfirstTidx, firstLoopIdx = ΣfirstKidx, innerLoopNum = oΣ)
        sigma = Parquet.sigma(sigmaPara, extK, true, name = :Σ)
        @assert all(x -> x[1] == ΣfirstTidx, sigma.extT) "all sigma should share the same in Tidx\n$sigma"

        #combine sigmas with the same out Tidx
        df = transform(sigma, :extT => ByRow(x -> [x[1], (x[2], extT[2]),]) => [:Tin, :GT,])
        allsame(df, :Tin)
        groups = mergeby(df, :GT, operator = Sum())

        #for each sigma group, attach a G, all pair ΣG share the same in and out Tidx
        append!(ΣGpairs, [ΣG(g, oG, GfirstTidx, GfirstKidx, ΣfirstTidx) for g in eachrow(groups)])
    end
    ΣGmerged = mergeby(ΣGpairs, operator = Sum(), name = :gΣG)
    @assert nrow(ΣGmerged) == 1
    return Diagram(GreenId(para, k = extK, t = extT), Prod(), [g0, ΣGmerged.diagram[1]], name = :G)
end

# function green(para, extK = DiagTree.getK(para.totalLoopNum, 1), extT = para.hasTau ? (1, 2) : (0, 0), subdiagram = false; name = :G)
#     @assert para.diagType == GreenDiag
#     @assert para.innerLoopNum >= 0
#     @assert length(extK) == para.totalLoopNum
#     @assert length(extT) == 2

#     (subdiagram == false) && uidreset()

#     tin, tout = extT[1], extT[2]
#     innerTauNum = para.innerLoopNum * para.interactionTauNum
#     tstart = para.firstTauIdx
#     tend = para.firstTauIdx + innerTauNum - 1

#     if isValidG(para) == false
#         return nothing
#     end

#     if para.innerLoopNum == 0
#         return Diagram(GreenId(para, k = extK, t = extT), name = name)
#     end

#     # ################# after this step, the Green's function must be nontrivial! ##################
#     if para.hasTau
#         @assert tin < tstart || tin > tend "external T index cann't be with in [$tstart, $tend]"
#         @assert tout < tstart || tout > tend "external T index cann't be with in [$tstart, $tend]"
#     end

#     if (para.extra isa ParquetBlocks) == false
#         parquetblocks = ParquetBlocks(phi = [PPr, PHEr], ppi = [PHr, PHEr], Γ4 = [PPr, PHr, PHEr])
#         para = reconstruct(para, extra = parquetblocks)
#     end

#     function ΣG(group, li, paraG, tleft)
#         #type: Instant or Dynamic
#         g = Diagram(GreenId(paraG, k = extK, t = group[:GT]), name = Symbol("g#$li")) #there is only one G diagram for a extT
#         # println(group)
#         pairT = (t = (tleft, (group[:GT][2])),)
#         return Diagram(GenericId(para, pairT), Prod(), [g, group[:diagram]], name = Symbol("Σg#$li"))
#     end

#     para0 = reconstruct(para, innerLoopNum = 0) #parameter for g0
#     gleft = Diagram(GreenId(para0, k = extK, t = (tin, tstart)), name = :gleft)

#     partition = []
#     for n = 1:para.innerLoopNum
#         #generate all possible ordered partition of the loops
#         append!(partition, orderedPartition(para.innerLoopNum, n))
#         #e.g., loopNum =5, n =2 ==> ordered = [[4, 1], [1, 4], [3, 2], [2, 3]]
#     end
#     # #e.g., loopNum =5 ==> partition = [[4, 1], [1, 4], [3, 2], [2, 3], [1, 1, 1, 1], ...]

#     Gall = []
#     for p in partition
#         #e.g., p = [4, 1]
#         firstTauIdx, maxTau = findFirstTauIdx(p, [SigmaDiag for i in p], para.firstTauIdx, para.interactionTauNum)
#         @assert maxTau <= para.totalTauNum "maxTau = $maxTau > $(para.totalTauNum)"
#         firstLoopIdx, maxLoop = findFirstLoopIdx(p, para.firstLoopIdx)
#         @assert maxLoop <= para.totalLoopNum "maxLoop = $maxLoop > $(para.totalLoopNum)"

#         if all([isValidSigma(para.filter, loop, true) for loop in p]) == false
#             continue
#         end

#         ΣGpair = []
#         #for each Σ, build ΣG, all of them share the same in and out Tidx
#         for (li, loop) in enumerate(p)
#             #tleft, tright are the in and out Tidx of the Green's function on the right hand side of the current sigma
#             tleft = firstTauIdx[li]
#             #if li == length(p), then the current sigma is the last one, tright must be extT[2] !
#             tright = li < length(p) ? firstTauIdx[li+1] : tout
#             #e.g., loop = 4 or 1
#             sigmaPara = reconstruct(para, diagType = SigmaDiag, firstTauIdx = tleft, firstLoopIdx = firstLoopIdx[li], innerLoopNum = loop)
#             sigma = Parquet.sigma(sigmaPara, extK, true, name = Symbol("Σ#$li"))
#             @assert all(x -> x[1] == tleft, sigma.extT) "all sigma should share the same in Tidx\n$sigma"

#             #combine sigmas with the same out Tidx
#             df = transform(sigma, :extT => ByRow(x -> [x[1], (x[2], tright),]) => [:Tin, :GT,])
#             allsame(df, :Tin)
#             groups = mergeby(df, :GT, operator = Sum())

#             #for each sigma group, attach a G, all pair ΣG share the same in and out Tidx
#             sigmaG = [ΣG(g, li, para0, tleft) for g in eachrow(groups)]
#             push!(ΣGpair, Diagram(GenericId(para), Sum(), sigmaG))
#         end
#         push!(Gall, Diagram(GenericId(para, (t = extT,)), Prod(), ΣGpair, name = Symbol("g(Σg)^$(length(ΣGpair))")))
#     end
#     # d = Diagram(GenericId(para), Sum(), Gall) # ΣgΣg...Σg all share the same in and out Tidx, thus can be combined into one diagram 
#     d = mergeby(Gall) # ΣgΣg...Σg all share the same in and out Tidx, thus can be combined into one diagram 
#     @assert nrow(d) == 1
#     return Diagram(GreenId(para, k = extK, t = extT), Prod(), [gleft, d.diagram[1]], name = :G)
# end