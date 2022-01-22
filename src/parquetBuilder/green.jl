"""
    function buildG(para, extK, extT, subdiagram = false; F = [I, U, S], V = [I, T, U], All = union(F, V), diag = newDiagTree(para, :G))
    
    Build composite Green's function.
    By definition, para.firstTauIdx is the first Tau index of the left most self-energy subdiagram.

- `extT`: [Tau index of the left leg, Tau index of the right leg]

"""
function buildG(para, extK, extT, subdiagram = false; name = :none)
    @assert para.diagType == GreenDiag
    tin, tout = extT[1], extT[2]
    extT = [tin, tout]
    innerTauNum = para.innerLoopNum * para.interactionTauNum
    tstart = para.firstTauIdx
    tend = para.firstTauIdx + innerTauNum - 1

    diags = Diagram{para.weightType}[]

    if isValidG(para) == false
        return nothing
    end

    if para.innerLoopNum == 0
        return Diagram(GreenId(para, k = extK, t = extT), name = name)
    end

    return nothing

    # ################# after this step, the Green's function must be nontrivial! ##################
    # @assert tin < tstart || tin > tend "external T index cann't be with in [$tstart, $tend]"
    # @assert tout < tstart || tout > tend "external T index cann't be with in [$tstart, $tend]"

    # gleft = DiagTree.addpropagator!(diag, :Gpool, 0, :gleft; site = [tin, tstart], loop = extK)

    # partition = []
    # for n = 1:para.innerLoopNum
    #     #generate all possible ordered partition of the loops
    #     append!(partition, orderedPartition(para.innerLoopNum, n))
    #     #e.g., loopNum =5, n =2 ==> ordered = [[4, 1], [1, 4], [3, 2], [2, 3]]
    # end
    # #e.g., loopNum =5 ==> partition = [[4, 1], [1, 4], [3, 2], [2, 3], [1, 1, 1, 1], ...]

    # Gall = [gleft,]
    # # println(partition)
    # for p in partition
    #     #e.g., p = [4, 1]
    #     firstTauIdx, maxTau = findFirstTauIdx(p, [SigmaDiag for i in p], para.firstTauIdx, para.interactionTauNum)
    #     firstLoopIdx, maxLoop = findFirstLoopIdx(p, para.firstLoopIdx)

    #     if all([isValidSigma(para.filter, loop, true) for loop in p]) == false
    #         continue
    #     end

    #     Gnode = []
    #     for (li, loop) in enumerate(p)
    #         tleft = firstTauIdx[li]
    #         tright = li < length(p) ? firstTauIdx[li+1] : tout
    #         #e.g., loop = 4 or 1
    #         sigmaPara = reconstruct(para, firstTauIdx = tleft, firstLoopIdx = firstLoopIdx[li], innerLoopNum = loop)
    #         # println("sigma loop=", loop)
    #         diag, instant, dynamic = buildSigma(sigmaPara, extK, true; F = F, V = V, All = All, diag = diag)
    #         root = vcat(instant, dynamic)
    #         @assert isempty(root) == false
    #         # push!(sigma, [(s, (s.object.para[1], s.object.para[2])) for s in root]) #external TauIdx of each component as a sigma

    #         nodes = []
    #         for r in root
    #             #build node Σ(tleft, t)*g(t, tright)
    #             t1, t2 = r.object.para
    #             # println(t1, ", ", t2, ", ", loop)
    #             @assert t1 == tleft "sigma left Tidx = $t1 is not equal to the firstTauIdx = $tleft"
    #             if li < length(p) #not the last g
    #                 @assert t2 < tright "sigma right Tidx = $t2 should be smaller than $tright"
    #             end
    #             g = DiagTree.addpropagator!(diag, :Gpool, 0, :g; site = [t2, tright], loop = extK)
    #             push!(nodes, DiagTree.addnode!(diag, MUL, :Σg, [r, g]; para = [tleft, tright]))
    #         end
    #         n = DiagTree.addnode!(diag, ADD, :Σgsum, nodes; para = [tleft, tright])
    #         #build node \sum_t Σ(tleft, t)*g(t, tright)
    #         @assert n.index > 0
    #         push!(Gnode, n)
    #     end
    #     n = DiagTree.addnode!(diag, MUL, :gΣg, Gnode; para = extT)
    #     @assert n.index > 0
    #     push!(Gall, n)
    # end
    # G = DiagTree.addnode!(diag, ADD, :BoldG, Gall; para = extT)
    # @assert G.index > 0
    # # DiagTree.showTree(diag, Gidx)
    # return diag, G
end