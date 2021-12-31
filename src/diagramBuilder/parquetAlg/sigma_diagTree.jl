
"""
    function buildSigma(para, externLoop; F = [I, U, S], V = [I, T, U], All = union(F, V), diag = newDiagTree(para, :Sigma), subdiagram = false)
    
    build sigma diagram. 
    When sigma is created as a subdiagram, then no Fock diagram is generated if para.filter contains NoFock, and no sigma diagram is generated if para.filter contains Girreducible

"""
function buildSigma(para, externLoop, subdiagram = false; F = [I, U, S], V = [I, T, U], All = union(F, V), diag = newDiagTree(para, :Sigma))
    @assert para.innerLoopNum >= 1
    @assert length(externLoop) == para.totalLoopNum
    tright = para.firstTauIdx - 1 + para.innerLoopNum * para.interactionTauNum
    @assert para.totalTauNum >= tright "totalTauNum = $(para.totalTauNum) is not enough, sigma requires $tright\npara=$para"
    @assert para.totalLoopNum >= para.firstLoopIdx -1 + para.innerLoopNum

    K = zero(externLoop)
    K[para.firstLoopIdx] = 1.0
    t0 = para.firstTauIdx

    function collapse!(dict, diag, nodes, factor)
        for nidx in nodes
            node = diag.nodePool.object[nidx]
            extT = node.para
            sigmaT = (extT[INL], extT[OUTR])
            if haskey(dict, sigmaT) == false
                dict[sigmaT] = []
            end
            push!(dict[sigmaT], (nidx, extT, factor))
        end
    end

    if subdiagram && (Girreducible in para.filter)
        return diag, []
    end

    root = []

    ############### Fock-type diagram #######################################
    #  /=== W ===\
    # /           \
    # ----- G -----
    ##########################################################################
    qe = K - externLoop

    if para.interactionTauNum == 1
        factor = 1 / (2π)^para.loopDim
        # g = DiagTree.addPropagator!(diag, :Gpool, 0, :Gsigma; loop = K, site = (t0, t0))
        paraG = reconstruct(para, innerLoopNum = para.innerLoopNum - 1, firstLoopIdx = para.firstLoopIdx + 1, firstTauIdx = t0 + 1)
        diag, g = buildG(paraG, K, [t0, t0]; diag = diag)
        if g != 0
            # @assert g != 0 "Get invalid G index in GW diagram with G para =\n $paraG"
            v = DiagTree.addPropagator!(diag, :Vpool, 1, :Vsigma; loop = qe)
            n = DiagTree.addNodeByName!(diag, DiagTree.MUL, :GV, factor; Gpool = g, Vpool = v, para = [t0, t0])
            push!(root, n)
        end
    elseif para.interactionTauNum == 2
        factor = 1 / (2π)^para.loopDim
        paraG = reconstruct(para, innerLoopNum = para.innerLoopNum - 1, firstLoopIdx = para.firstLoopIdx + 1, firstTauIdx = t0 + 2)
        diag, gv = buildG(paraG, K, [t0, t0]; diag = diag)
        if gv != 0
            # @assert gv != 0 "Get invalid G index in GW diagram with G para =\n $paraG"
            # gv = DiagTree.addPropagator!(diag, :Gpool, 0, :Gsigma; loop = K, site = (t0, t0))
            v = DiagTree.addPropagator!(diag, :Vpool, 1, :Vsigma; loop = qe, site = (t0, t0 + 1))
            nv = DiagTree.addNodeByName!(diag, DiagTree.MUL, :GV, factor; Gpool = gv, Vpool = v, para = (t0, t0))
            push!(root, nv)
        end

        diag, gw = buildG(paraG, K, [t0, t0 + 1]; diag = diag)
        if gw != 0
            # @assert gw != 0 "Get invalid G index in GW diagram"
            gw = DiagTree.addPropagator!(diag, :Gpool, 0, :Gsigma; loop = K, site = (t0, t0 + 1))
            w = DiagTree.addPropagator!(diag, :Wpool, 1, :Wsigma; loop = qe, site = (t0, t0 + 1))

            nw = DiagTree.addNodeByName!(diag, DiagTree.MUL, :GW, factor; Gpool = gw, Wpool = w, para = [t0, t0 + 1])
            push!(root, nw)
        end
    else
        error("not implemented!")
    end

    if para.innerLoopNum >= 2
        ################# Sigma beyond the Fock diagram #################################
        #all self-energy diagram will be countered twice, thus a factor 1/2 is needed.
        factor = 1 / (2π)^para.loopDim * (1 / 2)
        KinL, KoutR = externLoop, externLoop
        KinR, KoutL = K, K
        for ver4LoopNum in 1:para.innerLoopNum-1
            ver4Para = reconstruct(para, firstLoopIdx = para.firstLoopIdx + 1, innerLoopNum = ver4LoopNum)
            diag, ver4, dir, ex = buildVer4(ver4Para, [KinL, KoutL, KinR, KoutR],
                [T,], F, V, All; Fouter = [], Allouter = All, diag = diag)

            dict = Dict{Tuple{Int,Int},Vector{Any}}()
            collapse!(dict, diag, dir, 1.0)
            collapse!(dict, diag, ex, para.spin)
            #exchange ver4 has an additional Fermi loop compared to the direct counterpart when plugged into sigma

            for key in keys(dict)
                nodes = []
                for (nidx, extT, spinFactor) in dict[key]

                    tpair = [extT[INL], extT[OUTR]]
                    paraG = reconstruct(para, innerLoopNum = para.innerLoopNum - 1 - ver4LoopNum,
                        firstLoopIdx = para.firstLoopIdx + 1 + ver4LoopNum, firstTauIdx = maxTauIdx(ver4) + 1)
                    diag, g = buildG(paraG, K, tpair; diag = diag)
                    # println(g, ", from", paraG)
                    # @assert g != 0 "Get invalid G index in GW diagram"
                    # g = DiagTree.addPropagator!(diag, :Gpool, 0, :Gsigma; site = (extT[OUTL], extT[INR]), loop = K)
                    if g != 0
                        push!(nodes, DiagTree.addNodeByName!(diag, DiagTree.MUL, :sigma, factor * spinFactor;
                            child = nidx, Gpool = g, para = tpair))
                    end
                end

                if isempty(nodes) == false
                    rootidx = DiagTree.addNode!(diag, DiagTree.ADD, :sigma;
                        child = nodes, para = collect(key)) #key = (extT[INL], extT[OUTR])
                    # println(rootidx)
                    push!(root, rootidx)
                end
            end
        end
    end

    # make sure all incoming Tau idx is equal
    for nidx in root
        node = DiagTree.getNode(diag, nidx)
        @assert node.para[1] == para.firstTauIdx
        @assert node.para[2] <= tright
    end
    @assert isempty(root) == false "Sigma diagram doesn't exist for \n$para"
    return diag, root
end

function orderedPartition(total, n)
    unorderedPartition = collect(partitions(total, n))
    #e.g., loopNum =5, n =2 ==> unordered = [[4, 1], [3, 2]]
    orderedPartition = Vector{Vector{Int}}([])
    for p in unorderedPartition
        append!(orderedPartition, Set(permutations(p)))
    end
    #e.g., loopNum =5, n =2 ==> ordered = [[4, 1], [1, 4], [3, 2], [2, 3]]
    return orderedPartition
end

"""
    function buildG(para, externLoop, extT, subdiagram = false; F = [I, U, S], V = [I, T, U], All = union(F, V), diag = newDiagTree(para, :G))
    
    Build composite Green's function.
    By definition, para.firstTauIdx is the first Tau index of the left most self-energy subdiagram.

- `extT`: [Tau index of the left leg, Tau index of the right leg]

"""
function buildG(para, externLoop, extT, subdiagram = false; F = [I, U, S], V = [I, T, U], All = union(F, V), diag = newDiagTree(para, :G))
    tin, tout = extT[1], extT[2]
    innerTauNum = para.innerLoopNum * para.interactionTauNum
    tstart = para.firstTauIdx
    tend = para.firstTauIdx + innerTauNum - 1

    if (NoFock in para.filter) && para.innerLoopNum == 1
        return diag, 0
    end

    if (Girreducible in para.filter) && para.innerLoopNum > 0
        return diag, 0
    end

    if para.innerLoopNum == 0
        g = DiagTree.addPropagator!(diag, :Gpool, 0, :G; site = [tin, tout], loop = externLoop)
        return diag, g
    end

    ################# after this step, the Green's function must be nontrivial! ##################
    @assert tin < tstart || tin > tend "external T index cann't be with in [$tstart, $tend]"
    @assert tout < tstart || tout > tend "external T index cann't be with in [$tstart, $tend]"

    gleft = DiagTree.addPropagator!(diag, :Gpool, 0, :gleft; site = [tin, tstart], loop = externLoop)

    partition = []
    for n = 1:para.innerLoopNum
        #generate all possible ordered partition of the loops
        append!(partition, orderedPartition(para.innerLoopNum, n))
        #e.g., loopNum =5, n =2 ==> ordered = [[4, 1], [1, 4], [3, 2], [2, 3]]
    end
    #e.g., loopNum =5 ==> partition = [[4, 1], [1, 4], [3, 2], [2, 3], [1, 1, 1, 1], ...]

    Gall = []
    # println(partition)
    for p in partition
        #e.g., p = [4, 1]
        Gnode = []
        firstTauIdx = para.firstTauIdx
        firstLoopIdx = para.firstLoopIdx
        for (li, loop) in enumerate(p)
            #e.g., loop = 4 or 1
            sigmaPara = reconstruct(para, firstTauIdx = firstTauIdx, firstLoopIdx = firstLoopIdx, innerLoopNum = loop)
            diag, root = buildSigma(sigmaPara, externLoop, true; F = F, V = V, All = All, diag = diag)
            tleft = firstTauIdx
            tright = (li == length(p)) ? tout : tleft + loop * para.interactionTauNum

            nodes = []
            for r in root
                #build node Σ(tleft, t)*g(t, tright)
                n = DiagTree.getNode(diag, r)
                t1, t2 = n.para
                # println(t1, ", ", t2, ", ", loop)
                @assert t1 == firstTauIdx "sigma left Tidx = $t1 is not equal to the firstTauIdx = $firstTauIdx"
                if li < length(p) #not the last g
                    @assert t2 < tright "sigma right Tidx = $t2 should be smaller than $tright"
                end
                g = DiagTree.addPropagator!(diag, :Gpool, 0, :g; site = [t2, tright], loop = externLoop)
                push!(nodes, DiagTree.addNodeByName!(diag, DiagTree.MUL, :Σg, child = r, Gpool = g, para = [tleft, tright]))
            end

            if isempty(nodes) == false
                #build node \sum_t Σ(tleft, t)*g(t, tright)
                push!(Gnode, DiagTree.addNode!(diag, DiagTree.ADD, :Σgsum, child = nodes, para = [tleft, tright]))
            end

            firstLoopIdx += loop
            firstTauIdx += loop * para.interactionTauNum
        end
        if isempty(Gnode) == false
            push!(Gall, DiagTree.addNodeByName!(diag, DiagTree.MUL, :G; child = Gnode, Gpool = gleft, para = [tin, tout]))
        end
    end

    @assert isempty(Gall) == false

    Gidx = DiagTree.addNodeByName!(diag, DiagTree.ADD, :Gsum; child = Gall, para = [tin, tout])
    return diag, Gidx
end