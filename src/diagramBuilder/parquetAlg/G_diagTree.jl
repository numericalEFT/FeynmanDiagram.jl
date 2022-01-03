# check if G exist without creating objects in the pool
function isValidG(filter, innerLoopNum::Int)
    if (NoFock in filter) && innerLoopNum == 1
        return false
    end

    if (Girreducible in filter) && innerLoopNum > 0
        return false
    end

    return true
end

function isValidG(para::GenericPara)
    return isValidG(para.filter, para.innerLoopNum)
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

    if isValidG(para) == false
        return diag, zero(Component)
    end

    if para.innerLoopNum == 0
        g = DiagTree.addpropagator!(diag, :Gpool, 0, :G; site = [tin, tout], loop = externLoop)
        return diag, g
    end

    return diag, zero(Component)

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
            # println("sigma loop=", loop)
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
            push!(Gall, DiagTree.addNodeByName!(diag, DiagTree.MUL, :gΣg; child = Gnode, Gpool = gleft, para = [tin, tout]))
        end
    end

    @assert isempty(Gall) == false

    Gidx = DiagTree.addNodeByName!(diag, DiagTree.ADD, :gΣg_sum; child = Gall, para = [tin, tout])
    # DiagTree.showTree(diag, Gidx)
    return diag, Gidx, true
end