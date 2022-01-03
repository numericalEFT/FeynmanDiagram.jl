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

    # gleft = DiagTree.addPropagator!(diag, :Gpool, 0, :gleft; site = [tin, tstart], loop = externLoop)

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
        firstTauIdx = findFirstTauIdx(p, [SigmaDiag for i in p], para.firstTauIdx, para.interactionTauNum)
        firstLoopIdx = findFirstLoopIdx(p, para.firstLoopIdx)

        if all([isValidSigma(para.filter, loop, true) for loop in p]) == false
            continue
        end

        # Gnode = []
        sigma = []
        extT = []
        for (li, loop) in enumerate(p)
            #e.g., loop = 4 or 1
            sigmaPara = reconstruct(para, firstTauIdx = firstTauIdx[li], firstLoopIdx = firstLoopIdx[li], innerLoopNum = loop)
            # println("sigma loop=", loop)
            diag, root = buildSigma(sigmaPara, externLoop, true; F = F, V = V, All = All, diag = diag)
            @assert isempty(root) == false
            push!(sigma, root)
            push!(extT, [(s.object.para[1], s.object.para[2]) for s in root]) #external TauIdx of each component as a sigma
        end

    end

    @assert isempty(Gall) == false

    Gidx = DiagTree.addNodeByName!(diag, DiagTree.ADD, :gΣg_sum; child = Gall, para = [tin, tout])
    # DiagTree.showTree(diag, Gidx)
    return diag, Gidx, true
end