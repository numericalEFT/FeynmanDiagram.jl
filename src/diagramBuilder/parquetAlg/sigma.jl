
function buildSigma(para, externLoop; F = [I, U, S], V = [I, T, U], All = union(F, V), diag = newDiagTree(para, :Sigma), subdiagram = false)
    @assert para.innerLoopNum >= 1

    @assert length(externLoop) == para.totalLoopNum
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
        return diag
    end

    if para.innerLoopNum == 1
        # Fock diagram

        #if it is a Fock subdiagram, then check NoFock filter
        if subdiagram && (NoFock in para.filter)
            return diag
        end

        factor = 1 / (2π)^para.loopDim

        qe = K - externLoop
        if para.interactionTauNum == 1
            g = DiagTree.addPropagator!(diag, :Gpool, 0, :Gsigma; loop = K, site = (t0, t0))
            v = DiagTree.addPropagator!(diag, :Vpool, 1, :Vsigma; loop = qe)
            n = DiagTree.addNodeByName!(diag, DiagTree.MUL, :GV, factor; Gpool = g, Vpool = v, para = [t0, t0])
            return diag, [n,]
        elseif para.interactionTauNum == 2
            gv = DiagTree.addPropagator!(diag, :Gpool, 0, :Gsigma; loop = K, site = (t0, t0))
            v = DiagTree.addPropagator!(diag, :Vpool, 1, :Vsigma; loop = qe, site = (t0, t0 + 1))
            nv = DiagTree.addNodeByName!(diag, DiagTree.MUL, :GV, factor; Gpool = gv, Vpool = v, para = (t0, t0))

            gw = DiagTree.addPropagator!(diag, :Gpool, 0, :Gsigma; loop = K, site = (t0, t0 + 1))
            w = DiagTree.addPropagator!(diag, :Wpool, 1, :Wsigma; loop = qe, site = (t0, t0 + 1))

            nw = DiagTree.addNodeByName!(diag, DiagTree.MUL, :GW, factor; Gpool = gw, Wpool = w, para = [t0, t0 + 1])
            return diag, [nv, nw]
        else
            error("not implemented!")
        end
    elseif para.innerLoopNum >= 2
        #all self-energy diagram will be countered twice, thus a factor 1/2 is needed.
        factor = 1 / (2π)^para.loopDim * (1 / 2)
        KinL, KoutR = externLoop, externLoop
        KinR, KoutL = K, K
        ver4Para = reconstruct(para, firstLoopIdx = para.firstLoopIdx + 1, innerLoopNum = para.innerLoopNum - 1)
        diag, ver4, dir, ex = buildVer4(ver4Para, [KinL, KoutL, KinR, KoutR],
            [T,], F, V, All; Fouter = [], Allouter = All, diag = diag)

        dict = Dict{Tuple{Int,Int},Vector{Any}}()
        collapse!(dict, diag, dir, 1.0)
        collapse!(dict, diag, ex, para.spin)
        #exchange ver4 has an additional Fermi loop compared to the direct counterpart when plugged into sigma

        root = []
        for key in keys(dict)
            nodes = []
            for (nidx, extT, spinFactor) in dict[key]
                g = DiagTree.addPropagator!(diag, :Gpool, 0, :Gsigma; site = (extT[OUTL], extT[INR]), loop = K)
                push!(nodes, DiagTree.addNodeByName!(diag, DiagTree.MUL, :sigma, factor * spinFactor;
                    child = nidx, Gpool = g, para = [extT[INL], extT[OUTR]]))
            end

            rootidx = DiagTree.addNode!(diag, DiagTree.ADD, :sigma;
                child = nodes, para = collect(key)) #key = (extT[INL], extT[OUTR])
            # println(rootidx)
            push!(root, rootidx)
        end
        return diag, root
    end
end