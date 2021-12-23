################## Generate Expression Tree ########################
# mutable struct verWeight
#     isInteraction::Bool #is a propagator or a node
#     di::Int #index to the direct term
#     ex::Int #index to the exchange term
#     function NodeInfo(isInteraction, di = -1, ex = -1)
#         return new(isPropagator, di, ex)
#     end
# end

# function Base.zero(::Type{NodeInfo})
#     return NodeInfo(false, -1, -1)
# end

# function addNode(diag, node::NodeInfo, nidx, isDirect)
#     MUL, ADD = 1, 2
#     if isDirect
#         if node.di < 0
#             new = DiagTree.addNode!(diag, ADD, 1.0, [], [nidx,])
#             node.di = new
#         else
#             diagnode = diag.tree[node.di]
#             push!(diagnode.nodes, nidx)
#         end
#     else
#         if node.ex < 0
#             new = DiagTree.addNode!(diag, ADD, 1.0, [], [nidx,])
#             node.ex = new
#         else
#             diagnode = diag.tree[node.ex]
#             push!(diagnode.nodes, nidx)
#         end
#     end
# end

# function split(g0, gc, Lw, Rw, isLdirect, isRdirect)
#     propagators = [g0, gc]
#     nodes = []
#     if Lw.isPropagator
#         push!(propagators, isLdirect ? Lw.di : Lw.ex)
#     else
#         push!(nodes, isLdirect ? Lw.di : Lw.ex)
#     end
#     if Rw.isPropagator
#         push!(propagators, isRdirect ? Rw.di : Rw.ex)
#     else
#         push!(nodes, isRdirect ? Rw.di : Rw.ex)
#     end
#     return propagators, nodes
# end

function _newDiag(para::Para, legK, evalK::Function)
    Kbasis = Vector{Float64}
    Kpool = DiagTree.cachedPool(Kbasis, Vector{Float64})
    GTbasis = Tuple{Int,Int}
    GTpool = DiagTree.uncachedPool(GTbasis)
    if para.interactionTauNum == 2
        WTbasis = Tuple{Int,Int}
        WTpool = DiagTree.uncachedPool(WTbasis)
        Gpool = DiagTree.propagatorPool(para.greenType[1], para.greenType[2])
        Wpool = DiagTree.propagatorPool(para.wType[1], para.wType[2])
        return DiagTree.Diagrams((Kpool, GTpool, WTpool), (Gpool, Wpool), para.nodeType[1], para.nodeType[2])
    elseif para.interactionTauNum == 1
        Gpool = DiagTree.propagatorPool(para.greenType[1], para.greenType[2])
        Wpool = DiagTree.propagatorPool(para.wType[1], para.wType[2])
        return DiagTree.Diagrams((Kpool, GTpool), (Gpool, Wpool), para.nodeType[1], para.nodeType[2])
    else
        error("not implemented!")
    end
end

# addNode!(Td, para, diag, lver, rver, DI, EX)
function addNode!(nodes, para, diag, Lver, Rver, l, r, lc, rc, g0, gc, factor = para.nodeType[2](1))
    lLopNum, rLopNum = Lver.loopNum, Rver.loopNum
    Lw, Rw = Lver.weight[l][lc], Rver.weight[r][rc]
    if lLopNum == 0 && rLopNum == 0
        components, child = [[g0, gc], [Lw, Rw]], []
    elseif lLopNum == 0 && rLopNum > 0
        components, child = [[g0, gc], [Lw,]], [Rw,]
    elseif lLopNum > 0 && rLopNum == 0
        components, child = [[g0, gc], [Rw,]], [Lw,]
    else
        components, child = [[g0, gc], []], [Lw, Rw]
    end
    if (Lw != 0 && Rw != 0)
        push!(nodes, DiagTree.addNode(diag, DiagTree.MUL, components, child; factor = factor))
    end
    return nodes
end
# if (chan == T) {
#     W[DIR] = Lw[DIR] * Rw[DIR] * SPIN + Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
#     W[EX] = Lw[EX] * Rw[EX];
#   } else if (chan == U) {
#     W[EX] = Lw[DIR] * Rw[DIR] * SPIN + Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
#     W[DIR] = Lw[EX] * Rw[EX];
#   } else if (chan == S) {
#     // see the note "code convention"
#     W[DIR] = Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
#     W[EX] = Lw[DIR] * Rw[DIR] + Lw[EX] * Rw[EX];
#   }

# function addNode(diag::Diagrams, operator, components, childNodes; factor = 1.0, parent = 0, para = 0, currWeight = 0.0)
function node4Tbubble(para, diag, lver, rver, l, r, g0, gc)
    Factor = SymFactor[T] / (2π)^para.dim
    Td, Te = [], []
    addNode!(Td, para, diag, lver, rver, l, r, DI, DI, g0, gc, para.spin)
    addNode!(Td, para, diag, lver, rver, l, r, DI, EX, g0, gc)
    addNode!(Td, para, diag, lver, rver, l, r, EX, DI, g0, gc)
    nodeTd = DiagTree.addNode(diag, DiagTree.ADD, [[], []], Td; factor = para.nodeType[2](Factor))

    addNode!(Te, para, diag, lver, rver, l, r, EX, EX, g0, gc)
    if isempty(Te) == false
        nodeTe = Te[1]
    else
        nodeTe = 0
    end
    return @SVector [nodeTd, nodeTe]
end

# bubbletoDiagTree!(ver4Nodes, para, diag, ver4, b, legK, Kidx, Tidx, evalK, evalT, factor)
function bubbletoDiagTree!(ver4Nodes, para, diag, ver4, bubble, legK, Kidx::Int, Tidx::Int, evalK::Function, evalT::Function, factor = 1.0)
    KinL, KoutL, KinR, KoutR = legK[1], legK[2], legK[3], legK[4]
    @assert KinL + KinR ≈ KoutL + KoutR

    b = bubble
    c = b.chan
    G = ver4.G

    Gorder, Vorder, Worder = 0, 1, 1
    MUL, ADD = 1, 2

    K = zero(KinL)
    K[Kidx] = 1
    Kt = KoutL + K - KinL
    Ku = KoutR + K - KinL
    Ks = KinL + KinR - K


    # Factor = SymFactor[c] * PhaseFactor
    Llopidx = Kidx + 1
    Rlopidx = Kidx + 1 + b.Lver.loopNum
    Lver, Rver = b.Lver, b.Rver
    LLegK, RLegK = [], []
    if c == T
        LLegK = [KinL, KoutL, Kt, K]
        RLegK = [K, Kt, KinR, KoutR]
    elseif c == U
        LLegK = [KinL, KoutR, Ku, K]
        RLegK = [K, Ku, KinR, KoutL]
    else
        # S channel
        LLegK = [KinL, Ks, KinR, K]
        RLegK = [K, KoutL, Ks, KoutR]
    end
    ver4toDiagTree(para, Lver.loopNum, LLegK, Llopidx, Lver.Tidx, evalK, evalT, factor, diag, Lver)
    ver4toDiagTree(para, Rver.loopNum, RLegK, Rlopidx, Rver.Tidx, evalK, evalT, factor, diag, Rver)

    rN = length(b.Rver.weight)
    for (l, Lw) in enumerate(b.Lver.weight)
        for (r, Rw) in enumerate(b.Rver.weight)

            map = b.map[(l-1)*rN+r]
            # println("evalT: ", G[1].Tpair[map.G0])
            # evaTpair(Tpair) = Tuple([evalT(t) for t in Tpair])
            tpair = G[1].Tpair[map.G0]
            tbasis = (2, tpair, (evalT(tpair[1]), evalT(tpair[2])))
            g0 = DiagTree.addPropagator(diag, 1, Gorder, [(1, K, eval(K)), tbasis])
            if c == T
                Kc = Kt
            elseif c == U
                Kc = Ku
            elseif c == S
                Kc = Ks
            end
            tpair = G[c].Tpair[map.Gx]
            tbasis = (2, tpair, (evalT(tpair[1]), evalT(tpair[2])))
            gc = DiagTree.addPropagator(diag, 1, Gorder, [(1, Kc, eval(Kc)), tbasis])

            # w = (ver4.level == 1 && isFast) ? ver4.weight[ChanMap[c]] : ver4.weight[map.ver]
            w = ver4.weight[map.ver]

            if c == T
                Tde = node4Tbubble(para, diag, b.Lver, b.Rver, l, r, g0, gc)
                push!(ver4Nodes[map.ver], Tde)
                # if c == T || c == U
                #direct


                # nsum = []
                # ps, ns = split(g0, gc, Lw, Rw, true, true)
                # (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, spin * SymFactor[c], ps, ns))
                # ps, ns = split(g0, gc, Lw, Rw, true, false)
                # (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns))
                # ps, ns = split(g0, gc, Lw, Rw, false, true)
                # (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns))

                # if isempty(nsum) == false
                #     nt = DiagTree.addNode!(diag, ADD, 1.0, [], nsum)
                #     # DiagTree.showTree(diag, nt)
                #     addNode(diag, w, nt, c == T ? true : false) #direct for T, exchange for T
                # end

                # #exchange
                # ps, ns = split(g0, gc, Lw, Rw, false, false)
                # if (-1 in ps) == false
                #     nee = DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns)
                #     addNode(diag, w, nee, c == T ? false : true) #exchange for T, direct for U
                # end
                # DiagTree.showTree(diag)
            elseif c == S
                # nsum = []
                # ps, ns = split(g0, gc, Lw, Rw, true, false)
                # (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns))
                # ps, ns = split(g0, gc, Lw, Rw, false, true)
                # (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns))
                # if isempty(nsum) == false
                #     nd = DiagTree.addNode!(diag, ADD, 1.0, [], nsum)
                #     addNode(diag, w, nd, true)
                # end

                # nsum = []
                # ps, ns = split(g0, gc, Lw, Rw, true, true)
                # (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns))
                # ps, ns = split(g0, gc, Lw, Rw, false, false)
                # (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns))
                # if isempty(nsum) == false
                #     ne = DiagTree.addNode!(diag, ADD, 1.0, [], nsum)
                #     addNode(diag, w, ne, false)
                # end
            else
                error("not implemented!")
            end
        end
    end
    return diag
end

function ver4toDiagTree(para::Para, loopNum::Int, legK, Kidx::Int, Tidx::Int, evalK::Function, evalT::Function, factor = 1.0,
    diag = _newDiag(para, legK, evalK), ver4 = Ver4{SVector{2,Int}}(para, loopNum, Tidx))


    KinL, KoutL, KinR, KoutR = legK[1], legK[2], legK[3], legK[4]

    @assert KinL + KinR ≈ KoutL + KoutR
    # KoutR = KinL + KinR - KoutL
    GType, VType, WType = 1, 2, 3
    Gorder, Vorder, Worder = 0, 1, 1
    MUL, ADD = 1, 2

    evaTpair(Tpair) = Tuple([evalT(t) for t in Tpair])

    qd = KinL - KoutL
    qe = KinR - KoutL
    Tidx = ver4.Tidx

    if ver4.loopNum == 0
        qd, qe = [1, qd, evalK(qd)], [1, qe, evalK(qe)]
        if ver4.interactionTauNum == 2
            td = [2, (Tidx, Tidx + 1), evaTpair((Tidx, Tidx + 1))]
            te = td
            vd = DiagTree.addPropagator(diag, 2, Vorder, [qd, td])
            ve = DiagTree.addPropagator(diag, 2, Vorder, [qe, te])
            wd = DiagTree.addPropagator(diag, 2, Vorder, [qd, td])
            we = DiagTree.addPropagator(diag, 2, Vorder, [qe, te])
            ver4.weight[1] = @SVector [vd, ve]
            ver4.weight[2] = @SVector [wd, 0]
            ver4.weight[3] = @SVector [0, we]
        elseif ver4.interactionTauNum == 1
            vd = DiagTree.addPropagator(diag, 2, Vorder, [qd,])
            ve = DiagTree.addPropagator(diag, 2, Vorder, [qe,])
            ver4.weight[1] = @SVector [vd, ve]
        else
            error("not implemented!")
        end
        return diag, ver4
    end

    # LoopNum>=1
    ver4Nodes = []
    for i in 1:length(ver4.weight)
        ver4.weight[i] = @SVector [0, 0]
        push!(ver4Nodes, [])
    end


    for b in ver4.bubble
        bubbletoDiagTree!(ver4Nodes, para, diag, ver4, b, legK, Kidx, Tidx, evalK, evalT, factor)
    end

    # nodeTd = DiagTree.addNode(diag, DiagTree.ADD, [[], []], Td; factor = para.nodeFactorType(Factor))
    for i in 1:length(ver4.weight)
        ver4dNodes = []
        ver4eNodes = []
        for n in ver4Nodes[i]
            if n[DI] != 0
                push!(ver4dNodes, n[DI])
            end
            if n[EX] != 0
                push!(ver4eNodes, n[EX])
            end
        end
        nodeD = DiagTree.addNode(diag, DiagTree.ADD, [[], []], ver4dNodes; factor = para.nodeType[2](factor))
        nodeE = DiagTree.addNode(diag, DiagTree.ADD, [[], []], ver4eNodes; factor = para.nodeType[2](factor))
        ver4.weight[i] = @SVector [nodeD, nodeE]
    end

    return diag, ver4
end