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

    # print_tree(ver4)
    # println("testing ...")
    # println(ver4.Tpair)

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
    for i in 1:length(ver4.weight)
        ver4.weight[i] = @SVector [0, 0]
    end

    G = ver4.G

    K = zero(KinL)
    K[Kidx] = 1
    Kt = KoutL + K - KinL
    Ku = KoutR + K - KinL
    Ks = KinL + KinR - K

    for b in ver4.bubble
        c = b.chan
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

                # if c == T || c == U
                #     #direct
                #     nsum = []
                #     ps, ns = split(g0, gc, Lw, Rw, true, true)
                #     (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, spin * SymFactor[c], ps, ns))
                #     ps, ns = split(g0, gc, Lw, Rw, true, false)
                #     (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns))
                #     ps, ns = split(g0, gc, Lw, Rw, false, true)
                #     (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns))

                #     if isempty(nsum) == false
                #         nt = DiagTree.addNode!(diag, ADD, 1.0, [], nsum)
                #         # DiagTree.showTree(diag, nt)
                #         addNode(diag, w, nt, c == T ? true : false) #direct for T, exchange for T
                #     end

                #     #exchange
                #     ps, ns = split(g0, gc, Lw, Rw, false, false)
                #     if (-1 in ps) == false
                #         nee = DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns)
                #         addNode(diag, w, nee, c == T ? false : true) #exchange for T, direct for U
                #     end
                #     # DiagTree.showTree(diag)
                # elseif c == S
                #     nsum = []
                #     ps, ns = split(g0, gc, Lw, Rw, true, false)
                #     (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns))
                #     ps, ns = split(g0, gc, Lw, Rw, false, true)
                #     (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns))
                #     if isempty(nsum) == false
                #         nd = DiagTree.addNode!(diag, ADD, 1.0, [], nsum)
                #         addNode(diag, w, nd, true)
                #     end

                #     nsum = []
                #     ps, ns = split(g0, gc, Lw, Rw, true, true)
                #     (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns))
                #     ps, ns = split(g0, gc, Lw, Rw, false, false)
                #     (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns))
                #     if isempty(nsum) == false
                #         ne = DiagTree.addNode!(diag, ADD, 1.0, [], nsum)
                #         addNode(diag, w, ne, false)
                #     end
                # else
                #     error("not implemented!")
                # end
            end
        end
    end

    # if ver4.level == 1
    #     for w in ver4.weight
    #         w.di > 0 && push!(diag.root, w.di)
    #         w.ex > 0 && push!(diag.root, w.ex)
    #     end
    # end
    return diag, ver4
end

# function eval(ver4::Ver4, KinL, KoutL, KinR, KoutR, Kidx::Int, fast = false)
#     if ver4.loopNum == 0
#         ver4.weight[1] = interaction(KinL, KoutL, KinR, KoutR, ver4.inBox, norm(varK[0])) :
#         ver4.weight[1] = interaction(KinL, KoutL, KinR, KoutR, ver4.inBox)
#         return
#     end

#     # LoopNum>=1
#     for w in ver4.weight
#         w .= 0.0 # initialize all weights
#     end
#     G = ver4.G
#     K, Kt, Ku, Ks = (varK[Kidx], ver4.K[1], ver4.K[2], ver4.K[3])
#     eval(G[1], K, varT)
#     bubWeight = counterBubble(K)

#     for c in ver4.chan
#         if c == T || c == TC
#             Kt .= KoutL .+ K .- KinL
#             if (!ver4.inBox)
#                 eval(G[T], Kt)
#             end
#         elseif c == U || c == UC
#             # can not be in box!
#             Ku .= KoutR .+ K .- KinL
#             eval(G[U], Ku)
#         else
#             # S channel, and cann't be in box!
#             Ks .= KinL .+ KinR .- K
#             eval(G[S], Ks)
#         end
#     end
#     for b in ver4.bubble
#         c = b.chan
#         Factor = SymFactor[c] * PhaseFactor
#         Llopidx = Kidx + 1
#         Rlopidx = Kidx + 1 + b.Lver.loopNum

#         if c == T || c == TC
#             eval(b.Lver, KinL, KoutL, Kt, K, Llopidx)
#             eval(b.Rver, K, Kt, KinR, KoutR, Rlopidx)
#         elseif c == U || c == UC
#             eval(b.Lver, KinL, KoutR, Ku, K, Llopidx)
#             eval(b.Rver, K, Ku, KinR, KoutL, Rlopidx)
#         else
#             # S channel
#             eval(b.Lver, KinL, Ks, KinR, K, Llopidx)
#             eval(b.Rver, K, KoutL, Ks, KoutR, Rlopidx)
#         end

#         rN = length(b.Rver.weight)
#         gWeight = 0.0
#         for (l, Lw) in enumerate(b.Lver.weight)
#             for (r, Rw) in enumerate(b.Rver.weight)
#                 map = b.map[(l-1)*rN+r]

#                 if ver4.inBox || c == TC || c == UC
#                     gWeight = bubWeight * Factor
#                 else
#                     gWeight = G[1].weight[map.G] * G[c].weight[map.Gx] * Factor
#                 end

#                 if fast && ver4.level == 0
#                     pair = ver4.Tpair[map.ver]
#                     dT =
#                         varT[pair[INL]] - varT[pair[OUTL]] + varT[pair[INR]] -
#                         varT[pair[OUTR]]
#                     gWeight *= cos(2.0 * pi / Beta * dT)
#                     w = ver4.weight[ChanMap[c]]
#                 else
#                     w = ver4.weight[map.ver]
#                 end

#                 if c == T || c == TC
#                     w[DI] +=
#                         gWeight *
#                         (Lw[DI] * Rw[DI] * SPIN + Lw[DI] * Rw[EX] + Lw[EX] * Rw[DI])
#                     w[EX] += gWeight * Lw[EX] * Rw[EX]
#                 elseif c == U || c == UC
#                     w[DI] += gWeight * Lw[EX] * Rw[EX]
#                     w[EX] +=
#                         gWeight *
#                         (Lw[DI] * Rw[DI] * SPIN + Lw[DI] * Rw[EX] + Lw[EX] * Rw[DI])
#                 else
#                     # S channel,  see the note "code convention"
#                     w[DI] += gWeight * (Lw[DI] * Rw[EX] + Lw[EX] * Rw[DI])
#                     w[EX] += gWeight * (Lw[DI] * Rw[DI] + Lw[EX] * Rw[EX])
#                 end

#             end
#         end

#     end
# end