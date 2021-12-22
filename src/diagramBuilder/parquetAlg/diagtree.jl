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

function buildTree(para::Para, loopNum::Int, legK, Kidx::Int, Tidx::Int, evalK::Function, evalT::Function, factor = 1.0,
    diag = _newDiag(para, legK, evalK), ver4 = Ver4{Int}(para, loopNum, Tidx))


    KinL, KoutL, KinR, KoutR = legK[1], legK[2], legK[3], legK[4]

    @assert KinL + KinR ≈ KoutL + KoutR
    # KoutR = KinL + KinR - KoutL
    GType, VType, WType = 1, 2, 3
    Gorder, Vorder, Worder = 0, 1, 1
    MUL, ADD = 1, 2

    evaT4(Tpair) = Tuple([evalT[t] for t in Tpair])

    qd = KinL - KoutL
    qe = KinR - KoutL
    Tidx = ver4.Tidx

    # print_tree(ver4)
    # println("testing ...")
    # println(ver4.Tpair)

    if ver4.loopNum == 0
        # exit(0)
        if ver4.interactionTauNum == 2
            td = [2, ver4.Tpair[2], evaT4(ver4.Tpair[2])]
            te = [2, ver4.Tpair[3], evaT4(ver4.Tpair[3])]
            println("extT: ", ver4.Tpair[1])
            qd, qe = [1, qd, evalK(qd)], [1, qe, evalK(qe)]
            vd = DiagTree.addPropagator(diag, 2, Vorder, [qd, td])
            ve = DiagTree.addPropagator(diag, 2, Vorder, [qe, te])
            wd = DiagTree.addPropagator(diag, 2, Vorder, [qd, td])
            we = DiagTree.addPropagator(diag, 2, Vorder, [qe, te])
            # time-dependent interaction has different time configurations for the direct and exchange components
            # if 1 in ver4.interactionTauNum
            #     vd = DiagTree.addPropagator!(diag, VType, Vorder, qd, [Tidx, Tidx], Wsym, -1.0)[1]
            #     ve = DiagTree.addPropagator!(diag, VType, Vorder, qe, [Tidx, Tidx], Wsym, 1.0)[1]
            #     ver4.weight[1] = NodeInfo(true, vd, ve)
            # elseif 2 in ver4.interactionTauNum
            #     wd = DiagTree.addPropagator!(diag, WType, Worder, qd, [Tidx, Tidx + 1], Wsym, -1.0)[1]
            #     we = DiagTree.addPropagator!(diag, WType, Worder, qe, [Tidx, Tidx + 1], Wsym, 1.0)[1]
            #     ver4.weight[2] = NodeInfo(true, wd, -1)
            #     ver4.weight[3] = NodeInfo(true, -1, we)
        else
            error("not implemented!")
        end
        # return diag, ver4
    end

    # # LoopNum>=1
    # for w in ver4.weight
    #     w = NodeInfo(false)
    # end

    # K, Kt, Ku, Ks = similar(KinL), similar(KinL), similar(KinL), similar(KinL)
    # G = ver4.G

    # K = zero(KinL)
    # K[Kidx] = 1

    # for c in ver4.chan
    #     if c == T
    #         Kt = KoutL + K - KinL
    #     elseif c == U
    #         Ku = KoutR + K - KinL
    #     else
    #         Ks = KinL + KinR - K
    #     end
    # end

    # for b in ver4.bubble
    #     c = b.chan
    #     # Factor = SymFactor[c] * PhaseFactor
    #     Llopidx = Kidx + 1
    #     Rlopidx = Kidx + 1 + b.Lver.loopNum
    #     Lver, Rver = b.Lver, b.Rver
    #     LLegK, RLegK = [], []
    #     if c == T
    #         LLegK = [KinL, KoutL, Kt, K]
    #         RLegK = [K, Kt, KinR, KoutR]
    #     elseif c == U
    #         LLegK = [KinL, KoutR, Ku, K]
    #         RLegK = [K, Ku, KinR, KoutL]
    #     else
    #         # S channel
    #         LLegK = [KinL, Ks, KinR, K]
    #         RLegK = [K, KoutL, Ks, KoutR]
    #     end
    #     diagramTree(para, Lver.loopNum, LLegK, Llopidx, Lver.Tidx, WeightType, Gsym, Wsym, spin, 1.0, diag, Lver)
    #     diagramTree(para, Rver.loopNum, RLegK, Rlopidx, Rver.Tidx, WeightType, Gsym, Wsym, spin, 1.0, diag, Rver)

    #     rN = length(b.Rver.weight)
    #     for (l, Lw) in enumerate(b.Lver.weight)
    #         for (r, Rw) in enumerate(b.Rver.weight)

    #             map = b.map[(l-1)*rN+r]
    #             g0 = DiagTree.addPropagator!(diag, GType, Gorder, K, collect(G[1].Tpair[map.G0]), Gsym)[1]
    #             gc = DiagTree.addPropagator!(diag, GType, Gorder, K, collect(G[c].Tpair[map.Gx]), Gsym)[1]

    #             # w = (ver4.level == 1 && isFast) ? ver4.weight[ChanMap[c]] : ver4.weight[map.ver]
    #             w = ver4.weight[map.ver]

    #             if c == T || c == U
    #                 #direct
    #                 nsum = []
    #                 ps, ns = split(g0, gc, Lw, Rw, true, true)
    #                 (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, spin * SymFactor[c], ps, ns))
    #                 ps, ns = split(g0, gc, Lw, Rw, true, false)
    #                 (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns))
    #                 ps, ns = split(g0, gc, Lw, Rw, false, true)
    #                 (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns))

    #                 if isempty(nsum) == false
    #                     nt = DiagTree.addNode!(diag, ADD, 1.0, [], nsum)
    #                     # DiagTree.showTree(diag, nt)
    #                     addNode(diag, w, nt, c == T ? true : false) #direct for T, exchange for T
    #                 end

    #                 #exchange
    #                 ps, ns = split(g0, gc, Lw, Rw, false, false)
    #                 if (-1 in ps) == false
    #                     nee = DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns)
    #                     addNode(diag, w, nee, c == T ? false : true) #exchange for T, direct for U
    #                 end
    #                 # DiagTree.showTree(diag)
    #             elseif c == S
    #                 nsum = []
    #                 ps, ns = split(g0, gc, Lw, Rw, true, false)
    #                 (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns))
    #                 ps, ns = split(g0, gc, Lw, Rw, false, true)
    #                 (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns))
    #                 if isempty(nsum) == false
    #                     nd = DiagTree.addNode!(diag, ADD, 1.0, [], nsum)
    #                     addNode(diag, w, nd, true)
    #                 end

    #                 nsum = []
    #                 ps, ns = split(g0, gc, Lw, Rw, true, true)
    #                 (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns))
    #                 ps, ns = split(g0, gc, Lw, Rw, false, false)
    #                 (-1 in ps) || push!(nsum, DiagTree.addNode!(diag, MUL, SymFactor[c], ps, ns))
    #                 if isempty(nsum) == false
    #                     ne = DiagTree.addNode!(diag, ADD, 1.0, [], nsum)
    #                     addNode(diag, w, ne, false)
    #                 end
    #             else
    #                 error("not implemented!")
    #             end
    #         end
    #     end
    # end

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