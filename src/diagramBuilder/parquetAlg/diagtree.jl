################## Generate Expression Tree ########################
function _newDiag(para::Para, legK, evalK::Function)
    # function LoopPool(name::Symbol, dim::Int, N::Int, type::DataType)
    Kpool = DiagTree.LoopPool(:K, para.dim, para.loopNum, Float64)
    Tbasis = Tuple{Int,Int}
    Tpool = DiagTree.uncachedPool(Tbasis)
    weightType = para.weightType
    if para.interactionTauNum == 2
        Gpool = DiagTree.propagatorPool(:G, weightType; paraType = Symbol)
        Wpool = DiagTree.propagatorPool(:VW, weightType, paraType = Symbol)
        return DiagTree.Diagrams((Kpool, Tpool), (Gpool, Wpool), weightType, nodeParaType = Symbol)
    elseif para.interactionTauNum == 1
        Gpool = DiagTree.propagatorPool(:G, weightType, paraType = Symbol)
        Wpool = DiagTree.propagatorPool(:V, weightType, paraType = Symbol)
        return DiagTree.Diagrams((Kpool, Tpool), (Gpool, Wpool), weightType, nodeParaType = Symbol)
    else
        error("not implemented!")
    end
end

# addNode!(Td, para, diag, lver, rver, DI, EX)
function addNode!(name::Symbol, nodes, para, diag, Lver, Rver, l, r, lc, rc, g0, gc; factor = 1)
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
        push!(nodes, DiagTree.addNode(diag, DiagTree.MUL, components, child; factor = factor, para = name))
    end
    return nodes
end

function node4Tbubble(para, diag, lver, rver, l, r, g0, gc)
    #     W[DIR] = Lw[DIR] * Rw[DIR] * SPIN + Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
    #     W[EX] = Lw[EX] * Rw[EX];
    Factor = SymFactor[T] / (2π)^para.dim
    Td, Te = [], []
    addNode!(:dxd, Td, para, diag, lver, rver, l, r, DI, DI, g0, gc; factor = para.spin)
    addNode!(:dxe, Td, para, diag, lver, rver, l, r, DI, EX, g0, gc)
    addNode!(:exd, Td, para, diag, lver, rver, l, r, EX, DI, g0, gc)
    nodeTd = DiagTree.addNode(diag, DiagTree.ADD, [[], []], Td; factor = Factor, para = :Td)

    addNode!(:Te, Te, para, diag, lver, rver, l, r, EX, EX, g0, gc; factor = Factor)
    if isempty(Te) == false
        nodeTe = Te[1]
    else
        nodeTe = 0
    end
    return @SVector [nodeTd, nodeTe]
end

function node4Ububble(para, diag, lver, rver, l, r, g0, gc)
    #     W[EX] = Lw[DIR] * Rw[DIR] * SPIN + Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
    #     W[DIR] = Lw[EX] * Rw[EX];
    Factor = SymFactor[U] / (2π)^para.dim
    Ud, Ue = [], []
    addNode!(:dxd, Ue, para, diag, lver, rver, l, r, DI, DI, g0, gc; factor = para.spin)
    addNode!(:dxe, Ue, para, diag, lver, rver, l, r, DI, EX, g0, gc)
    addNode!(:exd, Ue, para, diag, lver, rver, l, r, EX, DI, g0, gc)
    nodeUe = DiagTree.addNode(diag, DiagTree.ADD, [[], []], Ue; factor = Factor, para = :Ue)

    addNode!(:Ud, Ud, para, diag, lver, rver, l, r, EX, EX, g0, gc; factor = Factor)
    if isempty(Ud) == false
        nodeUd = Ud[1]
    else
        nodeUd = 0
    end
    return @SVector [nodeUd, nodeUe]
end

function node4Sbubble(para, diag, lver, rver, l, r, g0, gc)
    #     W[DIR] = Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
    #     W[EX] = Lw[DIR] * Rw[DIR] + Lw[EX] * Rw[EX];
    Factor = SymFactor[S] / (2π)^para.dim
    Sd, Se = [], []
    addNode!(:dxe, Sd, para, diag, lver, rver, l, r, DI, EX, g0, gc)
    addNode!(:exd, Sd, para, diag, lver, rver, l, r, EX, DI, g0, gc)
    nodeSd = DiagTree.addNode(diag, DiagTree.ADD, [[], []], Sd; factor = Factor, para = :Sd)

    addNode!(:dxd, Se, para, diag, lver, rver, l, r, DI, DI, g0, gc)
    addNode!(:exe, Se, para, diag, lver, rver, l, r, EX, EX, g0, gc)
    nodeSe = DiagTree.addNode(diag, DiagTree.ADD, [[], []], Se; factor = Factor, para = :Se)
    return @SVector [nodeSd, nodeSe]
end

# bubbletoDiagTree!(ver4Nodes, para, diag, ver4, b, legK, Kidx, Tidx, evalK, evalT, factor)
function bubbletoDiagTree!(ver4Nodes, para, diag, ver4, bubble, legK, Kidx::Int, factor = 1.0)
    KinL, KoutL, KinR, KoutR = legK[1], legK[2], legK[3], legK[4]
    @assert KinL + KinR ≈ KoutL + KoutR

    b = bubble
    c = b.chan
    G = ver4.G

    K = zero(KinL)
    K[Kidx] = 1
    Kt = KoutL + K - KinL
    Ku = KoutR + K - KinL
    Ks = KinL + KinR - K
    # println("KinL: ", KinL)
    # println("KinR: ", KinR)
    # println("Ks: ", Ks)

    Gorder = 0
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
    ver4toDiagTree(para, LLegK, Llopidx, Lver.Tidx, Lver.loopNum, factor, diag, Lver)
    ver4toDiagTree(para, RLegK, Rlopidx, Rver.Tidx, Rver.loopNum, factor, diag, Rver)

    rN = length(b.Rver.weight)
    for (l, Lw) in enumerate(b.Lver.weight)
        for (r, Rw) in enumerate(b.Rver.weight)

            map = b.map[(l-1)*rN+r]
            g0 = DiagTree.addPropagator(diag, 1, Gorder; site = G[1].Tpair[map.G0], loop = K, para = :G0)
            if c == T
                Kc = Kt
            elseif c == U
                Kc = Ku
            elseif c == S
                Kc = Ks
            end
            gc = DiagTree.addPropagator(diag, 1, Gorder; site = G[c].Tpair[map.Gx], loop = Kc; para = :Gx)

            if c == T
                Tde = node4Tbubble(para, diag, b.Lver, b.Rver, l, r, g0, gc)
                push!(ver4Nodes[map.ver], Tde)
            elseif c == U
                Ude = node4Ububble(para, diag, b.Lver, b.Rver, l, r, g0, gc)
                push!(ver4Nodes[map.ver], Ude)
            elseif c == S
                Sde = node4Sbubble(para, diag, b.Lver, b.Rver, l, r, g0, gc)
                push!(ver4Nodes[map.ver], Sde)
            else
                error("not implemented!")
            end
        end
    end
    return diag
end

function ver4toDiagTree(para::Para, legK, Kidx::Int, Tidx::Int, loopNum = para.loopNum, factor = 1.0,
    diag = _newDiag(para, legK, evalK), ver4 = Ver4{SVector{2,Int}}(para, loopNum, Tidx))


    KinL, KoutL, KinR, KoutR = legK[1], legK[2], legK[3], legK[4]

    # KoutR = KinL + KinR - KoutL
    @assert KinL + KinR ≈ KoutL + KoutR

    evaTpair(Tpair) = Tuple([evalT(t) for t in Tpair])

    qd = KinL - KoutL
    qe = KinR - KoutL
    Tidx = ver4.Tidx

    Vorder = 1
    if ver4.loopNum == 0
        qd, qe = [1, qd, evalK(qd)], [1, qe, evalK(qe)]
        if ver4.interactionTauNum == 2
            td = [Tidx, Tidx + 1]
            te = td
            vd = DiagTree.addPropagator(diag, 2, Vorder; site = td, loop = qd, para = :Vd)
            ve = DiagTree.addPropagator(diag, 2, Vorder; site = te, loop = qe, para = :Ve)
            wd = DiagTree.addPropagator(diag, 2, Vorder; site = td, loop = qd, para = :Wd)
            we = DiagTree.addPropagator(diag, 2, Vorder; site = te, loop = qe, para = :We)
            ver4.weight[1] = @SVector [vd, ve]
            ver4.weight[2] = @SVector [wd, 0]
            ver4.weight[3] = @SVector [0, we]
        elseif ver4.interactionTauNum == 1
            vd = DiagTree.addPropagator(diag, 2, Vorder; loop = qd, para = :Wd)
            ve = DiagTree.addPropagator(diag, 2, Vorder; loop = qe, para = :We)
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

    function split(weightList)
        dir, ex = [], []
        for w in weightList
            (w[DI] != 0) && push!(dir, w[DI])
            (w[EX] != 0) && push!(ex, w[EX])
        end
        return dir, ex
    end

    for i in 1:length(ver4.weight)
        ver4dNodes, ver4eNodes = split(ver4Nodes[i])
        nodeD = DiagTree.addNode(diag, DiagTree.ADD, [[], []], ver4dNodes; factor = factor, para = :dir)
        nodeE = DiagTree.addNode(diag, DiagTree.ADD, [[], []], ver4eNodes; factor = factor, para = :ex)
        ver4.weight[i] = @SVector [nodeD, nodeE]
    end
    dir, ex = split(ver4.weight)

    return diag, ver4, dir, ex
end

# function splitWeight(ver4)
#     dir = []
#     ex = []
#     for w in ver4.weight
#         push!(dir, w[DI])
#         push!(ex, w[DI])
#     end
#     return dir, ex
# end
