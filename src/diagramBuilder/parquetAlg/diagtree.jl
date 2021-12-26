# struct ExternalSite
#     inL::Int
#     outL::Int
#     inR::Int
#     outR::Int
#     function ExternalSite(e::Tuple{Int,Int,Int,Int})
#         return new(e[1], e[2], e[3], e[4])
#     end
# end

################## Generate Expression Tree ########################
function _newDiag(para::Para, KbasisSample)
    # function LoopPool(name::Symbol, dim::Int, N::Int, type::DataType)
    loopBasisDim = length(KbasisSample)
    println("LoopBasis Dim derived from LegK: $loopBasisDim")
    Kpool = DiagTree.LoopPool(:K, para.dim, loopBasisDim, Float64)
    Tbasis = Tuple{Int,Int}
    # Tpool = DiagTree.uncachedPool(Tbasis)
    weightType = para.weightType
    if para.interactionTauNum == 2
        Gpool = DiagTree.propagatorPool(:Gpool, weightType)
        Vpool = DiagTree.propagatorPool(:Vpool, weightType)
        Wpool = DiagTree.propagatorPool(:Wpool, weightType)
        return DiagTree.Diagrams(Kpool, (Gpool, Vpool, Wpool), weightType, nodeParaType = Tuple{Int,Int,Int,Int})
    elseif para.interactionTauNum == 1
        Gpool = DiagTree.propagatorPool(:Gpool, weightType)
        Vpool = DiagTree.propagatorPool(:Vpool, weightType)
        return DiagTree.Diagrams(Kpool, (Gpool, Vpool), weightType, nodeParaType = Tuple{Int,Int,Int,Int})
    else
        error("not implemented!")
    end
end

# addNode!(Td, para, diag, lver, rver, DI, EX)
function addNode!(tauNum, name::Symbol, nodes, diag, Lver, Rver, map, lc, rc, g0, gc; factor = 1, para = Tuple{Int,Int,Int,Int}(zeros(0)))
    # function split(l, r, Lw, Rw)
    #     if l==1 && r==1
    #         return [], [Lw, Rw]
    #     elseif l==1 && r>1
    # end
    l, r = map.lv, map.rv
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
        push!(nodes, DiagTree.addNode(diag, DiagTree.MUL, components, child; factor = factor, name = name, para = para))
    end
    return nodes
end

function node4Tbubble(para, diag, ver4, lver, rver, g0, gc, map)
    #     W[DIR] = Lw[DIR] * Rw[DIR] * SPIN + Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
    #     W[EX] = Lw[EX] * Rw[EX];
    Factor = SymFactor[T] / (2π)^para.dim
    Td, Te = [], []
    extT = ver4.Tpair[map.ver]
    addNode!(para.interactionTauNum, :dxd, Td, diag, lver, rver, map, DI, DI, g0, gc; factor = para.spin, para = extT)
    addNode!(para.interactionTauNum, :dxe, Td, diag, lver, rver, map, DI, EX, g0, gc; para = extT)
    addNode!(para.interactionTauNum, :exd, Td, diag, lver, rver, map, EX, DI, g0, gc; para = extT)
    nodeTd = DiagTree.addNode(diag, DiagTree.ADD, [[], []], Td; factor = Factor, name = :Td, para = extT)

    addNode!(para.interactionTauNum, :Te, Te, diag, lver, rver, map, EX, EX, g0, gc; factor = Factor, para = extT)
    if isempty(Te) == false
        nodeTe = Te[1]
    else
        nodeTe = 0
    end
    return @SVector [nodeTd, nodeTe]
end

function node4Ububble(para, diag, ver4, lver, rver, g0, gc, map)
    #     W[EX] = Lw[DIR] * Rw[DIR] * SPIN + Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
    #     W[DIR] = Lw[EX] * Rw[EX];
    Factor = SymFactor[U] / (2π)^para.dim
    Ud, Ue = [], []
    extT = ver4.Tpair[map.ver]
    addNode!(para.interactionTauNum, :dxd, Ue, diag, lver, rver, map, DI, DI, g0, gc; factor = para.spin, para = extT)
    addNode!(para.interactionTauNum, :dxe, Ue, diag, lver, rver, map, DI, EX, g0, gc; para = extT)
    addNode!(para.interactionTauNum, :exd, Ue, diag, lver, rver, map, EX, DI, g0, gc; para = extT)
    nodeUe = DiagTree.addNode(diag, DiagTree.ADD, [[], []], Ue; factor = Factor, name = :Ue, para = extT)

    addNode!(para.interactionTauNum, :Ud, Ud, diag, lver, rver, map, EX, EX, g0, gc; factor = Factor, para = extT)
    if isempty(Ud) == false
        nodeUd = Ud[1]
    else
        nodeUd = 0
    end
    return @SVector [nodeUd, nodeUe]
end

function node4Sbubble(para, diag, ver4, lver, rver, g0, gc, map)
    #     W[DIR] = Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
    #     W[EX] = Lw[DIR] * Rw[DIR] + Lw[EX] * Rw[EX];
    Factor = SymFactor[S] / (2π)^para.dim
    Sd, Se = [], []
    extT = ver4.Tpair[map.ver]
    addNode!(para.interactionTauNum, :dxe, Sd, diag, lver, rver, map, DI, EX, g0, gc; para = extT)
    addNode!(para.interactionTauNum, :exd, Sd, diag, lver, rver, map, EX, DI, g0, gc; para = extT)
    nodeSd = DiagTree.addNode(diag, DiagTree.ADD, [[], []], Sd; factor = Factor, name = :Sd, para = extT)

    addNode!(para.interactionTauNum, :dxd, Se, diag, lver, rver, map, DI, DI, g0, gc; para = extT)
    addNode!(para.interactionTauNum, :exe, Se, diag, lver, rver, map, EX, EX, g0, gc; para = extT)
    nodeSe = DiagTree.addNode(diag, DiagTree.ADD, [[], []], Se; factor = Factor, name = :Se, para = extT)
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
            g0 = DiagTree.addPropagator(diag, 1, Gorder; site = G[1].Tpair[map.G0], loop = K, name = :G0)

            if c == T
                gc = DiagTree.addPropagator(diag, 1, Gorder; site = G[c].Tpair[map.Gx], loop = Kt, name = :Gt)
                Tde = node4Tbubble(para, diag, ver4, b.Lver, b.Rver, g0, gc, map)
                push!(ver4Nodes[map.ver], Tde)
            elseif c == U
                gc = DiagTree.addPropagator(diag, 1, Gorder; site = G[c].Tpair[map.Gx], loop = Ku, name = :Gu)
                Ude = node4Ububble(para, diag, ver4, b.Lver, b.Rver, g0, gc, map)
                push!(ver4Nodes[map.ver], Ude)
            elseif c == S
                gc = DiagTree.addPropagator(diag, 1, Gorder; site = G[c].Tpair[map.Gx], loop = Ks, name = :Gs)
                Sde = node4Sbubble(para, diag, ver4, b.Lver, b.Rver, g0, gc, map)
                push!(ver4Nodes[map.ver], Sde)
            else
                error("not implemented!")
            end
        end
    end
    return diag
end

function ver4toDiagTree(para::Para, legK, Kidx::Int, Tidx::Int, loopNum = para.loopNum, factor = 1.0,
    diag = _newDiag(para, legK[1]), ver4 = Ver4{SVector{2,Int}}(para, loopNum, Tidx))

    KinL, KoutL, KinR, KoutR = legK[1], legK[2], legK[3], legK[4]
    @assert KinL + KinR ≈ KoutL + KoutR

    qd = KinL - KoutL
    qe = KinR - KoutL
    Tidx = ver4.Tidx

    if ver4.loopNum == 0
        Vorder = 1
        if ver4.interactionTauNum == 2
            td = [Tidx, Tidx + 1]
            te = td
            vd = DiagTree.addPropagator(diag, :Vpool, Vorder; site = td, loop = qd, name = :Vd)
            ve = DiagTree.addPropagator(diag, :Vpool, Vorder; site = te, loop = qe, name = :Ve)
            wd = DiagTree.addPropagator(diag, :Wpool, Vorder; site = td, loop = qd, name = :Wd)
            we = DiagTree.addPropagator(diag, :Wpool, Vorder; site = te, loop = qe, name = :We)
            ver4.weight[1] = @SVector [vd, ve]
            ver4.weight[2] = @SVector [wd, 0]
            ver4.weight[3] = @SVector [0, we]
        elseif ver4.interactionTauNum == 1
            vd = DiagTree.addPropagator(diag, :Vpool, Vorder; loop = qd, name = :Vd)
            ve = DiagTree.addPropagator(diag, :Vpool, Vorder; loop = qe, name = :Ve)
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
    ver4Nodes = [[] for w in ver4.weight]

    for b in ver4.bubble
        bubbletoDiagTree!(ver4Nodes, para, diag, ver4, b, legK, Kidx, factor)
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
        nodeD = DiagTree.addNode(diag, DiagTree.ADD, [[], []], ver4dNodes; factor = factor, name = :dir, para = ver4.Tpair[i])
        nodeE = DiagTree.addNode(diag, DiagTree.ADD, [[], []], ver4eNodes; factor = factor, name = :ex, para = ver4.Tpair[i])
        ver4.weight[i] = @SVector [nodeD, nodeE]
    end
    dir, ex = split(ver4.weight)

    return diag, ver4, dir, ex
end
