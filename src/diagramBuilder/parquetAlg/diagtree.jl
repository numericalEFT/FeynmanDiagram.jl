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
function _newDiag(weightType::DataType, para)
    # function LoopPool(name::Symbol, dim::Int, N::Int, type::DataType)
    loopBasisDim = para.internalLoopNum + para.externalLoopNum
    println("LoopBasis Dim derived from LegK: $loopBasisDim")
    Kpool = DiagTree.LoopPool(:K, para.loopDim, loopBasisDim, Float64)
    Tbasis = Tuple{Int,Int}
    # Tpool = DiagTree.uncachedPool(Tbasis)
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
function addNode!(nodes, tauNum, name::Symbol, diag, Lver, Rver, map, lc, rc, g0, gc, para = (0, 0, 0, 0), factor = 1.0)
    # function split(l, r, Lw, Rw)
    #     if l==1 && r==1
    #         return [], [Lw, Rw]
    #     elseif l==1 && r>1
    # end
    isNode(ver) = (ver.loopNum > 0)
    isW(idx) = (idx > 1) #for tauNum=2 case, each interaction has three weight, 0, 1, and 2; the first 0 is for instant, the 1 and 2 for dynamic

    l, r = map.lv, map.rv
    Lw, Rw = Lver.weight[l][lc], Rver.weight[r][rc]

    if tauNum == 1
        components, child = [[g0, gc], []], []
        isNode(Lver) ? push!(child, Lw) : push!(components[2], Lw)
        isNode(Rver) ? push!(child, Rw) : push!(components[2], Rw)
    else
        components, child = [[g0, gc], [], []], []
        # components : [GPool, VPool, WPool]
        isNode(Lver) ? push!(child, Lw) : (isW(l) ? push!(components[3], Lw) : push!(components[2], Lw))
        isNode(Rver) ? push!(child, Rw) : (isW(r) ? push!(components[3], Rw) : push!(components[2], Rw))
    end

    if (Lw != 0 && Rw != 0)
        push!(nodes, DiagTree.addNode!(diag, DiagTree.MUL, name, factor; propagator = components, child = child, para = para))
    end
    return nodes
end

function node4Tbubble(para, diag, ver4, lver, rver, g0, gc, map)
    #     W[DIR] = Lw[DIR] * Rw[DIR] * SPIN + Lw[DIR] * Rw[EX] + Lw[EX] * Rw[DIR];
    #     W[EX] = Lw[EX] * Rw[EX];
    Factor = SymFactor[T] / (2π)^para.loopDim
    Td, Te = [], []
    extT = ver4.Tpair[map.ver]
    addNode!(Td, para.interactionTauNum, :dxd, diag, lver, rver, map, DI, DI, g0, gc, extT, para.spin)
    addNode!(Td, para.interactionTauNum, :dxe, diag, lver, rver, map, DI, EX, g0, gc, extT)
    addNode!(Td, para.interactionTauNum, :exd, diag, lver, rver, map, EX, DI, g0, gc, extT)
    nodeTd = DiagTree.addNode!(diag, DiagTree.ADD, :Td, Factor; child = Td, para = extT)

    addNode!(Te, para.interactionTauNum, :Te, diag, lver, rver, map, EX, EX, g0, gc, extT, Factor)
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
    Factor = SymFactor[U] / (2π)^para.loopDim
    Ud, Ue = [], []
    extT = ver4.Tpair[map.ver]
    addNode!(Ue, para.interactionTauNum, :dxd, diag, lver, rver, map, DI, DI, g0, gc, extT, para.spin)
    addNode!(Ue, para.interactionTauNum, :dxe, diag, lver, rver, map, DI, EX, g0, gc, extT)
    addNode!(Ue, para.interactionTauNum, :exd, diag, lver, rver, map, EX, DI, g0, gc, extT)
    nodeUe = DiagTree.addNode!(diag, DiagTree.ADD, :Ue, Factor; child = Ue, para = extT)

    addNode!(Ud, para.interactionTauNum, :Ud, diag, lver, rver, map, EX, EX, g0, gc, extT, Factor)
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
    Factor = SymFactor[S] / (2π)^para.loopDim
    Sd, Se = [], []
    extT = ver4.Tpair[map.ver]
    addNode!(Sd, para.interactionTauNum, :dxe, diag, lver, rver, map, DI, EX, g0, gc, extT)
    addNode!(Sd, para.interactionTauNum, :exd, diag, lver, rver, map, EX, DI, g0, gc, extT)
    nodeSd = DiagTree.addNode!(diag, DiagTree.ADD, :Sd, Factor; child = Sd, para = extT)

    addNode!(Se, para.interactionTauNum, :dxd, diag, lver, rver, map, DI, DI, g0, gc, extT)
    addNode!(Se, para.interactionTauNum, :exe, diag, lver, rver, map, EX, EX, g0, gc, extT)
    nodeSe = DiagTree.addNode!(diag, DiagTree.ADD, :Se, Factor; child = Se, para = extT)
    return @SVector [nodeSd, nodeSe]
end

# bubbletoDiagTree!(ver4Nodes, para, diag, ver4, b, legK, Kidx, Tidx, evalK, evalT, factor)
function bubbletoDiagTree!(weightType::DataType, ver4Nodes, para, diag, ver4, bubble, legK, Kidx::Int, factor = 1.0)
    # KinL, KoutL, KinR, KoutR = legK[1], legK[2], legK[3], legK[4]
    # @assert KinL + KinR ≈ KoutL + KoutR
    KinL, KoutL, KinR = legK[1], legK[2], legK[3]
    KoutR = KinL + KinR - KoutL

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
    ver4toDiagTree(weightType, para, LLegK, Llopidx, Lver.Tidx, Lver.loopNum, factor, diag, Lver)
    ver4toDiagTree(weightType, para, RLegK, Rlopidx, Rver.Tidx, Rver.loopNum, factor, diag, Rver)

    rN = length(b.Rver.weight)
    for (l, Lw) in enumerate(b.Lver.weight)
        for (r, Rw) in enumerate(b.Rver.weight)

            map = b.map[(l-1)*rN+r]
            g0 = DiagTree.addPropagator!(diag, :Gpool, Gorder, :G0; site = G[1].Tpair[map.G0], loop = K)

            if c == T
                gc = DiagTree.addPropagator!(diag, :Gpool, Gorder, :Gt; site = G[c].Tpair[map.Gx], loop = Kt)
                Tde = node4Tbubble(para, diag, ver4, b.Lver, b.Rver, g0, gc, map)
                push!(ver4Nodes[map.ver], Tde)
            elseif c == U
                gc = DiagTree.addPropagator!(diag, :Gpool, Gorder, :Gu; site = G[c].Tpair[map.Gx], loop = Ku)
                Ude = node4Ububble(para, diag, ver4, b.Lver, b.Rver, g0, gc, map)
                push!(ver4Nodes[map.ver], Ude)
            elseif c == S
                gc = DiagTree.addPropagator!(diag, :Gpool, Gorder, :Gs; site = G[c].Tpair[map.Gx], loop = Ks)
                Sde = node4Sbubble(para, diag, ver4, b.Lver, b.Rver, g0, gc, map)
                push!(ver4Nodes[map.ver], Sde)
            else
                error("not implemented!")
            end
        end
    end
    return diag
end

function ver4toDiagTree(para, diag, ver4, legK, factor = 1.0)

    totalLoopNum(para) = para.totalLoopNum

    KinL, KoutL, KinR = legK[1], legK[2], legK[3]
    KoutR = KinL + KinR - KoutL

    qd = KinL - KoutL
    qe = KinR - KoutL

    if ver4.loopNum == 0
        Tidx = para.firstTauIdx + ver4.TidxOffset
        Vorder = 1
        if ver4.para.interactionTauNum == 2
            td = [Tidx, Tidx + 1]
            te = td
            vd = DiagTree.addPropagator!(diag, :Vpool, Vorder, :Vd; site = td, loop = qd)
            ve = DiagTree.addPropagator!(diag, :Vpool, Vorder, :Ve; site = te, loop = qe)
            wd = DiagTree.addPropagator!(diag, :Wpool, Vorder, :Wd; site = td, loop = qd)
            we = DiagTree.addPropagator!(diag, :Wpool, Vorder, :We; site = te, loop = qe)
            ver4.weight[1] = @SVector [vd, ve]
            ver4.weight[2] = @SVector [wd, 0]
            ver4.weight[3] = @SVector [0, we]
        elseif ver4.para.interactionTauNum == 1
            vd = DiagTree.addPropagator!(diag, :Vpool, Vorder, :Vd; loop = qd)
            ve = DiagTree.addPropagator!(diag, :Vpool, Vorder, :Ve; loop = qe)
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
        bubbletoDiagTree!(weightType, ver4Nodes, para, diag, ver4, b, legK, Kidx, factor)
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
        nodeD = DiagTree.addNode!(diag, DiagTree.ADD, :dir, factor; child = ver4dNodes, para = ver4.Tpair[i])
        nodeE = DiagTree.addNode!(diag, DiagTree.ADD, :ex, factor; child = ver4eNodes, para = ver4.Tpair[i])
        ver4.weight[i] = @SVector [nodeD, nodeE]
    end
    dir, ex = split(ver4.weight)

    return diag, ver4, dir, ex
end

"""
    function build(weightType::DataType, para::Para, LegK)
    
    build DiagTree for the one-particle-irreducible 4-point vertex function using the parquet algorithm

# Arguments:
- weightType   : type of the weight of the propagators and the vertex functions
- para         : parameters to generate the diagram tree
- LegK         : momentum basis of external legs, only three of them are expected: [left in, left out, right in], the dimension of each legK is called loopBasis dimension.
"""
# function build(weightType::DataType, para, LegK)
#     return ver4toDiagTree(weightType, para, LegK)
# end

"""
    function build(weightType::DataType, para::Para, LegK)
    
    build DiagTree for the one-particle-irreducible 4-point vertex function using the parquet algorithm

# Arguments:
- weightType   : type of the weight of the propagators and the vertex functions
- para         : parameters to generate the diagram tree
- LegK         : momentum basis of external legs, only three of them are expected: [left in, left out, right in], the dimension of each legK is called loopBasis dimension.
"""
function buildVer4(para, LegK, chan, F, V, All = union(F, V); Fouter = F, Vouter = V, Allouter = All, factor = 1.0)

    @assert length(LegK[1]) == length(LegK[2]) == length(LegK[3]) == para.totalLoopNum

    # @assert 0 <= para.internalLoopNum "internal LoopNum $loopNum is not in [0, $(para.internalLoopNum)]"
    # if loopNum > 0
    #     @assert Kidx <= totalLoopNum(para) "Kidx $Kidx can't be larger than total loop number $(totalLoopNum(para))"
    # end

    weightType = para.weightType
    Kpool = DiagTree.LoopPool(:K, para.loopDim, para.totalLoopNum, Float64)

    if para.interactionTauNum == 2
        Gpool = DiagTree.propagatorPool(:Gpool, weightType)
        Vpool = DiagTree.propagatorPool(:Vpool, weightType)
        Wpool = DiagTree.propagatorPool(:Wpool, weightType)
        diag = DiagTree.Diagrams(Kpool, (Gpool, Vpool, Wpool), weightType, nodeParaType = Tuple{Int,Int,Int,Int})
    elseif para.interactionTauNum == 1
        Gpool = DiagTree.propagatorPool(:Gpool, weightType)
        Vpool = DiagTree.propagatorPool(:Vpool, weightType)
        diag = DiagTree.Diagrams(Kpool, (Gpool, Vpool), weightType, nodeParaType = Tuple{Int,Int,Int,Int})
    else
        error("not implemented!")
    end

    ver4 = Ver4{SVector{2,Int}}(para, chan, F, V, All; Fouter = Fouter, Vouter = Vouter, Allouter = Allouter)
    return ver4toDiagTree(para, diag, ver4, LegK, factor)
end
