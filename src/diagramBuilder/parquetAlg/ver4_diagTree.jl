# addNode!(Td, para, diag, lver, rver, DI, EX)
isPropagator(ver) = (ver.para.innerLoopNum == 0)

function addNode!(nodes, diag, map, name::Symbol, lc, rc, g0, gc, factor = 1.0)
    ver4, Lver, Rver = map.v, map.l, map.r
    para = ver4.para
    extT = collect(ver4.Tpair[map.vidx])
    tauNum = ver4.para.interactionTauNum
    isW(idx) = (idx > 1) #for tauNum=2 case, each interaction has three weight, 0, 1, and 2; the first 0 is for instant, the 1 and 2 for dynamic

    lidx, ridx = map.lidx, map.ridx
    Lw, Rw = Lver.weight[lidx][lc], Rver.weight[ridx][rc]

    Gpool, child = [], []
    g0[2] ? push!(child, g0[1]) : push!(Gpool, g0[1])
    gc[2] ? push!(child, gc[1]) : push!(Gpool, gc[1])

    # println(g0)
    # println(gc)

    if tauNum == 1 || tauNum == 0
        Vpool = []
        # components : [GPool, VPool]
        isPropagator(Lver) ? push!(Vpool, Lw) : push!(child, Lw)
        isPropagator(Rver) ? push!(Vpool, Rw) : push!(child, Rw)
        if (Lw != 0 && Rw != 0)
            push!(nodes, DiagTree.addNodeByName!(diag, DiagTree.MUL, name, factor;
                Gpool = Gpool, Vpool = Vpool, child = child, para = extT))
        end
    elseif tauNum == 2
        Vpool, Wpool = [], []
        # components : [GPool, VPool, WPool]
        if isPropagator(Lver)
            isW(lidx) ? push!(Wpool, Lw) : push!(Vpool, Lw)
        else
            push!(child, Lw)
        end
        if isPropagator(Rver)
            isW(ridx) ? push!(Wpool, Rw) : push!(Vpool, Rw)
        else
            push!(child, Rw)
        end
        if (Lw != 0 && Rw != 0)
            push!(nodes, DiagTree.addNodeByName!(diag, DiagTree.MUL, name, factor;
                Gpool = Gpool, Vpool = Vpool, Wpool = Wpool, child = child, para = extT))
        end
    else
        error("tauNum = $tauNum has not yet been implemented!")
    end

    return nodes
end

# bubbletoDiagTree!(ver4Nodes, para, diag, ver4, b, legK, Kidx, Tidx, evalK, evalT, factor)
function bubbletoDiagTree!(diag, ver4, bubble, factor = 1.0)
    para = ver4.para

    b = bubble
    c = b.chan
    G = ver4.G
    Gorder = 0
    Lver, Rver = b.Lver, b.Rver

    ver4toDiagTree!(diag, Lver, factor)
    ver4toDiagTree!(diag, Rver, factor)

    for map in b.map
        Factor = b.factor
        extT = collect(map.v.Tpair[map.vidx])

        if isValidG(map.G0.para) == false || isValidG(map.Gx.para) == false
            map.node = @SVector [0, 0]
            continue
        end

        # g0 = DiagTree.addPropagator!(diag, :Gpool, Gorder, :G0; site = map.G0.Tpair, loop = map.G0.loopBasis)
        # g0Para = reconstruct(para, innerLoopNum = G0.loopNum, firstLoopIdx = G0.loopIdx, firstTauIdx = G0.Tspan[1])
        diag, g0idx, isnode0 = buildG(map.G0.para, map.G0.loopBasis, map.G0.Tpair; diag = diag)
        diag, gxidx, isnodex = buildG(map.Gx.para, map.Gx.loopBasis, map.Gx.Tpair; diag = diag)

        g0 = (g0idx, isnode0)
        gc = (gxidx, isnodex)

        if c == T
            # gc = DiagTree.addPropagator!(diag, :Gpool, Gorder, :Gt; site = map.Gx.Tpair, loop = map.Gx.loopBasis)

            Td, Te = [], []
            if removeBubble(map, c, DI, DI) == false
                addNode!(Td, diag, map, :dxd, DI, DI, g0, gc, para.spin)
            end
            addNode!(Td, diag, map, :dxe, DI, EX, g0, gc)
            addNode!(Td, diag, map, :exd, EX, DI, g0, gc)
            nodeTd = DiagTree.addNode!(diag, DiagTree.ADD, :Td, Factor; child = Td, para = extT)

            addNode!(Te, diag, map, :Te, EX, EX, g0, gc, Factor)
            nodeTe = isempty(Te) ? 0 : Te[1]
            map.node = @SVector [nodeTd, nodeTe]

        elseif c == U
            # gc = DiagTree.addPropagator!(diag, :Gpool, Gorder, :Gu; site = map.Gx.Tpair, loop = map.Gx.loopBasis)

            Ud, Ue = [], []
            if removeBubble(map, c, DI, DI) == false
                addNode!(Ue, diag, map, :dxd, DI, DI, g0, gc, para.spin)
            end
            addNode!(Ue, diag, map, :dxe, DI, EX, g0, gc)
            addNode!(Ue, diag, map, :exd, EX, DI, g0, gc)
            nodeUe = DiagTree.addNode!(diag, DiagTree.ADD, :Ue, Factor; child = Ue, para = extT)

            addNode!(Ud, diag, map, :Ud, EX, EX, g0, gc, Factor)
            nodeUd = isempty(Ud) ? 0 : Ud[1]
            map.node = @SVector [nodeUd, nodeUe]

        elseif c == S
            # gc = DiagTree.addPropagator!(diag, :Gpool, Gorder, :Gs; site = map.Gx.Tpair, loop = map.Gx.loopBasis)

            Sd, Se = [], []
            addNode!(Sd, diag, map, :dxe, DI, EX, g0, gc)
            addNode!(Sd, diag, map, :exd, EX, DI, g0, gc)
            nodeSd = DiagTree.addNode!(diag, DiagTree.ADD, :Sd, Factor; child = Sd, para = extT)

            addNode!(Se, diag, map, :dxd, DI, DI, g0, gc)
            addNode!(Se, diag, map, :exe, EX, EX, g0, gc)
            nodeSe = DiagTree.addNode!(diag, DiagTree.ADD, :Se, Factor; child = Se, para = extT)
            map.node = @SVector [nodeSd, nodeSe]
        else
            error("not implemented!")
        end

    end
end

function ver4toDiagTree!(diag, ver4, factor = 1.0)
    para = ver4.para

    legK = ver4.legK
    KinL, KoutL, KinR = legK[1], legK[2], legK[3]

    qd = KinL - KoutL
    qe = KinR - KoutL

    if ver4.para.innerLoopNum == 0
        Tidx = para.firstTauIdx
        Vorder = 1
        if para.interactionTauNum == 2
            td = [Tidx, Tidx + 1]
            te = td
            vd = notProper(para, qd) ? 0 : DiagTree.addPropagator!(diag, :Vpool, Vorder, :Vd; site = td, loop = qd)
            ve = notProper(para, qe) ? 0 : DiagTree.addPropagator!(diag, :Vpool, Vorder, :Ve; site = te, loop = qe)
            wd = notProper(para, qd) ? 0 : DiagTree.addPropagator!(diag, :Wpool, Vorder, :Wd; site = td, loop = qd)
            we = notProper(para, qe) ? 0 : DiagTree.addPropagator!(diag, :Wpool, Vorder, :We; site = te, loop = qe)
            ver4.weight[1] = @SVector [vd, ve]
            ver4.weight[2] = @SVector [wd, 0]
            ver4.weight[3] = @SVector [0, we]
            return diag, ver4, [vd, wd], [ve, we]
        elseif para.interactionTauNum == 1 || para.interactionTauNum == 0
            vd = notProper(para, qd) ? 0 : DiagTree.addPropagator!(diag, :Vpool, Vorder, :Vd; loop = qd)
            ve = notProper(para, qe) ? 0 : DiagTree.addPropagator!(diag, :Vpool, Vorder, :Ve; loop = qe)
            ver4.weight[1] = @SVector [vd, ve]
            return diag, ver4, [vd,], [ve,]
        else
            error("not implemented!")
        end
    end

    # LoopNum>=1
    for i in 1:length(ver4.weight)
        ver4.weight[i] = @SVector [0, 0]
    end

    for b in ver4.bubble
        bubbletoDiagTree!(diag, ver4, b, factor)
    end

    dir, ex = [], []
    for i in 1:length(ver4.weight)
        childdir, childex = [], []
        for map in ver4.child[i]
            w = map.node
            (w[DI] != 0) && push!(childdir, w[DI])
            (w[EX] != 0) && push!(childex, w[EX])
        end
        nodeD = DiagTree.addNode!(diag, DiagTree.ADD, :dir, factor; child = childdir, para = collect(ver4.Tpair[i]))
        nodeE = DiagTree.addNode!(diag, DiagTree.ADD, :ex, factor; child = childex, para = collect(ver4.Tpair[i]))
        ver4.weight[i] = @SVector [nodeD, nodeE]
        (nodeD != 0) && push!(dir, nodeD)
        (nodeE != 0) && push!(ex, nodeE)
    end
    return diag, ver4, dir, ex
end

"""
    function buildVer4(para, LegK, chan, F, V, All = union(F, V);
        Fouter = F, Vouter = V, Allouter = All, factor = 1.0, diag = newDiagTree(para, Tuple{Int,Int,Int,Int}, :Ver4))
    
    build DiagTree for the one-particle-irreducible 4-point vertex function using the parquet algorithm

# Arguments:
- weightType   : type of the weight of the propagators and the vertex functions
- para         : parameters to generate the diagram tree
- LegK         : momentum basis of external legs, only three of them are expected: [left in, left out, right in], the dimension of each legK is called loopBasis dimension.
"""
function buildVer4(para, LegK, chan, F, V, All = union(F, V);
    Fouter = F, Vouter = V, Allouter = All, factor = 1.0, diag = newDiagTree(para, :Ver4))

    ver4 = Ver4{SVector{2,Int}}(para, LegK, chan, F, V, All; Fouter = Fouter, Vouter = Vouter, Allouter = Allouter)
    return ver4toDiagTree!(diag, ver4, factor)
end
