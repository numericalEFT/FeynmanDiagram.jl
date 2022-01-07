function bubbletoDiagTree!(diag, ver4, bubble)
    para = ver4.para

    b = bubble
    c = b.chan
    G = ver4.G
    Gorder = 0
    Lver, Rver = b.Lver, b.Rver

    ver4toDiagTree!(diag, Lver)
    ver4toDiagTree!(diag, Rver)

    for map in b.map
        factor = b.factor
        extT = collect(map.v.Tpair[map.vidx])
        Lw = map.l.weight[map.lidx]
        Rw = map.r.weight[map.ridx]

        # there is no need to add G0 or Gx and other nodes into diag if G0 and Gx don't exist
        if isValidG(map.G0.para) == false || isValidG(map.Gx.para) == false
            map.node = @MVector [zero(Component), zero(Component)]
            continue
        end

        # g0 = DiagTree.addpropagator!(diag, :Gpool, Gorder, :G0; site = map.G0.Tpair, loop = map.G0.loopBasis)
        # gc = DiagTree.addpropagator!(diag, :Gpool, Gorder, :Gx; site = map.Gx.Tpair, loop = map.Gx.loopBasis)

        diag, g0 = buildG(map.G0.para, map.G0.loopBasis, map.G0.Tpair; diag = diag)
        diag, gc = buildG(map.Gx.para, map.Gx.loopBasis, map.Gx.Tpair; diag = diag)

        if c == T || c == U
            if removeBubble(map, c, DI, DI)
                dd = zero(Component)
            else
                dd = DiagTree.addnode!(diag, DiagTree.MUL, Symbol("$(c)dd"), [g0, gc, Lw[DI], Rw[DI]], factor * para.spin; para = extT)
            end
        elseif c == S
            dd = DiagTree.addnode!(diag, DiagTree.MUL, Symbol("$(c)dd"), [g0, gc, Lw[DI], Rw[DI]], factor; para = extT)
        else
            error("not implemented!")
        end
        de = DiagTree.addnode!(diag, DiagTree.MUL, Symbol("$(c)de"), [g0, gc, Lw[DI], Rw[EX]], factor; para = extT)
        ed = DiagTree.addnode!(diag, DiagTree.MUL, Symbol("$(c)ed"), [g0, gc, Lw[EX], Rw[DI]], factor; para = extT)
        ee = DiagTree.addnode!(diag, DiagTree.MUL, Symbol("$(c)ee"), [g0, gc, Lw[EX], Rw[EX]], factor; para = extT)

        if c == T
            Td = DiagTree.addnode!(diag, DiagTree.ADD, :Td, [dd, de, ed]; para = extT)
            map.node = @MVector [Td, ee]
        elseif c == U
            Ue = DiagTree.addnode!(diag, DiagTree.ADD, :Ue, [dd, de, ed]; para = extT)
            map.node = @MVector [ee, Ue]

        elseif c == S
            Sd = DiagTree.addnode!(diag, DiagTree.ADD, :Sd, [de, ed]; para = extT)
            Se = DiagTree.addnode!(diag, DiagTree.ADD, :Se, [dd, ee]; para = extT)
            map.node = @MVector [Sd, Se]
        else
            error("not implemented!")
        end

    end
end

function ver4toDiagTree!(diag, ver4)
    para = ver4.para

    legK = ver4.legK
    KinL, KoutL, KinR = legK[1], legK[2], legK[3]

    qd = KinL - KoutL
    qe = KinR - KoutL

    for i in 1:length(ver4.weight)
        ver4.weight[i] = @MVector [zero(Component), zero(Component)]
    end

    if ver4.para.innerLoopNum == 0
        Tidx = para.firstTauIdx
        Vorder = 1
        if para.interactionTauNum == 2
            td = [Tidx, Tidx + 1]
            te = td
            vd = notProper(para, qd) ? zero(Component) : DiagTree.addpropagator!(diag, :Vpool, Vorder, :Vd; site = td, loop = qd)
            ve = notProper(para, qe) ? zero(Component) : DiagTree.addpropagator!(diag, :Vpool, Vorder, :Ve; site = te, loop = qe)
            wd = notProper(para, qd) ? zero(Component) : DiagTree.addpropagator!(diag, :Wpool, Vorder, :Wd; site = td, loop = qd)
            we = notProper(para, qe) ? zero(Component) : DiagTree.addpropagator!(diag, :Wpool, Vorder, :We; site = te, loop = qe)
            ver4.weight[1] = @MVector [vd, ve]
            ver4.weight[2][DI] = wd
            ver4.weight[3][EX] = we
            return diag, ver4, [vd, wd], [ve, we]
        elseif para.interactionTauNum == 1 || para.interactionTauNum == 0
            vd = notProper(para, qd) ? zero(Component) : DiagTree.addpropagator!(diag, :Vpool, Vorder, :Vd; loop = qd)
            ve = notProper(para, qe) ? zero(Component) : DiagTree.addpropagator!(diag, :Vpool, Vorder, :Ve; loop = qe)
            ver4.weight[1] = @MVector [vd, ve]
            return diag, ver4, [vd,], [ve,]
        else
            error("not implemented!")
        end
    end

    # LoopNum>=1
    for b in ver4.bubble
        bubbletoDiagTree!(diag, ver4, b)
    end

    dir, ex = [], []
    for i in 1:length(ver4.weight)
        childdir, childex = [], []
        for map in ver4.child[i]
            w = map.node
            (w[DI] != 0) && push!(childdir, w[DI])
            (w[EX] != 0) && push!(childex, w[EX])
        end
        nodeD = DiagTree.addnode!(diag, DiagTree.ADD, :dir, childdir; para = collect(ver4.Tpair[i]))
        nodeE = DiagTree.addnode!(diag, DiagTree.ADD, :ex, childex; para = collect(ver4.Tpair[i]))
        ver4.weight[i] = @MVector [nodeD, nodeE]
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
    Fouter = F, Vouter = V, Allouter = All, diag = newDiagTree(para, :Ver4))

    # ver4 = Ver4{SVector{2,DiagTree.Component}}(para, LegK, chan, F, V, All; Fouter = Fouter, Vouter = Vouter, Allouter = Allouter)
    ver4 = Ver4{MVector{2,Component}}(para, LegK, chan, F, V, All; Fouter = Fouter, Vouter = Vouter, Allouter = Allouter)
    return ver4toDiagTree!(diag, ver4)
end
