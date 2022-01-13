function bareVer4!(nodes, diag, para, legK)
    KinL, KoutL, KinR = legK[1], legK[2], legK[3]
    t0 = para.firstTauIdx

    q = [KinL - KoutL, KinR - KoutL]

    if para.hasTau
        extT_ins = [(t0, t0, t0, t0), (t0, t0, t0, t0)]
        extT_dyn = [(t0, t0, t0 + 1, t0 + 1), (t0, t0 + 1, t0 + 1, t0)]
        innerT_ins = [(0, 0), (0, 0)]
        innerT_dyn = [(t0, t0 + 1), (t0, t0 + 1)]
    else
        extT_ins = [(0, 0, 0, 0), (0, 0, 0, 0)]
        extT_dyn = extT_ins
        innerT_ins = [(0, 0), (0, 0)]
        innerT_dyn = innerT_ins
    end

    function addnode!(name, type, _extT, _innerT, merge::Bool)
        Vorder = 1
        poolName = symbol(name, type, "pool")


        # add!(nodes, ; children = [])

        vd, ve = [], []
        if notProper(para, q[DI]) == false
            sign = para.isFermi ? -1.0 : 1.0
            v = DiagTree.addpropagator!(diag, poolName, Vorder, symbol(name, type, "di"), sign; loop = q[DI], site = _innerT[DI])
            push!(vd, v)
        end
        if notProper(para, q[EX]) == false
            v = DiagTree.addpropagator!(diag, poolName, Vorder, symbol(name, type, "ex"), 1.0; loop = q[EX], site = _innerT[EX])
            push!(ve, v)
        end

        @assert isempty(vd) == false || isempty(ve) == false

        # "$_extT: external T of direct and exchange diagrams are different, impossible to merge!"
        # if DI and EX have the same external T, then it is possible to merge them into a same node
        if _extT[DI] == _extT[EX]
            id_diex = Vertex4(name, type, BOTH, legK, _extT[DI], para)
            # println("add both")
            add!(nodes, id_diex, children = vcat(vd, ve))
        else
            id_di = Vertex4(name, type, DI, legK, _extT[DI], para)
            id_ex = Vertex4(name, type, EX, legK, _extT[EX], para)
            add!(nodes, id_di, children = [vd,])
            add!(nodes, id_ex, children = [ve,])
        end
    end

    for interaction in para.interaction
        name = interaction.name
        type = interaction.type
        @assert name == UpUp || name == UpDown

        if Instant ∈ type && Dynamic ∉ type
            addnode!(name, Instant, extT_ins, innerT_ins, true)
        elseif Instant ∉ type && Dynamic ∈ type
            addnode!(name, Dynamic, extT_dyn, innerT_dyn, false)
        elseif Instant ∈ type && Dynamic ∈ type
            #if hasTau, instant interaction has an additional fake tau variable, making it similar to the dynamic interaction
            addnode!(name, Instant, extT_ins, innerT_dyn, true)
            addnode!(name, Dynamic, extT_dyn, innerT_dyn, false)
        end

        if D_Instant ∈ type && D_Dynamic ∉ type
            addnode!(name, D_Instant, extT_ins, innerT_ins, true)
        elseif D_Instant ∉ type && D_Dynamic ∈ type
            addnode!(name, D_Dynamic, extT_dyn, innerT_dyn, false)
        elseif D_Instant ∈ type && D_Dynamic ∈ type
            #if hasTau, instant interaction has an additional fake tau variable, making it similar to the dynamic interaction
            addnode!(name, D_Instant, extT_ins, innerT_dyn, true)
            addnode!(name, D_Dynamic, extT_dyn, innerT_dyn, false)
        end
    end

    return nodes
end
