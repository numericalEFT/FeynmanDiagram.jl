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

    function addnode!(interactionname, ver4name, type, _extT, _innerT, merge::Bool)
        Vorder = 1
        poolName = symbol(interactionname, type, "pool")

        vd, ve = [], []
        if notProper(para, q[DI]) == false
            sign = para.isFermi ? -1.0 : 1.0
            name = symbol(interactionname, type, "di")
            v = DiagTree.addpropagator!(diag, poolName, Vorder, name, sign; loop = q[DI], site = _innerT[DI])
            push!(vd, v)
        end
        if notProper(para, q[EX]) == false
            name = symbol(interactionname, type, "ex")
            v = DiagTree.addpropagator!(diag, poolName, Vorder, name, 1.0; loop = q[EX], site = _innerT[EX])
            push!(ve, v)
        end

        @assert isempty(vd) == false || isempty(ve) == false

        # "$_extT: external T of direct and exchange diagrams are different, impossible to merge!"
        # if DI and EX have the same external T, then it is possible to merge them into a same node
        if _extT[DI] == _extT[EX]
            id_diex = Vertex4(ver4name, type, BOTH, legK, _extT[DI], para)
            # println("add both")
            add!(nodes, id_diex, children = vcat(vd, ve))
        else
            id_di = Vertex4(ver4name, type, DI, legK, _extT[DI], para)
            id_ex = Vertex4(ver4name, type, EX, legK, _extT[EX], para)
            add!(nodes, id_di, children = [vd,])
            add!(nodes, id_ex, children = [ve,])
        end
    end

    function addUpDown!(interactionname, ver4name, type)
        @assert ver4name == UpUp || ver4name == UpDown

        if Instant ∈ type && Dynamic ∉ type
            addnode!(interactionname, ver4name, Instant, extT_ins, innerT_ins, true)
        elseif Instant ∉ type && Dynamic ∈ type
            addnode!(interactionname, ver4name, Dynamic, extT_dyn, innerT_dyn, false)
        elseif Instant ∈ type && Dynamic ∈ type
            #if hasTau, instant interaction has an additional fake tau variable, making it similar to the dynamic interaction
            addnode!(interactionname, ver4name, Instant, extT_ins, innerT_dyn, true)
            addnode!(interactionname, ver4name, Dynamic, extT_dyn, innerT_dyn, false)
        end

        if D_Instant ∈ type && D_Dynamic ∉ type
            addnode!(interactionname, ver4name, D_Instant, extT_ins, innerT_ins, true)
        elseif D_Instant ∉ type && D_Dynamic ∈ type
            addnode!(interactionname, ver4name, D_Dynamic, extT_dyn, innerT_dyn, false)
        elseif D_Instant ∈ type && D_Dynamic ∈ type
            #if hasTau, instant interaction has an additional fake tau variable, making it similar to the dynamic interaction
            addnode!(interactionname, ver4name, D_Instant, extT_ins, innerT_dyn, true)
            addnode!(interactionname, ver4name, D_Dynamic, extT_dyn, innerT_dyn, false)
        end
    end

    for interaction in para.interaction
        name = interaction.name
        type = interaction.type
        if name == UpUp || name == UpDown
            addUpDown!(name, name, type)
        elseif name == ChargeCharge
            addUpDown!(name, UpUp, type)
            addUpDown!(name, UpDown, type)
        else
            error("not implemented!")
        end
    end

    return nodes
end
