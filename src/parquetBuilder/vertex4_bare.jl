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

    function addbare!(responsename, type, diex, _innerT, _q)
        @assert diex == DI || diex == EX
        Vorder = 1
        poolName = symbol(responsename, type, "pool")
        name = (diex == DI) ?
               symbol(responsename, type, "D") :
               symbol(responsename, type, "E")

        sign = (diex == DI) ? -1.0 : 1.0
        if para.isFermi == false
            sign = abs(sign)
        end

        if notProper(para, _q) == false
            return DiagTree.addpropagator!(diag, poolName, Vorder, name, sign; loop = _q, site = _innerT)
        else
            return nothing
        end
    end

    function addver4!(responsename, type, _extT, vd, ve)
        # "$_extT: external T of direct and exchange diagrams are different, impossible to merge!"
        # if DI and EX have the same external T, then it is possible to merge them into a same node
        if isnothing(vd) && isnothing(ve)
            return
        end
        if _extT[DI] == _extT[EX] && isnothing(vd) == false && isnothing(ve) == false
            id_diex = Vertex4(responsename, type, BOTH, legK, _extT[DI], para)
            add!(nodes, id_diex, children = [vd, ve])
        else
            id_di = Vertex4(responsename, type, DI, legK, _extT[DI], para)
            id_ex = Vertex4(responsename, type, EX, legK, _extT[EX], para)
            (isnothing(vd) == false) && add!(nodes, id_di, children = [vd,])
            (isnothing(ve) == false) && add!(nodes, id_ex, children = [ve,])
        end
    end

    function addresponse!(responsename, type, _extT, _innerT)
        if responsename == UpUp || responsename == UpDown
            vd = addbare!(responsename, type, DI, _innerT[DI], q[DI])
            ve = addbare!(responsename, type, EX, _innerT[EX], q[EX])
            addver4!(responsename, type, _extT, vd, ve)
        elseif responsename == ChargeCharge
            # UpUp channel
            vuud = addbare!(ChargeCharge, type, DI, _innerT[DI], q[DI])
            vuue = addbare!(ChargeCharge, type, EX, _innerT[EX], q[EX])
            addver4!(UpUp, type, _extT, vuud, vuue)

            # UpDown channel
            vupd = addbare!(ChargeCharge, type, DI, _innerT[DI], q[DI])
            vupe = nothing
            # UpDown, exchange channel doesn't exist for the charge-charge interaction
            addver4!(UpDown, type, _extT, vupd, vupe)
        else
            error("not implemented!")
        end
    end

    for interaction in para.interaction
        responsename = interaction.name
        typeVec = interaction.type

        if Instant ∈ typeVec && Dynamic ∉ typeVec
            addresponse!(responsename, Instant, extT_ins, innerT_ins)
        elseif Instant ∉ typeVec && Dynamic ∈ typeVec
            addresponse!(responsename, Dynamic, extT_dyn, innerT_dyn)
        elseif Instant ∈ typeVec && Dynamic ∈ typeVec
            #if hasTau, instant interaction has an additional fake tau variable, making it similar to the dynamic interaction
            addresponse!(responsename, Instant, extT_ins, innerT_dyn)
            addresponse!(responsename, Dynamic, extT_dyn, innerT_dyn)
        end

        if D_Instant ∈ typeVec && D_Dynamic ∉ typeVec
            addresponse!(responsename, D_Instant, extT_ins, innerT_ins)
        elseif D_Instant ∉ typeVec && D_Dynamic ∈ typeVec
            addresponse!(responsename, D_Dynamic, extT_dyn, innerT_dyn)
        elseif D_Instant ∈ typeVec && D_Dynamic ∈ typeVec
            #if hasTau, instant interaction has an additional fake tau variable, making it similar to the dynamic interaction
            addresponse!(responsename, D_Instant, extT_ins, innerT_dyn)
            addresponse!(responsename, D_Dynamic, extT_dyn, innerT_dyn)
        end
    end

    return nodes
end
