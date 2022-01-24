function bareVer4!(diag, para::GenericPara, legK, diex = [DI, EX], nodeName = nothing)
    @assert para.diagType == Ver4Diag

    KinL, KoutL, KinR = legK[1], legK[2], legK[3]
    t0 = para.firstTauIdx

    nodes = Node{Vertex4}[]

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

    function addbare!(response::Response, type::AnalyticProperty, _diex::Int, _innerT, _q)
        @assert _diex == DI || _diex == EX
        Vorder = 1
        poolName = symbol(response, type, "pool")
        name = (_diex == DI) ?
               symbol(response, type, "D") :
               symbol(response, type, "E")

        sign = (_diex == DI) ? -1.0 : 1.0
        if para.isFermi == false
            sign = abs(sign)
        end

        if notProper(para, _q) == false && _diex in diex
            #create new bare ver4 only if _diex is required in the diex table 
            return DiagTree.addpropagator!(diag, poolName, Vorder, name, sign; loop = _q, site = _innerT)
        else
            return nothing
        end
    end

    function addver4!(response::Response, type, _extT, vd, ve)
        # "$_extT: external T of direct and exchange diagrams are different, impossible to merge!"
        # if DI and EX have the same external T, then it is possible to merge them into a same node
        if isnothing(vd) && isnothing(ve)
            return
        end
        if _extT[DI] == _extT[EX] && isnothing(vd) == false && isnothing(ve) == false
            id_diex = Vertex4(response, type, BOTH, legK, _extT[DI], para)
            add!(nodes, id_diex, children = [vd, ve])
        else
            id_di = Vertex4(response, type, DI, legK, _extT[DI], para)
            id_ex = Vertex4(response, type, EX, legK, _extT[EX], para)
            (isnothing(vd) == false) && add!(nodes, id_di, children = [vd,])
            (isnothing(ve) == false) && add!(nodes, id_ex, children = [ve,])
        end
    end

    function addresponse!(response::Response, type, _extT, _innerT)
        if response == UpUp || response == UpDown
            vd = addbare!(response, type, DI, _innerT[DI], q[DI])
            ve = addbare!(response, type, EX, _innerT[EX], q[EX])
            addver4!(response, type, _extT, vd, ve)
        elseif response == ChargeCharge
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
        response = interaction.response
        typeVec = interaction.type

        if Instant ∈ typeVec && Dynamic ∉ typeVec
            addresponse!(response, Instant, extT_ins, innerT_ins)
        elseif Instant ∉ typeVec && Dynamic ∈ typeVec
            addresponse!(response, Dynamic, extT_dyn, innerT_dyn)
        elseif Instant ∈ typeVec && Dynamic ∈ typeVec
            #if hasTau, instant interaction has an additional fake tau variable, making it similar to the dynamic interaction
            addresponse!(response, Instant, extT_ins, innerT_dyn)
            addresponse!(response, Dynamic, extT_dyn, innerT_dyn)
        end

        if D_Instant ∈ typeVec && D_Dynamic ∉ typeVec
            addresponse!(response, D_Instant, extT_ins, innerT_ins)
        elseif D_Instant ∉ typeVec && D_Dynamic ∈ typeVec
            addresponse!(response, D_Dynamic, extT_dyn, innerT_dyn)
        elseif D_Instant ∈ typeVec && D_Dynamic ∈ typeVec
            #if hasTau, instant interaction has an additional fake tau variable, making it similar to the dynamic interaction
            addresponse!(response, D_Instant, extT_ins, innerT_dyn)
            addresponse!(response, D_Dynamic, extT_dyn, innerT_dyn)
        end
    end

    for n in nodes
        if isnothing(nodeName)
            nodeName = symbol(n.id.response, n.id.type, "Int")
        end
        generate_node_from_children!(diag, n, ADD, 1.0, nodeName; para = n.id)
    end

    return nodes
end

function bareVer4(para::GenericPara, legK, diex::Vector{Permutation} = [Di, Ex])
    # @assert para.diagType == Ver4Diag

    KinL, KoutL, KinR = legK[1], legK[2], legK[3]
    t0 = para.firstTauIdx

    nodes = Diagram{para.weightType}[]

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

    function bare(response::Response, type::AnalyticProperty, _diex::Permutation, _innerT, _q)
        @assert _diex == Di || _diex == Ex

        sign = (para.isFermi && (_diex == Di)) ? -1.0 : 1.0

        if notProper(para, _q) == false && _diex in diex
            #create new bare ver4 only if _diex is required in the diex table 
            vid = InteractionId(para, response, type, k = _q, t = _innerT, permu = _diex)
            return Diagram(vid, factor = sign)
        else
            return nothing
        end
    end

    function addver4!(response::Response, type, _extT, vd, ve)
        id_di = Ver4Id(para, response, type, k = legK, t = _extT[DI])
        (isnothing(vd) == false) && push!(nodes, Diagram(id_di, Sum(), [vd,]))

        id_ex = Ver4Id(para, response, type, k = legK, t = _extT[EX])
        (isnothing(ve) == false) && push!(nodes, Diagram(id_ex, Sum(), [ve,]))
        # end
    end

    function addresponse!(response::Response, type, _extT, _innerT)
        if response == UpUp
            vd = bare(response, type, Di, _innerT[DI], q[DI])
            ve = bare(response, type, Ex, _innerT[EX], q[EX])
            addver4!(UpUp, type, _extT, vd, ve)
        elseif response == UpDown
            vd = bare(UpDown, type, Di, _innerT[DI], q[DI])
            ve = nothing
            addver4!(UpDown, type, _extT, vd, ve)
        elseif response == ChargeCharge
            # UpUp channel
            vuud = bare(ChargeCharge, type, Di, _innerT[DI], q[DI])
            vuue = bare(ChargeCharge, type, Ex, _innerT[EX], q[EX])
            addver4!(UpUp, type, _extT, vuud, vuue)

            # UpDown channel
            vupd = bare(ChargeCharge, type, Di, _innerT[DI], q[DI])
            vupe = nothing
            # UpDown, exchange channel doesn't exist for the charge-charge interaction
            addver4!(UpDown, type, _extT, vupd, vupe)
        else
            error("not implemented!")
        end
    end

    for interaction in para.interaction
        response = interaction.response
        typeVec = interaction.type

        if Instant ∈ typeVec && Dynamic ∉ typeVec
            addresponse!(response, Instant, extT_ins, innerT_ins)
        elseif Instant ∉ typeVec && Dynamic ∈ typeVec
            addresponse!(response, Dynamic, extT_dyn, innerT_dyn)
        elseif Instant ∈ typeVec && Dynamic ∈ typeVec
            #if hasTau, instant interaction has an additional fake tau variable, making it similar to the dynamic interaction
            addresponse!(response, Instant, extT_ins, innerT_dyn)
            addresponse!(response, Dynamic, extT_dyn, innerT_dyn)
        end

        if D_Instant ∈ typeVec && D_Dynamic ∉ typeVec
            addresponse!(response, D_Instant, extT_ins, innerT_ins)
        elseif D_Instant ∉ typeVec && D_Dynamic ∈ typeVec
            addresponse!(response, D_Dynamic, extT_dyn, innerT_dyn)
        elseif D_Instant ∈ typeVec && D_Dynamic ∈ typeVec
            #if hasTau, instant interaction has an additional fake tau variable, making it similar to the dynamic interaction
            addresponse!(response, D_Instant, extT_ins, innerT_dyn)
            addresponse!(response, D_Dynamic, extT_dyn, innerT_dyn)
        end
    end

    return nodes
end
