
function zeroLoopVer4Node!(nodes, diag, para, legK)
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

    function add!(name, type, _extT, _innerT)
        Vorder = 1
        poolName = symbol(name, type, "pool")

        vd = notProper(para, q[DI]) ? zero(Component) : DiagTree.addpropagator!(diag,
            poolName, Vorder, symbol(name, type, "di"); loop = q[DI], site = collect(_innerT[DI]))
        push!(nodes, Ver4identifier(name, type, DI, legK, _extT[DI], vd))

        ve = notProper(para, q[EX]) ? zero(Component) : DiagTree.addpropagator!(diag,
            poolName, Vorder, symbol(name, type, "ex"); loop = q[EX], site = collect(_innerT[EX]))
        push!(nodes, Ver4identifier(name, type, EX, legK, _extT[EX], ve))
    end

    for interaction in para.interaction
        name = interaction.name
        type = interaction.type

        if name == ChargeCharge || name == SpinSpin

            if Instant ∈ type && Dynamic ∉ type
                add!(name, Instant, extT_ins, innerT_ins)
            elseif Instant ∉ type && Dynamic ∈ type
                add!(name, Dynamic, extT_dyn, innerT_dyn)
            elseif Instant ∈ type && Dynamic ∈ type
                #if hasTau, instant interaction has an additional fake tau variable, making it similar to the dynamic interaction
                add!(name, Instant, extT_ins, innerT_dyn)
                add!(name, Dynamic, extT_dyn, innerT_dyn)
            end

            if D_Instant ∈ type && D_Dynamic ∉ type
                add!(name, D_Instant, extT_ins, innerT_ins)
            elseif D_Instant ∉ type && D_Dynamic ∈ type
                add!(name, D_Dynamic, extT_dyn, innerT_dyn)
            elseif D_Instant ∈ type && D_Dynamic ∈ type
                #if hasTau, instant interaction has an additional fake tau variable, making it similar to the dynamic interaction
                add!(name, D_Instant, extT_ins, innerT_dyn)
                add!(name, D_Dynamic, extT_dyn, innerT_dyn)
            end
        else
            @error("Interaction $name is not implemented!")
        end
    end

    return nodes
end

function legBasis(chan::Channel, legK, loopIdx)
    KinL, KoutL, KinR, KoutR = legK[1], legK[2], legK[3], legK[4]
    K = zero(KinL)
    K[loopIdx] = 1
    if chan == T
        Kx = KoutL + K - KinL
        LLegK = [KinL, KoutL, Kx, K]
        RLegK = [K, Kx, KinR, KoutR]
    elseif chan == U
        Kx = KoutR + K - KinL
        LLegK = [KinL, KoutR, Kx, K]
        RLegK = [K, Kx, KinR, KoutL]
    elseif chan == S
        Kx = KinL + KinR - K
        LLegK = [KinL, Kx, KinR, K]
        RLegK = [K, KoutL, Kx, KoutR]
    else
        @error("not implemented!")
    end
    return LLegK, K, RLegK, Kx
end


function mapBubbleT() end

