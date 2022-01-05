mutable struct Ver4identifier
    ######### properties that defines a unique ver4 ###################
    name::InteractionName #composite, chargecharge, spinspin, ...
    type::InteractionType
    DiEx::Int # 1 for direct, 2 for exchange
    legK::Vector{Vector{Float64}}
    extT::Tuple{Int,Int,Int,Int} #all possible extT from different interactionType
    node::Component
end


function Base.isequal(a::Ver4identifier, b::Ver4identifier)
    # only parent is allowed to be different
    if (a.name != b.name) || (a.type != b.type) || (a.extT != b.extT) || (a.DiEx != b.DiEx) || (a.legK != b.legK)
        return false
    else
        if a.node != b.node
            @warn("Two identical ver4 has different Component: \n$a and $b")
        end
        return true
    end
end
Base.:(==)(a::Ver4identifier, b::Ver4identifier) = Base.isequal(a, b)

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

function mapBubbleT() end

