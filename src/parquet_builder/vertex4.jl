"""
    function Ver4{W}(para::Para, loopNum = para.internalLoopNum, tidx = 1; chan = para.chan, F = para.F, V = para.V, level = 1, id = [1,]) where {W}

    Generate 4-vertex diagrams using Parquet Algorithm

#Arguments
- `para`: parameters. It should provide internalLoopNum, interactionTauNum, firstTauIdx
- `chan`: list of channels of the current 4-vertex. 
- `F`   : channels of left sub-vertex for the particle-hole and particle-hole-exchange bubbles
- `V`   : channels of left sub-vertex for the particle-particle bubble
- `All`   : channels of right sub-vertex of all channels
- `Fouter`   : channels of left sub-vertex for the particle-hole and particle-hole-exchange bubbles, only take effect for the outermost bubble
- `Vouter`   : channels of left sub-vertex for the particle-particle bubble, only take effect for the outermost bubble
- `Allouter`   : channels of right sub-vertex of all channels
- `loopNum`: momentum loop degrees of freedom of the 4-vertex diagrams
- `tidx`: the first τ variable index. It will be the τ variable of the left incoming electron for all 4-vertex diagrams
- `level`: level in the diagram tree
- `id`: the first element will be used as the id of the Ver4. All nodes in the tree will be labeled in preorder depth-first search
"""
function vertex4(para::GenericPara, legK, chan::AbstractVector, subdiagram = false; level = 1,
    phi_toplevel = para.extra.phi, ppi_toplevel = para.extra.ppi, Γ4_toplevel = para.extra.Γ4, name = :none)

    (subdiagram == false) && uidreset()

    @assert para.extra isa ParquetBlocks
    @assert para.totalTauNum >= maxVer4TauIdx(para) "Increase totalTauNum!\n$para"
    @assert para.totalLoopNum >= maxVer4LoopIdx(para) "Increase totalLoopNum\n$para"

    phi, ppi = para.extra.phi, para.extra.ppi

    @assert (PHr in phi) == false "PHi vertex is particle-hole irreducible, so that PHr channel is not allowed in $phi"
    @assert (PPr in ppi) == false "PPi vertex is particle-particle irreducible, so that PPr channel is not allowed in $ppi"
    @assert (PHr in phi_toplevel) == false "PHi vertex is particle-hole irreducible, so that PHr channel is not allowed in $phi_toplevel"
    @assert (PPr in ppi_toplevel) == false "PPi vertex is particle-particle irreducible, so that PPr channel is not allowed in $ppi_toplevel"

    @assert length(legK[1]) == length(legK[2]) == length(legK[3]) == para.totalLoopNum

    KinL, KoutL, KinR = legK[1], legK[2], legK[3]
    KoutR = (length(legK) > 3) ? legK[4] : KinL + KinR - KoutL
    @assert KoutR ≈ KinL + KinR - KoutL
    legK = [KinL, KoutL, KinR, KoutR]


    loopNum = para.innerLoopNum
    @assert loopNum >= 0

    diags = Diagram{para.weightType}[]

    if loopNum == 0
        append!(diags, bareVer4(para, legK, [Di, Ex]))
    else # loopNum>0
        for c in chan
            if c == Alli
                continue
            end

            partition = orderedPartition(loopNum - 1, 4, 0)

            for p in partition

                if c == PHr || c == PHEr || c == PPr
                    bub = bubble(para, legK, c, p, level, name, phi_toplevel, ppi_toplevel, Γ4_toplevel)
                    # println(bub)
                    append!(diags, bub)
                end
            end
        end
        # # TODO: add envolpe diagrams
    end
    @assert all(x -> x.id isa Ver4Id, diags) "not all id are Ver4Id! $diags"
    @assert all(x -> x.id.extK ≈ legK, diags) "not all extK are the same! $diags"

    # @assert isempty(diags) == false "got empty ver4! $chan with\n $para\n"
    if isempty(diags)
        return DataFrame(response = [], type = [], extT = [], diagram = [])
    end

    df = toDataFrame(diags, expand = true)
    # println(df[:, [:response, :type, :extT, :diagram]])
    groups = mergeby(df, [:response, :type, :extT], name = name,
        getid = g -> Ver4Id(para, g[1, :response], g[1, :type], k = legK, t = g[1, :extT]) #generate id from the dataframe
    )
    return groups
end

function bubble(para::GenericPara, legK, chan::TwoBodyChannel, partition::Vector{Int}, level::Int, name::Symbol,
    phi_toplevel, ppi_toplevel, Γ4_toplevel)

    diag = Diagram{para.weightType}[]

    TauNum = para.interactionTauNum # maximum tau number for each bare interaction
    oL, oG0, oR, oGx = partition[1], partition[2], partition[3], partition[4]
    if isValidG(para.filter, oG0) == false || isValidG(para.filter, oGx) == false
        return diag
    end

    #the first loop idx is the inner loop of the bubble!
    LoopIdx = para.firstLoopIdx
    idx, maxLoop = findFirstLoopIdx(partition, LoopIdx + 1)
    LfirstLoopIdx, G0firstLoopIdx, RfirstLoopIdx, GxfirstLoopIdx = idx
    @assert maxLoop == maxVer4LoopIdx(para)

    diagType = [Ver4Diag, GreenDiag, Ver4Diag, GreenDiag]
    idx, maxTau = findFirstTauIdx(partition, diagType, para.firstTauIdx, TauNum)
    LfirstTauIdx, G0firstTauIdx, RfirstTauIdx, GxfirstTauIdx = idx
    @assert maxTau == maxVer4TauIdx(para) "Partition $partition with tauNum configuration $idx. maxTau = $maxTau, yet $(maxTauIdx(para)) is expected!"

    lPara = reconstruct(para, innerLoopNum = oL, firstLoopIdx = LfirstLoopIdx, firstTauIdx = LfirstTauIdx)
    rPara = reconstruct(para, innerLoopNum = oR, firstLoopIdx = RfirstLoopIdx, firstTauIdx = RfirstTauIdx)
    gxPara = reconstruct(para, innerLoopNum = oGx, firstLoopIdx = GxfirstLoopIdx, firstTauIdx = GxfirstTauIdx)
    g0Para = reconstruct(para, innerLoopNum = oG0, firstLoopIdx = G0firstLoopIdx, firstTauIdx = G0firstTauIdx)

    phi, ppi, Γ4 = para.extra.phi, para.extra.ppi, para.extra.Γ4
    if chan == PHr || chan == PHEr
        Γi = (level == 1) ? phi_toplevel : phi
        Γf = (level == 1) ? Γ4_toplevel : Γ4
    elseif chan == PPr
        Γi = (level == 1) ? ppi_toplevel : ppi
        Γf = (level == 1) ? Γ4_toplevel : Γ4
    else
        error("chan $chan isn't implemented!")
    end

    LLegK, K, RLegK, Kx = legBasis(chan, legK, LoopIdx)
    # println(K, ", ", Kx)

    Lver = vertex4(lPara, LLegK, Γi, true; level = level + 1, name = :Γi)
    isempty(Lver) && return diag
    # println("Γf: ", Γf)
    Rver = vertex4(rPara, RLegK, Γf, true; level = level + 1, name = :Γf)
    isempty(Rver) && return diag

    for ldiag in Lver.diagram
        for rdiag in Rver.diagram
            extT, G0T, GxT = tauBasis(chan, ldiag.id.extT, rdiag.id.extT)
            # diag, g0 = buildG(bubble.g0, K, (LvT[OUTR], RvT[INL]); diag = diag)
            # diag, gc = buildG(bubble.gx, Kx, (RvT[OUTL], LvT[INR]); diag = diag)
            # g0 = DiagTree.addpropagator!(diag, :Gpool, 0, :G0; site = G0T, loop = K)
            # gc = DiagTree.addpropagator!(diag, :Gpool, 0, :Gx; site = GxT, loop = Kx)
            g0 = Diagram(GreenId(g0Para, k = K, t = G0T), name = :G0)
            gx = Diagram(GreenId(gxPara, k = Kx, t = GxT), name = :Gx)
            append!(diag, bubble2diag(para, chan, ldiag, rdiag, legK, g0, gx))
        end
    end
    return diag
end

function bubble2diag(para, chan, ldiag, rdiag, extK, g0, gx)
    lid, rid = ldiag.id, rdiag.id
    ln, rn = lid.response, rid.response
    lo, ro = lid.para.innerLoopNum, rid.para.innerLoopNum
    vtype = typeMap(lid.type, rid.type)

    extT, G0T, GxT = tauBasis(chan, lid.extT, rid.extT)
    Factor = factor(para, chan)
    spin(response) = (response == UpUp ? "↑↑" : "↑↓")

    diag = Diagram{para.weightType}[]

    function add(Lresponse::Response, Rresponse::Response, Vresponse::Response, factor = 1.0)
        if ln == Lresponse && rn == Rresponse
            nodeName = Symbol("$(spin(Lresponse))x$(spin(Rresponse)) → $chan,")
            id = Ver4Id(para, Vresponse, vtype, k = extK, t = extT, chan = chan)
            push!(diag, Diagram(id, Prod(), [g0, gx, ldiag, rdiag], factor = factor * Factor, name = nodeName))
        end
    end

    if chan == PHr
        add(UpUp, UpUp, UpUp, 1.0)
        add(UpDown, UpDown, UpUp, 1.0)
        add(UpUp, UpDown, UpDown, 1.0)
        add(UpDown, UpUp, UpDown, 1.0)
    elseif chan == PHEr
        add(UpUp, UpUp, UpUp, 1.0)
        add(UpUp, UpUp, UpDown, 1.0)
        add(UpDown, UpDown, UpUp, 1.0)
        add(UpDown, UpDown, UpDown, 1.0)
        #! the sign here is from the spin symmetry, not from the fermionic statistics
        add(UpUp, UpDown, UpDown, -1.0)
        #! the sign here is from the spin symmetry, not from the fermionic statistics
        add(UpDown, UpUp, UpDown, -1.0)
    elseif chan == PPr
        add(UpUp, UpUp, UpUp, 1.0)
        #! the sign here is from the spin symmetry, not from the fermionic statistics
        add(UpDown, UpDown, UpDown, -2.0)
        add(UpUp, UpDown, UpDown, 1.0)
        add(UpDown, UpUp, UpDown, 1.0)
    else
        error("chan $chan isn't implemented!")
    end

    return diag
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

        # there is an overall sign coming from Taylor expansion of exp(-S) depsite the statistics
        if _diex == Di
            sign = -1.0
        elseif _diex == Ex
            sign = para.isFermi ? 1.0 : -1.0
        else
            error("not implemented!")
        end

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

######################### utility functions ############################
maxVer4TauIdx(para) = (para.innerLoopNum + 1) * para.interactionTauNum + para.firstTauIdx - 1
maxVer4LoopIdx(para) = para.firstLoopIdx + para.innerLoopNum - 1

function legBasis(chan::TwoBodyChannel, legK, loopIdx)
    KinL, KoutL, KinR, KoutR = legK[1], legK[2], legK[3], legK[4]
    K = zero(KinL)
    K[loopIdx] = 1
    if chan == PHr
        Kx = KoutL + K - KinL
        LLegK = [KinL, KoutL, Kx, K]
        RLegK = [K, Kx, KinR, KoutR]
    elseif chan == PHEr
        Kx = KoutR + K - KinL
        LLegK = [KinL, KoutR, Kx, K]
        RLegK = [K, Kx, KinR, KoutL]
    elseif chan == PPr
        Kx = KinL + KinR - K
        LLegK = [KinL, Kx, KinR, K]
        RLegK = [K, KoutL, Kx, KoutR]
    else
        error("not implemented!")
    end

    # check conservation and momentum assignment
    @assert LLegK[INL] ≈ KinL
    @assert LLegK[INL] + LLegK[INR] ≈ LLegK[OUTL] + LLegK[OUTR]
    @assert RLegK[INL] + RLegK[INR] ≈ RLegK[OUTL] + RLegK[OUTR]

    return LLegK, K, RLegK, Kx
end

function tauBasis(chan::TwoBodyChannel, LvT, RvT)
    G0T = (LvT[OUTR], RvT[INL])
    if chan == PHr
        extT = (LvT[INL], LvT[OUTL], RvT[INR], RvT[OUTR])
        GxT = (RvT[OUTL], LvT[INR])
    elseif chan == PHEr
        extT = (LvT[INL], RvT[OUTR], RvT[INR], LvT[OUTL])
        GxT = (RvT[OUTL], LvT[INR])
    elseif chan == PPr
        extT = (LvT[INL], RvT[OUTL], LvT[INR], RvT[OUTR])
        GxT = (LvT[OUTL], RvT[INR])
    else
        error("not implemented!")
    end

    # make sure all tidx are used once and only once
    t1 = sort(vcat(collect(G0T), collect(GxT), collect(extT)))
    t2 = sort(vcat(collect(LvT), collect(RvT)))
    @assert t1 == t2 "chan $(chan): G0=$G0T, Gx=$GxT, external=$extT don't match with Lver4 $LvT and Rver4 $RvT"
    @assert extT[INL] == LvT[INL]
    return extT, G0T, GxT
end


function factor(para, chan)
    Factor = SymFactor[Int(chan)] / (2π)^para.loopDim
    if para.isFermi == false
        Factor = abs(Factor)
    end
    return Factor
end

function typeMap(ltype, rtype)
    if (ltype == Instant || ltype == Dynamic) && (rtype == Instant || rtype == Dynamic)
        return Dynamic
    elseif (ltype == D_Instant || ltype == D_Dynamic) && (rtype == Instant || rtype == Dynamic)
        return D_Dynamic
    elseif (ltype == Instant || ltype == Dynamic) && (rtype == D_Instant || rtype == D_Dynamic)
        return D_Dynamic
    else
        return nothing
    end
end