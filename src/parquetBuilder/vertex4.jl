function buildVer4(para, LegK, chan, subdiagram = false; F = [I, U, S], V = [I, T, U], All = union(F, V),
    Fouter = F, Vouter = V, Allouter = All, diag = newDiagTree(para, :Ver4))

    ver4 = Ver4(diag, para, LegK, chan, subdiagram; F = F, V = V, All = All, Fouter = Fouter, Vouter = Vouter, Allouter = Allouter)

    return diag, ver4.nodes
end

# function mapT(para, lnode, rnode)
#     lname, rname = lnode.name, rnode.name
#     lc, rc = lnode.DiEx, rnode.DiEx
#     if lname == ChargeCharge && rname == ChargeCharge
#         if lc == DI && rc == DI
#             return [(ChargeCharge, DI, para.spin),]
#         elseif (lc == DI && rc == EX) || (lc == EX && rc == DI)
#             return [(ChargeCharge, DI, 1.0),]
#         else
#             lc == EX && rc == EX
#             return [(ChargeCharge, EX, 1.0),]
#         end
#     else
#         error("not implemented!")
#     end
# end

function addT!(diag, bubble, lnode, rnode, K, Kx)
    ver4, lver, rver = bubble.parent, bubble.lver, bubble.rver
    lid, rid = lnode.id, rnode.id
    ln, rn = lid.name, rid.name
    lc, rc = lid.DiEx, rid.DiEx
    vtype = typeMap(lid.type, rid.type)

    LvT, RvT = lid.extT, rid.extT
    extT = (LvT[INL], LvT[OUTL], RvT[INR], RvT[OUTR])
    extK = ver4.extK

    # diag, g0 = buildG(bubble.g0, K, (LvT[OUTR], RvT[INL]); diag = diag)
    # diag, gc = buildG(bubble.gx, Kx, (RvT[OUTL], LvT[INR]); diag = diag)
    g0 = DiagTree.addpropagator!(diag, :Gpool, 0, :G0; site = (LvT[OUTR], RvT[INL]), loop = K)
    gc = DiagTree.addpropagator!(diag, :Gpool, 0, :Gx; site = (RvT[OUTL], LvT[INR]), loop = Kx)

    function add(nodeName, responseName, type, diex, factor = 1.0)
        id = Vertex4(responseName, type, diex, extK, extT)
        n = DiagTree.addnode!(diag, MUL, nodeName, [g0, gc, lnode.node, rnode.node], factor; para = id)
        add!(ver4.nodes, id, [n,])
    end

    function map(vname)
        if lc == DI && rc == DI
            if removeBubble(bubble, DI, DI) == false
                add(:Tdd, vname, vtype, DI, ver4.para.spin)
            end
        elseif lc == DI && rc == EX
            add(:Tde, vname, vtype, DI)
        elseif lc == EX && rc == DI
            add(:Ted, vname, vtype, DI)
        elseif lc == EX && rc == EX
            add(:Tee, vname, vtype, EX)
        else
            error("not exist!")
        end
    end

    if ln == ChargeCharge && rn == ChargeCharge
        map(ChargeCharge)
    elseif ln == SpinSpin && rn == SpinSpin
        map(SpinSpin)
    end

end

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
struct Ver4
    diag::Any
    para::GenericPara
    chan::Vector{Channel} # list of channels
    F::Vector{Channel}
    V::Vector{Channel}
    All::Vector{Channel}
    Fouter::Vector{Channel}
    Vouter::Vector{Channel}
    Allouter::Vector{Channel}

    ###### vertex topology information #####################
    level::Int

    ####### weight and tau table of the vertex  ###############
    extK::Vector{Vector{Float64}}

    nodes::Vector{Node{Vertex4}}

    function Ver4(diag, para, legK, chan, subdiagram = false; F = [I, U, S], V = [I, T, U], All = union(F, V),
        Fouter = F, Vouter = V, Allouter = All, level = 1)

        if level > 1
            @assert Set(F) == Set(Fouter)
            @assert Set(V) == Set(Vouter)
            @assert Set(All) == Set(Allouter)
        end

        @assert (T in F) == false "F vertex is particle-hole irreducible, so that T channel is not allowed in F=$F"
        @assert (S in V) == false "V vertex is particle-particle irreducible, so that S channel is not allowed in V=$V"
        @assert (T in Fouter) == false "F vertex is particle-hole irreducible, so that T channel is not allowed in F=$Fouter"
        @assert (S in Vouter) == false "V vertex is particle-particle irreducible, so that S channel is not allowed in V=$Vouter"

        @assert length(legK[1]) == length(legK[2]) == length(legK[3]) == para.totalLoopNum

        KinL, KoutL, KinR = legK[1], legK[2], legK[3]
        KoutR = (length(legK) > 3) ? legK[4] : KinL + KinR - KoutL
        @assert KoutR ≈ KinL + KinR - KoutL
        legK = [KinL, KoutL, KinR, KoutR]

        # g = @SVector [Vector{Green}([]) for i = 1:16]
        ver4 = new(diag, para, chan, F, V, All, Fouter, Vouter, Allouter, level, legK, [])

        @assert para.totalTauNum >= maxTauIdx(ver4) "Increase totalTauNum!\n$para"
        @assert para.totalLoopNum >= maxLoopIdx(ver4) "Increase totalLoopNum\n$para"

        loopNum = para.innerLoopNum
        @assert loopNum >= 0

        if loopNum == 0
            bareVer4!(ver4.nodes, diag, para, legK)
        else # loopNum>0
            for c in chan
                if c == I
                    continue
                end

                partition = orderedPartition(loopNum - 1, 4, 0)

                for p in partition

                    if c == T
                        # println(p)
                        addBubble!(ver4, c, p, level)
                        # if isnothing(bubble) == false && length(bubble.map) > 0  # if zero, bubble diagram doesn't exist
                        # push!(ver4.bubble, bubble)
                    end
                end
            end
            # # TODO: add envolpe diagrams
            # # for c in II
            # # end
            # test(ver4) # more test
        end

        for n in ver4.nodes
            if loopNum == 0 && subdiagram == true
                @assert length(n.children) == 1 "$(n.children)"
                n.node = n.children[1]
            else
                n.node = DiagTree.addnode!(diag, ADD, :sum, n.children; para = n.id)
            end
            @assert n.node != zero(Component)
        end

        return ver4
    end
end

struct Bubble
    channel::Channel
    g0::GenericPara
    gx::GenericPara
    lver::Ver4
    rver::Ver4
    parent::Ver4
end

function addBubble!(ver4::Ver4, chan::Channel, partition::Vector{Int}, level::Int)
    diag = ver4.diag
    para = ver4.para
    TauNum = para.interactionTauNum # maximum tau number for each bare interaction
    oL, oG0, oR, oGx = partition[1], partition[2], partition[3], partition[4]
    if isValidG(para.filter, oG0) == false || isValidG(para.filter, oGx) == false
        return
    end

    LoopIdx = para.firstLoopIdx
    idx, maxLoop = findFirstLoopIdx(partition, LoopIdx + 1)
    LfirstLoopIdx, G0firstLoopIdx, RfirstLoopIdx, GxfirstLoopIdx = idx
    @assert maxLoop == maxLoopIdx(ver4)

    diagType = [Ver4Diag, GreenDiag, Ver4Diag, GreenDiag]
    idx, maxTau = findFirstTauIdx(partition, diagType, para.firstTauIdx, TauNum)
    LfirstTauIdx, G0firstTauIdx, RfirstTauIdx, GxfirstTauIdx = idx
    @assert maxTau == maxTauIdx(ver4) "Partition $partition with tauNum configuration $idx. maxTau = $maxTau, yet $(maxTauIdx(ver4)) is expected!"

    if chan == T || chan == U
        LverChan = (level == 1) ? ver4.Fouter : ver4.F
        RverChan = (level == 1) ? ver4.Allouter : ver4.All
    elseif chan == S
        LverChan = (level == 1) ? ver4.Vouter : ver4.V
        RverChan = (level == 1) ? ver4.Allouter : ver4.All
    else
        error("chan $chan isn't implemented!")
    end

    LLegK, K, RLegK, Kx = legBasis(chan, ver4.extK, LoopIdx)

    lPara = reconstruct(para, innerLoopNum = oL, firstLoopIdx = LfirstLoopIdx, firstTauIdx = LfirstTauIdx)
    Lver = Ver4(diag, lPara, LLegK, LverChan, true; F = ver4.F, V = ver4.V, All = ver4.All, level = level + 1)

    rPara = reconstruct(para, innerLoopNum = oR, firstLoopIdx = RfirstLoopIdx, firstTauIdx = RfirstTauIdx)
    Rver = Ver4(diag, rPara, RLegK, RverChan, true; F = ver4.F, V = ver4.V, All = ver4.All, level = level + 1)

    gxPara = reconstruct(para, innerLoopNum = oGx, firstLoopIdx = GxfirstLoopIdx, firstTauIdx = GxfirstTauIdx)
    g0Para = reconstruct(para, innerLoopNum = oG0, firstLoopIdx = G0firstLoopIdx, firstTauIdx = G0firstTauIdx)

    bubble = Bubble(chan, g0Para, gxPara, Lver, Rver, ver4)

    for lnode in Lver.nodes
        for rnode in Rver.nodes

            if chan == T
                addT!(diag, bubble, lnode, rnode, K, Kx)
            elseif chan == U
            elseif chan == S
            else
                @error("Not implemented!")
            end
        end
    end

    Factor = SymFactor[Int(chan)] / (2π)^para.loopDim
    if para.isFermi == false
        Factor = abs(Factor)
    end

    return
end

function maxTauIdx(ver4::Ver4)
    para = ver4.para
    return (para.innerLoopNum + 1) * para.interactionTauNum + para.firstTauIdx - 1
end

function maxLoopIdx(ver4::Ver4)
    para = ver4.para
    return para.firstLoopIdx + para.innerLoopNum - 1
end

function compare(A, B)
    # check if the elements of XY are the same as Z
    XY, Z = copy(A), copy(B)
    for e in XY
        if (e in Z) == false
            return false
        end
        Z = (idx = findfirst(x -> x == e, Z)) > 0 ? deleteat!(Z, idx) : Z
    end
    return length(Z) == 0
end

function test(ver4)
    para = ver4.para
    if length(ver4.bubble) == 0
        return
    end

    @assert maxTauIdx(ver4) <= para.totalTauNum
    @assert maxLoopIdx(ver4) <= para.totalLoopNum

    G = ver4.G
    for bub in ver4.bubble
        Lver, Rver = bub.Lver, bub.Rver
        # G0 = G[1]
        # Gx = G[Int(bub.chan)]
        for map in bub.map
            if ver4.para.interactionTauNum > 0
                LverT, RverT = collect(Lver.Tpair[map.lidx]), collect(Rver.Tpair[map.ridx]) # 8 τ variables relevant for this bubble
                G1T, GxT = collect(map.G0.Tpair), collect(map.Gx.Tpair) # 4 internal variables
                ExtT = collect(ver4.Tpair[map.vidx]) # 4 external variables
                @assert compare(vcat(G1T, GxT, ExtT), vcat(LverT, RverT)) "chan $(bub.chan): G1=$G1T, Gx=$GxT, external=$ExtT don't match with Lver4 $LverT and Rver4 $RverT"

                tauSet = Set(vcat(G1T, GxT, ExtT))
                for t in tauSet
                    @assert t <= maxTauIdx(ver4) "Tauidx $t is too large! 
                    firstTauIdx = $(para.firstTauIdx), maxTauIdx =$(maxTauIdx(ver4)), loopNum=$(para.innerLoopNum)\n$para"
                end
            end
        end
    end
end
