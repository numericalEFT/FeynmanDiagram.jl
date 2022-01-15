function buildVer4(para, LegK, chan, subdiagram = false; F = [I, U, S], V = [I, T, U], All = union(F, V),
    Fouter = F, Vouter = V, Allouter = All, diag = newDiagTree(para, :Ver4))

    ver4 = Ver4(diag, para, LegK, chan, subdiagram; F = F, V = V, All = All, Fouter = Fouter, Vouter = Vouter, Allouter = Allouter)

    return diag, ver4.nodes
end

maxVer4TauIdx(para) = (para.innerLoopNum + 1) * para.interactionTauNum + para.firstTauIdx - 1
maxVer4LoopIdx(para) = para.firstLoopIdx + para.innerLoopNum - 1

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

    level::Int
    extK::Vector{Vector{Float64}}
    nodes::Vector{Node{Vertex4}}

    function Ver4(diag, para, legK, chan, subdiagram = false; F = [I, U, S], V = [I, T, U], All = union(F, V),
        Fouter = F, Vouter = V, Allouter = All, level = 1, name = :none)

        @assert para.totalTauNum >= maxVer4TauIdx(para) "Increase totalTauNum!\n$para"
        @assert para.totalLoopNum >= maxVer4LoopIdx(para) "Increase totalLoopNum\n$para"

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

                    if c == T || c == U || c == S
                        addBubble!(ver4, c, p, level)
                    end
                end
            end
            # # TODO: add envolpe diagrams
            # # for c in II
            # # end
            # test(ver4) # more test
        end

        for n in ver4.nodes
            generate_node_from_children!(diag, n, ADD, 1.0, name; para = n.id)
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

    #the first loop idx is the inner loop of the bubble!
    LoopIdx = para.firstLoopIdx
    idx, maxLoop = findFirstLoopIdx(partition, LoopIdx + 1)
    LfirstLoopIdx, G0firstLoopIdx, RfirstLoopIdx, GxfirstLoopIdx = idx
    @assert maxLoop == maxVer4LoopIdx(para)

    diagType = [Ver4Diag, GreenDiag, Ver4Diag, GreenDiag]
    idx, maxTau = findFirstTauIdx(partition, diagType, para.firstTauIdx, TauNum)
    LfirstTauIdx, G0firstTauIdx, RfirstTauIdx, GxfirstTauIdx = idx
    @assert maxTau == maxVer4TauIdx(para) "Partition $partition with tauNum configuration $idx. maxTau = $maxTau, yet $(maxTauIdx(para)) is expected!"

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
    # println(K, ", ", Kx)

    lPara = reconstruct(para, innerLoopNum = oL, firstLoopIdx = LfirstLoopIdx, firstTauIdx = LfirstTauIdx)
    Lver = Ver4(diag, lPara, LLegK, LverChan, true; F = ver4.F, V = ver4.V, All = ver4.All, level = level + 1, name = :F)

    rPara = reconstruct(para, innerLoopNum = oR, firstLoopIdx = RfirstLoopIdx, firstTauIdx = RfirstTauIdx)
    Rver = Ver4(diag, rPara, RLegK, RverChan, true; F = ver4.F, V = ver4.V, All = ver4.All, level = level + 1, name = :All)

    gxPara = reconstruct(para, innerLoopNum = oGx, firstLoopIdx = GxfirstLoopIdx, firstTauIdx = GxfirstTauIdx)
    g0Para = reconstruct(para, innerLoopNum = oG0, firstLoopIdx = G0firstLoopIdx, firstTauIdx = G0firstTauIdx)

    bubble = Bubble(chan, g0Para, gxPara, Lver, Rver, ver4)

    for lnode in Lver.nodes
        for rnode in Rver.nodes
            addBubble2Diag!(diag, bubble, lnode, rnode, K, Kx)
        end
    end

    return bubble
end

function addBubble2Diag!(diag, bubble, lnode, rnode, K, Kx)
    chan = bubble.channel
    ver4, lver, rver = bubble.parent, bubble.lver, bubble.rver
    lid, rid = lnode.id, rnode.id
    ln, rn = lid.name, rid.name
    lc, rc = lid.DiEx, rid.DiEx
    vtype = typeMap(lid.type, rid.type)

    extT, G0T, GxT = tauBasis(chan, lid.extT, rid.extT)
    extK = ver4.extK

    Factor = factor(ver4.para, chan)

    # diag, g0 = buildG(bubble.g0, K, (LvT[OUTR], RvT[INL]); diag = diag)
    # diag, gc = buildG(bubble.gx, Kx, (RvT[OUTL], LvT[INR]); diag = diag)
    g0 = DiagTree.addpropagator!(diag, :Gpool, 0, :G0; site = G0T, loop = K)
    gc = DiagTree.addpropagator!(diag, :Gpool, 0, :Gx; site = GxT, loop = Kx)

    spin(response) = (response == UpUp ? "↑↑" : "↑↓")

    function add(Lresponse, Rresponse, responseName, factor = 1.0)
        if ln == Lresponse && rn == Rresponse
            nodeName = Symbol("$(spin(Lresponse))x$(spin(Rresponse)) → $chan,")
            id = Vertex4(responseName, vtype, BOTH, extK, extT, ver4.para)
            n = DiagTree.addnode!(diag, MUL, nodeName, [g0, gc, lnode.node, rnode.node], factor * Factor; para = id)
            add!(ver4.nodes, id, children = [n,])
        end
    end

    if chan == T
        add(UpUp, UpUp, UpUp, 1.0)
        add(UpDown, UpDown, UpUp, 1.0)
        add(UpUp, UpDown, UpDown, 1.0)
        add(UpDown, UpUp, UpDown, 1.0)
    elseif chan == U
        add(UpUp, UpUp, UpUp, 1.0)
        add(UpUp, UpUp, UpDown, 1.0)
        add(UpDown, UpDown, UpUp, 1.0)
        add(UpDown, UpDown, UpDown, 1.0)
        add(UpUp, UpDown, UpDown, -1.0)
        add(UpDown, UpUp, UpDown, -1.0)
    elseif chan == S
        add(UpUp, UpUp, UpUp, 1.0)
        add(UpDown, UpDown, UpDown, -2.0)
        add(UpUp, UpDown, UpDown, 1.0)
        add(UpDown, UpUp, UpDown, 1.0)
    else
        error("chan $chan isn't implemented!")
    end
end
