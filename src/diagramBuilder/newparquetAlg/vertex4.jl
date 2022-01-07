mutable struct Weight{T} <: FieldVector{2,T}
    d::T
    e::T
    Weight{T}() where {T} = new{T}(T(0), T(0))
    Weight(d::T, e::T) where {T} = new{T}(d::T, e::T)
    Weight{T}(d, e) where {T} = new{T}(T(d), T(e))
end

Base.zero(::Type{Weight{T}}) where {T} = Weight(T(0), T(0))
Base.abs(w::Weight) = abs(w.d) + abs(w.e) # define abs(Weight)

mutable struct Green
    para::GenericPara
    loopBasis::Vector{Float64}
    Tpair::Tuple{Int,Int}
    weight::Float64
    function Green(para, loopBasis, tpair::Tuple{Int,Int})
        return new(para, loopBasis, tpair, 0.0)
    end
end

# add Tpairs to Green's function (in, out) or vertex4 (inL, outL, inR, outR)
function addTidx(obj, _Tidx)
    for (i, Tidx) in enumerate(obj.Tpair)
        if Tidx == _Tidx
            return i
        end
    end
    push!(obj.Tpair, _Tidx)
    push!(obj.weight, zero(eltype(obj.weight))) # add zero to the weight table of the object
    push!(obj.child, [])
    return length(obj.Tpair)
end

"""
    struct IdxMap{_Ver4}

    Map left vertex Tpair[lidx], right vertex Tpair[ridx], the shared Green's function G0 and the channel specific Green's function Gx to the top level 4-vertex Tpair[vidx]

# Arguments
- l::_Ver4 : left vertex 
- r::_Ver4 : right vertex
- v::_Ver4 : composte vertex
- lidx::Int : left sub-vertex index
- ridx::Int : right sub-vertex index
- vidx::Int : composite vertex index
- G0::Green : shared Green's function index
- Gx::Green : channel specific Green's function index
"""
mutable struct IdxMap{_Ver4}
    l::_Ver4 #left vertex 
    r::_Ver4 #right vertex
    v::_Ver4 #composte vertex
    lidx::Int # left sub-vertex index
    ridx::Int # right sub-vertex index
    vidx::Int # composite vertex index
    G0::Green # shared Green's function index
    Gx::Green # channel specific Green's function index
    node::Any # useful for constructing DiagTree
    function IdxMap(l::V, r::V, v::V, lidx, ridx, vidx, G0::Green, Gx::Green) where {V}
        return new{V}(l, r, v, lidx, ridx, vidx, G0, Gx, nothing)
    end
end

struct Bubble{_Ver4} # template Bubble to avoid mutually recursive struct
    id::Int
    chan::Channel
    loopIdx::Int # the index of the inner loop of the bubble
    factor::Float64
    Lver::_Ver4
    Rver::_Ver4
    map::Vector{IdxMap{_Ver4}}

    function Bubble(ver4::_Ver4, chan::Channel, partition::Vector{Int}, level::Int, _id::Vector{Int}) where {_Ver4}
        # @assert chan in para.chan "$chan isn't a bubble channels!"
        # @assert oL < ver4.loopNum "LVer loopNum must be smaller than the ver4 loopNum"
        para = ver4.para
        @assert sum(partition) == para.innerLoopNum - 1 "partition = $partition should sum to ver4 loopNum -1 =$(para.innerLoopNum-1)"
        for p in partition
            @assert p >= 0
        end

        TauNum = para.interactionTauNum # maximum tau number for each bare interaction

        idbub = _id[1] # id vector will be updated later, so store the current id as the bubble id
        _id[1] += 1

        oL, oG0, oR, oGx = partition[1], partition[2], partition[3], partition[4]

        if isValidG(para.filter, oG0) == false || isValidG(para.filter, oGx) == false
            return nothing
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

        ############# prepare LegK for left, right subvertex ################
        legK = ver4.legK
        KinL, KoutL, KinR, KoutR = legK[1], legK[2], legK[3], legK[4]
        K = zero(KinL)
        K[LoopIdx] = 1
        Kt = KoutL + K - KinL
        Ku = KoutR + K - KinL
        Ks = KinL + KinR - K

        LLegK, RLegK = [], []
        if chan == T
            LLegK = [KinL, KoutL, Kt, K]
            RLegK = [K, Kt, KinR, KoutR]
        elseif chan == U
            LLegK = [KinL, KoutR, Ku, K]
            RLegK = [K, Ku, KinR, KoutL]
        elseif chan == S
            LLegK = [KinL, Ks, KinR, K]
            RLegK = [K, KoutL, Ks, KoutR]
        else
            error("not implemented!")
        end

        lPara = reconstruct(para, innerLoopNum = oL, firstLoopIdx = LfirstLoopIdx, firstTauIdx = LfirstTauIdx)
        Lver = _Ver4(lPara, LLegK, LverChan, ver4.F, ver4.V, ver4.All;
            level = level + 1, id = _id)

        rPara = reconstruct(para, innerLoopNum = oR, firstLoopIdx = RfirstLoopIdx, firstTauIdx = RfirstTauIdx)
        Rver = _Ver4(rPara, RLegK, RverChan, ver4.F, ver4.V, ver4.All;
            level = level + 1, id = _id)

        @assert lPara.firstTauIdx == para.firstTauIdx "Lver Tidx must be equal to vertex4 Tidx! LoopNum: $(para.innerLoopNum), LverLoopNum: $(lPara.innerLoopNum), chan: $chan"

        ############## construct IdxMap ########################################
        map = []
        for (lt, LvT) in enumerate(Lver.Tpair)
            for (rt, RvT) in enumerate(Rver.Tpair)
                VerTidx = 0

                if chan == T
                    VerT = (LvT[INL], LvT[OUTL], RvT[INR], RvT[OUTR])
                    GTx = (RvT[OUTL], LvT[INR])
                    Kx = Kt
                elseif chan == U
                    VerT = (LvT[INL], RvT[OUTR], RvT[INR], LvT[OUTL])
                    GTx = (RvT[OUTL], LvT[INR])
                    Kx = Ku
                elseif chan == S
                    VerT = (LvT[INL], RvT[OUTL], LvT[INR], RvT[OUTR])
                    GTx = (LvT[OUTL], RvT[INR])
                    Kx = Ks
                else
                    error("This channel is invalid!")
                end

                gxPara = reconstruct(para, innerLoopNum = oGx, firstLoopIdx = GxfirstLoopIdx, firstTauIdx = GxfirstTauIdx)
                Gx = Green(gxPara, Kx, GTx)
                push!(ver4.G[Int(chan)], Gx)

                GT0 = (LvT[OUTR], RvT[INL])
                g0Para = reconstruct(para, innerLoopNum = oG0, firstLoopIdx = G0firstLoopIdx, firstTauIdx = G0firstTauIdx)
                G0 = Green(g0Para, K, GT0)
                push!(ver4.G[1], G0)

                VerTidx = addTidx(ver4, VerT)
                for tpair in ver4.Tpair
                    @assert tpair[1] == para.firstTauIdx "InL Tidx must be the same for all Tpairs in the vertex4"
                end

                ###### test if the internal + exteranl variables is equal to the total 8 variables of the left and right sub-vertices ############
                Total1 = vcat(collect(LvT), collect(RvT))
                Total2 = vcat(collect(GT0), collect(GTx), collect(VerT))
                # println(Total1)
                # println(Total2)
                @assert compare(Total1, Total2) "chan $(chan): G0=$GT0, Gx=$GTx, external=$VerT don't match with Lver4 $LvT and Rver4 $RvT"

                idxmap = IdxMap(Lver, Rver, ver4, lt, rt, VerTidx, G0, Gx)
                push!(map, idxmap)
                push!(ver4.child[VerTidx], idxmap)
            end
        end

        Factor = SymFactor[Int(chan)] / (2π)^para.loopDim
        if para.isFermi == false
            Factor = abs(Factor)
        end

        return new{_Ver4}(idbub, chan, LoopIdx, Factor, Lver, Rver, map)
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

#Remark:
- AbstractTrees interface is implemented for Ver4. So one can use the API in https://juliacollections.github.io/AbstractTrees.jl/stable/ to manipulate/print the tree structre of Ver4.
- There are three different methods to print/visualize the tree structre: 
1) `print_tree(ver4::Ver4)` or `print_tree(bub::Bubble)` to print the tree to terminal. This function is provided by AbstractTrees API. 
2) `newick(ver4::Ver4)` or `newick(bub::Bubble)` to serilize the tree to a newick format string. You may save the string to a text file, then visualize it with a newick format visualizer application. 
3) `showTree(ver4::Ver4)` to visualize the tree using the python package ete3. You have to install ete3 properly to use this function.
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

    ######  components of vertex  ##########################
    # bubble::Vector{Bubble{Ver4{W}}}

    ####### weight and tau table of the vertex  ###############
    legK::Vector{Vector{Float64}}

    node::Set{Ver4identifier}
    # child::Dict{Ver4identifier,IdxMap{Ver4}}

    function Ver4(para, legK, chan, F, V, All = union(F, V);
        Fouter = F, Vouter = V, Allouter = All,
        level = 1, diag = newDiagTree(para, :Ver4)
    )

        if level > 1
            @assert Set(F) == Set(Fouter)
            @assert Set(V) == Set(Vouter)
            @assert Set(All) == Set(Allouter)
        end

        @assert (T in F) == false "F vertex is particle-hole irreducible, so that T channel is not allowed in F"
        @assert (S in V) == false "V vertex is particle-particle irreducible, so that S channel is not allowed in V"
        @assert (T in Fouter) == false "F vertex is particle-hole irreducible, so that T channel is not allowed in F"
        @assert (S in Vouter) == false "V vertex is particle-particle irreducible, so that S channel is not allowed in V"

        @assert length(legK[1]) == length(legK[2]) == length(legK[3]) == para.totalLoopNum

        KinL, KoutL, KinR = legK[1], legK[2], legK[3]
        KoutR = (length(legK) > 3) ? legK[4] : KinL + KinR - KoutL
        @assert KoutR ≈ KinL + KinR - KoutL
        legK = [KinL, KoutL, KinR, KoutR]

        # g = @SVector [Vector{Green}([]) for i = 1:16]
        ver4 = new(diag, para, chan, F, V, All, Fouter, Vouter, Allouter, level, legK, Set([]))

        @assert para.totalTauNum >= maxTauIdx(ver4) "Increase totalTauNum!\n$para"
        @assert para.totalLoopNum >= maxLoopIdx(ver4) "Increase totalLoopNum\n$para"

        loopNum = para.innerLoopNum
        @assert loopNum >= 0

        if loopNum == 0
            zeroLoopVer4Node!(ver4.node, diag, para, legK)
        else # loopNum>0
            # for c in chan
            #     if c == I
            #         continue
            #     end

            #     partition = orderedPartition(loopNum - 1, 4, 0)

            #     for p in partition
            #         # println(p)
            #         bubble = Bubble(ver4, c, p, level, id)
            #         if isnothing(bubble) == false && length(bubble.map) > 0  # if zero, bubble diagram doesn't exist
            #             push!(ver4.bubble, bubble)
            #         end
            #     end
            # end
            # # TODO: add envolpe diagrams
            # # for c in II
            # # end
            # test(ver4) # more test
        end
        return diag, ver4.node
    end
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
