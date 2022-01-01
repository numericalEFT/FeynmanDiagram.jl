function orderedPartition(_total, n, lowerbound = 1)
    @assert lowerbound >= 0
    total = _total - n * (lowerbound - 1)
    @assert total >= n
    unorderedPartition = collect(partitions(total, n))
    #e.g., loopNum =5, n =2 ==> unordered = [[4, 1], [3, 2]]
    orderedPartition = Vector{Vector{Int}}([])
    for p in unorderedPartition
        p = p .+ (lowerbound - 1)
        @assert sum(p) == _total
        for i in p
            @assert i >= lowerbound
        end
        append!(orderedPartition, Set(permutations(p)))
    end
    #e.g., loopNum =5, n =2 ==> ordered = [[4, 1], [1, 4], [3, 2], [2, 3]]
    return orderedPartition
end

function findFirstLoopIdx(partition, isG, firstidx::Int)
    @assert length(partition) == length(isG)
    accumulated = accumulate(+, partition; init = firstidx) #  idx[i] = firstidx + p[1]+p[2]+...+p[i]
    firstLoopIdx = [firstidx,]
    append!(firstLoopIdx, accumulated[1:end-1])
    maxLoopIdx = accumulated[end] - 1
    return firstLoopIdx, maxLoopIdx
end

function findFirstTauIdx(partition, isG, firstidx) end

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

        oL, oR, oG0, oGx = partition[1], partition[2], partition[3], partition[4]
        # example 1: [1, 2, 1, 1] - total loop list [1, 2, 3, 4, 5, 6], tau list [1, 2, 3, 4, 5, 6]
        # example 2: [1, 2, 0, 0] - total loop list [1, 2, 3, 4], tau list [1, 2, 3, 4, 5]

        LoopIdx = para.firstLoopIdx
        # eg1: 1 - [1, ], eg2: 1 - [1, ]
        LfirstLoopIdx = para.firstLoopIdx + 1
        # eg1: 2 - [2, ], eg2: 2 - [2, ]
        G0firstLoopIdx = LfirstLoopIdx + oL #first inner loop index of G0 (may not exist if G0 loop number is zero) 
        # eg1: 3 - [3, ], eg2: 3 - [,]
        RfirstLoopIdx = LfirstLoopIdx + oL + oG0
        # eg1: 4 - [4, 5], eg2: 3 - [3, 4]
        GxfirstLoopIdx = RfirstLoopIdx + oR #first inner loop index of Gx (may not exist if Gx loop number is zero)
        # eg1: 6 - [6, ], eg2: 5 - [,]
        loopRight = RfirstLoopIdx + oR + oGx - 1
        # eg1: 4+2+1-1 == 6, eg2: 3+2+0-1 == 4
        @assert loopRight == maxLoopIdx(ver4)

        LfirstTauIdx = para.firstTauIdx
        # eg1: 1 - [1, 2], eg2: 1 - [1, 2]
        G0firstTauIdx = LfirstTauIdx + (oL + 1) * TauNum
        # eg1: 3 - [3, ], eg2: 3 - [, ]
        RfirstTauIdx = LfirstTauIdx + (oL + 1) * TauNum + oG0 * TauNum
        # eg1: 4 - [4, 5, 6], eg2: 3 - [3, 4, 5]
        GxfirstTauIdx = RfirstTauIdx + (oR + 1) * TauNum
        # eg1: 7 - [7, ], eg2: 6 - [, ]
        tauRight = RfirstTauIdx + (oR + 1) * TauNum + oGx * TauNum - 1
        # eg1: 4+3+1-1==7, eg2: 3+3+0-1==5
        @assert tauRight == maxTauIdx(ver4)

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
struct Ver4{W}
    para::GenericPara
    chan::Vector{Channel} # list of channels
    F::Vector{Channel}
    V::Vector{Channel}
    All::Vector{Channel}
    Fouter::Vector{Channel}
    Vouter::Vector{Channel}
    Allouter::Vector{Channel}

    ###### vertex topology information #####################
    id::Int
    level::Int

    ######  components of vertex  ##########################
    G::SVector{16,Vector{Green}}  # large enough to host all Green's function
    bubble::Vector{Bubble{Ver4{W}}}

    ####### weight and tau table of the vertex  ###############
    legK::Vector{Vector{Float64}}
    Tpair::Vector{Tuple{Int,Int,Int,Int}}
    child::Vector{Vector{IdxMap{Ver4{W}}}}
    weight::Vector{W}

    function Ver4{W}(para, legK, chan, F, V, All = union(F, V);
        Fouter = F, Vouter = V, Allouter = All,
        level = 1, id = [1,]
    ) where {W}

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

        g = @SVector [Vector{Green}([]) for i = 1:16]
        ver4 = new{W}(para, chan, F, V, All, Fouter, Vouter, Allouter, id[1], level, g, [], [KinL, KoutL, KinR, KoutR], [], [[],], [])

        @assert para.totalTauNum >= maxTauIdx(ver4) "Increase totalTauNum!\n$para"
        @assert para.totalLoopNum >= maxLoopIdx(ver4) "Increase totalLoopNum\n$para"

        id[1] += 1

        loopNum = para.innerLoopNum

        @assert loopNum >= 0


        if loopNum == 0
            tidx = para.firstTauIdx
            # bare interaction may have one, two or four independent tau variables
            if para.interactionTauNum == 1 || para.interactionTauNum == 0  # instantaneous interaction (interaction without time variable)
                addTidx(ver4, (tidx, tidx, tidx, tidx)) #direct instant intearction
                addTidx(ver4, (tidx, tidx, tidx, tidx)) #exchange instant interaction
            end
            if para.interactionTauNum == 2  # interaction with incoming and outing τ varibales
                addTidx(ver4, (tidx, tidx, tidx, tidx))  # direct and exchange instant interaction
                addTidx(ver4, (tidx, tidx, tidx + 1, tidx + 1))  # direct dynamic interaction
                addTidx(ver4, (tidx, tidx + 1, tidx + 1, tidx))  # exchange dynamic interaction
            end
            if para.interactionTauNum == 4  # interaction with incoming and outing τ varibales
                error("Not implemented!")
                # addTidx(ver4, (tidx, tidx + 1, tidx + 2, tidx + 3))  # direct dynamic interaction
                # addTidx(ver4, (tidx, tidx + 3, tidx + 2, tidx + 1))  # exchange dynamic interaction
            end
        else # loopNum>0
            for c in chan
                if c == I
                    continue
                end

                partition = orderedPartition(loopNum - 1, 4, 0)
                if Girreducible in para.filter
                    #if one-partitcle irreducible, then G0 at p[3] and Gx at p[4] must be loop 0
                    partition = [p for p in partition if p[3] == 0 && p[4] == 0]
                end

                for p in partition
                    # println(p)
                    bubble = Bubble(ver4, c, p, level, id)
                    if length(bubble.map) > 0  # if zero, bubble diagram doesn't exist
                        push!(ver4.bubble, bubble)
                    end
                end
            end
            # TODO: add envolpe diagrams
            # for c in II
            # end
            test(ver4) # more test
        end
        return ver4
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

function tpair(ver4, MaxT = 18)
    s = "\u001b[31m$(ver4.id):\u001b[0m"
    if ver4.para.innerLoopNum > 0
        s *= "$(ver4.para.innerLoopNum)lp, T$(length(ver4.Tpair))⨁ "
    else
        s *= "⨁ "
    end
    # if ver4.loopNum <= 1
    for (ti, T) in enumerate(ver4.Tpair)
        if ti <= MaxT
            s *= "($(T[1]),$(T[2]),$(T[3]),$(T[4]))"
        else
            s *= "..."
            break
        end
    end
    # end
    return s
end

##### pretty print of Bubble and Ver4  ##########################
Base.show(io::IO, bub::Bubble) = AbstractTrees.printnode(io::IO, bub)
Base.show(io::IO, ver4::Ver4) = AbstractTrees.printnode(io::IO, ver4)

################## implement AbstractTrees interface #######################
# refer to https://github.com/JuliaCollections/AbstractTrees.jl for more details
function AbstractTrees.children(ver4::Ver4)
    return ver4.bubble
end

function AbstractTrees.children(bubble::Bubble)
    return (bubble.Lver, bubble.Rver)
end

function iterate(ver4::Ver4{W}) where {W}
    if length(ver4.bubble) == 0
        return nothing
    else
        return (ver4.bubble[1], 1)
    end
end

function iterate(bub::Bubble)
    return (bub.Lver, false)
end

function iterate(ver4::Ver4{W}, state) where {W}
    if state >= length(ver4.bubble) || length(ver4.bubble) == 0
        return nothing
    else
        return (ver4.bubble[state+1], state + 1)
    end
end

function iterate(bub::Bubble, state::Bool)
    state && return nothing
    return (bub.Rver, true)
end

Base.IteratorSize(::Type{Ver4{W}}) where {W} = Base.SizeUnknown()
Base.eltype(::Type{Ver4{W}}) where {W} = Ver4{W}

Base.IteratorSize(::Type{Bubble{Ver4{W}}}) where {W} = Base.SizeUnknown()
Base.eltype(::Type{Bubble{Ver4{W}}}) where {W} = Bubble{Ver4{W}}

AbstractTrees.printnode(io::IO, ver4::Ver4) = print(io, tpair(ver4))
AbstractTrees.printnode(io::IO, bub::Bubble) = print(io,
    "\u001b[32m$(bub.id): $(bub.chan) $(bub.Lver.para.innerLoopNum)Ⓧ $(bub.Rver.para.innerLoopNum)\u001b[0m")
