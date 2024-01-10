mutable struct Green
    Tpair::Tuple{Int,Int}
    weight::Float64
    function Green(inT, outT)
        return new((inT, outT), 0.0)
    end
    function Green(tpair::Tuple{Int,Int})
        return new(tpair, 0.0)
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
    chan::Int
    Lver::_Ver4
    Rver::_Ver4
    map::Vector{IdxMap{_Ver4}}

    function Bubble(ver4::_Ver4, chan::Int, oL::Int, level::Int, _id::Vector{Int}) where {_Ver4}
        # @assert chan in para.chan "$chan isn't a bubble channels!"
        @assert oL < ver4.loopNum "LVer loopNum must be smaller than the ver4 loopNum"

        idbub = _id[1] # id vector will be updated later, so store the current id as the bubble id
        _id[1] += 1

        oR = ver4.loopNum - 1 - oL # loopNum of the right vertex
        lLpidxOffset = ver4.loopidxOffset + 1
        rLpidxOffset = lLpidxOffset + oL
        LTidx = ver4.TidxOffset  # the first τ index of the left vertex
        TauNum = interactionTauNum(ver4.para) # maximum tau number for each bare interaction
        RTidx = LTidx + (oL + 1) * TauNum   # the first τ index of the right sub-vertex

        if chan == T || chan == U
            LverChan = (level == 1) ? ver4.Fouter : ver4.F
            RverChan = (level == 1) ? ver4.Allouter : ver4.All
        elseif chan == S
            LverChan = (level == 1) ? ver4.Vouter : ver4.V
            RverChan = (level == 1) ? ver4.Allouter : ver4.All
        else
            error("chan $chan isn't implemented!")
        end

        # println("left ver chan: ", LsubVer, ", loop=", oL)
        # W = eltype(ver4.weight)
        Lver = _Ver4(ver4.para, LverChan, ver4.F, ver4.V, ver4.All;
            loopNum=oL, loopidxOffset=lLpidxOffset, tidxOffset=LTidx,
            level=level + 1, id=_id)

        Rver = _Ver4(ver4.para, RverChan, ver4.F, ver4.V, ver4.All;
            loopNum=oR, loopidxOffset=rLpidxOffset, tidxOffset=RTidx,
            level=level + 1, id=_id)

        @assert Lver.TidxOffset == ver4.TidxOffset "Lver Tidx must be equal to vertex4 Tidx! LoopNum: $(ver4.loopNum), LverLoopNum: $(Lver.loopNum), chan: $chan"

        ############## construct IdxMap ########################################
        map = []
        for (lt, LvT) in enumerate(Lver.Tpair)
            for (rt, RvT) in enumerate(Rver.Tpair)
                VerTidx = 0

                if chan == T
                    VerT = (LvT[INL], LvT[OUTL], RvT[INR], RvT[OUTR])
                    GTx = (RvT[OUTL], LvT[INR])
                elseif chan == U
                    VerT = (LvT[INL], RvT[OUTR], RvT[INR], LvT[OUTL])
                    GTx = (RvT[OUTL], LvT[INR])
                elseif chan == S
                    VerT = (LvT[INL], RvT[OUTL], LvT[INR], RvT[OUTR])
                    GTx = (LvT[OUTL], RvT[INR])
                else
                    throw("This channel is invalid!")
                end

                Gx = Green(GTx)
                push!(ver4.G[Int(chan)], Gx)

                GT0 = (LvT[OUTR], RvT[INL])
                G0 = Green(GT0)
                push!(ver4.G[1], G0)

                VerTidx = addTidx(ver4, VerT)
                for tpair in ver4.Tpair
                    @assert tpair[1] == ver4.TidxOffset "InL Tidx must be the same for all Tpairs in the vertex4"
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
        return new{_Ver4}(idbub, chan, Lver, Rver, map)
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
    para::Any
    chan::Vector{Int} # list of channels
    F::Vector{Int}
    V::Vector{Int}
    All::Vector{Int}
    Fouter::Vector{Int}
    Vouter::Vector{Int}
    Allouter::Vector{Int}

    ###### vertex topology information #####################
    id::Int
    level::Int

    #######  vertex properties   ###########################
    loopNum::Int
    loopidxOffset::Int # offset of the loop index from the first internal loop
    TidxOffset::Int # offset of the Tidx from the tau index of the left most incoming leg

    ######  components of vertex  ##########################
    G::SVector{16,Vector{Green}}  # large enough to host all Green's function
    bubble::Vector{Bubble{Ver4{W}}}

    ####### weight and tau table of the vertex  ###############
    Tpair::Vector{Tuple{Int,Int,Int,Int}}
    child::Vector{Vector{IdxMap{Ver4{W}}}}
    weight::Vector{W}

    function Ver4{W}(para, chan, F=[I, U, S], V=[I, T, U], All=union(F, V);
        loopNum=para.innerLoopNum, loopidxOffset=0, tidxOffset=0,
        Fouter=F, Vouter=V, Allouter=All,
        level=1, id=[1,]
    ) where {W}

        @assert para.totalTauNum >= (loopNum + 1) * interactionTauNum(para) "$para"

        if level > 1
            @assert Set(F) == Set(Fouter)
            @assert Set(V) == Set(Vouter)
            @assert Set(All) == Set(Allouter)
        end

        @assert (T in F) == false "F vertex is particle-hole irreducible, so that T channel is not allowed in F"
        @assert (S in V) == false "V vertex is particle-particle irreducible, so that S channel is not allowed in V"
        @assert (T in Fouter) == false "F vertex is particle-hole irreducible, so that T channel is not allowed in F"
        @assert (S in Vouter) == false "V vertex is particle-particle irreducible, so that S channel is not allowed in V"

        g = @SVector [Vector{Green}([]) for i = 1:16]
        ver4 = new{W}(para, chan, F, V, All, Fouter, Vouter, Allouter, id[1], level, loopNum, loopidxOffset, tidxOffset, g, [], [], [[],], [])
        id[1] += 1
        @assert loopNum >= 0


        if loopNum == 0
            tidx = tidxOffset
            # bare interaction may have one, two or four independent tau variables
            if interactionTauNum(para) == 1  # instantaneous interaction
                addTidx(ver4, (tidx, tidx, tidx, tidx)) #direct instant intearction
                addTidx(ver4, (tidx, tidx, tidx, tidx)) #exchange instant interaction
            end
            if interactionTauNum(para) == 2  # interaction with incoming and outing τ varibales
                addTidx(ver4, (tidx, tidx, tidx, tidx))  # direct and exchange instant interaction
                addTidx(ver4, (tidx, tidx, tidx + 1, tidx + 1))  # direct dynamic interaction
                addTidx(ver4, (tidx, tidx + 1, tidx + 1, tidx))  # exchange dynamic interaction
            end
            if interactionTauNum(para) == 4  # interaction with incoming and outing τ varibales
                error("Not implemented!")
                # addTidx(ver4, (tidx, tidx + 1, tidx + 2, tidx + 3))  # direct dynamic interaction
                # addTidx(ver4, (tidx, tidx + 3, tidx + 2, tidx + 1))  # exchange dynamic interaction
            end
        else # loopNum>0
            for c in chan
                for ol = 0:loopNum-1
                    bubble = Bubble(ver4, c, ol, level, id)
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
    if length(ver4.bubble) == 0
        return
    end

    G = ver4.G
    for bub in ver4.bubble
        Lver, Rver = bub.Lver, bub.Rver
        @assert Rver.loopidxOffset + Rver.loopNum == ver4.loopidxOffset + ver4.loopNum "Rver loopidx offset = $(Rver.loopidxOffset) + loopNum = $(Rver.loopNum) is not equal to ver4 loopidx offset = $(ver4.loopidxOffset) + loopNum = $(ver4.loopNum)"
        @assert Lver.loopNum + Rver.loopNum + 1 == ver4.loopNum
        for map in bub.map
            LverT, RverT = collect(Lver.Tpair[map.lidx]), collect(Rver.Tpair[map.ridx]) # 8 τ variables relevant for this bubble
            G1T, GxT = collect(map.G0.Tpair), collect(map.Gx.Tpair) # 4 internal variables
            ExtT = collect(ver4.Tpair[map.vidx]) # 4 external variables
            @assert compare(vcat(G1T, GxT, ExtT), vcat(LverT, RverT)) "chan $(bub.chan): G1=$G1T, Gx=$GxT, external=$ExtT don't match with Lver4 $LverT and Rver4 $RverT"
        end
    end
end