"""
Para(chan, interactionTauNum)

Parameters to generate diagrams using Parquet algorithm

#Arguments

- `chan`: list of channels of sub-vertices
- `interactionTauNum`: τ degrees of freedom of the bare interaction
"""
struct Para
    chan::Vector{Int}
    F::Vector{Int}
    V::Vector{Int}
    interactionTauNum::Vector{Int} # list of possible τ degrees of freedom of the bare interaction 0, 2, or 4
    function Para(chan, interactionTauNum)

        for tnum in interactionTauNum
            @assert tnum == 1 || tnum == 2 || tnum == 4
        end

        for c in chan
            @assert c in Allchan "$chan $c isn't implemented!"
        end
        F = intersect(chan, Fchan)
        V = intersect(chan, Vchan)

        return new(chan, F, V, interactionTauNum)
    end
end

struct Green{W}
    Tpair::Vector{Tuple{Int,Int}}
    weight::Vector{W}
    function Green{W}() where {W} 
        return new{W}([], [])
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
    return length(obj.Tpair)
end

"""
   struct IdxMap(lv, rv, G0, Gx, ver)

    Map left vertex Tpair[lv], right vertex Tpair[rv], the shared Green's function G0[G0] and the channel specific Green's function Gx[Gx] to the top level 4-vertex Tpair[ver]
"""
struct IdxMap
    lv::Int # left sub-vertex index
    rv::Int # right sub-vertex index
    G0::Int # shared Green's function index
    Gx::Int # channel specific Green's function index
    ver::Int # composite vertex index
end

struct Bubble{_Ver4,W} # template Bubble to avoid mutually recursive struct
    id::Int
    chan::Int
    Lver::_Ver4
    Rver::_Ver4
    map::Vector{IdxMap}

    function Bubble{_Ver4,W}(ver4::_Ver4, chan::Int, oL::Int, para::Para, level::Int, _id::Vector{Int}) where {_Ver4,W}
        @assert chan in para.chan "$chan isn't a bubble channels!"
        @assert oL < ver4.loopNum "LVer loopNum must be smaller than the ver4 loopNum"

        idbub = _id[1] # id vector will be updated later, so store the current id as the bubble id
        _id[1] += 1  

        oR = ver4.loopNum - 1 - oL # loopNum of the right vertex
        LTidx = ver4.Tidx  # the first τ index of the left vertex
        maxTauNum = maximum(para.interactionTauNum) # maximum tau number for each bare interaction
        RTidx = LTidx + (oL + 1) * maxTauNum + 1  # the first τ index of the right sub-vertex

        if chan == T || chan == U
            LsubVer = para.F
            RsubVer = para.chan
        elseif chan == S
            LsubVer = para.V
            RsubVer = para.chan
        else
            error("chan $chan isn't implemented!")
        end

        Lver = _Ver4{W}(oL, LTidx, para; chan=LsubVer, level=level + 1, id=_id)
        Rver = _Ver4{W}(oR, RTidx, para; chan=RsubVer, level=level + 1, id=_id)

        @assert Lver.Tidx == ver4.Tidx "Lver Tidx must be equal to vertex4 Tidx! LoopNum: $(ver4.loopNum), LverLoopNum: $(Lver.loopNum), chan: $chan"

        ############## construct IdxMap ########################################
        map = []
        G = ver4.G
        for (lt, LvT) in enumerate(Lver.Tpair)
            for (rt, RvT) in enumerate(Rver.Tpair)
                GT0 = (LvT[OUTR], RvT[INL])
                GT0idx = addTidx(G[1], GT0)
                GTxidx, VerTidx = 0, 0

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

                VerTidx = addTidx(ver4, VerT)
                GTxidx = addTidx(G[chan], GTx)

                for tpair in ver4.Tpair
                    @assert tpair[1] == ver4.Tidx "InL Tidx must be the same for all Tpairs in the vertex4"
                end

                ###### test if the internal + exteranl variables is equal to the total 8 variables of the left and right sub-vertices ############
                Total1 = vcat(collect(LvT), collect(RvT))
                Total2 = vcat(collect(GT0), collect(GTx), collect(VerT))
                @assert compare(Total1, Total2) "chan $(ChanName[chan]): G0=$GT0, Gx=$GTx, external=$VerT don't match with Lver4 $LvT and Rver4 $RvT" 

                push!(map, IdxMap(lt, rt, GT0idx, GTxidx, VerTidx))
            end
        end
        return new(idbub, chan, Lver, Rver, map)
    end
end

"""
    Ver4{W}(loopNum, tidx, para::Para; chan=nothing, level=1, id=[1, ]) where {W}

    Generate 4-vertex diagrams using Parquet Algorithm

#Arguments
- `loopNum`: momentum loop degrees of freedom of the 4-vertex diagrams
- `tidx`: the first τ variable index. It is also the τ variable of the left incoming electron for all 4-vertex diagrams
- `para`: parameters
- `chan`: list of channels of the current 4-vertex. If not specified, it is set to be `para.chan`
- `level`: level in the diagram tree
- `id`: the first element will be used as the id of the Ver4. All nodes in the tree will be labeled in preorder depth-first search

#Remark:
- The argument `chan` and `para.chan` are different. The former is the channels of current 4-vertex, while the later is the channels of the sub-vertices
- AbstractTrees interface is implemented for Ver4. So one can use the API in https://juliacollections.github.io/AbstractTrees.jl/stable/ to manipulate/print the tree structre of Ver4.
- There are three different methods to print/visualize the tree structre: 
1) `print_tree(ver4::Ver4)` or `print_tree(bub::Bubble)` to print the tree to terminal. This function is provided by AbstractTrees API. 
2) `newick(ver4::Ver4)` or `newick(bub::Bubble)` to serilize the tree to a newick format string. You may save the string to a text file, then visualize it with a newick format visualizer application. 
3) `showTree(ver4::Ver4)` to visualize the tree using the python package ete3. You have to install ete3 properly to use this function.
"""
struct Ver4{W}
    ###### vertex topology information #####################
    id::Int
    level::Int
    
    #######  vertex properties   ###########################
    loopNum::Int
    chan::Vector{Int} # list of channels
    Tidx::Int # inital Tidx

    ######  components of vertex  ##########################
    G::SVector{16,Green}  # large enough to host all Green's function
    bubble::Vector{Bubble{Ver4}}

    ####### weight and tau table of the vertex  ###############
    Tpair::Vector{Tuple{Int,Int,Int,Int}}
    weight::Vector{W}

    function Ver4{W}(loopNum, tidx, para::Para; chan=para.chan, level=1, id=[1, ]) where {W}
        g = @SVector [Green{W}() for i = 1:16]
        ver4 = new{W}(id[1], level, loopNum, chan, tidx, g, [], [], [])
        id[1] += 1
        @assert loopNum >= 0
        if loopNum == 0
            # bare interaction may have one, two or four independent tau variables
            if 1 in para.interactionTauNum  # instantaneous interaction
                addTidx(ver4, (tidx, tidx, tidx, tidx)) 
            end
            if 2 in para.interactionTauNum  # interaction with incoming and outing τ varibales
                addTidx(ver4, (tidx, tidx, tidx + 1, tidx + 1))  # direct dynamic interaction
                addTidx(ver4, (tidx, tidx + 1, tidx + 1, tidx))  # exchange dynamic interaction
            end
            if 4 in para.interactionTauNum  # interaction with incoming and outing τ varibales
                addTidx(ver4, (tidx, tidx + 1, tidx + 2, tidx + 3))  # direct dynamic interaction
                addTidx(ver4, (tidx, tidx + 3, tidx + 2, tidx + 1))  # exchange dynamic interaction
            end
        else # loopNum>0
            for c in para.chan
                for ol = 0:loopNum - 1
                    bubble = Bubble{Ver4,W}(ver4, c, ol, para, level, id)
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
        for map in bub.map
            LverT, RverT = collect(Lver.Tpair[map.lv]), collect(Rver.Tpair[map.rv]) # 8 τ variables relevant for this bubble
            G1T, GxT = collect(G[1].Tpair[map.G0]), collect(G[bub.chan].Tpair[map.Gx]) # 4 internal variables
            ExtT = collect(ver4.Tpair[map.ver]) # 4 external variables
            @assert compare(vcat(G1T, GxT, ExtT), vcat(LverT, RverT)) "chan $(ChanName[bub.chan]): G1=$G1T, Gx=$GxT, external=$ExtT don't match with Lver4 $LverT and Rver4 $RverT" 
end
end
end

function tpair(ver4, MaxT=18)
    s = "\u001b[31m$(ver4.id):\u001b[0m"
    if ver4.loopNum > 0
        s *= "$(ver4.loopNum)lp, T$(length(ver4.Tpair))⨁ "
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

function iterate(bub::Bubble{Ver4{W},W}) where {W}
    return (bub.Lver, false)
end

function iterate(ver4::Ver4{W}, state) where {W}
    if state >= length(ver4.bubble) || length(ver4.bubble) == 0
        return nothing
        else
        return (ver4.bubble[state + 1], state + 1)
    end
end

function iterate(bub::Bubble{Ver4{W},W}, state::Bool) where {W}
    state && return nothing
    return (bub.Rver, true)
end

Base.IteratorSize(::Type{Ver4{W}}) where W = Base.SizeUnknown()
Base.eltype(::Type{Ver4{W}}) where W = Ver4{W}

Base.IteratorSize(::Type{Bubble{Ver4{W},W}}) where W = Base.SizeUnknown()
Base.eltype(::Type{Bubble{Ver4{W},W}}) where W = Bubble{Ver4{W},W}

AbstractTrees.printnode(io::IO, ver4::Ver4) = print(io, tpair(ver4))
AbstractTrees.printnode(io::IO, bub::Bubble) = print(io, "\u001b[32m$(bub.id): $(ChanName[bub.chan]) $(bub.Lver.loopNum)Ⓧ $(bub.Rver.loopNum)\u001b[0m")


# function eval(ver4::Ver4, KinL, KoutL, KinR, KoutR, Kidx::Int, fast=false)
#     if ver4.loopNum == 0
#         DiagType == POLAR ?
#         ver4.weight[1] = interaction(KinL, KoutL, KinR, KoutR, ver4.inBox, norm(varK[0])) :
#         ver4.weight[1] = interaction(KinL, KoutL, KinR, KoutR, ver4.inBox)
#         return
#     end

#     # LoopNum>=1
#     for w in ver4.weight
#         w .= 0.0 # initialize all weights
#     end
#     G = ver4.G
#     K, Kt, Ku, Ks = (varK[Kidx], ver4.K[1], ver4.K[2], ver4.K[3])
#     eval(G[1], K, varT)
#     bubWeight = counterBubble(K)

#     for c in ver4.chan
#         if c == T || c == TC
#             Kt .= KoutL .+ K .- KinL
#             if (!ver4.inBox)
#                 eval(G[T], Kt)
#             end
#         elseif c == U || c == UC
#             # can not be in box!
#             Ku .= KoutR .+ K .- KinL
#             eval(G[U], Ku)
#         else
#             # S channel, and cann't be in box!
#             Ks .= KinL .+ KinR .- K
#             eval(G[S], Ks)
#         end
#     end
#     for b in ver4.bubble
#         c = b.chan
#         Factor = SymFactor[c] * PhaseFactor
#         Llopidx = Kidx + 1
#         Rlopidx = Kidx + 1 + b.Lver.loopNum

#         if c == T || c == TC
#             eval(b.Lver, KinL, KoutL, Kt, K, Llopidx)
#             eval(b.Rver, K, Kt, KinR, KoutR, Rlopidx)
#         elseif c == U || c == UC
#             eval(b.Lver, KinL, KoutR, Ku, K, Llopidx)
#             eval(b.Rver, K, Ku, KinR, KoutL, Rlopidx)
#         else
#             # S channel
#             eval(b.Lver, KinL, Ks, KinR, K, Llopidx)
#             eval(b.Rver, K, KoutL, Ks, KoutR, Rlopidx)
#         end

#         rN = length(b.Rver.weight)
#         gWeight = 0.0
#         for (l, Lw) in enumerate(b.Lver.weight)
#             for (r, Rw) in enumerate(b.Rver.weight)
#                 map = b.map[(l - 1) * rN + r]

#                     if ver4.inBox || c == TC || c == UC
#                     gWeight = bubWeight * Factor
#                 else
#                     gWeight = G[1].weight[map.G] * G[c].weight[map.Gx] * Factor
#                 end

#                 if fast && ver4.level == 0
#                     pair = ver4.Tpair[map.ver]
#                     dT =
#                         varT[pair[INL]] - varT[pair[OUTL]] + varT[pair[INR]] -
#                         varT[pair[OUTR]]
#                     gWeight *= cos(2.0 * pi / Beta * dT)
#                     w = ver4.weight[ChanMap[c]]
#                 else
#                     w = ver4.weight[map.ver]
#                 end

#                 if c == T || c == TC
#                     w[DI] +=
#                         gWeight *
#                         (Lw[DI] * Rw[DI] * SPIN + Lw[DI] * Rw[EX] + Lw[EX] * Rw[DI])
#                     w[EX] += gWeight * Lw[EX] * Rw[EX]
#                 elseif c == U || c == UC
#                     w[DI] += gWeight * Lw[EX] * Rw[EX]
#                     w[EX] +=
#                         gWeight *
#                         (Lw[DI] * Rw[DI] * SPIN + Lw[DI] * Rw[EX] + Lw[EX] * Rw[DI])
#                 else
#                     # S channel,  see the note "code convention"
#                     w[DI] += gWeight * (Lw[DI] * Rw[EX] + Lw[EX] * Rw[DI])
#                     w[EX] += gWeight * (Lw[DI] * Rw[DI] + Lw[EX] * Rw[EX])
#                 end

#             end
#         end

#     end
# end

# function _expandBubble(children, text, style, bub::Bubble, parent)
#     push!(children, zeros(Int, 0))
#     @assert parent == length(children)
#     dict = Dict(
#         I => ("I", "yellow"),
#         T => ("T", "red"),
#         TC => ("Tc", "pink"),
#         U => ("U", "blue"),
#         UC => ("Uc", "navy"),
#         S => ("S", "green"),
#     )
#     push!(text, "$(dict[bub.chan][1])\n$(bub.Lver.loopNum)-$(bub.Rver.loopNum)")
#     push!(style, "fill:$(dict[bub.chan][2])")

#     current = length(children) + 1
#     push!(children[parent], current)
#     _expandVer4(children, text, style, bub.Lver, current) # left vertex 

#     current = length(children) + 1
#     push!(children[parent], current)
#     _expandVer4(children, text, style, bub.Rver, current) # right vertex
# end

# function _expandVer4(children, text, style, ver4::Ver4, parent)
#     push!(children, zeros(Int, 0))
#     @assert parent == length(children)
#     # println("Ver4: $(ver4.level), Bubble: $(length(ver4.bubble))")
#     if ver4.loopNum > 0
#         info = "O$(ver4.loopNum)\nT[$(length(ver4.Tpair))]\n"
#         for t in ver4.Tpair
#             info *= "[$(t[1]) $(t[2]) $(t[3]) $(t[4])]\n"
#         end
#     else
#         info = "O$(ver4.loopNum)"
#     end
#     push!(text, info)

#     ver4.inBox ? push!(style, "stroke-dasharray:3,2") : push!(style, "")

#     for bub in ver4.bubble
#         # println(bub.chan)
#         current = length(children) + 1
#         push!(children[parent], current)
#         _expandBubble(children, text, style, bub, current)
#     end
# end

# function visualize(ver4::Ver4)
#     children, text, style = (Vector{Vector{Int}}(undef, 0), [], [])
#     _expandVer4(children, text, style, ver4, 1)

#     # text = ["one\n(second line)", "2", "III", "four"]
#     # style = ["", "fill:red", "r:14", "opacity:0.7"]
#     # link_style = ["", "stroke:blue", "", "stroke-width:10px"]
#     tooltip = ["pops", "up", "on", "hover"]
#     t = D3Trees.D3Tree(
#         children,
#         text=text,
#         style=style,
#         tooltip=tooltip,
#         # link_style = link_style,
#         title="Vertex4 Tree",
#         init_expand=2,
#     )

#     D3Trees.inchrome(t)
# end