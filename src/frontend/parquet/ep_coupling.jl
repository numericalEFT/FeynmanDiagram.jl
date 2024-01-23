"""
    function ep_coupling(para::DiagPara;
        extK=[getK(para.totalLoopNum, 1), getK(para.totalLoopNum, 2), getK(para.totalLoopNum, 3)],
        channels::AbstractVector=[PHr, PHEr, PPr, Alli],
        subdiagram=false,
        name=:none, resetuid=false,
        blocks::ParquetBlocks=ParquetBlocks()
    )

Generate electron-phonon 4-vertex diagrams using Parquet Algorithm. The right incoming Tau will be set to the last Tau for all diagrams
|         |
Γ3 -------|
|         |

# Arguments
- `para`            : parameters. It should provide internalLoopNum, interactionTauNum, firstTauIdx
- `extK`            : basis of external loops as a vector [left in, left out, right in, right out]. 
- `channels`        : vector of channels in the left Γ3 diagrams. 
- `subdiagram`      : a sub-vertex or not
- `name`            : name of the vertex
- `resetuid`        : restart uid count from 1
- `blocks`          : building blocks of the Parquet equation. See the struct ParquetBlocks for more details.

# Output
- A DataFrame with fields :response, :extT, :diagram, :hash. 

# Output
- A DataFrame with fields :response, :type, :extT, :diagram, :hash
"""
function ep_coupling(para::DiagPara;
    extK=[getK(para.totalLoopNum, 1), getK(para.totalLoopNum, 2), getK(para.totalLoopNum, 3)],
    channels::AbstractVector=[PHr, PHEr, PPr, Alli],
    subdiagram=false,
    name=:none, resetuid=false,
    blocks::ParquetBlocks=ParquetBlocks()
)

    @warn("ep vertex4 breaks SU(2) spin symmetry!")
    if NoBubble in para.filter
        @warn("the RPA chain counterterms from the renormalization of the outgoing interaction leg in the ep vertex4 have not yet been implemented!")
    end

    for k in extK
        @assert length(k) >= para.totalLoopNum "expect dim of extK>=$(para.totalLoopNum), got $(length(k))"
    end

    legK = [k[1:para.totalLoopNum] for k in extK[1:3]]
    push!(legK, legK[1] + legK[3] - legK[2])

    resetuid && IR.uidreset()

    @assert para.totalTauNum >= maxVer4TauIdx(para) "Increase totalTauNum!\n$para"
    @assert para.totalLoopNum >= maxVer4LoopIdx(para) "Increase totalLoopNum\n$para"

    loopNum = para.innerLoopNum
    # @assert loopNum >= 0

    ver4df = DataFrame(response=Response[], type=AnalyticProperty[], extT=Tuple{Int,Int,Int,Int}[], diagram=Graph{Ftype,Wtype}[])

    partition = orderedPartition(loopNum - 1, 4, 0)

    for p in partition
        if p[3] == 0 #oL, oG0, oR, oGx
            ep_bubble!(ver4df, para, legK, channels, p, name, blocks, 1.0)
        end
    end

    if NoBubble in para.filter
        # add RPA bubble counter-diagram to remove the bubble
        ep_RPA_chain!(ver4df, para, legK, name, -1.0)
    end
    # println(bub)

    diags = ver4df.diagram
    @assert all(x -> x.properties isa Ver4Id, diags) "not all id are Ver4Id! $diags"
    @assert all(x -> x.properties.extK ≈ legK, diags) "not all extK are the same! $diags"

    ver4df = merge_vertex4(para, ver4df, name, legK)
    @assert all(x -> x[1] == para.firstTauIdx, ver4df.extT) "not all extT[1] are equal to the first Tau index $(para.firstTauIdx)! $ver4df"
    # @assert all(x -> x[3] == para.totalTauNum, ver4df.extT) "not all extT[3] are equal to the first Tau index $(para.totalTauNum)! $ver4df"
    # @assert all(x -> x[4] == para.totalTauNum, ver4df.extT) "not all extT[4] are equal to the first Tau index $(para.totalTauNum)! $ver4df"
    # println(typeof(groups))
    return ver4df
end

function ep_bubble!(ver4df::DataFrame, para::DiagPara, legK, chans::Vector{TwoBodyChannel}, partition::Vector{Int}, name::Symbol, blocks::ParquetBlocks,
    extrafactor=1.0)

    @assert partition[3] == 0

    TauNum = interactionTauNum(para) # maximum tau number for each bare interaction
    oL, oG0, oR, oGx = partition[1], partition[2], partition[3], partition[4]
    if isValidG(para.filter, oG0) == false || isValidG(para.filter, oGx) == false
        return
    end

    #the first loop idx is the inner loop of the bubble!
    LoopIdx = para.firstLoopIdx
    idx, maxLoop = findFirstLoopIdx(partition, LoopIdx + 1)
    LfirstLoopIdx, G0firstLoopIdx, RfirstLoopIdx, GxfirstLoopIdx = idx
    @assert maxLoop == maxVer4LoopIdx(para)

    type = [Ver4Diag, GreenDiag, Ver4Diag, GreenDiag]
    idx, maxTau = findFirstTauIdx(partition, type, para.firstTauIdx, TauNum)
    LfirstTauIdx, G0firstTauIdx, RfirstTauIdx, GxfirstTauIdx = idx
    @assert maxTau == maxVer4TauIdx(para) "Partition $partition with tauNum configuration $idx. maxTau = $maxTau, yet $(maxTauIdx(para)) is expected!"

    lPara = reconstruct(para, type=Ver4Diag, innerLoopNum=oL, firstLoopIdx=LfirstLoopIdx, firstTauIdx=LfirstTauIdx)
    rPara = reconstruct(para, type=Ver4Diag, innerLoopNum=oR, firstLoopIdx=RfirstLoopIdx, firstTauIdx=RfirstTauIdx)
    gxPara = reconstruct(para, type=GreenDiag, innerLoopNum=oGx, firstLoopIdx=GxfirstLoopIdx, firstTauIdx=GxfirstTauIdx)
    g0Para = reconstruct(para, type=GreenDiag, innerLoopNum=oG0, firstLoopIdx=G0firstLoopIdx, firstTauIdx=G0firstTauIdx)

    LLegK, K, RLegK, Kx = legBasis(PHr, legK, LoopIdx)

    Lver = vertex4(lPara, LLegK, true; channels=chans, name=:Γf, blocks=blocks)
    isempty(Lver) && return

    Rver = DataFrame(response=Response[], type=AnalyticProperty[], extT=Tuple{Int,Int,Int,Int}[], diagram=Graph{Ftype,Wtype}[])
    bareVer4(Rver, rPara, RLegK, [Di,], false) # the external tau is right aligned
    Rver = merge_vertex4(rPara, Rver, :bare, RLegK)

    @assert isempty(Rver) == false

    for ldiag in Lver.diagram
        for rdiag in Rver.diagram
            LvT, RvT = ldiag.properties.extT, rdiag.properties.extT
            extT, G0T, GxT = tauBasis(PHr, ldiag.properties.extT, rdiag.properties.extT)
            g0 = green(g0Para, K, G0T, true, name=:G0, blocks=blocks)
            gx = green(gxPara, Kx, GxT, true, name=:Gx, blocks=blocks)
            @assert g0 isa Graph && gx isa Graph
            # append!(diag, bubble2diag(para, chan, ldiag, rdiag, legK, g0, gx, extrafactor))
            bubble2diag!(ver4df, para, PHr, ldiag, rdiag, legK, g0, gx, extrafactor)
        end
    end
    return
end

function ep_RPA_chain!(ver4df::DataFrame, para::DiagPara, legK, name::Symbol, extrafactor)
    new_filter = union(union(para.filter, Girreducible), DirectOnly)
    para_rpa = reconstruct(para, filter=new_filter)
    blocks = ParquetBlocks(; phi=[], ppi=[], Γ4=[PHr,])
    ep_bubble!(ver4df, para_rpa, legK, [PHr,], [0, 0, para.innerLoopNum - 1, 0], Symbol("$(name)_ep_RPA_CT"), blocks, extrafactor)
    return
end