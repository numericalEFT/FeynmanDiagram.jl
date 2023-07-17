"""
    function ep_vertex4(para, extK = [DiagTree.getK(para.totalLoopNum, 1), DiagTree.getK(para.totalLoopNum, 2)],
        subdiagram = false; name = :Γ3, chan = [PHr, PHEr, PPr, Alli], resetuid = false, 
        blocks::ParquetBlocks=ParquetBlocks()
        )

    Generate 3-vertex diagrams using Parquet Algorithm.
    With imaginary-time variables, all vertex3 generated has the same bosonic Tidx ``extT[1]=para.firstTauIdx`` and the incoming fermionic Tidx ``extT[2]=para.firstTauIdx+1``.

#Arguments
- `para`            : parameters. It should provide internalLoopNum, interactionTauNum, firstTauIdx
- `extK`            : basis of external loops as a vector [bosonic leg, fermionic in, fermionic out]. 
- `subdiagram`      : a sub-vertex or not
- `name`            : name of the vertex
- `chan`            : vector of channels of the current 4-vertex. 
- `resetuid`        : restart uid count from 1
- `blocks`          : building blocks of the Parquet equation. See the struct ParquetBlocks for more details.

# Output
- A DataFrame with fields :response, :extT, :diagram, :hash. 
"""
function ep_vertex4(para::DiagPara{WW};
    extK=[DiagTree.getK(para.totalLoopNum, 1), DiagTree.getK(para.totalLoopNum, 2), DiagTree.getK(para.totalLoopNum, 3)],
    subdiagram=false;
    name=:ep_Γ4,
    chan=[PHr, PHEr, PPr, Alli],
    resetuid=false,
    blocks::ParquetBlocks=ParquetBlocks()
) where {WW}

    resetuid && uidreset()
    @assert para.type == Ver4Diag
    @assert para.innerLoopNum >= 1 "Only generates vertex corrections with more than one internal loops."
    for k in extK
        @assert length(k) >= para.totalLoopNum "expect dim of extK>=$(para.totalLoopNum), got $(length(k))"
    end

    extK = [k[1:para.totalLoopNum] for k in extK[1:3]]
    push!(extK, extK[1] + extK[3] - extK[2])
    Kin, Kout = extK[1], extK[2]
    q = Kin - Kout
    @assert ((q ≈ KinL) == false) && ((q ≈ KoutL) == false) "The bosonic q cann't be same as the fermionic k. Ohterwise the proper diagram check will fail!"

    para = _properVer3Para(para, q)

    t0 = para.totalTauNum # the right incoming/outing Tau will be set to the last Tau index
    vertex3 = DataFrame(response=Response[], extT=Tuple{Int,Int,Int}[], diagram=Diagram{WW}[])

    K = zero(q)
    LoopIdx = para.firstLoopIdx
    K[LoopIdx] = 1.0
    # extT = (t0, t0 + 1)
    legK = [Kin, Kout, K, K .+ q]
    RlegK = [Kin, Kout, K .+ q, K]

    ######################## Π0 = GG #########################################
    for (oVer4, oGin, oGout, oRVer4) in orderedPartition(para.innerLoopNum - 1, 4, 0)
        # ! Vertex4 must be in the first place, because we want to make sure that the TinL of the vertex4 start with t0+1
        if oRVer4 != 0
            continue
        end

        idx, maxLoop = findFirstLoopIdx([oVer4, oGin, oGout, oRVer4], LoopIdx + 1)
        @assert maxLoop <= para.totalLoopNum "maxLoop = $maxLoop > $(para.totalLoopNum)"
        Ver4Kidx, GinKidx, GoutKidx, RVer4Kidx = idx

        ver4t0 = para.firstTauIdx
        idx, maxTau = findFirstTauIdx([oVer4, oGin, oGout, oRVer4], [Ver4Diag, GreenDiag, GreenDiag, Ver4Diag], ver4t0, interactionTauNum(para))
        @assert maxTau <= para.totalTauNum "maxTau = $maxTau > $(para.totalTauNum)"
        Ver4Tidx, GinTidx, GoutTidx, RVer4Tidx = idx

        if isValidG(para.filter, oGin) && isValidG(para.filter, oGout)
            paraGin = reconstruct(para, type=GreenDiag, innerLoopNum=oGin,
                firstLoopIdx=GinKidx, firstTauIdx=GinTidx)
            paraGout = reconstruct(para, type=GreenDiag, innerLoopNum=oGout,
                firstLoopIdx=GoutKidx, firstTauIdx=GoutTidx)
            paraVer4 = reconstruct(para, type=Ver4Diag, innerLoopNum=oVer4,
                firstLoopIdx=Ver4Kidx, firstTauIdx=Ver4Tidx)
            paraRVer4 = reconstruct(para, type=Ver4Diag, innerLoopNum=oRVer4,
                firstLoopIdx=RVer4Kidx, firstTauIdx=RVer4Tidx)

            ver4 = vertex4(paraVer4, legK, chan, true; blocks=blocks)

            rver4 = vertex4(paraRVer4, RlegK, [], true; blocks=blocks)

            if isnothing(ver4) || isempty(ver4)
                continue
            end

            @assert isnothing(rver4) == false
            @assert isempty(rver4) == false

            if para.hasTau
                @assert all(x -> x[INL] == ver4t0, ver4.extT) "The TinL of the inner Γ4 must be firstTauIdx+1"
            end

            if para.hasTau
                @assert all(x -> x[INL] == t0, rver4.extT) "The TinR of the right Γ4 must be the last Tau index"
            end

            #transform extT coloum into extT for Vertex4 and the extT for Gin and Gout
            df = transform(ver4, :extT => ByRow(x -> [(t0, x[INL], x[OUTL]), (t0, x[INR]), (x[OUTR], t0)]) => [:extT, :GinT, :GoutT])

            groups = mergeby(WW, df, [:response, :GinT, :GoutT, :extT], operator=Sum())

            for v4 in eachrow(groups)
                response = v4.response
                @assert response == UpUp || response == UpDown
                #type: Instant or Dynamic
                ver3id = Ver3Id(para, response, k=extK, t=v4.extT)
                gin = green(paraGin, K, v4.GinT, true, name=:Gin, blocks=blocks)
                gout = green(paraGout, K .+ q, v4.GoutT, true, name=:Gout, blocks=blocks)
                @assert gin isa Diagram && gout isa Diagram

                ver3diag = Diagram{WW}(ver3id, Prod(), [gin, gout, v4.diagram], name=name)
                push!(vertex3, (response=response, extT=v4.extT, diagram=ver3diag))
            end

        end
    end
    # println(vertex3)

    if isempty(vertex3) == false
        # Factor = 1 / (2π)^para.loopDim
        Factor = 1.0
        vertex3 = mergeby(WW, vertex3, [:response, :extT]; name=name, factor=Factor,
            getid=g -> Ver3Id(para, g[1, :response], k=extK, t=g[1, :extT])
        )
    end
    return vertex3
end

function _properVer3Para(p::DiagPara{W}, q) where {W}
    # ############ reset transferLoop to be q ################
    if Proper in p.filter
        if (length(p.transferLoop) != length(q)) || (!(p.transferLoop ≈ q)) #first check if the dimension is wrong
            return derivepara(p, transferLoop=q)
        end
    end
    return p
end
