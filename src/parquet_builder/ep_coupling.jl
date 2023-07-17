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


"""
    ep_vertex4(para::DiagPara,
        extK = [DiagTree.getK(para.totalLoopNum, 1), DiagTree.getK(para.totalLoopNum, 2), DiagTree.getK(para.totalLoopNum, 3)],
        chan::AbstractVector = [PHr, PHEr, PPr, Alli],
        subdiagram = false;
        level = 1, name = :none, resetuid = false,
        subchannel::Symbol=:All, #:All, :W, :Lver3, :Rver3, :RPA
        blocks::ParquetBlocks=ParquetBlocks(),
        blockstoplevel::ParquetBlocks=blocks
        )

    Generate 4-vertex diagrams using Parquet Algorithm

# Arguments
- `para`            : parameters. It should provide internalLoopNum, interactionTauNum, firstTauIdx
- `extK`            : basis of external loops as a vector [left in, left out, right in, right out]. 
- `chan`            : vector of channels of the current 4-vertex. 
- `subdiagram`      : a sub-vertex or not
- `name`            : name of the vertex
- `level`           : level in the diagram tree
- `resetuid`        : restart uid count from 1
- `subchannel`      : :All, :W, :Lver3, :Rver3, :RPA to select all, W-interaction, left-vertex-correction, right-vertex-correction or RPA-interaction diagrams. See the remarks for more details.
- `blocks`          : building blocks of the Parquet equation. See the struct ParquetBlocks for more details.
- `blockstoplevel`  : building blocks of the Parquet equation at the toplevel. See the struct ParquetBlocks for more details.

# Output
- A DataFrame with fields :response, :type, :extT, :diagram, :hash

# Remarks
- ep diagram (if loopNum>0, then the right incoming Tau will be set to the last Tau for all diagrams):
|         |
Γ3 -------|
|         |

"""
function vertex4(para::DiagPara{W};
    extK=[DiagTree.getK(para.totalLoopNum, 1), DiagTree.getK(para.totalLoopNum, 2), DiagTree.getK(para.totalLoopNum, 3)],
    chan::AbstractVector=[PHr, PHEr, PPr, Alli], subdiagram=false,
    level=1, name=:none, resetuid=false,
    blocks::ParquetBlocks=ParquetBlocks()
) where {W}

    for k in extK
        @assert length(k) >= para.totalLoopNum "expect dim of extK>=$(para.totalLoopNum), got $(length(k))"
    end

    legK = [k[1:para.totalLoopNum] for k in extK[1:3]]
    push!(legK, legK[1] + legK[3] - legK[2])

    resetuid && uidreset()

    @assert para.totalTauNum >= maxVer4TauIdx(para) "Increase totalTauNum!\n$para"
    @assert para.totalLoopNum >= maxVer4LoopIdx(para) "Increase totalLoopNum\n$para"

    phi, ppi = blocks.phi, blocks.ppi
    phi_toplevel, ppi_toplevel = blockstoplevel.phi, blockstoplevel.ppi
    blockstoplevel = ParquetBlocks(phi=phi_toplevel, ppi=ppi_toplevel)

    @assert (PHr in phi) == false "PHi vertex is particle-hole irreducible, so that PHr channel is not allowed in $phi"
    @assert (PPr in ppi) == false "PPi vertex is particle-particle irreducible, so that PPr channel is not allowed in $ppi"
    @assert (PHr in phi_toplevel) == false "PHi vertex is particle-hole irreducible, so that PHr channel is not allowed in $phi_toplevel"
    @assert (PPr in ppi_toplevel) == false "PPi vertex is particle-particle irreducible, so that PPr channel is not allowed in $ppi_toplevel"

    loopNum = para.innerLoopNum
    # @assert loopNum >= 0

    ver4df = DataFrame(response=Response[], type=AnalyticProperty[], extT=Tuple{Int,Int,Int,Int}[], diagram=Diagram{W}[])

    if loopNum == 0
        if DirectOnly in para.filter
            permutation = [Di,]
        else
            permutation = [Di, Ex]
        end
        bareVer4(ver4df, para, legK, permutation)
    else # loopNum>0
        for c in chan
            if c == Alli
                continue
            end

            partition = orderedPartition(loopNum - 1, 4, 0)

            for p in partition
                if c == PHr || c == PHEr || c == PPr
                    bubble!(ver4df, para, legK, c, p, level, name, blocks, blockstoplevel, 1.0, subchannel)
                end
            end

            if (NoBubble in para.filter) && (c == PHr || c == PHEr)
                # add RPA bubble counter-diagram to remove the bubble
                RPA_chain!(ver4df, para, legK, c, level, name, -1.0)
            end
            # println(bub)
        end
        # # TODO: add envolpe diagrams
    end
    diags = ver4df.diagram
    @assert all(x -> x.id isa Ver4Id, diags) "not all id are Ver4Id! $diags"
    @assert all(x -> x.id.extK ≈ legK, diags) "not all extK are the same! $diags"

    # @assert isempty(diags) == false "got empty ver4! $chan with\n $para\n"
    if isempty(ver4df) == false
        ver4df = mergeby(W, ver4df, [:response, :type, :extT], name=name,
            getid=g -> Ver4Id(para, g[1, :response], g[1, :type], k=legK, t=g[1, :extT]) #generate id from the dataframe
        )
    end
    @assert all(x -> x[1] == para.firstTauIdx, ver4df.extT) "not all extT[1] are equal to the first Tau index $(para.firstTauIdx)! $ver4df"
    # println(typeof(groups))
    return ver4df
end