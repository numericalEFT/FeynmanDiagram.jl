function Vertex3(para, extK, subdiagram = false; name = :Γ3, chan = [PHr, PHEr, PPr, Alli])
    (subdiagram == false) && uidreset()
    @assert para.diagType == Ver3Diag
    @assert para.innerLoopNum >= 1 "Only generates vertex corrections with more than one internal loops."
    for k in extK
        @assert length(k) == para.totalLoopNum
    end

    q, Kin = extK[1], extK[2]
    Kout = length(extK) == 3 ? extK[3] : Kin .- q
    @assert ((q ≈ Kin) == false) && ((q ≈ Kout) == false) "The bosonic q cann't be same as the fermionic k. Ohterwise the proper diagram check will fail!"
    extK = [q, Kin, Kout]

    if Proper in para.filter
        para = reconstruct(para, transferLoop = q)
    end

    if (para.extra isa ParquetBlocks) == false
        parquetblocks = ParquetBlocks(phi = [PPr, PHEr], ppi = [PHr, PHEr], Γ4 = [PPr, PHr, PHEr])
        para = reconstruct(para, extra = parquetblocks)
    end

    K = zero(q)
    LoopIdx = para.firstLoopIdx
    K[LoopIdx] = 1.0
    t0 = para.firstTauIdx
    # extT = (t0, t0 + 1)
    legK = [Kin, Kout, K, K .+ q]

    vertex3 = DataFrame()

    ######################## Π0 = GG #########################################
    for (oGin, oGout, oVer4) in orderedPartition(para.innerLoopNum - 1, 3, 0)

        idx, maxLoop = findFirstLoopIdx([oGin, oGout, oVer4], LoopIdx + 1)
        @assert maxLoop <= para.totalLoopNum "maxLoop = $maxLoop > $(para.totalLoopNum)"
        GinKidx, GoutKidx, Ver4Kidx = idx

        idx, maxTau = findFirstTauIdx([oGin, oGout, oVer4], [GreenDiag, GreenDiag, Ver4Diag], para.firstTauIdx, para.interactionTauNum)
        @assert maxTau <= para.totalTauNum
        GinTidx, GoutTidx, Ver4Tidx = idx

        paraGin = reconstruct(para, diagType = GreenDiag, innerLoopNum = oGin,
            firstLoopIdx = GinKidx, firstTauIdx = GinTidx)
        paraGout = reconstruct(para, diagType = GreenDiag, innerLoopNum = oGout,
            firstLoopIdx = GoutKidx, firstTauIdx = GoutTidx)
        paraVer4 = reconstruct(para, diagType = Ver4Diag, innerLoopNum = oVer4,
            firstLoopIdx = Ver4Kidx, firstTauIdx = Ver4Tidx)

        if isValidG(paraGin) && isValidG(paraGout)
            ver4 = buildVer4(paraVer4, legK, chan, true)

            #transform extT coloum into extT for Vertex4 and the extT for Gin and Gout
            df = transform(ver4, :extT => ByRow(x -> [(t0, x[INL], x[OUTL]), (t0, x[INR]), (x[OUTR], t0)]) => [:extT, :GinT, :GoutT])

            groups = mergeby(df, [:response, :GinT, :GoutT, :extT], operator = Sum())

            for v4 in eachrow(groups)
                response = v4[:response]
                @assert response == UpUp || response == UpDown
                #type: Instant or Dynamic
                ver3id = Ver3Id(para, response, k = extK, t = v4[:extT])
                gin = buildG(paraGin, K, v4[:GinT], true, name = :Gin)
                gout = buildG(paraGout, K .+ q, v4[:GoutT], true, name = :Gout)
                @assert gin isa Diagram && gout isa Diagram

                ver3diag = Diagram(ver3id, Prod(), [gin, gout, v4[:diagram]], name = name)
                push!(vertex3, (response = response, extT = v4[:extT], diagram = ver3diag))
            end
        end
    end

    ver3 = mergeby(vertex3, [:response, :extT]; name = name,
        getid = g -> Ver3Id(para, g[1, :response], k = extK, t = g[1, :extT])
    )
    return ver3
end