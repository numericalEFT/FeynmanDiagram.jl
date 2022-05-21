
"""
    function polarization(para, extK = DiagTree.getK(para.totalLoopNum, 1), subdiagram = false; name = :Π, resetuid = false)

    Generate polarization diagrams using Parquet Algorithm.

# Arguments
- `para`            : parameters. It should provide internalLoopNum, interactionTauNum, firstTauIdx
- `extK`            : basis of external loop. 
- `subdiagram`      : a sub-vertex or not
- `name`            : name of the vertex
- `resetuid`        : restart uid count from 1

# Output
- A DataFrame with fields `:response`, `:diagram`, `:hash`. 
- All polarization share the same external Tau index. With imaginary-time variables, they are extT = (para.firstTauIdx, para.firstTauIdx+1)
"""
function polarization(para, extK=DiagTree.getK(para.totalLoopNum, 1), subdiagram=false; name=:Π, resetuid=false)
    resetuid && uidreset()
    @assert para.diagType == PolarDiag
    @assert para.innerLoopNum >= 1
    # @assert length(extK) == para.totalLoopNum
    @assert length(extK) >= para.totalLoopNum "expect dim of extK>=$(para.totalLoopNum), got $(length(extK))"
    extK = extK[1:para.totalLoopNum]

    #polarization diagram is always proper
    para = reconstruct(para, filter=union(Proper, para.filter), transferLoop=extK)

    if (para.extra isa ParquetBlocks) == false
        para = reconstruct(para, extra=ParquetBlocks())
    end

    K = zero(extK)
    LoopIdx = para.firstLoopIdx
    K[LoopIdx] = 1.0
    @assert (K ≈ extK) == false
    t0 = para.firstTauIdx
    extT = para.hasTau ? (t0, t0 + 1) : (t0, t0)
    legK = [extK, K, K .- extK]

    polar = DataFrame()

    for (oVer3, oGin, oGout) in orderedPartition(para.innerLoopNum - 1, 3, 0)
        # ! Vertex3 must be in the first place, because we want to make sure that the bosonic extT of the vertex3 start with t0+1

        idx, maxLoop = findFirstLoopIdx([oVer3, oGin, oGout], LoopIdx + 1) # GGΓ3 consumes one internal loop
        @assert maxLoop <= para.totalLoopNum "maxLoop = $maxLoop > $(para.totalLoopNum)"
        Ver3Kidx, GinKidx, GoutKidx = idx

        if isValidG(para.filter, oGin) && isValidG(para.filter, oGout)
            if oVer3 == 0
                ######################## Π0 = GG #########################################
                gt0 = para.hasTau ? extT[2] + 1 : extT[1]
                idx, maxTau = findFirstTauIdx([oGin, oGout], [GreenDiag, GreenDiag], gt0, para.interactionTauNum)
                @assert maxTau <= para.totalTauNum "maxTau = $maxTau > $(para.totalTauNum)"
                GinTidx, GoutTidx = idx

                paraGin = reconstruct(para, diagType=GreenDiag, innerLoopNum=oGin,
                    firstLoopIdx=GinKidx, firstTauIdx=GinTidx)
                paraGout = reconstruct(para, diagType=GreenDiag, innerLoopNum=oGout,
                    firstLoopIdx=GoutKidx, firstTauIdx=GoutTidx)

                response = UpUp
                polarid = PolarId(para, response, k=extK, t=extT)
                gin = green(paraGin, K, (extT[1], extT[2]), true, name=:Gin)
                gout = green(paraGout, K .- extK, (extT[2], extT[1]), true, name=:Gout)
                @assert gin isa Diagram && gout isa Diagram "$gin or $gout is not a single diagram"

                sign = para.isFermi ? -1.0 : 1.0
                polardiag = Diagram(polarid, Prod(), [gin, gout], name=name, factor=sign)
                push!(polar, (response=response, extT=extT, diagram=polardiag))
            else
                ##################### composite polarization #####################################
                idx, maxTau = findFirstTauIdx([oVer3, oGin, oGout], [Ver3Diag, GreenDiag, GreenDiag], extT[2], para.interactionTauNum)
                @assert maxTau <= para.totalTauNum "maxTau = $maxTau > $(para.totalTauNum)"
                Ver3Tidx, GinTidx, GoutTidx = idx

                paraGin = reconstruct(para, diagType=GreenDiag, innerLoopNum=oGin,
                    firstLoopIdx=GinKidx, firstTauIdx=GinTidx)
                paraGout = reconstruct(para, diagType=GreenDiag, innerLoopNum=oGout,
                    firstLoopIdx=GoutKidx, firstTauIdx=GoutTidx)

                paraVer3 = reconstruct(para, diagType=Ver3Diag, innerLoopNum=oVer3,
                    firstLoopIdx=Ver3Kidx, firstTauIdx=Ver3Tidx)
                ver3 = vertex3(paraVer3, legK, true)
                if isnothing(ver3) || isempty(ver3)
                    continue
                end
                if para.hasTau
                    @assert all(x -> x[1] == extT[2], ver3[:, :extT]) "The bosonic T must be firstTauIdx+1 if hasTau\n$ver3"
                    @assert all(x -> x[2] == ver3[1, :extT][2], ver3[:, :extT]) "The TinL must be firstTauIdx+2 if hasTau\n$ver3"
                end

                #transform extT coloum into extT for Vertex4 and the extT for Gin and Gout
                df = transform(ver3, :extT => ByRow(x -> [extT, (extT[1], x[2]), (x[3], extT[1])]) => [:extT, :GinT, :GoutT])

                groups = mergeby(df, [:response, :GinT, :GoutT, :extT], operator=Sum())

                for v3 in eachrow(groups)
                    response = v3[:response]
                    @assert response == UpUp || response == UpDown
                    #type: Instant or Dynamic
                    polarid = PolarId(para, response, k=extK, t=v3[:extT])
                    gin = green(paraGin, K, v3[:GinT], true, name=:Gin)
                    gout = green(paraGout, K .- extK, v3[:GoutT], true, name=:Gout)
                    @assert gin isa Diagram && gout isa Diagram

                    polardiag = Diagram(polarid, Prod(), [gin, gout, v3[:diagram]], name=name)
                    push!(polar, (response=response, extT=v3[:extT], diagram=polardiag))
                end
            end
        end
    end

    # for (oGinL, oGoutL, oGinR, oGoutR, ver4) in orderedPartition(para.innerLoopNum - 1, 5, 0)
    # end
    if isempty(polar)
        return DataFrame(response=[], extT=[], diagram=[])
    end

    # legK = [extK, K, K, extK]
    Factor = 1 / (2π)^para.loopDim
    polar = mergeby(polar, [:response, :extT]; name=name, factor=Factor,
        getid=g -> PolarId(para, g[1, :response], k=extK, t=extT)
    )
    return polar
end