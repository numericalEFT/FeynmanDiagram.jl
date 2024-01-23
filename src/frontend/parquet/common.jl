
function build(para::DiagPara, extK=nothing, subdiagram=false; channels=[PHr, PHEr, PPr, Alli])
    if para.type == Ver4Diag
        if isnothing(extK)
            extK = [getK(para.totalLoopNum, 1), getK(para.totalLoopNum, 2), getK(para.totalLoopNum, 3)]
        end
        return vertex4(para, extK, subdiagram, channels=channels)
    elseif para.type == SigmaDiag
        if isnothing(extK)
            extK = getK(para.totalLoopNum, 1)
        end
        return sigma(para, extK, subdiagram)
    elseif para.type == PolarDiag
        if isnothing(extK)
            extK = getK(para.totalLoopNum, 1)
        end
        return polarization(para, extK, subdiagram)
    elseif para.type == Ver3Diag
        if isnothing(extK)
            extK = [getK(para.totalLoopNum, 1), getK(para.totalLoopNum, 2)]
        end
        return vertex3(para, extK, subdiagram, channels=channels)
    else
        error("not implemented!")
    end
end

function orderedPartition(_total, n, lowerbound=1)
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

"""
    function innerTauNum(type::DiagramType, innerLoopNum, interactionTauNum)
    
    internal imaginary-time degrees of freedom for a given diagram type and internal loop number.
    For the vertex functions (self-energy, polarization, vertex3, and vertex4), innerTauNum is equivalent to tauNum.
    For the Green function, tauNum = innerTauNum + external tauNum 
"""
function innerTauNum(type::DiagramType, innerLoopNum, interactionTauNum)
    if type == Ver4Diag
        return (innerLoopNum + 1) * interactionTauNum
    elseif type == SigmaDiag
        return innerLoopNum * interactionTauNum
    elseif type == GreenDiag
        return innerLoopNum * interactionTauNum
    elseif type == VacuumDiag
        return (innerLoopNum - 1) * interactionTauNum
    elseif type == PolarDiag
        return 1 + innerTauNum(Ver3Diag, innerLoopNum - 1, interactionTauNum)
    elseif type == Ver3Diag
        return 1 + innerTauNum(Ver4Diag, innerLoopNum - 1, interactionTauNum)
    else
        error("not implemented!")
    end
end

function interactionTauNum(hasTau::Bool, interactionSet)
    if hasTau == false
        return 0
    end
    for interaction in interactionSet
        if Dynamic in interaction.type
            return 2
        end
    end
    return 1
end

function firstTauIdx(type, offset::Int=0)
    if type == GreenDiag
        return 3 + offset
    elseif type == Ver3Diag
        return 1 + offset
    elseif type == PolarDiag
        return 1 + offset
    else
        return 1 + offset
    end
end

function firstLoopIdx(type, offset::Int=0)
    if type == Ver4Diag #three extK
        return 4 + offset
    elseif type == SigmaDiag #one extK
        return 2 + offset
    elseif type == GreenDiag #one extK
        return 2 + offset
    elseif type == PolarDiag #one extK
        return 2 + offset
    elseif type == Ver3Diag #two extK
        return 3 + offset
    elseif type == VacuumDiag #no extK
        return 1 + offset
    else
        error("not implemented!")
    end
end

function totalTauNum(type, innerLoopNum, interactionTauNum, offset::Int=0)
    return firstTauIdx(type, offset) + innerTauNum(type, innerLoopNum, interactionTauNum) - 1
end

function totalLoopNum(type, innerLoopNum, offset::Int=0)
    return firstLoopIdx(type, offset) + innerLoopNum - 1
end

function totalTauNum(para, type::Symbol=:none)
    return para.totalTauNum
    # if type == :Ver4
    #     return (para.internalLoopNum + 1) * para.interactionTauNum
    # else
    #     error("not implemented!")
    # end
end

function totalLoopNum(para, type::Symbol=:none)
    return para.totalLoopNum
end

function getK(loopNum::Int, loopIdx::Int)
    k = zeros(loopNum)
    k[loopIdx] = 1.0
    return k
end


function findFirstLoopIdx(partition, firstidx::Int)
    ## example: firstidx = 1
    # partition = [1, 1, 2, 1], then the loop partition = [1][2][34][5], thus firstTauIdx = [1, 2, 3, 5]
    # partition = [1, 0, 2, 0], then the loop partition = [1][][23][], thus firstTauIdx = [1, 2, 2, 4]
    # @assert length(partition) == length(isG)
    accumulated = accumulate(+, partition; init=firstidx) #  idx[i] = firstidx + p[1]+p[2]+...+p[i]
    firstLoopIdx = [firstidx,]
    append!(firstLoopIdx, accumulated[1:end-1])
    maxLoopIdx = accumulated[end] - 1
    return firstLoopIdx, maxLoopIdx
end

function findFirstTauIdx(partition::Vector{Int}, type::Vector{DiagramType}, firstidx::Int, _tauNum::Int)
    ## example: type =[Vertex4, GreenDiag, Vertex4, GreenDiag], firstidx = 1
    # n-loop G has n*_tauNum DOF, while n-loop ver4 has (n+1)*_tauNum DOF
    # partition = [1, 1, 2, 1], then the tau partition = [12][3][456][7], thus firstTauIdx = [1, 3, 4, 7]
    # partition = [1, 0, 2, 0], then the tau partition = [12][][345][], thus firstTauIdx = [1, 3, 3, 6]
    @assert length(partition) == length(type)
    @assert _tauNum >= 0
    taupartition = [innerTauNum(type[i], p, _tauNum) for (i, p) in enumerate(partition)]
    accumulated = accumulate(+, taupartition; init=firstidx) #  idx[i] = firstidx + p[1]+p[2]+...+p[i]
    firstTauidx = [firstidx,]
    append!(firstTauidx, accumulated[1:end-1])
    maxTauIdx = accumulated[end] - 1
    return firstTauidx, maxTauIdx
end

function allsame(df, name::Symbol)
    @assert all(x -> x == df[1, name], df[!, name]) "Not all rows of the $name field are the same.\n$df"
end
function allsame(df, names::Vector{Symbol})
    for name in names
        allsame(df, name)
    end
end
function allsametype(df, name::Symbol)
    @assert all(x -> typeof(x) == typeof(df[1, name]), df[!, name]) "Not all rows of the $name field are the same type.\n$df"
end
function allsametype(df, names::Vector{Symbol})
    for name in names
        allsametype(df, name)
    end
end

