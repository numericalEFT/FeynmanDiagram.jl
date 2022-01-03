import ..Filter
import ..Wirreducible  #remove all polarization subdiagrams
import ..Girreducible  #remove all self-energy inseration
import ..NoHatree
import ..NoFock
import ..NoBubble  # true to remove all bubble subdiagram
import ..Proper  #ver4, ver3, and polarization diagrams may require to be irreducible along the transfer momentum/frequency

import ..DiagramType
import ..GreenDiag
import ..SigmaDiag
import ..PolarDiag
import ..Ver3Diag
import ..Ver4Diag

import ..GenericPara

import ..innerTauNum

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

function findFirstLoopIdx(partition, firstidx::Int)
    ## example: firstidx = 1
    # partition = [1, 1, 2, 1], then the loop partition = [1][2][34][5], thus firstTauIdx = [1, 2, 3, 5]
    # partition = [1, 0, 2, 0], then the loop partition = [1][][23][], thus firstTauIdx = [1, 2, 2, 4]
    # @assert length(partition) == length(isG)
    accumulated = accumulate(+, partition; init = firstidx) #  idx[i] = firstidx + p[1]+p[2]+...+p[i]
    firstLoopIdx = [firstidx,]
    append!(firstLoopIdx, accumulated[1:end-1])
    maxLoopIdx = accumulated[end] - 1
    return firstLoopIdx, maxLoopIdx
end

function findFirstTauIdx(partition::Vector{Int}, diagType::Vector{DiagramType}, firstidx::Int, _tauNum::Int)
    ## example: diagType =[Vertex4, GreenDiagram, Vertex4, GreenDiagram], firstidx = 1
    # n-loop G has n*_tauNum DOF, while n-loop ver4 has (n+1)*_tauNum DOF
    # partition = [1, 1, 2, 1], then the tau partition = [12][3][456][7], thus firstTauIdx = [1, 3, 4, 7]
    # partition = [1, 0, 2, 0], then the tau partition = [12][][345][], thus firstTauIdx = [1, 3, 3, 6]
    @assert length(partition) == length(diagType)
    taupartition = [innerTauNum(diagType[i], p, _tauNum) for (i, p) in enumerate(partition)]
    accumulated = accumulate(+, taupartition; init = firstidx) #  idx[i] = firstidx + p[1]+p[2]+...+p[i]
    firstTauidx = [firstidx,]
    append!(firstTauidx, accumulated[1:end-1])
    maxTauIdx = accumulated[end] - 1
    return firstTauidx, maxTauIdx
end

function newDiagTree(para, name::Symbol = :none)
    weightType = para.weightType
    Kpool = DiagTree.LoopPool(:K, para.loopDim, para.totalLoopNum, Float64)
    nodeParaType = Vector{Int}

    if para.interactionTauNum == 2
        Gpool = DiagTree.propagatorPool(:Gpool, weightType)
        Vpool = DiagTree.propagatorPool(:Vpool, weightType)
        Wpool = DiagTree.propagatorPool(:Wpool, weightType)
        return DiagTree.Diagrams(Kpool, (Gpool, Vpool, Wpool), weightType, nodeParaType = nodeParaType, name = name)
    elseif para.interactionTauNum == 1 || para.interactionTauNum == 0
        Gpool = DiagTree.propagatorPool(:Gpool, weightType)
        Vpool = DiagTree.propagatorPool(:Vpool, weightType)
        return DiagTree.Diagrams(Kpool, (Gpool, Vpool), weightType, nodeParaType = nodeParaType, name = name)
    else
        error("not implemented!")
    end
end