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

import ..ChargeCharge
import ..SpinSpin
import ..InteractionName

import ..Instant
import ..Dynamic
import ..InteractionType

import ..Interaction

import ..GenericPara

import ..innerTauNum

# mutable struct Weight <: FieldVector{2,Float64}
#     d::Float64
#     e::Float64
#     Weight() = new(0.0, 0.0)
#     Weight(d, e) = new(d, e)
# end

# const Base.zero(::Type{Weight}) = Weight(0.0, 0.0)
# const Base.abs(w::Weight) = abs(w.d) + abs(w.e) # define abs(Weight)

mutable struct Weight{T} <: FieldVector{2,T}
    d::T
    e::T
    Weight{T}() where {T} = new{T}(T(0), T(0))
    Weight(d::T, e::T) where {T} = new{T}(d::T, e::T)
    Weight{T}(d, e) where {T} = new{T}(T(d), T(e))
end

const Base.zero(::Type{Weight{T}}) where {T} = Weight{T}(0, 0)
const Base.abs(w::Weight) = abs(w.d) + abs(w.e) # define abs(Weight)

mutable struct CompositeWeight{T}
    interaction::Interaction
    instant::Weight{T}
    dynamic::Weight{T}
    d_instant::Weight{T}
    d_dynamic::Weight{T}
end


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

# struct ParameterizedComponent
#     component::Component
#     para::Any
# end
# const ComponentExtT = Tuple{Component,Tuple{Int,Int}}

# function connectComponentsbyGreen(diag, originalPara::GenericPara, name, loopBasis, extT::Tuple{Int,Int}, loopNumofG::Vector{Int},
#     diagType::Vector{DiagramType}, componentsVector::Vector{Vector{ComponentExt}}, factor = 1.0; para = extT)
#     #each component.para must be (tin, tout) or [tin, tout]
#     @assert length(loopNumofG) == (length(componentsVector) + 1)

#     for loop in loopNumofG
#         if isValidG(para.filter, loop) == false
#             @error("Some of the Green's function doesn't exist in the loopNum list $loopNumofG")
#         end
#     end
#     for c in componentsVector
#         @assert isempty(c) "Some of the components are empty!$componentsVector"
#     end

#     nodes = []
#     for configuration in Iterators.product(Tuple(componentsVector)...)
#         components = [ct[1] for ct in configuration]
#         _extT = [ct[2] for ct in configuration]

#         ########## prepare G extT ##################
#         GextT = [(extT[1], _extT[1][1]),]
#         for i in 1:length(_extT)-1
#             push!(GextT, (_extT[i][2], _extT[i+1][1]))
#         end
#         push!(GextT, (_extT[end][2], extT[2]))

#         for tpair in GextT
#             push!(components, DiagTree.addpropagator!(diag, :Gpool, 0, :G; site = tpair, loop = loopBasis))
#         end
#         node = DiagTree.addnode!(diag, MUL, name, components, factor; para = para)
#         @assert node.index > 0
#         push!(nodes, node)
#     end
#     n = DiagTree.addnode!(diag, ADD, name, nodes, factor; para = para)
#     @assert n.index > 0
#     return n
# end
