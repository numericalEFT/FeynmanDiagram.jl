using FeynmanDiagram
import FeynmanDiagram.ComputationalGraphs as IR
using MCIntegration, Lehmann
using Random, LinearAlgebra

const TAU_CUTOFF = 1e-10
inds = [12, 33, 37, 82, 83, 88, 102, 123, 127, 172, 173, 178]

function main()
    dim = 3
    kF = 1.919
    β = 3.0
    # para = Parquet.DiagPara(type=Parquet.Ver4Diag, innerLoopNum=4)

    # partition = [(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0), (3, 1, 0), (3, 0, 1), (4, 0, 0)]
    partition = [(4, 0, 0)]
    randseed = 1234
    λ = 1.2
    # MaxLoopNum = maximum([p[1] for p in partition]) + 3
    MaxLoopNum = 7
    Random.seed!(randseed)

    # FeynGraphs = diagdict_parquet(:vertex4, partition)
    # ver4df = Parquet.vertex4(para)
    diags, extT, responses = GV.eachorder_ver4diag(4)

    # diags = ver4df.diagram
    IR.optimize!(diags)
    # extT_labels = Vector{Vector{Int}}[]
    # spin_conventions = Vector{FrontEnds.Response}[]
    leaf_maps = Vector{Dict{Int,Graph}}()

    # for p in partition
    # push!(extT_labels, FeynGraphs[p][2])
    # push!(spin_conventions, FeynGraphs[p][3])
    # end

    # for (i, key) in enumerate(partition)
    # _, leafmap = Compilers.compile(FeynGraphs[key][1])
    # push!(leaf_maps, leafmap)
    # end
    _, leafmap = Compilers.compile(diags)
    push!(leaf_maps, leafmap)

    leafStat, loopBasis, leafval_map = _leafstates(leaf_maps, MaxLoopNum)

    momLoopPool = FrontEnds.LoopPool(:K, dim, loopBasis)
    # root = zeros(Float64, maximum(length.(extT_labels)))
    varK = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=3)
    varK.data[:, 1] .= [kF, 0.0, 0.0]
    varK.data[:, 2] .= [kF, 0.0, 0.0]
    varK.data[:, 3] .= 0.0
    varK.data[:, 4:9] = rand(Float64, (dim, 6))
    varT = MCIntegration.Continuous(0.0, β, offset=1)
    varT.data[1] = 0.0
    # varT.data[2:end] .= rand(16) * β


    leafval, leafType, leafOrders, leafτ_i, leafτ_o, leafMomIdx = leafStat

    FrontEnds.update(momLoopPool, varK.data[:, 1:MaxLoopNum])

    for (idx, p) in enumerate(partition)
        for (i, lftype) in enumerate(leafType[idx])
            if lftype == 0
                continue
            elseif lftype == 1 #fermionic 
                τ = varT[leafτ_o[idx][i]] - varT[leafτ_i[idx][i]]
                # kq = varK.data[:, leafMomIdx[idx][i]]
                kq = FrontEnds.loop(momLoopPool, leafMomIdx[idx][i])
                ϵ = dot(kq, kq) - kF^2
                order = leafOrders[idx][i][1]
                leafval[idx][i] = green_derive(τ, ϵ, β, order)
            elseif lftype == 2 #bosonic 
                # kq = varK.data[:, leafMomIdx[idx][i]]
                kq = FrontEnds.loop(momLoopPool, leafMomIdx[idx][i])
                order = leafOrders[idx][i][2]
                # leafval[idx][i] = Propagator.interaction_derive(τ1, τ2, kq, para, idorder; idtype=Instant, tau_num=1)
                invK = 1.0 / (dot(kq, kq) + λ)
                leafval[idx][i] = 8π / invK * (λ * invK)^order
            else
                error("this leaftype $lftype not implemented!")
            end
        end
        # graphfuncs! = funcGraphs![idx]
        # graphfuncs!(root, leafval[idx])
        for g in diags
            IR.eval!(g, leafval_map[idx], leafval[idx])
        end
    end

    return diags
    # return [FeynGraphs[p][1] for p in partition]
end

function green_derive(τ, ϵ, β, order)
    if order == 0
        result = green(τ, ϵ, β)
    elseif order == 1
        result = -Spectral.kernelFermiT_dω(τ, ϵ, β)
    elseif order == 2
        result = Spectral.kernelFermiT_dω2(τ, ϵ, β) / 2.0
    elseif order == 3
        result = -Spectral.kernelFermiT_dω3(τ, ϵ, β) / 6.0
    elseif order == 4
        result = Spectral.kernelFermiT_dω4(τ, ϵ, β) / 24.0
    elseif order == 5
        result = -Spectral.kernelFermiT_dω5(τ, ϵ, β) / 120.0
    else
        error("not implemented!")
        # result = Propagator.green(τ, ϵ, β) * 0.0
    end
    return result
end

function green(τ::T, ω::T, β::T) where {T}
    #generate green function of fermion
    if τ ≈ T(0.0)
        τ = -TAU_CUTOFF
    end
    if τ > T(0.0)
        return ω > T(0.0) ?
               exp(-ω * τ) / (1 + exp(-ω * β)) :
               exp(ω * (β - τ)) / (1 + exp(ω * β))
    else
        return ω > T(0.0) ?
               -exp(-ω * (τ + β)) / (1 + exp(-ω * β)) :
               -exp(-ω * τ) / (1 + exp(ω * β))
    end
end

function _leafstates(leaf_maps::Vector{Dict{Int,G}}, maxloopNum::Int) where {G<:Graph}

    num_g = length(leaf_maps)
    leafType = [Vector{Int}() for _ in 1:num_g]
    leafOrders = [Vector{Vector{Int}}() for _ in 1:num_g]
    leafInTau = [Vector{Int}() for _ in 1:num_g]
    leafOutTau = [Vector{Int}() for _ in 1:num_g]
    leafLoopIndex = [Vector{Int}() for _ in 1:num_g]
    leafValue = [Vector{Float64}() for _ in 1:num_g]
    leafval_map = [Dict{Int,Int}() for _ in 1:num_g]

    loopbasis = Vector{Float64}[]
    # tau_labels = Vector{Int}[]
    for (ikey, leafmap) in enumerate(leaf_maps)
        len_leaves = length(keys(leafmap))
        sizehint!(leafType[ikey], len_leaves)
        sizehint!(leafOrders[ikey], len_leaves)
        sizehint!(leafInTau[ikey], len_leaves)
        sizehint!(leafOutTau[ikey], len_leaves)
        sizehint!(leafLoopIndex[ikey], len_leaves)
        leafValue[ikey] = ones(Float64, len_leaves)

        valIdx = 1
        for idx in 1:len_leaves
            leaf = leafmap[idx]
            @assert IR.isleaf(leaf)
            diagId, leaf_orders = leaf.properties, leaf.orders
            loopmom = copy(diagId.extK)
            len = length(loopmom)
            @assert maxloopNum >= len

            if maxloopNum > length(loopmom)
                Base.append!(loopmom, zeros(Float64, maxloopNum - len))
            end
            flag = true
            for bi in eachindex(loopbasis)
                if loopbasis[bi] ≈ loopmom
                    push!(leafLoopIndex[ikey], bi)
                    flag = false
                    break
                end
            end
            if flag
                push!(loopbasis, loopmom)
                push!(leafLoopIndex[ikey], length(loopbasis))
            end

            # push!(tau_labels, collect(diagId.extT))
            push!(leafInTau[ikey], diagId.extT[1])
            push!(leafOutTau[ikey], diagId.extT[2])

            push!(leafOrders[ikey], leaf_orders)
            push!(leafType[ikey], FrontEnds.index(typeof(diagId)))
            leafval_map[ikey][leaf.id] = valIdx
            valIdx += 1
        end
    end

    return (leafValue, leafType, leafOrders, leafInTau, leafOutTau, leafLoopIndex), loopbasis, leafval_map
end

main()