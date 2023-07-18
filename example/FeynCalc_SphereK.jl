# This example demonstrated how to calculate the diagrams of free electrons without counterterms
# in various orders using the FeynmanDiagram and MCIntegration module.
using FeynmanDiagram, MCIntegration, Lehmann
using LinearAlgebra, Random, Printf
using StaticArrays, AbstractTrees
using Profile

Steps = 1e8
Base.@kwdef struct Para
    rs::Float64 = 1.0
    beta::Float64 = 40.0
    spin::Int = 2
    Qsize::Int = 6
    n::Int = 0 # external Matsubara frequency
    dim::Int = 3
    me::Float64 = 0.5
    λ::Float64 = 1.0

    kF::Float64 = (dim == 3) ? (9π / (2spin))^(1 / 3) / rs : sqrt(4 / spin) / rs
    extQ::Vector{SVector{3,Float64}} = [@SVector [q, 0.0, 0.0] for q in LinRange(0.0 * kF, 2.5 * kF, Qsize)]
    β::Float64 = beta / (kF^2 / 2me)
end

function green(τ::T, ω::T, β::T) where {T}
    #generate green function of fermion
    if τ == T(0.0)
        τ = -1e-10
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

# macro apply_graphfunc(graphfuncs, idx, root, leaf)
#     quote
#         $(esc(graphfuncs))[$(esc(idx))]($(esc(root)), $(esc(leaf)))
#     end
# end

# function apply_graphfunc(graphfuncs, idx, root, leaf)
#     graphfuncs[idx](root, leaf)
#     return nothing
# end

function integrand(idx, vars, config) #for the mcmc algorithm
    X, T, Ext = vars
    para = config.userdata[1]
    MaxOrder = config.userdata[2]
    @assert idx in 1:MaxOrder "$(idx) is not a valid integrand"

    leaf, leafType, leafτ_i, leafτ_o, leafMomIdx = config.userdata[3]
    LoopPool = config.userdata[4]
    MomVar, kVar, root = config.userdata[5:7]
    graphfuncs! = config.userdata[8]
    # graphfunc1!, graphfunc2! = config.userdata[6:7]

    kF, β, me, λ = para.kF, para.β, para.me, para.λ

    MomVar[:, 1] = para.extQ[Ext[1]]
    Ri = view(X[1], 1:MaxOrder)
    kVar .= Ri ./ (1.0 .- Ri)
    θ = view(X[2], 1:MaxOrder)
    ϕ = view(X[3], 1:MaxOrder)
    MomVar[1, 2:end] .= kVar .* sin.(θ) .* cos.(ϕ)
    MomVar[2, 2:end] .= kVar .* sin.(θ) .* sin.(ϕ)
    MomVar[3, 2:end] .= kVar .* cos.(θ)

    FrontEnds.update(LoopPool, MomVar)

    for (i, lf) in enumerate(leafType[idx])
        if lf == 0
            continue
        elseif (lf == 1)
            τ_l = T[leafτ_o[idx][i]] - T[leafτ_i[idx][i]]
            kq = FrontEnds.loop(LoopPool, leafMomIdx[idx][i])
            ω = (dot(kq, kq) - kF^2) / (2me)
            # leaf[i] = Spectral.kernelFermiT(τ_l, ω, β) # green function of Fermion
            leaf[idx][i] = green(τ_l, ω, β) # green function of Fermion
        else
            kq = FrontEnds.loop(LoopPool, leafMomIdx[idx][i])
            leaf[idx][i] = 8π / (dot(kq, kq) + λ)
        end
    end

    graphfuncs![idx](root, leaf[idx])
    # @apply_graphfunc(graphfuncs!, idx, root, leaf[idx])
    # apply_graphfunc(graphfuncs!, idx, root, leaf[idx])
    # if idx == 1
    #     graphfunc1!(root, leaf[idx])
    # elseif idx == 2
    #     graphfunc2!(root, leaf[idx])
    # end
    for i in 1:idx
        root[1] *= Ri[i]^2 / (1.0 - Ri[i])^4 * sin(θ[i])
    end

    root[1] *= 1.0 / (2π)^(3idx)
    return root[1]
end

function integrand(vars, config) # for vegas or vegasmc algorithm
    X, T, Ext = vars
    para = config.userdata[1]
    MaxOrder = config.userdata[2]

    leaf, leafType, leafτ_i, leafτ_o, leafMomIdx = config.userdata[3]
    LoopPool = config.userdata[4]
    MomVar, kVar, root = config.userdata[5:7]
    graphfunc! = config.userdata[8]
    # graphfunc1!, graphfunc2! = config.userdata[6:7]

    kF, β, me, λ = para.kF, para.β, para.me, para.λ

    MomVar[:, 1] = para.extQ[Ext[1]]
    Ri = view(X[1], 1:MaxOrder)
    kVar .= Ri ./ (1.0 .- Ri)
    θ = view(X[2], 1:MaxOrder)
    ϕ = view(X[3], 1:MaxOrder)
    MomVar[1, 2:end] .= kVar .* sin.(θ) .* cos.(ϕ)
    MomVar[2, 2:end] .= kVar .* sin.(θ) .* sin.(ϕ)
    MomVar[3, 2:end] .= kVar .* cos.(θ)

    FrontEnds.update(LoopPool, MomVar)

    for (i, lf) in enumerate(leafType)
        if lf == 0
            continue
        elseif (lf == 1)
            τ_l = T[leafτ_o[i]] - T[leafτ_i[i]]
            kq = FrontEnds.loop(LoopPool, leafMomIdx[i])
            ω = (dot(kq, kq) - kF^2) / (2me)
            # leaf[i] = Spectral.kernelFermiT(τ_l, ω, β) # green function of Fermion
            leaf[i] = green(τ_l, ω, β) # green function of Fermion
        else
            kq = FrontEnds.loop(LoopPool, leafMomIdx[i])
            leaf[i] = 8π / (dot(kq, kq) + λ)
        end
    end
    graphfunc!(root, leaf)

    for order in 1:MaxOrder
        factor = Ri[order]^2 / (1.0 - Ri[order])^4 * sin(θ[order])
        for i in order:MaxOrder
            root[i] *= factor
        end
        root[order] *= 1.0 / (2π)^(3order)
    end
    return root
end

function LeafInfor(FeynGraph::Graph, FermiLabel::LabelProduct, BoseLabel::LabelProduct)
    #read information of each leaf from the generated graph and its LabelProduct, the information include type, loop momentum, imaginary time. For mcmc algorithm.
    num_g = length(FeynGraph.subgraphs)
    LeafType = [Vector{Int}() for _ in 1:num_g]
    LeafInTau = [Vector{Int}() for _ in 1:num_g]
    LeafOutTau = [Vector{Int}() for _ in 1:num_g]
    LeafLoopIndex = [Vector{Int}() for _ in 1:num_g]
    Leaf = [Vector{Float64}() for _ in 1:num_g]

    for ig in 1:num_g
        for g in Leaves(FeynGraph.subgraphs[ig])
            if (g.type == FeynmanDiagram.ComputationalGraphs.Interaction)
                push!(LeafType[ig], 0)
                In = Out = g.vertices[1][1].label
                push!(LeafLoopIndex[ig], 1)
            elseif (isfermionic(g.vertices[1]))
                push!(LeafType[ig], 1)
                In, Out = g.vertices[1][1].label, g.vertices[2][1].label
                push!(LeafLoopIndex[ig], FrontEnds.linear_to_index(FermiLabel, In)[end]) #the label of LoopPool for each fermionic leaf
            else
                push!(LeafType[ig], 2)
                In, Out = g.vertices[1][1].label, g.vertices[2][1].label
                push!(LeafLoopIndex[ig], FrontEnds.linear_to_index(BoseLabel, In)[end]) #the label of LoopPool for each bosonic leaf
            end
            push!(Leaf[ig], 1.0)
            push!(LeafInTau[ig], FermiLabel[In][1])
            push!(LeafOutTau[ig], FermiLabel[Out][1])
            # push!(LeafLoopIndex[ig], FrontEnds.linear_to_index(FermiLabel, In)[end]) #the label of LoopPool for each leaf
        end
    end
    return Leaf, LeafType, LeafInTau, LeafOutTau, LeafLoopIndex
end

function LeafInfor_total(FeynGraph::Graph, FermiLabel::LabelProduct, BoseLabel::LabelProduct)
    #read information of each leaf from the generated graph and its LabelProduct, the information include type, loop momentum, imaginary time. For vegas or vegasmc algorithm.
    LeafType = Vector{Int}()
    LeafInTau = Vector{Int}()
    LeafOutTau = Vector{Int}()
    LeafLoopIndex = Vector{Int}()
    Leaf = Vector{Float64}()
    for g in Leaves(FeynGraph)
        if (g.type == FeynmanDiagram.ComputationalGraphs.Interaction)
            push!(LeafType, 0)
            In = Out = g.vertices[1][1].label
        elseif (isfermionic(g.vertices[1]))
            push!(LeafType, 1)
            In, Out = g.vertices[1][1].label, g.vertices[2][1].label
        else
            push!(LeafType, 2)
            In, Out = g.vertices[1][1].label, g.vertices[2][1].label
        end
        push!(Leaf, 1.0)
        push!(LeafInTau, FermiLabel[In][1])
        push!(LeafOutTau, FermiLabel[Out][1])
        push!(LeafLoopIndex, FrontEnds.linear_to_index(FermiLabel, In)[end])
    end
    return Leaf, LeafType, LeafInTau, LeafOutTau, LeafLoopIndex
end

function measure(vars, obs, weight, config) # for vegas and vegasmc algorithms
    N = config.userdata[2]
    Ext = vars[end]
    for i in 1:N
        obs[i][Ext[1]] += weight[i]
    end
end

function measure(idx, vars, obs, weight, config) # for the mcmc algorithm
    Ext = vars[end]
    obs[idx][Ext[1]] += weight
end

@inline function green(str)
    return "\u001b[32m$str\u001b[0m"
end

function run(steps, alg, MaxOrder::Int)
    para = Para()
    extQ, Qsize, dim = para.extQ, para.Qsize, para.dim
    kF, β = para.kF, para.β
    FeynGraph, FermiLabel, BoseLabel = PolarDiagrams(:charge, MaxOrder)
    println(green("Diagrams with the largest order $MaxOrder has been read."))
    # SinGraph, FermiLabel, BoseLabel = PolarEachOrder(:charge, MaxOrder,0,0)
    # println(green("Diagram with order $MaxOrder has been read."))
    # FeynGraph = Graph(SinGraph,factor=1.0)
    #=
    function PolarEachOrder(type::Symbol, order::Int, VerOrder::Int=0, SigmaOrder::Int=0; loopPool::Union{LoopPool,Nothing}=nothing,
        # tau_labels::Union{Nothing,Vector{Int}}=nothing, GTypes::Union{Nothing,Vector{Int}}=nothing, VTypes::Union{Nothing,Vector{Int}}=nothing)

    Generates a `Graph`: the polarization diagrams with static interactions of a given order, where the actual order of diagrams equals to `order + VerOrder + 2 * SigmaOrder`.
    =#

    # Ps = Compilers.to_julia_str([FeynGraph,], name="eval_graph!")
    # Pexpr = Meta.parse(Ps)
    # println(Pexpr)
    # eval(Pexpr)
    # funcGraph!(x, y) = Base.invokelatest(eval_graph!, x, y)

    if alg == :mcmc
        # funcGraphs! = [Compilers.compile([FeynGraph.subgraphs[i],]) for i in 1:MaxOrder] #Compile graphs into a julia static function Vector. 
        funcGraphs! = Dict{Int,Function}(i => Compilers.compile([FeynGraph.subgraphs[i],]) for i in 1:MaxOrder)
        LeafStat = LeafInfor(FeynGraph, FermiLabel, BoseLabel)
    else
        funcGraphs! = Compilers.compile(FeynGraph.subgraphs)
        LeafStat = LeafInfor_total(FeynGraph, FermiLabel, BoseLabel)
    end
    # println(green("Julia static function from Graph has been compiled."))

    LoopPool = FermiLabel.labels[3]
    # println(green("Leaf information has been extracted."))
    root = zeros(Float64, MaxOrder)

    T = Continuous(0.0, β; alpha=3.0, adapt=true, offset=1)
    T.data[1] = 0.0
    R = Continuous(0.0, 1.0; alpha=3.0, adapt=true)
    θ = Continuous(0.0, 1π; alpha=3.0, adapt=true)
    ϕ = Continuous(0.0, 2π; alpha=3.0, adapt=true)
    X = CompositeVar(R, θ, ϕ)
    Ext = Discrete(1, length(extQ); adapt=false)

    MomVar = Matrix{Float64}(undef, dim, MaxOrder + 1)
    kVar = Vector{Float64}(undef, MaxOrder)

    dof = [[Order, Order, 1] for Order in 1:MaxOrder] # degrees of freedom of the diagram
    obs = [zeros(Float64, Qsize) for _ in 1:MaxOrder]

    reweight_goal = [4.0^(i - 1) for i in 1:MaxOrder]
    # reweight_goal = [2.0^(i - 1) for i in 1:MaxOrder]
    push!(reweight_goal, 1.0)

    println(green("Start computing integral with $alg:"))
    # result = integrate(integrand; measure=measure, userdata=(para, MaxOrder, LeafStat, LoopPool, MomVar, kVar, root, funcGraphs!),
    #     var=(X, T, Ext), dof=dof, obs=obs, solver=alg, reweight_goal=reweight_goal,
    #     neval=steps, print=-1, block=2) # gets compiled
    # Profile.clear_malloc_data() # clear allocations
    result = integrate(integrand; measure=measure, userdata=(para, MaxOrder, LeafStat, LoopPool, MomVar, kVar, root, funcGraphs!),
        var=(X, T, Ext), dof=dof, obs=obs, solver=alg, reweight_goal=reweight_goal,
        neval=steps, print=0, block=16)

    if isnothing(result) == false
        avg, std = result.mean, result.stdev

        for o in 1:MaxOrder
            println("Order:$o")
            @printf("%10s  %10s   %10s \n", "q/kF", "avg", "err")
            for (idx, q) in enumerate(extQ)
                q = q[1]
                if (MaxOrder == 1)
                    @printf("%10.6f  %10.6f ± %10.6f\n", q / kF, avg[idx], std[idx])
                else
                    @printf("%10.6f  %10.6f ± %10.6f\n", q / kF, avg[o][idx], std[o][idx])
                end
            end
        end
        report(result)
        report(result.config)
    end
end

# run(Steps, :mcmc, 1)
run(Steps, :mcmc, 3)
# run(Steps, :vegasmc, 3)

