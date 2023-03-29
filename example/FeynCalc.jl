# This example demonstrated how to calculate the diagrams of free electrons without counterterms
# in various orders using the FeynmanDiagram and MCIntegration module.
using FeynmanDiagram, MCIntegration, Lehmann
using LinearAlgebra, Random, Printf
using StaticArrays, AbstractTrees

Steps = 1e7

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

function integrand(vars, config)
    #generate the MC integrand function
    R, Theta, Phi, T, Ext = vars
    para = config.userdata[1]
    Order = config.userdata[2]
    leaf, leafType, leafτ_i, leafτ_o, leafMom = config.userdata[3]
    graphfunc! = config.userdata[4]

    kF, β, me, λ = para.kF, para.β, para.me, para.λ

    Ri = [R[i] for i in 1:Order]
    r = [R[i] / (1 - R[i]) for i in 1:Order]
    θ = [Theta[i] for i in 1:Order]
    ϕ = [Phi[i] for i in 1:Order]

    factor = 1.0 / (2π)^(para.dim)  #each momentum loop is ∫dkxdkydkz/(2π)^3
    factor *= r .^ 2 ./ (1 .- Ri) .^ 2 .* sin.(θ)

    τ = [(i == 0 ? 0.0 : T[i]) for i = 0:Order]
    extidx = Ext[1]
    q = para.extQ[extidx]
    k = [(i == 0 ? q : [r[i] * sin(θ[i]) * cos(ϕ[i]), r[i] * sin(θ[i]) * sin(ϕ[i]), r[i] * cos(θ[i])]) for i in 0:Order]
    root = [0.0,]
    for (i, lf) in enumerate(leafType)
        if (lf == 0)
            continue
        elseif (lf == 1)
            τ_l = τ[leafτ_o[i]] - τ[leafτ_i[i]]
            kq = sum(leafMom[i] .* k)
            ω = (dot(kq, kq) - kF^2) / (2me)
            # leaf[i] = Spectral.kernelFermiT(τ_l, ω, β) # green function of Fermion
            leaf[i] = green(τ_l, ω, β) # green function of Fermion
        else
            kq = sum(leafMom[i] .* k)
            leaf[i] = (8 * π) / (dot(kq, kq) + λ)
        end
    end
    return graphfunc!(root, leaf) * prod(factor)
end

function LeafInfor(FeynGraph::Graph, FermiLabel::LabelProduct, BoseLabel::LabelProduct)
    #read information of each leaf from the generated graph and its LabelProduct, the information include type, loop momentum, imaginary time.
    LeafType = Vector{Int}(undef, 0)
    LeafInTau = Vector{Int}(undef, 0)
    LeafOutTau = Vector{Int}(undef, 0)
    LeafLoopMom = Vector{Vector{Float64}}(undef, 0)
    Leaf = Vector{Float64}(undef, 0)
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
        push!(LeafLoopMom, FermiLabel[In][3])
    end
    return Leaf, LeafType, LeafInTau, LeafOutTau, LeafLoopMom
end

function integrand(idx, vars, config) #for the mcmc algorithm
    return integrand(vars, config)::Float64
end

function measure(vars, obs, weight, config) # for vegas and vegasmc algorithms
    Ext = vars[end]
    obs[1][Ext[1]] += weight[1]
end

function measure(idx, vars, obs, weight, config) # for the mcmc algorithm
    measure(vars, obs, weight, config)
end

@inline function green(str)
    return "\u001b[32m$str\u001b[0m"
end

function run(steps, Order::Int)
    para = Para()
    extQ, Qsize = para.extQ, para.Qsize
    kF, β = para.kF, para.β
    LoopNum = Order + 1
    FeynGraph, FermiLabel, BoseLabel = PolarEachOrder(:charge, Order, 0, 0)
    println(green("Diagram with order $Order has been read."))
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

    funcGraph! = Compilers.compile([FeynGraph,]) #Compile graphs into a julia static function. 
    println(green("Julia static function from Graph has been compiled."))

    LeafStat = LeafInfor(FeynGraph, FermiLabel, BoseLabel)
    println(green("Leaf information has been extracted."))

    T = Continuous(0.0, β; alpha=3.0, adapt=true)
    R = Continuous(0.0, 1.0; alpha=3.0, adapt=true)
    θ = Continuous(0.0, 1π; alpha=3.0, adapt=true)
    ϕ = Continuous(0.0, 2π; alpha=3.0, adapt=true)
    Ext = Discrete(1, length(extQ); adapt=false)

    dof = [[Order, Order, Order, Order, 1],] # degrees of freedom of the diagram
    obs = [zeros(Float64, Qsize),]

    println(green("Start computing integral:"))
    result = integrate(integrand; measure=measure, userdata=(para, Order, LeafStat, funcGraph!),
        var=(R, θ, ϕ, T, Ext), dof=dof, obs=obs, solver=:vegasmc,
        neval=steps, print=0, block=32, debug=true)

    if isnothing(result) == false
        avg, std = result.mean, result.stdev

        @printf("%10s  %10s   %10s \n", "q/kF", "avg", "err")
        for (idx, q) in enumerate(extQ)
            q = q[1]
            @printf("%10.6f  %10.6f ± %10.6f\n", q / kF, avg[idx], std[idx])
        end
        report(result)
    end
end

# run(Steps, 1)
# run(Steps, 2)
run(Steps, 3)

