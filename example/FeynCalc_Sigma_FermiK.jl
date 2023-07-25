# This example demonstrated how to calculate the self-energy diagrams of free electrons with μ and λ counterterms
# in various orders using the FeynmanDiagram and MCIntegration module.
using FeynmanDiagram, MCIntegration, Lehmann
using LinearAlgebra, Random, Printf, Measurements
using StaticArrays, AbstractTrees
# using Profile, ProfileView
import ElectronLiquid: UEG

diag_type = :sigma
datatype = ComplexF64
# datatype = Float64
# has_counterterm = true
has_counterterm = false
Steps = 1e7
# Steps = 2e8
Base.@kwdef struct Para
    rs::Float64 = 0.5
    beta::Float64 = 50.0
    spin::Int = 2
    # Qsize::Int = 5
    dim::Int = 2
    me::Float64 = 0.5
    λ::Float64 = 4.0

    kF::Float64 = (dim == 3) ? (9π / (2spin))^(1 / 3) / rs : sqrt(4 / spin) / rs
    # extQ::Vector{SVector{3,Float64}} = [@SVector [q, 0.0, 0.0] for q in LinRange(0.0 * kF, 2.5 * kF, Qsize)]
    extQ::Vector{SVector{2,Float64}} = [@SVector [(1.0 + q) * kF, 0.0] for q in [-0.1, -0.05, 0, 0.05, 0.1]]
    ngrid::Vector{Int} = [0,] # external Matsubara frequency
    β::Float64 = beta / (kF^2 / 2me)
    μ::Float64 = kF^2 / 2me
end

function green(τ::T, ω::T, β::T) where {T}
    #generate green function of fermion
    if τ ≈ T(0.0)
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

function phase(varT, extT, l, β)
    tin, tout = varT[extT[1]], varT[extT[2]]
    return exp(1im * π * (2l + 1) / β * (tout - tin))
end

function integrand(idx, vars, config) #for the mcmc algorithm
    varK, varT, varN, Ext = vars
    para = config.userdata[1]
    MaxLoopNum = config.userdata[2]
    leaf, leafType, leafτ_i, leafτ_o, leafMomIdx, extT_labels = config.userdata[3]
    LoopPool = config.userdata[4]
    root = config.userdata[5]
    graphfuncs! = config.userdata[6]

    dim, β, me, λ, μ, ngrid = para.dim, para.β, para.me, para.λ, para.μ, para.ngrid
    loopNum = config.dof[idx][1]

    extidx = Ext[1]
    varK.data[:, 1] .= para.extQ[extidx]
    FrontEnds.update(LoopPool, varK.data[:, 1:MaxLoopNum])
    for (i, lf) in enumerate(leafType[idx])
        if lf == 0
            continue
        elseif isodd(lf) #fermionic 
            τ = varT[leafτ_o[idx][i]] - varT[leafτ_i[idx][i]]
            kq = FrontEnds.loop(LoopPool, leafMomIdx[idx][i])
            ϵ = dot(kq, kq) / (2me) - μ
            order = (lf - 1) / 2
            if order == 0
                leaf[idx][i] = green(τ, ϵ, β)
                # if τ ≈ 0.0
                #     leaf[idx][i] = Spectral.kernelFermiT(-1e-10, ϵ, β)
                # else
                #     leaf[idx][i] = Spectral.kernelFermiT(τ, ϵ, β)
                # end
            elseif order == 1
                leaf[idx][i] = -Spectral.kernelFermiT_dω(τ, ϵ, β)
            elseif order == 2
                leaf[idx][i] = Spectral.kernelFermiT_dω2(τ, ϵ, β) / 2.0
            elseif order == 3
                leaf[idx][i] = -Spectral.kernelFermiT_dω3(τ, ϵ, β) / 6.0
            elseif order == 4
                leaf[idx][i] = Spectral.kernelFermiT_dω4(τ, ϵ, β) / 24.0
            elseif order == 5
                leaf[idx][i] = -Spectral.kernelFermiT_dω5(τ, ϵ, β) / 120.0
            else
                error("not implemented!")
            end
        else
            kq = FrontEnds.loop(LoopPool, leafMomIdx[idx][i])
            # leaf[idx][i] = 8π / (dot(kq, kq) + λ)
            order = lf / 2 - 1
            invK = 1.0 / (sqrt(dot(kq, kq)) + λ)
            leaf[idx][i] = 4π * invK * (λ * invK)^order
        end
    end
    graphfuncs![idx](root, leaf[idx])  # allocations due to run-time variable `idx`

    n = ngrid[varN[1]]
    factor = 1.0 / (2π)^(dim * loopNum)

    weight = sum(root[i] * phase(varT, extT, n, β) for (i, extT) in enumerate(extT_labels[idx]))
    return weight * factor
end

function LeafInfor(FeynGraphs::Dict{Tuple{Int,Int,Int},Tuple{Vector{G},Vector{Vector{Int}}}}, FermiLabel::LabelProduct, BoseLabel::LabelProduct,
    graph_keys) where {G<:Graph}
    #read information of each leaf from the generated graph and its LabelProduct, the information include type, loop momentum, imaginary time.
    num_g = length(graph_keys)
    LeafType = [Vector{Int}() for _ in 1:num_g]
    LeafInTau = [Vector{Int}() for _ in 1:num_g]
    LeafOutTau = [Vector{Int}() for _ in 1:num_g]
    LeafLoopIndex = [Vector{Int}() for _ in 1:num_g]
    Leaf = [Vector{Float64}() for _ in 1:num_g]
    ExtT_index = [Vector{Vector{Int}}() for _ in 1:num_g]

    for (ig, key) in enumerate(graph_keys)
        ExtT_index[ig] = FeynGraphs[key][2]  # external tau variables
        for j in eachindex(ExtT_index[ig])
            for g in Leaves(FeynGraphs[key][1][j])
                if g.type == FeynmanDiagram.ComputationalGraphs.GenericDiag
                    push!(LeafType[ig], 0)
                    In = Out = 1
                    push!(Leaf[ig], 0.0)
                    push!(LeafLoopIndex[ig], 1)
                else
                    if g.type == FeynmanDiagram.ComputationalGraphs.Interaction
                        push!(LeafType[ig], 0)
                        In = Out = g.vertices[1][1].label
                        push!(LeafLoopIndex[ig], 1)
                    elseif (isfermionic(g.vertices[1]))
                        In, Out = g.vertices[2][1].label, g.vertices[1][1].label
                        if FermiLabel[In][2] in [-2, -3]
                            push!(LeafType[ig], 0)
                            push!(LeafLoopIndex[ig], 1)
                        else
                            push!(LeafType[ig], FermiLabel[In][2] * 2 + 1)
                            push!(LeafLoopIndex[ig], FrontEnds.linear_to_index(FermiLabel, In)[end]) #the label of LoopPool for each fermionic leaf
                        end
                    else
                        In, Out = g.vertices[2][1].label, g.vertices[1][1].label
                        push!(LeafType[ig], BoseLabel[In][2] * 2 + 2)
                        push!(LeafLoopIndex[ig], FrontEnds.linear_to_index(BoseLabel, In)[end]) #the label of LoopPool for each bosonic leaf
                    end
                    push!(Leaf[ig], 1.0)
                end
                push!(LeafInTau[ig], FermiLabel[In][1])
                push!(LeafOutTau[ig], FermiLabel[Out][1])
            end
        end
    end
    return Leaf, LeafType, LeafInTau, LeafOutTau, LeafLoopIndex, ExtT_index
end

function measure(vars, obs, weight, config) # for vegas and vegasmc algorithms
    N = length(config.dof)
    n = vars[3][1]  #matsubara frequency
    k = vars[4][1]  #K
    for idx in 1:N
        obs[idx][n, k] += weight
    end
end

function measure(idx, vars, obs, weight, config) # for the mcmc algorithm
    n = vars[3][1]  #matsubara frequency
    k = vars[4][1]  #K
    obs[idx][n, k] += weight
end

@inline function green(str)
    return "\u001b[32m$str\u001b[0m"
end

function run(steps, gkeys::Vector{Tuple{Int,Int,Int}}; alpha=3.0)
    para = Para()
    extQ, ngrid = para.extQ, para.ngrid
    dim, kF, β = para.dim, para.kF, para.β
    MaxLoopNum = maximum([key[1] for key in gkeys]) + 2

    # FeynGraphs, FermiLabel, BoseLabel = SigmaDiagrams(MaxOrder, has_counterterm, para.dim)
    FeynGraphs, FermiLabel, BoseLabel = SigmaDiagrams(gkeys, para.dim)
    println(green("Diagrams have been read with keys: $gkeys."))

    # funcGraph!(i) = Compilers.compile([FeynGraph.subgraphs[i],]) #Compile graph i into a julia static function. 
    # println(green("Julia static function from Graph has been compiled."))
    # gkeys = keys(FeynGraphs)
    neighbor = UEG.neighbor(gkeys)
    funcGraphs! = Dict{Int,Function}(i => Compilers.compile(FeynGraphs[key][1]) for (i, key) in enumerate(gkeys))

    LoopPool = FermiLabel.labels[3]
    LeafStat = LeafInfor(FeynGraphs, FermiLabel, BoseLabel, gkeys)
    # println(green("Leaf information has been extracted."))
    # root = zeros(Float64, 2)
    root = zeros(Float64, 24)

    # T = MCIntegration.Continuous(0.0, β; grid=collect(LinRange(0.0, β, 1000)), offset=1, alpha=alpha)
    T = Continuous(0.0, β; alpha=alpha, adapt=true, offset=2)
    T.data[1], T.data[2] = 0.0, 0.0
    K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= extQ[1]
    X = MCIntegration.Discrete(1, length(ngrid), alpha=alpha)
    Ext = MCIntegration.Discrete(1, length(extQ); adapt=false)

    dof = Vector{Int}[]
    reweight_goal = Float64[]
    for (order, verOrder, SigmaOrder) in gkeys
        # print("($order, $verOrder, $SigmaOrder) ")
        # order, verOrder, SigmaOrder = [parse(Int, c) for c in g.name]
        push!(dof, [order, order, 1, 1])
        push!(reweight_goal, 4.0^(order + verOrder - 1))
    end
    # reweight_goal = [1.0, 50.0, 200.0]
    push!(reweight_goal, 2.0) # reweight for the normalization part
    obs = [zeros(datatype, length(ngrid), length(extQ)) for _ in gkeys]

    println(green("Start computing MCMC integral:"))
    # result = integrate(integrand; measure=measure, userdata=(para, MaxLoopNum, LeafStat, LoopPool, root, funcGraphs!),
    #     var=(K, T, Ext), dof=dof, obs=obs, solver=:mcmc, reweight_goal=reweight, gamma=0.5,
    #     neval=steps, print=-1, block=2) # gets compiled
    # Profile.clear_malloc_data() # clear allocations

    result = integrate(integrand; measure=measure, userdata=(para, MaxLoopNum, LeafStat, LoopPool, root, funcGraphs!),
        # var=(K, T, X, Ext), dof=dof, obs=obs, type=datatype, solver=:mcmc, niter=10, block=16, neval=steps,
        var=(K, T, X, Ext), dof=dof, obs=obs, type=datatype, solver=:mcmc, neval=steps,
        reweight_goal=reweight_goal, neighbor=neighbor, print=0, parallel=:thread)
    # ProfileView.view()

    if isnothing(result) == false
        datadict = Dict{eltype(gkeys),Any}()
        for (o, key) in enumerate(gkeys)
            if length(dof) == 1
                avg, std = result.mean, result.stdev
            else
                avg, std = result.mean[o], result.stdev[o]
            end
            # r = measurement.(real(avg), real(std))
            # i = measurement.(imag(avg), imag(std))
            r = measurement.(real(avg), real(std)) ./ β
            i = measurement.(imag(avg), imag(std)) ./ β
            data = Complex.(r, i)
            datadict[key] = data
            println("Group ", key)
            # println(datadict[key])
            @printf("%10s  %10s   %10s   %10s   %10s \n", "q/kF", "real(avg)", "err", "imag(avg)", "err")
            for (in, n) in enumerate(ngrid)
                println("n = $n")
                for (iq, q) in enumerate(extQ)
                    @printf("%10.6f  %10.6f ± %10.6f   %10.6f ± %10.6f\n", q[1] / kF, r[in, iq].val, r[in, iq].err, i[in, iq].val, i[in, iq].err)
                end
            end
        end
        report(result)
        report(result.config)
        return datadict, result
    else
        return nothing, nothing
    end
end

function run(steps, MaxOrder::Int; alpha=3.0)
    gkeys = UEG.partition(MaxOrder)
    run(steps, gkeys, alpha=alpha)
end

# gkeys = [(1, 0, 2), (2, 0, 2), (3, 0, 1)]
# gkeys = [(1, 0, 0), (2, 0, 0), (3, 0, 0)]
# run(Steps, gkeys)

# run(Steps, 2)
run(Steps, 3)
# run(Steps, 4)

