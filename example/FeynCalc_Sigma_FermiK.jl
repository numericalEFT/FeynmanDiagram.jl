# This example demonstrated how to calculate the self-energy diagrams of free electrons with μ and λ counterterms
# in various orders using the FeynmanDiagram and MCIntegration module.
using FeynmanDiagram, MCIntegration, Lehmann
using LinearAlgebra, Random, Printf, Measurements
using StaticArrays, AbstractTrees
using Profile, ProfileView
import ElectronLiquid: UEG

diag_type = :sigma
datatype = ComplexF64
# datatype = Float64
has_counterterm = true
Steps = 1e8
Base.@kwdef struct Para
    rs::Float64 = 0.5
    beta::Float64 = 50.0
    spin::Int = 2
    Qsize::Int = 5
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

# function _apply_graphfunc(graphfuncs, ::Val{I}, root, leaf) where {I}
#     graphfuncs[I](root, leaf)
#     return nothing
# end

# function _apply_graphfunc(graphfuncs, I, root, leaf)
#     graphfuncs[I](root, leaf)
#     return nothing
# end

function integrand(idx, vars, config) #for the mcmc algorithm
    K, T, N, Ext = vars
    para = config.userdata[1]
    MaxOrder = config.userdata[2]
    # @assert idx in 1:MaxOrder "$(idx) is not a valid integrand"
    leaf, leafType, leafτ_i, leafτ_o, leafMomIdx = config.userdata[3]
    LoopPool = config.userdata[4]
    root = config.userdata[5]
    graphfuncs! = config.userdata[6]
    # graphfunc1!, graphfunc2! = config.userdata[6:7]

    dim, β, me, λ, μ, ngrid = para.dim, para.β, para.me, para.λ, para.μ, para.ngrid
    loopNum = config.dof[idx][1]

    extidx = Ext[1]
    K.data[:, 1] .= para.extQ[extidx]
    # FrontEnds.update(LoopPool, K.data[:, 1:MaxOrder+1])
    FrontEnds.update(LoopPool, K.data[:, 1:MaxOrder+2])
    for (i, lf) in enumerate(leafType[idx])
        if lf == 0
            continue
        elseif isodd(lf)
            τ = T[leafτ_o[idx][i]] - T[leafτ_i[idx][i]]
            kq = FrontEnds.loop(LoopPool, leafMomIdx[idx][i])
            ϵ = dot(kq, kq) / (2me) - μ
            order = (lf - 1) / 2
            if order == 0
                if τ ≈ 0.0
                    leaf[idx][i] = Spectral.kernelFermiT(-1e-10, ϵ, β)
                else
                    leaf[idx][i] = Spectral.kernelFermiT(τ, ϵ, β)
                end
            elseif order == 1
                leaf[idx][i] = -Spectral.kernelFermiT_dω(τ, ϵ, β)
            elseif order == 2
                leaf[idx][i] = Spectral.kernelFermiT_dω2(τ, ϵ, β)
            elseif order == 3
                leaf[idx][i] = -Spectral.kernelFermiT_dω3(τ, ϵ, β)
            elseif order == 4
                leaf[idx][i] = Spectral.kernelFermiT_dω4(τ, ϵ, β)
            elseif order == 5
                leaf[idx][i] = -Spectral.kernelFermiT_dω5(τ, ϵ, β)
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
    # @apply_graphfunc(graphfuncs!, idx, root, leaf[idx])
    # _apply_graphfunc(graphfuncs!, Val(idx), root, leaf[idx])
    # _apply_graphfunc(graphfuncs!, idx, root, leaf[idx])
    # if idx == 1
    #     graphfunc1!(root, leaf[idx])
    # elseif idx == 2
    #     graphfunc2!(root, leaf[idx])
    # end
    n = ngrid[N[1]]

    phase1 = 1.0 / (2π)^(dim * loopNum)
    phase2 = phase1 * exp(1im * π * (2n + 1) / β * T[2])

    # root[1] *= phase
    return root[1] * phase1 + root[2] * phase2
end

# function LeafInfor(FeynGraph::Graph, FermiLabel::LabelProduct, BoseLabel::LabelProduct)
function LeafInfor(FeynGraphs::Dict{Tuple{Int,Int,Int},Vector{G}}, FermiLabel::LabelProduct, BoseLabel::LabelProduct,
    graph_keys) where {G<:Graph}
    #read information of each leaf from the generated graph and its LabelProduct, the information include type, loop momentum, imaginary time.
    num_g = length(graph_keys)
    LeafType = [Vector{Int}() for _ in 1:num_g]
    LeafInTau = [Vector{Int}() for _ in 1:num_g]
    LeafOutTau = [Vector{Int}() for _ in 1:num_g]
    LeafLoopIndex = [Vector{Int}() for _ in 1:num_g]
    Leaf = [Vector{Float64}() for _ in 1:num_g]

    for (ig, key) in enumerate(graph_keys)
        for j in 1:2
            for g in Leaves(FeynGraphs[key][j])
                if g.type == FeynmanDiagram.ComputationalGraphs.GenericDiag
                    push!(LeafType[ig], 0)
                    In = Out = 1
                    push!(Leaf[ig], 0.0)
                else
                    if g.type == FeynmanDiagram.ComputationalGraphs.Interaction
                        push!(LeafType[ig], 0)
                        In = Out = g.vertices[1][1].label
                    elseif (isfermionic(g.vertices[1]))
                        In, Out = g.vertices[1][1].label, g.vertices[2][1].label
                        if FermiLabel[In][2] in [-2, -3]
                            push!(LeafType[ig], 0)
                        else
                            push!(LeafType[ig], FermiLabel[In][2] * 2 + 1)
                        end
                    else
                        In, Out = g.vertices[1][1].label, g.vertices[2][1].label
                        push!(LeafType[ig], BoseLabel[In][2] * 2 + 2)
                    end
                    push!(Leaf[ig], 1.0)
                end
                push!(LeafInTau[ig], FermiLabel[In][1])
                push!(LeafOutTau[ig], FermiLabel[Out][1])
                push!(LeafLoopIndex[ig], FrontEnds.linear_to_index(FermiLabel, In)[end]) #the label of LoopPool for each leaf
            end
        end
        println(key)
        println(LeafType[ig])
        println(Leaf[ig])
    end
    return Leaf, LeafType, LeafInTau, LeafOutTau, LeafLoopIndex
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

function run(steps, MaxOrder::Int; alpha=3.0)
    para = Para()
    extQ, ngrid = para.extQ, para.ngrid
    dim, kF, β = para.dim, para.kF, para.β
    FeynGraphs, FermiLabel, BoseLabel = SigmaDiagrams(diag_type, MaxOrder, has_counterterm, para.dim)
    println(green("Diagrams with the largest order $MaxOrder has been read."))

    # Ps = Compilers.to_julia_str([FeynGraph,], name="eval_graph!")
    # Pexpr = Meta.parse(Ps)
    # println(Pexpr)
    # eval(Pexpr)
    # funcGraph!(x, y) = Base.invokelatest(eval_graph!, x, y)

    # funcGraphs! = Tuple([Compilers.compile([FeynGraph.subgraphs[i],]) for i in 1:MaxOrder]) #Compile graphs into a julia static function Vector. 
    # funcGraph!(i) = Compilers.compile([FeynGraph.subgraphs[i],]) #Compile graph i into a julia static function. 
    # println(green("Julia static function from Graph has been compiled."))
    # inds = eachindex(FeynGraph.subgraphs)
    # funcGraphs! = Dict{Int,Function}(i => Compilers.compile([FeynGraph.subgraphs[i],]) for i in inds)
    # gkeys = keys(FeynGraphs)
    gkeys = UEG.partition(MaxOrder)
    neighbor = UEG.neighbor(gkeys)
    funcGraphs! = Dict{Int,Function}(i => Compilers.compile(FeynGraphs[key]) for (i, key) in enumerate(gkeys))

    LoopPool = FermiLabel.labels[3]
    LeafStat = LeafInfor(FeynGraphs, FermiLabel, BoseLabel, gkeys)
    # println(green("Leaf information has been extracted."))
    root = zeros(Float64, 2)

    T = Continuous(0.0, β; alpha=alpha, adapt=true, offset=1)
    # T = Continuous(0.0, β; adapt=false, offset=1)
    T.data[1] = 0.0
    K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= extQ[1]
    X = MCIntegration.Discrete(1, length(ngrid), alpha=alpha)
    Ext = MCIntegration.Discrete(1, length(extQ); adapt=false)

    dof = Vector{Int}[]
    reweight_goal = Float64[]
    for (order, verOrder, SigmaOrder) in gkeys
        println("($order, $verOrder, $SigmaOrder)")
        # order, verOrder, SigmaOrder = [parse(Int, c) for c in g.name]
        push!(dof, [order, order - 1, 1, 1])
        push!(reweight_goal, 4.0^(order + verOrder - 1))
    end
    push!(reweight_goal, 2.0) # reweight for the normalization part
    obs = [zeros(datatype, length(ngrid), length(extQ)) for _ in gkeys]

    println(green("Start computing MCMC integral:"))
    # result = integrate(integrand; measure=measure, userdata=(para, MaxOrder, LeafStat, LoopPool, root, funcGraphs!),
    #     var=(K, T, Ext), dof=dof, obs=obs, solver=:mcmc, reweight_goal=reweight, gamma=0.5,
    #     neval=steps, print=-1, block=2) # gets compiled
    # Profile.clear_malloc_data() # clear allocations

    result = integrate(integrand; measure=measure, userdata=(para, MaxOrder, LeafStat, LoopPool, root, funcGraphs!),
        var=(K, T, X, Ext), dof=dof, obs=obs, type=datatype, solver=:mcmc, niter=10, block=16, neval=steps,
        reweight_goal=reweight_goal, neighbor=neighbor, print=0, parallel=:thread)
    # result = integrate(integrand; measure=measure, userdata=(para, MaxOrder, LeafStat, LoopPool, root, funcGraphs!),
    #     var=(K, T, Ext), dof=dof, obs=obs, solver=:mcmc, niter=10, reweight=reweight_goal, gamma=0,
    #     neval=steps, print=0, block=16, parallel=:thread)   # strictly setting reweight as reweight_goal (with gamma=0)
    # ProfileView.view()

    if isnothing(result) == false
        datadict = Dict{eltype(gkeys),Any}()
        for (o, key) in enumerate(gkeys)
            if length(dof) == 1
                avg, std = result.mean, result.stdev
            else
                avg, std = result.mean[o], result.stdev[o]
            end
            r = measurement.(real(avg), real(std))
            i = measurement.(imag(avg), imag(std))
            data = Complex.(r, i)
            datadict[key] = data
            println("Group ", key)
            println(datadict[key])
        end

        # for key in gkeys
        #     println("Group ", key)
        #     println(datadict[key])
        #     # @printf("%10s  %10s   %10s   %10s   %10s \n", "q/kF", "real(avg)", "err", "imag(avg)", "err")
        #     # for (in, n) in enumerate(ngrid)
        #     #     for (idx, q) in enumerate(extQ)
        #     #         q = q[1]
        #     #         println(avg)
        #     #         println(std)
        #     #         if (MaxOrder == 1)
        #     #             @printf("%10.6f  %10.6f ± %10.6f  %10.6f ± %10.6f\n", q / kF, real(avg[idx]), real(std[idx]), imag(avg[idx]), imag(std[idx]))
        #     #             @printf("%10.6f  %10.6f ± %10.6f  %10.6f ± %10.6f\n", q / kF, real(avg[idx]), real(std[idx]), imag(avg[idx]), imag(std[idx]))
        #     #             # println(q / kF, avg[idx], std[idx])
        #     #         else
        #     #             @printf("%10.6f  %10.6f ± %10.6f  %10.6f ± %10.6f\n", q / kF, real(avg[o][idx]), real(std[o][idx]), imag(avg[idx]), imag(std[idx]))
        #     #             @printf("%10.6f  %10.6f ± %10.6f  %10.6f ± %10.6f\n", q / kF, real(avg[o][idx]), real(std[o][idx]), imag(avg[idx]), imag(std[idx]))
        #     #             # println(q / kF, avg[o][idx], std[idx])
        #     #         end
        #     #     end
        #     # end
        # end
        report(result)
        report(result.config)
        return datadict, result
    else
        return nothing, nothing
    end
end

# run(Steps, 1)
run(Steps, 2)
# run(Steps, 3)
# run(Steps, 4)

