push!(LOAD_PATH, pwd())
using Atom
using Lehmann
using MCIntegration
using Printf
using Measurements
using JLD2
using DataStructures
# using LinearAlgebra: det
using LinearAlgebra
using Random
using Dates

include("generate_freeE_NLErec.jl")
include("dof_utils.jl")

struct ParaMC
    μ::Float64
    U::Float64
    t::Float64
    β::Float64
    n::Int
    Lx::Int
    Ly::Int
    lambda::Float64
    dμ::Float64
    order::Int
    ϵk::Array{Float64,4}
end

paraid(p::ParaMC) = Dict(
    "order" => p.order,
    "beta" => p.β,
    "lambda" => p.lambda,
    "mu" => p.μ,
    "dmu" => p.dμ,
    "U" => p.U,
    "Lx" => p.Lx,
    "Ly" => p.Ly,
)
short(p::ParaMC) = join(["$(k)_$(v)" for (k, v) in sort!(OrderedDict(paraid(p)))], "_")

function neighbor(partitions; order_diff=1)
    n = Vector{Tuple{Int,Int}}()
    Nnorm = length(partitions) + 1 # the index of the normalization diagram is the N+1
    for (ip, p) in enumerate(partitions)
        # if p[1] == 1 # if there is only one loop, then the diagram can be connected to the normalization diagram
        if p[1] in [0, 1, 2] # if there is only one loop, then the diagram can be connected to the normalization diagram
            push!(n, (ip, Nnorm))
        end
        for (idx, np) in enumerate(partitions)
            if idx >= ip
                continue
            end
            if (np[1] == p[1] && (np[2] == p[2] || np[2] == p[2] + 1 || np[2] == p[2] - 1)) ||
               ((np[1] == p[1] + 2 || np[1] == p[1] - 2) && np[2] == p[2]) ||
               ((np[1] == p[1] + order_diff || np[1] == p[1] - order_diff) && np[2] == p[2]) ||
               ((np[1] == p[1] + order_diff || np[1] == p[1] - order_diff) && (np[2] == p[2] + 1 || np[2] == p[2] - 1)) ||
               ((np[1] == p[1] + order_diff || np[1] == p[1] - order_diff) && (np[2] == p[2] + 2 || np[2] == p[2] - 2))
                #the first index is the number of loops; the second index is the number of space variables
                push!(n, (ip, idx))
            end
        end
    end
    # println(n)
    return n
end

function disperion_PBC(Lx, Ly, t, dμ=0.0)
    No = 2 # spin up/down
    ϵk = zeros(Float64, (No, No, Lx, Ly)) # julia column major, the first index is the major index

    for xi in 1:Lx
        for yi in 1:Ly
            k = [2π * (xi - 1) / Lx, 2π * (yi - 1) / Ly]
            ϵk[1, 1, xi, yi] = -2t * sum(cos.(k)) + dμ
            ϵk[2, 2, xi, yi] = -2t * sum(cos.(k)) + dμ
        end
    end
    return ϵk
end

function disperion_FBC(Lx, Ly, t, dμ=0.0)
    No = 2 # spin up/down
    ϵk = zeros(Float64, (No, No, Lx, Ly)) # julia column major, the first index is the major index

    for xi in 1:Lx
        for yi in 1:Ly
            k = [π * xi / (Lx + 1), π * yi / (Ly + 1)]
            ϵk[1, 1, xi, yi] = -2t * sum(cos.(k)) + dμ
            ϵk[2, 2, xi, yi] = -2t * sum(cos.(k)) + dμ
        end
    end
    return ϵk
end

@fastmath function propagator(τ::T, ω::T, β::T) where {T}
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

function propagator_derivative(τ, ϵ, β, order)
    if order == 0
        result = propagator(τ, ϵ, β)
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
    end
    return result
end

function hopping_counterterm_PBC(para::ParaMC, τ::T, r1::Vector{Int}, r2::Vector{Int}, orbital::Int, order::Int) where {T}
    β, Lx, Ly, ϵk = para.β, para.Lx, para.Ly, para.ϵk
    g2c = 0.0
    N = Lx * Ly
    for xi in 1:Lx
        for yi in 1:Ly
            k = [2π * (xi - 1) / Lx, 2π * (yi - 1) / Ly]
            ω = -1.0 / ϵk[orbital, orbital, xi, yi]
            lambda = sign(ω) * para.lambda
            ω /= lambda
            g2c_τ = 0.0
            for o in 0:order
                g2c_τ += propagator_derivative(τ, ω, β, o) * ω^o * binomial(order, o) * (-1)^o
            end
            g2c += cos(dot(k, (r1 - r2))) * g2c_τ / lambda / N
        end
    end
    return g2c
end

function hopping_counterterm_FBC(para::ParaMC, r1::Vector{Int}, r2::Vector{Int}, orbital::Int)
    β, Lx, Ly, ϵk = para.β, para.Lx, para.Ly, para.ϵk
    g2c = 0.0
    prefactor = 4 / (Lx + 1) / (Ly + 1)
    for xi in 1:Lx
        for yi in 1:Ly
            k = [π * xi / (Lx + 1), π * yi / (Ly + 1)]
            ω = -1.0 / ϵk[orbital, orbital, xi, yi]
            lambda = sign(ω) * para.lambda
            ω /= lambda

            g2c_τ = 0.0
            for o in 0:order
                g2c_τ += propagator_derivative(τ, ω, β, o) * ω^o * binomial(order, o) * (-1)^o
            end

            phi_r1 = prod(sin.(k .* r1))
            phi_r2 = prod(sin.(k .* r2))
            g2c += phi_r1 * phi_r2 * g2c_τ / lambda
        end
    end
    return g2c * prefactor
end

function integrand(idx, vars, config)
    para, root, graphfuncs! = config.userdata[1:3]
    leafval, leafType, leafOrders, leafSites, leafτ_i, leafτ_o, leaforbitals_i, leaforbitals_o = config.userdata[4]
    model = config.userdata[5]
    varT, (varRx, varRy), varT_D = vars
    τp = varT_D[1]

    num_varR = config.dof[idx][2] + 1
    # varR = collect(zip(varRx[1:num_varR], varRy[1:num_varR]))
    # if length(Set(varR)) != length(varR)
    #     return 0.0
    # end
    # 使用双重循环检查重叠 (O(N^2) 但无内存分配，对于小阶数 N 极快)
    has_overlap = false
    @inbounds for i in 1:num_varR
        xi, yi = varRx[i], varRy[i]
        for j in (i+1):num_varR
            if xi == varRx[j] && yi == varRy[j]
                has_overlap = true
                break
            end
        end
        if has_overlap
            break
        end
    end

    if has_overlap
        return 0.0
    end

    for (i, lftype) in enumerate(leafType[idx])
        if lftype == 0
            continue
        elseif lftype == 3  # BareGreenNId
            τi, τo = varT[leafτ_i[idx][i]], varT[leafτ_o[idx][i]]
            orbitals_i, orbitals_o = leaforbitals_i[idx][i], leaforbitals_o[idx][i]
            _gn = Green.GreenN(model, vcat(τi, τo), vcat(orbitals_i, orbitals_o))
            order = leafOrders[idx][i][2]

            if order == 0
                leafval[idx][i] = Green.Gn(model, _gn)
            elseif order == 1
                leafval[idx][i] = Green.dGn_dU_estimator(model, _gn, τp)
            end
        elseif lftype == 4  # BareHoppingId
            τ = varT[leafτ_o[idx][i][1]] - varT[leafτ_i[idx][i][1]]
            r1 = collect(varR[leafSites[idx][i][1]])
            r2 = collect(varR[leafSites[idx][i][2]])

            order = leafOrders[idx][i][1]
            orbital = leaforbitals_i[idx][i][1]
            leafval[idx][i] = hopping_counterterm_PBC(para, τ, r1, r2, orbital, order)
        else
            error("this leaftype $lftype not implemented!")
        end
    end

    graphfuncs![idx](root, leafval[idx])

    return root[1]
end

function double_occupancy(model, para::ParaMC, diagram, _neighbor; neval=1e6, print=0, dtype=ComplexF64, kwargs...)
    partition, diagpara, FeynGraphs = diagram

    println("Start compiling...", now())
    funcGraphs! = Dict{Int,Function}()
    leaf_maps = Vector{Dict{Int,Graph}}()
    for (i, key) in enumerate(partition)
        funcGraphs![i], leafmap = Compilers.compile(FeynGraphs[key])
        push!(leaf_maps, leafmap)
    end
    println("Compile finished.", now())

    leafStat = FeynmanDiagram.leafstates(leaf_maps, dtype=dtype)

    root = zeros(dtype, 1)
    T_doublon = Continuous(0.0, para.β; adapt=true, alpha=3.0)
    T = Continuous(0.0, para.β; offset=1, adapt=true, alpha=3.0)
    # T = Continuous(0.0, para.β; offset=1, adapt=false)
    T.data[1] = 0.0

    R = Discrete([(1, para.Lx), (1, para.Ly)]; offset=1, adapt=true, alpha=3.0) # the fixed R is [1, 1]

    dof = build_dof(diagpara; include_probe=true)
    obs = zeros(dtype, length(diagpara))

    println("dof: ", dof)

    config = Configuration(; var=(T, R, T_doublon), dof=dof, obs=obs, type=dtype, neighbor=_neighbor,
        userdata=(para, root, funcGraphs!, leafStat, model))
    result = integrate(integrand; config=config, neval=neval, thermal_ratio=0.2, print=print, solver=:mcmc, kwargs...)

    if isnothing(result) == false
        if print >= 0
            report(result.config)
            println(report(result, pick=o -> first(o)))
            println(result)
        end
        if print >= -2
            println(result)
        end

        datadict = Dict{eltype(partition),Any}()
        for (o, key) in enumerate(partition)
            avg, std = result.mean[o], result.stdev[o]
            datadict[key] = -measurement.(avg, std)
        end
        return datadict, result
    else
        return nothing, nothing
    end
end

function double_occupancy_MC(model, para::ParaMC; neval=1e6, partition=partition(para.order), reweight_goal=nothing,
    print=0, filename::Union{String,Nothing}=nothing, dtype=ComplexF64, _neighbor=nothing)

    println("Generating diagrams...", now())
    diagram = generate_Gnderiv1(partition)
    prinln("Generate finished.", now())

    partition = diagram[1]
    println("partition: ", partition)
    if isnothing(reweight_goal)
        reweight_goal = Float64[]
        for (order, sOrder) in partition
            if sOrder == 0
                push!(reweight_goal, 4.0)
            else
                push!(reweight_goal, 1.0)
            end
        end
        push!(reweight_goal, 2.0)
    end

    if isnothing(_neighbor)
        # _neighbor = neighbor(partition, order_diff=2)
        _neighbor = neighbor(partition, order_diff=1)
    end

    Dloc = Green.thermal_expectation(model, model.D)
    println("The local double occupancy (0-th order) is: ", Dloc)

    doublon, result = double_occupancy(model, para, diagram, _neighbor; neval=neval,
        reweight_goal=reweight_goal, dtype=dtype, print=print)

    if isnothing(doublon) == false
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (doublon,)
            end
        end
        for (ip, key) in enumerate(partition)
            println("Group ", key)
            # @printf("%10s   %10s \n", "avg", "err")
            @printf("%10s   %10s   %10s   %10s \n", "real(avg)", "err", "imag(avg)", "err")
            # @printf("%10.6f ± %10.6f\n", doublon[key].val, doublon[key].err)
            r, i = real(doublon[key]), imag(doublon[key])
            @printf("%10.6f ± %10.6f    %10.6f ± %10.6f\n", r.val, r.err, i.val, i.err)
        end
    end
    return doublon, result
end