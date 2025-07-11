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

include("free_energy_v1.jl")

struct ParaMC
    μ::Float64
    U::Float64
    β::Float64
    n::Int
    Lx::Int
    Ly::Int
    lambda::Float64
    order::Int
    ϵk::Array{Float64,4}
end

paraid(p::ParaMC) = Dict(
    "order" => p.order,
    "beta" => p.β,
    "lambda" => p.lambda,
    "mu" => p.μ,
    "U" => p.U,
    "Lx" => p.Lx,
    "Ly" => p.Ly,
)
short(p::ParaMC) = join(["$(k)_$(v)" for (k, v) in sort!(OrderedDict(paraid(p)))], "_")

function F0(para::ParaMC)
    return -log((1 + exp(para.β)) * (1 + exp(-para.β))) / para.β
end
# println("F0 = ", F0(para))

function disperion_PBC(Lx, Ly, t)
    No = 2 # spin up/down
    ϵk = zeros(Float64, (No, No, Lx, Ly)) # julia column major, the first index is the major index

    for xi in 1:Lx
        for yi in 1:Ly
            kx, ky = 2π * (xi - 1) / Lx, 2π * (yi - 1) / Ly
            ϵk[1, 1, xi, yi] = -2t * (cos(kx) + cos(ky))
            ϵk[2, 2, xi, yi] = -2t * (cos(kx) + cos(ky))
        end
    end
    return ϵk
end

function disperion_FBC(Lx, Ly, t)
    No = 2 # spin up/down
    ϵk = zeros(Float64, (No, No, Lx, Ly)) # julia column major, the first index is the major index

    for xi in 1:Lx
        for yi in 1:Ly
            k = [π * xi / (Lx + 1), π * yi / (Ly + 1)]
            ϵk[1, 1, xi, yi] = -2t * sum(cos.(k))
            ϵk[2, 2, xi, yi] = -2t * sum(cos.(k))
        end
    end
    return ϵk
end

function propagator(τ::T, ω::T, β::T) where {T}
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

# function green_counterterm_PBC(para::ParaMC, τ::T, rx::Int, ry::Int, orbital::Int, order::Int) where {T}
function green_counterterm_PBC(para::ParaMC, τ::T, r1::Vector{Int}, r2::Vector{Int}, orbital::Int, order::Int) where {T}
    β, Lx, Ly, ϵk = para.β, para.Lx, para.Ly, para.ϵk
    g2c = 0.0
    N = Lx * Ly
    for xi in 1:Lx
        for yi in 1:Ly
            # kx, ky = 2π * (xi - 1) / Lx, 2π * (yi - 1) / Ly
            k = [2π * (xi - 1) / Lx, 2π * (yi - 1) / Ly]
            ω = -1.0 / ϵk[orbital, orbital, xi, yi]
            lambda = sign(ω) * para.lambda
            ω /= lambda
            g2c_τ = 0.0
            for o in 0:order
                g2c_τ += propagator_derivative(τ, ω, β, o) * ω^o * binomial(order, o) * (-1)^o
            end

            # g2c += cos(kx * rx + ky * ry) * g2c_τ / lambda / N
            # g2c += cos(dot(k, (r1 + r2))) * g2c_τ / lambda / N
            # println(exp(im * dot(k, (r1 - r2))))
            # g2c += cos(dot(k, (r1 - r2))) * g2c_τ / lambda / N
            g2c += exp(im * dot(k, (r1 - r2))) * g2c_τ / lambda / N
        end
    end
    return g2c
end

function green_counterterm_FBC(para::ParaMC, τ::T, r1::Vector{Int}, r2::Vector{Int}, orbital::Int, order::Int) where {T}
    β, Lx, Ly, ϵk = para.β, para.Lx, para.Ly, para.ϵk
    g2c = 0.0
    prefactor = 4 / (Lx + 1) / (Ly + 1)
    for xi in 1:Lx
        for yi in 1:Ly
            k = [π * xi / (Lx + 1), π * yi / (Ly + 1)]
            # ω = -1.0 / ϵk[orbital, orbital, xi, yi] / para.lambda
            ω = -1.0 / ϵk[orbital, orbital, xi, yi]
            lambda = sign(ω) * para.lambda
            ω /= lambda

            # println("ω: ", ω)
            g2c_τ = 0.0
            for o in 0:order
                g2c_τ += propagator_derivative(τ, ω, β, o) * ω^o * binomial(order, o) * (-1)^o
                # g2c_τ += propagator_derivative(τ, ω, β, o) * ω^o * binomial(order, o)
            end

            phi_r1 = prod(sin.(k .* r1))
            phi_r2 = prod(sin.(k .* r2))
            # g2c += phi_r1 * phi_r2 * g2c_τ / para.lambda
            g2c += phi_r1 * phi_r2 * g2c_τ / lambda
            # println("phi: ", phi_r1 * phi_r2)
        end
    end
    return g2c * prefactor
end

function g2(para, τ)
    β, μ, U = para.β, para.μ, para.U
    Z = 1 + 2 * exp(β * μ) + exp(β * (2μ - U))

    if τ > 0
        return (exp(μ * τ) + exp(β * μ) * exp((μ - U) * τ)) / Z
    else
        return -(exp(β * μ) * exp(μ * τ) + exp(β * (2μ - U)) * exp((μ - U) * τ)) / Z
    end
end

@inline function find_duplicates_with_indices(ri::Vector{T}, ro::Vector{T}) where {T}
    element_indices = Dict{eltype(ri),Tuple{Vector{Int},Vector{Int}}}()
    # element_indices = Dict{eltype(ri),Vector{Int}}()

    sites = Set(vcat(ri, ro))
    for elem in sites
        idx_i = findall(x -> x == elem, ri)
        idx_o = findall(x -> x == elem, ro)
        if length(idx_i) != length(idx_o)
            return nothing
        end
        element_indices[elem] = (idx_i, idx_o)
        # element_indices[elem] = collect(Iterators.flatten(zip(idx_i, idx_o)))
    end
    # for (idx, elem) in pairs(v)
    #     indices = get!(() -> Int[], element_indices, elem)
    #     push!(indices, idx)
    # end
    return element_indices
end

@inline function wick_sign(v::Vector{Int})
    sign = 1
    for i in eachindex(v)
        for j in (i+1):length(v)
            if (v[i] > v[j]) || (v[i] == v[j] && iseven(i) && isodd(j))
                sign *= -1
            end
        end
    end
    return sign
end

function integrand(idx, vars, config)
    para, root, graphfuncs! = config.userdata[1:3]
    leafval, leafType, leafOrders, leafSites, leafτ_i, leafτ_o, leaforbitals_i, leaforbitals_o = config.userdata[4]
    model = config.userdata[5]
    varT, varR = vars

    varR_all = varR
    # varR_all = vcat(varR[1:numR], varR[1:numR])
    # shuffle!(config.rng, varR_all)

    # println("varT: ", varT[1], " ", varT[2])
    # println("idx: ", idx, " ", varR[1:4idx])

    for (i, lftype) in enumerate(leafType[idx])
        if lftype == 0
            continue
        elseif lftype == 4  # BareHoppingId
            # println("hopping t: ", leafτ_o[idx][i][1], " ", leafτ_i[idx][i][1])
            τ = varT[leafτ_o[idx][i][1]] - varT[leafτ_i[idx][i][1]]
            r1 = [varR_all[leafSites[idx][i][1]], 1]
            r2 = [varR_all[leafSites[idx][i][2]], 1]

            order = leafOrders[idx][i][1]
            orbital = leaforbitals_i[idx][i][1]
            # leafval[idx][i] = green_counterterm_FBC(para, τ, r1, r2, orbital, order)
            leafval[idx][i] = green_counterterm_PBC(para, τ, r1, r2, orbital, order)
        elseif lftype == 5  # BareGreenNId
            # println("bare green: ", leafτ_o[idx][i], " ", leafτ_i[idx][i])
            # τ = vcat(varT[leafτ_i[idx][i]], varT[leafτ_o[idx][i]])
            # orbitals = vcat(leaforbitals_i[idx][i], leaforbitals_o[idx][i])

            τi, τo = varT[leafτ_i[idx][i]], varT[leafτ_o[idx][i]]
            orbitals_i, orbitals_o = leaforbitals_i[idx][i], leaforbitals_o[idx][i]

            Np = Int(length(leafSites[idx][i]) / 2)
            sites_i = varR_all[leafSites[idx][i][1:Np]]
            sites_o = varR_all[leafSites[idx][i][Np+1:end]]
            r_dict = find_duplicates_with_indices(sites_i, sites_o)

            # println("r_dict: ", r_dict)
            if isnothing(r_dict)
                leafval[idx][i] = 0.0
                continue
            end

            # r_locs = values(r_dict)
            # if any(isodd(length(locs)) for locs in r_locs)
            #     leafval[idx][i] = 0.0
            #     continue
            # end
            leafval[idx][i] = wick_sign(collect(Iterators.flatten(zip(sites_i, sites_o))))
            # leafval[idx][i] = 1.0
            # println("r:, $(leafSites[idx][i])", "τ: $τ", r_dict)

            # println(collect(Iterators.flatten(zip(sites_i, sites_o))), leafval[idx][i])

            for (loc_i, loc_o) in values(r_dict)
                # println(idx, " loc_i: ", loc_i, " loc_o: ", loc_o)
                # loc = vcat(loc_i, loc_o)
                # _gn = Green.GreenN(model, τ[loc], orbitals[loc])
                τ = vcat(τi[loc_i], τo[loc_o])
                orbitals = vcat(orbitals_i[loc_i], orbitals_o[loc_o])
                _gn = Green.GreenN(model, τ, orbitals)
                leafval[idx][i] *= Green.Gn(model, _gn)
            end
        else
            error("this leaftype $lftype not implemented!")
        end
    end

    graphfuncs![idx](root, leafval[idx])

    # if idx >= 2
    # 	# if idx == 4
    # 	# println(idx, " ", root[1])
    # 	println(idx, " ", varR.data[1:4], " ", root[1])
    # end
    # idx == 2 && println("idx: ", idx, " ", varR_all[1:4idx], " ", floor.(Int, config.propose[2, :, :]), floor.(Int, config.accept[2, :, :]), " ", root[1])
    return root[1]
end

function freeE(model, para::ParaMC, diagram; neval=1e6, print=0, dtype=ComplexF64, kwargs...)
    partition, diagpara, FeynGraphs = diagram

    funcGraphs! = Dict{Int,Function}()
    leaf_maps = Vector{Dict{Int,Graph}}()
    for (i, key) in enumerate(partition)
        funcGraphs![i], leafmap = Compilers.compile(FeynGraphs[key])
        push!(leaf_maps, leafmap)
        # println(funcGraphs![i])
    end

    leafStat = FeynmanDiagram.leafstates(leaf_maps, dtype=dtype)

    root = zeros(dtype, 1)
    T = Continuous(0.0, para.β; offset=1, adapt=true)
    T.data[1] = 0.0
    # R = Discrete(1, para.Lx, adapt=false)
    R = Discrete(1, para.Lx, adapt=true)

    dof = [[p.totalTauNum - 1, p.innerLoopNum * 2] for p in diagpara]
    # dof = [[p.totalTauNum - 1, p.innerLoopNum] for p in diagpara]
    obs = zeros(dtype, length(diagpara))
    global_updates = [false, true]

    println("dof: ", dof)

    config = Configuration(; var=(T, R), dof=dof, obs=obs, type=dtype, global_updates=global_updates,
        userdata=(para, root, funcGraphs!, leafStat, model))
    result = integrate(integrand; config=config, neval=neval, print=print, solver=:mcmc, kwargs...)

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

function freeE_MC(model, para::ParaMC; neval=1e6, partition=partition(para.order), reweight_goal=nothing,
    print=0, filename::Union{String,Nothing}=nothing, dtype=ComplexF64, _neighbor=nothing)
    # partition, diagpara, FeynGraphs = free_energy(partition)
    diagram = free_energy(partition)

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
        _neighbor = neighbor(partition)
    end

    freeEnergy, result = freeE(model, para, diagram; neval=neval, neighbor=_neighbor,
        # freeEnergy, result = freeE(model, para, diagram; neval=neval,
        reweight_goal=reweight_goal, dtype=dtype, print=print)

    if isnothing(freeEnergy) == false
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (freeEnergy,)
            end
        end
        for (ip, key) in enumerate(partition)
            println("Group ", key)
            # @printf("%10s   %10s \n", "avg", "err")
            @printf("%10s   %10s   %10s   %10s \n", "real(avg)", "err", "imag(avg)", "err")
            # @printf("%10.6f ± %10.6f\n", freeEnergy[key].val, freeEnergy[key].err)
            r, i = real(freeEnergy[key]), imag(freeEnergy[key])
            @printf("%10.6f ± %10.6f    %10.6f ± %10.6f\n", r.val, r.err, i.val, i.err)
        end
    end
    return freeEnergy, result
end