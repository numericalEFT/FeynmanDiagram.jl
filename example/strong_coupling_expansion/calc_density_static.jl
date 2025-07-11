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

include("density_static.jl")

struct ParaMC
    μ::Float64
    U::Float64
    t::Float64
    β::Float64
    n::Int
    Lx::Int
    Ly::Int
    order::Int
end

paraid(p::ParaMC) = Dict(
    "order" => p.order,
    "beta" => p.β,
    "mu" => p.μ,
    "U" => p.U,
    "Lx" => p.Lx,
    "Ly" => p.Ly,
)
short(p::ParaMC) = join(["$(k)_$(v)" for (k, v) in sort!(OrderedDict(paraid(p)))], "_")

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

function green_counterterm_PBC(para::ParaMC, r1::Vector{Int}, r2::Vector{Int}, orbital::Int)
    L = [para.Lx, para.Ly]
    delta12 = abs.(r1 - r2)
    delta = min.(delta12, L .- delta12)

    if sum(delta) == 1
        return para.t
    else
        return 0.0
    end
end

function green_counterterm_FBC(para::ParaMC, r1::Vector{Int}, r2::Vector{Int}, orbital::Int)
    delta = abs.(r1 - r2)

    if sum(delta) == 1
        return para.t
    else
        return 0.0
    end
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

    sites = Set(vcat(ri, ro))
    for elem in sites
        idx_i = findall(x -> x == elem, ri)
        idx_o = findall(x -> x == elem, ro)
        if length(idx_i) != length(idx_o)
            return nothing
        end
        element_indices[elem] = (idx_i, idx_o)
    end
    return element_indices
end

@inline function permu_sign(v::Vector{Int})
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
    varT, varRx = vars

    for (i, lftype) in enumerate(leafType[idx])
        if lftype == 0
            continue
        elseif lftype == 4  # BareHoppingId
            τ = varT[leafτ_o[idx][i][1]] - varT[leafτ_i[idx][i][1]]
            r1 = [varRx[leafSites[idx][i][1]], 1]
            r2 = [varRx[leafSites[idx][i][2]], 1]
            orbital = leaforbitals_i[idx][i][1]
            leafval[idx][i] = green_counterterm_PBC(para, r1, r2, orbital) #/ para.β
        elseif lftype == 5  # BareGreenNId
            Np = Int(length(leafSites[idx][i]) / 2)
            sites_i = varRx[leafSites[idx][i][1:Np]]
            sites_o = varRx[leafSites[idx][i][Np+1:end]]
            r_dict = find_duplicates_with_indices(sites_i, sites_o)

            if isnothing(r_dict)
                leafval[idx][i] = 0.0
                continue
            end

            leafval[idx][i] = permu_sign(collect(Iterators.flatten(zip(sites_i, sites_o))))

            τi, τo = varT[leafτ_i[idx][i]], varT[leafτ_o[idx][i]]
            orbitals_i, orbitals_o = leaforbitals_i[idx][i], leaforbitals_o[idx][i]

            for (loc_i, loc_o) in values(r_dict)
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

    return root[1]
end

function density(model, para::ParaMC, diagram; neval=1e6, print=0, dtype=Float64, kwargs...)
    partition, diagpara, FeynGraphs = diagram

    funcGraphs! = Dict{Int,Function}()
    leaf_maps = Vector{Dict{Int,Graph}}()
    for (i, key) in enumerate(partition)
        funcGraphs![i], leafmap = Compilers.compile(FeynGraphs[key])
        push!(leaf_maps, leafmap)
    end

    leafStat = FeynmanDiagram.leafstates(leaf_maps, dtype=dtype)

    root = zeros(dtype, 1)
    T = Continuous(0.0, para.β; offset=2, adapt=true)
    # T = Continuous(0.0, para.β; offset=1, adapt=false)
    # T.data[1] = 0.0
    # T.data[2] = -1e-10
    # T.data[1] = 1e-10
    T.data[1] = 1e-5
    T.data[2] = 0.0
    # R = Discrete(1, para.Lx, adapt=false)
    Rx = Discrete(1, para.Lx, adapt=true)
    Ry = Discrete(1, para.Ly, adapt=true)

    dof = [[p.totalTauNum, p.innerLoopNum * 2 + 1] for p in diagpara]
    obs = zeros(dtype, length(diagpara))
    global_updates = [false, true]

    println("dof: ", dof)

    config = Configuration(; var=(T, Rx), dof=dof, obs=obs, type=dtype, global_updates=global_updates,
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
            datadict[key] = measurement.(avg, std)
        end
        return datadict, result
    else
        return nothing, nothing
    end
end

function density_MC(model, para::ParaMC; neval=1e6, partition=partition(para.order), reweight_goal=nothing,
    print=0, filename::Union{String,Nothing}=nothing, dtype=Float64, _neighbor=nothing)
    diagram = density_info(partition)

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

    _density, result = density(model, para, diagram; neval=neval, neighbor=_neighbor,
        reweight_goal=reweight_goal, dtype=dtype, print=print)

    if isnothing(_density) == false
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (_density,)
            end
        end
        for (ip, key) in enumerate(partition)
            println("Group ", key)
            # @printf("%10s   %10s \n", "avg", "err")
            @printf("%10s   %10s   %10s   %10s \n", "real(avg)", "err", "imag(avg)", "err")
            # @printf("%10.6f ± %10.6f\n", _density[key].val, _density[key].err)
            r, i = real(_density[key]), imag(_density[key])
            @printf("%10.6f ± %10.6f    %10.6f ± %10.6f\n", r.val, r.err, i.val, i.err)
        end
    end
    return _density, result
end