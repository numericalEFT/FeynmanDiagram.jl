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

include("generate_freeE_NLE.jl")

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

function hopping_PBC(para::ParaMC, r1::Vector{Int}, r2::Vector{Int}, orbital::Int)
    L = [para.Lx, para.Ly]
    delta12 = abs.(r1 - r2)
    delta = min.(delta12, L .- delta12)


    if sum(delta) == 1 || all(delta .== 1)
        return para.t
    else
        return 0.0
    end
end

function hopping_FBC(para::ParaMC, r1::Vector{Int}, r2::Vector{Int}, orbital::Int)
    delta = abs.(r1 - r2)

    if sum(delta) == 1 || all(delta .== 1)
        return para.t
    else
        return 0.0
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
            leafval[idx][i] = hopping_PBC(para, r1, r2, orbital)
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

function indices_to_triangular_lattice(indices::AbstractVector{Int}, L::Int)
    n = length(indices)
    coords = Vector{Vector{Int}}(undef, n)

    @inbounds for i in 1:n
        # Convert from 1-based index to 0-based for math
        idx = indices[i] - 1

        # Calculate the integer coordinates in the oblique basis
        u = (idx % L) + 1  # Coordinate along the first basis vector (horizontal)
        v = (idx ÷ L) + 1  # Coordinate along the second basis vector (skewed)

        coords[i] = [u, v]
    end
    return coords
end

function freeE(model, para::ParaMC, diagram; neval=1e6, print=0, dtype=ComplexF64, kwargs...)
    partition, diagpara, FeynGraphs = diagram

    funcGraphs! = Dict{Int,Function}()
    leaf_maps = Vector{Dict{Int,Graph}}()
    for (i, key) in enumerate(partition)
        funcGraphs![i], leafmap = Compilers.compile(FeynGraphs[key])
        push!(leaf_maps, leafmap)
    end

    leafStat = FeynmanDiagram.leafstates(leaf_maps, dtype=dtype)

    root = zeros(dtype, 1)
    T = Continuous(0.0, para.β; offset=1, adapt=true)
    # T = Continuous(0.0, para.β; offset=1, adapt=false)
    T.data[1] = 0.0
    # T = Continuous(0.0, para.β; adapt=true)
    R = Discrete(1, para.Lx * para.Ly, adapt=true)
    coords = indices_to_triangular_lattice(collect(1:para.Lx*para.Ly), para.Lx)

    dof = [[p.totalTauNum - 1, p.innerLoopNum * 2] for p in diagpara]
    # dof = [[p.totalTauNum, p.innerLoopNum * 2] for p in diagpara]
    obs = zeros(dtype, length(diagpara))
    global_updates = [false, true]

    println("dof: ", dof)

    config = Configuration(; var=(T, R), dof=dof, obs=obs, type=dtype, global_updates=global_updates,
        userdata=(para, root, funcGraphs!, leafStat, model, coords))
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
            datadict[key] = -measurement.(avg, std) / (para.Lx * para.Ly)
            # r = measurement.(real(avg), real(std))
            # i = measurement.(imag(avg), imag(std))
            # data = Complex.(r, i)
            # datadict[key] = data
        end
        return datadict, result
    else
        return nothing, nothing
    end
end

function freeE_MC(model, para::ParaMC; neval=1e6, partition=partition(para.order), reweight_goal=nothing,
    print=0, filename::Union{String,Nothing}=nothing, dtype=ComplexF64, _neighbor=nothing)
    diagram = free_energy(partition, dynamic_hop=false)

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