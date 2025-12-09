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
    model, coords = config.userdata[5:6]
    varT_D, varT, varRx = vars
    τp = varT_D[1]

    for (i, lftype) in enumerate(leafType[idx])
        if lftype == 0
            continue
        elseif lftype == 4  # BareHoppingId
            τ = varT[leafτ_o[idx][i][1]] - varT[leafτ_i[idx][i][1]]
            r1 = coords[varRx[leafSites[idx][i][1]]]
            r2 = coords[varRx[leafSites[idx][i][2]]]

            order = leafOrders[idx][i][1]
            orbital = leaforbitals_i[idx][i][1]
            leafval[idx][i] = hopping_counterterm_PBC(para, τ, r1, r2, orbital, order)
        elseif lftype == 5  # GreenNId
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

            order = leafOrders[idx][i][2]
            if order == 0
                for (loc_i, loc_o) in values(r_dict)
                    τ = vcat(τi[loc_i], τo[loc_o])
                    orbitals = vcat(orbitals_i[loc_i], orbitals_o[loc_o])

                    _gn = Green.GreenN(model, τ, orbitals)
                    leafval[idx][i] *= Green.Gn(model, _gn)
                end
            elseif order == 1
                len_greenN = length(r_dict)
                Gvec = Vector{Float64}(undef, len_greenN)
                dGvec_dU = Vector{Float64}(undef, len_greenN)
                if len_greenN == 1
                    τ = vcat(τi, τo)
                    orbitals = vcat(orbitals_i, orbitals_o)

                    _gn = Green.GreenN(model, τ, orbitals)
                    leafval[idx][i] *= Green.dGn_dU_estimator(model, _gn, τp)
                else
                    for (i, (loc_i, loc_o)) in enumerate(values(r_dict))
                        τ = vcat(τi[loc_i], τo[loc_o])
                        orbitals = vcat(orbitals_i[loc_i], orbitals_o[loc_o])

                        _gn = Green.GreenN(model, τ, orbitals)
                        Gvec[i] = Green.Gn(model, _gn)
                        dGvec_dU[i] = Green.dGn_dU_estimator(model, _gn, τp)
                    end
                    leafval[idx][i] *= deriv_prod(Gvec, dGvec_dU)
                end
            else
                error("this order $order not implemented!")
            end
        else
            error("this leaftype $lftype not implemented!")
        end
    end

    graphfuncs![idx](root, leafval[idx])

    return root[1]
end

function indices_to_lattice(indices::AbstractVector{Int}, L::Int)
    n = length(indices)
    coords = Vector{Vector{Int}}(undef, n)

    @inbounds for i in 1:n
        idx = indices[i] - 1
        rx = (idx ÷ L) + 1
        ry = (idx % L) + 1
        coords[i] = [rx, ry]
    end
    return coords
end

"""
Calculates the derivative of prod(A) with respect to U.

Arguments:
- A: The vector of values.
- dA_dU: The vector of derivatives of each element in A w.r.t. U.
"""
function deriv_prod(A::AbstractVector, dA_dU::AbstractVector)
    @assert length(A) == length(dA_dU) "Vectors must have the same length"

    # Find the indices of zero elements
    zero_indices = findall(iszero, A)
    num_zeros = length(zero_indices)

    if num_zeros == 0
        # --- Case 1: No zeros ---
        # Use the stable log-derivative formula
        return prod(A) * sum(dA_dU ./ A)
    elseif num_zeros == 1
        # --- Case 2: Exactly one zero ---
        k = zero_indices[1] # Get the index of the single zero

        # We need the product of all non-zero elements.
        # This is an efficient way to calculate it without allocating new arrays from slicing.
        p_others = 1.0
        for i in eachindex(A)
            if i != k
                p_others *= A[i]
            end
        end
        return p_others * dA_dU[k]
    else
        # --- Case 3: Two or more zeros ---
        return 0.0
    end
end

function double_occupancy(model, para::ParaMC, diagram; neval=1e6, print=0, dtype=ComplexF64, kwargs...)
    partition, diagpara, FeynGraphs = diagram

    funcGraphs! = Dict{Int,Function}()
    leaf_maps = Vector{Dict{Int,Graph}}()
    for (i, key) in enumerate(partition)
        funcGraphs![i], leafmap = Compilers.compile(FeynGraphs[key])
        push!(leaf_maps, leafmap)
    end

    leafStat = FeynmanDiagram.leafstates(leaf_maps, dtype=dtype)

    root = zeros(dtype, 1)
    T_doublon = Continuous(0.0, para.β; adapt=true, alpha=3.0)
    T = Continuous(0.0, para.β; offset=1, adapt=true, alpha=3.0)
    # T = Continuous(0.0, para.β; offset=1, adapt=false)
    T.data[1] = 0.0
    # T = Continuous(0.0, para.β; adapt=true)
    # R = Discrete(1, para.Lx, adapt=false)
    R = Discrete(1, para.Lx * para.Ly, offset=1, adapt=true, alpha=3.0)
    R.data[1] = 1
    coords = indices_to_lattice(collect(1:para.Lx*para.Ly), para.Ly)

    dof = [[1, p.totalTauNum - 1, p.innerLoopNum * 2 - 1] for p in diagpara]
    # dof = [[p.totalTauNum, p.innerLoopNum * 2] for p in diagpara]
    obs = zeros(dtype, length(diagpara))

    global_updates = [false, false, true]
    println("dof: ", dof)

    # config = Configuration(; var=(T_doublon, T, R), dof=dof, obs=obs, type=dtype,
    config = Configuration(; var=(T_doublon, T, R), dof=dof, obs=obs, type=dtype, global_updates=global_updates,
        userdata=(para, root, funcGraphs!, leafStat, model, coords))
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
            # datadict[key] = -measurement.(avg, std) / (para.Lx * para.Ly)
            datadict[key] = -measurement.(avg, std)
        end
        return datadict, result
    else
        return nothing, nothing
    end
end

function double_occupancy_MC(model, para::ParaMC; neval=1e6, partition=partition(para.order), reweight_goal=nothing,
    print=0, filename::Union{String,Nothing}=nothing, dtype=ComplexF64, _neighbor=nothing)
    diagram = generate_Gnderiv1(partition)

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

    Dloc = Green.thermal_expectation(model, model.D)
    println("The local double occupancy (0-th order) is: ", Dloc)

    doublon, result = double_occupancy(model, para, diagram; neval=neval, neighbor=_neighbor,
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