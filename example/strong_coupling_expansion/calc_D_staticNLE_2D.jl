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
    dμ::Float64
    order::Int
end

paraid(p::ParaMC) = Dict(
    "order" => p.order,
    "beta" => p.β,
    "mu" => p.μ,
    "dmu" => p.dμ,
    "U" => p.U,
    "Lx" => p.Lx,
    "Ly" => p.Ly,
)
short(p::ParaMC) = join(["$(k)_$(v)" for (k, v) in sort!(OrderedDict(paraid(p)))], "_")

function hopping_PBC(para::ParaMC, r1::Vector{Int}, r2::Vector{Int}, orbital::Int)
    L = [para.Lx, para.Ly]
    delta12 = abs.(r1 - r2)
    delta = min.(delta12, L .- delta12)

    sum_d = sum(delta)
    if sum_d == 1
        # return para.t - para.dμ
        return para.t
    elseif sum_d == 0
        # return -para.dμ
        return para.dμ
    else
        return 0.0
    end
end

function hopping_FBC(para::ParaMC, r1::Vector{Int}, r2::Vector{Int}, orbital::Int)
    delta = abs.(r1 - r2)

    sum_d = sum(delta)
    if sum_d == 1
        # return para.t - para.dμ
        return para.t
    elseif sum_d == 0
        return -para.dμ
        # return para.dμ
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

# minimum-image displacement on a periodic chain of length L
@inline function min_image(dx::Int, L::Int)
    half = div(L, 2)  # floor(L/2)
    if dx > half
        dx -= L
    elseif dx < -half
        dx += L
    end
    return dx
end

"""
Build an importance-sampling histogram p0 over site indices 1:Ns
for a finite Lx×Ly PBC lattice, given that site 1 is the pinned root.

Arguments:
  Lx, Ly :: Int
  α      :: Float64  (range parameter for exp(-α r))

Returns:
  p0 :: Vector{Float64} of length Ns, normalized to sum(p0)=1.
"""
function build_spatial_histogram(Lx::Int, Ly::Int; α::Float64=1.0)
    Ns = Lx * Ly

    # lattice coordinates for each site index
    site_indices = collect(1:Ns)
    coords = indices_to_lattice(site_indices, Ly)

    # root is index 1
    (x_root, y_root) = coords[1]

    # unnormalized weights
    w = zeros(Float64, Ns)

    @inbounds for s in 1:Ns
        (xs, ys) = coords[s]

        # displacement with periodic wrap (minimum image)
        dx = min_image(xs - x_root, Lx)
        dy = min_image(ys - y_root, Ly)

        # Euclidean distance on torus
        r = sqrt(dx * dx + dy * dy)

        # importance weight: decays with distance
        w[s] = exp(-α * r)
    end

    # normalize to get probabilities
    Z = sum(w)
    if Z == 0.0
        # fallback to uniform just in case
        return fill(1.0 / Ns, Ns)
    else
        return w ./ Z
    end
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

function double_occupancy(model, para::ParaMC, diagram, _neighbor; neval=1e6, print=0, dtype=ComplexF64, kwargs...)
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
    T.data[1] = 0.0
    R = Discrete(1, para.Lx * para.Ly; offset=1, adapt=true, alpha=3.0,
        distribution=build_spatial_histogram(para.Lx, para.Ly))
    R.data[1] = 1
    coords = indices_to_lattice(collect(1:para.Lx*para.Ly), para.Ly)

    dof = [[1, p.totalTauNum - 1, p.innerLoopNum * 2 - 1] for p in diagpara]
    # dof = [[p.totalTauNum, p.innerLoopNum * 2] for p in diagpara]
    obs = zeros(dtype, length(diagpara))
    global_updates = [false, false, true]

    println("dof: ", dof)

    config = Configuration(; var=(T_doublon, T, R), dof=dof, obs=obs, type=dtype,
        global_updates=global_updates, neighbor=_neighbor,
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
            datadict[key] = -measurement.(avg, std)
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

function double_occupancy_MC(model, para::ParaMC; neval=1e6, partition=partition(para.order), reweight_goal=nothing,
    print=0, filename::Union{String,Nothing}=nothing, dtype=ComplexF64, _neighbor=nothing)
    diagram = generate_Gnderiv1(partition, dynamic_hop=false)

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