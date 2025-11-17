push!(LOAD_PATH, pwd())
using Atom
using Lehmann
using MCIntegration
using Printf
using Measurements
using JLD2
using DataStructures
using LinearAlgebra
using Random
using StaticArrays

include("generate_freeE_NLE.jl")

#  Imaginary-time helpers  
@inline function wrap_tau_sign(Δτ::Float64, β::Float64)
    # fermionic antiperiodicity: G(Δτ<0) = -G(Δτ+β)
    if Δτ >= 0.0
        return Δτ, 1.0
    else
        return Δτ + β, -1.0
    end
end

@inline function lininterp(xgrid::AbstractVector{<:Real}, yvals::AbstractVector{<:Real}, x::Real)
    i = searchsortedlast(xgrid, x)
    if i == length(xgrid)
        return yvals[end]
    elseif i == 0
        return yvals[1]
    else
        xL = xgrid[i]
        xR = xgrid[i+1]
        yL = yvals[i]
        yR = yvals[i+1]
        w = (x - xL) / (xR - xL)
        return yL + w * (yR - yL)
    end
end

#  Bare propagator kernels 
function propagator(τ::T, ω::T, β::T) where {T}
    # τ=0⁺ regularization
    if τ ≈ T(0.0)
        τ = -1e-10
    end
    if τ > T(0.0)
        if ω > T(0.0)
            return exp(-ω * τ) / (1 + exp(-ω * β))
        else
            return exp(ω * (β - τ)) / (1 + exp(ω * β))
        end
    else
        # τ<0 branch fallback; normally handled in wrap_tau_sign
        if ω > T(0.0)
            return -exp(-ω * (τ + β)) / (1 + exp(-ω * β))
        else
            return -exp(-ω * τ) / (1 + exp(ω * β))
        end
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
        error("propagator_derivative: order>5 not implemented")
    end
    return result
end

#########################################
#  Momentum mesh / dispersion in k-space
#########################################

"Brillouin-zone sampling mesh: Lkx×Lky points for TL integration. dμ can shift chemical potential."
function dispersion_kmesh(Lkx::Int, Lky::Int, t::Float64, dμ::Float64)
    Norb = 2
    ϵk = zeros(Float64, Norb, Norb, Lkx, Lky)
    kx_arr = zeros(Float64, Lkx, Lky)
    ky_arr = zeros(Float64, Lkx, Lky)

    for ix in 1:Lkx
        for iy in 1:Lky
            kx = 2π * (ix - 1) / Lkx
            ky = 2π * (iy - 1) / Lky
            ek = -2t * (cos(kx) + cos(ky)) + dμ
            ϵk[1, 1, ix, iy] = ek
            ϵk[2, 2, ix, iy] = ek
            kx_arr[ix, iy] = kx
            ky_arr[ix, iy] = ky
        end
    end
    return ϵk, kx_arr, ky_arr
end

#########################################
#  Pre-tabulation in k, then FT to real #
#########################################

"Build g_m(k,τ)/λ_eff for m=0..max_order, all k, all τ."
function build_gmkτ(ϵk, β, λ::Float64, max_order::Int, τgrid::Vector{Float64})
    Norb, _, Lkx, Lky = size(ϵk)
    Ntau = length(τgrid)
    gmkτ = zeros(Float64, Norb, max_order + 1, Lkx, Lky, Ntau)

    # for ix in 1:Lkx, iy in 1:Lky, orb in 1:Norb
    for ix in 1:Lkx, iy in 1:Lky
        # ω = -1.0 / ϵk[orb, orb, ix, iy]      # matches previous convention
        ω = -1.0 / ϵk[1, 1, ix, iy]      # matches previous convention
        λeff = sign(ω) * λ
        ωscaled = ω / λeff
        for (itau, τval) in enumerate(τgrid)
            for m in 0:max_order
                acc = 0.0
                @inbounds for o in 0:m
                    acc += propagator_derivative(τval, ωscaled, β, o) *
                           (ωscaled^o) *
                           binomial(m, o) *
                           (-1)^o
                end
                gmkτ[1:Norb, m+1, ix, iy, itau] .= acc / λeff
            end
        end
    end
    return gmkτ
end

#############################################
#  Displacement bases (Rsample vs Rtable)   #
#############################################

"Generate list of integer displacement vectors in a square box [-R:R]^2."
function make_displacement_list(R::Int)
    ΔRlist = SVector{2,Int}[]
    for dx in -R:R, dy in -R:R
        push!(ΔRlist, SVector{2,Int}(dx, dy))
    end
    ΔRdict = Dict{Tuple{Int,Int},Int}()
    for (idx, ΔR) in enumerate(ΔRlist)
        ΔRdict[(ΔR[1], ΔR[2])] = idx
    end
    return ΔRlist, ΔRdict
end

#############################################
#  Fourier transform k→real for pretab      #
#############################################

"Fourier transform g_m(k,τ) to real-space displacements in ΔRlist_table.
Return Ctable[orb, m+1, ridx, itau]."
function k_to_realspace(gmkτ, ΔRlist_table, kx_arr, ky_arr)
    Norb, Mplus1, Lkx, Lky, Ntau = size(gmkτ)
    Nk = Lkx * Lky
    Nr = length(ΔRlist_table)

    Ctable = zeros(Float64, Norb, Mplus1, Nr, Ntau)

    for (ridx, ΔR) in enumerate(ΔRlist_table)
        dx = ΔR[1]
        dy = ΔR[2]
        # for orb in 1:Norb, m1 in 1:Mplus1, itau in 1:Ntau
        for m1 in 1:Mplus1, itau in 1:Ntau
            acc = 0.0
            @inbounds for ix in 1:Lkx, iy in 1:Lky
                phase = kx_arr[ix, iy] * dx + ky_arr[ix, iy] * dy
                # acc += cos(phase) * gmkτ[orb, m1, ix, iy, itau]
                acc += cos(phase) * gmkτ[1, m1, ix, iy, itau]
            end
            # Ctable[orb, m1, ridx, itau] .= acc / Nk
            Ctable[1:Norb, m1, ridx, itau] .= acc / Nk
        end
    end
    return Ctable
end

#############################################
#  Physics-aware importance distribution    #
#############################################

"Build importance distribution p0 over the *sampling* displacement list (ΔRlist_sample),
using |C_0(ΔR, τ)| integrated over τ.

We need:
- ΔRlist_sample (size NsampR) : the displacements we'll actually sample for vertices
- ΔRdict_table / Ctable : the big-table FT info
- τgrid, β

Returns:
- p0::Vector{Float64}, length NsampR, normalized
"
function build_spatial_importance_from_Ctable(ΔRlist_sample,
    ΔRdict_table::Dict{Tuple{Int,Int},Int},
    Ctable,
    τgrid::Vector{Float64};
    orbital::Int=1, m::Int=0)

    Ntau = length(τgrid)
    β = τgrid[end]  # assuming range starts at 0 and ends at β

    NsampR = length(ΔRlist_sample)
    W = zeros(Float64, NsampR)

    @inbounds for s in 1:NsampR
        dx = ΔRlist_sample[s][1]
        dy = ΔRlist_sample[s][2]
        ridx = get(ΔRdict_table, (dx, dy), nothing)
        if ridx === nothing
            W[s] = 0.0
            continue
        end

        acc = 0.0
        for itau in 1:Ntau
            valτ = Ctable[orbital, m+1, ridx, itau]
            acc += abs(valτ)
        end
        acc *= β / Ntau
        W[s] = acc
    end

    Z = sum(W)
    if Z == 0.0
        return fill(1.0 / NsampR, NsampR)
    else
        return W ./ Z
    end
end

#############################################
#  PreTab + ParaMC structs                  #
#############################################

"All tabulated physics objects needed at runtime."
struct RealSpacePreTab
    τgrid::Vector{Float64}                     # Ntau points in [0,β]
    Ctable::Array{Float64,4}                   # [orb, m+1, ridx_table, itau]
    ΔRlist_table::Vector{SVector{2,Int}}       # displacements for evaluation (big box, Rtable)
    ΔRdict_table::Dict{Tuple{Int,Int},Int}     # (dx,dy) -> ridx_table
end

"Proposal info for Monte Carlo sampling of vertex coordinates."
struct SpatialProposal
    ΔRlist_sample::Vector{SVector{2,Int}}      # allowed relative coords for MC (smaller box, Rsample)
    root_idx::Int
    p0::Vector{Float64}                        # importance distribution over ΔRlist_sample
end

"Simulation parameters for TL strong-coupling DiagMC (pinned root)."
struct ParaMC
    μ::Float64
    U::Float64
    t::Float64
    β::Float64
    lambda::Float64
    dμ::Float64
    order::Int
    pretabs::RealSpacePreTab
    proposal::SpatialProposal
end

paraid(p::ParaMC) = Dict(
    "beta" => p.β,
    "lambda" => p.lambda,
    "mu" => p.μ,
    "dmu" => p.dμ,
    "U" => p.U,
    "order" => p.order,
)

short(p::ParaMC) = join(["$(k)_$(v)" for (k, v) in sort!(OrderedDict(paraid(p)))], "_")

#############################################
#  Fast evaluation of a bare line           #
#############################################

"Lookup C_m(Δr, Δτ) with τ antiperiodicity and τ interpolation.
Δr = r_a - r_b (both relative to root), so Δr is in Rtable box by construction."
function bare_line_fast(para::ParaMC,
    Δr_ab::SVector{2,Int},
    Δτ_ab::Float64,
    orbital::Int,
    m::Int)

    pre = para.pretabs
    β = para.β

    ridx = get(pre.ΔRdict_table, (Δr_ab[1], Δr_ab[2]), nothing)
    ridx === nothing && return 0.0  # outside tabulation range -> negligible

    τ_mod, signfac = wrap_tau_sign(Δτ_ab, β)

    arrτ = @view pre.Ctable[orbital, m+1, ridx, :]
    valτ = lininterp(pre.τgrid, arrτ, τ_mod)

    return signfac * valτ
end

#############################################
#  Permutation sign and duplicate utility   #
#############################################

@inline function permu_sign(v::Vector{Int})
    sgn = 1
    for i in eachindex(v)
        for j in (i+1):length(v)
            if (v[i] > v[j]) || (v[i] == v[j] && iseven(i) && isodd(j))
                sgn *= -1
            end
        end
    end
    return sgn
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

#############################################
#  Coordinate handling in integrand         #
#############################################

# """
# We now treat spatial coordinates as:
# - vertex 1 (root): r_root = (0,0), τ_root = 0
# - all other vertices: sampled relative coords from ΔRlist_sample via MCIntegration
#   (their variable indices live in Discrete(...) domains over that list)

# We assume MCIntegration passes:
#   vars = (T, R)
# where
#   T is a Continuous() vector of times for non-root vertices.
#   R is a Discrete()  vector of indices into ΔRlist_sample for non-root vertices.

# leafStat:
#   leafType[idx][i] encodes what to build:
#     - lftype == 4 : BareHopping-like line with (site a, site b), times, order m, orbital
#     - lftype == 5 : local multi-point GreenN product (kept from your previous code)

# IMPORTANT:
# We must reconstruct absolute coords of each site label:
# - site label = 1 means the root (0,0) and τ=0
# - site label > 1 means pull from R[...] and T[...] appropriately
# """
# function get_vertex_space_time(site_label::Int,
#     varR,
#     varT,
#     ΔRlist_sample::Vector{SVector{2,Int}})
#     if site_label == 1
#         # pinned root
#         return SVector{2,Int}(0, 0), 0.0
#     else
#         idx = site_label - 1
#         r_rel = ΔRlist_sample[varR[idx]]
#         τ_rel = varT[idx]
#         return r_rel, τ_rel
#     end
# end

#############################################
#  integrand for MCIntegration.integrate    #
#############################################

function integrand(idx, vars, config)
    para, root, graphfuncs!, leafStat, model = config.userdata
    leafval, leafType, leafOrders, leafSites,
    leafτ_i, leafτ_o, leaforbitals_i, leaforbitals_o = leafStat

    varT, varR = vars

    ΔRlist_sample = para.proposal.ΔRlist_sample

    for (i, lftype) in enumerate(leafType[idx])
        if lftype == 0
            continue
        elseif lftype == 4            # BareHoppingId / counterterm-like line
            r_o = ΔRlist_sample[varR[leafSites[idx][i][1]]]
            r_i = ΔRlist_sample[varR[leafSites[idx][i][2]]]
            τ_o, τ_i = varT[leafτ_o[idx][i][1]], varT[leafτ_i[idx][i][1]]

            if Rsample in abs.(r_o) || Rsample in abs.(r_i)
                @warn "bare_line_fast: vertex displacement outside Rsample box!"
            end

            Δr = r_o .- r_i
            Δτ = τ_o - τ_i
            m_order = leafOrders[idx][i][1]
            orb = leaforbitals_i[idx][i][1]
            leafval[idx][i] = bare_line_fast(para, Δr, Δτ, orb, m_order)
        elseif lftype == 5      # BareGreenNId: product of local GreenN blocks
            # We must build τ arrays and orbital arrays for matched sites.
            Np = Int(length(leafSites[idx][i]) ÷ 2)
            sites_i = varR[leafSites[idx][i][1:Np]]
            sites_o = varR[leafSites[idx][i][Np+1:end]]

            # map each site label -> "index in varT/varR or root"
            # BUT: Green.GreenN expects absolute τ's, not differences. We defined τ_root=0, τ_other=varT[j] as absolute times already, so we can just read them.
            r_dict = find_duplicates_with_indices(sites_i, sites_o)

            if isnothing(r_dict)
                leafval[idx][i] = 0.0
                continue
            end

            leafval[idx][i] = permu_sign(collect(Iterators.flatten(zip(sites_i, sites_o))))

            τs_i, τs_o = varT[leafτ_i[idx][i]], varT[leafτ_o[idx][i]]
            orbitals_i, orbitals_o = leaforbitals_i[idx][i], leaforbitals_o[idx][i]

            for (loc_i, loc_o) in values(r_dict)
                τarr = vcat(τs_i[loc_i], τs_o[loc_o])
                orbs = vcat(orbitals_i[loc_i], orbitals_o[loc_o])
                _gn = Green.GreenN(model, τarr, orbs)
                leafval[idx][i] *= Green.Gn(model, _gn)
            end

        else
            error("integrand: leafType $lftype not implemented")
        end
    end

    graphfuncs![idx](root, leafval[idx])
    return root[1]
end

#############################################
#  Driver to evaluate free energy group     #
#############################################

function freeE(model, para::ParaMC, diagram;
    neval=1e6,
    print=0,
    dtype=ComplexF64,
    kwargs...)

    partition, diagpara, FeynGraphs = diagram

    # Compile diagram groups
    funcGraphs! = Dict{Int,Function}()
    leaf_maps = Vector{Dict{Int,Graph}}()
    for (i, key) in enumerate(partition)
        funcGraphs![i], leafmap = Compilers.compile(FeynGraphs[key])
        push!(leaf_maps, leafmap)
    end

    leafStat = FeynmanDiagram.leafstates(leaf_maps, dtype=dtype)

    root = zeros(dtype, 1)

    # ---- Coordinate/time variables for MCIntegration ----
    # This matches your previous pattern, except now "space variables" are *indices into ΔRlist_sample*, not full site indices.

    dof = [[p.totalTauNum - 1,   # number of non-root τ's
        p.innerLoopNum * 2 - 1]  # number of non-root spatial labels 
           for p in diagpara]

    # Continuous times for non-root vertices:
    T = Continuous(0.0, para.β; offset=1, adapt=true, alpha=3.0)
    T.data[1] = 0.0  # root τ=0 fixed analytically, not sampled;

    # Discrete spatial displacement indices for non-root vertices:
    # Domain is 1 : length(ΔRlist_sample)
    R = Discrete(1, length(para.proposal.ΔRlist_sample); offset=1, adapt=true, alpha=3.0,
        distribution=para.proposal.p0)
    R.data[1] = para.proposal.root_idx   # root r=0 fixed analytically, not sampled;

    obs = zeros(dtype, length(diagpara))

    config = Configuration(; var=(T, R), dof=dof, obs=obs, type=dtype,
        userdata=(para, root, funcGraphs!, leafStat, model))

    result = integrate(integrand;
        config=config,
        neval=neval,
        print=print,
        solver=:mcmc,
        kwargs...)

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

#############################################
#  Top-level convenience wrapper            #
#############################################

function freeE_MC(model, para::ParaMC;
    neval=1e6,
    partition=partition(para.order),
    reweight_goal=nothing,
    print=0,
    filename::Union{String,Nothing}=nothing,
    dtype=ComplexF64,
    _neighbor=nothing)

    diagram = free_energy(partition)
    partition = diagram[1]
    println("partition: ", partition)

    # Same logic as before for MCMC tuning defaults etc.
    if isnothing(reweight_goal)
        reweight_goal = Float64[]
        for (_order, sOrder) in partition
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

    freeEnergy, result = freeE(model, para, diagram;
        neval=neval,
        neighbor=_neighbor,
        reweight_goal=reweight_goal,
        dtype=dtype,
        print=print)

    if freeEnergy !== nothing
        if filename !== nothing
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
            @printf("%10s   %10s   %10s   %10s \n",
                "real(avg)", "err", "imag(avg)", "err")
            r, i = real(freeEnergy[key]), imag(freeEnergy[key])
            @printf("%10.6f ± %10.6f    %10.6f ± %10.6f\n",
                r.val, r.err, i.val, i.err)
        end
    end

    return freeEnergy, result
end

#############################################
#  Builder: construct ParaMC end-to-end     #
#############################################

"""
build_paraMC_TL_SC(
    μ, U, t, β;
    lambda, dμ, order,
    Lkx, Lky,
    Rsample, Rtable,
    Ntau
)

Creates:
- τ grid
- k-mesh dispersion
- g_m(k,τ)
- Fourier transform to Ctable on displacement box Rtable
- proposal distribution p0 on Rsample box
- ParaMC with pinned-root formulation
"""
function build_paraMC_TL_SC(; μ::Float64,
    U::Float64,
    t::Float64,
    β::Float64,
    lambda::Float64,
    dμ::Float64,
    order::Int,
    deriv_order::Int,
    Lkx::Int,
    Lky::Int,
    Rsample::Int,
    Rtable::Int,
    Ntau::Int)

    @assert Rtable >= 2 * Rsample "Need Rtable ≥ 2*Rsample so Δr_ab stays tabulated"

    # time grid
    τgrid = collect(range(0.0, β; length=Ntau))

    # dispersion & k-mesh
    ϵk, kx_arr, ky_arr = dispersion_kmesh(Lkx, Lky, t, dμ)

    # build g_m(k,τ)
    gmkτ = build_gmkτ(ϵk, β, lambda, deriv_order, τgrid)

    # displacement sets
    ΔRlist_table, ΔRdict_table = make_displacement_list(Rtable)
    ΔRlist_sample, ΔRdict_sample = make_displacement_list(Rsample)

    # FT to real space on full table range
    Ctable = k_to_realspace(gmkτ, ΔRlist_table, kx_arr, ky_arr)

    pretabs = RealSpacePreTab(τgrid, Ctable, ΔRlist_table, ΔRdict_table)

    # build spatial importance distribution over the *sampling* box
    p0 = build_spatial_importance_from_Ctable(ΔRlist_sample,
        ΔRdict_table,
        Ctable,
        τgrid;
        orbital=1,
        m=0)

    proposal = SpatialProposal(ΔRlist_sample, ΔRdict_sample[(0, 0)], p0)

    return ParaMC(μ, U, t, β, lambda, dμ, order, pretabs, proposal)
end
