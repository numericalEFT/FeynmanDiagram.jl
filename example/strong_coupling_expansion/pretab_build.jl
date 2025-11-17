############ pretab_build.jl ############
# Build and cache TL strong-coupling pretabulation for DDMC.

using FFTW
using Lehmann
using StaticArrays
using LinearAlgebra
using JLD2

########## Physics helper functions ##########

@inline function wrap_tau_sign(Δτ::Float64, β::Float64)
    if Δτ >= 0.0
        return Δτ, 1.0
    else
        return Δτ + β, -1.0
    end
end

@inline function propagator(τ::T, ω::T, β::T) where {T}
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

########## Displacement bases ##########

"Generate list of integer displacement vectors in [-R,R]^2 and the lookup dict."
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

########## k-space dispersion / mesh ##########

"Non-interacting dispersion ε_k on a uniform k-mesh, Hubbard-like."
function dispersion_kmesh(Lkx::Int, Lky::Int, t::Float64, dμ::Float64)
    Norb = 2
    ϵk = zeros(Float64, Norb, Norb, Lkx, Lky)
    kx_arr = zeros(Float64, Lkx, Lky)
    ky_arr = zeros(Float64, Lkx, Lky)

    for ix in 1:Lkx, iy in 1:Lky
        kx = 2π * (ix - 1) / Lkx
        ky = 2π * (iy - 1) / Lky
        ek = -2t * (cos(kx) + cos(ky)) + dμ
        ϵk[1, 1, ix, iy] = ek
        ϵk[2, 2, ix, iy] = ek
        kx_arr[ix, iy] = kx
        ky_arr[ix, iy] = ky
    end
    return ϵk, kx_arr, ky_arr
end

########## Build g_m(k,τ) ##########

"Compute g_m(k,τ)/λ_eff for m=0..max_order, all orbitals, all k, all τ."
function build_gmkτ(ϵk, β, λ::Float64, max_order::Int, τgrid::Vector{Float64})
    Norb, _, Lkx, Lky = size(ϵk)
    Ntau = length(τgrid)
    gmkτ = zeros(Float64, Norb, max_order + 1, Lkx, Lky, Ntau)

    for ix in 1:Lkx, iy in 1:Lky, orb in 1:Norb
        ω = -1.0 / ϵk[orb, orb, ix, iy]
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
                gmkτ[orb, m+1, ix, iy, itau] = acc / λeff
            end
        end
    end
    return gmkτ
end

########## k -> real-space via FFT ##########

"""
Compute real-space Ctable from gmkτ with a batched 2D inverse FFT.

Input:
  gmkτ[orb, m+1, kx, ky, itau]
Output:
  Ctable[orb, m+1, ridx, itau] with ridx indexing ΔRlist_table.

We FFT to get a full real-space grid G[orb,m,rx,ry,τ], then pick only
the displacements we care about via ΔRlist_table and periodic wrap.
"""
function k_to_realspace_fft(gmkτ,
    ΔRlist_table::Vector{SVector{2,Int}},
    ΔRdict_table::Dict{Tuple{Int,Int},Int})

    Norb, Mplus1, Lkx, Lky, Ntau = size(gmkτ)
    Nr = length(ΔRlist_table)

    G = ComplexF64.(gmkτ)  # copy so FFTW can work in-place
    plan = plan_ifft!(G, (3, 4))  # inverse FFT along kx,ky dims

    plan * G  # now G is (unnormalized? no, FFTW.ifft! includes 1/(Lkx*Lky))

    # Allocate compact table
    Ctable = zeros(Float64, Norb, Mplus1, Nr, Ntau)

    @inbounds for ridx in 1:Nr
        dx = ΔRlist_table[ridx][1]
        dy = ΔRlist_table[ridx][2]
        rx = mod(dx, Lkx) + 1
        ry = mod(dy, Lky) + 1
        for orb in 1:Norb, m1 in 1:Mplus1, itau in 1:Ntau
            Ctable[orb, m1, ridx, itau] = real(G[orb, m1, rx, ry, itau])
        end
    end

    return Ctable
end

########## Spatial importance distribution ##########

"""
Build importance distribution p0 over ΔRlist_sample, based on |C_0(ΔR,τ)| integrated over τ.
Returns Vector p0 normalized to sum=1.
"""
function build_spatial_importance_from_Ctable(ΔRlist_sample,
    ΔRdict_table::Dict{Tuple{Int,Int},Int},
    Ctable,
    τgrid::Vector{Float64};
    orbital::Int=1,
    m::Int=0)

    Ntau = length(τgrid)
    β = τgrid[end]
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
    return Z == 0.0 ? fill(1.0 / NsampR, NsampR) : (W ./ Z)
end

########## Structs we'll serialize ##########

struct RealSpacePreTab
    τgrid::Vector{Float64}
    Ctable::Array{Float64,4}                   # [orb, m+1, ridx_table, itau]
    ΔRlist_table::Vector{SVector{2,Int}}       # full "evaluation box"
    ΔRdict_table::Dict{Tuple{Int,Int},Int}     # (dx,dy)->ridx in table
end

struct SpatialProposal
    ΔRlist_sample::Vector{SVector{2,Int}}      # "sampling box"
    p0::Vector{Float64}                        # importance weights over ΔRlist_sample
end

"Parameter signature for cache safety / ParaMC reconstruction."
function pretab_signature(; μ, U, t, β, lambda, dμ,
    order, Lkx, Lky,
    Rsample, Rtable, Ntau)
    return Dict(
        "mu" => μ,
        "U" => U,
        "t" => t,
        "beta" => β,
        "lambda" => lambda,
        "dmu" => dμ,
        "order" => order,
        "Lkx" => Lkx,
        "Lky" => Lky,
        "Rsample" => Rsample,
        "Rtable" => Rtable,
        "Ntau" => Ntau,
    )
end

########## Build + Save pipeline ##########

"""
build_pretab_TL_SC:
Generate (sig, pretabs, proposal) for TL strong-coupling DiagMC.

Rtable is set internally to 2*Rsample.
"""
function build_pretab_TL_SC(; μ::Float64,
    U::Float64,
    t::Float64,
    β::Float64,
    lambda::Float64,
    dμ::Float64,
    order::Int,
    Lkx::Int,
    Lky::Int,
    Rsample::Int,
    Ntau::Int)

    Rtable = 2 * Rsample

    # time grid
    τgrid = collect(range(0.0, β; length=Ntau))

    # dispersion in k-space
    ϵk, kx_arr, ky_arr = dispersion_kmesh(Lkx, Lky, t, dμ)

    # g_m(k,τ)
    gmkτ = build_gmkτ(ϵk, β, lambda, order, τgrid)

    # displacement lists
    ΔRlist_table, ΔRdict_table = make_displacement_list(Rtable)
    ΔRlist_sample, _dummy = make_displacement_list(Rsample)

    # k -> real space (bare object tabulation)
    Ctable = k_to_realspace_fft(gmkτ, ΔRlist_table, ΔRdict_table)

    pretabs = RealSpacePreTab(τgrid, Ctable, ΔRlist_table, ΔRdict_table)

    # importance proposal over Rsample box
    p0 = build_spatial_importance_from_Ctable(
        ΔRlist_sample,
        ΔRdict_table,
        Ctable,
        τgrid;
        orbital=1,
        m=0,
    )
    proposal = SpatialProposal(ΔRlist_sample, p0)

    sig = pretab_signature(; μ, U, t, β, lambda, dμ,
        order, Lkx, Lky,
        Rsample, Rtable, Ntau)

    return sig, pretabs, proposal
end

"Save pretab + proposal + signature to a .jld2 file."
function save_pretab_jld2(filename::String,
    sig::Dict{String,Any},
    pretabs::RealSpacePreTab,
    proposal::SpatialProposal)

    jldopen(filename, "w") do f
        f["signature"] = sig
        f["τgrid"] = pretabs.τgrid
        f["Ctable"] = pretabs.Ctable
        f["ΔRlist_table"] = pretabs.ΔRlist_table
        f["ΔRdict_table"] = pretabs.ΔRdict_table
        f["ΔRlist_sample"] = proposal.ΔRlist_sample
        f["p0"] = proposal.p0
    end
end

########################################################
# EXAMPLE USAGE (run this file as a script standalone) #
########################################################

if abspath(PROGRAM_FILE) == @__FILE__
    # choose parameters
    μ = -0.5
    U = 4.0
    t = 1.0
    β = 5.0
    lambda = 1.0
    dμ = 0.0
    order = 4
    Lkx = 64
    Lky = 64
    Rsample = 4
    Ntau = 2000

    sig, pretabs, proposal = build_pretab_TL_SC(; μ, U, t, β,
        lambda, dμ,
        order,
        Lkx, Lky,
        Rsample,
        Ntau)

    outfile = @sprintf("pretab_mu%.3f_U%.3f_beta%.3f_R%d.jld2",
        μ, U, β, Rsample)

    println("Saving pretab to $outfile ...")
    save_pretab_jld2(outfile, sig, pretabs, proposal)
    println("Done.")
end
