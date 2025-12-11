using StaticArrays
using LinearAlgebra
using Printf
using Lehmann
using FFTW

# --- Basic Physics Kernels (Unchanged) ---

@inline function wrap_tau_sign(Δτ::Float64, β::Float64)
    if Δτ >= 0.0
        return Δτ, 1.0
    else
        return Δτ + β, -1.0
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

# --- Generalized Storage Structure ---

struct RealSpacePreTab
    is_uniform::Bool                   # Flag to enable fast path
    dτ::Float64                        # Grid spacing (used only if is_uniform)
    τgrid::Vector{Float64}             # Full grid (used if !is_uniform)

    β::Float64
    # Layout: [τ, orbital, order, spatial_index]
    Ctable::Array{Float64,4}
    Rtable::Int
    # width::Int
end

# --- Construction ---

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

# Mapping (dx, dy) -> Compact Index
# (dx, dy) -> 0 <= y <= x <= R 
@inline function get_compact_index(dx::Int, dy::Int)
    ax = abs(dx)
    ay = abs(dy)

    if ay > ax
        ax, ay = ay, ax
    end

    # Triangular Number Indexing
    return (ax * (ax + 1)) >> 1 + ay + 1
end

"""
    build_pretab(...)

Builds the lookup table for an arbitrary τgrid. 
Automatically detects if τgrid is uniform to enable O(1) lookups.
"""
function build_pretab(t, β, τgrid::Vector{Float64};
    lambda, dμ, deriv_order, Lkx, Lky, Rtable)

    println("Info: Starting memory-efficient pre-tabulation...")

    Ntau = length(τgrid)
    Norb = 2
    Mplus1 = deriv_order + 1
    # Sum_{x=0 to R} (x+1) = (R+1)(R+2)/2
    Nr = ((Rtable + 1) * (Rtable + 2)) ÷ 2

    mem_gb = Ntau * Norb * Mplus1 * Nr * 8 / 1024^3
    println("Info: Rtable=$Rtable. Storing $Nr spatial points (vs original $(2*Rtable+1)^2).")
    println("Info: Final Ctable Size: $(round(mem_gb, digits=2)) GB")

    # Check uniformity with a small tolerance
    dτ = τgrid[2] - τgrid[1]
    is_uniform = all(isapprox.(diff(τgrid), dτ; atol=1e-9))

    Ctable = zeros(Float64, Ntau, Norb, Mplus1, Nr)

    ϵk, _, _ = dispersion_kmesh(Lkx, Lky, t, dμ)
    Gk_buffer = zeros(ComplexF64, Lkx, Lky)
    p_ifft = plan_ifft(Gk_buffer)

    println("Info: Computing & FFT per tau point...")

    # pre-calculate ω and λeff
    ω_arr = zeros(Float64, Lkx, Lky)
    λeff_arr = zeros(Float64, Lkx, Lky)
    ωscaled_arr = zeros(Float64, Lkx, Lky)
    for ix in 1:Lkx, iy in 1:Lky
        ω = -1.0 / ϵk[1, 1, ix, iy]
        λ = sign(ω) * lambda
        ω_arr[ix, iy] = ω
        λeff_arr[ix, iy] = λ
        ωscaled_arr[ix, iy] = ω / λ
    end

    for (itau, τval) in enumerate(τgrid)
        if itau % 1000 == 0
            print("\rProgress: $itau / $Ntau")
        end

        for m1 in 1:Mplus1 # m1 = order + 1
            m = m1 - 1

            # Assume orb=1,1 and orb=2,2 symmetric
            @inbounds for iy in 1:Lky, ix in 1:Lkx
                acc = 0.0
                ωsc = ωscaled_arr[ix, iy]

                term = 1.0
                for o in 0:m
                    val = propagator_derivative(τval, ωsc, β, o)
                    # acc += val * (ωsc^o) * binomial(m, o) * (-1)^o
                    acc += val * term * binomial(m, o)
                    term *= -ωsc
                end
                Gk_buffer[ix, iy] = acc / λeff_arr[ix, iy]
            end

            # B. Perform FFT (k -> r)
            Gr_full = p_ifft * Gk_buffer

            # C. Truncate and fill Ctable
            # Ctable layout: [itau, orb, m1, ridx]
            idx_compact = 0
            for dx in 0:Rtable
                for dy in 0:dx
                    # Map dx, dy to FFT output indices (1-based, periodic)
                    idx_x = dx + 1
                    idx_y = dy + 1

                    val_r = real(Gr_full[idx_x, idx_y])

                    # map to linear index
                    idx_compact += 1

                    Ctable[itau, 1, m1, idx_compact] = val_r
                    Ctable[itau, 2, m1, idx_compact] = val_r
                end
            end
        end
    end
    println("\nInfo: Pre-tabulation finished.")

    return RealSpacePreTab(is_uniform, dτ, τgrid, β, Ctable, Rtable)
end

# --- The Generalized Hot Path ---

"""
    bare_line_fast(pretab, Δr, Δτ, orbital, m)

Computes hopping counterterm. Uses O(1) arithmetic if grid is uniform,
otherwise uses O(log N) binary search.
"""
@inline function bare_line_fast(pre::RealSpacePreTab,
    Δr::SVector{2,Int},
    Δτ::Float64,
    orbital::Int,
    m::Int)

    # 1. Spatial Check (Early Exit)
    R = pre.Rtable
    dx, dy = Δr[1], Δr[2]
    ax = abs(dx)
    ay = abs(dy)
    if ax > R || ay > R
        return 0.0
    end

    # Spatial Indexing
    if ay > ax
        ax, ay = ay, ax
    end
    ridx = (ax * (ax + 1)) >> 1 + ay + 1

    # 2. Time Wrapping
    τ_mod, signfac = wrap_tau_sign(Δτ, pre.β)

    # 3. Find Grid Index & Weight
    idx = 0
    w = 0.0
    Ntau = size(pre.Ctable, 1)

    if pre.is_uniform
        # --- Fast Path (Arithmetic) ---
        x = τ_mod / pre.dτ
        idx = floor(Int, x) + 1

        # Clamp upper bound (if τ_mod ≈ β)
        if idx >= Ntau
            idx = Ntau - 1
        end

        # Weight (x - floor(x))
        w = (τ_mod - (idx - 1) * pre.dτ) / pre.dτ
    else
        # --- General Path (Binary Search) ---
        # Returns index of the last value <= τ_mod
        idx = searchsortedlast(pre.τgrid, τ_mod)

        # Handle boundaries
        if idx == 0
            idx = 1 # Should not happen if τgrid[1] == 0.0 and τ_mod >= 0
        elseif idx >= Ntau
            idx = Ntau - 1
        end

        # Calculate weight based on specific grid points
        @inbounds begin
            t_L = pre.τgrid[idx]
            t_R = pre.τgrid[idx+1]
        end
        w = (τ_mod - t_L) / (t_R - t_L)
    end

    # 4. Interpolate
    # Ctable layout: [τ, orb, m, r]
    @inbounds begin
        yL = pre.Ctable[idx, orbital, m+1, ridx]
        yR = pre.Ctable[idx+1, orbital, m+1, ridx]
    end

    val = yL + w * (yR - yL)

    return signfac * val
end