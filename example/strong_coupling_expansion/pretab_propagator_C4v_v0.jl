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
        result = Spectral.kernelFermiT_dω(τ, ϵ, β)
    elseif order == 2
        result = Spectral.kernelFermiT_dω2(τ, ϵ, β) / 2.0
    elseif order == 3
        result = Spectral.kernelFermiT_dω3(τ, ϵ, β) / 6.0
    elseif order == 4
        result = Spectral.kernelFermiT_dω4(τ, ϵ, β) / 24.0
    elseif order == 5
        result = Spectral.kernelFermiT_dω5(τ, ϵ, β) / 120.0
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

@inline function propagator_lamderiv!(Gk_buffer, Lkx, Lky, m, τval, λeff_arr, ωscaled_arr, β)
    @inbounds for iy in 1:Lky, ix in 1:Lkx
        acc = 0.0
        ωsc = ωscaled_arr[ix, iy]

        term = 1.0
        for o in 0:m
            val = propagator_derivative(τval, ωsc, β, o)
            acc += val * term * binomial(m, o)
            term *= ωsc
        end
        Gk_buffer[ix, iy] = acc / λeff_arr[ix, iy]
    end
end

@inline function propagator_muderiv!(Gk_buffer, Lkx, Lky, m, τval, λeff_arr, ωscaled_arr, β)
    fact_m = factorial(m)
    @inbounds for iy in 1:Lky, ix in 1:Lkx
        acc = 0.0
        ωsc = ωscaled_arr[ix, iy]

        if m == 0
            Gk_buffer[ix, iy] = propagator_derivative(τval, ωsc, β, 0) / λeff_arr[ix, iy]
            continue
        end

        term = 1.0
        for o in 1:m
            val = propagator_derivative(τval, ωsc, β, o)
            term *= ωsc
            acc += val * term * binomial(m - 1, o - 1)
        end
        Gk_buffer[ix, iy] = acc * λeff_arr[ix, iy]^(m - 1) * (-ωsc)^m * fact_m
    end
end

# 辅助函数：计算 d^M/dμ^M [ (ωsc)^k * G^{(n)} ]
# 注意：返回值已经包含了适当的 λ 因子，使得在主函数最后除以 λ 时能给出正确量纲。
# 这里的缩放逻辑是：返回值的量级为 λ^M (相对于无量纲量)，
# 这样主函数除以 λ 后，M=0, L=0 的项就是 1/λ。
@inline function calc_mu_deriv_component(M, k, n, τval, ωsc, β, λ)
    # M: mu 导数阶数
    # k: ωsc 的幂次 (来自 lambda 导数展开)
    # n: G 的导数阶数

    # === 情况 1: M = 0 (无 mu 导数) ===
    if M == 0
        # 此时只需要返回 (ωsc)^k * G^{(n)}
        # 这里的 scaling 是 1.0 (即 λ^0)。
        # 主函数最后除以 λ，得到 λ^{-1}，符合 M=0, L=0 的 propagator 定义。
        return (ωsc^k) * propagator_derivative(τval, ωsc, β, n)
    end

    # === 情况 2: M > 0 (应用莱布尼茨法则) ===
    # 目标：计算 d^M/dμ^M [ (ωsc)^k * G^{(n)} ]

    res = 0.0

    # 莱布尼茨展开: sum_{j=0}^{M} binomial(M, j) * [d^j (ωsc)^k] * [d^{M-j} G^{(n)}]
    # 注意：只在 j <= k 时 d^j (ωsc)^k 才不为零
    for j in 0:min(M, k)

        # --- Part A: (ωsc)^k 关于 μ 的 j 阶导数 ---
        # 导数公式: k * (k-1) * ... * (k-j+1) * (ωsc)^(k-j) * (-1/λ)^j
        # 我们把 (-1/λ)^j 拆解为：(-1)^j * (1/λ)^j

        deriv_A_coeff = 1.0
        for i in 0:(j-1)
            deriv_A_coeff *= (k - i)
        end
        term_A = deriv_A_coeff * (ωsc)^(k - j)

        # --- Part B: G^{(n)} 关于 μ 的 P = M-j 阶导数 ---
        # 使用您 muderiv! 中的逻辑：
        # d^P G / dμ^P ~ acc * λ^(P-1) * (-ωsc)^P * P!
        # *关键修改*：为了配合主函数最后的 /λ，这里我们将 scaling 提升为 λ^P (即乘了一个 λ)
        # 这样 Part B 的 scaling 为 λ^P。

        P = M - j
        fact_P = factorial(P)

        # 计算 G 的展开 (muderiv 核心循环)
        acc_B = 0.0
        if P == 0
            acc_B = propagator_derivative(τval, ωsc, β, n)
        else
            term_inner = 1.0
            for p in 1:P
                val = propagator_derivative(τval, ωsc, β, n + p)
                term_inner *= ωsc
                acc_B += val * term_inner * binomial(P - 1, p - 1)
            end
            # 应用 muderiv 的系数，但注意 scaling调整：
            # 原公式: acc * λ^(P-1) * (-ωsc)^P * P!
            # 这里的 scaling_B 我们拆分处理
            acc_B *= (-ωsc)^P * fact_P
        end

        # --- 组合 Part A 和 Part B ---
        # 总项 = C(M,j) * [Part A] * [Part B]
        # 缩放因子分析：
        # Part A 带来了 (1/λ)^j
        # Part B 我们赋予 λ^P 的权重 (P = M-j)
        # 总 λ 因子: λ^P * (1/λ)^j = λ^{M-j} * λ^{-j} = λ^{M-2j}

        lambda_factor = λ^(M - 2 * j)

        # 符号 (-1)^j 来自 d(ωsc)/dμ
        sign_factor = (-1)^j

        res += binomial(M, j) * (term_A * sign_factor) * acc_B * lambda_factor
    end

    return res
end

@inline function propagator_mixderiv!(Gk_buffer, Lkx, Lky, M, L, τval, λeff_arr, ωscaled_arr, β)
    # M: Order of mu derivative
    # L: Order of lambda derivative

    @inbounds for iy in 1:Lky, ix in 1:Lkx
        λ = λeff_arr[ix, iy]
        ωsc = ωscaled_arr[ix, iy]

        total_acc = 0.0

        # sum C(L, o) * (ωsc)^o * G^{(o)} 

        for o in 0:L
            term_val = calc_mu_deriv_component(M, o, o, τval, ωsc, β, λ)
            total_acc += binomial(L, o) * term_val
        end
        Gk_buffer[ix, iy] = total_acc / λ
    end
end

"""
    build_pretab(...)

Builds the lookup table for an arbitrary τgrid. 
Automatically detects if τgrid is uniform to enable O(1) lookups.
"""
function build_pretab(t, β, τgrid::Vector{Float64}, propagator!::Function=propagator_lamderiv!;
    lambda, dμ, deriv_order, Lkx, Lky, Rtable)

    # println("Info: Starting memory-efficient pre-tabulation...")

    Ntau = length(τgrid)
    Norb = 2
    Mplus1 = deriv_order + 1
    # Sum_{x=0 to R} (x+1) = (R+1)(R+2)/2
    Nr = ((Rtable + 1) * (Rtable + 2)) ÷ 2

    mem_gb = Ntau * Norb * Mplus1 * Nr * 8 / 1024^3
    # println("Info: Rtable=$Rtable. Storing $Nr spatial points (vs original $(2*Rtable+1)^2).")
    # println("Info: Final Ctable Size: $(round(mem_gb, digits=2)) GB")

    # Check uniformity with a small tolerance
    dτ = τgrid[2] - τgrid[1]
    is_uniform = all(isapprox.(diff(τgrid), dτ; atol=1e-9))

    Ctable = zeros(Float64, Ntau, Norb, Mplus1, Nr)

    ϵk, _, _ = dispersion_kmesh(Lkx, Lky, t, dμ)
    Gk_buffer = zeros(ComplexF64, Lkx, Lky)
    p_ifft = plan_ifft(Gk_buffer)

    Gk_buffer1 = zeros(ComplexF64, Lkx, Lky)

    # println("Info: Computing & FFT per tau point...")

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
        # if itau % 1000 == 0
        #     print("\rProgress: $itau / $Ntau")
        # end

        for m1 in 1:Mplus1 # m1 = order + 1
            m = m1 - 1


            # Assume orb=1,1 and orb=2,2 symmetric
            # propagator!(Gk_buffer, Lkx, Lky, m, τval, λeff_arr, ωscaled_arr, β)
            propagator_mixderiv!(Gk_buffer, Lkx, Lky, m, 0, τval, λeff_arr, ωscaled_arr, β)
            # propagator_mixderiv!(Gk_buffer1, Lkx, Lky, 0, m, τval, λeff_arr, ωscaled_arr, β)


            # idx = findfirst(.!isapprox.(Gk_buffer, Gk_buffer1))
            # if idx !== nothing
            #     println(idx, " ", Gk_buffer[idx])
            #     println(Gk_buffer1[idx])
            #     @warn "Propagator mismatch for m=$m, tau=$τval"
            # end

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
    # println("\nInfo: Pre-tabulation finished.")

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