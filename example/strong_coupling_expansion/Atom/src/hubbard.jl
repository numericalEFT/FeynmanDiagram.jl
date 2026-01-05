module Hubbard
include("common.jl")
# include("hilbert.jl")
using ..Hilbert
# include("green.jl")
using ..Green
using LinearAlgebra
# import QuantumStatistics.Basis:tau2matfreq, tau2dlr

function fermiHubbard(t, U, μ, sites, bonds, h=0.0)
    Nsite = length(sites)
    Fock = Hilbert.BinaryFock(Nsite)
    cu⁺ = [Hilbert.creation(Fock, s, UP) for s in sites]
    cd⁺ = [Hilbert.creation(Fock, s, DOWN) for s in sites]
    cu⁻ = [c' for c in cu⁺]
    cd⁻ = [c' for c in cd⁺]
    nu = [cu⁺[s] * cu⁻[s] for s in sites]
    nd = [cd⁺[s] * cd⁻[s] for s in sites]

    K = sum([cu⁺[i] * cu⁻[j] + cd⁺[i] * cd⁻[j] for (i, j) in bonds])
    V = sum(nu[s] * nd[s] for s in sites)
    C = sum(nu[s] + nd[s] for s in sites)
    M = sum(nu[s] - nd[s] for s in sites)

    # show(stdout, "text/plain", Matrix(K))
    # println()

    H = -t * K + U * V - μ * C - h * M

    # show(stdout, "text/plain", Matrix(H))
    # println()

    return H, cu⁺, cd⁺
end

function hubbardAtom(type, U, μ, β, h=0.0)
    E = [0.0, -μ - h, -μ + h, U - 2μ]
    H = zeros(Float64, (4, 4))
    H[diagind(H)] = E

    # |0>=|00>=1, |↑>=|10>=2, |↓>=|01>=3, |↑↓>=|11>=4
    # the first is the forck state for ↑ spin, the second is forck state for ↓

    cpup = zeros(Float64, (4, 4))
    cpdown = zeros(Float64, (4, 4))

    cpup[2, 1], cpup[4, 3] = 1, 1
    cpdown[3, 1], cpdown[4, 2] = 1, -1
    cmup, cmdown = cpup', cpdown'

    @assert abs(tr(cpup * cmdown)) < 1e-16
    @assert abs(tr(cpdown * cmup)) < 1e-16

    m = Model(β, H, [cpup, cpdown])
    println("U=$U, μ=$μ, β=$β")
    println("Model Hilbert space: $(m.dim)")
    println("The lowest eigen energy per site:\n$(E)")
    return m
end

"""
    get_mu_from_n(target_n::Float64, beta, U, h=0.0)

Calculate chemical potential mu given density n for Hubbard atom model.
Uses a numerically stable quadratic solver to handle large beta/U.
"""
function get_mu_from_n(target_n::Float64, beta, U, h=0.0)
    if target_n <= 1e-14
        return -Inf
    elseif target_n >= 2.0 - 1e-14
        return Inf
    end

    # 2. Coefficients for quadratic equation: Ax^2 + Bx + C = 0
    # where x = exp(beta * mu)

    # Handle potential underflow for very large beta * U
    # exp(-beta * U) might become 0.0, which is fine, but we must handle it gracefully.
    factor_double = exp(-beta * U)
    factor_single = 2.0 * cosh(beta * h)

    a = (2.0 - target_n) * factor_double
    b = (1.0 - target_n) * factor_single
    c = -target_n

    discriminant = b^2 - 4 * a * c

    if discriminant < 0
        error("Discriminant < 0. This should not happen for physical n.")
    end

    sqrt_disc = sqrt(discriminant)
    x = 0.0

    # 3. Stable Quadratic Solver (Citardauq Formula)
    # To avoid cancellation errors when b is large and close to sqrt(b^2 - 4ac)

    if b > 0
        # Case: n < 1. b is positive and large.
        # Standard formula (-b + sqrt(...)) involves cancellation (large - large).
        # Use rationalized form: x = -2c / (b + sqrt(...))
        x = -2 * c / (b + sqrt_disc)
    else
        # Case: n >= 1. b is negative or zero.
        # -b is positive. No cancellation in (-b + sqrt(...)).
        numerator = -b + sqrt_disc
        denominator = 2 * a

        # Check for division by zero if A is extremely small (U very large)
        if abs(denominator) < 1e-100
            # Fallback to linear solution Bx + C = 0 -> x = -C/B
            x = -c / b
        else
            x = numerator / denominator
        end
    end

    if x <= 0
        error("Solver returned non-positive x: $x. Parameters might be too extreme.")
    end

    return log(x) / beta
end

function hubbardAtom2(type, t, U, μ, β)
    Nsite = 2
    bonds = [(1, 2), (2, 1)]
    sites = [s for s in 1:Nsite]

    if type == :fermi
        H, cu⁺, cd⁺ = fermiHubbard(t, U, μ, sites, bonds)
    else
        @error("Not implemented!")
    end

    m = Model(β, Matrix(H), vcat(Matrix.(cu⁺), Matrix.(cd⁺)))

    println("t=$t, U=$U, μ=$μ, β=$β")
    println("Model Hilbert space: $(m.dim)")
    println("The lowest eigen energy per site:\n$(m.E[1:8] / Nsite)")
    # E, U = eigen(Matrix(H))
    # println(E)
    return m
end

function hubbardAtom4(type, t, U, μ, β)
    Nsite = 4
    bonds = [(1, 2), (2, 1), (2, 3), (3, 2), (3, 4), (4, 3), (1, 4), (4, 1)]
    sites = [s for s in 1:Nsite]

    if type == :fermi
        H, cu⁺, cd⁺ = fermiHubbard(t, U, μ, sites, bonds)
    else
        @error("Not implemented!")
    end

    m = Model(β, Matrix(H), vcat(Matrix.(cu⁺), Matrix.(cd⁺)))

    println("t=$t, U=$U, μ=$μ, β=$β")
    println("Model Hilbert space: $(m.dim)")
    println("The lowest eigen energy per site:\n$(m.E[1:16] / Nsite)")
    # E, U = eigen(Matrix(H))
    # println(m.E[1])
    return m
end

# const m = Model(10.0, 5.0, β=10.0)
# propagator(m)
end