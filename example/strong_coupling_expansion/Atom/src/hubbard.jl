module Hubbard
include("common.jl")
# include("hilbert.jl")
using ..Hilbert
# include("green.jl")
using ..Green
using LinearAlgebra
# import QuantumStatistics.Basis:tau2matfreq, tau2dlr

function fermiHubbard(t, U, μ, sites, bonds)
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

    # show(stdout, "text/plain", Matrix(K))
    # println()

    H = -t * K + U * V - μ * C

    # show(stdout, "text/plain", Matrix(H))
    # println()

    return H, cu⁺, cd⁺
end

function hubbardAtom(type, U, μ, β)
    E = [0.0, -μ, -μ, U - 2μ]
    H = zeros(Float, (4, 4))
    H[diagind(H)] = E

    # |0>=|00>=1, |↑>=|10>=2, |↓>=|01>=3, |↑↓>=|11>=4
    # the first is the forck state for ↑ spin, the second is forck state for ↓

    cpup = zeros(Float, (4, 4)) 
    cpdown = zeros(Float, (4, 4)) 

    cpup[2, 1], cpup[4, 3] = 1, 1
    cpdown[3, 1], cpdown[4, 2] = 1, 1
    cmup, cmdown = cpup', cpdown'
    
    @assert abs(tr(cpup * cmdown)) < 1e-16
    @assert abs(tr(cpdown * cmup)) < 1e-16

    m = Model(β, H, [cpup, cpdown])
    println("U=$U, μ=$μ, β=$β")
    println("Model Hilbert space: $(m.dim)")
    println("The lowest eigen energy per site:\n$(E)")
    return m
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

    m = Model(β, H, hcat(cu⁺, cd⁺))

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

    m = Model(β, H, hcat(cu⁺, cd⁺))

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