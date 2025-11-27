module Green

include("common.jl")
using LinearAlgebra, Combinatorics

export Model, GreenN
export density, Heisenberg, thermalavg, parity
export Gn, G2, G_with_D, dGn_dU_estimator

# -----------------------------------------------------------------------------
# Model: local (atomic) problem in its eigenbasis
# -----------------------------------------------------------------------------
# Notes vs your original version (fileciteturn2file0):
#  1. We fix the eigen-sorting bug: we sort eigenvalues AND reorder eigenvectors.
#  2. We keep orbital support: c⁺[orb], c⁻[orb], n[orb].
#  3. We add D, the local double-occupancy operator n_up * n_dn in eigenbasis.
#     By default we take doublon_orbitals = (1,2), i.e. orbital 1 = ↑, 2 = ↓.
#  4. Field names/exports stay compatible with your code.
# -----------------------------------------------------------------------------

struct Model{N,No}
    isfermi::Bool
    β::Float64
    dim::Int           # Hilbert-space dimension
    Norbital::Int      # number of orbitals/spin flavors
    E::SVector{N,Float64}    # eigen-energies (sorted ascending)
    Z::Float64               # physical partition sum = Tr[e^{-β H_loc}]
    lnZ::Float64             # log(Z) for stable algebra
    w::SVector{N,Float64}    # normalized Boltzmann weights e^{-β Ê} with Σw=1
    Enorm::SVector{N,Float64}  # shifted energies E + lnZ/β
    Hdiag::Matrix{Float64}        # diagonal H in eigenbasis (physical)
    Hnorm::Matrix{Float64}        # shifted H with partition 1
    ΔE::Matrix{Float64}      # pairwise energy differences E[a]-E[b]
    c⁺::SVector{No,Matrix{Float64}}  # c† for each orbital, in eigenbasis
    c⁻::SVector{No,Matrix{Float64}}  # c  for each orbital, in eigenbasis
    n::SVector{No,Matrix{Float64}}   # n = c† c for each orbital, in eigenbasis
    Ntot::Matrix{Float64}            # total density operator ∑ n_orb
    M::Matrix{Float64}               # local magnetization operator in eigenbasis
    D::Matrix{Float64}               # local double-occupancy operator in eigenbasis

    function Model(β, H, c⁺_fock::Vector{Matrix{Float64}}, isfermi::Bool=true;
        doublon_orbitals::Tuple{Int,Int}=(UP, DOWN))
        dim = size(H, 1)
        @assert size(H) == size(c⁺_fock[1])
        @assert size(H) == (dim, dim)
        Norbital = length(c⁺_fock)

        # Diagonalize local Hamiltonian H (Fock basis -> eigenbasis)
        F = eigen(Float64.(Matrix(H)))
        E_unsorted = F.values
        U_unsorted = F.vectors
        p = sortperm(E_unsorted)
        E_sorted = E_unsorted[p]
        U_sorted = U_unsorted[:, p]

        # Partition function Z = Σ_a e^{-β E_a}, computed stably
        logw = -β .* E_sorted
        logw_max = maximum(logw)
        lnZ = logw_max + log(sum(exp.(logw .- logw_max)))
        Z = exp(lnZ)
        w = exp.(logw .- lnZ)  # normalized weights, sum(w)=1

        # Diagonal H in eigenbasis (physical and normalized copies)
        Hdiag = zeros(Float64, (dim, dim))
        Hnorm = zeros(Float64, (dim, dim))
        shift = lnZ / β
        Enorm = E_sorted .+ shift
        Hdiag[diagind(Hdiag)] = E_sorted
        Hnorm[diagind(Hnorm)] = Enorm
        ΔE = Matrix{Float64}(undef, dim, dim)
        @inbounds for a in 1:dim, b in 1:dim
            ΔE[a, b] = E_sorted[a] - E_sorted[b]
        end

        # Rotate creation operators into eigenbasis
        c⁺ = [U_sorted' * o * U_sorted for o in c⁺_fock]
        # Define annihilation operators as Hermitian adjoint of c† in that basis
        c⁻ = [adjoint(op) for op in c⁺]
        # Densities n_orb = c†_orb c_orb in eigenbasis
        n_ops = [c⁺[i] * c⁻[i] for i in 1:Norbital]
        Ntot = reduce(+, n_ops)

        # Local double-occupancy operator D = n_up * n_dn
        up, dn = doublon_orbitals
        Mop = n_ops[up] - n_ops[dn]
        Dop = n_ops[up] * n_ops[dn]

        return new{dim,Norbital}(isfermi,
            β,
            dim,
            Norbital,
            E_sorted,
            Z,
            lnZ,
            w,
            Enorm,
            Hdiag,
            Hnorm,
            ΔE,
            c⁺,
            c⁻,
            n_ops,
            Ntot,
            Mop,
            Dop)
    end
end

# -----------------------------------------------------------------------------
# Finite-T trace utilities in eigenbasis
# -----------------------------------------------------------------------------

"""
thermalavg(O::Operator, E, β, Z)

Return ⟨O⟩ = Tr[ O e^{-β H} ] / Z, assuming O is expressed in the eigenbasis of H
and E are the eigenvalues of H (so e^{-βH} is diagonal with entries e^{-βE[a]}).
Matches your original API. (fileciteturn2file0)
"""
function thermalavg(O::Matrix{Float64}, w, Z)
    @assert size(O, 1) == length(w)
    return dot(diag(O), w)
end

"""
Heisenberg(O::Operator, spectrum, τ)

Heisenberg-evolve a local operator O with the *atomic* Hamiltonian:
    O(τ) = e^{τ H} O e^{-τ H}
`spectrum` can be the eigenvalues `E` or a precomputed matrix ΔE with entries E[a]-E[b].
Matches (and fixes nothing) from your code. (fileciteturn2file0)
"""
function Heisenberg(O::Matrix{Float64}, spectrum, τ)
    n = size(O, 1)
    @assert size(O, 2) == n
    if abs(τ) < 1e-10
        return O
    end

    use_diff_matrix = spectrum isa AbstractMatrix
    if use_diff_matrix
        @assert size(spectrum) == (n, n)
    else
        @assert length(spectrum) == n
    end

    Out = similar(O)
    if use_diff_matrix
        ΔE = spectrum
        @inbounds for a in 1:n, b in 1:n
            Out[a, b] = O[a, b] * exp(ΔE[a, b] * τ)
        end
    else
        Evals = spectrum
        @inbounds for a in 1:n, b in 1:n
            Out[a, b] = O[a, b] * exp((Evals[a] - Evals[b]) * τ)
        end
    end
    return Out
end

const τ_tol = 1e-12

@inline function propagator(m::Model, Δτ::Real)
    if abs(Δτ) < τ_tol
        return nothing
    else
        return @. exp(-Δτ * m.Enorm)
    end
end

@inline function copy_scaled!(dest::Matrix{Float64}, src::Matrix{Float64}, factors)
    if factors === nothing
        copyto!(dest, src)
    else
        for j in axes(src, 2)
            @inbounds @simd for i in axes(src, 1)
                dest[i, j] = factors[i] * src[i, j]
            end
        end
    end
    return dest
end

function chain_trace(m::Model, τ_ord::Vector{Float64}, ops_ord::Vector{Matrix{Float64}})
    nops = length(ops_ord)
    @assert nops == length(τ_ord) && nops > 0
    dim = m.dim
    current = Matrix{Float64}(undef, dim, dim)
    scaled = Matrix{Float64}(undef, dim, dim)
    work = Matrix{Float64}(undef, dim, dim)

    scale = propagator(m, m.β + τ_ord[end] - τ_ord[1])
    copy_scaled!(current, ops_ord[1], scale)

    for i in 2:nops
        scale = propagator(m, τ_ord[i-1] - τ_ord[i])
        copy_scaled!(scaled, ops_ord[i], scale)
        mul!(work, current, scaled)
        current, work = work, current
    end
    return tr(current)
end

@inline function ordered_chain_trace(m::Model, τ::Vector{Float64}, ops::Vector{Matrix{Float64}}, perm::Vector{Int})
    return chain_trace(m, τ[perm], ops[perm])
end

function fermionic_sign(m::Model, perm::Vector{Int}, NF::Int)
    if !m.isfermi
        return 1
    end
    fermion_positions = Vector{Int}(undef, NF)
    for newpos in 1:length(perm)
        oldidx = perm[newpos]
        if oldidx <= NF
            fermion_positions[oldidx] = newpos
        end
    end
    return parity(sortperm(fermion_positions))
end

@inline thermal_expectation(m::Model, O::Matrix{Float64}) = chain_trace(m, Float64[0.0], Matrix{Float64}[O])


"""
parity(p)

Fermionic sign of a permutation p (same as your implementation).
"""
function parity(p)
    n = length(p)
    not_seen = Set{Int}(1:n)
    seen = Set{Int}()
    cycles = Array{Int,1}[]
    while !isempty(not_seen)
        cycle = Int[]
        x = pop!(not_seen)
        while !in(x, seen)
            push!(cycle, x)
            push!(seen, x)
            x = p[x]
            pop!(not_seen, x, 0)
        end
        push!(cycles, cycle)
    end
    cycle_lengths = map(length, cycles)
    even_cycles = filter(i -> i % 2 == 0, cycle_lengths)
    return (length(even_cycles) % 2 == 0) ? 1 : -1
end

# -----------------------------------------------------------------------------
# GreenN: N-body (2N-leg) atomic Green's function container
# -----------------------------------------------------------------------------
"""
GreenN(m::Model, τ, orbital, isfermi=true)

Construct struct to store the variables to evaulate N-body Green's functions.
The leg index is assumed to be [1, 2, 3, ...,2N], where the incoming legs are 1:N, and the outgoing legs are N+1:2N
The full Green's function is defined as,
```math
Gn = <Tτ c⁺(1)c⁺(2)...c⁺(N)c(N+1)...c(2N)>
```
e.g.,
1->------->-3
	| G4 |        
2->------->-4
All other Green's function are derived from the above full Green's function
"""
struct GreenN
    N::Int                  # n-body or 2n-point
    τ::Vector{Float64}        # times for each leg (length 2N)
    orbital::Vector{Int}    # orbital/spin index for each leg (length 2N)
    ops::Vector{Matrix{Float64}}   # bare operators for each leg (creation first)

    function GreenN(m::Model, τ, orbital)
        N = length(τ) ÷ 2
        @assert length(τ) == length(orbital) == 2N "Length of τ and orbital must be 2N."

        ops = Vector{Matrix{Float64}}(undef, 2N)
        # incoming: creation legs 1..N
        for i in 1:N
            ops[i] = m.c⁺[orbital[i]]
        end
        # outgoing: annihilation legs N+1..2N
        for i in (N+1):2N
            ops[i] = m.c⁻[orbital[i]]
        end

        return new(N, τ, orbital, ops)
    end
end

# density for a given orbital (unchanged API)
density(m::Model, orbital) = thermal_expectation(m, m.n[orbital])

# -----------------------------------------------------------------------------
# Bare Gn and G2 (unchanged semantics)
# -----------------------------------------------------------------------------

"""
Gn(m,g)
Return the full 2N-leg atomic Green's function
  G^N = ⟨ Tτ c†(1) c†(2) ... c†(N) c(N+1) ... c(2N) ⟩.
e.g.,
1->------->-3
	| G4 |        
2->------->-4
"""
function Gn(m::Model, g::GreenN)
    perm = sortperm(g.τ; rev=true)
    sign = fermionic_sign(m, perm, length(g.τ))
    Gval = ordered_chain_trace(m, g.τ, g.ops, perm)
    return sign * Gval
end

# Your explicit G2 left untouched (for debugging / benchmarking analytic form)
function G2(m::Model, g::GreenN)
    if g.N != 1
        throw(AssertionError("G2 function is for n=2"))
    end
    τi, τo = g.τ[1], g.τ[2]
    β, E = m.β, m.E
    c⁺, c⁻ = m.c⁺[g.orbital[1]], m.c⁻[g.orbital[2]]

    G = 0.0
    if (τi < τo)
        for j in 1:m.dim
            for k in 1:m.dim
                G += exp(-(β - τo + τi) * E[j] - (τo - τi) * E[k]) * c⁻[j, k] * c⁺[k, j]
            end
        end
    else
        for j in 1:m.dim
            for k in 1:m.dim
                G += -exp(-(β - τi + τo) * E[j] - (τi - τo) * E[k]) * c⁺[j, k] * c⁻[k, j]
            end
        end
    end
    return G / m.Z
end

# -----------------------------------------------------------------------------
# G_with_D: insert local doublon operator D(τp) = n↑n↓(τp)
# -----------------------------------------------------------------------------
"""
G_with_D(m, g, τp)

Return ⟨ Tτ[ (fermion legs in g) · D(τp) ] ⟩, where D = n↑ n↓.
We:
  1. Build D(τp) in Heisenberg picture.
  2. Merge it with the fermionic legs.
  3. Time-order all operators together.
  4. Compute fermionic sign from *only* the fermionic legs. D is even-parity,
     so moving it through does not add extra minus signs.
"""
function G_with_D(m::Model, g::GreenN, τp::Real)
    τ_ext = [g.τ; τp]
    op_ext = vcat(g.ops, [m.D])

    NF = length(g.τ)
    perm_all = sortperm(τ_ext; rev=true)
    fermion_sign = fermionic_sign(m, perm_all, NF)

    GwD = ordered_chain_trace(m, τ_ext, op_ext, perm_all)
    return fermion_sign * GwD
end

# -----------------------------------------------------------------------------
# dGn_dU_estimator: unbiased stochastic estimator for ∂G/∂U
# -----------------------------------------------------------------------------

"""
dGn_dU_estimator(m, g, τp)

One-shot estimator for ∂G/∂U using explicit operator insertion, *no τ-grid*.
Theory:
  ∂G/∂U = - ∫₀^β dτ' [ ⟨Tτ(legs · D(τ'))⟩ - G · ⟨D⟩ ].
If τp ~ Uniform(0,β), then
  E_{τp}[ β ( GwD(τp) - G Dloc ) ] = ∫₀^β dτ' [...],
so
  ∂G/∂U = E_{τp}[ -β ( GwD(τp) - G Dloc ) ].

Call this with a random τp ∈ [0,β) inside your MC measurement loop and
average its return value. That average → ∂G/∂U.
"""
function dGn_dU_estimator(m::Model, g::GreenN, τp::Real)
    Gval = Gn(m, g)                        # G = ⟨Tτ legs⟩
    Dloc = thermal_expectation(m, m.D)     # ⟨D⟩ = ⟨n↑n↓⟩_atom
    GwD = G_with_D(m, g, τp)              # ⟨Tτ(legs · D(τp))⟩

    # return -m.β * (GwD - Gval * Dloc)
    return -GwD + Gval * Dloc
end

function G_with_N(m::Model, g::GreenN, τp::Real)
    τ_ext = [g.τ; τp]
    op_ext = vcat(g.ops, [m.Ntot])

    NF = length(g.τ)
    perm_all = sortperm(τ_ext; rev=true)
    fermion_sign = fermionic_sign(m, perm_all, NF)

    GwN = ordered_chain_trace(m, τ_ext, op_ext, perm_all)
    return fermion_sign * GwN
end

function dGn_dμ_estimator(m::Model, g::GreenN, τp::Real)
    Gval = Gn(m, g)  # ⟨Tτ legs⟩
    # local density expectation ⟨N⟩
    Nloc = thermal_expectation(m, m.Ntot)
    # correlator with insertion
    GwN = G_with_N(m, g, τp)
    return GwN - Gval * Nloc
end


end # module Green
