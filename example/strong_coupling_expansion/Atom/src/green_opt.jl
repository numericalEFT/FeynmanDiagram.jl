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
    β::Float
    dim::Int           # Hilbert-space dimension
    Norbital::Int      # number of orbitals/spin flavors
    E::SVector{N,Float}  # eigen-energies (sorted ascending)
    Z::Float           # partition sum = Tr[e^{-β H_loc}]
    w::SVector{N,Float}    # Boltzmann weights e^{-β E}
    Hdiag::Operator    # diagonal H in eigenbasis
    c⁺::SVector{No,Operator}  # c† for each orbital, in eigenbasis
    c⁻::SVector{No,Operator}  # c  for each orbital, in eigenbasis
    n::SVector{No,Operator}   # n = c† c for each orbital, in eigenbasis
    M::Operator               # local magnetization operator in eigenbasis
    D::Operator               # local double-occupancy operator in eigenbasis

    function Model(β, H, c⁺_fock::Vector{Operator}, isfermi::Bool=true;
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

        # Partition function Z = Σ_a e^{-β E_a}
        w = exp.(-β .* E_sorted)
        Z = sum(w)

        # Diagonal H in eigenbasis
        Hdiag = zeros(Float, (dim, dim))
        Hdiag[diagind(Hdiag)] = E_sorted

        # Rotate creation operators into eigenbasis
        c⁺ = [U_sorted' * o * U_sorted for o in c⁺_fock]
        # Define annihilation operators as Hermitian adjoint of c† in that basis
        c⁻ = [adjoint(op) for op in c⁺]
        # Densities n_orb = c†_orb c_orb in eigenbasis
        n_ops = [c⁺[i] * c⁻[i] for i in 1:Norbital]

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
            Hdiag,
            c⁺,
            c⁻,
            n_ops,
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
function thermalavg(O::Operator, w, Z)
    @assert size(O, 1) == length(w)
    return dot(diag(O), w) / Z
end

"""
Heisenberg(O::Operator, E, τ)

Heisenberg-evolve a local operator O with the *atomic* Hamiltonian:
    O(τ) = e^{τ H} O e^{-τ H}
In eigenbasis this is just elementwise multiplication by exp(+τE[a]) and exp(-τE[b]).
Matches (and fixes nothing) from your code. (fileciteturn2file0)
"""
function Heisenberg(O::Operator, E, τ)
    n = length(E)
    @assert size(O) == (n, n)
    if abs(τ) < 1e-10
        return O
    end
    # allocate result same size/type as O
    Out = similar(O)
    @inbounds for a in 1:n, b in 1:n
        Out[a, b] = O[a, b] * exp((E[a] - E[b]) * τ)
    end
    return Out
end

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
    τ::Vector{Float}        # times for each leg (length 2N)
    orbital::Vector{Int}    # orbital/spin index for each leg (length 2N)
    hop::Vector{Operator}   # Heisenberg-evolved operators for each leg

    function GreenN(m::Model, τ, orbital)
        N = length(τ) ÷ 2
        @assert length(τ) == length(orbital) == 2N "Length of τ and orbital must be 2N."

        hop = Vector{Operator}(undef, 2N)
        # incoming: creation legs 1..N
        for i in 1:N
            hop[i] = Heisenberg(m.c⁺[orbital[i]], m.E, τ[i])
        end
        # outgoing: annihilation legs N+1..2N
        for i in (N+1):2N
            hop[i] = Heisenberg(m.c⁻[orbital[i]], m.E, τ[i])
        end

        return new(N, τ, orbital, hop)
    end
end

# density for a given orbital (unchanged API)
density(m::Model, orbital) = thermalavg(m.n[orbital], m.w, m.Z)

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
    τ = g.τ
    hop = g.hop

    perm = sortperm(τ)
    ordered_hop = hop[perm]

    # Time-ordered product: multiply from latest τ to earliest τ
    M = ordered_hop[end]
    for op in Iterators.reverse(ordered_hop[1:(end-1)])
        M *= op
    end

    Gval = thermalavg(M, m.w, m.Z)
    if m.isfermi
        return Gval * parity(perm)
    else
        return Gval
    end
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
    # Heisenberg-evolve D at τp
    Dτ = Heisenberg(m.D, m.E, τp)

    # Extended lists of times/operators
    τ_ext = [g.τ; τp]
    # op_ext = [g.hop; Dτ]
    op_ext = vcat(g.hop, [Dτ])

    NF = length(g.τ)          # number of fermionic legs = 2N
    perm_all = sortperm(τ_ext)

    # Figure out how the original fermionic legs were permuted inside perm_all
    fermion_positions = Vector{Int}(undef, NF)
    for newpos in eachindex(perm_all)
        oldidx = perm_all[newpos]
        if oldidx <= NF
            fermion_positions[oldidx] = newpos
        end
    end
    fermion_sign = m.isfermi ? parity(sortperm(fermion_positions)) : 1

    # Build fully time-ordered product (latest τ last in 'perm_all')
    op_ord = op_ext[perm_all]
    M = op_ord[end]
    for op in Iterators.reverse(op_ord[1:(end-1)])
        M *= op
    end

    GwD = thermalavg(M, m.w, m.Z)
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
    Dloc = thermalavg(m.D, m.w, m.Z)  # ⟨D⟩ = ⟨n↑n↓⟩_atom
    GwD = G_with_D(m, g, τp)              # ⟨Tτ(legs · D(τp))⟩

    # return -m.β * (GwD - Gval * Dloc)
    return -GwD + Gval * Dloc
end

function G_with_N(m::Model, g::GreenN, τp::Real)
    # Heisenberg-evolve D at τp
    Nτ = Heisenberg(sum(m.n), m.E, τp)

    # Extended lists of times/operators
    τ_ext = [g.τ; τp]
    op_ext = vcat(g.hop, [Nτ])

    NF = length(g.τ)          # number of fermionic legs = 2N
    perm_all = sortperm(τ_ext)

    # Figure out how the original fermionic legs were permuted inside perm_all
    fermion_positions = Vector{Int}(undef, NF)
    for newpos in eachindex(perm_all)
        oldidx = perm_all[newpos]
        if oldidx <= NF
            fermion_positions[oldidx] = newpos
        end
    end
    fermion_sign = m.isfermi ? parity(sortperm(fermion_positions)) : 1

    # Build fully time-ordered product (latest τ last in 'perm_all')
    op_ord = op_ext[perm_all]
    M = op_ord[end]
    for op in Iterators.reverse(op_ord[1:(end-1)])
        M *= op
    end

    GwN = thermalavg(M, m.w, m.Z)
    return fermion_sign * GwN
end

function dGn_dμ_estimator(m::Model, g::GreenN, τp::Real)
    Gval = Gn(m, g)  # ⟨Tτ legs⟩
    # total density operator N = sum over orbitals
    Nmat = zero(m.n[1])
    for orb in 1:m.Norbital
        Nmat .+= m.n[orb]
    end
    # local density expectation ⟨N⟩
    Nloc = thermalavg(Nmat, m.w, m.Z)
    # correlator with insertion
    GwN = G_with_N(m, g, τp)  # you'll implement same as GwD but using Nmat
    return GwN - Gval * Nloc
end


end # module Green
