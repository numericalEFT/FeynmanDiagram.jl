module Green

include("common.jl")
using LinearAlgebra, Combinatorics, StaticArrays

export Model, GreenN, GreenWorkspace
export density, Heisenberg, thermalavg, parity
export Gn, G2, dGn_dU_estimator

# -----------------------------------------------------------------------------
# Model: local (atomic) problem in its eigenbasis
# -----------------------------------------------------------------------------
struct Model{N,No}
    isfermi::Bool
    β::Float64
    dim::Int           # Hilbert-space dimension
    Norbital::Int      # number of orbitals/spin flavors
    E::SVector{N,Float64}    # eigen-energies (sorted ascending)
    Z::Float64               # physical partition sum
    lnZ::Float64             # log(Z)
    w::SVector{N,Float64}    # normalized Boltzmann weights
    Enorm::SVector{N,Float64}  # shifted energies E + lnZ/β

    # Operators in eigenbasis
    # Note: Even for N<=16, we keep Matrix here to avoid compilation overhead and flexibility, but we use non-allocating operations in traces.
    Hdiag::Matrix{Float64}
    Hnorm::Matrix{Float64}
    ΔE::Matrix{Float64}
    c⁺::SVector{No,Matrix{Float64}}
    c⁻::SVector{No,Matrix{Float64}}
    n::SVector{No,Matrix{Float64}}
    Ntot::Matrix{Float64}
    M::Matrix{Float64}
    D::Matrix{Float64}

    function Model(β, H, c⁺_fock::Vector{Matrix{Float64}}, isfermi::Bool=true;
        doublon_orbitals::Tuple{Int,Int}=(1, 2))
        dim = size(H, 1)
        @assert size(H) == size(c⁺_fock[1])
        Norbital = length(c⁺_fock)

        F = eigen(Float64.(Matrix(H)))
        E_unsorted = F.values
        U_unsorted = F.vectors
        p = sortperm(E_unsorted)
        E_sorted = E_unsorted[p]
        U_sorted = U_unsorted[:, p]

        logw = -β .* E_sorted
        logw_max = maximum(logw)
        lnZ = logw_max + log(sum(exp.(logw .- logw_max)))
        Z = exp(lnZ)
        w = exp.(logw .- lnZ)

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

        c⁺ = [U_sorted' * o * U_sorted for o in c⁺_fock]
        c⁻ = [adjoint(op) for op in c⁺]
        n_ops = [c⁺[i] * c⁻[i] for i in 1:Norbital]
        Ntot = reduce(+, n_ops)

        up, dn = doublon_orbitals
        Mop = n_ops[up] - n_ops[dn]
        Dop = n_ops[up] * n_ops[dn]

        return new{dim,Norbital}(isfermi, β, dim, Norbital,
            SVector{dim}(E_sorted), Z, lnZ, SVector{dim}(w), SVector{dim}(Enorm),
            Hdiag, Hnorm, ΔE,
            SVector{Norbital}(c⁺), SVector{Norbital}(c⁻), SVector{Norbital}(n_ops),
            Ntot, Mop, Dop)
    end
end

# -----------------------------------------------------------------------------
# GreenWorkspace: Memory pool for Zero-Allocation calculations
# -----------------------------------------------------------------------------
struct GreenWorkspace
    # Matrix buffers (Ping-Pong buffers for trace)
    mat_A::Matrix{Float64}
    mat_B::Matrix{Float64}
    mat_C::Matrix{Float64}

    # Data buffers
    τ_buffer::Vector{Float64}
    ops_buffer::Vector{Matrix{Float64}}
    perm_buffer::Vector{Int}

    function GreenWorkspace(model::Model, max_order::Int)
        dim = model.dim
        # max_len needs to hold standard legs (2*order) plus insertion (1)
        max_len = 2 * max_order + 2

        new(
            Matrix{Float64}(undef, dim, dim),
            Matrix{Float64}(undef, dim, dim),
            Matrix{Float64}(undef, dim, dim),
            Vector{Float64}(undef, max_len),
            Vector{Matrix{Float64}}(undef, max_len),
            Vector{Int}(undef, max_len)
        )
    end
end

# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------

function thermalavg(O::Matrix{Float64}, w)
    @assert size(O, 1) == length(w)
    return dot(diag(O), w)
end

function Heisenberg(O::Matrix{Float64}, spectrum, τ)
    if abs(τ) < 1e-10
        return O
    end
    n = size(O, 1)
    Out = similar(O)
    # Optimized for spectrum being ΔE matrix
    if spectrum isa AbstractMatrix
        @inbounds for b in 1:n, a in 1:n
            Out[a, b] = O[a, b] * exp(spectrum[a, b] * τ)
        end
    else
        # Fallback if vector passed
        Evals = spectrum
        @inbounds for b in 1:n, a in 1:n
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

# In-place scaling: dest = scale .* src (where scale acts on rows or cols depending on usage)
# Since propagator is diagonal exp(-beta*E), and we operate in eigenbasis:
# (exp(-tH) * O)[i,j] = exp(-t*E_i) * O[i,j] -> Row scaling
@inline function copy_scaled!(dest::Matrix{Float64}, src::Matrix{Float64}, factors)
    if factors === nothing
        copyto!(dest, src)
    else
        # Optimized loop order for column-major Julia
        @inbounds for j in axes(src, 2)
            @simd for i in axes(src, 1)
                dest[i, j] = factors[i] * src[i, j]
            end
        end
    end
    return dest
end

"""
parity(p)
Calculate the parity of a permutation using bitmask.
Zero allocation and very fast for N <= 64.
"""
function parity(p::AbstractVector{Int})
    n = length(p)
    seen = 0  # Int64 as a bitmask
    cycles = 0

    @inbounds for i in 1:n
        # Check if i-th bit is 0 (not seen)
        if (seen & (1 << (i - 1))) == 0
            cycles += 1
            j = i
            # Traverse the cycle
            while (seen & (1 << (j - 1))) == 0
                seen |= (1 << (j - 1)) # Mark j as seen
                j = p[j]
            end
        end
    end
    return (n - cycles) % 2 == 0 ? 1 : -1
end

function fermionic_sign(m::Model, perm::AbstractVector{Int}, NF::Int)
    if !m.isfermi
        return 1
    end

    # Calculate inversions for the subsequence of fermions
    # Only suitable for small NF (<= 16 is fine for O(N^2))
    inversions = 0
    @inbounds for i in 1:length(perm)
        idx_i = perm[i]
        if idx_i <= NF # It's a fermion
            for j in (i+1):length(perm)
                idx_j = perm[j]
                if idx_j <= NF
                    if idx_i > idx_j
                        inversions += 1
                    end
                end
            end
        end
    end

    return (inversions % 2 == 0) ? 1 : -1
end

# -----------------------------------------------------------------------------
# GreenN Structure
# -----------------------------------------------------------------------------
struct GreenN
    N::Int
    τ::Vector{Float64}
    orbital::Vector{Int}
    ops::Vector{Matrix{Float64}}

    function GreenN(m::Model, τ, orbital)
        N = length(τ) ÷ 2
        @assert length(τ) == length(orbital) == 2N
        ops = Vector{Matrix{Float64}}(undef, 2N)
        for i in 1:N
            ops[i] = m.c⁺[orbital[i]]
        end
        for i in (N+1):2N
            ops[i] = m.c⁻[orbital[i]]
        end
        return new(N, τ, orbital, ops)
    end
end

# -----------------------------------------------------------------------------
# Core Calculation Functions (Optimized)
# -----------------------------------------------------------------------------

"""
chain_trace_perm!(ws, m, perm, nops)

Compute Tr[ e^{-beta H} T_tau (Op1(t1) ... OpN(tn)) ] using the workspace buffers.
Zero allocation.
"""
function chain_trace_perm!(ws::GreenWorkspace, m::Model, perm::AbstractVector{Int}, nops::Int)
    # Pointers to buffers (Note: simple variable assignment, not copying data)
    current = ws.mat_A
    next_op = ws.mat_B
    product = ws.mat_C

    τs = ws.τ_buffer
    ops = ws.ops_buffer

    # 1. Start with the operator at the largest time (τ[perm[1]])
    # The trace formula is Tr [ ... U(τ2, τ1) Op1 U(τ1, β+τn) ... ] 
    # But usually we order them t1 > t2 > ... > tn.
    # Trace = Tr [ e^{-βH} e^{(β-t1)H} O1 e^{-(t1-t2)H} O2 ... e^{-tn H} ]
    #       = Tr [ e^{-(β+tn-t1)H} O1 e^{-(t1-t2)H} O2 ... ]

    idx1 = perm[1]
    idx_end = perm[nops]

    # Initial propagator: exp(-(β + τ_end - τ_start) * H)
    # Note: τs are raw times. 
    scale = propagator(m, m.β + τs[idx_end] - τs[idx1])

    # current = scale .* Ops[idx1]
    copy_scaled!(current, ops[idx1], scale)

    # 2. Chain multiply
    @inbounds for i in 2:nops
        p_curr = perm[i]
        p_prev = perm[i-1]

        # Propagator from prev time to current time: exp(-(τ_prev - τ_curr) * H)
        dt = τs[p_prev] - τs[p_curr]
        scale = propagator(m, dt)

        # Prepare next operator matrix: next_op = scale .* Op
        copy_scaled!(next_op, ops[p_curr], scale)

        # Multiply: product = current * next_op
        mul!(product, current, next_op)

        # Swap pointers: current now points to the result, product points to old buffer
        current, product = product, current
    end

    return tr(current)
end

"""
Gn(m, g, ws)
Compute Green's function with workspace (Zero allocation).
"""
function Gn(m::Model, g::GreenN, ws::GreenWorkspace)
    n_legs = 2 * g.N

    # Copy to buffer
    @inbounds for i in 1:n_legs
        ws.τ_buffer[i] = g.τ[i]
        ws.ops_buffer[i] = g.ops[i]
        ws.perm_buffer[i] = i
    end

    # Sort permutation based on time
    # We use a view to sort only the valid part of the buffer
    p_view = view(ws.perm_buffer, 1:n_legs)
    sort!(p_view, by=i -> ws.τ_buffer[i], rev=true)

    # Sign
    sign = fermionic_sign(m, p_view, n_legs) # n_legs is total legs, but sign needs N pairs? 
    # Actually fermionic_sign takes total length in your code logic.

    # Trace
    Gval = chain_trace_perm!(ws, m, p_view, n_legs)

    return sign * Gval
end

# Wrappers for compatibility or simple cases
Gn(m::Model, g::GreenN) = Gn(m, g, GreenWorkspace(m, g.N)) # Not recommended for loops

"""
dGn_dU_estimator(m, g, τp, ws)
Optimized estimator for ∂G/∂U.
"""
function dGn_dU_estimator(m::Model, g::GreenN, τp::Real, ws::GreenWorkspace)
    # -----------------------------------------------------
    # Part A: Compute G_with_D ( < T Op... D(τp) > )
    # -----------------------------------------------------
    N = 2 * g.N
    n_ext = N + 1

    # 1. Fill buffers (Op... then D)
    @inbounds for i in 1:N
        ws.τ_buffer[i] = g.τ[i]
        ws.ops_buffer[i] = g.ops[i]
        ws.perm_buffer[i] = i
    end
    ws.τ_buffer[n_ext] = τp
    ws.ops_buffer[n_ext] = m.D
    ws.perm_buffer[n_ext] = n_ext

    # 2. Sort permutation for time ordering
    p_view = view(ws.perm_buffer, 1:n_ext)
    sort!(p_view, by=i -> ws.τ_buffer[i], rev=true)

    # 3. Fermion Sign
    # The D operator is bosonic (even parity), so inserting it doesn't change 
    # the relative order of fermions. We pass the perm of the N+1 operators,
    # but tell fermionic_sign to only look for the original N legs (indices <= N).
    fsign = fermionic_sign(m, p_view, N)

    # 4. Trace
    GwD = chain_trace_perm!(ws, m, p_view, n_ext) * fsign

    # -----------------------------------------------------
    # Part B: Compute Gval * Dloc
    # -----------------------------------------------------
    # For maximum performance, we could cache Gval if calculated before.
    # Here we recalculate it efficiently using the same workspace (carefully).

    # Note: We must restore the perm buffer for N items, or just reuse the first N slots.
    # The data in τ_buffer and ops_buffer 1:N is still valid and untouched!
    # We just need to resort the first N indices.

    p_view_G = view(ws.perm_buffer, 1:N)
    @inbounds for i in 1:N
        p_view_G[i] = i # Reset indices
    end
    sort!(p_view_G, by=i -> ws.τ_buffer[i], rev=true)

    # Sign for G
    fsign_G = fermionic_sign(m, p_view_G, N)

    # Trace for G
    Gval = chain_trace_perm!(ws, m, p_view_G, N) * fsign_G

    # Local Double Occupancy
    Dloc = thermalavg(m.D, m.w)

    # Result
    return -GwD + Gval * Dloc
end

# Backward compatibility (slow, allocates)
function dGn_dU_estimator(m::Model, g::GreenN, τp::Real)
    ws = GreenWorkspace(m, g.N)
    return dGn_dU_estimator(m, g, τp, ws)
end

function G2(m::Model, g::GreenN)
    # Kept analytic G2 for testing/debugging
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

end # module