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
    Gn_1stderiv_estimator(m, g, Op, τp, ws)

    Generic estimator for ∂G/∂λ where H(λ) = H₀ - λ⋅Op.
Used for:
  - μ (Op = Ntot)
  - h (Op = M)
  - U (Op = D, if using Hubbard-U parameter)

Formula: ∂G/∂λ = ∫ dτ <T G Op(τ)> - <G><Op(τ)>
"""
function Gn_1stderiv_estimator(m::Model, g::GreenN, Op::Matrix{Float64}, τp::Real, ws::GreenWorkspace)
    # 1. Setup Buffers
    N = 2 * g.N
    n_ext = N + 1

    @inbounds for i in 1:N
        ws.τ_buffer[i] = g.τ[i]
        ws.ops_buffer[i] = g.ops[i]
        ws.perm_buffer[i] = i
    end
    ws.τ_buffer[n_ext] = τp
    ws.ops_buffer[n_ext] = Op
    ws.perm_buffer[n_ext] = n_ext

    # 2. Compute <T Legs... Op>
    # Sort permutation for time ordering
    p_view = view(ws.perm_buffer, 1:n_ext)
    sort!(p_view, by=i -> ws.τ_buffer[i], rev=true)

    # 3. Fermion Sign
    # Fermion sign depends only on the first N indices (the fermions)
    # The Op operator is bosonic (even parity), so inserting it doesn't change the relative order of fermions. 
    fsign = fermionic_sign(m, p_view, N)
    GwOp = chain_trace_perm!(ws, m, p_view, n_ext) * fsign

    # 4. Compute Gval (Legs only)
    # We re-sort the first N indices. Data in buffer 1:N is untouched.
    p_view_G = view(ws.perm_buffer, 1:N)
    @inbounds for i in 1:N
        p_view_G[i] = i
    end
    sort!(p_view_G, by=i -> ws.τ_buffer[i], rev=true)

    fsign_G = fermionic_sign(m, p_view_G, N)
    Gval = chain_trace_perm!(ws, m, p_view_G, N) * fsign_G

    # 5. Local average <Op>
    Oploc = thermalavg(Op, m.w)

    # Result: Connected part (with minus sign for derivative definition usually)
    # If just measuring correlator <T G Op>, remove the minus.
    # Assuming derivative definition:
    return GwOp - Gval * Oploc
end

"""
    Gn_2ndderiv_estimator(m, g, Op1, τ1, Op2, τ2, ws)

Efficient estimator for the second order derivative with respect to fields coupling to Op1 and Op2.
Computes the connected 3-point correlator <T Legs... Op1 Op2>_connected.
"""
function Gn_2ndderiv_estimator(m::Model, g::GreenN,
    Op1::Matrix{Float64}, τ1::Real,
    Op2::Matrix{Float64}, τ2::Real,
    ws::GreenWorkspace)

    N = 2 * g.N
    idx_1 = N + 1
    idx_2 = N + 2
    total_len = N + 2

    # =========================================================
    # 1. Fill Buffer ONCE (Legs + Op1 + Op2)
    # =========================================================
    @inbounds for i in 1:N
        ws.τ_buffer[i] = g.τ[i]
        ws.ops_buffer[i] = g.ops[i]
    end
    ws.τ_buffer[idx_1] = τ1
    ws.ops_buffer[idx_1] = Op1
    ws.τ_buffer[idx_2] = τ2
    ws.ops_buffer[idx_2] = Op2

    # Averages <Op> (Fast, scalar)
    O1_loc = thermalavg(Op1, m.w)
    O2_loc = thermalavg(Op2, m.w)

    # =========================================================
    # 2. Term A: < T Legs... Op1 Op2 >
    # =========================================================
    # Reset perm buffer for all items
    @inbounds for i in 1:total_len
        ws.perm_buffer[i] = i
    end

    p_view_all = view(ws.perm_buffer, 1:total_len)
    sort!(p_view_all, by=i -> ws.τ_buffer[i], rev=true)

    # Sign only depends on original legs (indices <= N)
    fsign = fermionic_sign(m, p_view_all, N)
    GwO1O2 = chain_trace_perm!(ws, m, p_view_all, total_len) * fsign

    # =========================================================
    # 3. Term B: < T Legs... Op1 >
    # =========================================================
    # We use the buffer slots 1:N and idx_1
    # We must construct a specific permutation list
    # (We cannot overwrite the main buffer, but we can overwrite perm_buffer)

    # Construct indices for this subset
    p_subset = view(ws.perm_buffer, 1:(N+1))
    @inbounds for i in 1:N
        p_subset[i] = i
    end
    p_subset[N+1] = idx_1

    sort!(p_subset, by=i -> ws.τ_buffer[i], rev=true)
    fsign_1 = fermionic_sign(m, p_subset, N)
    GwO1 = chain_trace_perm!(ws, m, p_subset, N + 1) * fsign_1

    # =========================================================
    # 4. Term C: < T Legs... Op2 >
    # =========================================================
    # Indices: 1:N and idx_2
    @inbounds for i in 1:N
        p_subset[i] = i
    end
    p_subset[N+1] = idx_2

    sort!(p_subset, by=i -> ws.τ_buffer[i], rev=true)
    fsign_2 = fermionic_sign(m, p_subset, N)
    GwO2 = chain_trace_perm!(ws, m, p_subset, N + 1) * fsign_2

    # =========================================================
    # 5. Term D: < T Legs... > (Gval)
    # =========================================================
    p_legs = view(ws.perm_buffer, 1:N)
    @inbounds for i in 1:N
        p_legs[i] = i
    end
    sort!(p_legs, by=i -> ws.τ_buffer[i], rev=true)

    fsign_G = fermionic_sign(m, p_legs, N)
    Gval = chain_trace_perm!(ws, m, p_legs, N) * fsign_G

    # =========================================================
    # 6. Term E: < T Op1 Op2 > (Bosonic Correlator)
    # =========================================================
    # Indices: idx_1, idx_2
    # We reuse the start of perm buffer
    p_corr = view(ws.perm_buffer, 1:2)
    p_corr[1] = idx_1
    p_corr[2] = idx_2
    sort!(p_corr, by=i -> ws.τ_buffer[i], rev=true)

    # No fermion sign for bosonic operators
    O1O2_corr = chain_trace_perm!(ws, m, p_corr, 2)

    # =========================================================
    # Final Combination (Connected Part)
    # =========================================================
    # Formula: <AB>_c = <AB> - <A><B>
    # Here A = GreenLegs, B = (Op1, Op2 insertions)
    # The full expanded form for ∂²G/∂U² is:

    # < T G O1 O2 >
    # - < T G O1 > <O2>
    # - < T G O2 > <O1>
    # + 2 < G > <O1> <O2>
    # - < G > ( < T O1 O2 > - <O1><O2> )

    term1 = GwO1O2
    term2 = GwO1 * O2_loc
    term3 = GwO2 * O1_loc
    term4 = 2 * Gval * O1_loc * O2_loc
    term5 = Gval * (O1O2_corr - O1_loc * O2_loc)

    return term1 - term2 - term3 + term4 - term5
end

"""
dGn_dU_estimator(m, g, τp, ws)

Correct estimator for ∂G/∂U = - ∫ dτ [ <T G D(τ)> - <G><D> ].
It inserts D at time τp (sampled by MC) and computes the trace.
"""
dGn_dU_estimator(m::Model, g::GreenN, τp::Real, ws::GreenWorkspace) = (-1) * Gn_1stderiv_estimator(m, g, m.D, τp, ws)
# Backward compatibility (slow, allocates)
function dGn_dU_estimator(m::Model, g::GreenN, τp::Real)
    ws = GreenWorkspace(m, g.N)
    return dGn_dU_estimator(m, g, τp, ws)
end

dGn_dμ_estimator(m::Model, g::GreenN, τp::Real, ws::GreenWorkspace) = Gn_1stderiv_estimator(m, g, m.Ntot, τp, ws)
function dGn_dμ_estimator(m::Model, g::GreenN, τp::Real)
    ws = GreenWorkspace(m, g.N)
    return dGn_dμ_estimator(m, g, τp, ws)
end

dGn_dh_estimator(m::Model, g::GreenN, τp::Real, ws::GreenWorkspace) = Gn_1stderiv_estimator(m, g, m.M, τp, ws)
function dGn_dh_estimator(m::Model, g::GreenN, τp::Real)
    ws = GreenWorkspace(m, g.N)
    return dGn_dh_estimator(m, g, τp, ws)
end


"""
    chain_trace_moments(ws, m, perm, nops, diag_Op)

Compute trace moments for a diagonal operator `diag_Op` (vector).
Returns (Tr[G], Tr[Op * G], Tr[Op^2 * G]).
Cost: One matrix multiplication chain + 3 dot products.
"""
function chain_trace_moments(ws::GreenWorkspace, m::Model, perm::AbstractVector{Int}, nops::Int, diag_Op::AbstractVector{Float64})
    # 1. Compute G
    current = ws.mat_A
    next_op = ws.mat_B
    product = ws.mat_C
    τs, ops = ws.τ_buffer, ws.ops_buffer

    idx1, idx_end = perm[1], perm[nops]

    scale = propagator(m, m.β + τs[idx_end] - τs[idx1])
    copy_scaled!(current, ops[idx1], scale)

    @inbounds for i in 2:nops
        p_curr, p_prev = perm[i], perm[i-1]
        dt = τs[p_prev] - τs[p_curr]
        scale = propagator(m, dt)
        copy_scaled!(next_op, ops[p_curr], scale)
        mul!(product, current, next_op)
        current, product = product, current
    end

    # 2. compute trace moments for diagonal operator
    # Tr[Op * M] = sum(Op[i] * M[i,i]) 
    tr_0 = 0.0 # <G>
    tr_1 = 0.0 # <Op G>
    tr_2 = 0.0 # <Op^2 G>

    @inbounds for i in 1:m.dim
        val = current[i, i]
        op_val = diag_Op[i]

        tr_0 += val
        tr_1 += val * op_val
        tr_2 += val * op_val * op_val
    end

    return tr_0, tr_1, tr_2
end

"""
    dGn_dβ_estimator(m, g, ws)

Estimator for ∂G/∂β = -(<HG> - <H><G>).
Efficient: No extra operator insertion, just weighted trace.
"""
function dGn_dβ_estimator(m::Model, g::GreenN, ws::GreenWorkspace)
    n_legs = 2 * g.N

    # 1. Setup Buffers & Sort
    @inbounds for i in 1:n_legs
        ws.τ_buffer[i] = g.τ[i]
        ws.ops_buffer[i] = g.ops[i]
        ws.perm_buffer[i] = i
    end
    p_view = view(ws.perm_buffer, 1:n_legs)
    sort!(p_view, by=i -> ws.τ_buffer[i], rev=true)
    fsign = fermionic_sign(m, p_view, n_legs)

    # 2. Compute <G> and <HG>
    tr_G, tr_HG, _ = chain_trace_moments(ws, m, p_view, n_legs, m.E)

    Gval = fsign * tr_G
    HGval = fsign * tr_HG

    # 3. Compute <H>
    E_avg = dot(m.E, m.w)

    # 4. Result - <HG> + <H><G>
    return -(HGval - E_avg * Gval)
end

function dGn_dβ_estimator(m::Model, g::GreenN)
    ws = GreenWorkspace(m, g.N)
    return dGn_dβ_estimator(m, g, ws)
end

"""
    d2Gn_dβ2_estimator(m, g, ws)

Efficient estimator for the second derivative ∂²G/∂β².
Formula: (<H²G> - <H²><G>) - 2<H>(<HG> - <H><G>)
Cost: Same as computing G (O(dim) overhead only). No τ_p sampling needed.
"""
function d2Gn_dβ2_estimator(m::Model, g::GreenN, ws::GreenWorkspace)
    n_legs = 2 * g.N

    # 1. Setup Buffers & Sort
    @inbounds for i in 1:n_legs
        ws.τ_buffer[i] = g.τ[i]
        ws.ops_buffer[i] = g.ops[i]
        ws.perm_buffer[i] = i
    end
    p_view = view(ws.perm_buffer, 1:n_legs)
    sort!(p_view, by=i -> ws.τ_buffer[i], rev=true)

    # 2. Fermion Sign
    fsign = fermionic_sign(m, p_view, n_legs)

    # 3. Compute Traces: <G>, <HG>, <H²G> in one pass
    # We pass 'm.E' (Eigenvalues) as the diagonal operator
    tr_G, tr_HG, tr_H2G = chain_trace_moments(ws, m, p_view, n_legs, m.E)

    # Apply sign
    Gval = fsign * tr_G
    HGval = fsign * tr_HG
    H2Gval = fsign * tr_H2G

    # 4. Compute Static Moments <H> and <H²>
    # m.w are the normalized Boltzmann weights
    E_avg = dot(m.E, m.w)
    E2_avg = dot(m.E .^ 2, m.w)

    # 5. Assemble Formula
    # Term A: <H²G> - <H²><G>
    term_A = H2Gval - E2_avg * Gval

    # Term B: <HG> - <H><G>
    term_B = HGval - E_avg * Gval

    # Result = Term A - 2 * <H> * Term B
    return term_A - 2 * E_avg * term_B
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