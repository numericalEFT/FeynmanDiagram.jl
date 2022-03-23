module Green
include("common.jl")
using LinearAlgebra, Combinatorics
export Model, GreenN
export density, Heisenberg, thermalavg
export Gn, Gnc
# include("hilbert.jl")
# using .Hilbert

"""
Model(H::Operator, c⁺::Vector{Operator}, c⁻::Vector{Operator})

Construct a model.  All input operators are in the Fock space. But all members of Model struct is in the eigenspace.
Operator should be a matrix type.

# Arguments:
- β: inverse temperature
- H: Hamiltonian in the Fock space
- c⁺, c⁻: creation/anniliation operators for each orbital, all in the Fock space
"""
struct Model{N,No}
    isfermi::Bool
    β::Float
    dim::Int # dim of Hilbert space
    Norbital::Int # number of orbitals
    E::SVector{N,Float} # eigen energy
    Z::Float # partition sum
    H::Operator # diagnoalized Hamitlonian
    c⁺::SVector{No,Operator} # creation operator in the eigenspace
    c⁻::SVector{No,Operator} # anniliation operator in the eigenspace
    n::SVector{No,Operator}  # density operator in the eigenspace

    function Model(β, H, _c⁺, isfermi = true)
        dim = size(H)[1]
        @assert size(H) == size(_c⁺[1])
        @assert size(H) == (dim, dim)
        Norbital = length(_c⁺)

        E, U = eigen(Float64.(Matrix(H)))  # U'*H*U will diagonalize the Hamiltonian
        Z = sum(exp.(-β * E))
        E = sort(E)

        Hdiag = zeros(Float, (dim, dim))
        Hdiag[diagind(Hdiag)] = E
        c⁺ = [U' * o * U for o in _c⁺]
        c⁻ = [U' * o' * U for o in _c⁺]
        n = [c⁺[i] * c⁻[i] for i in 1:Norbital]

        return new{dim,Norbital}(isfermi, β, dim, Norbital, E, Z, Hdiag, c⁺, c⁻, n)
    end
end

function thermalavg(o::Operator, E, β, Z)
    if !(size(o) == (length(E), length(E)))
        throw(AssertionError("Dimension of Operator[$(o.m), $(o.n)] doesn't match with $(length(E))"))
    end
    avg = 0.0
    for di in 1:length(E)
        avg += o[di, di] * exp(-β * E[di])
    end
    return avg / Z
end
"""
Heisenberg(o::Operator, E, τ)

   Transform operator o into Heisenberg picture
    exp(H * τ) * o * exp(-H * τ)
"""
function Heisenberg(o::Operator, E, τ)
    if !(size(o) == (length(E), length(E)))
        throw(AssertionError("Dimension of Operator[$(o.m), $(o.n)] doesn't match with $(length(E))"))
    end
    if abs(τ) < 1.0e-10
        return o  # if τ≈0, no evoluation is needed
    end
    oh = copy(o)
    for ri in 1:length(E)
        oh[ri, :] *= exp(E[ri] * τ)
    end
    for ci in 1:length(E)
        oh[:, ci] *= exp(-E[ci] * τ)
    end
    return oh
end

function parity(p)
    """
    calculate the parity of a given permutation of the array [1, 2, 3, ...]
    """
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
    length(even_cycles) % 2 == 0 ? 1 : -1
end

struct Partition
    """
    2-partition of 2N-legs of Green's functions
    """
    l::Vector{Vector{Int}}
    r::Vector{Vector{Int}}
    sign::Vector{Int}
    function Partition(order)
        N = 2^order - 2
        function addincoming(l, r)
            # map 1, 2, 3, ... to order+1, order+2, ...
            l = l .+ order
            r = r .+ order
            nl, nr = length(l), length(r)
            inl = [i for i in 1:nl]
            inr = [i for i in nl+1:order]
            append!(inl, l)
            append!(inr, r)
            return inl, inr
        end

        left = Vector{Vector{Int}}(undef, N)
        right = similar(left)
        sign = Vector{Int}(undef, N)
        for (si, s) in enumerate(partitions(1:order, 2))
            # generate 2-partition of the outgoing legs
            # all incoming legs are kept the same
            s1, s2 = addincoming(s[1], s[2])
            left[2si-1], right[2si-1] = s1, s2
            sign[2si-1] = parity(vcat(s1, s2))

            q1, q2 = addincoming(s[2], s[1])
            left[2si], right[2si] = q1, q2
            sign[2si] = parity(vcat(q1, q2))
            # println("($s1, $s2) -> $p1, ($q1, $q2) -> $p2")
        end
        return new(left, right, sign)
    end
end

"""
GreenN(m::Model, τ, orbital, isfermi=true)

Construct struct to store the variables to evaulate N-body Green's functions.
The leg index is assumed to be [1, 2, 3, ...,2N], where the incoming legs are 1:N, and the outgoing legs are N+1:2N
The full Green's function is defined as,
```math
Gn=<Tτ c(1)c(2)...c(N)c⁺(N+1)...c⁺(2N)>
```
e.g.,
1->------->-3
    | G4 |        
2->------->-4
All other Green's function are derived from the above full Green's function

#Arguments
- m: Model struct
- τ: array stores the imaginary-time of legs
- orbital: array stores the orbital/spin of legs
- N: particle number (Note that the leg numebr is 2N)
"""
struct GreenN
    N::Int # n-body or 2n-point
    τ::Vector{Float} # 1:2n, array of imaginary-time
    orbital::Vector{Int} # 1:2n, array of orbitals 
    hop::Vector{Operator}
    partition::Vector{Partition}
    function GreenN(m::Model, τ, orbital, N = Int(length(τ) / 2))
        # @assert length(τ) == length(orbital)
        # N = Int(length(τ) / 2)

        hop = Vector{Operator}(undef, 2N)
        for idx in 1:N
            hop[idx] = Heisenberg(m.c⁺[orbital[idx]], m.E, τ[idx])
        end
        for idx in N+1:2N
            hop[idx] = Heisenberg(m.c⁻[orbital[idx]], m.E, τ[idx])
        end

        partition = [Partition(o) for o in 2:N]
        return new(N, τ, orbital, hop, partition)
    end
end

function density(m::Model, orbital)
    n = 0.0
    for (ei, e) in enumerate(m.E)
        n += exp(-m.β * e) * m.n[orbital][ei, ei]
    end
    return n / m.Z
end

"""
Gn=<Tτ c(1)c(2)...c(N)c⁺(N+1)...c⁺(2N)>
e.g.,
1->------->-3
    | G4 |        
2->------->-4
"""
@inline function Gn(m::Model, g::GreenN)
    return Gn(m, g, [i for i in 1:g.N*2])
end

function Gn(m::Model, g::GreenN, idx, level = 1)
    # printstyled("...."^level, "start $level Gn-$idx\n", color=:green)
    if g.N == 1 # fast treatment of one-body Green's function
        i, o = idx[1], idx[2]
        println(g.hop[o])
        println(g.hop[i])
        if g.τ[i] < g.τ[o]
            return thermalavg(g.hop[o] * g.hop[i], m.E, m.β, m.Z)
        else
            G = thermalavg(g.hop[i] * g.hop[o], m.E, m.β, m.Z)
            m.isfermi && (G *= -1.0)
            # printstyled("...."^level, "get $level Gn-$idx -> $(idx[perm])*$(parity(perm))   $G\n", color=:green)
            return G
        end
    else
        perm = sortperm(g.τ[idx])
        # if length(idx) == 4
        #     println(g.τ[idx], "  and  ", g.orbital[idx], " --> ", perm)
        # end
        # M=c(τ_{2N})...c(2)c(1)
        M = prod(g.hop[reverse(idx[perm])])  # idx[perm] gives the idx array with the correct time ordering
        # note that the storage order of hop operators are the opposite of the phyical order of the operators
        G = thermalavg(M, m.E, m.β, m.Z) # G=tr(exp(-βH) M)/Z
        m.isfermi && (G *= parity(perm))
        # printstyled("...."^level, "get $level Gn-$idx -> $(idx[perm])*$(parity(perm))   $G\n", color=:green)
        return G
    end
end

"""
1->------->-4       1->------->-4     1->---   --->-4
    | G4 |      -                  +         X   
2->------->-3       2->------->-3     2->---   --->-3
"""
@inline function Gnc(m::Model, g::GreenN)
    return Gnc(m, g, [i for i in 1:g.N*2])
end

function Gnc(m::Model, g::GreenN, idx, level = 1)
    # printstyled("...."^level, "start $level Gc-$idx\n", color=:red)
    if length(idx) == 2 # case of G2
        return Gn(m, g, idx, level)
    end

    G = Gn(m, g, idx, level)

    order = Int(length(idx) / 2)
    p = g.partition[order-1]
    for (li, l) in enumerate(p.l)
        # printstyled("...."^level, "begin $(level): ", idx, " => ", idx[l], ", ", idx[p.r[li]], ": parity: ", 1, " --> ", p.sign[li], "\n", color=:yellow)
        disconnected = Gnc(m, g, idx[l], level + 1) * Gn(m, g, idx[p.r[li]], level + 1)
        m.isfermi && (disconnected *= p.sign[li])
        G -= disconnected
        # printstyled("...."^level, "end $(level)\n", color=:yellow)
    end
    # printstyled("...."^level, "end $level Gc-$idx $G\n", color=:red)
    return G
end

function G2(m::Model, g::GreenN)
    if g.N != 1
        throw(AssertionError("G2 function is for n=2"))
    end
    τi, τo = g.τ[1], g.τ[2]
    β, E = m.β, m.E
    c⁺, c⁻ = m.c⁺[g.orbital[1]], m.c⁻[g.orbital[2]] # in is creation operator, out is anniliation operator

    G = 0.0
    if (τi < τo)
        for j in 1:m.dim
            for k in 1:m.dim
                G += exp(-(β - τo + τi) * E[j] - (τo - τi) * E[k]) * c⁻[j, k] * c⁺[k, j]
            end
        end
    else # equal-time means 0^-
        for j in 1:m.dim
            for k in 1:m.dim
                G += -exp(-(β - τi + τo) * E[j] - (τi - τo) * E[k]) * c⁺[j, k] * c⁻[k, j]
            end
        end
    end
    return G / m.Z
end

end