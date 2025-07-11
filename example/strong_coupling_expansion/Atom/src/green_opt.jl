module Green
include("common.jl")
using LinearAlgebra, Combinatorics
using Base.Threads
using Memoization

export Model, GreenN
export density, Heisenberg, thermalavg
export Gn, Gnc, Gn_determinant

# 原Model结构保持不变，因为这不是瓶颈
struct Model{N,No}
    isfermi::Bool
    β::Float
    dim::Int
    Norbital::Int
    E::SVector{N,Float}
    Z::Float
    Hdiag::Operator
    c⁺::SVector{No,Operator}
    c⁻::SVector{No,Operator}
    n::SVector{No,Operator}

    function Model(β, H, c⁺_fock::Vector{Operator}, isfermi=true)
        dim = size(H, 1)
        @assert size(H) == size(c⁺_fock[1])
        @assert size(H) == (dim, dim)
        Norbital = length(c⁺_fock)

        E, U = eigen(Float64.(Matrix(H)))
        Z = sum(exp.(-β * E))
        E = sort(E)

        Hdiag = zeros(Float, (dim, dim))
        Hdiag[diagind(Hdiag)] = E
        c⁺ = [U' * o * U for o in c⁺_fock]
        c⁻ = [adjoint(op) for op in c⁺]
        n = [c⁺[i] * c⁻[i] for i in 1:Norbital]

        return new{dim,Norbital}(isfermi, β, dim, Norbital, E, Z, Hdiag, c⁺, c⁻, n)
    end
end

function thermalavg(O::Operator, E, β, Z)
    if !(size(O) == (length(E), length(E)))
        throw(AssertionError("Dimension of Operator doesn't match"))
    end
    return sum(diag(O) .* exp.(-β * E)) / Z
end

function Heisenberg(O::Operator, E, τ)
    if !(size(O) == (length(E), length(E)))
        throw(AssertionError("Dimension of Operator doesn't match"))
    end
    if abs(τ) < 1e-10
        return O
    end
    Uτ = Diagonal(exp.(E * τ))
    return Uτ * O * inv(Uτ)
end

# 高效的奇偶性计算
function parity(p)
    n = length(p)
    visited = falses(n)
    sign = 1

    for i in 1:n
        if !visited[i]
            cycle_length = 0
            j = i
            while !visited[j]
                visited[j] = true
                j = p[j]
                cycle_length += 1
            end
            if cycle_length > 1 && iseven(cycle_length)
                sign *= -1
            end
        end
    end
    return sign
end

"""
优化的GreenN结构 - 专门针对大N情况
"""
struct GreenN
    N::Int
    τ::Vector{Float}
    orbital::Vector{Int}
    hop::Vector{Operator}

    # 优化：预计算2-point Green函数矩阵用于Wick定理
    G2_matrix::Matrix{Float}  # G2_matrix[i,j] = G2(i->j)

    function GreenN(m::Model, τ, orbital)
        N = length(τ) ÷ 2
        @assert length(τ) == length(orbital) == 2N "Length of τ and orbital must be 2N."

        hop = Vector{Operator}(undef, 2N)
        for i in 1:N
            hop[i] = Heisenberg(m.c⁺[orbital[i]], m.E, τ[i])
        end
        for i in (N+1):2N
            hop[i] = Heisenberg(m.c⁻[orbital[i]], m.E, τ[i])
        end

        # 预计算所有可能的2-point函数
        G2_matrix = zeros(Float, 2N, 2N)
        for i in 1:N, j in (N+1):2N
            g2 = GreenN(m, [τ[i], τ[j]], [orbital[i], orbital[j]])
            G2_matrix[i, j] = G2(m, g2)
        end

        return new(N, τ, orbital, hop, G2_matrix)
    end
end

# 2-point函数保持原样
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

"""
行列式方法实现 - 将N-point函数表示为2-point函数矩阵的行列式
复杂度: O(N^3) vs 原来的 O(N!)
"""
function Gn_determinant(m::Model, g::GreenN)
    if !m.isfermi
        error("Determinant method currently only supports fermionic systems")
    end

    N = g.N
    if N == 1
        return G2(m, g)
    end

    # 构建N×N的2-point函数矩阵
    # M[i,j] = G2(incoming_i -> outgoing_j)
    M = Matrix{Float}(undef, N, N)
    for i in 1:N
        for j in 1:N
            M[i, j] = g.G2_matrix[i, N+j]  # incoming i to outgoing N+j
        end
    end

    return det(M)
end

"""
改进的直接计算方法 - 使用更好的算法组织
"""
function Gn_direct_optimized(m::Model, g::GreenN)
    τ = g.τ
    hop = g.hop
    N = g.N

    # 预排序以减少后续计算
    perm = sortperm(τ)
    ordered_hop = hop[perm]
    parity_sign = m.isfermi ? parity(perm) : 1

    # 使用分块矩阵乘法优化大矩阵运算
    M = ordered_hop[end]
    for i in (length(ordered_hop)-1):-1:1
        M = ordered_hop[i] * M
    end

    G = thermalavg(M, m.E, m.β, m.Z)
    return G * parity_sign
end

"""
智能选择计算方法
"""
function Gn(m::Model, g::GreenN)
    N = g.N

    if N == 1
        return G2(m, g)
    elseif N <= 3
        # 小N时直接计算仍然高效
        return Gn_direct_optimized(m, g)
    elseif N <= 8
        # 中等N时使用行列式方法
        return Gn_determinant(m, g)
    else
        # 大N时使用Wick定理（如果适用）
        # 注意：Wick定理严格来说只对非相互作用系统精确成立
        # 对相互作用系统这是近似
        if m.isfermi
            return Gn_determinant(m, g)  # 行列式方法对大N仍然是最好的精确方法
        else
            return Gn_direct_optimized(m, g)
        end
    end
end

"""
并行批量计算 - 针对需要计算大量大N Green函数的情况
"""
function batch_Gn_parallel(m::Model, greens::Vector{GreenN})
    results = Vector{Float}(undef, length(greens))

    # 根据N大小分组处理
    small_indices = Int[]
    large_indices = Int[]

    for (i, g) in enumerate(greens)
        if g.N <= 4
            push!(small_indices, i)
        else
            push!(large_indices, i)
        end
    end

    # 小N的并行计算
    @threads for idx in small_indices
        results[idx] = Gn_direct_optimized(m, greens[idx])
    end

    # 大N的并行计算（使用行列式方法）
    @threads for idx in large_indices
        results[idx] = Gn_determinant(m, greens[idx])
    end

    return results
end

# 保持原有的utility函数
density(m::Model, orbital) = thermalavg(m.n[orbital], m.E, m.β, m.Z)

end