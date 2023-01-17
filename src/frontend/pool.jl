"""
    struct LoopPool{T}

    Pool of loop basis. Each loop basis corresponds to a loop variable.
    A loop variable is a linear combination of N independent loops. The combination coefficients is what we call a loop basis.
    For example, if a loop is a momentum K, then

    varibale_i = K_1*basis[1, i] + K_2*basis[2, i] + K_3*basis[3, i] + ...

# Members
- name::Symbol : name of the pool
- dim::Int     : dimension of a loop variable (for example, the dimension of a momentum-frequency loop variable is (d+1) where d is the spatial dimension)
- basis::Matrix{T}    : Matrix of (N x Nb) that stores the loop basis, where Nb is the number of loop basis (or number of loop variables).
- current::Matrix{T}  : Matrix of (dim x Nb) that stores the loop variables, where Nb is the number of loop basis (or number of loop variables).
"""
mutable struct LoopPool{T} <: AbstractVector{T}
    name::Symbol
    dim::Int #dimension
    loopNum::Int #number of independent loops
    basis::Matrix{T} # loopNum x N
    current::Matrix{T} # dim x loopNum

    function LoopPool(name::Symbol, dim::Int, loopNum::Int, type::DataType=Float64)
        basis = Matrix{type}(undef, loopNum, 0) # Nx0 matrix
        current = Matrix{type}(undef, dim, 0) # dimx0 matrix
        return new{type}(name, dim, loopNum, basis, current)
    end
    function LoopPool(name::Symbol, dim::Int, basis::AbstractVector{Vector{T}}) where {T}
        @assert isempty(basis) == false
        loopNum = length(basis[1])
        N = length(basis) #number of basis
        @assert all(x -> length(x) == loopNum, basis)
        current = rands(T, (dim, N))
        return new{T}(name, dim, loopNum, basis, current)
    end
end

Base.size(pool::LoopPool) = (size(pool.basis)[2],)
Base.IndexStyle(::Type{<:LoopPool}) = IndexLinear()
Base.getindex(pool::LoopPool, i::Int) = pool.basis[:, i]
Base.setindex!(pool::LoopPool, v, i::Int) = setindex!(pool.basis, v, i)

function update(pool::LoopPool, variable=rand(eltype(pool.current), pool.dim, pool.loopNum))
    @assert size(variable)[1] == pool.dim
    # loopNum = size(pool.basis)[1]
    loopNum = pool.loopNum

    ############### higher level LinearAlgebra call, no allocation, takes ~0.8μs ################
    mul!(pool.current, view(variable, :, 1:loopNum), pool.basis) #use view will use less memory than variable[:, 1:loopNum]
end

current(pool::LoopPool, idx) = view(pool.current, :, idx)

hasloop(pool::LoopPool) = (pool.dim > 0) && (pool.loopNum > 0)

function append(pool::LoopPool, basis::AbstractVector)
    for bi in 1:length(pool)
        if pool.basis[:, bi] ≈ basis
            # if (pool.basis[:, bi] ≈ basis) || (pool.basis[:, bi] ≈ -basis)
            return bi
        end
    end

    pool.basis = hcat(pool.basis, basis)
    pool.current = hcat(pool.current, rand(eltype(pool.current), pool.dim))
    # pool.size += 1
    # @assert pool.size <= size(pool.basis)[2] "Too many loop basis!. Increase maxSize when creates the LoopPool!"
    # pool.basis[:, pool.size] = basis
    return length(pool)
end