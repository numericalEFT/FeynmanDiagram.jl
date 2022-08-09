"""
    struct CachedPool{O,T}

        Use this pool to host the objects that are heavy to evaluate so that one wants to cache their status.
        The user should defines a compare 

# Members
- name::Symbol : name of the pool
- object::O    : object
- current::T      : current status
- new::T       : the new status wants to assign later
- version::Int128 : the current version
- excited::Bool   : if set to excited, then the current status needs to be replaced with the new status
"""
mutable struct CachedPool{O,T}
    name::Symbol
    object::Vector{O}
    current::Vector{T}
    new::Vector{T}
    version::Vector{Int128}
    excited::Vector{Bool}

    function CachedPool(name::Symbol, objType::DataType, weightType::DataType)
        object = Vector{objType}(undef, 0)
        current = Vector{weightType}(undef, 0)
        _new = Vector{weightType}(undef, 0)
        version = Vector{Int128}(undef, 0)
        excited = Vector{Bool}(undef, 0)
        return new{objType,weightType}(name, object, current, _new, version, excited)
    end
    # function CachedPool{T}(obj::Vector{O}) where {O,T}
    #     weight = zeros(5, length(obj))
    #     return new{O,T}(obj, weight)
    # end
end

Base.length(pool::CachedPool) = length(pool.object)
Base.size(pool::CachedPool) = size(pool.object)
Base.show(io::IO, pool::CachedPool) = print(io, pool.object)
# Base.view(pool::Pool, inds...) = Base.view(pool.pool, inds...)

#index interface for Pool
Base.getindex(pool::CachedPool, i) = pool.object[i]
Base.setindex!(pool::CachedPool, v, i) = setindex!(pool.object, v, i)
Base.firstindex(pool::CachedPool) = 1
Base.lastindex(pool::CachedPool) = length(pool)

# iterator interface
function Base.iterate(pool::CachedPool)
    if length(pool) == 0
        return nothing
    else
        return (pool.object[1], 1)
    end
end

function Base.iterate(pool::CachedPool, state)
    if state >= length(pool) || length(pool) == 0
        return nothing
    else
        return (pool.object[state+1], state + 1)
    end
end


function append(pool::CachedPool, object)
    # @assert para isa eltype(pool.pool)
    # for (oi, o) in enumerate(pool.object)
    #     if o == object
    #         return oi #existing obj
    #     end
    # end

    id = length(pool.object) + 1
    push!(pool.object, object)

    return id #new momentum
end

function initialize!(pool::CachedPool{O,T}) where {O,T}
    N = length(pool)
    pool.current = zeros(T, N)
    pool.new = zeros(T, N)
    pool.version = ones(T, N)
    pool.excited = Vector{Bool}(zeros(Int, N))
end

function updateAll(pool::CachedPool, ignoreCache::Bool, eval::Function; kwargs...)
    N = length(pool.object)
    if ignoreCache
        T = eltype(pool.current)
        for (idx, o) in enumerate(pool.object)
            pool.current[idx] = T(eval(obj; kwargs...))
        end
    else
        error("not implemented!")
        # pool.version .= 1
        # pool.excited .= false
    end
end

"""
    struct LoopPool{T}

    Pool of loop basis. Each loop basis corresponds to a loop variable.
    A loop variable is a linear combination of N independent loops. The combination coefficients is what we call a loop basis.
    For example, if a loop is a momentum K, then

    varibale_i = K_1*basis[1, i] + K_2*basis[2, i] + K_3*basis[3, i] + ...

# Members
- name::Symbol : name of the pool
- dim::Int     : dimension of a loop variable (for example, the dimension of a momentum-frequency loop variable is (d+1) where d is the spatial dimension)
- N::Int       : number of independent loops (dimension of loop basis)
- basis::Matrix{T}    : Matrix of (N x Nb) that stores the loop basis, where Nb is the number of loop basis (or number of loop variables).
- current::Matrix{T}  : Matrix of (dim x Nb) that stores the loop variables, where Nb is the number of loop basis (or number of loop variables).
"""
mutable struct LoopPool{T}
    name::Symbol
    dim::Int #dimension
    N::Int #number of basis
    basis::Matrix{T}
    current::Matrix{T}

    function LoopPool(name::Symbol, dim::Int, N::Int, type::DataType=Float64)
        basis = Matrix{type}(undef, N, 0) # Nx0 matrix
        current = Matrix{type}(undef, dim, 0) # dimx0 matrix
        return new{type}(name, dim, N, basis, current)
    end
end

Base.length(pool::LoopPool) = size(pool.basis)[2]
Base.size(pool::LoopPool) = length(pool)
Base.show(io::IO, pool::LoopPool) = print(io, pool.basis)
# Base.view(pool::Pool, inds...) = Base.view(pool.pool, inds...)

#index interface for Pool
Base.getindex(pool::LoopPool, i) = pool.basis[:, i]
Base.setindex!(pool::LoopPool, v, i) = setindex!(pool.basis, v, i)
Base.firstindex(pool::LoopPool) = 1
Base.lastindex(pool::LoopPool) = length(pool)

# iterator interface
function Base.iterate(pool::LoopPool)
    if length(pool) == 0
        return nothing
    else
        return (pool.basis[:, 1], 1)
    end
end

function Base.iterate(pool::LoopPool, state)
    if state >= length(pool) || length(pool) == 0
        return nothing
    else
        return (pool.basis[:, state+1], state + 1)
    end
end

function update(pool::LoopPool, variable=rand(eltype(pool.current), pool.dim, pool.N))
    loopNum = size(pool.basis)[1]
    pool.current = view(variable, :, 1:loopNum) * pool.basis
    # B = view(pool.basis, 1:loopNum, :)
    # B = view(pool.basis, 1:loopNum, :)
    # C = pool.current[:, 1:length(pool)]
    # LinearAlgebra.mul!(pool.current, A, pool.basis)
    # LinearAlgebra.BLAS.gemm!('N', 'N', false, variable[:, 1:loopNum], pool.basis[1:loopNum, 1:length(pool)], false, pool.current[:, 1:length(pool)])
end

# current(pool::LoopPool, idx) = pool.current[:, idx]
current(pool::LoopPool, idx) = view(pool.current, :, idx)

hasloop(pool::LoopPool) = (pool.dim > 0) && (pool.N > 0)

function append(pool::LoopPool, basis::AbstractVector)
    for bi in 1:length(pool)
        if pool.basis[:, bi] ≈ basis
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
