"""
    mutable struct Cache{O,T}

        Cache an object that is so heavy that one wants to cache their status without evaluating them on the fly.
        The user should defines a compare 

# Members
- object::O    : parameters that are needed to evaluate the object
- id::Int      : index of the object
- factor::T    : reweight factor of the current object, defaut to be 1
- curr::T      : current status
- new::T       : the new status wants to assign later
- version::Int128 : the current version
- excited::Bool   : if set to excited, then the current status needs to be replaced with the new status
"""
mutable struct Cache{O,T}
    object::O
    id::Int
    curr::T
    new::T
    version::Int128
    excited::Bool #if set to excited, then the current weight needs to be replaced with the new weight
    function Cache(object::O, curr::T, id, version = 1, excited = false) where {O,T}
        return new{O,T}(object, id, curr, curr, version, excited)
    end
    function Cache{T}(object::O, id, version = 1, excited = false) where {O,T}
        return new{O,T}(object, id, zero(T), zero(T), version, excited)
    end
    function Cache{O,T}(object, curr, id, version = 1, excited = false) where {O,T}
        # println("object: ", object)
        return new{O,T}(object, id, curr, curr, version, excited)
    end
end

Base.show(io::IO, obj::Cache) = print(io, "id=$(obj.id): obj=$(obj.object) curr=$(obj.curr)")

"""
    struct Pool{O}

        Pool of (cached) objects.

# Members
- pool::Vector{O} : Vector that hosts the (cached) object
"""
struct Pool{O}
    pool::Vector{O}

    function Pool{O}() where {O}
        pool = Vector{O}(undef, 0)
        return new{O}(pool)
    end
    function Pool(obj::Vector{O}) where {O}
        return new{O}(obj)
    end
end

function cachedPool(objType::DataType, currType::DataType)
    COBJ = Cache{objType,currType}
    return Pool{COBJ}()
end

function uncachedPool(objType::DataType)
    return Pool{objType}()
end

Base.length(pool::Pool) = length(pool.pool)
Base.size(pool::Pool) = size(pool.pool)
Base.show(io::IO, pool::Pool) = print(io, pool.pool)
# Base.view(pool::Pool, inds...) = Base.view(pool.pool, inds...)

#index interface for Pool
Base.getindex(pool::Pool, i) = pool.pool[i]
Base.setindex!(pool::Pool, v, i) = setindex!(pool.pool, v, i)
Base.firstindex(pool::Pool) = 1
Base.lastindex(pool::Pool) = length(pool)

# iterator interface
function Base.iterate(pool::Pool)
    if length(pool) == 0
        return nothing
    else
        return (pool.pool[1], 1)
    end
end

function Base.iterate(pool::Pool, state)
    if state >= length(pool) || length(pool) == 0
        return nothing
    else
        return (pool[state+1], state + 1)
    end
end


"""
    struct PoolwithParameter{PARA,O}

        Pool of (cached) objects, support additional parameter.

# Members
- para::PARA  : additional user-defined parameter
- pool::Vector{O} : Vector that hosts the (cached) object
"""
struct PoolwithParameter{PARA,O}
    para::PARA
    pool::Vector{O}

    function PoolwithParameter{O}(para::PARA) where {PARA,O}
        pool = Vector{O}(undef, 0)
        return new{PARA,O}(para, pool)
    end
    function PoolwithParameter(_para::PARA, obj::Vector{O}) where {PARA,O}
        return new{PARA,O}(_para, obj)
    end
end

Base.length(pool::PoolwithParameter) = length(pool.pool)
Base.size(pool::PoolwithParameter) = size(pool.pool)
Base.show(io::IO, pool::PoolwithParameter) = print(io, "para: $(pool.para), pool: $(pool.pool)")
# Base.view(pool::PoolwithParameter, inds...) = Base.view(pool.pool, inds...)

#index interface for Pool
Base.getindex(pool::PoolwithParameter, i) = pool.pool[i]
Base.setindex!(pool::PoolwithParameter, v, i) = setindex!(pool.pool, v, i)
Base.firstindex(pool::PoolwithParameter) = 1
Base.lastindex(pool::PoolwithParameter) = length(pool)

# iterator interface
function Base.iterate(pool::PoolwithParameter)
    if length(pool) == 0
        return nothing
    else
        return (pool.pool[1], 1)
    end
end

function Base.iterate(pool::PoolwithParameter, state)
    if state >= length(pool) || length(pool) == 0
        return nothing
    else
        return (pool[state+1], state + 1)
    end
end

"""
    struct CachedPool{O}

        Pool of (cached) objects.

# Members
- pool::Vector{O} : Vector that hosts the (cached) object
"""
struct CachedPool{O,T}
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


function isCached(pool)
    ObjectInPoolType = eltype(pool.pool)
    return ObjectInPoolType <: Cache
end

function append(pool, object, curr, isCached)
    # @assert para isa eltype(pool.pool)
    for (oi, o) in enumerate(pool.pool)
        if isCached
            if o.object == object
                return oi #existing obj
            end
        else
            if o == object
                return oi #existing obj
            end
        end
    end

    id = length(pool) + 1
    if isCached
        CACHEDOBJECT = eltype(pool.pool)
        O = fieldtype(CACHEDOBJECT, :object)
        W = fieldtype(CACHEDOBJECT, :curr)
        obj = Cache{O,W}(object, curr, id)
    else
        obj = object
    end
    push!(pool.pool, obj)
    return id #new momentum
end

function append(pool::Pool, object, evaluate::Function, isCached)
    append(pool, object, evaluate(object), isCached)
end

function append(pool::PoolwithParameter, object, evaluate::Function, isCached)
    append(pool, object, evaluate(pool.para, object), isCached)
end

function append(pool::CachedPool, object)
    # @assert para isa eltype(pool.pool)
    for (oi, o) in enumerate(pool.object)
        if o == object
            return oi #existing obj
        end
    end

    id = length(pool.object) + 1
    push!(pool.object, object)
    return id #new momentum
end

function initialize(pool::CachedPool; eval::Function = nothing, value = 0)
    N = length(pool.object)
    pool.version = zeros(Int128, 0)
    pool.excited = zeros(Bool, 0)
    T = eltype(fieldtype(pool, :current))
    if isnothing(eval)
        pool.current = [T(value) for i in 1:N]
    else
        pool.current = [T(eval(obj)) for obj in pool.object]
    end
    pool.new = deepcopy(pool.current)
end

"""
    struct LoopPool{O}

        Pool of (cached) objects.

# Members
- pool::Vector{O} : Vector that hosts the (cached) object
"""
mutable struct LoopPool{T}
    name::Symbol
    dim::Int #dimension
    N::Int #number of basis
    basis::Matrix{T}
    current::Matrix{T}

    function LoopPool(name::Symbol, dim::Int, N::Int, type::DataType)
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

function initialize(pool::LoopPool, variable::AbstractVector = rand(eltype(fieldtype(pool, :current)), pool.N))
    @assert length(variable) == pool.N
    T = eltype(fieldtype(pool, :current))
    pool.current[:, 1:length(pool)] = T.(variable) * pool.basis[:, 1:length(pool)]
end

function append(pool::LoopPool, basis::AbstractVector)
    for bi in 1:length(pool)
        if pool.basis[:, bi] ≈ basis
            return bi
        end
    end

    pool.basis = hcat(basis, pool.basis)
    # pool.size += 1
    # @assert pool.size <= size(pool.basis)[2] "Too many loop basis!. Increase maxSize when creates the LoopPool!"
    # pool.basis[:, pool.size] = basis
    return length(pool)
end
