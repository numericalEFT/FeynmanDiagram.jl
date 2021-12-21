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

# Base.show(io::IO, obj::Cache) = print(io, "id $(obj.id): obj: $(obj.object) curr: $(obj.curr)")
Base.show(io::IO, obj::Cache) = print(io, "obj: $(obj.object) curr: $(obj.curr)")

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


# function append(pool::Pool, obj::CachedObject)
#     for (oi, o) in enumerate(pool)
#         if o.para == obj.para
#             return oi, false #existing obj
#         end
#     end
#     push!(pool, obj)
#     return length(pool), true #new momentum
# end

# function append(pool::Pool, para, curr)
#     @assert para isa eltype(pool.pool)

#     for (oi, o) in enumerate(pool)
#         if o.para == para
#             return oi, false #existing obj
#         end
#     end
#     id = length(pool)
#     push!(pool, CachedObject(para, curr, id))
#     return id, true #new momentum
# end

# struct SubPool{O,T}
#     pool::Pool{O,T}
#     idx::Vector{Int}
#     function SubPool(pool::Pool{O,T}, idx = []) where {O,T}
#         return new{O,T}(pool, idx)
#     end
# end
