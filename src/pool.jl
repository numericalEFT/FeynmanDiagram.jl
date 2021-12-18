"""
    mutable struct CachedObject{PARA,T}

        CachedOjects are too heavy so that one wants to cache their status without evaluating them on the fly.
        The user should defines a compare 

    para::PARA : para needs to evaluate the object
    id::Int : index of the object
    curr::T : current status
    new::T : the new status wants to assign later
    version::Int128 : the current version
    excited::Bool : if set to excited, then the current status needs to be replaced with the new status
"""
mutable struct CachedObject{PARA,T}
    para::PARA
    id::Int
    curr::T
    new::T
    version::Int128
    excited::Bool #if set to excited, then the current weight needs to be replaced with the new weight
    function CachedObject(para::P, curr::T, id, version = 1, excited = false) where {P,T}
        return new{P,T}(para, id, curr, curr, version, excited)
    end
end

struct Pool{O<:CachedObject}
    pool::Vector{O}

    function Pool{O}()
        return new{O}([])
    end
    function Pool(obj::Vector{O})
        return new{O}(obj)
    end
end

Base.length(pool::Pool) = length(pool.pool)
Base.size(pool::Pool) = size(pool.pool)

# function append(pool::Pool, obj::CachedObject)
#     for (oi, o) in enumerate(pool)
#         if o.para == obj.para
#             return oi, false #existing obj
#         end
#     end
#     push!(pool, obj)
#     return length(pool), true #new momentum
# end

function append(pool::Pool, para, curr)
    @assert para isa eltype(pool.pool)

    for (oi, o) in enumerate(pool)
        if o.para == para
            return oi, false #existing obj
        end
    end
    id = length(pool)
    push!(pool, CachedObject(para, curr, id))
    return id, true #new momentum
end

struct Collection{O<:CachedObject}
    pool::Pool{O}
    idx::Vector{Int}
    function Collection(pool::Pool{O}, idx = []) where {O<:CachedObject}
        return new(pool, idx)
    end
end
