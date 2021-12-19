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

Base.show(io::IO, obj::CachedObject) = print(io, "id: $(obj.id), para: $(obj.para) curr: $(obj.curr)")

struct Pool{P,T}
    pool::Vector{CachedObject{P,T}}

    function Pool{P,T}() where {P,T}
        pool = Vector{CachedObject{P,T}}(undef, 0)
        return new{P,T}(pool)
    end
    function Pool(obj::Vector{CachedObject{P,T}}) where {P,T}
        return new{P,T}(obj)
    end
end

Base.length(pool::Pool) = length(pool.pool)
Base.size(pool::Pool) = size(pool.pool)
Base.show(io::IO, pool::Pool) = print(io, pool.pool)
Base.view(pool::Pool, inds...) where {N} = Base.view(pool.pool, inds...)

#index interface for Pool
Base.getindex(pool::Pool, i) = pool.pool[i]
Base.setindex!(pool::Pool, v, i) = setindex!(pool.pool, v, i)
Base.firstindex(pool::Pool) = 1
Base.lastindex(pool::Pool) = length(pool)

function append(pool, para::P, curr::T) where {P,T}
    # @assert para isa eltype(pool.pool)
    for (oi, o) in enumerate(pool.pool)
        if o.para == para
            return oi, false #existing obj
        end
    end
    id = length(pool)
    push!(pool.pool, CachedObject(para, curr, id))
    return id, true #new momentum
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

# struct SubPool{O}
#     pool::Pool{O}
#     idx::Vector{Int}
#     function Collection(pool::Pool{O}, idx = []) where {O<:CachedObject}
#         return new(pool, idx)
#     end
# end
