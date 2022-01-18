abstract type DiagramId end

toDict(d::DiagramId) = error("toDict not implemented!")
Base.show(io::IO, d::DiagramId) = error("Base.show not implemented!")
Base.isequal(a::DiagramId, b::DiagramId) = error("Base.isequal not implemented!")
Base.(==)(a::DiagramId, b::DiagramId) = Base.iequal(a, b)

struct Diagram{T<:Diagram}
    id::T
    subdiagrams::Vector{Diagram}
    # parent::Diagram

    Diagram(id::T) where {T<:DiagramId} = new{T}(id, [])
end

function add_subdiagram!(parent::Diagram, child::Diagram)
    for c in parent.children
        if c.id == child.id
            return false
        end
    end
    push!(parent.children, child)
end

function toDataFrame(diagVec::Vector{Diagram})
    df = DataFrame()
    for diag in diagVec
        d = toDict(diag.id)
        d[:hash] = hash(diag.id)
        d[:subdiagram] = [hash(d) for d in diag.subdiagrams]
        push!(df, d)
    end
    return df
end