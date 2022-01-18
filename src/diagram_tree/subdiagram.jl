abstract type DiagramId end

toDict(d::DiagramId) = error("toDict not implemented!")
Base.show(io::IO, d::DiagramId) = error("Base.show not implemented!")
Base.isequal(a::DiagramId, b::DiagramId) = error("Base.isequal not implemented!")
Base.(==)(a::DiagramId, b::DiagramId) = Base.iequal(a, b)

mutable struct Diagram{T<:Diagram}
    id::T
    # parent::Diagram
    subdiagrams::Vector{Diagram}

    Diagram(id::T) where {T<:DiagramId} = new{T}(id, [])
end

function addsubdiagram!(parent::Diagram, child::Diagram)
    for c in parent.children
        if c.id == child.id
            return false
        end
    end
    push!(parent.children, child)
end

function toDataFrame(diagVec::Vector{Diagram})

end