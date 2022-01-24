function compile(diag::Diagram{W}) where {W}
    exprtree = newExprTree(diag.id.para, diag.name)

end

function newExprTree(para, name::Symbol = :none)
    weightType = para.weightType
    Kpool = LoopPool(:K, para.loopDim, para.totalLoopNum, Float64)
    # nodeParaType = Vector{Int}
    nodeParaType = Any
    propagatorPool = []
    push!(propagatorPool, propagatorPool(:Gpool, weightType))
    for interaction in para.interaction
        response = interaction.response
        for type in interaction.type
            push!(propagatorPool, propagatorPool(symbol(response, type, "pool"), weightType))
        end
    end
    return Diagrams(Kpool, Tuple(propagatorPool), weightType, nodeParaType = nodeParaType, name = name)
end