function readDiag(io::IO)

    return Diagram
end

function printDiag(diag::Union{Tuple,AbstractVector}, kwargs...)

    f = open("./DiagFiles/Diag.txt", "w")
    writedlm(f, [])
end
