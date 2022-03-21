function connectedGreen(para, site::AbstractVector, orbital::AbstractVector, extT::AbstractVector = collect(1:length(orbital)), subdiagram = false; name = Symbol("Gc$(length(site))"), resetuid = false, even = true)
    # @assert para.diagType == GreenNDiag
    @assert length(extT) == length(orbital) == length(site)
    N = length(extT)

    resetuid && uidreset()

    Gc = []

    # for paired Green's function, odd number of legs always leads to zero 
    if even && (length(site) % 2 == 1)
        return nothing
    end

    Gfull = fullGreen(para, site, orbital, extT, true; resetuid = false, even = even)
    if isnothing(Gfull)
        return nothing
    end
    push!(Gc, Gfull)

    for (lind, rind) in partitions(collect(1:N), 2)
        if even && (length(lind) % 2 == 1)
            continue
        end
        lR, lo, lt = site[lind], orbital[lind], extT[lind]
        rR, ro, rt = site[rind], orbital[rind], extT[rind]
        subGc = connectedGreen(para, lR, lo, lt, true; resetuid = false, even = even)
        subGn = fullGreen(para, rR, ro, rt, true; resetuid = false, even = even)
        if isnothing(subGn) || isnothing(subGc)
            continue
        end
        p = parity(vcat(lind, rind))
        push!(Gc, Diagram(GenericId(para), Prod(), [subGc, subGn], factor = -p)) #additional minus sign because Gc(s) = Gn(s) - \sum_o Gc(o)Gn(s-o)
    end

    return Diagram(ConnectedGreenNId(para, orbital, extT, site), Sum(), Gc, name = name)
end