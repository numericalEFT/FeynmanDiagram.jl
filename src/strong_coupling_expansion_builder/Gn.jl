function fullGreen(para, site::AbstractVector, orbital::AbstractVector, extT::AbstractVector = collect(1:length(orbital)), subdiagram = false; name = Symbol("Gn$(length(site))"), resetuid = false, even = true)
    # @assert para.diagType == GreenNDiag
    @assert length(extT) == length(orbital) == length(site)
    if even
        @assert length(extT) % 2 == 0
    end

    resetuid && uidreset()

    gn = []
    uniqueR = Set(site)
    permutation = [] # keep track of the permutation after the site index rearrangement
    for r in uniqueR
        ind = findall(x -> x == r, site)
        if even && (length(ind) % 2 == 1)
            continue
        end
        t = extT[ind]
        o = orbital[ind]
        bareGId = BareGreenNId(para, o, t, r)
        push!(gn, Diagram(bareGId, name = Symbol("gn$(length(t))")))
        append!(permutation, ind)
    end

    if isempty(gn)
        return nothing
    else
        return Diagram(GreenNId(para, orbital, extT, site), Prod(), gn, name = name, factor = parity(permutation))
    end
end