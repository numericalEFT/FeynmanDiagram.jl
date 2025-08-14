function connectedGreen(para, site::Vector{Int}, orbital::AbstractVector, extT::AbstractVector, creation::AbstractVector;
    # function connectedGreen(para, site, orbital, extT, creation;
    ext_site::Vector{Int}=Int[], ext_orbital::Vector{Int}=Int[], ext_T::Vector{Int}=Int[], ext_creation::Vector{Bool}=Bool[],
    name=Symbol("Gc$(length(site))"), resetuid=false, num_orbital::Int=2)

    @assert length(extT) == length(orbital) == length(site) == length(creation)
    @assert isdisjoint(ext_site, site)
    @assert length(ext_site) == length(ext_orbital) == length(ext_T) == length(ext_creation)

    # N = length(site)
    resetuid && IR.uidreset()
    Gc = []

    Gfull = fullGreen(para, site, orbital, extT, creation;
        ext_site=ext_site, ext_orbital=ext_orbital, ext_T=ext_T, ext_creation=ext_creation, resetuid=false, num_orbital=num_orbital)
    push!(Gc, Gfull)

    uniqueR = unique(site)
    N = length(uniqueR)
    # for mask in 1:(2^N-2)
    #     S_idx = Int[]
    #     for b in 1:N
    #         (mask & (1 << (b - 1))) != 0 && push!(S_idx, b)
    #     end
    #     lidx = findall(x -> x in uniqueR[S_idx], site)
    #     ridx = findall(x -> x ∉ uniqueR[S_idx], site)
    for (lind, rind) in partitions(collect(1:N), 2)
        lidx = findall(x -> x in uniqueR[lind], site)
        ridx = findall(x -> x in uniqueR[rind], site)
        subGc = connectedGreen(para, site[lidx], orbital[lidx], extT[lidx], creation[lidx];
            ext_site=ext_site, ext_orbital=ext_orbital, ext_T=ext_T, ext_creation=ext_creation, resetuid=false)
        subGn = fullGreen(para, site[ridx], orbital[ridx], extT[ridx], creation[ridx]; resetuid=false, num_orbital=num_orbital)

        push!(Gc, Graph([subGc, subGn], properties=GenericId(para), operator=Prod(), factor=-1.0))
    end

    if isempty(ext_site)
        property = VacuumId(para)
    else
        property = ConnectedGreenNId(para, orbital=ext_orbital, t=ext_T, r=ext_site, creation=ext_creation)
    end

    num_repeated_vertextype = count_repeated_occurrences(site)
    factor = 1.0 / prod(factorial.(num_repeated_vertextype))

    return Graph(Gc, properties=property, factor=factor, operator=Sum(), name=name)
end

@inline function count_repeated_occurrences(arr::Vector{Int})
    counts = Dict{Int,Int}()
    for num in arr
        counts[num] = get(counts, num, 0) + 1
    end

    occurrences = Dict{Int,Int}()
    for count in values(counts)
        occurrences[count] = get(occurrences, count, 0) + 1
    end

    return collect(values(occurrences))
end
