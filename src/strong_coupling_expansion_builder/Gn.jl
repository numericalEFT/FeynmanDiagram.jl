# function fullGreen(para, site::AbstractVector, orbital::AbstractVector, extT::AbstractVector = collect(1:length(orbital)), subdiagram = false; name = Symbol("Gn$(length(site))"), resetuid = false, even = true)
#     # @assert para.type == GreenNDiag
#     @assert length(extT) == length(orbital) == length(site)
#     if even
#         @assert length(extT) % 2 == 0
#     end

#     resetuid && uidreset()

#     gn = []
#     uniqueR = Set(site)
#     permutation = [] # keep track of the permutation after the site index rearrangement
#     for r in uniqueR
#         ind = findall(x -> x == r, site)
#         if even && (length(ind) % 2 == 1)
#             continue
#         end
#         t = extT[ind]
#         o = orbital[ind]
#         bareGId = BareGreenNId(para, o, t, r)
#         push!(gn, Diagram(bareGId, name = Symbol("gn$(length(t))")))
#         append!(permutation, ind)
#     end

#     if isempty(gn)
#         return nothing
#     else
#         return Diagram(GreenNId(para, orbital, extT, site), Prod(), gn, name = name, factor = parity(permutation))
#     end
# end

function fullGreen(para, hop::Vector{BareHoppingId}, subdiagram=false; name=Symbol("Gn$(length(hop)*2)"), resetuid=false, even=true)
    # @assert para.type == GreenNDiag
    # @assert length(extT) == length(orbital) == length(site)
    # if even
    #     @assert length(extT) % 2 == 0
    # end
    extT, orbital, site, _creation = [], [], [], []
    for h in hop
        append!(extT, h.extT)
        append!(site, h.site)
        append!(_creation, [true, false])
        append!(orbital, h.orbital)
    end
    # println("calculate: ", hop, " . site: ", site)

    resetuid && uidreset()

    gn = []
    uniqueR = Set(site)
    permutation = [] # keep track of the permutation after the site index rearrangement
    for r in uniqueR
        ind = findall(x -> x == r, site)
        if even && (length(ind) % 2 == 1)
            return nothing
        end
        t = extT[ind]
        o = orbital[ind]
        c = _creation[ind]
        bareGId = BareGreenNId(para, orbital=o, t=t, r=r, creation=c)
        push!(gn, Diagram(bareGId, name=Symbol("gn$(length(t))")))
        append!(permutation, ind)
    end

    if isempty(gn)
        return nothing
    else
        for h in hop
            push!(gn, Diagram(h, name=:hop))
        end
        # println(permutation)
        return Diagram(GreenNId(para, orbital=orbital, t=extT, r=site, creation=_creation), Prod(), gn, name=name, factor=parity(permutation))
    end
end