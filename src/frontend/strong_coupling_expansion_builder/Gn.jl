function fullGreen(para, hop::Vector{BareHoppingId}; name=Symbol("Gn$(length(hop)*2)"), resetuid=false, even=true)
    extT, orbital, site, _creation = [], [], [], []
    for h in hop
        append!(extT, h.extT)
        append!(site, h.site)
        append!(_creation, [false, true])  # for GreenN; [true, false] for hopping.
        append!(orbital, h.orbital)
    end
    # println("calculate: ", hop, " . site: ", site)

    resetuid && IR.uidreset()

    GnId = GreenNId(para, orbital=orbital, t=extT, r=site, creation=_creation)
    gn = [Graph([], properties=GnId, name=Symbol("gn$(length(extT))"))]

    for h in hop
        push!(gn, Graph([], properties=h, name=:hop))
    end
    # println(permutation)
    # return Graph(gn, properties=GreenNId(para, orbital=orbital, t=extT, r=site, creation=_creation),
    return Graph(gn, operator=Prod(), name=name)
end

function vacuum_order(para, hop::Vector{BareHoppingId}; name=Symbol("vac_o$(length(hop)*2)"),
    prefactor=1.0, resetuid=false, even=true)
    extT, orbital, site, _creation = [], [], [], []
    for h in hop
        append!(extT, h.extT)
        append!(site, h.site)
        append!(_creation, [false, true])  # for GreenN; [true, false] for hopping.
        append!(orbital, h.orbital)
    end
    # println("calculate: ", hop, " . site: ", site)

    resetuid && IR.uidreset()

    GnId = GreenNId(para, orbital=orbital, t=extT, r=site, creation=_creation)
    gn = [Graph([], properties=GnId, name=Symbol("gn$(length(extT))"))]

    for h in hop
        push!(gn, Graph([], properties=h, name=:hop))
    end
    # println(permutation)
    # return Graph(gn, properties=GreenNId(para, orbital=orbital, t=extT, r=site, creation=_creation),
    return Graph(gn, operator=Prod(), name=name, factor=prefactor)
end

function fullGreen_v1(para, hop::Vector{BareHoppingId}; external_T::Vector{Int}=[], external_sites::Vector{Int}=[],
    external_creation::Vector{Bool}=[], external_orbitals::Vector{Int}=[], prefactor=1.0,
    name=Symbol("Gn$(length(hop)*2)"), resetuid=false, even=true
)

    @assert length(external_T) == length(external_sites) == length(external_creation) == length(external_orbitals)

    extT, orbital, site, _creation = [], [], [], []
    gn = Graph[]
    for h in hop
        append!(extT, h.extT)
        append!(site, h.site)
        append!(_creation, [false, true])  # for GreenN; [true, false] for hopping.
        append!(orbital, h.orbital)
        push!(gn, Graph([], properties=h, name=:hop))
    end

    append!(extT, external_T)
    append!(site, external_sites)
    append!(_creation, external_creation)
    append!(orbital, external_orbitals)
    # println("calculate: ", hop, " . site: ", site)

    resetuid && IR.uidreset()

    GnId = GreenNId(para, orbital=orbital, t=extT, r=site, creation=_creation)
    push!(gn, Graph([], properties=GnId, name=Symbol("gn$(length(extT))")))

    return Graph(gn, operator=Prod(), name=name, factor=prefactor)
end