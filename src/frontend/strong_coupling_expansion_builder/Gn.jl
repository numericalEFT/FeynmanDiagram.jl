# function derangement(n::Integer)
#     if n < 0
#         throw(DomainError(n, "n must be non-negative."))
#     end
#     if n == 0
#         return BigInt(1)
#     end

#     d_prev = BigInt(1)
#     d_current = BigInt(0)

#     for i = 1:n
#         # D(i) = i * D(i-1) + (-1)^i
#         sign = iseven(i) ? 1 : -1
#         d_current = i * d_prev + sign
#         d_prev = d_current
#     end

#     return d_current
# end

function fullGreen(para, hop::Vector{BareHoppingId}; name=Symbol("Gn$(length(hop)*2)"),
    resetuid=false, even=true, has_root=false) #    prefactor=1.0,
    extT, orbital, site, _creation = [], [], [], []
    N = length(hop)
    for h in hop
        append!(extT, h.extT)
        append!(site, h.site)
        append!(_creation, [false, true])  # for GreenN; [true, false] for hopping.
        append!(orbital, h.orbital)
    end

    resetuid && IR.uidreset()

    GnId = GreenNId(para, orbital=orbital, t=extT, r=site, creation=_creation)
    gn = [Graph([], properties=GnId, name=Symbol("gn$(length(extT))"))]

    for h in hop
        push!(gn, Graph([], properties=h, name=:hop))
    end

    return Graph(gn, operator=Prod(), name=name)
end

function vacuum_order(para, hop::Vector{<:BareHoppingId}; name=Symbol("vac_o$(length(hop)*2)"),
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

function fullGreen_external(para, hop::Vector{BareHoppingId}; external_T::Vector{Int}=[], external_sites::Vector{Int}=[],
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

function fullGreen_topology(para, hop::Vector{BareHoppingId};
    name=Symbol("Gn$(length(hop)*2)"), is_local_Gn::Bool=true,
    resetuid=false, even=true, has_root=false)
    extT, orbital, site, _creation = [], [], Int[], []
    N = length(hop)
    for h in hop
        append!(extT, h.extT)
        append!(site, h.site)
        append!(_creation, [false, true])  # for GreenN; [true, false] for hopping.
        append!(orbital, h.orbital)
    end

    resetuid && IR.uidreset()

    if !is_local_Gn
        GnId = GreenNId(para, orbital=orbital, t=extT, r=site, creation=_creation)
        gn = [Graph([], properties=GnId, name=Symbol("gn$(length(extT))"))]
        for h in hop
            push!(gn, Graph([], properties=h, name=:hop))
        end
        return Graph(gn, operator=Prod(), name=name)
    else
        pos_map = Dict{Int,Vector{Int}}()
        for (i, x) in enumerate(site)
            push!(get!(pos_map, x, Int[]), i)
        end

        gn = Graph[]
        for (x, idx) in pos_map
            GnId = BareGreenNId(para, r=x, orbital=orbital[idx], t=extT[idx], creation=_creation[idx])
            push!(gn, Graph([], properties=GnId, name=Symbol("gn$(length(idx))")))
        end
        for h in hop
            push!(gn, Graph([], properties=h, name=:hop))
        end

        sign = permu_sign(site)
        return Graph(gn, operator=Prod(), name=name, factor=sign)

        # if is_connected(site) 
        # vec_gn = Graph[]
        # inds_conn = get_connected_component_indices(site)
        # for idx in inds_conn
        #     o = orbital[idx]
        #     t = extT[idx]
        #     r = site[idx]
        #     c = _creation[idx]
        #     GnId = GreenNId(para, orbital=o, t=t, r=r, creation=c)
        #     gn = [Graph([], properties=GnId, name=Symbol("gn$(length(t))"))]
        #     for h in hop[idx]
        #         push!(gn, Graph([], properties=h, name=:hop))
        #     end
        #     push!(vec_gn, Graph(gn, operator=Prod(), name=name))
        # end
        # return Graph(vec_gn, operator=Prod(), name=name)
    end
end

@inline function permu_sign(v::Vector{Int})
    # swap to GreenN c*cc*c... operator sequence
    n = length(v)
    for i in 1:2:min(n - 1, n)
        if i + 1 <= n
            v[i], v[i+1] = v[i+1], v[i]
        end
    end

    sign = 1
    for i in eachindex(v)
        for j in (i+1):length(v)
            if (v[i] > v[j]) || (v[i] == v[j] && iseven(i) && isodd(j))
                sign *= -1
            end
        end
    end
    return sign
end