# function fullGreen(para, operators::Vector{OperatorId}; name = Symbol("Gn$(length(hop)*2)"), resetuid = false, even = true)
# 	extT, orbital, site, _creation = [], [], [], []
# 	for op in operators
# 		push!(extT, op.extT)
# 		push!(orbital, op.orbital)
# 		push!(site, op.site)
# 		push!(_creation, op.creation)
# 	end
# 	# println("calculate: ", hop, " . site: ", site)

# 	resetuid && IR.uidreset()

# 	gn = Graph[]
# 	uniqueR = Set(site)
# 	permutation = [] # keep track of the permutation after the site index rearrangement
# 	for r in uniqueR
# 		ind = findall(x -> x == r, site)
# 		if even && (length(ind) % 2 == 1)
# 			return nothing
# 		end
# 		t = extT[ind]
# 		o = orbital[ind]
# 		c = _creation[ind]

# 		bareGId = BareGreenNId(para, orbital = o, t = t, r = r, creation = c)
# 		push!(gn, Graph([], properties = bareGId, name = Symbol("gn$(length(t))")))
# 		append!(permutation, ind)
# 	end

# 	# for h in hop
# 	# 	push!(gn, Graph([], properties = h, name = :hop))
# 	# end
# 	# println(permutation)
# 	return Graph(gn, properties = GreenNId(para, orbital = orbital, t = extT, r = site, creation = _creation), operator = Prod(), name = name, factor = parity(permutation))
# end


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