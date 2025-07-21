function fullGreen(para, site::Vector{Int}, orbital::AbstractVector, extT::AbstractVector, creation::AbstractVector;
	# bareGreenN::Union{Function, Nothing} = nothing,
	bareGreenN::Union{Function, Nothing} = Gnc,
	ext_site::Vector{Int} = Int[], ext_orbital::Vector{Int} = Int[], ext_T::Vector{Int} = Int[], ext_creation::Vector{Bool} = Bool[],
	name = Symbol("Gn$(length(site))"), resetuid = false, num_orbital::Int = 2)

	@assert length(extT) == length(orbital) == length(site) == length(creation)
	@assert isdisjoint(ext_site, site)
	@assert length(ext_site) == length(ext_orbital) == length(ext_T) == length(ext_creation)

	resetuid && IR.uidreset()

	all_sites = isempty(ext_site) ? site : vcat(ext_site, site)
	len_ext = length(ext_site)

	Gn = [[] for _ in 1:num_orbital]
	for (i, ri) in enumerate(all_sites)
		if i <= len_ext && ext_creation[i]
			oi = ext_orbital[i]
			ti = ext_T[i]
		else
			c = creation[i-len_ext]
			oi = orbital[i-len_ext][c][1]
			ti = extT[i-len_ext][c][1]
		end

		for (j, rj) in enumerate(all_sites)
			if j <= len_ext && !ext_creation[j]
				oj = ext_orbital[j]
				tj = ext_T[j]
			else
				a = .!creation[j-len_ext]
				oj = orbital[j-len_ext][a][1]
				tj = extT[j-len_ext][a][1]
			end
			oi != oj && continue

			# hopping operator is conjugate to the half operator of a given site.
			push!(Gn[oi], Graph([], properties = BareHoppingId(para, (rj, ri), (oj, oi), (tj, ti))))
		end
	end

	dets_g = []
	for i in 1:num_orbital
		isempty(Gn[i]) && continue
		push!(dets_g, Graph(Gn[i], operator = Det(), name = :det))
	end

	gn = [Graph(dets_g, operator = Prod())]


	uniqueR = Set(site)
	for r in uniqueR
		ind = findall(x -> x == r, site)
		t = collect(Iterators.flatten(extT[ind]))
		o = collect(Iterators.flatten(orbital[ind]))
		c = collect(Iterators.flatten(creation[ind]))

		bareGId = BareGreenNId(para, orbital = o, t = t, r = r, creation = c)
		if isnothing(bareGreenN)
			bareGN = Graph([], properties = bareGId, name = Symbol("gn$(length(t))"), factor = prefactor(o[[2m - 1 for m in 1:length(ind)]], num_orbital))
		else
			bareGN = bareGreenN(bareGId) * prefactor(o[[2m - 1 for m in 1:length(ind)]], num_orbital)
		end
		push!(gn, bareGN)

		# println("t: $t, o: $o, c: $c, factor: $(prefactor(o[[2m - 1 for m in 1:length(ind)]], num_orbital))")
	end

	if isempty(ext_site)
		property = VacuumId(para)
	else
		property = GreenNId(para, orbital = ext_orbital, t = ext_T, r = ext_site, creation = ext_creation)
	end
	return Graph(gn, properties = property, operator = Prod(), name = name)
end

function Gnc(GId::DiagramId)
	N = length(GId.extT)
	if N == 2
		return Graph([], properties = GId, name = Symbol("bareGnc$N"))
	end

	order = Int(N / 2)
	@assert order * 2 == N "BareGreenN must be even N"

	o = vcat(GId.orbital[GId.creation], GId.orbital[.!GId.creation])
	t = vcat(GId.extT[GId.creation], GId.extT[.!GId.creation])
	GId = BareGreenNId(GId.para, orbital = o, t = t, r = GId.site, creation = vcat([true for _ in 1:order], [false for _ in 1:order]))

	G = Graph([], properties = GId, name = Symbol("bareGn$N"))
	Gc = [G]

	p = Partition(Int(N / 2))
	# println(p.l)
	for (li, l) in enumerate(p.l)
		sub_gncId = BareGreenNId(GId.para, orbital = GId.orbital[l], t = GId.extT[l], r = GId.site, creation = GId.creation[l])
		sub_gnc = Gnc(sub_gncId)

		r = p.r[li]
		sub_gnId = BareGreenNId(GId.para, orbital = GId.orbital[r], t = GId.extT[r], r = GId.site, creation = GId.creation[r])
		sub_gn = Graph([], properties = sub_gnId, name = Symbol("gn$(length(r))"))

		# push!(Gc, Graph([sub_gnc, sub_gn], properties = GenericId(GId.para), operator = Prod(), factor = -1.0 * p.sign[li]))
		push!(Gc, Graph([sub_gnc, sub_gn], properties = GenericId(GId.para), operator = Prod(), factor = p.sign[li]))
	end

	return Graph(Gc, properties = GId, operator = Sum(), name = Symbol("bareGnc$N"))
end

function prefactor(orbitals, num_orbital::Int)
	m = length(orbitals)
	_factor = 1.0
	for i in 1:m
		for j in (i+1):m
			_factor *= (num_orbital - (orbitals[i] == orbitals[j] ? 1 : 0))
		end
	end
	# return _factor / (2m)
	return _factor * (-1)^m / factorial(m)^2
	# return _factor / factorial(m)^2
	# return (-1)^m / factorial(m)^2
end

function fullGreen(para, hop::Vector{BareHoppingId}, subdiagram = false; name = Symbol("Gn$(length(hop)*2)"), resetuid = false, even = true)
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

	resetuid && IR.uidreset()

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
		bareGId = BareGreenNId(para, orbital = o, t = t, r = r, creation = c)
		push!(gn, Graph([], properties = bareGId, name = Symbol("gn$(length(t))")))
		append!(permutation, ind)
	end

	if isempty(gn)
		return nothing
	else
		for h in hop
			push!(gn, Graph([], properties = h, name = :hop))
		end
		# println(permutation)
		return Graph(gn, properties = GreenNId(para, orbital = orbital, t = extT, r = site, creation = _creation), operator = Prod(), name = name, factor = parity(permutation))
	end
end
