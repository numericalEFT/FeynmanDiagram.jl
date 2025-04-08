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

function fullGreen(para, site::Vector{Int}, orbital::AbstractVector, extT::AbstractVector, creation::AbstractVector;
	ext_site::Vector{Int} = Int[], ext_orbital::Vector{Int} = Int[], ext_T::Vector{Int} = Int[], ext_creation::Vector{Bool} = Bool[],
	name = Symbol("Gn$(length(site))"), resetuid = false, num_orbital::Int = 2)

	@assert length(extT) == length(orbital) == length(site) == length(creation)
	@assert isdisjoint(ext_site, site)
	@assert length(ext_site) == length(ext_orbital) == length(ext_T) == length(ext_creation)

	resetuid && IR.uidreset()

	all_sites = isempty(ext_site) ? site : vcat(ext_site, site)
	len_ext = length(ext_site)

	Gn = []
	for (i, ri) in enumerate(all_sites)
		if i <= len_ext
			oi = ext_orbital[i]
			ti = ext_T[i]
		else
			oi = orbital[i-len_ext][1]
			ti = extT[i-len_ext][1]
		end

		for (j, rj) in enumerate(all_sites)
			if j <= len_ext
				oj = ext_orbital[j]
				tj = ext_T[j]
			else
				oj = orbital[j-len_ext][2]
				tj = extT[j-len_ext][2]
			end
			# push!(Gn, Graph([], properties = BareHoppingId(para, (ri, rj), (oi, oj), (ti, tj))))
			push!(Gn, Graph([], properties = BareHoppingId(para, (rj, ri), (oj, oi), (tj, ti))))
		end
	end

	gn = [Graph(Gn, operator = Det(), name = :det)]

	uniqueR = Set(site)
	for r in uniqueR
		ind = findall(x -> x == r, site)
		t = collect(Iterators.flatten(extT[ind]))
		o = collect(Iterators.flatten(orbital[ind]))
		c = collect(Iterators.flatten(creation[ind]))
		bareGId = BareGreenNId(para, orbital = o, t = t, r = r, creation = c)
		push!(gn, Graph([], properties = bareGId, name = Symbol("gn$(length(t))"),
			factor = prefactor(o[[2m - 1 for m in 1:length(ind)]], num_orbital)))
	end

	if isempty(ext_site)
		property = VacuumId(para)
	else
		property = GreenNId(para, orbital = ext_orbital, t = ext_T, r = ext_site, creation = ext_creation)
	end
	return Graph(gn, properties = property, operator = Prod(), name = name)
end

function prefactor(orbitals, num_orbital::Int)
	m = length(orbitals)
	_factor = 1.0
	for i in 1:m
		for j in (i+1):m
			_factor *= (num_orbital - (orbitals[i] == orbitals[j] ? 1 : 0))
		end
	end
	return _factor / (2m)
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
