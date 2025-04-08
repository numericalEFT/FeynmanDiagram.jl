function connectedGreen(para, site::Vector{Int}, orbital::AbstractVector, extT::AbstractVector, creation::AbstractVector;
	# function connectedGreen(para, site, orbital, extT, creation;
	ext_site::Vector{Int} = Int[], ext_orbital::Vector{Int} = Int[], ext_T::Vector{Int} = Int[], ext_creation::Vector{Bool} = Bool[],
	name = Symbol("Gc$(length(site))"), resetuid = false, num_orbital::Int = 2)

	@assert length(extT) == length(orbital) == length(site) == length(creation)
	@assert isdisjoint(ext_site, site)
	@assert length(ext_site) == length(ext_orbital) == length(ext_T) == length(ext_creation)

	# N = length(site)
	resetuid && IR.uidreset()
	Gc = []

	Gfull = fullGreen(para, site, orbital, extT, creation;
		ext_site = ext_site, ext_orbital = ext_orbital, ext_T = ext_T, ext_creation = ext_creation, resetuid = false, num_orbital = num_orbital)
	push!(Gc, Gfull)

	uniqueR = unique(site)
	N = length(uniqueR)
	for (lind, rind) in partitions(collect(1:N), 2)
		lidx = findall(x -> x in uniqueR[lind], site)
		ridx = findall(x -> x in uniqueR[rind], site)
		subGc = connectedGreen(para, site[lidx], orbital[lidx], extT[lidx], creation[lidx];
			ext_site = ext_site, ext_orbital = ext_orbital, ext_T = ext_T, ext_creation = ext_creation, resetuid = false)
		subGn = fullGreen(para, site[ridx], orbital[ridx], extT[ridx], creation[ridx]; resetuid = false, num_orbital = num_orbital)

		push!(Gc, Graph([subGc, subGn], properties = GenericId(para), operator = Prod(), factor = -1.0))
	end

	if isempty(ext_site)
		property = VacuumId(para)
	else
		property = ConnectedGreenNId(para, orbital = ext_orbital, t = ext_T, r = ext_site, creation = ext_creation)
	end
	return Graph(Gc, properties = property, operator = Sum(), name = name)
end

function connectedGreen(para, hop::Vector{BareHoppingId}, subdiagram = false; name = Symbol("Gc$(length(hop)*2)"), resetuid = false, even = true)
	# @assert para.type == GreenNDiag
	# @assert length(extT) == length(orbital) == length(site)
	# println("wip")
	# println(hop)
	N = length(hop)

	resetuid && IR.uidreset()

	Gc = []

	# for paired Green's function, odd number of legs always leads to zero 
	# if even && (length(site) % 2 == 1)
	#     return nothing
	# end

	Gfull = fullGreen(para, hop, true; resetuid = false, even = even)
	if isnothing(Gfull)
		return nothing
	end
	push!(Gc, Gfull)

	for (lind, rind) in partitions(collect(1:N), 2)
		#this partition will not generate fermionic sign because the hopping term is always a bosonic operator
		# if even && (length(lind) % 2 == 1)
		#     continue
		# end
		subGc = connectedGreen(para, hop[lind], true; resetuid = false, even = even)
		subGn = fullGreen(para, hop[rind], true; resetuid = false, even = even)
		if isnothing(subGn) || isnothing(subGc)
			continue
		end
		push!(Gc, Graph([subGc, subGn], properties = GenericId(para), operator = Prod(), factor = -1.0)) #additional minus sign because Gc(s) = Gn(s) - \sum_o Gc(o)Gn(s-o)
	end

	extT, orbital, site, creation = [], [], [], []
	for h in hop
		append!(extT, h.extT)
		append!(site, h.site)
		append!(creation, [true, false])
		append!(orbital, h.orbital)
	end
	return Graph(Gc, properties = ConnectedGreenNId(para, orbital = orbital, t = extT, r = site, creation = creation), operator = Sum(), name = name)
end

# function connectedGreen(para, site::AbstractVector, orbital::AbstractVector, extT::AbstractVector = collect(1:length(orbital)), subdiagram = false; name = Symbol("Gc$(length(site))"), resetuid = false, even = true)
#     # @assert para.type == GreenNDiag
#     @assert length(extT) == length(orbital) == length(site)
#     N = length(extT)

#     resetuid && uidreset()

#     Gc = []

#     # for paired Green's function, odd number of legs always leads to zero 
#     if even && (length(site) % 2 == 1)
#         return nothing
#     end

#     Gfull = fullGreen(para, site, orbital, extT, true; resetuid = false, even = even)
#     if isnothing(Gfull)
#         return nothing
#     end
#     push!(Gc, Gfull)

#     for (lind, rind) in partitions(collect(1:N), 2)
#         if even && (length(lind) % 2 == 1)
#             continue
#         end
#         lR, lo, lt = site[lind], orbital[lind], extT[lind]
#         rR, ro, rt = site[rind], orbital[rind], extT[rind]
#         subGc = connectedGreen(para, lR, lo, lt, true; resetuid = false, even = even)
#         subGn = fullGreen(para, rR, ro, rt, true; resetuid = false, even = even)
#         if isnothing(subGn) || isnothing(subGc)
#             continue
#         end
#         p = parity(vcat(lind, rind))
#         push!(Gc, Diagram(GenericId(para), Prod(), [subGc, subGn], factor = -p)) #additional minus sign because Gc(s) = Gn(s) - \sum_o Gc(o)Gn(s-o)
#     end

#     return Diagram(ConnectedGreenNId(para, orbital, extT, site), Sum(), Gc, name = name)
# end