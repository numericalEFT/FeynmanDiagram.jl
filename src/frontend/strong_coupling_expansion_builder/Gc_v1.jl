
function connectedVacuum(para, hop::Vector{BareHoppingId}; name=Symbol("vacuum_c$(length(hop)*2)"), resetuid=false, even=true)
    N = length(hop)
    all_compositions = find_compositions(N)

    vac_graph = Graph[]
    for compos in all_compositions
        num_compos = sum(compos) - 1
        factor = (-1)^num_compos * factorial(num_compos) / prod(factorial.(compos))
        # println("Composition: $compos", " Factor: $factor")
        idx_hop = 1
        vac_vec = Graph[]
        for (o, num) in enumerate(compos)
            num == 0 && continue
            gvec = Graph[]
            for _ in 1:num
                # println(idx_hop, " ", idx_hop + o - 1)
                push!(gvec, vacuum_order(para, hop[idx_hop:idx_hop+o-1]; prefactor=(-1)^o / factorial(o)))
                idx_hop += o
            end
            push!(vac_vec, Graph(gvec, operator=Prod()))
        end
        push!(vac_graph, Graph(vac_vec, operator=Prod(), factor=factor))
    end

    return Graph(vac_graph, properties=VacuumId(para), operator=Sum(), name=name)
end

@inline function find_compositions(n::Int)
    if n <= 0
        println("Input must be a positive integer.")
        return Vector{Vector{Int}}[]
    end

    @inline function partition_to_coeffs(partition::Vector{Int})
        if isempty(partition)
            return Int[]
        end
        max_part = maximum(partition)
        coeffs = zeros(Int, max_part)

        for part in partition
            coeffs[part] += 1
        end
        return coeffs
    end

    all_compositions = Vector{Int}[]
    for p in partitions(n)
        # Convert the partition format (e.g., [4, 1, 1]) to the
        # desired coefficient format (e.g., b_1=2, b_4=1 => [2, 0, 0, 1])
        coeffs = partition_to_coeffs(p)
        if coeffs[1] == 0
            push!(all_compositions, coeffs)
        end
    end
    return all_compositions
end

# function connectedGreen_v1(para, hop::Vector{BareHoppingId}; external_T::Vector{Int}=[], external_sites::Vector{Int}=[],
#     external_creation::Vector{Bool}=[], external_orbitals::Vector{Int}=[],
#     prefactor=1.0, name=Symbol("Gc$(length(hop)*2)"), resetuid=false, even=true
# )

#     @assert length(external_T) == length(external_sites) == length(external_creation) == length(external_orbitals)

#     N = length(hop)

#     extT, orbital, site, _creation = [], [], [], []
#     resetuid && IR.uidreset()

#     Gc = [fullGreen_v1(para, hop; external_T=external_T, external_sites=external_sites,
#         external_creation=external_creation, external_orbitals=external_orbitals, resetuid=false, even=even)]

#     for (lind, rind) in partitions(collect(1:N), 2)
#         #this partition will not generate fermionic sign because the hopping term is always a bosonic operator
#         # if even && (length(lind) % 2 == 1)
#         #     continue
#         # end
#         subGc = connectedGreen_v1(para, hop[lind]; resetuid=false, even=even)
#         subGn = fullGreen_v1(para, hop[rind]; resetuid=false, even=even)
#         if isnothing(subGn) || isnothing(subGc)
#             continue
#         end
#         push!(Gc, Graph([subGc, subGn], properties=GenericId(para), operator=Prod())) #additional minus sign because Gc(s) = Gn(s) - \sum_o Gc(o)Gn(s-o)
#     end
# end

function connectedGreen(para, hop::Vector{BareHoppingId}; prefactor=1.0, name=Symbol("Gc$(length(hop)*2)"), resetuid=false, even=true)
    # @assert para.type == GreenNDiag
    # @assert length(extT) == length(orbital) == length(site)
    N = length(hop)

    extT, orbital, site, _creation = [], [], [], []
    for h in hop
        append!(extT, h.extT)
        append!(site, h.site)
        append!(_creation, [true, false])
        append!(orbital, h.orbital)
    end

    resetuid && IR.uidreset()

    Gc = [fullGreen(para, hop; resetuid=false, even=even)]

    for (lind, rind) in partitions(collect(1:N), 2)
        #this partition will not generate fermionic sign because the hopping term is always a bosonic operator
        # if even && (length(lind) % 2 == 1)
        #     continue
        # end
        subGc = connectedGreen(para, hop[lind]; resetuid=false, even=even)
        subGn = fullGreen(para, hop[rind]; resetuid=false, even=even)
        if isnothing(subGn) || isnothing(subGc)
            continue
        end
        push!(Gc, Graph([subGc, subGn], properties=GenericId(para), operator=Prod())) #additional minus sign because Gc(s) = Gn(s) - \sum_o Gc(o)Gn(s-o)
    end

    # extT, orbital, site, creation = [], [], [], []
    ext_T, ext_orbital, ext_site, ext_creation = [], [], [], []
    # for h in hop
    # 	append!(extT, h.extT)
    # 	append!(site, h.site)
    # 	append!(creation, [true, false])
    # 	append!(orbital, h.orbital)
    # end
    if isempty(ext_site)
        property = VacuumId(para)
    else
        property = ConnectedGreenNId(para, orbital=ext_orbital, t=ext_T, r=ext_site, creation=ext_creation)
    end

    sg_factors = ones(length(Gc))
    sg_factors[2:end] .= -1.0
    return Graph(Gc, subgraph_factors=sg_factors, properties=property, operator=Sum(), name=name, factor=prefactor)
end

function connectedGreen_o2(para, hop::Vector{BareHoppingId}; prefactor=1.0, name=Symbol("Gc$(length(hop)*2)"), resetuid=false, even=true)
    # @assert para.type == GreenNDiag
    # @assert length(extT) == length(orbital) == length(site)
    N = length(hop)

    extT, orbital, site, _creation = [], [], [], []
    gn = Graph[]
    for h in hop
        append!(extT, h.extT)
        append!(site, h.site)
        append!(_creation, [false, true])
        append!(orbital, h.orbital)
        push!(gn, Graph([], properties=h, name=:hop))
    end

    gvec = Graph[]
    labels = [[1, 4], [3, 2]]
    for idx in labels
        GnId = GreenNId(para, orbital=orbital[idx], t=extT[idx], r=site[idx], creation=_creation[idx])
        push!(gvec, Graph([], properties=GnId, name=Symbol("gn2")))
    end

    push!(gn, Graph(gvec, operator=Prod()))

    return Graph(gn, properties=VacuumId(para), operator=Prod(), name=name, factor=prefactor * -1.0)
end

function connectedGreen_o3(para, hop::Vector{BareHoppingId}; prefactor=1.0, name=Symbol("Gc$(length(hop)*2)"), resetuid=false, even=true)
    # @assert para.type == GreenNDiag
    # @assert length(extT) == length(orbital) == length(site)
    N = length(hop)

    println(hop)

    extT, orbital, site, _creation = [], [], [], []
    gn = Graph[]
    for h in hop
        append!(extT, h.extT)
        append!(site, h.site)
        append!(_creation, [false, true])
        append!(orbital, h.orbital)
        push!(gn, Graph([], properties=h, name=:hop))
    end

    gvec = Graph[]
    labels = [[1, 4], [5, 2], [3, 6]]
    for idx in labels
        GnId = GreenNId(para, orbital=orbital[idx], t=extT[idx], r=site[idx], creation=_creation[idx])
        push!(gvec, Graph([], properties=GnId, name=Symbol("gn2")))
    end

    gn_o2 = [Graph(gvec, operator=Prod())]

    gvec = Graph[]
    labels = [[1, 6], [3, 2], [5, 4]]
    for idx in labels
        GnId = GreenNId(para, orbital=orbital[idx], t=extT[idx], r=site[idx], creation=_creation[idx])
        push!(gvec, Graph([], properties=GnId, name=Symbol("gn2")))
    end
    push!(gn_o2, Graph(gvec, operator=Prod()))

    push!(gn, Graph(gn_o2, operator=Sum()))
    # push!(gn, Graph(gvec, operator=Prod()))

    return Graph(gn, properties=VacuumId(para), operator=Prod(), name=name, factor=prefactor)
end

# function connectedGreen(para, operators::Vector{OperatorId}, hop_id = [[2i-1, 2i] for i in Int(length(operators)/2)];
# 	name = Symbol("Gc$(length(hop)*2)"), resetuid = false, even = true)
# 	# @assert para.type == GreenNDiag
# 	# @assert length(extT) == length(orbital) == length(site)

# 	extT, orbital, site, _creation = [], [], [], []
# 	for op in operators
# 		push!(extT, op.extT)
# 		push!(orbital, op.orbital)
# 		push!(site, op.site)
# 		push!(_creation, op.creation)
# 	end
# 	N = length(site)

# 	resetuid && IR.uidreset()

# 	Gc = []
# 	Gfull = fullGreen(para, operators; resetuid = false, even = even)
# 	if isnothing(Gfull)
# 		return nothing
# 	end
# 	push!(Gc, Gfull)

# 	hop_set = Set([Set(h) for h in hop_id])
# 	for (lind, rind) in partitions(collect(1:N), 2)
# 		#this partition will not generate fermionic sign because the hopping term is always a bosonic operator
# 		# if even && (length(lind) % 2 == 1)
# 		#     continue
# 		# end
# 		for (l, r) in Iterators.product(lind, rind)
# 			# the splitted sites should not be connected by hopping
# 			Set([l, r]) in hop_set && continue
# 		end
# 		subGc = connectedGreen(para, operators[lind]; resetuid = false, even = even)
# 		subGn = fullGreen(para, operators[rind]; resetuid = false, even = even)
# 		if isnothing(subGn) || isnothing(subGc)
# 			continue
# 		end
# 		push!(Gc, Graph([subGc, subGn], properties = GenericId(para), operator = Prod(), factor = -1.0)) #additional minus sign because Gc(s) = Gn(s) - \sum_o Gc(o)Gn(s-o)
# 	end

# 	# extT, orbital, site, creation = [], [], [], []
# 	ext_T, ext_orbital, ext_site, ext_creation = [], [], [], []
# 	if isempty(ext_site)
# 		property = VacuumId(para)
# 	else
# 		property = ConnectedGreenNId(para, orbital = ext_orbital, t = ext_T, r = ext_site, creation = ext_creation)
# 	end
# 	return Graph(Gc, properties = property, operator = Sum(), name = name)
# end
