
# function connectedVacuum(para, hop::Vector{BareHoppingId{T}}; name=Symbol("vacuum_c$(length(hop)*2)"), resetuid=false, even=true) where {T}
function connectedVacuum(para, hop::Vector{<:BareHoppingId}; name=Symbol("vacuum_c$(length(hop)*2)"), resetuid=false, even=true)
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
        # if coeffs[1] == 0  # no order-1 contributions
        push!(all_compositions, coeffs)
        # end
    end
    return all_compositions
end

function connectedGreen(para, hop::Vector{BareHoppingId}; name=Symbol("Gc$(length(hop)*2)"),
    resetuid=false, even=true, prefactor=1.0)
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

    Gc = [fullGreen(para, hop; resetuid=false, even=even, has_root=true)]

    for (lind, rind) in partitions(collect(1:N), 2)
        subGc = connectedGreen(para, hop[lind]; resetuid=false, even=even, name=Symbol("Gc$(length(lind)*2)"))
        subGn = fullGreen(para, hop[rind]; resetuid=false, even=even)
        if isnothing(subGn) || isnothing(subGc)
            continue
        end
        push!(Gc, Graph([subGc, subGn], properties=GenericId(para), operator=Prod())) #additional minus sign because Gc(s) = Gn(s) - \sum_o Gc(o)Gn(s-o)
    end

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

    return Graph(Gc, factor=prefactor, subgraph_factors=sg_factors, properties=property, operator=Sum(), name=name)

end


function connectedGreen!(para, hop_sets::Vector{Vector{BareHoppingId}}, Gc_pool=Dict(), Gn_pool=Dict();
    name=Symbol("Gc$(length(hop)*2)"), is_local_Gn::Bool=true,
    resetuid=false, even=true, prefactors=[1.0 for _ in 1:length(hop_sets)])
    @assert allequal(length.(hop_sets))
    N = length(hop_sets[1])
    graphs = Graph[]

    # println("$N hoppings ", prefactors)

    for hop in hop_sets
        site = Tuple{Int,Int}[]
        for h in hop
            push!(site, h.site)
        end

        if haskey(Gc_pool, hop)
            push!(graphs, Gc_pool[hop])
            continue
        end

        resetuid && IR.uidreset()

        if haskey(Gn_pool, hop)
            Gc = [Gn_pool[hop]]
        else
            Gc = [fullGreen_topology(para, hop; is_local_Gn=is_local_Gn, resetuid=false, even=even, has_root=true)]
        end

        # println(site)
        for (lind, rind) in partitions(collect(1:N), 2)
            if !(is_connected(site[lind]) && is_closed(site[rind]))
                continue
            end
            # println(site[lind], "  ", site[rind])
            if haskey(Gc_pool, hop[lind])
                subGc = Gc_pool[hop[lind]]
            else
                # site_unique = unique(collect(Iterators.flatten(site[lind])))
                # sym_factor = get_symmetry_factor(length(site_unique), site[lind])
                subGc, _, _ = connectedGreen!(para, [hop[lind]], Gc_pool, Gn_pool;
                    is_local_Gn=is_local_Gn, resetuid=false, even=even, name=Symbol("Gc$(length(lind)*2)"))
                # prefactors=[1.0 / sym_factor,], resetuid=false, even=even, name=Symbol("Gc$(length(lind)*2)"))
                Gc_pool[hop[lind]] = subGc
            end
            if haskey(Gn_pool, hop[rind])
                subGn = Gn_pool[hop[rind]]
            else
                subGn = fullGreen_topology(para, hop[rind]; is_local_Gn=is_local_Gn, resetuid=false, even=even)
                Gn_pool[hop[rind]] = subGn
            end
            push!(Gc, Graph([subGc, subGn], properties=GenericId(para), operator=Prod())) #additional minus sign because Gc(s) = Gn(s) - \sum_o Gc(o)Gn(s-o)
        end
        # ext_T, ext_orbital, ext_site, ext_creation = [], [], [], []
        # if isempty(ext_site)
        property = VacuumId(para)
        # else
        #     property = ConnectedGreenNId(para, orbital=ext_orbital, t=ext_T, r=ext_site, creation=ext_creation)
        # end
        sg_factors = ones(length(Gc))
        sg_factors[2:end] .= -1.0
        push!(graphs, Graph(Gc, subgraph_factors=sg_factors, properties=property, operator=Sum(), name=name))
    end

    total_Gc = Graph(graphs, subgraph_factors=prefactors, operator=Sum(), properties=VacuumId(para), name=name)
    return total_Gc, Gc_pool, Gn_pool
end