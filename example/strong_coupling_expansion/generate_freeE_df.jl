using FeynmanDiagram
import FeynmanDiagram.Parquet: DiagPara, Interaction, VacuumDiag
import FeynmanDiagram.ComputationalGraphs: Sum
import FeynmanDiagram.FrontEnds: ConnectedGreenNId, BareHoppingId, VacuumId, UpUp, UpDown, Dynamic
using Parameters

function generate_vectors(order::Int)
    function helper(current_vector::Vector{Int}, remaining_length::Int)
        if remaining_length == 0
            return [copy(current_vector)]
        end

        results = Vector{Vector{Int}}()
        last_value = current_vector[end]

        # Option 1: Repeat the last value
        push!(results, helper(vcat(current_vector, [last_value]), remaining_length - 1)...)

        # Option 2: Increment the last value
        push!(results, helper(vcat(current_vector, [last_value + 1]), remaining_length - 1)...)

        return results
    end

    return helper([1], order - 1)
end

function generate_topologies(order::Int)
    all_vectors = generate_vectors(order)
    groups = Dict{Vector{Int},Vector{Vector{Int}}}()

    @inline function get_block_lengths(vec::Vector{Int})::Vector{Int}
        if isempty(vec)
            return Int[]
        end
        blocks = [1]
        current = vec[1]
        for x in vec[2:end]
            if x == current
                blocks[end] += 1
            else
                push!(blocks, 1)
                current = x
            end
        end
        return blocks
    end

    for vec in all_vectors
        block_lengths = get_block_lengths(vec)
        key = sort(block_lengths)
        if haskey(groups, key)
            push!(groups[key], vec)
        else
            groups[key] = [vec]
        end
    end

    @inline function count_permutations(block_lengths::Vector{Int})::Int
        counts = Dict{Int,Int}()
        for len in block_lengths
            counts[len] = get(counts, len, 0) + 1
        end
        k = length(block_lengths)
        permutations = factorial(k)
        for cnt in values(counts)
            permutations ÷= factorial(cnt)
        end
        return permutations
    end

    result = Vector{Tuple{Vector{Int},Int}}()
    for (block_lengths, vecs) in groups
        # For the strong coupling expansion based on the Hubbard atom, local-site hoppings are not allowed
        if 2 * maximum(block_lengths) > sum(block_lengths)
            continue
        end

        representative = vecs[1]
        symmetry_factor = count_permutations(block_lengths)
        push!(result, (representative, symmetry_factor))
    end

    return result
end

function assign_orbitals(num_sites::Int, orbital_options=[[1, 1], [2, 2]])
    # Each site can have either [1,1] or [2,2]
    all_combinations = Iterators.product(ntuple(_ -> orbital_options, num_sites)...)
    return [collect(comb) for comb in all_combinations]
end

function partition(order::Int)
    par = [
        # order 1
        (1, 0),
        # order 2
        (2, 0), (1, 1),
        # order 3
        (3, 0), (2, 1), (1, 2),
        # order 4
        (4, 0), (3, 1), (2, 2), (1, 3),
        #order 5
        (5, 0), (4, 1), (3, 2), (2, 3), (1, 4),
        #order 6
        (6, 0), (5, 1), (4, 2), (3, 3), (2, 4), (1, 5),
    ]
    return sort([p for p in par if p[1] + p[2] <= order])
end

function free_energy(_partition::Vector{T}; filter=[], leaf_dep_funcs::Vector{Function}=Function[pr->pr isa BareHoppingId]) where {T}

    diagpara = []
    inter = [Interaction(UpDown, [Dynamic])]

    orders = union([p[1] for p in _partition])
    max_order = maximum(orders)
    min_order = minimum(orders)
    max_totalorder = maximum([sum(p) for p in _partition])
    dict_graphs = Dict{NTuple{2,Int},Vector{Graph}}()

    for order in orders
        para = DiagPara(type=VacuumDiag, innerLoopNum=order, hasTau=true, interaction=inter, totalTauNum=2order, filter=filter)
        push!(diagpara, para)
        println("Order: ", order)

        topologies = generate_topologies(order)

        extT = [[2 * i - 1, 2 * i] for i in 1:order]
        creations = [[true, false] for _ in 1:order]

        graphs_fE = Graph[]
        sub_factors = []
        orbitals = assign_orbitals(order)

        tops = []
        for (sites, factor) in topologies
            println("sites: ", sites)
            # println("factor: ", factor)
            for orbital in orbitals
                println(orbital, "t: ", extT)
                push!(graphs_fE, SCE.connectedGreen(para, sites, orbital, extT, creations))
                push!(sub_factors, factor)

                push!(tops, [sites, orbital])
            end
        end

        println("len of graphs: ", length(graphs_fE))

        property = VacuumId(para)
        graph_order = [Graph(graphs_fE, operator=Sum(), properties=property, name=Symbol("F_$order"))]

        optimize!(graph_order)
        optimize!(graph_order)

        renormalization_orders = [max_totalorder - order]

        dict_graph_order = taylorAD(graph_order, renormalization_orders, leaf_dep_funcs)
        for key in keys(dict_graph_order)
            p = (order, key...)
            if p in _partition
                dict_graphs[p] = dict_graph_order[key]
            end
        end
    end

    diagpara = Vector{DiagPara}()
    partitions = sort(collect(keys(dict_graphs)))
    for p in partitions
        push!(diagpara, DiagPara(type=VacuumDiag, innerLoopNum=p[1], hasTau=true, interaction=inter, totalTauNum=2 * p[1], filter=filter))
    end

    return (partitions, diagpara, dict_graphs)
end
