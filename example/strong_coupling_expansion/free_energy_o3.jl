using FeynmanDiagram
import FeynmanDiagram.Parquet: DiagPara, Interaction, VacuumDiag
import FeynmanDiagram.ComputationalGraphs: Sum
import FeynmanDiagram.FrontEnds: ConnectedGreenNId, BareHoppingId, VacuumId, UpUp, UpDown, Dynamic
using Parameters

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

function neighbor(partitions)
    n = Vector{Tuple{Int,Int}}()
    Nnorm = length(partitions) + 1 # the index of the normalization diagram is the N+1
    for (ip, p) in enumerate(partitions)
        # if p[1] == 1 # if there is only one loop, then the diagram can be connected to the normalization diagram
        if p[1] in [0, 1, 2] # if there is only one loop, then the diagram can be connected to the normalization diagram
            push!(n, (ip, Nnorm))
        end
        for (idx, np) in enumerate(partitions)
            if idx >= ip
                continue
            end
            # if np[1] == p[1] || np[1] == p[1] + 2 || np[1] == p[1] - 2 #the first index is the number of loops
            if np[1] == p[1] || np[1] == p[1] + 1 || np[1] == p[1] - 1 #the first index is the number of loops
                push!(n, (ip, idx))
            end
        end
    end
    # println(n)
    return n
end

function free_energy(_partition::Vector{T}; filter=[], leaf_dep_funcs::Vector{Function}=Function[pr->pr isa BareHoppingId], num_orbitals::Int=2) where {T}

    diagpara = []
    inter = [Interaction(UpDown, [Dynamic])]

    orders = union([p[1] for p in _partition])
    max_totalorder = maximum([sum(p) for p in _partition])
    dict_graphs = Dict{NTuple{2,Int},Vector{Graph}}()

    hop_orbitals = [[i, i] for i in 1:num_orbitals]

    for order in orders
        para = DiagPara(type=VacuumDiag, innerLoopNum=order, hasTau=true, interaction=inter, totalTauNum=2order, filter=filter)
        push!(diagpara, para)
        println("Order: ", order)

        graphs_fE = Graph[]
        orbitals_all = assign_orbitals(order, hop_orbitals)
        if order == 3
            for orbital in [(1, 1), (2, 2)]
                hoppings = BareHoppingId[]
                for hop_idx in 1:order
                    push!(hoppings, BareHoppingId(para, (2hop_idx - 1, 2hop_idx), orbital, (2hop_idx - 1, 2hop_idx)))
                end
                push!(graphs_fE, SCE.connectedGreen_o3(para, hoppings, prefactor=1.0 / factorial(order)))
            end
        elseif order == 2
            for orbital in [(1, 1), (2, 2)]
                hoppings = BareHoppingId[]
                for hop_idx in 1:order
                    push!(hoppings, BareHoppingId(para, (2hop_idx - 1, 2hop_idx), orbital, (2hop_idx - 1, 2hop_idx)))
                end
                push!(graphs_fE, SCE.connectedGreen_o2(para, hoppings, prefactor=1.0 / factorial(order)))
            end
        else
            for orbital in orbitals_all
                hoppings = BareHoppingId[]
                for hop_idx in 1:order
                    push!(hoppings, BareHoppingId(para, (2hop_idx - 1, 2hop_idx), Tuple(orbital[hop_idx]), (2hop_idx - 1, 2hop_idx)))
                end
                push!(graphs_fE, SCE.connectedGreen(para, hoppings, prefactor=1.0 / factorial(order)))
            end
        end

        println("len of graphs: ", length(graphs_fE))

        # order == 3 && exit()

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
