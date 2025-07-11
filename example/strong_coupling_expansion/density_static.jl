using FeynmanDiagram
import FeynmanDiagram.Parquet: DiagPara, Interaction, VacuumDiag
import FeynmanDiagram.ComputationalGraphs: Sum
import FeynmanDiagram.FrontEnds: ConnectedGreenNId, BareHoppingId, VacuumId, GenericId, UpUp, UpDown, Dynamic
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

function density_info(_partition::Vector{T}; filter=[], leaf_dep_funcs::Vector{Function}=Function[pr->pr isa BareHoppingId], num_orbitals::Int=2) where {T}

    diagpara = []
    inter = [Interaction(UpDown, [Dynamic])]

    orders = union([p[1] for p in _partition])
    max_totalorder = maximum([sum(p) for p in _partition])
    dict_graphs = Dict{NTuple{2,Int},Vector{Graph}}()

    hop_orbitals = [[i, i] for i in 1:num_orbitals]

    external_T = [1, 2]
    external_sites = [1, 1]
    external_creation = [true, false]
    external_orbitals = [1, 1]
    for order in orders
        para = DiagPara(type=VacuumDiag, innerLoopNum=order, hasTau=true, interaction=inter, totalTauNum=order, filter=filter)
        push!(diagpara, para)
        println("Order: ", order)

        graphs_fE = Graph[]
        orbitals_all = assign_orbitals(order, hop_orbitals)

        for orbital in orbitals_all
            hoppings = BareHoppingId[]
            for hop_idx in 1:order
                push!(hoppings, BareHoppingId(para, (2hop_idx, 2hop_idx + 1), Tuple(orbital[hop_idx]), (hop_idx + 2, hop_idx + 2)))
            end
            push!(graphs_fE, SCE.fullGreen_v1(para, hoppings; external_T=external_T, external_sites=external_sites,
                external_creation=external_creation, external_orbitals=external_orbitals,
                prefactor=(-1)^order / factorial(order)))
        end

        property = GenericId(para)

        graph_order = [Graph(graphs_fE, operator=Sum(), properties=property, name=Symbol("n_$order"))]
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
        push!(diagpara, DiagPara(type=VacuumDiag, innerLoopNum=p[1], hasTau=true, interaction=inter, totalTauNum=p[1], filter=filter))
    end

    return (partitions, diagpara, dict_graphs)
end