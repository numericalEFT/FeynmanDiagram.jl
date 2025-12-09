using FeynmanDiagram
import FeynmanDiagram.Parquet: DiagPara, Interaction, VacuumDiag, GreenDiag
import FeynmanDiagram.ComputationalGraphs: Sum
import FeynmanDiagram.FrontEnds: ConnectedGreenNId, BareHoppingId, BareGreenNId, GreenNId, VacuumId, UpUp, UpDown, Dynamic
using Parameters
using JLD2

current_dir = @__DIR__
file_path = joinpath(current_dir, "graphs_SCEo8.jld2")

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

function free_energy_recursion(_partition::Vector{T}; filter=[], leaf_dep_funcs::Vector{Function}=Function[pr->pr isa BareHoppingId],
    num_orbitals::Int=2, dynamic_hop=true) where {T}

    topologies = jldopen(file_path, "r")["free_energy"]
    order_numvar = collect(keys(topologies))

    diagpara = []
    inter = [Interaction(UpDown, [Dynamic])]

    orders = union([p[1] for p in _partition])
    max_totalorder = maximum([sum(p) for p in _partition])
    dict_graphs = Dict{NTuple{3,Int},Vector{Graph}}()
    gc_pool, gn_pool = Dict{Vector{BareHoppingId},Graph}(), Dict{Vector{BareHoppingId},Graph}()

    hop_orbitals = [[i, i] for i in 1:num_orbitals]

    for order in orders
        para = DiagPara(type=VacuumDiag, innerLoopNum=order, hasTau=true, interaction=inter, totalTauNum=order, filter=filter)
        para_hop = DiagPara(type=GreenDiag, innerLoopNum=0, hasTau=true, interaction=inter)
        push!(diagpara, para)

        orbitals_all = assign_orbitals(order, hop_orbitals)
        keys_topo = [k for k in order_numvar if k[1] == order]

        for key_topo in keys_topo  # given perturbation order and independent number of site variables
            sym_factors = Int[]
            vec_hop_set = Vector{Vector{BareHoppingId}}()
            for (hop_inds, sym_factor) in topologies[key_topo]
                for orbital in orbitals_all
                    hoppings = BareHoppingId[]
                    for hop_idx in 1:order
                        if dynamic_hop
                            push!(hoppings, BareHoppingId(para, hop_inds[hop_idx], Tuple(orbital[hop_idx]), (2hop_idx - 1, 2hop_idx)))
                        else
                            push!(hoppings, BareHoppingId(para, hop_inds[hop_idx], Tuple(orbital[hop_idx]), (hop_idx, hop_idx)))
                        end
                    end
                    push!(vec_hop_set, hoppings)
                    push!(sym_factors, sym_factor)
                end
            end
            # println("len of graphs $key_topo: ", length(graphs_fE))
            graphs_fE, _, _ = SCE.connectedGreen!(para, vec_hop_set, gc_pool, gn_pool; name=Symbol("F_$order"), even=true,
                prefactors=1.0 ./ sym_factors, is_local_Gn=true)
            # prefactors=(-1)^order ./ sym_factors, is_local_Gn=true)

            graph_order = [graphs_fE,]
            optimize!(graph_order)
            optimize!(graph_order)

            renormalization_orders = [max_totalorder - order]

            dict_graph_order = taylorAD(graph_order, renormalization_orders, leaf_dep_funcs)
            for key in keys(dict_graph_order)
                p = (key_topo..., key...)
                if (p[1], p[3]) in _partition
                    dict_graphs[p] = dict_graph_order[key]
                end
            end
        end
    end

    diagpara = Vector{DiagPara}()
    partitions = sort(collect(keys(dict_graphs)))
    for p in partitions
        if dynamic_hop
            push!(diagpara, DiagPara(type=VacuumDiag, innerLoopNum=p[2], hasTau=true, interaction=inter, totalTauNum=p[1] * 2, filter=filter))
        else
            push!(diagpara, DiagPara(type=VacuumDiag, innerLoopNum=p[2], hasTau=true, interaction=inter, totalTauNum=p[1], filter=filter))
        end
    end

    return (partitions, diagpara, dict_graphs)
end

function generate_Gnderiv1(_partition::Vector{T}; filter=[],
    leaf_dep_funcs::Vector{Function}=[pr -> pr isa BareHoppingId, pr -> pr isa BareGreenNId],
    num_orbitals::Int=2, dynamic_hop=true) where {T}

    topologies = jldopen(file_path, "r")["free_energy"]
    order_numvar = collect(keys(topologies))

    diagpara = []
    inter = [Interaction(UpDown, [Dynamic])]

    orders = union([p[1] for p in _partition])
    max_totalorder = maximum([sum(p) for p in _partition])
    dict_graphs = Dict{NTuple{3,Int},Vector{Graph}}()
    gc_pool, gn_pool = Dict{Vector{BareHoppingId},Graph}(), Dict{Vector{BareHoppingId},Graph}()

    hop_orbitals = [[i, i] for i in 1:num_orbitals]

    for order in orders
        para = DiagPara(type=VacuumDiag, innerLoopNum=order, hasTau=true, interaction=inter, totalTauNum=order, filter=filter)
        para_hop = DiagPara(type=GreenDiag, innerLoopNum=0, hasTau=true, interaction=inter)
        push!(diagpara, para)

        orbitals_all = assign_orbitals(order, hop_orbitals)
        keys_topo = [k for k in order_numvar if k[1] == order]

        for key_topo in keys_topo  # given perturbation order and independent number of site variables
            sym_factors = Int[]
            vec_hop_set = Vector{Vector{BareHoppingId}}()
            for (hop_inds, sym_factor) in topologies[key_topo]
                for orbital in orbitals_all
                    hoppings = BareHoppingId[]
                    for hop_idx in 1:order
                        if dynamic_hop
                            push!(hoppings, BareHoppingId(para, hop_inds[hop_idx], Tuple(orbital[hop_idx]), (2hop_idx - 1, 2hop_idx)))
                        else
                            push!(hoppings, BareHoppingId(para, hop_inds[hop_idx], Tuple(orbital[hop_idx]), (hop_idx, hop_idx)))
                        end
                    end
                    push!(vec_hop_set, hoppings)
                    push!(sym_factors, sym_factor)
                end
            end
            # println("len of graphs $key_topo: ", length(graphs_fE))
            graphs_fE, _, _ = SCE.connectedGreen!(para, vec_hop_set, gc_pool, gn_pool; name=Symbol("F_$order"), even=true,
                prefactors=1.0 ./ sym_factors, is_local_Gn=true)

            graph_order = [graphs_fE,]
            optimize!(graph_order)
            optimize!(graph_order)

            renormalization_orders = [max_totalorder - order, 1]

            # println("renormalization_orders: ", renormalization_orders)

            dict_graph_order = taylorAD(graph_order, renormalization_orders, leaf_dep_funcs)
            for key in keys(dict_graph_order)
                p = (key_topo..., key[1])
                if (p[1], p[3]) in _partition && key[2] == 1
                    dict_graphs[p] = dict_graph_order[key]
                end
            end
        end
    end

    diagpara = Vector{DiagPara}()
    partitions = sort(collect(keys(dict_graphs)))
    for p in partitions
        if dynamic_hop
            push!(diagpara, DiagPara(type=VacuumDiag, innerLoopNum=p[2], hasTau=true, interaction=inter, totalTauNum=p[1] * 2, filter=filter))
        else
            push!(diagpara, DiagPara(type=VacuumDiag, innerLoopNum=p[2], hasTau=true, interaction=inter, totalTauNum=p[1], filter=filter))
        end
    end

    return (partitions, diagpara, dict_graphs)
end
