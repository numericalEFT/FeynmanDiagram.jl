include("./input.jl")
include("./calc_free_energy_staticNLE.jl")

function neighbor(partitions)
    n = Vector{Tuple{Int,Int}}()
    Nnorm = length(partitions) + 1 # the index of the normalization diagram is the N+1
    for (ip, p) in enumerate(partitions)
        if p[1] in [0, 1] # if there is only one loop, then the diagram can be connected to the normalization diagram
            push!(n, (ip, Nnorm))
        end
        for (idx, np) in enumerate(partitions)
            if idx >= ip
                continue
            end
            if np[1] == p[1] || np[1] == p[1] + 1 || np[1] == p[1] - 1 #the first index is the number of loops
                # if np[1] == p[1] || np[1] == p[1] + 2 || np[1] == p[1] - 2 #the first index is the number of loops
                push!(n, (ip, idx))
            end
        end
    end
    # println(n)
    return n
end

for (_μ, _U, _β, order) in Iterators.product(μ, U, β, orders)
    para = ParaMC(_μ, _U, t, _β, 0, Lx, Ly, order)
    println(short(para))

    model = Hubbard.hubbardAtom(:fermi, _U, _μ, _β)

    _partition = [(2, 0), (3, 0), (4, 0)]
    # _partition = [(2, 0), (3, 0),]
    # _partition = [(2, 0), (4, 0),]
    # _partition = [(2, 0), (4, 0), (6, 0)]
    # reweight_goal = Float64[]
    # for (order, sOrder) in partition
    # 	reweight_factor = 2.0^(2order + 2sOrder - 2)
    # 	if (order, sOrder) == (1, 0)
    # 		reweight_factor = 4.0
    # 	end
    # 	push!(reweight_goal, reweight_factor)
    # end
    # push!(reweight_goal, 4.0)

    freeE_MC(model, para, partition=_partition, neval=neval, filename=freeE_filename,
        dtype=Float64, _neighbor=neighbor(_partition))#, print=1)
end
