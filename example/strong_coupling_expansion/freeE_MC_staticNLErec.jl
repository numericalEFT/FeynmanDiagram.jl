include("./input.jl")
include("./common.jl")
# include("./calc_free_energy_staticNLErec.jl")
include("./calc_free_energy_staticNLErec_2D.jl")

function partition(order::Int, hasodd::Bool=false)
    par = []
    dord = hasodd ? 1 : 2
    for i in 2:dord:order
        push!(par, (i, 0))
    end
    return par
end

for (_μ, _U, _β, _dμ, order) in Iterators.product(μ, U, β, dμ, orders)
    para = ParaMC(_μ, _U, t, _β, 0, Lx, Ly, _dμ, order)
    println(short(para))

    model = Hubbard.hubbardAtom(:fermi, _U, _μ + _dμ, _β)

    _partition = partition_static(order)
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
        dtype=Float64)#, _neighbor=neighbor(_partition)) print=1)
end
