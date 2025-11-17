include("./input.jl")
include("./calc_D_dynNLErec_2D.jl")

for (_μ, _U, _β, lam, _dμ, order) in Iterators.product(μ, U, β, lambdas, dμ, orders)
    ϵk = disperion_PBC(Lx, Ly, t, _dμ)

    println(ϵk)

    para = ParaMC(_μ, _U, t, _β, 0, Lx, Ly, lam, _dμ, order, ϵk)
    println(short(para))

    model = Hubbard.hubbardAtom(:fermi, _U, _μ + _dμ, _β)

    # _partition = [(2, 0), (2, 1), (2, 2), (2, 3), (2, 4), (3, 0), (3, 1)]
    # _partition = [(2, 0), (2, 1), (2, 2), (2, 3), (3, 0), (3, 1), (3, 2), (3, 3), (4, 0), (4, 1), (4, 2), (4, 3)]
    _partition = [(2, 0), (2, 1), (2, 2), (2, 3), (3, 0), (3, 1), (3, 2), (3, 3)]
    # _partition = [(2, 0), (2, 1), (2, 2), (2, 3), (2, 4), (4, 0), (4, 1), (4, 2), (4, 3), (4, 4)]

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

    double_occupancy_MC(model, para, partition=_partition, neval=neval, filename=D_filename,
        dtype=Float64)#, _neighbor=neighbor(_partition)), print=1)
end
