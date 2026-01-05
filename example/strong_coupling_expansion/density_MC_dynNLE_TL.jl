include("./input.jl")
# include("./calc_N_dynNLE_2DTL.jl")
include("./calc_N_dynNLErec_2D.jl")


for (_μ, _U, _β, lam, _dμ, order) in Iterators.product(μ, U, β, lambdas, dμ, orders)
    model = Hubbard.hubbardAtom(:fermi, _U, _μ + _dμ, _β)

    # _partition = [(2, 0), (2, 1), (2, 2), (3, 0), (3, 1), (3, 2)]
    # _partition = [(2, 0), (2, 1), (2, 2), (2, 3), (3, 0), (3, 1), (3, 2), (3, 3)]
    # _partition = [(1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2), (3, 0), (3, 1), (3, 2)]
    # deriv_order = maximum(p[2] for p in _partition)

    _partition = partition_dynmu(order)

    para = ParaMC(
        μ=_μ,
        U=_U,
        t=1.0,
        β=_β,
        lambda=lam,
        dμ=_dμ,
        order=order,
        Lkx=Lkx,
        Lky=Lky,
        Lx=Lx,
        Ly=Ly,
        Rmax=Rmax
    )

    # _partition = [(2, 0), (2, 1), (2, 2), (2, 3), (2, 4), (3, 0), (3, 1)]
    # _partition = [(2, 0), (2, 1), (2, 2), (2, 3), (2, 4), (4, 0), (4, 1), (4, 2), (4, 3), (4, 4)]

    # reweight_goal = Float64[]
    # for (order, sOrder) in partition
    # 	reweight_factor = 2.0^(2order + 2sOrder - 2)
    # 	if (order, sOrder) == (1, 0)
    # 		reweight_factor = 4.0
    # 	end
    # 	push!(reweight_goal, reweight_factor)
    # end
    # push!(reweight_goal, 4.0)

    density_MC(model, para, partition=_partition, neval=neval, filename=N_filename,
        dtype=Float64, Ntau=20000, file_pretab=pretab_filename)#, print=1)
end
