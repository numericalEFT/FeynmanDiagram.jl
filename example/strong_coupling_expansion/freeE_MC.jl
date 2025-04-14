include("./input.jl")
include("./calc_free_energy.jl")

for (_μ, _U, _β, lam, order) in Iterators.product(μ, U, β, lambdas, orders)
	ϵk = disperion_FBC(Lx, Ly, t)
	para = ParaMC(_μ, _U, _β, 0, Lx, Ly, lam, order, ϵk)
	println(short(para))

	model = Hubbard.hubbardAtom(:fermi, _U, _μ, _β)

	_partition = partition(order)
	# reweight_goal = Float64[]
	# for (order, sOrder) in partition
	# 	reweight_factor = 2.0^(2order + 2sOrder - 2)
	# 	if (order, sOrder) == (1, 0)
	# 		reweight_factor = 4.0
	# 	end
	# 	push!(reweight_goal, reweight_factor)
	# end
	# push!(reweight_goal, 4.0)

	freeE_MC(model, para, partition = _partition, neval = neval, filename = freeE_filename)

	# println(res)
	println("F0 = ", F0(para), "\n")
end
