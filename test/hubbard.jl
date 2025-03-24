
using Lehmann
using MCIntegration

@testset "Hubbard Atom" begin
	import FeynmanDiagram.FrontEnds: ChargeCharge, SpinSpin, UpUp, UpDown, Instant
	import FeynmanDiagram.Parquet: DiagPara, Interaction, SigmaDiag
	struct ParaMC
		μ::Float64
		U::Float64
		β::Float64
		n::Int
	end

	function compare(a, b, err, ratio = 3.0)
		@test abs(real(a) - real(b)) < ratio * real(err)
		@test abs(imag(a) - imag(b)) < ratio * imag(err)
	end

	function sigma(orders::AbstractVector, filter = [])
		inter = [Interaction(UpDown, [Instant])]
		sigma = Dict{Int, Vector{Graph}}()
		diagpara = []

		extT_labels = Vector{Vector{Int}}[]
		for o in orders
			p = DiagPara(type = SigmaDiag, innerLoopNum = o, hasTau = true, interaction = inter, filter = filter)
			push!(diagpara, p)
			graph_df = Parquet.build(p)
			optimize!(graph_df.diagram, level = 1)
			optimize!(graph_df.diagram, level = 1)
			sigma[o] = graph_df.diagram
			# plot_tree(_s, maxdepth=15)
			push!(extT_labels, [collect(g.properties.extT) for g in graph_df.diagram])
		end

		return (diagpara, sigma, extT_labels)
	end

	@inline function phase(varT, extT, l, β)
		tin, tout = varT[extT[1]], varT[extT[2]]
		return exp(1im * π * (2l + 1) / β * (tout - tin))
		# return cos(π * (2l + 1) / β * (tout - tin))
	end

	function propagator(τ, para)
		if τ ≈ 0.0
			τ = -1e-8
		end
		ϵ, β = -para.μ, para.β
		return Spectral.kernelFermiT(τ, ϵ, β)
	end


	function integrand(idx, vars, config)
		para, extT_labels, root, graphfuncs! = config.userdata[1:4]
		leafval, leafType, leafOrders, leafτ_i, leafτ_o, leafMomIdx = config.userdata[5]
		varT = vars

		for (i, lftype) in enumerate(leafType[idx])
			if lftype == 0
				continue
			elseif lftype == 1
				τ = varT[leafτ_o[idx][i]] - varT[leafτ_i[idx][i]]
				leafval[idx][i] = propagator(τ, para)
			elseif lftype == 2
				leafval[idx][i] = para.U
			else
				error("this leaftype $lftype not implemented!")
			end
		end

		graphfuncs![idx](root, leafval[idx])

		w = sum(root[i] * phase(varT, extT, para.n, para.β) for (i, extT) in enumerate(extT_labels[idx]))
		return w #the current implementation of sigma has an additional minus sign compared to the standard defintion
	end

	function sigmaMC(para::ParaMC, orders::Vector{Int}, neval; kwargs...)
		μ, U, β = para.μ, para.U, para.β
		T = Continuous(0.0, β; offset = 1)
		T.data[1] = 0.0

		diagpara, FeynGraphs, extT = sigma(orders)

		funcGraphs! = Dict{Int, Function}()
		leaf_maps = Vector{Dict{Int, Graph}}()
		for (i, key) in enumerate(orders)
			funcGraphs![i], leafmap = Compilers.compile(FeynGraphs[key])
			push!(leaf_maps, leafmap)
		end
		leafStat, loopbasis = FeynmanDiagram.leafstates(leaf_maps, maximum(orders) + 1)

		root = zeros(ComplexF64, maximum(length.(extT)))
		dof = [[diagpara[o].totalTauNum - 1] for o in 1:length(orders)]
		obs = zeros(ComplexF64, length(orders))

		config = Configuration(; var = (T,), dof = dof, obs = obs,
			userdata = (para, extT, root, funcGraphs!, leafStat), type = ComplexF64)
		result = integrate(integrand; config = config, neval = neval, solver = :mcmc, kwargs...)
		# ExprTree.showTree(diag[1], 3)
		return result
	end

	U = 1.0
	β = 2.3
	μ = 0.0
	n = 0
	para = ParaMC(μ, U, β, n)
	neval = 1e6
	print = 0
	orders = [1, 2, 3, 4]
	# orders = [2]
	result = sigmaMC(para, orders, neval; print = print)
	avg, std = result.mean, result.stdev
	expect = [-U / 2,
		(2im + π) * β * U^2 / 8 / π,
		(4 - π^2) * β^2 * U^3 / 32 / π^2,
		-(24im - 12π + 6im * π^2 + π^3) * β^3 * U^4 / 384 / π^3,
	]

	for (oi, o) in enumerate(orders)
		println("order $o :  ", avg[oi], " +- ", std[oi], "  ~  ", expect[o])
		compare(avg[oi], expect[o], std[oi])
	end
end
