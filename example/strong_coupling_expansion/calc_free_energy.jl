push!(LOAD_PATH, pwd())
using Atom
using Lehmann
using MCIntegration
using Printf
using Measurements
using JLD2
using DataStructures

include("free_energy.jl")

struct ParaMC
	μ::Float64
	U::Float64
	β::Float64
	n::Int
	Lx::Int
	Ly::Int
	lambda::Float64
	order::Int
end

paraid(p::ParaMC) = Dict(
	"order" => p.order,
	"beta" => p.β,
	"lambda" => p.lambda,
	"mu" => p.μ,
	"U" => p.U,
	"Lx" => p.Lx,
	"Ly" => p.Ly,
)
short(p::ParaMC) = join(["$(k)_$(v)" for (k, v) in sort!(OrderedDict(paraid(p)))], "_")

# t, U, μ, β, n = 1.0, 3.0, 1.0, 1.0, 0
# t, U, μ, β, n = 1.0, 8.0, 1.0, 1.0, 0
t, U, μ, β, n = 1.0, 4.0, 0.5, 2.0, 0
Lx, Ly = 2, 1
# lam = 0.2
# lam = 0.1
lam = 0.01
order = 2

para = ParaMC(μ, U, β, n, Lx, Ly, lam, order)
m = Hubbard.hubbardAtom(:fermi, U, μ, β)


function disperion_PBC(Lx, Ly, t)
	No = 2 # spin up/down
	ϵk = zeros(Float64, (No, No, Lx, Ly)) # julia column major, the first index is the major index

	for xi in 1:Lx
		for yi in 1:Ly
			kx, ky = 2π * (xi - 1) / Lx, 2π * (yi - 1) / Ly
			ϵk[1, 1, xi, yi] = -2 * t * (cos(kx) + cos(ky))
			ϵk[2, 2, xi, yi] = -2 * t * (cos(kx) + cos(ky))
		end
	end
	return ϵk
end

function disperion_FBC(Lx, Ly, t)
	No = 2 # spin up/down
	ϵk = zeros(Float64, (No, No, Lx, Ly)) # julia column major, the first index is the major index

	for xi in 1:Lx
		for yi in 1:Ly
			k = [π * xi / (Lx + 1), π * yi / (Ly + 1)]
			ϵk[1, 1, xi, yi] = -2t * sum(cos.(k))
			ϵk[2, 2, xi, yi] = -2t * sum(cos.(k))
		end
	end
	return ϵk
end

const ϵk = disperion_FBC(para.Lx, para.Ly, t)

function propagator(τ::T, ω::T, β::T) where {T}
	if τ ≈ T(0.0)
		τ = -1e-10
	end
	if τ > T(0.0)
		return ω > T(0.0) ?
			   exp(-ω * τ) / (1 + exp(-ω * β)) :
			   exp(ω * (β - τ)) / (1 + exp(ω * β))
	else
		return ω > T(0.0) ?
			   -exp(-ω * (τ + β)) / (1 + exp(-ω * β)) :
			   -exp(-ω * τ) / (1 + exp(ω * β))
	end
end

function propagator_derivative(τ, ϵ, β, order)
	if order == 0
		result = propagator(τ, ϵ, β)
	elseif order == 1
		result = -Spectral.kernelFermiT_dω(τ, ϵ, β)
	elseif order == 2
		result = Spectral.kernelFermiT_dω2(τ, ϵ, β) / 2.0
	elseif order == 3
		result = -Spectral.kernelFermiT_dω3(τ, ϵ, β) / 6.0
	elseif order == 4
		result = Spectral.kernelFermiT_dω4(τ, ϵ, β) / 24.0
	elseif order == 5
		result = -Spectral.kernelFermiT_dω5(τ, ϵ, β) / 120.0
	else
		error("not implemented!")
	end
	return result
end

function green_counterterm_PBC(para::ParaMC, τ::T, rx::Int, ry::Int, orbital::Int, order::Int) where {T}
	β, Lx, Ly = para.β, para.Lx, para.Ly
	g2c = 0.0
	N = Lx * Ly
	for xi in 1:Lx
		for yi in 1:Ly
			kx, ky = 2π * (xi - 1) / Lx, 2π * (yi - 1) / Ly
			ω = -1.0 / ϵk[orbital, orbital, xi, yi] / para.lambda
			g2c_τ = 0.0
			for o in 0:order
				g2c_τ += propagator_derivative(τ, ω, β, o) * ω^o * binomial(order, o) * (-1)^o
			end
			# g2c += exp(im * (kx * rx + ky * ry)) * g2c_τ / para.lambda / N
			g2c += cos(kx * rx + ky * ry) * g2c_τ / para.lambda / N
		end
	end
	return g2c
end

function green_counterterm_FBC(para::ParaMC, τ::T, r1::Vector{Int}, r2::Vector{Int}, orbital::Int, order::Int) where {T}
	β, Lx, Ly = para.β, para.Lx, para.Ly
	g2c = 0.0
	prefactor = 4 / (Lx + 1) / (Ly + 1)
	for xi in 1:Lx
		for yi in 1:Ly
			k = [π * xi / (Lx + 1), π * yi / (Ly + 1)]
			ω = -1.0 / ϵk[orbital, orbital, xi, yi] / para.lambda
			g2c_τ = 0.0
			for o in 0:order
				g2c_τ += propagator_derivative(τ, ω, β, o) * ω^o * binomial(order, o) * (-1)^o
				# g2c_τ += propagator_derivative(τ, ω, β, o) * ω^o * binomial(order, o)
			end

			phi_r1 = prod(sin.(k .* r1))
			phi_r2 = prod(sin.(k .* r2))
			g2c += phi_r1 * phi_r2 * g2c_τ / para.lambda
		end
	end
	return g2c * prefactor
end


function integrand(idx, vars, config)
	para, root, graphfuncs! = config.userdata[1:3]
	leafval, leafType, leafOrders, leafSites, leafτ_i, leafτ_o, leaforbitals_i, leaforbitals_o = config.userdata[4]
	model = config.userdata[5]
	varT, varR = vars

	for (i, lftype) in enumerate(leafType[idx])
		if lftype == 0
			continue
		elseif lftype == 3  # BareGreenNId
			# τ = vcat(varT[leafτ_i[idx][i]], varT[leafτ_o[idx][i]])
			τ = vcat(varT[leafτ_o[idx][i]], varT[leafτ_i[idx][i]])
			orbitals = vcat(leaforbitals_o[idx][i], leaforbitals_i[idx][i])
			_gn = Green.GreenN(model, τ, orbitals)
			# leafval[idx][i] = Green.Gnc(model, _gn)
			leafval[idx][i] = Green.Gn(model, _gn)
		elseif lftype == 4  # BareHoppingId
			# if leaforbitals_i[idx][i][1] != leaforbitals_o[idx][i][1]
			# 	leafval[idx][i] = 0.0
			# 	continue
			# end
			τ = varT[leafτ_o[idx][i][1]] - varT[leafτ_i[idx][i][1]]
			r1x, r2x = varR[leafSites[idx][i][1]], varR[leafSites[idx][i][2]]
			r1y, r2y = 1, 1
			# @assert leaforbitals[idx][i][1] == leaforbitals[idx][i][2]

			order = leafOrders[idx][i][1]
			orbital = leaforbitals_i[idx][i][1]
			leafval[idx][i] = green_counterterm_FBC(para, τ, [r1x, r1y], [r2x, r2y], orbital, order)
		else
			error("this leaftype $lftype not implemented!")
		end
	end

	graphfuncs![idx](root, leafval[idx])

	return root[1]
end

function freeE(model, para::ParaMC, diagram; neval = 1e6, print = 0, kwargs...)
	partition, diagpara, FeynGraphs = diagram

	funcGraphs! = Dict{Int, Function}()
	leaf_maps = Vector{Dict{Int, Graph}}()
	for (i, key) in enumerate(partition)
		funcGraphs![i], leafmap = Compilers.compile(FeynGraphs[key])
		push!(leaf_maps, leafmap)
		# println(funcGraphs![i])
	end

	leafStat = FeynmanDiagram.leafstates(leaf_maps)

	root = zeros(Float64, 1)
	# T = Continuous(0.0, para.β; offset = 1, adapt = true)
	# T.data[1] = 0.0
	T = Continuous(0.0, para.β; adapt = true)
	R = Discrete(1, 2)

	dof = [[p.totalTauNum - 1, p.innerLoopNum] for p in diagpara]
	obs = zeros(Float64, length(diagpara))

	config = Configuration(; var = (T, R), dof = dof, obs = obs, type = Float64,
		userdata = (para, root, funcGraphs!, leafStat, model))
	result = integrate(integrand; config = config, neval = neval, print = print, solver = :mcmc, kwargs...)

	if isnothing(result) == false
		# if print >= 0
		# 	report(result.config)
		# 	println(report(result, pick = o -> first(o)))
		# 	println(result)
		# end
		if print >= -2
			println(result)
		end

		datadict = Dict{eltype(partition), Any}()
		for (o, key) in enumerate(partition)
			avg, std = result.mean[o], result.stdev[o]
			# r = measurement.(real(avg), real(std))
			# i = measurement.(imag(avg), imag(std))
			# data = Complex.(r, i)
			datadict[key] = measurement.(avg, std) / (-para.β)
		end
		return datadict, result
	else
		return nothing, nothing
	end
end

function freeE_MC(model, para::ParaMC; neval = 1e6, partition = partition(para.order), reweight_goal = nothing,
	print = 0, filename::Union{String, Nothing} = nothing)
	# partition, diagpara, FeynGraphs = free_energy(partition)
	diagram = free_energy(partition)

	partition = diagram[1]
	println("partition: ", partition)
	if isnothing(reweight_goal)
		reweight_goal = Float64[]
		for (order, sOrder) in partition
			if sOrder == 0
				push!(reweight_goal, 2.0)
			else
				push!(reweight_goal, 1.0)
			end
		end
		push!(reweight_goal, 2.0)
	end

	_neighbor = neighbor(partition)

	freeEnergy, result = freeE(model, para, diagram; neval = neval, neighbor = _neighbor, reweight_goal = reweight_goal, print = print)

	if isnothing(freeEnergy) == false
		if isnothing(filename) == false
			jldopen(filename, "a+") do f
				key = "$(short(para))"
				if haskey(f, key)
					@warn("replacing existing data for $key")
					delete!(f, key)
				end
				f[key] = (freeEnergy,)
			end
		end
		for (ip, key) in enumerate(partition)
			println("Group ", key)
			@printf("%10s   %10s \n", "avg", "err")
			@printf("%10.6f ± %10.6f\n", freeEnergy[key].val, freeEnergy[key].err)
		end
	end
	return freeEnergy, result
end
# _partition = [(1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (2, 0), (2, 1), (2, 2), (2, 3), (2, 4), (3, 0), (3, 1), (3, 2)]
# _partition = [(1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (2, 0), (2, 1), (2, 2), (2, 3), (2, 4)]
_partition = [(1, 0), (1, 1), (1, 2), (1, 3), (2, 0), (2, 1), (2, 2)]
# _partition = partition(3)
# freeE_MC(m, para, partition = _partition, neval = 1e8, filename = "data_freeE.jld2")
# freeE_MC(m, para, partition = _partition, neval = 4e6, filename = "data_freeE.jld2")
freeE_MC(m, para, partition = _partition, neval = 2e6)

# println(res)
