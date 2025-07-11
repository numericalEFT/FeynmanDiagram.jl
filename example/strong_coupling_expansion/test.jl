
β, Lx, Ly = 0.5, 2, 1
lambda = 1.0

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

const ϵk = disperion_FBC(Lx, Ly, 1.0)

function green_counterterm_FBC(τ::T, r1::Vector{Int}, r2::Vector{Int}, orbital::Int; order::Int = 0) where {T}

	g2c = 0.0
	prefactor = 4 / (Lx + 1) / (Ly + 1)
	for xi in 1:Lx
		for yi in 1:Ly
			k = [π * xi / (Lx + 1), π * yi / (Ly + 1)]
			println(k ./ π)
			ω = -1.0 / ϵk[orbital, orbital, xi, yi] / lambda
			println("ω:", ω)
			g2c_τ = propagator(τ, ω, β)
			# for o in 0:order
			# 	g2c_τ += propagator(τ, ω, β) * ω^o * binomial(order, o) * (-1)^o
			# 	# g2c_τ += propagator_derivative(τ, ω, β, o) * ω^o * binomial(order, o)
			# end

			println("g2:", g2c_τ)
			phi_r1 = prod(sin.(k .* r1))
			phi_r2 = prod(sin.(k .* r2))
			println("phi_r1:", phi_r1, "phi_r2:", phi_r2)
			g2c += phi_r1 * phi_r2 * g2c_τ / lambda
		end
	end
	println("prefactor: ", prefactor)
	return g2c * prefactor
end

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

green_counterterm_FBC(0.0, [2, 1], [2, 1], 1)
