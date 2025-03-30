using FeynmanDiagram
import FeynmanDiagram.Parquet: DiagPara, VacuumDiag
import FeynmanDiagram.ComputationalGraphs: Sum
import FeynmanDiagram.FrontEnds: ConnectedGreenNId, VacuumId

function generate_vectors(order::Int)
	function helper(current_vector::Vector{Int}, remaining_length::Int)
		if remaining_length == 0
			return [copy(current_vector)]
		end

		results = Vector{Vector{Int}}()
		last_value = current_vector[end]

		# Option 1: Repeat the last value
		push!(results, helper(vcat(current_vector, [last_value]), remaining_length - 1)...)

		# Option 2: Increment the last value
		push!(results, helper(vcat(current_vector, [last_value + 1]), remaining_length - 1)...)

		return results
	end

	return helper([1], order - 1)
end

function generate_topologies(order::Int)
	all_vectors = generate_vectors(order)
	groups = Dict{Vector{Int}, Vector{Vector{Int}}}()

	@inline function get_block_lengths(vec::Vector{Int})::Vector{Int}
		if isempty(vec)
			return Int[]
		end
		blocks = [1]
		current = vec[1]
		for x in vec[2:end]
			if x == current
				blocks[end] += 1
			else
				push!(blocks, 1)
				current = x
			end
		end
		return blocks
	end

	for vec in all_vectors
		block_lengths = get_block_lengths(vec)
		key = sort(block_lengths)
		if haskey(groups, key)
			push!(groups[key], vec)
		else
			groups[key] = [vec]
		end
	end

	@inline function count_permutations(block_lengths::Vector{Int})::Int
		counts = Dict{Int, Int}()
		for len in block_lengths
			counts[len] = get(counts, len, 0) + 1
		end
		k = length(block_lengths)
		permutations = factorial(k)
		for cnt in values(counts)
			permutations ÷= factorial(cnt)
		end
		return permutations
	end

	result = Vector{Tuple{Vector{Int}, Int}}()
	for (block_lengths, vecs) in groups
		representative = vecs[1]
		symmetry_factor = count_permutations(block_lengths)
		push!(result, (representative, symmetry_factor))
	end

	return result
end

function free_energy(max_order::Int)

	fE = Dict{Int, Graph}()
	for order in 1:max_order
		para = DiagPara(type = VacuumDiag, innerLoopNum = order, hasTau = true)
		println("Order: ", order)

		topologies = generate_topologies(order)

		extT = [[2 * i - 1, 2 * i] for i in 1:order]
		creations = [(true, false) for _ in 1:order]

		graphs_fE = []
		sub_factors = []
		for (sites, factor) in topologies
			push!(graphs_fE, SCE.connectedGreen(para, sites, extT, extT, creations))
			push!(sub_factors, factor)
		end

		property = VacuumId(para)
		fE[order] = Graph(graphs_fE, subgraph_factors = sub_factors, operator = Sum(), properties = property, name = Symbol("F_$order"))
	end

	return fE
end


# fE = free_energy(3)
