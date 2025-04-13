
function parity(p)
	"""
	calculate the parity of a given permutation of the array [1, 2, 3, ...]
	"""
	n = length(p)
	not_seen = Set{Int}(1:n)
	seen = Set{Int}()
	cycles = Array{Int, 1}[]
	while !isempty(not_seen)
		cycle = Int[]
		x = pop!(not_seen)
		while !in(x, seen)
			push!(cycle, x)
			push!(seen, x)
			x = p[x]
			pop!(not_seen, x, 0)
		end
		push!(cycles, cycle)
	end
	cycle_lengths = map(length, cycles)
	even_cycles = filter(i -> i % 2 == 0, cycle_lengths)
	length(even_cycles) % 2 == 0 ? 1 : -1
end

struct Partition
	"""
	2-partition of 2N-legs of Green's functions
	"""
	l::Vector{Vector{Int}}
	r::Vector{Vector{Int}}
	sign::Vector{Int}
	function Partition(order)
		N = 2^order - 2
		function addincoming(l, r)
			# map 1, 2, 3, ... to order+1, order+2, ...
			l = l .+ order
			r = r .+ order
			nl, nr = length(l), length(r)
			inl = [i for i in 1:nl]
			inr = [i for i in nl+1:order]
			append!(inl, l)
			append!(inr, r)
			return inl, inr
		end

		left = Vector{Vector{Int}}(undef, N)
		right = similar(left)
		sign = Vector{Int}(undef, N)
		for (si, s) in enumerate(partitions(1:order, 2))
			# generate 2-partition of the outgoing legs
			# all incoming legs are kept the same
			s1, s2 = addincoming(s[1], s[2])
			left[2si-1], right[2si-1] = s1, s2
			sign[2si-1] = parity(vcat(s1, s2))

			q1, q2 = addincoming(s[2], s[1])
			left[2si], right[2si] = q1, q2
			sign[2si] = parity(vcat(q1, q2))
			# println("($s1, $s2) -> $p1, ($q1, $q2) -> $p2")
		end
		return new(left, right, sign)
	end
end