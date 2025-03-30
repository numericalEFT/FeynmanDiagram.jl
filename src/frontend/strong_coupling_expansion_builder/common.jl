
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
