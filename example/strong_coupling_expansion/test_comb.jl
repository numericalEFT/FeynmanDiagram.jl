#=
To run this version of the script, you may need the Combinatorics.jl package.
You can install it by opening the Julia REPL and running:
import Pkg
Pkg.add("Combinatorics")
=#
using Combinatorics

#=
This script finds all sets of non-negative integer coefficients {b_1, b_2, ..., b_k}
such that the "weighted" sum equals a target integer `n`:
    1*b_1 + 2*b_2 + 3*b_3 + ... + k*b_k = n

This is equivalent to finding all partitions of the integer `n`.
For example, if n = 4, the partitions are:
- 4
- 3 + 1
- 2 + 2
- 2 + 1 + 1
- 1 + 1 + 1 + 1

These correspond to the coefficient sets:
- b_4=1  => {0, 0, 0, 1}
- b_3=1, b_1=1 => {1, 0, 1, 0}
- b_2=2  => {0, 2, 0, 0}
- b_2=1, b_1=2 => {2, 1, 0, 0}
- b_1=4  => {4, 0, 0, 0}
Note: The sets are represented sparsely, with the vector length being the max part.
=#

"""
    find_integer_compositions(n::Int)

Finds all sets of coefficients {b_i} such that ∑ i*b_i = n, using the
Combinatorics.jl package. This is the main function for the calculation.

# Arguments
- `n::Int`: The positive integer to be composed.

# Returns
- `Vector{Vector{Int}}`: A vector of vectors, where each inner vector
  represents a set of coefficients {b_1, b_2, ...}.
"""
function find_integer_compositions(n::Int)
    if n <= 0
        println("Input must be a positive integer.")
        return Vector{Vector{Int}}[]
    end

    # results = Vector{Vector{Int}}[]
    results = Vector{Int}[]
    # The partitions(n) function from Combinatorics.jl returns an iterator
    # for all partitions of the integer n.
    for p in partitions(n)
        # Convert the partition format (e.g., [4, 1, 1]) to the
        # desired coefficient format (e.g., b_1=2, b_4=1 => [2, 0, 0, 1])
        push!(results, partition_to_coeffs(p))
    end
    return results
end


"""
    partition_to_coeffs(partition::Vector{Int})

Converts a partition representation into a coefficient `b_i` representation.

For example, the partition [4, 2, 2, 1] for n=9 means we used one 4, two 2s,
and one 1. This corresponds to b_1=1, b_2=2, b_4=1. The function would
return the vector [1, 2, 0, 1].

# Arguments
- `partition::Vector{Int}`: A vector representing a partition, e.g., [3, 1].

# Returns
- `Vector{Int}`: The corresponding coefficient vector {b_1, b_2, ...}.
"""
function partition_to_coeffs(partition::Vector{Int})
    if isempty(partition)
        return Int[]
    end
    # The size of the coefficient vector is determined by the largest number
    # in the partition.
    max_part = maximum(partition)
    coeffs = zeros(Int, max_part)

    for part in partition
        coeffs[part] += 1
    end
    return coeffs
end


# --- Main Execution ---
function main()
    println("--- Integer Composition Solver ---")
    print("Enter a positive integer n: ")
    try
        n_str = readline()
        n = parse(Int, n_str)

        if n > 0
            println("\nFinding compositions for n = $n")
            println("such that 1*b₁ + 2*b₂ + 3*b₃ + ... = $n")
            println("-"^30)

            all_compositions = find_integer_compositions(n)

            println(all_compositions)

            if isempty(all_compositions)
                println("No compositions found.")
            else
                println("Found $(length(all_compositions)) sets of coefficients {b₁, b₂, ...}:")
                # Sort results for consistent output, making it easier to read.
                # Sorting by length first, then by content.
                sort!(all_compositions, by=x -> (length(x), x))
                for (i, coeffs) in enumerate(all_compositions)
                    print("  Set $i: {")
                    print(join(coeffs, ", "))
                    println("}")
                end
            end
        else
            println("Please enter a positive integer.")
        end
    catch e
        if isa(e, ArgumentError)
            println("Invalid input. Please enter an integer.")
        else
            rethrow(e)
        end
    end
end

# Run the main function
main()


# --- Original Recursive Implementation (for reference) ---

"""
    find_integer_compositions_recursive(n::Int)

(Reference Implementation) Finds all sets of coefficients {b_i} such that ∑ i*b_i = n.
"""
function find_integer_compositions_recursive(n::Int)
    if n <= 0
        println("Input must be a positive integer.")
        return Vector{Vector{Int}}[]
    end
    results = Vector{Vector{Int}}[]
    find_compositions_recursive_helper(n, n, [], results)
    return results
end

"""
    find_compositions_recursive_helper(target, max_val, current_partition, results)

(Reference Implementation) A recursive helper function to find integer compositions.
"""
function find_compositions_recursive_helper(target::Int, max_val::Int, current_partition::Vector{Int}, results::Vector{Vector{Int}})
    if target == 0
        push!(results, partition_to_coeffs(current_partition))
        return
    end

    if max_val < 1 || max_val > target
        return
    end

    for i in min(max_val, target):-1:1
        # CORRECTED LINE: Use the splat operator '...' to create a new vector.
        # [current_partition; i] is invalid syntax for this operation.
        new_partition = [current_partition..., i]
        find_compositions_recursive_helper(target - i, i, new_partition, results)
    end
end
