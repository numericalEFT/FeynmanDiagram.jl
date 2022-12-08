"""
    struct OperatorProduct <: AbstractVector{QuantumOperator}

    struct of a quantum-operator product. It is a subtype of `AbstractVector` and
    inherits a large set of Vector behaviors including iteration and indexing. 

# Members:
- `operators::Vector{QuantumOperator}`  vector of quantum operators
"""
struct OperatorProduct <: AbstractVector{QuantumOperator}
    operators::Vector{QuantumOperator}
    function OperatorProduct(products::Vector{OperatorProduct})
        return new(union(products...))
    end
    function OperatorProduct(operators::Vector{QuantumOperator})
        return new(operators)
    end
    function OperatorProduct(operator::QuantumOperator)
        return new([operator,])
    end
    function OperatorProduct(operators::OperatorProduct)
        # return new([operators...,])
        return operators
    end
end

"""
    Create a OperatorProduct with one quantum operator from given label `i`.
    It supports the following abbreviated function form:

'''
const 𝑓⁻ = fermionic_annihilation
const 𝑓⁺ = fermionic_creation
const 𝑓 = majorana
const 𝑏⁻ = bosonic_annihilation
const 𝑏⁺ = bosonic_creation
const 𝜙 = real_classic
'''
"""
fermionic_annihilation(i) = OperatorProduct(QuantumOperator(:f⁻, i))
fermionic_creation(i) = OperatorProduct(QuantumOperator(:f⁺, i))
majorana(i) = OperatorProduct(QuantumOperator(:f, i))
bosonic_annihilation(i) = OperatorProduct(QuantumOperator(:b⁻, i))
bosonic_creation(i) = OperatorProduct(QuantumOperator(:b⁺, i))
real_classic(i) = OperatorProduct(QuantumOperator(:ϕ, i))
const 𝑓⁻ = fermionic_annihilation
const 𝑓⁺ = fermionic_creation
const 𝑓 = majorana
const 𝑏⁻ = bosonic_annihilation
const 𝑏⁺ = bosonic_creation
const 𝜙 = real_classic

Base.eltype(::Type{OperatorProduct}) = QuantumOperator
Base.getindex(o::OperatorProduct, i::Int) = o.operators[i]
Base.setindex!(o::OperatorProduct, v::QuantumOperator, i::Int) = o.operators[i] = v
Base.length(o::OperatorProduct) = length(o.operators)
Base.size(o::OperatorProduct) = size(o.operators)

Base.show(io::IO, o::OperatorProduct) = print(io, reduce(*, ["$o" for o in o.operators]))
Base.show(io::IO, ::MIME"text/plain", o::OperatorProduct) = Base.show(io, o)

"""
    Base.:*(o1::Union{QuantumOperator, OperatorProduct}, o2::Union{QuantumOperator, OperatorProduct})

    `o1 * o2` returns the quantum operator product of `o1` and `o2`
"""
function Base.:*(o1::QuantumOperator, o2::QuantumOperator)
    return OperatorProduct([o1, o2])
end

function Base.:*(o1::OperatorProduct, o2::QuantumOperator)
    return OperatorProduct([o1.operators; o2])
end

function Base.:*(o1::QuantumOperator, o2::OperatorProduct)
    return OperatorProduct([o1; o2.operators])
end

function Base.:*(o1::OperatorProduct, o2::OperatorProduct)
    return OperatorProduct([o1.operators; o2.operators])
end

"""
    Base.adjoint(o::OperatorProduct)

    Return the conjuated composite operator of `o`.
"""
function Base.adjoint(o::OperatorProduct)
    return OperatorProduct([op' for op in reverse(o)])
end

"""
    function isfermionic(o::OperatorProduct)

    Check if `o` is a fermionic composite operator.
"""
function isfermionic(o::OperatorProduct)
    numf = 0
    for op in o
        isfermionic(op) && (numf += 1)
    end
    isodd(numf) && return true
    return false
end

"""
    function normal_order(operator::OperatorProduct)

    Computes the permutation required to convert a OperatorProduct to normal-ordered form. 
    Returns the associated statistical sign and permutation.
"""
function normal_order(operator::OperatorProduct)
    sign = 1
    permutation = collect(eachindex(operator.operators))
    return sign, permutation
end

"""
    function normal_order(operator::OperatorProduct)

    Convert a OperatorProduct to correlator-ordered form. 
    Returns the associated statistical sign and ordered OperatorProduct.
"""
function correlator_order(operator::OperatorProduct)
    num = length(operator)
    ind_pair, ind_unpair = 0, num + 1
    ordering = Int[]
    for (i, op) in enumerate(operator)
        if op' in operator[i+1:end]
            ind_pair += 1
            push!(ordering, !iscreation(op) ? ind_pair : num + 1 - ind_pair)
        elseif op' in operator[1:i-1]
            push!(ordering, num + 1 - ordering[findlast(isequal(op'), operator[1:i-1])])
        else
            push!(ordering, !iscreation(op) ? ind_unpair : -ind_unpair)
        end
    end
    ind_ann, ind_cre = 0, 0
    for (i, value) in enumerate(ordering)
        if value == ind_unpair
            ind_ann += 1
            ordering[i] = ind_pair + ind_ann
        elseif value == -ind_unpair
            ind_cre += 1
            ordering[i] = num + 1 - ind_pair - ind_cre
        end
    end

    permutation = ordering[isfermionic.(operator)]
    sign = isempty(permutation) ? 1 : parity(sortperm(permutation))

    return sign, operator[sortperm(ordering)]
end

"""
The parity of a permutation P is +1 if the number of 2-cycles (swaps) in an n-cycle
decomposition with n ≤ 2 is even, and -1 if the number of 2-cycles is odd.
"""
function parity(p::W) where {W<:AbstractVector{Int}}
    count = 0
    p_swap = copy(p)
    for i in eachindex(p)
        while p_swap[i] != i
            count += 1
            j = p_swap[i]  # NOTE: we cannot swap in place with nested access in Julia
            p_swap[i], p_swap[j] = p_swap[j], p_swap[i]
        end
    end
    return 2 * iseven(count) - 1
end

"""
calculate the parity of a given permutation of the array [1, 2, 3, ...]
"""
function parity_old(p)
    n = length(p)
    not_seen = Set{Int}(1:n)
    seen = Set{Int}()
    cycles = Array{Int,1}[]
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