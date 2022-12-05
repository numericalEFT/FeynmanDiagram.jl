mutable struct QuantumExpr <: AbstractVector{QuantumOperator}
    operators::Vector{QuantumOperator}
    function QuantumExpr(operators::Vector{QuantumOperator})
        return new(operators)
    end
    function QuantumExpr(operator::QuantumOperator)
        return new([operator,])
    end
    function QuantumExpr(operators::QuantumExpr)
        # return new([operators...,])
        return operators
    end
end

fermionic_annihilation(i) = QuantumExpr(QuantumOperator(:f⁻, i))
fermionic_creation(i) = QuantumExpr(QuantumOperator(:f⁺, i))
majorana(i) = QuantumExpr(QuantumOperator(:f, i))
bosonic_annihilation(i) = QuantumExpr(QuantumOperator(:b⁻, i))
bosonic_creation(i) = QuantumExpr(QuantumOperator(:b⁺, i))
real_classic(i) = QuantumExpr(QuantumOperator(:ϕ, i))

const 𝑓⁻ = fermionic_annihilation
const 𝑓⁺ = fermionic_creation
const 𝑓 = majorana
const 𝑏⁻ = bosonic_annihilation
const 𝑏⁺ = bosonic_creation
const 𝜙 = real_classic


#TODO: make compositeoperator norm ordered when it is created

Base.eltype(::Type{QuantumExpr}) = QuantumOperator
Base.getindex(o::QuantumExpr, i::Int) = o.operators[i]
Base.setindex!(o::QuantumExpr, v::QuantumOperator, i::Int) = o.operators[i] = v
Base.length(o::QuantumExpr) = length(o.operators)
Base.size(o::QuantumExpr) = size(o.operators)

Base.show(io::IO, o::QuantumExpr) = print(io, reduce(*, ["$o" for o in o.operators]))
Base.show(io::IO, ::MIME"text/plain", o::QuantumExpr) = Base.show(io, o)

function Base.:*(o1::QuantumOperator, o2::QuantumOperator)
    return QuantumExpr([o1, o2])
end

function Base.:*(o1::QuantumExpr, o2::QuantumOperator)
    return QuantumExpr([o1.operators; o2])
end

function Base.:*(o1::QuantumOperator, o2::QuantumExpr)
    return QuantumExpr([o1; o2.operators])
end

function Base.:*(o1::QuantumExpr, o2::QuantumExpr)
    return QuantumExpr([o1.operators; o2.operators])
end

function Base.adjoint(o::QuantumExpr)
    return QuantumExpr([op' for op in reverse(o)])
end

function isfermionic(o::QuantumExpr)
    numf = 0
    for op in o
        isfermionic(op) && (numf += 1)
    end
    isodd(numf) && return true
    return false
end

"""
    Converts a QuantumExpr to normal-ordered form in place and returns the associated statistical sign.
"""
function normal_order!(operator::QuantumExpr)
    sign = 1
    return sign
end

"""
    Computes the permutation required to convert a QuantumExpr to normal-ordered form. 
    Returns the associated statistical sign and permutation.
"""
function normal_order(operator::QuantumExpr)
    sign = 1
    permutation = collect(eachindex(operator.operators))
    return sign, permutation
end

function correlator_order!(operator::QuantumExpr)
    sign, ordering = correlator_order(operator)
    operator.operators = operator[sortperm(ordering)]
    return sign
end

function correlator_order(operator::QuantumExpr)
    num = length(operator)
    ind_pair, ind_unpair = 0, num + 1
    ordering = Int[]
    for (i, op) in enumerate(operator)
        if op' in operator[i+1:end]
            ind_pair += 1
            push!(ordering, !iscreation(op) ? ind_pair : num + 1 - ind_pair)
        elseif op' in operator[1:i-1]
            push!(ordering, num + 1 - ordering[findlast(isequal(op'), operator)])
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

    return sign, ordering
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