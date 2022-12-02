import ..Filter
# import ..Wirreducible  #remove all polarization subdiagrams
# import ..Girreducible  #remove all self-energy inseration
import ..NoHartree
import ..NoFock
import ..NoBubble  # true to remove all bubble subdiagram
import ..Proper  #ver4, ver3, and polarization diagrams may require to be irreducible along the transfer momentum/frequency

# import ..DiagType
# import ..Vacuum
# import ..Tadpole
# import ..FermiPropagator
# import ..BosePropagator
# import ..FermiSelfEnergy
# import ..BoseSelfEnergy
# import ..VertexDiag
# import ..GncDiag
# import ..GndDiag

import ..Composite
import ..ChargeCharge
import ..SpinSpin
import ..UpUp
import ..UpDown
import ..Response

import ..Instant
import ..Dynamic
import ..D_Instant
import ..D_Dynamic
import ..AnalyticProperty

import ..symbol
import ..short

import ..Interaction

import ..DiagramPara

import ..innerTauNum

# unique id
# uid() = abs(rand(Int)) % 10000

# let z = 0
#     global function uid()
#         z += 1
#         return z
#     end
# end
const _counter = [0,]

mutable struct DType
    factor::DataType
    weight::DataType
end

const _dtype = DType(Float64, Float64)

function set_datatype(; factor=Float64, weight=Float64)
    _dtype.factor = factor
    _dtype.weight = weight
end

function uid()
    _counter[1] += 1
    return _counter[1]
end

function uidreset()
    _counter[1] = 0
end

const _labelcounter = [0,]

function label()
    _labelcounter[1] += 1
    return _labelcounter[1]
end

function labelreset()
    _labelcounter[1] = 0
end

function getK(loopNum::Int, loopIdx::Int)
    k = zeros(loopNum)
    k[loopIdx] = 1.0
    return k
end

"""
The parity of a permutation P of length n is even if P has an even number of even-length cycles, 
and odd otherwise, where a cycle's parity is -1 if it is even in length. 

The total parity of P is then given by: sgn(P) = (-1)^(n - n_cycles)
"""
function parity(p::W) where {W<:AbstractVector{Int}}
    n = length(p)
    visited = Set{Int}()
    unvisited = Set{Int}(1:n)
    # Continue searching until we visit every item
    n_cycles = 0
    while isempty(unvisited) == false
        cycle = Int[]
        x = pop!(unvisited)
        # Continue until we complete the cycle
        while in(x, visited) == false
            push!(cycle, x)
            push!(visited, x)
            x = p[x]
            pop!(unvisited, x, 0)
        end
        # Found a cycle
        n_cycles += 1
    end
    return (-1)^(n - n_cycles)
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
