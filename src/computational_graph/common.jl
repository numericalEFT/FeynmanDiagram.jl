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

