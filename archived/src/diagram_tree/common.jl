import ..Filter
import ..Wirreducible  #remove all polarization subdiagrams
import ..Girreducible  #remove all self-energy inseration
import ..NoHartree
import ..NoFock
import ..NoBubble  # true to remove all bubble subdiagram
import ..Proper  #ver4, ver3, and polarization diagrams may require to be irreducible along the transfer momentum/frequency

import ..DiagramType
import ..GreenDiag
import ..SigmaDiag
import ..PolarDiag
import ..Ver3Diag
import ..Ver4Diag
import ..GnDiag
import ..GcDiag

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

import ..DiagPara

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

function uid()
    _counter[1] += 1
    return _counter[1]
end

function uidreset()
    _counter[1] = 0
end

function getK(loopNum::Int, loopIdx::Int)
    k = zeros(loopNum)
    k[loopIdx] = 1.0
    return k
end
