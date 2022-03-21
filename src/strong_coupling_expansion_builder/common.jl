import ..Filter
import ..Wirreducible  #remove all polarization subdiagrams
import ..Girreducible  #remove all self-energy inseration
import ..NoHatree
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
import ..GenericPara
import ..innerTauNum

import ..Diagram

import ..DiagramId
import ..Ver4Id
import ..Ver3Id
import ..GreenId
import ..SigmaId
import ..PolarId
import ..BareInteractionId
import ..BareGreenId
import ..BareGreenNId
import ..GreenNId
import ..ConnectedGreenNId

import ..TwoBodyChannel
import ..Alli
import ..PHr
import ..PHEr
import ..PPr
import ..AnyChan

import ..Permutation
import ..Di
import ..Ex
import ..DiEx

import ..uidreset
import ..toDataFrame
import ..mergeby

# function build(para::GenericPara, extT = nothing, subdiagram = false)
#     if para.diagType == GcDiag
#         if isnothing(extK)
#             extK = [DiagTree.getK(para.totalLoopNum, 1), DiagTree.getK(para.totalLoopNum, 2), DiagTree.getK(para.totalLoopNum, 3)]
#         end
#         return vertex4(para, extK, [PHr, PHEr, PPr, Alli], subdiagram)
#     elseif para.diagType == SigmaDiag
#         if isnothing(extK)
#             extK = DiagTree.getK(para.totalLoopNum, 1)
#         end
#         return sigma(para, extK, subdiagram)
#     elseif para.diagType == PolarDiag
#         if isnothing(extK)
#             extK = DiagTree.getK(para.totalLoopNum, 1)
#         end
#         return polarization(para, extK, subdiagram)
#     elseif para.diagType == Ver3Diag
#         if isnothing(extK)
#             extK = [DiagTree.getK(para.totalLoopNum, 1), DiagTree.getK(para.totalLoopNum, 2)]
#         end
#         return vertex3(para, extK, subdiagram)
#     else
#         error("not implemented!")
#     end
# end
