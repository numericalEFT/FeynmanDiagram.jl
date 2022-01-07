module Builder
using Parameters

# const DI, EX = 1, 2
# const INL, OUTL, INR, OUTR = 1, 2, 3, 4
# # orginal diagrams T, U, S; particle-hole counterterm Ts, Us; and their counterterm Tc, Uc, Sc, Tsc, Usc 
# # const I, T, U, S, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc = 1:12
# @enum Channel I, T, U, S, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc
# const ChanName = ["I", "T", "U", "S", "Ts", "Us", "Ic", "Tc", "Uc", "Sc", "Tsc", "Usc"]
# const SymFactor = [1.0, -1.0, 1.0, -0.5, +1.0, -1.0]

# # const Fchan = [I, U, S, Ts, Us, Ic, Uc, Sc, Tsc, Usc]
# # const Vchan = [I, T, U, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc]
# const Allchan = [I, T, U, S, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc]

# @enum DiagramType begin
#     SigmaDiag          #self-energy
#     GreenDiag          #green's function
#     PolarDiag          #polarization
#     Ver3Diag           #3-point vertex function
#     Ver4Diag           #4-point vertex function
# end

# @enum Filter begin
#     Wirreducible  #remove all polarization subdiagrams
#     Girreducible  #remove all self-energy inseration
#     NoHatree
#     NoFock
#     NoBubble  # true to remove all bubble subdiagram
#     Proper  #ver4, ver3, and polarization diagrams may require to be irreducible along the transfer momentum/frequency
# end

# @enum InteractionName begin
#     Composite
#     ChargeCharge
#     SpinSpin
#     ProperChargeCharge
#     ProperSpinSpin
# end

# @enum InteractionType begin
#     Instant
#     Dynamic
#     D_Instant #derivative of instant interaction
#     D_Dynamic #derivative of the dynamic interaction
# end

# export SigmaDiag, PolarDiag, Ver3Diag, Ver4Diag
# export Wirreducible, Girreducible, NoBubble, NoHatree, Proper
# export InteractionName, ChargeCharge, SpinSpin
# export InteractionType, Instant, Dynamic, D_Instant, D_Dynamic

# include("parameter.jl")
# export GenericPara, Interaction

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

import ..Composite
import ..ChargeCharge
import ..SpinSpin
import ..ResponseName

import ..Instant
import ..Dynamic
import ..AnalyticProperty

import ..Interaction

import ..GenericPara

import ..innerTauNum

using ..DiagTree

include("parquetAlg/parquet.jl")
export Parquet

# include("newparquetAlg/parquet.jl")

end