module Parquet
using StaticArrays, PyCall
using AbstractTrees
using Parameters
using ..DiagTree

import ..Filter
import ..Wirreducible  #remove all polarization subdiagrams
import ..Girreducible  #remove all self-energy inseration
import ..NoHatree
import ..NoFock
import ..NoBubble  # true to remove all bubble subdiagram
import ..Proper  #ver4, ver3, and polarization diagrams may require to be irreducible along the transfer momentum/frequency

const DI, EX = 1, 2
const INL, OUTL, INR, OUTR = 1, 2, 3, 4
# orginal diagrams T, U, S; particle-hole counterterm Ts, Us; and their counterterm Tc, Uc, Sc, Tsc, Usc 
# const I, T, U, S, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc = 1:12
@enum Channel I = 1 T U S Ts Us Ic Tc Uc Sc Tsc Usc
# const ChanName = ["I", "T", "U", "S", "Ts", "Us", "Ic", "Tc", "Uc", "Sc", "Tsc", "Usc"]
const SymFactor = [1.0, -1.0, 1.0, -0.5, +1.0, -1.0]

# const Fchan = [I, U, S, Ts, Us, Ic, Uc, Sc, Tsc, Usc]
# const Vchan = [I, T, U, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc]
const Allchan = [I, T, U, S, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc]


include("vertex4.jl")
include("eval.jl")
include("io.jl")
include("filter.jl")
include("ver4toDiagTree.jl")
include("sigma.jl")
end