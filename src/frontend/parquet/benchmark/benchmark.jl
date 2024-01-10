module Benchmark
using StaticArrays, PyCall
using AbstractTrees
using Parameters, Combinatorics

import ..DiagPara
import ..interactionTauNum

const DI, EX = 1, 2
const INL, OUTL, INR, OUTR = 1, 2, 3, 4
# orginal diagrams T, U, S; particle-hole counterterm Ts, Us; and their counterterm Tc, Uc, Sc, Tsc, Usc 
const I, T, U, S, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc = 1:12
# const ChanName = ["I", "T", "U", "S", "Ts", "Us", "Ic", "Tc", "Uc", "Sc", "Tsc", "Usc"]
const SymFactor = [1.0, -1.0, 1.0, -0.5, +1.0, -1.0]

# const Fchan = [I, U, S, Ts, Us, Ic, Uc, Sc, Tsc, Usc]
# const Vchan = [I, T, U, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc]
# const Allchan = [I, T, U, S, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc]

include("vertex4.jl")
include("vertex4_eval.jl")
include("vertex4_io.jl")
include("diagram_count.jl")
end