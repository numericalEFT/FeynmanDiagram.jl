module Manual
using StaticArrays, PyCall
using AbstractTrees
using ..DiagTree

const INL, OUTL, INR, OUTR = 1, 2, 3, 4
# orginal diagrams T, U, S; particle-hole counterterm Ts, Us; and their counterterm Tc, Uc, Sc, Tsc, Usc 
const I, T, U, S, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc = 1:12
const ChanName = ["I", "T", "U", "S", "Ts", "Us", "Ic", "Tc", "Uc", "Sc", "Tsc", "Usc"]
const SymFactor = [1.0, -1.0, 1.0, -0.5, +1.0, -1.0]

include("vertex4.jl")
# include("io.jl")
end