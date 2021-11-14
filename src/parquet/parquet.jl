module Parquet
using StaticArrays, PyCall
using AbstractTrees

const INL, OUTL, INR, OUTR = 1, 2, 3, 4
# orginal diagrams T, U, S; particle-hole counterterm Ts, Us; and their counterterm Tc, Uc, Sc, Tsc, Usc 
const I, T, U, S, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc = 1:12
const ChanName = ["I", "T", "U", "S", "Ts", "Us", "Ic", "Tc", "Uc", "Sc", "Tsc", "Usc"]

const Fchan = [I, U, S, Ts, Us, Ic, Uc, Sc, Tsc, Usc]
const Vchan = [I, T, U, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc]
const Allchan = [I, T, U, S, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc]

include("vertex4.jl")
include("io.jl")
end