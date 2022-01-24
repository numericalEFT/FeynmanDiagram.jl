module Parquet
using StaticArrays, PyCall
using AbstractTrees
using Parameters, Combinatorics
using DataFrames
using ..DiagTree
using ..DiagTreeNew


const DI, EX, BOTH = 1, 2, 3
const INL, OUTL, INR, OUTR = 1, 2, 3, 4
# orginal diagrams T, U, S; particle-hole counterterm Ts, Us; and their counterterm Tc, Uc, Sc, Tsc, Usc 
# const I, T, U, S, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc = 1:12
@enum Channel I = 1 T U S Ts Us Ic Tc Uc Sc Tsc Usc

Base.length(r::Channel) = 1
Base.iterate(r::Channel) = (r, nothing)
function Base.iterate(r::Channel, ::Nothing) end

# const ChanName = ["I", "T", "U", "S", "Ts", "Us", "Ic", "Tc", "Uc", "Sc", "Tsc", "Usc"]
const SymFactor = [1.0, -1.0, 1.0, -0.5, +1.0, -1.0]

# const Fchan = [I, U, S, Ts, Us, Ic, Uc, Sc, Tsc, Usc]
# const Vchan = [I, T, U, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc]
# const Allchan = [I, T, U, S, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc]

include("common.jl")
export ParquetBlocks

include("filter.jl")
include("vertex4.jl")

include("sigma.jl")
include("green.jl")
include("vertex3.jl")
include("polarization.jl")

include("benchmark/benchmark.jl")
end