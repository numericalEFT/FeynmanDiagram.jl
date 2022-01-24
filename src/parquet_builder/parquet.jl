module Parquet
using StaticArrays, PyCall
using AbstractTrees
using Parameters, Combinatorics
using DataFrames
using ..DiagTree


const DI, EX, BOTH = 1, 2, 3
const INL, OUTL, INR, OUTR = 1, 2, 3, 4
# orginal diagrams T, U, S; particle-hole counterterm Ts, Us; and their counterterm Tc, Uc, Sc, Tsc, Usc 
# symmetry factor for Alli, PHr, PHEr, PPr, PHrc, PHErc 
const SymFactor = [1.0, -1.0, 1.0, -0.5, +1.0, -1.0]

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