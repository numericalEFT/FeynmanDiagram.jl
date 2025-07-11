module SCE
using StaticArrays, PyCall
using AbstractTrees
using Parameters, Combinatorics
using DataFrames
import ..ComputationalGraphs as IR
import ..ComputationalGraphs: Graph
import ..ComputationalGraphs: Sum, Prod, Det
using ..FrontEnds: BareHoppingId, BareGreenNId, ConnectedGreenNId, GreenNId, GenericId, VacuumId, DiagramId

# const DI, EX, BOTH = 1, 2, 3
# const INL, OUTL, INR, OUTR = 1, 2, 3, 4
# orginal diagrams T, U, S; particle-hole counterterm Ts, Us; and their counterterm Tc, Uc, Sc, Tsc, Usc 
# symmetry factor for Alli, PHr, PHEr, PPr, PHrc, PHErc 

include("common.jl")
# include("Gn.jl")
# include("Gc.jl")
include("Gn_v1.jl")
include("Gc_v1.jl")

end