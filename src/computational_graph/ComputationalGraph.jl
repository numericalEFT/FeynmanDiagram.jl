module ComputationalGraph
using AbstractTrees
using Printf, PyCall, DataFrames, Parameters

@enum TwoBodyChannel Alli = 1 PHr PHEr PPr AnyChan
@enum Permutation Di = 1 Ex DiEx

# export TwoBodyChannel, Alli, PHr, PHEr, PPr, AnyChan
# export Permutation, Di, Ex, DiEx

Base.length(r::TwoBodyChannel) = 1
Base.iterate(r::TwoBodyChannel) = (r, nothing)
function Base.iterate(r::TwoBodyChannel, ::Nothing) end

Base.length(r::Permutation) = 1
Base.iterate(r::Permutation) = (r, nothing)
function Base.iterate(r::Permutation, ::Permutation) end

include("common.jl")
include("GreenDiagram.jl")
# include("tree.jl")
# include("operation.jl")
# include("io.jl")
# include("eval.jl")
# include("optimize.jl")


const INL, OUTL, INR, OUTR = 1, 2, 3, 4

export GreenDiagram, ExternalVertice
# export addSubDiagram!
# export evalDiagTree!
# export evalDiagTreeKT!
# export Operator, Sum, Prod
# export uidreset
# export toDataFrame, mergeby, plot_tree

end