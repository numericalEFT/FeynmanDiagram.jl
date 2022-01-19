module DiagTreeNew
using AbstractTrees
using Printf, PyCall, DataFrames

uid() = abs(rand(Int)) % 10000

include("common.jl")
include("traits.jl")
include("diagram.jl")
include("diagram_io.jl")
export Diagram
export addSubDiagram!
export toDataFrame
export evalDiagNode!
export evalDiagTree!
export Operator, Sum, Prod
export DiagramId, Ver4Id, Ver3Id, GreenId, SigmaId, PolarId, InteractionId

# struct Cache and struct Pool
# include("pool.jl")

# struct Propagator, Node and Diagram
# include("object.jl")

# # IO operations
# include("io.jl")

# # diagram evaluation
# include("eval.jl")

end