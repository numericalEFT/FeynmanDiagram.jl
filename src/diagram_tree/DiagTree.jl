module DiagTreeNew
using AbstractTrees
using Printf, PyCall, DataFrames

include("common.jl")
include("traits.jl")
include("diagram.jl")
include("diagram_io.jl")
export DiagramId
export Diagram
export addSubDiagram!
export toDataFrame
export evalDiagNode!
export evalDiagTree!
export Operator, Sum, Prod

# struct Cache and struct Pool
# include("pool.jl")

# struct Propagator, Node and Diagram
# include("object.jl")

# # IO operations
# include("io.jl")

# # diagram evaluation
# include("eval.jl")

end