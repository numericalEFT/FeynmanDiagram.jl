module Parquet
using StaticArrays, PyCall
using AbstractTrees
using ..DiagTree

const DI, EX = 1, 2
const INL, OUTL, INR, OUTR = 1, 2, 3, 4
# orginal diagrams T, U, S; particle-hole counterterm Ts, Us; and their counterterm Tc, Uc, Sc, Tsc, Usc 
const I, T, U, S, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc = 1:12
const ChanName = ["I", "T", "U", "S", "Ts", "Us", "Ic", "Tc", "Uc", "Sc", "Tsc", "Usc"]
const SymFactor = [1.0, -1.0, 1.0, -0.5, +1.0, -1.0]

const Fchan = [I, U, S, Ts, Us, Ic, Uc, Sc, Tsc, Usc]
const Vchan = [I, T, U, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc]
const Allchan = [I, T, U, S, Ts, Us, Ic, Tc, Uc, Sc, Tsc, Usc]

# function Base.getproperty(obj::GenericPara, sym::Symbol)
#     if sym === :internalSiteNum
#         return 
#     elseif sym === :n
#         return (obj.dim == 3) ? (obj.EF*2*obj.me)^(3/2)/(6π^2)*obj.spin : obj.me*obj.EF/π
#     elseif sym === :Rs
#         return (obj.dim == 3) ? (3 / (4π*obj.n))^(1 / 3) : sqrt(1/(πobj.n))
#     elseif sym === :a0
#         return 4π*obj.ϵ0/(obj.me*obj.e0^2)
#     elseif sym === :rs
#         return obj.Rs/obj.a0
#     elseif sym === :kF
#         return sqrt(2*obj.me*obj.EF)
#     else # fallback to getfield
#         return getfield(obj, sym)
#     end
# end


include("vertex4.jl")
include("io.jl")
include("diagtree.jl")
end