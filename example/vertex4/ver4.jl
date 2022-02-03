using Printf, LinearAlgebra
using CompositeGrids
using ElectronGas
using Parameters, Random, DataFrames
using MCIntegration
using Lehmann

using FeynmanDiagram
using StaticArrays

include("parameter.jl")
include("interaction.jl")


const steps = 1e6
const Order = 1
const isF = false

const RefK = [kF, 0.0, 0.0]
const qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 6 * kF], [0.0, 2kF], 16, 0.01 * kF, 8)
const τgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 16, β * 1e-4, 8)
const lgrid = [1, 2]
const Nl = length(lgrid)
const Nk = 16
const θgrid = collect(LinRange(0.1, π, Nk)) # external angle grid
const ExtK = [[kF * cos(θ), kF * sin(θ), 0.0] for θ in θgrid]

vqinv = [(q^2 + mass2) / (4π * e0^2) for q in qgrid.grid]
const dW0 = TwoPoint.dWRPA(vqinv, qgrid.grid, τgrid.grid, dim, EF, kF, β, spin, me) # dynamic part of the effective interaction

KinL = KoutL = [1.0, 0, 0]
KinR = KoutR = [0, 1.0, 0]
legK = [KinL, KoutL, KinR, KoutR]

diagPara(order) = GenericPara(diagType = Ver4Diag, innerLoopNum = order, hasTau = true, loopDim = dim, spin = spin, firstLoopIdx = 3,
    interaction = [FeynmanDiagram.Interaction(ChargeCharge, [Instant, Dynamic]),],  #instant charge-charge interaction
    filter = [
        Girreducible,
        Proper,   #one interaction irreduble diagrams or not
        # NoBubble, #allow the bubble diagram or not
    ],
    transferLoop = KinL - KoutL
)

println("Build the diagrams into an experssion tree ...")
const para = [diagPara(o) for o in 1:Order]
ver4 = [Parquet.vertex4(para[i], legK, [PHr,]) for i in 1:Order]   #diagram of different orders
#different order has different set of K, T variables, thus must have different exprtrees
println(ver4)
"""
DataFrame[6×5 DataFrame
 Row │ response  type     extT          diagram                            hash  
     │ Any       Any      Any           Diagram…                           Int64 
─────┼───────────────────────────────────────────────────────────────────────────
   1 │ UpUp      Dynamic  (1, 1, 2, 2)  ↑↑Dyn,t(1, 1, 2, 2)=0.0=⨁ (17, 2…     81
   2 │ UpUp      Dynamic  (1, 2, 1, 2)  ↑↑x↑↑ → PPr, PPr ↑↑Dyn,t(1, 2, 1…     71
   3 │ UpUp      Dynamic  (1, 2, 2, 1)  ↑↑Dyn,t(1, 2, 2, 1)=0.0=⨁ (43, 5…     82
   4 │ UpDown    Dynamic  (1, 1, 2, 2)  ↑↓Dyn,t(1, 1, 2, 2)=0.0=⨁ (20, 2…     83
   5 │ UpDown    Dynamic  (1, 2, 1, 2)  ↑↓Dyn,t(1, 2, 1, 2)=0.0=⨁ (74, 7…     84
   6 │ UpDown    Dynamic  (1, 2, 2, 1)  ↑↓Dyn,t(1, 2, 2, 1)=0.0=⨁ (44, 4…     85]
"""
# plot_tree(ver4uu[1][1])
const diag = [ExprTree.build(ver4[o].diagram) for o in 1:Order]    #experssion tree representation of diagrams 
const rootuu = [[idx for idx in d.root if d.nodePool.object[idx].para.response == UpUp] for d in diag]
const rootud = [[idx for idx in d.root if d.nodePool.object[idx].para.response == UpDown] for d in diag]
const extTuu = [[diag[ri].nodePool.object[idx].para.extT for idx in root] for (ri, root) in enumerate(rootuu)]
const extTud = [[diag[ri].nodePool.object[idx].para.extT for idx in root] for (ri, root) in enumerate(rootud)]
println(rootuu)
println(extTuu)
println(rootud)
println(extTud)
ExprTree.showTree(diag[1], rootuu[1][1])
# ExprTree.showTree(diag[1], rootud[1][1])

# exit(0)

##################### propagator and interaction evaluation ##############
function eval(id::GreenId, K, varT)
    τin, τout = varT[id.extT[1]], varT[id.extT[2]]
    ϵ = dot(K, K) / (2me) - μ
    τ = τout - τin
    if τ ≈ 0.0
        return Spectral.kernelFermiT(-1e-8, ϵ, β)
    else
        return Spectral.kernelFermiT(τ, ϵ, β)
    end
end

eval(id::InteractionId, K, varT) = e0^2 / ϵ0 / (dot(K, K) + mass2)

@inline function phase(varT, extT)
    # println(extT)
    tInL, tOutL, tInR, tOutR = varT[extT[INL]], varT[extT[OUTL]], varT[extT[INR]],
    varT[extT[OUTR]]
    if (isF)
        return cos(π / β * ((tInL + tOutL) - (tInR + tOutR)))
    else
        return cos(π / β * ((tInL - tOutL) + (tInR - tOutR)))
    end
end

function integrand(config)
    order = config.curr
    x = config.var[3][1]
    varK, varT = config.var[1], config.var[2]
    varK.data[:, 2] .= ExtK[x]
    # @assert varK.data[:, 1] ≈ [kF, 0.0, 0.0]
    w = ExprTree.evalNaive!(diag[order], varK.data, varT.data, eval)
    # println(rootuu)
    # for (ri, root) in enumerate(rootuu)
    #     println(ri, ",  ", root)
    #     println(diag[order].nodePool.current[root] * phase(varT, extTuu[order][ri]))
    # end
    # @assert w[1] ≈ diag[order].nodePool.current[13]
    # @assert w[2] ≈ diag[order].nodePool.current[14]
    if !isempty(rootuu[order])
        wuu = sum(diag[order].nodePool.current[root] * phase(varT, extTuu[order][ri]) for (ri, root) in enumerate(rootuu[order]))
    else
        wuu = 0.0
    end
    if !isempty(rootud[order])
        wud = sum(diag[order].nodePool.current[root] * phase(varT, extTud[order][ri]) for (ri, root) in enumerate(rootud[order]))
    else
        wud = 0.0
    end
    # println(wuu, ",  ", wud)
    return Weight(wuu / β, wud / β)
end

function measure(config)
    factor = 1.0 / config.reweight[config.curr]
    x = config.var[3][1]
    # println(config.observable[1][1])
    if config.curr == 1
        weight = integrand(config)
        config.observable[x, 1] += weight.d / abs(weight) * factor
        config.observable[x, 2] += weight.e / abs(weight) * factor
    else
        return
    end
end

function MC()
    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset = 2)
    K.data[:, 1] .= RefK
    T = MCIntegration.Tau(β, β / 2.0)
    # X = MCIntegration.Continuous([-1.0, 1.0], 0.2) #x=cos(θ)
    X = MCIntegration.Discrete(1, Nk)

    # for (ti, t) in enumerate(T.data)
    #     t[1] = β * rand()
    #     t[2] = β * rand()
    # end

    dof = [[1, 2, 1],] # K, T, ExtKidx
    obs = zeros(Nk, 2) # observable for the Fock diagram 

    config = MCIntegration.Configuration(steps, (K, T, X), dof, obs)
    avg, std = MCIntegration.sample(config, integrand, measure; print = 0, Nblock = 16)

    function info(idx, di)
        return @sprintf("   %8.4f ±%8.4f", avg[idx, di], std[idx, di])
    end

    if isnothing(avg) == false
        avg *= NF
        std *= NF
        N = size(avg)[1]
        grid = θgrid
        println("UpUp ver4: ")
        for li in 1:N
            @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], avg[li, 1], std[li, 1])
        end
        println("UpDown ver4: ")
        for li in 1:N
            @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], avg[li, 2], std[li, 2])
        end

        println("S ver4: ")
        for li in 1:N
            @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], (avg[li, 1] + avg[li, 2]), (std[li, 1] + std[li, 2]))
        end
        println("A ver4: ")
        for li in 1:N
            @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], (avg[li, 1] - avg[li, 2]), (std[li, 1] - std[li, 2]))
        end

    end

end

MC()