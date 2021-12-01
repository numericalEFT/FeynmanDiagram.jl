using Printf, LinearAlgebra
using CompositeGrids
using ElectronGas
using Parameters, Random
using MCIntegration
using Lehmann

using ExpressionTree
using StaticArrays

Parquet = GWKT.Parquet
Manual = GWKT.Manual
DiagTree = GWKT.DiagTree

include("parameter.jl")
include("interaction.jl")


const steps = 1e7
const isF = false
const isProper = true #one interaction irreduble diagrams or not
const hasBubble = false #allow the bubble diagram or not
const Nk = 8
const θgrid = collect(LinRange(0.0, π, Nk)) # external angle grid
const ExtK = [@SVector [kF * cos(θ), kF * sin(θ), 0.0] for θ in θgrid]
const RefK = @SVector [kF, 0.0, 0.0]

struct Para{Q,T}
    dW0::Matrix{Float64}
    qgrid::Q
    τgrid::T # dedicated τgrid for dynamic interaction
    diag::Vector{DiagTree.Diagrams{Float64}}
    dir::Vector{Int}
    ex::Vector{Int}
    function Para(loopOrder::Int)
        qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 6 * kF], [0.0, 2kF], 16, 0.01 * kF, 8)
        τgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 16, β * 1e-4, 8)

        vqinv = [(q^2 + mass2) / (4π * e0^2) for q in qgrid.grid]
        dW0 = TwoPoint.dWRPA(vqinv, qgrid.grid, τgrid.grid, dim, EF, kF, β, spin, me) # dynamic part of the effective interaction

        # chan = [1, 2, 3]
        KinL = KoutL = [1, 0, 0]
        KinR = KoutR = [0, 1, 0]
        legK = [KinL, KoutL, KinR, KoutR]
        Gsym = [:mirror]
        Wsym = [:mirror, :timereversal]



        diag1, dir1, ex1 = Manual.build([1, 2, 3], legK, 3, spin, isProper, hasBubble, Gsym, Wsym)
        diag2, dir2, ex2 = Manual.build([2,], legK, 3, spin, isProper, hasBubble, Gsym, Wsym)
        diag3, dir3, ex3 = Manual.build([3,], legK, 3, spin, isProper, hasBubble, Gsym, Wsym)
        diag = [diag1, diag2, diag3]
        dir = [dir1, dir2, dir3]
        ex = [ex1, ex2, ex3]
        # DiagTree.showTree(diag1, dir1)
        # DiagTree.showTree(diag3, dir3)

        return new{typeof(qgrid),typeof(τgrid)}(dW0, qgrid, τgrid, diag, dir, ex)
    end
end

function integrand(config, _isF = isF)
    if config.curr == 1
        extKidx = config.var[3][1]
        KInL = RefK
        KOutL = ExtK[extKidx]
        return eval(config, KInL, KOutL, _isF)
    else
        error("impossible!")
    end
end

function phaseT(tInL, tOutL, _isF)
    if (_isF)
        return cos(π * (tInL + tOutL) / β)
    else
        return cos(π * (tInL - tOutL) / β)
    end
end

function eval(config, KInL, KOutL, _isF)
    para = config.para
    K, T = config.var[1], config.var[2]
    Qd = KInL - KOutL
    k1, k2 = K[1], K[1] - Qd
    t1, t2 = T[1], T[2] # t1, t2 both have two tau variables

    vle, wle = interactionDynamic(para, Vector(KInL) - k1, t1, t2)

    # possible green's functions on the top
    ϵ1, ϵ2 = (dot(k1, k1) - kF^2) / (2me), (dot(k2, k2) - kF^2) / (2me)

    g1 = Spectral.kernelFermiT(0.0 - t1, ϵ1, β)
    g2 = Spectral.kernelFermiT(t1 - 0.0, ϵ2, β)
    g3 = Spectral.kernelFermiT(t2 - 0.0, ϵ2, β)

    w = 0.0
    w += g1 * g2 * vle / (2π)^dim * phaseT(t1, t1, _isF)
    w += g1 * g3 * wle / (2π)^dim * phaseT(t1, t2, _isF)

    return w
end

function measure(config)
    factor = 1.0 / config.reweight[config.curr]
    extKidx = config.var[3][1]
    # println(config.observable[1][1])
    if config.curr == 1
        weight = integrand(config)
        wF = integrand(config, true)
        config.observable[extKidx, 1] += wF / abs(weight) * factor
        wA = integrand(config, false)
        config.observable[extKidx, 2] += wA / abs(weight) * factor
    else
        return
    end
end

function l0(data)
    avg = 0.0
    dθ = (θgrid[2] - θgrid[1]) / π
    avg += data[1] * sin(θgrid[1]) * dθ / 2
    dθ = (θgrid[end] - θgrid[end-1]) / π
    avg += data[end] * sin(θgrid[end]) * dθ / 2
    for i = 2:length(θgrid)-1
        dθ = (θgrid[i+1] - θgrid[i]) / π
        avg += data[i] * sin(θgrid[i]) * dθ
    end
    avg /= 2 # normalization factor
    return avg
end

function MC()
    para = Para(1)

    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF)
    T = MCIntegration.Tau(β, β / 2.0)
    ExtKidx = MCIntegration.Discrete(1, Nk)

    dof = [[1, 2, 1],] # K, T, ExtKidx
    obs = zeros(Nk, 2) # observable for the Fock diagram 

    config = MCIntegration.Configuration(steps, (K, T, ExtKidx), dof, obs; para = para)
    avg, std = MCIntegration.sample(config, integrand, measure; print = 10, Nblock = 16, reweight = 1e5)

    if isnothing(avg) == false
        # avg *= NF
        # std *= NF
        println("   theta            F                  A                  A-F")
        # println("Ver3: ")
        for (ki, theta) in enumerate(θgrid)
            println(@sprintf("%8.4f  %8.4f ±%8.4f  %8.4f ±%8.4f  %8.4f ±%8.4f", theta, avg[ki, 1], std[ki, 1], avg[ki, 2], std[ki, 2], avg[ki, 2] - avg[ki, 1], std[ki, 2] + std[ki, 1]))
        end
        # println("l=0 component: ")
        # println(@sprintf("averge   %8.4f ±%8.4f  %8.4f ±%8.4f  %8.4f ±%8.4f", l0(avg[:, 1]), l0(std[:, 1]), l0(avg[:, 2]), l0(std[:, 2]), l0(avg[:, 2] - avg[:, 1]), l0(std[:, 2] + std[:, 1])))
    end

end

MC()