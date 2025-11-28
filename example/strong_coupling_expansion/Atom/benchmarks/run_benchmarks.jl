#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using BenchmarkTools
using Atom
using Atom.Green: GreenN, Gn, dGn_dU_estimator, dGn_dμ_estimator
using LinearAlgebra: diagind
using Printf: @sprintf

const N_LIST = [1, 2, 3, 4]
const τ_PROBE = 0.37

function reference_model(; U=6.0, μ=3.0, β=25.0)
    E = [0.0, -μ, -μ, U - 2μ]
    H = zeros(Float64, 4, 4)
    H[diagind(H)] = E

    cpup = zeros(Float64, 4, 4)
    cpdown = zeros(Float64, 4, 4)
    cpup[2, 1], cpup[4, 3] = 1, 1
    cpdown[3, 1], cpdown[4, 2] = 1, 1

    return Green.Model(β, H, [cpup, cpdown], true)
end

function leg_times(N::Int)
    if N == 1
        return [0.15, 0.85]
    end
    creation = collect(range(0.05, 0.45; length=N))
    annihilation = collect(range(0.55, 0.95; length=N))
    return vcat(creation, annihilation)
end

function build_green_objects(m::Green.Model)
    greens = Dict{Int,GreenN}()
    for N in N_LIST
        τ = leg_times(N)
        orbital = fill(UP, 2N)
        greens[N] = GreenN(m, τ, orbital)
    end
    return greens
end

function benchmark_suite()
    m = reference_model()
    greens = build_green_objects(m)
    suite = BenchmarkGroup()
    labels = Dict{Int,Dict{Symbol,String}}()
    for N in N_LIST
        g = greens[N]
        label_gn = @sprintf("Gn:%2d-leg", 2N)
        label_dU = @sprintf("dGn/dU:%2d-leg", 2N)
        label_dμ = @sprintf("dGn/dμ:%2d-leg", 2N)
        labels[N] = Dict(:Gn => label_gn, :dU => label_dU, :dμ => label_dμ)
        suite[label_gn] = @benchmarkable Gn($m, $g)
        suite[label_dU] = @benchmarkable dGn_dU_estimator($m, $g, $τ_PROBE)
        suite[label_dμ] = @benchmarkable dGn_dμ_estimator($m, $g, $τ_PROBE)
    end
    return suite, labels
end

function summarise(results, labels)
    println("\nBenchmark summary (minimum times):")
    for N in N_LIST
        labs = labels[N]
        gn_time = minimum(results[labs[:Gn]]).time / 1e3
        dU_time = minimum(results[labs[:dU]]).time / 1e3
        dμ_time = minimum(results[labs[:dμ]]).time / 1e3
        println(@sprintf("N=%d (legs=%2d): Gn=%7.3f μs, dGn/dU=%7.3f μs, dGn/dμ=%7.3f μs",
                         N, 2N, gn_time, dU_time, dμ_time))
    end
end

function main()
    suite, labels = benchmark_suite()
    samples = parse(Int, get(ENV, "ATOM_BENCH_SAMPLES", "1000"))
    println("Running Atom benchmarks with $samples samples …")
    results = run(suite; verbose=true, samples=samples)
    summarise(results, labels)
end

isinteractive() || main()
