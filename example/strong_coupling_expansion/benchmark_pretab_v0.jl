include("pretab_propagator_C4v.jl")
# include("pretab_propagator_general.jl")
using BenchmarkTools
using Test
using Random
using Plots
gr()

function exact_hopping_value(t, β, dμ, lambda, deriv_order, Lx, Ly,
    τ::Float64, r)
    g2c = 0.0
    N = Lx * Ly

    for xi in 1:Lx
        for yi in 1:Ly
            kx = 2π * (xi - 1) / Lx
            ky = 2π * (yi - 1) / Ly

            ek = -2t * (cos(kx) + cos(ky)) + dμ

            ω = -1.0 / ek
            λeff = sign(ω) * lambda
            ωscaled = ω / λeff

            g2c_τ = 0.0
            for o in 0:deriv_order
                val = propagator_derivative(τ, ωscaled, β, o)
                g2c_τ += val * (ωscaled^o) * binomial(deriv_order, o) * (-1)^o
            end

            phase = kx * r[1] + ky * r[2]
            g2c += cos(phase) * g2c_τ / λeff
        end
    end
    return g2c / N
end

function exact_hopping_muderiv(t, β, dμ, lambda, deriv_order, Lx, Ly,
    τ::Float64, r)
    g2c = 0.0
    N = Lx * Ly
    fact_m = factorial(deriv_order)

    for xi in 1:Lx
        for yi in 1:Ly
            kx = 2π * (xi - 1) / Lx
            ky = 2π * (yi - 1) / Ly

            ek = -2t * (cos(kx) + cos(ky)) + dμ

            ω = -1.0 / ek
            λeff = sign(ω) * lambda
            ωscaled = ω / λeff

            g2c_τ = 0.0
            term = 1.0
            for o in 1:deriv_order
                val = propagator_derivative(τ, ωscaled, β, o)
                term *= ωscaled
                g2c_τ += val * term * binomial(deriv_order - 1, o - 1)
            end

            phase = kx * r[1] + ky * r[2]
            g2c += cos(phase) * g2c_τ * λeff^(deriv_order - 1) * (-ωscaled)^deriv_order * fact_m
        end
    end
    return g2c / N
end

function run_general_test()
    println("=== Setting up Parameters ===")
    # Lx, Ly = 32, 32
    Lx, Ly = 64, 64
    t = 1.0
    β = 10.0
    lambda = 0.01
    dμ = -0.2
    order = 1
    # Ntau = 2000
    Ntau = 20000
    # Rtable = 20
    Rtable = 4

    # Case A: Uniform Grid
    grid_uniform = collect(range(0.0, β; length=Ntau))
    grid_uniform[1] = 1e-10
    println("\n[Building Uniform Table...]")
    pretab_uni = build_pretab(t, β, grid_uniform, propagator_muderiv!;
        lambda=lambda, dμ=dμ, deriv_order=order, Lkx=Lx, Lky=Ly, Rtable=Rtable)

    # Case B: Non-Uniform Grid (e.g., Power Grid)
    # Dense near 0 and β, sparse in middle
    power = 3
    grid_power = [(i / (Ntau - 1))^power * β for i in 0:(Ntau-1)]
    # grid_power[1] = 1e-10
    println("[Building Non-Uniform Table...]")
    pretab_non = build_pretab(t, β, grid_power, propagator_muderiv!;
        lambda=lambda, dμ=dμ, deriv_order=order, Lkx=Lx, Lky=Ly, Rtable=Rtable)

    println("\n=== 1. 精度测试: 密集 Tau 扫描 ===")

    # 测试参数
    # r_test = SVector{2,Int}(1, 0)    # 测试相邻格点 hopping
    r_test = SVector{2,Int}(2, 0)
    # r_test = SVector{2,Int}(20, 0)    # 测试相邻格点 hopping
    m_test = 1         # 0阶导数

    # 生成测试点：比插值网格密得多的点，且故意避开网格点
    # test_points = range(0.001, β - 0.001, length=5000)
    test_points = range(0.001, β - 0.001, length=550)

    errors = Float64[]
    vals_exact = Float64[]
    vals_fast = Float64[]

    println("正在扫描 $(length(test_points)) 个测试点...")

    for τ in test_points
        # v_ex = exact_hopping_value(t, β, dμ, lambda, m_test, Lx, Ly, τ, r_test)
        v_ex = exact_hopping_muderiv(t, β, dμ, lambda, m_test, Lx, Ly, τ, r_test)

        # v_fast = bare_line_fast(pretab_uni, r_test, τ, 1, m_test)
        v_fast = bare_line_fast(pretab_non, r_test, τ, 1, m_test)

        push!(vals_exact, v_ex)
        push!(vals_fast, v_fast)
        push!(errors, abs(v_ex - v_fast))
    end

    rel_err = errors ./ vals_exact
    max_err = maximum(errors)
    avg_err = sum(errors) / length(errors)
    max_val = maximum(abs.(vals_exact))

    println("-"^50)
    @printf "测试点数量: %d\n" length(test_points)
    @printf "最大绝对值: %.6f\n" max_val
    @printf "最大绝对误差: %.2e\n" max_err
    @printf "最大相对误差: %.2e\n" maximum(rel_err)
    @printf "平均绝对误差: %.2e\n" avg_err
    @printf "相对误差(Avg): %.2e\n" (avg_err / max_val)
    println("-"^50)

    if max_err < 1e-4
        println("✅ 精度测试通过 (Error < 1e-4)")
    else
        println("⚠️ 注意：精度误差较大，建议增加 N_grid 大小")
    end

    num_pltpt = 550

    p1 = plot(test_points[1:num_pltpt], vals_exact[1:num_pltpt], label="Exact", lw=2, color=:black, alpha=0.6, title="Interpolation Check")
    plot!(p1, test_points[1:num_pltpt], vals_fast[1:num_pltpt], label="Interp", style=:dash, color=:red, lw=1.5)
    # scatter!(p1, grid_visual, [bare_line_fast(pretab_vis, r_test, g, 0) for g in grid_visual],
    #     label="Grid Points", color=:red, markersize=3)
    ylabel!(p1, "G(τ)")

    final_plot = plot(p1, layout=(1, 1), size=(800, 600), dpi=300)
    # savefig(final_plot, "benchmark_pretab_uni_m$m_test.pdf")
    # savefig(final_plot, "benchmark_pretab_nonuni_m$m_test.pdf")
    # savefig(final_plot, "benchmark_pretab_nonuni_r20_m$m_test.pdf")
    savefig(final_plot, "benchmark_pretab_muderiv_nonuni_r20_m$m_test.pdf")
    println("✅ 图片已保存为 'benchmark_pretab.pdf'")

    println("\n=== Performance Benchmark ===")
    r_bench = SVector{2,Int}(1, 0)
    τ_bench = 1.2345
    m_bench = 0

    println("-> Exact Calculation (Loop over k):")
    b_ex = @benchmark exact_hopping_value($t, $β, $dμ, $lambda, $m_bench, $Lx, $Ly, $τ_bench, $r_bench)
    display(b_ex)
    @btime exact_hopping_value($t, $β, $dμ, $lambda, $m_bench, $Lx, $Ly, $τ_bench, $r_bench)

    println("\n1. Uniform Grid (Should utilize O(1) path):")
    # This should be very fast (~5-10ns)
    @btime bare_line_fast($pretab_uni, $r_bench, $τ_bench, 1, 0)

    println("\n2. Non-Uniform Grid (Should utilize O(log N) path):")
    # This should be slower due to searchsortedlast (~40-80ns depending on Ntau)
    @btime bare_line_fast($pretab_non, $r_bench, $τ_bench, 1, 0)
end

run_general_test()