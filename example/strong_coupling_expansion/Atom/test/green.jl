@testset "2-point Green's function" begin
    using LinearAlgebra
    using Atom.Green:GreenN, Gn, Gnc
    
    # construct model
    # U, μ, β = 1.0, 0.4, 5.0
    U, μ, β = 10.0, 5.0, 100.0
    rtol = 1e-6

    E = [0.0, -μ, -μ, U - 2μ]
    H = zeros(Float64, (4, 4))
    H[diagind(H)] = E

    # |0>=|00>=1, |↑>=|10>=2, |↓>=|01>=3, |↑↓>=|11>=4
    # the first is the forck state for ↑ spin, the second is forck state for ↓

    cpup = zeros(Float64, (4, 4)) 
    cpdown = zeros(Float64, (4, 4)) 

    cpup[2, 1], cpup[4, 3] = 1, 1
    cpdown[3, 1], cpdown[4, 2] = 1, 1
    cmup, cmdown = cpup', cpdown'
    
    @assert abs(tr(cpup * cmdown)) < 1e-16
    @assert abs(tr(cpdown * cmup)) < 1e-16

    m = Green.Model(β, H, [cpup, cpdown], true)

    G(τ) = (exp(μ * τ) + exp(μ * β - (U - μ) * τ)) / (1 + 2 * exp(μ * β) + exp(-(U - 2μ) * β))

    @test Green.density(m, UP) ≈ G(β)
    # @test thermalavg(m.n[UP], m.E, m.β, m.Z) ≈ G(β)

    τi, τo = 0.1, 0.8
    g2 = GreenN(m, [τi, τo], [UP, DOWN])
    @test abs(Gn(m, g2)) < 1e-16
    g2 = GreenN(m, [τi, τo], [UP, UP])
    @test Gn(m, g2) ≈ G(τo - τi)

    # G2(m, g2)
    Gn(m, g2)
    # @time G2(m, g2)
    @time Gn(m, g2)

    τi, τo = 0.8, 0.1
    g2 = GreenN(m, [τi, τo], [UP, DOWN])
    @test abs(Gn(m, g2)) < 1e-16
    g2 = GreenN(m, [τi, τo], [UP, UP])
    @test Gn(m, g2) ≈ -G(τo - τi + β)

    τ = [0.1, 0.8, 0.8 - 1e-10, 0.1 - 1e-10]
    orbital = [UP, UP, UP, UP]
    g4 = GreenN(m, τ, [UP, UP, UP, UP])
    @test Gn(m, g4) ≈ G(β) # density up or down is conserved, and n^2=m
    gc = 0.5 - 0.5^2 - G(τ[2] - τ[1])^2 # <n^2>-<n>^2=<n>-<n>^2
    @test Gnc(m, g4) ≈ gc
    # exit(0)

    Gn(m, g4)
    Gnc(m, g4)
    @time Gn(m, g4)
    @time Gnc(m, g4)

    δ = 1e-6
    # τ = [t for t in LinRange(3δ, β, 16)]
    for (ti, t) in enumerate(LinRange(3δ, β, 16))
    t4 = [0.0, 2δ, δ, t]
    g4 = GreenN(m, t4, [UP, UP, UP, UP])
        # println(Gnc(m, g4, [1,2,3,4]), " vs ", G(t))

        # <Tτ a(t)a(δ)a⁺(2δ)a⁺(0)> ≈ -<Tτ a(t)a(0)>
    @test abs(Gn(m, g4, [1,2,3,4]) + G(t)) < 10δ

        # <Tτ a(t)a(δ)a⁺(2δ)a⁺(0)>_{connected} ≈ 0.0
    @test abs(Gnc(m, g4, [1,2,3,4])) < 10δ
end

    δ = 1e-6
    for (ti, t) in enumerate(LinRange(5δ, β, 16))
    t6 = [0.0, 2δ, 4δ, δ, 3δ, t]
    g6 = GreenN(m, t6, [UP, UP, UP, UP, UP, UP])
    @test abs(Gn(m, g6, [1, 2,3,4,5,6]) + G(t)) < 10δ
    @test abs(Gnc(m, g6, [1, 2,3,4,5,6])) < 10δ
end

    t6 = [0.0, 2δ, 4δ, δ, 3δ, 0.5]
    g6 = GreenN(m, t6, [UP, UP, UP, UP, UP, UP])
    Gn(m, g6, [1,2,3,4,5,6])
    Gnc(m, g6, [1,2,3,4,5,6])
    @time Gn(m, g6, [1,2,3,4,5,6])
    @time Gnc(m, g6, [1,2,3,4,5,6])

    t12 = [0.0, 2δ, 4δ, δ, 3δ, 0.1, 0.15, 0.13, 0.11, 0.12, 0.05, 0.01]
    g12 = GreenN(m, t12, [UP for i in 1:12])
    Gn(m, g12)
    Gnc(m, g12)
    @time Gn(m, g12)
    @time Gnc(m, g12)
    println(Gnc(m, g12))

    # t14 = [0.0, 2δ, 4δ, δ, 3δ, 0.1, 0.15, 0.13, 0.11, 0.12, 0.05, 0.01, 0.06, 0.03]
    # g14 = Green(m, t14, [UP for i in 1:14])
    # Gn(m, g14)
    # Gnc(m, g14)
    # @time Gn(m, g14)
    # @time Gnc(m, g14)
    # println(Gnc(m, g12))

end