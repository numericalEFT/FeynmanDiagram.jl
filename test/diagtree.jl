@testset "DiagTree" begin
    DiagTree = GWKT.DiagTree
    # Write your tests here.
    Gsymmetry = [:mirror, :particlehole]
    Wsymmetry = [:mirror, :timereversal]
    kbasis = [1, 0, 1, 1]
    Tidx = [(1, 2), (1, 2), (2, 1), (2, 1)]
    Kidx = [kbasis, -kbasis, kbasis, -kbasis]
    N = length(Tidx)
    Gtype, Wtype = 1, 2

    ############# Test G with symmetry ############################
    diag = DiagTree.Diagrams{Float64}()
    idx = [DiagTree.addPropagator!(diag, Gtype, 1, Kidx[i], Tidx[i], Gsymmetry)[1] for i = 1:N]
    @test idx == [1, 1, 1, 1]
    ############# Test G without symmetry ############################
    diag = DiagTree.Diagrams{Float64}()
    idx = [DiagTree.addPropagator!(diag, Gtype, 1, Kidx[i], Tidx[i], [])[1] for i = 1:N]
    @test idx == [1, 2, 3, 4]
    ############# Test W with symmetry ############################
    diag = DiagTree.Diagrams{Float64}()
    idx = [DiagTree.addPropagator!(diag, Wtype, 1, Kidx[i], Tidx[i], Wsymmetry)[1] for i = 1:N]
    @test idx == [1, 1, 1, 1]
    ############# Test W without symmetry ############################
    diag = DiagTree.Diagrams{Float64}()
    idx = [DiagTree.addPropagator!(diag, Wtype, 1, Kidx[i], Tidx[i], [])[1] for i = 1:N]
    @test idx == [1, 2, 3, 4]
end

@testset "Diagrams" begin

    """
        k1-k3                     k2+k3 
        |                         | 
    t1.L ↑     t1.L       t2.L     ↑ t2.L
        |-------------->----------|
        |       |  k3+k4   |      |
        |   v   |          |  v   |
        |       |    k4    |      |
        |--------------<----------|
    t1.L ↑    t1.L        t2.L     ↑ t2.L
        |                         | 
        k1                        k2
    """
    # We only consider the direct part of the above diagram
    DiagTree = GWKT.DiagTree
    Gsym = [:mirror]
    Wsym = [:mirror, :timereversal]
    Gtype, Wtype = 1, 2
    spin = 2.0
    kF, β, mass2 = 1.919, 0.5, 1.0
    K1 = K2 = [kF, 0.0, 0.0]
    K3 = [0.5 * kF, 0.5 * kF, 0.0]
    K4 = [0.0, kF, 0.0]
    varK = [K1, K2, K3, K4]
    varT = [0.2 * β, 0.3 * β]

    diag = DiagTree.Diagrams{Float64}()
    gK = [[0, 0, 1, 1], [0, 0, 0, 1]]
    gT = [[1, 2], [2, 1]]
    g = [DiagTree.addPropagator!(diag, Gtype, 0, gK[i], gT[i], Gsym)[1] for i = 1:2]

    vdK = [[0, 0, 1, 0], [0, 0, 1, 0]]
    vdT = [[1, 1], [2, 2]]
    vd = [DiagTree.addPropagator!(diag, Wtype, 0, vdK[i], vdT[i], Wsym)[1] for i = 1:2]

    veK = [[1, 0, -1, -1], [0, 1, 0, -1]]
    veT = [[1, 1], [2, 2]]
    ve = [DiagTree.addPropagator!(diag, Wtype, 0, veK[i], veT[i], Wsym)[1] for i = 1:2]

    MUL, ADD = 1, 2

    gg_n = DiagTree.addNode!(diag, MUL, 1.0, [g[1], g[2]], [])
    vdd = DiagTree.addNode!(diag, MUL, spin, [vd[1], vd[2]], [])
    vde = DiagTree.addNode!(diag, MUL, -1.0, [vd[1], ve[2]], [])
    ved = DiagTree.addNode!(diag, MUL, -1.0, [ve[1], vd[2]], [])
    vsum = DiagTree.addNode!(diag, ADD, 1.0, [], [vdd, vde, ved])
    root = DiagTree.addNode!(diag, MUL, 1.0, [], [gg_n, vsum])

    #make sure the total number of diagrams are correct
    evalPropagator1(type, K, Tidx, varT) = 1.0
    @test DiagTree.evalNaive(diag, evalPropagator1, varK, varT) ≈ -2 + 1 * spin

    function evalPropagator2(type, K, Tidx, varT)
        if type == Gtype
            ϵ = dot(K, K) / 2 - kF^2
            τ = varT[Tidx[2]] - varT[Tidx[1]]
            return Spectral.kernelFermiT(τ, ϵ, β)
        elseif type == Wtype
            return 8π / (dot(K, K) + mass2)
        else
            error("not implemented")
        end
    end

    getK(basis, varK) = sum([basis[i] * K for (i, K) in enumerate(varK)])

    gw = [evalPropagator2(Gtype, getK(gK[i], varK), gT[i], varT) for i = 1:2]
    vdw = [evalPropagator2(Wtype, getK(vdK[i], varK), vdT[i], varT) for i = 1:2]
    vew = [evalPropagator2(Wtype, getK(veK[i], varK), veT[i], varT) for i = 1:2]

    Vweight = spin * vdw[1] * vdw[2] - vdw[1] * vew[2] - vew[1] * vdw[2]
    Weight = gw[1] * gw[2] * Vweight

    # println(DiagTree.evalNaive(diag, evalPropagator2, varK, varT))
    # println(Weight)
    @test DiagTree.evalNaive(diag, evalPropagator2, varK, varT) ≈ Weight

    DiagTree.showTree(diag)

end