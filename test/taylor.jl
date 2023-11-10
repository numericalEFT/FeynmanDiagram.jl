using FeynmanDiagram: Taylor as Taylor

@testset verbose = true "TaylorSeries" begin
    using FeynmanDiagram.Taylor:
        getcoeff, set_variables, taylor_factorial
    a, b, c, d, e = set_variables("a b c d e", orders=[3, 3, 3, 3, 3])
    F1 = (a + b) * (a + b) * (a + b)
    print("$(F1)")
    @test getcoeff(F1, [2, 1, 0, 0, 0]) == 3.0
    @test getcoeff(F1, [1, 2, 0, 0, 0]) == 3.0
    @test getcoeff(F1, [3, 0, 0, 0, 0]) == 1.0
    @test getcoeff(F1, [0, 3, 0, 0, 0]) == 1.0
    F2 = (1 + a) * (3 + 2c)
    @test getcoeff(F2, [0, 0, 0, 0, 0]) == 3.0
    @test getcoeff(F2, [1, 0, 0, 0, 0]) == 3.0
    @test getcoeff(F2, [0, 0, 1, 0, 0]) == 2.0
    @test getcoeff(F2, [1, 0, 1, 0, 0]) == 2.0
    F3 = (a + b)^3
    @test getcoeff(F1, [2, 1, 0, 0, 0]) == 3.0
    @test getcoeff(F1, [1, 2, 0, 0, 0]) == 3.0
    @test getcoeff(F1, [3, 0, 0, 0, 0]) == 1.0
    @test getcoeff(F1, [0, 3, 0, 0, 0]) == 1.0
    using FeynmanDiagram.ComputationalGraphs:
        eval!, forwardAD, node_derivative, backAD, build_all_leaf_derivative, count_operation
    using FeynmanDiagram.Utility:
        taylorexpansion!, build_derivative_backAD!
    g1 = Graph([])
    g2 = Graph([])
    g3 = Graph([], factor=2.0)
    G3 = g1
    G4 = 1.0 * g1 * g1
    G5 = 1.0 * (3.0 * G3 + 0.5 * G4)
    G6 = (1.0 * g1 + 2.0 * g2) * (g1 + g3)

    set_variables("x y z", orders=[2, 3, 2])
    for G in [G3, G4, G5, G6]
        T, taylormap = taylorexpansion!(G)
        T_compare = build_derivative_backAD!(G)
        for (order, coeff) in T_compare.coeffs
            @test eval!(coeff) == eval!(taylor_factorial(order) * T.coeffs[order])
        end
    end

end


function getdiagram(spin=2.0, D=3, Nk=4, Nt=2)
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

    DiagTree.uidreset()
    # We only consider the direct part of the above diagram

    paraG = DiagParaF64(type=GreenDiag,
        innerLoopNum=0, totalLoopNum=Nk, loopDim=D,
        hasTau=true, totalTauNum=Nt)
    paraV = paraG

    # #construct the propagator table
    gK = [[0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 0.0, 1.0]]
    gT = [(1, 2), (2, 1)]
    g = [Diagram{Float64}(BareGreenId(paraG, k=gK[i], t=gT[i]), name=:G) for i in 1:2]

    vdK = [[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 1.0, 0.0]]
    # vdT = [[1, 1], [2, 2]]
    vd = [Diagram{Float64}(BareInteractionId(paraV, ChargeCharge, k=vdK[i], permu=Di), name=:Vd) for i in 1:2]

    veK = [[1, 0, -1, -1], [0, 1, 0, -1]]
    # veT = [[1, 1], [2, 2]]
    ve = [Diagram{Float64}(BareInteractionId(paraV, ChargeCharge, k=veK[i], permu=Ex), name=:Ve) for i in 1:2]

    Id = GenericId(paraV)
    # contruct the tree
    ggn = Diagram{Float64}(Id, Prod(), [g[1], g[2]])
    vdd = Diagram{Float64}(Id, Prod(), [vd[1], vd[2]], factor=spin)
    vde = Diagram{Float64}(Id, Prod(), [vd[1], ve[2]], factor=-1.0)
    ved = Diagram{Float64}(Id, Prod(), [ve[1], vd[2]], factor=-1.0)
    vsum = Diagram{Float64}(Id, Sum(), [vdd, vde, ved])
    root = Diagram{Float64}(Id, Prod(), [vsum, ggn], factor=1 / (2π)^D, name=:root)

    return root, gK, gT, vdK, veK
end

@testset "Taylor AD of DiagTree" begin

    DiagTree.uidreset()
    # We only consider the direct part of the above diagram
    spin = 0.5
    D = 3
    kF, β, mass2 = 1.919, 0.5, 1.0
    Nk, Nt = 4, 2

    root, gK, gT, vdK, veK = getdiagram(spin, D, Nk, Nt)

    #optimize the diagram
    DiagTree.optimize!([root,])

    # autodiff
    droot_dg = DiagTree.derivative([root,], BareGreenId)[1]
    droot_dv = DiagTree.derivative([root,], BareInteractionId)[1]
    # plot_tree(droot_dg)
    factor = 1 / (2π)^D
    DiagTree.eval!(root; eval=(x -> 1.0))
    @test root.weight ≈ (-2 + spin) * factor

    DiagTree.eval!(droot_dg; eval=(x -> 1.0))
    @test droot_dg.weight ≈ (-2 + spin) * 2 * factor

    DiagTree.eval!(droot_dv; eval=(x -> 1.0))
    @test droot_dv.weight ≈ (-2 + spin) * 2 * factor

    set_variables("x"; orders=[2])
    g, map = FrontEnds.Graph!(root)
    var_dependence = FrontEnds.extract_var_dependence(map, BareGreenId)
    t, taylormap = taylorexpansion!(g, var_dependence)
    order = [0]
    @test eval!(taylormap[g.id].coeffs[order]) ≈ (-2 + spin) * factor

    order = [1]
    @test eval!(taylormap[g.id].coeffs[order]) ≈ (-2 + spin) * factor * 2 * taylor_factorial(order)

    var_dependence = FrontEnds.extract_var_dependence(map, BareInteractionId)

    t, taylormap = taylorexpansion!(g, var_dependence)
    order = [0]
    @test eval!(taylormap[g.id].coeffs[order]) ≈ (-2 + spin) * factor

    order = [1]
    @test eval!(taylormap[g.id].coeffs[order]) ≈ (-2 + spin) * factor * 2 * taylor_factorial(order)
end