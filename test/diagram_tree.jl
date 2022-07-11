@testset "Diagram" begin
    # Diagram = DiagTreeNew.Diagram
    # DiagramId = DiagTreeNew.DiagramId
    # add_subdiagram! = DiagTreeNew.add_subdiagram!

    struct ID <: DiagramId
        uid::Int
    end
    Base.show(io::IO, d::ID) = print(io, d.uid)
    # Base.isequal(a::ID, b::ID) = (a.index == b.index)
    # Base.Dict(d::ID) = Dict(:id => d.index)

    DiagTree.uidreset()

    W = Int
    ll = Diagram{W}(ID(3))
    l = Diagram{W}(ID(1), Sum(), [ll,])
    r = Diagram{W}(ID(2))
    root = Diagram{W}(ID(0), Sum(), [l, r])
    print_tree(root)
    """
    4 : 0=0=⨁ (2, 3)
    ├─ 2 : 1=0=⨁ (1)
    │  └─ 1 : 3=0
    └─ 3 : 2=0
    """

    collect(PostOrderDFS(root))
    @test [node.id.uid for node in PostOrderDFS(root)] == [3, 1, 2, 0]
    @test [node.id.uid for node in PreOrderDFS(root)] == [0, 1, 3, 2]
    @test [node.id.uid for node in Leaves(root)] == [3, 2]

    # eval(d::ID, vargs...) = d.uid
    @test DiagTree.eval!(root; eval=(d -> d.uid)) == sum(node.id.uid for node in Leaves(root))

    print_tree(root)
    # DiagTreeNew.plot_tree(root)

    println(toDataFrame([root,]))
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

    paraG = GenericPara(diagType=GreenDiag,
        innerLoopNum=0, totalLoopNum=Nk, loopDim=D,
        hasTau=true, totalTauNum=Nt)
    paraV = paraG

    # #construct the propagator table
    gK = [[0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 0.0, 1.0]]
    gT = [(1, 2), (2, 1)]
    g = [Diagram(BareGreenId(paraG, k=gK[i], t=gT[i]), name=:G) for i in 1:2]

    vdK = [[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 1.0, 0.0]]
    # vdT = [[1, 1], [2, 2]]
    vd = [Diagram(BareInteractionId(paraV, ChargeCharge, k=vdK[i], permu=Di), name=:Vd) for i in 1:2]

    veK = [[1, 0, -1, -1], [0, 1, 0, -1]]
    # veT = [[1, 1], [2, 2]]
    ve = [Diagram(BareInteractionId(paraV, ChargeCharge, k=veK[i], permu=Ex), name=:Ve) for i in 1:2]

    Id = GenericId(paraV)
    # contruct the tree
    ggn = Diagram(Id, Prod(), [g[1], g[2]])
    vdd = Diagram(Id, Prod(), [vd[1], vd[2]], factor=spin)
    vde = Diagram(Id, Prod(), [vd[1], ve[2]], factor=-1.0)
    ved = Diagram(Id, Prod(), [ve[1], vd[2]], factor=-1.0)
    vsum = Diagram(Id, Sum(), [vdd, vde, ved])
    root = Diagram(Id, Prod(), [vsum, ggn], factor=1 / (2π)^D, name=:root)

    return root, gK, gT, vdK, veK
end

@testset "Generic Diagrams" begin

    DiagTree.uidreset()
    # We only consider the direct part of the above diagram
    spin = 2.0
    D = 3
    kF, β, mass2 = 1.919, 0.5, 1.0
    Nk, Nt = 4, 2

    root, gK, gT, vdK, veK = getdiagram(spin, D, Nk, Nt)

    #optimize the diagram
    root = DiagTree.optimize(root)

    # autodiff
    droot_dg = DiagTree.derivative(root, BareGreenId)
    droot_dv = DiagTree.derivative(root, BareInteractionId)
    # plot_tree(droot_dg)

    DiagTree.eval!(root; eval=(x -> 1.0))
    @test root.weight ≈ -2 + spin

    DiagTree.eval!(droot_dg; eval=(x -> 1.0))
    @test root.weight ≈ (-2 + spin) * 2

    DiagTree.eval!(droot_dv; eval=(x -> 1.0))
    @test root.weight ≈ (-2 + spin) * 2

    # #more sophisticated test of the weight evaluation
    varK = rand(D, Nk)
    varT = [rand() * β for t in 1:Nt]

    function evalG(K, τBasis, varT, order=0)
        ϵ = dot(K, K) / 2 - kF^2
        τ = varT[τBasis[2]] - varT[τBasis[1]]
        if order == 0
            return Spectral.kernelFermiT(τ, ϵ, β)
        elseif order == 1
            return Spectral.kernelFermiT(τ, ϵ, β) * 3.1415
        else
            error("not implemented!")
        end
    end

    function evalV(K, order=0)
        if order == 0
            return 8π / (dot(K, K) + mass2)
        elseif order == 1
            return 8π / (dot(K, K) + mass2) * 3.1415
        else
            error("not implemented!")
        end
    end

    # # getK(basis, varK) = sum([basis[i] * K for (i, K) in enumerate(varK)])
    getK(basis, varK) = varK * basis

    eval(id::BareGreenId, varK, varT) = evalG(getK(id.extK, varK), id.extT, varT, id.order[1])
    eval(id::BareInteractionId, varK, varT) = evalV(getK(id.extK, varK), id.order[2])

    gw = [evalG(getK(gK[i], varK), gT[i], varT) for i = 1:2]
    vdw = [evalV(getK(vdK[i], varK)) for i = 1:2]
    vew = [evalV(getK(veK[i], varK)) for i = 1:2]

    dgw = [evalG(getK(gK[i], varK), gT[i], varT, 1) for i = 1:2]
    dvdw = [evalV(getK(vdK[i], varK), 1) for i = 1:2]
    dvew = [evalV(getK(veK[i], varK), 1) for i = 1:2]

    Vweight = spin * vdw[1] * vdw[2] - vdw[1] * vew[2] - vew[1] * vdw[2]
    Gweight = gw[1] * gw[2]
    Weight = Gweight * Vweight / (2π)^D

    dVweight = spin * (dvdw[1] * vdw[2] + vdw[1] * dvdw[2]) -
               (dvdw[1] * vew[2] + vdw[1] * dvew[2]) -
               (dvew[1] * vdw[2] + vew[1] * dvdw[2])

    dGweight = dgw[1] * gw[2] + gw[1] * dgw[2]
    dWeight_dg = dGweight * Vweight / (2π)^D
    dWeight_dv = Gweight * dVweight / (2π)^D

    # print_tree(root)
    DiagTree.eval!(root, varK, varT; eval=eval)
    @test root.weight ≈ Weight

    DiagTree.eval!(droot_dg, varK, varT; eval=eval)
    @test droot_dg.weight ≈ dWeight_dg

    DiagTree.eval!(droot_dv, varK, varT; eval=eval)
    @test droot_dv.weight ≈ dWeight_dv

    ############### test diagram optimization #################
    uniqueG, uniqueInt = DiagTree.removeDuplicatedLeaves!([root,], verbose=1)
    @test length(uniqueG) == 2
    @test length(uniqueInt) == 3
    DiagTree.eval!(root, varK, varT; eval=eval)
    @test root.weight ≈ Weight
end

@testset "optimize" begin
    DiagTree.uidreset()

    W = Int
    lll = Diagram{W}(ID(5))
    ll = Diagram{W}(ID(3), Prod(), [lll,])
    l = Diagram{W}(ID(1), Sum(), [ll,])
    r = Diagram{W}(ID(2))
    root = Diagram{W}(ID(0), Sum(), [l, r])
    # print_tree(root)
    """
    5 : 0=0=⨁ (3, 4)
    ├─ 3 : 1=0=⨁ (2)
    │  └─ 2 : 3=0=Ⓧ (1)
    │     └─ 1 : 5=0
    └─ 4 : 2=0
    """

    #remove the 2, which only has one child
    DiagTree.removeOneChildParent!([root,])
    """
    4 : 0=0=⨁ (1, 3)
    ├─ 1 : 3=0
    └─ 3 : 2=0
    """

    # print_tree(root)

    @test root.subdiagram[1].hash == 1
    # print_tree(root)

    # root, gK, gT, vdK, veK = getdiagram()
    # uniqueG, uniqueInt = DiagTree.removeDuplicatedLeaves!(root, verbose = 1)
    # @test length(uniqueG) == 2
    # @test length(uniqueInt) == 3
end