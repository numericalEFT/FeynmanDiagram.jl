@testset "Diagram" begin
    # Diagram = DiagTreeNew.Diagram
    # DiagramId = DiagTreeNew.DiagramId
    # add_subdiagram! = DiagTreeNew.add_subdiagram!

    struct ID <: DiagramId
        index::Int
    end
    Base.show(io::IO, d::ID) = print(io, d.index)
    # Base.isequal(a::ID, b::ID) = (a.index == b.index)
    # Base.Dict(d::ID) = Dict(:id => d.index)
    DiagTreeNew.eval(d::ID) = d.index

    root = Diagram(ID(0), Sum())
    l = Diagram(ID(1), Sum())
    r = Diagram(ID(2), Sum())
    addSubDiagram!(root, [l, r])
    addSubDiagram!(l, Diagram(ID(3), Sum()))

    collect(PostOrderDFS(root))
    @test [node.id.index for node in PostOrderDFS(root)] == [3, 1, 2, 0]
    @test [node.id.index for node in PreOrderDFS(root)] == [0, 1, 3, 2]
    @test [node.id.index for node in Leaves(root)] == [3, 2]
    @test evalDiagTree!(root) == sum(node.id.index for node in Leaves(root))

    print_tree(root)
    DiagTreeNew.plot_tree(root)

    println(toDataFrame([root,]))
end

@testset "Generic Diagrams" begin

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
    spin = 2.0
    D = 3
    kF, β, mass2 = 1.919, 0.5, 1.0

    # varK = [rand(D) for i = 1:4] #k1, k2, k3, k4

    varK = rand(D, 4)
    varT = [rand() * β, rand() * β]

    K0 = [0.0, 0.0, 0.0]
    T0 = 0.0

    calcK(para, basis) = sum([para[i] .* basis[i] for i = 1:length(para)])
    calcT(para, basis) = para[basis[2]] - para[basis[1]]

    gorder, vorder = 0, 1

    weightType = Float64

    # function LoopPool(name::Symbol, dim::Int, N::Int, type::DataType)
    MomPool = DiagTree.LoopPool(:K, D, 4)

    GPool = DiagTree.propagatorPool(:Gpool, weightType)
    VPool = DiagTree.propagatorPool(:Vpool, weightType)

    diag = DiagTree.Diagrams(MomPool, (GPool, VPool), weightType)

    # #construct the propagator table
    gK = [[0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 0.0, 1.0]]
    gT = [(1, 2), (2, 1)]
    g = [DiagTree.addPropagator!(diag, :Gpool, gorder, :G; site = gT[i], loop = gK[i]) for i = 1:2]

    vdK = [[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 1.0, 0.0]]
    vdT = [[1, 1], [2, 2]]
    vd = [DiagTree.addPropagator!(diag, :Vpool, vorder, :Vd; loop = vdK[i]) for i = 1:2]

    veK = [[1, 0, -1, -1], [0, 1, 0, -1]]
    veT = [[1, 1], [2, 2]]
    ve = [DiagTree.addPropagator!(diag, :Vpool, vorder, :Ve; loop = veK[i]) for i = 1:2]
    # ve = [DiagTree.addPropagator!(diag, Wtype, 1, veK[i], veT[i], Wsym)[1] for i = 1:2]
    # # W order is 1

    # # contruct the tree
    MUL, ADD = DiagTree.MUL, DiagTree.ADD
    ggn = DiagTree.addNodeByName!(diag, MUL, :gxg, 1.0; Gpool = [g[1], g[2]])
    vdd = DiagTree.addNodeByName!(diag, MUL, :dxd, spin; Vpool = [vd[1], vd[2]])
    vde = DiagTree.addNodeByName!(diag, MUL, :dxe, -1.0; Vpool = [vd[1], ve[2]])
    ved = DiagTree.addNodeByName!(diag, MUL, :exd, -1.0; Vpool = [ve[1], vd[2]])
    vsum = DiagTree.addNodeByName!(diag, ADD, :sum, 1.0; child = [vdd, vde, ved])
    root = DiagTree.addNode!(diag, MUL, :root, 1.0; child = [ggn, vsum])
    push!(diag.root, root)

    # printBasisPool(diag)
    # printPropagator(diag)
    # printNodes(diag)
    # DiagTree.showTree(diag, diag.root[1])

    # #make sure the total number of diagrams are correct

    evalPropagator1(idx, object, K, varT, diag) = 1.0
    @test DiagTree.evalNaive(diag, varK, varT, evalPropagator1)[1] ≈ -2 + 1 * spin

    # #more sophisticated test of the weight evaluation

    function evalG(K, τBasis, varT)
        ϵ = dot(K, K) / 2 - kF^2
        τ = varT[τBasis[2]] - varT[τBasis[1]]
        return Spectral.kernelFermiT(τ, ϵ, β)
    end

    evalV(K) = 8π / (dot(K, K) + mass2)

    function evalPropagator2(idx, object, K, varT, diag)
        if idx == 1
            return evalG(K, object.siteBasis, varT)
        elseif idx == 2
            return evalV(K)
        else
            error("not implemented")
        end
    end

    # getK(basis, varK) = sum([basis[i] * K for (i, K) in enumerate(varK)])
    getK(basis, varK) = varK * basis

    gw = [evalG(getK(gK[i], varK), gT[i], varT) for i = 1:2]
    vdw = [evalV(getK(vdK[i], varK)) for i = 1:2]
    vew = [evalV(getK(veK[i], varK)) for i = 1:2]

    Vweight = spin * vdw[1] * vdw[2] - vdw[1] * vew[2] - vew[1] * vdw[2]
    Weight = gw[1] * gw[2] * Vweight

    # println(DiagTree.printPropagator(diag))
    @test DiagTree.evalNaive(diag, varK, varT, evalPropagator2)[1] ≈ Weight

end