@testset "LoopPool" begin
    dim, N = 3, 4
    loopPool = DiagTree.LoopPool(:K, dim, N, Float64)
    idx1 = DiagTree.append(loopPool, [1.0, 0.0, 0.0, 0.0])
    idx2 = DiagTree.append(loopPool, [1.0, 1.0, 0.0, 0.0])
    idx3 = DiagTree.append(loopPool, [1.0, 1.0, 0.0, 0.0])
    idx4 = DiagTree.append(loopPool, [1.0, 0.0, 0.0, 0.0])
    @test length(loopPool) == 2
    @test idx1 == idx4
    @test idx2 == idx3
end

@testset "CachedPool" begin
    objType, weightType = Int, Float64
    pool = DiagTree.CachedPool(:P, objType, weightType)
    idx1 = DiagTree.append(pool, 1)
    idx2 = DiagTree.append(pool, 2)
    idx3 = DiagTree.append(pool, 2)
    idx4 = DiagTree.append(pool, 1)
    @test length(pool) == 2
    @test idx1 == idx4
    @test idx2 == idx3

    # test symmetry operator
    # symmetry = Var.refection(Float64, 3)
    # a = Var.VectorVariable([1.0, 2.0, 2.0])
    # b = Var.VectorVariable([-1.0, -2.0, -2.0])
    # c = Var.VectorVariable([1.0, 2.0, 2.0])
    # @test isequal(a, c)
    # @test isequal(a, b) == false #two vectors are different if the symmetries are different, regardless of the basis 

    #test Node
    # node1 = DiagTree.Node(1; components = [[1, 2], [3, 4]], child = [1, 2])
    # node2 = DiagTree.Node(1; components = [[1, 2], [3, 4]], child = [1, 2])
    # @test (node1 != node2) == false

    #test diagram
end

@testset "Propagator" begin
    order1, order2 = 1, 2
    para1, para2 = 5, 6
    factor1, factor2 = 1.1, 1.2
    loopidx1, loopidx2 = 1, 2
    sitebasis1, sitebasis2 = [1, 2], [2, 3]
    function propagator(; order = order1, para = para1, factor = factor1, loopidx = loopidx1, sitebasis = sitebasis1)
        return DiagTree.Propagator{Float64}(:none, order, para, factor, loopidx, sitebasis)
    end
    @test propagator() != propagator(order = order2)
    @test propagator() != propagator(para = para2)
    @test propagator() != propagator(factor = factor2)
    @test propagator() != propagator(loopidx = loopidx2)
    @test propagator() != propagator(sitebasis = sitebasis2)
end

@testset "Node" begin
    operation1, operation2 = DiagTree.ADD, DiagTree.MUL
    para1, para2 = 5, 6
    factor1, factor2 = 1.1, 1.2
    components1, components2 = [[1,], [2,]], [[1,],]
    child1, child2 = [3, 4], [3, 5]
    function node(; operation = operation1, para = para1, factor = factor1, components = components1, child = child1)
        return DiagTree.Node{Float64}(:none, operation, para, components, child, factor)
    end
    @test node() != node(operation = operation2)
    @test node() != node(para = para2)
    @test node() != node(factor = factor2)
    @test node() != node(components = components2)
    @test node() != node(child = child2)
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
    Gtype, Wtype = 1, 2
    spin = 2.0
    D = 3
    kF, β, mass2 = 1.919, 0.5, 1.0

    varK = [rand(D) for i = 1:4] #k1, k2, k3, k4
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
    g = [DiagTree.addPropagator(diag, :Gpool, gorder, :G; site = gT[i], loop = gK[i]) for i = 1:2]

    vdK = [[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 1.0, 0.0]]
    vdT = [[1, 1], [2, 2]]
    vd = [DiagTree.addPropagator(diag, :Vpool, vorder, :Vd; loop = vdK[i]) for i = 1:2]

    veK = [[1, 0, -1, -1], [0, 1, 0, -1]]
    veT = [[1, 1], [2, 2]]
    ve = [DiagTree.addPropagator(diag, :Vpool, vorder, :Ve; loop = veK[i]) for i = 1:2]
    # ve = [DiagTree.addPropagator!(diag, Wtype, 1, veK[i], veT[i], Wsym)[1] for i = 1:2]
    # # W order is 1

    # # contruct the tree
    MUL, ADD = DiagTree.MUL, DiagTree.ADD
    ggn = DiagTree.addNodeByName(diag, MUL, :gxg; Gpool = [g[1], g[2]], factor = 1.0)
    vdd = DiagTree.addNodeByName(diag, MUL, :dxd; Vpool = [vd[1], vd[2]], factor = spin)
    vde = DiagTree.addNodeByName(diag, MUL, :dxe; Vpool = [vd[1], ve[2]], factor = -1.0)
    ved = DiagTree.addNodeByName(diag, MUL, :exd; Vpool = [ve[1], vd[2]], factor = -1.0)
    vsum = DiagTree.addNodeByName(diag, ADD, :sum; child = [vdd, vde, ved], factor = 1.0)
    root = DiagTree.addNode(diag, MUL, :root; child = [ggn, vsum], factor = 1.0)
    push!(diag.root, root)

    printBasisPool(diag)
    printPropagator(diag)
    printNodes(diag)
    DiagTree.showTree(diag, diag.root[1])
    # DiagTree.printBasisPool(diag, "test.txt")

    # #make sure the total number of diagrams are correct
    # evalPropagator1(type, K, Tidx, varT, factor = 1.0, para = nothing) = 1.0
    # @test DiagTree.evalNaive(diag, evalPropagator1, varK, varT) ≈ -2 + 1 * spin

    # #more sophisticated test of the weight evaluation
    # function evalPropagator2(type, K, Tidx, varT, factor = 1.0, para = nothing)
    #     if type == Gtype
    #         ϵ = dot(K, K) / 2 - kF^2
    #         τ = varT[Tidx[2]] - varT[Tidx[1]]
    #         return Spectral.kernelFermiT(τ, ϵ, β)
    #     elseif type == Wtype
    #         return 8π / (dot(K, K) + mass2)
    #     else
    #         error("not implemented")
    #     end
    # end

    # getK(basis, varK) = sum([basis[i] * K for (i, K) in enumerate(varK)])

    # gw = [evalPropagator2(Gtype, getK(gK[i], varK), gT[i], varT) for i = 1:2]
    # vdw = [evalPropagator2(Wtype, getK(vdK[i], varK), vdT[i], varT) for i = 1:2]
    # vew = [evalPropagator2(Wtype, getK(veK[i], varK), veT[i], varT) for i = 1:2]

    # Vweight = spin * vdw[1] * vdw[2] - vdw[1] * vew[2] - vew[1] * vdw[2]
    # Weight = gw[1] * gw[2] * Vweight

    # # println(DiagTree.evalNaive(diag, evalPropagator2, varK, varT))
    # # println(Weight)
    # @test DiagTree.evalNaive(diag, evalPropagator2, varK, varT) ≈ Weight


end

# @testset "DiagTree" begin
#     DiagTree = GWKT.DiagTree
#     # Write your tests here.
#     Gsymmetry = [:mirror, :particlehole]
#     Wsymmetry = [:mirror, :timereversal]
#     kbasis = [1, 0, 1, 1]
#     Tidx = [(1, 2), (1, 2), (2, 1), (2, 1)]
#     Kidx = [kbasis, -kbasis, kbasis, -kbasis]
#     N = length(Tidx)
#     Gtype, Wtype = 1, 2

#     ############# Test G with symmetry ############################
#     diag = DiagTree.Diagrams{Float64}()
#     idx = [DiagTree.addPropagator!(diag, Gtype, 0, Kidx[i], Tidx[i], Gsymmetry)[1] for i = 1:N]
#     #G has an order 0
#     @test idx == [1, 1, 1, 1]
#     ############# Test G without symmetry ############################
#     diag = DiagTree.Diagrams{Float64}()
#     idx = [DiagTree.addPropagator!(diag, Gtype, 0, Kidx[i], Tidx[i], [])[1] for i = 1:N]
#     #G has an order 0
#     @test idx == [1, 2, 3, 4]
#     ############# Test W with symmetry ############################
#     diag = DiagTree.Diagrams{Float64}()
#     idx = [DiagTree.addPropagator!(diag, Wtype, 1, Kidx[i], Tidx[i], Wsymmetry)[1] for i = 1:N]
#     #W has an order 1
#     @test idx == [1, 1, 1, 1]
#     ############# Test W without symmetry ############################
#     diag = DiagTree.Diagrams{Float64}()
#     idx = [DiagTree.addPropagator!(diag, Wtype, 1, Kidx[i], Tidx[i], [])[1] for i = 1:N]
#     #W has an order 1
#     @test idx == [1, 2, 3, 4]
# end

# @testset "Diagrams" begin

#     """
#         k1-k3                     k2+k3 
#         |                         | 
#     t1.L ↑     t1.L       t2.L     ↑ t2.L
#         |-------------->----------|
#         |       |  k3+k4   |      |
#         |   v   |          |  v   |
#         |       |    k4    |      |
#         |--------------<----------|
#     t1.L ↑    t1.L        t2.L     ↑ t2.L
#         |                         | 
#         k1                        k2
#     """
#     # We only consider the direct part of the above diagram
#     DiagTree = GWKT.DiagTree
#     Gsym = [:mirror]
#     Wsym = [:mirror, :timereversal]
#     Gtype, Wtype = 1, 2
#     spin = 2.0
#     kF, β, mass2 = 1.919, 0.5, 1.0
#     varK = [rand(3) for i = 1:4] #k1, k2, k3, k4
#     varT = [rand() * β, rand() * β]

#     diag = DiagTree.Diagrams{Float64}()

#     #construct the propagator table
#     gK = [[0, 0, 1, 1], [0, 0, 0, 1]]
#     gT = [[1, 2], [2, 1]]
#     g = [DiagTree.addPropagator!(diag, Gtype, 0, gK[i], gT[i], Gsym)[1] for i = 1:2]
#     # G order is 0

#     vdK = [[0, 0, 1, 0], [0, 0, 1, 0]]
#     vdT = [[1, 1], [2, 2]]
#     vd = [DiagTree.addPropagator!(diag, Wtype, 1, vdK[i], vdT[i], Wsym)[1] for i = 1:2]
#     # W order is 1

#     veK = [[1, 0, -1, -1], [0, 1, 0, -1]]
#     veT = [[1, 1], [2, 2]]
#     ve = [DiagTree.addPropagator!(diag, Wtype, 1, veK[i], veT[i], Wsym)[1] for i = 1:2]
#     # W order is 1


#     # contruct the tree
#     MUL, ADD = 1, 2
#     gg_n = DiagTree.addNode!(diag, MUL, 1.0, [g[1], g[2]], [])
#     vdd = DiagTree.addNode!(diag, MUL, spin, [vd[1], vd[2]], [])
#     vde = DiagTree.addNode!(diag, MUL, -1.0, [vd[1], ve[2]], [])
#     ved = DiagTree.addNode!(diag, MUL, -1.0, [ve[1], vd[2]], [])
#     vsum = DiagTree.addNode!(diag, ADD, 1.0, [], [vdd, vde, ved])
#     root = DiagTree.addNode!(diag, MUL, 1.0, [], [gg_n, vsum])

#     # DiagTree.showTree(diag)

#     #make sure the total number of diagrams are correct
#     evalPropagator1(type, K, Tidx, varT, factor = 1.0, para = nothing) = 1.0
#     @test DiagTree.evalNaive(diag, evalPropagator1, varK, varT) ≈ -2 + 1 * spin

#     #more sophisticated test of the weight evaluation
#     function evalPropagator2(type, K, Tidx, varT, factor = 1.0, para = nothing)
#         if type == Gtype
#             ϵ = dot(K, K) / 2 - kF^2
#             τ = varT[Tidx[2]] - varT[Tidx[1]]
#             return Spectral.kernelFermiT(τ, ϵ, β)
#         elseif type == Wtype
#             return 8π / (dot(K, K) + mass2)
#         else
#             error("not implemented")
#         end
#     end

#     getK(basis, varK) = sum([basis[i] * K for (i, K) in enumerate(varK)])

#     gw = [evalPropagator2(Gtype, getK(gK[i], varK), gT[i], varT) for i = 1:2]
#     vdw = [evalPropagator2(Wtype, getK(vdK[i], varK), vdT[i], varT) for i = 1:2]
#     vew = [evalPropagator2(Wtype, getK(veK[i], varK), veT[i], varT) for i = 1:2]

#     Vweight = spin * vdw[1] * vdw[2] - vdw[1] * vew[2] - vew[1] * vdw[2]
#     Weight = gw[1] * gw[2] * Vweight

#     # println(DiagTree.evalNaive(diag, evalPropagator2, varK, varT))
#     # println(Weight)
#     @test DiagTree.evalNaive(diag, evalPropagator2, varK, varT) ≈ Weight


# end