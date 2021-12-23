@testset "Pool" begin
    # pool = DiagTree.Pool{Int64,Float64}()
    # DiagTree.append(pool, 1, 2.0)
    # DiagTree.append(pool, 2, 3.0)
    # DiagTree.append(pool, 3, 4.0)

    # v = DiagTree.SubPool(pool, [1, 3])
    # # v = view(pool, [1, 3])
    # new = 1.5
    # v.pool[v.idx[1]].curr = new
    # # println(typeof(v))


    # @test pool[1].curr ≈ new #test view only returns the reference


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

    # reflection = Var.refection(Float64, D)

    # println(typeof(varK))
    # Mom = Var.VectorVariable{Vector{Vector{Float64}},Float64}
    # Tpair = Var.ScalarVariable{Vector{Float64},Tuple{Int,Int}}
    # Base.isequal(a::Mom, b::Mom) = (Mom.basis ≈ Mom.basis) || (Mom.basis ≈ -Mom.basis)


    K0 = [0.0, 0.0, 0.0]
    T0 = 0.0

    calcK(para, basis) = sum([para[i] .* basis[i] for i = 1:length(para)])
    calcT(para, basis) = para[basis[2]] - para[basis[1]]

    gorder, vorder = 0, 1

    Cache = DiagTree.Cache

    # MomPool = DiagTree.Pool{Cache{MomBasis,Vector{Float64}}}()
    # TpairPool = DiagTree.Pool{TpairBasis}()
    MomBasis = Vector{Float64}
    TpairBasis = Tuple{Int,Int}

    MomPool = DiagTree.cachedPool(MomBasis, Vector{Float64})
    TpairPool = DiagTree.uncachedPool(TpairBasis)

    GPool = DiagTree.propagatorPool(Float64, Float64)
    VPool = DiagTree.propagatorPool(Float64, Float64)

    diag = DiagTree.Diagrams((MomPool, TpairPool), (GPool, VPool), Float64, Float64)

    # #construct the propagator table
    gK = [[0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 0.0, 1.0]]
    gT = [(1, 2), (2, 1)]
    g = [DiagTree.addPropagator(diag, 1, gorder, [[1, gK[i], K0], [2, gT[i], T0]]) for i = 1:2]

    vdK = [[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 1.0, 0.0]]
    vdT = [[1, 1], [2, 2]]
    vd = [DiagTree.addPropagator(diag, 2, vorder, [(1, vdK[i], K0),]) for i = 1:2]

    veK = [[1, 0, -1, -1], [0, 1, 0, -1]]
    veT = [[1, 1], [2, 2]]
    ve = [DiagTree.addPropagator(diag, 2, vorder, [(1, veK[i], K0),]) for i = 1:2]
    # ve = [DiagTree.addPropagator!(diag, Wtype, 1, veK[i], veT[i], Wsym)[1] for i = 1:2]
    # # W order is 1

    # # contruct the tree
    MUL, ADD = DiagTree.MUL, DiagTree.ADD
    ggn = DiagTree.addNode(diag, MUL, [[g[1], g[2]], []], [], factor = 1.0)
    vdd = DiagTree.addNode(diag, MUL, [[], [vd[1], vd[2]]], [], factor = spin)
    vde = DiagTree.addNode(diag, MUL, [[], [vd[1], ve[2]]], [], factor = -1.0)
    ved = DiagTree.addNode(diag, MUL, [[], [ve[1], vd[2]]], [], factor = -1.0)
    vsum = DiagTree.addNode(diag, ADD, [[], []], [vdd, vde, ved], factor = 1.0)
    root = DiagTree.addNode(diag, MUL, [[], []], [ggn, vsum], factor = 1.0)
    push!(diag.root, root)

    DiagTree.showTree(diag)

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