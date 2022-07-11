@testset "LoopPool" begin
    dim, N = 3, 4
    loopPool = ExprTree.LoopPool(:K, dim, N, Float64)
    basis1 = [1.0, 0.0, 0.0, 1.0]
    basis2 = [1.0, 1.0, 0.0, 0.0]
    basis3 = [1.0, 0.0, -1.0, 1.0]
    idx1 = ExprTree.append(loopPool, basis1)
    idx2 = ExprTree.append(loopPool, basis2)
    idx3 = ExprTree.append(loopPool, basis2)
    idx4 = ExprTree.append(loopPool, basis1)
    idx5 = ExprTree.append(loopPool, basis3)
    @test length(loopPool) == 3
    @test idx1 == idx4
    @test idx2 == idx3

    varK = rand(dim, N)
    ExprTree.update(loopPool, varK)
    @test ExprTree.current(loopPool, 1) ≈ varK * basis1
    @test ExprTree.current(loopPool, 2) ≈ varK * basis2
    @test ExprTree.current(loopPool, 3) ≈ varK * basis3

end

# @testset "CachedPool" begin
#     objType, weightType = Int, Float64
#     pool = ExprTree.CachedPool(:P, objType, weightType)
#     idx1 = ExprTree.append(pool, 1)
#     idx2 = ExprTree.append(pool, 2)
#     idx3 = ExprTree.append(pool, 2)
#     idx4 = ExprTree.append(pool, 1)
#     @test length(pool) == 2
#     @test idx1 == idx4
#     @test idx2 == idx3
# end

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
    MomPool = ExprTree.LoopPool(:K, D, 4)

    diag = ExprTree.ExpressionTree(loopBasis=MomPool, nodePara=Int, weight=weightType)

    # #construct the propagator table
    gK = [[0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 0.0, 1.0]]
    gT = [(1, 2), (2, 1)]
    g = [ExprTree.addpropagator!(diag, :G; site=gT[i], loop=gK[i], para=1) for i = 1:2]

    vdK = [[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 1.0, 0.0]]
    vdT = [[1, 1], [2, 2]]
    vd = [ExprTree.addpropagator!(diag, :Vd; loop=vdK[i], para=2) for i = 1:2]

    veK = [[1, 0, -1, -1], [0, 1, 0, -1]]
    veT = [[1, 1], [2, 2]]
    ve = [ExprTree.addpropagator!(diag, :Ve; loop=veK[i], para=2) for i = 1:2]
    # ve = [ExprTree.addPropagator!(diag, Wtype, 1, veK[i], veT[i], Wsym)[1] for i = 1:2]
    # # W order is 1

    # # contruct the tree
    MUL, ADD = ExprTree.MUL, ExprTree.ADD
    ggn = ExprTree.addnode!(diag, MUL, :gxg, [g[1], g[2]], 1.0, para=0)
    vdd = ExprTree.addnode!(diag, MUL, :dxd, [vd[1], vd[2]], spin, para=0)
    vde = ExprTree.addnode!(diag, MUL, :dxe, [vd[1], ve[2]], -1.0, para=0)
    ved = ExprTree.addnode!(diag, MUL, :exd, [ve[1], vd[2]], -1.0, para=0)
    vsum = ExprTree.addnode!(diag, ADD, :sum, [vdd, vde, ved], 1.0, para=0)
    root = ExprTree.addnode!(diag, MUL, :root, [ggn, vsum], 1.0, para=0)
    push!(diag.root, root)
    ExprTree.initialize!(diag.node)

    # printBasisPool(diag)
    # printPropagator(diag)
    # ExprTree.printNodes(diag)
    # ExprTree.showTree(diag, diag.root[1])

    # #make sure the total number of diagrams are correct
    let
        DiagTree.eval(para, K, Tbasis, varT) = 1.0
        ExprTree.evalKT!(diag, varK, varT)
        @test diag[1] ≈ -2 + 1 * spin
    end

    # #more sophisticated test of the weight evaluation
    let
        function evalG(K, τBasis, varT)
            ϵ = dot(K, K) / 2 - kF^2
            τ = varT[τBasis[2]] - varT[τBasis[1]]
            return Spectral.kernelFermiT(τ, ϵ, β)
        end

        evalV(K) = 8π / (dot(K, K) + mass2)

        function DiagTree.eval(para, K, Tbasis, varT)
            if para[1] == 1
                return evalG(K, Tbasis, varT)
            elseif para[1] == 2
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

        # println(ExprTree.printPropagator(diag))
        ExprTree.evalKT!(diag, varK, varT)
        @test diag[1] ≈ Weight
    end

end