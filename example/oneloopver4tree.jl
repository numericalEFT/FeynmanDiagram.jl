using ExpressionTree

DiagTree = GWKT.DiagTree
Parquet = GWKT.Parquet
Gsym = [:mirror]
Wsym = [:mirror, :timereversal]
Gtype, Wtype = 1, 2

"""
Build tree with KinL = KoutL = [1, 0, 0, 0], KinR = KoutR = [0, 1, 0]
"""
function build(spin)
    diag = DiagTree.Diagrams{Float64}()
    KinL = KoutL = [1, 0, 0, 0]
    KinR = KoutR = [0, 1, 0, 0]
    legK = [KinL, KoutL, KinR, KoutR]
    Td, Te = buildT(diag, legK, spin)
    DiagTree.showTree(diag, Td[1])
    DiagTree.showTree(diag, Te[1])
end


function buildT(diag, legK, spin)
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
    Sym = -1.0
    KinL, KinR, KoutL, KoutR = legK
    qd = KinL - KoutL
    K = [0, 0, 0, 1]
    Td, Te = [], []
    #construct the propagator table
    gK = [KinL + K - KoutL, K]
    gT = [[1, 3], [3, 1]]
    g = [DiagTree.addPropagator!(diag, Gtype, 0, gK[i], gT[i], Gsym)[1] for i = 1:2]
    # G order is 0

    vdK = [qd, qd]
    vdT = [[1, 1], [1, 1]]
    vd = [DiagTree.addPropagator!(diag, Wtype, 1, vdK[i], vdT[i], Wsym)[1] for i = 1:2]
    # W order is 1

    veK = [K - KoutL, KinR - K]
    veT = [[1, 1], [1, 1]]
    ve = [DiagTree.addPropagator!(diag, Wtype, 1, veK[i], veT[i], Wsym)[1] for i = 1:2]
    # W order is 1


    # contruct the tree
    MUL, ADD = 1, 2
    vdd = DiagTree.addNode!(diag, MUL, spin * Sym, [g[1], g[2], vd[1], vd[2]], [])
    vde = DiagTree.addNode!(diag, MUL, Sym, [g[1], g[2], vd[1], ve[2]], [])
    ved = DiagTree.addNode!(diag, MUL, Sym, [g[1], g[2], ve[1], vd[2]], [])
    push!(Td, DiagTree.addNode!(diag, ADD, 1.0, [], [vdd, vde, ved]; extT = [1, 1, 3, 3]))
    push!(Te, DiagTree.addNode!(diag, MUL, Sym, [g[1], g[2], ve[1], ve[2]], []; extT = [1, 1, 3, 3]))
    return Td, Te
end

build(2)