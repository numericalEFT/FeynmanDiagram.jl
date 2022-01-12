using FeynmanDiagram
using AbstractTrees
# using NewickTree
using StaticArrays
const Weight = SVector{2,Float64}

# Parquet = Builder.Parquet
Parquet = ParquetNew

chan = [Parquet.T, Parquet.U, Parquet.S]

F = [Parquet.U, Parquet.S]
V = [Parquet.T, Parquet.U]

###################### ver4 to DiagTree ###########################################
para = GenericPara(
    loopDim = 3,
    innerLoopNum = 1,
    totalLoopNum = 4,
    totalTauNum = 3,
    spin = 2,
    hasTau = true,
    weightType = Float64,
    firstLoopIdx = 2,
    firstTauIdx = 1,
    filter = [Builder.NoFock,],
    # interaction = [Interaction(ChargeCharge, Instant),]
    interaction = [Interaction(UpUp, Instant), Interaction(UpDown, Instant),]
)

K0 = zeros(para.totalLoopNum)
KinL, KoutL, KinR, KoutR = deepcopy(K0), deepcopy(K0), deepcopy(K0), deepcopy(K0)
KinL[1] = KoutL[1] = 1
KinR[2] = KoutR[2] = 1
legK = [KinL, KoutL, KinR, KoutR]

varK = [rand(para.loopDim) for i in 1:para.totalLoopNum]
varT = [rand() for i in 1:para.totalTauNum]
evalK(basis) = sum([basis[i] * varK[i] for i in 1:para.totalLoopNum])
evalT(Tidx) = varT[Tidx]

diag, nodes = Parquet.buildVer4(para, legK, chan, F = F, V = V)
output = Parquet.classify!(diag, nodes, :name)
println(output)
uu, ud = output
# println(dir)
# dir, ex = Parquet.merge(diag, nodes, :DiEx)
# dir = [node.nidx for ]
# rootDir = DiagTree.addnode!(diag, DiagTree.ADD, :dir, dir.children, para = [0, 0, 0, 0])
# rootEx = DiagTree.addnode!(diag, DiagTree.ADD, :ex, ex.children, para = [0, 0, 0, 0])
# diag.root = [rootDir.index, rootEx.index]
# DiagTree.showTree(diag, rootDir.index)
# diag.root = [dir[2].index, ex[2].index]

DiagTree.showTree(diag, uu[2].index)
DiagTree.showTree(diag, ud[2].index)

##################### lower level subroutines  #######################################
# ver4 = Parquet.Ver4{Float64}(para, legK, chan, F, V)

########## use AbstractTrees interface to print/manipulate the tree
# print_tree(ver4)

# [println(node) for node in Leaves(ver4)]  #print all loopNum=0 ver4

# println("Iterate the tree use the AbstractTrees interface: ")
# [println(node) for node in PostOrderDFS(ver4)]  # iterator ver4 in depth-first search (children before parents)
# [println(node) for node in PreOrderDFS(ver4)]  # iterator ver4 (parents before children)

########## print tree to a newick format file  ##############
# io = open("./test.newick", "w")
# write(io, Parquet.newick(ver4))
# close(io)

########## use ete3 package to visualize tree
# Parquet.showTree(ver4, verbose = 1, depth = 3)  # visualize tree using python3 package ete3

######################################## self-energy  ################################################

para = Builder.GenericPara(
    loopDim = 3,
    innerLoopNum = 2,
    totalLoopNum = 4,
    totalTauNum = 3,
    spin = 2,
    # interactionTauNum = 1,
    firstLoopIdx = 2,
    firstTauIdx = 1,
    # weightType = Float64,
    hasTau = true,
    filter = [Builder.NoHatree,]
)

K0 = zeros(para.totalLoopNum)
K0[1] = 1.0
sigma, instant, dynamic = Parquet.buildSigma(para, K0)
# println(root)
rootidx = DiagTree.addnode!(sigma, DiagTree.ADD, :sum, vcat(instant, dynamic); para = [0, 0])
DiagTree.showTree(sigma, rootidx)


exit(0)

##################################### vertex 3   #################################################
para = Builder.GenericPara(
    loopDim = 3,
    innerLoopNum = 1,
    totalLoopNum = 4,
    totalTauNum = 3,
    spin = 2,
    interactionTauNum = 1,
    firstLoopIdx = 2,
    weightType = Float64,
    filter = [Builder.Proper,]
)

K0 = zeros(para.totalLoopNum)
Kin, Kout = deepcopy(K0), deepcopy(K0)
Kin[1] = 1
Kout[2] = 1
legK = [Kin, Kout]
Parquet.buildVer3(para, legK)