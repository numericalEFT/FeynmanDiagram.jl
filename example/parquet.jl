using ExpressionTree
using AbstractTrees
# using NewickTree
using StaticArrays
const Weight = SVector{2,Float64}

Parquet = Builder.Parquet

chan = [Parquet.T, Parquet.U, Parquet.S]

F = [Parquet.U, Parquet.S]
V = [Parquet.T, Parquet.U]

###################### ver4 to DiagTree ###########################################
para = Builder.GenericPara(
    loopDim = 3,
    innerLoopNum = 1,
    totalLoopNum = 3,
    totalTauNum = 2,
    spin = 2,
    interactionTauNum = 1,
    weightType = Float64,
    filter = [Builder.NoBubble,]
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

# diag, ver4, dir, ex = Parquet.buildVer4(para, legK, chan, F, V)
# rootDir = DiagTree.addNode!(diag, DiagTree.ADD, :dir; child = dir, para = [0, 0, 0, 0])
# rootEx = DiagTree.addNode!(diag, DiagTree.ADD, :ex; child = ex, para = [0, 0, 0, 0])
# diag.root = [rootDir, rootEx]

# DiagTree.showTree(diag, rootDir)

##################### lower level subroutines  #######################################
ver4 = Parquet.Ver4{Float64}(para, chan, legK, F, V)

########## use AbstractTrees interface to print/manipulate the tree
print_tree(ver4)

# [println(node) for node in Leaves(ver4)]  #print all loopNum=0 ver4

println("Iterate the tree use the AbstractTrees interface: ")
[println(node) for node in PostOrderDFS(ver4)]  # iterator ver4 in depth-first search (children before parents)
# [println(node) for node in PreOrderDFS(ver4)]  # iterator ver4 (parents before children)

########## print tree to a newick format file  ##############
# io = open("./test.newick", "w")
# write(io, Parquet.newick(ver4))
# close(io)

########## use ete3 package to visualize tree
Parquet.showTree(ver4, verbose = 1, depth = 3)  # visualize tree using python3 package ete3

exit(0)

######################################## self-energy  ################################################

para = Builder.GenericPara(
    loopDim = 3,
    innerLoopNum = 3,
    totalLoopNum = 4,
    totalTauNum = 3,
    spin = 2,
    interactionTauNum = 1,
    firstLoopIdx = 2,
    weightType = Float64,
    filter = [Builder.NoHatree,]
)

K0 = zeros(para.totalLoopNum)
K0[1] = 1.0
sigma, root = Parquet.buildSigma(para, K0)
# println(root)
rootidx = DiagTree.addNode!(sigma, DiagTree.ADD, :sum; child = root, para = [0, 0])
DiagTree.showTree(sigma, rootidx)