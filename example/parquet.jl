using ExpressionTree
using AbstractTrees
# using NewickTree
using StaticArrays
const Weight = SVector{2,Float64}

Parquet = Builder.Parquet

chan = [Parquet.T, Parquet.U, Parquet.S]

F = [Parquet.U, Parquet.S]
V = [Parquet.T, Parquet.U]

interactionTauNum = 1
loopNum = 3
spin = 2
Kdim = 3

# para = Parquet.Para(chan, F, V, loopNum, 2, Kdim, interactionTauNum, spin)
para = Builder.GenericPara(
    loopDim = 3,
    internalLoopNum = 1,
    totalLoopNum = 3,
    spin = 2,
    interactionTauNum = 1,
    weightType = Float64
)

K0 = zeros(2 + para.totalLoopNum)
KinL, KoutL, KinR, KoutR = deepcopy(K0), deepcopy(K0), deepcopy(K0), deepcopy(K0)
KinL[1] = KoutL[1] = 1
KinR[2] = KoutR[2] = 1
legK = [KinL, KoutL, KinR, KoutR]

varK = [rand(para.loopDim) for i in 1:2+para.totalLoopNum]
varT = [rand() for i in 1:2*(para.totalLoopNum+1)]
evalK(basis) = sum([basis[i] * varK[i] for i in 1:para.totalLoopNum])
evalT(Tidx) = varT[Tidx]


# diag, ver4, dir, ex = Parquet.build(Float64, para, legK)
# rootDir = DiagTree.addNode!(diag, DiagTree.ADD, :dir; child = dir, para = (0, 0, 0, 0))
# rootEx = DiagTree.addNode!(diag, DiagTree.ADD, :ex; child = ex, para = (0, 0, 0, 0))
# diag.root = [rootDir, rootEx]

# DiagTree.showTree(diag, rootDir)
# print_tree(ver4)

##################### lower level subroutines  #######################################
ver4 = Parquet.Ver4{Float64}(para, chan, F, V)

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
# Parquet.showTree(ver4, verbose = 1, depth = 3)  # visualize tree using python3 package ete3