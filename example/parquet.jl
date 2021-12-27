using ExpressionTree
using AbstractTrees
# using NewickTree
using StaticArrays
const Weight = SVector{2,Float64}

Parquet = Builder.Parquet

chan = [Parquet.T, Parquet.U, Parquet.S]
interactionTauNum = 1
loopNum = 1
spin = 2
Kdim = 3

K0 = zeros(2 + loopNum)
KinL, KoutL, KinR, KoutR = deepcopy(K0), deepcopy(K0), deepcopy(K0), deepcopy(K0)
KinL[1] = KoutL[1] = 1
KinR[2] = KoutR[2] = 1
legK = [KinL, KoutL, KinR, KoutR]

varK = [rand(Kdim) for i in 1:2+loopNum]
varT = [rand() for i in 1:2*(loopNum+1)]
evalK(basis) = sum([basis[i] * varK[i] for i in 1:3])
evalT(Tidx) = varT[Tidx]

diag, ver4, dir, ex = Parquet.build(Float64, chan, loopNum, legK, Kdim, 3, interactionTauNum, spin)
rootDir = DiagTree.addNode!(diag, DiagTree.ADD, :dir; child = dir, para = (0, 0, 0, 0))
rootEx = DiagTree.addNode!(diag, DiagTree.ADD, :ex; child = ex, para = (0, 0, 0, 0))
diag.root = [rootDir, rootEx]

DiagTree.showTree(diag, rootDir)
# print_tree(ver4)

##################### lower level subroutines  #######################################
para = Parquet.Para(Float64, Kdim, loopNum, loopNum + 2, chan, interactionTauNum, spin)
ver4 = Parquet.Ver4{Float64}(para, loopNum)

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
# Parquet.showTree(ver4, para, verbose = 1, depth = 3)  # visualize tree using python3 package ete3