using ExpressionTree
using AbstractTrees
# using NewickTree
using StaticArrays
const Weight = SVector{2,Float64}

Parquet = Builder.Parquet

chan = [Parquet.T, Parquet.U, Parquet.S]
interactionTauNum = 1
loopNum = 2
spin = 2

para = Parquet.Para(3, chan, interactionTauNum, (Float64, Float64), (Float64, Float64), (Weight, Float64), spin)

ver4 = Parquet.Ver4{Float64}(para, loopNum, 1)

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

para = Parquet.Para(3, chan, interactionTauNum, (Float64, Float64), (Float64, Float64), (Float64, Float64), spin)
K0 = zeros(2 + loopNum)
KinL = KoutL = KinR = KoutR = deepcopy(K0)
KinL[1] = KoutL[1] = 1
KinR[2] = KoutR[2] = 2
legK = [KinL, KoutL, KinR, KoutR]

varK = [rand(3) for i in 1:2+loopNum]
varT = [rand() for i in 1:2*(loopNum+1)]
evalK(basis) = sum([basis[i] * varK[i] for i in 1:3])
evalT(Tidx) = varT[Tidx]
# diag, ver4 = Parquet.diagramTree(para, 1, legK, 3, 1, Float64, Gsym, Wsym, spin)
# ver4 = Parquet.Ver4{Int}(para, 1, 1)
# println("ver4444...")
# println(ver4)
diag, ver4, dir, ex = Parquet.ver4toDiagTree(para, loopNum, legK, 3, 1, evalK, evalT, 1.0)
rootDir = DiagTree.addNode(diag, DiagTree.ADD, [[], []], dir; para = "dir")
rootEx = DiagTree.addNode(diag, DiagTree.ADD, [[], []], ex; para = "ex")
diag.root = [rootDir, rootEx]

DiagTree.showTree(diag, rootDir)

# print_tree(ver4)
# # nodeDi = length(diag.tree) - 2 #the last second node is for ve*vec
# # nodeEx = length(diag.tree)
# # GWKT.DiagTree.showTree(diag, nodeDi)
# # GWKT.DiagTree.showTree(diag, nodeEx)
# println(diag.root)
# for root in diag.root
#     GWKT.DiagTree.showTree(diag, root)
# end
