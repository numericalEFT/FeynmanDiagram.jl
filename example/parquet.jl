using ExpressionTree
using AbstractTrees
# using NewickTree
using StaticArrays
const Weight = SVector{2,Float64}

chan = [Parquet.T, Parquet.U, Parquet.S]

para = Parquet.Para(chan, [1, 2])

ver4 = Parquet.Ver4{Weight}(1, 1, para)

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

Parquet.diagramTree(ver4)
