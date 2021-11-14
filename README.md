# ExpressionTree

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://numericalEFT.github.io/ExpressionTree.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://numericalEFT.github.io/ExpressionTree.jl/dev)
[![Build Status](https://github.com/numericalEFT/ExpressionTree.jl/workflows/CI/badge.svg)](https://github.com/numericalEFT/ExpressionTree.jl/actions)
[![Coverage](https://codecov.io/gh/numericalEFT/ExpressionTree.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/numericalEFT/ExpressionTree.jl)
<!-- [![Coverage](https://codecov.io/gh/houpc/ExpressionTree.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/houpc/ExpressionTree.jl) -->

This library includes several algorithms to generate the optimized expression tree representations of high-order Feynman diagrams. Here we give a tutorial on the algorithms.

## Parquet Algorithm

The ExpressionTree.Parquet model uses parquet equations to generate high-order Feynman diagrams. In this representation, a 4-point one-particle-irreducible (1PI) diagram of loop order N is a product of a left and a right 4-point 1PI subdiagrams (loop order M and loop order N-M where M=0, 1, ..., N). Note that there are three ways to glue the left and right subdiagrams: the particle-hole channel (Parquet.T), the particle-hole exchange channel (Parquet.U) and the particle-particle channel (Parquet.S). Each of them could have a counter-term (Parquet.TC, Parquet.UC and Parquet.SC).

```julia
using AbstractTrees, StaticArrays, ExpressionTree

# Define the type of the weight of a 4-vertex diagram
const Weight = SVector{2,Float64} 

# Define the possibilities to build a 4-vertex diagram from the left and right 4-vertex subdiagrams.
chan = [Parquet.T, Parquet.U, Parquet.S] 

# Generate the parameter for the expression tree. The second argument lists the possible number of imaginary-time variables in the bare 4-vertex (namely, the bare interaction of your model). For example, the instaneous Coulomb interaction only has one time variable, while the retared effective interaction has two time variables.
para = Parquet.Para(chan, [1, 2]) 

# Generate an expression tree for a set of 4-vertex diagrams with loop order 1, initial imaginary-time index 1, and the parameter set para.
ver4 = Parquet.Ver4{Weight}(1, 1, para) 

# use AbstractTrees interface to print/manipulate the tree
print_tree(ver4)

# [println(node) for node in Leaves(ver4)]  #print all loopNum=0 ver4
# [println(node) for node in PostOrderDFS(ver4)]  # iterator ver4 in depth-first search (children before parents)
# [println(node) for node in PreOrderDFS(ver4)]  # iterator ver4 (parents before children)

# You can also use ete3 python3 package to visualize tree
Parquet.showTree(ver4, para, verbose = 1, depth = 3)

# You can also print tree to a newick format file, there are many visualization software for the newick format
io = open("./test.newick", "w")
write(io, Parquet.newick(ver4))
close(io)
```

Run the above script will print out the tree on the terminal,
![terminal](docs/figures/terminal_example.png?raw=true "Terminal Ouptut")
where the 4-pair (t1, t2, t3, t4) gives the imaginary-time indices of the four legs (left in, left out, right in, right out). For each channel, the loop order of the left and the right subdiagrams is given.

and generate a visualization of tree using the ete3 python3 package,
![ete](docs/figures/ete_example.png?raw=true "Ete3 visualization")
where "lp" stands for the loop order of the current set of the 4-vertex diagrams.

