# FeynmanDiagram

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://numericalEFT.github.io/FeynmanDiagram.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://numericalEFT.github.io/FeynmanDiagram.jl/dev)
[![Build Status](https://github.com/numericalEFT/FeynmanDiagram.jl/workflows/CI/badge.svg)](https://github.com/numericalEFT/FeynmanDiagram.jl/actions)
[![Coverage](https://codecov.io/gh/numericalEFT/FeynmanDiagram.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/numericalEFT/FeynmanDiagram.jl)

This package implements a mini-compiler that compiles generic Feynman diagrams into expression tree representations for fast computation. 

## Infrastructure

In general, Feynman diagrams represents high-order integral. The integrand are propagators/interactions composed by the basis arithmetic operations (multiplication, addition). The sequence of calculating the integrand by combining the propagators/interactions with the arithmetic operatos can be represented as an algebraic expression tree. In this sense, the expression tree provides an intermediate representation (IR) for Feynman diagrams that completely independent of the diagram type. 

![infrastructure](docs/figures/diagram_compiler.jpeg?raw=true "Compiler Infrastructure")
Base on this observation, we develop a package to compile the integrand of Feynman diagrams into machine code so that one can evaluate the it efficiently. The infrastructure of this package is similar to the modern compiler LLVM for generic programming language. There are three layers: a front-end translates a source code into an IR as an expression tree, then a mid-end optimizes and transforms the IR, and a back-end to compiles the IR to machine code. 

- The front-end supports Feynman diagrams from weak coupling expansion or strong coupling expansion. The user can incorprate new types of diagrams by writing their own front-end.

- The mid-end performs universal optimizations and transformations of one expression tree to another. The possible optimizations of the expression tree includes: remove common nodes/leaves, remove zero-valued nodes/leaves, merge small nodes into a large one. The possible transformations include automatic differentiation (which can be useful to derive the diagrams for the specific heat, RG flow equation, etc.), renormalization of the propagators and the interactions, and analytic Matsubara-frequency integration (work in progress).

- The back-end provides a universal subroutine to evalue the expression tree efficiently. 

## Supported Front-end

### 1. Parquet Algorithm for generic Weak Coupling Expansion

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

### 2. Generic Strong Coupling Expansion (work in progress)
### 3. Hand-drawing Feynman diagrams (work in progress)

## Expression Tree visualization
If you have installed ete3 python3 package (http://etetoolkit.org/), you will get a visualization of tree by running the above script,
![ete](docs/figures/ete_example.png?raw=true "Ete3 visualization")
where "lp" stands for the loop order of the current set of the 4-vertex diagrams.

Note that we rely on "PyCall.jl" to call the ete3 python functions from julia. Therefore, you have to install ete3 python package use the python distribution associated with PyCall. According to the tutorial of PyCall (https://github.com/JuliaPy/PyCall.jl), "by default on Mac and Windows systems, Pkg.add("PyCall") or Pkg.build("PyCall") will use the Conda.jl package to install a minimal Python distribution (via Miniconda) that is private to Julia (not in your PATH). You can use the Conda Julia package to install more Python packages, and import Conda to print the Conda.PYTHONDIR directory where python was installed. On GNU/Linux systems, PyCall will default to using the python3 program (if any, otherwise python) in your PATH."

