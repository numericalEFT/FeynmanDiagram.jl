# FeynmanDiagram

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://numericalEFT.github.io/FeynmanDiagram.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://numericalEFT.github.io/FeynmanDiagram.jl/dev)
[![Build Status](https://github.com/numericalEFT/FeynmanDiagram.jl/workflows/CI/badge.svg)](https://github.com/numericalEFT/FeynmanDiagram.jl/actions)
[![Coverage](https://codecov.io/gh/numericalEFT/FeynmanDiagram.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/numericalEFT/FeynmanDiagram.jl)

This package implements a mini-compiler that compiles generic Feynman diagrams into expression tree representations for fast computation. 

## Infrastructure

In general, Feynman diagrams represents high-order integral. The integrand are propagators/interactions composed by the basis arithmetic operations (multiplication, addition). The sequence of calculating the integrand by combining the propagators/interactions with the arithmetic operatos can be represented as an algebraic expression tree. In this sense, the expression tree provides an intermediate representation (IR) for Feynman diagrams that completely independent of the diagram type. 

![infrastructure](assets/diagram_compiler.jpeg?raw=true "Compiler Infrastructure")

Base on this observation, we develop a package to compile the integrand of Feynman diagrams into machine code so that one can evaluate the it efficiently. The infrastructure of this package is similar to the modern compiler LLVM for generic programming language. There are three layers: a front-end translates a source code into an IR as an expression tree, then a mid-end optimizes and transforms the IR, and a back-end to compiles the IR to machine code. 

- The front-end supports Feynman diagrams from weak coupling expansion or strong coupling expansion. The user can incorprate new types of diagrams by writing their own front-end.

- The mid-end performs universal optimizations and transformations of one expression tree to another. The possible optimizations of the expression tree includes: remove common nodes/leaves, remove zero-valued nodes/leaves, merge small nodes into a large one. The possible transformations include automatic differentiation (which can be useful to derive the diagrams for the specific heat, RG flow equation, etc.), renormalization of the propagators and the interactions, and analytic Matsubara-frequency integration (work in progress).

- The back-end provides a universal subroutine to evalue the expression tree efficiently. 

## Supported Front-end

### 1. Generic Weak Coupling Expansion based on the Parquet Algorithm

This algorithm generates the Feynman diagrams of weak coupling expansion. It supports the diagrams of self-energy, polarization, 3-point vertex function and 4-point vertex function. The internal degrees of freedom can be either the loop variables (e.g., momentum or frequency) or the site variables (e.g., imaginary-time or lattice site).

The main idea of the algorithm is to use the parquet equation to build high-order-vertex-function diagrams from the lower order sub-diagrams. 

The following code is a simple example to generate the one-loop 4-point vertex function diagrams, then visualize the expression tree.

```julia
using FeynmanDiagram

# Define a parameter structure for the 4-vertex diagram with one-loop, in the momentum and the imaginary-time representation. Require the diagrams to be green's function irreducible.
para = GenericPara(diagType =Ver4Diag, innerLoopNum = 1,hasTau = true, filter=[Girreducible,])

ver4=Parquet.build(para) #build the diagram tree with the parquet algorithm.

plot_tree(ver4) # visualize the generated diagram tree

tree=ExprTree.build(ver4.diagram) #optimize the diagram tree to get an optimized expression tree
```

The generated diagram tree is as shown in the following figure. The leaves of the tree are the propagators (labeled with `G`) and the interactions (labeled with `Ins`). By default, the interactions is assumed to spin-symmetric. A typical example is the Coulomb interaction.
![tree](assets/ver4tree.png?raw=true "Diagram Tree")


### 2. Generic Strong Coupling Expansion (work in progress)
### 3. Hand-drawing Feynman diagrams (work in progress)

## Expression Tree visualization
To visualize the diagram tree, you need to install the ete3 python3 package (http://etetoolkit.org/).

Note that we rely on "PyCall.jl" to call the ete3 python functions from julia. Therefore, you have to install ete3 python package use the python distribution associated with PyCall. According to the tutorial of PyCall (https://github.com/JuliaPy/PyCall.jl), "by default on Mac and Windows systems, Pkg.add("PyCall") or Pkg.build("PyCall") will use the Conda.jl package to install a minimal Python distribution (via Miniconda) that is private to Julia (not in your PATH). You can use the Conda Julia package to install more Python packages, and import Conda to print the Conda.PYTHONDIR directory where python was installed. On GNU/Linux systems, PyCall will default to using the python3 program (if any, otherwise python) in your PATH."

