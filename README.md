# FeynmanDiagram

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://numericalEFT.github.io/FeynmanDiagram.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://numericalEFT.github.io/FeynmanDiagram.jl/dev)
[![Build Status](https://github.com/numericalEFT/FeynmanDiagram.jl/workflows/CI/badge.svg)](https://github.com/numericalEFT/FeynmanDiagram.jl/actions)
[![Coverage](https://codecov.io/gh/numericalEFT/FeynmanDiagram.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/numericalEFT/FeynmanDiagram.jl)

`FeynmanDiagram.jl` is a Julia package that compiles Feynman diagrams in Quantum Field Theory (QFT) into a compact computational graph representation for fast evalution. It utilizes Taylor-mode Automatic Differentiation (AD) for field-theoretic renormalization, showcasing the synergy between QFT and AI tech stack to address computational challenges in QFT.

## Key Features

- **Computational Graphs for Feynman Diagrams**: Utilizes computational graphs analogous to those in Machine Learning (ML) architectures, where nodes represent mathematical operations and edges denote the tensor flow.
  
- **Compiler Architecture**: Implements a three-stage compiler process, analogous to modern programming language compilers, to transform Feynman diagrams into executable code, significantly enhancing the adaptability and computational efficiency of QFT calculations.

- **Interdisciplinary Approach**: Bridges QFT and AI tech stack by adapting ML algorithms, such as Taylor-mode AD, for QFT calculations, enhancing computational efficiency with tools like JAX, TensorFlow, PyTorch, and Mindspore.

## Compiler Architecture Overview

In general, Feynman diagrams represents high-order integral. The integrand are propagators/interactions composed by the basis arithmetic operations (multiplication, addition, power, etc). The sequence of calculating the integrand by combining the propagators/interactions with the arithmetic operatos can be represented as an algebraic computational graph. In this sense, the computational graph serves as an intermediate representation that standardizes the treatment of various diagram types, ensuring a consistent approach across different QFT calculations.

![infrastructure](assets/diagram_compiler.svg?raw=true "Compiler Infrastructure")

Base on this observation, the compiler architecture is designed to process Feynman diagrams through a structured, three-stage procedure, drawing parallels with the architectures of the modern compiler LLVM. This approach not only enhances the adaptability of the compiler for various QFT calculations but also significantly boosts computational efficiency. The three-stage procedure is as follows:

- **Front End**: Generates the specific Feynman diagrams and maps them into a unified intermediate representation as a computational graph. It identifies diagram elements, mapping them into corresponding nodes. Users can also incorporate new diagram types by extending the front end.

- **Intermediate Representation**: Applies optimizations and automatic differentiation to the static computational graph, enabling comprehensive analysis and optimization. Optimizations include removing redundant nodes/leaves, flattens chains, merges linear combinations, removing zero-valued nodes, and so on. AD is used to derive the diagrams for the specific heat, RG flow equation, etc, or the renormalized Feynman diagrams. The incorporated Taylor-mode AD is utilized for efficient high-order derivative graph computations.

- **Back End**: Translates the optimized graph into executable code for a variety of computing platforms, supporting multiple programming languages and integration with software and hardware ecosystems.

## Usage

### Installation

To install the package, you can simply use the following command in the Julia package manager:

```julia
julia> using Pkg
julia> Pkg.add("FeynmanDiagram.jl")
```

### Example: Weak Coupling Expansion with Parquet Algorithm

This algorithm generates Feynman diagrams for weak coupling expansions, supporting various diagram types, such as self-energy, polarization, 3-point vertex function, and 4-point vertex function. The internal degrees of freedom can be either the loop variables (e.g., momentum or frequency) or the site variables (e.g., imaginary-time or lattice site). The main idea of the algorithm is to use the parquet equation to build high-order-vertex-function diagrams from the lower order sub-diagrams. 

The example below demonstrates generating one-loop 4-point vertex function diagrams and visualizing their computational graph.

```julia
using FeynmanDiagram
import FeynmanDiagram.FrontEnds: NoHartree, Girreducible

# Define a parameter structure for a 4-vertex diagram with one-loop in the momentum and the imaginary-time representation. Require the diagrams to be green's function irreducible.
para = Parquet.DiagPara(type = Parquet.Ver4Diag, innerLoopNum = 1, hasTau = true, filter=[NoHartree, Girreducible,])

# Generate the Feynman diagrams in a DataFrame using the parquet algorithm. `ver4df` is a DataFrame containing fields :response, :type, :extT, :diagram, and :hash.
ver4df = Parquet.build(para) 

ver4df_extT = Parquet.mergeby(ver4df, [:extT], name=:vertex4) # merge the Feynman diagrams in a DataFrame with the same `extT` field.

optimize!(ver4df_extT.diagram) # optimize the Graph for the given Feynman diagrams.

ver4 = ver4df_extT.diagram[1]
plot_tree(ver4) # visualize the generated first Feynman diagram and its computational graph using ete3.
```

The generated computational graph is as shown in the following figure. The leaves of the tree are the propagators (labeled with `G`) and the interactions (labeled with `Ins`). By default, the interactions is assumed to spin-symmetric. A typical example is the Coulomb interaction.
![tree](assets/ver4_ete3.svg?raw=true "Diagram Tree")

### Computational Graph visualization
Install the [ete3](http://etetoolkit.org/) Python package for visualization.

Note that we rely on [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) to call the ete3 python functions from julia. Therefore, you have to install ete3 python package use the python distribution associated with PyCall. According to the tutorial of PyCall, "by default on Mac and Windows systems, Pkg.add("PyCall") or Pkg.build("PyCall") will use the Conda.jl package to install a minimal Python distribution (via Miniconda) that is private to Julia (not in your PATH). You can use the Conda Julia package to install more Python packages, and import Conda to print the Conda.PYTHONDIR directory where python was installed. On GNU/Linux systems, PyCall will default to using the python3 program (if any, otherwise python) in your PATH."

## License
`FeynmanDiagram.jl` is open-source, available under the [MIT License](https://opensource.org/licenses/MIT). For more details, see the `license.md` file in the repository.