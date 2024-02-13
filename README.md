# FeynmanDiagram

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://numericalEFT.github.io/FeynmanDiagram.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://numericalEFT.github.io/FeynmanDiagram.jl/dev)
[![Build Status](https://github.com/numericalEFT/FeynmanDiagram.jl/workflows/CI/badge.svg)](https://github.com/numericalEFT/FeynmanDiagram.jl/actions)
[![Coverage](https://codecov.io/gh/numericalEFT/FeynmanDiagram.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/numericalEFT/FeynmanDiagram.jl)

`FeynmanDiagram.jl` is a Julia package designed to efficiently encode Feynman diagrams --- essential elements of Quantum Field Theory (QFT) —-- into compact computational graphs for fast evaluation. It employs Taylor-mode Automatic Differentiation (AD) specifically to implement field-theoretic renormalization schemes, a pivotal technique in QFT that significantly improves the convergence of Feynman diagrammatic series. This approach underscores the synergy between QFT and AI technologies, effectively addressing the sophisticated computational challenges in QFT.

## Key Features

- **Computational Graphs for Feynman Diagrams**: Utilizes computational graphs analogous to those in Machine Learning (ML) architectures, where nodes represent mathematical operations and edges denote the tensor flow.
  
- **Compiler Architecture**: Implements a three-stage compiler process, analogous to modern programming language compilers, to transform Feynman diagrams into executable code, significantly enhancing the adaptability and computational efficiency of QFT calculations.

- **Interdisciplinary Approach**: Bridges QFT and AI tech stack by adapting ML algorithms, such as Taylor-mode AD, for QFT calculations, enhancing computational efficiency with tools like JAX, TensorFlow, PyTorch, and Mindspore.

## Compiler Architecture Overview

In general, Feynman diagrams represents high-order integral. The integrand are propagators/interactions composed by the basis arithmetic operations (multiplication, addition, power, etc). The sequence of calculating the integrand by combining the propagators/interactions with the arithmetic operatos can be represented as an algebraic computational graph. In this sense, the computational graph serves as an intermediate representation that standardizes the treatment of various diagram types, ensuring a consistent approach across different QFT calculations.

![infrastructure](assets/diagram_compiler.svg?raw=true "Compiler Infrastructure")

Drawing from these insights, the architecture of our compiler is intentionally crafted to process Feynman diagrams via a strategic, three-stage process that reflects the advanced design principles of the LLVM compiler architecture. This approach enhances the compiler's flexibility for a wide array of QFT computations and significantly improves its computational efficiency. The workflow is organized into three critical stages:

- **Front End**: In this initial phase, Feynman diagrams are generated and transformed into a standardized intermediate representation, taking the form of static computational graphs. This stage features two key algorithms for generating generic Feynman diagrams for weak-coupling expansions:
   - `Parquet` module systematically organizes higher-order Feynman diagrams into a concise hierarchical structure of sub-diagrams, enabling efficient evaluation of repeated sub-diagrams. This module introduces algorithms for constructing streamlined computational graphs for two-, three-, and four-point vertex functions by exploiting the perturbative representations of the Dyson-Schwinger and parquet equations. This approach facilitates the analysis of a wide array of observables in quantum many-body problems.
   - `GV` module, specifically aimed at many-electron problems with Coulomb interactions, incorporates the algorithm proposed in [Nat Commun 10, 3725 (2019)](https://doi.org/10.1038/s41467-019-11708-6).
   - The front end also allows for the integration of new diagram types by users, enhancing its adaptability.

- **Intermediate Representation**:  At this stage, the compiler applies optimizations and Automatic Differentiation (AD) to the static computational graph. This process is geared towards refining the graph for thorough analysis and optimization. The optimizations are aimed at streamlining the graph by removing redundant elements, flattening nested chains, and eliminating zero-valued nodes, among other strategies. The incorporation of Taylor-mode AD is critical for efficiently calculating high-order derivatives, essential for deriving diagrams for specific heat, RG flow equations, or renormalized Feynman diagrams.

- **Back End**: This final phase is responsible for translating the optimized graph into executable code that is compatible with a broad spectrum of computing environments. It supports various programming languages and facilitates seamless integration with different software and hardware ecosystems, significantly extending the compiler's utility across multiple platforms.

## Usage

### Installation

To install the package, you can simply use the following command in the Julia package manager:

```julia
julia> using Pkg
julia> Pkg.add("FeynmanDiagram.jl")
```

### Example: Self-energy in FrontEnds

The algorithms in `FrontEnds` generates Feynman diagrams for weak coupling expansions, supporting various diagram types, such as self-energy, polarization, 3-point vertex function, and 4-point vertex function. The internal degrees of freedom can be either the loop variables (e.g., momentum or frequency) or the site variables (e.g., imaginary-time or lattice site).

The example below demonstrates generating two-loop self-energy diagrams by `Parquet` and optimizing their computational graphs.

#### Generate Diagram generation with the `Parquet` algorithm

```julia
julia> using FeynmanDiagram
julia> import FeynmanDiagram.FrontEnds: NoHartree

# Define a parameter for two-loop self-energy diagrams in the momentum and the imaginary-time representation. Exclude any diagrams containing Hartree subdiagrams. 
julia> para = Parquet.DiagPara(type = Parquet.SigmaDiag, innerLoopNum = 2, hasTau = true, filter=[NoHartree,]);

# Construct Feynman diagrams within a DataFrame utilizing the parquet algorithm. The resulting sigmadf DataFrame comprises two components: the instantaneous part and the dynamic part of the self-energy.
julia> sigmadf = Parquet.build(para) 
2×4 DataFrame
 Row │ diagram                            extT    hash   type
     │ Graph…                             Tuple…  Int64  Analytic…
─────┼─────────────────────────────────────────────────────────────
   1 │ 18721-ΣIns, k[1.0, 0.0, 0.0], t(…  (1, 1)  18721  Instant
   2 │ 18722-ΣDyn, k[1.0, 0.0, 0.0], t(…  (1, 2)  18722  Dynamic

# Optimize the Graph for the given Feynman diagrams.
julia> optimize!(sigmadf.diagram); 
```

#### Construct renormalized Feynman diagrams using the Taylor-mode AD

The example code below demonstrates how to build the renormalized Feynman diagrams for the self-energy with the Green's function counterterms and the interaction counterterms using Taylor-mode AD.

```julia
# Set the renormalization orders. The first element is the order of the Green's function counterterms, and the second element is the order of the interaction counterterms.
julia> renormalization_orders = [2, 3];

# Generate the Dict of Graph for the renormalized self-energy diagrams with the Green's function counterterms and the interaction counterterms.
julia> dict_sigma = taylorAD(sigmadf.diagram, [para.innerLoopNum, para.innerLoopNum], renormalization_orders);
```

### Example: Compile Feynman diagrams to different programming languages 
The Back End architecture enables the compiler to output source code in a range of other programming languages and machine learning frameworks. The example code below demonstrates how to use the `Compilers` to generate the source code for the self-energy diagrams in Julia, C, and Python.

```julia
#Access the two-loop self-energy data for the configuration with 1st-order Green's function counterterms and 1st-order interaction counterterms.
julia> g_o211 = dict_sigma[[2,1,1]]; 

# Compile the selected self-energy into a Julia RuntimeGeneratedFunction `func` and a `leafmap`.
# The `leafmap` associates vector indices of leaf values with their corresponding leaves (propagators and interactions). 
julia> func, leafmap = Compilers.compile(g_o211);

# Export the self-energy's source code to a Julia file.
julia> Compilers.compile_Julia(g_o211, "func_g211.jl");

# Export the self-energy's source code to a C file.
julia> Compilers.compile_C(g_o211, "func_g211.c");

# Export the self-energy's source code to a Python file for use with the JAX machine learning framework.
julia> Compilers.compile_Python(g_o211, :jax, "func_g211.py");
```

### Computational Graph visualization
#### Tree-type visualization using ETE
To visualize tree-type structures of the self-energy in the Parquet example, install the [ete3](http://etetoolkit.org/) Python package, a powerful toolkit for tree visualizations.

Execute the following command in Julia for tree-type visualization of the self-energy generated in the above `Parquet` example:
```julia
julia> plot_tree(sigmadf.diagram, depth = 3)
```
For installation instructions on using ete3 with [PyCall.jl](https://github.com/JuliaPy/PyCall.jl), please refer to the PyCall.jl documentation on how to configure and use external Python packages within Julia.

#### Graph visualization using DOT
The Back-End's `Compilers` module includes functionality to translate computational graphs into .dot files, which can be visualized using the `dot` command from Graphviz software. Below is an example code snippet that illustrates how to visualize the computational graph of the self-energy from the previously mentioned `Parquet` example.

```julia
# Convert the self-energy graph from the `Parquet` example into a dot file.
julia> Compilers.compile_dot(sigmadf.diagram, "sigma_o2.dot");
# Use Graphviz's dot command to create an SVG visualization of the computational graph.
julia> run(`dot -Tsvg sigma_o2.dot -o sigma_o2.svg`);
```
The resulting computational graphs are depicted as directed acyclic graphs by the dot compiler. In these visualizations, the leaves are indicated by green oval nodes, while the intermediate nodes take on a rectangular shape. On the penultimate level of the graph, the left blue node with the Sum operator signifies the dynamic part of the self-energy, and the right blue node denotes the instantaneous part of the self-energy.
![graph](assets/sigma_o2.svg?raw=true "Graph")

## License
`FeynmanDiagram.jl` is open-source, available under the [MIT License](https://opensource.org/licenses/MIT). For more details, see the `license.md` file in the repository.