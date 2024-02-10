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

Based on these insights, the compiler's architecture is designed to process Feynman diagrams via a structured, three-stage procedure, mirroring the design principles of the advanced LLVM compiler architecture. This strategy not only broadens the compiler's versatility for diverse QFT computations but also markedly enhances its computational efficiency. The procedure unfolds in three distinct stages:

- **Front End**: This initial phase generates Feynman diagrams and converts them into a unified intermediate representation, shaping them as static computational graphs. The current implementation includes two diagram-generated types: `GV` and `Parquet`.
    - `GV` module is used to generate Feynman diagrams, which are regrouped into a
much smaller number of sign-blessed groups to boost the computational efficiency. The grouping is achieved by leveraging the fermionic crossing symmetry and the free-energy generating functional. Details see [Nat Commun 10, 3725 (2019)](https://doi.org/10.1038/s41467-019-11708-6).
    - `Parquet` module is used to build the 4-point vertex function from a bottom-up approach using the parquet equations and generate other Feynman diagrams coupling with the Dyson-Schwinger equations.
  
  Users can also incorporate new diagram types by extending the front end.

- **Intermediate Representation**: At this stage, the compiler applies optimizations and AD to the static computational graph, facilitating thorough analysis and refinement. The optimizations aim to streamline the graph by eliminating redundant elements, flattening nested chains, and removing zero-valued nodes, among other enhancements. The Taylor-mode AD --- an algorithm for efficient high-order derivatives calculations --- is crucial for deriving specific heat diagrams, RG flow equations, or renormalized Feynman diagrams.

- **Back End**: The final phase focuses on converting the optimized graph into executable code, compatible with a wide range of computing environments. This stage supports various programming languages and ensures seamless integration with different software and hardware ecosystems, thereby extending the compiler's utility across multiple platforms.

## Usage

### Installation

To install the package, you can simply use the following command in the Julia package manager:

```julia
julia> using Pkg
julia> Pkg.add("FeynmanDiagram.jl")
```
### Example: Self-energy in FrontEnds

The algorithms in `FrontEnds` generates Feynman diagrams for weak coupling expansions, supporting various diagram types, such as self-energy, polarization, 3-point vertex function, and 4-point vertex function. The internal degrees of freedom can be either the loop variables (e.g., momentum or frequency) or the site variables (e.g., imaginary-time or lattice site).

The example below demonstrates generating two-loop self-energy diagrams by `Parquet`, and optimizing and visualizing their computational graphs.

#### Generated by `Parquet`:
```julia
julia> using FeynmanDiagram
julia> import FeynmanDiagram.FrontEnds: NoHartree
# Define a parameter structure for two-loop self-energy diagrams in the momentum and the imaginary-time representation. Require the diagrams to be green's function irreducible.
julia> para = Parquet.DiagPara(type = Parquet.SigmaDiag, innerLoopNum = 2 hasTau = true, filter=[NoHartree,]);

# Generate the Feynman diagrams in a DataFrame using the parquet algorithm. `sigmadf` is a DataFrame containing fields :type, :extT, :diagram, and :hash.
julia> sigmadf = Parquet.build(para) 
2×4 DataFrame
 Row │ diagram                            extT    hash   type
     │ Graph…                             Tuple…  Int64  Analytic…
─────┼─────────────────────────────────────────────────────────────
   1 │ 18721-ΣIns, k[1.0, 0.0, 0.0], t(…  (1, 1)  18721  Instant
   2 │ 18722-ΣDyn, k[1.0, 0.0, 0.0], t(…  (1, 2)  18722  Dynamic

# Optimize the Graph for the given Feynman diagrams.
julia> optimize!(sigmadf.diagram); 

# Visualize the computational graph as a tree using ETE.
julia> plot_tree(sigmadf.diagram) 
```

The generated computational graphs are shown in the following figures as trees, in which the first one is the instant self-energy diagram, and the second one is the dynamic self-energy diagram. The leaves of the tree are the propagators (labeled with `G`) and the charge-charge interactions (labeled with `ccIns`).
![tree](assets/sigmaIns_ete.svg?raw=true "Diagram Tree")
![tree](assets/sigma_ete.svg?raw=true "Diagram Tree")

#### Construct renormalized Feynman diagrams by Taylor-mode AD

The example code below demonstrates how to build the renormalized Feynman diagrams for the self-energy with the Green's function counterterms and the interaction counterterms using Taylor-mode AD.

```julia
julia> using FeynmanDiagram
# Set the renormalization orders. The first element is the order of Feynman diagrams, the second element is the order of the Green's function counterterms, and the second element is the order of the interaction counterterms.
julia> renormalization_orders = [(2,0,0), (2,1,0), (2,0,1), (2,2,0), (2,1,1), (2,0,2)];

# Generate the Dict of Graph for the renormalized self-energy diagrams with the Green's function counterterms and the interaction counterterms.
julia> dict_sigma = Parquet.diagdict_parquet(Parquet.SigmaDiag, renormalization_orders, filter=[FrontEnds.NoHartree])
Dict{Tuple{Int64, Int64, Int64}, Tuple{Vector{Graph}, Vector{Vector{Int64}}}} with 6 entries:
  (2, 0, 1) => ([19215[0, 1]=0.0=⨁ , 19581[0, 1]=0.0=⨁ ], [[1, 1], [1, 2]])
  (2, 0, 0) => ([19207=0.0=⨁ , 19573=0.0=⨁ ], [[1, 1], [1, 2]])
  (2, 0, 2) => ([19209[0, 2]=0.0=⨁ , 19576[0, 2]=0.0=⨁ ], [[1, 1], [1, 2]])
  (2, 2, 0) => ([19210[2]=0.0=⨁ , 19575[2]=0.0=⨁ ], [[1, 1], [1, 2]])
  (2, 1, 0) => ([19211[1]=0.0=⨁ , 19577[1]=0.0=⨁ ], [[1, 1], [1, 2]])
  (2, 1, 1) => ([19212[1, 1]=0.0=⨁ , 19578[1, 1]=0.0=⨁ ], [[1, 1], [1, 2]])
```

### Example: BackEnds's `Compilers` for self-energy
The Back End architecture enables the compiler to output source code in a range of other programming languages and machine learning frameworks. The example code below demonstrates how to use the `Compilers` to generate the source code for the self-energy diagrams in Julia, C, and Python.

```julia
julia> g_o211 = dict_sigma[(2,1,1)] # select the self-energy with 2nd-order Green's function counterterms and 1st-order interaction counterterms.

# Compile the self-energy to Julia RuntimeGeneratedFunction `func` and the `leafmap`, which maps the indices in the vector of leaf values to the corresponding leafs (propagators and interactions). 
julia> func, leafmap = Compilers.compile(g_o211);

# Generate the source code for the self-energy in Julia and save it to a file.
julia> Compilers.compile_Julia(g_o211, "func_g211.jl");

# Generate the source code for the self-energy in C-language and save it to a file.
julia> Compilers.compile_C(g_o211, "func_g211.c");

# Generate the source code for the self-energy in Python and JAX machine learning framework and save it to a file.
julia> Compilers.compile_Python(g_o211, :jax, "func_g211.py");
```

### Computational Graph visualization
#### Tree-type visualization using ETE
Install the [ete3](http://etetoolkit.org/) Python package for tree-type visualization.

Note that we rely on [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) to call the ete3 python functions from julia. Therefore, you have to install ete3 python package use the python distribution associated with PyCall. According to the tutorial of PyCall, "by default on Mac and Windows systems, Pkg.add("PyCall") or Pkg.build("PyCall") will use the Conda.jl package to install a minimal Python distribution (via Miniconda) that is private to Julia (not in your PATH). You can use the Conda Julia package to install more Python packages, and import Conda to print the Conda.PYTHONDIR directory where python was installed. On GNU/Linux systems, PyCall will default to using the python3 program (if any, otherwise python) in your PATH."

#### Graph visualization using DOT
Back End's `Compilers` support compile the computataional graphs to the `dot` file, which can be visualized using the `dot` command in Graphviz programs. The example code below demonstrates how to visualize the computational graph of the self-energy in the above `Parquet` example.

```julia
# Generate the dot file for the Graph of the self-energy in the above `Parquet` example.
julia> Compilers.compile_dot(sigmadf.diagram, "sigma_o2.dot");
julia> run(`dot -Tsvg sigma_o2.dot -o sigma_o2.svg`);
```

The generated computational graphs are shown in the following figure as a directed acyclic graph, by the `dot` compiler. The leaves are represented by the green oval nodes and the intermediate nodes are represented by the rectangular nodes. In the penultimate floor of the graph, the left blue node with Sum operator represents the dynamic self-energy, and the right blue node represents the instant self-energy. 

![graph](assets/sigma_o2.svg?raw=true "Graph")

## License
`FeynmanDiagram.jl` is open-source, available under the [MIT License](https://opensource.org/licenses/MIT). For more details, see the `license.md` file in the repository.