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

### Example: 4-point vertex function in FrontEnds

The algorithms in `FrontEnds` generates Feynman diagrams for weak coupling expansions, supporting various diagram types, such as self-energy, polarization, 3-point vertex function, and 4-point vertex function. The internal degrees of freedom can be either the loop variables (e.g., momentum or frequency) or the site variables (e.g., imaginary-time or lattice site).

The example below demonstrates generating one-loop 4-point vertex function diagrams by `Parquet` and `GV`, and optimizing and visualizing their computational graphs, respectively.

#### Generated by `Parquet`:
```julia
julia> using FeynmanDiagram
julia> import FeynmanDiagram.FrontEnds: NoHartree
# Define a parameter structure for a 4-vertex diagram with one-loop in the momentum and the imaginary-time representation. Require the diagrams to be green's function irreducible.
julia> para = Parquet.DiagPara(type = Parquet.Ver4Diag, innerLoopNum = 1, hasTau = true, filter=[NoHartree,]);

# Generate the Feynman diagrams in a DataFrame using the parquet algorithm. `ver4df` is a DataFrame containing fields :response, :type, :extT, :diagram, and :hash.
julia> ver4df = Parquet.build(para) 
6×5 DataFrame
 Row │ diagram                            extT          hash   response  type
     │ Graph…                             Tuple…        Int64  Response  Analytic…
─────┼─────────────────────────────────────────────────────────────────────────────
   1 │ 18671PHr ↑↑Dyn,t(1, 1, 2, 2)=0.0…  (1, 1, 2, 2)  18671  UpUp      Dynamic
   2 │ 18753PPr ↑↑Dyn,t(1, 2, 1, 2)=0.0…  (1, 2, 1, 2)  18753  UpUp      Dynamic
   3 │ 18713PHEr ↑↑Dyn,t(1, 2, 2, 1)=0.…  (1, 2, 2, 1)  18713  UpUp      Dynamic
   4 │ 18667PHr ↑↓Dyn,t(1, 1, 2, 2)=0.0…  (1, 1, 2, 2)  18667  UpDown    Dynamic
   5 │ 18750PPr ↑↓Dyn,t(1, 2, 1, 2)=0.0…  (1, 2, 1, 2)  18750  UpDown    Dynamic
   6 │ 18709PHEr ↑↓Dyn,t(1, 2, 2, 1)=0.…  (1, 2, 2, 1)  18709  UpDown    Dynamic

# Merge the Feynman diagrams in a DataFrame with the same :extT field.
julia> ver4df_extT = Parquet.mergeby(ver4df, [:extT], name=:vertex4)
3×3 DataFrame
 Row │ diagram                            extT          hash
     │ Graph…                             Tuple…        Int64
─────┼────────────────────────────────────────────────────────
   1 │ 18754-vertex4((1, 1, 2, 2),)=0.0…  (1, 1, 2, 2)  18754
   2 │ 18755-vertex4((1, 2, 1, 2),)=0.0…  (1, 2, 1, 2)  18755
   3 │ 18756-vertex4((1, 2, 2, 1),)=0.0…  (1, 2, 2, 1)  18756

# Optimize the Graph for the given Feynman diagrams.
julia> optimize!(ver4df_extT.diagram); 
julia> ver4_parquet = ver4df_extT.diagram[1]; # extract the first Feynman diagram with (1, 1, 2, 2) extT from the DataFrame.

# Visualize the generated first Feynman diagram and its computational graph using ete3.
julia> plot_tree(ver4_parquet)
```

The generated computational graph is shown in the following figure. The leaves of the tree are the propagators (labeled with `G`) and the charge-charge interactions (labeled with `ccIns`).
![tree](assets/ver4_parquet.svg?raw=true "Diagram Tree")

#### Generated by `GV`:
```julia
julia> using FeynmanDiagram
julia> import FeynmanDiagram.FrontEnds: NoHartree
# Generate the 4-vertex diagram with one-loop using the `GV` module. 
julia> ver4_vec, extTs, responses = GV.eachorder_ver4diag(1, filter=[NoHartree])
(Graph{Float64, Float64}[144PPr ↑↑Dyn,t(1, 2, 1, 2)=0.0=⨁ , 142PPr ↑↓Dyn,t(1, 2, 1, 2)=0.0=⨁ , 147PHEr ↑↑Dyn,t(1, 2, 2, 1)=0.0=⨁ , 145PHEr ↑↓Dyn,t(1, 2, 2, 1)=0.0=⨁ , 150PHr ↑↑Dyn,t(1, 1, 2, 2)=0.0=⨁ , 148PHr ↑↓Dyn,t(1, 1, 2, 2)=0.0=⨁ ], [(1, 2, 1, 2), (1, 2, 1, 2), (1, 2, 2, 1), (1, 2, 2, 1), (1, 1, 2, 2), (1, 1, 2, 2)], FeynmanDiagram.FrontEnds.Response[FeynmanDiagram.FrontEnds.UpUp, FeynmanDiagram.FrontEnds.UpDown, FeynmanDiagram.FrontEnds.UpUp, FeynmanDiagram.FrontEnds.UpDown, FeynmanDiagram.FrontEnds.UpUp, FeynmanDiagram.FrontEnds.UpDown])

julia> optimize!(ver4_vec); # Optimize the compuatational graphs.

julia> ver4_GV = ver4_vec[5] + ver4_vec[6] # `+` operator sums up two graphs with the same extT (1, 1, 2, 2).
76=0.0=⨁

julia> plot_tree(ver4_GV)
```

The generated computational graph is shown in the following figure. The leaves of the tree are the propagators and the charge-charge interactions (labeled with `ccIns`).
![tree](assets/ver4_GV.svg?raw=true "Diagram Tree")

#### Construct renormalized Feynman diagrams by Taylor-mode AD

The example code below demonstrates how to build the renormalized Feynman diagrams for the 4-vertex function with the Green's function counterterms and the interaction counterterms using Taylor-mode AD.

```julia
julia> using FeynmanDiagram
# Set the renormalization orders. The first element is the order of Feynman diagrams, the second element is the order of the Green's function counterterms, and the second element is the order of the interaction counterterms.
julia> renormalization_orders = [(1,0,0), (1,1,0), (1,0,1), (1,2,0), (1,1,1), (1,0,2)];

# Generate the Dict of Graph for the renormalized 4-vertex function with the Green's function counterterms and the interaction counterterms.
julia> dict_ver4 = Parquet.diagdict_parquet(Parquet.Ver4Diag, renormalization_orders, filter=[FrontEnds.NoHartree])
Dict{Tuple{Int64, Int64, Int64}, Tuple{Vector{Graph}, Vector{Vector{Int64}}, Vector{FeynmanDiagram.FrontEnds.Response}}} with 6 entries:
  (1, 2, 0) => ([18893[2]=0.0=⨁ , 19004[2]=0.0=⨁ , 19152[2]=0.0=⨁ , 19246[2]=0.0=⨁ , 19377[2]=0.0=⨁ , 19513[2]=0.0=⨁ ], [[1, 1, 2, 2], [1, 2, 1, 2], [1, 2,…
  (1, 0, 2) => ([18887[0, 2]=0.0=Ⓧ , 18998[0, 2]=0.0=Ⓧ , 19146[0, 2]=0.0=Ⓧ , 19240[0, 2]=0.0=Ⓧ , 19371[0, 2]=0.0=Ⓧ , 19507[0, 2]=0.0=Ⓧ ], [[1, 1, 2, 2], [1…
  (1, 1, 0) => ([18891[1]=0.0=⨁ , 19002[1]=0.0=⨁ , 19150[1]=0.0=⨁ , 19244[1]=0.0=⨁ , 19375[1]=0.0=⨁ , 19511[1]=0.0=⨁ ], [[1, 1, 2, 2], [1, 2, 1, 2], [1, 2,…
  (1, 1, 1) => ([18907[1, 1]=0.0=⨁ , 19018[1, 1]=0.0=⨁ , 19166[1, 1]=0.0=⨁ , 19260[1, 1]=0.0=⨁ , 19391[1, 1]=0.0=⨁ , 19527[1, 1]=0.0=⨁ ], [[1, 1, 2, 2], [1…
  (1, 0, 1) => ([18903[0, 1]=0.0=Ⓧ , 19014[0, 1]=0.0=Ⓧ , 19162[0, 1]=0.0=Ⓧ , 19256[0, 1]=0.0=Ⓧ , 19387[0, 1]=0.0=Ⓧ , 19523[0, 1]=0.0=Ⓧ ], [[1, 1, 2, 2], [1…
  (1, 0, 0) => ([18881=0.0=Ⓧ , 18992=0.0=Ⓧ , 19140=0.0=Ⓧ , 19234=0.0=Ⓧ , 19365=0.0=Ⓧ , 19501=0.0=Ⓧ ], [[1, 1, 2, 2], [1, 2, 1, 2], [1, 2, 2, 1], [1, 1, 2, …
```

### Example: BackEnds's `Compilers` for 4-point vertex function
The Back End architecture enables the compiler to output source code in a range of other programming languages and machine learning frameworks. The example code below demonstrates how to use the `Compilers` to generate the source code for the 4-point vertex function in Julia, C, and Python.

```julia
julia> g_o111 = dict_ver4[(1,1,1)] # select the 4-vertex function with both 1st-order Green's function and interaction counterterms.

# Compile the 4-vertex function to Julia RuntimeGeneratedFunction `func` and the `leafmap`, which maps the indices in the vector of leaf values to the corresponding leafs (propagators and interactions). 
julia> func, leafmap = Compilers.compile(g_o111);

# Generate the source code for the 4-vertex function in Julia and save it to a file.
julia> Compilers.compile_Julia(g_o111, "func_g111.jl");

# Generate the source code for the 4-vertex function in C-language and save it to a file.
julia> Compilers.compile_C(g_o111, "func_g111.c");

# Generate the source code for the 4-vertex function in Python and JAX machine learning framework and save it to a file.
julia> Compilers.compile_Python(g_o111, :jax, "func_g111.py");
```

### Computational Graph visualization
#### Tree visualization using ETE
Install the [ete3](http://etetoolkit.org/) Python package for tree-type visualization.

Note that we rely on [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) to call the ete3 python functions from julia. Therefore, you have to install ete3 python package use the python distribution associated with PyCall. According to the tutorial of PyCall, "by default on Mac and Windows systems, Pkg.add("PyCall") or Pkg.build("PyCall") will use the Conda.jl package to install a minimal Python distribution (via Miniconda) that is private to Julia (not in your PATH). You can use the Conda Julia package to install more Python packages, and import Conda to print the Conda.PYTHONDIR directory where python was installed. On GNU/Linux systems, PyCall will default to using the python3 program (if any, otherwise python) in your PATH."

#### Graph visualization using DOT
Back End's `Compilers` support compile the computataional graphs to the `dot` file, which can be visualized using the `dot` command in Graphviz programs. The example code below demonstrates how to visualize the computational graph of the 4-point vertex function.

```julia
# Generate the dot file for the Graph of the 4-point vertex function in the above `Parquet` example.
julia> Compilers.compile_dot([ver4_parquet], "ver4_parquet.dot");
julia> run(`dot -Tsvg ver4_parquet.dot -o ver4_parquet.svg`);

# Generate the dot file for the Graph of the 4-point vertex function in the above `GV` example.
julia> Compilers.compile_dot([ver4_GV], "ver4_GV.dot");
julia> run(`dot -Tsvg ver4_GV.dot -o ver4_GV.svg`);
```

The generated computational graphs are shown in the following figures by the `dot` compiler. The left figure is the graph generated by `Parquet`, and the right figure is the graph generated by `GV`.

<p align="center">
  <img src="assets/ver4_parquetdot.svg?raw=true" alt="parquet_dot" width="300"/>
  <img src="assets/ver4_GVdot.svg?raw=true" alt="GV_dot" width="300"/>
</p>


## License
`FeynmanDiagram.jl` is open-source, available under the [MIT License](https://opensource.org/licenses/MIT). For more details, see the `license.md` file in the repository.