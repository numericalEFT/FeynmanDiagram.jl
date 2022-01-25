var documenterSearchIndex = {"docs":
[{"location":"lib/builder/#Diagram-Tree-Builder","page":"Diagram Tree Builder","title":"Diagram Tree Builder","text":"","category":"section"},{"location":"lib/builder/","page":"Diagram Tree Builder","title":"Diagram Tree Builder","text":"Modules = [FeynmanDiagram.Builder]","category":"page"},{"location":"lib/parquet/#Diagram-Tree-Builder-based-on-the-Parquet-Algorithm","page":"Diagram Tree Builder based on the Parquet Algorithm","title":"Diagram Tree Builder based on the Parquet Algorithm","text":"","category":"section"},{"location":"lib/parquet/","page":"Diagram Tree Builder based on the Parquet Algorithm","title":"Diagram Tree Builder based on the Parquet Algorithm","text":"Modules = [FeynmanDiagram.Builder.Parquet]","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = FeynmanDiagram","category":"page"},{"location":"#FeynmanDiagram","page":"Home","title":"FeynmanDiagram","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for FeynmanDiagram.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [FeynmanDiagram]","category":"page"},{"location":"#FeynmanDiagram.innerTauNum-Tuple{FeynmanDiagram.DiagramType, Any, Any}","page":"Home","title":"FeynmanDiagram.innerTauNum","text":"function innerTauNum(diagType::DiagramType, innerLoopNum, interactionTauNum)\n\ninternal imaginary-time degrees of freedom for a given diagram type and internal loop number.\nFor the vertex functions (self-energy, polarization, vertex3, and vertex4), innerTauNum is equivalent to tauNum.\nFor the Green function, tauNum = innerTauNum + external tauNum\n\n\n\n\n\n","category":"method"},{"location":"#Library-Outline","page":"Home","title":"Library Outline","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n    \"lib/parquet.md\",\n    \"lib/builder.md\",\n    \"lib/diagtree.md\"\n]\nDepth = 2","category":"page"},{"location":"lib/diagtree/#Expression-Tree-Representation-of-Feynman-Diagrams","page":"Expression Tree Representation of Feynman Diagrams","title":"Expression Tree Representation of Feynman Diagrams","text":"","category":"section"},{"location":"lib/diagtree/","page":"Expression Tree Representation of Feynman Diagrams","title":"Expression Tree Representation of Feynman Diagrams","text":"Modules = [FeynmanDiagram.DiagTree]","category":"page"},{"location":"lib/diagtree/#FeynmanDiagram.DiagTree.plot_tree-Tuple{Diagram}","page":"Expression Tree Representation of Feynman Diagrams","title":"FeynmanDiagram.DiagTree.plot_tree","text":"showTree(diag::Diagrams, _root = diag.root[end]; verbose = 0, depth = 999)\n\nVisualize the diagram tree using ete3 python package\n\n#Arguments\n\ndiag: the Diagrams struct to visualize\n_root: the index of the root node to visualize\nverbose=0: the amount of information to show\ndepth=999: deepest level of the diagram tree to show\n\n\n\n\n\n","category":"method"}]
}
