using FeynmanDiagram
using Documenter

DocMeta.setdocmeta!(FeynmanDiagram, :DocTestSetup, :(using FeynmanDiagram); recursive=true)

makedocs(;
    modules=[FeynmanDiagram],
    authors="Kun Chen, Pengcheng Hou",
    repo="https://github.com/numericalEFT/FeynmanDiagram.jl/blob/{commit}{path}#{line}",
    sitename="FeynmanDiagram.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://numericaleft.github.io/FeynmanDiagram.jl",
        assets=["assets/custom.css"]
    ),
    checkdocs=:exports, # check only exported names within the modules
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "manual/counterterms.md",
            "manual/feynman_rule.md",
            "manual/interaction.md",
            "manual/hubbard_atom.md"
        ],
        "API reference" => Any[
            "lib/operator.md",
            "lib/computgraph.md",
            "lib/taylorseries.md",
            "lib/frontend.md",
            "lib/GV.md",
            "lib/parquet.md",
            "lib/backend.md",
            "lib/utility.md"
        ]
    ]
)

deploydocs(;
    repo="github.com/numericalEFT/FeynmanDiagram.jl",
    push_preview=true
)
