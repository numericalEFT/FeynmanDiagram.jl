using ExpressionTree
using Documenter

DocMeta.setdocmeta!(ExpressionTree, :DocTestSetup, :(using ExpressionTree); recursive = true)

makedocs(;
    modules = [ExpressionTree],
    authors = "Kun Chen, Pengcheng Hou",
    repo = "https://github.com/numericalEFT/ExpressionTree.jl/blob/{commit}{path}#{line}",
    sitename = "ExpressionTree.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://numericaleft.github.io/ExpressionTree.jl",
        assets = String[]
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
        ],
        "API reference" => Any[
            "lib/diagtree.md",
            "lib/parquet.md",
        ]
    ]
)

deploydocs(;
    repo = "github.com/numericalEFT/ExpressionTree.jl"
)
