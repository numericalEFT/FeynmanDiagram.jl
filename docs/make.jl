using ExpressionTree
using Documenter

DocMeta.setdocmeta!(ExpressionTree, :DocTestSetup, :(using ExpressionTree); recursive=true)

makedocs(;
    modules=[ExpressionTree],
    authors="Pengcheng Hou",
    repo="https://github.com/houpc/ExpressionTree.jl/blob/{commit}{path}#{line}",
    sitename="ExpressionTree.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://houpc.github.io/ExpressionTree.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/numericalEFT/ExpressionTree.jl",
)
