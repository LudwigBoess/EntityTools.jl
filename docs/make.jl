using EntityTools
using Documenter

DocMeta.setdocmeta!(EntityTools, :DocTestSetup, :(using EntityTools); recursive=true)

makedocs(;
    modules=[EntityTools],
    authors="Ludwig Böss <ludwig.boess@outlook.de> and contributors",
    sitename="EntityTools.jl",
    format=Documenter.HTML(;
        canonical="https://LudwigBoess.github.io/EntityTools.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/LudwigBoess/EntityTools.jl",
    devbranch="main",
)
