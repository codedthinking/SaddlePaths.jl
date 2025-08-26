using SaddlePaths
using Documenter

DocMeta.setdocmeta!(SaddlePaths, :DocTestSetup, :(using SaddlePaths); recursive=true)

makedocs(;
    modules=[SaddlePaths],
    authors="Koren",
    sitename="SaddlePaths.jl",
    format=Documenter.HTML(;
        canonical="https://koren.github.io/SaddlePaths.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/koren/SaddlePaths.jl",
    devbranch="main",
)
