using Documenter
using SaddlePaths

makedocs(
    sitename = "SaddlePaths.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://yourusername.github.io/SaddlePaths.jl",
        assets = String[],
    ),
    modules = [SaddlePaths],
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "DSL Reference" => "dsl.md",
        "Examples" => "examples.md",
        "API Reference" => "api.md",
        "Theory" => "theory.md",
    ]
)

deploydocs(
    repo = "github.com/yourusername/SaddlePaths.jl.git",
)