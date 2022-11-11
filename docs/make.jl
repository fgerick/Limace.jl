using Limace
using Documenter

DocMeta.setdocmeta!(Limace, :DocTestSetup, :(using Limace); recursive=true)

makedocs(;
    modules=[Limace],
    authors="Felix <felixgerick@gmail.com> and contributors",
    repo="https://github.com/fgerick/Limace.jl/blob/{commit}{path}#{line}",
    sitename="Limace.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fgerick.github.io/Limace.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fgerick/Limace.jl",
    devbranch="main",
)
