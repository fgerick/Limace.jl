using Limace
using Documenter, DocStringExtensions
using Documenter: MathJax3, DocMeta
using PlutoStaticHTML

const NOTEBOOK_DIR = joinpath(@__DIR__, "src", "notebooks")


DocMeta.setdocmeta!(Limace, :DocTestSetup, :(using Limace); recursive=true)

"""
    build()

Run all Pluto notebooks (".jl" files) in `NOTEBOOK_DIR`.
"""
function build()
    println("Building notebooks in $NOTEBOOK_DIR")
    oopts = OutputOptions(; append_build_context=false)
    output_format = documenter_output
    bopts = BuildOptions(NOTEBOOK_DIR; output_format)
    build_notebooks(bopts, oopts)
    return nothing
end

# Build the notebooks; defaults to true.
if get(ENV, "BUILD_DOCS_NOTEBOOKS", "true") == "true"
    build()
end

pages= [
    "Home" => "index.md",
    "Examples" => [
        "Inviscid inertial modes" => "notebooks/inertialmodes.md",
        "Torsional Alfvén modes" => "notebooks/torsionalmodes.md"
    ],
    "Bases" => "bases.md",
    "Forces" => "forces.md",
    "Poly" => "poly.md",
    # "Misc" => "misc.md"
]


makedocs(;
    modules=[Limace],
    repo="https://github.com/fgerick/Limace.jl/blob/{commit}{path}#{line}",
    sitename="Limace.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        mathengine = MathJax3(),
        canonical="https://fgerick.github.io/Limace.jl",
        edit_link="main",
        assets=String[],
        size_threshold_ignore = collect(values(Dict(Dict(pages)["Examples"])))
    ),
    pages,
)

deploydocs(;
    repo="github.com/fgerick/Limace.jl",
    devbranch="main",
)
