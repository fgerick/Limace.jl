using Limace
using Documenter, DocStringExtensions
using Documenter: MathJax3, DocMeta
using DocumenterCitations
using Literate

DocMeta.setdocmeta!(Limace, :DocTestSetup, :(using Limace); recursive=true)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:authoryear
)

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src", "examples")

examples = [
    "inertialmodes.jl",
    "torsionalmodes.jl"
]

if !(@isdefined LiveServer)
    for example in examples
        example_filepath = joinpath(EXAMPLES_DIR, example)
        Literate.markdown(example_filepath, OUTPUT_DIR)
        # Literate.notebook(example_filepath, OUTPUT_DIR)
    end
end

pages= [
    "Home" => "index.md",
    "Theoretical background" => "theory.md",
    "Examples" => [
        "Inviscid inertial modes" => "examples/inertialmodes.md",
        "Torsional Alfvén modes" => "examples/torsionalmodes.md"
    ],
    "Bases" => "bases.md",
    "Forces" => "forces.md",
    "Poly" => "poly.md",
    "Postprocessing" => "processing.md",
    "References" => "references.md",
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
    plugins=[bib],
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/fgerick/Limace.jl",
    devbranch="main",
)
