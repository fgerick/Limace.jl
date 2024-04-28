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

makedocs(;
    modules=[Limace],
    authors="Felix <felixgerick@gmail.com> and contributors",
    repo="https://github.com/fgerick/Limace.jl/blob/{commit}{path}#{line}",
    
    sitename="Limace.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        mathengine = MathJax3(),
        canonical="https://fgerick.github.io/Limace.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Poly" => "poly.md",
        "Bases" => ["Basis definition" => "bases/bases.md", 
                    "InsulatingBasis" => "bases/insulating.md",
                    "InviscidBasis" => "bases/inviscid.md",
                    "ViscousBasis" => "bases/viscous.md",
        ],
        "Forces" => ["Summary" => "forces/index.md", 
                    "Coriolis" => "forces/coriolis.md",
                    "Diffusion" => "forces/diffusion.md",
                    "Induction" => "forces/induction.md",
                    "Inertia" => "forces/inertial.md",
                    "Lorentz" => "forces/lorentz.md",
        ],
        "Examples" => ["Inviscid inertial modes" => "notebooks/inertialmodes.md",
        ],
        "Misc" => "misc.md"
    ],
)

deploydocs(;
    repo="github.com/fgerick/Limace.jl",
    devbranch="main",
)
