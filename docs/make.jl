using Limace
using Documenter, DocStringExtensions

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
        "Poly" => "poly.md",
        "Bases" => ["Basis definition" => "bases/bases.md", 
                    "InviscidBasis" => "bases/inviscid.md",
                    "InsulatingBasis" => "bases/insulating.md",
        ],
        "Forces" => ["Summary" => "forces/index.md", 
                    "Coriolis" => "forces/coriolis.md",
                    "Diffusion" => "forces/diffusion.md",
                    "Induction" => "forces/induction.md",
                    "Inertia" => "forces/inertial.md",
                    "Lorentz" => "forces/lorentz.md",
        ],
        "Misc" => "misc.md"
    ],
)

deploydocs(;
    repo="github.com/fgerick/Limace.jl",
    devbranch="main",
)
