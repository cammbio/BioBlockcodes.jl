using Documenter
using BioBlockcodes

DocMeta.setdocmeta!(BioBlockcodes, :DocTestSetup, :(using BioBlockcodes); recursive = true)

makedocs(;
    modules = [BioBlockcodes],
    authors = "Filip Cavar, Markus Gumbel",
    sitename = "BioBlockcodes.jl",
    format = Documenter.HTML(;
        # canonical = "https://fcavar.github.io/BioBlockcodes.jl",
        edit_link = "master",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/cammbio/BioBlockcodes.jl", devbranch = "master")
