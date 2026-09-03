using Documenter
using BioBlockcodes

# Load the plotting stack so the BioBlockcodesGraphMakie extension is active and
# the `@example` blocks in the "Plotting graphs" section can render figures.
using CairoMakie, GraphMakie, NetworkLayout

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
