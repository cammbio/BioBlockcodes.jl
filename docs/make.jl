using Documenter
using GCATCodes

DocMeta.setdocmeta!(GCATCodes, :DocTestSetup, :(using GCATCodes); recursive = true)

makedocs(;
    modules = [GCATCodes],
    authors = "Filip Cavar",
    sitename = "GCATCodes.jl",
    format = Documenter.HTML(;
        canonical = "https://fcavar.github.io/GCATCodes.jl",
        edit_link = "master",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/fcavar/GCATCodes.jl", devbranch = "master")
