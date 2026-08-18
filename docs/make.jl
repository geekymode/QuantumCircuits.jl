using Documenter
using QuantumCircuits
using CairoMakie          # loads the plotting extension so its docs build

DocMeta.setdocmeta!(QuantumCircuits, :DocTestSetup, :(using QuantumCircuits); recursive=true)

makedocs(
    sitename = "QuantumCircuits.jl",
    modules  = [QuantumCircuits],
    authors  = "Rethna Pulikkoonattu",
    repo     = Documenter.Remotes.GitHub("geekymode", "QuantumCircuits.jl"),
    # prettyurls=false outside CI emits page.html instead of page/index.html, so
    # a local build browses correctly straight off the filesystem (file:// URLs).
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical  = "https://geekymode.github.io/QuantumCircuits.jl",
        collapselevel = 1,
    ),
    pages = [
        "Home"          => "index.md",
        "Gray coding"   => "graycode.md",
        "Linear algebra" => "math.md",
        "Illustrations" => "plots.md",
        "API reference" => "api.md",
    ],
    checkdocs = :exports,
    doctest   = true,
)

deploydocs(
    repo      = "github.com/geekymode/QuantumCircuits.jl.git",
    devbranch = "main",
)
