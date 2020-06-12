using Documenter, Qwind

makedocs(;
    modules=[Qwind],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/arnauqb/Qwind.jl/blob/{commit}{path}#L{line}",
    sitename="Qwind.jl",
    authors="Arnau Quera-Bofarull",
    assets=String[],
)

deploydocs(;
    repo="github.com/arnauqb/Qwind.jl",
)
