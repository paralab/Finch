using Documenter

if !@isdefined(Finch)
    include("../src/Finch.jl");
    using .Finch
end

makedocs(
    sitename="Finch",
    pages = [
        "index.md",
        "start.md",
        "Usage" => [
            "configuration.md",
            "mesh.md",
            "entities.md",
            "conditions.md",
            "equation.md",
            "solution.md",
            "datatypes.md",
            "reorder.md",
            "misc.md"
        ],
        "examples.md" 
    ]
)

deploydocs(
    repo = "github.com/paralab/Finch.git",
)