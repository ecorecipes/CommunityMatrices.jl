using Documenter
using CommunityMatrices

tutorials = sort(filter(f -> endswith(f, ".md"), readdir(joinpath(@__DIR__, "src", "tutorials"))))
titles = Dict("01_plant_communities.md" => "Plant communities: structural stability and sign patterns",
              "02_gut_microbiota.md" => "A gut microbiota under antibiotics")

# Source and edit links need a git remote; outside a git checkout, build without them.
ingit = success(pipeline(`git -C $(@__DIR__) rev-parse`; stdout = devnull, stderr = devnull))

makedocs(;
         modules = [CommunityMatrices],
         authors = "Simon Frost",
         sitename = "CommunityMatrices.jl",
         checkdocs = :exports,
         remotes = ingit ? Dict() : nothing,
         format = Documenter.HTML(;
                                  edit_link = ingit ? "main" : nothing,
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://ecorecipes.github.io/CommunityMatrices.jl",
                                  size_threshold = 400 * 1024),
         pages = ["Home" => "index.md",
                  "Tutorials" => [get(titles, f, f) => "tutorials/$f" for f in tutorials],
                  "API Reference" => ["Interaction matrices and Lotka–Volterra" => "api/matrices.md",
                                      "Structural stability and sign patterns" => "api/structural.md",
                                      "Time series" => "api/timeseries.md",
                                      "Empirical data" => "api/data.md",
                                      "Metrics from ResilienceCore" => "api/metrics.md"]])

if get(ENV, "CI", "false") == "true"
    deploydocs(; repo = "github.com/ecorecipes/CommunityMatrices.jl.git", push_preview = true)
end
