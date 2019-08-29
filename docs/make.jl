using Documenter, OceanPlot

makedocs(modules = [OceanPlot], sitename = "OceanPlot.jl")

deploydocs(
    repo = "github.com/gher-ulg/OceanPlot.jl.git",
)
