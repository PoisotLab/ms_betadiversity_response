# Dependencies for the simulations
using StatsPlots
using PlotThemes
using Random
using StatsBase
using Base.Threads
using Statistics

theme(:mute)
default(; frame=:box, dpi=600, size=(500, 500), c=:batlow)

# Path to store the figures
_FIGPATH = joinpath("figures")
ispath(_FIGPATH) || mkpath(_FIGPATH)

# Path to store the experiments
_EXPPATH = joinpath("code", "experiments")

# Run every experiment in the folder
for f in filter(f -> endswith(f, ".jl"), readdir(_EXPPATH))
    @info f
    include(joinpath(_EXPPATH, f))
end

