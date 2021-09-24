# Dependencies for the simulations
using StatsPlots
using PlotThemes
using Random
using StatsBase
using Base.Threads
using Statistics

# This ensures the reproducibility of the results
Random.seed!(90210)
theme(:mute)
default(; frame=:box, dpi=600, size=(500, 500), c=:lapaz)

# β-diversity function we use
βₜ = (a,b,c) -> (b+c)/(2a+b+c)

_FIGPATH = joinpath("figures")
ispath(_FIGPATH) || mkpath(_FIGPATH)

# Numerical experiment 1
include(joinpath("code", "experiments", "01_sources_of_variation.jl"))

# Numerical experiment 2 - at some constant connectance, we examine the effect
# of changing the species overlap between two networks, as well as the
# proportion that overlapping species keep their interactions. This will provide
# evidence that the components of the beta diversity partition are responding to
# the correct ecological processes. Note that, as explained in the manuscript,
# this calculation is independent of network richness, and does not require any
# simulations.

Co = 0.4
p_values = LinRange(0.0, 1.0, 200) # Fraction of shared species
q_values = LinRange(0.0, 1.0, 202) # Fraction of rewired links

# Stores
OS = zeros(Float64, length(p_values), length(q_values))
WN = similar(OS)

for (i, p) in enumerate(p_values)
    for (j, q) in enumerate(q_values)
        A = p * (1 - q) * Co
        S = p * q * Co
        U = Co * (1 - p) * 2
        OS[i,j] = S / (2A + S)
        WN[i,j] = (S + U) / (2A + S + U)
    end
end

ST = WN .- OS

p1 = contour(q_values, p_values, OS, clim=(0, 1), aspectratio=1, fill=true, lc=:white, levels=4)
xaxis!(p1, "Rewiring probability", (0, 1))
yaxis!(p1, "Species sharing probability", (0, 1))
title!(p1, "βos")

p2 = contour(q_values, p_values, ST, clim=(0, 1), aspectratio=1, fill=true, lc=:white, levels=4)
xaxis!(p2, "Rewiring probability", (0, 1))
yaxis!(p2, "Species sharing probability", (0, 1))
title!(p2, "βst")

p3 = contour(q_values, p_values, WN, clim=(0, 1), aspectratio=1, fill=true, lc=:white, levels=4)
xaxis!(p3, "Rewiring probability", (0, 1))
yaxis!(p3, "Species sharing probability", (0, 1))
title!(p3, "βwn")

p4 = contour(q_values, p_values, ST./WN, clim=(0, 1), aspectratio=1, fill=true, lc=:white, levels=4)
xaxis!(p4, "Rewiring probability", (0, 1))
yaxis!(p4, "Species sharing probability", (0, 1))
title!(p4, "βst / βwn")

plot(p1, p2, p3, p4, size=(800, 800), dpi=600)
savefig("numexp2.png")

# Numerical experiment 3 - we vary the connectance of one of the networks
# relative to that of the other

Co1 = 0.25
Co2 = LinRange(0.0, 1.0, 200)

p = LinRange(0.0, 1.0, 204)'
q = 0.15

# The expressions are all calculated in main text
A = p.*(1.0.-q).*min.(Co1,Co2)
U = (1.0.-p).*Co1 .+ (1.0.-p).*Co2
S = Co1 .- (A .+ (1.0.-p).*Co1) .+ Co2 .- (A .+ (1.0.-p).*Co2)

OS = S ./ (2A .+ S)
WN = (S.+U) ./ (2A .+ S .+ U)
ST = WN .- OS

p1 = contour(Co2, p', OS', aspectratio=1, fill=true, lc=:white, levels=4, clim=(0,1))
vline!([Co1], c=:black, lab="", lw=3, ls=:dash)
xaxis!(p1, "Connectance of N₂", (0, 1))
yaxis!(p1, "Species sharing probability", (0, 1))
title!(p1, "βos")

p2 = contour(Co2, p', ST', aspectratio=1, fill=true, lc=:white, levels=4, clim=(0,1))
vline!([Co1], c=:black, lab="", lw=3, ls=:dash)
xaxis!(p2, "Connectance of N₂", (0, 1))
yaxis!(p2, "Species sharing probability", (0, 1))
title!(p2, "βst")

p3 = contour(Co2, p', WN', aspectratio=1, fill=true, lc=:white, levels=4, clim=(0,1))
vline!([Co1], c=:black, lab="", lw=3, ls=:dash)
xaxis!(p3, "Connectance of N₂", (0, 1))
yaxis!(p3, "Species sharing probability", (0, 1))
title!(p3, "βwn")

p4 = contour(Co2, p', (ST./WN)', aspectratio=1, fill=true, lc=:white, levels=4, clim=(0,1))
vline!([Co1], c=:black, lab="", lw=3, ls=:dash)
xaxis!(p4, "Connectance of N₂", (0, 1))
yaxis!(p4, "Species sharing probability", (0, 1))
title!(p4, "βst / βwn")

plot(p1, p2, p3, p4, size=(800, 800), dpi=600)
savefig("numexp3.png")
