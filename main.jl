# Dependencies for the simulations
using EcologicalNetworks
using StatsPlots
using PlotThemes
using Random
using StatsBase
using Base.Threads
using Statistics

# This ensures the reproducibility of the results
Random.seed!(90210)
theme(:mute)
default(frame=:box, dpi=600, size=(500, 500))

# Shortcut notation for the betadiv measure we use
βₜ = EcologicalNetworks.KGL08

# Numerical experiment 1: L1 = L2, we vary the number of shared, unique, and
# non-unique links, and examine the consequences on the values of the OS, ST,
# and WN components of beta diversity. The goal of this experiment is to show
# that, without considering connectance or species turnover, we can tie the
# values of these components to the sources of network dissimilarity.

# The values of L1 and L2 are entirely arbitrary -- indeed, the manuscript shows
# how the values are independant from these values, as we can collect the
# non-shared links in those that are rewired, and those that are connecting
# non-shared species. Using different values for L1 and L2 complexifies the
# problem a little, by introducing some combinatorics questions that would just
# stand in the way.
L1, L2 = 200, 200

# These are proportions - first, we pick a proportion of links (out of L1 or L2)
# that are shared, and then we divide the rest of the links in those that are
# not-shared and rewired (S in the manuscript), and those that are not-shared
# and involve non-overlapping species (U in the paper). Note that the species
# richness does not matter here.
shared_links = LinRange(0.0, 1.0, 200)
links_rewired = LinRange(0.0, 1.0, 202)

# We prepare matrices to fill the calculated values
OS = zeros(Float64, length(shared_links), length(links_rewired))
WN = similar(OS)

# Finally, we measure OS and WN - we will get the rest by substraction
for (i, sl) in enumerate(shared_links)
    for (j, lr) in enumerate(links_rewired)
        common = L1 * sl
        rewired = (L1 - common) * lr
        turned = L1 - (common + rewired)
        # We can assume that all links are from the same network, since the
        # measure is symetrical
        OS[i,j] = βₜ((a = common, b = 0.0, c = rewired))
        WN[i,j] = βₜ((a = common, b = turned, c = rewired))
    end
end

# We get the values of ST
ST = WN .- OS

# The following lines are plots of the results of numerical experiment 1

p1 = contour(links_rewired, shared_links, OS, clim=(0, 1), aspectratio=1, fill=true, lc=:white, levels=4)
xaxis!(p1, "Proportion of rewired links", (0, 1))
yaxis!(p1, "Proportion of shared links", (0, 1))
title!(p1, "βos")

p2 = contour(links_rewired, shared_links, ST, clim=(0, 1), aspectratio=1, fill=true, lc=:white, levels=4)
xaxis!(p2, "Proportion of rewired links", (0, 1))
yaxis!(p2, "Proportion of shared links", (0, 1))
title!(p2, "βst")

p3 = contour(links_rewired, shared_links, WN, clim=(0, 1), aspectratio=1, fill=true, lc=:white, levels=4)
xaxis!(p3, "Proportion of rewired links", (0, 1))
yaxis!(p3, "Proportion of shared links", (0, 1))
title!(p3, "βwn")

p4 = contour(links_rewired, shared_links, ST ./ WN, clim=(0, 1), aspectratio=1, fill=true, lc=:white, levels=4)
xaxis!(p4, "Proportion of rewired links", (0, 1))
yaxis!(p4, "Proportion of shared links", (0, 1))
title!(p4, "βst / βwn")

plot(p1, p2, p3, p4, size=(800, 800), dpi=600)
savefig("numexp1.png")

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
