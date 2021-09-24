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

p1 = contour(q_values, p_values, OS, clim=(0, 1), aspectratio=1, fill=true, lc=:white)
xaxis!(p1, "Rewiring probability", (0, 1))
yaxis!(p1, "Species sharing probability", (0, 1))
title!(p1, "βos")

p2 = contour(q_values, p_values, ST, clim=(0, 1), aspectratio=1, fill=true, lc=:white)
xaxis!(p2, "Rewiring probability", (0, 1))
yaxis!(p2, "Species sharing probability", (0, 1))
title!(p2, "βst")

p3 = contour(q_values, p_values, WN, clim=(0, 1), aspectratio=1, fill=true, lc=:white)
xaxis!(p3, "Rewiring probability", (0, 1))
yaxis!(p3, "Species sharing probability", (0, 1))
title!(p3, "βwn")

p4 = contour(q_values, p_values, ST./WN, clim=(0, 1), aspectratio=1, fill=true, lc=:white)
xaxis!(p4, "Rewiring probability", (0, 1))
yaxis!(p4, "Species sharing probability", (0, 1))
title!(p4, "βst / βwn")

plot(p1, p2, p3, p4, size=(800, 800), dpi=600)

savefig(joinpath(_FIGPATH, "numexp2.png"))