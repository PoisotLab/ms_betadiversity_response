# Numerical experiment 4 - difference in rewiring value

steps = 150

p_values = LinRange(0.0, 1.0, steps) # Fraction of shared species
q_values = LinRange(0.0, 1.0, steps) # Fraction of rewired links

# Stores
OS = zeros(Float64, length(p_values), length(q_values))
FR = similar(OS)
WN = similar(OS)

for (i, p) in enumerate(p_values)
    for (j, q) in enumerate(q_values)
        A = p * (1 - q)
        S = p * q
        U = (1 - p) * 2
        OS[i,j] = S / (2A + S)
        FR[i,j] = S / (2A + S + U)
        WN[i,j] = (S+U) / (2A + S + U)
    end
end

ST = WN .- OS
FT = WN .- FR

Δ = OS .- FR

p1 = heatmap(q_values, p_values, OS, clim=(0, 1), aspectratio=1, fill=true, lc=:white)
xaxis!(p1, "Rewiring probability", (0, 1))
yaxis!(p1, "Species sharing probability", (0, 1))
title!(p1, "βos")

p2 = heatmap(q_values, p_values, FR, clim=(0, 1), aspectratio=1, fill=true, lc=:white)
xaxis!(p2, "Rewiring probability", (0, 1))
yaxis!(p2, "Species sharing probability", (0, 1))
title!(p2, "βos (common denominator)")

p3 = heatmap(q_values, p_values, FT, clim=(0, 1), aspectratio=1, fill=true, lc=:white)
xaxis!(p3, "Rewiring probability", (0, 1))
yaxis!(p3, "Species sharing probability", (0, 1))
title!(p3, "βst (common denominator)")

p4 = heatmap(q_values, p_values, OS .- FR, clim=(0, 1), aspectratio=1, fill=true, lc=:white, c=:grayC)
xaxis!(p4, "Rewiring probability", (0, 1))
yaxis!(p4, "Species sharing probability", (0, 1))
title!(p4, "Difference in rewiring estimate")

plot(p1, p2, p3, p4, size=(800, 800), dpi=600)

savefig(joinpath(_FIGPATH, "numexp4.png"))