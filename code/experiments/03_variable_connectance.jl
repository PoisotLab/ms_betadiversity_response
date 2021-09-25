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

p1 = contour(Co2, p', OS', aspectratio=1, fill=true, levels=9, lc=:white, clim=(0,1))
vline!([Co1], c=:black, lab="", lw=3, ls=:dash)
xaxis!(p1, "Connectance of N₂", (0, 1))
yaxis!(p1, "Species sharing probability", (0, 1))
title!(p1, "βos")

p2 = contour(Co2, p', ST', aspectratio=1, fill=true, levels=9, lc=:white, clim=(0,1))
vline!([Co1], c=:black, lab="", lw=3, ls=:dash)
xaxis!(p2, "Connectance of N₂", (0, 1))
yaxis!(p2, "Species sharing probability", (0, 1))
title!(p2, "βst")

p3 = contour(Co2, p', WN', aspectratio=1, fill=true, levels=9, lc=:white, clim=(0,1))
vline!([Co1], c=:black, lab="", lw=3, ls=:dash)
xaxis!(p3, "Connectance of N₂", (0, 1))
yaxis!(p3, "Species sharing probability", (0, 1))
title!(p3, "βwn")

p4 = contour(Co2, p', (ST./WN)', aspectratio=1, fill=true, levels=9, lc=:white, clim=(0,1))
vline!([Co1], c=:black, lab="", lw=3, ls=:dash)
xaxis!(p4, "Connectance of N₂", (0, 1))
yaxis!(p4, "Species sharing probability", (0, 1))
title!(p4, "βst / βwn")

plot(p1, p2, p3, p4, size=(800, 800), dpi=600)

savefig(joinpath(_FIGPATH, "numexp3.png"))
