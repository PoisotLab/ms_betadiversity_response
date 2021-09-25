# Numerical experiment 5 - we vary the connectance of one of the networks
# relative to that of the other

q = LinRange(0.0, 1.0, 120)
p = 0.75
a = 10.0.^LinRange(0.0, 1.0, 97)'

# The expressions are all calculated in main text
A = p.^2.0 .* (1.0.-q)
U = (1.0.+a) .- p.^2.0 .- a.*p.^2.0
S = (1.0.+a).*(p.^2.0) .- 2.0.*(p.^2.0).*(1.0.-q)

OS = S ./ (2A .+ S)
WN = (S.+U) ./ (2A .+ S .+ U)

ST = WN .- OS

p1 = contour(log10.(vec(a)), q, OS, aspectratio=1, fill=true, levels=9, lc=:white, clim=(0,1))
xaxis!(p1, "Connectance ratio (log10)", (0, 1))
yaxis!(p1, "Rewiring probability", (0, 1))
title!(p1, "βos")

p2 = contour(log10.(vec(a)), q, ST, aspectratio=1, fill=true, levels=9, lc=:white, clim=(0,1))
xaxis!(p2, "Connectance ratio (log10)", (0, 1))
yaxis!(p2, "Rewiring probability", (0, 1))
title!(p2, "βst")

p3 = contour(log10.(vec(a)), q, WN, aspectratio=1, fill=true, levels=9, lc=:white, clim=(0,1))
xaxis!(p3, "Connectance ratio (log10)", (0, 1))
yaxis!(p3, "Rewiring probability", (0, 1))
title!(p3, "βwn")

p4 = contour(log10.(vec(a)), q, ST./WN, aspectratio=1, fill=true, levels=9, lc=:white, clim=(0,1))
xaxis!(p4, "Connectance ratio (log10)", (0, 1))
yaxis!(p4, "Rewiring probability", (0, 1))
title!(p4, "βst / βwn ")

plot(p1, p2, p3, p4, size=(800, 800), dpi=600)

savefig(joinpath(_FIGPATH, "numexp5.png"))










q = 0.15
p = LinRange(0.0, 1.0, 120)
a = 10.0.^LinRange(0.0, 1.0, 97)'

# The expressions are all calculated in main text
A = p.^2.0 .* (1.0.-q)
U = (1.0.+a) .- p.^2.0 .- a.*p.^2.0
S = (1.0.+a).*(p.^2.0) .- 2.0.*(p.^2.0).*(1.0.-q)

OS = S ./ (2A .+ S)
WN = (S.+U) ./ (2A .+ S .+ U)

ST = WN .- OS

p1 = contour(log10.(vec(a)), p, OS, aspectratio=1, fill=true, levels=9, lc=:white, clim=(0,1))
xaxis!(p1, "Connectance ratio (log10)", (0, 1))
yaxis!(p1, "Species sharing probability", (0, 1))
title!(p1, "βos")

p2 = contour(log10.(vec(a)), p, ST, aspectratio=1, fill=true, levels=9, lc=:white, clim=(0,1))
xaxis!(p2, "Connectance ratio (log10)", (0, 1))
yaxis!(p2, "Species sharing probability", (0, 1))
title!(p2, "βst")

p3 = contour(log10.(vec(a)), p, WN, aspectratio=1, fill=true, levels=9, lc=:white, clim=(0,1))
xaxis!(p3, "Connectance ratio (log10)", (0, 1))
yaxis!(p3, "Species sharing probability", (0, 1))
title!(p3, "βwn")

p4 = contour(log10.(vec(a)), p, ST./WN, aspectratio=1, fill=true, levels=9, lc=:white, clim=(0,1))
xaxis!(p4, "Connectance ratio (log10)", (0, 1))
yaxis!(p4, "Species sharing probability", (0, 1))
title!(p4, "βst / βwn ")

plot(p1, p2, p3, p4, size=(800, 800), dpi=600)

savefig(joinpath(_FIGPATH, "numexp6.png"))
