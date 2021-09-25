using StatsPlots
default(; frame=:box, dpi=500, size=(500,500), aspectratio=1, c=:lapaz)

# C - common links
# R - rewired links
# T - turnover links

# Number of steps at which we evaluate r and d
steps = 200

# Values of d (proportion of different links)
d = LinRange(0.0, 1.0, steps)

# Values of r (proportion of unshared links due to rewiring)
r = LinRange(0.0, 1.0, steps)

# Components as a function of r and d
C = 1.0.-d
R = d .* r'
T = d .* (1.0 .- r')

OS = R ./ (2.0.*C .+ R)
WN = (R .+ T) ./ (2.0.*C .+R .+T)
ST = WN .- OS

contour(r, d, OS, levels=9, fill=true)
xaxis!((0,1), "Proportion of rewired links")
yaxis!((0,1), "Proportion of unshared links")
