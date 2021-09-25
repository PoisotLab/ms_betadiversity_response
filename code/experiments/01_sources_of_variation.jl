
# Numerical experiment 1: L1 = L2, we vary the number of shared, unique, and
# non-unique links, and examine the consequences on the values of the OS, ST,
# and WN components of beta diversity. The goal of this experiment is to show
# that, without considering connectance or species turnover, we can tie the
# values of these components to the sources of network dissimilarity.

steps = 100

# These are proportions - first, we pick a proportion of links (out of L1 or L2)
# that are shared, and then we divide the rest of the links in those that are
# not-shared and rewired (S in the manuscript), and those that are not-shared
# and involve non-overlapping species (U in the paper). Note that the species
# richness does not matter here.

shared_links = LinRange(0.0, 1.0, steps)
links_rewired = LinRange(0.0, 1.0, steps)

# We prepare matrices to fill the calculated values
OS = zeros(Float64, steps, steps)
WN = similar(OS)

# Finally, we measure OS and WN - we will get the rest by substraction
for (i, sl) in enumerate(shared_links)
    for (j, lr) in enumerate(links_rewired)
        common = sl
        rewired = (1.0 - common) * lr
        turned = 1.0 - (common + rewired)
        # We can assume that all links are from the same network, since the
        # measure is symetrical
        OS[i,j] = beta_t(common, 0.0, rewired)
        WN[i,j] = beta_t(common, turned, rewired)
    end
end

# We get the values of ST
ST = WN .- OS

# The following lines are plots of the results of numerical experiment 1

p1 = contour(links_rewired, shared_links, OS, clim=(0, 1), aspectratio=1, fill=true, levels=9, lc=:white)
xaxis!(p1, "Proportion of rewired links", (0, 1))
yaxis!(p1, "Proportion of shared links", (0, 1))
title!(p1, "βos")

p2 = contour(links_rewired, shared_links, ST, clim=(0, 1), aspectratio=1, fill=true, levels=9, lc=:white)
xaxis!(p2, "Proportion of rewired links", (0, 1))
yaxis!(p2, "Proportion of shared links", (0, 1))
title!(p2, "βst")

p3 = contour(links_rewired, shared_links, WN, clim=(0, 1), aspectratio=1, fill=true, levels=9, lc=:white)
xaxis!(p3, "Proportion of rewired links", (0, 1))
yaxis!(p3, "Proportion of shared links", (0, 1))
title!(p3, "βwn")

p4 = contour(links_rewired, shared_links, ST ./ WN, clim=(0, 1), aspectratio=1, fill=true, levels=9, lc=:white)
xaxis!(p4, "Proportion of rewired links", (0, 1))
yaxis!(p4, "Proportion of shared links", (0, 1))
title!(p4, "βst / βwn")

plot(p1, p2, p3, p4, size=(800, 800), dpi=600)

savefig(joinpath(_FIGPATH, "numexp1.png"))