# Dependencies for the simulations
using EcologicalNetworks
using StatsPlots
using PlotThemes
using Random
using StatsBase
using Base.Threads

# This ensures the reproducibility of the results
Random.seed!(90210)
theme(:mute)
default(frame=:box, dpi=600, size=(500,500), c=:lapaz)

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
links_rewired = LinRange(0.0, 1.0, 220)

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

p1 = contour(links_rewired, shared_links, OS, clim=(0, 1), aspectratio=1, fill=true, lc=:white, levels=5)
xaxis!(p1, "Proportion of rewired links", (0, 1))
yaxis!(p1, "Proportion of shared links", (0, 1))
title!(p1, "βos")

p2 = contour(links_rewired, shared_links, ST, clim=(0, 1), aspectratio=1, fill=true, lc=:white, levels=5)
xaxis!(p2, "Proportion of rewired links", (0, 1))
yaxis!(p2, "Proportion of shared links", (0, 1))
title!(p2, "βst")

p3 = contour(links_rewired, shared_links, WN, clim=(0, 1), aspectratio=1, fill=true, lc=:white, levels=5)
xaxis!(p3, "Proportion of rewired links", (0, 1))
yaxis!(p3, "Proportion of shared links", (0, 1))
title!(p3, "βwn")

p4 = contour(links_rewired, shared_links, ST ./ WN, clim=(0, 1), aspectratio=1, fill=true, lc=:white, levels=5)
xaxis!(p4, "Proportion of rewired links", (0, 1))
yaxis!(p4, "Proportion of shared links", (0, 1))
title!(p4, "βst / βwn")

plot(p1, p2, p3, p4, size=(800,800), dpi=600)
savefig("numexp1.png")

# Numerical experiment 2 - at some constant connectance, we examine the effect
# of changing the species overlap between two networks, as well as the
# proportion that overlapping species keep their interactions. This will provide
# evidence that the components of the beta diversity partition are responding to
# the correct ecological processes.

# We will set a network of fixed size and connectance
M = simplify(rand(BipartiteProbabilisticNetwork(rand(Float64, 35, 35) .* 0.3)))
r = richness(M; dims=1)

# The second network will be generated with a number of shared species ranging
# from 0 (maximal dissimilarity) to r (no species turnover), and values of
# rewiring probability ranging from 0 to 1 - because we will rely on random
# draws, there will be a number of replicates to ensure some smooth results.
shared_nodes = 0:r
proba_rewiring = LinRange(0.0, 1.0, 35)
replicates = 20

# The storage objects are tensors, with a third dimension added for replicate
OS = zeros(Float64, length(shared_nodes), length(proba_rewiring), replicates)
WN = similar(OS)

# We will run the number of overlapping species sequentially
for (sni, sn) in enumerate(shared_nodes)
    # But the probability of rewiring is threaded
    Threads.@threads for sli in 1:length(proba_rewiring)
        sl = proba_rewiring[sli]
        # And then the replicates are ran sequentially within a thread
        for ri in 1:replicates
            shared_species = sample(species(M; dims=1), sn; replace=false)
            turned_species = "t" .* string.(((r + 1):(r + (r - length(shared_species)))).+4)
            new_species = vcat(shared_species, turned_species)
            N = BipartiteNetwork(zeros(Bool, size(M)), new_species, species(M; dims=2))
            sgraph = M[shared_species,:]
            # Rewiring
            newint = eltype(M)[]
            oldint = interactions(simplify(sgraph))
            for oint in oldint
                if rand() < sl
                    s1, s2 = rand(species(sgraph; dims=1)), rand(species(sgraph; dims=2))
                    tpl = (from = s1, to = s2)
                    while tpl in vcat(newint, oldint)
                        s1, s2 = rand(species(sgraph; dims=1)), rand(species(sgraph; dims=2))
                        tpl = (from = s1, to = s2)
                    end
                    push!(newint, tpl)
                else
                    push!(newint, oint)
                end
            end
            for nint in newint
                N[nint.from, nint.to] = true
            end
            # Other interactions
            for spfr in turned_species
                for spto in species(M; dims=2)
                    if rand() < connectance(M)
                        N[spfr, spto] = true
                    end
                end
            end
            # Simplify by removing species from N that have no interactions --
            # we could theoretically wrap this in a loop that would not accept a
            # random network unless it had the exact same number of species, but
            # the results are extremely robust, and so the extra computing time
            # is not worth it.
            simplify!(N)
            # Add the correct values to the storage objects
            OS[sni,sli,ri] = βₜ(βos(M, N))
            WN[sni,sli,ri] = βₜ(βwn(M, N))
        end
    end
end

# We still get ST by element-wise substraction of the 3D arrays, so the
# replicates are preserved -- not that it matters, since we will flatten then by
# the next step...
ST = WN.-OS

# The next lines are about visualisation of these results

p1 = heatmap(proba_rewiring, shared_nodes./r, mean(OS; dims=(3))[:,:,1], aspectratio=1, clim=(0,1))
xaxis!(p1, "Rewiring probability", (0,1))
yaxis!(p1, "Proportion of shared species", (0,1))
title!(p1, "βos")

p2 = heatmap(proba_rewiring, shared_nodes./r, mean(ST; dims=(3))[:,:,1], aspectratio=1, clim=(0,1))
xaxis!(p2, "Rewiring probability", (0,1))
yaxis!(p2, "Proportion of shared species", (0,1))
title!(p2, "βst")

p3 = heatmap(proba_rewiring, shared_nodes./r, mean(WN; dims=(3))[:,:,1], aspectratio=1, clim=(0,1))
xaxis!(p3, "Rewiring probability", (0,1))
yaxis!(p3, "Proportion of shared species", (0,1))
title!(p3, "βwn")

p4 = heatmap(proba_rewiring, shared_nodes./r, mean(ST./WN; dims=(3))[:,:,1], aspectratio=1, clim=(0,1))
xaxis!(p4, "Rewiring probability", (0,1))
yaxis!(p4, "Proportion of shared species", (0,1))
title!(p4, "βst / βwn")

plot(p1, p2, p3, p4, size=(800,800), dpi=600)
savefig("numexp2.png")