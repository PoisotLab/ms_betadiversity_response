using EcologicalNetworks:length
using EcologicalNetworks
using StatsPlots
using PlotThemes
using Random
using StatsBase

Random.seed!(90210)
theme(:mute)
default(frame=:box)

βₜ = EcologicalNetworks.KGL08

# Simulation 1 - L1 = L2, we vary the number of shared, unique, and non-unique links

L1, L2 = 200, 200

shared_links = LinRange(0.0, 1.0, 200)
links_rewired = LinRange(0.0, 1.0, 220)

OS = zeros(Float64, length(shared_links), length(links_rewired))
ST = similar(OS)
WN = similar(OS)

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

ST = WN .- OS

contour(links_rewired, shared_links, OS, clim=(0, 1), aspectratio=1, fill=true, lc=:white, levels=5)
xaxis!("Proportion of rewired links", (0, 1))
yaxis!("Proportion of shared links", (0, 1))
title!("βos")

contour(links_rewired, shared_links, ST, clim=(0, 1), aspectratio=1, fill=true, lc=:white, levels=5)
xaxis!("Proportion of rewired links", (0, 1))
yaxis!("Proportion of shared links", (0, 1))
title!("βst")

contour(links_rewired, shared_links, WN, clim=(0, 1), aspectratio=1, fill=true, lc=:white, levels=5)
xaxis!("Proportion of rewired links", (0, 1))
yaxis!("Proportion of shared links", (0, 1))
title!("βwn")

contour(links_rewired, shared_links, ST ./ WN, clim=(0, 1), aspectratio=1, fill=true, lc=:white, levels=5)
xaxis!("Proportion of rewired links", (0, 1))
yaxis!("Proportion of shared links", (0, 1))
title!("βst / βwn")

# Simulation 2 - simulate networks and measure the turnover in species
M = simplify(rand(BipartiteProbabilisticNetwork(rand(Float64, 25, 25) .* 0.3)))
r = richness(M; dims=1)

shared_nodes = 0:r
proba_rewiring = LinRange(0.0, 1.0, 20)
replicates = 12

OS = zeros(Float64, length(shared_nodes), length(proba_rewiring), replicates)
WN = similar(OS)

for (sni, sn) in enumerate(shared_nodes)
    for (sli, sl) in enumerate(proba_rewiring)
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
            # Simplify
            simplify!(N)
            OS[sni,sli,ri] = βₜ(βos(M, N))
            WN[sni,sli,ri] = βₜ(βwn(M, N))
        end
    end
end

ST = WN.-OS

heatmap(shared_links, shared_nodes./r, mean(WN; dims=(3))[:,:,1], aspectratio=1)
xaxis!("Proportion of rewired links", (0,1))
yaxis!("Proportion of shared species", (0,1))

heatmap(shared_links, shared_nodes./r, mean(OS; dims=(3))[:,:,1], aspectratio=1)
xaxis!("Proportion of rewired links", (0,1))
yaxis!("Proportion of shared species", (0,1))

heatmap(shared_links, shared_nodes./r, mean(ST; dims=(3))[:,:,1], aspectratio=1)
xaxis!("Proportion of rewired links", (0,1))
yaxis!("Proportion of shared species", (0,1))

heatmap(shared_links, shared_nodes./r, mean(ST./WN; dims=(3))[:,:,1], aspectratio=1)
xaxis!("Proportion of rewired links", (0,1))
yaxis!("Proportion of shared species", (0,1))