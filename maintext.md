---
author: Timothée Poisot
date: \today
title: "Dissimilarity of species interaction networks: how to quantify the impact of species turnover?"
---

Ecological networks are variable both in time and space [@Poisot2015SpeWhy;
@Trojelsgaard2016EcoNet] - this variability motivated the emergence of
methodology to compare ecological networks, in a way that meshes with the usual
approaches of comparison of ecological communities, *i.e.* $\beta$-diversity;
although the definiton of $\beta$-diversity is a contentious topic amongst
community ecologists [see *e.g.* @Tuomisto2010DivBet], the need to understand
network variability is motivated by the fact that species that make up the
networks do not react to their environment in the same way, and therefore the
$\beta$-diversity of networks may behave in complex ways.

@Poisot2012DisSpe and @Canard2014EmpEva have suggested an approach to
$\beta$-diversity for ecological networks which is based on the comparison of
shared and unique links among species, and differentiate this sharing of links
between common and unique species. This framework has been widely adopted since
its publication, with recent examples using it to understand the effect of fire
on pollination systems [@Baronio2021NatFir]; the impact of rewiring on
spatio-temporal network dynamics [@Campos-Moreno2021ImpInt]; the effects of
farming on rural and urban landscapes on species interactions
[@Olsson2021IntPla]; and as a tool to estimate the sampling completeness of
networks [@Souza2021PlaSam]. It has, similarly, received a number of extensions,
including the ability to account for interaction strength [@Magrach2017PlaNet],
the ability to handle probabilistic ecological networks [@Poisot2016StrPro], and
the integration into the Local Contribution to Beta Diversity
[@Legendre2013BetDiv] approach to understand how environment changes drive
network dissimilarity [@Poisot2017HosPar].

In a  recent contribution, @Frund2021DisSpe argues that the calculation of
network dissimilarity terms as outlined by @Poisot2012DisSpe is incorrect, as it
can lead to over-estimating the role of interactions between shared species in a
network ("rewiring"), and therefore underestimate the importance of species
turnover across networks. Here, I present a more thorough justification of the
methodological choices for the @Poisot2021DisSpe method, explain how information
about species turnover can be extracted from its decomposition, and conduct
numerical experiments to guide the interpretation of the $\beta$-diversity
values thus obtained.

## Partitioning network dissimilarity

The approach to quantifying the difference between pairs of networks suggested
by @Poisot2012DisSpe is a simple extension of the overall method by
@Koleff2003MeaBet for species dissimilarity baed on presence-absence data. The
objects to compare, $X_1$ and $X_2$, are partitioned into three sets, $a = |X_1
\cup X_2|$, $b = |X_2 \setminus X_1|$, and $c = |X_1 \setminus X_2|$, where
$|x|$ is the cardinality of set $x$, and $\setminus$ is the set substraction
operation.

In order to apply this framework to networks, we first need to establish a few
definitions. Although ecologists tend to think of networks as their adjacency
matrix, this representation muddies what the terms used in calculation of
dissimilarity ought to be. For this reason, we will use the definition of
networks as a pair of sets, whereing $\mathcal{G} = (V, E)$, where $V$ is a set
containing the vertices (species) of $\mathcal{G}$, and $E$ is a set of ordered
pairs, in which every pair is two elemenets of $V$, which representes the edges
(interactions). In the context of networks comparison (assuming the networks to
compare are $\mathcal{M}$ and $\mathcal{N}$), we can further decompose the
contents of these sets as 

$$\mathcal{M} = (V_c \cup V_m, E_c \cup E_{sm} \cup E_{um}) \,,$$

and 

$$\mathcal{M} = (V_c \cup V_n, E_c \cup E_{sn} \cup E_{un}) \,,$$

where $V_c$ is the set of shared species, $V_k$ are the species belonging only
to network $k$, $E_c$ are the shared edges, and $E_sk$ and $E_uk$ are the
interactions unique to $k$ involving, respectively, only species in $V_c$, and
at least one species from $V_k$.

The metaweb [@Dunne2006NetStr], which is to say the entire regional species pool
and their interaction, can be defined as $\mathcal{M} \cup \mathcal{N}$ (this
operation is commutative), which is to say

$$\mathcal{M} \cup \mathcal{N} = (V_c \cup V_m \cup V_n, E_c \cup E_{sm} \cup E_{um} \cup E_{sn} \cup E_{un}) \,.$$

Making the usual assumption that only species with at least one interaction are
present in the set of vertices, all elements of the set of vertices are used at
least once in the set of edges.

We can similarly define the intersection (similarly commutative) of two
networks:
 
$$\mathcal{M} \cap \mathcal{N} = (V_c, E_c \cup E_{sm} \cup E_{sn})\,,$$

namely all interactions that are established between species that are shared by
both networks.

The decomposition of $\beta$-diversity from @Poisot2012DisSpe uses these
components to measure $\beta_{os}$ (the interaction dissimilarity between shared
species, which @Frund2021DisSpe terms "rewiring"), and $\beta_{wn}$ (the overall
dissimilarity including non-shared species). We can express the components $a$,
$b$, and $c$ of @Koleff2003MeaBet as the cardinality of the following sets:

| Component    | $a$   | $b$                  | $c$                  |
| ------------ | ----- | -------------------- | -------------------- |
| $\beta_{os}$ | $E_c$ | $E_{sn}$             | $E_{sm}$             |
| $\beta_{wn}$ | $E_c$ | $E_{sn} \cup E_{un}$ | $E_{sm} \cup E_{um}$ |

Note that in *all* instances $a+b+c$ is equal to the total number of links in
either $\mathcal{M}\cup\mathcal{N}$ (for $\beta_{wn}$) or
$\mathcal{M}\cap\mathcal{N}$ (for $\beta_{os}$), and that the nodes of the
networks can correctly be induced from the edgelist.

## Quantifying the importance of species turnover

The difference between $\beta_{os}$ and $\beta_{wn}$ stems from the species
dissimilarity between $\mathcal{M}$ and $\mathcal{N}$, and it is easier to
understand the effect of turnover by picking a dissimilarity measure to work as
an exemplar. At this point, @Frund2021DisSpe introduce a substantial mistake in
their manuscript by (i) stating that S\o rensen's and Whittaker's measures of
dissimilarity are the same in the @Koleff2003MeaBet framework (they are not),
and (ii) noting Whittaker's measure as $(b+c)/(2a+b+c)$, which in the
@Koleff2003MeaBet framework is, in fact, $\beta_t$ [@Wilson1984MeaBet]. This
error, although serious, does not change the overall conclusions as these
measures behave in the same monotonously increasing way with regard to
dissimilarity, with a minimum of 0 and a maximum of 1. For the sake of
comparison, I will use $\beta_t$ moving forward, as @Koleff2003MeaBet also
recommend it as a default measure.

Based on a partition between three sets of cardinality $a$, $b$, and $c$, 

$$\beta_t = \frac{b+c}{2a+b+c}\,.$$

So as to simplify the notation of the following section, we will introduce a
series of new variables. Let $A = |E_c|$ be the number of links that are
identical between networks; $S = |E_{sn} \cup E_{sm}|$ be the number of links
that are not shared, but only involve shared species; and $U = |E_{un} \cup
E_{um}|$ the number of links that are not shared, and involve at least one
unique species. Adopting the perspective developed in the previous section,
wherein networks are sets and the measure of $\beta$-diversity operates on this
sets, highlights the conceptual issue in the @Frund2021DisSpe alternative
normalization: they are using components of the networks that are *not* in the
networks being compared.

There are two important points to note here. First, the number or proportion of
species that are shared is not involved in the calculation. Second, the
connectance of either network is not involved in the calculation. That all links
counted in *e.g.* $U$ come from $\mathcal{M}$, or that they are evenly
distributed between $\mathcal{M}$ and $\mathcal{N}$, has no impact on the
result. This is a desirable property of the approach: whatever quantitative
value on the components of dissimilarity can be interpreted in the light of the
connectance and species turnover *without* any risk of circularity. Therefore
the argument of @Frund2021DisSpe, whereby the $\beta_{os}$ component should
decrease with turnover, and be invariant to connectance, does not hold: the
point of the approach is to provide measures that can be interpreted in the
light of connectance and species turnover.

The final component of network dissimilarity in @Poisot2012DisSpe [which is
itself an extension of the initial approach by @Canard2014EmpEva] is
$\beta_{st}$, *i.e.* the part of $\beta_{wn}$ that is not explained by changes
in interactions between shared species ($\beta_{os}$), and therefore stems from
species turnover. This fraction is defined as $\beta_{st} =
\beta_{wn}-\beta_{os}$.

The expression of $\beta_{st}$ does not involve a partition into sets that can
be plugged into the framework of @Koleff2003MeaBet, because the part of
$\mathcal{M}$ and $\mathcal{N}$ that are composed of their unique species
cannot, by definition, share interactions. One could, theoretically, express
these as $\mathcal{M} \setminus \mathcal{N} = (V_m, E_{um})$ and $\mathcal{N}
\setminus \mathcal{M} = (V_v, E_{un})$ (note the non-commutativity here), but
the dissimilarity between these networks is trivially maximal for the measure
considered.

Using the $\beta_w$ measure of dissimilarity, we can write

$$\beta_{os} = \frac{S}{2A+S}\,,$$

and

$$\beta_{wn} = \frac{S+U}{2A+S+U}\,.$$

As $\beta_{os}$ has the form $x/y$, and $\beta_{wn}$ has the form $(x+k)/(y+k)$,
with $x = S$, $y = 2A+S$, and $k = U$, as long as $k \ge 0$, it is guaranteed
that $\beta_{wn} \ge \beta_{os}$, and therefore that $0 \ge \beta_{st} \ge 1$.

Although it holds little information, we can get an expression for $\beta_{st}$,
by bringing $\beta_{os}$ and $\beta_{wn}$ to a common denominator:

$$\beta_{st} = \frac{2AU}{(2A+S)(2A+S+U)}\,.$$

Note that this measure varies in a non-monotonic way with regards to the number
of shared interactions -- this is obvious when developing the denominator into

$$4A^2 + 2AS + 2AU + 2AS + S^2 + SU = 4A^2 + S^2 + 4AS + 2AU + SU\,,$$

which reveals why this value tends towards 0 whenever one component
($\beta_{wn}$ or $\beta_{os}$) explains the entire dissimilarity (*i.e.* its
value tends towards 1). Note further that for this reason, @Poisot2012DisSpe
suggest that $\beta_{st}/\beta{wn}$ is a better indicator of the *relative*
importance of turnover processes on network dissimilarity. This can be
calculated as 

$$\frac{\beta_{st}}{\beta_{wn}} = \frac{2AU}{(2A+S)(2A+S+U)}\times\frac{S+U}{2A+S+U}\,,$$

which reduces to 

$$\frac{\beta_{st}}{\beta_{wn}} = \frac{2AU}{(2A+S)(S+U)}\,.$$

This expression admits roots for $A=0$ (the turnover of species has no
contribution to the difference between $\beta_{wn}$ and $\beta_{os}$ if there
are no shared species), and  for $U = 0$ (the turnover of species has no
contribution if all species are shared).

## Is this decomposition over-estimating the effect of "rewiring"?

One of the arguments put forth by @Frund2021DisSpe is that the decomposition
outlined above will overestimate the effect of rewiring; I argue that this is
based on a misunderstanding of what $\beta_{st}$ achieves. It is paramount to
clarify that $\beta_{st}$ is not a direct measure of the importance of turnover:
it is a quantification of the relative impact of rewiring to overall
dissimilarity, which, all non-turnover mechanisms being accounted for in the
decomposition, can be explained by turnover mechanisms.

We can re-calculate the illustration of @Frund2021DisSpe, wherein a pair of
networks with two shared interactions ($A = 2$) receive either an interaction in
$S$, in $U$, or in both:

| $A$ | $S$ | $U$ | $\beta_{os}$ | $\beta_{wn}$ | $\beta_{st}$ | $\beta_{st}/\beta_{wn}$ |
| --- | --- | --- | ------------ | ------------ | ------------ | ----------------------- |
| 2   | 0   | 0   | 0            | 0            | 0            |                         |
| 2   | 1   | 0   | $1/5$        | $1/5$        | 0            | 0                       |
| 2   | 0   | 1   | 0            | $1/5$        | $1/5$        | 0                       |
| 2   | 1   | 1   | $1/5$        | $1/3$        | $2/15$       | $2/5$                   |

The over-estimation argument hinges on the fact that $\beta_{st} < \beta_{os}$
in the last situation (one interaction as rewiring, one as turnover). Reaching
the conclusion of an overestimation from this is based on a mis-interpretation
of what $\beta_{st}$ means. The correct interpretation is that, out of the
entire network dissimilarity, only three-fifths are explained by re-wiring. The
fact that this fraction is not exactly one-half comes from the fact that the
@Wilson1984MeaBet measure counts shared interactions *twice* (*i.e.* it has a
$2A$ term), which over-amplifies the effect of shared interactions as the
network is really small. Running the same calculations with $A = 10$ gives a
relative importance of the turnover processes of 47%, and $\beta_{st}$ goes to
$1/2$ as $A/(S+U)$ increases. As an additional caveat, the value of $\beta_{st}$ will depend on the measure of beta-diversity used. Measures that do not count the shared interaction twice are not going to amplify the effect of rewiring.

## Does the partition of network dissimilarity needs a new normalization?

Based on the arguments presented above, I do not think the suggestion of
@Frund2021DisSpe to change the denominator of $\beta_{os}$ makes sense as a
default; the strength of the original approach by @Poisot2012DisSpe is indeed
that the effect of turnover is based on a rigorous definition of networks as
graphs (as opposed to networks as matrices), in which the induction of vertices
from the edgelist being compared gives rise to biologically meaningful
denominators. The advantage of this approach is that at no time does the
turnover of species itself, or the connectance of the network, enter into the
calculation. As such, it is possible to use $\beta_{os}$ and $\beta_{wn}$ in
relationship to these terms, calculated externally [as was recently done by
*e.g.* @Higino2021BetPhy] without creating circularities.

The choice of changing the denominator hinges on what one admits as a definition
for $\beta_{st}$. If the point of $\beta_{st}$ is to be a component of overall
$\beta$-diversity as advocated by @Frund2021DisSpe and @Novotny2009BetDiv, a
change of numerator *might* be acceptable. Nevertheless, this change of
numerator contributes to blurring the frontier between a measure of interaction
dissimilarity and a measure of community dissimilarity, and may warrant a full
methodological assessment. Conversely, if as we argue in @Poisot2012DisSpe,
$\beta_{st}$ is to be meant as a *guide* to the interpretation of $\beta_{wn}$
and $\beta_{os}$, and related to actual measures of species turnover and network
connectance, one must not change the denominator. It is central to recognize
that there are multiple reasons to calculate network dissimilarity, and it is
our opinion that the arguments levied by @Frund2021DisSpe against the original
partition stem from a misunderstanding of what it intends to do (and does,
indeed, do well), not from intrinsic methodological issues in the partition
itself.

## References