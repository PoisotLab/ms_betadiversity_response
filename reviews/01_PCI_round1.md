# Editor's comments

The most salient point is that the preprint is a kind of response/comment to a
recent paper of Fründ in Ecosphere.

> I have de-emphasized this component of the manuscript in the revision - it is
> now framed as a "best practices" recommendations, with a discussion of the
> consequences of changing the denominator as suggested by Fründ towards the
> end. This is made clearer in the introduction, by explaining first the
> difficulties in interpreting the components, then the importance of choices
> about the denominator.

In addition, the ms does not follow a classical
Introduction/M&M/Results/Discussion section, arguably because of the
response/comment nature of this contribution.

> I do not think the "classical" structure would be appropriate, as it would
> essentially disconnect the explanations of how components are calculated from
> their use in the numerical experiments. To guide the readers I have added an
> expanded outline to the introduction, and strengthened the rationale for all
> of the numerical experiments in the text.

As suggested by reviewer 2, the ms would gain clarity if,

- an Introduction section more explicitly states the nature of the diverging
views between Fründ and you, and thus more clearly exposing the motivation of
challenging the recent paper of Fründ,

> As mentionned in a previous comment, I have moved this part of the manuscript
> towards the end, with more explicit recommendations about the use of the
> proposed alternative.

- a Discussion section synthesizes the pros and cons of both approaches, as it
seems that each method can be justified and used in an appropriate context.

> Same answer as above. In addition, I have clarified the mathematical notation
> throughout, so that it will hopefully be easier to follow.

# Reviewer 1

## Main comments

While there is no doubt that partitioning interaction turnover is a needed tool
for ecologists, there has been recent debate on how to best perform this
partitioning. I assess this manuscript from the perspective of the interested
user in applying such indexes and doing it correctly. This means that while I
mostly follow the maths decomposing the different indexes, I focus my review on
their interpretation, and I can't fully assess the more complex mathematical
derivations.

> I have attempted to clarify the notation *and* the way the derivations are
> presented, so that they are easier to follow by interested users as well.

I appreciate the detailed explanation provided and the rationale behind each
index. I think the calculation and interpretation of BetaOS are clear. This one
can be interpreted unambiguously as interaction rewiring among shared species.
To me, this is the key index to interpret for most ecological questions and can
be interpreted as a single index (probably along with the proportion of shared
species for context).

Next, Beta WN is also intuitive, but as depicted in Fig 2 it needs to be clearly
stated that it depends on both interaction rewiring and the proportion of shared
species. For most ecological questions it may be of secondary importance.

The most problematic term for me is BetaST. This is also stated in the
manuscript, and I agree that it has caused a larger degree of confusion.
Interpreting this term (beyond an error term, which is the simplest
interpretation) is complex. In fact, the manuscript often describes their
behaviour as "as expected" given its mathematical formulation, but looking at
the literature it is clear that most researchers (including myself) were not
expecting some of those behaviours. However, even when completely understanding
how it behaves, the interpretation is still too complex for me in order to be
useful. Fig 2 B is ilustrative to me. Note that at q = 1; BetaST = 0 regardless
of the proportion of species sharing! Hence, BetaST tells you nothing about the
contribution of species turnover when rewiring is high (ST stands for species
turnover, so it's normal people get confused). I had this discussion with
several researchers, and I can tell this is hard to grasp at first. Following
with figure 2B, you can see that the same value can mean two very different
things. BetaST can be low if there is high rewiring, or if rewiring is very low,
but they share most species. This is less accentuated in the relative importance
but is still the case. Hence, I would suggest giving clear recommendations on
not interpreting BetaST as a primary index, but only in cases when is pertinent,
and giving a clear context (i.e. the proportion of shared species, which I think
is much more useful for ecological questions). I know this is suggested in some
places of the manuscript, but in my opinion, it can be said stronger. In fact, I
would love to see a section on how to interpret each component, recommending the
interpretation of BetaOS as the most straightforward, and cautioning that
interpretation of beta WN, and especially ST need to be done in context and
can't be interpreted alone.

## Minor comments

line 13: Currently it is done in several ways, so maybe better cite here Poisot
et al 2012.

line 15-21: This is a very long sentence with "i.e", ";", "-", which I had to
read twice. Consider splitting it up. At least a point after "Tuomisto 2010)" I
think is needed.

lines 72-74: letters used in the text do not match those in the equation. (m =
k, I think)

line 85: The second "similarly" can be replaced by "also".

line 98: I would name between brackets the abbreviations of Sorensen and
Wittaker indexes when first mentioned a couple of lines above.

line 111-112: I think the reader needs first to be introduced briefly to what
Frund did differently to understand this sentence. This is not done until line
160.

line 150: I understand we do not need species yet, but this first numerical
exercise is quite abstract. Maybe guide the reader on the purpose of the
experiment. Which I believe is to describe the behaviour of the components.

line 165: I agree with this interpretation "it is a quantification of the
relative impact of rewiring to overall dissimilarity". My critique is that this
measure, as defined here, is hard to interpret ecologically and I may dare to
say that irrelevant for most ecological questions.

Line 165: But I do not follow the next sentence. What do you mean by "all
non-turnover mechanisms being accounted for in the decomposition, can be
explained by turnover mechanisms". See my point above on the ambiguity of
interpreting its values, which can emerge from very different ecological
situations.

line 171: Adding the illustration may help a lot. I had to write it myself (and
it was not easy for me).

line 228: formatting error when designing the title.

Fig 3: Maybe expressing the x-axes as connectance dissimilarity would help the
reader? I would also appreciate a comment in the text on further partitioning of
the components on changes due to differences in link number (related to
connectance) and true link turnover. Those are discussed by Frund, and I think
it can help interpret those values depending on your ecological question.

line 260: Here a caution on interpreting BetaST may be good as stated above.

# Reviewer 2

## Main comments

The paper is dedicated to an outstanding question, how to compare interaction
networks teasing apart the effect of species turnover and interaction rewiring.
This question is central in many studies, due to the recent development of
efficient data acquisition techniques in the "multiple networks era". 

The author did not write a cover letter. It would have been important in
particular regarding the following observations. From the last paragraph of the
first section, we early understand that this manuscript is actually responding
another, Fründ 2021, that gave some objective criticisms about the renowed
method developed by the author of the present manuscript in Poisot et al 2011.
In my opinion, the abstract is not clear about this fact because it does not
clearly state the debate which is central in the paper, i.e. the different views
proposed by Fründ and Poisot respectively. Indeed, having read the paper by
Fründ seems to be a prerequisite to read this manuscript. The paper does not
have a classical structure with the usual sections (e.g. no
"Introduction/Discussion" keywords) but is linear in the sense that the
arguments are displayed one after the other. Again, this is not clear enough to
easily see 1/ what is exactly the debate 2/ what is the author's strategy to
convince the reader in this debate.  Finally, the notations adopted in this
manuscript are completely different from those of Fründ, which does not help the
reader in evaluating the different conclusions raised by both authors.
 
About the numerical experiments, they are well conducted and clear but could be
better motivated.  The first one consists in redoing the small example proposed
in Fründ (here, this is clear). The conclusion is that actually both authors
have two different understanding of the betast but in fact the actual
incoherence between the two authors comes from the fact that Poisot counts the
shared interactions twice, which is actually not an issue as long as the network
is large enough. The next two experiments are dedicated to bipartite networks :
what is the reason for this choice? In fact, one can question the
generalisability of the conclusions to unipartite networks. The first one
explores the variation of the different beta measures with rewiring and
turnover, and actually shows these measures are responding as expected. The
second one  is investigating the possible link between the beta measures and
connectance. In fact, it is clear that connectance difference can induce a large
dissimilarity between two networks (again, not only bipartite networks). The
author shows his results for two connectance values in the bipartite case, but
one can expect a broader exploration of the link between intrinsic networks
properties (for a range of variation, not two values) and the beta measures.
Connectance is indeed important, but degree distribution is paramount as well
(at a fixed connectance, for instance). Finally the author gives his opinion
about the need for a new denominator as suggested by Fründ. However, the
previous experiments were not performed with an alternative denominator, and
then the reader can not compare different proposals. 
 
The paper is, again, dealing with an outstanding question using very important
methodology. But the reader could expect more systematic experiments (unipartite
case?) and a better organized manuscript.

## Minor comments

L.32 Maybe citing Ohlman et al, Eco. Lett. as another extension.

> Done.

L.53 Typo "baed"

> Changed.

L.54 Using "x" here and then using "x" L.56-57 could be misleading. Maybe just
say "|.|" is the cardinality operator ?

> Done.

L.66-68 This assumes that edges are directed. This can not be the case, for
mutualistic interaction (btw, replacing an undirected edge by two directed edges
is still possible)

L.70 Well, classically, Vm are the vertices of m (all of them)...

L.70 Is Ec with "c" for "in common"? If yes, tell it because it can help
memorizing.

L.71 Typo: capital N instead of capital M.

L.87 In my opinion, that's a good point to call it "rewiring".

L.89 a is supposed to be the union, and Ec is the intersection. Is there a
problem here in the analogy? Maybe I missed something, my apologies in advance.
Also, a/b/c are supposed to be cardinality, and this is not the case in the
table. Am I right?

L.98 betaSor and betaw, betat not defined before.

L.98 Typo : where is (i) ?

L.93-102 This paragraph is difficult to understand without being an expert of
the mentioned papers. What is the message here?

L.104 OK, here betat is defined. This is related to Dice (1945) or Sorensen
(1948)?

L.111-112 "they are using components of the networks that are not part of the
networks being compared" is not clear.

L.114 Which calculation? 

L.113-121 This part of the paper looks like an answer to Fründ. Is this an
"answer" paper ?

L.120 Typo : "very point"

L.130-143 Straightforward mathematics, details could go in appendix.

L.151 Once again, the choice of these capital letters for shared (A), rewired
(S), and unique (U) does not help.

L.151 In fact, ps is a proportion but pr is not, whereas pr X (1-ps) is a
proportion... (pr is a proportion of the proportion of non-shared links).

L.153 How can "shared links [be] rewired" if they are shared?

L.164-166 Maybe split this 3 lines sentence, to facilitate understanding ?

L.176 Precise what are these three-fifth (this is because betast/betawn=2/5)

L.178 This "2A" is the key point. This is enough to understand the full paragraph.

L.182 Then, should we move to such a measure that do not amplify the effect of
rewiring?

L.184 Why do we switch to bipartite networks here ?

L.184 Again, the choice of capital R for a number of species could be discussed.

L.196-199 Are these conclusions not obvious ? From the definitions of the beta indices.

L.197 Previously, we have "proportions" and now we have "probabilities". It
would be better to harmonize.

Figure2 Please harmonize the axis labels with Figure 1.

L.208-209 Typo: bad copy/paste here. 

L.210 Again, why bipartite networks? Also, why varying connectance? It is a
relevant idea but one can make the degree distribution vary as well. This
distribution can have huge impacts on beta indices. 

Also, it is normal to expect huge dissimilarity between network of different
connectance because the number of un-shared links is expected to be high. This
seems straighforward. Here, the reader can expect a formulation of the
hypothesis that will be tested : what is this numerical experiment for?

L.243 The formulation "I do not think" is not appropriate because the reader
will not expect an opinion but conclusions drawn for the previous experiments.

L.244 It would be necessary to recall the proposal of Fründ somewhere in the
paper. Did the author perform the same numerical experiments with the other
numerator?

L.246 "rigorous definition of networks as graphs (as opposed to networks as
matrices)". It is unclear why we have to oppose the matrix and the graph view,
since graphs and adjacency matrices belong to the same conceptual context. 

L.252-254 In fact, one can think that the debate is solved by choising a
definition. A matter of perspective, at some point.