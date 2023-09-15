# Network Motifs: Simple Building Blocks of Complex Networks

https://www.cs.cornell.edu/courses/cs6241/2020sp/readings/Milo-2002-motifs.pdf

\*\*Note: Currently, noone is signed up for the hacker role. I may decide to implement this simple idea myself. - https://github.com/aplbrain/dotmotif - https://github.com/aplbrain/grandiso-networkx - https://github.com/p-koo/learning_sequence_motifs - https://github.com/gtremper/Network-Motif

Check this out: https://networkrepository.com/graphvis.php

"To go beyond these global features would require an understanding of the basic structural elements particular to each class of networks (9).
To do this, we developed an algorithm for detecting network motifs: recurring, significant patterns of interconnections."

Key definition for network motifs: Recurring, significant patterns of interconnections.

"Here we generalize this approach to virtually any type of connectivity graph and find the striking appearance of motifs in networks representing a broad range of natural phenomena." - Is it really "striking"? Take a poll of the class beforehand - which of the 3-node motifs to we believe will be most common in e.g. the transcription network.

"Patterns that are functionally important but not statistically significant could exist, which would be missed by our approach."

"The two transcription networks show the same motifs: a three-node motif termed “feedforward loop” (11) and a four-node motif termed “bi-fan."

## Concurrent paper by same authors

[Network motifs in the transcriptional regulation network of Escherichia coli](https://www.nature.com/articles/ng881z.pdf)

- "The feedforward loop motif often occurs where an external
  signal causes a rapid response of many systems (such repression
  of sugar utilization systems in response to glucose, shift to anaerobic metabolism). The abundance of coherent feedforward
  loops, as opposed to incoherent ones, suggests a functional
  design (Table 1)."

## Questions

Why did they choose to focus on common connectivity patterns within 3/4-node subgroups?
Why do they think this measure is representative of "recurring, significant patterns of interconnections"?
What would other measures be?

- Increasing n in n-node subgroups.
- what about n-edge subgroups?
- use direction in selecting subgroups? E.g. to take regard of causal neighborhoods, like looking at connectivity patterns within markov blankets?

Are 3-node subgraphs expressive enough that finding a statistically highly occuring pattern within it is nontrivial/meaningful/surprising? E.g. they find the "feedforward loop" motif in transcription networks. This is significant if we consider the frequency of this connection pattern amonst _all possible_ ones. But looking at Fig. 1B, we can also see that the feedforward loop (#5) is also only one of the two fully connected graphs without bidirectional edges, and the other one is a loop. So, perhaps we could arrive at the same prediction (feedforward loop is most common 3-node motif) by hypothesising, say, two tendencies of natural networks:

1. Favor dense over spare connections (reasoning: nodes in natural networks occur in isolation less often than with neighbors)
2. Avoid tight loops (reasoning: loops are bad at propogating information, and the recursion muddies/saturates signals and obfuscates signal source)

These two principles would also lead us to #5 (feedforward loop).
So, I don't think we should reasonably assign a uniform prior across all possible connections, and find it slightly misleading to frame the findings with surprise at how non-uniform the actual distribution over motifs is.

Also, it could be that e.g. many 3-node network motifs are present in a >3-node motif that isn't studied here. If this were the case, positing explanations like the following may not be appropriate:

"The feedforward loop motif common to both types of networks may play a functional role in information processing. One possible function of this
circuit is to activate output only if the input signal is persistent and to allow a rapid deactivation when the input goes off (11)."

For example, when I see the 3-node network in the context of neurons, I think of skip-connections, which both have a somewhat different purpose (removing the burden of intermediate nodes adding the identity of their parent signals to their child nodes), and are better understood when we consider the intermediate layers that are being skipped.

- Actually, this reminds me of papers in convolution architecture - in particular in parameterized architectures like the alpha in [MobileNets](https://arxiv.org/pdf/1704.04861.pdf) / [MobileNetV2](https://arxiv.org/pdf/1801.04381.pdf). These provide, in my opinion, a more promising direction for a framework of thinking about common network connectivity patterns - particularly, efficiency vs. accuracy tradeoff with simultaneous representational capacity goals and energy resource constraints.

Put more generally, I think we need to look at higher-order properties of networks than frequency of specific N-node connectivity patterns, and I don't think we can/should avoid considering the functional properties of the network.

"This suggests that motifs can define broad classes of networks, each with specific types of elementary structures." - The suggestion with this wording, and throughout the paper, seems to be that these diferrent systems may be combining these "motifs" together in various waves, to acheive their goals, whereas I don't see why we should conclude they are fundamental building blocks rather than common byproducts of other mechanisms. I find the milder wording they use later in the same paragraph more appropriate: "Information processing seems to give rise to significantly different structures than does energy flow."

Again, in the last paragraph, the authors interpret network motifs in two different ways - as "elementary computational circuits" (which I think is unjustified and likely incorrect), and another milder interpretation as "structures that arise because of the special constrains under which the network has evolved".

**I can crisply demonstrate my issue here by finding two or more networks with vastly different properties/goals/contexts, with intuitively obvious structure when wee look at larger neighborhoods than 3/4, but with similar distributions of 3/4 node motifs.**
The idea is to challenge the general conclusion/interpretation of the paper that if we find that graphs in different domains share similar motif distributions, this may reflect similar underlying organizational principles. "... motifs can define broad classes of networks, each with specific types of elementary structures."

Ideas for networks to investigate:

- Ladders
- CNNs
- Spider web

## Need further understanding

"Edges are directed from a gene that encodes for a transcription factor protein to a gene transcriptionally regulated by that transcription factor."

Need to fully understand "gene transcriptionally regulated by that transcription factor" - I think I get it but be sure.

## Potential future research

A natural practical question is how graph neural networks perform with these motifs imposed as priors over the node connectivity.

## Presentation

Call back to lecture - triangles are important: friends of friends are my friends. This is the feedforward loop. Also, same lecture we talked about clustering coefficients, and how rand. networks have smaller cluster coeff. This is related. ALSO in that same lecture, we saw Erdos/Renyi random graphs.

## Slides

### Definition

Network Motifs: Patterns of interconnections occurring in complex networks at numbers that are significantly higher than those in randomized networks.

### Goals/approach:

Rather than looking at global statistical features, e.g.

- "Small world" property
- Scale-free networks

The authors aim to understand the "basic structural elements" of network classes.

Here we generalize this approach to virtually any type of connectivity graph and find the striking appearance of motifs in networks representing a broad range of natural phenomena.

### Methodology

- Scan directed graph for all possible n-node subgraphs (n = 3 and 4), and count the occurrences.
- Compare the real network to randomized networks.
- Find patterns appearing in the real network significantly more than in randomized network.

That's it!

### Methodology: Randomized networks

- Each node in the randomized network has the same number of incoming and outgoing edges as the corresponding node in the real network.
  - Accounts for patterns that appear only because of _single-node_ characteristics of the network.
    - Example (with picture): presence of nodes with many edges.
- Preserve the same number of appearances of all (n - 1)-node subgraphs as in the real network. (Only for n = 4 motif detection, not n = 3.)
  - Ensures a high significance is not assigned because the pattern has a highly significant subpattern.
    - Show an example.

"Network motifs" are patterns with probability P of appearing in a randomized network an equal or greater number of times than in the real network is lower than P = 0.01.

Limitation:

> Patterns that are functionally important but not statistically significant could exist, which would be missed by our approach.

### Methodology: 3-node motifs

Show all 3-node motifes.

### Methodology: 4-node motifs

Show important 4-node motifs.

### Findings

Distinct motifs in different classes of networks

- Food webs, genetic networks, internet

"Information processing" networks share similar motifs

- Biomolecules within a cell
- Synaptic connections between neurons in C. elegans

### Claims/interpretations/hypotheses

- **Motifs may define universal classes of networks.**
- **Motifs may uncover the basic building blocks of most networks.**
