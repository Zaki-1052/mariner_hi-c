1
2
3
4
5
6
7
8
9
10
11
12
The histone deubiquitinase Bap1 maintains the distinction between euchromatin and
heterochromatin during neurodevelopment
Kelly C Wang1#, Challana E. Tea1#, Yixing Du2, Addison C. Bullard1, Karen Tsai1, Donald Le1, Aditya
Parmar1, Zhendong Song1, Jonathan Tawil1, Amanda Roberts3, Kimberly Dore2, Cole J. Ferguson1*
1Department of Pathology, University of California San Diego, La Jolla, CA 92093, USA
2Department of Neurosciences, University of California San Diego, La Jolla, CA 92093, USA
3The Scripps Research Institute, La Jolla, CA 92093, USA
#co-first authors
*Corresponding author: cferguson@health.ucsd.edu
1
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
Abstract
Neurodevelopment unfolds through a complex sequence of genes being expressed and repressed,
much of which is controlled at the level of histone modifications that demarcate euchromatic
(transcriptionally active) from heterochromatic (inactive) loci. Complicating this picture, however, is our
recent work showing that, in neurons, ubiquitination of H2AK119 (H2Aub), which is carried out by
Polycomb repressive complex 1 (PRC1) and generally considered repressive, appears within both
euchromatin and heterochromatin in cerebellar neurons. To better understand the effects of H2Aub in
these compartments, we removed the H2Aub deubiquitinase Bap1 from excitatory granule neurons of the
mouse cerebellum. Chromatin mapping at early and adult timepoints revealed enrichment of H2Aub at
euchromatic loci, accompanied by loss of activating H3K27ac at enhancers, shifting these regions toward
a heterochromatic state. To our surprise, however, H2Aub was depleted from heterochromatin, which
along with other histone modification effectively euchromatized these regions. The direction of
transcriptional change at affected genes was consistent with the particular combination of perturbed
histone modifications within each locus, leading to dysregulation of genes required for neuronal growth,
synapse formation and function, cell adhesion pathways and lineage-defining transcription factors—
paralleling Bap1 mutations that increase metastatic risk in certain malignancies. To explore the effects of
gene dysregulation on behavior, we generated mice lacking Bap1 in the forebrain. These mice were
developmentally delayed and exhibited profound behavioral abnormalities along with hyperexcitability
and seizures. Our data reveal that histone ubiquitination is essential to maintaining distinct chromatin
states and lend insight into how Bap1 promotes proper neurodevelopment.
Teaser: Bap1 loss weakens the distinction between heterochromatin and euchromatin in mammalian
neurons.
2
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52
53
54
55
56
57
58
59
60
61
62
63
64
65
66
67
68
69
70
INTRODUCTION
Neurodevelopment proceeds through an ordered series of transitions between distinct chromatin
states. These states are specified in large part by post-translational modifications of histone proteins,
hundreds of which have been catalogued across the genome (1, 2). Among the most abundant and broadly
distributed modifications, H2AK119ub (Histone 2A ubiquitinated at Lys119, hereafter H2Aub) (3), is
deposited by non-canonical (also known as ‘variant’) Polycomb-Repressive Complex 1 (PRC1) (4).
H2Aub recruits Polycomb-Repressive Complex 2 (PRC2) to deposit H3K27me3 (5–7), the epigenetic
hallmark of facultative heterochromatin, thereby condensing chromatin to make it inaccessible to
transcription factors. At least this is the traditional view: recent work in mammalian neurons has
complicated this picture by showing that H2Aub often appears in the absence of H3K27me3 (8). In fact,
over the course of cerebellar neurodevelopment H2Aub becomes depleted from H3K27me3-positive
heterochromatin while becoming enriched within H3K27ac-positive active enhancers and promoters (8).
Studies in nematodes have also shown that H2Aub often appears outside of facultative heterochromatin
within enhancers that activate neuronal genes (9).
If the nature of H2Aub's activities remains uncertain, it is nonetheless clear that its proper
deposition and removal are essential to neurodevelopment. Mutations in genes encoding either of PRC1's
E3 ligases that deposit H2Aub (RING1 or RNF2) cause severe neurodevelopmental syndromes (10–12).
Removal of H2Aub is carried out by the histone deubiquitinase BAP1 (BRCA1-associated Protein 1), and
global Bap1 knockout mice die in utero (13). Because inactivating somatic mutations in BAP1 drive uveal
melanoma (14, 15), mesothelioma (16), and other malignancies (13), BAP1 has been studied primarily as
a tumor suppressor gene (17, 18). Indeed, germline loss-of-function BAP1 mutations cause a tumor
diathesis syndrome (19). However, a distinct set of germline BAP1 mutations underlie Küry–Isidor
Syndrome (KIS), which causes varying degrees of developmental delay, autistic features, and/or seizures
(20, 21). (These 11 documented individuals were young at ascertainment, so it is not clear whether they
will eventually prove to be cancer-prone.) Mutations in ASXL1 and ASXL3, genes that encode binding
partners in BAP1-containing protein complexes (22), also cause syndromic neurodevelopmental disorders
(Bohring-Opitz (23) and Bainbridge-Ropers (24), respectively).
To better understand the role of Bap1 in neurodevelopment, we generated two conditional mouse
mutants. To examine chromatin changes, we deleted Bap1 in the cerebellum, whose cellular homogeneity
(being mostly granule neurons) lends itself to chromatin mapping. We observed unexpectedly different
effects on euchromatic and heterochromatic loci, particularly affecting the expression of genes involved
in neuronal adhesion, calcium influx, synapse function and transcriptional regulation. To test the effect of
gene dysregulation on cortical function, we deleted Bap1 mainly in excitatory neurons of the forebrain.
The forebrain contains excitatory and inhibitory neurons alongside abundant glia; although this
3
71
72
73
74
75
76
77
78
79
80
81
82
83
84
85
86
87
88
89
90
91
92
93
94
95
96
97
98
99
100
101
102
103
104
heterogeneity makes chromatin mapping of Bap1 function infeasible, it allowed us to characterize
behavioral and neurophysiological alterations resulting from Bap1 loss. Our findings demonstrate that
Bap1 is essential for maintaining both euchromatin and heterochromatin identity to direct proper
neurodevelopment.
RESULTS
Loss of neuronal Bap1 increases H2Aub within euchromatin and depletes H2Aub from
heterochromatin
Because global Bap1 knockout is embyronic-lethal (13), we created a conditional deletion by
crossing Bap1f/f mice (13) with mice expressing cre recombinase under the control of the Math1/Atoh1
promoter, which is active in most granule neuron progenitors (25). Cerebellar granule neurons are highly
abundant excitatory glutamatergic neurons: they constitute >85% of cerebellar cells and >99% of its
neurons (26), imparting unmatched neuronal homogeneity as bulk input for in vivo chromatin mapping
and gene expression studies.
Before examining chromatin, we asked whether deletion of Bap1 from granule neurons produces a
deficit in cerebellar function. We evaluated Bap1f/f,Math1-cre mice at 4-7 months of age on the
accelerating rotarod. The mutant mice showed a shorter latency to fall, with no significant improvement
within daily trials (Figure 1A) or over consecutive days of testing (Figure S1A). Bap1 expression in
granule neurons is therefore required for cerebellar function, i.e., the chromatin changes described below
occur in a neuronal population whose physiology requires this deubiquitinase.
To examine how the loss of Bap1 affects the redistribution of H2Aub over the course of cerebellar
development, we chose postnatal day 12 (P12) as the ‘early’ timepoint in cerebellar development to
maximize the representation of immature granule neurons, which are produced in huge numbers from
postnatal days 4-16 (P4 to P16) in mice (8, 26, 27). Age 3 months served as the ‘late’ timepoint, when the
cerebellum achieves epigenetic, transcriptomic and functional maturity (8, 28). RNA-seq analysis of
nuclei isolated from the cerebellum of control mice confirmed that Bap1 and other components of the PR-
DUB complex were expressed at these timepoints (Figure S1B). Immunoblot showed reduction of Bap1
in the cerebellum of Bap1f/f,Math1-cre mice at both ages (Figure S1C). Interestingly, the total abundance
of Bap1’s substrate H2Aub was unchanged in these same lysates at both ages, indicting that Bap1 loss
does not alter the global level of H2Aub in neurons. This result also establishes that the percentile-based
normalization (8) we apply to chromatin mapping data below is not confounded by a genotype-dependent
change in overall H2Aub.
To examine the genomic distribution of H2Aub, we isolated cerebellar nuclei and mapped H2Aub
using CUT&RUN (Cleavage Under Targets and Release Using Nuclease), an antibody-based next-
4
105
106
107
108
109
110
111
112
113
114
115
116
117
118
119
120
121
122
123
124
125
126
127
128
129
130
131
132
133
134
135
136
137
138
generation sequencing method known for high resolution and low background (8, 29, 30). Principal
component analysis (PCA) of control and Bap1f/f,Math1-cre H2Aub CUT&RUN data at either the early
or late timepoint segregated biological replicates along PC1, with genotype accounting for most of the
observed variability (Figure S1D). To identify loci where H2Aub differed significantly between control
and Bap1f/f,Math1-cre, we used SEACR (Sparse Enrichment Analysis for CUT&RUN) (31) to call peak
regions before employing DiffBind (32) to identify regions where H2Aub is differentially abundant.
Volcano plots and minus-average plots, which depict the ratio (mutant/control) of H2Aub versus
significance or H2Aub abundance (concentration), respectively, demonstrated that there were far fewer
significantly dysregulated loci at P12 (Figure 1B and S1E) than in the adult (Figure 1C and S1F). The
predicted genotype effect predominated at both ages, with loss of the deubiquitinase Bap1 mostly
increasing its substrate H2Aub. The data also revealed, however, that 4% of regions at the adult timepoint
had undergone a paradoxical reduction in H2Aub (Figure 1C).
To understand the genomic contexts where H2Aub was perturbed and the effect of that perturbance,
we used ChIPseeker (33) to functionally annotate differentially ubiquitinated loci. In the mutant adult
cerebellum, H2Aub overabundance was most common within intronic, and especially distal intergenic,
sites—loci often associated with enhancer activity (Figure 1D)—while promoters were more commonly
depleted of H2Aub. We then input dysregulated H2Aub regions into the Cistrome repository (34, 35),
which compares user-defined regions of interest with thousands of chromatin mapping datasets to
produce a statistical measure of similarity known as the GIGGLE (Genomics Search Engine for Large-
Scale Integrated Genome) score. Loci abnormally enriched in H2Aub bore the greatest similarity to loci
(in the cerebellum or other brain regions) that are bound by the neuronal lineage-defining transcription
factors (Atoh1, Zic1/2, Tbr2, Neurod2, Wiz), the chromatin remodeler Chd4, or the RNA polymerase II
subunit Polr2a (Figure 1E). The highest GIGGLE scores belonged to acetylated modifications that are
associated with active cis-regulatory elements, such as H3K27ac and the acetylated version of the histone
variant H2A.z, and enhancer-associated H3K4me1. These results indicate that Bap1 loss leads to H2Aub
accumulation within euchromatic loci with features of active cis-regulatory elements.
A completely different collection of chromatin regulatory factors marked loci in the mutant
cerebellum that lost H2Aub, with components of PRC1 (Cbx7, Rnf2/Ring1b, C17orf96 and Pcgf6) or
PRC2 (Ezh2, Suz12, Jarid2) ranking the highest. These sites of H2Aub depletion most closely resembled
sites bound by PRC1 and PRC2 subunits in embryonic stem cells, i.e., stably polycomb repressed loci
(Figure 1F). The only histone modification to turn up was the PRC2 product H3K27me3 (especially
when detected in the cerebellum), confirming the heterochromatic nature of loci where Bap1 loss reduced
H2Aub.
These results suggest that the effects of Bap1 loss on H2Aub differ between euchromatic and
5
139
140
141
142
143
144
145
146
147
148
149
150
151
152
153
154
155
156
157
158
159
160
161
162
163
164
165
166
167
168
169
170
171
172
heterochromatic contexts. We therefore measured H2Aub within these distinct compartments. We first
detected H3K27ac (signifying euchromatic enhancers and promoters) and H3K27me3 (signifying
facultative heterochromatin) in the cerebellum of adult wild-type mice using CUT&RUN before calling
peak regions for each modification. Using a previously described percentile-based method (8) to
normalize CUT&RUN data across replicates of control and mutant genotypes, we quantified H2Aub
using edgeR (36) and limma (37), tools that permit the input of normalized data. In the mutant mice, the
H2Aub signal was intensified within H3K27ac peak regions (as detected by SEACR) (Figure 1G);
aggregate heatmaps depicting normalized H2Aub CUT&RUN data centered around H3K27ac peak
summits (as identified using the MACS2 narrow algorithm (38)) showed increased H2Aub genome-wide
(Figure 1H). H2Aub was also increased by Bap1 loss within active enhancers, as previously identified in
the adult cerebellum using hidden Markov modeling of chromatin states with the ChromHMM tool (8,
39) (Figure 1I). On the other hand, depicting H2Aub in relation to peak regions detected for H3K27me3
showed that Bap1 loss reduced H2Aub within facultative heterochromatin (Figure 1J and 1K) and
within bivalent loci identified using ChromHMM (Figure 1L). Differential analysis and metagene
heatmaps within protein-coding genes showed H2Aub depletion (Figure 1M and 1N), for example at the
heterochromatic Anxa2r1 locus (Figure 1O).
Seeking an explanation for the unexpected reduction of heterochromatic H2Aub, we used
CUT&RUN to examine the abundance of the PRC1 ubiquitin ligase Ring1B. Within the HoxA-HoxD
loci, which become repressed across tissues after embryonic development, both H2Aub and Ring1B were
reduced (Figure 2A). Other heterochromatic protein-coding genes also exhibited H2Aub depletion
accompanied by reduced Ring1B binding (Figure 2B). Aggregate heatmaps depicting Ring1B
CUT&RUN data within loci exhibiting significant reduction of H2Aub revealed concurrent reduction of
Ring1B at these sites (Figure 2C).
Bap1 loss in cerebellar granule neurons therefore bidirectionally perturbs H2Aub, with the direction
of perturbation dependent upon genomic context: H2Aub was enriched primarily within euchromatic cis-
regulatory loci while being depleted within facultative heterochromatin. In the next sections, we examine
the effects of Bap1 loss on the compartment-defining histone modifications of euchromatin and
heterochromatin.
Inverse relationship between euchromatin-associated H3K27ac and H2Aub in Bap1 mutant
Along with accessibility and transcriptional activity, histone acetylation is one of the defining
features of euchromatin. Having observed H2Aub gains within H3K27ac peak regions, we sought to
determine whether Bap1 loss affected H3K27ac per se. Detection of H3K27ac via CUT&RUN in control
and Bap1f/f,Math1-cre cerebellum segregated replicates according to genotype in PCA, whether data were
6
173
174
175
176
177
178
179
180
181
182
183
184
185
186
187
188
189
190
191
192
193
194
195
196
197
198
199
200
201
202
203
204
205
206
acquired at the early or late stage (Figure S2A). Using SEACR and Diffbind, we found that a greater
number of H3K27ac peak regions lost this mark than gained it at both time points, though more loci were
affected in adult neurons (Figure 3A, 3B, S2B and S2C). Plotting CUT&RUN data within ChromHMM-
defined active enhancers identified at the late timepoint (8) demonstrated that H3K27ac was generally
reduced by Bap1 loss (Figure 3C). At both ages, regions exhibiting significant reductions in H3K27ac
were much more likely to reside within distal intergenic loci, whereas regions with increased or invariant
H3K27ac were much more likely to reside within promoters or introns (Figure 3D and S2D). Comparing
up- and down-regulated H3K27ac loci with the Cistrome database demonstrated a striking inversion of
the pattern for H2Aub: regions that gained H3K27ac most resembled binding sites of Polycomb factors in
embryonic stem cells, whereas regions that lost H3K27ac were most similar to lineage-defining
transcription factors (as detected in the cerebellum) (Figure 3E, S2E and S2F). Comparisons with
histone modifications broke down less clearly, although regions in which Bap1 deletion reduced H3K27ac
were most similar to acetylated, euchromatin-associated modifications, as mapped in the cerebellum
and/or other brain regions (Figure 3F , S2E and S2F).
We next asked whether these euchromatic H2Aub gains and H3K27ac losses occurred at the same
loci. The intergenic enhancer upstream of the Mctp1 gene (encoding a calcium sensor that modulates
synaptic transmission) exhibited simultaneous increase in H2Aub and decreased H3K27ac, as did the
enhancer upstream of Rims1 (which encodes a presynaptic active zone protein that regulates
neurotransmitter release, synaptic vesicle docking, and calcium-dependent plasticity) (Figure 3G).
Importantly, these regions lacked heterochromatin-associated H3K27me3, confirming that Bap1 operates
within sites where PRC1 and PRC2 act independently. To determine if this inverse relationship between
H2Aub and H3K27ac held up genome-wide, we measured each modification within H3K27ac peak
regions identified in the cerebellum of adult controls and quantified the log₂ ratio of each modification,
mutant/control. This confirmed a strong inverse relationship: sites that gained H2Aub tended to lose
H3K27ac and vice versa (Figure 3H). In a lengthy intergenic enhancer cluster, increased H2Aub was
associated with reductions in H3K27ac at several loci (Figure 3I). Interestingly, the enhancer-associated
histone variant H2A.Z and its acetylated version H2A.Zac, both of which appeared in the Cistrome
analysis among upregulated H2Aub and downregulated H3K27ac sites, were decreased within these same
loci. Thus, Bap1 loss perturbs histone modifications beyond H2Aub and even histone variants at
enhancers, the combination of which is predicted to weaken their activity.
Neuronal Bap1 loss perturbs PRC2-dependent histone modifications
PRC1 deposition of H2Aub recruits and activates PRC2 to methylate H3K27 and concentrate
H3K27me3 within facultative heterochromatin (40); PRC2 can also generate H3K27me1 or H3K27me2.
7
207
208
209
210
211
212
213
214
215
216
217
218
219
220
221
222
223
224
225
226
227
228
229
230
231
232
233
234
235
236
237
238
239
240
Thee three methylation states show distinct genomic distributions (41). We previously showed that
H3K27me1 becomes enriched within gene bodies of transcribed genes over the course of normal
neurodevelopment and undergoes more remodeling than heterochromatin-associated H3K27me3 (8).
PRC2, then, like PRC1, regulates both euchromatin and heterochromatin, so we used CUT&RUN to
understand how PRC2-dependent histone modifications are affected by neuronal deletion of Bap1.
PCA segregated samples according to genotype along PC1, whether H3K27me3 CUT&RUN data
were acquired at the early or late timepoints (Figure S3A). H3K27me3 became more perturbed over time
and, at both ages, slightly more H3K27me3 peak regions were significantly downregulated than
upregulated by Bap1 loss (Figure 4A,B and S3B,C). In the early cerebellum, upregulated regions were
more likely to fall within promoters, whereas downregulated regions were more likely to be intronic or
intergenic (Figure 4C). In the late cerebellum, however, downregulated H3K27me3 peak regions became
similarly distributed to invariant regions. The most informative comparisons were across ages:
downregulated regions in the late cerebellum were distributed quite differently from downregulated
regions in the early cerebellum, whereas upregulated regions lost most of their promoter representation by
adulthood, while gaining distal intergenic representation.
Bap1-dependent regulation of H3K27me1 was also evident in PCA of CUT&RUN data, which
segregated samples by genotype along PC1 (Figure S4A). Quantification and differential analysis within
20519 protein coding genes identified 344 loci where H3K27me1 was downregulated (1.7%) and 458
where it was upregulated (2.2%) (Figure S4B and S4C). Interestingly, loci in which Bap1 loss increased
H3K27me1 harbored very little of the other assayed histone modifications (H3K27me3, H2Aub or
H3K27ac) (Figure S4D). In contrast, loci with depletion of H3K27me1 did have alterations in other
Bap1-dependent histone marks. For instance, Lsamp, which encodes an adhesion molecule crucial for
axon guidance, exhibited concurrent decreases in H3K27me1 and H3K27ac accompanied by increases in
H2Aub and H3K27me3 (Figure 4D). This collection of chromatin phenotypes reduced mRNA levels
(Figure 4E).
Observing increased H3K27me3 in the same regions as H2Aub prompted us to compare these
modifications within heterochromatin (H3K27me3 peak regions detected in the adult control cerebellum)
using log2-log2 plots, which showed a slight positive correlation upon loss of Bap1 (Figure 4F). This
suggests that while H2Aub and PRC2-dependent H3K27me3 can coincide at the genomic scale, the
relationship is weak and regulation at individual loci is largely independent. After observing genes where
H3K27me1 depletion was accompanied by increased H3K27me3, we also compared these modifications,
this time within protein-coding genes. Bap1-dependent changes in these modifications were weakly anti-
correlated, but also with a low R2 value (Figure 4G). In the absence of Bap1, then, loci that lose
monomethylation tend to gain trimethylation, and vice versa. But, once more, the mutual dependence
8
241
242
243
244
245
246
247
248
249
250
251
252
253
254
255
256
257
258
259
260
261
262
263
264
265
266
267
268
269
270
271
272
273
274
between these modifications was detectable only at the genomic scale: changes in one cannot predict the
effect on the other at individual loci. Moreover, because the mutant Bap1 allele completely removes
Bap1, we cannot separate changes that follow dysregulation of H2Aub from changes that are
compensatory or a secondary consequence of Bap1’s absence from chromatin-regulatory complexes.
Changes in histone modifications upon Bap1 deletion correspond to changes in chromatin
accessibility
The relationship between H2Aub and chromatin accessibility is complex. In one study employing
reconstituted nucleosome arrays and Drosophila lacking PR-DUB, the authors concluded that H2Aub
antagonizes compaction and increases accessibility at Polycomb targets, casting the modification as
opposing the condensed state (42). Work in mammalian cells has supported the opposite interpretation:
loss of Bap1 in embryonic stem cells caused H2Aub to accumulate diffusely across the genome, titrating
PRC2 away from its targets and driving genome-wide H3K27me3 deposition and general chromatin
compaction (43), consistent with independent evidence that Bap1 constrains pervasive H2Aub to preserve
the transcriptional potential of the genome (44). A recent study provides a way to reconcile these
opposing results by demonstrating that H2Aub acts in both directions depending upon its binding partner,
interfering with canonical PRC1 association with nucleosomes while facilitating H1-dependent
condensation (45). Because histone acetylation independently promotes accessibility, it remains difficult
to predict how the opposing changes in colocalized H2Aub and H3K27ac that we observed would
ultimately influence accessibility in the setting of Bap1 loss.
To determine how Bap1 loss affects chromatin accessibility, we performed ATAC-seq (Assay for
Transposase-Accessible Chromatin using sequencing) (46) in cerebellar neurons from adult
Bap1f/f,Math1-cre mice. PCA segregated replicates by genotype along PC1 (Figure S5A). DiffBind
identified >11,000 differentially accessible regions out of more than 80,000 consensus peak regions (as
detected by MACS2 narrow (38)), 10% of which became more accessible and 5% less so (Figures S5B
and S5C). Functional annotation demonstrated that regions with increased accessibility were more likely
to fall within promoters, whereas regions with decreased accessibility were more likely to reside within
intronic and intergenic regions (Figure S5D). Cistrome comparisons paralleled findings from chromatin
mapping experiments: loci with significantly increased accessibility showed similar genomic distributions
to those of PRC components mapped in ESCs, and modifications associated with repressed enhancers
(H3K4me1 and H3K27me3) (Figure S5E). In contrast, loci with significantly decreased accessibility
overlapped with binding sites of lineage-defining transcription factors and loci harboring active enhancer-
associated histone modifications (H3K27ac and acetylated H2A.z), as detected in cerebellar datasets
(Figure S5F). These results suggest that Bap1 regulates chromatin accessibility differently in
9
275
276
277
278
279
280
281
282
283
284
285
286
287
288
289
290
291
292
293
294
295
296
297
298
299
300
301
302
303
304
305
306
307
308
heterochromatin versus euchromatin, prompting us to examine ATAC-seq data in relation to H3K27ac and
H3K27me3. Quantification of the ATAC-seq signal within H3K27ac peak regions from adult control
cerebellum showed a predominant decrease, indicating that Bap1 loss reduces accessibility within active
enhancers (Figure S5G). For example, intronic enhancers within Rph3al (a.k.a. Noc2, which regulates
calcium ion-dependent exocytosis and synaptic vesicle fusion) showed reductions in ATAC-seq seq signal
accompanied by reduced H3K27ac (Figure S5H). In contrast, differential analysis of ATAC-seq data
within H3K27me3 peak regions showed an overall increase (Figure S5I), and we identified
heterochromatic loci where increases in chromatin accessibility were accompanied by simultaneous
reductions of H3K27me3 and H2Aub with increases in H3K27ac (Figure S5J). Thus, the changes in
chromatin state mediated by loss of Bap1 evoke the expected changes in chromatin accessibility within
both heterochromatin and euchromatin.
Loss of Bap1 weakens both euchromatin and heterochromatin identity
The genes Lsamp and Zbtb20 (encoding a transcriptional repressor that regulates neurogenesis)
each harbor multiple sites where decreases in H3K27ac colocalized with regions exhibiting increased
H3K27me3 and H2Aub, all while H3K27me1 was depleted throughout the gene body (Figure 4D and
5A). In these cases, excessive H2Aub in mutant neurons did stimulate ectopic deposition of H3K27me3.
This collection of changes effectively swaps activating euchromatin-associated modifications (H3K27ac
and H3K27me1) for repressive heterochromatic modifications (H2Aub and H3K27me3), resulting in
reduction of the corresponding mRNAs (Figure 4E and 5B). To determine if Bap1 loss caused
euchromatin to become more like heterochromatin on a genome-wide basis, we examined H3K27me3
specifically within H3K27ac peak regions in the cerebellum of control mice. PCA separated samples by
genotype (Figure S6A), while differential peak analysis demonstrated that many more sites experienced a
significant increase in H3K27me3 than experienced a decrease (Figure 5C and S6B). Aggregate plots
depicting normalized H3K27me3 CUT&RUN data centered around H3K27ac peak summits (Figure 5D),
ChromHMM-defined active enhancers (Figure 5E) or active promoters (Figure 5F) all demonstrated an
overall increase in H3K27me3 in the Bap1f/f,Math1-cre cerebellum.
We wondered whether Bap1 loss would also alter the state of heterochromatin so that it comes to
more closely resemble euchromatin. This appears to be the case: we observed typically heterochromatic
loci such as the Myc promoter gain the activating modifications H3K4me3 and H3K27ac while losing the
repressive modifications H2Aub and H3K27me3 (Figure 5G). Genome-wide quantitation of H3K27ac
within H3K27me3 peak regions produced PCA plots that separated samples by genotype (Figure S6C),
while differential peak analysis showed an overall increase in H3K27ac (Figure 5H and S6D). Similarly,
aggregate heatmaps depicting H3K27ac CUT&RUN data centered around H3K27me3 peak summits or
10
309
310
311
312
313
314
315
316
317
318
319
320
321
322
323
324
325
326
327
328
329
330
331
332
333
334
335
336
337
338
339
340
341
342
bivalent promoters showed an overall increase in H3K27ac in the Bap1 mutant cerebellum (Figure 5I
and 5J). Thus, heterochromatic loci became euchromatinized by Bap1 loss.
Given that several histone modifications controlled by Bap1 and/or H2Aub regulate chromatin
accessibility, we examined Bap1-mediated perturbations of epigenetic state in relation to changes in
chromatin accessibility by generating aggregate heatmaps of H3K27ac, H3K27me3 and H2Aub
CUT&RUN data within the two classes of differentially accessible loci. Loci where Bap1 loss increased
accessibility also exhibited increased H3K27ac and decreased H3K27me3 (Figure 5K). Conversely, loci
with decreased accessibility underwent changes in the opposite direction: reduction of H3K27ac with
increased H3K27me3 (Figure 5L). Thus, the direction of accessibility change at each class of loci
matches the direction predicted by the combination of changes in histone modifications.
After observing multiple aspects of chromatin state change in complementary ways at specific loci,
we used permutation analysis, which measures the non-randomness of the overlap between peak sets in
comparison to randomized genomic distributions, to gain a more generalized and quantitative
understanding of the relationship between loci whose attributes are altered by Bap1 loss. After 104
permutations, this analysis verified that increases in accessibility in ATAC-seq experiments were more
likely to coincide with increases in H3K27ac and decreases in H3K27me3, whereas decreases in
accessibility were more likely to coincide with sites where H3K27ac declined and H3K27me3 increased
(Figure 5M). The relation between accessibility and H2Aub was subtler: increased accessibility occurred
equally often within sites where H2Aub increased as where it decreased; and sites where accessibility
decreased were slightly more enriched for loci that experienced significant elevation of H2Aub. This
result accords with aggregate plots showing that H2Aub was increased by Bap1 loss at loci with
decreased accessibility (Figure 5L), whereas no overall genotype difference in H2Aub occurred at loci
where accessibility increased (Figure 5K). Increases in H2Aub were also associated with loci showing
gains in H3K27me3 and loss of H3K27ac , whereas decreases in H2Aub were accompanied by reduction
in HK27me3 and increases in H3K27ac. These results suggest that in postmitotic neurons the
consequences of H2Aub for chromatin accessibility are exerted largely through the complexes that
deposit and read it rather than through the modification on its own.
Collectively, these observations demonstrate that Bap1 controls a multifaceted chromatin regulatory
program that is required for proper formation of both euchromatin and heterochromatin. Remarkably,
neuronal loss of Bap1 weakens the identity of both compartments so that euchromatin becomes more like
heterochromatin while heterochromatin comes to more closely resemble euchromatin.
Bap1-dependent changes in histone modifications drive differential expression of specific classes of
neuronal genes
11
343
344
345
346
347
348
349
350
351
352
353
354
355
356
357
358
359
360
361
362
363
364
365
366
367
368
369
370
371
372
373
374
375
376
To confirm that these changes in chromatin state perturb the transcriptional program of neurons, we
sequenced RNA isolated from the cerebellum of Bap1f/f,Math1-cre mice at both early and late timepoints.
To focus on gene expression and minimize the effect of post-transcriptional regulation, we input RNA
extracted from cerebellar nuclei. PCA segregated replicates by age along PC1 followed by genotype along
PC2 (Figure 6A). PC1 accounted for 91% of the variation between samples, demonstrating the extent to
which the neuronal transcriptome is developmentally remodeled. Samples became increasingly
segregated along PC2 over time. Accordingly, more transcripts were significantly dysregulated at the late
compared to the early timepoint (Figure 6B and 6C). At both ages, about twice as many genes were
downregulated as upregulated. For example, the Syt1 locus, which encodes the calcium sensor
synaptotagmin-1 that resides in synaptic vesicles and triggers neurotransmitter release via membrane
fusion, became heterochromatinized (increased H2Aub and H3K27me3 but reduced H3K27ac, H3K4me3
and H3K27me1) (Figure 6D and 6E). Conversely, the normally repressed gene Iqsec3, which encodes an
organizer of GABAergic synapses, instead became euchromatinized, losing H2Aub and H3K27me3 while
gaining H3K27ac, resulting in overexpression of its mRNA (Figure 6F and 6G).
To compare changes in chromatin state to changes in transcription more broadly, we quantified
these histone modifications within significantly overexpressed or underexpressed genes before computing
the log₂ ratio of each modification’s abundance (mutant over control). At the adult stage, overexpressed
genes tended to harbor more H3K27ac and H3K4me3 and less H2Aub (Figure 6H), while
underexpressed genes gained H3K27me3 as they lost the euchromatin-associated modifications
H3K27ac, H3K4me3 and H3K27me1. Foreshadowing the adult pattern, genes that were overexpressed in
the early cerebellum showed gains in H3K27ac and loss of H2Aub, while underexpression was associated
with increased H3K27me3 (Figure 6I). These results suggest that transcriptional dysregulation in neurons
of Bap1f/f,Math1-cre mice stems from the combined effects of changes in multiple histone modifications
whose progressive perturbation intensifies transcriptomic dysregulation as granule neurons mature.
To understand the pathways most influenced by Bap1 loss, we examined the ontology of
significantly dysregulated genes at the late stage. Biological processes, cellular components and
molecular functions were all highly enriched for terms related to specific aspects of neuronal function
(Figure 6J). Differentially abundant transcripts encoded major classes of cell adhesion molecules
including several members of the cadherin, proto-cadherin, contactin, collagen, ephrin, laminin and
semaphorin families (Figure S7A). With the notable exception of Onecut2, which was overexpressed,
transcriptional regulators responsible for differentiation and maintenance of cell identity were strongly
underexpressed, including 4 members of the zinc finger of the cerebellum (ZIC) family (Figure S7B).
Transcripts encoding regulators of synaptic vesicle exocytosis, including the synaptotagmins Syt1, Syt2,
Syt4, Syt6, Syt7 and Syt9 (membrane proteins on synaptic vesicles that act as the primary Ca²⁺ sensors
12
377
378
379
380
381
382
383
384
385
386
387
388
389
390
391
392
393
394
395
396
397
398
399
400
401
402
403
404
405
406
407
408
409
410
for neurotransmitter release) were dysregulated (Figures 6E and S7C). Multiple transcripts encoding
calcium channels or associated proteins were also dysregulated, as was expression of calcium-permeable
glutamate receptors. Expression of multiple inositol triphosphate and ryanodine receptors was increased,
the effect of which would be to amplify calcium currents. Meanwhile, underexpression of mitochondrial
calcium importers would prolong calcium signaling in synapses (Figure S7D). Collectively, these
transcriptomic perturbances are very likely to disrupt synaptic activity, connectivity and neuronal identity.
Bap1 loss in excitatory neurons of the forebrain causes behavioral abnormalities and seizures
The preceding observations establish that Bap1 controls chromatin state in neurons to direct a gene
expression program that should be required for the development and function of neural circuits. This
requirement is consistent with work in Xenopus, where Bap1 loss prevents both the silencing of
pluripotency genes and the activation of lineage-commitment genes in ectoderm, mesoderm and neural
crest (47). Although the conditional deletion of Bap1 in cerebellar granule cells was enough to impair
motor learning (Figures 1A and S1A), we wanted to examine a greater variety of behaviors. We therefore
crossed Bap1f/f mice with mice expressing cre recombinase under the control of the Emx1 promoter,
which is expressed in ~90% of excitatory glutamatergic neurons in the forebrain, as well as in cortical glia
arising from the same lineage (48). Bap1f/f,Emx1-cre mice exhibited obvious postnatal growth delay, a
phenotype commonly associated with atypical neurodevelopment in mice (49, 50) (Figures 7A and 7B).
We then subjected 5- to 9-month-old male and female Bap1f/f,Emx1-cre mice to various behavioral assays
(see Methods). In the open field, the mutant mice showed striking and highly penetrant abnormalities in
exploratory behavior: instead of investigating the entire arena, mutants hugged its perimeter and corners
(Figure 7C), resulting in fewer entries to both the center and the border zones (Figure 7D) and less time
spent within the center (Figure 7E). This behavior was not due to locomotor deficits, as mutants’ average
and maximum speeds were, if anything, slightly greater than those of their wild-type littermates (Figure
7F). Neurotypical mice instinctively bury exposed marbles in bedding, but this goal-directed behavior
was essentially absent (Figure 7G). When suspended by the tail, the mice clasped their hindlimbs, a
nonspecific sign of neurological dysfunction (Figure 7H).
Careful observation of Bap1f/f,Emx1-cre mice suggested that they experienced episodes resembling
absence seizures. Some KIS patients with BAP1 mutation also experience seizures (20). To examine the
neurophysiology of the forebrain of Bap1f/f,Emx1-cre mice, we performed continuous
electroencephalography (EEG) recordings over three consecutive days in adult mice (age range 7 to 11
months). Spike electrical discharges (rapid, paroxysmal, high-amplitude waveforms that signify
hyperexcitability) were exceptionally rare in control mice but frequent in mutants (Figure 7I), which
experienced thousands of these daily epileptiform discharges (Figure 7J).
13
411
412
413
414
415
416
417
418
419
420
421
422
423
424
425
426
427
428
429
430
431
432
433
434
435
436
437
438
439
440
441
442
443
444
Frequency-response curves of EEG data demonstrated increased power within the 3-8 Hz range
across all vigilance states, but especially during wakefulness and rapid eye movement (REM) sleep
(Figure 7K). Decomposition of EEG data (i.e., spectral analysis) into delta (0.5–4 Hz), theta (4–8 Hz),
alpha (8–13 Hz), beta (13–30 Hz) and gamma (>30 Hz) bands revealed alterations in cortical oscillatory
dynamics in the mutant brain. During wakefulness and REM sleep, mutants exhibited increases in the
power density of low-frequency delta and theta bands accompanied by reductions in higher frequency
bands (Figure 7L). Similar changes were present to a lesser degree during slow wave sleep. These results
indicate that Bap1 loss shifts cortical oscillatory power toward lower frequencies and blunts the spectral
differences that normally distinguish vigilance states.
To determine whether these network-level changes reflect altered excitability at the single-neuron
level, we performed whole-cell patch clamp recordings of CA1 pyramidal neurons in acute hippocampal
slices from one-month-old animals. Increased resting membrane potential (VM) of mutant neurons
indicated depolarization at baseline (Figure 8A), and mutant neurons fired more spontaneous action
potentials (APs) (Figures 8B and 8C). The rheobase (i.e., the minimum injected current needed to
stimulate an AP) was lower in mutant neurons (Figure 8D and 8E), and injection of 50, 100 or 150 pA
(but not higher currents) evoked greater numbers of APs in mutant neurons (Figures 8F and 8G) without
affecting their frequency (Figure 8H) or amplitude (Figure 8I).
Excitatory postsynaptic potentials (EPSPs) within dendrites are generated by the combined influx
of Na⁺ through voltage-gated Sodium channels and Ca²⁺ through AMPA and NMDA glutamate receptors.
These inward currents depolarize the postsynaptic membrane, increasing the probability that the sum of
postsynaptic currents exceeds the threshold required to generate an AP. To determine whether the
hyperexcitability phenotype observed in mutant neurons arises from alterations in synaptic transmission,
we recorded extracellular field EPSPs (fEPSPs) in response to electrical stimulation of Schaffer collateral
axons originating from CA3 and terminating in the stratum radiatum of CA1. The amplitude of triggered
fiber volleys was comparable between the control and mutant across increasing stimulus intensities
(Figures 8J and 8K), indicating that Bap1 is not required for the integrity of glutamatergic projections
from CA3 onto CA1. However, the input-output curves of fEPSPs demonstrated that mutant neurons have
exaggerated synaptic responses, as evidenced by the trend toward greater fEPSP amplitude (Figure 8L)
and significantly increased slope (Figure 8M). Thus, the hippocampus of Bap1f/f,Emx1-cre mice
experiences higher levels of basal synaptic transmission (responses to single APs rather than high-
frequency stimulation) between CA3-CA1 synapses.
Finally, to determine if high basal synaptic transmission influences short-term plasticity in the
mutant hippocampus, we performed paired-pulse facilitation (PPF), a paradigm in which two APs are
triggered in rapid succession. The amplitude of both the first and second fEPSPs were increased in the
14
445
446
447
448
449
450
451
452
453
454
455
456
457
458
459
460
461
462
463
464
465
466
467
468
469
470
471
472
473
474
475
476
477
478
hippocampus of Bap1f/f,Emx1-cre mice, with the amplitude increasing during the second stimulation
(Figure 8N). The ratio of the initial slope of the second deflection to the initial slope of the first (i.e. the
paired-pulse ratio, PPR) was decreased in the mutant (Figure 8O). A decrease in PPR indicates that
enhanced presynaptic neurotransmitter release contributes to the observed increase in basal synaptic
transmission and affects short-term synaptic plasticity. These observations reflect profound changes in
synaptic transmission upon Bap1 loss.
DISCUSSION
Our examination of Bap1 deletion in neurons revealed two fundamental roles for this
deubiquitinating enzyme in the mammalian brain. First, Bap1 controls expression of essential neuronal
genes via a multifaceted program of chromatin regulation. Bap1 seems to be especially important for
regulating genes required for synaptic function, including genes regulating calcium signaling and/or
calcium responsiveness within the pre-synapse (e.g., synaptotagmins), post-synapse (voltage-gated
calcium channels), ER (ryanodine receptors) and mitochondria (Mcu and Micu1). As a result, Bap1-
mutant mice display profound behavioral abnormalities and epileptiform activity, while their hippocampal
pyramidal neurons exhibit intrinsic hyperexcitability. Second, and more surprisingly, loss of Bap1
bidirectionally perturbs histone modifications beyond its direct substrate H2Aub, weakening the
distinctions between euchromatin and heterochromatin.
Comparisons with previous studies in embryonic stem cells (44) highlight both conserved and
neuron-specific aspects of Bap1 function. In both types of cells, Bap1 regulates H2Aub beyond
facultative heterochromatin, as its loss in ESCs also increases H2Aub at euchromatic cis-regulatory
elements (44). This coordinated shift reflects a broader epigenetic reprograming that depletes
euchromatin-associated modifications to promote inappropriate heterochromatinization and suppress gene
expression. In neurons, we find that these changes occur at sites bound by lineage-defining transcription
factors, indicating that Bap1 is critical for maintaining enhancer and promoter landscapes required for
mature neuronal gene expression. Unlike Bap1-deficient ESCs, however, where H2Aub increases
genome-wide, including within heterochromatin (44), mutant neurons had a subset of heterochromatic
loci in mutant neurons that showed reduction of H2Aub, which is the opposite of what one would expect
after deletion of a deubiquitinase. We observed depletion of the ubiquitin ligase Ring1B from facultative
heterochromatin, as also occurs in ESCs (44), which provides a plausible explanation for the reduced
H2Aub we observe at these sites in neurons. Although this remains to be tested, we speculate that the
effects of Bap1 loss on heterochromatin reflect a redistribution of PRC1 activity between chromatin
compartments. PRC1 is present at very low abundance relative to its potential substrate. Even at the most
highly occupied sites in ESCs, fewer than 0.3 Ring1B molecules are estimated per kilobase of enriched
15
479
480
481
482
483
484
485
486
487
488
489
490
491
492
493
494
495
496
497
498
499
500
501
502
503
504
505
506
507
508
509
510
511
512
chromatin (51). But why, then, should heterochromatic H2Aub be preserved in Bap1-deficient ESCs but
lost in Bap1-deficient neurons? The two settings differ in important ways. ESCs divide rapidly and
express PRC1 and PRC2 components at high levels, so the Ring1B pool available to canonical PRC1 may
never become limiting. Neurons are postmitotic and express these components at far lower levels. The
compartments themselves also differ: facultative heterochromatin in mature neurons has been
consolidated over the course of development, while Polycomb repression in ESCs is poised, turns over
quickly, and is largely focal, with canonical PRC1 compacting chromatin at bivalent CpG-island
promoters and Hox clusters rather than across extended territories (52, 53). Finally, facultative
heterochromatin is particularly abundant in cerebellar granule neurons, which accumulate H3K27me3
more than ten-fold over broad genomic domains during postnatal maturation (54). The demand on
Ring1B scales with the extent of Polycomb domains, implying less capacity to safeguard against the
consequences of Ring1B redistribution to non-canonical complexes. A depletion that is buffered in one
context may be unmasked in the other.
Ring1B forms the shared catalytic core of both canonical PRC1 (which is recruited to H3K27me3
through Cbx-family chromodomain subunits (55)), and non-canonical PRC1 (which lacks Cbx subunits,
engages chromatin independently of H3K27me3, and is responsible for most genome-wide catalytic
H2Aub deposition (56)). Given these facts, we favor a model in which unrestrained non-canonical PRC1
activity across euchromatin in Bap1-deleted neurons draws upon the same limited Ring1B pool that
sustains canonical PRC1 within facultative heterochromatin, producing the coordinate loss of Ring1B and
H2Aub that we observe there. There are several possible reasons we did not detect a corresponding
accumulation of Ring1B at euchromatic sites by CUT&RUN. Non-canonical PRC1 constitutes the
mobile, low-residence population, whose transient engagement dispersed across a large genomic territory
yields little signal at any individual locus; it is much easier to capture the clustered, stable binding
characteristic of canonical PRC1 within Polycomb domains. Live-cell single-particle tracking of
endogenously tagged Ring1B demonstrates that only a small minority of molecules is chromatin-bound at
any moment, that the majority of binding events are brief, and that stable, long-lived association is largely
restricted to Polycomb bodies (51). The remaining molecules sample chromatin transiently, and it is this
mobile population that deposits H2Aub across most of the genome (56), including euchromatin, and
represents the fraction that Bap1 normally erases (44). In addition, euchromatic loci in mutant neurons
remain partly acetylated at H3K27, and because acetylation and trimethylation compete for the same
residue (57), residual acetylation constrains the substrate available for Cbx-mediated recruitment of
Ring1B-containing canonical PRC1. Since stable canonical PRC1 residence requires multivalent
engagement across a contiguous H3K27me3 domain, the ectopic H3K27me3 we detect is likely
insufficient to nucleate a detectable Ring1B footprint.
16
513
514
515
516
517
518
519
520
521
522
523
524
525
526
527
528
529
530
531
532
533
534
535
536
537
538
539
540
541
542
543
544
545
546
It is unclear why only a subset of euchromatic loci that gain H2Aub in Bap1-deficient neurons go
on to acquire H3K27me3. Although we observe clear colocalization of pathologically enriched H2Aub
and H3K27me3 in Bap1 mutant neurons at genes like Lsamp and Zbtb20 (Figures 4D and 5A), loci like
Mctp1 and Rims1 (Figure 3G) demonstrate that H2Aub accumulation alone is insufficient to direct PRC2
activity to a given locus. Residual H3K27 acetylation could restrict trimethylation at some sites, but other
locus-specific features (e.g., CpG dinucleotide density, DNA methylation state or the preexisting
abundance of H3K27me1 and H3K27me2 available for further methylation), could also determine
whether changes in H2Aub translates into H3K27me3. Resolving the relative contributions of canonical
and non-canonical PRC1 to this redistribution will require mapping the two complexes in postmitotic
neurons.
Chromatin accessibility is the most direct functional readout of the distinction between euchromatin
and facultative heterochromatin (58): euchromatin adopts an open configuration in which regulatory
elements remain available to transcription factors and RNA polymerase, whereas facultative
heterochromatin is compacted and refractory to binding (59). These differences emerge as a consequence
of the histone modifications that define each compartment. Acetylation of H3K27 neutralizes the positive
charge of histone tails, weakening their contacts with DNA and with adjacent nucleosomes (59), while
recruiting bromodomain-containing readers and remodelers that slide or evict nucleosomes (60).
H3K27me3, by contrast, is bound by the chromodomain subunits of canonical PRC1 (55), which bridges
nucleosomes to drive the compaction characteristic of the repressed state (53). H2Aub fits less neatly into
this scheme. Reconstituted nucleosome arrays and PR-DUB-deficient Drosophila indicate that H2Aub
antagonizes compaction at Polycomb targets (42), whereas work in mammalian cells supports the
opposite reading, with diffuse H2Aub accumulation in Bap1-null embryonic stem cells titrating PRC2
away from its targets and compacting chromatin genome-wide (43). A recent study helps reconcile these
results by showing that H2Aub acts in both directions depending on its binding partner, interfering with
canonical PRC1’s association with nucleosomes while facilitating H1-dependent condensation (45). Our
own data are consistent with this dual behavior: accessibility in the Bap1 mutant tracked H3K27ac and
H3K27me3 closely, whereas its relationship to H2Aub was weak and bidirectional.
At the transcriptional level, changes in chromatin state following neuronal Bap1 loss produce a
program that bears features of the Class 2 state in uveal melanoma, which, when it carries a BAP1
deletion, is associated with increased metastatic risk (61). In both neurons and uveal melanoma cells, loss
of Bap1 impairs expression of differentiation-associated transcription factors and structural genes, while
inappropriately enhancing expression of developmental regulators such as the tumor suppressor Onecut2
(62). In both systems, these changes are consistent with a shift toward a more progenitor-like state of
dedifferentiation (63). Neuronal Bap1 loss also leads to misexpression of genes involved in axon
17
547
548
549
550
551
552
553
554
555
556
557
558
559
560
561
562
563
564
565
566
guidance and cellular navigation, such as semaphorins, ephrins, and their receptors (e.g., Robo1),
paralleling cell adhesion pathways that melanoma co-opts to promote invasion and metastasis.
One of the major unanswered questions is whether the changes in epigenetic state result from
changes in H2Aub and its recruitment of factors to chromatin, or from losing a structural function of
Bap1, which interacts with many chromatin-regulatory protein complexes. Distinguishing these
possibilities will require separating Bap1's catalytic and scaffolding activities, which the null allele
employed here cannot achieve. Introducing individual patient variants, either as knock-in alleles in mice
or by editing human induced neurons, would allow such a distinction, since Küry-Isidor syndrome arises
from missense substitutions that fall mainly within the catalytic UCH domain rather than regions
mediating association with ASXL proteins and other chromatin-regulatory interactors (20). Comparing an
allelic series across chromatin, transcriptomic and physiological readouts would establish which of the
changes we describe depend upon deubiquitinase activity itself, and which follow from the loss of Bap1
as a platform for complex assembly.
In summary, Bap1 loss in mammalian neurons induces widespread reprogramming of
euchromatin and heterochromatin, disrupting the balance between enhancer activation and Polycomb-
mediated repression that is required for proper neuronal gene expression. We conclude that Bap1
preserves neuronal identity, supports circuit stability, and prevents hyperexcitable and developmentally
dysregulated neurologic states. These findings expand our understanding of ubiquitin-dependent
chromatin regulation in the brain and provide insight into the molecular basis of atypical
neurodevelopment.
18
567
568
569
570
571
572
573
574
575
576
577
578
579
580
581
582
583
584
585
586
587
588
589
590
591
592
593
594
595
596
597
598
599
600
MATERIALS & METHODS
Ethics statements
Mice were cared for in accordance with NIH guidelines. All experimental methods were approved
by the UCSD Institutional Committee on the Use and Care of Animals under the protocol number
S20121.
Mice
Sex matched-littermate controls on the Jackson Labs C57BL/6 background were used in all
experiments, and the ages at which animals were used is reported in figure legends. Animals were housed
in a 12:12 light:dark cycle. Biological replicates were sex-matched littermates, and both male and female
animals were examined during the major experiments. None of the main findings in this work varied by
sex in the tissues we examined.
Nuclear isolation
Mice were deeply anesthetized with isofluorane before decapitation and extraction of the whole
cerebellum. Tissue was finely minced in lysis buffer (Nuclei EZ Lysis Buffer supplemented with 1x Halt
combined protease and phosphatase Inhibitor cocktail, 10 mM sodium butyrate and 1.5 mM
iodoacetamide) prior to mechanical homogenization with Dounce A pestle in 2 mL of Lysis Buffer. This
suspension was incubated on ice for 5 minutes before centrifugation at 500 x g for 5 minutes, with
acceleration and deceleration on the lowest setting, at 4 degrees. The supernatant was aspirated and nuclei
were resuspended in EZ lysis buffer and incubated on ice for an additional 5 minutes prior to
centrifugation. Examination of nuclei under DIC demonstrated abundant intact nuclei. Intact nuclei were
nevertheless stained by trypan blue, indicating permeabilization by this procedure. Nuclei were
resuspended in CUT&RUN wash buffer and filtered through a 40 µm mesh. This process routinely
yielded ~6 million nuclei from a single mouse cerebellum. In our experience, detergent-based mechanical
homogenization disrupts the nuclei of larger neurons such as pyramidal and Purkinje cells, further
enriching cerebellar nuclear preparations for the small nuclei of granule neurons that serve as input for
chromatin mapping and transcriptomic studies.
Native CUT&RUN
For CUT&RUN performed under native conditions, nuclei were resuspended in wash buffer (20
mM HEPES pH 7.5, 150 mM NaCl, 0.5 mM spermidine supplemented with Roche complete EDTA-free
protease inhibitor tablet). Pelleted nuclei were resuspended twice in wash buffer and subjected to
19
601
602
603
604
605
606
607
608
609
610
611
612
613
614
615
616
617
618
619
620
621
622
623
624
625
626
627
628
629
630
631
632
633
634
centrifugation at 500 x g for 5 minutes. Nuclei were resuspended in 1 mL, filtered by gravity through a 40
µm mesh, and counted with a Countess automated cell counter (Invitrogen). In general, between 250,000-
500,000 nuclei were used in each CUT&RUN experiment, and the cellular input was kept the same across
all replicates. Nuclei were bound to Concanavalin A-coated magnetic beads. Following bead binding, the
supernatant was discarded, and the bead-nuclei slurry was resuspended in 75 uL of antibody buffer (wash
buffer supplemented with 0.005% digitonin and 2 mM EDTA) containing primary antibodies at a
concentration of 1:50. Antibody incubation was carried out overnight at 4°C on a nutator. Samples were
washed with cold cell permeabilization buffer (Wash buffer with 0.005% digitonin) to remove unbound
antibody and then incubated with 1x Protein A/G-MNase in cell permeabilization buffer for 15 minutes at
room temperature on a shaker. Samples were washed with cold cell permeabilization buffer and then
resuspended in 50 µL cell permeabilization buffer containing 2 mM CaCl2. Samples were digested for 2
hours at 4°C on a nutator. The digestion reactions were quenched with 33 µL of STOP buffer (340 mM
NaCl, 20 mM EDTA pH 8, 4 mM EGTA pH 7.7, RNase A 0.05 mg/mL, Glycogen 0.05 mg/mL). Samples
were then incubated for 10 minutes at 37°C to release digested chromatin fragments. The supernatant
containing CUT&RUN fragments was transferred to fresh tubes and DNA was purified via column
purification.
Fixed CUT&RUN
CUT&RUN detection of Ring1B under cross-linked conditions was carried out using a similar
protocol overall, but with important differences. After tissue homogenization, the supernatant was
removed and nuclei were briefly fixed in a 1% formaldehyde in PBS for 1 minute at room temperature
before quenching in 125 mM glycine. Subsequently, wash, antibody and permeabilization buffers were
supplemented with 1% Triton X-100 and 0.05% SDS. After stopping the chromatin digestions and
retrieving liberated chromatin complexes, decrosslinking was performed 0.8 µL of 10% SDS and 1 µL of
20 µg/µL of Proteinase K were added to each sample before incubation at 55°C overnight.
CUT&RUN library preparation and sequencing
Library preparation for CUT&RUN experiments was performed according to the manufacturer’s
specifications using SPRIselect beads. The concentration of DNA yielded from CUT&RUN was
measured on a Qubit device and between 25 and 50 ng of input DNA was used for library preparation,
with the amount of input DNA kept constant for all samples using the same antibody. 14 PCR cycles were
used in the final amplification with Illumina-compatible adaptors, universal i5 primer and barcoded i7
primers. After measuring the concentration of the amplified product on Nanodrop, all samples were
diluted to the same concentration for a given antibody, ranging from 10 ng/µL to 40 ng/µL. Tapestation
20
635
636
637
638
639
640
641
642
643
644
645
646
647
648
649
650
651
652
653
654
655
656
657
658
659
660
661
662
663
664
665
666
667
668
was performed to determine the precise concentration and fragment distribution range and to ensure there
were not significant amount adaptor dimers. Samples were pooled for sequencing 100 bp reads in paired-
end configuration on an Illumina NovaSeqX platform at the UCSD Institute for Genomic Medicine. We
aimed for >25 million reads for all samples except for H3K4me3, for which 15 million reads were
adequate given the narrow distribution of this histone modification.
CUT&RUN data processing
We used FastQC to evaluate read quality, including base quality scores and adapter contamination.
We first trimmed compressed sequenced reads in fastq.gz format were first trimmed using Trimmomatic
v.0.39 (64). We then aligned samples to the mm10 genome using Bowtie2 v.2.3.4.1 (65) and sorted BAM
files, removing unmapped fragments using SAMtools v.1.14 (66). During alignment we set a minimum
Phred score of 33 and employed the --dovetail option. We then converted alignment files to .bigwig files
using the bamCoverge tool in deepTools v.3.5.0 (67) with 50 bp bins. Duplicate reads were not removed.
Peak calling was performed using SEACR (Sparse Enrichment Analysis for CUT&RUN) v.1.1 (31) and
MACS2 (Model-Based Analysis of ChIP-seq) (38) in broad and narrow formats. Samples of different
replicates collected in the same experiment were normalized to each other to correct for different loading
between samples. .bigwig and .bed files were examined in Integrated Genome Viewer and visualized
using KaryoplotR v.1.22.0 (68). We obtained regional overlap between different samples using the
‘intersect’ functionality of BEDTools v.2.29.2 (69).
Normalization of CUT&RUN data
After parsing the genome into 50 bp bins, we identified local maxima and generated a blacklist of
271 intergenic peaks representing spurious alignment artifacts. We quantified the height of the leftover
authentic peaks and recorded the value of the height for the 99th percentile peak recorded. The ratio of the
values of the 99th percentile peaks between samples served as a scaling factor by which the height of
every bin in one sample could be multiplied to normalize samples to one another. Across the datasets
reported here, these scaling factors were small, generally requiring corrections of less than 10%. We
chose this data-driven approach over exogenous spike-in normalization for three reasons specific to
chromatin profiling of neuronal nuclei. First, spike-in strategies assume proportional recovery of
endogenous and exogenous chromatin, an assumption that is difficult to satisfy when neuronal nuclei lyse
and aggregate over the course of a two-day CUT&RUN protocol, because sample loss is not experienced
equally by the two chromatin pools. Second, the number of biological replicates required here exceeded
what could be processed in a single batch, and spike-in ratios are themselves sensitive to batch-to-batch
variation in the amount and quality of added chromatin. Third, spike-in normalization performs least well
21
669
670
671
672
673
674
675
676
677
678
679
680
681
682
683
684
685
686
687
688
689
690
691
692
693
694
695
696
697
698
699
700
701
702
in precisely the situation examined here, in which a perturbation redistributes a modification
bidirectionally between genomic compartments without changing its global abundance; under these
conditions a single global scaling factor derived from an exogenous standard can obscure locus-specific
effects. Consistent with this reasoning, immunoblot analysis of nuclear lysates demonstrated that total
H2Aub was unchanged between genotypes at both ages, confirming that the within-sample scaling
applied here is not compensating for a global shift in modification abundance. Differential analyses were
additionally performed with DiffBind using non-normalized BAM files as input, and these analyses
recovered the same direction of change at the loci reported.
DiffBind analysis
Differential peak analysis was performed on SEACR-called peaks using DiffBind (v. 3.14) (32)
with default parameters (summit width=400 bp). Non-normalized BAM files were used as input.
Significantly differential peaks were identified using DESeq2 (70) within Diffbind and only peaks with
FDR<0.05 were used for further analysis. We functionally annotated significant peaks using the
“annotatePeak” function of ChIPseeker (v.1.4) (33) and the TxDb.Mmusculus.UCSC.mm10.knownGene
(v. 3.10) database.
Heatmaps of signal grouped by developmental H2Aub dynamics
To visualize chromatin state changes across distinct categories of H2Aub regulation, peak regions
were grouped based on differential analysis of H2Aub signal using the DiffBind R package. Peaks were
categorized as no change, significantly increased, or significantly decreased over neurodevelopment
(adjusted p < 0.05, fold change > 1). For each group, CUT&RUN signal for H2Aub, H3K27me3,
H3K4me3 and H3K27ac across both timepoints was extracted from bigWig files using computeMatrix
with the reference-point setting. plotHeatmap was then used to visualize the signal centered around peak
summits and sorted by decreasing histone modification abundance.
Genome-wide aggregate heatmaps
Genome-wide heatmaps were generated using deepTools v.3.5. Plots centered around peak summits
were produced in reference to .bed files output from peak calling analysis of CUT&RUN data from
mouse tissue. computeMatrix was used to represent each locus from bigwig CUT&RUN data, where the
rows were loci and the columns were 50 bp bins. Matrices were subsequently plotted as heatmaps using
the plotHeatmap tool. For metagene plots, in which the signal in bigwig tracks was scaled to gene length,
the metagene flag was used during computeMatrix. plotHeatmap was then used to visualize the
abundance of different histone PTMs.
22
703
704
705
706
707
708
709
710
711
712
713
714
715
716
717
718
719
720
721
722
723
724
725
726
727
728
729
730
731
732
733
734
735
736
ATAC-seq
ATAC-seq experiments were using the commercial kit (Active Motif, #53150) under native
conditions. Following cell lysis in an ice-cold lysis buffer, nuclei were resuspended in 1 mL, filtered by
gravity through a 40 µm mesh, and counted with a Countess automated cell counter (Invitrogen). In
general, 200,000 nuclei were used in each ATAC-seq experiment. The nuclei were incubated in the
tagmentation master mix for 60 min at 37 °C, and the resulting DNA was purified using DNA purification
columns. PCR amplification of the tagmented DNA was then performed to generate libraries with
appropriately indexed primers. We performed a 2 step bead clean-up using SPRI select beads. The first
cleanup used a 0.55x ratio to remove large fragments, retaining unbound fragments in the supernatant.
The second used a 1.5x ratio to bind small fragment on the beads, which were retained and into library
preparation. DNA libraries were sequenced on the Illumina NovaSeq X Plus platform.
Permutation test
Pairwise enrichment between genomic region sets was assessed using permutation testing
implemented in the regioneReloaded R package. All genomic coordinates were mapped to the mm10
(GRCm38) reference genome. Region sets were loaded as GRanges objects using rtracklayer, restricted to
standard autosomes, and sorted prior to analysis. Permutation tests were performed with number of
overlapping regions (numOverlaps) as the test statistic and 10,000 permutation iterations using the
resampleRegions randomization function to generate the null distribution.
RNA-seq
Mouse cerebellar nuclei were isolated using the same nuclear isolation protocol for CUT&RUN.
Isolated nuclei were homogenized in Trizol. Total nuclear RNA was isolated using the column-based
Qiagen RNeasy kit. Library prep was performed using 1 ug of input RNA with the Illumina Stranded
Total RNA Prep Ligation with Ribo-Zero Plus kit. After depleting ribosomal RNA and performing cDNA
synthesis, the DNA was quantitated using Qubit and equalized across samples. Illumina-compatible dual
index primers were used during 15 PCR cycles. The prepared libraries were subsequently sequencing
using the Illumina NovaSeq S4 platform to generate 100 bp paired-end reads. Each sample was sequenced
to a depth of at least 30 million reads. RNA-seq libraries were aligned to mm10 using STAR v.2.7.9a (71).
Alignment files were then converted to bigwig format using deepTools bamCoverage v.3.5.0. DESeq2
(70) was then used to compute normalized reads and adjusted P-values.
Gene ontology analysis
23
737
738
739
740
741
742
743
744
745
746
747
748
749
750
751
752
753
754
755
756
757
758
759
760
761
762
763
764
765
766
767
768
769
770
A list of the chromosomal coordinates of all known genes in the mm10 genome along with their
Refseq IDs was obtained from the UCSC Genome Browser. The closest RefSeq genes within 10 kb of
each .bed file locus was obtained using BEDTools v.2.29.2. The gene IDs left after this filtering was the
gene list used as input to ShinyGO v.0.77 (72) for GO analysis with an FDR cutoff of 0.05. Rather than
testing against the whole-genome default, enrichment was evaluated against a custom background
comprising all genes detected as expressed in our cerebellar nuclear RNA-seq datasets. Restricting the
background in this way ensures that the reported enrichments reflect overrepresentation relative to the
neuronal transcriptome rather than relative to the genome as a whole, and therefore that terms describing
neuronal structure and function are not recovered simply as a consequence of the tissue examined.
Behavior analysis – general features of cohort that included Bap1f/f,Emx1-cre mice
Mice ranged from 5.1 to 9.2 months at the start of behavioral testing, with a mean of age of 7.1 and
a median age of 7.3 months. Behavioral testing took 2 months to complete. Mutants and their littermate
controls were included with similar numbers of males and females.
Behavior analysis of Bap1f/f,Math1-cre mice – Accelerating Rota-Rod
Mice ranged from 4.0 to 6.8 months at the start of testing, with a mean of age of 5.4 and a median
age of 5.2 months. We used the Accelerated Rota-Rod 7650 (UGO Basile, Italy), which has a diameter of
3 cm. All mice were placed in testing room 1 hour prior to experiment for habituation. On the first day of
testing, they were allowed to habituate on the stationary Rota-Rod for 5 minutes. During testing,
acceleration to a final speed of 40 rpm occurred over 5 minutes. Latency to fall was tested three times
each day with at least 30 minutes between trials. The trial was ended if a mouse stayed on for 15 minutes.
Testing occurred on 4 consecutive days.
Behavior analysis – open field test
The open field tests were conducted in a dedicated sound-attenuated room and in a square box of
composed of white reflective panels 60 × 60 × 40 cm (length × width × height). The center zone was
defined as the 40 x 40 cm central square of the arena, whereas the border zone was the 10 cm perimeter.
Testing was conducted in a dedicated behavioral room under dim illumination to minimize stress related
to lighting. An overhead video camera was positioned above the center of the arena and connected to a
computer running ANY -maze software for automated tracking and analysis. Before experimental trials, an
experiment file was created in ANY -maze for automated video-based tracking of the animal’s center
point/center of mass. All subjects were entered in the order they were tested. The camera view and arena
boundaries were aligned within the software, and tracking performance was verified using a mock animal
24
771
772
773
774
775
776
777
778
779
780
781
782
783
784
785
786
787
788
789
790
791
792
793
794
795
796
797
798
799
800
801
802
803
804
placed in the arena. If detection was suboptimal, video settings (e.g., contrast) were adjusted until the
software reliably tracked the mock animal. Animals were brought to the testing room and allowed to
acclimate for at least 30 min prior to testing. Testing was conducted in a dedicated behavioral room under
dim illumination (low-wattage lighting) to minimize stress. For each trial, a single animal was placed
gently in the center of the arena facing away from the experimenter, and recording/tracking was initiated
immediately. Animals were allowed to freely explore the arena for 15 min. At the end of the session, the
animal was removed and returned to its home cage. The arena was thoroughly cleaned between subjects
using 70% ethanol to reduce olfactory cues and prevent carryover effects, and the procedure was repeated
for all remaining animals. Locomotor activity was quantified as total distance traveled and line crossings
using a 16-grid overlay. Spatial distribution measures were used to evaluate center versus periphery
exploration and thigmotaxis, with the center zone defined as all area outside the periphery. Heatmaps and
Track plots were then generated using ANY -maze software the compare the activity of the animals.
Session videos and representative movement traces were saved for visualization and quality control.
Behavior analysis - Paw clasping
Mice were picked up by the distal third of their tails and observed for 10 s. They were given a
rating based on clasping of the front and/or back paws: 0 – no paw clasping, 1 – occasional clasping of
front paws, and 2 – almost constant clasping of front paws and occasional clasping of paws (73).
Behavior analysis - Marble burying test
Mice were placed individually in standard mouse cages containing bedding 5 cm in depth with 20
small marbles arranged in 4 evenly spaced rows of 5 on top of the bedding material. The number of
marbles buried (at least 2/3 covered by bedding) was determined in a 30 min trial (74).
Electroencephalography
The same behavioral cohort was subjected to electroencephalography (EEG) analysis after the
completion of behavioral experiments, as EEG requires surgery. Mice ranged from 7.0 to 11.0 months of
age, with a mean age of 9.0 months and a median age of 9.1 months. EEG was performed as previously
described (75–78). Mice ranged EEG data were recorded from stainless steel screw electrodes implanted
on the frontal and parietal bone over the hippocampus (coordinates: 2.0 mm posterior and 2.0 mm lateral
to bregma), and under general anesthesia (1-1.5% isoflurane). A fourth EEG electrode was implanted over
the cerebellum, used to ground the animal to reduce signal artifacts. Insulated leads from the EEG
electrodes were crimped to male pins (220-P02) and then cemented to the skull with dental acrylic.
Following surgical implantation mice were allowed 1-2 weeks to recover prior to the study. To record
25
805
806
807
808
809
810
811
812
813
814
815
816
817
818
819
820
821
822
823
824
825
826
827
828
829
830
831
832
833
834
835
836
837
838
EEG, mice were connected to commutators (PlasticOne) with flexible recording cables allowing their
unrestricted movements within the cage and habituated to the recording cages for 24h before 72h of test
recording. The EEG and EMG signals were amplified in a Grass Model 7D polygraph in a frequency
range of 0.3 to 10 KHz. The EEG and EMG were displayed on a computer monitor and stored with a
resolution of 128 Hz in the hard drive of a computer for the off-line analysis of the vigilance states and
spectral analysis using software supplied by Kissei Comptec. The polygraphic results were analyzed
semiautomatically by 15-second epochs and classified as W, SWS and REM sleep. The total time of these
vigilance states was calculated across the first & third day. In addition to standard sleep analysis, EEG
spectral power analysis in the different states of vigilance was performed by Fourier fast transformer
(FFT) analysis using 4-second epochs. The polygraphic results also were analyzed for the total number of
spikes, latency to a possible first seizure, as well as duration and frequency of any detected seizures.
Artifacts in the raw EEG recordings (electrical noise), induced by exploratory behavior and grooming,
were identified and excluded from the analysis. To detect epileptiform spikes, the following parameters
were used: Threshold amplitude value=200000 µV; Maximum amplitude value 1500000 µV; Minimum
spike duration=5 ms; and Maximum spike duration=60 ms.
Acute brain slice preparation
Hippocampal slices were prepared as described previously with some modifications (79). Briefly,
after anesthesia with isoflurane inhalation, mice were decapitated, and brains were rapidly removed from
skulls and placed into ice-cold oxygenated (95 % O2/5 % CO2) sucrose cutting solution with reduced Na+
and Ca2+ and increased Mg2+ for one minute. The sucrose cutting solution contained (in mM): 80 NaCl,
3.5 KCl, 25 NaHCO3, 1.25 NaH2PO4, 0.5 CaCl2, 4.5 MgCl2, 100 sucrose, and 10 glucose (osmolarity 330
– 340 mOsm; pH 7.3 – 7.4). The brain was cut at the midline, and the two hemispheres were mounted
onto the slicing stage of the vibratome at a 30-degree angle. The slicing chamber was filled with pre-
cooled sucrose cutting solution, and oxygenation was continuously provided throughout the procedure.
Transverse hippocampal slices (400 μm for whole-cell and field potential recording) were cut at 4 °C with
a vibratome (Leica VT1000 S) and transferred to oxygenated sucrose cutting solution pre-warmed at 34
°C. Slices were recovered in oxygenated sucrose cutting solution at 34 °C for 35 min and then in
oxygenated standard artificial cerebrospinal fluid (aCSF) at room temperature (RT) for 1 hour (hr) before
recording. The standard aCSF contained (in mM): 125 NaCl, 25 NaHCO3, 1.25 NaH2PO4, 3.5 KCl, 2
CaCl2, 1 MgCl2, and 10 glucose (osmolarity 300 – 310 mOsm; pH 7.3 – 7.4).
Whole-cell patch-clamp recording
Brain slices were transferred into the recording chamber with a continuous flow of oxygenated
26
839
840
841
842
843
844
845
846
847
848
849
850
851
852
853
854
855
856
857
858
859
860
861
862
863
864
865
866
867
868
869
870
871
872
standard artificial cerebrospinal fluid (aCSF). The bath aCSF was set to perfuse at a rate of 1.5-2.0 ml/min
and was heated to 30-32 °C with a temperature controller (Warner Instrument, TC-324B). A slice holder
was gently put on top of the slice, and then the slices were allowed to rest in the chamber for about 5 min
before all recordings. CA1 pyramidal neuron layer was identified under the microscope of the recording
rig and neurons with clear cell body morphology without swelling or atrophy were selected for patch-
clamp recording. Borosilicate glass pipettes (Warner, Hamden, CT) were pulled from a micropipette
puller (Sutter, Model P-97). The recording electrodes had a resistance of 2-5 MΩ when filled with an
internal solution containing (in mM) 130 K-gluconate, 1.0 CaCl2, 1.0 MgCl2, 5 EGTA, 10 HEPES, 4 Mg-
ATP, and 0.3 Na2-GTP (osmolarity 290±5 mOsm, pH 7.25 -7.30). A MultiClamp 700B amplifier, an
Axon Digidata 1550B, and Clampex 11 software (Molecular Devices, San Jose, CA, USA) were used for
data acquisition. Signals were digitized at 2-10 kHz, and filtered at 2 kHz. The liquid junction potential
was compensated for prior to forming the cell-attached mode for all recordings. A minimum of 2 GΩ seal
resistance was required before rupturing the membrane for whole-cell configuration. For each recorded
cell, membrane test in Clampex software was used to read the membrane capacitance (Cm), membrane
resistance (Rm), and access resistance (Ra) immediately after whole-cell configuration. Then current
clamp mode I=0 or IC without holding current was used to record membrane potential (VM) for 2 - 3
minutes. Rheobase (10 ms duration of minimal current injection from 0 pA with 10 pA increment until the
first AP appears), and excitability (500 ms duration of -50 to +500 pA current injection steps with 50 pA
increment) were also recorded in current clamp mode. Then cells were held at -70 mV in voltage clamp
mode for current recordings. After the holding current stabilized, a voltage step protocol was applied to
record the Na+ currents and K+ currents. Voltage steps were from -180 to +40 mV with 10 mV increment
and 25 ms duration, and with the baseline holding at -70 mV. Only the cells with Ra <= Rm x 20% and
less than 20% fluctuation of Ra were included for voltage clamp recordings. Each slice was kept in the
recording chamber for less than one hour to ensure that the cells were recorded in a healthy state.
Extracellular field potential recording
Brain slices were put in the recording chamber with the same perfusion rate and temperature as
described above. Extracellular field excitatory postsynaptic potentials (fEPSPs) were recorded in the
stratum radiatum layer of the CA1 hippocampal area. A gain of 100 was used in current clamp mode and
the signal was digitized at 2-10 kHz and filtered at 1 kHz. The stimulating electrode (FHC, cluster bipolar
microelectrodes) was placed in the stratum radiatum about 200 μm away from the cell body layer of CA1
pyramidal neurons. The recording electrode had a resistance of 2-5 MΩ when filled with standard aCSF
and was placed about 300 μm away from stimulating electrode in stratum radiatum. For each slice, 0 -10
V stimulating intensity with 1 V interval was used for the input-output test to define the baseline synaptic
27
873
874
875
876
877
878
879
880
881
882
883
884
885
886
887
888
889
890
891
892
893
894
895
896
897
898
899
900
901
902
903
transmission. ~50% of the stimulation intensity inducing the maximum fEPSP response in input-output
test was used for paired-pulse facilitation (PPF) recordings in the same slice. The paired-pulse facilitation
was performed by the delivery of two identical intensity stimuli at an interval of 50 ms, 3 sweeps were
recorded with a 20 s interval. Slices with fEPSP amplitude less than 0.4 mV at the maximum stimulation
(10 V) in the input-output test were not included for data analysis.
Data analysis and statistics (Electrophysiology)
All the data were analyzed in software ClampFit 11 (Molecular Devices San Jose, CA, USA). For
neuronal whole-cell recordings, the most negative membrane potential was used for resting VM
quantification, and the cell was counted as spontaneous firing if any spontaneous AP was observed during
the membrane potential recording period. The minimum current to induce the first AP was used to
quantify the rheobase of each cell. For excitability analysis, “Threshold search” function in ClampFit
software was used, the depolarization peaks with overshoot above 0 mV were identified as APs, and only
APs occurring during the current injection period were included for quantification of induced APs. For the
cells showing spontaneous firing at rest, APs occurring outside the range of current injection period were
not included for induced AP analysis. For fEPSP recordings, the peak amplitude and the slope (30% -
70%) of each sweep were plotted against the stimulating intensity accordingly for input-output (I-O)
plots. The ratio of paired-pulse facilitation (PPR) was calculated using the slope values of the two
fEPSPs: the slope of the second fEPSP divided by the slope of the first fEPSP.
Statistical analysis was performed using software GraphPad Prism 10. Means between two groups
were compared using an unpaired, two-tailed t-test. Repeated measures one-way ANOVA was used to
compare the number, frequency and amplitude of induced APs across current injection steps in whole-cell
recordings, as well as the input-output curves of fEPSP recordings between two groups. Fisher's exact test
was used to compare the incidence of spontaneous firing in the recorded neurons between two groups.
Data (summary) were presented as mean ± SEM. For all analysis, p ≤ 0.05 was considered statistically
significant.
Quantification And Statistical Analysis
For details of the quantitative methods used in CUT&RUN, RNA-seq and ATAC-seq, please refer
to the Methods sections and Figure Legends. We used GraphPad Prism 10 software to perform other
statistical tests described in Figure legends.
28
904
905
906
907
908
909
910
911
912
913
914
915
916
917
918
919
920
921
922
923
924
925
926
927
928
929
930
931
932
933
934
935
936
937
938
939
940
941
942
REFERENCES
1. 2. 3. 4. 5. 6. 7. 8. 9. 10. 11. 12. B. D. Strahl, C. D. Allis, The language of covalent histone modifications. 403 (2000).
Y. Zhao, B. A. Garcia, Comprehensive Catalog of Currently Documented Histone Modifications.
Cold Spring Harb. Perspect. Biol. 7, a025064 (2015).
F. Mattiroli, L. Penengo, Histone Ubiquitination: An Integrative Signaling Platform in Genome
Stability. Trends Genet. 37, 566–581 (2021).
N. P. Blackledge, A. M. Farcas, T. Kondo, H. W. King, J. F. McGouran, L. L. P. Hanssen, S. Ito, S.
Cooper, K. Kondo, Y. Koseki, T. Ishikura, H. K. Long, T. W. Sheahan, N. Brockdorff, B. M.
Kessler, H. Koseki, R. J. Klose, Variant PRC1 Complex-Dependent H2A Ubiquitylation Drives
PRC2 Recruitment and Polycomb Domain Formation. Cell 157, 1445–1459 (2014).
A. Laugesen, J. W. Højfeldt, K. Helin, Molecular Mechanisms Directing PRC2 Recruitment and
H3K27 Methylation. Mol. Cell 74, 8–18 (2019).
V. Kasinath, C. Beck, P. Sauer, S. Poepsel, J. Kosmatka, M. Faini, D. Toso, R. Aebersold, E.
Nogales, JARID2 and AEBP2 regulate PRC2 in the presence of H2AK119ub1 and other histone
modifications. Science 371, eabc3393 (2021).
S. Cooper, M. Dienstbier, R. Hassan, L. Schermelleh, J. Sharif, N. P. Blackledge, V. De Marco, S.
Elderkin, H. Koseki, R. Klose, A. Heger, N. Brockdorff, Targeting Polycomb to Pericentric
Heterochromatin in Embryonic Stem Cells Reveals a Role for H2AK119u1 in PRC2 Recruitment.
Cell Rep. 7, 1456–1470 (2014).
A. Parmar, A. Srinivasan, L. Krockenberger, A. Augustine, O. Gong, A. C. Bullard, R. Kalra, L.
Ledvin, D. Pilz, J. Tawil, C. E. Tea, K. C. Wang, O. Urso, L. M. Kaube, Y. Sun, R. Sasik, K. J.
Gaulton, K. M. Fisch, C. J. Ferguson, Polycomb repressive complexes 1 and 2 independently and
dynamically regulate euchromatin during cerebellar neurodevelopment. PLOS Genet. 21,
e1011843 (2025).
K. MacGillivray, D. Fusca, L. Tan, R. Aram, A. L. Saltzman, Widespread association of Polycomb
complex–deposited histone H2A monoubiquitylation with enhancers and neuronal gene regulation.
Sci. Adv. (2025).
S. B. Pierce, M. D. Stewart, S. Gulsuner, T. Walsh, A. Dhall, J. M. McClellan, R. E. Klevit, M.-C.
King, De novo mutation in RING1 with epigenetic effects on neurodevelopment. Proc. Natl. Acad.
Sci. 115, 1558–1563 (2018).
R. L. Borges, G. González-Blanco, H. Arigela, Y. Huang, L. D. Caeiro, N. Fattakhov, S. Lepore, L.
Garcia-Martinez, M. Maurice, P. D. Mehta, E. J. Park, K. MacGillivray, J. Nehru, M. Chau, M. C.
Robayo, C. Abad, A. Bilbao-Martinez, F. Monteiro, X. Luo, S. Tan, D. Bilbao, S. Sidoli, B. Di
Stefano, K. Walz, A. L. Saltzman, R. E. Verdun, R. Shiekhattar, L. Morey, Unbalanced chromatin
binding of Polycomb complexes drives neurodevelopmental disorders. Mol. Cell 86, 604-624.e16
(2026).
X. Luo, K. Schoch, S. V. Jangam, V. H. Bhavana, H. K. Graves, S. Kansagra, J. M. Jasien, N. Stong,
B. Keren, C. Mignot, C. Ravelli, Undiagnosed Diseases Network, H. J. Bellen, M. F. Wangler, V.
Shashi, S. Yamamoto, Rare deleterious de novo missense variants in Rnf2/Ring2 are associated
29
943
944
945
946
947
948
949
950
951
952
953
954
955
956
957
958
959
960
961
962
963
964
965
966
967
968
969
970
971
972
973
974
975
976
977
978
979
980
981
982
983
984
985
13. 14. 15. 16. 17. 18. 19. 20. with a neurodevelopmental disorder with unique clinical features. Hum. Mol. Genet. 30, 1283–
1292 (2021).
A. Dey, D. Seshasayee, R. Noubade, D. M. French, J. Liu, M. S. Chaurushiya, D. S. Kirkpatrick, V.
C. Pham, J. R. Lill, C. E. Bakalarski, J. Wu, L. Phu, P. Katavolos, L. M. LaFave, O. Abdel-Wahab,
Z. Modrusan, S. Seshagiri, K. Dong, Z. Lin, M. Balazs, R. Suriben, K. Newton, S. Hymowitz, G.
Garcia-Manero, F. Martin, R. L. Levine, V. M. Dixit, Loss of the Tumor Suppressor BAP1 Causes
Myeloid Transformation. Science 337, 1541–1546 (2012).
K. N. Smit, M. J. Jager, A. De Klein, E. Kiliҫ, Uveal melanoma: Towards a molecular
understanding. Prog. Retin. Eye Res. 75, 100800 (2020).
J. W. Harbour, M. D. Onken, E. D. O. Roberson, S. Duan, L. Cao, L. A. Worley, M. L. Council, K.
A. Matatall, C. Helms, A. M. Bowcock, Frequent Mutation of BAP1 in Metastasizing Uveal
Melanomas. Science 330 (2010).
J. Xu, Y. Kadariya, M. Cheung, J. Pei, J. Talarchek, E. Sementino, Y. Tan, C. W. Menges, K. Q. Cai,
S. Litwin, H. Peng, J. Karar, F. J. Rauscher, J. R. Testa, Germline Mutation of Bap1 Accelerates
Development of Asbestos-Induced Malignant Mesothelioma. Cancer Res. 74, 4388–4397 (2014).
M. Carbone, H. Yang, H. I. Pass, T. Krausz, J. R. Testa, G. Gaudino, BAP1 and cancer. Nat. Rev.
Cancer 13, 153–159 (2013).
Y. Kadariya, M. Cheung, J. Xu, J. Pei, E. Sementino, C. W. Menges, K. Q. Cai, F. J. Rauscher, A. J.
Klein-Szanto, J. R. Testa, Bap1 Is a Bona Fide Tumor Suppressor: Genetic Evidence from Mouse
Models Carrying Heterozygous Germline Bap1 Mutations. Cancer Res. 76, 2836–2844 (2016).
S. Walpole, A. L. Pritchard, C. M. Cebulla, R. Pilarski, M. Stautberg, F. H. Davidorf, A. De La
Fouchardière, O. Cabaret, L. Golmard, D. Stoppa-Lyonnet, E. Garfield, C.-N. Njauw, M. Cheung,
J. A. Turunen, P. Repo, R.-S. Järvinen, R. Van Doorn, M. J. Jager, G. P. M. Luyten, M. Marinkovic,
C. Chau, M. Potrony, V. Höiom, H. Helgadottir, L. Pastorino, W. Bruno, V. Andreotti, B.
Dalmasso, G. Ciccarese, P. Queirolo, L. Mastracci, K. Wadt, J. F. Kiilgaard, M. R. Speicher, N.
Van Poppelen, E. Kilic, R. T. Al-Jamal, I. Dianzani, M. Betti, C. Bergmann, S. Santagata, S.
Dahiya, S. Taibjee, J. Burke, N. Poplawski, S. J. O’Shea, J. Newton-Bishop, J. Adlard, D. J.
Adams, A.-M. Lane, I. Kim, S. Klebe, H. Racher, J. W. Harbour, M. L. Nickerson, R. Murali, J. M.
Palmer, M. Howlie, J. Symmons, H. Hamilton, S. Warrier, W. Glasson, P. Johansson, C. D. Robles-
Espinoza, R. Ossio, A. De Klein, S. Puig, P. Ghiorzo, M. Nielsen, T. T. Kivelä, H. Tsao, J. R. Testa,
P. Gerami, M.-H. Stern, B. B. Paillerets, M. H. Abdel-Rahman, N. K. Hayward, Comprehensive
Study of the Clinical Phenotype of Germline BAP1 Variant-Carrying Families Worldwide. JNCI J.
Natl. Cancer Inst. 110, 1328–1341 (2018).
S. Küry, F. Ebstein, A. Mollé, T. Besnard, M.-K. Lee, V. Vignard, T. Hery, M. Nizon, G. M. S.
Mancini, J. C. Giltay, B. Cogné, K. McWalter, W. Deb, H. Mor-Shaked, H. Li, R. E. Schnur, I. M.
Wentzensen, A.-S. Denommé-Pichon, C. Fourgeux, F. W. Verheijen, E. Faurie, R. Schot, C. A.
Stevens, D. J. Smits, E. Barr, R. Sheffer, J. A. Bernstein, C. L. Stimach, E. Kovitch, V. Shashi, K.
Schoch, W. Smith, R. H. Van Jaarsveld, A. C. E. Hurst, K. Smith, E. H. Baugh, S. G. Bohm, E.
Vyhnálková, L. Ryba, C. Delnatte, J. Neira, D. Bonneau, A. Toutain, J. A. Rosenfeld, S. Audebert-
Bellanger, B. Gilbert-Dussardier, S. Odent, F. Laumonnier, S. I. Berger, A. C. M. Smith, F.
Bourdeaut, M.-H. Stern, R. Redon, E. Krüger, R. Margueron, S. Bézieau, J. Poschmann, B. Isidor,
Rare germline heterozygous missense variants in BRCA1-associated protein 1, BAP1, cause a
syndromic neurodevelopmental disorder. Am. J. Hum. Genet. 109, 361–372 (2022).
30
986
987
988
989
990
991
992
993
994
995
996
997
998
999
1000
1001
1002
1003
1004
1005
1006
1007
1008
1009
1010
1011
1012
1013
1014
1015
1016
1017
1018
1019
1020
1021
21. 22. 23. 24. 25. 26. 27. 28. 29. 30. 31. 32. 33. E. C. West, M. Chiappetta, A. A. Mattingly, M. T. Congedo, J. Evangelista, A. Campanella, C.
Sassorossi, S. Flamini, T. Rossi, M. Pistoni, L. Abenavoli, S. Margaritora, F. Lococo, L. Boccuto,
BRCA1‐associated protein 1: Tumor predisposition syndrome and Kury‐Isidor syndrome, from
genotype–phenotype correlation to clinical management. Clin. Genet. 105, 589–595 (2024).
J. C. Scheuermann, A. G. De Ayala Alonso, K. Oktaba, N. Ly-Hartig, R. K. McGinty, S. Fraterman,
M. Wilm, T. W. Muir, J. Müller, Histone H2A deubiquitinase activity of the Polycomb repressive
complex PR-DUB. Nature 465, 243–247 (2010).
A. Hoischen, B. W. M. Van Bon, B. Rodríguez-Santiago, C. Gilissen, L. E. L. M. Vissers, P. De
Vries, I. Janssen, B. Van Lier, R. Hastings, S. F. Smithson, R. Newbury-Ecob, S. Kjaergaard, J.
Goodship, R. McGowan, D. Bartholdi, A. Rauch, M. Peippo, J. M. Cobben, D. Wieczorek, G.
Gillessen-Kaesbach, J. A. Veltman, H. G. Brunner, B. B. B. A. De Vries, De novo nonsense
mutations in ASXL1 cause Bohring-Opitz syndrome. Nat. Genet. 43, 729–731 (2011).
A. Srivastava, K. C. Ritesh, Y.-C. Tsan, R. Liao, F. Su, X. Cao, M. C. Hannibal, C. E. Keegan, A. M.
Chinnaiyan, D. M. Martin, S. L. Bielas, De novo dominant ASXL3 mutations alter H2A
deubiquitination and transcription in Bainbridge–Ropers syndrome. Hum. Mol. Genet. 25, 597–
608 (2016).
V. Matei, S. Pauley, S. Kaing, D. Rowitch, K. W. Beisel, K. Morris, F. Feng, K. Jones, J. Lee, B.
Fritzsch, Smaller inner ear sensory epithelia in Neurog1 null mice are related to earlier hair cell
cycle exit. Dev. Dyn. 234, 633–650 (2005).
G. G. Consalez, D. Goldowitz, F. Casoni, R. Hawkes, Origins, Development, and
Compartmentation of the Granule Cells of the Cerebellum. Front. Neural Circuits 14, 611841
(2021).
M. Hibi, T. Shimizu, Development of the cerebellum and cerebellar neural circuits. Dev. Neurobiol.
72, 282–301 (2012).
L. Hammelrath, S. Škokić, A. Khmelinskii, A. Hess, N. Van Der Knaap, M. Staring, B. P. F.
Lelieveldt, D. Wiedermann, M. Hoehn, Morphological maturation of the mouse brain: An in vivo
MRI and histology investigation. NeuroImage 125, 144–152 (2016).
M. P. Meers, T. D. Bryson, J. G. Henikoff, S. Henikoff, Improved CUT&RUN chromatin profiling
tools. eLife 8, 1–16 (2019).
P. J. Skene, S. Henikoff, An efficient targeted nuclease strategy for high-resolution mapping of DNA
binding sites. eLife 6, 1–35 (2017).
M. P. Meers, D. Tenenbaum, S. Henikoff, Peak calling by Sparse Enrichment Analysis for
CUT&RUN chromatin profiling. Epigenetics Chromatin 12, 42 (2019).
R. Stark, G. Brown, DiffBind: Differential binding analysis of ChIP-Seq peak data.
G. Yu, L.-G. Wang, Q.-Y. He, ChIPseeker: an R/Bioconductor package for ChIP peak annotation,
comparison and visualization. Bioinformatics 31, 2382–2383 (2015).
31
1022
1023
1024
1025
1026
1027
1028
1029
1030
1031
1032
1033
1034
1035
1036
1037
1038
1039
1040
1041
1042
1043
1044
1045
1046
1047
1048
1049
1050
1051
1052
1053
1054
1055
1056
1057
1058
34. 35. 36. 37. 38. 39. 40. 41. 42. 43. 44. 45. 46. S. Mei, Q. Qin, Q. Wu, H. Sun, R. Zheng, C. Zang, M. Zhu, J. Wu, X. Shi, L. Taing, T. Liu, M.
Brown, C. A. Meyer, X. S. Liu, Cistrome Data Browser: a data portal for ChIP-Seq and chromatin
accessibility data in human and mouse. Nucleic Acids Res. 45, D658–D662 (2017).
R. Zheng, C. Wan, S. Mei, Q. Qin, Q. Wu, H. Sun, C.-H. Chen, M. Brown, X. Zhang, C. A. Meyer,
X. S. Liu, Cistrome Data Browser: expanded datasets and new tools for gene regulatory analysis.
Nucleic Acids Res. 47, D729–D735 (2019).
M. D. Robinson, D. J. McCarthy, G. K. Smyth, edgeR : a Bioconductor package for differential
expression analysis of digital gene expression data. Bioinformatics 26, 139–140 (2010).
M. E. Ritchie, B. Phipson, D. Wu, Y. Hu, C. W. Law, W. Shi, G. K. Smyth, limma powers
differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res.
43, e47–e47 (2015).
Y. Zhang, T. Liu, C. A. Meyer, J. Eeckhoute, D. S. Johnson, B. E. Bernstein, C. Nusbaum, R. M.
Myers, M. Brown, W. Li, X. S. Liu, Model-based Analysis of ChIP-Seq (MACS). Genome Biol. 9,
R137 (2008).
J. Ernst, M. Kellis, Chromatin-state discovery and genome annotation with ChromHMM. Nat.
Protoc. 12, 2478–2492 (2017).
N. P. Blackledge, N. R. Rose, R. J. Klose, Targeting Polycomb systems to regulate gene expression:
modifications to a complex story. Nat. Rev. Mol. Cell Biol. 16, 643–649 (2015).
K. J. Ferrari, A. Scelfo, S. Jammula, A. Cuomo, I. Barozzi, A. Stützer, W. Fischle, T. Bonaldi, D.
Pasini, Polycomb-Dependent H3K27me1 and H3K27me2 Regulate Active Transcription and
Enhancer Fidelity. Mol. Cell 53, 49–62 (2014).
J. Bonnet, I. Boichenko, R. Kalb, M. Le Jeune, S. Maltseva, M. Pieropan, K. Finkl, B. Fierz, J.
Müller, PR-DUB preserves Polycomb repression by preventing excessive accumulation of
H2Aub1, an antagonist of chromatin compaction. Genes Dev. 36, 1046–1061 (2022).
E. Conway, F. Rossi, D. Fernandez-Perez, E. Ponzo, K. J. Ferrari, M. Zanotti, D. Manganaro, S.
Rodighiero, S. Tamburri, D. Pasini, BAP1 enhances Polycomb repression by counteracting
widespread H2AK119ub1 deposition and chromatin condensation. Mol. Cell 81, 3526-3541.e8
(2021).
N. A. Fursova, A. H. Turberfield, N. P. Blackledge, E. L. Findlater, A. Lastuvkova, M. K. Huseyin,
P. Dobrinić, R. J. Klose, BAP1 constrains pervasive H2AK119ub1 to control the transcriptional
potential of the genome. Genes Dev. 35, 749–770 (2021).
J. Zhao, J. Lan, M. Wang, C. Liu, Z. Fang, A. Song, T. Zhang, L. Wang, B. Zhu, P. Chen, J. Yu, G.
Li, H2AK119ub1 differentially fine-tunes gene expression by modulating canonical PRC1- and
H1-dependent chromatin compaction. Mol. Cell 84, 1191-1205.e7 (2024).
J. D. Buenrostro, P. G. Giresi, L. C. Zaba, H. Y. Chang, W. J. Greenleaf, Transposition of native
chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins
and nucleosome position. Nat. Methods 10, 1213–1218 (2013).
32
1059
1060
1061
1062
1063
1064
1065
1066
1067
1068
1069
1070
1071
1072
1073
1074
1075
1076
1077
1078
1079
1080
1081
1082
1083
1084
1085
1086
1087
1088
1089
1090
1091
1092
1093
1094
1095
1096
47. 48. 49. 50. 51. 52. 53. 54. 55. 56. 57. 58. 59. 60. J. N. Kuznetsov, T. H. Aguero, D. A. Owens, S. Kurtenbach, M. G. Field, M. A. Durante, D. A.
Rodriguez, M. L. King, J. W. Harbour, BAP1 regulates epigenetic switch from pluripotency to
differentiation in developmental lineages giving rise to BAP1-mutant cancers. Sci. Adv. 5,
eaax1738 (2019).
J. A. Gorski, T. Talley, M. Qiu, L. Puelles, J. L. R. Rubenstein, K. R. Jones, Cortical Excitatory
Neurons and Glia, But Not GABAergic Neurons, Are Produced in the Emx1-Expressing Lineage.
J. Neurosci. 22, 6309–6314 (2002).
W. Wang, W. N. Frankel, Overlaps, gaps, and complexities of mouse models of Developmental and
Epileptic Encephalopathy. Neurobiol. Dis. 148, 105220 (2021).
E. Damianidou, L. Mouratidou, C. Kyrousi, Research models of neurodevelopmental disorders: The
right model in the right place. Front. Neurosci. 16, 1031075 (2022).
M. K. Huseyin, R. J. Klose, Live-cell single particle tracking of PRC1 reveals a highly dynamic
system with low target site occupancy. Nat. Commun. 12, 887 (2021).
S. Kundu, F. Ji, H. Sunwoo, G. Jain, J. T. Lee, R. I. Sadreyev, J. Dekker, R. E. Kingston, Polycomb
Repressive Complex 1 Generates Discrete Compacted Domains that Change during
Differentiation. Mol. Cell 65, 432-446.e5 (2017).
R. Eskeland, M. Leeb, G. R. Grimes, C. Kress, S. Boyle, D. Sproul, N. Gilbert, Y. Fan, A. I.
Skoultchi, A. Wutz, W. A. Bickmore, Ring1B Compacts Chromatin Structure and Represses Gene
Expression Independent of Histone Ubiquitination. Mol. Cell 38, 452–464 (2010).
I. Laas, M. R. Paul, N. Bhanu, L. Feng, E.-E. Govek, T. S. Carroll, C. D. Allis, M. E. Hatten, K.
Mätlik, Long-term maintenance of H3K27me3 in postmitotic neurons is dispensable for gene
expression regulation.
E. Bernstein, E. M. Duncan, O. Masui, J. Gil, E. Heard, C. D. Allis, Mouse Polycomb Proteins Bind
Differentially to Methylated Histone H3 and RNA and Are Enriched in Facultative
Heterochromatin. Mol. Cell. Biol. 26, 2560–2569 (2006).
N. A. Fursova, N. P. Blackledge, M. Nakayama, S. Ito, Y. Koseki, A. M. Farcas, H. W. King, H.
Koseki, R. J. Klose, Synergy between Variant PRC1 Complexes Defines Polycomb-Mediated
Gene Repression. Mol. Cell 74, 1020-1036.e8 (2019).
D. Pasini, M. Malatesta, H. R. Jung, J. Walfridsson, A. Willer, L. Olsson, J. Skotte, A. Wutz, B.
Porse, O. N. Jensen, K. Helin, Characterization of an antagonistic switch between histone H3
lysine 27 methylation and acetylation in the transcriptional regulation of Polycomb group target
genes. Nucleic Acids Res. 38, 4958–4969 (2010).
P. Trojer, D. Reinberg, Facultative Heterochromatin: Is There a Distinctive Molecular Signature?
Mol. Cell 28, 1–13 (2007).
K. A. Gelato, W. Fischle, Role of histone modifications in defining chromatin structure and
function. Biol. Chem. 389, 353–363 (2008).
M. Yun, J. Wu, J. L. Workman, B. Li, Readers of histone modifications. Cell Res. 21, 564–578
(2011).
33
1097
1098
1099
1100
1101
1102
1103
1104
1105
1106
1107
1108
1109
1110
1111
1112
1113
1114
1115
1116
1117
1118
1119
1120
1121
1122
1123
1124
1125
1126
1127
1128
1129
1130
1131
1132
1133
61. 62. 63. 64. 65. 66. 67. 68. 69. 70. 71. 72. 73. 74. 75. T. R. See, G. Stålhammar, S. Phillips, H. E. Grossniklaus, BAP1 Immunoreactivity Correlates with
Gene Expression Class in Uveal Melanoma. Ocul. Oncol. Pathol. 6, 129–137 (2020).
M. Rotinen, S. You, J. Yang, S. G. Coetzee, M. Reis-Sobreiro, W.-C. Huang, F. Huang, X. Pan, A.
Yáñez, D. J. Hazelett, C.-Y. Chu, K. Steadman, C. M. Morrissey, P. S. Nelson, E. Corey, L. W. K.
Chung, S. J. Freedland, D. Di Vizio, I. P. Garraway, R. Murali, B. S. Knudsen, M. R. Freeman,
ONECUT2 is a targetable master regulator of lethal prostate cancer that suppresses the androgen
axis. Nat. Med. 24, 1887–1898 (2018).
K. A. Matatall, O. A. Agapova, M. D. Onken, L. A. Worley, A. M. Bowcock, J. W. Harbour, BAP1
deficiency causes loss of melanocytic cell identity in uveal melanoma. BMC Cancer 13, 371
(2013).
A. M. Bolger, M. Lohse, B. Usadel, Trimmomatic: a flexible trimmer for Illumina sequence data.
Bioinformatics 30, 2114–2120 (2014).
B. Langmead, S. L. Salzberg, Fast gapped-read alignment with Bowtie 2. Nat. Methods 9, 357–359
(2012).
H. Li, B. Handsaker, A. Wysoker, T. Fennell, J. Ruan, N. Homer, G. Marth, G. Abecasis, R. Durbin,
1000 Genome Project Data Processing Subgroup, The Sequence Alignment/Map format and
SAMtools. Bioinformatics 25, 2078–2079 (2009).
F. Ramírez, D. P. Ryan, B. Grüning, V. Bhardwaj, F. Kilpert, A. S. Richter, S. Heyne, F. Dündar, T.
Manke, deepTools2: a next generation web server for deep-sequencing data analysis. Nucleic Acids
Res. 44, W160–W165 (2016).
B. Gel, E. Serra, karyoploteR: an R/Bioconductor package to plot customizable genomes displaying
arbitrary data. Bioinformatics 33, 3088–3090 (2017).
A. R. Quinlan, I. M. Hall, BEDTools: A flexible suite of utilities for comparing genomic features.
Bioinformatics 26, 841–842 (2010).
M. I. Love, W. Huber, S. Anders, Moderated estimation of fold change and dispersion for RNA-seq
data with DESeq2. Genome Biol. 15, 550 (2014).
A. Dobin, C. A. Davis, F. Schlesinger, J. Drenkow, C. Zaleski, S. Jha, P. Batut, M. Chaisson, T. R.
Gingeras, STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15–21 (2013).
S. X. Ge, D. Jung, R. Yao, ShinyGO: a graphical gene-set enrichment tool for animals and plants.
Bioinformatics 36, 2628–2629 (2020).
R. Lalonde, C. Strazielle, Brain regions and genes affecting limb-clasping responses. Brain Res.
Rev. 67, 252–259 (2011).
R. M. J. Deacon, Digging and marble burying in mice: simple methods for in vivo identification of
biological impacts. Nat. Protoc. 1, 122–124 (2006).
P. B. Hedlund, S. Huitron-Resendiz, S. J. Henriksen, J. G. Sutcliffe, 5-HT7 Receptor Inhibition and
Inactivation Induce Antidepressantlike Behavior and Sleep Pattern. Biol. Psychiatry 58, 831–837
(2005).
34
1134
1135
1136
1137
1138
1139
1140
1141
1142
1143
1144
1145
1146
1147
1148
1149
1150
1151
1152
1153
1154
76. 77. 78. 79. 80. 81. 82. 83. S. Huitron-Resendiz, M. Sanchez-Alavez, D. N. Wills, B. F. Cravatt, S. J. Henriksen,
Characterization of the Sleep-wake Patterns in Mice Lacking Fatty Acid Amide Hydrolase. Sleep
27, 857–865 (2004).
S. Huitron-Resendiz, Urotensin II Modulates Rapid Eye Movement Sleep through Activation of
Brainstem Cholinergic Neurons. J. Neurosci. 25, 5465–5474 (2005).
P. Bourgin, S. Huitrón-Reséndiz, A. D. Spier, V. Fabre, B. Morte, J. R. Criado, J. G. Sutcliffe, S. J.
Henriksen, L. De Lecea, Hypocretin-1 Modulates Rapid Eye Movement Sleep through Activation
of Locus Coeruleus Neurons. J. Neurosci. 20, 7760–7765 (2000).
Y. Du, F. H. Brennan, P. G. Popovich, M. Zhou, Microglia maintain the normal structure and
function of the hippocampal astrocyte network. Glia 70, 1359–1379 (2022).
A. R. Quinlan, I. M. Hall, BEDTools: a flexible suite of utilities for comparing genomic features.
Bioinformatics 26, 841–842 (2010).
S. Heinz, C. Benner, N. Spann, E. Bertolino, Y. C. Lin, P. Laslo, J. X. Cheng, C. Murre, H. Singh, C.
K. Glass, Simple Combinations of Lineage-Determining Transcription Factors Prime cis-
Regulatory Elements Required for Macrophage and B Cell Identities. Mol. Cell 38, 576–589
(2010).
G. Yu, L.-G. Wang, Y. Han, Q.-Y. He, clusterProfiler: an R Package for Comparing Biological
Themes Among Gene Clusters. OMICS J. Integr. Biol. 16, 284–287 (2012).
R. Malinverni, D. Corujo, B. Gel, M. Buschbeck, regioneReloaded: evaluating the association of
multiple genomic region sets. Bioinformatics 39, btad704 (2023).
35
1155
1156
1157
1158
1159
1160
1161
1162
1163
1164
1165
1166
1167
1168
1169
1170
1171
1172
1173
1174
1175
1176
1177
1178
1179
1180
1181
1182
1183
1184
1185
1186
1187
ACKNOWLEDGMENTS
Funding: This work was supported by the NIH grant R01NS138924 (C.J.F.). We thank Vicky Brandt for
helpful discussions and comments on the manuscript. We also thank Kristen Jepsen and the UCSD IGM
for performing and advising on sequencing experiments. Challana Tea was supported in part by an
institutional Ruth L. Kirschstein National Research Service Award from the National Institute of General
Medical Sciences, T32GM145427. This publication includes data generated at the UC San Diego IGM
Genomics Center utilizing an Illumina NovaSeq X Plus that was purchased with funding from a National
Institutes of Health SIG grant (#S10 OD026929).
Author contributions:
Conceptualization: K.C.W., C.E.T., C.J.F.
Data curation: K.C.W., C.E.T., A.P., C.J.F.
Formal analysis: K.C.W., C.E.T., A.P.,
Funding Acquisition: C.E.T., K.D., C.J.F.
Investigation: K.C.W., C.E.T., A.C.B., K.T., Z.S., J.T., A.R., C.J.F.
Methodology: C.E.T., Y.D., K.T., Z.S., A.R.,
Software: C.E.T., A.C.B., D.L., A.P.,
Supervision: C.E.T., K.D., C.J.F.
Validation: K.C.W., C.E.T., K.T., A.R., C.J.F.
Visualization: K.C.W., C.E.T., Y.D., A.C.B., D.L., A.P., A.R., K.D., C.J.F.
Writing – original draft: K.C.W., Y.D.,
Writing – review & editing: K.C.W., C.E.T., C.J.F.
Competing Interests: The authors declare no competing interests.
Data, code, and materials availability: Further information and requests for resources and/or reagents
should be directed to the Lead Contact, Cole Ferguson (cferguson@health.ucsd.edu). Any other non-
commercially available reagent or input biomaterial generated or used during this study will be made
freely available to interested researchers. Raw and analyzed data is available under the NCBI BioProject
PRJNA1150596. All original code is available at: https://github.com/FergusonLab. Any additional
information required to reanalyze the data reported in this paper is available from the lead contact upon
request.
36
1188
1189
1190
1191
1192
1193
1194
1195
1196
1197
1198
1199
1200
1201
1202
1203
1204
1205
1206
1207
1208
1209
1210
1211
1212
1213
1214
1215
1216
1217
1218
1219
1220
1221
MAIN FIGURE LEGENDS
Figure 1. Loss of Bap1 in cerebellar granule neurons of adult Bap1f/f,Math1-cre mice increases
H2Aub in euchromatin and depletes H2Aub from heterochromatin and gene bodies.
(A) Performance in the accelerating Rota-Rod test of cerebellar-dependent motor performance and
learning. Three daily trials are over four consecutive days of testing. Animals ranged in age from 4 to
7 months. Error bars SEM, p-values by two-way ANOVA with Sidak test.
(B) DiffBind (32) differential analysis of H2Aub abundance in CUT&RUN data from the P12 cerebellum
(n=4 control, n=4 mutant), as analyzed within consensus peak regions identified using SEACR (31).
Correcting for multiple comparisons, the false discovery rate (FDR) was set to <0.05.
(C) Differential analysis of H2Aub abundance in CUT&RUN data from the adult (age 3 months)
cerebellum (n=6 control, n=6 mutant).
(D) ChIPseeker (33) functional annotation of different classes of H2Aub peak regions based upon how
they behave in Bap1 mutant neurons.
(E) GIGGLE score comparisons between peak regions with increased H2Aub and the Cistrome (34, 35)
database of chromatin factors (left) and modified histones (right). The source of experimental data is
indicated by the color key.
(F) Comparisons between peak regions with decreased H2Aub and the Cistrome database.
(G) Volcano plot depicting H2Aub abundance, as detected in normalized CUT&RUN data from control
and mutant cerebellum, within H3K27ac peak regions identified by SEACR in the control
cerebellum. Quantitation and adjusted p-values were computed using edgeR (36) and limma (37)
(n=6 control, 6 mutant).
(H) Genome-wide aggregate heatmaps depicting normalized H2Aub CUT&RUN data from the adult
cerebellum within H3K27ac peak summits, as detected using MACS2 narrow (38).
(I) Genome-wide aggregate heatmaps depicting normalized H2Aub CUT&RUN data from the adult
cerebellum within strong active enhancers, as detected using Chromatin Hidden Markov
(ChromHMM) Modeling (39) of CUT&RUN data previously acquired from the adult cerebellum (8).
(J) Volcano plot depicting H2Aub abundance in normalized CUT&RUN data from late cerebellum,
within H3K27me3 peak regions. Quantitation and adjusted p-values were computed using edgeR and
limma (n=6 control, 6 mutant).
(K) Genome-wide aggregate heatmaps depicting normalized H2Aub CUT&RUN data from the adult
cerebellum within H3K27me3 peak summits, as detected using MACS2 narrow.
(L) Genome-wide aggregate heatmaps depicting normalized H2Aub CUT&RUN data from the adult
cerebellum within bivalent loci, as identified using ChromHMM analysis of CUT&RUN data from
37
1222
1223
1224
1225
1226
1227
1228
1229
the adult cerebellum.
(M) Volcano plot depicting H2Aub abundance in normalized CUT&RUN data from the late cerebellum
within protein coding genes. Quantitation and adjusted p-values were computed using edgeR and
limma (n=6 control, 6 mutant).
(N) Metagene heatmaps depicting normalized H2Aub CUT&RUN data from the adult cerebellum across
protein-coding genes.
(O) Normalized CUT&RUN data for the indicated histone modifications in the adult control and mutant
cerebellum at the heterochromatic Anxa2r1 locus.
38
1230
1231
1232
1233
1234
1235
1236
1237
1238
1239
Figure 2: Reduced association of Ring1B from heterochromatic loci where neuronal Bap1 loss
depletes H2Aub in the adult mouse cerebellum.
(A) Normalized CUT&RUN tracks depicting the Hox loci in mouse cerebellar granule neurons. Each
locus exhibits H2Aub depletion accompanied by reduced association of Ring1B. H2Aub and
H3K27me3 were detected under native conditions, Ring1B under fixed.
(B) Similar data within other protein-coding genes which exhibit statistically significant reduction in
H2Aub.
(C) Aggregate heatmap depicting normalized Ring1B CUT&RUN data within loci which DiffBind
determined had significant H2Aub reduction (FDR<0.05).
39
1240
1241
1242
1243
1244
1245
1246
1247
1248
1249
1250
1251
1252
1253
1254
1255
1256
1257
1258
1259
1260
1261
1262
1263
Figure 3. Neuronal loss of Bap1 primarily decreases H3K27ac, particularly in zones where H2Aub
increases.
(A) DiffBind differential analysis of H3K27ac abundance in CUT&RUN data from the P12 cerebellum,
as analyzed within consensus peak regions identified using SEACR (n=4 control, n=4 mutant).
(B) DiffBind differential analysis of H3K27ac abundance in CUT&RUN data from the adult cerebellum
(n=5 control, n=5 mutant).
(C) Genome-wide aggregate heatmaps depicting normalized H3K27ac CUT&RUN data from the adult
cerebellum within strong active enhancers, as identified using ChromHMM analysis of CUT&RUN
data from the adult cerebellum.
(D) ChIPseeker functional annotation of different classes of H3K27ac peak regions based upon how they
change in the cerebellum of adult Bap1 mutants.
(E) GIGGLE score comparisons of chromatin factors within peak regions exhibiting increased (left) or
decreased (right) H3K27ac in the adult cerebellum. The data source is indicated by the color key.
(F) GIGGLE score comparisons of histone modifications within peak regions exhibiting increased (left)
or decreased (right) H3K27ac. The data source is indicated by the color key.
(G) Normalized CUT&RUN detection of the indicated histone modifications in the adult control and
mutant at intergenic enhancers located upstream of the Mctp1 (left) or Rims1 genes. Shading showed
regions affected in mutant neurons.
(H) Log2-log2 scatterplot comparing the ratio of H3K27ac (mutant/control) to the ratio of H2Aub
(mutant/control) within H3K27ac peak regions identified in the adult cerebellum (n=6 for
H2AK119ub, n=5 for H3K27ac). Quantitation using Limma and EdgeR. Genome-wide
anticorrelation indicates that gain of H2Aub is associated with loss of H3K27ac and vice versa.
(I) Normalized CUT&RUN detection of the indicated histone modifications in adult neurons within a
large intergenic enhancer cluster. H2A.z is a H2A variant and H2A.zAc is its acetylated version.
40
1264
1265
1266
1267
1268
1269
1270
1271
1272
1273
1274
1275
1276
1277
1278
1279
1280
1281
1282
1283
1284
Figure 4. Loss of Bap1 in cerebellar granule neurons of adult Bap1f/f,Math1-cre mice bidirectionally
perturbs PRC2-dependent modifications.
(A) Diffbind differential analysis of H3K27me3 abundance in CUT&RUN data from the P12 cerebellum,
as analyzed within consensus peak regions identified using SEACR (n=4 control, n=4 mutant).
(B) Differential analysis of H3K27me3 abundance in CUT&RUN data from the adult cerebellum (n=4
control, n=4 mutant).
(C) ChIPseeker functional annotation of different classes of H3K27me3 peak regions, segregated by how
they are perturbed by Bap1 loss at each timepoint.
(D) Normalized CUT&RUN detection of the indicated histone modifications in the adult control and
mutant at the Lsamp gene shown alongside RNA-seq analysis of nuclear RNA.
(E) Transcripts per million (TPM) in RNA-seq analysis performed on nuclei isolated from the adult
cerebellum. FDR computed by DESeq2 (70).
(F) Log2-log2 scatterplot comparing the ratio of H3K27me3 (mutant/control) to the ratio of H2Aub
(mutant/control) within H3K27me3 peak regions identified in the late cerebellum. Quantitation using
Limma and EdgeR. Most points lie in the bottom left quadrant indicating loss of both H2Aub and
H3K27me3 in the Bap1 mutant. Weak positive correlation indicates that these modifications tend to
change together.
(G) Log2-log2 scatterplot comparing the ratio of H3K27me3 (mutant/control) to the ratio of H3K27me1
(mutant/control) within all protein coding genes. Quantitation using Limma and EdgeR. Anti-
correlation indicates that these PRC2-dependent modifications tend to change in opposite directions
when Bap1 is lost.
41
1285
1286
1287
1288
1289
1290
1291
1292
1293
1294
1295
1296
1297
1298
1299
1300
1301
1302
1303
1304
1305
1306
1307
1308
1309
1310
1311
1312
1313
1314
Figure 5. Neuronal loss of Bap1 makes heterochromatin resemble euchromatin and euchromatin
resemble heterochromatin.
(A) Normalized CUT&RUN detection of the indicated histone modifications in the adult control and
mutant at Zbtb20 gene.
(B) TPM and adjusted p-value for Zbtb20, as computed using DESeq2 on RNA-seq, analysis of nuclear
RNA from the adult cerebellum.
(C) DiffBind differential analysis of H3K27me3 abundance in CUT&RUN data from the adult cerebellum
analyzed within H3K27ac peak regions identified using SEACR.
(D) Genome-wide aggregate heatmaps depicting normalized H3K27me3 CUT&RUN data from the adult
cerebellum within H3K27ac peak summits, as detected using MACS2 narrow.
(E) Aggregate heatmaps depicting normalized H3K27me3 CUT&RUN data from the adult cerebellum
within ChromHMM-defined strong active enhancers.
(F) Aggregate heatmaps depicting normalized H3K27me3 CUT&RUN data from the adult cerebellum
within ChromHMM-defined active promoters.
(G) Normalized CUT&RUN detection of the indicated histone modifications in the typically
heterochromatic Myc gene.
(H) DiffBind differential analysis of H3K27ac abundance in CUT&RUN data from the adult cerebellum,
as analyzed within H3K27me3 peak regions identified using SEACR.
(I) Aggregate heatmaps depicting normalized H3K27ac CUT&RUN data from the adult cerebellum
within H3K27me3 peak summits, as detected using MACS2 narrow.
(J) Heatmaps depicting normalized H3K27ac CUT&RUN data from the adult cerebellum within
ChromHMM-defined bivalent promoters.
(K) Heatmaps depicting normalized H3K27ac, H3K27me3 and H2Aub CUT&RUN data from the adult
cerebellum within sites with increased chromatin accessibility in the mutant cerebellum, as
determined using ATAC-seq.
(L) Heatmaps depicting normalized H3K27ac, H3K27me3 and H2Aub CUT&RUN data from the adult
cerebellum within sites with decreased chromatin accessibility in the mutant cerebellum, as
determined using ATAC-seq.
(M) Permutation analysis of peak datasets for ATAC-seq data and the indicated histone modifications.
10,000 permutations were performed with the whole genome set as the universe.
42
1315
1316
1317
1318
1319
1320
1321
1322
1323
1324
1325
1326
1327
1328
1329
1330
1331
1332
1333
1334
Figure 6. Bap1 loss progressively perturbs the transcriptome of neuronal nuclei.
(A) PCA of RNA-seq data from the P12 and adult timepoints (n=7 adult controls, 7 adult mutants, 5
controls at P12, 5 mutants at P12). RNA was extracted from nuclei isolated from the cerebellum.
(B) Volcano plot of RNA-seq data from the P12 cerebellum. Adjusted p-values computed using DESeq2
(70). Significant transcripts had adjusted p-values <0.05 with fold-change >1.5.
(C) RNA-seq data from adult cerebellum. Subsets of these plots are shown to illustrate the overall
direction of changes.
(D) Reduced abundance of mRNA derived from the Syt1 gene, which loses activating histone
modifications (H3K27ac, H3K4me3, H3K27me1) while gaining repressive modifications (H2Aub,
H3K27me3).
(E) RNA-seq analysis of Syt1 at the adult timepoint.
(F) Increased abundance of mRNA derived from the Iqsec3 gene, which gains H3K27ac while losing
H2Aub and H3K27me3.
(G) RNA-seq analysis of Iqsec3 at the adult timepoint.
(H) Ratio of the abundance of the indicted histone modifications in mutant versus control neurons (at the
adult stage) within genes whose transcripts are increased or decreased by Bap1 loss.
(I) (J) Similar analysis at age P12.
Combined gene ontology analysis of all significantly upregulated or downregulated loci, as
determined using RNA-seq analysis of nuclei isolated from the cerebellum of adult mice.
43
1335
1336
1337
1338
1339
1340
1341
1342
1343
1344
1345
1346
1347
1348
1349
1350
1351
1352
1353
1354
1355
1356
1357
1358
1359
Figure 7. Loss of Bap1 in excitatory neurons of the forebrain of Bap1f/f,Emx1-cre mice perturbs
behavior and causes epileptiform neural activity.
(A) Photo of littermates.
(B) Body weight at the indicated ages (n= 96 controls, n=42 mutants, error bars SEM).
(C) Tracking data for individual mice during 15-minute open field test (n=8 control males, 9 control
females, 6 mutant males and 8 mutant females). For this and subsequent behavior analysis, mice
ranged in age from 5 to 9 months with a median age of 7 months.
(D) Number of entries to the border (10 cm wide perimeter) and center (central 40 x 40 cm square) of the
arena during the open field test. P-values by Kruskal-Wallis (K-W) test followed by Dunn’s multiple
comparisons test.
(E) Time spent in the border and center zones during the open field test. P-values by K-W test with
Dunn’s test.
(F) Average and maximum speed during the open field test.
(G) Number of marbles buried (out of 20 total) over 30 minutes. P-values by K-W test with Dunn’s test.
(H) Number of paw clasping episodes during 10 seconds of tail suspension. P-values by K-W test with
Dunn’s test.
(I) Number of hourly spike discharges over three days of continuous electroencephalography (EEG)
recordings (n=2 control males, 2 mutant males, 2 control females and 2 mutant females). EEG was
performed between 7 to 11 months with a median age of 9 months.
(J) Number of daily spike discharges during continuous EEG recordings.
(K) Frequency response curve of EEG power density during different vigilance states. Traces represent
mean power (µV²/Hz) ± SEM averaged across animals within each state (n=4 controls, 4 mutants).
(L) Power spectral density of EEG recordings was computed using fast Fourier transform across vigilance
states. Canonical frequency bands are indicated as follows: delta (0.5–4 Hz), theta (4–8 Hz), alpha (8–
13 Hz), beta (13–30 Hz), and gamma (>30 Hz) (n=4 controls, 4 mutants).
44
1360
1361
1362
1363
1364
1365
1366
1367
1368
1369
1370
1371
1372
1373
1374
1375
1376
1377
1378
1379
1380
1381
1382
1383
1384
1385
1386
1387
Figure 8. Hyperexcitability and increased basal synaptic transmission of hippocampal CA1
pyramidal neurons from Bap1f/f,Emx1-cre mice.
(A) Resting membraine potential (VM) in hippocampal pyramidal neurons in acutely prepared slices from
mice at age one month (n=20 control neurons, 23 mutant neurons, p-value by unpaired t-test)
(B) Representative traces of resting membrane potential of control and mutant hippocampal neurons, the
latter of which fires spontaneous action potentials.
(C) Fraction of neurons exhibiting spontaneous action potentials (APs).
(D) Representative rheobase recordings, with the current injection protocol shown on the left and the
output on the right.
(E) Rheobase currents (n=21 control neurons, 25 mutant neurons, p-value by unpaired t-test)
(F) Representative recordings of induced APs. Left: input current steps from -100 pA to 500 pA in 50 pA
increments, 50 ms duration. Right: Resulting voltage sweeps.
(G) Number of APs generated for a given current injection (n=20 control neurons, 25 mutant neurons,
*indicates p-value <0.05 as computed by Mixed-effects modeling with Geisser-Greenhouse
correction followed by Fisher’s least significance difference [LSD] test).
(H) Frequency of APs in response to the current injection protocol in (G).
(I) Amplitude of induced APs in response to the current injection protocol in (G).
(J) Representative input-output recordings of fiber volley amplitude and field excitatory postsynaptic
potential (fEPSP).
(K) Fiber volley amplitude (n=16 control, 15 mutant slices, not significant by mixed-effects modeling
with Geisser-Greenhouse correction followed by Fisher’s LSD test).
(L) fEPSP amplitude (n=16 control, 15 mutant slices, not significant by mixed-effects modeling with
Geisser-Greenhouse correction followed by Fisher’s LSD test).
(M) fEPSP slope (n=16 control, 15 mutant slices, *indicates p-value <0.02 as computed by mixed-effects
modeling with Geisser-Greenhouse correction followed by Fisher’s LSD test).
(N) Representative paired-pulse facilitation (PPF) traces: two consecutive fEPSPs evoked by Schaffer
Collateral stimulation with a 50 ms interstimulus interval.
(O) Summary of paired-pulse ratios (n=16 control, 15 mutant slices, p-value by unpaired t-test).
45
1388
1389
1390
1391
1392
1393
1394
1395
1396
1397
1398
1399
1400
1401
1402
1403
1404
1405
1406
1407
1408
Supplementary Materials for
The histone deubiquitinase Bap1 maintains the distinction between euchromatin and
heterochromatin during neurodevelopment
Wang et al.
*Corresponding author. Email: cferguson@health.ucsd.edu
This PDF file includes:
Figs. S1 to S7
Table S1
46
1409
1410
1411
1412
1413
1414
1415
1416
1417
1418
1419
1420
1421
1422
1423
1424
1425
1426
SUPPLEMENTAL FIGURE LEGENDS
Figure S1: Bap1 loss perturbs H2Aub and cerebellar neurodevelopment.
(A) Mean latency to fall in the accelerating Rota-Rod test of cerebellar-dependent motor performance and
learning. Daily trials were averaged among controls and mutants. Animals ranged in age from 4 to 7
months. Error bars SEM, p-values by two-way ANOVA with Sidak test.
(B) Abundance of transcripts from the PR-DUB complex, as determined by RNA-seq analysis of nuclei
isolated from the cerebellum at the indicated ages. The transcripts per million (TPM) metric was
computed using DESeq2 (70).
(C) Immunoblot detection of Bap1, H2Aub and H3K27ac in nuclei isolated from the mouse cerebellum at
the indicated ages. Lysates were prepared from nuclei subjected to benzonase digestion to include
chromatin.
(D) Principal Component Analysis (PCA) of H2Aub CUT&RUN data from replicates of the early and
adult cerebellum.
(E) Minus-average plot depicting H2Aub CUT&RUN data from nuclei isolated from the developing
cerebellum (age P12) of control and Bap1f/f,Math1-cre mutant mice, as analyzed using DiffBind (32).
(F) Minus-average plot of H2Aub CUT&RUN data from the adult cerebellum.
47
1427
1428
1429
1430
1431
1432
1433
1434
1435
1436
1437
1438
Figure S2: Bap1 loss perturbs H3K27ac in the mouse cerebellum.
(A) Principal Component Analysis (PCA) of H3K27ac CUT&RUN data from the early or late
cerebellum.
(B) Minus-average plot depicting H3K27ac CUT&RUN data from nuclei isolated from developing
cerebellum (age P12) of control and Bap1f/f,Math1-cre mutant mice, as analyzed using DiffBind (32).
(C) Minus-average plot of H3K27ac CUT&RUN data from the adult (3 month) cerebellum .
(D) ChIPseeker (33) functional annotation of different classes of H3K27ac peak regions based upon how
they behave in Bap1 mutant neurons.
(E) GIGGLE score comparisons between peak regions with increased H3K27ac in the developing
cerebellum and the Cistrome (34, 35) database of chromatin factors (left) and modified histones
(right). The source of experimental data is indicated by the color key.
(F) Comparisons between peak regions with decreased H3K27ac and the Cistrome database.
48
1439
1440
1441
1442
1443
1444
1445
Figure S3: Bap1 loss perturbs H3K27me3 in the developing cerebellum.
(A) Principal Component Analysis (PCA) of H3K27me3 CUT&RUN data from the P12 and adult
cerebellum.
(B) Minus-average plot depicting H3K27me3 CUT&RUN data from nuclei isolated from the P12
cerebellum of control and Bap1f/f,Math1-cre mutant mice, as analyzed using DiffBind.
(C) Minus-average plot depicting H3K27me3 CUT&RUN data from nuclei isolated from the adult
cerebellum.
49
1446
1447
1448
1449
1450
1451
1452
1453
1454
Figure S4: Bap1 loss bidirectionally perturbs PRC2-dependent H3K27me1.
(A) Principal Component Analysis (PCA) of replicate H3K27me1 CUT&RUN data from the adult
cerebellum.
(B) Volcano plot depicting H3K27me1 abundance within protein-coding genes, as analyzed using
DiffBind.
(C) Minus-average plot of H3K27me1 within protein-coding genes.
(D) Three genes with significant increases in H3K27me1. These genes essentially lack any significant
amounts of the other modifications.
50
1455
1456
1457
1458
1459
1460
1461
1462
1463
1464
1465
1466
1467
1468
1469
1470
1471
1472
1473
1474
Figure S5: Bap1 loss bidirectionally perturbs chromatin accessibility in ways that complement
changes in histone modifications.
(A) Principal Component Analysis (PCA) of ATAC-seq data from the adult cerebellum (n=6 control, 6
mutant).
(B) Volcano plot depicting ATAC-seq data from nuclei isolated from the adult (3 month), as quantified
using DiffBind, within peak regions identified by MACS2 (38) narrow.
(C) Minus-average plot.
(D) ChIPseeker (33) functional annotation of different classes of ATAC-seq peak regions based upon how
they behave in Bap1 mutant neurons.
(E) GIGGLE score comparisons between regions with increased accessibility and the Cistrome (34, 35)
database of chromatin factors (left) and modified histones (right). The source of experimental data is
indicated by the color key.
(F) Comparisons between regions with decreased accessibility and the Cistrome database.
(G) Volcano plot depicting ATAC-seq data within H3K27ac peak regions from control cerebellum as
computed using DiffBind.
(H) Locus with multiple sites where reduced accessibility is accompanied by reductions in H3K27ac.
(I) Volcano plot depicting ATAC-seq data within H3K27me3 peak regions from control cerebellum as
computed using DiffBind.
(J) Locus with multiple sites where increased accessibility is accompanied by reduced
heterochromatinization (reduced H2Aub and H3K27me3, increased H3K27ac).
51
1475
1476
1477
1478
1479
1480
1481
1482
1483
1484
Figure S6: Bap1 loss promotes deposition of H3K27me3 within euchromatin and H3K27ac within
heterochromatin.
(A) PCA of H3K27me3 CUT&RUN data evaluated within H3K27ac peak regions, as detected in the adult
control cerebellum.
(B) Minus-average plot depicting H3K27me3 CUT&RUN data from nuclei isolated from the adult
cerebellum as analyzed using DiffBind specifically within H3K27ac peak regions.
(C) PCA of H3K27ac CUT&RUN data evaluated within H3K27me3 peak regions, as detected in the adult
control cerebellum.
(D) Minus-average plot depicting H3K27ac CUT&RUN data from nuclei isolated from the adult
cerebellum as analyzed using DiffBind specifically within H3K27me3 peak regions.
52
1485
1486
1487
1488
1489
1490
1491
1492
1493
1494
1495
1496
1497
1498
1499
1500
1501
1502
1503
1504
1505
1506
1507
1508
1509
1510
1511
1512
1513
1514
1515
1516
1517
1518
Figure S7: Bap1 loss perturbs the abundance of genes required for neuronal connectivity, synaptic
communication and identity. Significantly dysregulated transcripts encoding cell adhesion molecules.
Transcripts per million (TPM) and adjusted p-values derived from RNA-seq experiments using cerebellar
nuclei isolated from the adult cerebellum. Data were analyzed using DESeq2 (70). Many dysregulated
transcripts fell into shared categories.
(A) Genes which are essential for axonal targeting and therefore the establishment and maintenance of
neuronal connectivity. Cadherins are Ca²⁺-dependent adhesion proteins that establish and maintain
precise synaptic connections. Protocadherins provide molecular diversity for neuron–neuron
recognition, specifying synaptic connectivity by enabling neurons to distinguish self from non-self.
Contactins regulate axon guidance, synapse formation, and axon–glia interactions. Collagen proteins
constitute part of the perisynaptic extracellular matrix that stabilize and modulate pre- and
postsynaptic interactions. Ephrins are membrane-bound signaling proteins that interact with Eph
receptors to control axon guidance and synapse development by mediating contact-dependent
repulsive and attractive cues. Lama5 and Lamb3 encode extracellular matrix laminin components that
assemble into heterotrimeric basement membrane proteins that regulate migration, axon guidance,
and synapse stabilization by interacting with the perisynaptic matrix. Immunoglobulin superfamily
cell adhesion molecules synaptic recognition to specify synaptic partners and stabilize synaptic
contacts. Robo1 encodes a guidance receptor that regulates axon pathfinding and synaptic circuit
wiring by mediating repulsive responses to Slit ligands. Semaphorins are secreted and membrane-
associated guidance cues that shape synapse formation and plasticity by regulating axon targeting,
dendritic patterning, and contact-dependent repulsive or stabilizing signaling at synapses.
(B) Dysregulated transcripts encoding transcriptional regulators and key transcription factors required for
neuron-specific transcriptomic regulation. Zic proteins are key regulators of cerebellar development
that control granule neuron progenitor proliferation, differentiation, migration, and maintenance of
identity.
(C) Dysregulated transcripts encoding regulators of synaptic transmission. Synaptotagmins act as the
principal Ca²⁺ sensors for synaptic vesicle fusion, triggering neurotransmitter release by coupling
calcium influx to the activation of the SNARE-mediated membrane fusion machinery.
(D) Dysregulated transcripts regulating membrane calcium flux. Cacna1B encodes the α1B subunit of an
N-type voltage-gated calcium channel that mediates rapid, activity-dependent Ca²⁺ influx at
presynaptic terminals to trigger neurotransmitter release. Cacna1c encodes a voltage-gated L-type
calcium channel subunit that mediates sustained Ca²⁺ influx in response to strong depolarization,
linking neuronal activity to gene expression and plasticity. Cacna1i encodes a T-type voltage-gated
calcium channel subunit that mediates low-threshold, transient Ca²⁺ currents important for burst
53
1519
1520
1521
1522
1523
1524
1525
1526
1527
1528
1529
1530
1531
1532
1533
1534
1535
firing, rhythmic activity, and subthreshold excitability. Cacnb1 encodes the β1 auxiliary subunit of
voltage-gated calcium channels, which regulates channel trafficking to the membrane, modulates
channel gating properties, and fine-tunes the amplitude and kinetics of calcium influx in neurons.
Cacng8 encodes an auxiliary transmembrane protein that modulates AMPA receptor trafficking,
synaptic targeting, and gating properties to shape excitatory synaptic transmission and plasticity.
Camk4 encodes a nuclear Ca²⁺/calmodulin-dependent kinase that decodes neuronal activity into
changes in gene transcription, thereby regulating synaptic plasticity, learning-related gene programs,
and long-term neuronal adaptations. Grik4 and Grik5 encode subunits of ionotropic kainate-type
glutamate receptors that assemble into heteromeric channels mediating excitatory synaptic
transmission and modulating synaptic strength, neuronal excitability, and presynaptic
neurotransmitter release. Itpr2 and Itpr3 encode ER membrane calcium-release channels that open in
response to IP₃ to mobilize intracellular Ca²⁺ stores. Mcu encodes the pore-forming subunit of the
mitochondrial calcium uniporter that enables Ca²⁺ entry into the mitochondrial matrix. Micu1 encodes
a Ca²⁺-sensing regulatory protein that sets a threshold for Mcu opening, preventing mitochondrial
calcium uptake at low cytosolic Ca²⁺ and allowing it during strong neuronal activity. Ryr2 and Ryr3
encode endoplasmic reticulum calcium-release channels that mediate Ca²⁺-induced Ca²⁺ release in
neurons, amplifying intracellular calcium signals.
54
1536
Table S1: Reagents used in this study
REAGENT or RESOURCE SOURCE Antibodies
H2AK119ub CST H3K27me3 CST H3K27ac Active Motif H3K4me3 CST H3K27me1 Active Motif Ring1B CST H2A.z Active Motif H2A.zAc (Lys4/Lys7) CST Chemicals, Peptides, and Recombinant Proteins
Nuclei EZ Lysis Buffer Sigma-Aldrich Protease and Phosphatase
inhibitor cocktail
Thermofischer IDENTIFIER
CAT# 8240S
RRID: AB_10891618
CAT# 9733S
RRID: AB_2616029
CAT# 39133
RRID: AB_2561016
CAT# 9751S
RRID: AB_2616028
CAT# 61016
RRID: AB_2715573
CAT# 5694S
RRID: AB_10705604
CAT# 39013
RRID: AB_2615081
CAT# 75336
RRID: AB_2799867
CAT# N3408-200ML
CAT# 1891280
55
Spermidine Sigma-Aldrich Sodium Butyrate Sigma-Aldrich Complete Mini, EDTA-free Roche Iodoacetamide Millipore-Sigma Digitonin Sigma-Aldrich pAG-MNase Epicypher Trizol Invitrogen RNase A Promega Glycogen Invitrogen EDTA VWR EGTA EMD Millipore RNasin® Ribonuclease
Inhibitor
Promega DMSO Sigma-Aldrich Chloroform Fischer Chemical Deposited Data
Raw data and analyzed data This paper Experimental Models: Organisms/Strains
Bap1tm2c(EUCOMM)Hmgu/BrugJ, The Jackson CAT# 85558-1G
CAT# 19-137
CAT# 11836170001
CAT# A3221-1VL
CAT# D141
CAT# 15-1016
CAT# 15596018
CAT# A7973
CAT# 10814010
CAT# 0322
CAT# 324626
CAT# N261B
CAT# D2650-100ML
CAT# C606-1
NCBI BioProject PRJNA1418631
Cat #:031874
56
Bap1 f/f B6.Cg-Tg(Atoh1-cre)1Bfri/J,
Math1/Atoh1-cre mice
B6.129S2-
Emx1t m 1 ( cre) K r j /J, cre , Emx1-Cre mice
Emx1IRES
C57BL/6J Software and Algorithms
DiffBind chIPSeeker SEACR MACS2 SAMtools v.1.14 BEDTools ShinyGO ChromHMM KaryoploteR v.1.22.0 deepTools STAR Laboratory RRID:IMSR_JAX:031874
The Jackson
Laboratory
Cat# 011104
RRID:IMSR_JAX:011104
The Jackson
Laboratory
Cat #:005628
RRID:IMSR_JAX:005628
The Jackson
Laboratory
Cat# 000664
RRID:IMSR_JAX:000664
(32)
(33)
(31)
(38)
(66)
(80)
(72)
(39)
(68)
(67)
(71)
57
Limma EdgeR HOMER clusterProfiler MultiClamp 700B Commander Clampex 11 RegioneReloaded DESeq2 GraphPad Prism10
Adobe Illustrator 2021
Equipment
Vibratome Micropipette puller Temperature controller MultiClamp 700B amplifier Accelerated Rota-Rod (Jones
and Roberts) For Mice 7650
Axon Digidata 1550B Other
(37)
(36)
(81)
(82)
Molecular Devices
Molecular Devices
(83)
(70)
Leica, VT1000 S
Sutter, Model P-97
Warner Instrument,
TC-324B
Molecular Devices
Ugo Basile
Molecular Devices
58
Concanavalin A-coated
Epicypher magnetic beads
RNeasy Mini Kit Qiagen NEBNext Ultra II DNA
New England Biolabs Library Prep Kit
Illumina Stranded Total RNA
Illumina Prep Ligation with Ribo-Zero
Plus Kit
Qubit™ dsDNA Quantification
Invitrogen Assay Kits
Qubit™ RNA Broad Range
Invitrogen (BR) Assay Kit
40 µM filters Corning ATAC-seq kit Active Motif CAT# 21-1401
CAT# 74104
CAT# E7645L
CAT# 20040525
CAT# Q32854
CAT# Q10210
CAT# 431750
CAT# 53150
1537
59