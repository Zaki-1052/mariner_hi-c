[3:35 PM]5 files [3:35 PM]just to be clear these are absolute levels
[3:36 PM]and yeah neuronal genes are constituitively enriched in ctrl,
[3:37 PM]but the enrichment is weakening in the mut? the gained genes are only 15.4% neuronal vs 23.5% genome-wide and no sig terms
[3:39 PM]remember this was the fold change from earlier btw
2 files [3:40 PM]49/49 in this latest section for gsea terms with positive nes (the ctrl having k119ub neuronal genes)
Jai Ramchandra  [3:41 PM]
i think that makes sense. like in the ctrl, ub controls more neuronal genes, but then in the mutant, those sites are prob remodeled more to heterochrom (we saw increase in K27me3 and decrease in K27ac there). so mecp2 is going there. but overall, Ub change by bap1 is broad and not specific to just neuronal genes. it just seems that maybe they are the ones that attract mecp2 the most.
Zakir Alibhai  [3:42 PM]
(also in the logs from R it suggested to add scoreType = "pos", since all levels are above 1 which will improve power...u can assume p vals will get a bit better but ill say if it changes significantly)
Jai Ramchandra  [3:42 PM]
ahhh okok
[3:43 PM]i guess to confirm the theory, we would need to overlap the neuronal sites in both ctrl and mut with the ATAC data? and then see if the neuronal genes are disproportionally shifting to heterochrom?
[3:44 PM]i feel like the answer is prob yes since we saw that the neuronal genes that mecp2 bound had the decrease in K27ac and increase in K27me3