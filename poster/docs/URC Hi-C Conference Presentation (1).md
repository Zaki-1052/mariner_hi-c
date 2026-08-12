By: Zakir Alibhai
Regulation of chromatin conformation by the histone deubiquitinase BAP1 in the brain
Cole Ferguson Lab
Neuropathology
Undergraduate Research Conference
Hi, my name is zakir, im a 2nd year bioinformatics researcher in the cole ferguson lab, and my talk will be on the 
15 seconds

Chromatin folding
Long-range contacts

The 3D Genome
Chromatin folding
Long-range contacts
Gene regulation
Hi-C detection
Looping:
So, before we talk about the 3D genome, some quick background biology. Click.
Now: you may remember we have the loong dna helix here on the right, but it needs to fit inside the nucleus of a cell, so it wraps around spools called histones, and one bead of that is called a nucleosome. Then you take that string and coil it together into chromatin, which then makes up the classic chromosome. click.
So with that brief review, I’ll now start out by explaining the first term from my title, click - 
which deals with a misconception about dna: the genome as we saw isn't just a linear sequence or a flat string. DNA is packaged into chromatin and folded into specific 3D shapes inside the nucleus. Click. This is true of every cell, and we can measure the structures that form at different scales, as can be seen here in this diagram. Click. 
And this folding brings regions that are really far apart on the DNA strand into physical contact. Click. Click. Click.These contacts are called loops, which I will mainly be focusing on. [fig] click. Click. Click. 
Now, these contacts have an actual purpose. Enhancers, which are gene switches usually far away from the genes they control, have to physically touch their target promoters to turn genes on. If you disrupt the folding, you disrupt which genes are on or off. click.
Hi-C is the technique that lets us measure these contacts genome-wide. [Point to figure] On top is the contact matrix, where the bright spots off the diagonal show where two regions are touching. And to the side, we can see a physical loop structure with the essential elements.
2:00

Histone modifications
BAP1 (deubiquitinase) 
Neurodevelopment
Mouse cerebellum
A Tumor Suppressor in the Brain
BAP1 removes the epigenetic modification H2AK119ub
Parmar A, Srinivasan A, Krockenberger L, Augustine A, Gong O, et al. (2025) Polycomb repressive complexes 1 and 2 independently and dynamically regulate euchromatin during cerebellar neurodevelopment. PLOS Genetics 21(9): e1011843. https://doi.org/10.1371/journal.pgen.1011843


Click. Now, chromatin folding, or conformation, is shaped by chemical tags on histones, which are the proteins that DNA wraps around. Click. Click. One of those tags is called H2AK119ub: a ubiquitin molecule attached to histone H2A. click. 
This is traditionally understood as a gene-silencing mark, but recent studies from my lab have shown otherwise. Click. 
It's placed on chromatin by a protein repressive complex called Polycomb, or PRC1 here. [diagram] Click.
BAP1 is the enzyme that erases this histone modification. It's called a deubiquitinase, because it literally removes a ubiquitin, which is usually repressive. Going forward, when I say BAP1, I mean this enzyme, and I’ll be calling H2AK119ub ubiquitinated histone. [show erasing direction] Click.
BAP1 was originally studied as a tumor suppressor, where mutations were linked to cancer. But when BAP1 is lost in the brain, both patients and mouse models show severe developmental delay, epilepsy, and motor and speech problems. Click. 
Despite these stakes, BAP1's molecular role in the brain is still relatively new territory. That's where my lab comes in: specifically, we knock out BAP1 in the mouse cerebellum, the part of the brain that coordinates movement, to study its importance.
3:00

Does BAP1 loss reshape 3D folding in developing neurons?
BAP1-KO
Epigenetic Dysregulation
???
Click. So, our lab has previously characterized what BAP1 loss does to the epigenome – and as you'd expect, ubiquitinated histone piles up when BAP1 is gone. 
But we’ll see this is just the start of a cascade of regulatory changes. Click.
What hadn't been asked, and what we wanted to test, is whether all of this epigenetic disruption also reshapes the physical folding of the genome. Click.
So: does BAP1 loss reshape 3D folding in the brain? That’s our ultimate question
3:30
H3K27me3 redistributes, active marks

Asking the 3D Question
Isolate nuclei from cerebella


Sequence DNA








3 mice per condition
~2 mo. old
Genotype comparison
Two timepoints
Hi-C profiling
Computational pipelines

Hi-C 
Contact
Matrices:



Durand, N. C., Robinson, J. T., Shamim, M. S., Machol, I., Mesirov, J. P., Lander, E. S., & Aiden, E. L. (2016). Juicebox provides a visualization system for Hi-C contact maps with unlimited zoom. Cell Systems, 3(1), 99–101. https://doi.org/10.1016/j.cels.2015.07.012


Click. So to test this… click.
We compared control vs. mutant mice, 3 mice per group, at two ages: early postnatal and adult. click.
We used Hi-C, or chromatin conformation capture, which is the method described earlier. Click. Click. 
Here you can see both our process and what we call “contact maps”. Click. Click. Click. 
These are the different scales at which we can measure how close regions of the DNA are touching, and how much. Click.
As a bioinformatics student, my role was computational: I built our differential analysis pipelines you’ll see the results of. Click.
Initially, we just wanted to tell at a broad resolution whether those compartments you saw earlier changed. Click. 
But the changes at a smaller scale - what I started on - ended up being even more interesting. 
4:00


Zooming In on Loops


Loops = regulation
Strength change
Gained / lost
Genome-wide
Loops control regulatory chromatin interactions
Click. Now, That smaller scale is loops: [click] we established that these are where distant regulatory elements physically meet. click.
They’re basically our measure of gene regulation (whether it is “off” or “on”) that we care about. click.
So at each detected loop, we asked: is the contact strength different between the control (or wiltype) and mutant (our knockout model)? click.
We call a loop that's weaker in the mutant 'lost.' And Stronger in the mutant is 'gained.' I’ll keep using these terms next slide. click.
When we looked across the genome, thousands of these loops were different, and in specifically meaningful ways.
4:30-45

Thousands of Loops Dysregulated


12-step pipeline
Mariner + edgeR
Replicates/resolutions
~2,800 differential loops
Click. So now I'll focus on our findings - ultimately how loops change, at what distances, and what that tells us. click.
The first question just asked what this loop dysregulation actually looked like. So click.
I built a 12-step pipeline in R, using the mariner Hi-C package and the edgeR framework - click. identifying the loops, aggregating them in binned matrices, and calculating their differential strength. click.
With three mice per condition, and over a billion “reads” or DNA fragments, I could use statistics to detect really small differences.
That code is what generated the volcano plot on this slide. Now, to read this plot, we can see that each dot is a loop, classified by their effect size and importance. The x-axis is the change in strength between the two conditions, and the y-axis is the false discovery rate, or how “real” the loop is. The top left are the ones that weakened the most, and the  top right are the ones that strengthened the most. click.
Now, we see that in the adult cerebellum, almost 3000 loops were significantly differential, or up to one in five of all detected loops at the largest resolution.
And they split into two groups: about 1,200 lost — weaker in the mutant — and about 1,600 gained — stronger in the mutant. 
So to recap, this means when we knock out bap1, a great number of contacts change in our mutant, possibly related to neuronal dysregulation.
Then, the next obvious step was to look at what kinds of loops were on either side of this very uneven split.
5:30-6

Preferential Loss of Long Loops


Length divide
Lost = long
Gained = short
Switched anchors
Click. The clearest pattern came from looking at the distances between the anchors, of each loop. In other words, how far apart either end sat on the DNA. click.
This first plot shows these lengths for our lost and gained loops. click.
The interesting part is the massive separation between the two curves, where the red curve of lost loops sits at longer distances, and the blue curve, where we gain loops, shifts shorter. click.
These long range contacts are about three times more likely to be lost than gaind, where the median is about switched in half. click.
From the volcano plot on the last slide, we know the genome isn't losing contacts overall, and it also isn't gaining them overall. click.
We believe from this plot [click] that It's trading long-range contacts for shorter-range ones.
To test this, I looked at regions where the same anchor participated in both a lost loop and a gained loop.
This shows over 200 loci, or genes, where this phenomenon occurred, and again, while there are more gained overall, lost loops are longer.
So we have a pattern, and it’s happening consistently in the same place. The next question to ask is: what kind of chromatin is driving it?
7:00

Polycomb Drives the Asymmetry


Loop directions
Polycomb-enriched
Structural role
So here [click] I've taken all the differential loops and split them by direction — strengthened in the mutant on the left, weakened on the right. Click.
Each bar is colored by the chromatin state, or type of loop, at its anchors, using histone modifications to classify the regions. click.
This is essentially asking: when a loop goes up versus when a loop goes down, what kind of chromatin is at its anchor? click.
The green here is Polycomb, which is repressed. About 22 percent of our gained loops vs about 5 percent of the lost. 
On the other side — the blues and reds are active enhancers and active promoters. Together they make up about 40 percent of the lost loops, and less than ten opposite. 
So gains happen at Polycomb, and losses happen at active regulatory elements. click.
From this, we believe that BAP1, through adjusting repressive marks, is affecting which regions contact which. When it’s gone, Polycomb regions collapse inward onto short-range contacts — the gained. Active regions lose their distant partners — the lost. 
If we can link all of these structural events to the same enzyme, where the gap between distances gets much larger when we look at repressive polycomb regions, well. That would be a complete picture.
8:00

BAP1 Loss Rewires the 3D Genome


H2AK119ub accumulation
Dynamic cycling
Progressive damage
18x
Click. Here’s the mechanistic part I want to connect, [click] linking how much ubiquitinated histone is changing at the anchors, between the lost and gained loops. 
Click. On the left, for lost loops, ubiquitinated histone is significantly elevated in the mutant with a positive median. 
At gained loops, it’s depleted. Click. So the loops that die are specifically the ones at anchors where it’s piling up. click.
Click. As we remember from that first molecule slide, PRC1 writes ubiquitinated histone, and BAP1 erases it. Click. 
There’s some complex regulatory cycle we don’t fully understand yet at play here, but we know that when its mark gets stuck at a loop region, it fails.
Click. Lastly, we have a neurodevelopmental component, where at an early timepoint when the mice are only 13 days old, we see a couple hundred dysregulated loops. But by adulthood, that increases [click] almost 20 times, to the effects ive been talking about. The failures we see when the mice are young start at polycomb sites, and spread from there. click.
So what starts as one molecule accumulating completely reorganizes the 3D structure of the genome over development.
And this is how we’ve established a mechanistic basis for deficits in neurodevlopment when bap1 is perturbed.
11:00


Thank You!
Any questions?
Acknowledgements
PI: Dr. Cole Ferguson
Graduate Mentor:
Challana Tea
Funding:
RO1 (NINDS) 
NIH Grant
Thank you!
11:30
