library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(rtracklayer)

TP53.region <- toGRanges("chr3:36,664,677-38,372,860")
kp <- plotKaryotype(genome = "hg38", zoom = TP53.region,cex=1)

genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
genes.data <- mergeTranscripts(genes.data)
genes.data <- addGeneNames(genes.data,org.Hs.eg.db)
kpPlotGenes(kp, data=genes.data,r0=0, r1=0.15, gene.name.cex = 0.5)

source("16_Track_BigWig.R")








#BIGWIG + BED FILE 1
at <- autotrack(current.track = 2, total.tracks = 2, r0=0.35, r1=1)
pclgRNA3n1H3K27ac.bw<-"pclgRNA3n1H3K27ac.bw"
kpPlotBigWig(kp, data="pclgRNA3n1H3K27ac.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="limegreen", border = NA)
kpAddLabels(kp, labels = "pclgRNA3n1H3K27ac", r0=at$r0, r1=at$r0, cex=0.5, label.margin = 0.015,col="limegreen")

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K27me1P107ctrl1.bw- "H2AK119ubP12ctrl1.bed"
kpPlotRegions(kp, data=H2AK119ubP12ctrl1.bed, col="dodgerblue4", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 2
at <- autotrack(current.track = 1, total.tracks = 2, r0=0.35, r1=1)
pclgRNA3n2H3K27me3.bw<-"pclgRNA3n2H3K27me3.bw"
kpPlotBigWig(kp, data="pclgRNA3n2H3K27me3.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="tomato3", border = NA)
kpAddLabels(kp, labels = "pclgRNA3n2H3K27me3", r0=at$r0, r1=at$r1, cex=0.5, label.margin = 0.015,col="tomato3")

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K27me3ctrl1peaks.bed<- "H3K27me3ctrl1peaks.bed"
kpPlotRegions(kp, data=H3K27me3ctrl1peaks.bed,col="tomato3", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 3
at <- autotrack(current.track = 3, total.tracks = 5, r0=0.35, r1=1)
H2AK119ubCerebellumLate1.bw<-"H2AK119ubCerebellumLate1.bw"
kpPlotBigWig(kp, data="H2AK119ubCerebellumLate1.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="tomato3", border = NA)
kpAddLabels(kp, labels = "H2AK119ubCerebellumLate1", cex=0.5, label.margin = 0.015,col="tomato3",r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K27acCtrl2peaks.bed <- "H3K27acCtrl2peaks.bed"
kpPlotRegions(kp, data=H3K27acCtrl2peaks.bed, col="limegreen", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 4
at <- autotrack(current.track = 2, total.tracks = 5, r0=0.35, r1=1)
H2AK119ubCerebellumLate2.bw<- "H2AK119ubCerebellumLate2.bw"
kpPlotBigWig(kp, data="H2AK119ubCerebellumLate2.bw", ymin=0, ymax=100, r0=at$r0, r1=at$r1,col="tomato3", border = NA)
kpAddLabels(kp, labels = "H2AK119ubCerebellumLate2", cex=0.5, label.margin = 0.015,col="tomato3",r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 1, total.tracks = 5, r0=0.35, r1=1)
CpGmm10.bed<- "CpGmm10.bed"
kpPlotRegions(kp, data=CpGmm10.bed, col="limegreen", r0=at$r0, r1=at$r1)
kpAddLabels(kp, labels = "CpGmm10", cex=0.5, label.margin = 0.015,col="limegreen",r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 5
at <- autotrack(current.track = 1, total.tracks = 5, r0=0.35, r1=1)
H3K27ac.bw <- "H3K27ac.bw"
kpPlotBigWig(kp, data="H3K27ac.bw", ymin=0, ymax=60, r0=at$r0, r1=at$r1,col="darkgreen", border = NA)
kpAddLabels(kp, labels = "H3K27ac", cex=0.5, label.margin = 0.015,col="darkgreen",r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K27acCtrl2peaks.bed<- "H3K27acCtrl2peaks.bed"
kpPlotRegions(kp, data=H3K27acCtrl2peaks.bed, col="darkgreen", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 6
at <- autotrack(current.track = 11, total.tracks = 16, r0=0.35, r1=1)
sorted_230927_i_12_BAP1_P56_mut_H3K4me3_S14_aligned_reads.bam_norm_rnorm.bw <- "sorted_230927_i_12_BAP1_P56_mut_H3K4me3_S14_aligned_reads.bam_norm_rnorm.bw"
kpPlotBigWig(kp, data="sorted_230927_i_12_BAP1_P56_mut_H3K4me3_S14_aligned_reads.bam_norm_rnorm.bw", ymin=0, ymax=300, r0=at$r0, r1=at$r1,col="red", border = NA)
kpAddLabels(kp, labels = "sorted_230927_i_12_BAP1_P56_mut_H3K4me3_S14_aligned_reads.bam_norm_rnorm", cex=0.5, label.margin = 0.015,col="red",r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K4me1ctrl1peaks.bed<- "H3K4me1ctrl1peaks.bed"
kpPlotRegions(kp, data=H3K4me1ctrl1peaks.bed, col="purple", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 7
at <- autotrack(current.track = 10, total.tracks = 16, r0=0.35, r1=1)
sorted_230927_i_7_BAP1_P56_ctrl_H3K27ac_S9_aligned_reads.bam_norm_rnorm.bw <- "sorted_230927_i_7_BAP1_P56_ctrl_H3K27ac_S9_aligned_reads.bam_norm_rnorm.bw"
kpPlotBigWig(kp, data="sorted_230927_i_7_BAP1_P56_ctrl_H3K27ac_S9_aligned_reads.bam_norm_rnorm.bw", ymin=0, ymax=150, r0=at$r0, r1=at$r1,col="navy", border = NA)
kpAddLabels(kp, labels = "sorted_230927_i_7_BAP1_P56_ctrl_H3K27ac_S9_aligned_reads.bam_norm_rnorm", cex=0.5, label.margin = 0.015,col="navy",r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K4me1ctrl1peaks.bed<- "H3K4me1ctrl1peaks.bed"
kpPlotRegions(kp, data=H3K4me1ctrl1peaks.bed, col="purple", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 8
at <- autotrack(current.track = 9, total.tracks = 16, r0=0.35, r1=1)
sorted_230927_i_8_BAP1_P56_mut_H3K27ac_S10_aligned_reads.bam_norm_rnorm.bw <- "sorted_230927_i_8_BAP1_P56_mut_H3K27ac_S10_aligned_reads.bam_norm_rnorm.bw"
kpPlotBigWig(kp, data="sorted_230927_i_8_BAP1_P56_mut_H3K27ac_S10_aligned_reads.bam_norm_rnorm.bw", ymin=0, ymax=150, r0=at$r0, r1=at$r1,col="red", border = NA)
kpAddLabels(kp, labels = "sorted_230927_i_8_BAP1_P56_mut_H3K27ac_S10_aligned_reads.bam_norm_rnorm", cex=0.5, label.margin = 0.015,col="red",r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K4me1ctrl1peaks.bed<- "H3K4me1ctrl1peaks.bed"
kpPlotRegions(kp, data=H3K4me1ctrl1peaks.bed, col="purple", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 9
at <- autotrack(current.track = 8, total.tracks = 16, r0=0.35, r1=1)
H3K27me1P12ctrl1.bw<-"H3K27me1P12ctrl1.bw"
kpPlotBigWig(kp, data="H3K27me1P12ctrl1.bw", ymin = 0, ymax=25, r0=at$r0, r1=at$r1,col="navy", border = NA)
kpAddLabels(kp, labels = "H3K27me1P12ctrl1", r0=at$r0, r1=at$r0, cex=0.5, label.margin = 0.015,col="navy")

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K27me1P107ctrl1.bw- "H2AK119ubP12ctrl1.bed"
kpPlotRegions(kp, data=H2AK119ubP12ctrl1.bed, col="navy", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 10
at <- autotrack(current.track = 7, total.tracks = 16, r0=0.35, r1=1)
H3K27me1P12ctrl2.bw<-"H3K27me1P12ctrl2.bw"
kpPlotBigWig(kp, data="H3K27me1P12ctrl2.bw", ymin = 0, ymax=25, r0=at$r0, r1=at$r1,col="navy", border = NA)
kpAddLabels(kp, labels = "H3K27me1P12ctrl2", r0=at$r0, r1=at$r1, cex=0.5, label.margin = 0.015,col="navy")

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K27me3ctrl1peaks.bed<- "H3K27me3ctrl1peaks.bed"
kpPlotRegions(kp, data=H3K27me3ctrl1peaks.bed,col="red", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 11
at <- autotrack(current.track = 6, total.tracks = 16, r0=0.35, r1=1)
H3K27me1P107ctrl1.bw<-"H3K27me1P107ctrl1.bw"
kpPlotBigWig(kp, data="H3K27me1P107ctrl1.bw", ymin = 0, ymax=25, r0=at$r0, r1=at$r1,col="red", border = NA)
kpAddLabels(kp, labels = "H3K27me1P107ctrl1", cex=0.5, label.margin = 0.015,col="red",r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K27acCtrl2peaks.bed <- "H3K27acCtrl2peaks.bed"
kpPlotRegions(kp, data=H3K27acCtrl2peaks.bed, col="darkgreen", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 12
at <- autotrack(current.track = 5, total.tracks = 16, r0=0.35, r1=1)
H3K27me1P107ctrl2.bw<- "H3K27me1P107ctrl2.bw"
kpPlotBigWig(kp, data="H3K27me1P107ctrl2.bw", ymin=0, ymax=25, r0=at$r0, r1=at$r1,col="red", border = NA)
kpAddLabels(kp, labels = "H3K27me1P107ctrl2", cex=0.5, label.margin = 0.015,col="red",r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K27acCtrl2peaks.bed<- "H3K27acCtrl2peaks.bed"
kpPlotRegions(kp, data=H3K27acCtrl2peaks.bed, col="magenta", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 13
at <- autotrack(current.track = 4, total.tracks = 16, r0=0.35, r1=1)
H3K27ac.bw <- "H3K27ac.bw"
kpPlotBigWig(kp, data="H3K27ac.bw", ymin=0, ymax=60, r0=at$r0, r1=at$r1,col="darkgreen", border = NA)
kpAddLabels(kp, labels = "H3K27ac", cex=0.5, label.margin = 0.015,col="darkgreen",r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K27acCtrl2peaks.bed<- "H3K27acCtrl2peaks.bed"
kpPlotRegions(kp, data=H3K27acCtrl2peaks.bed, col="darkgreen", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 14
at <- autotrack(current.track = 3, total.tracks = 16, r0=0.35, r1=1)
sorted_230927_i_12_BAP1_P56_mut_H3K4me3_S14_aligned_reads.bam_norm_rnorm.bw <- "sorted_230927_i_12_BAP1_P56_mut_H3K4me3_S14_aligned_reads.bam_norm_rnorm.bw"
kpPlotBigWig(kp, data="sorted_230927_i_12_BAP1_P56_mut_H3K4me3_S14_aligned_reads.bam_norm_rnorm.bw", ymin=0, ymax=300, r0=at$r0, r1=at$r1,col="red", border = NA)
kpAddLabels(kp, labels = "sorted_230927_i_12_BAP1_P56_mut_H3K4me3_S14_aligned_reads.bam_norm_rnorm", cex=0.5, label.margin = 0.015,col="red",r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K4me1ctrl1peaks.bed<- "H3K4me1ctrl1peaks.bed"
kpPlotRegions(kp, data=H3K4me1ctrl1peaks.bed, col="purple", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 15
at <- autotrack(current.track = 2, total.tracks = 16, r0=0.35, r1=1)
sorted_230927_i_7_BAP1_P56_ctrl_H3K27ac_S9_aligned_reads.bam_norm_rnorm.bw <- "sorted_230927_i_7_BAP1_P56_ctrl_H3K27ac_S9_aligned_reads.bam_norm_rnorm.bw"
kpPlotBigWig(kp, data="sorted_230927_i_7_BAP1_P56_ctrl_H3K27ac_S9_aligned_reads.bam_norm_rnorm.bw", ymin=0, ymax=150, r0=at$r0, r1=at$r1,col="navy", border = NA)
kpAddLabels(kp, labels = "sorted_230927_i_7_BAP1_P56_ctrl_H3K27ac_S9_aligned_reads.bam_norm_rnorm", cex=0.5, label.margin = 0.015,col="navy",r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K4me1ctrl1peaks.bed<- "H3K4me1ctrl1peaks.bed"
kpPlotRegions(kp, data=H3K4me1ctrl1peaks.bed, col="purple", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 16
at <- autotrack(current.track = 1, total.tracks = 16, r0=0.35, r1=1)
sorted_230927_i_8_BAP1_P56_mut_H3K27ac_S10_aligned_reads.bam_norm_rnorm.bw <- "sorted_230927_i_8_BAP1_P56_mut_H3K27ac_S10_aligned_reads.bam_norm_rnorm.bw"
kpPlotBigWig(kp, data="sorted_230927_i_8_BAP1_P56_mut_H3K27ac_S10_aligned_reads.bam_norm_rnorm.bw", ymin=0, ymax=150, r0=at$r0, r1=at$r1,col="red", border = NA)
kpAddLabels(kp, labels = "sorted_230927_i_8_BAP1_P56_mut_H3K27ac_S10_aligned_reads.bam_norm_rnorm", cex=0.5, label.margin = 0.015,col="red",r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K4me1ctrl1peaks.bed<- "H3K4me1ctrl1peaks.bed"
kpPlotRegions(kp, data=H3K4me1ctrl1peaks.bed, col="purple", r0=at$r0, r1=at$r1)

#Broad Peak BED Files

at <- autotrack(current.track = 9, total.tracks = 10, r0=0.35, r1=1)
H2AK119ub.bed <- "H2AK119ub.bed"
extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric",
                         qValue = "numeric")
gr_broadPeak <- import.bed(H2AK119ub.bed, extraCols = extraCols_broadPeak)
kpPlotRegions(kp, data=gr_broadPeak, col="blue", r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 7, total.tracks = 10, r0=0.35, r1=1)
H3K27me3ctrl1peaks.bed <- "H3K27me3ctrl1peaks.bed"
extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric",
                         qValue = "numeric")
gr_broadPeak <- import.bed(H3K27me3ctrl1peaks.bed, extraCols = extraCols_broadPeak)
kpPlotRegions(kp, data=gr_broadPeak, col="red", r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 5, total.tracks = 10, r0=0.35, r1=1)
SUZ12.bed <- "SUZ12.bed"
extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric",
                         qValue = "numeric")
gr_broadPeak <- import.bed(SUZ12.bed, extraCols = extraCols_broadPeak)
kpPlotRegions(kp, data=gr_broadPeak, col="green", r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 3, total.tracks = 10, r0=0.35, r1=1)
EZH2.bed <- "EZH2.bed"
extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric",
                         qValue = "numeric")
gr_broadPeak <- import.bed(EZH2.bed, extraCols = extraCols_broadPeak)
kpPlotRegions(kp, data=gr_broadPeak, col="magenta", r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 1, total.tracks = 10, r0=0.35, r1=1)
CpGmm10.bed <- "CpGmm10.bed"
extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric",
                         qValue = "numeric")
gr_broadPeak <- import.bed(CpGmm10.bed, extraCols = extraCols_broadPeak)
kpPlotRegions(kp, data=gr_broadPeak, col="orange", r0=at$r0, r1=at$r1)





#BIGWIG + BED FILE 7
at <- autotrack(current.track = 10, total.tracks = 16, r0=0.35, r1=1)
liverH3K9acEARLY1.bw <- "liverH3K9acEARLY1.bw"
kpPlotBigWig(kp, data="liverH3K9acEARLY1.bw", ymin=0, ymax=100, r0=at$r0, r1=at$r1,col="limegreen", border = NA)
kpAddLabels(kp, labels = "liverH3K9acEARLY1.bw", cex=0.5, label.margin = 0.015,col="limegreen",r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 8
at <- autotrack(current.track = 9, total.tracks = 16, r0=0.35, r1=1)
liverH3K9acLATE1.bw <- "liverH3K9acLATE1.bw"
kpPlotBigWig(kp, data="liverH3K9acLATE1.bw", ymin=0, ymax=100, r0=at$r0, r1=at$r1,col="limegreen", border = NA)
kpAddLabels(kp, labels = "liverH3K9acLATE1.bw", cex=0.5, label.margin = 0.015,col="limegreen",r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 9
at <- autotrack(current.track = 8, total.tracks = 16, r0=0.35, r1=1)
kidneyH3K27acEARLY1.bw<-"kidneyH3K27acEARLY1.bw"
kpPlotBigWig(kp, data="kidneyH3K27acEARLY1.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="darkgreen", border = NA)
kpAddLabels(kp, labels = "kidneyH3K27acEARLY1", r0=at$r0, r1=at$r0, cex=0.5, label.margin = 0.015,col="darkgreen")

#BIGWIG + BED FILE 10
at <- autotrack(current.track = 7, total.tracks = 16, r0=0.35, r1=1)
kidneyH3K27acLATE1.bw<-"kidneyH3K27acLATE1.bw"
kpPlotBigWig(kp, data="kidneyH3K27acLATE1.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="darkgreen", border = NA)
kpAddLabels(kp, labels = "kidneyH3K27acLATE1", r0=at$r0, r1=at$r1, cex=0.5, label.margin = 0.015,col="darkgreen")

#BIGWIG + BED FILE 11
at <- autotrack(current.track = 6, total.tracks = 16, r0=0.35, r1=1)
liverH3K27acEARLY1.bw<-"liverH3K27acEARLY1.bw"
kpPlotBigWig(kp, data="liverH3K27acEARLY1.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="darkgreen", border = NA)
kpAddLabels(kp, labels = "liverH3K27acEARLY1", cex=0.5, label.margin = 0.015,col="darkgreen",r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 12
at <- autotrack(current.track = 5, total.tracks = 16, r0=0.35, r1=1)
liverH3K27acLATE1.bw<- "liverH3K27acLATE1.bw"
kpPlotBigWig(kp, data="liverH3K27acLATE1.bw", ymin=0, ymax=100, r0=at$r0, r1=at$r1,col="darkgreen", border = NA)
kpAddLabels(kp, labels = "liverH3K27acLATE1", cex=0.5, label.margin = 0.015,col="darkgreen",r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 13
at <- autotrack(current.track = 4, total.tracks = 16, r0=0.35, r1=1)
kidneyH3K27me3EARLY1.bw <- "kidneyH3K27me3EARLY1.bw"
kpPlotBigWig(kp, data="kidneyH3K27me3EARLY1.bw", ymin=0, ymax=100, r0=at$r0, r1=at$r1,col="red3", border = NA)
kpAddLabels(kp, labels = "kidneyH3K27me3EARLY1", cex=0.5, label.margin = 0.015,col="red3",r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 14
at <- autotrack(current.track = 3, total.tracks = 16, r0=0.35, r1=1)
kidneyH3K27me3LATE1.bw <- "kidneyH3K27me3LATE1.bw"
kpPlotBigWig(kp, data="kidneyH3K27me3LATE1.bw", ymin=0, ymax=100, r0=at$r0, r1=at$r1,col="red3", border = NA)
kpAddLabels(kp, labels = "kidneyH3K27me3LATE1.bw", cex=0.5, label.margin = 0.015,col="red3",r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 15
at <- autotrack(current.track = 2, total.tracks = 16, r0=0.35, r1=1)
liverH3K27me3EARLY1.bw <- "liverH3K27me3EARLY1.bw"
kpPlotBigWig(kp, data="liverH3K27me3EARLY1.bw", ymin=0, ymax=100, r0=at$r0, r1=at$r1,col="red3", border = NA)
kpAddLabels(kp, labels = "liverH3K27me3EARLY1.bw", cex=0.5, label.margin = 0.015,col="red3",r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 16
at <- autotrack(current.track = 1, total.tracks = 16, r0=0.35, r1=1)
liverH3K27me3LATE1.bw <- "liverH3K27me3LATE1.bw"
kpPlotBigWig(kp, data="liverH3K27me3LATE1.bw", ymin=0, ymax=100, r0=at$r0, r1=at$r1,col="red3", border = NA)
kpAddLabels(kp, labels = "liverH3K27me3LATE1.bw", cex=0.5, label.margin = 0.015,col="red3",r0=at$r0, r1=at$r1)




