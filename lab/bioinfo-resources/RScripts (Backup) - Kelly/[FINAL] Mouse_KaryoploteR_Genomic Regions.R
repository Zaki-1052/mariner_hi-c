library(karyoploteR)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(rtracklayer)

setwd("/Users/kellywang/RStudio_Scripts")

TP53.region <- toGRanges("chr8:20,321,012-20,366,504")
kp <- plotKaryotype(genome = "mm10", zoom = TP53.region,cex=1)

genes.data <- makeGenesDataFromTxDb(TxDb.Mmusculus.UCSC.mm10.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
genes.data <- mergeTranscripts(genes.data)
genes.data <- addGeneNames(genes.data,org.Mm.eg.db)
kpPlotGenes(kp, data=genes.data,r0=0, r1=0.15, gene.name.cex = 0.5)

source("[FINAL] 16_Track_BigWig.R")



#BIGWIG + BED FILE 1
at <- autotrack(current.track = 8, total.tracks = 8, r0=0.35, r1=1)
H2AK119ubCerebellumEarly1.bw<-"H2AK119ubCerebellumEarly1.bw"
kpPlotBigWig(kp, data="H2AK119ubCerebellumEarly1.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="dodgerblue4", border = NA)
kpAddLabels(kp, labels = "H2AK119ubCerebellumEarly1", r0=at$r0, r1=at$r0, cex=0.5, label.margin = 0.015,col="dodgerblue4")

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K27me1P107ctrl1.bw- "H2AK119ubP12ctrl1.bed"
kpPlotRegions(kp, data=H2AK119ubP12ctrl1.bed, col="dodgerblue4", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 2
at <- autotrack(current.track = 7, total.tracks = 8, r0=0.35, r1=1)
H2AK119ubCerebellumLate1.bw<-"H2AK119ubCerebellumLate1.bw"
kpPlotBigWig(kp, data="H2AK119ubCerebellumLate1.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="dodgerblue4", border = NA)
kpAddLabels(kp, labels = "H2AK119ubCerebellumLate1", r0=at$r0, r1=at$r1, cex=0.5, label.margin = 0.015,col="dodgerblue4")

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K27me3ctrl1peaks.bed<- "H3K27me3ctrl1peaks.bed"
kpPlotRegions(kp, data=H3K27me3ctrl1peaks.bed,col="tomato3", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 3
at <- autotrack(current.track = 6, total.tracks = 8, r0=0.35, r1=1)
H3K27acCerebellumEarly1.bw<-"H3K27acCerebellumEarly1.bw"
kpPlotBigWig(kp, data="H3K27acCerebellumEarly1.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="limegreen", border = NA)
kpAddLabels(kp, labels = "H3K27acCerebellumEarly1", cex=0.5, label.margin = 0.015,col="limegreen",r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K27acCtrl2peaks.bed <- "H3K27acCtrl2peaks.bed"
kpPlotRegions(kp, data=H3K27acCtrl2peaks.bed, col="limegreen", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 4
at <- autotrack(current.track = 5, total.tracks = 8, r0=0.35, r1=1)
H3K27acCerebellumLate1.bw<- "H3K27acCerebellumLate1.bw"
kpPlotBigWig(kp, data="H3K27acCerebellumLate1.bw", ymin=0, ymax=100, r0=at$r0, r1=at$r1,col="limegreen", border = NA)
kpAddLabels(kp, labels = "H3K27acCerebellumLate1", cex=0.5, label.margin = 0.015,col="limegreen",r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 1, total.tracks = 5, r0=0.35, r1=1)
CpGmm10.bed<- "CpGmm10.bed"
kpPlotRegions(kp, data=CpGmm10.bed, col="limegreen", r0=at$r0, r1=at$r1)
kpAddLabels(kp, labels = "CpGmm10", cex=0.5, label.margin = 0.015,col="limegreen",r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 5
at <- autotrack(current.track = 4, total.tracks = 8, r0=0.35, r1=1)
H3K27me3CerebellumEarly1.bw <- "H3K27me3CerebellumEarly1"
kpPlotBigWig(kp, data= "H3K27me3CerebellumEarly1.bw", ymin=0, ymax=100, r0=at$r0, r1=at$r1,col="tomato3", border = NA)
kpAddLabels(kp, labels = "H3K27me3CerebellumEarly1", cex=0.5, label.margin = 0.015,col="tomato3",r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K27acCtrl2peaks.bed<- "H3K27acCtrl2peaks.bed"
kpPlotRegions(kp, data=H3K27acCtrl2peaks.bed, col="darkgreen", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 6
at <- autotrack(current.track = 3, total.tracks = 8, r0=0.35, r1=1)
H3K27me3CerebellumLate1.bw <- "H3K27me3CerebellumLate1.bw"
kpPlotBigWig(kp, data="H3K27me3CerebellumLate1.bw", ymin=0, ymax=100, r0=at$r0, r1=at$r1,col="tomato3", border = NA)
kpAddLabels(kp, labels = "H3K27me3CerebellumLate1", cex=0.5, label.margin = 0.015,col="tomato3",r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K4me1ctrl1peaks.bed<- "H3K4me1ctrl1peaks.bed"
kpPlotRegions(kp, data=H3K4me1ctrl1peaks.bed, col="purple", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 7
at <- autotrack(current.track = 2, total.tracks = 8, r0=0.35, r1=1)
RNAcerebellumEarly1.bw <- "RNAcerebellumEarly1.bw"
kpPlotBigWig(kp, data="RNAcerebellumEarly1.bw", ymin=0, ymax=100, r0=at$r0, r1=at$r1,col="orange", border = NA)
kpAddLabels(kp, labels = "RNAcerebellumEarly1", cex=0.5, label.margin = 0.015,col="orange",r0=at$r0, r1=at$r1)

at <- autotrack(current.track = 0, total.tracks = 0, r0=0.35, r1=1)
H3K4me1ctrl1peaks.bed<- "H3K4me1ctrl1peaks.bed"
kpPlotRegions(kp, data=H3K4me1ctrl1peaks.bed, col="purple", r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 8
at <- autotrack(current.track = 1, total.tracks = 8, r0=0.35, r1=1)
RNAcerebellumLat11.bw <- "RNAcerebellumLate1.bw"
kpPlotBigWig(kp, data="RNAcerebellumLate1.bw", ymin=0, ymax=100, r0=at$r0, r1=at$r1,col="orange", border = NA)
kpAddLabels(kp, labels = "RNAcerebellumLate1", cex=0.5, label.margin = 0.015,col="orange",r0=at$r0, r1=at$r1)

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




