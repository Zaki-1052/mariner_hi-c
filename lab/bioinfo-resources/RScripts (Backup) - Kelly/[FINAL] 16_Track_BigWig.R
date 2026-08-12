#BIGWIG + BED FILE 1
at <- autotrack(current.track = 18, total.tracks = 18, r0=0.35, r1=1)
H2AK119ubCtrl.bw<-"H2AK119ubCtrl.bw"
kpPlotBigWig(kp, data="H2AK119ubCtrl.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="royalblue2", border = NA)
kpAddLabels(kp, labels = "H2AK119ubCtrl", r0=at$r0, r1=at$r1, cex=0.5, label.margin = 0.015,col="royalblue2")

#BIGWIG + BED FILE 2
at <- autotrack(current.track = 17, total.tracks = 18, r0=0.35, r1=1)
H2AK119ubMut.bw<-"H2AK119ubMut.bw"
kpPlotBigWig(kp, data="H2AK119ubMut.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="royalblue2", border = NA)
kpAddLabels(kp, labels = "H2AK119ubMut", r0=at$r0, r1=at$r1, cex=0.5, label.margin = 0.015,col="royalblue2")

#BIGWIG + BED FILE 3
at <- autotrack(current.track = 16, total.tracks = 18, r0=0.35, r1=1)
H3K27acCtrl.bw<-"H3K27acCtrl.bw"
kpPlotBigWig(kp, data="H3K27acCtrl.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="limegreen", border = NA)
kpAddLabels(kp, labels = "H3K27acCtrl", cex=0.5, label.margin = 0.015,col="limegreen",r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 4
at <- autotrack(current.track = 15, total.tracks = 18, r0=0.35, r1=1)
H3K27acMut.bw<- "H3K27acMut.bw"
kpPlotBigWig(kp, data="H3K27acMut.bw", ymin=0, ymax=100, r0=at$r0, r1=at$r1,col="limegreen", border = NA)
kpAddLabels(kp, labels = "H3K27acMut", cex=0.5, label.margin = 0.015,col="limegreen",r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 5
at <- autotrack(current.track = 14, total.tracks = 18, r0=0.35, r1=1)
H3K27me3Ctrl.bw<-"H3K27me3Ctrl.bw"
kpPlotBigWig(kp, data="H3K27me3Ctrl.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="tomato", border = NA)
kpAddLabels(kp, labels = "H3K27me3Ctrl", r0=at$r0, r1=at$r1, cex=0.5, label.margin = 0.015,col="tomato")

#BIGWIG + BED FILE 6
at <- autotrack(current.track = 13, total.tracks = 18, r0=0.35, r1=1)
H3K27me3Mut.bw<-"H3K27me3Mut.bw"
kpPlotBigWig(kp, data="H3K27me3Mut.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="tomato", border = NA)
kpAddLabels(kp, labels = "H3K27me3Mut", r0=at$r0, r1=at$r1, cex=0.5, label.margin = 0.015,col="tomato")

#BIGWIG + BED FILE 7
at <- autotrack(current.track = 12, total.tracks = 18, r0=0.35, r1=1)
H3K27me1Ctrl.bw<-"H3K27me1Ctrl.bw"
kpPlotBigWig(kp, data="H3K27me1Ctrl.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="orange", border = NA)
kpAddLabels(kp, labels = "H3K27me1Ctrl", r0=at$r0, r1=at$r1, cex=0.5, label.margin = 0.015,col="orange")

#BIGWIG + BED FILE 8
at <- autotrack(current.track = 11, total.tracks = 18, r0=0.35, r1=1)
H3K27me1Mut.bw<-"H3K27me1Mut.bw"
kpPlotBigWig(kp, data="H3K27me1Mut.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="orange", border = NA)
kpAddLabels(kp, labels = "H3K27me1Mut", r0=at$r0, r1=at$r1, cex=0.5, label.margin = 0.015,col="orange")

#BIGWIG + BED FILE 9
at <- autotrack(current.track = 10, total.tracks = 18, r0=0.35, r1=1)
H3K4me3Ctrl.bw<-"H3K4me3Ctrl.bw"
kpPlotBigWig(kp, data="H3K4me3Ctrl.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="aquamarine2", border = NA)
kpAddLabels(kp, labels = "H3K4me3Ctrl", r0=at$r0, r1=at$r1, cex=0.5, label.margin = 0.015,col="aquamarine2")

#BIGWIG + BED FILE 10
at <- autotrack(current.track = 9, total.tracks = 18, r0=0.35, r1=1)
H3K4me3Mut.bw<-"H3K4me3Mut.bw"
kpPlotBigWig(kp, data="H3K4me3Mut.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="aquamarine2", border = NA)
kpAddLabels(kp, labels = "H3K4me3Mut", r0=at$r0, r1=at$r1, cex=0.5, label.margin = 0.015,col="aquamarine2")

#BIGWIG + BED FILE 11
at <- autotrack(current.track = 8, total.tracks = 18, r0=0.35, r1=1)
H2AzCtrl.bw<-"H2AzCtrl.bw"
kpPlotBigWig(kp, data="H2AzCtrl.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="orchid2", border = NA)
kpAddLabels(kp, labels = "H2AzCtrl", cex=0.5, label.margin = 0.015,col="orchid2",r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 12
at <- autotrack(current.track = 7, total.tracks = 18, r0=0.35, r1=1)
H2AzMut.bw<- "H2AzMut.bw"
kpPlotBigWig(kp, data="H2AzMut.bw", ymin=0, ymax=100, r0=at$r0, r1=at$r1,col="orchid2", border = NA)
kpAddLabels(kp, labels = "H2AzMut", cex=0.5, label.margin = 0.015,col="orchid2",r0=at$r0, r1=at$r1)

#BIGWIG + BED FILE 13
at <- autotrack(current.track = 6, total.tracks = 18, r0=0.35, r1=1)
H3K9acCtrl.bw<-"H3K9acCtrl.bw"
kpPlotBigWig(kp, data="H3K9acCtrl.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="seagreen4", border = NA)
kpAddLabels(kp, labels = "H3K9acCtrl", r0=at$r0, r1=at$r1, cex=0.5, label.margin = 0.015,col="seagreen4")

#BIGWIG + BED FILE 14
at <- autotrack(current.track = 5, total.tracks = 18, r0=0.35, r1=1)
H3K9acMut.bw<-"H3K9acMut.bw"
kpPlotBigWig(kp, data="H3K9acMut.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="seagreen4", border = NA)
kpAddLabels(kp, labels = "H3K9acMut", r0=at$r0, r1=at$r1, cex=0.5, label.margin = 0.015,col="seagreen4")

#BIGWIG + BED FILE 15
at <- autotrack(current.track = 4, total.tracks = 18, r0=0.35, r1=1)
H4K12acCtrl.bw<-"H4K12acCtrl.bw"
kpPlotBigWig(kp, data="H4K12acCtrl.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="slateblue3", border = NA)
kpAddLabels(kp, labels = "H4K12acCtrl", r0=at$r0, r1=at$r1, cex=0.5, label.margin = 0.015,col="slateblue3")

#BIGWIG + BED FILE 16
at <- autotrack(current.track = 3, total.tracks = 18, r0=0.35, r1=1)
H4K12acMut.bw<-"H4K12acMut.bw"
kpPlotBigWig(kp, data="H4K12acMut.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="slateblue3", border = NA)
kpAddLabels(kp, labels = "H4K12acMut", r0=at$r0, r1=at$r1, cex=0.5, label.margin = 0.015,col="slateblue3")

#BIGWIG + BED FILE 17
at <- autotrack(current.track = 2, total.tracks = 18, r0=0.35, r1=1)
RNActrl.bw<-"RNActrl.bw"
kpPlotBigWig(kp, data="RNActrl.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="deeppink2", border = NA)
kpAddLabels(kp, labels = "RNActrl", r0=at$r0, r1=at$r1, cex=0.5, label.margin = 0.015,col="deeppink2")

#BIGWIG + BED FILE 18
at <- autotrack(current.track = 1, total.tracks = 18, r0=0.35, r1=1)
RNAmut.bw<-"RNAmut.bw"
kpPlotBigWig(kp, data="RNAmut.bw", ymin = 0, ymax=100, r0=at$r0, r1=at$r1,col="deeppink2", border = NA)
kpAddLabels(kp, labels = "RNAmut", r0=at$r0, r1=at$r1, cex=0.5, label.margin = 0.015,col="deeppink2")





