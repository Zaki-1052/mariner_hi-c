#!/usr/bin/env Rscript
# /expanse/lustre/projects/csd940/zalibhai/mariner/scripts/prepare_loops.R

library(mariner)
library(InteractionSet)
library(GenomicRanges)

ctrlFile <- "/expanse/lustre/projects/csd940/ctea/nf-hic/hiccups/hiccups_results/resorted_ctrl/merged_loops.bedpe"
mutFile <- "/expanse/lustre/projects/csd940/ctea/nf-hic/hiccups/hiccups_results/resorted_mut/merged_loops.bedpe"

bedpe_colnames <- c(
		      "chr1", "x1", "x2", "chr2", "y1", "y2",
		        "name", "score", "strand1", "strand2", "color",
		        "observed", "expectedBL", "expectedDonut", "expectedH", "expectedV",
			  "fdrBL", "fdrDonut", "fdrH", "fdrV",
			  "numCollapsed", "centroid1", "centroid2", "radius"
			  )

read_bedpe <- function(filepath) {
	  bedpe <- read.table(filepath, header = FALSE, skip = 2, 
			                            stringsAsFactors = FALSE, sep = "\t",
						                          comment.char = "")
  colnames(bedpe) <- bedpe_colnames
    bedpe$resolution <- bedpe$x2 - bedpe$x1
    bedpe_5kb <- bedpe[bedpe$resolution == 5000, ]
      bedpe_subset <- bedpe_5kb[1:nrow(bedpe_5kb), ]
      bedpe_subset <- bedpe_subset[bedpe_subset$observed > 0, ]
        return(bedpe_subset)
}

ctrlLoops <- read_bedpe(ctrlFile)
mutLoops <- read_bedpe(mutFile)

cat(sprintf("%d ctrl, %d mut\n", nrow(ctrlLoops), nrow(mutLoops)))

gi_ctrl <- as_ginteractions(ctrlLoops, 
			                                keep.extra.columns = TRUE,
							                            starts.in.df.are.0based = TRUE)

gi_mut <- as_ginteractions(mutLoops,
			                              keep.extra.columns = TRUE,
						                                 starts.in.df.are.0based = TRUE)

gi_list <- list(ctrl = gi_ctrl, mut = gi_mut)
saveRDS(gi_list, "outputs/full/01_ginteractions.rds")

merged <- mergePairs(
		       x = gi_list,
		         radius = 10e3,
		         column = "observed",
			   selectMax = TRUE,
			   method = "manhattan"
			   )

cluster_sizes <- sapply(clusters(merged), nrow)
cat(sprintf("%d merged from %d input\n", 
	                length(merged), length(gi_ctrl) + length(gi_mut)))
cat(sprintf("clusters: %d single, %d paired\n", 
	                sum(cluster_sizes == 1), sum(cluster_sizes == 2)))

saveRDS(merged, "outputs/full/02_merged.rds")

binned <- assignToBins(
		         x = merged,
			   binSize = 5e3,
			   pos1 = "center",
			     pos2 = "center"
			   )

cat(sprintf("%d binned\n", length(binned)))
saveRDS(binned, "outputs/full/03_binned.rds")

buffered <- pixelsToMatrices(
			       x = binned,
			         buffer = 2
			       )

anchor1_width <- unique(width(anchors(buffered, type = "first")))
anchor2_width <- unique(width(anchors(buffered, type = "second")))
cat(sprintf("%d loops, %.0fx%.0f pixels\n", 
	                length(buffered), (anchor1_width-1)/5000, (anchor2_width-1)/5000))

saveRDS(buffered, "outputs/full/04_buffered.rds")

cat(sprintf("Ready for extraction: %d loops x %d pixels = %d total\n",
	                length(buffered), 25, length(buffered) * 25))
