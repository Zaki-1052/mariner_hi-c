# For Juicebox visualization (if you want to compare with original .hic files)
counts <- readRDS("outputs/full/06_counts_matrix.rds")
y <- readRDS("outputs/full/06_edgeR_input.rds")
write.table(
  data.frame(
    chr1 = y$genes$chr1,
    start1 = y$genes$start1,
    end1 = y$genes$end1,
    chr2 = y$genes$chr2,
    start2 = y$genes$start2,
    end2 = y$genes$end2,
    ctrl = counts[,1],
    mut = counts[,2]
  ),
  file = "loops_with_counts.bedpe",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# For IGV/UCSC browser (BED format for each sample)
for (sample in colnames(counts)) {
  bed_data <- data.frame(
    chr = y$genes$chr1,
    start = y$genes$start1,
    end = y$genes$end1,
    name = paste0("loop_", 1:nrow(counts)),
    score = counts[,sample]
  )
  write.table(bed_data,
              file = paste0("loops_", sample, ".bed"),
              sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
}
