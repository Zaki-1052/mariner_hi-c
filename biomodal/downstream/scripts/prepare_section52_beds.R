# biomodal/downstream/scripts/prepare_section52_beds.R
# Generate sub-gene-feature and metagene-bin BED files for Section 52 analysis
# Parses GENCODE vM25 GTF to extract exon/intron/UTR/splice-site intervals
# for the 124 max-significance genes, formatted for modality regional-frac.
#
# Run from downstream/ directory:
#   Rscript scripts/prepare_section52_beds.R
#
# Then: git add, commit, push. On HPC: git pull.

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomicFeatures)
})

REPO_ROOT <- normalizePath(file.path(getwd(), "../.."))

# ---------------------------------------------------------------------------
# Input paths
# ---------------------------------------------------------------------------
GTF_PATH <- file.path(REPO_ROOT, "abc/reference/gencode.vM25.annotation.gtf.gz")
GENE_TABLE_PATH <- file.path(getwd(), "plots/visualizations/tables/max_significance_genes_merged.tsv")

stopifnot(file.exists(GTF_PATH))
stopifnot(file.exists(GENE_TABLE_PATH))

# Output directory
OUT_DIR <- file.path(getwd(), "modality/section52_beds")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load gene list (124 merged max-significance genes)
# ---------------------------------------------------------------------------
cat("Loading max-significance gene list...\n")
gene_table <- read.table(GENE_TABLE_PATH, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
target_genes <- gene_table$gene
cat(sprintf("  %d target genes loaded.\n", length(target_genes)))

# ---------------------------------------------------------------------------
# Load GENCODE vM25 GTF
# ---------------------------------------------------------------------------
cat("Loading GENCODE vM25 GTF (this takes ~1-2 min)...\n")
gtf <- import(GTF_PATH)

# Filter to target genes only (drastically reduces working set)
gtf_sub <- gtf[gtf$gene_name %in% target_genes]
cat(sprintf("  GTF subset: %d records for %d target genes.\n",
            length(gtf_sub),
            length(unique(gtf_sub$gene_name))))

# Check which target genes are missing from GTF
missing_genes <- setdiff(target_genes, unique(gtf_sub$gene_name))
if (length(missing_genes) > 0) {
  cat(sprintf("  WARNING: %d genes not found in GTF: %s\n",
              length(missing_genes), paste(missing_genes, collapse = ", ")))
}

# ---------------------------------------------------------------------------
# Select canonical transcript per gene
# ---------------------------------------------------------------------------
select_canonical_transcript <- function(gene_gtf) {
  gene_name <- unique(gene_gtf$gene_name)
  txs <- gene_gtf[which(gene_gtf$type == "transcript")]

  # Prefer protein_coding transcripts
  pc_txs <- txs[which(!is.na(txs$transcript_type) & txs$transcript_type == "protein_coding")]
  if (length(pc_txs) == 0) pc_txs <- txs

  # Look for appris_principal_1 tag
  has_tag <- sapply(pc_txs$tag, function(t) {
    if (is.null(t) || all(is.na(t))) return(FALSE)
    "appris_principal_1" %in% t
  })

  if (any(has_tag)) {
    return(pc_txs$transcript_id[which(has_tag)[1]])
  }

  # Fallback: longest transcript by exon span
  tx_ids <- pc_txs$transcript_id
  exons <- gene_gtf[which(gene_gtf$type == "exon" & gene_gtf$transcript_id %in% tx_ids)]

  if (length(exons) == 0) return(NA_character_)

  total_width <- tapply(width(exons), exons$transcript_id, sum)
  names(total_width)[which.max(total_width)]
}

cat("Selecting canonical transcripts...\n")
canonical_tx <- sapply(split(gtf_sub, gtf_sub$gene_name), select_canonical_transcript)
canonical_tx <- canonical_tx[!is.na(canonical_tx)]
cat(sprintf("  %d canonical transcripts selected.\n", length(canonical_tx)))

# ---------------------------------------------------------------------------
# Extract sub-gene features per gene
# ---------------------------------------------------------------------------
SPLICE_WINDOW <- 50  # bp on each side of exon-intron boundary

extract_sub_gene_features <- function(gene_name, tx_id, gtf_sub) {
  tx_gtf <- gtf_sub[which(gtf_sub$transcript_id == tx_id)]
  strand_val <- as.character(strand(tx_gtf))[1]

  # Get exons sorted by genomic position
  exons <- tx_gtf[which(tx_gtf$type == "exon")]
  exons <- sort(exons)
  if (length(exons) == 0) return(NULL)

  # Get UTRs
  utrs <- tx_gtf[which(tx_gtf$type == "UTR")]

  # Classify UTRs as 5' or 3' based on strand and position relative to CDS
  utr_features <- GRanges()
  cds_records <- tx_gtf[which(tx_gtf$type == "CDS")]
  if (length(utrs) > 0 && length(cds_records) > 0) {
    first_cds <- min(start(cds_records))
    last_cds <- max(end(cds_records))

    for (i in seq_along(utrs)) {
      u <- utrs[i]
      if (strand_val == "+") {
        utr_type <- if (end(u) <= first_cds) "5UTR" else "3UTR"
      } else {
        utr_type <- if (start(u) >= last_cds) "5UTR" else "3UTR"
      }
      mcols(u)$feature_type <- utr_type
      utr_features <- c(utr_features, u)
    }
  }

  # Coding exon regions = exons minus UTRs
  if (length(utr_features) > 0) {
    coding_exons <- setdiff(reduce(exons), reduce(utr_features))
  } else {
    coding_exons <- reduce(exons)
  }

  # Derive introns from gaps between exons
  all_exons_reduced <- reduce(exons)
  gene_range <- range(all_exons_reduced)
  introns <- setdiff(gene_range, all_exons_reduced)

  # Create splice site windows at internal exon-intron boundaries
  splice_sites <- GRanges()
  if (length(all_exons_reduced) > 1) {
    chr <- as.character(seqnames(all_exons_reduced))[1]
    for (i in seq_len(length(all_exons_reduced) - 1)) {
      exon_end <- end(all_exons_reduced[i])
      next_exon_start <- start(all_exons_reduced[i + 1])

      # Donor site: exon_end ± SPLICE_WINDOW
      donor_start <- max(exon_end - SPLICE_WINDOW + 1, start(gene_range))
      donor_end <- min(exon_end + SPLICE_WINDOW, end(gene_range))
      donor <- GRanges(seqnames = chr,
                       ranges = IRanges(donor_start, donor_end),
                       strand = strand_val)
      mcols(donor)$feature_type <- "SpliceSite_Donor"
      mcols(donor)$rank <- i

      # Acceptor site: next_exon_start ± SPLICE_WINDOW
      acceptor_start <- max(next_exon_start - SPLICE_WINDOW, start(gene_range))
      acceptor_end <- min(next_exon_start + SPLICE_WINDOW - 1, end(gene_range))
      acceptor <- GRanges(seqnames = chr,
                          ranges = IRanges(acceptor_start, acceptor_end),
                          strand = strand_val)
      mcols(acceptor)$feature_type <- "SpliceSite_Acceptor"
      mcols(acceptor)$rank <- i

      splice_sites <- c(splice_sites, donor, acceptor)
    }
  }

  # Subtract splice sites from introns to avoid double-counting
  splice_reduced <- reduce(splice_sites)
  if (length(introns) > 0 && length(splice_reduced) > 0) {
    introns_trimmed <- setdiff(introns, splice_reduced)
  } else {
    introns_trimmed <- introns
  }

  # Build output data.frame
  result <- data.frame(
    Chromosome = character(0), Start = integer(0), End = integer(0),
    Annotation = character(0), Name = character(0),
    stringsAsFactors = FALSE
  )

  # Modality BED format: both start and end are GTF_coord - 1
  # (matches gencode.vM25.mouse.genes.annotation.bed.gz convention)
  add_features <- function(gr, feature_type, result) {
    if (length(gr) == 0) return(result)
    gr <- sort(gr)
    n <- length(gr)
    rbind(result, data.frame(
      Chromosome = as.character(seqnames(gr)),
      Start = start(gr) - 1,
      End = end(gr) - 1,
      Annotation = feature_type,
      Name = paste0(gene_name, "_", feature_type, "_",
                    seq_len(n)),
      stringsAsFactors = FALSE
    ))
  }

  add_splice_features <- function(gr, result) {
    if (length(gr) == 0) return(result)
    gr <- sort(gr)
    rbind(result, data.frame(
      Chromosome = as.character(seqnames(gr)),
      Start = start(gr) - 1,
      End = end(gr) - 1,
      Annotation = mcols(gr)$feature_type,
      Name = paste0(gene_name, "_", mcols(gr)$feature_type, "_",
                    mcols(gr)$rank),
      stringsAsFactors = FALSE
    ))
  }

  # UTRs
  if (length(utr_features) > 0) {
    for (ft in c("5UTR", "3UTR")) {
      sub <- utr_features[which(mcols(utr_features)$feature_type == ft)]
      result <- add_features(sub, ft, result)
    }
  }

  # Coding exons
  result <- add_features(coding_exons, "Exon", result)

  # Splice sites (preserve donor/acceptor labels)
  result <- add_splice_features(splice_sites, result)

  # Trimmed introns
  result <- add_features(introns_trimmed, "Intron", result)

  # Filter out any negative-width intervals
  result <- result[result$End >= result$Start, ]
  result
}

cat("Extracting sub-gene features...\n")
all_features <- list()
strand_lookup <- data.frame(gene = character(0), strand = character(0),
                            stringsAsFactors = FALSE)

for (gene_name in names(canonical_tx)) {
  tx_id <- canonical_tx[gene_name]
  features <- extract_sub_gene_features(gene_name, tx_id, gtf_sub)
  if (!is.null(features) && nrow(features) > 0) {
    all_features[[gene_name]] <- features
  }

  tx_records <- gtf_sub[which(gtf_sub$transcript_id == tx_id)]
  strand_val <- as.character(strand(tx_records))[1]
  strand_lookup <- rbind(strand_lookup,
                         data.frame(gene = gene_name, strand = strand_val,
                                    stringsAsFactors = FALSE))
}

features_df <- do.call(rbind, all_features)
rownames(features_df) <- NULL

# Strip chr prefix for modality format
features_df$Chromosome <- gsub("^chr", "", features_df$Chromosome)

# Sort by chromosome and start
features_df <- features_df[order(features_df$Chromosome, features_df$Start), ]

cat(sprintf("  %d sub-gene feature intervals across %d genes.\n",
            nrow(features_df), length(all_features)))

# Feature-type summary
feature_counts <- table(features_df$Annotation)
cat("  Feature breakdown:\n")
for (ft in names(feature_counts)) {
  cat(sprintf("    %s: %d\n", ft, feature_counts[ft]))
}

# ---------------------------------------------------------------------------
# Generate metagene bins (100 equal-width bins per gene)
# ---------------------------------------------------------------------------
cat("\nGenerating metagene bins...\n")
N_BINS <- 100

metagene_rows <- list()
for (gene_name in names(canonical_tx)) {
  gene_info <- gene_table[gene_table$gene == gene_name, ]
  if (nrow(gene_info) == 0) next

  chr <- gsub("^chr", "", gene_info$chr[1])
  gene_start <- gene_info$start[1]
  gene_end <- gene_info$end[1]
  gene_len <- gene_end - gene_start

  bin_width <- gene_len / N_BINS
  bins <- data.frame(
    Chromosome = chr,
    Start = as.integer(gene_start + floor((0:(N_BINS - 1)) * bin_width)),
    End = as.integer(gene_start + floor((1:N_BINS) * bin_width)),
    Annotation = "Bin",
    Name = paste0(gene_name, "_bin_", sprintf("%03d", 1:N_BINS)),
    stringsAsFactors = FALSE
  )
  # Ensure last bin reaches gene_end exactly
  bins$End[N_BINS] <- gene_end

  metagene_rows[[gene_name]] <- bins
}

metagene_df <- do.call(rbind, metagene_rows)
rownames(metagene_df) <- NULL
metagene_df <- metagene_df[order(metagene_df$Chromosome, metagene_df$Start), ]

cat(sprintf("  %d metagene bins across %d genes.\n",
            nrow(metagene_df), length(metagene_rows)))

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------
cat("\nWriting BED files...\n")

# BED 1: sub-gene features
features_path <- file.path(OUT_DIR, "sub_gene_features.annotation.bed.gz")
gz <- gzfile(features_path, "w")
write.table(features_df, gz, sep = "\t", row.names = FALSE, quote = FALSE)
close(gz)
cat(sprintf("  Wrote: %s (%d rows)\n", features_path, nrow(features_df)))

# BED 2: metagene bins
metagene_path <- file.path(OUT_DIR, "metagene_bins.annotation.bed.gz")
gz <- gzfile(metagene_path, "w")
write.table(metagene_df, gz, sep = "\t", row.names = FALSE, quote = FALSE)
close(gz)
cat(sprintf("  Wrote: %s (%d rows)\n", metagene_path, nrow(metagene_df)))

# Strand lookup for metagene bin reversal in section_52 R script
strand_path <- file.path(OUT_DIR, "gene_strand_lookup.tsv")
write.table(strand_lookup, strand_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Wrote: %s (%d genes)\n", strand_path, nrow(strand_lookup)))

# ---------------------------------------------------------------------------
# Verification checks
# ---------------------------------------------------------------------------
cat("\n========================================\n")
cat("VERIFICATION\n")
cat("========================================\n")
cat("Run these checks before committing:\n\n")
cat("  # Confirm chr prefix stripped\n")
cat("  gunzip -c modality/section52_beds/sub_gene_features.annotation.bed.gz | head -3\n\n")
cat("  # Confirm Syt1 exons present\n")
cat("  gunzip -c modality/section52_beds/sub_gene_features.annotation.bed.gz | grep 'Syt1_Exon' | head\n\n")
cat("  # Row counts\n")
cat(sprintf("  # Expected: ~%d sub-gene features, %d metagene bins\n",
            nrow(features_df), nrow(metagene_df)))
cat("  gunzip -c modality/section52_beds/sub_gene_features.annotation.bed.gz | wc -l\n")
cat("  gunzip -c modality/section52_beds/metagene_bins.annotation.bed.gz | wc -l\n\n")
cat("Next steps:\n")
cat("  1. git add modality/section52_beds/ && git commit && git push\n")
cat("  2. On HPC: cd /expanse/lustre/projects/csd940/zalibhai/mariner_hi-c && git pull\n")
cat("  3. On HPC: sbatch biomodal/downstream/modality/run_section52_hpc.sb\n")
cat("\nDone.\n")
