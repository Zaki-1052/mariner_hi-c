# ML/cmpts/scripts/A1_run_calder2.R
# Stage A1: Run CALDER2 subcompartment calling on a single Hi-C sample.
#
# CALDER writes large intermediate .Rds files (>100MB) to DATA_ROOT,
# then this script copies the small result files (<1MB) into CODE_ROOT
# for local rsync.
#
# Usage:
#   Rscript --vanilla A1_run_calder2.R <timepoint> <sample> <data_root> <code_root>
#     <timepoint> : 250402 | 250831
#     <sample>    : ctrl_merged | mut_merged
#     <data_root> : HPC data directory (e.g. /expanse/.../sniper)
#     <code_root> : repo directory (e.g. /expanse/.../mariner_hi-c/ML/cmpts)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript A1_run_calder2.R <timepoint> <sample> <data_root> <code_root>")
}

TP        <- args[1]
SAMPLE    <- args[2]
DATA_ROOT <- args[3]
CODE_ROOT <- args[4]

HIC_ROOT <- "/expanse/lustre/projects/csd940/zalibhai/stripes/StripeCaller/data/hic"
hic_path <- file.path(HIC_ROOT, TP, paste0(SAMPLE, ".hic"))
work_dir <- file.path(DATA_ROOT, "calder2", TP, SAMPLE)
repo_dir <- file.path(CODE_ROOT, "outputs", "calder2", TP, SAMPLE)

# ── Pre-flight validation ──

if (!file.exists(hic_path)) stop(paste("Input .hic not found:", hic_path))
if (file.info(hic_path)$size == 0) stop(paste("Input .hic is empty:", hic_path))
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(repo_dir, recursive = TRUE, showWarnings = FALSE)

cat("===========================================\n")
cat("A1: CALDER2 Subcompartment Calling\n")
cat("===========================================\n")
cat(sprintf("Timepoint:  %s\n", TP))
cat(sprintf("Sample:     %s\n", SAMPLE))
cat(sprintf("Input:      %s\n", hic_path))
cat(sprintf("Work dir:   %s\n", work_dir))
cat(sprintf("Repo dir:   %s\n", repo_dir))
cat(sprintf("Start:      %s\n", date()))
cat("===========================================\n\n")

# ── Run CALDER2 ──

library(CALDER)

CALDER(
  contact_file_hic       = hic_path,
  chrs                   = as.character(1:19),
  bin_size               = 50000,
  save_dir               = work_dir,
  save_intermediate_data = TRUE,
  genome                 = "mm10",
  n_cores                = 8,
  sub_domains            = FALSE,
  single_binsize_only    = FALSE
)

# ── Verify outputs ──

tsv_file <- file.path(work_dir, "sub_compartments", "all_sub_compartments.tsv")
bed_file <- file.path(work_dir, "sub_compartments", "all_sub_compartments.bed")
cor_file <- file.path(work_dir, "sub_compartments", "cor_with_ref.txt")

if (!file.exists(tsv_file)) stop(paste("Missing output:", tsv_file))
if (!file.exists(bed_file)) stop(paste("Missing output:", bed_file))
if (!file.exists(cor_file)) stop(paste("Missing output:", cor_file))

# ── Copy results to repo ──

src_dir  <- file.path(work_dir, "sub_compartments")
dest_dir <- file.path(repo_dir, "sub_compartments")
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

src_files <- list.files(src_dir, full.names = TRUE)
for (sf in src_files) {
  file.copy(sf, dest_dir, overwrite = TRUE)
}
cat(sprintf("\nCopied %d result files to repo: %s\n", length(src_files), dest_dir))

# ── Summary statistics ──

comp <- data.table::fread(tsv_file)

expected_chroms <- paste0("chr", 1:19)
chroms_found    <- unique(comp$chr)
missing_chroms  <- setdiff(expected_chroms, chroms_found)
if (length(missing_chroms) > 0) {
  warning(paste("Missing chromosomes in output:", paste(missing_chroms, collapse = ", ")))
}

comp$label_d2 <- sub("^([AB]\\.[12]).*", "\\1", comp$comp_name)
label_dist    <- table(comp$label_d2)
label_pct     <- round(100 * prop.table(label_dist), 1)

cat("\n--- Summary Statistics ---\n")
cat(sprintf("Total rows:         %d\n", nrow(comp)))
cat(sprintf("Chromosomes found:  %d / 19\n", length(intersect(expected_chroms, chroms_found))))
cat("\nDepth-2 label distribution:\n")
for (lbl in names(label_dist)) {
  cat(sprintf("  %-5s  %5d bins  (%5.1f%%)\n", lbl, label_dist[[lbl]], label_pct[[lbl]]))
}

cor_tab <- data.table::fread(cor_file)
cat("\nReference correlations (cor_with_ref.txt):\n")
cat(sprintf("  Median cor: %.3f\n", median(cor_tab$cor, na.rm = TRUE)))
cat(sprintf("  Min cor:    %.3f (chr %s)\n",
    min(cor_tab$cor, na.rm = TRUE),
    cor_tab$chr[which.min(cor_tab$cor)]))
n_low <- sum(cor_tab$cor < 0.3, na.rm = TRUE)
if (n_low > 0) warning(sprintf("%d chromosome(s) have cor < 0.3 — inspect output", n_low))

cat("\n===========================================\n")
cat(sprintf("A1 COMPLETE: %s / %s\n", TP, SAMPLE))
cat(sprintf("Finished:    %s\n", date()))
cat("===========================================\n")
