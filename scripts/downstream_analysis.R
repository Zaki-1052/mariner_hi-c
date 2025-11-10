#!/usr/bin/env Rscript
# scripts/downstream_analysis.R
# Downstream Analysis: Merged Loops & Anchor Characterization
# Author: Zakir Alibhai
# Date: 2025-10-30
#
# Purpose:
#   1. Merge final results from 5kb, 10kb, 25kb with overlap removal
#   2. Validate all loops files for volcano plots
#   3. Characterize loop anchors with genomic features and gene annotations
#   4. Export to multiple formats (TSV, RDS, BEDPE)
#
# Usage:
#   Rscript scripts/downstream_analysis.R [--gtf PATH] [--tolerance KB]
#
# Arguments:
#   --gtf PATH       Path to GTF file for gene annotations (optional)
#   --tolerance KB   Overlap tolerance in kb (default: 10)

# =============================================================================
# SETUP
# =============================================================================

cat("\n")
cat("========================================\n")
cat("Downstream Analysis: Merged Loops\n")
cat("========================================\n\n")

suppressPackageStartupMessages({
		  library(mariner)
		  	    library(InteractionSet)
		  	    library(GenomicRanges)
			    	      library(tidyverse)
			    	      library(ggplot2)
				      	        library(patchwork)
})

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
gtf_file <- NULL
tolerance_kb <- 10

if (length(args) > 0) {
		  for (i in seq_along(args)) {
			  		      if (args[i] == "--gtf" && i < length(args)) {
						      			            gtf_file <- args[i + 1]
    } else if (args[i] == "--tolerance" && i < length(args)) {
	    	          tolerance_kb <- as.numeric(args[i + 1])
            }
  }
}

cat(sprintf("Configuration:\n"))
cat(sprintf("  Overlap tolerance: %d kb\n", tolerance_kb))
if (!is.null(gtf_file)) {
		  cat(sprintf("  GTF file: %s\n", gtf_file))
} else {
		  cat("  GTF file: Using TxDb.Mmusculus.UCSC.mm10.knownGene\n")
}
cat("\n")

# Set working directory (adjust if needed)
if (dir.exists("/expanse/lustre/projects/csd940/zalibhai/mariner")) {
		  setwd("/expanse/lustre/projects/csd940/zalibhai/mariner")
} else if (basename(getwd()) == "mariner-final") {
		  # Already in project directory
} else {
		  warning("Working directory may not be correct. Current: ", getwd())
}

# Create output directory
output_dir <- "outputs/merged_loops"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Output directory: %s\n\n", output_dir))

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Match loops across GInteractions with tolerance
#'
#' @param coords1 GInteractions object
#' @param coords2 GInteractions object
#' @param tolerance_bp Tolerance in base pairs
#' @return data.frame with index1, index2, distance
match_loops <- function(coords1, coords2, tolerance_bp) {
		  # Convert to data.frame with coordinates
		  df1 <- data.frame(
				    			        chr1 = as.character(seqnames(anchors(coords1, "first"))),
												    start1 = start(anchors(coords1, "first")),
												    end1 = end(anchors(coords1, "first")),
												    				        chr2 = as.character(seqnames(anchors(coords1, "second"))),
												    				        start2 = start(anchors(coords1, "second")),
																						    end2 = end(anchors(coords1, "second")),
																						    index1 = 1:length(coords1)
																						    					      )

  df2 <- data.frame(
		    		        chr1 = as.character(seqnames(anchors(coords2, "first"))),
								    start1 = start(anchors(coords2, "first")),
								    end1 = end(anchors(coords2, "first")),
								    			        chr2 = as.character(seqnames(anchors(coords2, "second"))),
								    			        start2 = start(anchors(coords2, "second")),
																    end2 = end(anchors(coords2, "second")),
																    index2 = 1:length(coords2)
																    				      )

      # Find overlaps
      matches <- data.frame(index1 = integer(), index2 = integer(), distance = numeric())

      for (i in 1:nrow(df1)) {
	      	        # Same chromosome check
	      	        same_chr <- df2$chr1 == df1$chr1[i] & df2$chr2 == df1$chr2[i]

              # Distance check (within tolerance for both anchors)
              dist1 <- abs(df2$start1 - df1$start1[i])
	      	    dist2 <- abs(df2$start2 - df1$start2[i])

	      	    close_enough <- same_chr & dist1 <= tolerance_bp & dist2 <= tolerance_bp

		    	        if (any(close_enough)) {
								      # Take closest match
								      distances <- dist1 + dist2
		    	          best_match <- which(close_enough)[which.min(distances[close_enough])]
				  		        matches <- rbind(matches, data.frame(
											       							             index1 = i,
																			     								             index2 = best_match,
																			     								             distance = distances[best_match]
																												     									           ))
				  		      }
		    	      }

            return(matches)
}

# =============================================================================
# SECTION 1: LOAD FINAL RESULTS FROM ALL RESOLUTIONS
# =============================================================================

cat("\n========================================\n")
cat("SECTION 1: Loading Final Results\n")
cat("\n========================================\n")

resolutions <- c(5000, 10000, 25000)
final_results_list <- list()
coords_list <- list()

for (res in resolutions) {
		  res_kb <- res / 1000
  cat(sprintf("Loading %dkb results...\n", res_kb))

      # Load final results TSV
      final_file <- sprintf("outputs/edgeR_results_res_%dkb/primary_analysis/final_results.tsv", res_kb)
      coords_file <- sprintf("outputs/res_%dkb/03_binned.rds", res_kb)

            if (!file.exists(final_file)) {
		    	          warning(sprintf("  Final results file not found: %s", final_file))
              next
	      	  }

            if (!file.exists(coords_file)) {
		    	          warning(sprintf("  Coordinates file not found: %s", coords_file))
	              next
		      	    }

	            # Load data
	            final_df <- read.table(final_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
		            coords <- readRDS(coords_file)

		    	  cat(sprintf("  Final results: %d loops\n", nrow(final_df)))
			  	  cat(sprintf("  Coordinates: %d interactions\n", length(coords)))

			  	    # Match final results to coordinates
			  	    coords_df <- data.frame(
							      				        chr1 = as.character(seqnames(anchors(coords, "first"))),
																	    start1 = start(anchors(coords, "first")),
																	    end1 = end(anchors(coords, "first")),
																	    					        chr2 = as.character(seqnames(anchors(coords, "second"))),
																	    					        start2 = start(anchors(coords, "second")),
																													    end2 = end(anchors(coords, "second"))
																													  )

			  	    # Find matching indices
			  	    matches <- numeric(nrow(final_df))
				    	      for (i in 1:nrow(final_df)) {
						      		          idx <- which(
											     				             coords_df$chr1 == final_df$chr1[i] &
																	     						           coords_df$start1 == final_df$start1[i] &
																								   							         coords_df$end1 == final_df$end1[i] &
																																 								       coords_df$chr2 == final_df$chr2[i] &
																																								       								             coords_df$start2 == final_df$start2[i] &
																																																	     									           coords_df$end2 == final_df$end2[i]
																																																										   									       )
				    	        if (length(idx) > 0) {
										      matches[i] <- idx[1]
								    }
								  }

				    	      # Filter to matched loops
				    	      matches <- matches[matches > 0]

					      	        if (length(matches) == 0) {
											    warning(sprintf("  No matches found between final results and coordinates"))
					      	          next
							  		    }

					      	        # Store filtered coordinates and results
					      	        coords_filtered <- coords[matches]
									  final_df_filtered <- final_df[matches > 0, ]

									  # Add resolution column
									  final_df_filtered$resolution <- res
									  		    final_df_filtered$resolution_kb <- res_kb

									  		    final_results_list[[as.character(res)]] <- final_df_filtered
											    		      coords_list[[as.character(res)]] <- coords_filtered

											    		      cat(sprintf("  Matched: %d loops\n\n", length(matches)))
}

if (length(final_results_list) == 0) {
		  stop("ERROR: No final results loaded. Run generate_final_results.py first.")
}

cat(sprintf("✓ Loaded final results from %d resolutions\n", length(final_results_list)))
cat(sprintf("  Total loops before merging: %d\n\n", sum(sapply(final_results_list, nrow))))

# =============================================================================
# SECTION 1B: MERGE ALL_RESULTS (INCLUDING NON-SIGNIFICANT LOOPS)
# =============================================================================

cat("\n========================================\n")
cat("SECTION 1B: Merging All Results (Complete Dataset)\n")
cat("\n========================================\n")

cat("Loading all_results_primary.tsv from all resolutions...\n")
cat("(This includes ALL loops regardless of significance)\n\n")

all_results_list <- list()
all_coords_list <- list()

for (res in resolutions) {
	  res_kb <- res / 1000
  cat(sprintf("Loading %dkb all_results...\n", res_kb))

    # Load all results TSV
    all_results_file <- sprintf("outputs/edgeR_results_res_%dkb/primary_analysis/all_results_primary.tsv", res_kb)
    coords_file <- sprintf("outputs/res_%dkb/03_binned.rds", res_kb)

      if (!file.exists(all_results_file)) {
	          warning(sprintf("  All results file not found: %s", all_results_file))
        next
	  }

      if (!file.exists(coords_file)) {
	          warning(sprintf("  Coordinates file not found: %s", coords_file))
          next
	    }

        # Load data
        all_res_df <- read.table(all_results_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        coords <- readRDS(coords_file)

	  cat(sprintf("  All results: %d loops\n", nrow(all_res_df)))
	  cat(sprintf("  Coordinates: %d interactions\n", length(coords)))

	    # Match all results to coordinates
	    coords_df <- data.frame(
				        chr1 = as.character(seqnames(anchors(coords, "first"))),
					    start1 = start(anchors(coords, "first")),
					    end1 = end(anchors(coords, "first")),
					        chr2 = as.character(seqnames(anchors(coords, "second"))),
					        start2 = start(anchors(coords, "second")),
						    end2 = end(anchors(coords, "second"))
						  )

	    # Find matching indices
	    matches <- numeric(nrow(all_res_df))
	      for (i in 1:nrow(all_res_df)) {
		          idx <- which(
				             coords_df$chr1 == all_res_df$chr1[i] &
						           coords_df$start1 == all_res_df$start1[i] &
							         coords_df$end1 == all_res_df$end1[i] &
								       coords_df$chr2 == all_res_df$chr2[i] &
								             coords_df$start2 == all_res_df$start2[i] &
									           coords_df$end2 == all_res_df$end2[i]
									       )
	        if (length(idx) > 0) {
			      matches[i] <- idx[1]
		    }
		  }

	      # Filter to matched loops
	      matches <- matches[matches > 0]

	        if (length(matches) == 0) {
			    warning(sprintf("  No matches found between all results and coordinates"))
	          next
		    }

	        # Store filtered coordinates and results
	        coords_filtered <- coords[matches]
		  all_res_df_filtered <- all_res_df[matches > 0, ]

		  # Add resolution column
		  all_res_df_filtered$resolution <- res
		    all_res_df_filtered$resolution_kb <- res_kb

		    all_results_list[[as.character(res)]] <- all_res_df_filtered
		      all_coords_list[[as.character(res)]] <- coords_filtered

		      cat(sprintf("  Matched: %d loops\n\n", length(matches)))
}

if (length(all_results_list) == 0) {
	  stop("ERROR: No all_results loaded. Cannot proceed.")
}

cat(sprintf("✓ Loaded all_results from %d resolutions\n", length(all_results_list)))
cat(sprintf("  Total loops before merging: %d\n\n", sum(sapply(all_results_list, nrow))))

# Merge all loops with overlap removal
cat("Merging all_results with overlap removal...\n")
cat(sprintf("Overlap tolerance: %d kb\n", tolerance_kb))
cat("Priority: Lowest FDR (most significant)\n\n")

# Combine all loops
all_loops_all_df <- bind_rows(all_results_list)
all_loops_all_coords <- Reduce(c, all_coords_list)

cat(sprintf("Total loops to process: %d\n", nrow(all_loops_all_df)))
cat(sprintf("  5kb: %d loops\n", sum(all_loops_all_df$resolution == 5000)))
cat(sprintf("  10kb: %d loops\n", sum(all_loops_all_df$resolution == 10000)))
cat(sprintf("  25kb: %d loops\n\n", sum(all_loops_all_df$resolution == 25000)))

# Build overlap graph (same algorithm as SECTION 2)
cat("Detecting overlaps...\n")

tolerance_bp_all <- tolerance_kb * 1000
n_loops_all <- length(all_loops_all_coords)

# Create adjacency matrix
overlap_graph_all <- matrix(FALSE, nrow = n_loops_all, ncol = n_loops_all)

# Pre-extract coordinates
cat("  Pre-extracting coordinates...\n")
start_time <- Sys.time()

chr1_all_loops <- as.character(seqnames(anchors(all_loops_all_coords, "first")))
start1_all_loops <- start(anchors(all_loops_all_coords, "first"))
chr2_all_loops <- as.character(seqnames(anchors(all_loops_all_coords, "second")))
start2_all_loops <- start(anchors(all_loops_all_coords, "second"))

cat(sprintf("    Extracted coordinates for %d loops (%.1fs)\n",
	                n_loops_all, as.numeric(Sys.time() - start_time)))

# Group by chromosome pairs
cat("  Grouping by chromosome pairs...\n")
chr_pair_key_all <- paste(chr1_all_loops, chr2_all_loops, sep = "_")
chr_pair_groups_all <- split(1:n_loops_all, chr_pair_key_all)

cat(sprintf("    Found %d unique chromosome pairs\n", length(chr_pair_groups_all)))
cat(sprintf("    Starting pairwise comparisons...\n"))

# Vectorized comparisons
start_time <- Sys.time()
comparisons_done <- 0
loops_processed <- 0

for (group_name in names(chr_pair_groups_all)) {
	  group_indices <- chr_pair_groups_all[[group_name]]
  n_group <- length(group_indices)

    if (n_group < 2) {
	        loops_processed <- loops_processed + 1
      next
        }

    for (i_idx in 1:(n_group - 1)) {
	        i <- group_indices[i_idx]
        j_indices <- group_indices[(i_idx + 1):n_group]

	    # Vectorized distance calculation
	    dist1 <- abs(start1_all_loops[i] - start1_all_loops[j_indices])
	    dist2 <- abs(start2_all_loops[i] - start2_all_loops[j_indices])

	        # Find overlaps
	        overlaps <- (dist1 <= tolerance_bp_all) & (dist2 <= tolerance_bp_all)

	        # Update adjacency matrix
	        if (any(overlaps)) {
			      overlap_j <- j_indices[overlaps]
		      overlap_graph_all[i, overlap_j] <- TRUE
		            overlap_graph_all[overlap_j, i] <- TRUE
		          }

		    comparisons_done <- comparisons_done + length(j_indices)
		    loops_processed <- loops_processed + 1

		        # Progress reporting
		        if (loops_processed %% 100 == 0) {
				      elapsed <- as.numeric(Sys.time() - start_time)
		          rate <- loops_processed / elapsed
			        remaining <- (n_loops_all - loops_processed) / rate
			        cat(sprintf("    Processed %d/%d loops (%d comparisons, %.1fs elapsed, ~%.1fs remaining)\r",
					                      loops_processed, n_loops_all, comparisons_done, elapsed, remaining))
				    }
		      }

      loops_processed <- loops_processed + 1
}

elapsed <- as.numeric(Sys.time() - start_time)
cat(sprintf("    Processed %d/%d loops (%d comparisons, %.1fs total)      \n",
	                n_loops_all, n_loops_all, comparisons_done, elapsed))

# Find connected components
cat("\nFinding overlap clusters...\n")

visited_all <- rep(FALSE, n_loops_all)
clusters_all <- list()
cluster_id_all <- 0

for (i in 1:n_loops_all) {
	  if (!visited_all[i]) {
		      cluster_id_all <- cluster_id_all + 1
    current_cluster <- c(i)
        visited_all[i] <- TRUE

        # BFS
        queue <- c(i)

	    while (length(queue) > 0) {
		          current <- queue[1]
	      queue <- queue[-1]

	            neighbors <- which(overlap_graph_all[current, ] & !visited_all)

	            if (length(neighbors) > 0) {
			            visited_all[neighbors] <- TRUE
		            current_cluster <- c(current_cluster, neighbors)
			            queue <- c(queue, neighbors)
			          }
		        }

	    clusters_all[[cluster_id_all]] <- current_cluster
	      }
}

n_singleton_all <- sum(sapply(clusters_all, length) == 1)
n_merged_all <- sum(sapply(clusters_all, length) > 1)

cat(sprintf("  Found %d clusters\n", length(clusters_all)))
cat(sprintf("    Singleton clusters (no overlaps): %d\n", n_singleton_all))
cat(sprintf("    Multi-loop clusters (will merge): %d\n", n_merged_all))
cat(sprintf("    Total loops to merge: %d\n\n", sum(sapply(clusters_all[sapply(clusters_all, length) > 1], length))))

# Select representative loop (lowest FDR) and track comprehensive metadata
cat("Selecting representative loops and tracking metadata...\n")

keep_indices_all <- numeric(length(clusters_all))
overlap_report_all <- data.frame(
				   cluster_id = integer(),
				     n_overlaps = integer(),
				     kept_index = integer(),
				       kept_resolution = character(),
				       kept_fdr = numeric(),
				         n_resolutions = integer(),
				         resolutions_list = character(),
					   merged_indices = character(),
					   merged_resolutions = character(),
					     stringsAsFactors = FALSE
					   )

# For each loop, track multi-resolution metadata
all_loops_all_df$n_resolutions_detected <- NA
all_loops_all_df$resolutions_list <- NA
all_loops_all_df$is_multi_resolution <- FALSE
all_loops_all_df$FDR_5kb <- NA_real_
all_loops_all_df$FDR_10kb <- NA_real_
all_loops_all_df$FDR_25kb <- NA_real_
all_loops_all_df$logFC_5kb <- NA_real_
all_loops_all_df$logFC_10kb <- NA_real_
all_loops_all_df$logFC_25kb <- NA_real_
all_loops_all_df$kept_from_resolution <- NA

for (i in seq_along(clusters_all)) {
	  cluster_indices <- clusters_all[[i]]
  cluster_fdr <- all_loops_all_df$FDR[cluster_indices]

    # Select loop with lowest FDR
    best_idx_in_cluster <- which.min(cluster_fdr)
    best_idx_global <- cluster_indices[best_idx_in_cluster]
      keep_indices_all[i] <- best_idx_global

      # Track which resolutions detected this loop cluster
      cluster_resolutions <- unique(all_loops_all_df$resolution[cluster_indices])
        n_res <- length(cluster_resolutions)
        res_list <- paste(sort(cluster_resolutions / 1000), collapse = ",")

	  # For the kept loop, store multi-resolution metadata
	  all_loops_all_df$n_resolutions_detected[best_idx_global] <- n_res
	  all_loops_all_df$resolutions_list[best_idx_global] <- paste0(sort(cluster_resolutions / 1000), "kb", collapse = ",")
	    all_loops_all_df$is_multi_resolution[best_idx_global] <- n_res > 1
	    all_loops_all_df$kept_from_resolution[best_idx_global] <- all_loops_all_df$resolution_kb[best_idx_global]

	      # Store FDR and logFC from each resolution (if present in cluster)
	      for (idx in cluster_indices) {
		          res <- all_loops_all_df$resolution[idx]
	        if (res == 5000) {
			      all_loops_all_df$FDR_5kb[best_idx_global] <- all_loops_all_df$FDR[idx]
		      all_loops_all_df$logFC_5kb[best_idx_global] <- all_loops_all_df$logFC[idx]
		          } else if (res == 10000) {
				        all_loops_all_df$FDR_10kb[best_idx_global] <- all_loops_all_df$FDR[idx]
		            all_loops_all_df$logFC_10kb[best_idx_global] <- all_loops_all_df$logFC[idx]
			        } else if (res == 25000) {
					      all_loops_all_df$FDR_25kb[best_idx_global] <- all_loops_all_df$FDR[idx]
			          all_loops_all_df$logFC_25kb[best_idx_global] <- all_loops_all_df$logFC[idx]
				      }
		  }

	      # Record in overlap report
	      if (length(cluster_indices) > 1) {
		          overlap_report_all <- rbind(overlap_report_all, data.frame(
										           cluster_id = i,
											         n_overlaps = length(cluster_indices),
											         kept_index = best_idx_global,
												       kept_resolution = sprintf("%dkb", all_loops_all_df$resolution_kb[best_idx_global]),
												       kept_fdr = all_loops_all_df$FDR[best_idx_global],
												             n_resolutions = n_res,
												             resolutions_list = res_list,
													           merged_indices = paste(cluster_indices, collapse = ","),
													           merged_resolutions = paste(sprintf("%dkb", all_loops_all_df$resolution_kb[cluster_indices]), collapse = ","),
														         stringsAsFactors = FALSE
														       ))
	        }
}

cat(sprintf("  Selected %d representative loops\n", length(keep_indices_all)))
cat(sprintf("  Removed %d overlapping loops\n\n", n_loops_all - length(keep_indices_all)))

# Create merged all_results dataset
merged_all_results_df <- all_loops_all_df[keep_indices_all, ]
merged_all_results_coords <- all_loops_all_coords[keep_indices_all]

# Add overlap annotation
merged_all_results_df$n_overlaps <- sapply(clusters_all, length)
merged_all_results_df$is_merged <- merged_all_results_df$n_overlaps > 1

cat(sprintf("✓ Non-redundant all_results set: %d loops\n", nrow(merged_all_results_df)))
cat(sprintf("  From 5kb: %d loops\n", sum(merged_all_results_df$resolution == 5000)))
cat(sprintf("  From 10kb: %d loops\n", sum(merged_all_results_df$resolution == 10000)))
cat(sprintf("  From 25kb: %d loops\n", sum(merged_all_results_df$resolution == 25000)))
cat(sprintf("  Multi-resolution loops: %d (%.1f%%)\n\n",
	                sum(merged_all_results_df$is_multi_resolution),
			            100 * mean(merged_all_results_df$is_multi_resolution)))

# Save merged all_results
write.table(merged_all_results_df,
	                file.path(output_dir, "merged_all_results.tsv"),
			            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("✓ Merged all_results saved: %s\n", file.path(output_dir, "merged_all_results.tsv")))

# Save overlap report
if (nrow(overlap_report_all) > 0) {
	  write.table(overlap_report_all,
		                    file.path(output_dir, "merged_all_results_overlap_report.tsv"),
				                  sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("✓ Overlap report saved: %s\n\n", file.path(output_dir, "merged_all_results_overlap_report.tsv")))
}

# Save RDS object
merged_all_results_coords_with_metadata <- merged_all_results_coords
mcols(merged_all_results_coords_with_metadata) <- merged_all_results_df
saveRDS(merged_all_results_coords_with_metadata,
	        file.path(output_dir, "merged_all_results.rds"))
cat(sprintf("✓ RDS object saved: %s\n\n", file.path(output_dir, "merged_all_results.rds")))

# =============================================================================
# SECTION 2: MERGE LOOPS WITH OVERLAP REMOVAL
# =============================================================================

cat("\n========================================\n")
cat("SECTION 2: Merging Loops (Overlap Removal)\n")
cat("\n========================================\n")

cat(sprintf("Overlap tolerance: %d kb\n", tolerance_kb))
cat("Priority: Lowest FDR (most significant)\n\n")

# Combine all loops into a single data structure
all_loops_df <- bind_rows(final_results_list)
# Use Reduce instead of do.call to properly concatenate GInteractions
# (do.call on named lists creates a list instead of concatenating)
all_loops_coords <- Reduce(c, coords_list)

cat(sprintf("Total loops to process: %d\n", nrow(all_loops_df)))
cat(sprintf("  5kb: %d loops\n", sum(all_loops_df$resolution == 5000)))
cat(sprintf("  10kb: %d loops\n", sum(all_loops_df$resolution == 10000)))
cat(sprintf("  25kb: %d loops\n\n", sum(all_loops_df$resolution == 25000)))

# Build overlap graph
cat("Detecting overlaps...\n")

tolerance_bp <- tolerance_kb * 1000
n_loops <- length(all_loops_coords)

# Create adjacency matrix for overlaps
overlap_graph <- matrix(FALSE, nrow = n_loops, ncol = n_loops)

# OPTIMIZATION 1: Pre-extract all coordinates once (instead of millions of times)
cat("  Pre-extracting coordinates...\n")
start_time <- Sys.time()

chr1_all <- as.character(seqnames(anchors(all_loops_coords, "first")))
start1_all <- start(anchors(all_loops_coords, "first"))
chr2_all <- as.character(seqnames(anchors(all_loops_coords, "second")))
start2_all <- start(anchors(all_loops_coords, "second"))

cat(sprintf("    Extracted coordinates for %d loops (%.1fs)\n",
	    	                n_loops, as.numeric(Sys.time() - start_time)))

# OPTIMIZATION 2: Group by chromosome pairs to reduce comparison space
cat("  Grouping by chromosome pairs...\n")
chr_pair_key <- paste(chr1_all, chr2_all, sep = "_")
chr_pair_groups <- split(1:n_loops, chr_pair_key)

cat(sprintf("    Found %d unique chromosome pairs\n", length(chr_pair_groups)))
cat(sprintf("    Starting pairwise comparisons...\n"))

# OPTIMIZATION 3: Vectorized comparisons within chromosome groups
start_time <- Sys.time()
comparisons_done <- 0
loops_processed <- 0

for (group_name in names(chr_pair_groups)) {
		  group_indices <- chr_pair_groups[[group_name]]
  n_group <- length(group_indices)

      # Skip single-loop groups (no comparisons needed)
      if (n_group < 2) {
	      	        loops_processed <- loops_processed + 1
        next
	        }

      # Compare all pairs within this chromosome group
      for (i_idx in 1:(n_group - 1)) {
	      	        i <- group_indices[i_idx]

              # Get remaining indices in group
              j_indices <- group_indices[(i_idx + 1):n_group]

	      	    # Vectorized distance calculation
	      	    dist1 <- abs(start1_all[i] - start1_all[j_indices])
	      	    dist2 <- abs(start2_all[i] - start2_all[j_indices])

		    	        # Find overlaps (within tolerance on both anchors)
		    	        overlaps <- (dist1 <= tolerance_bp) & (dist2 <= tolerance_bp)

		    	        # Update adjacency matrix (symmetric)
		    	        if (any(overlaps)) {
								      overlap_j <- j_indices[overlaps]
						      overlap_graph[i, overlap_j] <- TRUE
						      		            overlap_graph[overlap_j, i] <- TRUE
						      		          }

						    comparisons_done <- comparisons_done + length(j_indices)
						    loops_processed <- loops_processed + 1

						    		        # Progress reporting
						    		        if (loops_processed %% 100 == 0) {
														      elapsed <- as.numeric(Sys.time() - start_time)
						    		          rate <- loops_processed / elapsed
									  			        remaining <- (n_loops - loops_processed) / rate
									  			        cat(sprintf("    Processed %d/%d loops (%d comparisons, %.1fs elapsed, ~%.1fs remaining)\r",
														      					                      loops_processed, n_loops, comparisons_done, elapsed, remaining))
																	    }
						    		      }

            # Account for last loop in group
            loops_processed <- loops_processed + 1
}

elapsed <- as.numeric(Sys.time() - start_time)
cat(sprintf("    Processed %d/%d loops (%d comparisons, %.1fs total)      \n",
	    	                n_loops, n_loops, comparisons_done, elapsed))

# Find connected components (clusters of overlapping loops)
cat("\nFinding overlap clusters...\n")

visited <- rep(FALSE, n_loops)
clusters <- list()
cluster_id <- 0

for (i in 1:n_loops) {
		  if (!visited[i]) {
			  		      # Start new cluster
			  		      cluster_id <- cluster_id + 1
    current_cluster <- c(i)
            visited[i] <- TRUE

            # Find all connected loops (BFS)
            queue <- c(i)

	    	    while (length(queue) > 0) {
			    		          current <- queue[1]
	    	      queue <- queue[-1]

		      	            # Find neighbors
		      	            neighbors <- which(overlap_graph[current, ] & !visited)

		      	            if (length(neighbors) > 0) {
					    			            visited[neighbors] <- TRUE
				    		            current_cluster <- c(current_cluster, neighbors)
							    			            queue <- c(queue, neighbors)
							    			          }
				    		        }

	    	    clusters[[cluster_id]] <- current_cluster
		    	      }
}

n_singleton <- sum(sapply(clusters, length) == 1)
n_merged <- sum(sapply(clusters, length) > 1)

cat(sprintf("  Found %d clusters\n", length(clusters)))
cat(sprintf("    Singleton clusters (no overlaps): %d\n", n_singleton))
cat(sprintf("    Multi-loop clusters (will merge): %d\n", n_merged))
cat(sprintf("    Total loops to merge: %d\n\n", sum(sapply(clusters[sapply(clusters, length) > 1], length))))

# Select representative loop from each cluster (lowest FDR)
cat("Selecting representative loops (lowest FDR)...\n")

keep_indices <- numeric(length(clusters))
overlap_report <- data.frame(
			     			       cluster_id = integer(),
						       			         n_overlaps = integer(),
						       			         kept_index = integer(),
										 				   kept_resolution = character(),
										 				   kept_fdr = numeric(),
														   				     merged_indices = character(),
														   				     merged_resolutions = character(),
																		     				       stringsAsFactors = FALSE
																		     				     )

for (i in seq_along(clusters)) {
		  cluster_indices <- clusters[[i]]
  cluster_fdr <- all_loops_df$FDR[cluster_indices]

      # Select loop with lowest FDR
      best_idx_in_cluster <- which.min(cluster_fdr)
      best_idx_global <- cluster_indices[best_idx_in_cluster]
            keep_indices[i] <- best_idx_global

            # Record in overlap report
            if (length(cluster_indices) > 1) {
		    	          overlap_report <- rbind(overlap_report, data.frame(
											 								           cluster_id = i,
																				   									         n_overlaps = length(cluster_indices),
																				   									         kept_index = best_idx_global,
																														 										       kept_resolution = sprintf("%dkb", all_loops_df$resolution_kb[best_idx_global]),
																														 										       kept_fdr = all_loops_df$FDR[best_idx_global],
																																								       										             merged_indices = paste(cluster_indices, collapse = ","),
																																								       										             merged_resolutions = paste(sprintf("%dkb", all_loops_df$resolution_kb[cluster_indices]), collapse = ","),
																																																			     											           stringsAsFactors = FALSE
																																																			     											         ))
	            }
}

cat(sprintf("  Selected %d representative loops\n", length(keep_indices)))
cat(sprintf("  Removed %d overlapping loops\n\n", n_loops - length(keep_indices)))

# Create merged dataset
merged_loops_df <- all_loops_df[keep_indices, ]
merged_loops_coords <- all_loops_coords[keep_indices]

# Add overlap annotation
merged_loops_df$n_overlaps <- sapply(clusters, length)
merged_loops_df$is_merged <- merged_loops_df$n_overlaps > 1

cat(sprintf("✓ Non-redundant loop set: %d loops\n", nrow(merged_loops_df)))
cat(sprintf("  From 5kb: %d loops\n", sum(merged_loops_df$resolution == 5000)))
cat(sprintf("  From 10kb: %d loops\n", sum(merged_loops_df$resolution == 10000)))
cat(sprintf("  From 25kb: %d loops\n\n", sum(merged_loops_df$resolution == 25000)))

# Save overlap report
if (nrow(overlap_report) > 0) {
		  write.table(overlap_report,
			      		                    file.path(output_dir, "overlap_report.tsv"),
							    				                  sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("✓ Overlap report saved: %s\n\n", file.path(output_dir, "overlap_report.tsv")))
}

# =============================================================================
# SECTION 3: VALIDATE ALL LOOPS FILES FOR VOLCANO PLOTS
# =============================================================================

cat("\n========================================\n")
cat("SECTION 3: Validating All Loops Files\n")
cat("\n========================================\n")

volcano_status <- data.frame(
			     			       resolution = character(),
						       			         file_exists = logical(),
						       			         n_loops = integer(),
										 				   has_logFC = logical(),
										 				   has_FDR = logical(),
														   				     has_PValue = logical(),
														   				     stringsAsFactors = FALSE
																		     				     )

for (res_kb in c(5, 10, 25)) {
		  all_results_file <- sprintf("outputs/edgeR_results_res_%dkb/primary_analysis/all_results_primary.tsv", res_kb)

  file_exists <- file.exists(all_results_file)

      if (file_exists) {
	      	        df <- read.table(all_results_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, nrows = 1)

        volcano_status <- rbind(volcano_status, data.frame(
							   							       resolution = sprintf("%dkb", res_kb),
														       							             file_exists = TRUE,
														       							             n_loops = NA,  # Will count later
																						     								           has_logFC = "logFC" %in% colnames(df),
																						     								           has_FDR = "FDR" %in% colnames(df),
																															   									         has_PValue = "PValue" %in% colnames(df),
																															   									         stringsAsFactors = FALSE
																																									 										     ))

	          # Count actual loops
	          df_full <- read.table(all_results_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	          volcano_status$n_loops[nrow(volcano_status)] <- nrow(df_full)

		  	      cat(sprintf("%s: ✓ Valid (%d loops)\n", sprintf("%dkb", res_kb), nrow(df_full)))
		  	    } else {
				    		        volcano_status <- rbind(volcano_status, data.frame(
													       									         resolution = sprintf("%dkb", res_kb),
																							 										       file_exists = FALSE,
																							 										       n_loops = NA,
																																	       										             has_logFC = FALSE,
																																	       										             has_FDR = FALSE,
																																												     											           has_PValue = FALSE,
																																												     											           stringsAsFactors = FALSE
																																																								   												       ))

			    	        cat(sprintf("%s: ✗ File not found\n", sprintf("%dkb", res_kb)))
							  }
}

cat("\n")

# Save validation report
write.table(volcano_status,
	    	                file.path(output_dir, "volcano_files_validation.tsv"),
							            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("✓ Validation complete\n"))
cat(sprintf("  Files saved in: outputs/edgeR_results_res_{RES}kb/primary_analysis/all_results_primary.tsv\n\n"))

# =============================================================================
# SECTION 4: LOOP ANCHOR CHARACTERIZATIONS
# =============================================================================

cat("\n========================================\n")
cat("SECTION 4: Loop Anchor Characterizations\n")
cat("\n========================================\n")

# --- 4A. Basic Genomic Features ---

cat("4A. Computing basic genomic features...\n")

# Extract anchor coordinates
anchor1_gr <- anchors(merged_loops_coords, "first")
anchor2_gr <- anchors(merged_loops_coords, "second")

# Calculate midpoints
anchor1_mid <- start(anchor1_gr) + (end(anchor1_gr) - start(anchor1_gr)) / 2
anchor2_mid <- start(anchor2_gr) + (end(anchor2_gr) - start(anchor2_gr)) / 2

# Basic features
merged_loops_df$anchor1_chr <- as.character(seqnames(anchor1_gr))
merged_loops_df$anchor1_start <- start(anchor1_gr)
merged_loops_df$anchor1_end <- end(anchor1_gr)
merged_loops_df$anchor1_midpoint <- anchor1_mid
merged_loops_df$anchor1_size <- end(anchor1_gr) - start(anchor1_gr)

merged_loops_df$anchor2_chr <- as.character(seqnames(anchor2_gr))
merged_loops_df$anchor2_start <- start(anchor2_gr)
merged_loops_df$anchor2_end <- end(anchor2_gr)
merged_loops_df$anchor2_midpoint <- anchor2_mid
merged_loops_df$anchor2_size <- end(anchor2_gr) - start(anchor2_gr)

# Loop distance
merged_loops_df$loop_distance <- anchor2_mid - anchor1_mid

# Distance categories
merged_loops_df$distance_category <- cut(
					 					   merged_loops_df$loop_distance,
										   					     breaks = c(0, 100000, 500000, 1000000, Inf),
										   					     labels = c("<100kb", "100-500kb", "500kb-1Mb", ">1Mb"),
															     					       include.lowest = TRUE
															     					     )

# Chromosome type
merged_loops_df$is_intrachromosomal <- merged_loops_df$anchor1_chr == merged_loops_df$anchor2_chr

cat(sprintf("  Computed features for %d loops\n", nrow(merged_loops_df)))
cat(sprintf("  Distance range: %d bp - %d bp\n",
	    	                min(merged_loops_df$loop_distance),
							            max(merged_loops_df$loop_distance)))
cat(sprintf("  All loops are intrachromosomal: %s\n\n",
	    	                ifelse(all(merged_loops_df$is_intrachromosomal), "Yes", "No")))

# --- 4B. Gene Annotations ---

cat("4B. Annotating with gene information...\n")

# Load gene annotations
if (!is.null(gtf_file) && file.exists(gtf_file)) {
		  cat(sprintf("  Loading GTF: %s\n", gtf_file))
  suppressPackageStartupMessages({
	  	      library(GenomicFeatures)
		      	        })
      txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
} else {
		  cat("  Loading TxDb.Mmusculus.UCSC.mm10.knownGene\n")
  suppressPackageStartupMessages({
	  	      library(TxDb.Mmusculus.UCSC.mm10.knownGene)
		      	          library(org.Mm.eg.db)
		      	        })
      txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
}

# Extract genes
genes <- genes(txdb)

cat(sprintf("  Loaded %d genes\n", length(genes)))

# Get TSS (transcription start sites)
tss <- resize(genes, width = 1, fix = "start")

# Find nearest gene for each anchor
cat("  Finding nearest genes for anchor1...\n")
nearest1 <- distanceToNearest(anchor1_gr, tss)
merged_loops_df$anchor1_nearest_gene_idx <- NA
merged_loops_df$anchor1_distance_to_tss <- NA

if (length(nearest1) > 0) {
		  merged_loops_df$anchor1_nearest_gene_idx[queryHits(nearest1)] <- subjectHits(nearest1)
  merged_loops_df$anchor1_distance_to_tss[queryHits(nearest1)] <- mcols(nearest1)$distance
}

cat("  Finding nearest genes for anchor2...\n")
nearest2 <- distanceToNearest(anchor2_gr, tss)
merged_loops_df$anchor2_nearest_gene_idx <- NA
merged_loops_df$anchor2_distance_to_tss <- NA

if (length(nearest2) > 0) {
		  merged_loops_df$anchor2_nearest_gene_idx[queryHits(nearest2)] <- subjectHits(nearest2)
  merged_loops_df$anchor2_distance_to_tss[queryHits(nearest2)] <- mcols(nearest2)$distance
}

# Add gene IDs
merged_loops_df$anchor1_nearest_gene <- NA
merged_loops_df$anchor2_nearest_gene <- NA

for (i in 1:nrow(merged_loops_df)) {
		  if (!is.na(merged_loops_df$anchor1_nearest_gene_idx[i])) {
			  		      gene_idx <- merged_loops_df$anchor1_nearest_gene_idx[i]
    merged_loops_df$anchor1_nearest_gene[i] <- names(genes)[gene_idx]
          }

  if (!is.na(merged_loops_df$anchor2_nearest_gene_idx[i])) {
	  	      gene_idx <- merged_loops_df$anchor2_nearest_gene_idx[i]
        merged_loops_df$anchor2_nearest_gene[i] <- names(genes)[gene_idx]
	        }
}

# Classify anchors by proximity to TSS
# Promoter: within 2kb of TSS
merged_loops_df$anchor1_is_promoter <- !is.na(merged_loops_df$anchor1_distance_to_tss) &
		                                        merged_loops_df$anchor1_distance_to_tss <= 2000
											merged_loops_df$anchor2_is_promoter <- !is.na(merged_loops_df$anchor2_distance_to_tss) &
																		                                        merged_loops_df$anchor2_distance_to_tss <= 2000

																																# Classify loop type
																																merged_loops_df$loop_type <- case_when(
																																				       														         merged_loops_df$anchor1_is_promoter & merged_loops_df$anchor2_is_promoter ~ "Promoter-Promoter",
																																																			 															   merged_loops_df$anchor1_is_promoter | merged_loops_df$anchor2_is_promoter ~ "Promoter-Enhancer",
																																																			 															   TRUE ~ "Enhancer-Enhancer"
																																																																		   															   )

										cat(sprintf("\n  Gene annotation complete\n"))
																				cat(sprintf("  Loops with genes near both anchors: %d (%.1f%%)\n",
																					    											                sum(!is.na(merged_loops_df$anchor1_nearest_gene) & !is.na(merged_loops_df$anchor2_nearest_gene)),
																																															            100 * mean(!is.na(merged_loops_df$anchor1_nearest_gene) & !is.na(merged_loops_df$anchor2_nearest_gene))))
																				cat(sprintf("\n  Loop type classification:\n"))
																														cat(sprintf("    Promoter-Promoter: %d (%.1f%%)\n",
																															    											                sum(merged_loops_df$loop_type == "Promoter-Promoter"),
																																																									            100 * mean(merged_loops_df$loop_type == "Promoter-Promoter")))
																														cat(sprintf("    Promoter-Enhancer: %d (%.1f%%)\n",
																															    											                sum(merged_loops_df$loop_type == "Promoter-Enhancer"),
																																																									            100 * mean(merged_loops_df$loop_type == "Promoter-Enhancer")))
																																								cat(sprintf("    Enhancer-Enhancer: %d (%.1f%%)\n\n",
																																									    											                sum(merged_loops_df$loop_type == "Enhancer-Enhancer"),
																																																																			            100 * mean(merged_loops_df$loop_type == "Enhancer-Enhancer")))

																																								# --- 4C. Characterization Summaries ---

																																								cat("4C. Generating characterization summaries and plots...\n")

																																																		# Summary statistics
																																																		summary_text <- c(
																																																				  												    "=========================================",
																																																																    												      "LOOP ANCHOR CHARACTERIZATION SUMMARY",
																																																																    												      "=========================================",
																																																																												      												        sprintf("Total non-redundant loops: %d", nrow(merged_loops_df)),
																																																																												      												        sprintf("  From 5kb resolution: %d", sum(merged_loops_df$resolution == 5000)),
																																																																																																						  sprintf("  From 10kb resolution: %d", sum(merged_loops_df$resolution == 10000)),
																																																																																																						  sprintf("  From 25kb resolution: %d", sum(merged_loops_df$resolution == 25000)),
																																																																																																						  													    "",
																																																																																																						  													    "Loop Distance Distribution:",
																																																																																																																			    													      sprintf("  <100kb: %d (%.1f%%)",
																																																																																																																																		  														                sum(merged_loops_df$distance_category == "<100kb", na.rm = TRUE),
																																																																																																																																																																		          100 * mean(merged_loops_df$distance_category == "<100kb", na.rm = TRUE)),
																																																				  												    sprintf("  100-500kb: %d (%.1f%%)",
																																																																	      													              sum(merged_loops_df$distance_category == "100-500kb", na.rm = TRUE),
																																																																															      														                100 * mean(merged_loops_df$distance_category == "100-500kb", na.rm = TRUE)),
																																																				  												    sprintf("  500kb-1Mb: %d (%.1f%%)",
																																																																	      													              sum(merged_loops_df$distance_category == "500kb-1Mb", na.rm = TRUE),
																																																																															      														                100 * mean(merged_loops_df$distance_category == "500kb-1Mb", na.rm = TRUE)),
																																																				  												    sprintf("  >1Mb: %d (%.1f%%)",
																																																																	      													              sum(merged_loops_df$distance_category == ">1Mb", na.rm = TRUE),
																																																																															      														                100 * mean(merged_loops_df$distance_category == ">1Mb", na.rm = TRUE)),
																																																				  												    "",
																																																																    												      "Gene Association:",
																																																																    												      sprintf("  Anchor1 with nearby gene (<50kb): %d (%.1f%%)",
																																																																														  													                sum(!is.na(merged_loops_df$anchor1_distance_to_tss) &
																																																																																													      															                  merged_loops_df$anchor1_distance_to_tss < 50000),
																																																																														  													                100 * mean(!is.na(merged_loops_df$anchor1_distance_to_tss) &
																																																																																														     																                        merged_loops_df$anchor1_distance_to_tss < 50000)),
																																																				  												    sprintf("  Anchor2 with nearby gene (<50kb): %d (%.1f%%)",
																																																																	      													              sum(!is.na(merged_loops_df$anchor2_distance_to_tss) &
																																																																																																                merged_loops_df$anchor2_distance_to_tss < 50000),
																																																																	      													              100 * mean(!is.na(merged_loops_df$anchor2_distance_to_tss) &
																																																																																	       																                      merged_loops_df$anchor2_distance_to_tss < 50000)),
																																																				  												    "",
																																																																    												      "Loop Type Classification:",
																																																																    												      sprintf("  Promoter-Promoter: %d (%.1f%%)",
																																																																														  													                sum(merged_loops_df$loop_type == "Promoter-Promoter"),
																																																																																																												          100 * mean(merged_loops_df$loop_type == "Promoter-Promoter")),
																																																				  												    sprintf("  Promoter-Enhancer: %d (%.1f%%)",
																																																																	      													              sum(merged_loops_df$loop_type == "Promoter-Enhancer"),
																																																																															      														                100 * mean(merged_loops_df$loop_type == "Promoter-Enhancer")),
																																																				  												    sprintf("  Enhancer-Enhancer: %d (%.1f%%)",
																																																																	      													              sum(merged_loops_df$loop_type == "Enhancer-Enhancer"),
																																																																															      														                100 * mean(merged_loops_df$loop_type == "Enhancer-Enhancer")),
																																																				  												    "",
																																																																    												      "Differential Loops:",
																																																																    												      sprintf("  Up in mutant: %d (%.1f%%)",
																																																																														  													                sum(merged_loops_df$logFC > 0),
																																																																																																												          100 * mean(merged_loops_df$logFC > 0)),
																																																				  												    sprintf("  Down in mutant: %d (%.1f%%)",
																																																																	      													              sum(merged_loops_df$logFC < 0),
																																																																															      														                100 * mean(merged_loops_df$logFC < 0)),
																																																				  												    "",
																																																																    												      "========================================="
																																																																    												    )

																																																		writeLines(summary_text, file.path(output_dir, "characterization_summary.txt"))

																																																												cat("  ✓ Summary saved\n")

																																																												# Generate plots

																																																												# Plot 1: Distance distribution
																																																												p1 <- ggplot(merged_loops_df, aes(x = loop_distance / 1e6)) +
																																																																								  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
																																																																								  											    scale_x_continuous(trans = "log10") +
																																																																																			    											      labs(
																																																																																															       												       title = "Distribution of Loop Distances",
																																																																																																											       												           x = "Loop Distance (Mb, log10)",
																																																																																																											       												           y = "Count"
																																																																																																																								   													     ) +
  theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "plots", "distance_distribution.pdf"), p1, width = 8, height = 6)
cat("  ✓ Distance distribution plot saved\n")

# Plot 2: Chromosome distribution
chr_counts <- merged_loops_df %>%
		  group_by(anchor1_chr) %>%
		  	    summarise(n = n()) %>%
			    	      arrange(desc(n))

			            p2 <- ggplot(chr_counts, aes(x = reorder(anchor1_chr, -n), y = n)) +
					    	        geom_bar(stat = "identity", fill = "steelblue", color = "black", alpha = 0.7) +
									  labs(
									       		           title = "Loop Distribution by Chromosome",
												   			       x = "Chromosome",
												   			       y = "Number of Loops"
															       			         ) +
  theme_minimal() +
      theme(
	    	      plot.title = element_text(hjust = 0.5, face = "bold"),
		      	          axis.text.x = element_text(angle = 45, hjust = 1)
		      	        )

    ggsave(file.path(output_dir, "plots", "chromosome_distribution.pdf"), p2, width = 10, height = 6)
      cat("  ✓ Chromosome distribution plot saved\n")

      # Plot 3: Gene proximity distribution
      gene_proximity_df <- data.frame(
				      				    anchor = rep(c("Anchor1", "Anchor2"), each = nrow(merged_loops_df)),
								    				      distance_to_tss = c(merged_loops_df$anchor1_distance_to_tss,
															      							                        merged_loops_df$anchor2_distance_to_tss)
								    				    ) %>%
          filter(!is.na(distance_to_tss))

  p3 <- ggplot(gene_proximity_df, aes(x = distance_to_tss / 1000, fill = anchor)) +
	  	  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
		  	    scale_x_log10() +
			    	      labs(
					       		       title = "Distance to Nearest Gene TSS",
							       		           x = "Distance to TSS (kb, log10)",
							       		           y = "Count",
										   			       fill = "Anchor"
										   			     ) +
  theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "plots", "gene_proximity_distribution.pdf"), p3, width = 8, height = 6)
cat("  ✓ Gene proximity distribution plot saved\n")

# Plot 4: Loop type classification
loop_type_counts <- merged_loops_df %>%
		  group_by(loop_type, direction) %>%
		  	    summarise(n = n(), .groups = "drop")

		        p4 <- ggplot(loop_type_counts, aes(x = loop_type, y = n, fill = direction)) +
					      geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.7) +
					      	        scale_fill_manual(values = c("up_in_mutant" = "#d73027", "down_in_mutant" = "#4575b4")) +
									  labs(
									       		           title = "Loop Type Classification by Direction",
												   			       x = "Loop Type",
												   			       y = "Number of Loops",
															       			           fill = "Direction"
															       			         ) +
  theme_minimal() +
      theme(
	    	      plot.title = element_text(hjust = 0.5, face = "bold"),
		      	          axis.text.x = element_text(angle = 45, hjust = 1)
		      	        )

    ggsave(file.path(output_dir, "plots", "loop_type_classification.pdf"), p4, width = 8, height = 6)
      cat("  ✓ Loop type classification plot saved\n\n")

      # =============================================================================
      # SECTION 5: EXPORT TO MULTIPLE FORMATS
      # =============================================================================

      cat("\n========================================\n")
        cat("SECTION 5: Exporting Results\n")
        cat("\n========================================\n")

	  # --- 5A. TSV Export ---

	  cat("5A. Exporting TSV files...\n")

	  # Main results file with all features
	  write.table(merged_loops_df,
		      	                  file.path(output_dir, "characterized_loops.tsv"),
					  			              sep = "\t", quote = FALSE, row.names = FALSE)
	    cat(sprintf("  ✓ %s (%d loops)\n",
				                  "characterized_loops.tsv", nrow(merged_loops_df)))

	    # Non-redundant loops (basic info)
	    non_redundant_df <- merged_loops_df %>%
		    	    dplyr::select(loop_id, chr1, start1, end1, chr2, start2, end2,
					      			                  logFC, logCPM, FDR, PValue, significant,
										  					                  resolution_kb, n_overlaps, loop_distance, loop_type,
										  					                  direction, category)

	      write.table(non_redundant_df,
			  	                  file.path(output_dir, "non_redundant_loops.tsv"),
						  			              sep = "\t", quote = FALSE, row.names = FALSE)
	        cat(sprintf("  ✓ %s (%d loops)\n",
			    	                  "non_redundant_loops.tsv", nrow(non_redundant_df)))

	        # Anchor annotations (per-anchor detail)
	        anchor1_df <- merged_loops_df %>%
				    dplyr::select(loop_id, anchor = anchor1_chr,
						  			                  start = anchor1_start, end = anchor1_end, midpoint = anchor1_midpoint,
											  					                  nearest_gene = anchor1_nearest_gene,
											  					                  distance_to_tss = anchor1_distance_to_tss,
																		  							                  is_promoter = anchor1_is_promoter) %>%
		    mutate(anchor_num = 1)

	    anchor2_df <- merged_loops_df %>%
		    	  dplyr::select(loop_id, anchor = anchor2_chr,
					    			                start = anchor2_start, end = anchor2_end, midpoint = anchor2_midpoint,
															                nearest_gene = anchor2_nearest_gene,
															                distance_to_tss = anchor2_distance_to_tss,
																								                is_promoter = anchor2_is_promoter) %>%
	      mutate(anchor_num = 2)

      anchor_annotations <- bind_rows(anchor1_df, anchor2_df)

      write.table(anchor_annotations,
		  	                file.path(output_dir, "anchor_annotations.tsv"),
								            sep = "\t", quote = FALSE, row.names = FALSE)
      cat(sprintf("  ✓ %s (%d anchors)\n\n",
		  	                "anchor_annotations.tsv", nrow(anchor_annotations)))

      # --- 5B. RDS Export ---

      cat("5B. Exporting RDS files...\n")

      # Add metadata to GInteractions
      mcols(merged_loops_coords) <- merged_loops_df

      saveRDS(merged_loops_coords,
	      	        file.path(output_dir, "non_redundant_loops.rds"))
      cat(sprintf("  ✓ %s (GInteractions object)\n\n",
		  	                "non_redundant_loops.rds"))

      # --- 5C. BEDPE Export ---

      cat("5C. Exporting BEDPE files...\n")

      # Create BEDPE format
      bedpe_df <- data.frame(
			     		         chr1 = merged_loops_df$chr1,
						 			   start1 = merged_loops_df$start1,
						 			   end1 = merged_loops_df$end1,
									   			     chr2 = merged_loops_df$chr2,
									   			     start2 = merged_loops_df$start2,
												     			       end2 = merged_loops_df$end2,
												     			       name = merged_loops_df$loop_id,
															       			         score = -log10(merged_loops_df$FDR),  # For visualization
															       			         strand1 = ".",
																			 				   strand2 = ".",
																			 				   color = ifelse(merged_loops_df$logFC > 0, "255,0,0", "0,0,255")  # Red=up, Blue=down
																							   				   )

      write.table(bedpe_df,
		  	                file.path(output_dir, "non_redundant_loops.bedpe"),
								            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      cat(sprintf("  ✓ %s (for Juicebox/IGV)\n\n",
		  	                "non_redundant_loops.bedpe"))

      # =============================================================================
      # FINAL SUMMARY
      # =============================================================================

      cat("\n========================================\n")
      cat("DOWNSTREAM ANALYSIS COMPLETE\n")
      cat("\n========================================\n")

      cat("Output directory: outputs/merged_loops/\n\n")

      cat("Files generated:\n\n")

      cat("Main results:\n")
      cat("  1. non_redundant_loops.tsv\n")
      cat("     → Merged loop list with basic info (%d loops)\n", nrow(non_redundant_df))
      cat("\n")
      cat("  2. characterized_loops.tsv\n")
      cat("     → Complete characterization with all features (%d loops)\n", nrow(merged_loops_df))
      cat("\n")
      cat("  3. anchor_annotations.tsv\n")
      cat("     → Per-anchor gene annotations (%d anchors)\n", nrow(anchor_annotations))
      cat("\n")

      if (nrow(overlap_report) > 0) {
	      	  cat("  4. overlap_report.tsv\n")
        cat(sprintf("     → Details of merged overlapping loops (%d clusters)\n", nrow(overlap_report)))
	    cat("\n")
      }

      cat("  5. non_redundant_loops.rds\n")
      cat("     → GInteractions object for R analysis\n")
      cat("\n")

      cat("  6. non_redundant_loops.bedpe\n")
      cat("     → BEDPE format for genome browser visualization\n")
      cat("\n")

      cat("Quality control:\n")
      cat("  7. volcano_files_validation.tsv\n")
      cat("     → Validation of all_results files for volcano plots\n")
      cat("\n")

      cat("  8. characterization_summary.txt\n")
      cat("     → Summary statistics report\n")
      cat("\n")

      cat("Plots (outputs/merged_loops/plots/):\n")
      cat("  9. distance_distribution.pdf\n")
      cat("  10. chromosome_distribution.pdf\n")
      cat("  11. gene_proximity_distribution.pdf\n")
      cat("  12. loop_type_classification.pdf\n")
      cat("\n")

      cat("Summary statistics:\n")
      cat(sprintf("  Input loops: %d (from 3 resolutions)\n", n_loops))
      cat(sprintf("  Non-redundant loops: %d\n", nrow(merged_loops_df)))
      cat(sprintf("  Loops removed (overlaps): %d\n", n_loops - nrow(merged_loops_df)))
      cat(sprintf("  Promoter-involved loops: %d (%.1f%%)\n",
		  	                sum(merged_loops_df$loop_type != "Enhancer-Enhancer"),
								            100 * mean(merged_loops_df$loop_type != "Enhancer-Enhancer")))
      cat("\n")

      cat("\n========================================\n")
      cat("\n")


