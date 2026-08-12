# matrix_exploration.R

library(mariner)
library(InteractionSet)
library(HDF5Array)

# Load the HDF5-backed extracted matrices
pixels <- loadHDF5SummarizedExperiment(
  dir = "outputs/full",
  prefix = "05_extracted"
)

# Examine specific loops in detail
examine_loop <- function(loop_id, pixels) {
  # Get the 5x5 matrices for both samples
  count_array <- counts(pixels)
  
  ctrl_matrix <- as.matrix(count_array[,,loop_id, 1])
  mut_matrix <- as.matrix(count_array[,,loop_id, 2])
  
  # Plot both matrices side by side
  par(mfrow = c(1, 3), mar = c(2,2,3,1))
  
  # Control matrix
  image(ctrl_matrix, main = paste("Loop", loop_id, "- Control"),
        col = viridis(100), axes = FALSE)
  text(expand.grid(seq(0,1,0.25), seq(0,1,0.25)),
       labels = round(as.vector(ctrl_matrix), 1))
  
  # Mutant matrix
  image(mut_matrix, main = paste("Loop", loop_id, "- Mutant"),
        col = viridis(100), axes = FALSE)
  text(expand.grid(seq(0,1,0.25), seq(0,1,0.25)),
       labels = round(as.vector(mut_matrix), 1))
  
  # Difference matrix
  diff_matrix <- mut_matrix - ctrl_matrix
  image(diff_matrix, main = "Difference (Mut - Ctrl)",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        axes = FALSE)
  text(expand.grid(seq(0,1,0.25), seq(0,1,0.25)),
       labels = round(as.vector(diff_matrix), 1))
  
  # Return statistics
  return(list(
    ctrl_sum = sum(ctrl_matrix, na.rm = TRUE),
    mut_sum = sum(mut_matrix, na.rm = TRUE),
    ctrl_center = ctrl_matrix[3,3],
    mut_center = mut_matrix[3,3],
    shift_detected = which.max(mut_matrix) != which.max(ctrl_matrix)
  ))
}

# Examine loops with biggest changes
top_changes <- order(abs(M), decreasing = TRUE)[1:5]
for (i in top_changes) {
  stats <- examine_loop(i, pixels)
  print(paste("Loop", i, "- Shift detected:", stats$shift_detected))
}
