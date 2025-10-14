# filter_genes.R - Simple version with the two key parameters
library(edgeR)

filter_genes_edger <- function(count_matrix, metadata, condition_var,
                               min_count = 10, min_prop = 0.7) {
  
  # Ensure parameters are numeric (in case from command line)
  min_count <- as.numeric(min_count)
  min_prop <- as.numeric(min_prop)
  
  # Align samples
  common_samples <- intersect(colnames(count_matrix), as.character(metadata[,1]))
  count_matrix <- count_matrix[, common_samples, drop = FALSE]
  metadata_aligned <- metadata[match(common_samples, as.character(metadata[,1])), , drop = FALSE]
  
  # Create group factor
  group <- factor(metadata_aligned[[condition_var]])
  
  # Create DGEList
  y <- DGEList(counts = count_matrix, group = group)
  
  # Apply filterByExpr with user-specified min.count and min.prop
  # Let edgeR handle min.total.count and large.n automatically
  keep <- filterByExpr(y, group = group, min.count = min_count, min.prop = min_prop)
  
  # Filter
  filtered_matrix <- count_matrix[keep, , drop = FALSE]
  
  return(list(
    filtered_matrix = filtered_matrix,
    keep_genes = keep,
    n_original = nrow(count_matrix),
    n_filtered = nrow(filtered_matrix),
    n_removed = sum(!keep),
    prop_retained = sum(keep) / nrow(count_matrix),
    filter_params = list(
      min_count = min_count,
      min_prop = min_prop
    ),
    method = "edgeR filterByExpr"
  ))
}