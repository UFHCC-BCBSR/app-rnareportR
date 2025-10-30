# filter_genes.R - Updated with NOISeq option
library(edgeR)
library(NOISeq)

filter_genes <- function(count_matrix, metadata, condition_var,
                         method = "edgeR",
                         # edgeR parameters
                         min_count = 10, min_prop = 0.7,
                         # NOISeq parameters
                         noiseq_method = 1, cv_cutoff = 100, cpm = 1, p_adj = "fdr") {
  
  # Ensure parameters are numeric (in case from command line)
  min_count <- as.numeric(min_count)
  min_prop <- as.numeric(min_prop)
  noiseq_method <- as.numeric(noiseq_method)
  cv_cutoff <- as.numeric(cv_cutoff)
  cpm <- as.numeric(cpm)
  
  # Align samples
  common_samples <- intersect(colnames(count_matrix), as.character(metadata[,1]))
  count_matrix <- count_matrix[, common_samples, drop = FALSE]
  metadata_aligned <- metadata[match(common_samples, as.character(metadata[,1])), , drop = FALSE]
  
  # Create group factor
  group <- factor(metadata_aligned[[condition_var]])
  
  # Create DGEList
  y <- DGEList(counts = count_matrix, group = group)
  
  if (method == "edgeR") {
    # edgeR filtering
    keep <- filterByExpr(y, group = group, min.count = min_count, min.prop = min_prop)
    filtered_matrix <- count_matrix[keep, , drop = FALSE]
    
    filter_params <- list(
      min_count = min_count,
      min_prop = min_prop
    )
    method_name <- "edgeR filterByExpr"
    
  } else if (method == "NOISeq") {
    # NOISeq filtering
    filtered_counts <- filtered.data(
      count_matrix,
      factor = group,
      norm = FALSE,
      depth = NULL,
      method = noiseq_method,
      cv.cutoff = cv_cutoff,
      cpm = cpm,
      p.adj = p_adj
    )
    
    # Determine which genes were kept
    keep <- rownames(count_matrix) %in% rownames(filtered_counts)
    filtered_matrix <- filtered_counts
    
    filter_params <- list(
      noiseq_method = noiseq_method,
      cv_cutoff = cv_cutoff,
      cpm = cpm,
      p_adj = p_adj
    )
    method_name <- "NOISeq filtered.data"
    
  } else {
    stop("Method must be either 'edgeR' or 'NOISeq'")
  }
  
  return(list(
    filtered_matrix = filtered_matrix,
    keep_genes = keep,
    n_original = nrow(count_matrix),
    n_filtered = nrow(filtered_matrix),
    n_removed = sum(!keep),
    prop_retained = sum(keep) / nrow(count_matrix),
    filter_params = filter_params,
    method = method_name
  ))
}