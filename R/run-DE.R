# run_DE.R - Differential Expression Analysis Functions
library(DESeq2)
library(edgeR)
library(limma)
library(dplyr)

#' Standardize differential expression results columns
standardize_deg_columns <- function(res_df, method, gene_annotations = NULL) {
  # Ensure Gene column exists and is first
  if (!"Gene" %in% colnames(res_df)) {
    if (!is.null(rownames(res_df))) {
      res_df$Gene <- rownames(res_df)
    } else {
      res_df$Gene <- paste0("Gene_", 1:nrow(res_df))
    }
  }
  
  # Add gene_symbol from annotations if available
  if (!is.null(gene_annotations) && "gene_symbol" %in% colnames(gene_annotations)) {
    res_df$gene_symbol <- gene_annotations[res_df$Gene, "gene_symbol"]
  }

  # Create method-specific output formats
  if (method == "deseq2") {
    final_df <- data.frame(
      SYMBOL = if("gene_symbol" %in% colnames(res_df)) res_df$gene_symbol else rep(NA, nrow(res_df)),
      log2FoldChange = res_df$log2FoldChange,
      baseMean = res_df$baseMean,
      lfcSE = res_df$lfcSE,
      stat = if("stat" %in% colnames(res_df)) res_df$stat else res_df$log2FoldChange / (res_df$lfcSE + 0.001),
      pvalue = res_df$pvalue,
      adj.P.value = res_df$padj,
      contrast = res_df$contrast,
      ensembleID = res_df$Gene,
      stringsAsFactors = FALSE
    )
  } else if (method == "edger_GLM") {
    final_df <- data.frame(
      SYMBOL = if("gene_symbol" %in% colnames(res_df)) res_df$gene_symbol else rep(NA, nrow(res_df)),
      logFC = res_df$logFC,
      logCPM = res_df$logCPM,
      P.Value = res_df$PValue,
      adj.P.value = res_df$FDR,
      contrast = res_df$contrast,
      ensembleID = res_df$Gene,
      stringsAsFactors = FALSE
    )
  } else {
    # limma-voom
    final_df <- data.frame(
      SYMBOL = if("gene_symbol" %in% colnames(res_df)) res_df$gene_symbol else rep(NA, nrow(res_df)),
      logFC = res_df$logFC,
      AveExpr = res_df$AveExpr,
      t = res_df$t,
      P.value = res_df$P.Value,
      adj.P.value = res_df$adj.P.Val,
      B = res_df$B,
      contrast = res_df$contrast,
      ensembleID = res_df$Gene,
      stringsAsFactors = FALSE
    )
  }
  
  return(final_df)
}

#' Run DESeq2 differential expression analysis
run_deseq2_analysis <- function(count_matrix, metadata, condition_var, batch_var = NULL,
                                alpha = 0.05, lfc_threshold = 0, shrinkage = TRUE, 
                                gene_annotations = NULL, paired_design = FALSE) {
  tryCatch({
    # Ensure count matrix is integer
    count_matrix <- round(count_matrix)
    
    # Align samples
    common_samples <- intersect(colnames(count_matrix), as.character(metadata[,1]))
    count_matrix <- count_matrix[, common_samples, drop = FALSE]
    metadata <- metadata[match(common_samples, as.character(metadata[,1])), , drop = FALSE]
    rownames(metadata) <- as.character(metadata[,1])
    
    condition <- factor(metadata[[condition_var]])
    condition_levels <- levels(condition)
    n_levels <- length(condition_levels)
    
    # Create design formula
    if (!is.null(batch_var) && batch_var %in% colnames(metadata)) {
      batch <- factor(metadata[[batch_var]])
      design_formula <- as.formula("~ batch + condition")
      colData <- DataFrame(condition = condition, batch = batch)
      if (paired_design) {
        cat("DESeq2: Using paired design with", batch_var, "as blocking factor\n")
      }
    } else {
      design_formula <- as.formula("~ condition")
      colData <- DataFrame(condition = condition)
    }
    
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(
      countData = count_matrix,
      colData = colData,
      design = design_formula
    )
    
    # Run DESeq2
    dds <- DESeq(dds, quiet = TRUE)
    
    if (n_levels == 2) {
      # Simple pairwise comparison
      res <- results(dds, alpha = alpha, lfcThreshold = lfc_threshold)
      
      if (shrinkage) {
        tryCatch({
          res <- lfcShrink(dds, coef = 2, res = res, type = "apeglm", quiet = TRUE)
        }, error = function(e) {
          cat("Shrinkage failed, using unshrunk results:", e$message, "\n")
        })
      }
      
      res_df <- as.data.frame(res)
      res_df <- res_df[complete.cases(res_df[, c("pvalue", "padj")]), ]
      res_df$contrast <- paste0(condition_levels[2], "_vs_", condition_levels[1])
      res_df <- standardize_deg_columns(res_df, "deseq2", gene_annotations)
      res_df <- res_df[order(res_df$adj.P.value), ]
      comparison <- paste(condition_levels, collapse = " vs ")
      
      method_name <- if (paired_design) "DESeq2 (paired design)" else "DESeq2"
    } else {
      # Multiple contrasts - all pairwise
      all_results <- list()
      for(i in 1:(n_levels-1)) {
        for(j in (i+1):n_levels) {
          contrast_name <- paste0(condition_levels[j], "_vs_", condition_levels[i])
          res <- results(dds, contrast = c("condition", condition_levels[j], condition_levels[i]),
                         alpha = alpha, lfcThreshold = lfc_threshold)
          
          if (shrinkage) {
            tryCatch({
              res <- lfcShrink(dds, contrast = c("condition", condition_levels[j], condition_levels[i]),
                               res = res, type = "apeglm", quiet = TRUE)
            }, error = function(e) {
              cat("Shrinkage failed for", contrast_name, ":", e$message, "\n")
            })
          }
          
          res_df <- as.data.frame(res)
          res_df <- res_df[complete.cases(res_df[, c("pvalue", "padj")]), ]
          res_df$contrast <- contrast_name
          res_df <- standardize_deg_columns(res_df, "deseq2", gene_annotations)
          all_results[[contrast_name]] <- res_df
        }
      }
      res_df <- do.call(rbind, all_results)
      rownames(res_df) <- NULL
      comparison <- paste(names(all_results), collapse = ", ")
      method_name <- "DESeq2 (multiple contrasts)"
    }
    
    return(list(
      deg_table = res_df,
      method = method_name,
      n_total = nrow(count_matrix),
      n_significant = sum(res_df$adj.P.value < alpha, na.rm = TRUE),
      dds_object = dds,
      design_formula = design_formula,
      comparison = comparison,
      error = FALSE
    ))
  }, error = function(e) {
    return(list(
      error = TRUE,
      message = paste("DESeq2 analysis failed:", e$message),
      method = "DESeq2"
    ))
  })
}

#' Run edgeR GLM differential expression analysis
run_edger_glm_analysis <- function(count_matrix, metadata, condition_var, batch_var = NULL,
                                   alpha = 0.05, lfc_threshold = 0, gene_annotations = NULL,
                                   paired_design = FALSE) {
  # Ensure count matrix is integer
  count_matrix <- round(count_matrix)
  
  # Align samples
  common_samples <- intersect(colnames(count_matrix), as.character(metadata[,1]))
  count_matrix <- count_matrix[, common_samples, drop = FALSE]
  metadata <- metadata[match(common_samples, as.character(metadata[,1])), , drop = FALSE]
  
  condition <- factor(metadata[[condition_var]])
  condition_levels <- levels(condition)
  n_levels <- length(condition_levels)
  
  # Create design matrix
  if (!is.null(batch_var) && batch_var %in% colnames(metadata)) {
    if (paired_design) {
      cat("edgeR: Using paired design with", batch_var, "as blocking factor\n")
    }
    if (n_levels == 2) {
      design <- model.matrix(as.formula(paste("~", batch_var, "+", condition_var)), data = metadata)
    } else {
      design <- model.matrix(as.formula(paste("~ 0 +", condition_var, "+", batch_var)), data = metadata)
      colnames(design) <- gsub(condition_var, "", colnames(design))
    }
  } else {
    if (n_levels == 2) {
      design <- model.matrix(as.formula(paste("~", condition_var)), data = metadata)
    } else {
      design <- model.matrix(as.formula(paste("~ 0 +", condition_var)), data = metadata)
      colnames(design) <- gsub(condition_var, "", colnames(design))
    }
  }
  
  # Create DGEList
  y <- DGEList(counts = count_matrix)
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  
  # Fit GLM
  fit <- glmQLFit(y, design)
  
  if (n_levels == 2) {
    # Simple comparison
    qlf <- glmQLFTest(fit, coef = ncol(design))
    res_df <- topTags(qlf, n = Inf)$table
    res_df$Gene <- rownames(res_df)
    res_df$contrast <- paste0(condition_levels[2], "_vs_", condition_levels[1])
    res_df <- standardize_deg_columns(res_df, "edger_GLM", gene_annotations)
    comparison <- paste(condition_levels, collapse = " vs ")
    
    method_name <- if (paired_design) "edgeR GLM (paired design)" else "edgeR (GLM)"
  } else {
    # Multiple contrasts
    contrasts_list <- list()
    for(i in 1:(n_levels-1)) {
      for(j in (i+1):n_levels) {
        contrast_name <- paste0(condition_levels[j], "_vs_", condition_levels[i])
        contrast_string <- paste0(condition_levels[j], "-", condition_levels[i])
        contrasts_list[[contrast_name]] <- contrast_string
      }
    }
    
    contrast_matrix <- makeContrasts(contrasts = contrasts_list, levels = design)
    all_results <- list()
    
    for(i in 1:ncol(contrast_matrix)) {
      contrast_name <- colnames(contrast_matrix)[i]
      qlf <- glmQLFTest(fit, contrast = contrast_matrix[,i])
      res <- topTags(qlf, n = Inf)$table
      res$contrast <- contrast_name
      res$Gene <- rownames(res)
      res <- standardize_deg_columns(res, "edger_GLM", gene_annotations)
      all_results[[contrast_name]] <- res
    }
    
    res_df <- do.call(rbind, all_results)
    rownames(res_df) <- NULL
    comparison <- paste(names(contrasts_list), collapse = ", ")
    method_name <- "edgeR (GLM, multiple contrasts)"
  }
  
  return(list(
    deg_table = res_df,
    method = method_name,
    n_total = nrow(count_matrix),
    n_significant = sum(res_df$adj.P.value < alpha, na.rm = TRUE),
    dge_object = y,
    design_matrix = design,
    comparison = comparison
  ))
}

#' Run limma-voom differential expression analysis with optional paired design
run_limma_voom_analysis <- function(count_matrix, metadata, condition_var, batch_var = NULL,
                                    alpha = 0.05, lfc_threshold = 0, gene_annotations = NULL,
                                    paired_design = FALSE) {
  # Ensure count matrix is integer
  count_matrix <- round(count_matrix)
  
  # Align samples
  common_samples <- intersect(colnames(count_matrix), as.character(metadata[,1]))
  count_matrix <- count_matrix[, common_samples, drop = FALSE]
  metadata <- metadata[match(common_samples, as.character(metadata[,1])), , drop = FALSE]
  
  condition <- factor(metadata[[condition_var]])
  condition_levels <- levels(condition)
  n_levels <- length(condition_levels)
  
  # Determine if we should use duplicateCorrelation
  use_duplicate_correlation <- FALSE
  block <- NULL
  
  if (paired_design && !is.null(batch_var) && batch_var %in% colnames(metadata)) {
    block <- factor(metadata[[batch_var]])
    use_duplicate_correlation <- TRUE
    cat("limma-voom: Using duplicateCorrelation for paired design with", batch_var, "as random effect\n")
  }
  
  # Create design matrix
  if (use_duplicate_correlation) {
    # For paired data: simple design, use block in duplicateCorrelation
    if (n_levels == 2) {
      design <- model.matrix(as.formula(paste("~", condition_var)), data = metadata)
    } else {
      design <- model.matrix(as.formula(paste("~ 0 +", condition_var)), data = metadata)
      colnames(design) <- gsub(condition_var, "", colnames(design))
    }
  } else if (!is.null(batch_var) && batch_var %in% colnames(metadata)) {
    # For fixed batch effect: include batch in design
    if (n_levels == 2) {
      design <- model.matrix(as.formula(paste("~", batch_var, "+", condition_var)), data = metadata)
    } else {
      design <- model.matrix(as.formula(paste("~ 0 +", condition_var, "+", batch_var)), data = metadata)
      colnames(design) <- gsub(condition_var, "", colnames(design))
    }
  } else {
    # No batch variable
    if (n_levels == 2) {
      design <- model.matrix(as.formula(paste("~", condition_var)), data = metadata)
    } else {
      design <- model.matrix(as.formula(paste("~ 0 +", condition_var)), data = metadata)
      colnames(design) <- gsub(condition_var, "", colnames(design))
    }
  }
  
  # Create DGEList for filtering and normalization
  y <- DGEList(counts = count_matrix)
  y <- calcNormFactors(y)
  
  # Voom transformation
  v <- voom(y, design, plot = FALSE)
  
  # Fit linear model with duplicateCorrelation if paired
  if (use_duplicate_correlation) {
    cat("Estimating within-subject correlation...\n")
    dupcor <- duplicateCorrelation(v, design, block = block)
    cat("Consensus correlation:", round(dupcor$consensus.correlation, 3), "\n")
    
    # Re-run voom with correlation
    v <- voom(y, design, plot = FALSE, block = block, correlation = dupcor$consensus.correlation)
    
    # Fit with correlation
    fit <- lmFit(v, design, block = block, correlation = dupcor$consensus.correlation)
  } else {
    # Standard fit
    fit <- lmFit(v, design)
  }
  
  if (n_levels == 2) {
    # Simple comparison
    fit <- eBayes(fit)
    coef_to_use <- if (!use_duplicate_correlation && !is.null(batch_var)) ncol(design) else 2
    res <- topTable(fit, coef = coef_to_use, n = Inf, sort.by = "none")
    res$Gene <- rownames(res)
    res$contrast <- paste0(condition_levels[2], "_vs_", condition_levels[1])
    res_df <- standardize_deg_columns(res, "limma_voom", gene_annotations)
    comparison <- paste(condition_levels, collapse = " vs ")
    
    method_name <- if (use_duplicate_correlation) {
      "limma-voom (paired design with duplicateCorrelation)"
    } else {
      "limma-voom"
    }
  } else {
    # Multiple contrasts
    contrasts_list <- list()
    for(i in 1:(n_levels-1)) {
      for(j in (i+1):n_levels) {
        contrast_name <- paste0(condition_levels[j], "_vs_", condition_levels[i])
        contrast_string <- paste0(condition_levels[j], "-", condition_levels[i])
        contrasts_list[[contrast_name]] <- contrast_string
      }
    }
    
    contrast_matrix <- makeContrasts(contrasts = contrasts_list, levels = design)
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    all_results <- list()
    for(i in 1:ncol(contrast_matrix)) {
      contrast_name <- colnames(contrast_matrix)[i]
      res <- topTable(fit2, coef = i, n = Inf, sort.by = "none")
      res$contrast <- contrast_name
      res$Gene <- rownames(res)
      res <- standardize_deg_columns(res, "limma_voom", gene_annotations)
      all_results[[contrast_name]] <- res
    }
    
    res_df <- do.call(rbind, all_results)
    rownames(res_df) <- NULL
    comparison <- paste(names(contrasts_list), collapse = ", ")
    method_name <- "limma-voom (multiple contrasts)"
  }
  
  result <- list(
    deg_table = res_df,
    method = method_name,
    n_total = nrow(count_matrix),
    n_significant = sum(res_df$adj.P.value < alpha, na.rm = TRUE),
    voom_object = v,
    design_matrix = design,
    comparison = comparison
  )
  
  if (use_duplicate_correlation) {
    result$correlation <- dupcor$consensus.correlation
    result$paired_design <- TRUE
  }
  
  return(result)
}

#' Main wrapper function to run differential expression analysis
run_differential_expression <- function(method, count_matrix, metadata, condition_var,
                                        batch_var = NULL, alpha = 0.05, lfc_threshold = 0,
                                        shrinkage = TRUE, gene_annotations = NULL,
                                        paired_design = FALSE) {
  tryCatch({
    switch(method,
           "deseq2" = run_deseq2_analysis(count_matrix, metadata, condition_var, batch_var,
                                          alpha, lfc_threshold, shrinkage, gene_annotations,
                                          paired_design),
           "edger_GLM" = run_edger_glm_analysis(count_matrix, metadata, condition_var, batch_var,
                                                alpha, lfc_threshold, gene_annotations,
                                                paired_design),
           "limma_voom" = run_limma_voom_analysis(count_matrix, metadata, condition_var, batch_var,
                                                  alpha, lfc_threshold, gene_annotations,
                                                  paired_design),
           stop("Unknown method: ", method)
    )
  }, error = function(e) {
    return(list(
      error = TRUE,
      message = paste("Analysis failed:", e$message),
      method = method
    ))
  })
}
