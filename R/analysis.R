# R/analysis.R
run_rnaseq_analysis <- function(report_params) {
  # Load required libraries
  library(limma)
  library(edgeR)
  library(tidyverse)
  library(openxlsx)
  library(pheatmap)
  library(data.table)
  library(org.Hs.eg.db)
  library(Homo.sapiens)
  library(RColorBrewer)
  library(ggrepel)
  library(hues)
  library(Biobase)
  library(stringr)
  library(ggfortify)
  library(qs)
  library(DT)
  library(dplyr)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(glue)
  library(DESeq2)  # NEW: Added for DESeq2 support
  
  source("R/helper-funcs.R")
  source("R/filter-genes.R")      # NEW: Source filtering functions
  source("R/run-DE.R")             # NEW: Source DE analysis functions
  
  select <- dplyr::select
  
  if (!requireNamespace("statmod", quietly = TRUE)) {
    install.packages("statmod", repos = "https://cloud.r-project.org")
  }
  
  #############################################
  #########  SETUP DGE OBJECT   ###############
  #############################################
  
  annotation_obj <- get(report_params$annotation_db, envir = asNamespace(report_params$annotation_db))
  
  # Treat empty or missing batch_var as NULL
  if (!("batch_var" %in% names(report_params)) || is.null(report_params$batch_var) ||
      report_params$batch_var == "" || report_params$batch_var == "None") {
    report_params$batch_var <- NULL
  }
  
  # Set default filtering method if not provided
  if (!("filter_method" %in% names(report_params))) report_params$filter_method <- "edgeR"
  
  # Set default edgeR filtering parameters if not provided
  if (!("filter_min_count" %in% names(report_params))) report_params$filter_min_count <- 10
  if (!("filter_min_prop" %in% names(report_params))) report_params$filter_min_prop <- 0.7
  
  # Set default NOISeq filtering parameters if not provided
  if (!("noiseq_method" %in% names(report_params))) report_params$noiseq_method <- 1
  if (!("cv_cutoff" %in% names(report_params))) report_params$cv_cutoff <- 100
  if (!("cpm" %in% names(report_params))) report_params$cpm <- 1
  if (!("p_adj" %in% names(report_params))) report_params$p_adj <- "fdr"
  
  # Ensure numeric conversion (in case they came from command line as strings)
  report_params$filter_min_count <- as.numeric(report_params$filter_min_count)
  report_params$filter_min_prop <- as.numeric(report_params$filter_min_prop)
  report_params$noiseq_method <- as.numeric(report_params$noiseq_method)
  report_params$cv_cutoff <- as.numeric(report_params$cv_cutoff)
  report_params$cpm <- as.numeric(report_params$cpm)
  # p_adj stays as character (it's "fdr", "bonferroni", etc.)

  # Set default DE tool if not provided
  if (!("DE_tool" %in% names(report_params))) report_params$DE_tool <- "limma_voom"
  
  # Use the permanent output directory (not temp)
  output_dir <- report_params$output_path
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Step 1: List RSEM files
  files.list <- list.files(path = report_params$rsem_dir, pattern = 'genes.results$', full.names = TRUE)
  
  # Step 2: Extract sample names from filenames (remove suffix)
  file_sample_names <- gsub(".genes.results$", "", basename(files.list))
  
  # Step 3: Create data frame mapping sample names to files
  file_df_raw <- data.frame(SampleName = file_sample_names, File = files.list, stringsAsFactors = FALSE)
  
  # Step 4: Load sample keys
  sample.keys <- read.csv(report_params$sample_data)
  colnames(sample.keys)[1] <- "SampleName"
  
  # Step 5: Flexible matching - find one matching file per sample
  matched_files <- character(nrow(sample.keys))
  matched_names <- character(nrow(sample.keys))
  
  for (i in seq_along(sample.keys$SampleName)) {
    sample_name <- sample.keys$SampleName[i]
    matches <- grep(sample_name, file_sample_names, value = TRUE)
    
    if (length(matches) == 0) {
      stop(paste0("No file matched sample name: ", sample_name))
    } else if (length(matches) > 1) {
      stop(paste0("Multiple files matched sample name: ", sample_name, "\nMatches: ", paste(matches, collapse = ", ")))
    } else {
      matched_names[i] <- matches
      matched_files[i] <- files.list[which(file_sample_names == matches)]
    }
  }
  
  file_df <- data.frame(SampleName = sample.keys$SampleName, File = matched_files, stringsAsFactors = FALSE)
  
  # Step 6: Read only matched files
  DGE <- edgeR::readDGE(files = file_df$File, columns = c(1, 5))
  
  # Step 7: Rename samples
  colnames(DGE$counts) <- file_df$SampleName
  rownames(DGE$samples) <- file_df$SampleName
  DGE$samples$SampleName <- file_df$SampleName
  
  # Step 8: Merge metadata (safe way that preserves order)
  DGE$samples <- dplyr::left_join(DGE$samples, sample.keys, by = "SampleName")
  rownames(DGE$samples) <- DGE$samples$SampleName
  
  # Step 9: Get current gene identifiers from DGE (could be SYMBOL or ENSEMBL)
  gene_ids <- trimws(rownames(DGE))
  
  # Step 10: Detect type — assume SYMBOL if fewer than half look like Ensembl IDs
  n_ensembl <- sum(grepl("^ENS", gene_ids))
  n_total <- length(gene_ids)
  is_symbol <- (n_ensembl / n_total) < 0.5
  
  # Step 11: Retrieve gene annotations from org.Hs.eg.db or similar
  if (is_symbol) {
    # Input rownames are SYMBOLs, we want to convert to ENSEMBL IDs
    genes <- AnnotationDbi::select(annotation_obj,
                                   keys = gene_ids,
                                   columns = c("SYMBOL", "ENSEMBL"),
                                   keytype = "SYMBOL")
    genes <- genes[!duplicated(genes$SYMBOL), ]
    matched_genes <- genes[match(gene_ids, genes$SYMBOL), ]
    
    # Update rownames to Ensembl IDs
    valid_ensembl <- !is.na(matched_genes$ENSEMBL)
    rownames(DGE) <- matched_genes$ENSEMBL
    DGE <- DGE[valid_ensembl, , keep.lib.sizes = FALSE]
    matched_genes <- matched_genes[valid_ensembl, ]
  } else {
    # Input rownames are already ENSEMBL
    genes <- AnnotationDbi::select(annotation_obj,
                                   keys = gene_ids,
                                   columns = c("ENSEMBL", "SYMBOL"),
                                   keytype = "ENSEMBL")
    genes <- genes[!duplicated(genes$ENSEMBL), ]
    matched_genes <- genes[match(gene_ids, genes$ENSEMBL), ]
  }
  
  # Step 12: Attach the annotation table
  DGE$genes <- matched_genes
  
  #############################################
  #########  FILTER GENES  ####################
  #############################################
  
  # Choose method: "edgeR" or "NOISeq"
  filter_method <- report_params$filter_method  # Set this in your params
  
  message(sprintf("Filtering genes with %s...", filter_method))
  
  # Define treatment groups
  treatment.all <- as.factor(DGE$samples[[report_params$group_var]])
  
  # Call filtering
  if (filter_method == "edgeR") {
    filter_result <- filter_genes(
      count_matrix = DGE$counts,
      metadata = sample.keys,
      condition_var = report_params$group_var,
      method = "edgeR",
      min_count = report_params$filter_min_count,
      min_prop = report_params$filter_min_prop
    )
  } else if (filter_method == "NOISeq") {
    filter_result <- filter_genes(
      count_matrix = DGE$counts,
      metadata = sample.keys,
      condition_var = report_params$group_var,
      method = "NOISeq",
      noiseq_method = report_params$noiseq_method,  # typically 1
      cv_cutoff = report_params$cv_cutoff,          # typically 100
      cpm = report_params$cpm,                      # typically 1
      p_adj = report_params$p_adj                   # typically "fdr"
    )
  }
  
  # Create filtered DGE object (same for both methods)
  DGE.filtered <- DGE[filter_result$keep_genes, , keep.lib.sizes = FALSE]
  
  # Normalize using TMM method (same for both methods)
  DGE.filtered <- calcNormFactors(DGE.filtered, method = 'TMM')
  
  message(sprintf("Genes retained: %d of %d (%.1f%%)",
                  filter_result$n_filtered,
                  filter_result$n_original,
                  filter_result$prop_retained * 100))
  
  #############################################
  ########## BEGIN DE ANALYSIS ################
  #############################################
  
  message(paste("Running differential expression analysis with", report_params$DE_tool))
  
  # Extract contrasts
  if (!is.null(report_params$contrasts) && file.exists(report_params$contrasts)) {
    contrast_strings <- readLines(report_params$contrasts)
    contrast_strings <- contrast_strings[!grepl("^\\s*#", contrast_strings)]
    contrast_strings <- contrast_strings[contrast_strings != ""]
    contrast_strings <- str_replace_all(contrast_strings, "_vs_", "-")
  } else {
    contrast_strings <- unique(unlist(strsplit(sample.keys$Contrast, ";")))
  }
  
  # Prepare count matrix and metadata for DE analysis
  count_matrix_for_de <- DGE.filtered$counts
  metadata_for_de <- DGE.filtered$samples
  
  # Create gene annotations data frame for integration
  gene_annotations <- data.frame(
    gene_id = rownames(DGE.filtered),
    gene_symbol = DGE.filtered$genes$SYMBOL,
    stringsAsFactors = FALSE
  )
  rownames(gene_annotations) <- gene_annotations$gene_id
  
  # NEW: Run DE analysis based on selected method
  if (report_params$DE_tool == "limma_voom") {
    
    # ===== LIMMA-VOOM ANALYSIS =====
    
    # Create design matrix
    design.mat <- model.matrix(~ 0 + DGE.filtered$samples[[report_params$group_var]])
    rownames(design.mat) <- DGE.filtered$samples$SampleName
    colnames(design.mat) <- make.names(levels(as.factor(DGE.filtered$samples[[report_params$group_var]])))
    
    # Create voom plot and save voom object
    png(filename = file.path(output_dir, "voom_plot.png"))
    v <- voom(DGE.filtered, design.mat, plot = TRUE)
    dev.off()
    
    # Optional batch correction
    if (!is.null(report_params$batch_var) && report_params$batch_var %in% colnames(sample.keys)) {
      message("Batch correction: Running duplicateCorrelation and removeBatchEffect...")
      batch <- sample.keys[[report_params$batch_var]]
      design <- model.matrix(~ 0 + sample.keys[[report_params$group_var]])
      
      corfit <- duplicateCorrelation(v, design, block = batch)
      corrected_matrix <- removeBatchEffect(
        v$E,
        batch = batch,
        design = design,
        correlation = corfit$consensus.correlation,
        block = batch
      )
      
      v$E_corrected <- corrected_matrix
      DGE.filtered$E_corrected <- corrected_matrix
      message("Batch correction complete.")
    } else {
      message("No batch correction applied.")
      v$E_corrected <- v$E
      DGE.filtered$E_corrected <- v$E
    }
    
    # Fit linear model
    if (!is.null(report_params$batch_var) && report_params$batch_var %in% colnames(sample.keys)) {
      vfit <- lmFit(v, design.mat, block = batch, correlation = corfit$consensus.correlation)
    } else {
      vfit <- lmFit(v, design.mat)
    }
    
    # Generate contrast list
    contrasts_list <- setNames(
      lapply(contrast_strings, function(contrast) {
        groups <- trimws(unlist(strsplit(contrast, "-")))
        group1 <- make.names(groups[1])
        group2 <- make.names(groups[2])
        print(paste("Processing contrast:", group1, "-", group2))
        makeContrasts(
          contrasts = paste0("`", group1, "` - `", group2, "`"),
          levels = colnames(vfit$design)
        )
      }),
      paste0(gsub("-", "_vs_", contrast_strings))
    )
    
    # Compute contrasts and apply eBayes
    efit_list <- lapply(contrasts_list, function(contr) {
      vfit_contr <- contrasts.fit(vfit, contr)
      eBayes(vfit_contr)
    })
    
    # Create result data frames dynamically
    efit_results_list <- setNames(
      lapply(seq_along(efit_list), function(i) {
        create_efit_results_df(efit_list[[i]])
      }),
      paste0("efit_", names(efit_list), "_results_df")
    )
    
  } else {
    
    # ===== DESEQ2 OR EDGER GLM ANALYSIS =====
    
    # Run the appropriate DE method
    de_results <- run_differential_expression(
      method = report_params$DE_tool,
      count_matrix = count_matrix_for_de,
      metadata = metadata_for_de,
      condition_var = report_params$group_var,
      batch_var = report_params$batch_var,
      alpha = 0.05,
      lfc_threshold = 0,
      shrinkage = TRUE,
      gene_annotations = gene_annotations
    )
    
    # Check for errors
    if (!is.null(de_results$error) && de_results$error) {
      stop(paste("DE analysis failed:", de_results$message))
    }
    
    # Convert to efit_results_list format for compatibility with downstream code
    # Split by contrast if multiple contrasts exist
    if ("contrast" %in% colnames(de_results$deg_table)) {
      unique_contrasts <- unique(de_results$deg_table$contrast)
      
      efit_results_list <- setNames(
        lapply(unique_contrasts, function(contrast_name) {
          contrast_df <- de_results$deg_table[de_results$deg_table$contrast == contrast_name, ]
          # Rename columns to match limma-voom format
          contrast_df <- contrast_df %>%
            rename(
              SYMBOL = if("SYMBOL" %in% colnames(.)) SYMBOL else gene_symbol,
              ensembleID = if("ensembleID" %in% colnames(.)) ensembleID else Gene
            )
          # Add contrast column
          contrast_df$contrast <- contrast_name
          contrast_df
        }),
        paste0("efit_", gsub("_vs_", "_", unique_contrasts), "_results_df")
      )
    } else {
      # Single contrast
      contrast_name <- gsub(" vs ", "_vs_", de_results$comparison)
      efit_results_list <- list()
      efit_results_list[[paste0("efit_", contrast_name, "_results_df")]] <- de_results$deg_table %>%
        rename(
          SYMBOL = if("SYMBOL" %in% colnames(.)) SYMBOL else gene_symbol,
          ensembleID = if("ensembleID" %in% colnames(.)) ensembleID else Gene
        ) %>%
        mutate(contrast = contrast_name)
    }
    
    # For compatibility, create a placeholder for voom objects
    # (won't be used in plotting but needed for function signatures)
    v <- NULL
    efit_list <- NULL
    
    # Generate logCPM matrix for visualization
    DGE.filtered$E_corrected <- cpm(DGE.filtered, log = TRUE)
  }
  
  #############################################
  ######### ENRICHMENT ANALYSIS ###############
  #############################################
  
  # Helper function for extracting top DE genes
  # In analysis.R, update the top_DE_entrezIDs function:
  top_DE_entrezIDs <- function(df, direction = "up", min_genes = 5) {
    msg <- NULL
    used_relaxed_filter <- FALSE  # Track if we used relaxed criteria
    
    # Determine logFC column name
    logfc_col <- if("logFC" %in% colnames(df)) "logFC" else if("log2FoldChange" %in% colnames(df)) "log2FoldChange" else "logFC"
    
    # Initial filter
    if (direction == "up") {
      filtered_genes <- df %>% filter(!!sym(logfc_col) > 0.58 & adj.P.value < 0.05)
    } else if (direction == "down") {
      filtered_genes <- df %>% filter(!!sym(logfc_col) < -0.58 & adj.P.value < 0.05)
    } else {
      stop("Invalid direction. Choose 'up' or 'down'")
    }
    
    # If too few genes, relax criteria
    if (nrow(filtered_genes) < min_genes) {
      pval_col <- if("P.value" %in% colnames(df)) "P.value" else if("pvalue" %in% colnames(df)) "pvalue" else "P.Value"
      
      if (direction == "up") {
        filtered_genes <- df %>% filter(!!sym(logfc_col) > 0 & !!sym(pval_col) < 0.05)
      } else {
        filtered_genes <- df %>% filter(!!sym(logfc_col) < 0 & !!sym(pval_col) < 0.05)
      }
      
      # If relaxed filter found genes, keep a warning
      if (nrow(filtered_genes) >= min_genes) {
        used_relaxed_filter <- TRUE
        msg <- glue::glue(
          "⚠️ Using relaxed thresholds (raw P < 0.05, any logFC) for '{direction}' genes. ",
          "Found {nrow(filtered_genes)} genes."
        )
      } else {
        # Still no genes even with relaxed filter
        msg <- glue::glue(
          "❌ No DE genes found for direction '{direction}', even after relaxing thresholds."
        )
        return(list(entrez_ids = character(0), message = msg))
      }
    }
    
    # Try mapping to Entrez IDs
    entrez_ids <- tryCatch({
      mapIds(
        annotation_obj,
        keys = filtered_genes$ensembleID,
        column = "ENTREZID",
        keytype = "ENSEMBL",
        multiVals = "first"
      )
    }, error = function(e) {
      msg <<- glue::glue(
        "❌ Error mapping to Entrez IDs for direction '{direction}': {conditionMessage(e)}"
      )
      return(character(0))
    })
    
    # If mapping failed
    if (length(na.omit(entrez_ids)) == 0) {
      return(list(entrez_ids = character(0), message = msg))
    }
    
    # Return genes with message (NULL if conservative filter worked, warning if relaxed)
    return(list(entrez_ids = na.omit(entrez_ids), message = msg))
  }
  
  # Apply function dynamically for each contrast
  top_DE_entrezIDs_list <- setNames(
    lapply(efit_results_list, function(df) {
      list(
        up = top_DE_entrezIDs(df, "up"),
        down = top_DE_entrezIDs(df, "down")
      )
    }),
    paste0("top_DE_entrezIDs_", names(efit_results_list))
  )
  
  # Universe for enrichment
  universe_entrez <- mapIds(
    annotation_obj,
    keys = rownames(DGE.filtered),
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  #############################################
  ######### RETURN RESULTS ####################
  #############################################
  
  # Return all necessary objects in a list
  return(list(
    # Raw and filtered data (names match Rmd expectations)
    dge_list_raw = DGE,
    #dge_list_NOISeq = DGE.filtered,  # Primary filtered object used downstream
    dge_list_filt = DGE.filtered,   # Duplicate for backward compatibility
    
    # DE results
    efit_list = if(exists("efit_list")) efit_list else NULL,  # NULL for non-limma
    efit_results_dfs = efit_results_list,
    
    # Enrichment data
    entrez_ids_list = top_DE_entrezIDs_list,
    universe_entrez = universe_entrez,
    
    # Metadata and visualization
    sample.keys = sample.keys,
    lcpm_matrix = cpm(DGE.filtered, log = TRUE),
    
    # Method-specific outputs
    voom_plot_path = if(report_params$DE_tool == "limma_voom") {
      file.path(output_dir, "voom_plot.png")
    } else {
      NULL
    },
    de_method = report_params$DE_tool,
    filter_stats = filter_result
  ))
}