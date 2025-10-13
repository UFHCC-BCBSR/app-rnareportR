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
  library(NOISeq)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(glue)
  
  source("R/HRK_funcs.R")
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
  
  tmp_out_dir <- tempfile("rmd_tmpdir_")
  dir.create(tmp_out_dir)
  
  # Save it in the list of parameters
  report_params$out_dir <- tmp_out_dir
  
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
  
  # Step 1: Get current gene identifiers from DGE (could be SYMBOL or ENSEMBL)
  gene_ids <- trimws(rownames(DGE))
  
  # Step 2: Detect type — assume SYMBOL if fewer than half look like Ensembl IDs
  n_ensembl <- sum(grepl("^ENS", gene_ids))
  n_total <- length(gene_ids)
  is_symbol <- (n_ensembl / n_total) < 0.5
  
  # Step 3: Retrieve gene annotations from org.Hs.eg.db or similar
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
    DGE <- DGE[valid_ensembl, , keep.lib.sizes = FALSE]  # Remove rows without valid ENSEMBL
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
  
  # Step 4: Attach the annotation table
  DGE$genes <- matched_genes
  
  #############################################
  #########  FILTER GENES  ####################
  #############################################
  
  # Define treatments
  treatment.all <- as.factor(DGE$samples[[report_params$group_var]])
  
  # Filter genes based on expression level (edgeR)
  keep.exprs <- filterByExpr(DGE, group = treatment.all)
  DGE.edgeRfilt <- DGE[keep.exprs, , keep.lib.sizes = FALSE]
  
  # Normalize using TMM method (edgeR)
  DGE.edgeRfilt <- calcNormFactors(DGE.edgeRfilt, method = 'TMM')
  
  # Filter genes based on NOISeq
  DGE.NOIseqfilt <- DGE
  DGE.NOIseqfilt$counts <- filtered.data(
    DGE$counts,
    factor = DGE$samples[[report_params$group_var]],
    norm = FALSE,
    depth = NULL,
    method = 1,
    cv.cutoff = 100,
    cpm = 1,
    p.adj = "fdr"
  )
  
  DGE.NOIseqfilt$samples$lib.size <- apply(DGE.NOIseqfilt$counts, 2, sum)
  DGE.NOIseqfilt <- calcNormFactors(DGE.NOIseqfilt, method = 'TMM')
  
  #############################################
  ########## BEGIN DE TEST ####################
  #############################################
  
  # Create design matrix
  design.mat <- model.matrix(~ 0 + DGE$samples[[report_params$group_var]])
  rownames(design.mat) <- DGE$samples$SampleName
  colnames(design.mat) <- make.names(levels(as.factor(DGE$samples[[report_params$group_var]])))
  
  # Perform voom transformation
  # Create voom plot and save voom object
  png(filename = paste0(report_params$out_dir, "/voom_plot.png"))
  v <- voom(DGE.NOIseqfilt, design.mat, plot = TRUE)
  dev.off()
  
  # ----- Optional batch correction -----
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
    
    # Add corrected matrix to object
    v$E_corrected <- corrected_matrix
    DGE.NOIseqfilt$E_corrected <- corrected_matrix
    
    message("Batch correction complete.")
  } else {
    message("No batch correction applied.")
    v$E_corrected <- v$E
    DGE.NOIseqfilt$E_corrected <- v$E
  }
  
  # Fit linear model
  if (!is.null(report_params$batch_var) && report_params$batch_var %in% colnames(sample.keys)) {
    vfit <- lmFit(v, design.mat, block = batch, correlation = corfit$consensus.correlation)
  } else {
    vfit <- lmFit(v, design.mat)
  }
  
  # Extract contrasts
  if (!is.null(report_params$contrasts) && file.exists(report_params$contrasts)) {
    contrast_strings <- readLines(report_params$contrasts)
    contrast_strings <- contrast_strings[!grepl("^\\s*#", contrast_strings)]
    contrast_strings <- contrast_strings[contrast_strings != ""]
    contrast_strings <- str_replace_all(contrast_strings, "_vs_", "-")
  } else {
    contrast_strings <- unique(unlist(strsplit(sample.keys$Contrast, ";")))
  }
  
  # Generate contrast list
  contrasts_list <- setNames(
    lapply(contrast_strings, function(contrast) {
      groups <- trimws(unlist(strsplit(contrast, "-")))
      
      # Make group names syntactically valid
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
  
  # Helper function for extracting top DE genes
  top_DE_entrezIDs <- function(df, direction = "up", min_genes = 5) {
    msg <- NULL
    
    # Initial filter
    if (direction == "up") {
      filtered_genes <- df %>% filter(logFC > 0.58 & adj.P.value < 0.05)
    } else if (direction == "down") {
      filtered_genes <- df %>% filter(logFC < -0.58 & adj.P.value < 0.05)
    } else {
      stop("Invalid direction. Choose 'up' or 'down'")
    }
    
    # If too few genes, relax criteria
    if (nrow(filtered_genes) < min_genes) {
      msg <- glue::glue(
        "⚠️ {nrow(filtered_genes)} genes passed conservative DE filter for direction '{direction}'. ",
        "Relaxed filter used (logFC > 0 & adj.P < 0.1)."
      )
      
      if (direction == "up") {
        filtered_genes <- df %>% filter(logFC > 0 & adj.P.value < 0.1)
      } else {
        filtered_genes <- df %>% filter(logFC < 0 & adj.P.value < 0.1)
      }
    }
    
    # Still empty?
    if (nrow(filtered_genes) == 0) {
      msg <- glue::glue(
        "❌ No DE genes found for direction '{direction}', even after relaxing thresholds."
      )
      return(list(entrez_ids = character(0), message = msg))
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
        "❌ No DE genes found mapped to Entrez IDs for direction '{direction}': {conditionMessage(e)}"
      )
      return(character(0))
    })
    
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
    keys = rownames(DGE.NOIseqfilt),
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  # Return all necessary objects in a list
  return(list(
    dge_list_raw = DGE,
    dge_list_edgeR = DGE.edgeRfilt,
    dge_list_NOISeq = DGE.NOIseqfilt,
    efit_list = efit_list,
    efit_results_dfs = efit_results_list,
    entrez_ids_list = top_DE_entrezIDs_list,
    universe_entrez = universe_entrez,
    sample.keys = sample.keys,
    lcpm_matrix = cpm(DGE.NOIseqfilt, log = TRUE)
  ))
}