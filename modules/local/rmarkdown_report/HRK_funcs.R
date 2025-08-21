library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(plotly)
library(heatmaply)
library(org.Mm.eg.db)

generate_enrichment_plot <- function(gene_lists, de_results_df, universe_entrez, ont_category, significance_threshold = 0.05, top_n = 10, annotation_db) {
  pvalue_cutoff <- if (isTRUE(report_params$download_nonsig_enrich)) 1 else significance_threshold
  annotation_obj <- get(annotation_db, envir = asNamespace(annotation_db)) 
  names(gene_lists) <- gsub("efit_|_results_df", "", names(gene_lists))
  # Ensure gene lists are named and define contrast order
  if (is.null(names(gene_lists))) stop("Each gene list must be named!")
  contrast_order <- names(gene_lists)
  
  # Convert DE results dataframe to a lookup table for mapping Entrez IDs
  de_results_df <- de_results_df %>%
    mutate(ENTREZID = mapIds(annotation_obj, keys = ensembleID, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first"))
  
  # Run GO enrichment using compareCluster()
  formula_res <- compareCluster(
    Entrez ~ Contrast,
    data = bind_rows(lapply(contrast_order, function(contrast) {
      data.frame(
        Entrez = gene_lists[[contrast]],
        Contrast = contrast
      )
    })),
    fun = "enrichGO",
    universe = na.omit(universe_entrez),
    OrgDb = annotation_obj,
    keyType = "ENTREZID",
    ont = ont_category,
    pvalueCutoff = pvalue_cutoff
  )
  
  # Handle no results case
  if (is.null(formula_res) || nrow(formula_res@compareClusterResult) == 0) {
    formula_res <- new("compareClusterResult",
                       compareClusterResult = data.frame(
                         Cluster = factor(),
                         ID = character(),
                         Description = character(),
                         GeneRatio = character(),
                         BgRatio = character(),
                         pvalue = numeric(),
                         p.adjust = numeric(),
                         qvalue = numeric(),
                         geneID = character(),
                         Count = integer(),
                         stringsAsFactors = FALSE
                       ))
  }
  
  # Ensure clusters are ordered correctly
  formula_res@compareClusterResult$Cluster <- factor(
    formula_res@compareClusterResult$Cluster,
    levels = contrast_order
  )
  
  # **Filter only by significance threshold FIRST (to preserve all significant results)**
  filtered_results <- subset(
    formula_res@compareClusterResult,
    p.adjust <= significance_threshold
  )
  
  # Handle the special case: no enrichment found
  if (nrow(filtered_results) == 0) {
    message_plot <- ggplot() +
      annotate("text", x = 1, y = 1, label = paste0("No significant GO enrichment found\n(", ont_category, ")"), size = 6, hjust = 0.5) +
      theme_void() +
      ggtitle(paste("GO Term Enrichment (", ont_category, ")", sep = ""))
    
    interactive_plot <- ggplotly(message_plot)
    static_plot <- message_plot
    
    return(list(
      interactive_plot = interactive_plot,
      static_plot = static_plot,
      go_results = NULL
    ))
  }
  
  filtered_results$GeneSymbols <- sapply(seq_len(nrow(filtered_results)), function(i) {
    gene_list <- filtered_results$geneID[i]
    contrast_full <- as.character(filtered_results$Cluster[i]) 
    
    # Split into base contrast and direction
    contrast_base <- sub("\\.(up|down)$", "", contrast_full)
    direction <- sub("^.*\\.", "", contrast_full)
    
    entrez_ids <- unlist(strsplit(gene_list, "/"))
    
    # Match to relevant DE results
    de_sub <- de_results_df %>%
      filter(grepl(contrast_base, contrast),  # Match the base name
             ENTREZID %in% entrez_ids,
             case_when(
               direction == "up" ~ logFC > 0,
               direction == "down" ~ logFC < 0
             ))
    
    # Sort by p-value and get top 20 Entrez IDs
    top_genes <- de_sub %>%
      arrange(adj.P.value) %>%
      slice_head(n = 20) %>%
      pull(ENTREZID)
    
    gene_symbols <- mapIds(annotation_obj,
                           keys = top_genes,
                           column = "SYMBOL",
                           keytype = "ENTREZID",
                           multiVals = "first") %>%
      na.omit()
    
    if (length(gene_symbols) > 0) {
      paste(gene_symbols, collapse = "<br>")
    } else {
      NA_character_
    }
  })
  
  # **Save full filtered results for downloading**
  # Store all raw results before filtering (optional for download)
  all_results <- formula_res@compareClusterResult
  
  # Recalculate GeneSymbols for ALL terms if download_nonsig_enrich is TRUE
  if (isTRUE(report_params$download_nonsig_enrich)) {
    all_results$GeneSymbols <- sapply(seq_len(nrow(all_results)), function(i) {
      gene_list <- all_results$geneID[i]
      contrast_full <- as.character(all_results$Cluster[i]) 
      
      contrast_base <- sub("\\.(up|down)$", "", contrast_full)
      direction <- sub("^.*\\.", "", contrast_full)
      
      entrez_ids <- unlist(strsplit(gene_list, "/"))
      
      de_sub <- de_results_df %>%
        filter(grepl(contrast_base, contrast),
               ENTREZID %in% entrez_ids,
               case_when(
                 direction == "up" ~ logFC > 0,
                 direction == "down" ~ logFC < 0
               ))
      
      top_genes <- de_sub %>%
        arrange(adj.P.value) %>%
        slice_head(n = 20) %>%
        pull(ENTREZID)
      
      gene_symbols <- mapIds(annotation_obj,
                             keys = top_genes,
                             column = "SYMBOL",
                             keytype = "ENTREZID",
                             multiVals = "first") %>%
        na.omit()
      
      if (length(gene_symbols) > 0) {
        paste(gene_symbols, collapse = "<br>")
      } else {
        NA_character_
      }
    })
  }
  
  # Choose download results based on flag
  download_go_results <- if (isTRUE(report_params$download_nonsig_enrich)) {
    all_results %>%
      select(Cluster, Description, p.adjust, GeneSymbols, everything())
  } else {
    filtered_results %>%
      select(Cluster, Description, p.adjust, GeneSymbols, everything())
  }
  
  
  # **Identify top `n` GO terms across all clusters for plotting**
  top_GO_terms <- filtered_results %>%
    group_by(Cluster) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()
  
  # **Filter for plotting (only keeping top GO terms)**
  formula_res@compareClusterResult <- filtered_results %>%
    filter(Description %in% top_GO_terms)
  
  # Convert GeneRatio to numeric
  formula_res@compareClusterResult <- formula_res@compareClusterResult %>%
    mutate(GeneRatio = sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
  
  # **Restore GO Term Ordering Using Hierarchical Clustering**
  reorder_GO_terms <- function(df) {
    term_matrix <- table(df$Description, df$Cluster)  
    
    if (nrow(term_matrix) > 1 && length(unique(df$Description)) > 1) {
      term_dist <- dist(term_matrix, method = "binary")  
      term_hclust <- hclust(term_dist, method = "ward.D2")  
      term_order <- rownames(term_matrix)[term_hclust$order]
    } else {
      term_order <- unique(df$Description)  # Ensure a valid factor level
    }
    
    # Safely convert to a factor with the correct order
    df$Description <- factor(df$Description, levels = term_order)  
    return(df)
  }
  
  formula_res@compareClusterResult <- reorder_GO_terms(formula_res@compareClusterResult)
  
  # **Define GeneRatio bins**
  bin_breaks <- c(0, 0.01, 0.05, 0.10, max(formula_res@compareClusterResult$GeneRatio, na.rm = TRUE) + 0.01)  
  bin_labels <- c("≤0.01", "0.01 - 0.05", "0.05 - 0.10", "≥0.10")
  
  formula_res@compareClusterResult <- formula_res@compareClusterResult %>%
    mutate(GeneRatioCategory = cut(GeneRatio, breaks = bin_breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE))
  
  # **Define reference dot sizes**
  size_mapping <- c("≤0.01" = 2, "0.01 - 0.05" = 4, "0.05 - 0.10" = 6, "≥0.10" = 8)
  
  # **Ensure all necessary variables are characters/numeric**
  formula_res@compareClusterResult <- formula_res@compareClusterResult %>%
    mutate(
      p.adjust = as.numeric(as.character(p.adjust)),
      GeneRatio = as.numeric(as.character(GeneRatio)),
      plot_label = ifelse(
        sapply(strsplit(as.character(Description), " "), length) > 6,
        sapply(strsplit(as.character(Description), " "), function(words) {
          paste(c(words[1:3], "...", tail(words, 3)), collapse = " ")
        }),
        as.character(Description)
      ),
      tooltip_text = paste(
        "Cluster: ", as.character(Cluster), "<br>",  
        "GO Term: ", as.character(Description), "<br>",  
        "p.adjust: ", signif(p.adjust, 3), "<br>",
        "GeneRatio: ", signif(GeneRatio, 3), "<br>",
        "Top Genes:<br>", as.character(GeneSymbols)  
      )
    )
  
  # **Create ggplot object**
  p <- ggplot(formula_res@compareClusterResult, aes(
    x = Cluster, 
    y = plot_label, 
    size = GeneRatioCategory,  
    color = p.adjust,
    text = tooltip_text  
  )) +
    geom_point(alpha = 0.8) +  
    scale_size_manual(name = "Gene Ratio", values = size_mapping) +  
    scale_color_gradient(low = "red", high = "blue", 
                         limits = c(min(formula_res@compareClusterResult$p.adjust, na.rm = TRUE), 
                                    max(formula_res@compareClusterResult$p.adjust, na.rm = TRUE)),
                         name = "p.adjust") +  
    guides(color = guide_colorbar(title = "p.adjust")) +  
    ggtitle(paste("GO Term Enrichment (", ont_category, ")", sep = "")) +
    xlab("DE Gene list") +
    ylab("GO Term") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # **Convert to interactive plot (ensure tooltips work)**
  interactive_plot <- ggplotly(p, tooltip = "text") %>%  
    layout(legend = list(title = list(text = "Gene Ratio")))
  
  # **Static High-Resolution Plot (for manuscript)**
  static_plot <- dotplot(formula_res, showCategory = top_n) +
    ggtitle(paste("GO Term Enrichment (", ont_category, ")", sep = "")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # **Return interactive plot, static plot, and full GO results**
  return(list(
    interactive_plot = interactive_plot, 
    static_plot = static_plot, 
    go_results = download_go_results
  ))
}

# Run the function
#GO_BP_results <- generate_enrichment_plot(
#  gene_lists = gene_lists, 
#  de_results_df = bind_rows(efit.1.results_df, efit.2.results_df, efit.3.results_df), 
#  universe_entrez = universe_entrez,
#  ont_category = "BP"
#)
# Run the function
#GO_MF_results <- generate_enrichment_plot(
#  gene_lists = gene_lists, 
#  de_results_df = bind_rows(efit.1.results_df, efit.2.results_df, efit.3.results_df), 
#  universe_entrez = universe_entrez,
#  ont_category = "MF"
#)

# Save high-resolution PNG for manuscript
#ggsave("GO_Enrichment_BP.png", plots$static_plot, width = 10, height = 12, dpi = 300,bg = "white")

generate_kegg_enrichment_plot <- function(gene_lists, de_results_df, universe_entrez, significance_threshold = 0.05, top_n = 10,annotation_db) {
  pvalue_cutoff <- if (isTRUE(report_params$download_nonsig_enrich)) 1 else significance_threshold
  
  annotation_obj <- get(annotation_db, envir = asNamespace(annotation_db)) 
  names(gene_lists) <- gsub("efit_|_results_df","",names(gene_lists))
  
  # Ensure gene lists are named and define contrast order
  if (is.null(names(gene_lists))) stop("Each gene list must be named!")
  contrast_order <- names(gene_lists)
  
  # Convert DE results dataframe to a lookup table for mapping Entrez IDs
  de_results_df <- de_results_df %>%
    mutate(ENTREZID = mapIds(annotation_obj, keys = ensembleID, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first"))
  
  # Run KEGG enrichment using compareCluster()
  kegg_res <- compareCluster(
    Entrez ~ Contrast,
    data = bind_rows(lapply(contrast_order, function(contrast) {
      data.frame(
        Entrez = gene_lists[[contrast]],
        Contrast = contrast
      )
    })),
    fun = "enrichKEGG",
    universe = na.omit(universe_entrez),
    organism = report_params$organism,
    keyType = "kegg",
    pvalueCutoff = pvalue_cutoff
  )
  
  # Ensure clusters are ordered correctly
  kegg_res@compareClusterResult$Cluster <- factor(
    kegg_res@compareClusterResult$Cluster,
    levels = contrast_order
  )
  
  # **Filter only by significance threshold FIRST (to preserve all significant results)**
  filtered_results <- subset(
    kegg_res@compareClusterResult,
    p.adjust <= significance_threshold
  )
  
  # Convert Entrez IDs to Gene Symbols and select top 20 by p-value **within the correct cluster & direction**
  filtered_results$GeneSymbols <- sapply(seq_len(nrow(filtered_results)), function(i) {
    gene_list <- filtered_results$geneID[i]
    contrast_full <- as.character(filtered_results$Cluster[i]) 
    
    # Split into base contrast and direction
    contrast_base <- sub("\\.(up|down)$", "", contrast_full)
    direction <- sub("^.*\\.", "", contrast_full)
    
    entrez_ids <- unlist(strsplit(gene_list, "/"))
    
    # Match to relevant DE results
    de_sub <- de_results_df %>%
      filter(grepl(contrast_base, contrast),  # Match the base name
             ENTREZID %in% entrez_ids,
             case_when(
               direction == "up" ~ logFC > 0,
               direction == "down" ~ logFC < 0
             ))
    
    # Sort by p-value and get top 20 Entrez IDs
    top_genes <- de_sub %>%
      arrange(adj.P.value) %>%
      slice_head(n = 20) %>%
      pull(ENTREZID)
    
    gene_symbols <- mapIds(annotation_obj,
                           keys = top_genes,
                           column = "SYMBOL",
                           keytype = "ENTREZID",
                           multiVals = "first") %>%
      na.omit()
    
    if (length(gene_symbols) > 0) {
      paste(gene_symbols, collapse = "<br>")
    } else {
      NA_character_
    }
  })
  
  
  
  # Store all raw results before filtering (optional for download)
  all_results <- kegg_res@compareClusterResult
  
  # Recalculate GeneSymbols for ALL terms if download_nonsig_enrich is TRUE
  if (isTRUE(report_params$download_nonsig_enrich)) {
    all_results$GeneSymbols <- sapply(seq_len(nrow(all_results)), function(i) {
      gene_list <- all_results$geneID[i]
      contrast_full <- as.character(all_results$Cluster[i]) 
      
      contrast_base <- sub("\\.(up|down)$", "", contrast_full)
      direction <- sub("^.*\\.", "", contrast_full)
      
      entrez_ids <- unlist(strsplit(gene_list, "/"))
      
      de_sub <- de_results_df %>%
        filter(grepl(contrast_base, contrast),
               ENTREZID %in% entrez_ids,
               case_when(
                 direction == "up" ~ logFC > 0,
                 direction == "down" ~ logFC < 0
               ))
      
      top_genes <- de_sub %>%
        arrange(adj.P.value) %>%
        slice_head(n = 20) %>%
        pull(ENTREZID)
      
      gene_symbols <- mapIds(annotation_obj,
                             keys = top_genes,
                             column = "SYMBOL",
                             keytype = "ENTREZID",
                             multiVals = "first") %>%
        na.omit()
      
      if (length(gene_symbols) > 0) {
        paste(gene_symbols, collapse = "<br>")
      } else {
        NA_character_
      }
    })
  }
  
  # Choose download results based on flag
  download_kegg_results <- if (isTRUE(report_params$download_nonsig_enrich)) {
    all_results %>%
      select(Cluster, Description, p.adjust, GeneSymbols, everything())
  } else {
    filtered_results %>%
      select(Cluster, Description, p.adjust, GeneSymbols, everything())
  }
  
  
  # **Identify top `n` KEGG pathways across all clusters for plotting**
  top_KEGG_terms <- filtered_results %>%
    group_by(Cluster) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()
  
  # **Filter for plotting (only keeping top KEGG pathways)**
  kegg_res@compareClusterResult <- filtered_results %>%
    filter(Description %in% top_KEGG_terms)
  
  # Convert GeneRatio to numeric
  kegg_res@compareClusterResult <- kegg_res@compareClusterResult %>%
    mutate(GeneRatio = sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
  
  # **Restore KEGG Pathway Ordering Using Hierarchical Clustering**
  reorder_KEGG_terms <- function(df) {
    term_matrix <- table(df$Description, df$Cluster)  
    
    if (nrow(term_matrix) > 1) {
      term_dist <- dist(term_matrix, method = "binary")  
      term_hclust <- hclust(term_dist, method = "ward.D2")  
      term_order <- unique(rownames(term_matrix)[term_hclust$order])
    } else {
      term_order <- unique(df$Description)  # <--- this was the bug
    }
    
    df$Description <- factor(df$Description, levels = term_order)  
    return(df)
  }
  
  
  kegg_res@compareClusterResult <- reorder_KEGG_terms(kegg_res@compareClusterResult)
  
  # **Define GeneRatio bins**
  bin_breaks <- c(0, 0.01, 0.05, 0.10, max(kegg_res@compareClusterResult$GeneRatio, na.rm = TRUE) + 0.01)  
  bin_labels <- c("≤0.01", "0.01 - 0.05", "0.05 - 0.10", "≥0.10")
  
  kegg_res@compareClusterResult <- kegg_res@compareClusterResult %>%
    mutate(GeneRatioCategory = cut(GeneRatio, breaks = bin_breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE))
  
  # **Define reference dot sizes**
  size_mapping <- c("≤0.01" = 2, "0.01 - 0.05" = 4, "0.05 - 0.10" = 6, "≥0.10" = 8)
  
  # **Ensure all necessary variables are characters/numeric**
  kegg_res@compareClusterResult <- kegg_res@compareClusterResult %>%
    mutate(
      p.adjust = as.numeric(as.character(p.adjust)),
      GeneRatio = as.numeric(as.character(GeneRatio)),
      tooltip_text = paste(
        "Cluster: ", as.character(Cluster), "<br>",  
        "KEGG Pathway: ", as.character(Description), "<br>",  
        "p.adjust: ", signif(p.adjust, 3), "<br>",
        "GeneRatio: ", signif(GeneRatio, 3), "<br>",
        "Top Genes:<br>", as.character(GeneSymbols)  
      )
    )
  
  # **Create ggplot object**
  p <- ggplot(kegg_res@compareClusterResult, aes(
    x = Cluster, 
    y = Description, 
    size = GeneRatioCategory, 
    color = p.adjust,
    text = tooltip_text
  )) +
    geom_point(alpha = 0.8) +  
    scale_size_manual(name = "Gene Ratio", values = size_mapping) +  
    scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +  
    guides(color = guide_colorbar(title = "p.adjust")) +  
    ggtitle("KEGG Pathway Enrichment") +
    xlab("DE Gene list") +
    ylab("KEGG Pathway") +
    theme_minimal()+
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) 
  
  interactive_plot <- ggplotly(p, tooltip = "text") %>%  
    layout(legend = list(
      title = list(text = "Gene Ratio"),  # Ensure Gene Ratio is labeled
      colorbar = list(title = "p.adjust") # Ensure p.adjust is labeled for color
    ))
  
  
  
  # **Static KEGG Plot**
  static_plot <- dotplot(kegg_res, showCategory = top_n) +
    ggtitle("KEGG Pathway Enrichment") +
    theme_minimal()
  
  return(list(interactive_plot = interactive_plot, static_plot = static_plot, kegg_results = download_kegg_results))
}

# Run the function
#kegg_results <- generate_kegg_enrichment_plot(
#  gene_lists = gene_lists, 
#  de_results_df = bind_rows(efit.1.results_df, efit.2.results_df, efit.3.results_df), 
#  universe_entrez = universe_entrez
#)

#kegg_results$interactive_plot

generate_volcano_plot <- function(efit_results_df, contrast_name) {
  library(plotly)
  
  volcano.res <- efit_results_df
  volcano.res$group <- "NotSignificant"
  volcano.res[which(volcano.res$adj.P.value < 0.05 & volcano.res['logFC'] > 0.58), "group"] <- "Upregulated"
  volcano.res[which(volcano.res$adj.P.value < 0.05 & volcano.res['logFC'] < -0.58), "group"] <- "Downregulated"
  
  # Apply small offset to avoid log(1) = 0 issues
  volcano.res$neglog10p <- -log10(volcano.res$adj.P.value + 1e-10)  # Ensure no log(1)
  # Check if all y-values are small (<1) and set range accordingly
  if (max(volcano.res$neglog10p, na.rm = TRUE) < 1) {
    y_axis_range <- c(0, 1)  # Force range to 0-1
  } else {
    y_axis_range <- c(0, max(volcano.res$neglog10p, na.rm = TRUE) + 0.5)  # Dynamic range otherwise
  }
  plot_ly(
    data = volcano.res,
    x = volcano.res$logFC, 
    y = volcano.res$neglog10p,  # Use modified neg log10 p-values
    hoverinfo = "text", 
    text = paste(
      "Gene: ", volcano.res$SYMBOL, 
      "<br>Ensemble: ",volcano.res$ensembleID,
      "<br>logFC: ", round(volcano.res$logFC, 2),
      "<br>adj.P.val: ", round(volcano.res$adj.P.value, 2)
    ),
    mode = "markers", 
    color = volcano.res$group,
    colors = c(NotSignificant = "gray", Upregulated = "red", Downregulated = "blue")
  ) %>% 
    layout(
      title = paste(contrast_name, "(Adj. P-value vs. Fold-Change)"),
      xaxis = list(title = 'logFC'), 
      yaxis = list(
        title = '-log10(adj.PVal)',
        exponentformat = "none",  # Disable scientific notation
        separatethousands = TRUE,  # Keep numbers readable
        tickmode = "linear",  # Ensure even tick spacing
        range = y_axis_range 
      ),
      margin = list(l = 75, t = 150)
    )
  
  
}

generate_heatmap <- function(efit_results_df, lcpm_matrix, dge_list_NOISeq, title, num_genes = 50, fontsize_row = 12) {
  
  # Extract contrast name dynamically
  contrast_name <- unique(efit_results_df$contrast)[1]  # Ensure it's a single contrast
  contrast_groups <- unlist(strsplit(contrast_name, " - "))
  contrast_groups <- gsub("`", "", contrast_groups)                # remove backticks
  contrast_groups <- gsub("^X(?=\\d)", "", contrast_groups, perl = TRUE)  # remove leading 'X' before a digit
  
  
  # Filter samples that belong to these groups
  selected_samples <- dge_list_NOISeq$samples %>%
    filter(!!sym(report_params$group_var) %in% contrast_groups) %>%  # Select samples in the contrast
    select(SampleName, all_of(report_params$group_var))
  
  # Step 1: Select top 50 DE genes
  top50_sigOE_genes <- efit_results_df %>%
    filter(!is.na(SYMBOL)) %>%
    arrange(adj.P.value) %>%
    slice_head(n = num_genes) %>%
    select(ensembleID, SYMBOL)
  
  
  # Step 2: Subset the lcpm matrix for the selected samples
  lcpm_filtered <- lcpm_matrix[, colnames(lcpm_matrix) %in% selected_samples$SampleName]
  
  # Step 3: Ensure genes are properly filtered
  lcpm_filtered <- lcpm_filtered[rownames(lcpm_filtered) %in% top50_sigOE_genes$ensembleID, ]
  
  # Step 4: Transpose `lcpm_filtered`
  lcpm_filtered <- as.data.frame(t(lcpm_filtered))
  
  # Step 5: Ensure sample metadata matches
  sidecols <- selected_samples %>% select(SampleName, all_of(report_params$group_var))
  rownames(sidecols) <- sidecols$SampleName
  
  # Step 6: Reorder `lcpm_filtered` based on sorted `sidecols`
  sidecols <- sidecols %>% arrange(all_of(report_params$group_var))  
  lcpm_filtered <- t(lcpm_filtered)
  
  # Ensure order matches
  matching_indices <- match(rownames(sidecols), colnames(lcpm_filtered))  
  lcpm_filtered <- lcpm_filtered[, matching_indices, drop = FALSE]  
  rownames(lcpm_filtered) <- top50_sigOE_genes$SYMBOL  
  
  # Step 7: Define dynamic colors for the annotation
  group_levels <- unique(selected_samples[[report_params$group_var]])  
  color_palette <- brewer.pal(min(length(group_levels), 8), "Set2")  
  col_side_palette <- colorRampPalette(color_palette)  # Create a color palette function
  
  # Step 8: Create heatmap
  hm <- heatmaply(
    plot_method = "plotly",
    lcpm_filtered,
    show_dendrogram = FALSE,
    showticklabels = c(FALSE, TRUE),
    scale = "row",
    colors = c("darkblue", "white", "darkred"),
    col_side_colors = sidecols %>% select(all_of(report_params$group_var)),  
    xlab = "Samples",
    ylab = "Genes",
    main = title,
    margins = c(0, 0, 60, 0),
    fontsize_row = fontsize_row,
    cex.main = 1,
    col_side_palette = col_side_palette,  # Corrected palette function
    labRow = rownames(lcpm_filtered),
    side_color_colorbar_len = 0.2,   
    side_color_width = 0.05,
    height=700
  )
  
  return(hm)
}


#source("rnaseq_functions.R")
library(gtable)

# Function to generate MDS plots for multiple DGELists
generate_mds_plots <- function(dge_list, titles, top = 1000, group = "treatment") {
  if (length(dge_list) != length(titles)) {
    stop("Error: Length of DGEList and titles must be the same")
  }
  
  plots <- lapply(seq_along(dge_list), function(i) {
    dge <- dge_list[[i]]
    title <- titles[i]
    
    plot_MDS_2(lcpm = cpm(dge, log = TRUE), dge = dge, top = top, group = group) +
      theme(legend.position = "none") + 
      geom_point(size = 3) +
      ggtitle(title)
  })
  
  return(plots)
}

# Function to extract and adjust legend
get_legend <- function(p) {
  g <- ggplotGrob(
    p + theme(
      legend.key.width = unit(2, "cm"),  # Increase width of legend keys
      legend.spacing.x = unit(0.5, "cm"), # Increase spacing between legend items
      legend.text = element_text(size = 12),  # Adjust text size
      legend.title = element_text(size = 14, face = "bold") # Adjust title size
    )
  )
  legend <- gtable_filter(g, "guide-box")
  return(legend)
}

library(ggplot2)
library(base64enc)
library(htmltools)

download_button_png <- function(plot_object, output_name = "plot", width = 12, height = 6, dpi = 300) {
  # Step 1: Save the ggplot object as a temporary PNG file
  temp_file <- tempfile(fileext = ".png")
  ggsave(temp_file, plot = plot_object, width = width, height = height, dpi = dpi)
  
  # Step 2: Convert the image to Base64 encoding
  encoded_img <- base64enc::dataURI(file = temp_file, mime = "image/png")
  
  # Step 3: Create a styled HTML download button
  download_button <- HTML(paste0(
    '<a href="', encoded_img, '" download="', output_name, '.png" ',
    'class="btn btn-primary" style="padding:10px; font-size:16px; text-decoration:none; ',
    'color:white; background-color:#007BFF; border-radius:5px;">',
    '📥 Download ', output_name, '</a>'
  ))
  
  # Return the HTML download button
  return(download_button)
}

library(DT)
library(dplyr)

display_de_result_table <- function(efit, contrast_name = "Contrast") {
  # Step 1: Extract relevant columns from the efit object
  efit_results_df <- data.frame(
    ensembleID = rownames(efit$coefficients),  # Extract rownames as ensemble IDs
    logFC = efit$coefficients[, 1],           # Log fold change
    AveExpr = efit$Amean,                      # Average expression
    t = efit$t[, 1],                           # t-statistic
    P.value = efit$p.value[, 1],               # Raw p-value
    adj.P.value = p.adjust(efit$p.value[, 1], method = "BH"), # Adjusted p-value (BH)
    B = efit$lods[, 1],                         # Log-odds of differential expression
    contrast = contrast_name
  )
  
  # Step 2: Merge with gene symbols from efit$genes
  if (!is.null(efit$genes)) {
    efit_results_df <- merge(efit_results_df, efit$genes, by.x = "ensembleID", by.y = "ENSEMBL", all.x = TRUE)
  }
  
  # Step 3: Relocate SYMBOL and ensembleID for better display
  efit_results_df <- efit_results_df %>%
    relocate(SYMBOL) %>%
    relocate(ensembleID, .after = contrast)
  
  # Step 4: Create an interactive datatable
  datatable(
    efit_results_df %>% arrange(adj.P.value),  # Sort by adjusted p-value
    filter = "top",
    extensions = c('Buttons','Scroller'),
    caption = paste("Differential expression results for", contrast_name),
    options = list(
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel'),
      #pageLength = 10,
      autoWidth = TRUE,
      deferRender = TRUE,
      scrollY = 200,
      scroller = TRUE,
      order = list(list(which(names(efit_results_df) == "adj.P.value") - 1, "asc")),  # Sort by adj.P.value column
      columnDefs = list(list(className = 'dt-center', targets = "_all"))
    )
  ) %>%
    formatRound(columns = c("logFC", "AveExpr", "t", "P.value", "adj.P.value", "B"), digits = 3) %>%
    formatStyle(
      'logFC',
      background = styleInterval(c(-0.58, 0.58), c('lightblue', 'white', 'lightpink'))
    )
}

library(DT)
library(limma)
library(dplyr)
library(tidyr)

display_de_summary <- function(efit, lfc = 0.58, pval = 0.05, adjust_method = "BH", contrast) {
  # Step 1: Get summary table
  de_summary <- summary(decideTests(efit, lfc = lfc, adjust.method = adjust_method, p.value = pval))
  
  # Step 2: Convert to data frame
  de_summary_df <- as.data.frame(de_summary, responseName = "Number of Genes")
  
  # Step 3: Rename columns
  colnames(de_summary_df) <- c("Direction", "Contrast", "Number of Genes")
  
  # Step 4: Clean Contrast names
  de_summary_df$Contrast <- gsub("`", "", de_summary_df$Contrast)                # remove backticks
  de_summary_df$Contrast <- gsub("\\bX(?=[0-9])", "", de_summary_df$Contrast, perl = TRUE)  # remove X before a digit
  
  # Optional: Move 'Contrast' column to front
  de_summary_df <- de_summary_df %>% relocate(Contrast)
  
  # Step 5: Display datatable
  datatable(
    de_summary_df,
    rownames = FALSE,
    caption = paste("Summary of Differential Expression Results (abs(lfc) > 0.58, adj.p.val < 0.05) for contrast", contrast),
    options = list(
      dom = 'Bfrtip',
      buttons = c('copy'),
      pageLength = 10,
      autoWidth = TRUE
    )
  )
}


library(dplyr)

create_efit_results_df <- function(efit) {
  # Step 1: Extract contrast name dynamically
  contrast_name <- gsub("^X| X", "", colnames(efit$contrasts))
  
  # Step 2: Construct results data frame
  efit_results_df <- data.frame(
    ensembleID = rownames(efit$coefficients),  # Extract rownames as ensemble IDs
    logFC = efit$coefficients[, 1],           # Log fold change
    AveExpr = efit$Amean,                     # Average expression
    t = efit$t[, 1],                          # t-statistic
    P.value = efit$p.value[, 1],              # Raw p-value
    adj.P.value = p.adjust(efit$p.value[, 1], method = "BH"), # Adjusted p-value (BH)
    B = efit$lods[, 1],                       # Log-odds of differential expression
    contrast = contrast_name                  # Use extracted contrast name
  )
  
  # Step 3: Merge with gene symbols from efit$genes
  if (!is.null(efit$genes)) {
    efit_results_df <- merge(efit_results_df, efit$genes, by.x = "ensembleID", by.y = "ENSEMBL", all.x = TRUE)
  }
  
  # Step 4: Reorder columns for readability
  efit_results_df <- efit_results_df %>% relocate(SYMBOL, .after = ensembleID)
  
  return(efit_results_df)
}

library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)

# Function to add Entrez IDs to a single efit results df
add_entrez_ids <- function(results_df) {
  # Extract Ensembl IDs
  ensembl_ids <- results_df$ensembleID
  
  # Convert to Entrez ID
  entrez_ids <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")  # Use "first" to handle multiple mappings
  
  # Add Entrez IDs to DataFrame
  results_df <- results_df %>%
    mutate(EntrezID = entrez_ids)
  
  return(results_df)
}

library(plotly)
library(RColorBrewer)
library(car)

plot_pca <- function(dge, title, grp_var=report_params$group_var, show_legend = TRUE, combine_plots = FALSE) {
  # Extract log-transformed CPM values
  PCA_DATA <- t(get_log_matrix(dge))
  rownames(PCA_DATA) <- dge$samples$sample_id
  numsonly <- as.data.frame(PCA_DATA)
  numsonly <- as.data.frame(lapply(numsonly, as.numeric))  # Ensure numeric values
  
  # Perform PCA
  pca_res <- prcomp(numsonly, center = TRUE, scale. = FALSE, rank. = 2)
  scores <- pca_res$x
  
  # Extract group information
  group_labels <- dge$samples[[grp_var]]  
  colors <- RColorBrewer::brewer.pal(length(unique(group_labels)), "Set1")
  
  # Initialize plotly object
  fig <- plot_ly()
  
  # Function to plot ellipses
  plot_ellipse_function <- function(target) {
    target_data <- as.data.frame(scores[group_labels == target, 1:2])
    target_color <- colors[which(unique(group_labels) == target)]
    
    ellipse_coords <- car::dataEllipse(target_data$PC1, target_data$PC2, 
                                       levels = 0.68, plot.points = FALSE, add = TRUE, draw = FALSE)
    
    fig <<- fig %>%
      add_polygons(x = ellipse_coords[,1], y = ellipse_coords[,2],
                   line = list(color = target_color, dash = "dot"),
                   fillcolor = target_color, opacity = 0.3,
                   showlegend = FALSE, hoverinfo = "skip")
  }
  
  # Function to plot points
  plot_points_function <- function(target) {
    target_data <- as.data.frame(scores[group_labels == target, 1:2])
    rownames(target_data) <- rownames(PCA_DATA)[group_labels == target]
    target_color <- colors[which(unique(group_labels) == target)]
    
    text_data <- as.character(rownames(target_data))
    PC1 <- as.numeric(target_data$PC1)
    PC2 <- as.numeric(target_data$PC2)
    
    fig <<- fig %>%
      add_trace(data = target_data, x = ~PC1, y = ~PC2, 
                type = "scatter", mode = "markers", name = target,
                marker = list(color = target_color),
                text = ~text_data, hoverinfo = "text", showlegend = show_legend)
  }
  
  # Add ellipses and points
  invisible(lapply(unique(group_labels), plot_ellipse_function))
  invisible(lapply(unique(group_labels), plot_points_function))
  
  # Add title
  # Get % variance explained
  percentVar <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
  
  # Add axis titles with % variance
  fig <- fig %>% layout(
    title = title,
    xaxis = list(title = paste0("PC1 (", percentVar[1], "%)")),
    yaxis = list(title = paste0("PC2 (", percentVar[2], "%)"))
  )
  
  return(fig)
}
# Define function to subset samples by contrast and run PCA
plot_pca_by_contrast <- function(dge_list, contrast_name, group_var) {
  # Extract groups from contrast name
  groups <- unlist(strsplit(gsub("efit_|_results_df", "", contrast_name), "_vs_"))
  
  # Subset samples in the DGEList based on the groups in the contrast
  filtered_dge_list <- lapply(dge_list, function(dge) {
    sample_subset <- dge$samples[dge$samples[[report_params$group_var]] %in% groups, ]
    dge_filtered <- dge
    dge_filtered$counts <- dge$counts[, colnames(dge$counts) %in% sample_subset$SampleName]
    dge_filtered$samples <- sample_subset
    return(dge_filtered)
  })
  
  # Generate PCA plots for raw, NOISeq-filtered, and edgeR-filtered data
  pca_plots <- mapply(plot_pca, filtered_dge_list, 
                      title = paste0("PCA for ", gsub("efit_", "", contrast_name)),  
                      grp_var = report_params$group_var, 
                      show_legend = c(TRUE, FALSE, FALSE), 
                      SIMPLIFY = FALSE)
  
  # Combine plots into one row
  combined_plot <- subplot(pca_plots[[1]], pca_plots[[2]], pca_plots[[3]], 
                           nrows = 1, shareY = FALSE)
  
  # Add annotations
  annotations = list( 
    list( 
      x = 0.15, y = 0.95,  
      text = "PCA Before Filtering", xref = "paper", yref = "paper", 
      xanchor = "center", yanchor = "bottom", showarrow = FALSE 
    ),  
    list( 
      x = 0.5, y = 0.95,  
      text = "PCA After Filtering: NOISeq", xref = "paper", yref = "paper",  
      xanchor = "center", yanchor = "bottom", showarrow = FALSE
    ),
    list( 
      x = 0.85, y = 0.95,  
      text = "PCA After Filtering: edgeR", xref = "paper", yref = "paper",  
      xanchor = "center", yanchor = "bottom", showarrow = FALSE
    )
  )
  
  # Apply annotations
  combined_plot <- combined_plot %>% layout(
    annotations = annotations,
    title  = list(
      text = paste0("PCA for ", gsub("efit_", "", contrast_name)),
      y = 0.99,          # move title slightly lower (1 is top)
      yanchor = "top"    # anchor from the top so y behaves as expected
    )
  )
  
  return(combined_plot)
}

# Function to extract top DE genes and convert to Entrez IDs
top_DE_entrezIDs <- function(df, direction = "up") {
  if (direction == "up") {
    filtered_genes <- df %>% filter(logFC > 0.58 & adj.P.value < 0.05)
  } else if (direction == "down") {
    filtered_genes <- df %>% filter(logFC < -0.58 & adj.P.value < 0.05)
  } else {
    stop("Invalid direction. Choose either 'up' or 'down'.")
  }
  
  # Return empty vector if no significant genes
  if (nrow(filtered_genes) == 0) {
    return(character(0))  # Empty list instead of failure
  }
  
  entrez_ids <- mapIds(
    annotation_obj,
    keys = filtered_genes$ensembleID,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  return(entrez_ids)
}

# Function to extract top DE genes and convert to Entrez IDs
non_DE_entrezIDs <- function(df) {
  filtered_genes <- df %>% filter(adj.P.value > 0.05|abs(logFC) < 0.58)
  
  # Return empty vector if no significant genes
  if (nrow(filtered_genes) == 0) {
    return(character(0))  # Empty list instead of failure
  }
  
  entrez_ids <- mapIds(
    annotation_obj,
    keys = filtered_genes$ensembleID,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  return(entrez_ids)
}

download_de_result_table <- function(efit, contrast_name = "Contrast") {
  # Step 1: Extract relevant columns from the efit object
  efit_results_df <- data.frame(
    ensembleID = rownames(efit$coefficients),  # Extract rownames as ensemble IDs
    logFC = efit$coefficients[, 1],           # Log fold change
    AveExpr = efit$Amean,                      # Average expression
    t = efit$t[, 1],                           # t-statistic
    P.value = efit$p.value[, 1],               # Raw p-value
    adj.P.value = p.adjust(efit$p.value[, 1], method = "BH"), # Adjusted p-value (BH)
    B = efit$lods[, 1],                         # Log-odds of differential expression
    contrast = contrast_name
  )
  
  # Step 2: Merge with gene symbols from efit$genes
  if (!is.null(efit$genes)) {
    efit_results_df <- merge(efit_results_df, efit$genes, by.x = "ensembleID", by.y = "ENSEMBL", all.x = TRUE)
  }
  
  # Step 3: Relocate SYMBOL and ensembleID for better display
  efit_results_df <- efit_results_df %>%
    relocate(SYMBOL) %>%
    relocate(ensembleID, .after = contrast)
  return(efit_results_df)
}

plot_one_pca_by_contrast <- function(dge_list, contrast_name, group_var) {
  # Extract groups from contrast name
  groups <- unlist(strsplit(gsub("efit_|_results_df", "", contrast_name), "_vs_"))
  
  # Subset samples based on contrast groups
  sample_subset <- dge_list$samples[dge_list$samples[[group_var]] %in% groups, ]
  dge_filtered <- dge_list
  selected_samples <- sample_subset$SampleName
  
  dge_filtered$counts <- dge_list$counts[, colnames(dge_list$counts) %in% selected_samples]
  dge_filtered$samples <- sample_subset
  
  if (!is.null(dge_list$E_corrected)) {
    dge_filtered$E_corrected <- dge_list$E_corrected[, selected_samples]
  }
  
  # Generate PCA plot
  pca_plot <- plot_pca(dge_filtered,
                       title = paste0("PCA for ", gsub("efit_", "", contrast_name)),
                       grp_var = group_var,
                       show_legend = TRUE)
  
  return(pca_plot)
}


plot_pca_combined <- function(dge_list, grp_var = report_params$group_var, annotation_labels = NULL) {
  show_legend_flags <- c(TRUE, rep(FALSE, length(dge_list) - 1))
  pca_plots <- mapply(plot_pca, dge_list, grp_var = grp_var, show_legend = show_legend_flags, title = "", SIMPLIFY = FALSE)
  
  combined_plot <- subplot(pca_plots, nrows = 1, shareY = FALSE)
  
  if (!is.null(annotation_labels)) {
    annotations <- lapply(seq_along(annotation_labels), function(i) {
      list(
        x = (i - 0.5) / length(annotation_labels),
        y = 1,
        text = annotation_labels[i],
        xref = "paper",
        yref = "paper",
        xanchor = "center",
        yanchor = "bottom",
        showarrow = FALSE
      )
    })
    combined_plot <- combined_plot %>% layout(annotations = annotations)
  }
  
  return(combined_plot)
}


get_log_matrix <- function(dge) {
  if (!is.null(dge$E_corrected)) {
    return(dge$E_corrected)
  } else if (!is.null(dge$counts)) {
    return(cpm(dge, log = TRUE))
  } else {
    stop("Input DGE object must contain either 'counts' or 'E_corrected'.")
  }
}

