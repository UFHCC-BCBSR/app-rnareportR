library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(plotly)
library(heatmaply)
library(org.Mm.eg.db)

generate_enrichment_plot <- function(gene_lists, de_results_df, universe_entrez, ont_category, 
                                     significance_threshold = 0.05, top_n = 10, annotation_db) {
  
  library(BiocParallel)
  
  # Get n_cores from report_params if available, otherwise default to 1
  n_cores <- if (exists("report_params") && "n_cores" %in% names(report_params)) {
    as.numeric(report_params$n_cores)
  } else {
    1
  }
  
  pvalue_cutoff <- if (isTRUE(report_params$download_nonsig_enrich)) 1 else significance_threshold
  annotation_obj <- get(annotation_db, envir = asNamespace(annotation_db))
  names(gene_lists) <- gsub("efit_|_results_df", "", names(gene_lists))
  
  if (is.null(names(gene_lists))) stop("Each gene list must be named!")
  contrast_order <- names(gene_lists)
  
  # PRE-COMPUTE ENTREZ IDs for ALL DE results (do this ONCE)
  de_results_df <- de_results_df %>%
    mutate(ENTREZID = mapIds(annotation_obj, keys = ensembleID, 
                             column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first"))
  
  # PARALLELIZED: Run GO enrichment using compareCluster() with BiocParallel
  message(sprintf("Running GO enrichment with %d cores...", n_cores))
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
    pvalueCutoff = pvalue_cutoff,
  )
  
  # Handle no results case
  if (is.null(formula_res) || nrow(formula_res@compareClusterResult) == 0) {
    message_plot <- ggplot() +
      annotate("text", x = 1, y = 1, 
               label = paste0("No significant GO enrichment found\n(", ont_category, ")"), 
               size = 6, hjust = 0.5) +
      theme_void() +
      ggtitle(paste("GO Term Enrichment (", ont_category, ")", sep = ""))
    
    return(list(
      interactive_plot = ggplotly(message_plot),
      static_plot = message_plot,
      go_results = NULL
    ))
  }
  
  formula_res@compareClusterResult$Cluster <- factor(
    formula_res@compareClusterResult$Cluster,
    levels = contrast_order
  )
  
  filtered_results <- subset(
    formula_res@compareClusterResult,
    p.adjust <= significance_threshold
  )
  
  if (nrow(filtered_results) == 0) {
    message_plot <- ggplot() +
      annotate("text", x = 1, y = 1, 
               label = paste0("No significant GO enrichment found\n(", ont_category, ")"), 
               size = 6, hjust = 0.5) +
      theme_void() +
      ggtitle(paste("GO Term Enrichment (", ont_category, ")", sep = ""))
    
    return(list(
      interactive_plot = ggplotly(message_plot),
      static_plot = message_plot,
      go_results = NULL
    ))
  }
  
  # PRE-BUILD lookup table: ENTREZID -> SYMBOL (do this ONCE)
  unique_entrez <- unique(unlist(strsplit(filtered_results$geneID, "/")))
  entrez_to_symbol <- setNames(
    mapIds(annotation_obj, keys = unique_entrez, column = "SYMBOL", 
           keytype = "ENTREZID", multiVals = "first"),
    unique_entrez
  )
  
  # Function to get gene symbols (uses pre-computed lookups)
  get_gene_symbols <- function(i, results_df, de_df, entrez_symbol_map) {
    gene_list <- results_df$geneID[i]
    contrast_full <- as.character(results_df$Cluster[i])
    contrast_base <- sub("\\.(up|down)$", "", contrast_full)
    direction <- sub("^.*\\.", "", contrast_full)
    entrez_ids <- unlist(strsplit(gene_list, "/"))
    
    de_sub <- de_df %>%
      filter(grepl(contrast_base, contrast),
             ENTREZID %in% entrez_ids,
             case_when(
               direction == "up" ~ logFC > 0,
               direction == "down" ~ logFC < 0
             ))
    
    top_entrez <- de_sub %>%
      arrange(adj.P.value) %>%
      slice_head(n = 20) %>%
      pull(ENTREZID)
    
    gene_symbols <- na.omit(entrez_symbol_map[top_entrez])
    
    if (length(gene_symbols) > 0) {
      paste(gene_symbols, collapse = "<br>")
    } else {
      NA_character_
    }
  }
  
  # PARALLELIZED: Get gene symbols using BiocParallel (works on HPC)
  if (n_cores > 1) {
    message(sprintf("Computing gene symbols with %d cores...", n_cores))
    filtered_results$GeneSymbols <- bplapply(
      seq_len(nrow(filtered_results)),
      get_gene_symbols,
      results_df = filtered_results,
      de_df = de_results_df,
      entrez_symbol_map = entrez_to_symbol,
      BPPARAM = MulticoreParam(workers = n_cores)
    ) %>% unlist()
  } else {
    filtered_results$GeneSymbols <- sapply(
      seq_len(nrow(filtered_results)),
      get_gene_symbols,
      results_df = filtered_results,
      de_df = de_results_df,
      entrez_symbol_map = entrez_to_symbol
    )
  }
  
  all_results <- formula_res@compareClusterResult
  
  # If download_nonsig_enrich, compute symbols for ALL results
  if (isTRUE(report_params$download_nonsig_enrich)) {
    all_unique_entrez <- unique(unlist(strsplit(all_results$geneID, "/")))
    all_entrez_to_symbol <- setNames(
      mapIds(annotation_obj, keys = all_unique_entrez, column = "SYMBOL", 
             keytype = "ENTREZID", multiVals = "first"),
      all_unique_entrez
    )
    
    if (n_cores > 1) {
      all_results$GeneSymbols <- bplapply(
        seq_len(nrow(all_results)),
        get_gene_symbols,
        results_df = all_results,
        de_df = de_results_df,
        entrez_symbol_map = all_entrez_to_symbol,
        BPPARAM = MulticoreParam(workers = n_cores)
      ) %>% unlist()
    } else {
      all_results$GeneSymbols <- sapply(
        seq_len(nrow(all_results)),
        get_gene_symbols,
        results_df = all_results,
        de_df = de_results_df,
        entrez_symbol_map = all_entrez_to_symbol
      )
    }
  }
  
  download_go_results <- if (isTRUE(report_params$download_nonsig_enrich)) {
    all_results %>% select(Cluster, Description, p.adjust, GeneSymbols, everything())
  } else {
    filtered_results %>% select(Cluster, Description, p.adjust, GeneSymbols, everything())
  }
  
  # ... rest of plotting code unchanged ...
  top_GO_terms <- filtered_results %>%
    group_by(Cluster) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    pull(Description) %>%
    unique()
  
  formula_res@compareClusterResult <- filtered_results %>%
    filter(Description %in% top_GO_terms)
  
  formula_res@compareClusterResult <- formula_res@compareClusterResult %>%
    mutate(GeneRatio = sapply(strsplit(as.character(GeneRatio), "/"), 
                              function(x) as.numeric(x[1]) / as.numeric(x[2])))
  
  reorder_GO_terms <- function(df) {
    term_matrix <- table(df$Description, df$Cluster)
    if (nrow(term_matrix) > 1 && length(unique(df$Description)) > 1) {
      term_dist <- dist(term_matrix, method = "binary")
      term_hclust <- hclust(term_dist, method = "ward.D2")
      term_order <- rownames(term_matrix)[term_hclust$order]
    } else {
      term_order <- unique(df$Description)
    }
    df$Description <- factor(df$Description, levels = term_order)
    return(df)
  }
  
  formula_res@compareClusterResult <- reorder_GO_terms(formula_res@compareClusterResult)
  
  bin_breaks <- c(0, 0.01, 0.05, 0.10, 
                  max(formula_res@compareClusterResult$GeneRatio, na.rm = TRUE) + 0.01)
  bin_labels <- c("≤0.01", "0.01 - 0.05", "0.05 - 0.10", "≥0.10")
  
  formula_res@compareClusterResult <- formula_res@compareClusterResult %>%
    mutate(GeneRatioCategory = cut(GeneRatio, breaks = bin_breaks, labels = bin_labels, 
                                   include.lowest = TRUE, right = FALSE))
  
  size_mapping <- c("≤0.01" = 2, "0.01 - 0.05" = 4, "0.05 - 0.10" = 6, "≥0.10" = 8)
  
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
  
  p <- ggplot(formula_res@compareClusterResult, aes(
    x = Cluster, y = plot_label, size = GeneRatioCategory,
    color = p.adjust, text = tooltip_text
  )) +
    geom_point(alpha = 0.8) +
    scale_size_manual(name = "Gene Ratio", values = size_mapping) +
    scale_color_gradient(low = "red", high = "blue",
                         limits = c(min(formula_res@compareClusterResult$p.adjust, na.rm = TRUE),
                                    max(formula_res@compareClusterResult$p.adjust, na.rm = TRUE)),
                         name = "p.adjust") +
    guides(color = guide_colorbar(title = "p.adjust")) +
    ggtitle(paste("GO Term Enrichment (", ont_category, ")", sep = "")) +
    xlab("DE Gene list") + ylab("GO Term") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.title.x = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  interactive_plot <- ggplotly(p, tooltip = "text") %>%
    layout(legend = list(title = list(text = "Gene Ratio")))
  
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

if (is.null(kegg_res) || nrow(kegg_res@compareClusterResult) == 0) {
  message("No significant KEGG pathways found.")
  return(list(
    interactive_plot = NULL, 
    static_plot = NULL, 
    kegg_results = data.frame()
  ))
}

# Ensure clusters are ordered correctly
kegg_res@compareClusterResult$Cluster <- factor(
  kegg_res@compareClusterResult$Cluster,
  levels = contrast_order
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
  
  # Calculate -log10(p-value)
  volcano.res$neglog10p <- -log10(volcano.res$adj.P.value + 1e-300)  # Avoid log(0)
  
  # Better y-axis range - always show significance threshold
  sig_threshold <- -log10(0.05)  # = 1.3
  max_y <- max(volcano.res$neglog10p, na.rm = TRUE)
  y_axis_range <- c(0, max(sig_threshold + 0.5, max_y + 0.5))  # Always show above threshold
  
  # Count significant genes
  n_up <- sum(volcano.res$group == "Upregulated")
  n_down <- sum(volcano.res$group == "Downregulated")
  n_total <- nrow(volcano.res)
  
  # Create plot
  p <- plot_ly(
    data = volcano.res,
    x = ~logFC,
    y = ~neglog10p,
    type = "scatter",
    mode = "markers",
    color = ~group,
    colors = c(NotSignificant = "gray", Upregulated = "red", Downregulated = "blue"),
    marker = list(size = 5, opacity = 0.6),
    hoverinfo = "text",
    text = ~paste0(
      "<b>", SYMBOL, "</b>",
      "<br>Ensembl: ", ensembleID,
      "<br>logFC: ", round(logFC, 3),
      "<br>adj.P.value: ", formatC(adj.P.value, format = "e", digits = 2),
      "<br>-log10(adj.P): ", round(neglog10p, 2)
    )
  ) %>%
    layout(
      title = list(
        text = paste0(gsub("_", " ", contrast_name), 
                      "<br><sub>", n_up, " upregulated | ", 
                      n_down, " downregulated | ",
                      n_total - n_up - n_down, " not significant</sub>"),
        font = list(size = 16)
      ),
      xaxis = list(
        title = "log<sub>2</sub> Fold Change",
        zeroline = TRUE,
        zerolinewidth = 1,
        zerolinecolor = "black"
      ),
      yaxis = list(
        title = "-log<sub>10</sub>(Adjusted P-value)",
        range = y_axis_range
      ),
      shapes = list(
        # Horizontal line at significance threshold
        list(
          type = "line",
          x0 = min(volcano.res$logFC, na.rm = TRUE) - 1,
          x1 = max(volcano.res$logFC, na.rm = TRUE) + 1,
          y0 = sig_threshold,
          y1 = sig_threshold,
          line = list(color = "red", width = 1, dash = "dash")
        ),
        # Vertical line at logFC = 0.58
        list(
          type = "line",
          x0 = 0.58,
          x1 = 0.58,
          y0 = 0,
          y1 = max(y_axis_range),
          line = list(color = "blue", width = 1, dash = "dash")
        ),
        # Vertical line at logFC = -0.58
        list(
          type = "line",
          x0 = -0.58,
          x1 = -0.58,
          y0 = 0,
          y1 = max(y_axis_range),
          line = list(color = "blue", width = 1, dash = "dash")
        )
      ),
      annotations = if (n_up == 0 && n_down == 0) {
        list(
          list(
            x = 0.5,
            y = 0.95,
            xref = "paper",
            yref = "paper",
            text = "<b>No significant genes</b><br>(adj.P < 0.05 & |logFC| > 0.58)",
            showarrow = FALSE,
            font = list(size = 14, color = "darkred"),
            bgcolor = "rgba(255, 200, 200, 0.8)",
            bordercolor = "red",
            borderwidth = 2
          )
        )
      } else {
        list()
      },
      margin = list(l = 75, t = 100, b = 75)
    )
  
  return(p)
}

#' Generate delta plots showing per-subject changes for paired data
#' @param efit_results_dfs List of efit results dataframes
#' @param lcpm_matrix Log-CPM expression matrix
#' @param dge_list_filt Filtered DGEList object with sample metadata
#' @param top_n Number of top genes to plot (default 10)
generate_delta_plots <- function(efit_results_dfs, lcpm_matrix, dge_list_filt, top_n = 10) {
  
  delta_plot_list <- lapply(seq_along(efit_results_dfs), function(i) {
    
    # Extract contrast name and groups dynamically
    contrast_name <- unique(efit_results_dfs[[i]]$contrast)[1]
    contrast_groups <- unlist(strsplit(contrast_name, " - "))
    contrast_groups <- gsub("`", "", contrast_groups)
    contrast_groups <- gsub("^X(?=\\d)", "", contrast_groups, perl = TRUE)
    
    # Reference group is on the right side of " - " (index 2), comparison is on left (index 1)
    ref_group <- contrast_groups[2]
    comp_group <- contrast_groups[1]
    
    top_genes <- efit_results_dfs[[i]] %>%
      arrange(P.value) %>%
      head(top_n)
    
    delta_data <- do.call(rbind, lapply(1:nrow(top_genes), function(g) {
      gene_id <- top_genes$ensembleID[g]
      gene_name <- if(!is.na(top_genes$SYMBOL[g])) top_genes$SYMBOL[g] else gene_id
      
      expr <- lcpm_matrix[gene_id, ]
      subjects <- unique(dge_list_filt$samples$subject)
      
      deltas <- sapply(subjects, function(subj) {
        samples <- dge_list_filt$samples
        ref_sample <- samples$SampleName[samples$subject == subj & samples[[report_params$group_var]] == ref_group]
        comp_sample <- samples$SampleName[samples$subject == subj & samples[[report_params$group_var]] == comp_group]
        
        if (length(ref_sample) > 0 && length(comp_sample) > 0) {
          expr[comp_sample] - expr[ref_sample]
        } else {
          NA
        }
      })
      
      data.frame(
        gene = paste0(gene_name),
        subject = subjects,
        delta = deltas
      )
    }))
    
    ggplot(delta_data, aes(x = factor(subject), y = delta)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_point(size = 3, alpha = 0.7, color = "#3498db") +
      geom_segment(aes(x = factor(subject), xend = factor(subject), 
                       y = 0, yend = delta), 
                   color = "#3498db", alpha = 0.5) +
      facet_wrap(~ gene, scales = "free_y", ncol = 2) +
      labs(
        title = paste0("Per-subject changes (top ", top_n, " genes)\n",
                       gsub("efit_|_results_df","",names(efit_results_dfs)[i])),
        x = "Subject",
        y = paste0("Log2 fold change (", comp_group, " vs ", ref_group, ")")
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 8)
      )
  })
  
  return(delta_plot_list)
}

#' Generate interactive trajectory plots with searchable dropdown
#' @param efit_results_dfs List of efit results dataframes
#' @param lcpm_matrix Log-CPM expression matrix
#' @param dge_list_filt Filtered DGEList object with sample metadata
#' @param max_genes Maximum genes to include (default 100, supports type-to-filter)
generate_trajectory_plots <- function(efit_results_dfs, lcpm_matrix, dge_list_filt, max_genes = 100) {
  
  trajectory_plot_list <- lapply(seq_along(efit_results_dfs), function(i) {
    
    # Extract contrast name and groups dynamically
    contrast_name <- unique(efit_results_dfs[[i]]$contrast)[1]
    contrast_groups <- unlist(strsplit(contrast_name, " - "))
    contrast_groups <- gsub("`", "", contrast_groups)
    contrast_groups <- gsub("^X(?=\\d)", "", contrast_groups, perl = TRUE)
    
    # Reference group is on the right side of " - " (index 2), comparison is on left (index 1)
    ref_group <- contrast_groups[2]  # Left side of plot
    comp_group <- contrast_groups[1]  # Right side of plot
    timepoint_order <- c(ref_group, comp_group)
    
    # Get more genes for better search capability
    top_genes <- efit_results_dfs[[i]] %>%
      filter(!is.na(SYMBOL)) %>%
      arrange(P.value) %>%
      head(max_genes)
    
    all_traces <- list()
    buttons <- list()
    
    for(g in 1:nrow(top_genes)) {
      gene_id <- top_genes$ensembleID[g]
      gene_name <- top_genes$SYMBOL[g]
      
      expr <- lcpm_matrix[gene_id, ]
      subjects <- unique(dge_list_filt$samples$subject)
      
      traj_data <- do.call(rbind, lapply(subjects, function(subj) {
        samples <- dge_list_filt$samples
        ref_sample <- samples$SampleName[samples$subject == subj & samples[[report_params$group_var]] == ref_group]
        comp_sample <- samples$SampleName[samples$subject == subj & samples[[report_params$group_var]] == comp_group]
        
        if (length(ref_sample) > 0 && length(comp_sample) > 0) {
          data.frame(
            subject = subj,
            timepoint = factor(c(ref_group, comp_group), levels = timepoint_order),
            expression = c(expr[ref_sample], expr[comp_sample]),
            stringsAsFactors = FALSE
          )
        }
      }))
      
      # Create traces for this gene (one per subject)
      for(subj in subjects) {
        subj_data <- traj_data[traj_data$subject == subj, ]
        
        if(nrow(subj_data) > 0) {
          trace <- list(
            x = as.character(subj_data$timepoint),
            y = subj_data$expression,
            type = "scatter",
            mode = "lines+markers",
            name = paste("Subject", subj),
            visible = (g == 1),  # Only first gene visible initially
            line = list(width = 2),
            marker = list(size = 8),
            hovertemplate = paste0(
              "Subject: ", subj, "<br>",
              "Timepoint: %{x}<br>",
              "Expression: %{y:.2f}<br>",
              "<extra></extra>"
            )
          )
          
          all_traces[[length(all_traces) + 1]] <- trace
        }
      }
      
      # Create button for this gene
      n_subjects <- length(subjects)
      visible_vec <- rep(FALSE, length(subjects) * nrow(top_genes))
      start_idx <- (g - 1) * n_subjects + 1
      end_idx <- g * n_subjects
      visible_vec[start_idx:end_idx] <- TRUE
      
      buttons[[g]] <- list(
        method = "update",
        args = list(
          list(visible = visible_vec),
          list(title = paste0(gene_name, " - Paired trajectories"))
        ),
        label = gene_name
      )
    }
    
    # Create plotly figure
    fig <- plot_ly()
    
    # Add all traces
    for(trace in all_traces) {
      fig <- fig %>% add_trace(
        x = trace$x,
        y = trace$y,
        type = trace$type,
        mode = trace$mode,
        name = trace$name,
        visible = trace$visible,
        line = trace$line,
        marker = trace$marker,
        hovertemplate = trace$hovertemplate
      )
    }
    
    # Add dropdown menu
    fig <- fig %>%
      layout(
        title = list(
          text = paste0(top_genes$SYMBOL[1], " - Paired trajectories<br>",
                        "<sub>", gsub("efit_|_results_df","",names(efit_results_dfs)[i]), 
                        " (Type to search ", nrow(top_genes), " genes)</sub>"),
          font = list(size = 16)
        ),
        xaxis = list(
          title = "Timepoint",
          categoryorder = "array",
          categoryarray = timepoint_order
        ),
        yaxis = list(title = "Expression (log-CPM)"),
        updatemenus = list(
          list(
            active = 0,
            type = "dropdown",
            direction = "down",
            x = 0.01,
            y = 1.15,
            xanchor = "left",
            yanchor = "top",
            buttons = buttons
          )
        ),
        showlegend = TRUE,
        hovermode = "closest"
      )
    
    return(fig)
  })
  
  return(trajectory_plot_list)
}

#' Generate subject-ordered heatmap for paired data
#' @param efit_results_dfs List of efit results dataframes
#' @param lcpm_matrix Log-CPM expression matrix
#' @param dge_list_filt Filtered DGEList object with sample metadata
#' @param num_genes Number of genes to show (default 50)
generate_subject_heatmaps <- function(efit_results_dfs, lcpm_matrix, dge_list_filt, num_genes = 50) {
  
  heatmap_list <- lapply(seq_along(efit_results_dfs), function(i) {
    
    # Extract contrast name and groups dynamically
    contrast_name <- unique(efit_results_dfs[[i]]$contrast)[1]
    contrast_groups <- unlist(strsplit(contrast_name, " - "))
    contrast_groups <- gsub("`", "", contrast_groups)
    contrast_groups <- gsub("^X(?=\\d)", "", contrast_groups, perl = TRUE)
    
    # Reference group is on the right side of " - " (index 2), comparison is on left (index 1)
    ref_group <- contrast_groups[2]
    comp_group <- contrast_groups[1]
    group_order <- c(ref_group, comp_group)
    
    # Get top genes with gene symbols
    top_genes_df <- efit_results_dfs[[i]] %>%
      filter(!is.na(SYMBOL)) %>%
      arrange(P.value) %>%
      head(num_genes)
    
    top_genes <- top_genes_df$ensembleID
    gene_symbols <- top_genes_df$SYMBOL
    
    lcpm_top <- lcpm_matrix[top_genes, ]
    
    # Replace Ensemble IDs with gene symbols
    rownames(lcpm_top) <- gene_symbols
    
    # Filter samples to only those in the contrast groups
    samples_in_contrast <- dge_list_filt$samples %>%
      filter(!!sym(report_params$group_var) %in% contrast_groups)
    
    # Subset lcpm matrix to only samples in contrast
    lcpm_top <- lcpm_top[, colnames(lcpm_top) %in% samples_in_contrast$SampleName]
    
    annotation_col <- data.frame(
      subject = factor(samples_in_contrast$subject),
      group = factor(samples_in_contrast[[report_params$group_var]], levels = group_order)
    )
    rownames(annotation_col) <- samples_in_contrast$SampleName
    
    # Ensure lcpm_top columns match annotation_col rows
    lcpm_top <- lcpm_top[, rownames(annotation_col)]
    
    sample_order <- order(annotation_col$subject, annotation_col$group)
    
    n_subjects <- length(unique(annotation_col$subject))
    subject_colors <- setNames(
      rainbow(n_subjects, s = 0.6, v = 0.8), 
      levels(annotation_col$subject)
    )
    
    # Dynamic group colors
    group_colors <- setNames(
      c("#3498db", "#e74c3c"),
      group_order
    )
    
    ann_colors <- list(
      subject = subject_colors,
      group = group_colors
    )
    
    pheatmap(
      lcpm_top[, sample_order],
      scale = "row",
      annotation_col = annotation_col[sample_order, ],
      annotation_colors = ann_colors,
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      show_colnames = TRUE,
      show_rownames = TRUE,
      fontsize_row = 6,
      fontsize_col = 8,
      main = paste0("Top ", num_genes, " genes - grouped by subject\n", 
                    gsub("efit_|_results_df","",names(efit_results_dfs)[i])),
      color = colorRampPalette(c("blue", "white", "red"))(100)
    )
  })
  
  return(heatmap_list)
}

generate_heatmap <- function(efit_results_df, lcpm_matrix, dge_list_filt, title, num_genes = 50, fontsize_row = 12) {
  
  # Extract contrast name dynamically
  contrast_name <- unique(efit_results_df$contrast)[1]  # Ensure it's a single contrast
  contrast_groups <- unlist(strsplit(contrast_name, " - "))
  contrast_groups <- gsub("`", "", contrast_groups)                # remove backticks
  contrast_groups <- gsub("^X(?=\\d)", "", contrast_groups, perl = TRUE)  # remove leading 'X' before a digit
  
  
  # Filter samples that belong to these groups
  selected_samples <- dge_list_filt$samples %>%
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
  
  # Step 6: Arrange `lcpm_filtered` based on sorted `sidecols`
  sidecols <- sidecols %>% arrange(all_of(report_params$group_var))
  
  # transpose so that now samples are columns and genes are rows
  lcpm_filtered <- t(lcpm_filtered)
  
  # Ensure sample order matches data
  matching_indices <- match(rownames(sidecols), colnames(lcpm_filtered))  
  lcpm_filtered <- lcpm_filtered[, matching_indices, drop = FALSE]  
  
  # Add rownames to lcpm matrix subset ensuring order matches
  gene_indices <- match(top50_sigOE_genes$ensembleID, rownames(lcpm_filtered))
  valid_genes <- !is.na(gene_indices)
  lcpm_filtered <- lcpm_filtered[gene_indices[valid_genes], ]
  gene_symbols <- top50_sigOE_genes$SYMBOL[valid_genes]
  rownames(lcpm_filtered) <- gene_symbols

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
      autoWidth = TRUE,
      deferRender = TRUE,
      scrollY = 200,
      scroller = TRUE,
      order = list(list(which(names(efit_results_df) == "adj.P.value") - 1, "asc")),  # Sort by adj.P.value column
      columnDefs = list(list(className = 'dt-center', targets = "_all"))
    )
  ) %>%
    # Format numeric columns (not p-values) with 3 decimal places
    formatRound(columns = c("logFC", "AveExpr", "t", "B"), digits = 3) %>%
    # Format p-values with scientific notation (3 significant figures)
    formatSignif(columns = c("P.value", "adj.P.value"), digits = 3) %>%
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

# Improved plot_pca - Faster and cleaner
plot_pca <- function(dge, title = "", grp_var = report_params$group_var, show_legend = TRUE) {
  
  # Extract log-CPM efficiently
  lcpm <- get_log_matrix(dge)
  
  # Perform PCA directly on transposed matrix (samples in rows)
  pca_res <- prcomp(t(lcpm), center = TRUE, scale. = FALSE, rank. = 2)
  
  # Get variance explained
  percentVar <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
  
  # Create data frame for plotting
  plot_data <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    Group = factor(dge$samples[[grp_var]]),
    SampleID = rownames(pca_res$x),
    stringsAsFactors = FALSE
  )
  
  # Get unique groups and colors
  groups <- levels(plot_data$Group)
  n_groups <- length(groups)
  colors <- if (n_groups <= 8) {
    RColorBrewer::brewer.pal(max(3, n_groups), "Set1")[1:n_groups]
  } else {
    rainbow(n_groups)
  }
  names(colors) <- groups
  
  # Initialize plotly figure
  fig <- plot_ly()
  
  # Add ellipses and points for each group
  for (grp in groups) {
    grp_data <- plot_data[plot_data$Group == grp, ]
    grp_color <- colors[as.character(grp)]
    
    # Add ellipse (only if >2 points)
    if (nrow(grp_data) > 2) {
      ellipse_coords <- car::dataEllipse(
        grp_data$PC1, grp_data$PC2,
        levels = 0.68, 
        plot.points = FALSE, 
        draw = FALSE
      )
      
      fig <- fig %>%
        add_trace(
          x = ellipse_coords[, 1], 
          y = ellipse_coords[, 2],
          type = "scatter",
          mode = "lines",
          fill = "toself",
          fillcolor = grp_color,
          opacity = 0.3,
          line = list(color = grp_color, dash = "dot", width = 1),
          showlegend = FALSE,
          hoverinfo = "skip",
          name = paste0(grp, "_ellipse")
        )
    }
    
    # Add points
    fig <- fig %>%
      add_trace(
        data = grp_data,
        x = ~PC1, 
        y = ~PC2,
        type = "scatter", 
        mode = "markers",
        name = as.character(grp),
        marker = list(color = grp_color, size = 10),
        text = ~SampleID,
        hoverinfo = "text",
        showlegend = show_legend
      )
  }
  
  # Add layout with variance %
  fig <- fig %>% layout(
    title = title,
    xaxis = list(title = paste0("PC1 (", percentVar[1], "%)")),
    yaxis = list(title = paste0("PC2 (", percentVar[2], "%)"))
  )
  
  return(fig)
}

# Fixed plot_pca_combined - Preserves axis titles
plot_pca_combined <- function(dge_list, grp_var = report_params$group_var, annotation_labels = NULL) {
  
  # Generate individual PCA plots
  n_plots <- length(dge_list)
  show_legend_flags <- c(TRUE, rep(FALSE, n_plots - 1))
  
  pca_plots <- lapply(1:n_plots, function(i) {
    plot_pca(
      dge = dge_list[[i]], 
      title = "", 
      grp_var = grp_var, 
      show_legend = show_legend_flags[i]
    )
  })
  
  # Create subplot with preserved axis titles
  combined_plot <- subplot(
    pca_plots, 
    nrows = 1, 
    shareY = FALSE,
    shareX = FALSE,
    titleX = TRUE,  # KEY: Keep x-axis titles with variance %
    titleY = TRUE   # KEY: Keep y-axis titles
  )
  
  # Add subplot annotations if provided
  if (!is.null(annotation_labels)) {
    annotations <- lapply(seq_along(annotation_labels), function(i) {
      list(
        x = (i - 0.5) / length(annotation_labels),
        y = 1.05,  # Position above plot to avoid overlap
        text = paste0("<b>", annotation_labels[i], "</b>"),
        xref = "paper",
        yref = "paper",
        xanchor = "center",
        yanchor = "bottom",
        showarrow = FALSE,
        font = list(size = 14)
      )
    })
    combined_plot <- combined_plot %>% layout(
      annotations = annotations,
      margin = list(t = 60)  # Add top margin for annotations
    )
  }
  
  return(combined_plot)
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


get_log_matrix <- function(dge) {
  if (!is.null(dge$E_corrected)) {
    return(dge$E_corrected)
  } else if (!is.null(dge$counts)) {
    return(cpm(dge, log = TRUE))
  } else {
    stop("Input DGE object must contain either 'counts' or 'E_corrected'.")
  }
}

