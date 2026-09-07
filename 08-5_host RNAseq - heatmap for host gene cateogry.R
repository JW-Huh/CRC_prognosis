#-----------------------------------------------------------------#
# Batch heatmaps for selected host response pathways
#-----------------------------------------------------------------#


rm(list = ls())
options(stringsAsFactors = FALSE)
setwd("D:/2-연구/2-CRC metagenomics/")


load("host_RNAseq/results_clean_metabolite_host/pathway_gene_heatmaps/pathway_gene_heatmap_inputs.RData")

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("ComplexHeatmap")
}

if (!requireNamespace("circlize", quietly = TRUE)) {
  install.packages("circlize", type = "binary")
}

if (!requireNamespace("svglite", quietly = TRUE)) {
  install.packages("svglite", type = "binary")
}

library(dplyr)
library(tibble)
library(stringr)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(grid)

dir.create("figures", showWarnings = FALSE)
dir.create("host_RNAseq/results_clean_metabolite_host/pathway_gene_heatmaps",
           recursive = TRUE, showWarnings = FALSE)

#-----------------------------------------------------------------#
# Function
#-----------------------------------------------------------------#

plot_pathway_gene_heatmap <- function(
    term_regex,
    output_stub,
    plot_title = NULL,
    gene_heatmap_n = 35,
    gene_p_cutoff = 0.10,
    min_genes_to_plot = 5,
    use_leading_edge = TRUE
) {
  
  pathway_annot <- gsea_plot_df %>%
    dplyr::filter(
      stringr::str_detect(
        display_label,
        stringr::regex(term_regex, ignore_case = TRUE)
      )
    ) %>%
    dplyr::slice_head(n = 1)
  
  if (nrow(pathway_annot) == 0) {
    warning("No pathway matched term_regex: ", term_regex)
    return(NULL)
  }
  
  pathway_id <- pathway_annot$pathway_id[1]
  gs_name <- pathway_annot$gs_name[1]
  display_label <- pathway_annot$display_label[1]
  
  if (is.null(plot_title)) {
    plot_title <- paste0(display_label, " genes")
  }
  
  message("Selected pathway for ", output_stub, ": ", display_label)
  message("pathway_id: ", pathway_id)
  message("gs_name: ", gs_name)
  
  genes_all <- gene_sets_use_long %>%
    dplyr::filter(
      pathway_id == !!pathway_id,
      gs_name == !!gs_name
    ) %>%
    dplyr::pull(Gene) %>%
    unique()
  
  leading_edge <- fgsea_res %>%
    dplyr::filter(
      pathway_id == !!pathway_id,
      gs_name == !!gs_name
    ) %>%
    dplyr::pull(leadingEdge)
  
  if (length(leading_edge) > 0) {
    leading_edge <- leading_edge[[1]]
  } else {
    leading_edge <- character(0)
  }
  
  if (use_leading_edge) {
    genes_for_plot <- intersect(leading_edge, rownames(vst_int))
  } else {
    genes_for_plot <- character(0)
  }
  
  if (length(genes_for_plot) < min_genes_to_plot) {
    genes_for_plot <- intersect(genes_all, rownames(vst_int))
  }
  
  if (length(genes_for_plot) < min_genes_to_plot) {
    warning("Too few genes found for: ", display_label)
    return(NULL)
  }
  
  sample_annot <- col_int[colnames(vst_int), , drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::mutate(
      TRG_plot = factor(TRG_plot, levels = c("non_pCR", "pCR"))
    )
  
  #----------------------------#
  # Gene-level statistics
  #----------------------------#
  
  gene_stat_all <- purrr::map_dfr(
    genes_for_plot,
    function(gene) {
      x <- tibble::tibble(
        Sample = colnames(vst_int),
        expr = as.numeric(vst_int[gene, colnames(vst_int)]),
        Response = sample_annot$TRG_plot[
          match(colnames(vst_int), sample_annot$Sample)
        ]
      ) %>%
        dplyr::filter(!is.na(expr), !is.na(Response))
      
      wt <- tryCatch(
        stats::wilcox.test(expr ~ Response, data = x, exact = FALSE),
        error = function(e) NULL
      )
      
      tibble::tibble(
        Gene = gene,
        median_non_pCR = median(x$expr[x$Response == "non_pCR"], na.rm = TRUE),
        median_pCR = median(x$expr[x$Response == "pCR"], na.rm = TRUE),
        logFC_non_pCR_minus_pCR =
          median(x$expr[x$Response == "non_pCR"], na.rm = TRUE) -
          median(x$expr[x$Response == "pCR"], na.rm = TRUE),
        pval_wilcox_gene = ifelse(is.null(wt), NA_real_, wt$p.value)
      )
    }
  ) %>%
    dplyr::mutate(
      FDR_wilcox_gene = p.adjust(pval_wilcox_gene, method = "BH"),
      neglog10_wilcox_p = -log10(pval_wilcox_gene + 1e-300),
      neglog10_wilcox_p_capped = pmin(neglog10_wilcox_p, 5),
      p_symbol = dplyr::case_when(
        pval_wilcox_gene < 0.001 ~ "***",
        pval_wilcox_gene < 0.01 ~ "**",
        pval_wilcox_gene < 0.05 ~ "*",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::arrange(
      pval_wilcox_gene,
      dplyr::desc(abs(logFC_non_pCR_minus_pCR))
    )
  
  readr::write_csv(
    gene_stat_all,
    paste0(
      "host_RNAseq/results_clean_metabolite_host/pathway_gene_heatmaps/",
      output_stub,
      "_gene_logFC_wilcoxon_all.csv"
    )
  )
  
  gene_stat_tile <- gene_stat_all %>%
    dplyr::filter(
      !is.na(pval_wilcox_gene),
      pval_wilcox_gene < gene_p_cutoff
    ) %>%
    dplyr::slice_head(n = gene_heatmap_n)
  
  if (nrow(gene_stat_tile) < min_genes_to_plot) {
    message(
      "Fewer than ", min_genes_to_plot,
      " genes passed P < ", gene_p_cutoff,
      ". Falling back to top genes by absolute logFC."
    )
    
    gene_stat_tile <- gene_stat_all %>%
      dplyr::filter(!is.na(logFC_non_pCR_minus_pCR)) %>%
      dplyr::arrange(
        dplyr::desc(abs(logFC_non_pCR_minus_pCR))
      ) %>%
      dplyr::slice_head(n = gene_heatmap_n)
  }
  
  readr::write_csv(
    gene_stat_tile,
    paste0(
      "host_RNAseq/results_clean_metabolite_host/pathway_gene_heatmaps/",
      output_stub,
      "_gene_logFC_wilcoxon_plotted.csv"
    )
  )
  
  plot_genes <- gene_stat_tile %>%
    dplyr::pull(Gene)
  
  expr_z <- vst_int[plot_genes, , drop = FALSE]
  expr_z <- t(scale(t(expr_z)))
  expr_z[expr_z > 2.5] <- 2.5
  expr_z[expr_z < -2.5] <- -2.5
  
  gene_stat_tile_for_heatmap <- gene_stat_tile %>%
    dplyr::filter(Gene %in% rownames(expr_z)) %>%
    dplyr::select(
      Gene,
      logFC_non_pCR_minus_pCR,
      neglog10_wilcox_p_capped,
      pval_wilcox_gene,
      FDR_wilcox_gene,
      p_symbol
    ) %>%
    tibble::column_to_rownames("Gene")
  
  gene_stat_tile_for_heatmap <- gene_stat_tile_for_heatmap[
    rownames(expr_z),
    ,
    drop = FALSE
  ]
  
  #----------------------------#
  # Sample annotation
  #----------------------------#
  
  col_annot <- sample_annot %>%
    dplyr::select(Sample, Response = TRG_plot) %>%
    tibble::column_to_rownames("Sample")
  
  col_annot <- col_annot[colnames(expr_z), , drop = FALSE]
  
  #----------------------------#
  # Clustering
  #----------------------------#
  
  row_cor <- stats::cor(
    t(expr_z),
    use = "pairwise.complete.obs"
  )
  row_cor[is.na(row_cor)] <- 0
  diag(row_cor) <- 1
  
  col_cor <- stats::cor(
    expr_z,
    use = "pairwise.complete.obs"
  )
  col_cor[is.na(col_cor)] <- 0
  diag(col_cor) <- 1
  
  row_hc <- stats::hclust(
    stats::as.dist(1 - row_cor),
    method = "complete"
  )
  
  col_hc <- stats::hclust(
    stats::as.dist(1 - col_cor),
    method = "complete"
  )
  
  row_reorder_weight <- gene_stat_tile_for_heatmap[
    row_hc$labels,
    "logFC_non_pCR_minus_pCR"
  ]
  
  row_dend <- stats::reorder(
    stats::as.dendrogram(row_hc),
    wts = row_reorder_weight,
    agglo.FUN = mean
  )
  
  sample_score <- colMeans(expr_z, na.rm = TRUE)
  
  col_reorder_weight <- col_annot %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::mutate(
      Response_weight = dplyr::case_when(
        Response == "pCR" ~ 0,
        Response == "non_pCR" ~ 100,
        TRUE ~ 50
      ),
      Pathway_score = sample_score[Sample],
      reorder_weight =
        Response_weight +
        rank(Pathway_score, ties.method = "average") / 1000
    ) %>%
    dplyr::select(Sample, reorder_weight) %>%
    tibble::deframe()
  
  col_dend <- stats::reorder(
    stats::as.dendrogram(col_hc),
    wts = col_reorder_weight[col_hc$labels],
    agglo.FUN = mean
  )
  
  #----------------------------#
  # Colors
  #----------------------------#
  
  expr_col_fun <- circlize::colorRamp2(
    c(-2.5, 0, 2.5),
    c("#3B6EA8", "white", "#B14A4A")
  )
  
  logFC_cap <- max(
    0.5,
    stats::quantile(
      abs(gene_stat_tile_for_heatmap$logFC_non_pCR_minus_pCR),
      probs = 0.95,
      na.rm = TRUE
    )
  )
  
  logFC_col_fun <- circlize::colorRamp2(
    c(-logFC_cap, 0, logFC_cap),
    c(group_cols[["pCR"]], "white", group_cols[["non_pCR"]])
  )
  
  p_col_fun <- circlize::colorRamp2(
    c(1.3, 2, 5),
    c("#E8EDF8", "#8EA6D9", "#08519C")
  )
  
  p_text_col <- dplyr::case_when(
    gene_stat_tile_for_heatmap$neglog10_wilcox_p_capped >= 2 ~ "white",
    TRUE ~ "grey20"
  )
  
  ha_col <- ComplexHeatmap::HeatmapAnnotation(
    Response = col_annot$Response,
    col = list(
      Response = c(
        "non_pCR" = group_cols[["non_pCR"]],
        "pCR" = group_cols[["pCR"]]
      )
    ),
    annotation_name_gp = grid::gpar(fontsize = 8),
    simple_anno_size = grid::unit(3.2, "mm")
  )
  
  ha_row_left <- ComplexHeatmap::rowAnnotation(
    Gene = ComplexHeatmap::anno_text(
      rownames(gene_stat_tile_for_heatmap),
      just = "right",
      location = grid::unit(1, "npc"),
      gp = grid::gpar(fontsize = 6.4)
    ),
    `logFC` = gene_stat_tile_for_heatmap$logFC_non_pCR_minus_pCR,
    `P` = ComplexHeatmap::anno_simple(
      gene_stat_tile_for_heatmap$neglog10_wilcox_p_capped,
      col = p_col_fun,
      pch = gene_stat_tile_for_heatmap$p_symbol,
      pt_gp = grid::gpar(
        col = p_text_col,
        fontsize = 6.2,
        fontface = "bold"
      ),
      pt_size = grid::unit(2.5, "mm")
    ),
    col = list(
      `logFC` = logFC_col_fun
    ),
    annotation_name_side = "top",
    annotation_name_rot = 0,
    annotation_name_gp = grid::gpar(fontsize = 7),
    simple_anno_size = grid::unit(4.2, "mm"),
    annotation_width = grid::unit.c(
      grid::unit(38, "mm"),
      grid::unit(4.2, "mm"),
      grid::unit(4.2, "mm")
    )
  )
  
  ht <- ComplexHeatmap::Heatmap(
    expr_z,
    name = "z",
    col = expr_col_fun,
    cluster_rows = row_dend,
    cluster_columns = col_dend,
    left_annotation = ha_row_left,
    top_annotation = ha_col,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_dend_side = "right",
    column_dend_side = "top",
    column_title = plot_title,
    column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Row z-score",
      at = c(-2, -1, 0, 1, 2)
    ),
    use_raster = FALSE
  )
  
  #----------------------------#
  # Draw and save
  #----------------------------#
  
  grid::grid.newpage()
  
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = FALSE
  )
  
  svglite::svglite(
    paste0(
      "figures/host_RNAseq_",
      output_stub,
      "_gene_heatmap_clustered_left_gene_logFC_wilcox.svg"
    ),
    width = 7,
    height = max(5.2, 0.18 * nrow(expr_z) + 1.2)
  )
  
  grid::grid.newpage()
  
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = FALSE
  )
  
  grDevices::dev.off()
  
  return(
    list(
      pathway = pathway_annot,
      genes_all = genes_all,
      genes_used = plot_genes,
      gene_stat_all = gene_stat_all,
      gene_stat_tile = gene_stat_tile,
      heatmap = ht
    )
  )
}

#-----------------------------------------------------------------#
# Check matching pathway names before plotting
#-----------------------------------------------------------------#

pathway_terms <- tibble::tribble(
  ~output_stub, ~term_regex, ~plot_title,
  "emt", "epithelial.*mesenchymal|\\bEMT\\b", "Epithelial–mesenchymal transition genes",
  "ros_response", "oxidative.*stress|ROS", "Oxidative stress / ROS response genes",
  "nfkb_il6_jak_stat", "TNF.*NF|NF.*kB|IL-?6|JAK.*STAT", "TNF–NFκB / IL-6–JAK–STAT signaling genes",
  "g2m_checkpoint", "G2.?M", "G2/M checkpoint genes",
  "nad_niacin", "NAD|niacin", "NAD / niacin metabolism genes",
  "mucus_goblet", "mucus|goblet", "Mucus / goblet-cell program genes",
  "tcell_activation", "T[- ]cell activation", "T-cell activation genes"
)

pathway_match_check <- purrr::map_dfr(
  seq_len(nrow(pathway_terms)),
  function(i) {
    gsea_plot_df %>%
      dplyr::filter(
        stringr::str_detect(
          display_label,
          stringr::regex(pathway_terms$term_regex[i], ignore_case = TRUE)
        )
      ) %>%
      dplyr::mutate(
        output_stub = pathway_terms$output_stub[i],
        query_regex = pathway_terms$term_regex[i]
      ) %>%
      dplyr::select(
        output_stub,
        query_regex,
        display_label,
        pathway_id,
        gs_name,
        NES,
        padj
      )
  }
)

print(pathway_match_check)

readr::write_csv(
  pathway_match_check,
  "host_RNAseq/results_clean_metabolite_host/pathway_gene_heatmaps/pathway_match_check.csv"
)

#-----------------------------------------------------------------#
# Run batch heatmaps
#-----------------------------------------------------------------#

heatmap_results <- purrr::map(
  seq_len(nrow(pathway_terms)),
  function(i) {
    plot_pathway_gene_heatmap(
      term_regex = pathway_terms$term_regex[i],
      output_stub = pathway_terms$output_stub[i],
      plot_title = pathway_terms$plot_title[i],
      gene_heatmap_n = 35,
      gene_p_cutoff = 0.10,
      min_genes_to_plot = 5,
      use_leading_edge = TRUE
    )
  }
)

names(heatmap_results) <- pathway_terms$output_stub

save(
  heatmap_results,
  pathway_match_check,
  file = "host_RNAseq/results_clean_metabolite_host/pathway_gene_heatmaps/pathway_gene_heatmap_batch_results.RData"
)