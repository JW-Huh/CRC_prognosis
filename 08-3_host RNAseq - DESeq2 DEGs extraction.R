#-----------------------------------------------------------------#

#

# Host RNA-seq differential expression analysis

#

# Baseline tumor samples, pCR vs non-pCR

#

#-----------------------------------------------------------------#

rm(list = ls())
options(stringsAsFactors = FALSE)
setwd("D:/2-연구/2-CRC metagenomics/")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggbeeswarm)
  library(DESeq2)
  library(SummarizedExperiment)
})

dir.create("host_RNAseq/results_clean_metabolite_host", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

#-----------------------------------------------------------------#

# 1. Load integrated input object

#-----------------------------------------------------------------#

load("host_RNAseq/results_clean_metabolite_host/host_met_rna_matched_input.RData")

stopifnot(exists("cnt_rna"))
stopifnot(exists("col_rna"))
stopifnot(exists("cnt_int"))
stopifnot(exists("col_int"))
stopifnot(exists("met_int"))

col_rna <- as.data.frame(col_rna)
col_rna$TRG_plot <- factor(col_rna$TRG_plot, levels = c("non_pCR", "pCR"))
rownames(col_rna) <- col_rna$RNA_sample_id

cnt_rna <- cnt_rna[, rownames(col_rna), drop = FALSE]

stopifnot(identical(colnames(cnt_rna), rownames(col_rna)))
stopifnot(all(col_rna$TRG_plot %in% c("non_pCR", "pCR")))

message("RNA count matrix: ", nrow(cnt_rna), " genes x ", ncol(cnt_rna), " samples")

print(
  col_rna %>%
    dplyr::count(TRG_plot, name = "n")
)

#-----------------------------------------------------------------#

# 2. DESeq2 differential expression analysis

#-----------------------------------------------------------------#

dds_rna <- DESeq2::DESeqDataSetFromMatrix(
  countData = cnt_rna,
  colData = col_rna,
  design = ~ TRG_plot
)

dds_rna <- dds_rna[rowSums(DESeq2::counts(dds_rna) >= 10) >= 3, ]

message("Genes retained after low-count filtering: ", nrow(dds_rna))

dds_rna <- DESeq2::DESeq(dds_rna)

deg_rna <- DESeq2::results(
  dds_rna,
  contrast = c("TRG_plot", "pCR", "non_pCR"),
  alpha = 0.10
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  tibble::as_tibble() %>%
  dplyr::rename(
    pval_DESeq2 = pvalue,
    FDR_DESeq2 = padj,
    stat_DESeq2 = stat
  ) %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      log2FoldChange > 0 ~ "pCR_high",
      log2FoldChange < 0 ~ "non_pCR_high",
      TRUE ~ "no_difference"
    ),
    DEG_sig = dplyr::case_when(
      !is.na(FDR_DESeq2) & FDR_DESeq2 < 0.10 ~ "DESeq2_FDR<0.1",
      !is.na(pval_DESeq2) & pval_DESeq2 < 0.05 ~ "DESeq2_nominal_P<0.05",
      TRUE ~ "not_significant"
    ),
    neglog10_pval_DESeq2 = -log10(pval_DESeq2 + 1e-300),
    neglog10_FDR_DESeq2 = -log10(FDR_DESeq2 + 1e-300),
    rank_stat_DESeq2 = stat_DESeq2,
    rank_signed_logp_DESeq2 = sign(log2FoldChange) * neglog10_pval_DESeq2
  ) %>%
  dplyr::arrange(pval_DESeq2)

#-----------------------------------------------------------------#

# 3. VST-normalized expression matrix

#-----------------------------------------------------------------#

vst_rna <- SummarizedExperiment::assay(DESeq2::vst(dds_rna, blind = FALSE))

stopifnot(identical(colnames(vst_rna), rownames(col_rna)))

message("VST matrix: ", nrow(vst_rna), " genes x ", ncol(vst_rna), " samples")

#-----------------------------------------------------------------#

# 4. Wilcoxon rank-sum test using VST-normalized expression

#-----------------------------------------------------------------#

pval_wilcox <- apply(
  vst_rna,
  1,
  function(x) {
    suppressWarnings(
      wilcox.test(
        x[col_rna$TRG_plot == "pCR"],
        x[col_rna$TRG_plot == "non_pCR"],
        exact = FALSE
      )$p.value
    )
  }
)

deg_rna <- deg_rna %>%
  dplyr::left_join(
    tibble::tibble(
      Gene = names(pval_wilcox),
      pval_wilcox = as.numeric(pval_wilcox),
      FDR_wilcox = p.adjust(as.numeric(pval_wilcox), method = "BH"),
      neglog10_pval_wilcox = -log10(as.numeric(pval_wilcox) + 1e-300),
      neglog10_FDR_wilcox = -log10(p.adjust(as.numeric(pval_wilcox), method = "BH") + 1e-300)
    ),
    by = "Gene"
  ) %>%
  dplyr::mutate(
    DEG_sig_combined = dplyr::case_when(
      !is.na(FDR_DESeq2) & FDR_DESeq2 < 0.10 &
        !is.na(FDR_wilcox) & FDR_wilcox < 0.10 ~ "DESeq2_FDR<0.1_and_Wilcoxon_FDR<0.1",
      !is.na(FDR_DESeq2) & FDR_DESeq2 < 0.10 ~ "DESeq2_FDR<0.1_only",
      !is.na(FDR_wilcox) & FDR_wilcox < 0.10 ~ "Wilcoxon_FDR<0.1_only",
      !is.na(pval_DESeq2) & pval_DESeq2 < 0.05 |
        !is.na(pval_wilcox) & pval_wilcox < 0.05 ~ "nominal_P<0.05_only",
      TRUE ~ "not_significant"
    )
  ) %>%
  dplyr::arrange(pval_DESeq2)

#-----------------------------------------------------------------#

# 5. Prepare matched VST matrix for later metabolite correlation

#-----------------------------------------------------------------#

col_int <- as.data.frame(col_int)
col_int$TRG_plot <- factor(col_int$TRG_plot, levels = c("non_pCR", "pCR"))
rownames(col_int) <- col_int$RNA_sample_id

vst_int <- vst_rna[, rownames(col_int), drop = FALSE]

stopifnot(identical(colnames(vst_int), rownames(col_int)))
stopifnot(identical(rownames(met_int), colnames(vst_int)))

message("Matched VST matrix for metabolite correlation: ", nrow(vst_int), " genes x ", ncol(vst_int), " samples")
message("Matched metabolite matrix: ", nrow(met_int), " samples x ", ncol(met_int), " metabolites")




#-----------------------------------------------------------------#

# 6. Plot color setting

#-----------------------------------------------------------------#

group_cols <- c(
  pCR = "#4FAE9A",
  non_pCR = "#DE7872"
)

deg_cols <- c(
  "Not significant" = "#D9D9D9",
  "pCR high, nominal P<0.05" = "#B7DED6",
  "pCR high, FDR<0.1" = unname(group_cols["pCR"]),
  "non-pCR high, nominal P<0.05" = "#EDBBB8",
  "non-pCR high, FDR<0.1" = unname(group_cols["non_pCR"])
)

deg_rna <- deg_rna %>%
  dplyr::mutate(
    DEG_plot_sig = dplyr::case_when(
      !is.na(FDR_DESeq2) & FDR_DESeq2 < 0.10 & log2FoldChange > 0 ~ "pCR high, FDR<0.1",
      !is.na(FDR_DESeq2) & FDR_DESeq2 < 0.10 & log2FoldChange < 0 ~ "non-pCR high, FDR<0.1",
      !is.na(pval_DESeq2) & pval_DESeq2 < 0.05 & log2FoldChange > 0 ~ "pCR high, nominal P<0.05",
      !is.na(pval_DESeq2) & pval_DESeq2 < 0.05 & log2FoldChange < 0 ~ "non-pCR high, nominal P<0.05",
      TRUE ~ "Not significant"
    ),
    DEG_plot_sig = factor(
      DEG_plot_sig,
      levels = c(
        "Not significant",
        "pCR high, nominal P<0.05",
        "pCR high, FDR<0.1",
        "non-pCR high, nominal P<0.05",
        "non-pCR high, FDR<0.1"
      )
    )
  )

#-----------------------------------------------------------------#

# 7. MA plot

#-----------------------------------------------------------------#

library(ggrepel)

label_ma <- deg_rna %>%
  dplyr::filter(!is.na(FDR_DESeq2), FDR_DESeq2 < 0.10, !is.na(baseMean), !is.na(log2FoldChange)) %>%
  dplyr::arrange(FDR_DESeq2, dplyr::desc(abs(log2FoldChange))) %>%
  dplyr::slice_head(n = 20)


p_MA <- deg_rna %>%
  dplyr::filter(!is.na(baseMean), !is.na(log2FoldChange)) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = baseMean + 1,
      y = log2FoldChange,
      color = DEG_plot_sig
    )
  ) +
  ggplot2::geom_point(size = 1.1, alpha = 0.75) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  ggrepel::geom_text_repel(
    data = label_ma,
    ggplot2::aes(label = Gene),
    size = 3,
    max.overlaps = 100,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey50",
    show.legend = FALSE
  ) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_color_manual(values = deg_cols, drop = FALSE) +
  ggplot2::labs(
    x = "Mean normalized count + 1",
    y = "log2 fold change (pCR / non-pCR)",
    color = "DESeq2 significance",
    title = "MA plot: pCR vs non-pCR"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank()
  )

p_MA

ggplot2::ggsave(
  "figures/host_RNAseq_MA_plot_pCR_vs_non_pCR.svg",
  p_MA,
  width = 7,
  height = 5.5,
  device = "svg"
)


#-----------------------------------------------------------------#

# 8. Volcano plot

#-----------------------------------------------------------------#

label_volcano <- deg_rna %>%
  dplyr::filter(!is.na(FDR_DESeq2), FDR_DESeq2 < 0.10, !is.na(log2FoldChange), !is.na(pval_DESeq2)) %>%
  dplyr::arrange(FDR_DESeq2, dplyr::desc(abs(log2FoldChange))) %>%
  dplyr::slice_head(n = 20)

p_volcano <- deg_rna %>%
  dplyr::filter(!is.na(log2FoldChange), !is.na(pval_DESeq2)) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = log2FoldChange,
      y = neglog10_pval_DESeq2,
      color = DEG_plot_sig
    )
  ) +
  ggplot2::geom_point(size = 1.1, alpha = 0.75) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
  ggrepel::geom_text_repel(
    data = label_volcano,
    ggplot2::aes(label = Gene),
    size = 3,
    max.overlaps = 100,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey50",
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(values = deg_cols, drop = FALSE) +
  ggplot2::labs(
    x = "log2 fold change (pCR / non-pCR)",
    y = "-log10(DESeq2 p-value)",
    color = "DESeq2 significance",
    title = "Volcano plot: pCR vs non-pCR"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank()
  )

p_volcano

ggplot2::ggsave(
  "figures/host_RNAseq_volcano_pCR_vs_non_pCR.svg",
  p_volcano,
  width = 7,
  height = 5.5,
  device = "svg"
)


#-----------------------------------------------------------------#

# 9. Boxplot + beeswarm for significant genes

#-----------------------------------------------------------------#

top_boxplot_n <- 20

sig_genes <- deg_rna %>%
  dplyr::filter(
    (!is.na(FDR_DESeq2) & FDR_DESeq2 < 0.10) |
      (!is.na(pval_DESeq2) & pval_DESeq2 < 0.05)
  ) %>%
  dplyr::arrange(
    dplyr::if_else(is.na(FDR_DESeq2), Inf, FDR_DESeq2),
    pval_DESeq2
  ) %>%
  dplyr::slice_head(n = top_boxplot_n) %>%
  dplyr::pull(Gene)

if (length(sig_genes) > 0) {
  
  box_df <- vst_rna[sig_genes, , drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Gene") %>%
    tidyr::pivot_longer(
      cols = -Gene,
      names_to = "RNA_sample_id",
      values_to = "vst_expr"
    ) %>%
    dplyr::left_join(
      col_rna %>%
        dplyr::select(RNA_sample_id, TRG_plot),
      by = "RNA_sample_id"
    ) %>%
    dplyr::left_join(
      deg_rna %>%
        dplyr::select(Gene, log2FoldChange, pval_DESeq2, FDR_DESeq2, DEG_plot_sig),
      by = "Gene"
    ) %>%
    dplyr::mutate(
      Gene = factor(Gene, levels = rev(sig_genes)),
      TRG_plot = factor(TRG_plot, levels = c("non_pCR", "pCR"))
    )
  
  p_box_bee <- box_df %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = TRG_plot,
        y = vst_expr,
        color = TRG_plot,
        fill = TRG_plot
      )
    ) +
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      width = 0.55,
      alpha = 0.25,
      linewidth = 0.35
    ) +
    ggbeeswarm::geom_quasirandom(
      width = 0.15,
      size = 1.7,
      alpha = 0.85
    ) +
    ggplot2::facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
    ggplot2::scale_color_manual(values = group_cols, drop = FALSE) +
    ggplot2::scale_fill_manual(values = group_cols, drop = FALSE) +
    ggplot2::labs(
      x = NULL,
      y = "VST-normalized expression",
      color = "Response",
      fill = "Response",
      title = "Top significant genes: pCR vs non-pCR"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey95"),
      strip.text = ggplot2::element_text(size = 9)
    )
  
  p_box_bee
  
  ggplot2::ggsave(
    "figures/host_RNAseq_top_significant_genes_boxplot_beeswarm.svg",
    p_box_bee,
    width = 10,
    height = max(5, 2.2 * ceiling(length(sig_genes) / 4)),
    device = "svg"
  )
  
} else {
  message("No DESeq2 FDR<0.1 or nominal P<0.05 genes found. Boxplot + beeswarm figure was skipped.")
}

#-----------------------------------------------------------------#

# 10. DESeq2 p-value vs Wilcoxon p-value scatter plot

#-----------------------------------------------------------------#

label_pval_scatter <- deg_rna %>%
  dplyr::filter(
    !is.na(pval_DESeq2),
    !is.na(pval_wilcox),
    (
      (!is.na(FDR_DESeq2) & FDR_DESeq2 < 0.10) |
        (!is.na(FDR_wilcox) & FDR_wilcox < 0.10)
    )
  ) %>%
  dplyr::arrange(
    dplyr::if_else(is.na(FDR_DESeq2), Inf, FDR_DESeq2),
    dplyr::if_else(is.na(FDR_wilcox), Inf, FDR_wilcox),
    pval_DESeq2
  ) %>%
  dplyr::slice_head(n = 25)

p_pval_scatter <- deg_rna %>%
  dplyr::filter(
    !is.na(pval_DESeq2),
    !is.na(pval_wilcox)
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = neglog10_pval_DESeq2,
      y = neglog10_pval_wilcox,
      color = DEG_plot_sig
    )
  ) +
  ggplot2::geom_point(size = 1.1, alpha = 0.75) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.3
  ) +
  ggplot2::geom_hline(
    yintercept = -log10(0.05),
    linetype = "dotted",
    linewidth = 0.3
  ) +
  ggplot2::geom_vline(
    xintercept = -log10(0.05),
    linetype = "dotted",
    linewidth = 0.3
  ) +
  ggrepel::geom_text_repel(
    data = label_pval_scatter,
    ggplot2::aes(label = Gene),
    size = 3,
    max.overlaps = 100,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.color = "grey50",
    segment.linewidth = 0.25,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(values = deg_cols, drop = FALSE) +
  ggplot2::labs(
    x = "-log10(DESeq2 p-value)",
    y = "-log10(Wilcoxon p-value)",
    color = "DESeq2 significance",
    title = "DESeq2 vs Wilcoxon evidence"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank()
  )

p_pval_scatter

ggplot2::ggsave(
  "figures/host_RNAseq_DESeq2_vs_Wilcoxon_pvalue_scatter.svg",
  p_pval_scatter,
  width = 7,
  height = 5.5,
  device = "svg"
)



#-----------------------------------------------------------------#

# 11. Save result tables and objects

#-----------------------------------------------------------------#

write.csv(
  deg_rna,
  "host_RNAseq/results_clean_metabolite_host/deg_rna_pCR_vs_non_pCR_with_wilcoxon.csv",
  row.names = FALSE
)

save(
  dds_rna,
  deg_rna,
  vst_rna,
  vst_int,
  col_rna,
  col_int,
  met_int,
  met_vars,

  file = "host_RNAseq/results_clean_metabolite_host/host_rna_DESeq2_pCR_vs_non_pCR.RData"
)

#-----------------------------------------------------------------#

# 12. Quick result check

#-----------------------------------------------------------------#

print(
  deg_rna %>%
    dplyr::count(DEG_sig, direction, name = "n")
)

print(
  deg_rna %>%
    dplyr::count(DEG_sig_combined, direction, name = "n")
)

print(
  deg_rna %>%
    dplyr::select(
      Gene,
      baseMean,
      log2FoldChange,
      pval_DESeq2,
      FDR_DESeq2,
      pval_wilcox,
      FDR_wilcox,
      direction,
      DEG_sig,
      DEG_sig_combined
    ) %>%
    head(30)
)
