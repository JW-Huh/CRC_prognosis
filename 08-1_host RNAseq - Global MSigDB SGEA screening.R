#-----------------------------------------------------------------#
# 17-2d. Global MSigDB GSEA screening
#-----------------------------------------------------------------#

setwd("D:/2-연구/2-CRC metagenomics/")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(readr)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(forcats)
})

dir.create(
  "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/global_gsea",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create("figures", showWarnings = FALSE)

load("host_RNAseq/results_clean_metabolite_host/pathway_metabolite/global_gsea_inputs.RData")

#-----------------------------------------------------------------#
# 1. Prepare MSigDB table
#-----------------------------------------------------------------#

if (!exists("msig")) {
  msig <- msigdbr::msigdbr(species = "Homo sapiens")
}

if (!"gs_collection" %in% colnames(msig) && "gs_cat" %in% colnames(msig)) {
  msig <- msig %>% dplyr::rename(gs_collection = gs_cat)
}

if (!"gs_subcollection" %in% colnames(msig) && "gs_subcat" %in% colnames(msig)) {
  msig <- msig %>% dplyr::rename(gs_subcollection = gs_subcat)
}

if (!"gs_description" %in% colnames(msig)) {
  msig$gs_description <- NA_character_
}

if (!"gs_exact_source" %in% colnames(msig)) {
  msig$gs_exact_source <- NA_character_
}

colnames(msig)
# [1] "gs_name"          "Gene"             "gs_collection"    "gs_subcollection"
# [5] "gs_description"   "gs_exact_source"  "db_version" 

msig_all_long <- msig %>%
  dplyr::select(
    gs_name,
    Gene,
    gs_collection,
    gs_subcollection,
    gs_description,
    gs_exact_source
  ) %>%
  dplyr::filter(
    !is.na(Gene),
    Gene != "",
    Gene %in% rownames(vst_rna)
  ) %>%
  dplyr::distinct()

#-----------------------------------------------------------------#
# 2. Choose global screening collections
#-----------------------------------------------------------------#

msig_global_long <- msig_all_long %>%
  dplyr::filter(
    gs_collection == "H" |
      (
        gs_collection == "C2" &
          stringr::str_detect(gs_subcollection, "^CP")
      ) |
      (
        gs_collection == "C5" &
          gs_subcollection %in% c("GO:BP", "GO:MF", "GO:CC")
      )
  )

gene_set_size_global <- msig_global_long %>%
  dplyr::group_by(
    gs_name,
    gs_collection,
    gs_subcollection,
    gs_description,
    gs_exact_source
  ) %>%
  dplyr::summarise(
    n_genes_total = dplyr::n_distinct(Gene),
    .groups = "drop"
  )

min_gs_size_global <- 10
max_gs_size_global <- 500

use_gs_global <- gene_set_size_global %>%
  dplyr::filter(
    n_genes_total >= min_gs_size_global,
    n_genes_total <= max_gs_size_global
  ) %>%
  dplyr::pull(gs_name)

msig_global_use_long <- msig_global_long %>%
  dplyr::filter(gs_name %in% use_gs_global)

gene_sets_global <- split(
  msig_global_use_long$Gene,
  msig_global_use_long$gs_name
)

gene_sets_global <- lapply(gene_sets_global, unique)

message("Global gene sets used: ", length(gene_sets_global))

readr::write_csv(
  gene_set_size_global,
  "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/global_gsea/global_MSigDB_gene_set_size_table.csv"
)

#-----------------------------------------------------------------#
# 3. Prepare preranked DESeq2 statistic
#-----------------------------------------------------------------#

rank_tbl <- deg_rna %>%
  dplyr::filter(
    !is.na(Gene),
    Gene %in% rownames(vst_rna),
    !is.na(stat_DESeq2)
  ) %>%
  dplyr::group_by(Gene) %>%
  dplyr::arrange(dplyr::desc(abs(stat_DESeq2)), .by_group = TRUE) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(dplyr::desc(stat_DESeq2))

rank_stat <- rank_tbl$stat_DESeq2
names(rank_stat) <- rank_tbl$Gene
rank_stat <- sort(rank_stat, decreasing = TRUE)
rank_stat <- rank_stat[is.finite(rank_stat)]

#-----------------------------------------------------------------#
# 4. Run global preranked GSEA
#-----------------------------------------------------------------#

fgsea_global_raw <- fgsea::fgsea(
  pathways = gene_sets_global,
  stats = rank_stat,
  minSize = min_gs_size_global,
  maxSize = max_gs_size_global,
  eps = 0
)

fgsea_global_res <- fgsea_global_raw %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    gs_name = pathway,
    leadingEdge_string = vapply(
      leadingEdge,
      paste,
      collapse = ";",
      FUN.VALUE = character(1)
    ),
    Direction = dplyr::case_when(
      NES > 0 ~ "pCR_high",
      NES < 0 ~ "non_pCR_high",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::left_join(
    gene_set_size_global,
    by = "gs_name",
    relationship = "many-to-one"
  ) %>%
  dplyr::arrange(padj, pval, dplyr::desc(abs(NES)))

readr::write_csv(
  fgsea_global_res %>% dplyr::select(-leadingEdge),
  "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/global_gsea/GSEA_global_MSigDB_H_C2CP_C5GO.csv"
)

#-----------------------------------------------------------------#
# 5. Top term tables
#-----------------------------------------------------------------#

fgsea_global_top_pCR <- fgsea_global_res %>%
  dplyr::filter(
    Direction == "pCR_high",
    !is.na(padj)
  ) %>%
  dplyr::arrange(padj, dplyr::desc(NES)) %>%
  dplyr::slice_head(n = 100)

fgsea_global_top_non_pCR <- fgsea_global_res %>%
  dplyr::filter(
    Direction == "non_pCR_high",
    !is.na(padj)
  ) %>%
  dplyr::arrange(padj, NES) %>%
  dplyr::slice_head(n = 100)

readr::write_csv(
  fgsea_global_top_pCR %>% dplyr::select(-leadingEdge),
  "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/global_gsea/GSEA_global_top100_pCR_high.csv"
)

readr::write_csv(
  fgsea_global_top_non_pCR %>% dplyr::select(-leadingEdge),
  "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/global_gsea/GSEA_global_top100_non_pCR_high.csv"
)

#-----------------------------------------------------------------#
# 6. Simple dot plot: top global terms by direction
#-----------------------------------------------------------------#

format_gsea_label <- function(x) {
  x %>%
    stringr::str_remove("^HALLMARK_") %>%
    stringr::str_remove("^REACTOME_") %>%
    stringr::str_remove("^GOBP_") %>%
    stringr::str_remove("^GOCC_") %>%
    stringr::str_remove("^GOMF_") %>%
    stringr::str_remove("^WP_") %>%
    stringr::str_remove("^KEGG_") %>%
    stringr::str_remove("_WP[0-9]+$") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_lower() %>%
    stringr::str_to_sentence() %>%
    stringr::str_replace_all("\\bCrc\\b", "CRC") %>%
    stringr::str_replace_all("\\bAhr\\b", "AhR") %>%
    stringr::str_replace_all("\\bNfkb\\b", "NF-kB") %>%
    stringr::str_replace_all("\\bTnf\\b", "TNF") %>%
    stringr::str_replace_all("\\bIfn\\b", "IFN") %>%
    stringr::str_replace_all("\\bMyc\\b", "MYC") %>%
    stringr::str_replace_all("\\bEmt\\b", "EMT") %>%
    stringr::str_replace_all("\\bG2m\\b", "G2/M") %>%
    stringr::str_squish()
}

global_dot_df <- dplyr::bind_rows(
  fgsea_global_res %>%
    dplyr::filter(Direction == "pCR_high", !is.na(padj)) %>%
    dplyr::arrange(padj, dplyr::desc(NES)) %>%
    dplyr::slice_head(n = 20),
  fgsea_global_res %>%
    dplyr::filter(Direction == "non_pCR_high", !is.na(padj)) %>%
    dplyr::arrange(padj, NES) %>%
    dplyr::slice_head(n = 20)
) %>%
  dplyr::mutate(
    display_label = format_gsea_label(gs_name),
    neglog10_FDR = -log10(padj + 1e-300),
    neglog10_FDR_capped = pmin(neglog10_FDR, 10),
    display_label = factor(display_label, levels = rev(unique(display_label)))
  )

p_global_gsea <- global_dot_df %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = NES,
      y = display_label,
      color = Direction,
      size = neglog10_FDR_capped
    )
  ) +
  ggplot2::geom_vline(xintercept = 0, color = "grey65", linewidth = 0.3) +
  ggplot2::geom_point(alpha = 0.9) +
  ggplot2::scale_color_manual(
    values = c(
      "pCR_high" = group_cols[["pCR"]],
      "non_pCR_high" = group_cols[["non_pCR"]]
    )
  ) +
  ggplot2::scale_size_continuous(name = "-Log10(FDR)", range = c(1.8, 5.5)) +
  ggplot2::labs(
    x = "NES",
    y = NULL,
    color = "Enriched in",
    title = "Global MSigDB GSEA screening"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 7),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
  )

print(p_global_gsea)

ggplot2::ggsave(
  "figures/host_RNAseq_global_MSigDB_GSEA_top_terms.svg",
  p_global_gsea,
  width = 8.2,
  height = max(5.5, 0.16 * nrow(global_dot_df) + 1.8),
  device = "svg"
)

save(
  fgsea_global_raw,
  fgsea_global_res,
  fgsea_global_top_pCR,
  fgsea_global_top_non_pCR,
  global_dot_df,
  file = "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/global_gsea/GSEA_global_MSigDB_results.RData"
)
