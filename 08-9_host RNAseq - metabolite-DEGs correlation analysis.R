#-----------------------------------------------------------------#
#
# Fig. 5D-02. Candidate review, immune-module summary, and manual scatter panel
#
# Purpose:
#   Script 1에서 저장한 metabolite–host DEG screening 결과를 이용하여
#   1. 품질이 낮은 metabolite/gene 및 해석하기 어려운 gene을 제외하고,
#   2. DEG 여부와 무관하게 expression-QC를 통과한 면역 유전자를 포함하여
#      폭넓은 metabolite–host gene candidate table을 정리하고,
#   3. 비중복 면역 기능 모듈 score와 metabolite의 세 가지 상관을 계산하고,
#   4. Pearson/Spearman/partial Spearman을 3분할 원으로 표시하고,
#   5. RORC–indolepropionic acid 연관성을 별도 확인하며,
#   6. 사용자가 CSV에서 직접 선택한 pair만 최종 scatter panel로 그립니다.
#
# Important:
#   - Script 1의 correlation 계산은 다시 수행하지 않습니다.
#   - correlation candidate cutoff:
#       p < 0.20 OR absolute correlation >= 0.30
#   - main review는 rank-based association인 Spearman과
#     partial Spearman을 중심으로 합니다.
#   - Pearson 결과는 candidate table에 보존하지만 main bubble/scatter의
#     핵심 통계량으로 반복 표시하지 않습니다.
#   - host DEG status는 annotation으로 보존하지만 필터로 사용하지 않습니다.
#   - 모든 candidate의 개별 scatter는 생성하지 않습니다.
#   - 전체 candidate bubble은 기본적으로 생성하지 않으며, 요청 시 상위 N개만 표시합니다.
#   - pseudogene은 gene biotype 정보가 없으므로 보수적인 symbol pattern과
#     manual exclusion list만 사용합니다.
#
# Required input:
#   host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/
#     Fig5D_metabolite_host_priority_screening_results.RData
#
# Main outputs:
#   tables/Fig5D_candidate_all_methods_<analysis_set>.csv
#   tables/Fig5D_candidate_rank_based_<analysis_set>.csv
#   tables/Fig5D_manual_scatter_pairs_<analysis_set>.csv
#   tables/Fig5D_immune_module_gene_coverage_<analysis_set>.csv
#   tables/Fig5D_immune_module_coverage_summary_<analysis_set>.csv
#   tables/Fig5D_immune_module_all_tested_pairs_<analysis_set>.csv
#   tables/Fig5D_immune_module_gene_trends_<analysis_set>.csv
#   tables/Fig5D_immune_module_metabolite_summary_<analysis_set>.csv
#   tables/Fig5D_immune_module_coordinated_metabolites_<analysis_set>.csv
#   Fig5D_candidate_review_<analysis_set>.xlsx
#   figures/Fig5D_candidate_bubble_top_<N>_<analysis_set>.svg  [optional]
#   figures/Fig5D_immune_module_metabolite_pies_<analysis_set>.svg
#   figures/RORC_indolepropionic_acid_<analysis_set>.svg
#
# After manual selection:
#   figures/Fig5D_manual_scatter_panel_<analysis_set>.svg
#   figures/Fig5D_manual_scatter_vertical_<analysis_set>.svg
#   tables/Fig5D_manual_selected_pair_summary_<analysis_set>.csv
#
#-----------------------------------------------------------------#

rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("D:/2-연구/2-CRC metagenomics/")

for (pkg in c(
  "dplyr", "tidyr", "tibble", "stringr",
  "ggplot2", "svglite", "openxlsx", "ggforce", "patchwork"
)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, type = "binary")
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(svglite)
  library(openxlsx)
  library(ggforce)
  library(patchwork)
})

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create(
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

full_screening_rdata <- paste0(
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
  "Fig5D_microbe_metabolite_host_axis_screening_results.RData"
)

compatibility_rdata <- paste0(
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
  "Fig5D_metabolite_host_priority_screening_results.RData"
)

if (file.exists(full_screening_rdata)) {
  load(full_screening_rdata)
} else if (file.exists(compatibility_rdata)) {
  load(compatibility_rdata)
} else {
  stop(
    "Neither the full screening RData nor compatibility RData exists.",
    call. = FALSE
  )
}

#-----------------------------------------------------------------#
# 1. User-adjustable settings
#-----------------------------------------------------------------#

analysis_set <- "before"
# analysis_set <- "all_available"

correlation_p_cutoff <- 0.20
minimum_abs_correlation <- 0.30

# Metabolite IQR filtering is retained. Host genes are not subjected to an
# additional IQR quantile cutoff here because expression QC was already applied
# in Script 1. This preserves non-DEG immune genes with modest variability.
metabolite_iqr_quantile_cutoff <- 0.10
gene_iqr_quantile_cutoff <- 0.00

# Only strongly floor-dominated host genes are excluded.
gene_floor_fraction_cutoff <- 0.80

# Metabolite floor filters must match Script 1. If an older RData file lacks
# basal_cluster_fraction, the corresponding check is skipped conservatively.
met_exact_floor_fraction_cutoff <- 0.65
met_basal_cluster_fraction_cutoff <- 0.50

# Module-score and supplementary-plot settings.
minimum_module_genes_for_score <- 2

# Choose one curated module definition. Both definitions are included below.
#   "strict": small, specific, non-overlapping marker sets.
#   "broad" : wider, still non-overlapping sets containing genes whose higher
#             expression consistently supports the named cell state/program.
# module_gene_set_version <- "strict"
module_gene_set_version <- "broad"

# A coordinated trend is defined only by directional agreement of Pearson,
# Spearman, and partial Spearman. P values do not determine the outline;
# they determine circle radius and method-specific symbols.
# Circle size represents -log10 of the descriptive geometric-mean P value.

# One-character significance symbols are assigned independently to each
# correlation-method sector. These symbols do not represent a combined P value.
sector_symbol_p05 <- "*"
sector_symbol_p01 <- "#"
sector_symbol_p001 <- "$"

draw_immune_support_scatter <- TRUE
immune_support_scatter_gmean_p_cutoff <- 0.10
immune_support_scatter_max_pairs <- 24
immune_support_scatter_ncol <- 4
plot_metabolite_superclasses <- c(
  "SCFA / related",
  "Bile acid",
  "Tryptophan / indole"
)

# The all-candidate bubble plot becomes unreadable when many pairs pass the
# exploratory cutoff. It is therefore disabled by default. When enabled,
# only the highest-ranked candidates are displayed; all pairs remain in the
# CSV and Excel review tables.
draw_candidate_bubble <- FALSE
candidate_bubble_top_n <- 20

# Number of columns only for the manually selected combined scatter panel.
# This does not limit the number of selected pairs.
scatter_panel_ncol <- 3

group_cols <- c(
  "non_pCR" = "#E07A73",
  "pCR" = "#5AB49B"
)

gene_annotation_cols <- c(
  "Radiation / DNA damage / apoptosis" = "#D73027",
  "Chemoresistance / EMT / stroma" = "#7B3294",
  "Immune / ICI / inflammation" = "#F46D43",
  "Barrier / epithelial / AHR" = "#1B9E77",
  "Metabolite receptor / transport" = "#66A61E",
  "Cell cycle / proliferation" = "#E6AB02",
  "Other annotated" = "#BDBDBD"
)

# Add genes here only when they are unsuitable for interpretation or are
# confirmed pseudogene/poorly characterized candidates.
manual_gene_exclude <- c(
  "LRRC36"
)

manual_pair_file <- paste0(
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
  "Fig5D_manual_scatter_pairs_",
  analysis_set,
  ".csv"
)


#-----------------------------------------------------------------#
# 1B. Strict and broad non-overlapping immune functional modules
#-----------------------------------------------------------------#
#
# Each gene is assigned to only one displayed module within each version.
# All included genes are positive markers/regulators of the named state:
# higher expression contributes positively to the module score. The broad
# version expands coverage but deliberately avoids genes assigned elsewhere.
# NK is displayed first, followed by merged CD8 identity/cytotoxicity.
# The previous Tumor_killing_response module is removed because its mixed
# ligand/receptor composition was not readily distinguishable from NK/CD8.
#-----------------------------------------------------------------#

immune_gene_sets_strict <- list(
  NK_core = c(
    "KLRD1", "NCR1", "KLRK1", "KLRC1", "GNLY", "XCL1"
  ),
  
  CD8_effector = c(
    "CD3D", "CD3E", "TRAC", "CD8A", "CD8B", "RUNX3",
    "PRF1", "GZMB", "CTSW", "CCL5"
  ),
  
  CD8_exhaustion = c(
    "PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT", "ENTPD1", "LAYN"
  ),
  
  Th1_core = c(
    "TBX21", "IFNG", "CXCR3", "IL12RB2"
  ),
  
  Th2_core = c(
    "GATA3", "IL4", "IL5", "IL13", "PTGDR2"
  ),
  
  Th17_core = c(
    "RORC", "IL17A", "IL17F", "CCR6", "IL23R"
  ),
  
  Treg_core = c(
    "FOXP3", "IL2RA", "CTLA4", "IKZF2", "CCR8", "TNFRSF18"
  ),
  
  Neutrophil_core = c(
    "FCGR3B", "CEACAM8", "CSF3R", "S100A12", "FPR1", "CXCR1"
  ),
  
  MDSC_core = c(
    "OLR1", "ARG1", "S100A8", "S100A9", "CEBPB", "STAT3"
  ),
  
  cDC1 = c(
    "CLEC9A", "XCR1", "BATF3", "WDFY4", "IRF8"
  ),
  
  TLS_Bcell = c(
    "CXCL13", "CCL19", "CCL21", "MS4A1", "CD79A"
  )
)

immune_gene_sets_broad <- list(
  NK_core = c(
    "KLRD1", "NCR1", "KLRK1", "KLRC1", "KLRC2", "KLRF1",
    "GNLY", "FGFBP2", "XCL1", "XCL2", "FCGR3A", "FCER1G", "TYROBP"
  ),
  
  CD8_effector = c(
    "CD3D", "CD3E", "TRAC", "CD8A", "CD8B", "RUNX3",
    "PRF1", "GZMA", "GZMB", "GZMH", "GZMK", "CTSW", "CCL5",
    "FASLG", "TNFSF10", "IL2", "TNF", "LTA"
  ),
  
  CD8_exhaustion = c(
    "PDCD1", "TOX", "TOX2", "LAG3", "HAVCR2", "TIGIT",
    "ENTPD1", "LAYN", "CD244", "CD160", "NR4A1", "NR4A2",
    "NR4A3", "EOMES", "PRDM1"
  ),
  
  Th1_core = c(
    "TBX21", "IFNG", "CXCR3", "IL12RB2", "IL18R1",
    "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11"
  ),
  
  Th2_core = c(
    "GATA3", "IL4", "IL5", "IL13", "PTGDR2", "CCR4",
    "IL4R", "STAT6"
  ),
  
  Th17_core = c(
    "RORC", "IL17A", "IL17F", "CCR6", "IL23R", "KLRB1",
    "CCL20", "RORA", "IL21", "IL22", "IL26", "CSF2", "BHLHE40", "AHR"
  ),
  
  Treg_core = c(
    "FOXP3", "IL2RA", "CTLA4", "IKZF2", "CCR8", "TNFRSF18",
    "ICOS", "LRRC32", "FGL2", "IL10", "EBI3", "NT5E"
  ),
  
  Neutrophil_core = c(
    "FCGR3B", "CEACAM8", "CSF3R", "S100A12", "FPR1", "FPR2",
    "CXCR1", "MMP8", "MMP9", "LTF", "CAMP", "LCN2", "OLFM4",
    "MNDA", "SELL"
  ),
  
  MDSC_core = c(
    "OLR1", "ARG1", "S100A8", "S100A9", "CEBPB", "STAT3",
    "CXCR2", "PTGS2", "CYBB", "NCF1", "NCF2", "NCF4", "IDO1",
    "SLC7A2", "IL1B", "VCAN", "FCN1", "LILRB1", "CD14"
  ),
  
  cDC1 = c(
    "CLEC9A", "XCR1", "BATF3", "WDFY4", "IRF8", "CADM1",
    "THBD", "DNASE1L3", "SNX22", "CPVL", "CLNK"
  ),
  
  TLS_Bcell = c(
    "CXCL13", "CCL19", "CCL21", "MS4A1", "CD79A", "CD79B",
    "CD74", "IGKC", "BCL6", "IL21R", "CD37", "BANK1", "MZB1", "JCHAIN"
  )
)

if (!module_gene_set_version %in% c("strict", "broad")) {
  stop(
    "module_gene_set_version must be either 'strict' or 'broad'.",
    call. = FALSE
  )
}

immune_gene_sets <- switch(
  module_gene_set_version,
  strict = immune_gene_sets_strict,
  broad = immune_gene_sets_broad
)

module_output_tag <- paste0(
  analysis_set,
  "_",
  module_gene_set_version
)

immune_module_groups <- c(
  "NK_core" = "NK",
  "CD8_effector" = "CD8 T cell",
  "CD8_exhaustion" = "CD8 T cell",
  "Th1_core" = "Th1",
  "Th2_core" = "Th2",
  "Th17_core" = "Th17",
  "Treg_core" = "Treg",
  "Neutrophil_core" = "Neutrophil",
  "MDSC_core" = "MDSC",
  "cDC1" = "APC",
  "TLS_Bcell" = "TLS / B cell"
)

immune_master_regulators <- list(
  CD8_effector = c("CD8A", "CD8B", "RUNX3"),
  Th1_core = "TBX21",
  Th2_core = "GATA3",
  Th17_core = "RORC",
  Treg_core = "FOXP3",
  MDSC_core = c("CEBPB", "STAT3")
)

module_gene_duplicates <- enframe(
  immune_gene_sets,
  name = "Immune_module",
  value = "Gene"
) %>%
  unnest_longer(Gene) %>%
  count(Gene, name = "n_modules") %>%
  filter(n_modules > 1)

if (nrow(module_gene_duplicates) > 0) {
  stop(
    "Displayed immune modules contain duplicated genes in the ",
    module_gene_set_version,
    " definition: ",
    paste(module_gene_duplicates$Gene, collapse = ", "),
    call. = FALSE
  )
}

for (module_i in names(immune_master_regulators)) {
  if (!all(
    immune_master_regulators[[module_i]] %in%
    immune_gene_sets[[module_i]]
  )) {
    stop(
      "Required lineage regulator is missing from ",
      module_i,
      call. = FALSE
    )
  }
}

#-----------------------------------------------------------------#
# 2. Check input objects
#-----------------------------------------------------------------#

if (
  !exists("met_gene_all_results") ||
  is.null(met_gene_all_results[[analysis_set]]) ||
  nrow(met_gene_all_results[[analysis_set]]) == 0
) {
  stop(
    "met_gene_all_results for the selected analysis_set is missing or empty.",
    call. = FALSE
  )
}

if (
  !exists("sample_metadata_results") ||
  is.null(sample_metadata_results[[analysis_set]]) ||
  nrow(sample_metadata_results[[analysis_set]]) == 0
) {
  stop(
    "sample_metadata_results for the selected analysis_set is missing or empty.",
    call. = FALSE
  )
}

if (
  !exists("met_use_log_results") ||
  is.null(met_use_log_results[[analysis_set]]) ||
  ncol(met_use_log_results[[analysis_set]]) == 0
) {
  stop(
    "met_use_log_results for the selected analysis_set is missing or empty.",
    call. = FALSE
  )
}

if (
  !exists("vst_deg_results") ||
  is.null(vst_deg_results[[analysis_set]]) ||
  ncol(vst_deg_results[[analysis_set]]) == 0
) {
  stop(
    "vst_deg_results for the selected analysis_set is missing or empty.",
    call. = FALSE
  )
}

col_use <- sample_metadata_results[[analysis_set]]
met_use_log <- met_use_log_results[[analysis_set]]
vst_deg <- vst_deg_results[[analysis_set]]
candidate_input <- met_gene_all_results[[analysis_set]]

triad_input_object_report <- tibble(
  object_name = c(
    "sample_metadata_results",
    "met_use_log_results",
    "met_use_raw_results",
    "vst_deg_results",
    "species_use_log_results",
    "species_use_raw_results"
  ),
  available = c(
    exists("sample_metadata_results") &&
      !is.null(sample_metadata_results[[analysis_set]]),
    exists("met_use_log_results") &&
      !is.null(met_use_log_results[[analysis_set]]),
    exists("met_use_raw_results") &&
      !is.null(met_use_raw_results[[analysis_set]]),
    exists("vst_deg_results") &&
      !is.null(vst_deg_results[[analysis_set]]),
    exists("species_use_log_results") &&
      !is.null(species_use_log_results[[analysis_set]]) &&
      ncol(species_use_log_results[[analysis_set]]) > 0,
    exists("species_use_raw_results") &&
      !is.null(species_use_raw_results[[analysis_set]]) &&
      ncol(species_use_raw_results[[analysis_set]]) > 0
  )
)

write.csv(
  triad_input_object_report,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/diagnostics/",
    "Fig5D_9.8.2_triad_input_object_report_",
    analysis_set,
    ".csv"
  ),
  row.names = FALSE
)

if (!exists("met_use_raw_results")) {
  met_use_raw_results <- list()
}
if (!exists("species_use_log_results")) {
  species_use_log_results <- list()
}
if (!exists("species_use_raw_results")) {
  species_use_raw_results <- list()
}
if (!exists("species_qc_results")) {
  species_qc_results <- list()
}
if (!exists("species_met_results")) {
  species_met_results <- list()
}

if (!all(triad_input_object_report$available)) {
  warning(
    "Some full triad-recalculation matrices are unavailable. ",
    "Rerun the revised 7-9.8.1 script before 9.8.3. Missing: ",
    paste(
      triad_input_object_report$object_name[
        !triad_input_object_report$available
      ],
      collapse = ", "
    )
  )
}

if (
  !any(unique(unlist(immune_gene_sets)) %in% candidate_input$Gene)
) {
  warning(
    "No curated immune-module gene is present in candidate_input. ",
    "Rerun the revised Script 1 so that non-DEG immune genes are correlated."
  )
}

if (!all(c(
  "Gene", "Metabolite",
  "pearson_r", "pearson_p",
  "spearman_rho", "spearman_p",
  "partial_spearman_rho", "partial_spearman_p"
) %in% colnames(candidate_input))) {
  stop(
    "Required metabolite–host correlation columns are missing from met_gene_all_results.",
    call. = FALSE
  )
}

# Optional columns are given conservative defaults when absent.
if (!"host_DEG_selected" %in% colnames(candidate_input)) {
  candidate_input$host_DEG_selected <- NA
}
if (!"host_immune_target" %in% colnames(candidate_input)) {
  candidate_input$host_immune_target <-
    candidate_input$Gene %in% unique(unlist(immune_gene_sets))
}
if (!"host_housekeeping_control" %in% colnames(candidate_input)) {
  candidate_input$host_housekeeping_control <- FALSE
}
if (!"keep_gene_for_correlation" %in% colnames(candidate_input)) {
  candidate_input$keep_gene_for_correlation <- TRUE
}
if (!"host_gene_selected" %in% colnames(candidate_input)) {
  candidate_input$host_gene_selected <-
    candidate_input$keep_gene_for_correlation
}
if (!"host_selection_source" %in% colnames(candidate_input)) {
  candidate_input$host_selection_source <- case_when(
    coalesce(candidate_input$host_DEG_selected, FALSE) &
      coalesce(candidate_input$host_immune_target, FALSE) ~
      "DEG + curated immune module",
    coalesce(candidate_input$host_DEG_selected, FALSE) ~ "DEG",
    coalesce(candidate_input$host_immune_target, FALSE) ~
      "Curated immune module; DEG status unavailable or non-DEG",
    TRUE ~ "Expression-QC-passed host gene"
  )
}
if (!"met_DE_selected" %in% colnames(candidate_input)) {
  candidate_input$met_DE_selected <- FALSE
}
if (!"met_iqr" %in% colnames(candidate_input)) {
  candidate_input$met_iqr <- NA_real_
}
if (!"met_sd" %in% colnames(candidate_input)) {
  candidate_input$met_sd <- NA_real_
}
if (!"min_floor_fraction" %in% colnames(candidate_input)) {
  candidate_input$min_floor_fraction <- NA_real_
}
if (!"basal_cluster_fraction" %in% colnames(candidate_input)) {
  candidate_input$basal_cluster_fraction <- NA_real_
}
if (!"gene_iqr" %in% colnames(candidate_input)) {
  candidate_input$gene_iqr <- NA_real_
}
if (!"gene_sd" %in% colnames(candidate_input)) {
  candidate_input$gene_sd <- NA_real_
}
if (!"gene_min_floor_fraction" %in% colnames(candidate_input)) {
  candidate_input$gene_min_floor_fraction <- NA_real_
}
if (!"host_logFC" %in% colnames(candidate_input)) {
  candidate_input$host_logFC <- NA_real_
}
if (!"host_DESeq2_p" %in% colnames(candidate_input)) {
  candidate_input$host_DESeq2_p <- NA_real_
}
if (!"met_logFC" %in% colnames(candidate_input)) {
  candidate_input$met_logFC <- NA_real_
}
if (!"met_p" %in% colnames(candidate_input)) {
  candidate_input$met_p <- NA_real_
}
if (!"Metabolite_class" %in% colnames(candidate_input)) {
  candidate_input$Metabolite_class <- "Other metabolite"
}
if (!"host_oncology_category" %in% colnames(candidate_input)) {
  candidate_input$host_oncology_category <- "Other host gene"
}
if (!"Host_axis_summary" %in% colnames(candidate_input)) {
  candidate_input$Host_axis_summary <- "Not assigned"
}
if (!"pathway_label_summary" %in% colnames(candidate_input)) {
  candidate_input$pathway_label_summary <- "Not assigned"
}
if (!"Host_relevance_note" %in% colnames(candidate_input)) {
  candidate_input$Host_relevance_note <- "No specific note"
}
if (!"Metabolite_relevance_note" %in% colnames(candidate_input)) {
  candidate_input$Metabolite_relevance_note <- "No specific note"
}
if (!"top_microbe_support" %in% colnames(candidate_input)) {
  candidate_input$top_microbe_support <- "No weak species-metabolite support"
}
if (!"microbe_support_score" %in% colnames(candidate_input)) {
  candidate_input$microbe_support_score <- 0
}
if (!"priority_score" %in% colnames(candidate_input)) {
  candidate_input$priority_score <- 0
}
if (!"host_oncology_score" %in% colnames(candidate_input)) {
  candidate_input$host_oncology_score <- 0
}
if (!"metabolite_response_score" %in% colnames(candidate_input)) {
  candidate_input$metabolite_response_score <- 0
}
if (!"pathway_match_score" %in% colnames(candidate_input)) {
  candidate_input$pathway_match_score <- 0
}
if (!"direction_consistency_score" %in% colnames(candidate_input)) {
  candidate_input$direction_consistency_score <- 0
}
if (!"lowest_point_sensitivity" %in% colnames(candidate_input)) {
  candidate_input$lowest_point_sensitivity <- "not tested"
}
if (!"spearman_rho_no_lowest" %in% colnames(candidate_input)) {
  candidate_input$spearman_rho_no_lowest <- NA_real_
}
if (!"spearman_p_no_lowest" %in% colnames(candidate_input)) {
  candidate_input$spearman_p_no_lowest <- NA_real_
}

#-----------------------------------------------------------------#
# 3. Determine defensive IQR thresholds
#-----------------------------------------------------------------#
#
# These thresholds are calculated at the feature level, not from duplicated
# pair rows. A cutoff of 0 means that no additional IQR filtering is applied.
#-----------------------------------------------------------------#

metabolite_iqr_cutoff <- 0
gene_iqr_cutoff <- 0

if (
  exists("met_qc_results") &&
  !is.null(met_qc_results[[analysis_set]]) &&
  "met_iqr" %in% colnames(met_qc_results[[analysis_set]])
) {
  metabolite_iqr_cutoff <- quantile(
    met_qc_results[[analysis_set]]$met_iqr[
      is.finite(met_qc_results[[analysis_set]]$met_iqr) &
        met_qc_results[[analysis_set]]$met_iqr > 0
    ],
    metabolite_iqr_quantile_cutoff,
    na.rm = TRUE
  )
}

if (
  exists("host_expr_qc_results") &&
  !is.null(host_expr_qc_results[[analysis_set]]) &&
  "gene_iqr" %in% colnames(host_expr_qc_results[[analysis_set]])
) {
  gene_iqr_cutoff <- quantile(
    host_expr_qc_results[[analysis_set]]$gene_iqr[
      is.finite(host_expr_qc_results[[analysis_set]]$gene_iqr) &
        host_expr_qc_results[[analysis_set]]$gene_iqr > 0
    ],
    gene_iqr_quantile_cutoff,
    na.rm = TRUE
  )
}

if (!is.finite(metabolite_iqr_cutoff)) {
  metabolite_iqr_cutoff <- 0
}
if (!is.finite(gene_iqr_cutoff)) {
  gene_iqr_cutoff <- 0
}

message(
  "Additional visualization-stage IQR cutoffs: metabolite = ",
  signif(metabolite_iqr_cutoff, 4),
  "; gene = ",
  signif(gene_iqr_cutoff, 4)
)

housekeeping_control_pairs <- candidate_input %>%
  filter(coalesce(host_housekeeping_control, FALSE)) %>%
  arrange(Metabolite, Gene)

write.csv(
  housekeeping_control_pairs,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_housekeeping_control_pairs_",
    analysis_set,
    ".csv"
  ),
  row.names = FALSE
)

#-----------------------------------------------------------------#
# 4. Build quality-controlled candidate tables
#-----------------------------------------------------------------#
#
# Two candidate definitions are retained:
#
# all_method_pass:
#   Pearson, Spearman, or partial Spearman meets p/effect-size cutoff.
#
# rank_based_pass:
#   Spearman or partial Spearman meets p/effect-size cutoff.
#
# The rank-based table is used for the main bubble and scatter review because
# it is less dependent on linearity and less redundant with Panel C.
# Pearson statistics remain available in the exported review table.
#-----------------------------------------------------------------#

candidate_all_methods <- candidate_input %>%
  filter(
    Gene %in% colnames(vst_deg),
    Metabolite %in% colnames(met_use_log),
    coalesce(host_gene_selected, TRUE),
    coalesce(keep_gene_for_correlation, TRUE),
    !is.na(Gene),
    !is.na(Metabolite),
    Gene != "",
    Metabolite != "",
    !Gene %in% manual_gene_exclude,
    !coalesce(host_housekeeping_control, FALSE),
    is.na(met_iqr) | met_iqr >= metabolite_iqr_cutoff,
    is.na(gene_iqr) | gene_iqr >= gene_iqr_cutoff,
    is.na(met_sd) | met_sd > 0,
    is.na(gene_sd) | gene_sd > 0,
    is.na(min_floor_fraction) |
      min_floor_fraction < met_exact_floor_fraction_cutoff,
    is.na(basal_cluster_fraction) |
      basal_cluster_fraction < met_basal_cluster_fraction_cutoff,
    is.na(gene_min_floor_fraction) |
      coalesce(host_immune_target, FALSE) |
      gene_min_floor_fraction < gene_floor_fraction_cutoff
  ) %>%
  mutate(
    Gene = as.character(Gene),
    Metabolite = as.character(Metabolite),
    
    # Conservative pattern only. A generic P[0-9]+ suffix is intentionally
    # not used because it incorrectly removes real genes such as TNFAIP1.
    probable_pseudogene_or_low_annotation = str_detect(
      Gene,
      regex(
        paste(
          c(
            "^LOC[0-9]",
            "^LINC[0-9]",
            "^MIR[0-9]",
            "^SNORA",
            "^SNORD",
            "^RNU",
            "^AC[0-9]",
            "^AL[0-9]",
            "^AP[0-9]",
            "^RP[0-9]",
            "-AS[0-9]$"
          ),
          collapse = "|"
        ),
        ignore_case = FALSE
      )
    ),
    
    max_abs_rank_cor = pmax(
      abs(spearman_rho),
      abs(partial_spearman_rho),
      na.rm = TRUE
    ),
    max_abs_rank_cor = ifelse(
      is.infinite(max_abs_rank_cor),
      NA_real_,
      max_abs_rank_cor
    ),
    
    min_rank_cor_p = pmin(
      spearman_p,
      partial_spearman_p,
      na.rm = TRUE
    ),
    min_rank_cor_p = ifelse(
      is.infinite(min_rank_cor_p),
      NA_real_,
      min_rank_cor_p
    ),
    
    all_method_pass =
      coalesce(pearson_p < correlation_p_cutoff, FALSE) |
      coalesce(spearman_p < correlation_p_cutoff, FALSE) |
      coalesce(partial_spearman_p < correlation_p_cutoff, FALSE) |
      coalesce(abs(pearson_r) >= minimum_abs_correlation, FALSE) |
      coalesce(abs(spearman_rho) >= minimum_abs_correlation, FALSE) |
      coalesce(abs(partial_spearman_rho) >= minimum_abs_correlation, FALSE),
    
    rank_based_pass =
      coalesce(spearman_p < correlation_p_cutoff, FALSE) |
      coalesce(partial_spearman_p < correlation_p_cutoff, FALSE) |
      coalesce(abs(spearman_rho) >= minimum_abs_correlation, FALSE) |
      coalesce(abs(partial_spearman_rho) >= minimum_abs_correlation, FALSE),
    
    rank_direction_consistent =
      !is.na(spearman_rho) &
      !is.na(partial_spearman_rho) &
      sign(spearman_rho) == sign(partial_spearman_rho),
    
    adjusted_effect_retained =
      !is.na(spearman_rho) &
      !is.na(partial_spearman_rho) &
      sign(spearman_rho) == sign(partial_spearman_rho) &
      abs(partial_spearman_rho) >= 0.70 * abs(spearman_rho),
    
    lowest_point_stable = !lowest_point_sensitivity %in% c(
      "sign changed after lowest-met point removal",
      "sensitive to lowest-met point",
      "p-value sensitive to lowest-met point"
    ),
    
    annotation_text = paste(
      host_oncology_category,
      Host_axis_summary,
      pathway_label_summary,
      Host_relevance_note,
      sep = "; "
    ),
    
    gene_annotation = case_when(
      str_detect(
        annotation_text,
        regex(
          "radiation|radiotherapy|DNA damage|GADD45|p53|apoptosis",
          ignore_case = TRUE
        )
      ) ~ "Radiation / DNA damage / apoptosis",
      
      str_detect(
        annotation_text,
        regex(
          "chemotherapy|chemo|drug resistance|chemoresistance|EMT|matrix|ECM|stromal|invasion|metastasis",
          ignore_case = TRUE
        )
      ) ~ "Chemoresistance / EMT / stroma",
      
      str_detect(
        annotation_text,
        regex(
          "immune checkpoint|ICI|PD-1|PDL1|PD-L1|CTLA|interferon|T cell|lymphocyte|JAK|STAT|TNF|NF|cytokine|immune|leukocyte",
          ignore_case = TRUE
        )
      ) ~ "Immune / ICI / inflammation",
      
      str_detect(
        annotation_text,
        regex(
          "barrier|mucus|goblet|epithelial|AHR|IL22",
          ignore_case = TRUE
        )
      ) ~ "Barrier / epithelial / AHR",
      
      str_detect(
        annotation_text,
        regex(
          "metabolism|glycolysis|transport|receptor|SCFA|niacin|bile acid|HCAR|FFAR|FXR|TGR5|hypoxia",
          ignore_case = TRUE
        )
      ) ~ "Metabolite receptor / transport",
      
      str_detect(
        annotation_text,
        regex(
          "cell cycle|G2/M|proliferation|MYC",
          ignore_case = TRUE
        )
      ) ~ "Cell cycle / proliferation",
      
      TRUE ~ "Other annotated"
    ),
    
    radiation_relevance = str_detect(
      annotation_text,
      regex(
        "radiation|radiotherapy|DNA damage|GADD45|p53|apoptosis",
        ignore_case = TRUE
      )
    ),
    
    chemoresistance_relevance = str_detect(
      annotation_text,
      regex(
        "chemotherapy|chemo|drug resistance|chemoresistance|EMT|matrix|ECM|stromal|hypoxia",
        ignore_case = TRUE
      )
    ),
    
    ici_relevance = str_detect(
      annotation_text,
      regex(
        "immune checkpoint|ICI|PD-1|PDL1|PD-L1|CTLA|interferon|T cell|lymphocyte|immune",
        ignore_case = TRUE
      )
    ),
    
    crc_relevance =
      coalesce(host_oncology_score > 0, FALSE) |
      radiation_relevance |
      chemoresistance_relevance |
      ici_relevance,
    
    microbe_supported =
      coalesce(microbe_support_score > 0, FALSE) |
      (
        !is.na(top_microbe_support) &
          top_microbe_support != "" &
          top_microbe_support != "No weak species-metabolite support"
      ),
    
    correlation_review_score = case_when(
      min_rank_cor_p < 0.01 & max_abs_rank_cor >= 0.50 ~ 4,
      min_rank_cor_p < 0.05 & max_abs_rank_cor >= 0.40 ~ 3,
      min_rank_cor_p < 0.10 | max_abs_rank_cor >= 0.40 ~ 2,
      min_rank_cor_p < correlation_p_cutoff |
        max_abs_rank_cor >= minimum_abs_correlation ~ 1,
      TRUE ~ 0
    ),
    
    biological_review_score =
      correlation_review_score +
      coalesce(host_oncology_score, 0) +
      coalesce(metabolite_response_score, 0) +
      coalesce(pathway_match_score, 0) +
      ifelse(rank_direction_consistent, 1, 0) +
      ifelse(adjusted_effect_retained, 1, 0) +
      ifelse(lowest_point_stable, 1, -1) +
      ifelse(radiation_relevance, 2, 0) +
      ifelse(chemoresistance_relevance, 2, 0) +
      ifelse(ici_relevance, 2, 0) +
      ifelse(microbe_supported, 1, 0),
    
    review_tier = case_when(
      biological_review_score >= 12 ~ "High",
      biological_review_score >= 8 ~ "Moderate",
      biological_review_score >= 4 ~ "Putative",
      TRUE ~ "Exploratory"
    ),
    
    pair_label = paste0(
      str_replace_all(Metabolite, "_", " "),
      " | ",
      Gene
    )
  ) %>%
  filter(
    all_method_pass,
    !probable_pseudogene_or_low_annotation
  ) %>%
  distinct(Metabolite, Gene, .keep_all = TRUE) %>%
  arrange(
    desc(biological_review_score),
    desc(crc_relevance),
    desc(microbe_supported),
    desc(rank_direction_consistent),
    desc(lowest_point_stable),
    min_rank_cor_p,
    desc(max_abs_rank_cor),
    host_DESeq2_p,
    met_p
  ) %>%
  mutate(
    all_method_rank = row_number()
  )

candidate_rank_based <- candidate_all_methods %>%
  filter(rank_based_pass) %>%
  arrange(
    desc(biological_review_score),
    desc(crc_relevance),
    desc(microbe_supported),
    desc(rank_direction_consistent),
    desc(lowest_point_stable),
    min_rank_cor_p,
    desc(max_abs_rank_cor),
    host_DESeq2_p,
    met_p
  ) %>%
  mutate(
    rank_based_rank = row_number()
  )


#-----------------------------------------------------------------#
# 5. Immune module-score correlations and supporting genes
#-----------------------------------------------------------------#
#
# The primary module statistic is now a sample-level mean z-score across the
# available genes in each non-overlapping module. Each metabolite is correlated
# directly with this module score using Pearson, Spearman, and covariate-
# adjusted partial Spearman. This is more interpretable than taking a median
# of preselected gene-wise correlations and avoids selection-induced direction
# bias. Gene-wise associations are retained as supporting evidence and for
# selecting supplementary scatter plots.
#-----------------------------------------------------------------#

geometric_mean_p <- function(x) {
  x <- x[is.finite(x) & x >= 0]
  if (length(x) == 0) {
    return(NA_real_)
  }
  exp(mean(log(pmax(x, .Machine$double.xmin))))
}

classify_metabolite_superclass <- function(x) {
  x_std <- tolower(str_replace_all(x, "_", " "))
  
  case_when(
    x_std %in% c(
      "acetate", "propionate", "butyrate", "isobutyrate",
      "isovalerate", "valerate", "gamma aminobutyric acid"
    ) ~ "SCFA / related",
    
    str_detect(
      x_std,
      "cholic|deoxycholic|ursodeoxycholic|lithocholic|tauro|glyco"
    ) ~ "Bile acid",
    
    x_std %in% c(
      "indolepropionic acid", "indole acetic acid",
      "indole lactic acid", "indole", "tryptamine",
      "tryptophan", "kynurenic acid", "xanthurenic acid",
      "nicotinic acid"
    ) ~ "Tryptophan / indole",
    
    TRUE ~ "Other"
  )
}

calculate_module_correlations <- function(x, y, metadata, covariates) {
  ok <- is.finite(x) & is.finite(y)
  
  result <- list(
    n = sum(ok),
    pearson_r = NA_real_,
    pearson_p = NA_real_,
    spearman_rho = NA_real_,
    spearman_p = NA_real_,
    partial_spearman_rho = NA_real_,
    partial_spearman_p = NA_real_,
    partial_covariates = NA_character_
  )
  
  if (
    sum(ok) < 8 ||
    sd(x[ok], na.rm = TRUE) <= 0 ||
    sd(y[ok], na.rm = TRUE) <= 0
  ) {
    return(result)
  }
  
  pearson_result <- suppressWarnings(
    cor.test(x[ok], y[ok], method = "pearson")
  )
  spearman_result <- suppressWarnings(
    cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  )
  
  result$pearson_r <- unname(pearson_result$estimate)
  result$pearson_p <- pearson_result$p.value
  result$spearman_rho <- unname(spearman_result$estimate)
  result$spearman_p <- spearman_result$p.value
  
  if (length(covariates) > 0) {
    partial_df <- data.frame(
      metabolite_rank = rank(x[ok], ties.method = "average"),
      module_rank = rank(y[ok], ties.method = "average"),
      metadata[ok, covariates, drop = FALSE],
      check.names = FALSE
    )
    
    partial_df <- partial_df[
      complete.cases(partial_df),
      ,
      drop = FALSE
    ]
    
    for (covariate_i in covariates) {
      if (
        covariate_i %in% colnames(partial_df) &&
        (
          is.character(partial_df[[covariate_i]]) ||
          is.factor(partial_df[[covariate_i]])
        )
      ) {
        partial_df[[covariate_i]] <- factor(
          partial_df[[covariate_i]]
        )
      }
    }
    
    covariates_use <- covariates[
      covariates %in% colnames(partial_df) &
        vapply(
          covariates,
          function(z) {
            length(unique(na.omit(partial_df[[z]]))) > 1
          },
          logical(1)
        )
    ]
    
    if (
      length(covariates_use) > 0 &&
      nrow(partial_df) > length(covariates_use) + 3
    ) {
      formula_metabolite <- as.formula(
        paste(
          "metabolite_rank ~",
          paste(covariates_use, collapse = " + ")
        )
      )
      formula_module <- as.formula(
        paste(
          "module_rank ~",
          paste(covariates_use, collapse = " + ")
        )
      )
      
      partial_r <- suppressWarnings(
        cor(
          residuals(lm(formula_metabolite, data = partial_df)),
          residuals(lm(formula_module, data = partial_df)),
          method = "pearson"
        )
      )
      
      if (is.finite(partial_r) && abs(partial_r) < 1) {
        partial_df_residual <-
          nrow(partial_df) - length(covariates_use) - 2
        
        partial_p <- 2 * pt(
          -abs(
            partial_r * sqrt(
              partial_df_residual / (1 - partial_r^2)
            )
          ),
          df = partial_df_residual
        )
        
        result$partial_spearman_rho <- partial_r
        result$partial_spearman_p <- partial_p
        result$partial_covariates <- paste(
          covariates_use,
          collapse = ";"
        )
      }
    }
  }
  
  result
}

immune_module_gene_map <- enframe(
  immune_gene_sets,
  name = "Immune_module",
  value = "Gene"
) %>%
  unnest_longer(Gene) %>%
  mutate(
    Gene = as.character(Gene),
    Module_group = unname(immune_module_groups[Immune_module]),
    Master_regulator = mapply(
      function(module_i, gene_i) {
        gene_i %in% immune_master_regulators[[module_i]]
      },
      Immune_module,
      Gene
    )
  ) %>%
  distinct(Module_group, Immune_module, Gene, .keep_all = TRUE)

immune_module_gene_coverage <- immune_module_gene_map %>%
  mutate(
    in_vst_deg = Gene %in% colnames(vst_deg),
    in_correlation_table = Gene %in% candidate_input$Gene,
    host_DEG_selected_in_input = Gene %in% unique(
      candidate_input$Gene[
        coalesce(candidate_input$host_DEG_selected, FALSE)
      ]
    ),
    expression_QC_passed_in_input = Gene %in% unique(
      candidate_input$Gene[
        coalesce(candidate_input$keep_gene_for_correlation, TRUE)
      ]
    ),
    available_for_module_score =
      in_vst_deg & expression_QC_passed_in_input
  ) %>%
  arrange(
    Module_group,
    Immune_module,
    desc(Master_regulator),
    desc(available_for_module_score),
    Gene
  )

immune_module_coverage_summary <- immune_module_gene_coverage %>%
  group_by(Module_group, Immune_module) %>%
  summarise(
    n_module_genes = n_distinct(Gene),
    n_genes_available = n_distinct(
      Gene[available_for_module_score]
    ),
    coverage_fraction = n_genes_available / n_module_genes,
    available_genes = paste(
      sort(unique(Gene[available_for_module_score])),
      collapse = "; "
    ),
    unavailable_genes = paste(
      sort(unique(Gene[!available_for_module_score])),
      collapse = "; "
    ),
    master_regulators = paste(
      sort(unique(Gene[Master_regulator])),
      collapse = "; "
    ),
    master_regulator_available = all(
      Gene[Master_regulator] %in%
        Gene[available_for_module_score]
    ),
    .groups = "drop"
  ) %>%
  mutate(
    module_score_eligible =
      n_genes_available >= minimum_module_genes_for_score
  ) %>%
  arrange(Module_group, Immune_module)

immune_gene_module_labels <- immune_module_gene_map %>%
  group_by(Gene) %>%
  summarise(
    Immune_modules = paste(
      sort(unique(Immune_module)),
      collapse = "; "
    ),
    Immune_module_groups = paste(
      sort(unique(Module_group)),
      collapse = "; "
    ),
    .groups = "drop"
  )

candidate_all_methods <- candidate_all_methods %>%
  left_join(immune_gene_module_labels, by = "Gene") %>%
  mutate(
    Immune_modules = coalesce(
      Immune_modules,
      "Not in curated immune module"
    ),
    Immune_module_groups = coalesce(
      Immune_module_groups,
      "Not in curated immune module"
    )
  )

candidate_rank_based <- candidate_rank_based %>%
  left_join(immune_gene_module_labels, by = "Gene") %>%
  mutate(
    Immune_modules = coalesce(
      Immune_modules,
      "Not in curated immune module"
    ),
    Immune_module_groups = coalesce(
      Immune_module_groups,
      "Not in curated immune module"
    )
  )

module_samples <- Reduce(
  intersect,
  list(
    rownames(col_use),
    rownames(vst_deg),
    rownames(met_use_log)
  )
)

if (length(module_samples) < 8) {
  stop(
    "Too few common samples for immune module-score analysis.",
    call. = FALSE
  )
}

immune_module_score_matrix <- data.frame(
  row.names = module_samples
)

for (module_i in names(immune_gene_sets)) {
  module_genes <- intersect(
    immune_gene_sets[[module_i]],
    colnames(vst_deg)
  )
  
  module_genes <- module_genes[
    module_genes %in%
      immune_module_gene_coverage$Gene[
        immune_module_gene_coverage$available_for_module_score
      ]
  ]
  
  if (length(module_genes) < minimum_module_genes_for_score) {
    next
  }
  
  module_expression <- as.matrix(
    vst_deg[module_samples, module_genes, drop = FALSE]
  )
  
  module_expression_z <- scale(module_expression)
  module_expression_z[!is.finite(module_expression_z)] <- NA_real_
  
  module_score <- rowMeans(
    module_expression_z,
    na.rm = TRUE
  )
  module_score[!is.finite(module_score)] <- NA_real_
  
  immune_module_score_matrix[[module_i]] <- module_score
}

if (ncol(immune_module_score_matrix) == 0) {
  stop(
    "No immune module had enough expression-QC-passing genes.",
    call. = FALSE
  )
}

covar_cols <- intersect(
  c("Age", "Sex", "BMI"),
  colnames(col_use)
)

covar_cols <- covar_cols[
  vapply(
    covar_cols,
    function(covariate_i) {
      sum(!is.na(col_use[module_samples, covariate_i])) >= 8 &&
        length(
          unique(
            na.omit(col_use[module_samples, covariate_i])
          )
        ) > 1
    },
    logical(1)
  )
]

module_correlation_rows <- vector(
  "list",
  ncol(immune_module_score_matrix) * ncol(met_use_log)
)
row_i <- 1

for (module_i in colnames(immune_module_score_matrix)) {
  for (metabolite_i in colnames(met_use_log)) {
    result_i <- calculate_module_correlations(
      x = as.numeric(
        met_use_log[module_samples, metabolite_i]
      ),
      y = as.numeric(
        immune_module_score_matrix[
          module_samples,
          module_i
        ]
      ),
      metadata = col_use[module_samples, , drop = FALSE],
      covariates = covar_cols
    )
    
    module_correlation_rows[[row_i]] <- tibble(
      analysis_set = analysis_set,
      Module_group = unname(immune_module_groups[module_i]),
      Immune_module = module_i,
      Metabolite = metabolite_i,
      n = result_i$n,
      pearson_r = result_i$pearson_r,
      pearson_p = result_i$pearson_p,
      spearman_rho = result_i$spearman_rho,
      spearman_p = result_i$spearman_p,
      partial_spearman_rho = result_i$partial_spearman_rho,
      partial_spearman_p = result_i$partial_spearman_p,
      partial_covariates = result_i$partial_covariates
    )
    
    row_i <- row_i + 1
  }
}

immune_module_score_correlations <- bind_rows(
  module_correlation_rows[seq_len(row_i - 1)]
) %>%
  mutate(
    pearson_FDR = p.adjust(pearson_p, method = "BH"),
    spearman_FDR = p.adjust(spearman_p, method = "BH"),
    partial_spearman_FDR = p.adjust(
      partial_spearman_p,
      method = "BH"
    ),
    module_consensus_r = apply(
      cbind(
        pearson_r,
        spearman_rho,
        partial_spearman_rho
      ),
      1,
      function(z) {
        z <- z[is.finite(z)]
        if (length(z) == 0) NA_real_ else median(z)
      }
    ),
    module_direction_consistent = apply(
      cbind(
        pearson_r,
        spearman_rho,
        partial_spearman_rho
      ),
      1,
      function(z) {
        z <- z[is.finite(z) & z != 0]
        length(z) == 3 && length(unique(sign(z))) == 1
      }
    ),
    all_three_methods_available =
      is.finite(pearson_r) &
      is.finite(spearman_rho) &
      is.finite(partial_spearman_rho) &
      is.finite(pearson_p) &
      is.finite(spearman_p) &
      is.finite(partial_spearman_p),
    geometric_mean_method_p = ifelse(
      all_three_methods_available,
      mapply(
        function(p1, p2, p3) {
          geometric_mean_p(c(p1, p2, p3))
        },
        pearson_p,
        spearman_p,
        partial_spearman_p
      ),
      NA_real_
    ),
    maximum_abs_method_correlation = pmax(
      abs(pearson_r),
      abs(spearman_rho),
      abs(partial_spearman_rho),
      na.rm = TRUE
    ),
    minimum_method_p = pmin(
      pearson_p,
      spearman_p,
      partial_spearman_p,
      na.rm = TRUE
    ),
    maximum_abs_method_correlation = ifelse(
      is.infinite(maximum_abs_method_correlation),
      NA_real_,
      maximum_abs_method_correlation
    ),
    minimum_method_p = ifelse(
      is.infinite(minimum_method_p),
      NA_real_,
      minimum_method_p
    ),
    Metabolite_superclass = classify_metabolite_superclass(
      Metabolite
    )
  )

# Gene-level results are used only as supporting evidence, not to construct
# the module-level correlation coefficient.
immune_module_all_pairs <- candidate_input %>%
  inner_join(immune_module_gene_map, by = "Gene") %>%
  filter(
    Gene %in% colnames(vst_deg),
    Metabolite %in% colnames(met_use_log),
    coalesce(host_gene_selected, TRUE),
    coalesce(keep_gene_for_correlation, TRUE),
    !Gene %in% manual_gene_exclude,
    is.na(met_iqr) | met_iqr >= metabolite_iqr_cutoff,
    is.na(gene_iqr) | gene_iqr >= gene_iqr_cutoff,
    is.na(met_sd) | met_sd > 0,
    is.na(gene_sd) | gene_sd > 0,
    is.na(min_floor_fraction) |
      min_floor_fraction < met_exact_floor_fraction_cutoff,
    is.na(basal_cluster_fraction) |
      basal_cluster_fraction < met_basal_cluster_fraction_cutoff,
    is.na(gene_min_floor_fraction) |
      coalesce(host_immune_target, TRUE) |
      gene_min_floor_fraction < gene_floor_fraction_cutoff
  ) %>%
  mutate(
    average_rank_rho = case_when(
      is.finite(spearman_rho) &
        is.finite(partial_spearman_rho) ~
        (spearman_rho + partial_spearman_rho) / 2,
      is.finite(spearman_rho) ~ spearman_rho,
      is.finite(partial_spearman_rho) ~ partial_spearman_rho,
      TRUE ~ NA_real_
    ),
    rank_based_trend =
      coalesce(spearman_p < correlation_p_cutoff, FALSE) |
      coalesce(
        partial_spearman_p < correlation_p_cutoff,
        FALSE
      ) |
      coalesce(
        abs(spearman_rho) >= minimum_abs_correlation,
        FALSE
      ) |
      coalesce(
        abs(partial_spearman_rho) >=
          minimum_abs_correlation,
        FALSE
      ),
    all_method_direction_consistent = apply(
      cbind(
        pearson_r,
        spearman_rho,
        partial_spearman_rho
      ),
      1,
      function(z) {
        z <- z[is.finite(z) & z != 0]
        length(z) == 3 && length(unique(sign(z))) == 1
      }
    ),
    gene_all_three_methods_available =
      is.finite(pearson_r) &
      is.finite(spearman_rho) &
      is.finite(partial_spearman_rho) &
      is.finite(pearson_p) &
      is.finite(spearman_p) &
      is.finite(partial_spearman_p),
    gene_geometric_mean_method_p = ifelse(
      gene_all_three_methods_available,
      mapply(
        function(p1, p2, p3) {
          geometric_mean_p(c(p1, p2, p3))
        },
        pearson_p,
        spearman_p,
        partial_spearman_p
      ),
      NA_real_
    ),
    gene_minimum_method_p = pmin(
      pearson_p,
      spearman_p,
      partial_spearman_p,
      na.rm = TRUE
    ),
    gene_minimum_method_p = ifelse(
      is.infinite(gene_minimum_method_p),
      NA_real_,
      gene_minimum_method_p
    )
  ) %>%
  distinct(
    Module_group,
    Immune_module,
    Metabolite,
    Gene,
    .keep_all = TRUE
  )

immune_module_gene_trends <- immune_module_all_pairs %>%
  filter(rank_based_trend) %>%
  arrange(
    Module_group,
    Immune_module,
    Metabolite,
    gene_geometric_mean_method_p,
    desc(abs(average_rank_rho)),
    Gene
  )

immune_gene_support_summary <- immune_module_all_pairs %>%
  left_join(
    immune_module_score_correlations %>%
      select(
        Module_group,
        Immune_module,
        Metabolite,
        module_consensus_r,
        module_direction_consistent
      ),
    by = c(
      "Module_group",
      "Immune_module",
      "Metabolite"
    )
  ) %>%
  mutate(
    supports_module_direction =
      rank_based_trend &
      all_method_direction_consistent &
      is.finite(average_rank_rho) &
      is.finite(module_consensus_r) &
      sign(average_rank_rho) == sign(module_consensus_r)
  ) %>%
  group_by(Module_group, Immune_module, Metabolite) %>%
  summarise(
    n_genes_tested = n_distinct(Gene),
    n_trend_genes = n_distinct(Gene[rank_based_trend]),
    n_supporting_genes = n_distinct(
      Gene[supports_module_direction]
    ),
    trend_genes = paste(
      sort(unique(Gene[rank_based_trend])),
      collapse = "; "
    ),
    supporting_genes = paste(
      sort(unique(Gene[supports_module_direction])),
      collapse = "; "
    ),
    supporting_gene_fraction =
      n_supporting_genes / n_genes_tested,
    minimum_supporting_gene_p = ifelse(
      any(
        is.finite(
          gene_geometric_mean_method_p[
            supports_module_direction
          ]
        )
      ),
      min(
        gene_geometric_mean_method_p[
          supports_module_direction
        ],
        na.rm = TRUE
      ),
      NA_real_
    ),
    .groups = "drop"
  )

immune_module_metabolite_summary <-
  immune_module_score_correlations %>%
  left_join(
    immune_gene_support_summary,
    by = c(
      "Module_group",
      "Immune_module",
      "Metabolite"
    )
  ) %>%
  left_join(
    immune_module_coverage_summary %>%
      select(
        Module_group,
        Immune_module,
        n_module_genes,
        n_genes_available,
        coverage_fraction,
        master_regulators,
        master_regulator_available
      ),
    by = c("Module_group", "Immune_module")
  ) %>%
  mutate(
    across(
      c(n_genes_tested, n_trend_genes, n_supporting_genes),
      ~ coalesce(.x, 0L)
    ),
    supporting_gene_fraction = coalesce(
      supporting_gene_fraction,
      0
    ),
    # Coordinated trend is defined only by direction at the module-score level.
    # P values determine circle radius and method-specific symbols, not outline.
    coordinated_module_trend =
      all_three_methods_available &
      module_direction_consistent,
    pearson_symbol = case_when(
      is.finite(pearson_p) & pearson_p < 0.001 ~ sector_symbol_p001,
      is.finite(pearson_p) & pearson_p < 0.01 ~ sector_symbol_p01,
      is.finite(pearson_p) & pearson_p < 0.05 ~ sector_symbol_p05,
      TRUE ~ ""
    ),
    spearman_symbol = case_when(
      is.finite(spearman_p) & spearman_p < 0.001 ~ sector_symbol_p001,
      is.finite(spearman_p) & spearman_p < 0.01 ~ sector_symbol_p01,
      is.finite(spearman_p) & spearman_p < 0.05 ~ sector_symbol_p05,
      TRUE ~ ""
    ),
    partial_spearman_symbol = case_when(
      is.finite(partial_spearman_p) & partial_spearman_p < 0.001 ~ sector_symbol_p001,
      is.finite(partial_spearman_p) & partial_spearman_p < 0.01 ~ sector_symbol_p01,
      is.finite(partial_spearman_p) & partial_spearman_p < 0.05 ~ sector_symbol_p05,
      TRUE ~ ""
    )
  ) %>%
  arrange(
    Module_group,
    Immune_module,
    Metabolite_superclass,
    geometric_mean_method_p,
    desc(maximum_abs_method_correlation)
  )

immune_module_coordinated_metabolites <-
  immune_module_metabolite_summary %>%
  filter(coordinated_module_trend)

# One strongest supporting gene per significant module-metabolite cell is
# retained for a compact supplementary scatter panel.
immune_module_scatter_candidates <- immune_module_all_pairs %>%
  left_join(
    immune_module_metabolite_summary %>%
      select(
        Module_group,
        Immune_module,
        Metabolite,
        module_consensus_r,
        module_direction_consistent,
        coordinated_module_trend,
        geometric_mean_method_p,
        pearson_symbol,
        spearman_symbol,
        partial_spearman_symbol
      ),
    by = c(
      "Module_group",
      "Immune_module",
      "Metabolite"
    )
  ) %>%
  filter(
    coordinated_module_trend,
    is.finite(geometric_mean_method_p),
    geometric_mean_method_p < immune_support_scatter_gmean_p_cutoff,
    rank_based_trend,
    all_method_direction_consistent,
    is.finite(average_rank_rho),
    is.finite(module_consensus_r),
    sign(average_rank_rho) == sign(module_consensus_r)
  ) %>%
  group_by(Module_group, Immune_module, Metabolite) %>%
  arrange(
    gene_geometric_mean_method_p,
    desc(abs(average_rank_rho)),
    Gene,
    .by_group = TRUE
  ) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(
    geometric_mean_method_p,
    gene_geometric_mean_method_p,
    Module_group,
    Immune_module,
    Metabolite
  )

write.csv(
  immune_module_gene_coverage,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_immune_module_gene_coverage_",
    module_output_tag,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  immune_module_coverage_summary,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_immune_module_coverage_summary_",
    module_output_tag,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  cbind(
    Sample = rownames(immune_module_score_matrix),
    immune_module_score_matrix
  ),
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_immune_module_scores_",
    module_output_tag,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  immune_module_all_pairs,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_immune_module_all_tested_gene_pairs_",
    module_output_tag,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  immune_module_gene_trends,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_immune_module_gene_trends_",
    module_output_tag,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  immune_module_metabolite_summary,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_immune_module_score_correlations_",
    module_output_tag,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  immune_module_scatter_candidates,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_immune_module_scatter_candidates_",
    module_output_tag,
    ".csv"
  ),
  row.names = FALSE
)

# Display all curated modules in a fixed biologically interpretable order.
# NK is adjacent to the merged CD8 effector module. Neutrophil and MDSC are
# retained as separate, non-overlapping transcriptional programs.
plot_modules <- c(
  "NK_core",
  "CD8_effector",
  "CD8_exhaustion",
  "Th1_core",
  "Th2_core",
  "Th17_core",
  "Treg_core",
  "Neutrophil_core",
  "MDSC_core",
  "cDC1",
  "TLS_Bcell"
)

module_plot_order <- tibble(
  Immune_module = plot_modules,
  Module_group = unname(immune_module_groups[plot_modules]),
  x = seq_along(plot_modules)
) %>%
  left_join(
    immune_module_coverage_summary %>%
      select(
        Module_group,
        Immune_module,
        n_genes_available
      ),
    by = c("Module_group", "Immune_module")
  ) %>%
  mutate(
    n_genes_available = coalesce(n_genes_available, 0L),
    Module_display_label = paste0(
      str_replace_all(Immune_module, "_", " "),
      " (",
      n_genes_available,
      ")"
    )
  )

immune_module_plot_df <- immune_module_metabolite_summary %>%
  filter(
    Immune_module %in% plot_modules,
    Metabolite_superclass %in% plot_metabolite_superclasses
  ) %>%
  inner_join(
    module_plot_order %>%
      select(
        Immune_module,
        Module_group,
        x,
        Module_display_label
      ),
    by = c("Immune_module", "Module_group")
  )

p_immune_module_bubble <- NULL
p_immune_module_main <- NULL
p_immune_module_key <- NULL

if (nrow(immune_module_plot_df) > 0) {
  clustering_matrix <- immune_module_plot_df %>%
    select(
      Metabolite_superclass,
      Metabolite,
      Immune_module,
      pearson_r,
      spearman_rho,
      partial_spearman_rho
    ) %>%
    pivot_longer(
      cols = c(
        pearson_r,
        spearman_rho,
        partial_spearman_rho
      ),
      names_to = "Correlation_method",
      values_to = "Correlation_value"
    ) %>%
    unite(
      "Module_method",
      Immune_module,
      Correlation_method,
      remove = TRUE
    ) %>%
    pivot_wider(
      names_from = Module_method,
      values_from = Correlation_value,
      values_fill = 0
    )
  
  metabolite_order <- character()
  
  for (superclass_i in plot_metabolite_superclasses) {
    superclass_df <- clustering_matrix %>%
      filter(Metabolite_superclass == superclass_i)
    
    if (nrow(superclass_df) == 0) {
      next
    }
    
    if (nrow(superclass_df) == 1) {
      metabolite_order <- c(
        metabolite_order,
        superclass_df$Metabolite
      )
    } else {
      clustering_values <- as.matrix(
        superclass_df[
          ,
          setdiff(
            colnames(superclass_df),
            c("Metabolite_superclass", "Metabolite")
          ),
          drop = FALSE
        ]
      )
      clustering_values[!is.finite(clustering_values)] <- 0
      
      clustering_result <- hclust(
        dist(clustering_values),
        method = "ward.D2"
      )
      
      metabolite_order <- c(
        metabolite_order,
        superclass_df$Metabolite[clustering_result$order]
      )
    }
  }
  
  metabolite_labels <- str_replace_all(
    metabolite_order,
    "_",
    " "
  )
  
  metabolite_y <- tibble(
    Metabolite = metabolite_order,
    Metabolite_label = metabolite_labels,
    y = rev(seq_along(metabolite_order))
  )
  
  plot_base <- immune_module_plot_df %>%
    inner_join(
      metabolite_y %>%
        select(Metabolite, Metabolite_label, y),
      by = "Metabolite"
    ) %>%
    mutate(
      consensus_neglog10p = ifelse(
        is.finite(geometric_mean_method_p),
        -log10(
          pmax(
            geometric_mean_method_p,
            .Machine$double.xmin
          )
        ),
        0
      ),
      # Radius is fixed to an interpretable P-value scale. Values stronger
      # than P = 0.001 are capped to avoid overlap between adjacent cells.
      consensus_neglog10p_capped = pmin(
        consensus_neglog10p,
        3
      ),
      pie_radius = 0.10 +
        0.30 * sqrt(consensus_neglog10p_capped / 3),
      outline = ifelse(
        coordinated_module_trend,
        "Coordinated trend",
        "Single/mixed trend"
      )
    )
  
  # Angle convention used by ggforce::geom_arc_bar:
  # 0 radians is at 12 o'clock and angles increase clockwise.
  # Fixed sectors: upper-right = partial Spearman; bottom = Pearson;
  # upper-left = Spearman.
  slice_template <- tibble(
    Correlation_method = c(
      "Partial Spearman",
      "Pearson",
      "Spearman"
    ),
    start = c(0, 2 * pi / 3, 4 * pi / 3),
    end = c(2 * pi / 3, 4 * pi / 3, 2 * pi)
  )
  
  pie_df <- bind_rows(
    lapply(
      seq_len(nrow(plot_base)),
      function(i) {
        row_i <- plot_base[i, , drop = FALSE]
        
        bind_cols(
          row_i[rep(1, nrow(slice_template)), , drop = FALSE],
          slice_template
        ) %>%
          mutate(
            Correlation_value = c(
              row_i$partial_spearman_rho,
              row_i$pearson_r,
              row_i$spearman_rho
            ),
            Method_p = c(
              row_i$partial_spearman_p,
              row_i$pearson_p,
              row_i$spearman_p
            ),
            Method_symbol = c(
              row_i$partial_spearman_symbol,
              row_i$pearson_symbol,
              row_i$spearman_symbol
            ),
            r0 = 0,
            r = row_i$pie_radius,
            middle_angle = (start + end) / 2,
            symbol_radius = 0.58 * r,
            symbol_x = x + symbol_radius * sin(middle_angle),
            symbol_y = y + symbol_radius * cos(middle_angle)
          )
      }
    )
  )
  
  module_group_labels <- module_plot_order %>%
    group_by(Module_group) %>%
    summarise(
      xmid = mean(range(x)),
      .groups = "drop"
    )
  
  superclass_labels <- plot_base %>%
    distinct(Metabolite_superclass, Metabolite, y) %>%
    group_by(Metabolite_superclass) %>%
    summarise(
      ymid = mean(range(y)),
      ymin = min(y) - 0.5,
      .groups = "drop"
    )
  
  p_immune_module_main <- ggplot() +
    geom_hline(
      data = superclass_labels,
      aes(yintercept = ymin),
      color = "grey85",
      linewidth = 0.30
    ) +
    ggforce::geom_arc_bar(
      data = pie_df,
      aes(
        x0 = x,
        y0 = y,
        r0 = r0,
        r = r,
        start = start,
        end = end,
        fill = Correlation_value,
        color = outline
      ),
      linewidth = 0.45,
      alpha = 0.95
    ) +
    geom_text(
      data = pie_df %>% filter(Method_symbol != ""),
      aes(
        x = symbol_x,
        y = symbol_y,
        label = Method_symbol
      ),
      inherit.aes = FALSE,
      color = "black",
      size = 2.15,
      fontface = "bold"
    ) +
    geom_text(
      data = module_group_labels,
      aes(
        x = xmid,
        y = max(plot_base$y) + 1.15,
        label = Module_group
      ),
      fontface = "bold",
      size = 3.0
    ) +
    geom_text(
      data = superclass_labels,
      aes(
        x = min(plot_base$x) - 1.10,
        y = ymid,
        label = Metabolite_superclass
      ),
      hjust = 1,
      fontface = "bold",
      size = 3.0
    ) +
    scale_x_continuous(
      breaks = module_plot_order$x,
      labels = module_plot_order$Module_display_label,
      limits = c(
        min(plot_base$x) - 1.45,
        max(plot_base$x) + 0.65
      ),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = metabolite_y$y,
      labels = metabolite_y$Metabolite_label,
      limits = c(0.5, max(plot_base$y) + 1.55),
      expand = c(0, 0)
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Correlation"
    ) +
    scale_color_manual(
      values = c(
        "Single/mixed trend" = "grey65",
        "Coordinated trend" = "black"
      ),
      name = NULL
    ) +
    coord_fixed(clip = "off") +
    labs(
      x = NULL,
      y = NULL,
      caption = paste0(
        "Circle size represents -log10 of the descriptive geometric-mean ",
        "P value from Pearson, Spearman, and partial Spearman analyses ",
        "(capped at 3). A black outline indicates that all three correlation ",
        "coefficients have the same sign, irrespective of P value. Numbers ",
        "in parentheses are the genes available and used for ",
        "each module score. Symbols are assigned independently within each ",
        "sector: * P < 0.05, # P < 0.01, and $ P < 0.001 for that method."
      )
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 8
      ),
      axis.text.y = element_text(size = 8),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      plot.caption = element_text(
        hjust = 0,
        size = 7.5
      ),
      plot.margin = margin(8, 8, 8, 55)
    )
  
  # A dedicated key makes the fixed sector orientation explicit.
  # To make the circle-size legend visually match the actual circles in the
  # main panel, we compute a display-scaling factor that accounts for the
  # different coordinate ranges and patchwork panel widths.
  pie_plot_width <- max(
    12.5,
    5.4 + 0.68 * length(plot_modules)
  )
  
  pie_plot_height <- max(
    7.0,
    3.2 + 0.28 * n_distinct(plot_base$Metabolite)
  )
  
  patchwork_width_main <- 1
  patchwork_width_key <- 0.28
  
  main_x_range <- diff(c(
    min(plot_base$x) - 1.45,
    max(plot_base$x) + 0.65
  ))
  main_y_range <- diff(c(
    0.5,
    max(plot_base$y) + 1.55
  ))
  
  key_x_limits <- c(-1.45, 1.55)
  key_y_limits <- c(-5.95, 1.30)
  key_x_range <- diff(key_x_limits)
  key_y_range <- diff(key_y_limits)
  
  main_panel_width <- pie_plot_width * patchwork_width_main /
    (patchwork_width_main + patchwork_width_key)
  key_panel_width <- pie_plot_width * patchwork_width_key /
    (patchwork_width_main + patchwork_width_key)
  
  main_unit_scale <- min(
    main_panel_width / main_x_range,
    pie_plot_height / main_y_range
  )
  
  key_unit_scale_raw <- min(
    key_panel_width / key_x_range,
    pie_plot_height / key_y_range
  )
  
  legend_radius_scaling_factor <- main_unit_scale / key_unit_scale_raw
  
  sector_key_df <- tibble(
    start = c(0, 2 * pi / 3, 4 * pi / 3),
    end = c(2 * pi / 3, 4 * pi / 3, 2 * pi),
    fill_key = c("#F4A582", "#F7F7F7", "#92C5DE")
  )
  
  size_key_df <- tibble(
    consensus_neglog10p = c(
      -log10(0.10),
      -log10(0.05),
      -log10(0.01),
      -log10(0.001)
    ),
    p_label = c(
      "1.00  (P = 0.10)",
      "1.30  (P = 0.05)",
      "2.00  (P = 0.01)",
      "3.00  (P = 0.001)"
    ),
    y = c(-1.55, -2.30, -3.15, -4.10)
  ) %>%
    mutate(
      radius = 0.10 +
        0.30 * sqrt(
          pmin(consensus_neglog10p, 3) / 3
        ),
      radius_display = radius * legend_radius_scaling_factor
    )
  
  p_immune_module_key <- ggplot() +
    ggforce::geom_arc_bar(
      data = sector_key_df,
      aes(
        x0 = 0,
        y0 = 0,
        r0 = 0,
        r = 0.58,
        start = start,
        end = end,
        fill = fill_key
      ),
      color = "grey35",
      linewidth = 0.45,
      show.legend = FALSE
    ) +
    scale_fill_identity() +
    annotate(
      "text",
      x = 0,
      y = 1.10,
      label = "Sector orientation",
      fontface = "bold",
      size = 3.2
    ) +
    annotate(
      "text",
      x = 0,
      y = 0.87,
      label = "Clockwise from 12 o'clock",
      size = 2.6
    ) +
    annotate(
      "text",
      x = 0.93,
      y = 0.35,
      label = "Partial\nSpearman",
      hjust = 0,
      size = 2.8
    ) +
    annotate(
      "text",
      x = 0,
      y = -0.83,
      label = "Pearson",
      size = 2.8
    ) +
    annotate(
      "text",
      x = -0.93,
      y = 0.35,
      label = "Spearman",
      hjust = 1,
      lineheight = 0.9,
      size = 2.8
    ) +
    annotate(
      "text",
      x = 0,
      y = -1.05,
      label = "Circle size",
      fontface = "bold",
      size = 3.2
    ) +
    ggforce::geom_circle(
      data = size_key_df,
      aes(
        x0 = -0.48,
        y0 = y,
        r = radius_display
      ),
      fill = "white",
      color = "black",
      linewidth = 0.45
    ) +
    geom_text(
      data = size_key_df,
      aes(
        x = 0.10,
        y = y,
        label = p_label
      ),
      hjust = 0,
      size = 2.7
    ) +
    annotate(
      "text",
      x = 0,
      y = -4.72,
      label = expression(-log[10](geometric~mean~P)),
      size = 2.8
    ) +
    annotate(
      "text",
      x = 0,
      y = -5.25,
      label = "Sector-specific P symbols",
      fontface = "bold",
      size = 3.0
    ) +
    annotate(
      "text",
      x = 0,
      y = -5.63,
      label = "*  P < 0.05     #  P < 0.01     $  P < 0.001",
      size = 2.7
    ) +
    coord_fixed(
      xlim = key_x_limits,
      ylim = key_y_limits,
      clip = "off"
    ) +
    theme_void() +
    theme(
      plot.margin = margin(8, 8, 8, 8)
    )
  
  p_immune_module_bubble <- (
    p_immune_module_main |
      p_immune_module_key
  ) +
    patchwork::plot_layout(
      widths = c(1, 0.28),
      guides = "collect"
    ) &
    theme(legend.position = "right")
  
  p_immune_module_bubble
  
  ggsave(
    paste0(
      "figures/Fig5D_immune_module_metabolite_pies_",
      module_output_tag,
      ".svg"
    ),
    p_immune_module_bubble,
    width = pie_plot_width,
    height = pie_plot_height,
    device = "svg"
  )
}

p_immune_support_scatter <- NULL

if (
  draw_immune_support_scatter &&
  nrow(immune_module_scatter_candidates) > 0
) {
  immune_support_scatter_pairs <-
    immune_module_scatter_candidates %>%
    slice_head(n = immune_support_scatter_max_pairs) %>%
    mutate(
      Pair_label = paste0(
        str_replace_all(Metabolite, "_", " "),
        " | ",
        Immune_module,
        " | ",
        Gene
      ),
      Correlation_label = paste0(
        "Pearson r = ",
        sprintf("%.2f", pearson_r),
        ", p = ",
        formatC(pearson_p, format = "e", digits = 1),
        "\nSpearman rho = ",
        sprintf("%.2f", spearman_rho),
        ", p = ",
        formatC(spearman_p, format = "e", digits = 1),
        "\nPartial rho = ",
        sprintf("%.2f", partial_spearman_rho),
        ", p = ",
        formatC(
          partial_spearman_p,
          format = "e",
          digits = 1
        )
      )
    )
  
  immune_support_scatter_df <- tibble()
  
  for (i in seq_len(nrow(immune_support_scatter_pairs))) {
    immune_support_scatter_df <- bind_rows(
      immune_support_scatter_df,
      tibble(
        Sample = module_samples,
        TRG_plot = as.character(
          col_use[module_samples, "TRG_plot"]
        ),
        Pair_label = immune_support_scatter_pairs$Pair_label[i],
        Correlation_label =
          immune_support_scatter_pairs$Correlation_label[i],
        Metabolite_value = as.numeric(
          met_use_log[
            module_samples,
            immune_support_scatter_pairs$Metabolite[i]
          ]
        ),
        Gene_expression = as.numeric(
          vst_deg[
            module_samples,
            immune_support_scatter_pairs$Gene[i]
          ]
        )
      )
    )
  }
  
  immune_support_scatter_df <- immune_support_scatter_df %>%
    filter(
      is.finite(Metabolite_value),
      is.finite(Gene_expression)
    ) %>%
    mutate(
      TRG_plot = factor(
        TRG_plot,
        levels = c("pCR", "non_pCR")
      ),
      Pair_label = factor(
        Pair_label,
        levels = immune_support_scatter_pairs$Pair_label
      )
    )
  
  immune_support_label_df <- immune_support_scatter_df %>%
    group_by(Pair_label) %>%
    summarise(
      Metabolite_value = quantile(
        Metabolite_value,
        0.04,
        na.rm = TRUE
      ),
      Gene_expression = quantile(
        Gene_expression,
        0.96,
        na.rm = TRUE
      ),
      Correlation_label = first(Correlation_label),
      .groups = "drop"
    )
  
  p_immune_support_scatter <- ggplot(
    immune_support_scatter_df,
    aes(x = Metabolite_value, y = Gene_expression)
  ) +
    geom_point(
      aes(color = TRG_plot),
      size = 2.0,
      alpha = 0.85
    ) +
    geom_smooth(
      aes(group = 1),
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "grey35",
      linewidth = 0.50
    ) +
    geom_label(
      data = immune_support_label_df,
      aes(
        x = Metabolite_value,
        y = Gene_expression,
        label = Correlation_label
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.0,
      lineheight = 0.90,
      label.size = 0.12,
      fill = "white"
    ) +
    facet_wrap(
      ~ Pair_label,
      scales = "free",
      ncol = immune_support_scatter_ncol
    ) +
    scale_color_manual(
      values = group_cols,
      breaks = c("pCR", "non_pCR"),
      labels = c("pCR", "non-pCR"),
      name = NULL,
      drop = FALSE
    ) +
    labs(
      x = "Metabolite abundance, log10(x + 1e-06)",
      y = "Host gene expression, VST"
    ) +
    theme_classic(base_size = 9) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 7.5),
      legend.position = "top",
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 8)
    )
  
  ggsave(
    paste0(
      "figures/Fig5D_immune_support_scatter_panel_",
      module_output_tag,
      ".svg"
    ),
    p_immune_support_scatter,
    width = 3.6 * immune_support_scatter_ncol,
    height = max(
      4.0,
      3.0 * ceiling(
        nrow(immune_support_scatter_pairs) /
          immune_support_scatter_ncol
      )
    ),
    device = "svg"
  )
}

print(immune_module_coverage_summary)
print(immune_module_coordinated_metabolites)

#-----------------------------------------------------------------#
# 6. Create manual selection template
#-----------------------------------------------------------------#
#
# The file is created only when it does not already exist.
# Edit Include and Display_order, save the CSV, and rerun this script.
#
# Valid Include values:
#   TRUE, T, 1, YES, Y
#-----------------------------------------------------------------#

if (!file.exists(manual_pair_file)) {
  write.csv(
    candidate_rank_based %>%
      transmute(
        Include = FALSE,
        Display_order = NA_integer_,
        Metabolite,
        Gene,
        Metabolite_class,
        gene_annotation,
        Immune_module_groups,
        Immune_modules,
        review_tier,
        biological_review_score,
        spearman_rho,
        spearman_p,
        partial_spearman_rho,
        partial_spearman_p,
        rank_direction_consistent,
        adjusted_effect_retained,
        lowest_point_sensitivity,
        met_logFC,
        met_p,
        host_logFC,
        host_DESeq2_p,
        radiation_relevance,
        chemoresistance_relevance,
        ici_relevance,
        microbe_supported,
        top_microbe_support,
        Manual_note = ""
      ),
    manual_pair_file,
    row.names = FALSE
  )
  
  message(
    "Manual selection template created: ",
    manual_pair_file
  )
} else {
  
  manual_pair_update <- read.csv(
    manual_pair_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>%
    select(
      -any_of(c(
        "Immune_module_groups",
        "Immune_modules"
      ))
    ) %>%
    left_join(
      candidate_rank_based %>%
        select(
          Metabolite,
          Gene,
          Immune_module_groups,
          Immune_modules
        ),
      by = c("Metabolite", "Gene")
    )
  
  write.csv(
    manual_pair_update,
    manual_pair_file,
    row.names = FALSE
  )
  
  message(
    "Existing manual selections were preserved and immune-module columns were updated: ",
    manual_pair_file
  )
}

#-----------------------------------------------------------------#
# 7. Export candidate-review workbook
#-----------------------------------------------------------------#

wb <- createWorkbook()

addWorksheet(wb, "rank_based_candidates")
writeData(wb, "rank_based_candidates", candidate_rank_based)
freezePane(wb, "rank_based_candidates", firstRow = TRUE)

addWorksheet(wb, "all_method_candidates")
writeData(wb, "all_method_candidates", candidate_all_methods)
freezePane(wb, "all_method_candidates", firstRow = TRUE)

addWorksheet(wb, "housekeeping_controls")
writeData(wb, "housekeeping_controls", housekeeping_control_pairs)
freezePane(wb, "housekeeping_controls", firstRow = TRUE)

addWorksheet(wb, "manual_selection")
writeData(
  wb,
  "manual_selection",
  read.csv(
    manual_pair_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
)
freezePane(wb, "manual_selection", firstRow = TRUE)

addWorksheet(wb, "immune_module_genes")
writeData(wb, "immune_module_genes", immune_module_gene_coverage)
freezePane(wb, "immune_module_genes", firstRow = TRUE)

addWorksheet(wb, "immune_module_coverage")
writeData(
  wb,
  "immune_module_coverage",
  immune_module_coverage_summary
)
freezePane(wb, "immune_module_coverage", firstRow = TRUE)

addWorksheet(wb, "immune_module_scores")
writeData(
  wb,
  "immune_module_scores",
  cbind(
    Sample = rownames(immune_module_score_matrix),
    immune_module_score_matrix
  )
)
freezePane(wb, "immune_module_scores", firstRow = TRUE)

addWorksheet(wb, "module_correlations")
writeData(
  wb,
  "module_correlations",
  immune_module_metabolite_summary
)
freezePane(wb, "module_correlations", firstRow = TRUE)

addWorksheet(wb, "immune_gene_all_pairs")
writeData(wb, "immune_gene_all_pairs", immune_module_all_pairs)
freezePane(wb, "immune_gene_all_pairs", firstRow = TRUE)

addWorksheet(wb, "immune_gene_trends")
writeData(wb, "immune_gene_trends", immune_module_gene_trends)
freezePane(wb, "immune_gene_trends", firstRow = TRUE)

addWorksheet(wb, "scatter_candidates")
writeData(
  wb,
  "scatter_candidates",
  immune_module_scatter_candidates
)
freezePane(wb, "scatter_candidates", firstRow = TRUE)

addWorksheet(wb, "settings")
writeData(
  wb,
  "settings",
  tibble(
    setting = c(
      "analysis_set",
      "module_gene_set_version",
      "correlation_p_cutoff",
      "minimum_abs_correlation",
      "metabolite_iqr_cutoff",
      "gene_iqr_cutoff",
      "gene_floor_fraction_cutoff",
      "met_exact_floor_fraction_cutoff",
      "met_basal_cluster_fraction_cutoff",
      "minimum_module_genes_for_score",
      "immune_support_scatter_gmean_p_cutoff",
      "coordinated_trend_definition",
      "module_score_definition",
      "module_plot_methods",
      "significance_symbol_definition"
    ),
    value = c(
      analysis_set,
      module_gene_set_version,
      correlation_p_cutoff,
      minimum_abs_correlation,
      metabolite_iqr_cutoff,
      gene_iqr_cutoff,
      gene_floor_fraction_cutoff,
      met_exact_floor_fraction_cutoff,
      met_basal_cluster_fraction_cutoff,
      minimum_module_genes_for_score,
      immune_support_scatter_gmean_p_cutoff,
      "All three module-score correlations are available and have the same sign",
      "Mean of gene-wise z-scored VST expression",
      "Pearson; Spearman; partial Spearman",
      paste0(
        "Symbols are assigned separately to each method sector: ",
        sector_symbol_p05, " for P < 0.05; ",
        sector_symbol_p01, " for P < 0.01; ",
        sector_symbol_p001, " for P < 0.001. The circle outline is based ",
        "only on equal signs across the three methods."
      )
    )
  )
)

for (sheet_i in names(wb)) {
  first_row <- readWorkbook(
    wb,
    sheet = sheet_i,
    rows = 1
  )
  
  if (ncol(first_row) > 0) {
    addStyle(
      wb,
      sheet_i,
      style = createStyle(
        textDecoration = "bold",
        fgFill = "#E6E6E6",
        border = "Bottom"
      ),
      rows = 1,
      cols = seq_len(ncol(first_row)),
      gridExpand = TRUE,
      stack = TRUE
    )
    
    setColWidths(
      wb,
      sheet_i,
      cols = seq_len(ncol(first_row)),
      widths = "auto"
    )
  }
}

saveWorkbook(
  wb,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
    "Fig5D_candidate_review_",
    module_output_tag,
    ".xlsx"
  ),
  overwrite = TRUE
)

#-----------------------------------------------------------------#
# 8. Optional concise candidate bubble overview
#-----------------------------------------------------------------#
#
# The complete rank-based candidate list is retained in CSV/Excel. A single
# figure containing every pair is not generated because dense labels obscure
# both the text and correlation bubbles. Set draw_candidate_bubble to TRUE
# only when a compact overview of the highest-ranked candidates is needed.
#-----------------------------------------------------------------#

# Defensive defaults allow this section to run even if only this block was
# pasted into an older version of the script.
if (!exists("draw_candidate_bubble")) {
  draw_candidate_bubble <- FALSE
}
if (!exists("candidate_bubble_top_n")) {
  candidate_bubble_top_n <- 20
}

p_candidate_bubble <- NULL

if (draw_candidate_bubble) {
  
  bubble_pair_tbl <- candidate_rank_based %>%
    arrange(
      desc(biological_review_score),
      min_rank_cor_p,
      desc(max_abs_rank_cor),
      Metabolite,
      Gene
    ) %>%
    slice_head(n = candidate_bubble_top_n) %>%
    mutate(
      pair_label = paste0(
        str_replace_all(Metabolite, "_", " "),
        " | ",
        Gene
      ),
      pair_label = factor(
        pair_label,
        levels = rev(unique(pair_label))
      )
    )
  
  bubble_long <- bind_rows(
    bubble_pair_tbl %>%
      transmute(
        pair_label,
        method = "Spearman",
        rho = spearman_rho,
        p_value = spearman_p
      ),
    bubble_pair_tbl %>%
      transmute(
        pair_label,
        method = "Partial Spearman",
        rho = partial_spearman_rho,
        p_value = partial_spearman_p
      )
  ) %>%
    mutate(
      method = factor(
        method,
        levels = c("Spearman", "Partial Spearman")
      ),
      neglog10p = ifelse(
        is.finite(p_value),
        -log10(pmax(p_value, .Machine$double.xmin)),
        NA_real_
      )
    )
  
  p_candidate_bubble <- ggplot(
    bubble_long,
    aes(
      x = method,
      y = pair_label,
      color = rho,
      size = neglog10p
    )
  ) +
    geom_point(alpha = 0.90) +
    scale_color_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Correlation"
    ) +
    scale_size_continuous(
      range = c(1.8, 5.5),
      name = "-log10(p)"
    ) +
    scale_x_discrete(position = "top") +
    labs(
      x = NULL,
      y = NULL,
      title = paste0(
        "Top ",
        nrow(bubble_pair_tbl),
        " rank-based metabolite–gene candidates"
      ),
      subtitle = "All candidates are retained in the review tables"
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(face = "bold", size = 9),
      axis.text.y = element_text(size = 8),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9)
    )
  
  p_candidate_bubble
  
  ggsave(
    paste0(
      "figures/Fig5D_candidate_bubble_top_",
      nrow(bubble_pair_tbl),
      "_",
      analysis_set,
      ".svg"
    ),
    p_candidate_bubble,
    width = 7.2,
    height = max(5.0, 1.8 + 0.30 * nrow(bubble_pair_tbl)),
    device = "svg"
  )
  
} else {
  message(
    "Optional candidate bubble was intentionally skipped; ",
    "p_candidate_bubble remains NULL because draw_candidate_bubble is FALSE. ",
    "Set draw_candidate_bubble <- TRUE only to plot the top-ranked pairs."
  )
}

#-----------------------------------------------------------------#
# 9. Individual candidate scatter export omitted
#-----------------------------------------------------------------#
#
# All rank-based candidates remain available in the CSV and Excel tables.
# Scatter plots are generated only for pairs marked Include = TRUE in the
# manual selection file below.
#-----------------------------------------------------------------#

message(
  "Individual scatter export for every candidate was skipped. ",
  "Only manually selected pairs will be plotted."
)


#-----------------------------------------------------------------#
# 9B. Manual RORC–indolepropionic acid association
#-----------------------------------------------------------------#
#
# This focused plot is retained as a biological follow-up to the Th17 analysis.
# Nicotinic acid is not included. Pearson, Spearman, and Age/Sex/BMI-adjusted
# partial Spearman statistics are displayed in a single scatter plot.
#-----------------------------------------------------------------#

if (!"RORC" %in% colnames(vst_deg)) {
  stop(
    "RORC is not available in vst_deg. Check Script 1 gene selection and expression QC.",
    call. = FALSE
  )
}

ipa_name <- colnames(met_use_log)[
  tolower(colnames(met_use_log)) == "indolepropionic_acid"
]

if (length(ipa_name) != 1) {
  stop(
    "Indolepropionic acid column을 하나로 특정할 수 없습니다: ",
    paste(
      grep(
        "indole.*prop",
        colnames(met_use_log),
        ignore.case = TRUE,
        value = TRUE
      ),
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (!all(c("Age", "Sex", "BMI", "TRG_plot") %in% colnames(col_use))) {
  stop(
    "RORC–IPA partial correlation requires Age, Sex, BMI, and TRG_plot in sample metadata.",
    call. = FALSE
  )
}

common_samples <- Reduce(
  intersect,
  list(
    rownames(col_use),
    rownames(met_use_log),
    rownames(vst_deg)
  )
)

plot_IPA <- data.frame(
  Sample = common_samples,
  RORC = as.numeric(
    vst_deg[common_samples, "RORC", drop = TRUE]
  ),
  Indolepropionic_acid = as.numeric(
    met_use_log[common_samples, ipa_name, drop = TRUE]
  ),
  Age = as.numeric(
    col_use[common_samples, "Age", drop = TRUE]
  ),
  Sex = factor(
    col_use[common_samples, "Sex", drop = TRUE]
  ),
  BMI = as.numeric(
    col_use[common_samples, "BMI", drop = TRUE]
  ),
  TRG_plot = factor(
    col_use[common_samples, "TRG_plot", drop = TRUE],
    levels = names(group_cols)
  )
)

plot_IPA <- plot_IPA[
  complete.cases(plot_IPA),
  ,
  drop = FALSE
]

if (
  nrow(plot_IPA) < 8 ||
  sd(plot_IPA$Indolepropionic_acid) == 0 ||
  sd(plot_IPA$RORC) == 0
) {
  stop(
    "RORC–IPA correlation requires at least 8 complete samples and non-zero variation.",
    call. = FALSE
  )
}

pearson_IPA <- cor.test(
  plot_IPA$Indolepropionic_acid,
  plot_IPA$RORC,
  method = "pearson"
)

spearman_IPA <- cor.test(
  plot_IPA$Indolepropionic_acid,
  plot_IPA$RORC,
  method = "spearman",
  exact = FALSE
)

plot_IPA$RORC_adj <- residuals(
  lm(
    rank(RORC, ties.method = "average") ~
      rank(Age, ties.method = "average") +
      Sex +
      rank(BMI, ties.method = "average"),
    data = plot_IPA
  )
)

plot_IPA$IPA_adj <- residuals(
  lm(
    rank(Indolepropionic_acid, ties.method = "average") ~
      rank(Age, ties.method = "average") +
      Sex +
      rank(BMI, ties.method = "average"),
    data = plot_IPA
  )
)

partial_IPA <- cor.test(
  plot_IPA$IPA_adj,
  plot_IPA$RORC_adj,
  method = "pearson"
)

label_IPA <- paste0(
  "Pearson r = ",
  sprintf("%.2f", unname(pearson_IPA$estimate)),
  ", P = ",
  format.pval(
    pearson_IPA$p.value,
    digits = 2,
    eps = 0.001
  ),
  "\nSpearman \u03c1 = ",
  sprintf("%.2f", unname(spearman_IPA$estimate)),
  ", P = ",
  format.pval(
    spearman_IPA$p.value,
    digits = 2,
    eps = 0.001
  ),
  "\nPartial Spearman \u03c1 = ",
  sprintf("%.2f", unname(partial_IPA$estimate)),
  ", P = ",
  format.pval(
    partial_IPA$p.value,
    digits = 2,
    eps = 0.001
  )
)

p_RORC_IPA <- ggplot(
  plot_IPA,
  aes(
    x = Indolepropionic_acid,
    y = RORC,
    color = TRG_plot
  )
) +
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "black",
    linewidth = 0.8
  ) +
  scale_color_manual(
    values = group_cols,
    breaks = names(group_cols),
    drop = FALSE
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label_IPA,
    hjust = 1.05,
    vjust = 1.15,
    size = 3.4
  ) +
  labs(
    x = "Indolepropionic acid",
    y = "RORC expression (VST)",
    color = NULL,
    title = "RORC vs indolepropionic acid"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

p_RORC_IPA

ggsave(
  paste0(
    "figures/RORC_indolepropionic_acid_",
    analysis_set,
    ".svg"
  ),
  p_RORC_IPA,
  width = 5.2,
  height = 4.5,
  device = "svg"
)


#-----------------------------------------------------------------#
# 10. Draw manually selected scatter panels
#-----------------------------------------------------------------#
#
# This block runs only when Include is set to TRUE/T/1/YES/Y in the manual
# selection CSV. The number of selected pairs is not restricted.
#-----------------------------------------------------------------#

manual_pairs <- read.csv(
  manual_pair_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c(
  "Include", "Display_order", "Metabolite", "Gene"
) %in% colnames(manual_pairs))) {
  stop(
    "Manual selection CSV must contain Include, Display_order, Metabolite, and Gene.",
    call. = FALSE
  )
}

if (!"Manual_note" %in% colnames(manual_pairs)) {
  manual_pairs$Manual_note <- ""
}

manual_pairs <- manual_pairs %>%
  mutate(
    Include = toupper(str_trim(as.character(Include))) %in% c(
      "TRUE", "T", "1", "YES", "Y"
    ),
    Display_order = suppressWarnings(
      as.numeric(Display_order)
    )
  )

selected_pairs <- candidate_rank_based %>%
  inner_join(
    manual_pairs %>%
      filter(Include) %>%
      select(
        Metabolite,
        Gene,
        Display_order,
        Manual_note
      ),
    by = c("Metabolite", "Gene")
  ) %>%
  arrange(
    is.na(Display_order),
    Display_order,
    desc(biological_review_score),
    min_rank_cor_p,
    desc(max_abs_rank_cor)
  )

if (nrow(selected_pairs) > 0) {
  
  selected_pairs <- selected_pairs %>%
    mutate(
      pair_label = paste0(
        str_replace_all(Metabolite, "_", " "),
        "\n",
        Gene
      ),
      correlation_label = paste0(
        "rho = ",
        sprintf("%.2f", spearman_rho),
        ", p = ",
        formatC(spearman_p, format = "e", digits = 1),
        "\nadj. rho = ",
        sprintf("%.2f", partial_spearman_rho),
        ", p = ",
        formatC(
          partial_spearman_p,
          format = "e",
          digits = 1
        )
      )
    )
  
  plot_df <- tibble()
  
  for (i in seq_len(nrow(selected_pairs))) {
    
    plot_df <- bind_rows(
      plot_df,
      tibble(
        Sample = rownames(col_use),
        TRG_plot = as.character(col_use$TRG_plot),
        Metabolite = selected_pairs$Metabolite[i],
        Gene = selected_pairs$Gene[i],
        pair_label = selected_pairs$pair_label[i],
        correlation_label =
          selected_pairs$correlation_label[i],
        met_value = met_use_log[
          rownames(col_use),
          selected_pairs$Metabolite[i]
        ],
        gene_expr = vst_deg[
          rownames(col_use),
          selected_pairs$Gene[i]
        ]
      )
    )
  }
  
  plot_df <- plot_df %>%
    filter(
      is.finite(met_value),
      is.finite(gene_expr)
    ) %>%
    mutate(
      TRG_plot = str_trim(TRG_plot),
      TRG_plot = case_when(
        TRG_plot %in% c("pCR", "CR") ~ "pCR",
        TRG_plot %in% c(
          "non_pCR", "non-pCR", "nonCR", "non_CR"
        ) ~ "non_pCR",
        TRUE ~ TRG_plot
      ),
      TRG_plot = factor(
        TRG_plot,
        levels = c("pCR", "non_pCR")
      ),
      pair_label = factor(
        pair_label,
        levels = selected_pairs$pair_label
      )
    )
  
  label_df <- plot_df %>%
    group_by(pair_label) %>%
    summarise(
      met_value = quantile(
        met_value,
        0.04,
        na.rm = TRUE
      ),
      gene_expr = quantile(
        gene_expr,
        0.96,
        na.rm = TRUE
      ),
      correlation_label = first(correlation_label),
      .groups = "drop"
    )
  
  p_manual_scatter_panel <- ggplot(
    plot_df,
    aes(x = met_value, y = gene_expr)
  ) +
    geom_point(
      aes(color = TRG_plot),
      size = 2.4,
      alpha = 0.88
    ) +
    geom_smooth(
      aes(group = 1),
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "grey35",
      linewidth = 0.55
    ) +
    geom_label(
      data = label_df,
      aes(
        x = met_value,
        y = gene_expr,
        label = correlation_label
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.4,
      lineheight = 0.95,
      label.size = 0.15,
      fill = "white"
    ) +
    facet_wrap(
      ~ pair_label,
      scales = "free",
      ncol = scatter_panel_ncol
    ) +
    scale_color_manual(
      values = group_cols,
      breaks = c("pCR", "non_pCR"),
      labels = c("pCR", "non-pCR"),
      name = NULL,
      drop = FALSE
    ) +
    labs(
      x = "Metabolite abundance, log10(x + 1e-06)",
      y = "Host gene expression, VST"
    ) +
    theme_classic(base_size = 10) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        size = 9
      ),
      legend.position = "top",
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    )
  
  p_manual_scatter_panel
  
  ggsave(
    paste0(
      "figures/Fig5D_manual_scatter_panel_",
      analysis_set,
      ".svg"
    ),
    p_manual_scatter_panel,
    width = 3.7 * scatter_panel_ncol,
    height = max(
      4.2,
      3.2 * ceiling(
        nrow(selected_pairs) / scatter_panel_ncol
      )
    ),
    device = "svg"
  )
  
  p_manual_scatter_vertical <- ggplot(
    plot_df,
    aes(x = met_value, y = gene_expr)
  ) +
    geom_point(
      aes(color = TRG_plot),
      size = 2.4,
      alpha = 0.88
    ) +
    geom_smooth(
      aes(group = 1),
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "grey35",
      linewidth = 0.55
    ) +
    geom_label(
      data = label_df,
      aes(
        x = met_value,
        y = gene_expr,
        label = correlation_label
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.4,
      lineheight = 0.95,
      label.size = 0.15,
      fill = "white"
    ) +
    facet_wrap(
      ~ pair_label,
      scales = "free",
      ncol = 1
    ) +
    scale_color_manual(
      values = group_cols,
      breaks = c("pCR", "non_pCR"),
      labels = c("pCR", "non-pCR"),
      name = NULL,
      drop = FALSE
    ) +
    labs(
      x = "Metabolite abundance, log10(x + 1e-06)",
      y = "Host gene expression, VST"
    ) +
    theme_classic(base_size = 10) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        size = 9
      ),
      legend.position = "top",
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9)
    )
  
  ggsave(
    paste0(
      "figures/Fig5D_manual_scatter_vertical_",
      analysis_set,
      ".svg"
    ),
    p_manual_scatter_vertical,
    width = 4.8,
    height = max(
      4.2,
      3.1 * nrow(selected_pairs)
    ),
    device = "svg"
  )
  
  write.csv(
    selected_pairs,
    paste0(
      "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
      "Fig5D_manual_selected_pair_summary_",
      analysis_set,
      ".csv"
    ),
    row.names = FALSE
  )
  
  save(
    analysis_set,
    triad_input_object_report,
    met_use_log_results,
    met_use_raw_results,
    vst_deg_results,
    sample_metadata_results,
    species_use_log_results,
    species_use_raw_results,
    species_qc_results,
    species_met_results,
    candidate_all_methods,
    candidate_rank_based,
    housekeeping_control_pairs,
    module_gene_set_version,
    immune_gene_sets_strict,
    immune_gene_sets_broad,
    immune_gene_sets,
    immune_support_scatter_gmean_p_cutoff,
    immune_module_gene_coverage,
    immune_module_coverage_summary,
    immune_module_all_pairs,
    immune_module_gene_trends,
    immune_module_metabolite_summary,
    immune_module_coordinated_metabolites,
    p_immune_module_bubble,
    p_immune_support_scatter,
    immune_module_score_matrix,
    immune_module_score_correlations,
    immune_module_scatter_candidates,
    p_RORC_IPA,
    selected_pairs,
    plot_df,
    label_df,
    p_candidate_bubble,
    p_manual_scatter_panel,
    p_manual_scatter_vertical,
    group_cols,
    gene_annotation_cols,
    file = paste0(
      "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
      "Fig5D_candidate_review_and_manual_scatter_",
      module_output_tag,
      ".RData"
    )
  )
  
  message(
    "Manual scatter figures generated for ",
    nrow(selected_pairs),
    " selected pairs."
  )
  
} else {
  
  save(
    analysis_set,
    triad_input_object_report,
    met_use_log_results,
    met_use_raw_results,
    vst_deg_results,
    sample_metadata_results,
    species_use_log_results,
    species_use_raw_results,
    species_qc_results,
    species_met_results,
    candidate_all_methods,
    candidate_rank_based,
    housekeeping_control_pairs,
    module_gene_set_version,
    immune_gene_sets_strict,
    immune_gene_sets_broad,
    immune_gene_sets,
    immune_support_scatter_gmean_p_cutoff,
    immune_module_gene_coverage,
    immune_module_coverage_summary,
    immune_module_all_pairs,
    immune_module_gene_trends,
    immune_module_metabolite_summary,
    immune_module_coordinated_metabolites,
    p_immune_module_bubble,
    p_immune_support_scatter,
    immune_module_score_matrix,
    immune_module_score_correlations,
    immune_module_scatter_candidates,
    p_RORC_IPA,
    p_candidate_bubble,
    group_cols,
    gene_annotation_cols,
    file = paste0(
      "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
      "Fig5D_candidate_review_",
      module_output_tag,
      ".RData"
    )
  )
  
  message(
    "No pair is marked Include = TRUE yet. ",
    "Edit the manual selection CSV and rerun this script."
  )
}

message(
  "Done: module-score screening, tri-slice correlation plot, and supplementary scatter export completed."
)
