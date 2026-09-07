#-----------------------------------------------------------------#
#
# Fig. 5D. Block-based microbe–metabolite–host gene axis screening
#
# Purpose:
#   This script screens candidate axes in the following form:
#
#       microbe  --  metabolite  --  host gene expression
#
#   The script is intentionally divided into executable blocks so that
#   each intermediate object can be inspected before moving to the next
#   step. You can run one block at a time in RStudio.
#
# Main analysis logic:
#   1. Match RNA-seq, metabolite, metadata, and species abundance samples.
#   2. Remove metabolites with low variation or exact/near-minimum basal clustering.
#   3. Select generous host DEGs plus a comprehensive curated immune panel.
#   4. Correlate metabolites with host gene expression.
#   5. Correlate microbes with metabolites.
#   6. Save matched species, metabolite, host-expression, and metadata matrices.
#   7. Join species–metabolite and metabolite–host results into indirect axes.
#   8. Rank preliminary two-edge axes; the downstream triad script recalculates
#      all three edges using both Spearman and partial Spearman.
#   9. Draw diagnostic plots and final heatmaps.
#
# Important note about many-to-many joins:
#   The triad-building step intentionally joins species–metabolite pairs
#   with metabolite–gene pairs using:
#
#       by = c("analysis_set", "Metabolite")
#
#   If one metabolite is linked to multiple microbes and multiple host genes,
#   this creates all candidate species–metabolite–gene combinations. This is
#   expected for triad screening, but the degree of row expansion is diagnosed
#   and saved in:
#
#       diagnostics/Fig5D_join_key_expansion_<analysis_set>.csv
#       diagnostics/Fig5D_join_summary_<analysis_set>.csv
#
# Output folder:
#   host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/
#
# Major output files:
#   tables/Fig5D_metabolite_qc_<analysis_set>.csv
#     - One row per metabolite.
#     - Contains finite sample count, SD, IQR, floor fraction, and QC decision.
#
#   tables/Fig5D_host_DEG_generous_<analysis_set>.csv
#     - One row per selected host gene.
#     - Contains DESeq2 statistics, Wilcoxon p-value, combined minimum p-value,
#       and host direction.
#
#   tables/Fig5D_metabolite_host_all_pairs_<analysis_set>.csv
#     - One row per metabolite–host gene pair.
#     - Contains Pearson, Spearman, and partial Spearman correlations.
#
#   tables/Fig5D_metabolite_host_screened_pairs_<analysis_set>.csv
#     - Filtered metabolite–host gene pairs with at least weak trend:
#       Pearson p < 0.1 OR Spearman p < 0.1 OR partial Spearman p < 0.1.
#
#   tables/Fig5D_species_metabolite_pairs_<analysis_set>.csv
#     - Species–metabolite Spearman pairs with p < 0.1.
#
#   tables/Fig5D_candidate_indirect_axes_<analysis_set>.csv
#     - One row per candidate species–metabolite–host gene axis.
#     - Uses only the species–metabolite and metabolite–host edges.
#     - Direct species–host correlations are intentionally not recalculated
#       because they were already evaluated in Panel C.
#
#   Fig5D_microbe_metabolite_host_axis_screening_results.RData
#     - Essential R objects for downstream analyses and visualization.
#
# Suggested use:
#   1. Run Blocks 0–1 once.
#   2. Set analysis_set in Block 2 to either "before" or "all_available".
#   3. Run Blocks 2–7 stepwise, inspecting outputs after each block.
#   4. Run Block 8 for figures.
#   5. Change analysis_set and repeat Blocks 2–8 if needed.
#   6. Run Block 9 to save the combined RData object.
#
#-----------------------------------------------------------------#


# Block 0: package, folder 생성
# Block 1: input loading 및 입력 객체 점검
# Block 2: analysis_set 선택 후 RNA-seq / metabolite / metadata sample matching
# Block 3: metabolite QC 및 QC 그림 저장
# Block 4: generous DEG + curated immune-gene selection 및 QC
# Block 5: metabolite–host gene correlation
# Block 6: species–metabolite correlation only
# Block 7: many-to-many join 진단 및 indirect axis assembly
# Block 8: 최종 heatmap visualization
# Block 9: essential objects 저장
# Optional Block 10: 나중에 RData reload용


#=================================================================#
# Block 0. Fresh session, packages, folders
#=================================================================#

rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("D:/2-연구/2-CRC metagenomics/")

dir.create(
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/diagnostics",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/intermediate_RData",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

for (pkg in c(
  "dplyr", "tidyr", "tibble", "stringr", "ggplot2",
  "svglite", "circlize"
)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, type = "binary")
  }
}

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", type = "binary")
  }
  BiocManager::install("ComplexHeatmap", update = FALSE, ask = FALSE)
}


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(svglite)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
})


#-----------------------------------------------------------------#
# Colors
#-----------------------------------------------------------------#

if (!exists("group_cols")) {
  group_cols <- c(
    "pCR" = "#5AB4AC",
    "non_pCR" = "#D97A6C"
  )
}


#=================================================================#
# Block 1. Load inputs and inspect available objects
#=================================================================#

load("input/coherence_data.RData")

load(
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5_host_RNAseq_inputs.RData"
)

load(
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5B_C_downstream_inputs.RData"
)

input_object_report <- tibble(
  object_name = c(
    "coherence_data",
    "coherence_data$species_raw",
    "vst_int",
    "vst_rna",
    "col_int",
    "met_int",
    "met_rna",
    "deg_rna",
    "group_cols"
  ),
  exists_in_workspace = c(
    exists("coherence_data"),
    exists("coherence_data") && !is.null(coherence_data$species_raw),
    exists("vst_int"),
    exists("vst_rna"),
    exists("col_int"),
    exists("met_int"),
    exists("met_rna"),
    exists("deg_rna"),
    exists("group_cols")
  ),
  object_class = NA_character_,
  n_row = NA_integer_,
  n_col = NA_integer_
)

for (i in seq_len(nrow(input_object_report))) {
  object_name_i <- input_object_report$object_name[i]
  
  if (object_name_i == "coherence_data$species_raw") {
    if (exists("coherence_data") && !is.null(coherence_data$species_raw)) {
      input_object_report$object_class[i] <- paste(class(coherence_data$species_raw), collapse = ";")
      input_object_report$n_row[i] <- nrow(coherence_data$species_raw)
      input_object_report$n_col[i] <- ncol(coherence_data$species_raw)
    }
  } else if (exists(object_name_i)) {
    object_i <- get(object_name_i)
    input_object_report$object_class[i] <- paste(class(object_i), collapse = ";")
    if (!is.null(dim(object_i))) {
      input_object_report$n_row[i] <- dim(object_i)[1]
      input_object_report$n_col[i] <- dim(object_i)[2]
    }
  }
}

write.csv(
  input_object_report,
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/diagnostics/Fig5D_input_object_report.csv",
  row.names = FALSE
)

print(input_object_report)

# Result containers.
# These are lists indexed by analysis_set, e.g. met_gene_results[["before"]].
met_gene_results <- list()
met_gene_all_results <- list()
met_qc_results <- list()
host_deg_results <- list()
triad_results <- list()
sample_metadata_results <- list()
species_qc_results <- list()
species_met_results <- list()
species_gene_results <- list()
species_gene_target_results <- list()
join_diagnostics_results <- list()
analysis_summary_results <- list()
host_expr_qc_results <- list()
host_gene_screen_results <- list()
met_use_log_results <- list()
met_use_raw_results <- list()
vst_deg_results <- list()
species_use_log_results <- list()
species_use_raw_results <- list()
species_met_all_results <- list()


#=================================================================#
# Block 2. Choose one analysis set and match RNA-seq/metabolite data
#=================================================================#
#
# Run this block once for one analysis_set.
# After completing Blocks 2–8, change analysis_set and repeat if needed.
#
# Available options:
#   "before"        : pre-treatment / before TNT samples only
#   "all_available" : all samples available after matching
#
#-----------------------------------------------------------------#

analysis_set <- "before"
# analysis_set <- "all_available"

# Host-gene screening options:
#   "DEG_plus_immune": generous DEGs plus all curated immune-module genes
#                      that pass expression QC. Recommended.
#   "all_expressed"  : every host gene that passes expression QC.
#   "DEG_only"       : reproduce the previous DEG-restricted analysis.
host_gene_screen_mode <- "DEG_plus_immune"

# Standard genes use the ordinary expression-QC cutoff. Curated immune genes
# are retained with a more permissive cutoff because low-abundance cytokines,
# lineage regulators, and checkpoint genes may have zero IQR or many values
# near the lower expression floor. No DEG p-value cutoff is applied to them.
gene_floor_fraction_cutoff <- 0.80
curated_gene_floor_fraction_cutoff <- 0.95
curated_gene_min_distinct_values <- 3

# Metabolite QC additionally identifies a near-minimum basal cluster.
# A cutoff of 0.50 requires more than half of samples to lie above the basal
# window before a metabolite is treated as a continuous feature.
met_exact_floor_fraction_cutoff <- 0.65
met_basal_window_fraction_of_iqr <- 0.05
met_basal_cluster_fraction_cutoff <- 0.50

# Species are used only for species–metabolite analysis in this script.
# Direct species–host gene correlations are not recalculated.
species_top_n <- 300

# A small housekeeping panel is retained as a technical/global-expression
# control in the metabolite–host table, but not treated as an immune module.
housekeeping_control_genes <- c(
  "ACTB", "GAPDH", "HPRT1", "PPIA", "RPLP0"
)

immune_target_genes <- unique(c(
  "CD3D", "CD3E", "TRAC", "CD8A", "CD8B", "PRF1", "GZMA", "GZMB",
  "GZMH", "NKG7", "CTSW", "CCL5", "IFNG", "TNF", "IL2", "LTA",
  "CCL3", "CCL4", "XCL1", "XCL2", "CD69", "CD38", "TNFRSF9", "IL2RA",
  "PDCD1", "LAG3", "HAVCR2", "TIGIT", "ENTPD1", "TOX", "CXCL13", "LAYN",
  "TCF7", "SLAMF6", "CXCR5", "IL7R", "CCR7", "LEF1", "ITGAE", "CXCR6",
  "ITGA1", "ZNF683", "RUNX3", "MKI67", "TOP2A", "STMN1", "TYMS", "PCNA",
  "CDK1", "CCNB1", "CCNB2", "FASLG", "TNFSF10", "FAS", "TNFRSF10A", "TNFRSF10B",
  "TNFRSF1A", "CD274", "PDCD1LG2", "CD80", "CD86", "PVR", "NECTIN2", "LGALS9",
  "FGL1", "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "CXCR3", "FOXP3",
  "CTLA4", "IKZF2", "TNFRSF18", "CCR8", "ICOS", "IL10", "TGFB1", "EBI3",
  "IL12A", "NT5E", "LRRC32", "FGL2", "OLR1", "S100A8", "S100A9", "FCGR3B",
  "CSF3R", "CXCR2", "ARG1", "CEACAM8", "CD14", "CCR2", "VCAN", "FCN1",
  "LILRB1", "IL1B", "IDO1", "PTGS2", "CYBB", "NCF1", "NCF2", "NCF4",
  "RORC", "IL17A", "IL17F", "CCR6", "KLRB1", "IL23R", "RORA", "CCL20",
  "IL1R1", "CSF2", "TBX21", "BHLHE40", "IL17RA", "IL17RC", "TRAF3IP2", "NFKBIZ",
  "CXCL1", "CXCL2", "CXCL5", "CXCL6", "CXCL8", "CSF3", "KLRD1", "NCR1",
  "KLRK1", "FCGR3A", "TYROBP", "FCER1G", "GNLY", "CLEC9A", "XCR1", "BATF3",
  "IRF8", "WDFY4", "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2",
  "TAPBP", "NLRC5", "PSMB8", "PSMB9", "MS4A1", "CD79A", "CD79B", "CD74",
  "CCL19", "CCL21", "MZB1", "JCHAIN", "IGKC", "C1QA", "C1QB", "C1QC",
  "APOE", "TREM2", "SPP1", "CD163", "MRC1", "MARCO", "CSF1R",
  "GATA3", "IL4", "IL5", "IL13", "PTGDR2", "CEBPB", "STAT3",
  "KLRC1", "KLRC2", "KLRF1", "FGFBP2", "CD244", "CD160", "TOX2", "NR4A1", "NR4A2", "NR4A3", "EOMES", "PRDM1", "IL12RB2", "IL18R1", "IL4R", "STAT6", "IL21", "IL22", "IL26", "AHR", "FPR1", "FPR2", "CXCR1", "S100A12", "MMP8", "MMP9", "LTF", "CAMP", "LCN2", "OLFM4", "MNDA", "SELL", "SLC7A2", "CADM1", "THBD", "DNASE1L3", "SNX22", "CPVL", "CLNK", "BCL6", "IL21R", "CD37", "BANK1"
))

message("Block 2: matching samples for analysis_set = ", analysis_set)

if (analysis_set == "all_available" && exists("vst_rna")) {
  vst_use <- vst_rna
} else if (exists("vst_int")) {
  vst_use <- vst_int
} else if (exists("vst_rna")) {
  vst_use <- vst_rna
} else {
  stop("Neither vst_int nor vst_rna exists.", call. = FALSE)
}

vst_use <- as.matrix(vst_use)

storage.mode(vst_use) <- "numeric"

if (!exists("col_int")) {
  stop("col_int does not exist.", call. = FALSE)
}

col_use <- as.data.frame(col_int, check.names = FALSE)

if ("RNA_sample_id" %in% colnames(col_use)) {
  rownames(col_use) <- col_use$RNA_sample_id
}

if (!"TRG_plot" %in% colnames(col_use)) {
  stop("TRG_plot column is required in col_int.", call. = FALSE)
}

col_use <- col_use[
  rownames(col_use) %in% colnames(vst_use) &
    col_use$TRG_plot %in% c("pCR", "non_pCR"),
  ,
  drop = FALSE
]

if (analysis_set == "before") {
  if (!all(c("TNT", "timepoint") %in% colnames(col_use))) {
    stop("For analysis_set = 'before', TNT and timepoint columns are required.", call. = FALSE)
  }
  
  col_use <- col_use[
    col_use$TNT == "Before" &
      col_use$timepoint == "pre",
    ,
    drop = FALSE
  ]
}

col_use <- col_use[!duplicated(rownames(col_use)), , drop = FALSE]

if (nrow(col_use) < 8) {
  stop(
    analysis_set,
    ": too few RNA-seq samples after metadata matching. Check sample IDs.",
    call. = FALSE
  )
}

vst_use <- vst_use[, rownames(col_use), drop = FALSE]

# Metabolite matrix priority:
#   1. met_rna
#   2. met_int
#   3. numeric columns embedded in col_int
if (exists("met_rna")) {
  met_use <- as.data.frame(met_rna, check.names = FALSE)
  metabolite_source <- "met_rna"
} else if (exists("met_int")) {
  met_use <- as.data.frame(met_int, check.names = FALSE)
  metabolite_source <- "met_int"
} else {
  met_use <- col_use
  metabolite_source <- "numeric columns in col_int"
}

if (sum(rownames(met_use) %in% rownames(col_use)) >= 8) {
  col_use <- col_use[
    rownames(col_use) %in% rownames(met_use),
    ,
    drop = FALSE
  ]
  met_use <- met_use[rownames(col_use), , drop = FALSE]
} else if (
  "SampleID" %in% colnames(col_use) &&
  sum(col_use$SampleID %in% rownames(met_use), na.rm = TRUE) >= 8
) {
  col_use <- col_use[
    !is.na(col_use$SampleID) &
      col_use$SampleID %in% rownames(met_use),
    ,
    drop = FALSE
  ]
  met_use <- met_use[col_use$SampleID, , drop = FALSE]
  rownames(met_use) <- rownames(col_use)
} else {
  stop(
    analysis_set,
    ": metabolite matrix could not be matched to RNA-seq metadata.",
    call. = FALSE
  )
}

col_use <- col_use[
  rownames(col_use) %in% colnames(vst_use),
  ,
  drop = FALSE
]

vst_use <- vst_use[, rownames(col_use), drop = FALSE]
met_use <- met_use[rownames(col_use), , drop = FALSE]

# Keep numeric metabolite-like columns only.
# The exclusion list removes clinical variables and diversity indices that may
# appear in the same metadata table.
met_use <- met_use[
  ,
  setdiff(
    colnames(met_use)[vapply(met_use, is.numeric, logical(1))],
    c(
      "Sample_weight_g", "DW_ul",
      "Age", "Height", "Weight", "BMI", "TRG_score",
      "ASA", "Pre_Op_Tstage", "Pre_Op_Nstage", "AJCCstage",
      "CEA", "ObservedStrain", "Shannon", "InvSimpson",
      "ObservedStrain.x", "Shannon.x", "InvSimpson.x",
      "ObservedStrain.y", "Shannon.y", "InvSimpson.y",
      "Chart_numb"
    )
  ),
  drop = FALSE
]

if (ncol(met_use) == 0) {
  stop(analysis_set, ": no numeric metabolite columns remained.", call. = FALSE)
}

met_use[] <- lapply(met_use, as.numeric)

sample_summary <- tibble(
  analysis_set = analysis_set,
  metabolite_source = metabolite_source,
  n_matched_samples = nrow(col_use),
  n_pCR = sum(col_use$TRG_plot == "pCR", na.rm = TRUE),
  n_non_pCR = sum(col_use$TRG_plot == "non_pCR", na.rm = TRUE),
  n_host_genes_in_vst = nrow(vst_use),
  n_metabolite_columns_before_QC = ncol(met_use)
)

sample_metadata_results[[analysis_set]] <- col_use
analysis_summary_results[[analysis_set]] <- sample_summary

# write.csv(
#   sample_summary,
#   paste0(
#     "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/diagnostics/Fig5D_sample_matching_summary_",
#     analysis_set,
#     ".csv"
#   ),
#   row.names = FALSE
# )

# save(
#   analysis_set,
#   col_use,
#   vst_use,
#   met_use,
#   sample_summary,
#   file = paste0(
#     "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/intermediate_RData/Fig5D_block02_matched_inputs_",
#     analysis_set,
#     ".RData"
#   )
# )

print(sample_summary)


#=================================================================#
# Block 3. Metabolite QC and optional QC plots
#=================================================================#

message("Block 3: metabolite QC for analysis_set = ", analysis_set)

met_qc <- tibble(
  Metabolite = colnames(met_use),
  n_finite = NA_integer_,
  n_distinct_rounded = NA_integer_,
  met_sd = NA_real_,
  met_iqr = NA_real_,
  min_value = NA_real_,
  min_floor_fraction = NA_real_,
  basal_window_upper = NA_real_,
  basal_cluster_fraction = NA_real_,
  keep_metabolite = NA
)

for (i in seq_len(nrow(met_qc))) {
  met_x <- met_use[, met_qc$Metabolite[i]]
  met_x <- met_x[is.finite(met_x)]
  
  met_qc$n_finite[i] <- length(met_x)
  met_qc$n_distinct_rounded[i] <- length(unique(round(met_x, 8)))
  met_qc$met_sd[i] <- sd(met_x, na.rm = TRUE)
  met_qc$met_iqr[i] <- IQR(met_x, na.rm = TRUE)
  met_qc$min_value[i] <- min(met_x, na.rm = TRUE)
  met_qc$min_floor_fraction[i] <- mean(
    round(met_x, 8) == round(min(met_x, na.rm = TRUE), 8),
    na.rm = TRUE
  )
  
  # A positive technical floor may contain several slightly different values.
  # The basal cluster therefore includes values within 5% of the feature IQR
  # above the observed minimum, rather than requiring exact equality.
  met_qc$basal_window_upper[i] <- met_qc$min_value[i] +
    max(
      1e-08,
      met_basal_window_fraction_of_iqr * met_qc$met_iqr[i]
    )
  
  met_qc$basal_cluster_fraction[i] <- mean(
    met_x <= met_qc$basal_window_upper[i],
    na.rm = TRUE
  )
  
  met_qc$keep_metabolite[i] <-
    met_qc$n_finite[i] >= 8 &&
    met_qc$n_distinct_rounded[i] >= 4 &&
    is.finite(met_qc$met_sd[i]) &&
    met_qc$met_sd[i] > 0 &&
    is.finite(met_qc$met_iqr[i]) &&
    met_qc$met_iqr[i] > 0 &&
    met_qc$min_floor_fraction[i] < met_exact_floor_fraction_cutoff &&
    met_qc$basal_cluster_fraction[i] < met_basal_cluster_fraction_cutoff
}

met_qc_results[[analysis_set]] <- met_qc

# write.csv(
#   met_qc,
#   paste0(
#     "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/Fig5D_metabolite_qc_",
#     analysis_set,
#     ".csv"
#   ),
#   row.names = FALSE
# )

met_use <- met_use[
  ,
  met_qc$Metabolite[met_qc$keep_metabolite],
  drop = FALSE
]

if (ncol(met_use) == 0) {
  stop(analysis_set, ": all metabolites were removed by QC.", call. = FALSE)
}

met_use_raw_results[[analysis_set]] <- as.data.frame(
  met_use,
  check.names = FALSE
)

met_use_log <- log10(met_use + 1e-06)
met_use_log_results[[analysis_set]] <- met_use_log

met_qc_summary <- met_qc %>%
  summarise(
    analysis_set = analysis_set,
    n_metabolites_before_QC = n(),
    n_metabolites_after_QC = sum(keep_metabolite),
    n_removed_by_QC = sum(!keep_metabolite),
    n_exact_floor_dominated = sum(
      min_floor_fraction >= met_exact_floor_fraction_cutoff,
      na.rm = TRUE
    ),
    n_basal_cluster_dominated = sum(
      basal_cluster_fraction >= met_basal_cluster_fraction_cutoff,
      na.rm = TRUE
    ),
    median_exact_floor_fraction = median(
      min_floor_fraction,
      na.rm = TRUE
    ),
    median_basal_cluster_fraction = median(
      basal_cluster_fraction,
      na.rm = TRUE
    )
  )

# write.csv(
#   met_qc_summary,
#   paste0(
#     "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/diagnostics/Fig5D_metabolite_qc_summary_",
#     analysis_set,
#     ".csv"
#   ),
#   row.names = FALSE
# )

print(met_qc_summary)

# QC plot 1: floor-dominated metabolites.
p_met_qc_floor <- ggplot(
  met_qc,
  aes(x = basal_cluster_fraction, fill = keep_metabolite)
) +
  geom_histogram(bins = 30, color = "white", linewidth = 0.2) +
  labs(
    title = paste0("Metabolite basal-cluster QC: ", analysis_set),
    x = "Fraction of samples near the observed minimum",
    y = "Number of metabolites",
    fill = "Kept"
  ) +
  theme_classic(base_size = 11)

p_met_qc_floor

# ggsave(
#   paste0("figures/Fig5D_metabolite_QC_floor_fraction_", analysis_set, ".svg"),
#   p_met_qc_floor,
#   width = 6,
#   height = 4,
#   device = "svg"
# )

# QC plot 2: variance/IQR behavior.
p_met_qc_variance <- ggplot(
  met_qc,
  aes(x = met_sd, y = met_iqr, color = keep_metabolite)
) +
  geom_point(size = 2, alpha = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = paste0("Metabolite variance QC: ", analysis_set),
    x = "Metabolite SD",
    y = "Metabolite IQR",
    color = "Kept"
  ) +
  theme_classic(base_size = 11)

p_met_qc_variance

# ggsave(
#   paste0("figures/Fig5D_metabolite_QC_variance_", analysis_set, ".svg"),
#   p_met_qc_variance,
#   width = 6,
#   height = 4,
#   device = "svg"
# )

# save(
#   analysis_set,
#   met_qc,
#   met_use,
#   met_use_log,
#   file = paste0(
#     "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/intermediate_RData/Fig5D_block03_metabolite_QC_",
#     analysis_set,
#     ".RData"
#   )
# )


#=================================================================#
# Block 4. Host-gene selection using expression QC
#=================================================================#
#
# The original script restricted correlation testing to genes with
# DESeq2 p < 0.1 or Wilcoxon p < 0.1 and then retained only the top 2,000.
# That design excludes biologically important immune genes that are not
# differentially expressed between pCR and non_pCR.
#
# The default "DEG_plus_immune" mode now retains:
#   1. all generous DEGs that pass expression QC, and
#   2. every curated immune-module gene that passes expression QC,
#      regardless of DEG status.
#
# Expression QC is intentionally permissive. A gene is removed only when it
# has too few observations, essentially no variation, or is strongly
# floor-dominated across samples.
#=================================================================#

message(
  "Block 4: host-gene expression QC and selection for analysis_set = ",
  analysis_set,
  "; mode = ",
  host_gene_screen_mode
)

if (!exists("deg_rna")) {
  stop("deg_rna does not exist.", call. = FALSE)
}

required_deg_cols <- c(
  "Gene",
  "log2FoldChange",
  "pval_DESeq2",
  "FDR_DESeq2",
  "direction"
)

if (!all(required_deg_cols %in% colnames(deg_rna))) {
  stop(
    "deg_rna must contain these columns: ",
    paste(required_deg_cols, collapse = ", "),
    call. = FALSE
  )
}

if (!host_gene_screen_mode %in% c(
  "DEG_plus_immune",
  "all_expressed",
  "DEG_only"
)) {
  stop(
    "host_gene_screen_mode must be DEG_plus_immune, all_expressed, or DEG_only.",
    call. = FALSE
  )
}

host_stat <- tibble(
  Gene = rownames(vst_use)
) %>%
  left_join(
    deg_rna %>%
      transmute(
        Gene,
        host_logFC = log2FoldChange,
        host_DESeq2_p = pval_DESeq2,
        host_DESeq2_FDR = FDR_DESeq2,
        host_direction = direction
      ) %>%
      distinct(Gene, .keep_all = TRUE),
    by = "Gene"
  )

host_stat$host_wilcox_p <- NA_real_

for (i in seq_len(nrow(host_stat))) {
  if (
    sum(col_use$TRG_plot == "pCR") >= 2 &&
    sum(col_use$TRG_plot == "non_pCR") >= 2
  ) {
    host_stat$host_wilcox_p[i] <- suppressWarnings(
      wilcox.test(
        vst_use[
          host_stat$Gene[i],
          col_use$TRG_plot == "pCR"
        ],
        vst_use[
          host_stat$Gene[i],
          col_use$TRG_plot == "non_pCR"
        ],
        exact = FALSE
      )$p.value
    )
  }
}

host_expr_qc <- tibble(
  Gene = rownames(vst_use),
  gene_n_finite = NA_integer_,
  gene_n_distinct_rounded = NA_integer_,
  gene_mean = NA_real_,
  gene_median = NA_real_,
  gene_sd = NA_real_,
  gene_iqr = NA_real_,
  gene_min_value = NA_real_,
  gene_min_floor_fraction = NA_real_
)

for (i in seq_len(nrow(host_expr_qc))) {
  gene_x <- as.numeric(
    vst_use[
      host_expr_qc$Gene[i],
      rownames(col_use),
      drop = TRUE
    ]
  )
  gene_x <- gene_x[is.finite(gene_x)]
  
  host_expr_qc$gene_n_finite[i] <- length(gene_x)
  host_expr_qc$gene_n_distinct_rounded[i] <- length(
    unique(round(gene_x, 6))
  )
  host_expr_qc$gene_mean[i] <- mean(gene_x, na.rm = TRUE)
  host_expr_qc$gene_median[i] <- median(gene_x, na.rm = TRUE)
  host_expr_qc$gene_sd[i] <- sd(gene_x, na.rm = TRUE)
  host_expr_qc$gene_iqr[i] <- IQR(gene_x, na.rm = TRUE)
  host_expr_qc$gene_min_value[i] <- min(gene_x, na.rm = TRUE)
  host_expr_qc$gene_min_floor_fraction[i] <- mean(
    round(gene_x, 6) ==
      round(min(gene_x, na.rm = TRUE), 6),
    na.rm = TRUE
  )
}

host_expr_qc <- host_expr_qc %>%
  mutate(
    is_curated_immune_gene = Gene %in% immune_target_genes,
    
    keep_gene_general =
      gene_n_finite >= 8 &
      gene_n_distinct_rounded >= 4 &
      is.finite(gene_sd) &
      gene_sd > 0 &
      is.finite(gene_iqr) &
      gene_iqr > 0 &
      gene_min_floor_fraction < gene_floor_fraction_cutoff,
    
    keep_gene_curated =
      gene_n_finite >= 8 &
      gene_n_distinct_rounded >= curated_gene_min_distinct_values &
      is.finite(gene_sd) &
      gene_sd > 0 &
      gene_min_floor_fraction < curated_gene_floor_fraction_cutoff,
    
    # Curated immune genes use relaxed expression QC and are never filtered
    # by DEG p-value. Constant or virtually unmeasurable genes cannot support
    # a correlation and therefore remain excluded.
    keep_gene_for_correlation = if_else(
      is_curated_immune_gene,
      keep_gene_curated,
      keep_gene_general
    ),
    
    gene_qc_reason = case_when(
      keep_gene_for_correlation & is_curated_immune_gene ~
        "kept: curated immune gene with relaxed expression QC",
      keep_gene_for_correlation ~ "kept: standard expression QC",
      gene_n_finite < 8 ~ "n_finite < 8",
      !is.finite(gene_sd) | gene_sd <= 0 ~ "SD <= 0",
      is_curated_immune_gene &
        gene_n_distinct_rounded < curated_gene_min_distinct_values ~
        "curated gene has too few distinct values",
      !is_curated_immune_gene & gene_n_distinct_rounded < 4 ~
        "< 4 distinct values",
      !is_curated_immune_gene & (!is.finite(gene_iqr) | gene_iqr <= 0) ~
        "IQR <= 0",
      is_curated_immune_gene &
        gene_min_floor_fraction >= curated_gene_floor_fraction_cutoff ~
        "curated gene is >95% floor-dominated",
      !is_curated_immune_gene &
        gene_min_floor_fraction >= gene_floor_fraction_cutoff ~
        "strongly floor-dominated",
      TRUE ~ "removed by expression QC"
    )
  )

host_stat <- host_stat %>%
  left_join(
    host_expr_qc,
    by = "Gene",
    relationship = "one-to-one"
  ) %>%
  mutate(
    host_wilcox_FDR = p.adjust(host_wilcox_p, method = "BH"),
    host_min_p = pmin(
      host_DESeq2_p,
      host_wilcox_p,
      na.rm = TRUE
    ),
    host_min_p = ifelse(
      is.infinite(host_min_p),
      NA_real_,
      host_min_p
    ),
    host_neglog10p = -log10(
      pmax(host_min_p, .Machine$double.xmin)
    ),
    
    host_DEG_selected =
      coalesce(host_DESeq2_p < 0.1, FALSE) |
      coalesce(host_wilcox_p < 0.1, FALSE),
    
    host_immune_target = Gene %in% immune_target_genes,
    host_housekeeping_control = Gene %in% housekeeping_control_genes,
    
    host_gene_selected = case_when(
      host_gene_screen_mode == "DEG_plus_immune" ~
        (host_immune_target & keep_gene_curated) |
        (host_DEG_selected & keep_gene_general) |
        (host_housekeeping_control & keep_gene_general),
      
      host_gene_screen_mode == "all_expressed" ~
        keep_gene_for_correlation,
      
      host_gene_screen_mode == "DEG_only" ~
        keep_gene_general & host_DEG_selected,
      
      TRUE ~ FALSE
    ),
    
    host_selection_source = case_when(
      host_gene_selected &
        host_housekeeping_control ~ "Housekeeping control",
      host_gene_selected &
        host_DEG_selected &
        host_immune_target ~ "DEG + curated immune module",
      host_gene_selected &
        host_DEG_selected ~ "DEG",
      host_gene_selected &
        host_immune_target ~ "Curated immune module; non-DEG",
      host_gene_selected ~ "Expression-QC-passed host gene",
      TRUE ~ "Not selected"
    )
  )

host_deg <- host_stat %>%
  filter(host_DEG_selected) %>%
  arrange(
    host_min_p,
    desc(abs(host_logFC))
  )

host_gene_screen <- host_stat %>%
  filter(host_gene_selected) %>%
  arrange(
    desc(host_immune_target),
    desc(host_DEG_selected),
    host_min_p,
    desc(abs(host_logFC))
  )

if (nrow(host_gene_screen) == 0) {
  stop(
    analysis_set,
    ": no host genes passed the selected screening mode.",
    call. = FALSE
  )
}

immune_target_coverage <- tibble(Gene = immune_target_genes) %>%
  left_join(
    host_stat %>%
      select(
        Gene,
        host_DEG_selected,
        host_gene_selected,
        keep_gene_curated,
        keep_gene_for_correlation,
        gene_n_finite,
        gene_n_distinct_rounded,
        gene_mean,
        gene_median,
        gene_sd,
        gene_iqr,
        gene_min_floor_fraction,
        gene_qc_reason
      ),
    by = "Gene"
  ) %>%
  mutate(
    present_in_vst = Gene %in% rownames(vst_use),
    included_regardless_of_DEG = coalesce(host_gene_selected, FALSE),
    inclusion_status = case_when(
      !present_in_vst ~ "Not present in VST matrix",
      coalesce(host_gene_selected, FALSE) ~
        "Included: curated immune target; DEG p-value not required",
      TRUE ~ coalesce(gene_qc_reason, "Not selected")
    )
  ) %>%
  arrange(desc(included_regardless_of_DEG), Gene)

write.csv(
  immune_target_coverage,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_curated_immune_gene_coverage_",
    analysis_set,
    ".csv"
  ),
  row.names = FALSE
)

critical_immune_markers <- c(
  "CD3D", "CD8A", "CD8B", "PRF1", "GZMB", "PDCD1",
  "FOXP3", "RORC", "TBX21", "GATA3", "IL2", "CEBPB", "STAT3"
)

message(
  "Critical immune markers included: ",
  paste(
    critical_immune_markers[
      critical_immune_markers %in% host_gene_screen$Gene
    ],
    collapse = ", "
  )
)

missing_critical_markers <- setdiff(
  critical_immune_markers,
  host_gene_screen$Gene
)

if (length(missing_critical_markers) > 0) {
  warning(
    "Critical immune markers not usable after relaxed expression QC or absent from VST: ",
    paste(missing_critical_markers, collapse = ", ")
  )
}

host_deg_results[[analysis_set]] <- host_deg
host_expr_qc_results[[analysis_set]] <- host_expr_qc
host_gene_screen_results[[analysis_set]] <- host_gene_screen

write.csv(
  host_deg,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_host_DEG_generous_",
    analysis_set,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  host_expr_qc,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_host_expression_QC_",
    analysis_set,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  host_gene_screen,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_host_gene_screen_",
    analysis_set,
    ".csv"
  ),
  row.names = FALSE
)

# Keep the historical object name vst_deg for downstream compatibility,
# although it now contains every selected host gene, not only DEGs.
vst_deg <- t(
  vst_use[
    host_gene_screen$Gene,
    rownames(col_use),
    drop = FALSE
  ]
)

vst_deg_results[[analysis_set]] <- vst_deg

host_gene_screen_summary <- tibble(
  analysis_set = analysis_set,
  host_gene_screen_mode = host_gene_screen_mode,
  gene_floor_fraction_cutoff = gene_floor_fraction_cutoff,
  curated_gene_floor_fraction_cutoff = curated_gene_floor_fraction_cutoff,
  curated_gene_min_distinct_values = curated_gene_min_distinct_values,
  n_host_genes_in_vst = nrow(vst_use),
  n_expression_QC_pass = sum(
    host_stat$keep_gene_for_correlation,
    na.rm = TRUE
  ),
  n_generous_DEGs = sum(
    host_stat$host_DEG_selected,
    na.rm = TRUE
  ),
  n_curated_immune_genes_in_vst = sum(
    host_stat$host_immune_target,
    na.rm = TRUE
  ),
  n_curated_immune_genes_QC_pass = sum(
    host_stat$host_immune_target &
      host_stat$keep_gene_curated,
    na.rm = TRUE
  ),
  n_non_DEG_immune_genes_added = sum(
    host_stat$host_gene_selected &
      host_stat$host_immune_target &
      !host_stat$host_DEG_selected,
    na.rm = TRUE
  ),
  n_housekeeping_controls_added = sum(
    host_stat$host_gene_selected &
      host_stat$host_housekeeping_control,
    na.rm = TRUE
  ),
  n_host_genes_used_for_correlation = nrow(host_gene_screen)
)

write.csv(
  host_gene_screen_summary,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/diagnostics/",
    "Fig5D_host_gene_screen_summary_",
    analysis_set,
    ".csv"
  ),
  row.names = FALSE
)

print(host_gene_screen_summary)

p_host_deg_volcano <- ggplot(
  host_stat,
  aes(x = host_logFC, y = host_neglog10p)
) +
  geom_point(alpha = 0.25, size = 1.0) +
  geom_point(
    data = host_stat %>% filter(host_DEG_selected),
    alpha = 0.65,
    size = 1.3
  ) +
  geom_point(
    data = host_stat %>%
      filter(
        host_gene_selected,
        host_immune_target,
        !host_DEG_selected
      ),
    shape = 1,
    alpha = 0.90,
    size = 1.8,
    stroke = 0.65
  ) +
  labs(
    title = paste0(
      "Host-gene screening: ",
      analysis_set
    ),
    subtitle = paste0(
      "Filled: generous DEG; open: non-DEG curated immune gene; mode = ",
      host_gene_screen_mode
    ),
    x = "Host log2 fold-change",
    y = "-log10(minimum p-value)"
  ) +
  theme_classic(base_size = 11)

p_host_deg_volcano

#=================================================================#
# Block 5. Metabolite–host gene correlation
#=================================================================#

message("Block 5: metabolite–host gene correlations for analysis_set = ", analysis_set)

covar_cols <- intersect(c("Age", "Sex", "BMI"), colnames(col_use))

if (length(covar_cols) > 0) {
  covar_cols <- covar_cols[
    vapply(
      covar_cols,
      function(covar_i) {
        sum(!is.na(col_use[[covar_i]])) >= 8 &&
          length(unique(na.omit(col_use[[covar_i]]))) > 1
      },
      logical(1)
    )
  ]
}

message(
  "Partial Spearman covariates used: ",
  ifelse(length(covar_cols) == 0, "none", paste(covar_cols, collapse = ", "))
)

message(
  "Testing ",
  ncol(vst_deg) * ncol(met_use_log),
  " metabolite–gene pairs."
)

cor_rows <- vector(
  "list",
  ncol(vst_deg) * ncol(met_use_log)
)

row_i <- 1

for (g in colnames(vst_deg)) {
  for (m in colnames(met_use_log)) {
    met_x <- met_use_log[, m]
    gene_y <- vst_deg[, g]
    ok <- is.finite(met_x) & is.finite(gene_y)
    
    pearson_r <- NA_real_
    pearson_p <- NA_real_
    spearman_rho <- NA_real_
    spearman_p <- NA_real_
    partial_spearman_rho <- NA_real_
    partial_spearman_p <- NA_real_
    partial_covariates <- NA_character_
    
    if (
      sum(ok) >= 8 &&
      sd(met_x[ok], na.rm = TRUE) > 0 &&
      sd(gene_y[ok], na.rm = TRUE) > 0
    ) {
      pearson_result <- suppressWarnings(
        cor.test(
          met_x[ok],
          gene_y[ok],
          method = "pearson"
        )
      )
      
      spearman_result <- suppressWarnings(
        cor.test(
          met_x[ok],
          gene_y[ok],
          method = "spearman",
          exact = FALSE
        )
      )
      
      pearson_r <- unname(pearson_result$estimate)
      pearson_p <- pearson_result$p.value
      spearman_rho <- unname(spearman_result$estimate)
      spearman_p <- spearman_result$p.value
      
      if (length(covar_cols) > 0) {
        partial_df <- data.frame(
          met_rank = rank(met_x[ok], ties.method = "average"),
          gene_rank = rank(gene_y[ok], ties.method = "average"),
          col_use[ok, covar_cols, drop = FALSE],
          check.names = FALSE
        )
        
        partial_df <- partial_df[complete.cases(partial_df), , drop = FALSE]
        
        if (nrow(partial_df) >= 8) {
          for (cc in setdiff(colnames(partial_df), c("met_rank", "gene_rank"))) {
            if (is.character(partial_df[[cc]]) || is.factor(partial_df[[cc]])) {
              partial_df[[cc]] <- factor(partial_df[[cc]])
            }
          }
          
          covar_cols_use <- setdiff(colnames(partial_df), c("met_rank", "gene_rank"))
          
          covar_cols_use <- covar_cols_use[
            vapply(
              covar_cols_use,
              function(z) length(unique(na.omit(partial_df[[z]]))) > 1,
              logical(1)
            )
          ]
          
          if (length(covar_cols_use) > 0 && nrow(partial_df) > length(covar_cols_use) + 3) {
            partial_spearman_rho <- suppressWarnings(
              cor(
                resid(
                  lm(
                    as.formula(
                      paste("met_rank ~", paste(covar_cols_use, collapse = " + "))
                    ),
                    data = partial_df
                  )
                ),
                resid(
                  lm(
                    as.formula(
                      paste("gene_rank ~", paste(covar_cols_use, collapse = " + "))
                    ),
                    data = partial_df
                  )
                ),
                method = "pearson"
              )
            )
            
            if (is.finite(partial_spearman_rho) && abs(partial_spearman_rho) < 1) {
              partial_spearman_p <- 2 * pt(
                -abs(
                  partial_spearman_rho *
                    sqrt(
                      (nrow(partial_df) - length(covar_cols_use) - 2) /
                        (1 - partial_spearman_rho^2)
                    )
                ),
                df = nrow(partial_df) - length(covar_cols_use) - 2
              )
            }
            
            partial_covariates <- paste(covar_cols_use, collapse = ";")
          }
        }
      }
    }
    
    cor_rows[[row_i]] <- tibble(
      analysis_set = analysis_set,
      Gene = g,
      Metabolite = m,
      n = sum(ok),
      pearson_r = pearson_r,
      pearson_p = pearson_p,
      spearman_rho = spearman_rho,
      spearman_p = spearman_p,
      partial_spearman_rho = partial_spearman_rho,
      partial_spearman_p = partial_spearman_p,
      partial_covariates = partial_covariates
    )
    
    row_i <- row_i + 1
  }
}

met_gene_pair_tbl <- bind_rows(cor_rows[seq_len(row_i - 1)]) %>%
  mutate(
    pearson_FDR = p.adjust(pearson_p, method = "BH"),
    spearman_FDR = p.adjust(spearman_p, method = "BH"),
    partial_spearman_FDR = p.adjust(partial_spearman_p, method = "BH"),
    min_cor_p = pmin(
      pearson_p,
      spearman_p,
      partial_spearman_p,
      na.rm = TRUE
    ),
    max_abs_cor = pmax(
      abs(pearson_r),
      abs(spearman_rho),
      abs(partial_spearman_rho),
      na.rm = TRUE
    ),
    min_cor_p = ifelse(is.infinite(min_cor_p), NA_real_, min_cor_p),
    max_abs_cor = ifelse(is.infinite(max_abs_cor), NA_real_, max_abs_cor),
    best_cor_method = case_when(
      !is.na(partial_spearman_p) &
        partial_spearman_p <= pmin(pearson_p, spearman_p, partial_spearman_p, na.rm = TRUE) + 1e-15 ~
        "partial_spearman",
      !is.na(spearman_p) &
        spearman_p <= pmin(pearson_p, spearman_p, partial_spearman_p, na.rm = TRUE) + 1e-15 ~
        "spearman",
      !is.na(pearson_p) ~
        "pearson",
      TRUE ~
        NA_character_
    ),
    best_cor_r = case_when(
      best_cor_method == "partial_spearman" ~ partial_spearman_rho,
      best_cor_method == "spearman" ~ spearman_rho,
      best_cor_method == "pearson" ~ pearson_r,
      TRUE ~ NA_real_
    ),
    best_cor_p = case_when(
      best_cor_method == "partial_spearman" ~ partial_spearman_p,
      best_cor_method == "spearman" ~ spearman_p,
      best_cor_method == "pearson" ~ pearson_p,
      TRUE ~ NA_real_
    )
  ) %>%
  left_join(
    host_gene_screen %>%
      select(
        Gene,
        host_logFC,
        host_DESeq2_p,
        host_DESeq2_FDR,
        host_wilcox_p,
        host_wilcox_FDR,
        host_min_p,
        host_direction,
        host_DEG_selected,
        host_immune_target,
        host_housekeeping_control,
        host_gene_selected,
        host_selection_source,
        keep_gene_general,
        keep_gene_curated,
        keep_gene_for_correlation,
        gene_n_finite,
        gene_n_distinct_rounded,
        gene_mean,
        gene_median,
        gene_sd,
        gene_iqr,
        gene_min_value,
        gene_min_floor_fraction,
        gene_qc_reason
      ),
    by = "Gene",
    relationship = "many-to-one"
  ) %>%
  left_join(
    met_qc %>%
      select(
        Metabolite,
        n_distinct_rounded,
        met_sd,
        met_iqr,
        min_value,
        min_floor_fraction,
        basal_window_upper,
        basal_cluster_fraction
      ),
    by = "Metabolite",
    relationship = "many-to-one"
  )

met_gene_all_results[[analysis_set]] <- met_gene_pair_tbl

met_gene_results[[analysis_set]] <- met_gene_pair_tbl %>%
  filter(
    coalesce(pearson_p < 0.1, FALSE) |
      coalesce(spearman_p < 0.1, FALSE) |
      coalesce(partial_spearman_p < 0.1, FALSE)
  ) %>%
  arrange(
    best_cor_p,
    desc(abs(best_cor_r)),
    host_min_p
  )

write.csv(
  met_gene_all_results[[analysis_set]],
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/Fig5D_metabolite_host_all_pairs_",
    analysis_set,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  met_gene_results[[analysis_set]],
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/Fig5D_metabolite_host_screened_pairs_",
    analysis_set,
    ".csv"
  ),
  row.names = FALSE
)

met_gene_summary <- tibble(
  analysis_set = analysis_set,
  n_all_metabolite_gene_pairs = nrow(met_gene_all_results[[analysis_set]]),
  n_screened_pairs_p_lt_0_1 = nrow(met_gene_results[[analysis_set]]),
  n_screened_pairs_p_lt_0_05 = sum(
    met_gene_all_results[[analysis_set]]$best_cor_p < 0.05,
    na.rm = TRUE
  ),
  n_unique_metabolites_in_screened_pairs = n_distinct(met_gene_results[[analysis_set]]$Metabolite),
  n_unique_genes_in_screened_pairs = n_distinct(met_gene_results[[analysis_set]]$Gene)
)

write.csv(
  met_gene_summary,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/diagnostics/Fig5D_metabolite_host_summary_",
    analysis_set,
    ".csv"
  ),
  row.names = FALSE
)

print(met_gene_summary)

if (nrow(met_gene_results[[analysis_set]]) > 0) {
  top_met_gene_examples <- met_gene_results[[analysis_set]] %>%
    slice_head(n = 6) %>%
    mutate(pair_label = paste0(Metabolite, " ~ ", Gene))
  
  plot_df_met_gene <- bind_rows(
    lapply(
      seq_len(nrow(top_met_gene_examples)),
      function(i) {
        tibble(
          Sample = rownames(col_use),
          TRG_plot = col_use$TRG_plot,
          Metabolite = top_met_gene_examples$Metabolite[i],
          Gene = top_met_gene_examples$Gene[i],
          pair_label = top_met_gene_examples$pair_label[i],
          metabolite_log10 = met_use_log[
            rownames(col_use),
            top_met_gene_examples$Metabolite[i]
          ],
          gene_expression = vst_deg[
            rownames(col_use),
            top_met_gene_examples$Gene[i]
          ]
        )
      }
    )
  )
  
  p_met_gene_examples <- ggplot(
    plot_df_met_gene,
    aes(x = metabolite_log10, y = gene_expression)
  ) +
    geom_point(aes(shape = TRG_plot), size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
    facet_wrap(~ pair_label, scales = "free", ncol = 3) +
    labs(
      title = paste0("Top metabolite–host gene pairs: ", analysis_set),
      x = "Metabolite abundance, log10(x + 1e-06)",
      y = "Host gene expression, VST"
    ) +
    theme_classic(base_size = 10)
  
  ggsave(
    paste0("figures/Fig5D_top_metabolite_host_gene_scatter_", analysis_set, ".svg"),
    p_met_gene_examples,
    width = 10,
    height = 6,
    device = "svg"
  )
}

save(
  analysis_set,
  met_gene_pair_tbl,
  met_gene_results,
  met_gene_all_results,
  file = paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/intermediate_RData/Fig5D_block05_metabolite_host_correlations_",
    analysis_set,
    ".RData"
  )
)


#=================================================================#
# Block 6. Species–metabolite correlations only
#=================================================================#
#
# Direct species–host gene correlations are intentionally omitted here.
# They were already evaluated in Panel C and do not need to be recalculated
# for this metabolite-centered analysis.
#=================================================================#

message(
  "Block 6: species–metabolite correlations only for analysis_set = ",
  analysis_set
)

species_met_pairs <- tibble()
species_gene_pairs <- tibble()
species_qc <- tibble()
species_gene_results[[analysis_set]] <- tibble()
species_gene_target_results[[analysis_set]] <- tibble()
species_use_log_results[[analysis_set]] <- matrix(
  numeric(0),
  nrow = 0,
  ncol = 0
)
species_use_raw_results[[analysis_set]] <- matrix(
  numeric(0),
  nrow = 0,
  ncol = 0
)
species_met_all_results[[analysis_set]] <- tibble()

if (
  "SampleID" %in% colnames(col_use) &&
  exists("coherence_data") &&
  !is.null(coherence_data$species_raw)
) {
  col_species <- col_use[
    !is.na(col_use$SampleID) &
      col_use$SampleID %in% rownames(coherence_data$species_raw),
    ,
    drop = FALSE
  ]
  
  if (nrow(col_species) >= 8) {
    species_use <- as.data.frame(
      coherence_data$species_raw[col_species$SampleID, , drop = FALSE],
      check.names = FALSE
    )
    
    rownames(species_use) <- rownames(col_species)
    species_use[] <- lapply(species_use, as.numeric)
    
    met_species <- met_use_log[rownames(col_species), , drop = FALSE]
    
    species_qc <- tibble(
      Species = colnames(species_use),
      species_prevalence = NA_real_,
      species_mean_abundance = NA_real_,
      species_sd = NA_real_,
      species_unclear_name = NA
    )
    
    for (i in seq_len(nrow(species_qc))) {
      species_qc$species_prevalence[i] <- mean(
        species_use[, species_qc$Species[i]] > 0,
        na.rm = TRUE
      )
      species_qc$species_mean_abundance[i] <- mean(
        species_use[, species_qc$Species[i]],
        na.rm = TRUE
      )
      species_qc$species_sd[i] <- sd(
        species_use[, species_qc$Species[i]],
        na.rm = TRUE
      )
      species_qc$species_unclear_name[i] <- str_detect(
        species_qc$Species[i],
        regex(
          "GGB|SGB|CAG|UBA|MAG|uncultured|metagenome|_sp_|bacterium|oral_taxon",
          ignore_case = TRUE
        )
      )
    }
    
    species_qc <- species_qc %>%
      mutate(
        keep_species =
          species_prevalence >= 0.2 &
          is.finite(species_sd) &
          species_sd > 0 &
          !species_unclear_name
      ) %>%
      arrange(
        desc(keep_species),
        desc(species_prevalence),
        desc(species_mean_abundance)
      )
    
    species_qc_results[[analysis_set]] <- species_qc
    
    write.csv(
      species_qc,
      paste0(
        "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
        "Fig5D_species_qc_",
        analysis_set,
        ".csv"
      ),
      row.names = FALSE
    )
    
    species_keep <- species_qc %>%
      filter(keep_species) %>%
      arrange(
        desc(species_prevalence),
        desc(species_mean_abundance)
      ) %>%
      slice_head(n = species_top_n) %>%
      pull(Species)
    
    if (length(species_keep) > 0) {
      species_use_raw_selected <- as.matrix(
        species_use[, species_keep, drop = FALSE]
      )
      storage.mode(species_use_raw_selected) <- "numeric"
      
      species_use_log <- log10(
        species_use_raw_selected + 1e-06
      )
      
      species_use_raw_results[[analysis_set]] <-
        species_use_raw_selected
      species_use_log_results[[analysis_set]] <-
        species_use_log
      
      message(
        "Testing ",
        ncol(species_use_log) * ncol(met_species),
        " species–metabolite pairs (",
        ncol(species_use_log),
        " species × ",
        ncol(met_species),
        " metabolites)."
      )
      
      species_met_rows <- vector(
        "list",
        ncol(species_use_log) * ncol(met_species)
      )
      
      row_i <- 1
      
      for (s in colnames(species_use_log)) {
        for (m in colnames(met_species)) {
          species_x <- species_use_log[, s]
          met_y <- met_species[, m]
          ok <- is.finite(species_x) & is.finite(met_y)
          
          species_met_rho <- NA_real_
          species_met_p <- NA_real_
          
          if (
            sum(ok) >= 8 &&
            sd(species_x[ok], na.rm = TRUE) > 0 &&
            sd(met_y[ok], na.rm = TRUE) > 0
          ) {
            species_met_result <- suppressWarnings(
              cor.test(
                species_x[ok],
                met_y[ok],
                method = "spearman",
                exact = FALSE
              )
            )
            species_met_rho <- unname(species_met_result$estimate)
            species_met_p <- species_met_result$p.value
          }
          
          species_met_rows[[row_i]] <- tibble(
            analysis_set = analysis_set,
            Species = s,
            Metabolite = m,
            n_species_met = sum(ok),
            species_met_spearman_rho = species_met_rho,
            species_met_spearman_p = species_met_p
          )
          
          row_i <- row_i + 1
        }
      }
      
      species_met_pairs_all <- bind_rows(
        species_met_rows[seq_len(row_i - 1)]
      ) %>%
        mutate(
          species_met_spearman_FDR = p.adjust(
            species_met_spearman_p,
            method = "BH"
          )
        )
      
      species_met_all_results[[analysis_set]] <-
        species_met_pairs_all
      
      species_met_pairs <- species_met_pairs_all %>%
        filter(species_met_spearman_p < 0.1) %>%
        arrange(
          species_met_spearman_p,
          desc(abs(species_met_spearman_rho))
        )
      
      species_met_results[[analysis_set]] <- species_met_pairs
      
      write.csv(
        species_met_pairs_all,
        paste0(
          "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
          "Fig5D_species_metabolite_all_pairs_",
          analysis_set,
          ".csv"
        ),
        row.names = FALSE
      )
      
      write.csv(
        species_met_pairs,
        paste0(
          "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
          "Fig5D_species_metabolite_pairs_",
          analysis_set,
          ".csv"
        ),
        row.names = FALSE
      )
      
      microbe_summary <- tibble(
        analysis_set = analysis_set,
        n_species_samples = nrow(col_species),
        n_species_before_QC = ncol(species_use),
        n_species_after_QC_selected = length(species_keep),
        n_species_metabolite_pairs_tested = nrow(species_met_pairs_all),
        n_species_metabolite_pairs_p_lt_0_1 = nrow(species_met_pairs),
        species_host_gene_correlations_recalculated = FALSE
      )
      
      write.csv(
        microbe_summary,
        paste0(
          "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/diagnostics/",
          "Fig5D_microbe_correlation_summary_",
          analysis_set,
          ".csv"
        ),
        row.names = FALSE
      )
      
      print(microbe_summary)
      
      if (nrow(species_met_pairs) > 0) {
        top_species_met_examples <- species_met_pairs %>%
          slice_head(n = 6) %>%
          mutate(pair_label = paste0(Species, " ~ ", Metabolite))
        
        plot_df_species_met <- bind_rows(
          lapply(
            seq_len(nrow(top_species_met_examples)),
            function(i) {
              tibble(
                Sample = rownames(col_species),
                TRG_plot = col_species$TRG_plot,
                Species = top_species_met_examples$Species[i],
                Metabolite = top_species_met_examples$Metabolite[i],
                pair_label = top_species_met_examples$pair_label[i],
                species_log10 = species_use_log[
                  rownames(col_species),
                  top_species_met_examples$Species[i]
                ],
                metabolite_log10 = met_species[
                  rownames(col_species),
                  top_species_met_examples$Metabolite[i]
                ]
              )
            }
          )
        )
        
        p_species_met_examples <- ggplot(
          plot_df_species_met,
          aes(x = species_log10, y = metabolite_log10)
        ) +
          geom_point(aes(shape = TRG_plot), size = 2, alpha = 0.8) +
          geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
          facet_wrap(~ pair_label, scales = "free", ncol = 3) +
          labs(
            title = paste0("Top species–metabolite pairs: ", analysis_set),
            x = "Species abundance, log10(x + 1e-06)",
            y = "Metabolite abundance, log10(x + 1e-06)"
          ) +
          theme_classic(base_size = 10)
        
        ggsave(
          paste0(
            "figures/Fig5D_top_species_metabolite_scatter_",
            analysis_set,
            ".svg"
          ),
          p_species_met_examples,
          width = 10,
          height = 6,
          device = "svg"
        )
      }
    } else {
      warning(
        analysis_set,
        ": no species passed QC. Species–metabolite analysis skipped."
      )
    }
  } else {
    warning(
      analysis_set,
      ": fewer than 8 samples had valid species data. Analysis skipped."
    )
  }
} else {
  warning(
    analysis_set,
    ": SampleID or coherence_data$species_raw is unavailable. Analysis skipped."
  )
}

save(
  analysis_set,
  species_qc,
  species_met_pairs,
  species_gene_pairs,
  species_qc_results,
  species_met_results,
  species_gene_results,
  species_use_log_results,
  species_use_raw_results,
  species_met_all_results,
  met_use_raw_results,
  met_use_log_results,
  vst_deg_results,
  sample_metadata_results,
  file = paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
    "intermediate_RData/Fig5D_block06_species_metabolite_correlations_",
    analysis_set,
    ".RData"
  )
)


#=================================================================#
# Block 7. Indirect species–metabolite–host axis assembly
#=================================================================#
#
# These are two-edge candidate axes:
#
#   species -- metabolite -- host gene
#
# No direct species–host gene correlation is calculated or used here.
# Panel C can be consulted separately when direct support is needed.
#=================================================================#

message(
  "Block 7: indirect two-edge axis assembly for analysis_set = ",
  analysis_set
)

triad_results[[analysis_set]] <- tibble()
join_diagnostics_results[[analysis_set]] <- list()

met_gene_axis_results <- met_gene_results[[analysis_set]] %>%
  filter(!coalesce(host_housekeeping_control, FALSE))

if (
  nrow(species_met_pairs) > 0 &&
  nrow(met_gene_axis_results) > 0
) {
  species_met_key_count <- species_met_pairs %>%
    count(analysis_set, Metabolite, name = "n_species_rows")
  
  met_gene_key_count <- met_gene_axis_results %>%
    count(analysis_set, Metabolite, name = "n_host_gene_rows")
  
  join_key_expansion <- full_join(
    species_met_key_count,
    met_gene_key_count,
    by = c("analysis_set", "Metabolite")
  ) %>%
    mutate(
      n_species_rows = coalesce(n_species_rows, 0L),
      n_host_gene_rows = coalesce(n_host_gene_rows, 0L),
      expected_axis_rows = n_species_rows * n_host_gene_rows
    ) %>%
    arrange(desc(expected_axis_rows), Metabolite)
  
  join_summary <- tibble(
    analysis_set = analysis_set,
    n_species_metabolite_rows = nrow(species_met_pairs),
    n_metabolite_host_gene_rows = nrow(met_gene_axis_results),
    n_shared_metabolites = sum(
      join_key_expansion$n_species_rows > 0 &
        join_key_expansion$n_host_gene_rows > 0
    ),
    expected_axis_rows_after_join = sum(
      join_key_expansion$expected_axis_rows
    ),
    direct_species_host_edge_used = FALSE
  )
  
  join_diagnostics_results[[analysis_set]] <- list(
    join_key_expansion = join_key_expansion,
    join_summary = join_summary
  )
  
  write.csv(
    join_key_expansion,
    paste0(
      "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/diagnostics/",
      "Fig5D_join_key_expansion_",
      analysis_set,
      ".csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    join_summary,
    paste0(
      "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/diagnostics/",
      "Fig5D_join_summary_",
      analysis_set,
      ".csv"
    ),
    row.names = FALSE
  )
  
  triad_results[[analysis_set]] <- species_met_pairs %>%
    inner_join(
      met_gene_axis_results,
      by = c("analysis_set", "Metabolite"),
      relationship = "many-to-many"
    ) %>%
    mutate(
      species_met_sign = sign(species_met_spearman_rho),
      host_met_sign = sign(best_cor_r),
      expected_species_gene_sign = sign(
        species_met_spearman_rho * best_cor_r
      ),
      axis_direction = case_when(
        species_met_sign > 0 & host_met_sign > 0 ~
          "Species higher - metabolite higher - host gene higher",
        species_met_sign > 0 & host_met_sign < 0 ~
          "Species higher - metabolite higher - host gene lower",
        species_met_sign < 0 & host_met_sign > 0 ~
          "Species higher - metabolite lower - host gene lower",
        species_met_sign < 0 & host_met_sign < 0 ~
          "Species higher - metabolite lower - host gene higher",
        TRUE ~ "Direction unavailable"
      ),
      axis_evidence_level = case_when(
        species_met_spearman_p < 0.01 & best_cor_p < 0.01 ~
          "A. Both edges p < 0.01",
        species_met_spearman_p < 0.05 & best_cor_p < 0.05 ~
          "B. Both edges p < 0.05",
        species_met_spearman_p < 0.10 & best_cor_p < 0.10 ~
          "C. Both edges show trend",
        TRUE ~ "D. Exploratory"
      ),
      indirect_axis_rank_score =
        -log10(
          pmax(species_met_spearman_p, .Machine$double.xmin)
        ) * abs(species_met_spearman_rho) +
        -log10(
          pmax(best_cor_p, .Machine$double.xmin)
        ) * abs(best_cor_r),
      direct_species_host_correlation_recalculated = FALSE
    ) %>%
    arrange(
      axis_evidence_level,
      desc(indirect_axis_rank_score),
      species_met_spearman_p,
      best_cor_p
    )
  
  write.csv(
    triad_results[[analysis_set]],
    paste0(
      "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
      "Fig5D_candidate_indirect_axes_",
      analysis_set,
      ".csv"
    ),
    row.names = FALSE
  )
  
  # Historical filename retained for downstream compatibility.
  write.csv(
    triad_results[[analysis_set]],
    paste0(
      "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
      "Fig5D_candidate_triad_axes_",
      analysis_set,
      ".csv"
    ),
    row.names = FALSE
  )
  
  axis_summary <- triad_results[[analysis_set]] %>%
    count(axis_evidence_level, name = "n_axes")
  
  write.csv(
    axis_summary,
    paste0(
      "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/diagnostics/",
      "Fig5D_indirect_axis_summary_",
      analysis_set,
      ".csv"
    ),
    row.names = FALSE
  )
  
  print(axis_summary)
} else {
  warning(
    analysis_set,
    ": indirect axis assembly skipped because species–metabolite or ",
    "metabolite–host candidate pairs are empty."
  )
}

save(
  analysis_set,
  triad_results,
  join_diagnostics_results,
  file = paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
    "intermediate_RData/Fig5D_block07_indirect_axis_assembly_",
    analysis_set,
    ".RData"
  )
)


#=================================================================#
# Block 8. Final metabolite–host gene heatmap visualization
#=================================================================#
#
# This block is now inside the same script, but it is deliberately kept
# after the screening blocks. Run it after Block 5 or after Block 7.
#
# Output:
#   figures/Fig5D_metabolite_DEG_correlation_heatmap_<analysis_set>.svg
#
#-----------------------------------------------------------------#

message("Block 8: metabolite–host gene heatmap for analysis_set = ", analysis_set)

if (
  !is.null(met_gene_results[[analysis_set]]) &&
  nrow(met_gene_results[[analysis_set]]) > 0
) {
  plot_genes <- met_gene_results[[analysis_set]] %>%
    group_by(Gene) %>%
    summarise(
      best_p = min(best_cor_p, na.rm = TRUE),
      best_abs_cor = max(abs(best_cor_r), na.rm = TRUE),
      host_logFC = first(host_logFC),
      host_min_p = first(host_min_p),
      .groups = "drop"
    ) %>%
    arrange(
      best_p,
      desc(best_abs_cor),
      desc(abs(host_logFC))
    ) %>%
    slice_head(n = 35) %>%
    pull(Gene)
  
  plot_metabolites <- met_gene_results[[analysis_set]] %>%
    filter(Gene %in% plot_genes) %>%
    group_by(Metabolite) %>%
    summarise(
      best_p = min(best_cor_p, na.rm = TRUE),
      best_abs_cor = max(abs(best_cor_r), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(
      best_p,
      desc(best_abs_cor)
    ) %>%
    slice_head(n = 25) %>%
    pull(Metabolite)
  
  met_gene_rho_mat <- met_gene_all_results[[analysis_set]] %>%
    filter(
      Gene %in% plot_genes,
      Metabolite %in% plot_metabolites
    ) %>%
    select(Metabolite, Gene, spearman_rho) %>%
    pivot_wider(
      names_from = Gene,
      values_from = spearman_rho,
      values_fn = list(spearman_rho = mean)
    ) %>%
    as.data.frame(check.names = FALSE)
  
  rownames(met_gene_rho_mat) <- met_gene_rho_mat$Metabolite
  met_gene_rho_mat$Metabolite <- NULL
  met_gene_rho_mat <- as.matrix(
    met_gene_rho_mat[plot_metabolites, plot_genes, drop = FALSE]
  )
  
  met_gene_pearson_sig_mat <- met_gene_all_results[[analysis_set]] %>%
    filter(
      Gene %in% plot_genes,
      Metabolite %in% plot_metabolites
    ) %>%
    mutate(sig = coalesce(pearson_p < 0.05, FALSE)) %>%
    select(Metabolite, Gene, sig) %>%
    pivot_wider(
      names_from = Gene,
      values_from = sig,
      values_fill = FALSE,
      values_fn = list(sig = any)
    ) %>%
    as.data.frame(check.names = FALSE)
  
  rownames(met_gene_pearson_sig_mat) <- met_gene_pearson_sig_mat$Metabolite
  met_gene_pearson_sig_mat$Metabolite <- NULL
  met_gene_pearson_sig_mat <- as.matrix(
    met_gene_pearson_sig_mat[plot_metabolites, plot_genes, drop = FALSE]
  )
  
  met_gene_spearman_sig_mat <- met_gene_all_results[[analysis_set]] %>%
    filter(
      Gene %in% plot_genes,
      Metabolite %in% plot_metabolites
    ) %>%
    mutate(sig = coalesce(spearman_p < 0.05, FALSE)) %>%
    select(Metabolite, Gene, sig) %>%
    pivot_wider(
      names_from = Gene,
      values_from = sig,
      values_fill = FALSE,
      values_fn = list(sig = any)
    ) %>%
    as.data.frame(check.names = FALSE)
  
  rownames(met_gene_spearman_sig_mat) <- met_gene_spearman_sig_mat$Metabolite
  met_gene_spearman_sig_mat$Metabolite <- NULL
  met_gene_spearman_sig_mat <- as.matrix(
    met_gene_spearman_sig_mat[plot_metabolites, plot_genes, drop = FALSE]
  )
  
  met_gene_partial_sig_mat <- met_gene_all_results[[analysis_set]] %>%
    filter(
      Gene %in% plot_genes,
      Metabolite %in% plot_metabolites
    ) %>%
    mutate(sig = coalesce(partial_spearman_p < 0.05, FALSE)) %>%
    select(Metabolite, Gene, sig) %>%
    pivot_wider(
      names_from = Gene,
      values_from = sig,
      values_fill = FALSE,
      values_fn = list(sig = any)
    ) %>%
    as.data.frame(check.names = FALSE)
  
  rownames(met_gene_partial_sig_mat) <- met_gene_partial_sig_mat$Metabolite
  met_gene_partial_sig_mat$Metabolite <- NULL
  met_gene_partial_sig_mat <- as.matrix(
    met_gene_partial_sig_mat[plot_metabolites, plot_genes, drop = FALSE]
  )
  
  host_logfc_plot <- met_gene_all_results[[analysis_set]] %>%
    filter(Gene %in% plot_genes) %>%
    distinct(Gene, .keep_all = TRUE) %>%
    mutate(Gene = factor(Gene, levels = plot_genes)) %>%
    arrange(Gene)
  
  host_logfc_limit <- max(
    0.5,
    quantile(abs(host_logfc_plot$host_logFC), 0.95, na.rm = TRUE)
  )
  
  host_neglog10p_limit <- max(
    1,
    quantile(
      -log10(pmax(host_logfc_plot$host_min_p, .Machine$double.xmin)),
      0.95,
      na.rm = TRUE
    )
  )
  
  ht_met_gene <- Heatmap(
    met_gene_rho_mat,
    name = "Spearman rho",
    col = colorRamp2(
      c(-1, 0, 1),
      c("#2166AC", "#FFFFFF", "#B2182B")
    ),
    na_col = "#F2F2F2",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 90,
    column_names_max_height = unit(60, "mm"),
    top_annotation = HeatmapAnnotation(
      Host_LogFC = host_logfc_plot$host_logFC,
      Host_neglog10P = -log10(
        pmax(host_logfc_plot$host_min_p, .Machine$double.xmin)
      ),
      col = list(
        Host_LogFC = colorRamp2(
          c(-host_logfc_limit, 0, host_logfc_limit),
          c(group_cols["non_pCR"], "#F7F7F7", group_cols["pCR"])
        ),
        Host_neglog10P = colorRamp2(
          c(0, host_neglog10p_limit),
          c("#F7FBFF", "#2C5AA0")
        )
      ),
      annotation_name_gp = gpar(fontsize = 8)
    ),
    rect_gp = gpar(col = "white", lwd = 0.2),
    width = unit(max(120, ncol(met_gene_rho_mat) * 4.8), "mm"),
    height = unit(max(80, nrow(met_gene_rho_mat) * 4.8), "mm"),
    layer_fun = function(j, i, x, y, width, height, fill) {
      if (any(met_gene_pearson_sig_mat[cbind(i, j)], na.rm = TRUE)) {
        grid.points(
          x[met_gene_pearson_sig_mat[cbind(i, j)]],
          y[met_gene_pearson_sig_mat[cbind(i, j)]],
          pch = 1,
          size = unit(3.5, "mm"),
          gp = gpar(col = "black", lwd = 0.9)
        )
      }
      
      if (any(met_gene_spearman_sig_mat[cbind(i, j)], na.rm = TRUE)) {
        grid.points(
          x[met_gene_spearman_sig_mat[cbind(i, j)]],
          y[met_gene_spearman_sig_mat[cbind(i, j)]],
          pch = 2,
          size = unit(2.5, "mm"),
          gp = gpar(col = "black", lwd = 0.9)
        )
      }
      
      if (any(met_gene_partial_sig_mat[cbind(i, j)], na.rm = TRUE)) {
        grid.points(
          x[met_gene_partial_sig_mat[cbind(i, j)]],
          y[met_gene_partial_sig_mat[cbind(i, j)]],
          pch = 4,
          size = unit(2.2, "mm"),
          gp = gpar(col = "black", lwd = 0.9)
        )
      }
    }
  )
  
  lgd_cor <- Legend(
    title = "Correlation p < 0.05",
    labels = c("Pearson", "Spearman", "Partial Spearman"),
    type = "points",
    pch = c(1, 2, 4),
    legend_gp = gpar(col = "black", lwd = 0.9),
    size = unit(c(3.5, 2.5, 2.2), "mm"),
    labels_gp = gpar(fontsize = 8),
    title_gp = gpar(fontsize = 9, fontface = "bold")
  )
  
  svglite(
    paste0(
      "figures/Fig5D_metabolite_DEG_correlation_heatmap_",
      analysis_set,
      ".svg"
    ),
    width = max(10, min(22, 6 + ncol(met_gene_rho_mat) * 0.30)),
    height = max(7, min(18, 4 + nrow(met_gene_rho_mat) * 0.28))
  )
  
  draw(
    ht_met_gene,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    annotation_legend_list = list(lgd_cor),
    padding = unit(c(2, 4, 2, 2), "mm")
  )
  
  dev.off()
  
} else {
  message(analysis_set, ": no screened metabolite–gene pairs. Heatmap skipped.")
}


#=================================================================#
# Block 9. Save essential combined objects
#=================================================================#

save(
  met_gene_results,
  met_gene_all_results,
  met_qc_results,
  host_deg_results,
  host_expr_qc_results,
  host_gene_screen_results,
  met_use_log_results,
  met_use_raw_results,
  vst_deg_results,
  species_use_log_results,
  species_use_raw_results,
  species_met_all_results,
  triad_results,
  sample_metadata_results,
  species_qc_results,
  species_met_results,
  species_gene_results,
  species_gene_target_results,
  join_diagnostics_results,
  analysis_summary_results,
  file = paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
    "Fig5D_microbe_metabolite_host_axis_screening_results.RData"
  )
)

# This compatibility file also preserves the complete matched matrices needed
# for downstream three-edge Spearman and partial-Spearman recalculation.
save(
  met_gene_results,
  met_gene_all_results,
  met_qc_results,
  host_deg_results,
  host_expr_qc_results,
  host_gene_screen_results,
  met_use_log_results,
  met_use_raw_results,
  vst_deg_results,
  species_use_log_results,
  species_use_raw_results,
  species_met_all_results,
  species_qc_results,
  species_met_results,
  sample_metadata_results,
  analysis_summary_results,
  file = paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
    "Fig5D_metabolite_host_priority_screening_results.RData"
  )
)

message(
  "Saved full matched matrices and downstream-compatible RData files for triad recalculation."
)


#=================================================================#
# Optional Block 10. Reload saved objects later
#=================================================================#
#
# Use this block only in a new R session when you want to draw additional
# figures or inspect existing results without recomputing all correlations.
#
# load(
#   "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5D_microbe_metabolite_host_axis_screening_results.RData"
# )
#
# names(triad_results)
# head(triad_results[["before"]])
# head(met_gene_results[["before"]])
#
#=================================================================#
