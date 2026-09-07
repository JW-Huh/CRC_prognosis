#-----------------------------------------------------------------#
#
# Fig. 4C. Species-level microbe–host DEG Spearman-correlation heatmap
# with phylum- and family-level microbial taxonomy annotations
#
# Analysis design:
#   1. All before-TNT RNA-microbiome matched samples are used
#      (vst_rna / col_rna; expected n = 20), without metabolite matching.
#   2. The host panel is fixed to the 18 genes specified below.
#   3. All species require prevalence >= 0.30 (at least 6 of 20 samples).
#      Microbial LogFC and Wilcoxon p are retained in the analysis objects but
#      are not displayed in this taxonomy-annotation version.
#   4. A species is eligible when it has at least one Spearman p < 0.05
#      correlation or at least two Spearman p < 0.10 correlations with the
#      fixed 18-gene panel.
#   5. Weissella confusa is explicitly excluded. Bacteroides fragilis is
#      screened as a priority taxon and displayed in parallel included and
#      excluded figure versions.
#   6. Eligible species are ranked primarily by their correlation evidence
#      across the predefined host panel; priority status is a later tie-breaker.
#   7. Species are clustered globally by their 18-gene Spearman-rho profiles.
#      Phylum and family are displayed as side annotations but do not determine
#      row order or row blocks. Host genes remain in four fixed biological
#      blocks, while genes within each block are clustered by microbial
#      correlations.
#   8. Near-identical abundance-rank profiles are collapsed only among
#      non-priority species. Prespecified species are retained independently.
#   9. Eubacterium rectale and Roseburia inulinivorans are retained in both
#      figure variants when they pass prevalence >= 0.30, even if they do not
#      meet the correlation-based row-selection rule. Anaerostipes hadrus is
#      excluded because it showed no significant host-gene correlation.
#  10. Target-species group means and feature-name mappings are exported to
#      audit microbial LogFC direction and row-label alignment.
#  11. Phylum and family annotations are parsed from input/species.csv by
#      matching the s__, p__, and f__ fields in clade_name to the species rows.
#  12. Families represented by only one displayed species across the two
#      variants are grouped as Others, except Fusobacteriaceae,
#      Bifidobacteriaceae, Oscillospiraceae, and Veillonellaceae.
#
# Input:
#   input/species.csv
#   host_RNAseq/results_clean_metabolite_host/
#     host_met_rna_matched_input.RData
#   host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/
#     Fig5_host_RNAseq_inputs.RData
#   host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/
#     Fig5B_C_downstream_inputs.RData
#
# Output:
#   figures/Fig4C_v11_species_host_DEG_Spearman_taxonomy_Bfragilis_included.svg
#   figures/Fig4C_v11_species_host_DEG_Spearman_taxonomy_Bfragilis_excluded.svg
#   host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/
#     Fig5_priority_species_screen.csv
#     Fig5_priority_species_significant_pairs.csv
#     Fig5_pCR_enriched_species_with_multiple_host_correlations.csv
#     Fig5_species_block_coherence_ranking.csv
#     Fig5_species_weak_pattern_review.csv
#     Fig4C_v11_species_host_heatmap_final_species_Bfragilis_included.csv
#     Fig4C_v11_species_host_heatmap_final_species_Bfragilis_excluded.csv
#     Fig4C_v11_target_species_direction_audit.csv
#     Fig4C_v11_species_feature_name_map.csv
#     Fig4C_v11_species_taxonomy_map.csv
#     Fig4C_v11_forced_species_prevalence_audit.csv
#     Fig5_species_host_heatmap_final_genes.csv
#     Fig4C_v11_species_host_heatmap_Spearman_taxonomy_Bfragilis_variants.RData
#
#-----------------------------------------------------------------#

rm(list = ls())
options(stringsAsFactors = FALSE)

message(
  "Running Fig4C microbe-host script v11 ",
  paste0(
    "(matched n = 20; prevalence >= 30%; correlation-pattern clustering; ",
    "phylum/family parsed from input/species.csv; ",
    "direct species-name mapping; ",
    "forced butyrate producers; ",
    "B. fragilis variants and direction audit)"
  )
)

setwd("D:/2-연구/2-CRC metagenomics/")

#-----------------------------------------------------------------#
# 0. Packages
#-----------------------------------------------------------------#

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in c(
  "dplyr", "tidyr", "tibble", "stringr", "readr",
  "svglite", "circlize"
)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, type = "binary")
  }
}

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(svglite)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
})

#-----------------------------------------------------------------#
# 1. Output folders and input loading
#-----------------------------------------------------------------#

dir.create("figures", showWarnings = FALSE)

dir.create(
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq",
  recursive = TRUE,
  showWarnings = FALSE
)

# Load the sample-matching object created by
# "7-9.1 host RNAseq - data matching and DESeq2.R" into an isolated
# environment. Its col_rna and s objects originate from the same
# microbiome metadata backbone and must be used together for matching.
matching_env <- new.env(parent = baseenv())
load(
  "host_RNAseq/results_clean_metabolite_host/host_met_rna_matched_input.RData",
  envir = matching_env
)

required_matching_objects <- c("m", "s", "col_rna")
missing_matching_objects <- required_matching_objects[
  !vapply(required_matching_objects, exists, logical(1), envir = matching_env)
]

if (length(missing_matching_objects) > 0) {
  stop(
    "The sample-matching RData is missing: ",
    paste(missing_matching_objects, collapse = ", "),
    ". Re-run the data-matching script before this figure script.",
    call. = FALSE
  )
}

load(
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5_host_RNAseq_inputs.RData"
)

load(
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5B_C_downstream_inputs.RData"
)

#-----------------------------------------------------------------#
# 2. Colors and p-value labels
#-----------------------------------------------------------------#

if (!exists("group_cols")) {
  group_cols <- c(
    "pCR" = "#5AB4AC",
    "non_pCR" = "#D97A6C"
  )
}

if (!"pCR" %in% names(group_cols) && "CR" %in% names(group_cols)) {
  group_cols["pCR"] <- group_cols["CR"]
}

if (!"non_pCR" %in% names(group_cols) && "nonCR" %in% names(group_cols)) {
  group_cols["non_pCR"] <- group_cols["nonCR"]
}

if (!"pCR" %in% names(group_cols)) {
  group_cols["pCR"] <- "#5AB4AC"
}

if (!"non_pCR" %in% names(group_cols)) {
  group_cols["non_pCR"] <- "#D97A6C"
}

p_to_star <- function(p) {
  case_when(
    !is.na(p) & p < 0.001 ~ "***",
    !is.na(p) & p < 0.01 ~ "**",
    !is.na(p) & p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

#-----------------------------------------------------------------#
# 3. RNA-microbiome sample matching: before TNT only
#-----------------------------------------------------------------#

# Use the authoritative 20-sample RNA metadata table created by the matching
# script. Species abundances are subsequently selected from s using these
# exact microbiome SampleID values, without any second sample filter.
col_rna <- as.data.frame(matching_env$col_rna, stringsAsFactors = FALSE)

required_metadata_cols <- c(
  "RNA_sample_id", "SampleID", "TNT", "timepoint", "TRG_plot"
)

if (!all(required_metadata_cols %in% colnames(col_rna))) {
  stop(
    "Matched col_rna is missing required columns: ",
    paste(setdiff(required_metadata_cols, colnames(col_rna)), collapse = ", "),
    call. = FALSE
  )
}

rownames(col_rna) <- as.character(col_rna$RNA_sample_id)

# The matching script already restricts m_rna/col_rna to before-TNT tumor
# biopsies and excludes normal tissue. These checks reproduce that design
# without imposing an additional literal timepoint == "pre" requirement.
baseline_keep <-
  as.character(col_rna$TNT) == "Before" &
  tolower(as.character(col_rna$timepoint)) != "normal" &
  as.character(col_rna$TRG_plot) %in% c("pCR", "non_pCR")

col_rna <- col_rna[baseline_keep, , drop = FALSE]

if (nrow(col_rna) != 20) {
  stop(
    "The authoritative matched col_rna contains ", nrow(col_rna),
    " eligible baseline tumor samples rather than 20. Re-run and inspect ",
    "the data-matching script before continuing.",
    call. = FALSE
  )
}

missing_rna_expression <- setdiff(rownames(col_rna), colnames(vst_rna))

if (length(missing_rna_expression) > 0) {
  stop(
    "Matched RNA samples absent from vst_rna: ",
    paste(missing_rna_expression, collapse = ", "),
    call. = FALSE
  )
}

if (anyDuplicated(rownames(col_rna))) {
  stop("Duplicated RNA_sample_id values remained after matching.", call. = FALSE)
}

if (anyDuplicated(col_rna$SampleID)) {
  stop(
    "One microbiome SampleID is linked to more than one RNA sample.",
    call. = FALSE
  )
}

vst_rna <- vst_rna[, rownames(col_rna), drop = FALSE]

# Reconstruct the microbiome matrix exactly as in the diversity and coherence
# input-processing scripts. The saved object s is feature x sample: the
# Species column contains biological names and the Sample_* columns contain
# abundances. Species names must therefore come directly from s$Species; they
# must never be inferred from V### positions or another table's column order.
s_raw <- as.data.frame(
  matching_env$s,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"Species" %in% colnames(s_raw)) {
  stop(
    "The microbiome table s does not contain the required Species column. ",
    "Recreate m_s_t_only.RData with the diversity-analysis script.",
    call. = FALSE
  )
}

# The saved object s contains only species labels and abundance values. Recover
# higher taxonomy from the original species-level MetaPhlAn table, following
# the clade_name parsing used in the diversity-analysis script. family.csv
# provides family-level abundance but cannot map an individual species row back
# to its family, so the species-level lineage is used for this annotation.
if (!file.exists("input/species.csv")) {
  stop(
    "input/species.csv was not found under the project directory.",
    call. = FALSE
  )
}

species_taxonomy_source <- read.csv(
  "input/species.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"clade_name" %in% colnames(species_taxonomy_source)) {
  stop(
    "input/species.csv does not contain the required clade_name column.",
    call. = FALSE
  )
}

species_taxonomy_raw <- tibble(
  Species_source = stringr::str_match(
    as.character(species_taxonomy_source$clade_name),
    "(?:^|\\|)s__([^|]+)"
  )[, 2],
  Phylum = stringr::str_match(
    as.character(species_taxonomy_source$clade_name),
    "(?:^|\\|)p__([^|]+)"
  )[, 2],
  Family = stringr::str_match(
    as.character(species_taxonomy_source$clade_name),
    "(?:^|\\|)f__([^|]+)"
  )[, 2]
) %>%
  filter(!is.na(Species_source), str_squish(Species_source) != "") %>%
  mutate(
    Species_key = str_to_lower(
      str_replace_all(str_squish(Species_source), "_", " ")
    ),
    across(
      c(Phylum, Family),
      ~ str_squish(as.character(.x))
    ),
    across(
      c(Phylum, Family),
      ~ ifelse(
        is.na(.x) | .x == "",
        "Unclassified",
        .x
      )
    )
  )

taxonomy_conflicts <- species_taxonomy_raw %>%
  group_by(Species_key) %>%
  summarise(
    n_phyla = n_distinct(Phylum[Phylum != "Unclassified"]),
    n_families = n_distinct(Family[Family != "Unclassified"]),
    .groups = "drop"
  ) %>%
  filter(n_phyla > 1 | n_families > 1)

if (nrow(taxonomy_conflicts) > 0) {
  stop(
    "Conflicting phylum or family assignments were found for: ",
    paste(taxonomy_conflicts$Species_key, collapse = ", "),
    call. = FALSE
  )
}

species_taxonomy <- species_taxonomy_raw %>%
  group_by(Species_key) %>%
  summarise(
    Species_source = first(Species_source),
    Phylum = first(c(unique(Phylum[Phylum != "Unclassified"]), "Unclassified")),
    Family = first(c(unique(Family[Family != "Unclassified"]), "Unclassified")),
    .groups = "drop"
  )

matched_microbiome_ids <- as.character(col_rna$SampleID)
missing_microbiome_ids <- setdiff(matched_microbiome_ids, colnames(s_raw))

if (length(missing_microbiome_ids) > 0) {
  stop(
    "Matched microbiome SampleID columns absent from s: ",
    paste(missing_microbiome_ids, collapse = ", "),
    call. = FALSE
  )
}

species_feature_by_sample <- s_raw %>%
  dplyr::select(Species, dplyr::all_of(matched_microbiome_ids)) %>%
  dplyr::filter(
    !is.na(Species),
    stringr::str_squish(as.character(Species)) != ""
  ) %>%
  dplyr::mutate(
    Species = as.character(Species),
    dplyr::across(
      dplyr::all_of(matched_microbiome_ids),
      ~ suppressWarnings(as.numeric(as.character(.x)))
    )
  )

if (anyNA(species_feature_by_sample[, matched_microbiome_ids, drop = FALSE])) {
  stop(
    "Missing or non-numeric abundance values were detected in the 20 matched ",
    "SampleID columns of s.",
    call. = FALSE
  )
}

# Follow the upstream coherence preprocessing exactly: duplicate Species labels
# are summed before the feature-by-sample table is transposed.
species_feature_by_sample <- species_feature_by_sample %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(matched_microbiome_ids),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

species_feature_matrix <- as.data.frame(
  species_feature_by_sample,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(species_feature_matrix) <- species_feature_matrix$Species
species_feature_matrix$Species <- NULL

species_all <- base::t(as.matrix(species_feature_matrix))
mode(species_all) <- "numeric"
species_all <- species_all[matched_microbiome_ids, , drop = FALSE]

if (
  !identical(rownames(species_all), matched_microbiome_ids) ||
  anyNA(colnames(species_all)) ||
  anyDuplicated(colnames(species_all)) ||
  any(grepl("^V[0-9]+$", colnames(species_all)))
) {
  stop(
    "The directly reconstructed species matrix failed its SampleID or ",
    "species-name integrity check.",
    call. = FALSE
  )
}

species_taxonomy <- tibble(
  Species = colnames(species_all),
  Species_key = str_to_lower(
    str_replace_all(str_squish(colnames(species_all)), "_", " ")
  )
) %>%
  left_join(species_taxonomy, by = "Species_key") %>%
  mutate(
    Phylum = coalesce(Phylum, "Unclassified"),
    Family = coalesce(Family, "Unclassified")
  ) %>%
  mutate(Family_source = Family) %>%
  select(Species, Species_source, Phylum, Family_source, Family)

species_name_map <- tibble(
  source_feature = colnames(species_all),
  source_feature_index = seq_along(colnames(species_all)),
  Species = colnames(species_all)
)

write_csv(
  species_name_map,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
    "Fig4C_v11_species_feature_name_map.csv"
  )
)

species_int <- as.data.frame(species_all, check.names = FALSE)
rownames(species_int) <- rownames(col_rna)

message(
  "Reconstructed species matrix directly from s$Species and the matched ",
  "SampleID columns: ", nrow(species_int), " samples x ",
  ncol(species_int), " named species."
)

if (anyNA(species_int)) {
  stop(
    "Missing or non-numeric values were detected in the species table.",
    call. = FALSE
  )
}

if (nrow(species_int) < 8 || ncol(vst_rna) < 8) {
  stop(
    "Too few before-TNT matched RNA-seq and microbiome samples remained.",
    call. = FALSE
  )
}

if (!identical(rownames(species_int), colnames(vst_rna))) {
  stop("Sample order mismatch between species_int and vst_rna.", call. = FALSE)
}

if (length(unique(col_rna$TRG_plot)) < 2) {
  stop("TRG_plot must contain both pCR and non_pCR groups.", call. = FALSE)
}

if (nrow(col_rna) != 20) {
  stop(
    "Expected 20 baseline RNA-microbiome matched samples, but retained ",
    nrow(col_rna),
    ". Check SampleID mapping and baseline metadata filters.",
    call. = FALSE
  )
}

matching_audit <- col_rna %>%
  tibble::rownames_to_column("RNA_matrix_id") %>%
  dplyr::select(
    dplyr::any_of(c(
      "RNA_matrix_id", "RNA_sample_id", "SampleID", "SNU_ID", "SNU_tp",
      "SubjectID", "Chart_numb", "TNT", "timepoint", "TRG_plot"
    ))
  ) %>%
  dplyr::mutate(
    species_table_match = SampleID %in% rownames(species_all),
    RNA_matrix_match = RNA_matrix_id %in% colnames(vst_rna)
  )

readr::write_csv(
  matching_audit,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
    "Fig5_RNA_microbiome_20sample_matching_audit.csv"
  )
)

message(
  "Species-table matching method: direct s$Species x matched SampleID ",
  "reconstruction; duplicate labels summed before transposition."
)
message("RNA-microbiome matched baseline samples: ", nrow(col_rna))
print(col_rna %>% count(TRG_plot, name = "n"))

#-----------------------------------------------------------------#
# 4. Fixed host panel and priority microbial species
#-----------------------------------------------------------------#
#
# The host panel is fixed. Genes are retained whenever they are available
# in vst_rna and deg_rna, even when no individual microbe correlation
# reaches p < 0.05.
#
# Species screening:
#   - All species: prevalence >= 0.30, corresponding to at least six nonzero
#     samples among the 20 matched baseline samples.
#   - Microbial LogFC and Wilcoxon p are retained for audit and ranking
#     tie-breakers only; they do not determine entry into correlation screening.
#
# Final heatmap:
#   - A microbial species must have at least one Spearman p < 0.05 result or
#     at least two Spearman p < 0.10 results with the fixed host panel.
#   - Bacteroides fragilis and Clostridium symbiosum are screened as priority
#     taxa when they satisfy prevalence and correlation criteria.
#   - One output retains B. fragilis and the other excludes it, backfilling the
#     next-ranked eligible species to preserve the same heatmap dimensions.
#-----------------------------------------------------------------#

# Keep the four biological blocks fixed in the requested sequence. The vector
# below defines block membership; the display order within each block is later
# determined by hierarchical clustering of the microbial-correlation profiles.
target_genes <- c(
  "MSX1", "NOTUM", "UCA1", "NKD1", "HTRA1", "PLOD1",
  "NINJ1", "CX3CL1", "ADAM8", "ICAM1",
  "LIG4", "CLPTM1", "METTL14", "ZBTB7B",
  "PIK3CB", "CRK", "KLRK1", "RAB27A"
)

gene_block_tbl <- tibble(
  Gene = target_genes,
  Host_block = factor(
    c(
      rep("Adenoma signature + EMT", 6),
      rep("Gut homing + leukocyte trafficking", 4),
      rep("T-cell differentiation", 4),
      rep("Lymphocyte immunity + TNF signaling", 4)
    ),
    levels = c(
      "Adenoma signature + EMT",
      "Gut homing + leukocyte trafficking",
      "T-cell differentiation",
      "Lymphocyte immunity + TNF signaling"
    )
  )
)

if (any(!target_genes %in% rownames(vst_rna))) {
  stop(
    "Target genes absent from vst_rna: ",
    paste(setdiff(target_genes, rownames(vst_rna)), collapse = ", "),
    call. = FALSE
  )
}

target_genes <- target_genes[target_genes %in% rownames(vst_rna)]
gene_block_tbl <- gene_block_tbl %>%
  filter(Gene %in% target_genes)

if (length(target_genes) != 18) {
  stop("The fixed host panel must contain all 18 genes.", call. = FALSE)
}

forced_species_labels <- c(
  "Eubacterium rectale",
  "Roseburia inulinivorans"
)

forced_species_patterns <- c(
  "Eubacterium[_ ]rectale",
  "Roseburia[_ ]inulinivorans"
)

forced_species_regex <- regex(
  paste(forced_species_patterns, collapse = "|"),
  ignore_case = TRUE
)

#-----------------------------------------------------------------#
# 5. Species statistics and inclusive prevalence filtering
#-----------------------------------------------------------------#

species_stat <- tibble(
  Species = colnames(species_int),
  species_mean_pCR = NA_real_,
  species_mean_non_pCR = NA_real_,
  species_detected_pCR_n = NA_integer_,
  species_detected_non_pCR_n = NA_integer_,
  species_logFC = NA_real_,
  species_p = NA_real_,
  species_prevalence = NA_real_,
  species_nonzero_n = NA_integer_,
  species_mean_abundance = NA_real_
)

for (i in seq_len(nrow(species_stat))) {
  species_stat$species_mean_pCR[i] <- mean(
    species_int[
      col_rna$TRG_plot == "pCR",
      species_stat$Species[i]
    ],
    na.rm = TRUE
  )

  species_stat$species_mean_non_pCR[i] <- mean(
    species_int[
      col_rna$TRG_plot == "non_pCR",
      species_stat$Species[i]
    ],
    na.rm = TRUE
  )

  species_stat$species_detected_pCR_n[i] <- sum(
    species_int[
      col_rna$TRG_plot == "pCR",
      species_stat$Species[i]
    ] > 0,
    na.rm = TRUE
  )

  species_stat$species_detected_non_pCR_n[i] <- sum(
    species_int[
      col_rna$TRG_plot == "non_pCR",
      species_stat$Species[i]
    ] > 0,
    na.rm = TRUE
  )

  species_stat$species_logFC[i] <- log2(
    (species_stat$species_mean_pCR[i] + 1e-06) /
      (species_stat$species_mean_non_pCR[i] + 1e-06)
  )
  
  species_stat$species_p[i] <- suppressWarnings(
    wilcox.test(
      species_int[
        col_rna$TRG_plot == "pCR",
        species_stat$Species[i]
      ],
      species_int[
        col_rna$TRG_plot == "non_pCR",
        species_stat$Species[i]
      ],
      exact = FALSE
    )$p.value
  )
  
  species_stat$species_prevalence[i] <- mean(
    species_int[, species_stat$Species[i]] > 0,
    na.rm = TRUE
  )
  
  species_stat$species_nonzero_n[i] <- sum(
    species_int[, species_stat$Species[i]] > 0,
    na.rm = TRUE
  )
  
  species_stat$species_mean_abundance[i] <- mean(
    species_int[, species_stat$Species[i]],
    na.rm = TRUE
  )
}

species_stat <- species_stat %>%
  mutate(
    species_FDR = p.adjust(species_p, method = "BH"),
    species_neglog10p = -log10(pmax(species_p, .Machine$double.xmin)),
    species_neglog10p_capped = pmin(species_neglog10p, 10),
    species_p_star = p_to_star(species_p),
    species_logFC_direction = case_when(
      species_logFC > 0 ~ "pCR higher",
      species_logFC < 0 ~ "non-pCR higher",
      TRUE ~ "equal"
    ),
    
    # Priority taxa were selected a priori for CRC, oral-gut transmission,
    # immune regulation, or intestinal-homeostasis relevance. They remain
    # independent during duplicate-profile collapse and serve as tie-breakers.
    species_priority = str_detect(
      Species,
      regex(
        paste(
          c(
            "Eubacterium[_ ]rectale",
            "Roseburia[_ ]inulinivorans",
            "Faecalibacterium[_ ]prausnitzii",
            "Coprococcus[_ ]catus",
            "Clostridium[_ ]bolteae",
            "^Fusobacterium([_ ]|$)",
            "Parvimonas[_ ]micra",
            "Hungatella[_ ]hathewayi",
            "Gemella[_ ]morbillorum",
            "Peptostreptococcus[_ ]stomatis",
            "(Segatella|Prevotella)[_ ]copri",
            "Bacteroides[_ ]fragilis",
            "Clostridium[_ ]symbiosum",
            "Escherichia[_ ]coli",
            "^Campylobacter([_ ]|$)",
            "Ruminococcus[_ ]torques",
            "(Ruminococcus|Mediterraneibacter)[_ ]gnavus",
            "Prevotella[_ ]stercorea",
            "Bacteroides[_ ]thetaiotaomicron",
            "Weissella[_ ]cibaria",
            "Fusicatenibacter[_ ]saccharivorans",
            "^Porphyromonas([_ ]|$)",
            "^Bifidobacterium([_ ]|$)",
            "Streptococcus[_ ]salivarius",
            "Streptococcus[_ ]sanguinis",
            "^Veillonella([_ ]|$)"
          ),
          collapse = "|"
        ),
        ignore_case = TRUE
      )
    ),
    
    # These taxa are excluded before correlation testing and therefore
    # cannot re-enter through either the priority or general-species route.
    species_excluded = str_detect(
      Species,
      regex(
        "Weissella[_ ]confusa|Anaerostipes[_ ]hadrus",
        ignore_case = TRUE
      )
    ),
    
    species_unclear_name = str_detect(
      Species,
      regex(
        "GGB|SGB|CAG|UBA|MAG|uncultured|metagenome|_sp_|bacterium|oral_taxon",
        ignore_case = TRUE
      )
    )
  )

# Confirm the two forced taxa against the unfiltered, directly reconstructed
# 20-sample matrix. This separates a true missing taxonomy label from failure
# of the prevalence criterion and prevents misleading post-filter errors.
forced_species_present_in_input <- vapply(
  forced_species_patterns,
  function(x) any(str_detect(
    species_stat$Species,
    regex(paste0("^", x, "$"), ignore_case = TRUE)
  )),
  logical(1)
)

if (any(!forced_species_present_in_input)) {
  stop(
    "The following expected biological species names were absent from ",
    "s$Species before prevalence filtering: ",
    paste(
      forced_species_labels[!forced_species_present_in_input],
      collapse = ", "
    ),
    ". Check the upstream species.csv taxonomy labels.",
    call. = FALSE
  )
}

forced_species_prevalence_audit <- species_stat %>%
  filter(str_detect(Species, forced_species_regex)) %>%
  mutate(
    forced_display_order = match(
      str_replace_all(Species, "_", " "),
      forced_species_labels
    ),
    total_matched_n = nrow(col_rna),
    pCR_group_n = sum(as.character(col_rna$TRG_plot) == "pCR"),
    non_pCR_group_n = sum(as.character(col_rna$TRG_plot) == "non_pCR"),
    passes_prevalence_0.30 =
      species_prevalence >= 0.30 & species_nonzero_n >= 6
  ) %>%
  arrange(forced_display_order) %>%
  select(
    Species,
    total_matched_n,
    pCR_group_n,
    non_pCR_group_n,
    species_nonzero_n,
    species_detected_pCR_n,
    species_detected_non_pCR_n,
    species_prevalence,
    species_mean_pCR,
    species_mean_non_pCR,
    species_mean_abundance,
    species_logFC,
    species_logFC_direction,
    species_p,
    passes_prevalence_0.30
  )

message("Forced-species prevalence audit before filtering:")
print(forced_species_prevalence_audit, n = Inf)

write_csv(
  forced_species_prevalence_audit,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
    "Fig4C_v11_forced_species_prevalence_audit.csv"
  )
)

if (any(!forced_species_prevalence_audit$passes_prevalence_0.30)) {
  warning(
    "At least one forced butyrate producer had prevalence < 0.30 in the ",
    "directly reconstructed 20-sample matrix. See the forced-species audit ",
    "CSV. The script will continue without substituting species identities."
  )
}

species_stat <- species_stat %>%
  filter(
    !species_excluded,
    species_prevalence >= 0.30,
    species_nonzero_n >= 6,
    !species_unclear_name | species_priority
  ) %>%
  arrange(
    desc(species_priority),
    species_p,
    desc(species_prevalence),
    desc(species_mean_abundance),
    desc(abs(species_logFC))
  )

if (nrow(species_stat) == 0) {
  stop(
    "No species remained after the prevalence and name-quality filters.",
    call. = FALSE
  )
}

target_species_direction_audit <- species_stat %>%
  filter(
    str_detect(
      Species,
      regex(
        paste(
          c("Bacteroides[_ ]fragilis", forced_species_patterns),
          collapse = "|"
        ),
        ignore_case = TRUE
      )
    )
  ) %>%
  left_join(species_name_map, by = "Species") %>%
  mutate(
    pCR_group_n = sum(as.character(col_rna$TRG_plot) == "pCR"),
    non_pCR_group_n = sum(as.character(col_rna$TRG_plot) == "non_pCR")
  ) %>%
  select(
    source_feature,
    source_feature_index,
    Species,
    pCR_group_n,
    non_pCR_group_n,
    species_detected_pCR_n,
    species_detected_non_pCR_n,
    species_mean_pCR,
    species_mean_non_pCR,
    species_logFC,
    species_logFC_direction,
    species_p,
    species_FDR,
    species_prevalence,
    species_mean_abundance
  )

print(target_species_direction_audit, n = Inf)

write_csv(
  target_species_direction_audit,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
    "Fig4C_v11_target_species_direction_audit.csv"
  )
)

# Collapse near-identical abundance-rank profiles among non-priority species.
# Priority taxa remain separate because species-level biological interpretation
# is retained only when each taxon independently shows a significant host-gene
# association under the p < 0.05 / repeated p < 0.10 selection rule.
species_rank_mat <- vapply(
  species_stat$Species,
  function(sp) rank(species_int[, sp], ties.method = "average"),
  numeric(nrow(species_int))
)

if (is.null(dim(species_rank_mat))) {
  species_rank_mat <- matrix(
    species_rank_mat,
    ncol = 1,
    dimnames = list(rownames(species_int), species_stat$Species)
  )
} else {
  rownames(species_rank_mat) <- rownames(species_int)
  colnames(species_rank_mat) <- species_stat$Species
}

species_keep <- character(0)
species_rank_duplicate_log <- tibble(
  removed_species = character(0),
  retained_species = character(0),
  rank_similarity = numeric(0)
)

for (sp in species_stat$Species) {
  if (sd(species_rank_mat[, sp], na.rm = TRUE) == 0) {
    next
  }
  
  if (
    species_stat$species_priority[match(sp, species_stat$Species)] ||
    length(species_keep) == 0
  ) {
    species_keep <- c(species_keep, sp)
    next
  }
  
  rank_similarity <- vapply(
    species_keep,
    function(ref_sp) {
      suppressWarnings(
        cor(
          species_rank_mat[, sp],
          species_rank_mat[, ref_sp],
          method = "spearman",
          use = "complete.obs"
        )
      )
    },
    numeric(1)
  )
  
  rank_similarity[!is.finite(rank_similarity)] <- 0
  
  if (max(rank_similarity) >= 0.995) {
    species_rank_duplicate_log <- bind_rows(
      species_rank_duplicate_log,
      tibble(
        removed_species = sp,
        retained_species = species_keep[which.max(rank_similarity)],
        rank_similarity = max(rank_similarity)
      )
    )
  } else {
    species_keep <- c(species_keep, sp)
  }
}

species_stat <- species_stat %>%
  filter(Species %in% species_keep) %>%
  slice_head(n = 200)

if (nrow(species_stat) == 0) {
  stop("No nonredundant species remained for correlation analysis.", call. = FALSE)
}

#-----------------------------------------------------------------#
# 6. Host pathway annotation for the fixed gene panel
#-----------------------------------------------------------------#
#
# GSEA categories are used only to annotate and group the fixed host genes.
# Genes not represented in the selected GSEA categories remain in the panel
# under "Other selected host response".
#-----------------------------------------------------------------#

gsea_plot_df <- gsea_plot_df %>%
  filter(
    pass_display_filter,
    padj < 0.10
  ) %>%
  mutate(leadingEdge = as.character(leadingEdge)) %>%
  rowwise() %>%
  mutate(
    n_leadingEdge_in_vst = sum(
      str_split(ifelse(is.na(leadingEdge), "", leadingEdge), ";")[[1]] %in%
        rownames(vst_rna)
    )
  ) %>%
  ungroup() %>%
  filter(n_leadingEdge_in_vst >= 5) %>%
  arrange(
    padj,
    desc(abs(NES)),
    desc(n_leadingEdge_in_vst)
  )

if (nrow(gsea_plot_df) == 0) {
  stop("No usable GSEA category remained for host annotation.", call. = FALSE)
}

gsea_plot_df$select_for_heatmap <- FALSE
gsea_plot_df$new_gene_fraction <- NA_real_

for (i in seq_len(nrow(gsea_plot_df))) {
  if (sum(gsea_plot_df$select_for_heatmap) >= 10) {
    next
  }
  
  if (sum(gsea_plot_df$select_for_heatmap) == 0) {
    gsea_plot_df$select_for_heatmap[i] <- TRUE
    gsea_plot_df$new_gene_fraction[i] <- 1
  } else {
    gsea_plot_df$new_gene_fraction[i] <- length(
      setdiff(
        str_split(gsea_plot_df$leadingEdge[i], ";")[[1]],
        unique(
          unlist(
            str_split(
              gsea_plot_df$leadingEdge[gsea_plot_df$select_for_heatmap],
              ";"
            )
          )
        )
      )
    ) / max(
      1,
      length(str_split(gsea_plot_df$leadingEdge[i], ";")[[1]])
    )
    
    if (
      is.na(gsea_plot_df$new_gene_fraction[i]) ||
      gsea_plot_df$new_gene_fraction[i] >= 0.30 ||
      sum(gsea_plot_df$select_for_heatmap) < 6
    ) {
      gsea_plot_df$select_for_heatmap[i] <- TRUE
    }
  }
}

gsea_plot_df <- gsea_plot_df %>%
  filter(select_for_heatmap) %>%
  slice_head(n = 10)

selected_category_gene <- gsea_plot_df %>%
  select(
    Host_axis,
    pathway_label,
    pathway_id,
    gs_name,
    NES,
    padj,
    leadingEdge
  ) %>%
  separate_rows(leadingEdge, sep = ";") %>%
  transmute(
    Host_axis,
    pathway_label,
    pathway_id,
    gs_name,
    pathway_NES = NES,
    pathway_padj = padj,
    Gene = str_squish(leadingEdge)
  ) %>%
  filter(
    Gene %in% target_genes,
    Gene != ""
  ) %>%
  distinct()

if (nrow(selected_category_gene) == 0) {
  selected_category_gene <- gene_sets_use_long %>%
    inner_join(
      gsea_plot_df %>%
        select(
          Host_axis,
          pathway_label,
          pathway_id,
          gs_name,
          NES,
          padj
        ) %>%
        distinct(),
      by = c("pathway_id", "gs_name", "Host_axis")
    ) %>%
    transmute(
      Host_axis,
      pathway_label,
      pathway_id,
      gs_name,
      pathway_NES = NES,
      pathway_padj = padj,
      Gene
    ) %>%
    filter(Gene %in% target_genes) %>%
    distinct()
}

host_axis_order <- unique(as.character(gsea_plot_df$Host_axis))

if (
  length(host_axis_order) > 1 &&
  any(
    c(
      "Inflammatory leukocyte trafficking",
      "Leukocyte trafficking",
      "Integrin-mediated gut homing"
    ) %in% host_axis_order
  )
) {
  host_axis_order <- c(
    host_axis_order[1],
    intersect(
      c(
        "Inflammatory leukocyte trafficking",
        "Leukocyte trafficking",
        "Integrin-mediated gut homing"
      ),
      host_axis_order[-1]
    ),
    setdiff(
      host_axis_order[-1],
      c(
        "Inflammatory leukocyte trafficking",
        "Leukocyte trafficking",
        "Integrin-mediated gut homing"
      )
    )
  )
}

host_category_summary <- selected_category_gene %>%
  arrange(
    pathway_padj,
    desc(abs(pathway_NES))
  ) %>%
  group_by(Gene) %>%
  summarise(
    Host_axis = first(as.character(Host_axis)),
    primary_pathway_label = first(pathway_label),
    n_selected_pathways = n_distinct(pathway_id),
    best_pathway_padj = min(pathway_padj, na.rm = TRUE),
    max_abs_pathway_NES = max(abs(pathway_NES), na.rm = TRUE),
    .groups = "drop"
  )

host_stat <- deg_rna %>%
  filter(Gene %in% target_genes) %>%
  arrange(match(Gene, target_genes)) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  transmute(
    Gene,
    host_logFC = log2FoldChange,
    host_DESeq2_p = pval_DESeq2,
    host_DESeq2_FDR = FDR_DESeq2,
    host_direction = direction
  ) %>%
  left_join(host_category_summary, by = "Gene") %>%
  mutate(
    Host_axis = ifelse(
      is.na(Host_axis),
      "Other selected host response",
      Host_axis
    ),
    Host_axis = factor(
      Host_axis,
      levels = c(host_axis_order, "Other selected host response")
    ),
    n_selected_pathways = replace_na(n_selected_pathways, 0L),
    best_pathway_padj = replace_na(best_pathway_padj, 1),
    max_abs_pathway_NES = replace_na(max_abs_pathway_NES, 0),
    host_neglog10p = -log10(
      pmax(host_DESeq2_p, .Machine$double.xmin)
    ),
    host_neglog10p_capped = pmin(host_neglog10p, 10),
    host_DESeq2_p_star = p_to_star(host_DESeq2_p)
  )

if (nrow(host_stat) != 18) {
  stop(
    "Some fixed target genes were absent from deg_rna: ",
    paste(setdiff(target_genes, host_stat$Gene), collapse = ", "),
    call. = FALSE
  )
}

target_genes <- host_stat$Gene
gene_block_tbl <- gene_block_tbl %>%
  filter(Gene %in% target_genes)

#-----------------------------------------------------------------#
# 7. Spearman correlations between candidate species and target genes
#-----------------------------------------------------------------#

species_int <- log10(
  species_int[, species_stat$Species, drop = FALSE] + 1e-06
)

vst_rna <- t(
  vst_rna[
    target_genes,
    rownames(col_rna),
    drop = FALSE
  ]
)

rho_mat <- matrix(
  NA_real_,
  nrow = ncol(species_int),
  ncol = ncol(vst_rna),
  dimnames = list(colnames(species_int), colnames(vst_rna))
)

spearman_p_mat <- rho_mat
n_mat <- rho_mat

for (i in seq_len(nrow(rho_mat))) {
  for (j in seq_len(ncol(rho_mat))) {
    x <- species_int[, rownames(rho_mat)[i]]
    y <- vst_rna[, colnames(rho_mat)[j]]
    ok <- is.finite(x) & is.finite(y)
    
    n_mat[i, j] <- sum(ok)
    
    if (
      sum(ok) >= 8 &&
      sd(x[ok], na.rm = TRUE) > 0 &&
      sd(y[ok], na.rm = TRUE) > 0
    ) {
      spearman_result <- suppressWarnings(
        cor.test(
          x[ok],
          y[ok],
          method = "spearman",
          exact = FALSE
        )
      )
      
      rho_mat[i, j] <- unname(spearman_result$estimate)
      spearman_p_mat[i, j] <- spearman_result$p.value
    }
  }
}

cor_df_all <- expand_grid(
  Species = rownames(rho_mat),
  Gene = colnames(rho_mat)
) %>%
  mutate(
    n = n_mat[cbind(Species, Gene)],
    spearman_rho = rho_mat[cbind(Species, Gene)],
    spearman_p = spearman_p_mat[cbind(Species, Gene)],
    spearman_FDR = p.adjust(spearman_p, method = "BH"),
    spearman_p_star = p_to_star(spearman_p),
    spearman_sig = !is.na(spearman_p) & spearman_p < 0.05,
    spearman_trend = !is.na(spearman_p) & spearman_p < 0.10
  )

cor_df_sig <- cor_df_all %>%
  filter(spearman_sig) %>%
  arrange(
    spearman_p,
    desc(abs(spearman_rho))
  )

cor_df_screen <- cor_df_all %>%
  filter(spearman_trend) %>%
  arrange(
    spearman_p,
    desc(abs(spearman_rho))
  )

if (nrow(cor_df_screen) == 0) {
  stop(
    "No candidate species had Spearman p < 0.10 with the target genes.",
    call. = FALSE
  )
}

#-----------------------------------------------------------------#
# 8. Final species selection and matrix preparation
#-----------------------------------------------------------------#
#
# Selection and display use different criteria:
#
#   Selection:
#     1. Require at least one Spearman p < 0.05 correlation or at least two
#        Spearman p < 0.10 correlations across the fixed host panel.
#     2. Retain the two prespecified butyrate producers whenever they pass
#        prevalence >= 0.30, even when criterion 1 is not met.
#     3. Rank remaining species first by correlation evidence and block-level pattern
#        consistency. Priority status, microbial Wilcoxon p, prevalence, and
#        abundance are used only as later tie-breakers.
#
#   Display order:
#     Species are clustered globally from their Spearman-rho profiles across
#     all 18 host genes. Taxonomy is annotated after ordering and does not
#     define row blocks.
#-----------------------------------------------------------------#

host_stat_plot <- cor_df_sig %>%
  group_by(Gene) %>%
  summarise(
    gene_best_p = min(spearman_p, na.rm = TRUE),
    n_sig_pairs = n(),
    n_sig_species = n_distinct(Species),
    max_abs_spearman_rho = max(abs(spearman_rho), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  right_join(host_stat, by = "Gene") %>%
  mutate(
    gene_best_p = replace_na(gene_best_p, 1),
    n_sig_pairs = replace_na(n_sig_pairs, 0L),
    n_sig_species = replace_na(n_sig_species, 0L),
    max_abs_spearman_rho = replace_na(max_abs_spearman_rho, 0)
  )

species_pair_rank <- cor_df_all %>%
  group_by(Species) %>%
  summarise(
    species_best_p = min(c(spearman_p, 1), na.rm = TRUE),
    n_sig_pairs = sum(spearman_p < 0.05, na.rm = TRUE),
    n_p10_pairs = sum(spearman_p < 0.10, na.rm = TRUE),
    n_sig_genes = n_distinct(
      Gene[!is.na(spearman_p) & spearman_p < 0.05]
    ),
    n_p10_genes = n_distinct(
      Gene[!is.na(spearman_p) & spearman_p < 0.10]
    ),
    max_abs_spearman_rho = max(c(abs(spearman_rho), 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    cor_df_all %>%
      left_join(gene_block_tbl, by = "Gene") %>%
      filter(!is.na(Host_block)) %>%
      group_by(Species, Host_block) %>%
      summarise(
        block_sig_genes = sum(spearman_p < 0.05, na.rm = TRUE),
        block_p10_genes = sum(spearman_p < 0.10, na.rm = TRUE),
        block_mean_rho = mean(spearman_rho, na.rm = TRUE),
        block_sign_consistency = max(
          mean(spearman_rho >= 0, na.rm = TRUE),
          mean(spearman_rho <= 0, na.rm = TRUE)
        ),
        .groups = "drop"
      ) %>%
      mutate(
        block_abs_mean_rho = abs(block_mean_rho)
      ) %>%
      group_by(Species) %>%
      arrange(
        desc(block_sig_genes),
        desc(block_p10_genes),
        desc(block_sign_consistency),
        desc(block_abs_mean_rho),
        .by_group = TRUE
      ) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      transmute(
        Species,
        best_host_block = as.character(Host_block),
        max_sig_genes_in_block = block_sig_genes,
        max_p10_genes_in_block = block_p10_genes,
        block_mean_rho,
        block_abs_mean_rho,
        block_sign_consistency
      ),
    by = "Species"
  ) %>%
  left_join(species_stat, by = "Species") %>%
  mutate(
    forced_species = str_detect(Species, forced_species_regex),
    meets_correlation_rule =
      n_sig_genes >= 1 | n_p10_genes >= 2,
    selection_tier = case_when(
      forced_species ~ 0L,
      n_sig_genes >= 2 ~ 1L,
      n_sig_genes >= 1 ~ 2L,
      n_p10_genes >= 3 ~ 3L,
      TRUE ~ 4L
    )
  ) %>%
  filter(
    meets_correlation_rule | forced_species
  ) %>%
  arrange(
    desc(forced_species),
    selection_tier,
    desc(n_sig_genes),
    desc(n_p10_genes),
    desc(max_sig_genes_in_block),
    desc(max_p10_genes_in_block),
    desc(block_sign_consistency),
    desc(block_abs_mean_rho),
    species_best_p,
    desc(species_priority),
    species_p,
    desc(species_prevalence),
    desc(species_mean_abundance)
  )

# Both output versions are capped at 31 taxa. The included version places an
# eligible B. fragilis first so that it cannot be displaced solely by the row
# cap. The excluded version removes it and independently backfills the next
# ranked eligible species.
species_stat_plot_with_bfragilis <- species_pair_rank %>%
  mutate(
    Bfragilis_priority = str_detect(
      Species,
      regex("Bacteroides[_ ]fragilis", ignore_case = TRUE)
    )
  ) %>%
  arrange(
    desc(Bfragilis_priority),
    desc(forced_species),
    selection_tier,
    desc(n_sig_genes),
    desc(n_p10_genes),
    species_best_p
  ) %>%
  slice_head(n = 31) %>%
  select(-Bfragilis_priority)

species_stat_plot_without_bfragilis <- species_pair_rank %>%
  filter(
    !str_detect(
      Species,
      regex("Bacteroides[_ ]fragilis", ignore_case = TRUE)
    )
  ) %>%
  arrange(
    desc(forced_species),
    selection_tier,
    desc(n_sig_genes),
    desc(n_p10_genes),
    species_best_p
  ) %>%
  slice_head(n = 31)

forced_species_final <- species_stat$Species[
  str_detect(species_stat$Species, forced_species_regex)
]

if (
  any(!forced_species_final %in% species_stat_plot_with_bfragilis$Species) ||
    any(!forced_species_final %in% species_stat_plot_without_bfragilis$Species)
) {
  stop(
    paste(
      "At least one forced butyrate producer was lost during the 31-row cap.",
      "Check final species ranking."
    ),
    call. = FALSE
  )
}

if (
  nrow(species_stat_plot_with_bfragilis) == 0 ||
    nrow(species_stat_plot_without_bfragilis) == 0
) {
  stop(
    paste(
      "No species met the final correlation rule:",
      "at least one p < 0.05 or at least two p < 0.10."
    ),
    call. = FALSE
  )
}

if (any(str_detect(
  c(
    species_stat_plot_with_bfragilis$Species,
    species_stat_plot_without_bfragilis$Species
  ),
  regex("Weissella[_ ]confusa", ignore_case = TRUE)
))) {
  stop(
    "An explicitly excluded species entered the final heatmap.",
    call. = FALSE
  )
}

if (!any(str_detect(
  species_stat_plot_with_bfragilis$Species,
  regex("Bacteroides[_ ]fragilis", ignore_case = TRUE)
))) {
  warning(
    "Bacteroides fragilis did not meet the host-correlation inclusion rule. ",
    "The included version therefore cannot display it without overriding the ",
    "prespecified prevalence and correlation criteria."
  )
}

if (
  nrow(species_stat_plot_with_bfragilis) !=
    nrow(species_stat_plot_without_bfragilis)
) {
  warning(
    "The two variants contain different species counts because fewer than 31 ",
    "eligible replacement taxa were available."
  )
}

write_csv(
  species_pair_rank,
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5_species_block_coherence_ranking.csv"
)

write_csv(
  species_pair_rank %>%
    filter(
      !species_priority,
      n_sig_genes == 0,
      n_p10_genes >= 2,
      block_sign_consistency < 0.75
    ),
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5_species_weak_pattern_review.csv"
)

priority_species_check <- species_stat %>%
  filter(species_priority) %>%
  left_join(
    species_pair_rank %>%
      select(
        Species,
        species_best_p,
        n_sig_genes,
        n_p10_genes,
        forced_species,
        meets_correlation_rule,
        max_abs_spearman_rho
      ),
    by = "Species"
  ) %>%
  mutate(
    n_sig_genes = replace_na(n_sig_genes, 0L),
    n_p10_genes = replace_na(n_p10_genes, 0L),
    forced_species = replace_na(forced_species, FALSE),
    meets_correlation_rule = replace_na(meets_correlation_rule, FALSE),
    eligible_with_target_panel =
      meets_correlation_rule | forced_species,
    included_with_B_fragilis =
      Species %in% species_stat_plot_with_bfragilis$Species,
    included_without_B_fragilis =
      Species %in% species_stat_plot_without_bfragilis$Species
  ) %>%
  arrange(
    desc(eligible_with_target_panel),
    desc(n_sig_genes),
    desc(n_p10_genes),
    species_best_p,
    desc(species_prevalence),
    desc(species_mean_abundance)
  )

print(priority_species_check, n = 100)

priority_species_check %>%
  filter(
    str_detect(
      Species,
      regex(
        paste(
          c(
            "Weissella[_ ]cibaria",
            "Fusicatenibacter[_ ]saccharivorans",
            "^Porphyromonas([_ ]|$)",
            "^Bifidobacterium([_ ]|$)",
            "Streptococcus[_ ]salivarius",
            "Streptococcus[_ ]sanguinis",
            "^Veillonella([_ ]|$)"
          ),
          collapse = "|"
        ),
        ignore_case = TRUE
      )
    )
  ) %>%
  select(
    Species,
    species_prevalence,
    species_nonzero_n,
    species_mean_abundance,
    species_logFC,
    species_p,
    n_sig_genes,
    species_best_p,
    included_with_B_fragilis,
    included_without_B_fragilis
  ) %>%
  print(n = 100)

print(
  cor_df_sig %>%
    filter(Species %in% priority_species_check$Species) %>%
    left_join(
      species_stat %>%
        select(
          Species,
          species_prevalence,
          species_nonzero_n,
          species_mean_abundance,
          species_logFC,
          species_p
        ),
      by = "Species"
    ) %>%
    arrange(Species, spearman_p),
  n = 200
)

write_csv(
  priority_species_check,
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5_priority_species_screen.csv"
)

write_csv(
  cor_df_sig %>%
    filter(Species %in% priority_species_check$Species) %>%
    left_join(
      species_stat %>%
        select(
          Species,
          species_prevalence,
          species_nonzero_n,
          species_mean_abundance,
          species_logFC,
          species_p,
          species_FDR
        ),
      by = "Species"
    ) %>%
    arrange(Species, spearman_p),
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5_priority_species_significant_pairs.csv"
)

write_csv(
  species_pair_rank %>%
    filter(
      species_logFC > 0,
      n_sig_genes >= 2
    ) %>%
    arrange(
      desc(n_sig_genes),
      species_best_p,
      species_p,
      desc(species_prevalence),
      desc(species_mean_abundance)
    ),
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5_pCR_enriched_species_with_multiple_host_correlations.csv"
)

write_csv(
  species_stat_plot_with_bfragilis,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
    "Fig4C_v11_species_host_heatmap_final_species_Bfragilis_included.csv"
  )
)

write_csv(
  species_stat_plot_without_bfragilis,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
    "Fig4C_v11_species_host_heatmap_final_species_Bfragilis_excluded.csv"
  )
)

write_csv(
  host_stat_plot,
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5_species_host_heatmap_final_genes.csv"
)

message(
  "Final heatmaps: ",
  nrow(species_stat_plot_with_bfragilis),
  " and ",
  nrow(species_stat_plot_without_bfragilis),
  " species × ",
  nrow(host_stat_plot),
  " fixed genes"
)

rho_mat_all <- rho_mat
spearman_p_mat_all <- spearman_p_mat
n_mat_all <- n_mat

species_plot_variants <- list(
  Bfragilis_included = species_stat_plot_with_bfragilis,
  Bfragilis_excluded = species_stat_plot_without_bfragilis
)

figure_files <- c(
  Bfragilis_included =
    "figures/Fig4C_v11_species_host_DEG_Spearman_taxonomy_Bfragilis_included.svg",
  Bfragilis_excluded =
    "figures/Fig4C_v11_species_host_DEG_Spearman_taxonomy_Bfragilis_excluded.svg"
)

taxonomy_display_reference <- species_taxonomy %>%
  filter(
    Species %in% unique(
      unlist(
        lapply(species_plot_variants, function(x) x$Species),
        use.names = FALSE
      )
    )
  )

if (all(taxonomy_display_reference$Family == "Unclassified")) {
  stop(
    "No displayed species could be matched to an f__ family in ",
    "input/species.csv. Check whether the species labels and clade_name ",
    "strings were generated from the same taxonomy table.",
    call. = FALSE
  )
}

if (any(taxonomy_display_reference$Family == "Unclassified")) {
  warning(
    "Family was not resolved for: ",
    paste(
      taxonomy_display_reference$Species[
        taxonomy_display_reference$Family == "Unclassified"
      ],
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (any(taxonomy_display_reference$Phylum == "Unclassified")) {
  warning(
    "Phylum was not resolved for: ",
    paste(
      taxonomy_display_reference$Species[
        taxonomy_display_reference$Phylum == "Unclassified"
      ],
      collapse = ", "
    ),
    call. = FALSE
  )
}

family_display_counts <- taxonomy_display_reference %>%
  count(Family, name = "n_species")

singleton_family_exceptions <- c(
  "Fusobacteriaceae",
  "Bifidobacteriaceae",
  "Oscillospiraceae",
  "Veillonellaceae"
)

singleton_families <- family_display_counts %>%
  filter(
    n_species == 1,
    !Family %in% c(singleton_family_exceptions, "Unclassified")
  ) %>%
  pull(Family)

species_taxonomy <- species_taxonomy %>%
  mutate(
    Family = if_else(
      Family_source %in% singleton_families,
      "Others",
      Family_source
    )
  )

taxonomy_display_reference <- taxonomy_display_reference %>%
  mutate(
    Family = if_else(
      Family_source %in% singleton_families,
      "Others",
      Family_source
    )
  )

write_csv(
  species_taxonomy,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
    "Fig4C_v11_species_taxonomy_map.csv"
  )
)

phylum_levels_all <- c(
  sort(setdiff(unique(taxonomy_display_reference$Phylum), "Unclassified")),
  intersect("Unclassified", unique(taxonomy_display_reference$Phylum))
)

family_levels_all <- c(
  sort(
    setdiff(
      unique(taxonomy_display_reference$Family),
      c("Unclassified", "Others")
    )
  ),
  intersect("Unclassified", unique(taxonomy_display_reference$Family)),
  intersect("Others", unique(taxonomy_display_reference$Family))
)

phylum_cols_all <- c(
  "Actinobacteria" = "#C2D3EA",
  "Actinomycetota" = "#C2D3EA",
  "Bacteroidota" = "#F1CCAD",
  "Bacteroidetes" = "#F1CCAD",
  "Firmicutes" = "#C8DFC0",
  "Bacillota" = "#C8DFC0",
  "Fusobacteria" = "#E8BBB8",
  "Fusobacteriota" = "#E8BBB8",
  "Proteobacteria" = "#D6C7E7",
  "Pseudomonadota" = "#D6C7E7",
  "Unclassified" = "#EEEEEE"
)[phylum_levels_all]

if (anyNA(phylum_cols_all)) {
  phylum_cols_all[is.na(phylum_cols_all)] <- rep(
    c("#C7D8E6", "#E9D2BB", "#CADFCA", "#DDCAE3", "#E7D9AC"),
    length.out = sum(is.na(phylum_cols_all))
  )
}

family_cols_all <- c(
  "Lachnospiraceae" = "#7FA7C9",
  "Bacteroidaceae" = "#E5A36B",
  "Streptococcaceae" = "#8FBE88",
  "Actinomycetaceae" = "#B89ABD",
  "Fusobacteriaceae" = "#D77C78",
  "Oscillospiraceae" = "#C0A27A",
  "Bifidobacteriaceae" = "#78B7B1",
  "Veillonellaceae" = "#A99BD2",
  "Prevotellaceae" = "#C9BE72",
  "Unclassified" = "#F1F1F1",
  "Others" = "#E3E3E3"
)[family_levels_all]

if (anyNA(family_cols_all)) {
  family_cols_all[is.na(family_cols_all)] <- rep(
    c("#91B6CF", "#E5B486", "#9FC59A", "#C3A7C5", "#D99A96"),
    length.out = sum(is.na(family_cols_all))
  )
}

plot_results <- vector("list", length(species_plot_variants))
names(plot_results) <- names(species_plot_variants)

for (variant_name in names(species_plot_variants)) {
  species_stat_plot <- species_plot_variants[[variant_name]]
  rho_mat <- rho_mat_all
  spearman_p_mat <- spearman_p_mat_all
  n_mat <- n_mat_all

rho_mat <- rho_mat[
  species_stat_plot$Species,
  host_stat_plot$Gene,
  drop = FALSE
]

spearman_p_mat <- spearman_p_mat[
  species_stat_plot$Species,
  host_stat_plot$Gene,
  drop = FALSE
]

n_mat <- n_mat[
  species_stat_plot$Species,
  host_stat_plot$Gene,
  drop = FALSE
]

spearman_star_mat <- matrix(
  p_to_star(as.vector(spearman_p_mat)),
  nrow = nrow(spearman_p_mat),
  ncol = ncol(spearman_p_mat),
  dimnames = dimnames(spearman_p_mat)
)

cor_df <- expand_grid(
  Species = rownames(rho_mat),
  Gene = colnames(rho_mat)
) %>%
  mutate(
    n = n_mat[cbind(Species, Gene)],
    spearman_rho = rho_mat[cbind(Species, Gene)],
    spearman_p = spearman_p_mat[cbind(Species, Gene)],
    spearman_FDR = p.adjust(spearman_p, method = "BH"),
    spearman_p_star = p_to_star(spearman_p),
    spearman_sig = !is.na(spearman_p) & spearman_p < 0.05,
    spearman_trend = !is.na(spearman_p) & spearman_p < 0.10
  )

#-----------------------------------------------------------------#
# 9. Row and column ordering
#-----------------------------------------------------------------#
#
# Microbial rows:
#   Cluster all selected taxa together from their Spearman-rho profiles across
#   the predefined host panel. Phylum and family remain side annotations but
#   do not define row order or row blocks.
#
# Host columns:
#   Keep the requested four-block biological sequence and cluster genes within
#   each block according to their microbial-correlation profiles.
#-----------------------------------------------------------------#

species_stat_plot <- species_stat_plot %>%
  mutate(
    Response_enrichment = factor(
      if_else(
        species_logFC > 0,
        "pCR enriched",
        "non-pCR enriched"
      ),
      levels = c("non-pCR enriched", "pCR enriched")
    )
  )

rho_for_clustering <- rho_mat
rho_for_clustering[is.na(rho_for_clustering)] <- 0

if (nrow(rho_for_clustering) <= 1) {
  row_order <- rownames(rho_for_clustering)
} else {
  row_order <- rownames(rho_for_clustering)[
    hclust(dist(rho_for_clustering, method = "euclidean"))$order
  ]
}

col_order <- unlist(
  lapply(
    levels(gene_block_tbl$Host_block),
    function(host_block) {
      genes_in_block <- gene_block_tbl$Gene[
        gene_block_tbl$Host_block == host_block &
          gene_block_tbl$Gene %in% host_stat_plot$Gene
      ]

      if (length(genes_in_block) <= 1) {
        return(genes_in_block)
      }

      gene_correlation_profiles <- t(
        rho_mat[, genes_in_block, drop = FALSE]
      )
      gene_correlation_profiles[is.na(gene_correlation_profiles)] <- 0

      genes_in_block[
        hclust(dist(gene_correlation_profiles))$order
      ]
    }
  ),
  use.names = FALSE
)

rho_mat <- rho_mat[row_order, col_order, drop = FALSE]
spearman_p_mat <- spearman_p_mat[row_order, col_order, drop = FALSE]
spearman_star_mat <- spearman_star_mat[row_order, col_order, drop = FALSE]

species_stat_plot <- species_stat_plot %>%
  mutate(
    Species = factor(Species, levels = row_order)
  ) %>%
  arrange(Species) %>%
  mutate(
    Species = as.character(Species)
  )

if (!identical(species_stat_plot$Species, rownames(rho_mat))) {
  stop(
    "Microbial row order differs between the correlation matrix and annotations.",
    call. = FALSE
  )
}

species_taxonomy_plot <- species_taxonomy %>%
  filter(Species %in% species_stat_plot$Species) %>%
  arrange(match(Species, species_stat_plot$Species))

if (!identical(species_taxonomy_plot$Species, species_stat_plot$Species)) {
  stop(
    "Microbial taxonomy annotation is misaligned with heatmap rows.",
    call. = FALSE
  )
}

host_stat_plot <- host_stat_plot %>%
  mutate(
    Gene = factor(Gene, levels = col_order)
  ) %>%
  arrange(Gene) %>%
  mutate(
    Gene = as.character(Gene)
  )

column_split <- factor(
  gene_block_tbl$Host_block[
    match(col_order, gene_block_tbl$Gene)
  ],
  levels = levels(gene_block_tbl$Host_block)
)

row_split <- NULL

#-----------------------------------------------------------------#
# 10. Host selected-category annotation matrix
#-----------------------------------------------------------------#

selected_category_gene_plot <- selected_category_gene %>%
  filter(Gene %in% col_order) %>%
  mutate(
    pathway_label_clean = pathway_label,
    pathway_label_clean = str_replace_all(
      pathway_label_clean,
      "TNF-NF-kB / IL6-JAK-STAT signaling",
      "TNF/NF-kB & IL6/JAK-STAT"
    ),
    pathway_label_clean = str_replace_all(
      pathway_label_clean,
      "Mucus/goblet-cell program",
      "Mucus/goblet-cell"
    ),
    pathway_label_clean = str_replace_all(
      pathway_label_clean,
      "Th1 / Th2 / Th17 / Treg lineage",
      "T-cell lineage"
    ),
    pathway_label_clean = str_replace_all(
      pathway_label_clean,
      "Lymphocyte-mediated immunity",
      "Lymphocyte immunity"
    ),
    pathway_label_clean = str_replace_all(
      pathway_label_clean,
      "Inflammatory leukocyte trafficking",
      "Leukocyte trafficking"
    ),
    pathway_label_clean = str_replace_all(
      pathway_label_clean,
      "Oxidative stress/ROS response",
      "Oxidative stress / ROS"
    ),
    pathway_label_clean = str_replace_all(
      pathway_label_clean,
      "Colorectal adenoma signature",
      "Adenoma signature"
    ),
    pathway_label_clean = ifelse(
      str_detect(
        pathway_label_clean,
        regex("^Integrin-mediated gut homing$", ignore_case = TRUE)
      ),
      "Integrin-mediated gut homing",
      str_wrap(pathway_label_clean, width = 26)
    )
  ) %>%
  distinct(
    Host_axis,
    pathway_label,
    pathway_label_clean,
    pathway_id,
    pathway_padj,
    Gene
  )

if (nrow(selected_category_gene_plot) == 0) {
  category_mat <- matrix(
    "Included",
    nrow = 1,
    ncol = length(col_order),
    dimnames = list("Selected host genes", col_order)
  )
} else {
  category_mat <- selected_category_gene_plot %>%
    mutate(included = "Included") %>%
    dplyr::select(pathway_label_clean, Gene, included) %>%
    distinct() %>%
    pivot_wider(
      id_cols = pathway_label_clean,
      names_from = Gene,
      values_from = included,
      values_fill = "Not included"
    ) %>%
    as.data.frame(check.names = FALSE)
  
  rownames(category_mat) <- category_mat$pathway_label_clean
  category_mat$pathway_label_clean <- NULL
  
  for (g in setdiff(col_order, colnames(category_mat))) {
    category_mat[, g] <- "Not included"
  }
  
  category_mat <- as.matrix(category_mat[, col_order, drop = FALSE])
  category_mat <- category_mat[
    rowSums(category_mat == "Included") > 0,
    ,
    drop = FALSE
  ]
}

category_signature <- apply(
  category_mat == "Included",
  1,
  paste0,
  collapse = ""
)

if (any(duplicated(category_signature))) {
  category_dup_tbl <- tibble(
    pathway_label_clean = rownames(category_mat),
    category_signature = category_signature
  ) %>%
    left_join(
      selected_category_gene_plot %>%
        distinct(pathway_label_clean, pathway_padj),
      by = "pathway_label_clean"
    ) %>%
    group_by(category_signature) %>%
    arrange(pathway_padj, pathway_label_clean, .by_group = TRUE) %>%
    mutate(keep = row_number() == 1) %>%
    ungroup()
  
  category_mat <- category_mat[
    category_dup_tbl$pathway_label_clean[category_dup_tbl$keep],
    ,
    drop = FALSE
  ]
}

#-----------------------------------------------------------------#
# 11. Heatmap colors
#-----------------------------------------------------------------#

max_abs_rho <- max(abs(rho_mat), na.rm = TRUE)
max_abs_rho <- max(0.30, min(1, max_abs_rho))

cor_col_fun <- colorRamp2(
  c(-max_abs_rho, 0, max_abs_rho),
  c("#2166AC", "#FFFFFF", "#B2182B")
)

host_logfc_limit <- max(
  0.50,
  quantile(abs(host_stat_plot$host_logFC), 0.95, na.rm = TRUE)
)

host_logfc_col_fun <- colorRamp2(
  c(-host_logfc_limit, 0, host_logfc_limit),
  c(group_cols["non_pCR"], "#F7F7F7", group_cols["pCR"])
)

host_p_limit <- max(
  1,
  min(
    10,
    quantile(host_stat_plot$host_neglog10p_capped, 0.95, na.rm = TRUE)
  )
)

host_p_col_fun <- colorRamp2(
  c(0, host_p_limit / 2, host_p_limit),
  c("#F7FBFF", "#7AA6D1", "#2C5AA0")
)

phylum_levels <- phylum_levels_all[
  phylum_levels_all %in% species_taxonomy_plot$Phylum
]

family_levels <- family_levels_all[
  family_levels_all %in% species_taxonomy_plot$Family
]

phylum_cols <- phylum_cols_all[phylum_levels]
family_cols <- family_cols_all[family_levels]

#-----------------------------------------------------------------#
# 12. Heatmap annotations
#-----------------------------------------------------------------#

# The two taxonomy strips retain the same 4.5-mm width and 0.8-mm gap as the
# two microbial response-statistic strips in v8, preserving heatmap geometry.
# All phyla and families represented in a variant are shown in its legends.
ha_microbe <- rowAnnotation(
  Phylum = factor(species_taxonomy_plot$Phylum, levels = phylum_levels),
  Family = factor(species_taxonomy_plot$Family, levels = family_levels),
  col = list(
    Phylum = phylum_cols,
    Family = family_cols
  ),
  na_col = "#F2F2F2",
  simple_anno_size = unit(4.5, "mm"),
  show_annotation_name = TRUE,
  annotation_name_side = "bottom",
  # Rotate the microbial annotation labels vertically, matching the
  # orientation of the host-gene names below the heatmap.
  annotation_name_rot = 90,
  annotation_name_gp = gpar(fontsize = 8.2),
  gap = unit(0.8, "mm"),
  border = FALSE,
  gp = gpar(col = "white", lwd = 0.2),
  annotation_legend_param = list(
    Phylum = list(
      title = "Phylum",
      at = phylum_levels,
      labels = phylum_levels,
      labels_gp = gpar(fontsize = 9),
      title_gp = gpar(fontsize = 10, fontface = "bold")
    ),
    Family = list(
      title = "Family",
      at = family_levels,
      labels = family_levels,
      labels_gp = gpar(fontsize = 9),
      title_gp = gpar(fontsize = 10, fontface = "bold")
    )
  )
)

ht_category <- Heatmap(
  category_mat,
  name = "Host category",
  col = c("Not included" = "#F2F2F2", "Included" = "#303030"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  column_split = column_split,
  column_title = NULL,
  column_gap = unit(1.5, "mm"),
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 9),
  row_names_max_width = unit(58, "mm"),
  # Compress the pathway-membership block vertically so that the main
  # correlation matrix occupies more of the panel height.
  height = unit(max(15, nrow(category_mat) * 4.0), "mm"),
  rect_gp = gpar(col = "white", lwd = 0.2),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 9),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  )
)

ht_host_logfc <- Heatmap(
  matrix(
    host_stat_plot$host_logFC,
    nrow = 1,
    dimnames = list("Host LogFC", host_stat_plot$Gene)
  )[, col_order, drop = FALSE],
  name = "Host LogFC",
  col = host_logfc_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  column_split = column_split,
  column_title = NULL,
  column_gap = unit(1.5, "mm"),
  show_column_names = FALSE,
  # Match the strip height to the 4.5-mm microbial side annotations.
  # Keep the host annotation label horizontal for readability.
  row_names_rot = 0,
  row_names_gp = gpar(fontsize = 8.2),
  height = unit(4.5, "mm"),
  rect_gp = gpar(col = "white", lwd = 0.2),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 9),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  )
)

ht_host_p <- Heatmap(
  matrix(
    host_stat_plot$host_neglog10p_capped,
    nrow = 1,
    dimnames = list("Host DESeq2 -log10P", host_stat_plot$Gene)
  )[, col_order, drop = FALSE],
  name = "Host DESeq2 -log10P",
  col = host_p_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  column_split = column_split,
  column_title = NULL,
  column_gap = unit(1.5, "mm"),
  show_column_names = FALSE,
  # Match the strip height to the 4.5-mm microbial side annotations.
  # Keep the host annotation label horizontal for readability.
  row_names_rot = 0,
  row_names_gp = gpar(fontsize = 8.2),
  height = unit(4.5, "mm"),
  rect_gp = gpar(col = "white", lwd = 0.2),
  cell_fun = function(j, i, x, y, width, height, fill) {
    star <- host_stat_plot$host_DESeq2_p_star[match(col_order[j], host_stat_plot$Gene)]
    p_value_tile <- host_stat_plot$host_neglog10p_capped[
      match(col_order[j], host_stat_plot$Gene)
    ]
    
    if (nzchar(star)) {
      grid.text(
        star,
        x,
        y,
        gp = gpar(
          col = ifelse(p_value_tile >= 1.5, "white", "black"),
          fontsize = 7.5,
          fontface = "bold"
        )
      )
    }
  },
  heatmap_legend_param = list(
    title = "Host DESeq2 -Log10 P",
    labels_gp = gpar(fontsize = 9),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  )
)

ht_main <- Heatmap(
  rho_mat,
  name = "Spearman rho",
  col = cor_col_fun,
  na_col = "#F2F2F2",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  row_split = row_split,
  row_gap = unit(1.5, "mm"),
  row_title = NULL,
  column_split = column_split,
  column_gap = unit(1.5, "mm"),
  right_annotation = ha_microbe,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "right",
  row_labels = stringr::str_replace_all(rownames(rho_mat), "_", " "),
  row_names_gp = gpar(fontsize = 8.2, fontface = "italic"),
  row_names_max_width = unit(50, "mm"),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 90,
  column_names_max_height = unit(62, "mm"),
  rect_gp = gpar(col = "white", lwd = 0.2),
  # Compress only the horizontal tile dimension to 80% of the previous width;
  # vertical tile size and all legend dimensions remain unchanged.
  width = unit(ncol(rho_mat) * 3.68, "mm"),
  height = unit(nrow(rho_mat) * 4.6 + 1.5, "mm"),
  column_title_gp = gpar(fontsize = 8.5, fontface = "bold"),
  column_title_rot = 0,
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 9),
    title_gp = gpar(fontsize = 10, fontface = "bold")
  ),
  layer_fun = function(j, i, x, y, width, height, fill) {
    stars <- spearman_star_mat[cbind(i, j)]
    rho_values <- rho_mat[cbind(i, j)]
    star_index <- !is.na(stars) & stars != ""
    
    if (any(star_index)) {
      grid.text(
        stars[star_index],
        x[star_index],
        y[star_index],
        gp = gpar(
          col = ifelse(abs(rho_values[star_index]) >= 0.65, "white", "black"),
          fontsize = 7.5,
          fontface = "bold"
        )
      )
    }
  }
)

p_star_lgd <- Legend(
  title = "P value",
  labels = c(
    "*  p < 0.05",
    "**  p < 0.01",
    "***  p < 0.001"
  ),
  graphics = list(
    function(x, y, w, h) {
      grid.text("*", x, y, gp = gpar(fontsize = 9, fontface = "bold"))
    },
    function(x, y, w, h) {
      grid.text("**", x, y, gp = gpar(fontsize = 9, fontface = "bold"))
    },
    function(x, y, w, h) {
      grid.text("***", x, y, gp = gpar(fontsize = 9, fontface = "bold"))
    }
  ),
  labels_gp = gpar(fontsize = 9),
  title_gp = gpar(fontsize = 10, fontface = "bold")
)

#-----------------------------------------------------------------#
# 13. Save SVG figure
#-----------------------------------------------------------------#

svglite(
  figure_files[[variant_name]],
  # Reduce the canvas by the same absolute width removed from the 18 heatmap
  # columns, preserving the available space for the right-side legends.
  width = 7.45,
  height = 12.6
)

draw(
  ht_category %v% ht_host_logfc %v% ht_host_p %v% ht_main,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(p_star_lgd),
  merge_legends = TRUE,
  padding = unit(c(1, 2, 1, 2), "mm")
)

dev.off()

plot_results[[variant_name]] <- list(
  species_stat_plot = species_stat_plot,
  species_taxonomy_plot = species_taxonomy_plot,
  host_stat_plot = host_stat_plot,
  selected_category_gene_plot = selected_category_gene_plot,
  category_mat = category_mat,
  cor_df = cor_df,
  rho_mat = rho_mat,
  spearman_p_mat = spearman_p_mat,
  spearman_star_mat = spearman_star_mat,
  row_order = row_order,
  col_order = col_order,
  row_split = row_split,
  column_split = column_split,
  phylum_cols = phylum_cols,
  family_cols = family_cols,
  output_file = figure_files[[variant_name]]
)
}

#-----------------------------------------------------------------#
# 14. Save essential objects for downstream scripts
#-----------------------------------------------------------------#

setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)

base::save(
  col_rna,
  species_int,
  species_name_map,
  species_taxonomy,
  taxonomy_display_reference,
  family_display_counts,
  singleton_family_exceptions,
  singleton_families,
  vst_rna,
  species_stat,
  forced_species_prevalence_audit,
  target_species_direction_audit,
  species_rank_duplicate_log,
  host_stat_plot,
  gsea_plot_df,
  selected_category_gene,
  priority_species_check,
  species_pair_rank,
  cor_df_all,
  cor_df_sig,
  cor_df_screen,
  rho_mat_all,
  spearman_p_mat_all,
  n_mat_all,
  species_stat_plot_with_bfragilis,
  species_stat_plot_without_bfragilis,
  plot_results,
  figure_files,
  phylum_cols_all,
  family_cols_all,
  group_cols,
  target_genes,
  gene_block_tbl,
  forced_species_labels,
  forced_species_patterns,
  host_axis_order,
  file = paste0(
    "host_RNAseq/results_clean_metabolite_host/",
    "Figure5_host_RNAseq/",
    "Fig4C_v11_species_host_heatmap_Spearman_taxonomy_Bfragilis_variants.RData"
  ),
  compress = FALSE
)

message(
  "Final matched RNA-microbiome sample check: n = ", nrow(col_rna),
  " (pCR = ", sum(as.character(col_rna$TRG_plot) == "pCR"),
  ", non-pCR = ", sum(as.character(col_rna$TRG_plot) == "non_pCR"),
  ")."
)

message(
  "Forced-species prevalence in the matched baseline subset: ",
  paste0(
    stringr::str_replace_all(
      forced_species_prevalence_audit$Species,
      "_",
      " "
    ),
    " = ",
    sprintf(
      "%.0f%%",
      100 * forced_species_prevalence_audit$species_prevalence
    ),
    " (", forced_species_prevalence_audit$species_nonzero_n, "/",
    forced_species_prevalence_audit$total_matched_n, ")",
    collapse = "; "
  ),
  "."
)

message(
  "Final species inclusion rule: >=1 correlation with p < 0.05 or ",
  ">=2 correlations with p < 0.10."
)

message(
  "Displayed correlation stars: * p < 0.05; ** p < 0.01; ",
  "*** p < 0.001."
)

message(
  "Displayed microbial taxonomy: ",
  length(unique(taxonomy_display_reference$Phylum)), " phyla and ",
  length(unique(taxonomy_display_reference$Family)),
  " families across the two heatmap variants; each legend is restricted to ",
  "the taxa present in that variant."
)

message(
  "Families represented by one displayed species were grouped as Others: ",
  ifelse(
    length(singleton_families) == 0,
    "none",
    paste(singleton_families, collapse = ", ")
  ),
  ". Singleton-family exceptions retained separately: ",
  paste(singleton_family_exceptions, collapse = ", "),
  "."
)

message(
  "Forced in both variants after prevalence filtering: ",
  paste(forced_species_final, collapse = ", "), "."
)
