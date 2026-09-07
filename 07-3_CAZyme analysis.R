# ==================== CAZyme analysis v3 ====================
# Verified changes retained from v2:
#   1) mechanistic family selection uses nominal P < 0.1 without depending on
#      an undefined object;
#   2) GH32 is excluded from the starch/alpha-glucan-associated capacity score;
#   3) GH32 is evaluated separately and within the fructan/inulin-associated set;
#   4) Panel F has a separate GH32-excluded starch/alpha-glucan version;
#   5) modification times of saved SVG files are refreshed after ggsave().
#
# Fixed in v3:
#   6) subfamily diagnostic output uses the existing median_TPM_pCR and
#      median_TPM_non_pCR columns rather than undefined median_pCR and
#      median_non_pCR columns.


# ==================== 0. Setting ====================

# ---------- 0-1. Directory setting ----------

# ---------- 0-1. Directory setting ----------

rm(list = ls())
setwd("D:/2-연구/2-CRC metagenomics/")
options(stringsAsFactors = FALSE)
set.seed(123)

library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggforce)
library(purrr)
library(broom)
library(ggrepel)
library(vegan)
library(ComplexHeatmap)
library(circlize)
library(ggbeeswarm)



# Group color used consistently throughout all CAZyme analyses.
# CR in TRG_1 is renamed as pCR, and nonCR is renamed as non_pCR for plotting.

group_colors <- c("pCR" = "#4FAE9A", "non_pCR" = "#DE7872")

# Display labels for x-axis.
group_labels <- c("pCR" = "pCR", "non_pCR" = "non-pCR")

# Nominal significance cutoff used for feature selection.
# For this exploratory CAZyme analysis, p < 0.05 is used to select features for visualization.
sig_p_cutoff <- 0.05

# Nominal cutoff used to retain mechanistically linked CAZyme families.
mechanistic_p_cutoff <- 0.1

# ---------- 0-2. Rdata ----------
load("input/metabolite_metadata.RData")


# ==================== 1. Importing ====================

# ---------- 1-1. Metadata: Before only ----------
# Primary analysis is restricted to pretreatment baseline samples.
# Original response labels in m$TRG_1 are CR / nonCR.
# These are converted to pCR / non_pCR for cleaner downstream handling.

m_before <- m %>%
  dplyr::filter(
    TNT == "Before",
    TRG_1 %in% c("CR", "nonCR")
  ) %>%
  dplyr::mutate(
    Response = dplyr::case_when(
      TRG_1 == "CR" ~ "pCR",
      TRG_1 == "nonCR" ~ "non_pCR"
    ),
    Response = factor(Response, levels = c("pCR", "non_pCR"))
  ) %>%
  dplyr::arrange(Response, SampleID)

table(m_before$Response)

m_before %>%
  dplyr::select(SampleID, SNU_ID, SubjectID, TNT, TRG, TRG_score, TRG_1, Response) %>%
  head()


# ---------- 1-2. KO import ----------
# Kept as-is for future KO-linked carbohydrate metabolism analyses.

ko <- read_tsv("260224 final script/Input/merged_genefamilies_KO_named.tsv") %>% 
  dplyr::rename(GeneFamily_raw = `# Gene Family`) %>%
  tidyr::separate(
    GeneFamily_raw,
    into = c("KO", "Taxon"),
    sep = "\\|",
    fill = "right",
    extra = "merge",
    remove = FALSE
  ) %>%
  dplyr::mutate(
    Taxon = dplyr::na_if(Taxon, ""),
    Feature_level = dplyr::case_when(
      is.na(Taxon) ~ "KO_total",
      TRUE ~ "KO_stratified"
    ),
    Genus_raw = stringr::str_extract(Taxon, "(?<=g__)[^\\.]+"),
    Species_raw = stringr::str_extract(Taxon, "(?<=s__)[^\\.]+"),
    Genus = dplyr::case_when(
      Feature_level == "KO_total" ~ "TOTAL",
      is.na(Genus_raw) ~ "UNCLASSIFIED_GENUS",
      TRUE ~ Genus_raw
    ),
    Species = dplyr::case_when(
      Feature_level == "KO_total" ~ "TOTAL",
      is.na(Species_raw) ~ "UNCLASSIFIED_SPECIES",
      TRUE ~ Species_raw
    ),
    Taxon_status = dplyr::case_when(
      Feature_level == "KO_total" ~ "unstratified_total",
      is.na(Genus_raw) & is.na(Species_raw) ~ "taxonomically_unclassified",
      is.na(Species_raw) ~ "species_unclassified",
      TRUE ~ "species_classified"
    )
  ) %>%
  dplyr::relocate(
    KO, Feature_level, Taxon_status, Genus, Species, Taxon,
    .after = GeneFamily_raw
  ) %>%
  dplyr::select(-Genus_raw, -Species_raw)


# ---------- 1-3. CAZyme import ----------

cazy_class <- read_tsv("260224 final Input file/CAZyme/CAZy_class_TPM_matrix.tsv") %>%
  dplyr::rename(Feature = CAZy_class)

cazy_family <- read_tsv("260224 final Input file/CAZyme/CAZy_family_TPM_matrix.tsv") %>%
  dplyr::rename(Feature = CAZy_family)

cazy_subfamily <- read_tsv("260224 final Input file/CAZyme/CAZy_subfamily_TPM_matrix.tsv") %>%
  dplyr::rename(Feature = CAZy_subfamily_or_id)

cazy_module <- read_tsv("260224 final Input file/CAZyme/CAZy_module_scores.tsv")

gene_to_cazy <- read_tsv("260224 final Input file/CAZyme/gene_to_cazy_overview.tsv")

dim(cazy_class)
dim(cazy_family)
dim(cazy_subfamily)
dim(cazy_module)
dim(gene_to_cazy)

setdiff(m_before$SampleID, colnames(cazy_family))
setdiff(m_before$SampleID, cazy_module$sample_id)





# ==================== 2. CAZyme module score comparison ====================

# CAZyme module names used for plotting and tables.
# The original column name "sucrose_starch_score" is retained in the input file,
# but it is displayed as "Starch/α-glucan score" because the family/subfamily-level
# signal is mainly supported by GH13, GH77, CBM34, and CBM48, which are more
# consistent with starch- or α-glucan-associated carbohydrate processing than
# sucrose-specific degradation.

cazy_module_features <- c(
  "sucrose_starch_score",
  "fiber_degradation_score",
  "mucin_glycan_score",
  "rhamnose_pectin_score"
)

cazy_module_before <- cazy_module %>%
  dplyr::rename(SampleID = sample_id) %>%
  dplyr::inner_join(
    m_before %>%
      dplyr::select(
        SampleID, SNU_ID, SubjectID, Sex, Age, BMI,
        TNT, TRG, TRG_score, TRG_1, Response,
        ASA, Pre_Op_Tstage, Pre_Op_Nstage, AJCCstage, CEA
      ),
    by = "SampleID"
  ) %>%
  dplyr::arrange(Response, SampleID)

head(cazy_module_before)


# ---------- 2-1. Module score statistics ----------
# Wilcoxon rank-sum test is used for group comparison.
# Effect size for easy manual annotation is reported as:
#   median_diff = median(pCR) - median(non_pCR)

# ---------- Module score statistics ----------

cazy_module_result <- cazy_module_before %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(cazy_module_features),
    names_to = "Feature",
    values_to = "Score"
  ) %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(
    n_pCR = sum(Response == "pCR" & !is.na(Score)),
    n_non_pCR = sum(Response == "non_pCR" & !is.na(Score)),
    mean_pCR = mean(Score[Response == "pCR"], na.rm = TRUE),
    mean_non_pCR = mean(Score[Response == "non_pCR"], na.rm = TRUE),
    median_pCR = median(Score[Response == "pCR"], na.rm = TRUE),
    median_non_pCR = median(Score[Response == "non_pCR"], na.rm = TRUE),
    mean_diff_pCR_minus_non_pCR = mean_pCR - mean_non_pCR,
    median_diff_pCR_minus_non_pCR = median_pCR - median_non_pCR,
    p_value = wilcox.test(Score ~ Response, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    q_value = p.adjust(p_value, method = "BH"),
    Feature_label = dplyr::case_when(
      Feature == "sucrose_starch_score" ~ "Starch/α-glucan score",
      Feature == "fiber_degradation_score" ~ "Fiber degradation score",
      Feature == "mucin_glycan_score" ~ "Mucin glycan score",
      Feature == "rhamnose_pectin_score" ~ "Rhamnose/pectin score",
      TRUE ~ Feature
    )
  ) %>%
  dplyr::arrange(p_value); cazy_module_result

# write_tsv(cazy_module_result, 
#           "260224 final Input file/CAZyme/CR_nonCR_CAZy_module_result.tsv")

cazy_module_result %>%
  dplyr::select(
    Feature,
    median_pCR,
    median_non_pCR,
    median_diff_pCR_minus_non_pCR,
    p_value,
    q_value
  ) %>% 
  as.data.frame()


# ---------- 2-2. Module score plot ----------
# Violin shows the overall distribution.
# Boxplot is overlaid to show median and IQR.
# coef = 0 is used so that whiskers are effectively removed.
# Jittered points show individual samples.

p_cazy_module <- cazy_module_before %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(cazy_module_features),
    names_to = "Feature",
    values_to = "Score"
  ) %>%
  dplyr::mutate(
    Feature = factor(
      Feature,
      levels = c(
        "sucrose_starch_score",
        "fiber_degradation_score",
        "mucin_glycan_score",
        "rhamnose_pectin_score"
      ),
      labels = c(
        "Starch/α-glucan score",
        "Fiber degradation score",
        "Mucin glycan score",
        "Rhamnose/pectin score"
      )
    )
  ) %>%
  ggplot(aes(x = Response, y = Score, fill = Response, color = Response)) +
  geom_violin(
    width = 0.9,
    alpha = 0.30,
    trim = FALSE,
    linewidth = 0
  ) +
  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    coef = 0,
    fill = NA,
    linewidth = 0.55
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.3,
    alpha = 0.9
  ) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.y.npc = 0.98
  ) +
  facet_wrap(~ Feature, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_x_discrete(labels = group_labels) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(y = "CAZyme module score")

p_cazy_module

# ggsave(
#   "figures/CAZy_module_score_pCR_non_pCR_violin.svg",
#   p_cazy_module,
#   width = 8,
#   height = 5,
#   device = "svg"
# )

# ---------- CAZyme module: Starch/α-glucan score only ----------

cazy_module <- read_tsv("260224 final Input file/CAZyme/CAZy_module_scores.tsv") %>%
  dplyr::rename(SampleID = 1)

m_before <- m %>%
  dplyr::filter(
    TNT == "Before",
    TRG_1 %in% c("CR", "nonCR")
  ) %>%
  dplyr::mutate(
    Response = dplyr::case_when(
      TRG_1 == "CR" ~ "pCR",
      TRG_1 == "nonCR" ~ "non_pCR"
    ),
    Response = factor(Response, levels = c("pCR", "non_pCR"))
  ) %>%
  dplyr::arrange(Response, SampleID)

group_colors <- c("pCR" = "#4FAE9A", "non_pCR" = "#DE7872")
group_labels <- c("pCR" = "pCR", "non_pCR" = "non-pCR")

starch_module <- cazy_module %>%
  dplyr::select(SampleID, sucrose_starch_score) %>%
  dplyr::inner_join(
    m_before %>%
      dplyr::select(SampleID, Response, TRG_1, TRG_score),
    by = "SampleID"
  ) %>%
  dplyr::mutate(
    Module = "Starch/α-glucan score",
    Response = factor(Response, levels = c("pCR", "non_pCR"))
  )

# ---------- Box statistics without whiskers ----------
# We draw only the IQR box and median line manually.
# This avoids any whisker line from geom_boxplot().

starch_module_box <- starch_module %>%
  dplyr::group_by(Response) %>%
  dplyr::summarise(
    x = as.numeric(Response),
    q1 = quantile(sucrose_starch_score, 0.25, na.rm = TRUE),
    median = median(sucrose_starch_score, na.rm = TRUE),
    q3 = quantile(sucrose_starch_score, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    xmin = x - 0.13,
    xmax = x + 0.13
  )


# ---------- Publication-quality Starch/α-glucan module plot ----------
# Violin: distribution
# Box: IQR only, transparent fill, response-colored outline
# Median: response-colored horizontal line
# Jitter: individual samples
# No whiskers are drawn.

p_starch_module <- ggplot(
  starch_module,
  aes(x = Response, y = sucrose_starch_score, fill = Response, color = Response)
) +
  geom_violin(
    width = 0.92,
    trim = FALSE,
    alpha = 0.28,
    linewidth = 0.45
  ) +
  geom_rect(
    data = starch_module_box,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = q1,
      ymax = q3,
      color = Response
    ),
    inherit.aes = FALSE,
    fill = NA,
    linewidth = 0.75
  ) +
  geom_segment(
    data = starch_module_box,
    aes(
      x = xmin,
      xend = xmax,
      y = median,
      yend = median,
      color = Response
    ),
    inherit.aes = FALSE,
    linewidth = 0.9
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.6,
    alpha = 0.9
  ) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.y.npc = 0.94,
    size = 5
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_x_discrete(labels = group_labels) +
  labs(
    x = NULL,
    y = "CAZyme module score",
    title = "Starch/α-glucan score"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 13, hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

p_starch_module

# ggsave(
#   "figures/CAZy_module_starch_alpha_glucan_only.svg",
#   p_starch_module,
#   width = 2.5,
#   height = 4.3,
#   device = "svg"
# )


# ==================== 2-3. Targeted carbohydrate-processing capacity ====================

# The starch/alpha-glucan score excludes GH32 and summarizes families associated
# with starch, glycogen, maltodextrin, and alpha-glucan processing.
# The fructan/inulin score is kept separate because GH32 and GH68 include
# fructan-active enzymes, whereas CBM38 is a fructan-binding module.
# For each score, log10(TPM + 1) abundance is z-scored within each CAZy family,
# and the sample-level score is the mean z-score across the specified families.

dir.create("figures", showWarnings = FALSE)

cazy_capacity_family_sets <- tibble::tribble(
  ~Capacity,          ~Feature,
  "Starch/α-glucan",  "GH13",
  "Starch/α-glucan",  "GH65",
  "Starch/α-glucan",  "GH77",
  "Starch/α-glucan",  "GH97",
  "Starch/α-glucan",  "CBM20",
  "Starch/α-glucan",  "CBM25",
  "Starch/α-glucan",  "CBM26",
  "Starch/α-glucan",  "CBM34",
  "Starch/α-glucan",  "CBM48",
  "Fructan/inulin",   "GH32",
  "Fructan/inulin",   "GH68",
  "Fructan/inulin",   "CBM38"
)

# Internal checks for the substrate-specific score definitions.
stopifnot(
  !"GH32" %in% cazy_capacity_family_sets$Feature[
    cazy_capacity_family_sets$Capacity == "Starch/α-glucan"
  ],
  all(
    c("GH32", "GH68", "CBM38") %in%
      cazy_capacity_family_sets$Feature[
        cazy_capacity_family_sets$Capacity == "Fructan/inulin"
      ]
  )
)

if (!all(unique(cazy_capacity_family_sets$Feature) %in% cazy_family$Feature)) {
  warning(
    "CAZy families absent from the input matrix and omitted from score calculation: ",
    paste(
      setdiff(unique(cazy_capacity_family_sets$Feature), cazy_family$Feature),
      collapse = ", "
    )
  )
}

cazy_capacity_family_long <- cazy_family %>%
  dplyr::filter(Feature %in% unique(cazy_capacity_family_sets$Feature)) %>%
  dplyr::select(Feature, dplyr::all_of(m_before$SampleID)) %>%
  tidyr::pivot_longer(
    cols = -Feature,
    names_to = "SampleID",
    values_to = "TPM"
  ) %>%
  dplyr::inner_join(cazy_capacity_family_sets, by = "Feature") %>%
  dplyr::mutate(log_TPM = log10(TPM + 1)) %>%
  dplyr::group_by(Capacity, Feature) %>%
  dplyr::mutate(
    z_log_TPM = dplyr::case_when(
      is.na(sd(log_TPM, na.rm = TRUE)) | sd(log_TPM, na.rm = TRUE) == 0 ~ 0,
      TRUE ~ (log_TPM - mean(log_TPM, na.rm = TRUE)) /
        sd(log_TPM, na.rm = TRUE)
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(
    m_before %>%
      dplyr::select(SampleID, Response, TRG_1, TRG_score),
    by = "SampleID"
  ) %>%
  dplyr::mutate(
    Capacity = factor(
      Capacity,
      levels = c("Starch/α-glucan", "Fructan/inulin")
    ),
    Response = factor(Response, levels = c("pCR", "non_pCR"))
  )

cazy_capacity_scores <- cazy_capacity_family_long %>%
  dplyr::group_by(Capacity, SampleID, Response, TRG_1, TRG_score) %>%
  dplyr::summarise(
    CAZyme_capacity_score = mean(z_log_TPM, na.rm = TRUE),
    n_detected_families = sum(TPM > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Capacity, Response, SampleID)

cazy_capacity_results <- cazy_capacity_scores %>%
  dplyr::group_by(Capacity) %>%
  dplyr::summarise(
    n_pCR = sum(Response == "pCR"),
    n_non_pCR = sum(Response == "non_pCR"),
    median_pCR = median(
      CAZyme_capacity_score[Response == "pCR"],
      na.rm = TRUE
    ),
    median_non_pCR = median(
      CAZyme_capacity_score[Response == "non_pCR"],
      na.rm = TRUE
    ),
    median_diff_pCR_minus_non_pCR = median_pCR - median_non_pCR,
    p_value = wilcox.test(
      CAZyme_capacity_score ~ Response,
      exact = FALSE
    )$p.value,
    .groups = "drop"
  )

cazy_capacity_results


# ---------- Starch/alpha-glucan capacity excluding GH32 ----------

p_cazy_starch_alpha_glucan <- cazy_capacity_scores %>%
  dplyr::filter(Capacity == "Starch/α-glucan") %>%
  ggplot(
    aes(
      x = Response,
      y = CAZyme_capacity_score,
      fill = Response,
      color = Response
    )
  ) +
  geom_violin(
    width = 0.92,
    trim = FALSE,
    alpha = 0.28,
    linewidth = 0.45
  ) +
  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    coef = 0,
    fill = NA,
    linewidth = 0.65
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.6,
    alpha = 0.9
  ) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.y.npc = 0.94,
    size = 5
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_x_discrete(labels = group_labels) +
  labs(
    x = NULL,
    y = "Standardized CAZyme capacity score",
    title = "Starch/α-glucan-associated capacity"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 13, hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

p_cazy_starch_alpha_glucan

ggsave(
  "figures/CAZy_starch_alpha_glucan_capacity_excluding_GH32.svg",
  p_cazy_starch_alpha_glucan,
  width = 3.7,
  height = 4.6,
  device = "svg"
)
Sys.setFileTime(
  "figures/CAZy_starch_alpha_glucan_capacity_excluding_GH32.svg",
  Sys.time()
)


# ---------- GH32 abundance ----------

cazy_gh32 <- cazy_capacity_family_long %>%
  dplyr::filter(Feature == "GH32") %>%
  dplyr::select(
    SampleID, Response, TRG_1, TRG_score,
    GH32_TPM = TPM,
    log_GH32_TPM = log_TPM
  ) %>%
  dplyr::arrange(Response, SampleID)

cazy_gh32_result <- cazy_gh32 %>%
  dplyr::summarise(
    n_pCR = sum(Response == "pCR"),
    n_non_pCR = sum(Response == "non_pCR"),
    median_TPM_pCR = median(GH32_TPM[Response == "pCR"], na.rm = TRUE),
    median_TPM_non_pCR = median(
      GH32_TPM[Response == "non_pCR"],
      na.rm = TRUE
    ),
    median_diff_logTPM_pCR_minus_non_pCR =
      median(log_GH32_TPM[Response == "pCR"], na.rm = TRUE) -
      median(log_GH32_TPM[Response == "non_pCR"], na.rm = TRUE),
    p_value = wilcox.test(
      log_GH32_TPM ~ Response,
      exact = FALSE
    )$p.value
  )

cazy_gh32_result

p_cazy_gh32 <- ggplot(
  cazy_gh32,
  aes(
    x = Response,
    y = log_GH32_TPM,
    fill = Response,
    color = Response
  )
) +
  geom_violin(
    width = 0.92,
    trim = FALSE,
    alpha = 0.28,
    linewidth = 0.45
  ) +
  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    coef = 0,
    fill = NA,
    linewidth = 0.65
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.6,
    alpha = 0.9
  ) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.y.npc = 0.94,
    size = 5
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_x_discrete(labels = group_labels) +
  labs(
    x = NULL,
    y = "log10(TPM + 1)",
    title = "GH32 abundance"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 13, hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

p_cazy_gh32

ggsave(
  "figures/CAZy_GH32_abundance_pCR_non_pCR.svg",
  p_cazy_gh32,
  width = 3.7,
  height = 4.6,
  device = "svg"
)
Sys.setFileTime(
  "figures/CAZy_GH32_abundance_pCR_non_pCR.svg",
  Sys.time()
)


# ---------- Fructan/inulin-associated capacity ----------

p_cazy_fructan_inulin <- cazy_capacity_scores %>%
  dplyr::filter(Capacity == "Fructan/inulin") %>%
  ggplot(
    aes(
      x = Response,
      y = CAZyme_capacity_score,
      fill = Response,
      color = Response
    )
  ) +
  geom_violin(
    width = 0.92,
    trim = FALSE,
    alpha = 0.28,
    linewidth = 0.45
  ) +
  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    coef = 0,
    fill = NA,
    linewidth = 0.65
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.6,
    alpha = 0.9
  ) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.y.npc = 0.94,
    size = 5
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_x_discrete(labels = group_labels) +
  labs(
    x = NULL,
    y = "Standardized CAZyme capacity score",
    title = "Fructan/inulin-associated capacity",
    subtitle = "GH32, GH68, and CBM38"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

p_cazy_fructan_inulin

ggsave(
  "figures/CAZy_fructan_inulin_capacity_pCR_non_pCR.svg",
  p_cazy_fructan_inulin,
  width = 3.7,
  height = 4.6,
  device = "svg"
)
Sys.setFileTime(
  "figures/CAZy_fructan_inulin_capacity_pCR_non_pCR.svg",
  Sys.time()
)


# ==================== 3. CAZy family-level comparison ====================

cazy_family_before <- cazy_family %>%
  dplyr::select(Feature, dplyr::all_of(m_before$SampleID)) %>%
  tidyr::pivot_longer(
    cols = -Feature,
    names_to = "SampleID",
    values_to = "TPM"
  ) %>%
  dplyr::inner_join(
    m_before %>%
      dplyr::select(
        SampleID, SNU_ID, SubjectID, Sex, Age, BMI,
        TNT, TRG, TRG_score, TRG_1, Response,
        ASA, Pre_Op_Tstage, Pre_Op_Nstage, AJCCstage, CEA
      ),
    by = "SampleID"
  ) %>%
  dplyr::mutate(
    log_TPM = log10(TPM + 1)
  ) %>%
  dplyr::arrange(Response, SampleID, Feature)

head(cazy_family_before)

# ---------- 3-1. CAZy family statistics ----------
# Family-level tests are run on log10(TPM + 1).
# Sparse families (prevalence < 0.2) are retained in the table but not tested.
# Median difference is also reported for manual annotation if needed.

cazy_family_result <- cazy_family_before %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(
    n_pCR = sum(Response == "pCR"),
    n_non_pCR = sum(Response == "non_pCR"),
    prevalence = mean(TPM > 0),
    mean_TPM_pCR = mean(TPM[Response == "pCR"], na.rm = TRUE),
    mean_TPM_non_pCR = mean(TPM[Response == "non_pCR"], na.rm = TRUE),
    median_TPM_pCR = median(TPM[Response == "pCR"], na.rm = TRUE),
    median_TPM_non_pCR = median(TPM[Response == "non_pCR"], na.rm = TRUE),
    log2FC_pCR_vs_non_pCR = log2((mean_TPM_pCR + 1) / (mean_TPM_non_pCR + 1)),
    median_diff_logTPM_pCR_minus_non_pCR =
      median(log_TPM[Response == "pCR"], na.rm = TRUE) -
      median(log_TPM[Response == "non_pCR"], na.rm = TRUE),
    p_value = dplyr::case_when(
      prevalence >= 0.2 ~ wilcox.test(log_TPM ~ Response, exact = FALSE)$p.value,
      TRUE ~ NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    q_value = p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::arrange(p_value)

head(cazy_family_result, 30)

# write_tsv(cazy_family_result, "260224 final Input file/CAZyme/CR_nonCR_CAZy_family_result.tsv")


# ---------- 3-2. Targeted carbohydrate-related CAZy families ----------

target_cazy_families <- c(
  "GH13", "GH32", "GH65", "GH77", "GH97",
  "CBM20", "CBM25", "CBM26", "CBM48",
  "GH28", "GH78", "GH105", "GH106",
  "PL1", "PL9", "PL11",
  "CE8", "CE12",
  "GH2", "GH20", "GH29", "GH33", "GH35",
  "GH84", "GH89", "GH95", "GH101",
  "GH109", "GH110", "GH112", "GH123",
  "CE4",
  "GH3", "GH5", "GH10", "GH11", "GH16",
  "GH26", "GH30", "GH43", "GH51", "GH53",
  "GH67", "GH74",
  "PL7", "PL15",
  "CBM6", "CBM13", "CBM32"
)

cazy_family_target_result <- cazy_family_result %>%
  dplyr::filter(Feature %in% target_cazy_families) %>%
  dplyr::arrange(p_value)

cazy_family_target_result

cazy_family_target_result %>%
  dplyr::select(
    Feature,
    median_TPM_pCR,
    median_TPM_non_pCR,
    median_diff_logTPM_pCR_minus_non_pCR,
    p_value,
    q_value
  ) %>% 
  as.data.frame() %>% 
  filter(p_value <0.1)

# write_tsv(cazy_family_target_result, "260224 final Input file/CAZyme/CR_nonCR_CAZy_family_target_result.tsv")


# ---------- 3-3. Targeted CAZy family plot ----------

p_cazy_target_family <- cazy_family_before %>%
  dplyr::filter(Feature %in% target_cazy_families) %>%
  dplyr::mutate(
    Feature = factor(
      Feature,
      levels = cazy_family_target_result$Feature
    )
  ) %>%
  ggplot(aes(x = Response, y = log_TPM, fill = Response, color = Response)) +
  geom_violin(
    width = 0.9,
    alpha = 0.35,
    trim = FALSE,
    linewidth = 0
  ) +
  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    coef = 0,
    fill = "transparent",
    linewidth = 0.4
  ) +
  geom_jitter(
    width = 0.08,
    size = 1.7,
    alpha = 0.9
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  facet_wrap(~ Feature, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_x_discrete(labels = c("pCR" = "pCR", "non_pCR" = "non-pCR")) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(y = "log10(TPM + 1)")

p_cazy_target_family

# ggsave(
#   "figures/CAZy_target_family_pCR_non_pCR_violin.svg",
#   p_cazy_target_family,
#   width = 12,
#   height = 9
# )


# ==================== 4. CAZy class-level comparison ====================

cazy_class_before <- cazy_class %>%
  dplyr::select(Feature, dplyr::all_of(m_before$SampleID)) %>%
  tidyr::pivot_longer(
    cols = -Feature,
    names_to = "SampleID",
    values_to = "TPM"
  ) %>%
  dplyr::inner_join(
    m_before %>%
      dplyr::select(SampleID, Response, TRG_1, TNT, TRG, TRG_score),
    by = "SampleID"
  ) %>%
  dplyr::mutate(
    log_TPM = log10(TPM + 1)
  ) %>%
  dplyr::arrange(Response, SampleID, Feature)

# ---------- 4-1. CAZy class statistics ----------

cazy_class_result <- cazy_class_before %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(
    n_pCR = sum(Response == "pCR"),
    n_non_pCR = sum(Response == "non_pCR"),
    mean_TPM_pCR = mean(TPM[Response == "pCR"], na.rm = TRUE),
    mean_TPM_non_pCR = mean(TPM[Response == "non_pCR"], na.rm = TRUE),
    median_logTPM_pCR = median(log_TPM[Response == "pCR"], na.rm = TRUE),
    median_logTPM_non_pCR = median(log_TPM[Response == "non_pCR"], na.rm = TRUE),
    median_diff_logTPM_pCR_minus_non_pCR = median_logTPM_pCR - median_logTPM_non_pCR,
    log2FC_pCR_vs_non_pCR = log2((mean_TPM_pCR + 1) / (mean_TPM_non_pCR + 1)),
    p_value = wilcox.test(log_TPM ~ Response, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    q_value = p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::arrange(p_value)

cazy_class_result

cazy_class_result %>%
  dplyr::select(
    Feature,
    median_logTPM_pCR,
    median_logTPM_non_pCR,
    median_diff_logTPM_pCR_minus_non_pCR,
    p_value,
    q_value
  ) %>% 
  as.data.frame() %>% 
  filter(p_value < 0.1)

# write_tsv(cazy_class_result, "260224 final Input file/CAZyme/CR_nonCR_CAZy_class_result.tsv")



# ---------- 4-2. CAZy class plot ----------

p_cazy_class <- cazy_class_before %>%
  ggplot(aes(x = Response, y = log_TPM, fill = Response, color = Response)) +
  geom_violin(
    width = 0.9,
    alpha = 0.35,
    trim = FALSE,
    linewidth = 0
  ) +
  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    coef = 0,
    fill = "transparent",
    linewidth = 0.45
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.2,
    alpha = 0.9
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  facet_wrap(~ Feature, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_x_discrete(labels = c("pCR" = "pCR", "non_pCR" = "non-pCR")) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(y = "log10(TPM + 1)")

p_cazy_class

# ggsave(
#   "figures/CAZy_class_pCR_non_pCR_violin.svg",
#   p_cazy_class,
#   width = 7,
#   height = 5
# )


# ==================== AA class abundance plot ====================

# AA is a CAZyme class-level feature, not a module score.
# Therefore, its y-axis is log10(TPM + 1), not module score.

group_colors <- c("pCR" = "#4FAE9A", "non_pCR" = "#DE7872")
group_labels <- c("pCR" = "pCR", "non_pCR" = "non-pCR")

dir.create("figures", showWarnings = FALSE)

# ---------- 1. Before sample metadata ----------

m_before <- m %>%
  dplyr::filter(
    TNT == "Before",
    TRG_1 %in% c("CR", "nonCR")
  ) %>%
  dplyr::mutate(
    Response = dplyr::case_when(
      TRG_1 == "CR" ~ "pCR",
      TRG_1 == "nonCR" ~ "non_pCR"
    ),
    Response = factor(Response, levels = c("pCR", "non_pCR"))
  ) %>%
  dplyr::arrange(Response, SampleID)


# ---------- 2. AA class abundance table ----------

aa_class <- cazy_class %>%
  dplyr::filter(Feature == "AA") %>%
  dplyr::select(Feature, dplyr::all_of(m_before$SampleID)) %>%
  tidyr::pivot_longer(
    cols = -Feature,
    names_to = "SampleID",
    values_to = "TPM"
  ) %>%
  dplyr::inner_join(
    m_before %>%
      dplyr::select(SampleID, Response, TRG_1, TRG_score),
    by = "SampleID"
  ) %>%
  dplyr::mutate(
    log_TPM = log10(TPM + 1),
    Feature_label = "AA: Auxiliary activities"
  )


# ---------- 3. AA class statistics ----------

aa_class_result <- aa_class %>%
  dplyr::summarise(
    n_pCR = sum(Response == "pCR" & !is.na(log_TPM)),
    n_non_pCR = sum(Response == "non_pCR" & !is.na(log_TPM)),
    median_pCR = median(log_TPM[Response == "pCR"], na.rm = TRUE),
    median_non_pCR = median(log_TPM[Response == "non_pCR"], na.rm = TRUE),
    median_diff_pCR_minus_non_pCR = median_pCR - median_non_pCR,
    mean_pCR = mean(log_TPM[Response == "pCR"], na.rm = TRUE),
    mean_non_pCR = mean(log_TPM[Response == "non_pCR"], na.rm = TRUE),
    mean_diff_pCR_minus_non_pCR = mean_pCR - mean_non_pCR,
    p_value = wilcox.test(log_TPM ~ Response, exact = FALSE)$p.value
  )

aa_class_result

# write_tsv(
#   aa_class_result,
#   "260224 final Input file/CAZyme/AA_class_pCR_non_pCR_result.tsv"
# )


# ---------- 4. Manual IQR box without whiskers ----------

aa_class_box <- aa_class %>%
  dplyr::group_by(Response) %>%
  dplyr::summarise(
    x = as.numeric(Response),
    q1 = quantile(log_TPM, 0.25, na.rm = TRUE),
    median = median(log_TPM, na.rm = TRUE),
    q3 = quantile(log_TPM, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    xmin = x - 0.13,
    xmax = x + 0.13
  )


# ---------- 5. Publication-quality AA plot ----------

p_aa_class <- ggplot(
  aa_class,
  aes(x = Response, y = log_TPM, fill = Response, color = Response)
) +
  geom_violin(
    width = 0.92,
    trim = FALSE,
    alpha = 0.28,
    linewidth = 0.45
  ) +
  geom_rect(
    data = aa_class_box,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = q1,
      ymax = q3,
      color = Response
    ),
    inherit.aes = FALSE,
    fill = NA,
    linewidth = 0.75
  ) +
  geom_segment(
    data = aa_class_box,
    aes(
      x = xmin,
      xend = xmax,
      y = median,
      yend = median,
      color = Response
    ),
    inherit.aes = FALSE,
    linewidth = 0.9
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.6,
    alpha = 0.9
  ) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.y.npc = 0.94,
    size = 5
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_x_discrete(labels = group_labels) +
  labs(
    x = NULL,
    y = "CAZyme class abundance, log10(TPM + 1)",
    title = "AA: Auxiliary activities"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 13, hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

p_aa_class

# ggsave(
#   "figures/CAZy_AA_class_pCR_non_pCR.svg",
#   p_aa_class,
#   width = 2.5,
#   height = 4.3,
#   device = "svg"
# )





# ==================== 5. CAZy subfamily-level comparison ====================

# Subfamily-level features are high-dimensional.
# Use this as exploratory analysis with stricter prevalence filtering.

cazy_subfamily_before <- cazy_subfamily %>%
  dplyr::select(Feature, dplyr::all_of(m_before$SampleID)) %>%
  tidyr::pivot_longer(
    cols = -Feature,
    names_to = "SampleID",
    values_to = "TPM"
  ) %>%
  dplyr::inner_join(
    m_before %>%
      dplyr::select(SampleID, TRG_1, TNT, TRG, TRG_score),
    by = "SampleID"
  ) %>%
  dplyr::mutate(
    log_TPM = log10(TPM + 1)
  )

cazy_subfamily_result <- cazy_subfamily_before %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(
    n_nonCR = sum(TRG_1 == "nonCR"),
    n_CR = sum(TRG_1 == "CR"),
    prevalence = mean(TPM > 0),
    mean_TPM_nonCR = mean(TPM[TRG_1 == "nonCR"], na.rm = TRUE),
    mean_TPM_CR = mean(TPM[TRG_1 == "CR"], na.rm = TRUE),
    log2FC_CR_vs_nonCR = log2((mean_TPM_CR + 1) / (mean_TPM_nonCR + 1)),
    p_value = dplyr::case_when(
      prevalence >= 0.3 ~ wilcox.test(log_TPM ~ TRG_1, exact = FALSE)$p.value,
      TRUE ~ NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    q_value = p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::arrange(p_value)

cazy_subfamily_result %>% 
  as.data.frame() %>% 
  filter(p_value < 0.01)

  head(cazy_subfamily_result, 30)

# write_tsv(cazy_subfamily_result, 
#           "260224 final Input file/CAZyme/CR_nonCR_CAZy_subfamily_result.tsv")





# ====================  6. CAZyme summary analysis ====================

# This block replaces the previous duplicated CAZyme significant-feature sections.
# Main purpose:
#   1) reduce redundant feature display
#   2) keep interpretable Module / Class / Family features
#   3) add functional annotation for CAZyme interpretation
#   4) generate violin plot, median-difference table, and heatmap
#
# Subfamily-level results should remain separate as exploratory/supplementary analysis.

# ---------- Plotting and feature-selection options ----------


# ---------- 2. Ensure log2FC columns ----------

# Positive log2FC means higher in pCR.
# Negative log2FC means higher in non-pCR.

if (!"log2FC_pCR_vs_non_pCR" %in% colnames(cazy_class_result) &&
    "log2FC_CR_vs_nonCR" %in% colnames(cazy_class_result)) {
  
  cazy_class_result <- cazy_class_result %>%
    dplyr::mutate(
      log2FC_pCR_vs_non_pCR = log2FC_CR_vs_nonCR
    )
}

if (!"log2FC_pCR_vs_non_pCR" %in% colnames(cazy_class_result) &&
    all(c("mean_TPM_pCR", "mean_TPM_non_pCR") %in% colnames(cazy_class_result))) {
  
  cazy_class_result <- cazy_class_result %>%
    dplyr::mutate(
      log2FC_pCR_vs_non_pCR = log2((mean_TPM_pCR + 1) / (mean_TPM_non_pCR + 1))
    )
}

if (!"log2FC_pCR_vs_non_pCR" %in% colnames(cazy_class_result) &&
    all(c("mean_TPM_CR", "mean_TPM_nonCR") %in% colnames(cazy_class_result))) {
  
  cazy_class_result <- cazy_class_result %>%
    dplyr::mutate(
      log2FC_pCR_vs_non_pCR = log2((mean_TPM_CR + 1) / (mean_TPM_nonCR + 1))
    )
}

if (!"log2FC_pCR_vs_non_pCR" %in% colnames(cazy_family_result) &&
    "log2FC_CR_vs_nonCR" %in% colnames(cazy_family_result)) {
  
  cazy_family_result <- cazy_family_result %>%
    dplyr::mutate(
      log2FC_pCR_vs_non_pCR = log2FC_CR_vs_nonCR
    )
}

if (!"log2FC_pCR_vs_non_pCR" %in% colnames(cazy_family_result) &&
    all(c("mean_TPM_pCR", "mean_TPM_non_pCR") %in% colnames(cazy_family_result))) {
  
  cazy_family_result <- cazy_family_result %>%
    dplyr::mutate(
      log2FC_pCR_vs_non_pCR = log2((mean_TPM_pCR + 1) / (mean_TPM_non_pCR + 1))
    )
}

if (!"log2FC_pCR_vs_non_pCR" %in% colnames(cazy_family_result) &&
    all(c("mean_TPM_CR", "mean_TPM_nonCR") %in% colnames(cazy_family_result))) {
  
  cazy_family_result <- cazy_family_result %>%
    dplyr::mutate(
      log2FC_pCR_vs_non_pCR = log2((mean_TPM_CR + 1) / (mean_TPM_nonCR + 1))
    )
}

stopifnot("log2FC_pCR_vs_non_pCR" %in% colnames(cazy_class_result))
stopifnot("log2FC_pCR_vs_non_pCR" %in% colnames(cazy_family_result))


# ---------- 3. Module labels ----------

cazy_module_features <- c(
  "sucrose_starch_score",
  "fiber_degradation_score",
  "mucin_glycan_score",
  "rhamnose_pectin_score"
)

cazy_module_result <- cazy_module_result %>%
  dplyr::mutate(
    Feature_label = dplyr::case_when(
      Feature == "sucrose_starch_score" ~ "Starch/α-glucan score",
      Feature == "fiber_degradation_score" ~ "Fiber degradation score",
      Feature == "mucin_glycan_score" ~ "Mucin glycan score",
      Feature == "rhamnose_pectin_score" ~ "Rhamnose/pectin score",
      TRUE ~ Feature
    ),
    Functional_group = dplyr::case_when(
      Feature == "sucrose_starch_score" ~ "Starch/α-glucan utilization potential",
      Feature == "fiber_degradation_score" ~ "Plant fiber degradation potential",
      Feature == "mucin_glycan_score" ~ "Mucin/host glycan utilization potential",
      Feature == "rhamnose_pectin_score" ~ "Rhamnose/pectin degradation potential",
      TRUE ~ "Other CAZyme module"
    ),
    Functional_annotation = dplyr::case_when(
      Feature == "sucrose_starch_score" ~ "Module-level summary of CAZy families related to starch, maltodextrin, glycogen, and α-glucan processing; interpreted as starch/α-glucan potential rather than sucrose-specific degradation.",
      Feature == "fiber_degradation_score" ~ "Module-level summary of CAZy families related to complex plant polysaccharide degradation.",
      Feature == "mucin_glycan_score" ~ "Module-level summary of CAZy families related to host glycan and mucin-associated carbohydrate utilization.",
      Feature == "rhamnose_pectin_score" ~ "Module-level summary of CAZy families related to pectin, rhamnogalacturonan, and rhamnose-containing glycan degradation.",
      TRUE ~ "Other CAZyme module."
    )
  )


# ---------- 4. Class-level functional annotation ----------

cazy_class_result <- cazy_class_result %>%
  dplyr::mutate(
    Feature_label = dplyr::case_when(
      Feature == "AA" ~ "AA: Auxiliary activities",
      Feature == "GH" ~ "GH: Glycoside hydrolases",
      Feature == "GT" ~ "GT: Glycosyltransferases",
      Feature == "CE" ~ "CE: Carbohydrate esterases",
      Feature == "CBM" ~ "CBM: Carbohydrate-binding modules",
      Feature == "PL" ~ "PL: Polysaccharide lyases",
      TRUE ~ Feature
    ),
    Functional_group = dplyr::case_when(
      Feature == "AA" ~ "Auxiliary redox carbohydrate processing",
      Feature == "GH" ~ "Glycosidic bond hydrolysis",
      Feature == "GT" ~ "Glycan biosynthesis/modification",
      Feature == "CE" ~ "Carbohydrate de-esterification",
      Feature == "CBM" ~ "Carbohydrate binding/substrate targeting",
      Feature == "PL" ~ "Polysaccharide lyase-mediated degradation",
      TRUE ~ "Other CAZyme class"
    ),
    Functional_annotation = dplyr::case_when(
      Feature == "AA" ~ "Auxiliary activity class; redox-associated enzymes that can support oxidative processing of complex carbohydrates or biomass-associated compounds.",
      Feature == "GH" ~ "Glycoside hydrolase class; enzymes involved in hydrolysis or rearrangement of glycosidic bonds.",
      Feature == "GT" ~ "Glycosyltransferase class; enzymes involved in glycan biosynthesis or transfer of sugar moieties.",
      Feature == "CE" ~ "Carbohydrate esterase class; enzymes that remove ester-linked modifications from carbohydrates.",
      Feature == "CBM" ~ "Non-catalytic carbohydrate-binding modules that help target enzymes to carbohydrate substrates.",
      Feature == "PL" ~ "Polysaccharide lyase class; enzymes involved in non-hydrolytic cleavage of polysaccharides.",
      TRUE ~ "Other CAZyme class."
    )
  )


# ---------- 5. Family-level functional annotation ----------

cazy_family_result <- cazy_family_result %>%
  dplyr::mutate(
    Feature_label = dplyr::case_when(
      Feature == "GH13" ~ "GH13: α-amylase family",
      Feature == "GH77" ~ "GH77: α-glucanotransferase-like",
      Feature == "GH5" ~ "GH5: broad β-polysaccharide hydrolases",
      Feature == "CBM13" ~ "CBM13: carbohydrate-binding module",
      Feature == "CBM48" ~ "CBM48: α-glucan/glycogen-binding module",
      Feature == "CBM34" ~ "CBM34: starch-binding module",
      Feature == "GH53" ~ "GH53: β-1,4-galactanase-like",
      Feature == "GH43" ~ "GH43: xylan/arabinose-glycan related",
      Feature == "AA4" ~ "AA4: vanillyl-alcohol oxidase-like",
      Feature == "GT111" ~ "GT111: β-1,3-galactofuranosyltransferase",
      TRUE ~ Feature
    ),
    Functional_group = dplyr::case_when(
      Feature %in% c("GH13", "GH77", "CBM48", "CBM34") ~ "Starch/α-glucan-associated CAZymes",
      Feature %in% c("GH5", "CBM13", "GH53", "GH43") ~ "Plant polysaccharide/fiber-associated CAZymes",
      Feature == "AA4" ~ "Auxiliary redox/phenolic oxidation",
      Feature == "GT111" ~ "Glycan biosynthesis/modification",
      stringr::str_detect(Feature, "^AA") ~ "Auxiliary redox carbohydrate processing",
      stringr::str_detect(Feature, "^GT") ~ "Glycan biosynthesis/modification",
      stringr::str_detect(Feature, "^CBM") ~ "Carbohydrate binding/substrate targeting",
      stringr::str_detect(Feature, "^GH") ~ "Glycoside hydrolase activity",
      stringr::str_detect(Feature, "^CE") ~ "Carbohydrate de-esterification",
      stringr::str_detect(Feature, "^PL") ~ "Polysaccharide lyase-mediated degradation",
      TRUE ~ "Other CAZyme family"
    ),
    Functional_annotation = dplyr::case_when(
      Feature == "GH13" ~ "α-amylase family; interpreted here as starch, maltodextrin, glycogen, and α-glucan processing potential.",
      Feature == "GH77" ~ "4-α-glucanotransferase/amylomaltase-like family; supports starch-derived α-glucan remodeling rather than sucrose-specific degradation.",
      Feature == "GH5" ~ "Broad β-linked polysaccharide hydrolase family; may reflect complex plant polysaccharide degradation potential.",
      Feature == "CBM13" ~ "Carbohydrate-binding module; interpreted as substrate-binding or targeting capacity, potentially linked to plant glycan/xylan/arabinose-rich substrates.",
      Feature == "CBM48" ~ "α-glucan/glycogen-binding module; supports starch/glycogen/branched α-glucan-associated interpretation.",
      Feature == "CBM34" ~ "Starch-binding module; supports starch-active multidomain CAZyme interpretation.",
      Feature == "GH53" ~ "β-1,4-galactanase-like family; linked to galactan or pectin side-chain degradation.",
      Feature == "GH43" ~ "Family including xylan/arabinoxylan/arabinofuranose-related activities; interpreted as hemicellulose or plant glycan utilization potential.",
      Feature == "AA4" ~ "Auxiliary activity family 4; vanillyl-alcohol oxidase-like oxidative processing of phenolic/aromatic compounds.",
      Feature == "GT111" ~ "Glycosyltransferase family 111; UDP-Galf β-1,3-galactofuranosyltransferase-like glycan biosynthesis/modification.",
      stringr::str_detect(Feature, "^AA") ~ "Auxiliary activity family; broad redox-associated carbohydrate or biomass-related processing.",
      stringr::str_detect(Feature, "^GT") ~ "Glycosyltransferase family; glycan biosynthesis or glycan modification rather than carbohydrate degradation.",
      stringr::str_detect(Feature, "^CBM") ~ "Carbohydrate-binding module; substrate targeting rather than direct catalytic degradation.",
      stringr::str_detect(Feature, "^GH") ~ "Glycoside hydrolase family; broad carbohydrate degradation or glycosidic bond processing.",
      stringr::str_detect(Feature, "^CE") ~ "Carbohydrate esterase family; removal of ester-linked modifications from carbohydrate substrates.",
      stringr::str_detect(Feature, "^PL") ~ "Polysaccharide lyase family; non-hydrolytic cleavage of polysaccharides.",
      TRUE ~ "Other CAZyme family."
    )
  )


# ---------- 6. Select a compact family set ----------

# Family selection logic:
#   1) top 3 pCR-enriched families with p < 0.05
#   2) top non-pCR-enriched families with p < 0.05
#   3) mechanistically linked families with p < 0.1
#
# AA4 is removed from the main figure if the broader AA class is significant,
# because AA class and AA4 can be visually redundant in a compact overview.

mechanistic_p_cutoff <- 0.1

cazy_family_top_pcr <- cazy_family_result %>%
  dplyr::filter(
    !is.na(p_value),
    p_value < sig_p_cutoff,
    log2FC_pCR_vs_non_pCR > 0
  ) %>%
  dplyr::arrange(p_value) %>%
  dplyr::slice_head(n = 3) %>%
  dplyr::mutate(
    Selection_reason = "Top pCR-enriched family"
  )

cazy_family_top_non_pcr <- cazy_family_result %>%
  dplyr::filter(
    !is.na(p_value),
    p_value < sig_p_cutoff,
    log2FC_pCR_vs_non_pCR < 0
  ) %>%
  dplyr::arrange(p_value) %>%
  dplyr::slice_head(n = 3) %>%
  dplyr::mutate(
    Selection_reason = "Top non-pCR-enriched family"
  )

cazy_family_mechanistic <- cazy_family_result %>%
  dplyr::filter(
    !is.na(p_value),
    p_value < 0.1,
    Feature %in% c("GH13", "GH77", "GH5", "CBM13", "CBM48", "CBM34", "GH53", "GH43")
  ) %>%
  dplyr::arrange(p_value) %>%
  dplyr::mutate(
    Selection_reason = "Mechanistically linked family"
  )

cazy_family_selected_result <- dplyr::bind_rows(
  cazy_family_top_pcr,
  cazy_family_top_non_pcr,
  cazy_family_mechanistic
) %>%
  dplyr::distinct(Feature, .keep_all = TRUE) %>%
  dplyr::filter(
    !(Feature == "AA4" & any(cazy_class_result$Feature == "AA" & cazy_class_result$p_value < sig_p_cutoff))
  ) %>%
  dplyr::arrange(p_value)

cazy_family_interpretation_table <- cazy_family_result %>%
  dplyr::filter(
    p_value < 0.1 |
      Feature %in% c("AA4", "GT111", "GH13", "GH77", "GH5", "CBM13", "CBM48", "CBM34")
  ) %>%
  dplyr::arrange(p_value) %>%
  dplyr::select(
    Feature,
    Feature_label,
    Functional_group,
    Functional_annotation,
    prevalence,
    log2FC_pCR_vs_non_pCR,
    p_value,
    q_value
  )

write_tsv(
  cazy_family_selected_result,
  "260224 final Input file/CAZyme/CAZy_selected_family_for_main_figure.tsv"
)

write_tsv(
  cazy_family_interpretation_table,
  "260224 final Input file/CAZyme/CAZy_family_interpretation_table_p0.1.tsv"
)


# ---------- 7. Build compact Module + Class + Family table ----------

cazy_sig_result <- dplyr::bind_rows(
  cazy_module_result %>%
    dplyr::filter(p_value < sig_p_cutoff) %>%
    dplyr::transmute(
      Feature_level = "Module",
      Feature = Feature,
      Feature_label = Feature_label,
      Functional_group = Functional_group,
      Functional_annotation = Functional_annotation,
      p_value = p_value,
      q_value = q_value,
      log2FC_pCR_vs_non_pCR = NA_real_,
      Selection_reason = "Significant module"
    ),
  
  cazy_class_result %>%
    dplyr::filter(p_value < sig_p_cutoff) %>%
    dplyr::transmute(
      Feature_level = "Class",
      Feature = Feature,
      Feature_label = Feature_label,
      Functional_group = Functional_group,
      Functional_annotation = Functional_annotation,
      p_value = p_value,
      q_value = q_value,
      log2FC_pCR_vs_non_pCR = log2FC_pCR_vs_non_pCR,
      Selection_reason = "Significant CAZyme class"
    ),
  
  cazy_family_selected_result %>%
    dplyr::transmute(
      Feature_level = "Family",
      Feature = Feature,
      Feature_label = Feature_label,
      Functional_group = Functional_group,
      Functional_annotation = Functional_annotation,
      p_value = p_value,
      q_value = q_value,
      log2FC_pCR_vs_non_pCR = log2FC_pCR_vs_non_pCR,
      Selection_reason = Selection_reason
    )
) %>%
  dplyr::mutate(
    Feature_level = factor(Feature_level, levels = c("Module", "Class", "Family")),
    Facet_label = paste0(Feature_label, "\n", Functional_group)
  ) %>%
  dplyr::arrange(Feature_level, p_value)

# write_tsv(
#   cazy_sig_result,
#   "260224 final Input file/CAZyme/CAZy_selected_module_class_family_features.tsv"
# )


# ---------- 8. Plotting table ----------

cazy_sig_plot_data <- dplyr::bind_rows(
  cazy_module_before %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(cazy_module_features),
      names_to = "Feature",
      values_to = "Value"
    ) %>%
    dplyr::transmute(
      SampleID = SampleID,
      Response = Response,
      Feature_level = "Module",
      Feature = Feature,
      Value = Value,
      Value_scale = "Module score"
    ),
  
  cazy_class_before %>%
    dplyr::transmute(
      SampleID = SampleID,
      Response = Response,
      Feature_level = "Class",
      Feature = Feature,
      Value = log_TPM,
      Value_scale = "log10(TPM + 1)"
    ),
  
  cazy_family_before %>%
    dplyr::transmute(
      SampleID = SampleID,
      Response = Response,
      Feature_level = "Family",
      Feature = Feature,
      Value = log_TPM,
      Value_scale = "log10(TPM + 1)"
    )
) %>%
  dplyr::mutate(
    Response = factor(Response, levels = c("pCR", "non_pCR")),
    Feature_level = factor(Feature_level, levels = c("Module", "Class", "Family"))
  ) %>%
  dplyr::inner_join(
    cazy_sig_result,
    by = c("Feature_level", "Feature")
  ) %>%
  dplyr::mutate(
    Facet_label = factor(Facet_label, levels = cazy_sig_result$Facet_label)
  )


# ---------- 9. Median difference summary ----------

cazy_sig_median_difference <- cazy_sig_plot_data %>%
  dplyr::group_by(
    Feature_level,
    Feature,
    Feature_label,
    Functional_group,
    Functional_annotation,
    Facet_label,
    Value_scale,
    p_value,
    q_value,
    log2FC_pCR_vs_non_pCR,
    Selection_reason
  ) %>%
  dplyr::summarise(
    n_pCR = sum(Response == "pCR" & !is.na(Value)),
    n_non_pCR = sum(Response == "non_pCR" & !is.na(Value)),
    median_pCR = median(Value[Response == "pCR"], na.rm = TRUE),
    median_non_pCR = median(Value[Response == "non_pCR"], na.rm = TRUE),
    mean_pCR = mean(Value[Response == "pCR"], na.rm = TRUE),
    mean_non_pCR = mean(Value[Response == "non_pCR"], na.rm = TRUE),
    median_diff_pCR_minus_non_pCR = median_pCR - median_non_pCR,
    mean_diff_pCR_minus_non_pCR = mean_pCR - mean_non_pCR,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Direction = dplyr::case_when(
      median_diff_pCR_minus_non_pCR > 0 ~ "Higher in pCR",
      median_diff_pCR_minus_non_pCR < 0 ~ "Higher in non-pCR",
      TRUE ~ "No median difference"
    )
  ) %>%
  dplyr::arrange(Feature_level, p_value)

# write_tsv(
#   cazy_sig_median_difference,
#   "260224 final Input file/CAZyme/CAZy_selected_module_class_family_median_difference.tsv"
# )


# ---------- 10. Compact violin plot ----------

p_cazy_sig_all <- cazy_sig_plot_data %>%
  ggplot(aes(x = Response, y = Value, fill = Response, color = Response)) +
  geom_violin(
    width = 0.9,
    alpha = 0.30,
    trim = FALSE,
    linewidth = 0
  ) +
  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    coef = 0,
    fill = NA,
    linewidth = 0.55
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.2,
    alpha = 0.9
  ) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.y.npc = 0.98
  ) +
  facet_wrap(~ Facet_label, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_x_discrete(labels = group_labels) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 9),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(y = "Feature value")

p_cazy_sig_all

# ggsave(
#   "figures/CAZy_selected_module_class_family_violin.svg",
#   p_cazy_sig_all,
#   width = 10,
#   height = 6,
#   device = "svg"
# )


# ---------- 11. Median-difference dot plot ----------

p_cazy_sig_median_diff <- cazy_sig_median_difference %>%
  dplyr::mutate(
    Feature_label = factor(Feature_label, levels = rev(cazy_sig_result$Feature_label)),
    neg_log10_p = -log10(pmax(p_value, .Machine$double.xmin))
  ) %>%
  ggplot(aes(x = median_diff_pCR_minus_non_pCR, y = Feature_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_point(aes(size = neg_log10_p), color = "black", alpha = 0.85) +
  facet_grid(
    Feature_level ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Median difference: pCR - non-pCR",
    y = NULL,
    size = "-log10(p)"
  )

p_cazy_sig_median_diff

# ggsave(
#   "figures/CAZy_selected_module_class_family_median_difference.svg",
#   p_cazy_sig_median_diff,
#   width = 6.5,
#   height = 5,
#   device = "svg"
# )





# ==================== CAZyme profile-level analysis ====================

group_colors <- c("pCR" = "#4FAE9A", "non_pCR" = "#DE7872")
group_labels <- c("pCR" = "pCR", "non_pCR" = "non-pCR")

dir.create("figures", showWarnings = FALSE)

m_before <- m_before %>%
  dplyr::mutate(
    Response = factor(Response, levels = c("pCR", "non_pCR"))
  ) %>%
  dplyr::arrange(Response, SampleID)

# ---------- 1. CAZy family matrix ----------
# Rows are samples and columns are CAZy families.
# Only baseline samples are used.

cazy_family_before_mat <- cazy_family %>%
  dplyr::select(Feature, dplyr::all_of(m_before$SampleID))

rownames(cazy_family_before_mat) <- cazy_family_before_mat$Feature
cazy_family_before_mat$Feature <- NULL

cazy_family_before_mat <- t(as.matrix(cazy_family_before_mat))
cazy_family_before_mat[is.na(cazy_family_before_mat)] <- 0

cazy_family_before_mat <- cazy_family_before_mat[m_before$SampleID, ]

# Relative abundance version for distance-based analyses.
# This reduces the influence of total CAZyme load per sample.

cazy_family_before_rel <- sweep(
  cazy_family_before_mat,
  1,
  rowSums(cazy_family_before_mat),
  "/"
)

cazy_family_before_rel[is.na(cazy_family_before_rel)] <- 0
cazy_family_before_rel[is.infinite(cazy_family_before_rel)] <- 0


# ---------- 2-1. CAZyme PCoA: Bray-Curtis ----------
# Bray-Curtis is used as the main analysis because it is commonly used
# for microbiome compositional profiles and is less dominated by rare features
# than Canberra distance.

cazy_dist_bray <- vegan::vegdist(cazy_family_before_rel, method = "bray")

cazy_permanova_bray <- vegan::adonis2(
  cazy_dist_bray ~ Response,
  data = m_before,
  permutations = 9999
)

cazy_permanova_bray

cazy_pcoa_bray <- cmdscale(
  cazy_dist_bray,
  eig = TRUE,
  k = 2
)

cazy_pcoa_bray_df <- data.frame(
  SampleID = rownames(cazy_pcoa_bray$points),
  PCo1 = cazy_pcoa_bray$points[, 1],
  PCo2 = cazy_pcoa_bray$points[, 2]
) %>%
  dplyr::left_join(
    m_before %>%
      dplyr::select(SampleID, Response, TRG_score, TRG_1),
    by = "SampleID"
  )

cazy_pcoa_bray_var <- round(
  100 * cazy_pcoa_bray$eig[1:2] / sum(cazy_pcoa_bray$eig[cazy_pcoa_bray$eig > 0]),
  1
)

cazy_bray_R2 <- round(cazy_permanova_bray$R2[1], 3)
cazy_bray_p <- signif(cazy_permanova_bray$`Pr(>F)`[1], 3)

p_cazy_pcoa_bray <- cazy_pcoa_bray_df %>%
  ggplot(aes(x = PCo1, y = PCo2, color = Response)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(linewidth = 0.6, alpha = 0.8) +
  scale_color_manual(values = group_colors, labels = group_labels) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  labs(
    x = paste0("PCo1 (", cazy_pcoa_bray_var[1], "%)"),
    y = paste0("PCo2 (", cazy_pcoa_bray_var[2], "%)"),
    color = NULL,
    title = "CAZyme profile PCoA",
    subtitle = paste0("Bray-Curtis PERMANOVA: R² = ", cazy_bray_R2, ", p = ", cazy_bray_p)
  )

p_cazy_pcoa_bray

# ggsave(
#   "figures/CAZy_PCoA_BrayCurtis_pCR_non_pCR.svg",
#   p_cazy_pcoa_bray,
#   width = 6,
#   height = 5,
#   device = "svg"
# )

# # ---------- 2-2. CAZyme PCoA: Canberra sensitivity ----------
# # Canberra distance is included as a sensitivity analysis because the Cayman paper
# # used Canberra distance for CAZyme profile ordination.
# 
# cazy_dist_canberra <- vegan::vegdist(cazy_family_before_rel, method = "canberra")
# 
# cazy_permanova_canberra <- vegan::adonis2(
#   cazy_dist_canberra ~ Response,
#   data = m_before,
#   permutations = 9999
# )
# 
# cazy_permanova_canberra
# 
# cazy_pcoa_canberra <- cmdscale(
#   cazy_dist_canberra,
#   eig = TRUE,
#   k = 2
# )
# 
# cazy_pcoa_canberra_df <- data.frame(
#   SampleID = rownames(cazy_pcoa_canberra$points),
#   PCo1 = cazy_pcoa_canberra$points[, 1],
#   PCo2 = cazy_pcoa_canberra$points[, 2]
# ) %>%
#   dplyr::left_join(
#     m_before %>%
#       dplyr::select(SampleID, Response, TRG_score, TRG_1),
#     by = "SampleID"
#   )
# 
# cazy_pcoa_canberra_var <- round(
#   100 * cazy_pcoa_canberra$eig[1:2] / sum(cazy_pcoa_canberra$eig[cazy_pcoa_canberra$eig > 0]),
#   1
# )
# 
# cazy_canberra_R2 <- round(cazy_permanova_canberra$R2[1], 3)
# cazy_canberra_p <- signif(cazy_permanova_canberra$`Pr(>F)`[1], 3)
# 
# p_cazy_pcoa_canberra <- cazy_pcoa_canberra_df %>%
#   ggplot(aes(x = PCo1, y = PCo2, color = Response)) +
#   geom_point(size = 3, alpha = 0.9) +
#   stat_ellipse(linewidth = 0.6, alpha = 0.8) +
#   scale_color_manual(values = group_colors, labels = group_labels) +
#   theme_classic(base_size = 12) +
#   theme(
#     legend.position = "right",
#     panel.grid = element_blank()
#   ) +
#   labs(
#     x = paste0("PCo1 (", cazy_pcoa_canberra_var[1], "%)"),
#     y = paste0("PCo2 (", cazy_pcoa_canberra_var[2], "%)"),
#     color = NULL,
#     title = "CAZyme profile PCoA",
#     subtitle = paste0("Canberra PERMANOVA: R² = ", cazy_canberra_R2, ", p = ", cazy_canberra_p)
#   )
# 
# p_cazy_pcoa_canberra
# 
# ggsave(
#   "figures/CAZy_PCoA_Canberra_pCR_non_pCR.svg",
#   p_cazy_pcoa_canberra,
#   width = 6,
#   height = 5,
#   device = "svg"
# )



# ==================== CAZyme richness ====================

# Richness is defined as the number of CAZy families detected above a TPM threshold.
# The Cayman paper used >1 RPKM for richness; here TPM is used, so the threshold
# should be treated as an analysis choice rather than an identical replication.

cazy_richness_threshold <- 0.01

cazy_richness <- data.frame(
  SampleID = rownames(cazy_family_before_mat),
  CAZyme_richness = rowSums(cazy_family_before_mat > cazy_richness_threshold)
) %>%
  dplyr::left_join(
    m_before %>%
      dplyr::select(SampleID, Response, TRG_1, TRG_score),
    by = "SampleID"
  )

cazy_richness_result <- cazy_richness %>%
  dplyr::summarise(
    n_pCR = sum(Response == "pCR"),
    n_non_pCR = sum(Response == "non_pCR"),
    median_pCR = median(CAZyme_richness[Response == "pCR"], na.rm = TRUE),
    median_non_pCR = median(CAZyme_richness[Response == "non_pCR"], na.rm = TRUE),
    median_diff_pCR_minus_non_pCR = median_pCR - median_non_pCR,
    p_value = wilcox.test(CAZyme_richness ~ Response, exact = FALSE)$p.value
  )

cazy_richness_result

write_tsv(
  cazy_richness,
  "260224 final Input file/CAZyme/CAZy_family_richness_before.tsv"
)

write_tsv(
  cazy_richness_result,
  "260224 final Input file/CAZyme/CAZy_family_richness_CR_nonCR_result.tsv"
)

p_cazy_richness <- cazy_richness %>%
  ggplot(aes(x = Response, y = CAZyme_richness, fill = Response, color = Response)) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    fill = NA,
    linewidth = 0.6
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.4,
    alpha = 0.9
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_x_discrete(labels = group_labels) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(y = "Number of detected CAZy families")

p_cazy_richness

# ggsave(
#   "figures/CAZy_family_richness_before_pCR_non_pCR.svg",
#   p_cazy_richness,
#   width = 4,
#   height = 4.5,
#   device = "svg"
# )


# ==================== CAZyme substrate-ratio analysis ====================

# These family sets are proxy definitions.
# Replace them with Cayman substrate annotation if available.

mucin_targeting_families <- c(
  "GH2", "GH20", "GH29", "GH33", "GH35",
  "GH84", "GH89", "GH95", "GH101",
  "GH109", "GH110", "GH112", "GH123",
  "CE4"
)

dietary_fiber_targeting_families <- c(
  "GH3", "GH5", "GH10", "GH11", "GH16",
  "GH26", "GH28", "GH30", "GH43", "GH51",
  "GH53", "GH67", "GH74", "GH78",
  "GH105", "GH106",
  "PL1", "PL7", "PL9", "PL11", "PL15",
  "CE8", "CE12",
  "CBM6", "CBM13", "CBM32"
)

gag_targeting_families <- c(
  "PL8", "PL12", "PL13", "PL16", "PL21",
  "GH88", "GH89", "GH105",
  "CBM32", "CBM40"
)

cazy_substrate_ratio <- data.frame(
  SampleID = rownames(cazy_family_before_mat),
  dietary_fiber_TPM = rowSums(
    cazy_family_before_mat[, intersect(dietary_fiber_targeting_families, colnames(cazy_family_before_mat)), drop = FALSE]
  ),
  mucin_TPM = rowSums(
    cazy_family_before_mat[, intersect(mucin_targeting_families, colnames(cazy_family_before_mat)), drop = FALSE]
  ),
  GAG_TPM = rowSums(
    cazy_family_before_mat[, intersect(gag_targeting_families, colnames(cazy_family_before_mat)), drop = FALSE]
  )
) %>%
  dplyr::mutate(
    mucin_to_fiber_ratio = (mucin_TPM + 1) / (dietary_fiber_TPM + 1),
    GAG_to_fiber_ratio = (GAG_TPM + 1) / (dietary_fiber_TPM + 1),
    mucin_to_fiber_log2_ratio = log2(mucin_to_fiber_ratio),
    GAG_to_fiber_log2_ratio = log2(GAG_to_fiber_ratio)
  ) %>%
  dplyr::left_join(
    m_before %>%
      dplyr::select(SampleID, Response, TRG_1, TRG_score),
    by = "SampleID"
  )

# write_tsv(
#   cazy_substrate_ratio,
#   "260224 final Input file/CAZyme/CAZy_substrate_ratio_before.tsv"
# )

cazy_substrate_ratio_result <- cazy_substrate_ratio %>%
  tidyr::pivot_longer(
    cols = c(mucin_to_fiber_log2_ratio, GAG_to_fiber_log2_ratio),
    names_to = "Ratio",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Ratio = factor(
      Ratio,
      levels = c("mucin_to_fiber_log2_ratio", "GAG_to_fiber_log2_ratio"),
      labels = c("Mucin / dietary fiber", "GAG / dietary fiber")
    )
  ) %>%
  dplyr::group_by(Ratio) %>%
  dplyr::summarise(
    n_pCR = sum(Response == "pCR" & !is.na(Value)),
    n_non_pCR = sum(Response == "non_pCR" & !is.na(Value)),
    median_pCR = median(Value[Response == "pCR"], na.rm = TRUE),
    median_non_pCR = median(Value[Response == "non_pCR"], na.rm = TRUE),
    median_diff_pCR_minus_non_pCR = median_pCR - median_non_pCR,
    p_value = wilcox.test(Value ~ Response, exact = FALSE)$p.value,
    .groups = "drop"
  )

cazy_substrate_ratio_result

# write_tsv(
#   cazy_substrate_ratio_result,
#   "260224 final Input file/CAZyme/CAZy_substrate_ratio_CR_nonCR_result.tsv"
# )

p_cazy_substrate_ratio <- cazy_substrate_ratio %>%
  tidyr::pivot_longer(
    cols = c(mucin_to_fiber_log2_ratio, GAG_to_fiber_log2_ratio),
    names_to = "Ratio",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Ratio = factor(
      Ratio,
      levels = c("mucin_to_fiber_log2_ratio", "GAG_to_fiber_log2_ratio"),
      labels = c("Mucin / dietary fiber", "GAG / dietary fiber")
    )
  ) %>%
  ggplot(aes(x = Response, y = Value, fill = Response, color = Response)) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    fill = NA,
    linewidth = 0.6
  ) +
  geom_jitter(
    width = 0.08,
    size = 2.4,
    alpha = 0.9
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  facet_wrap(~ Ratio, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_x_discrete(labels = group_labels) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    strip.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank()
  ) +
  labs(y = "log2 ratio")

p_cazy_substrate_ratio

# ggsave(
#   "figures/CAZy_mucin_GAG_to_fiber_ratio_before.svg",
#   p_cazy_substrate_ratio,
#   width = 6.2,
#   height = 4.5,
#   device = "svg"
# )





# ==================== CAZy subfamily significance + substrate target figure ====================

# Before samples only
m_before <- m_before %>%
  dplyr::mutate(
    Response = factor(Response, levels = c("pCR", "non_pCR"))
  ) %>%
  dplyr::arrange(Response, SampleID)

group_colors <- c("pCR" = "#4FAE9A", "non_pCR" = "#DE7872")
group_labels <- c("pCR" = "pCR", "non_pCR" = "non-pCR")

subfamily_prevalence_cutoff <- 0.30
subfamily_mean_tpm_cutoff <- 1
subfamily_p_cutoff <- 0.05
top_n_each_direction <- 10

dir.create("figures", showWarnings = FALSE)


# ---------- 1. Long-format subfamily table ----------

cazy_subfamily_before <- cazy_subfamily %>%
  dplyr::select(Feature, dplyr::all_of(m_before$SampleID)) %>%
  tidyr::pivot_longer(
    cols = -Feature,
    names_to = "SampleID",
    values_to = "TPM"
  ) %>%
  dplyr::inner_join(
    m_before %>%
      dplyr::select(SampleID, Response, TRG_1, TRG_score),
    by = "SampleID"
  ) %>%
  dplyr::mutate(
    log_TPM = log2(TPM + 1),
    Parent_family = stringr::str_remove(Feature, "_.*$")
  ) %>%
  dplyr::arrange(Feature, Response, SampleID)


# ---------- 2. Subfamily-level statistics ----------
# Ranking is based on p-value after prevalence and abundance filtering.
# Effect size is retained only for direction and optional reporting,
# not shown as a separate tile in the final figure.

cazy_subfamily_result <- cazy_subfamily_before %>%
  dplyr::group_by(Feature, Parent_family) %>%
  dplyr::summarise(
    n_pCR = sum(Response == "pCR" & !is.na(TPM)),
    n_non_pCR = sum(Response == "non_pCR" & !is.na(TPM)),
    prevalence = mean(TPM > 0),
    mean_TPM_all = mean(TPM, na.rm = TRUE),
    median_TPM_all = median(TPM, na.rm = TRUE),
    mean_TPM_pCR = mean(TPM[Response == "pCR"], na.rm = TRUE),
    mean_TPM_non_pCR = mean(TPM[Response == "non_pCR"], na.rm = TRUE),
    median_TPM_pCR = median(TPM[Response == "pCR"], na.rm = TRUE),
    median_TPM_non_pCR = median(TPM[Response == "non_pCR"], na.rm = TRUE),
    median_pairwise_log2_diff_pCR_vs_non_pCR = median(
      as.vector(
        outer(
          log_TPM[Response == "pCR"],
          log_TPM[Response == "non_pCR"],
          "-"
        )
      ),
      na.rm = TRUE
    ),
    p_value = dplyr::case_when(
      prevalence >= subfamily_prevalence_cutoff &
        mean_TPM_all >= subfamily_mean_tpm_cutoff ~
        wilcox.test(log_TPM ~ Response, exact = FALSE)$p.value,
      TRUE ~ NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    q_value = p.adjust(p_value, method = "BH"),
    neg_log10_p = -log10(pmax(p_value, .Machine$double.xmin)),
    Association = dplyr::case_when(
      median_pairwise_log2_diff_pCR_vs_non_pCR > 0 ~ "pCR-enriched",
      median_pairwise_log2_diff_pCR_vs_non_pCR < 0 ~ "non-pCR-enriched",
      TRUE ~ "Neutral"
    )
  ) %>%
  dplyr::arrange(p_value)

cazy_subfamily_result %>% 
  as.data.frame() %>% 
  filter(p_value < 0.1)


# write_tsv(
#   cazy_subfamily_result,
#   "260224 final Input file/CAZyme/CAZy_subfamily_result_before_filtered.tsv"
# )


# ---------- 3. Functional labels ----------
# These labels are broad parent-family-level annotations.
# They are intended for figure readability, not strict subfamily-level enzymatic proof.

cazy_subfamily_result <- cazy_subfamily_result %>%
  dplyr::mutate(
    Functional_label = dplyr::case_when(
      Parent_family == "GH1" ~ "β-glucosidase-like",
      Parent_family == "GH2" ~ "β-galactosidase / β-mannosidase-like",
      Parent_family == "GH13" ~ "α-amylase / starch-active",
      Parent_family == "GH25" ~ "lysozyme / muramidase-like",
      Parent_family == "GH31" ~ "α-glucosidase-like",
      Parent_family == "GH32" ~ "fructan hydrolase-like",
      Parent_family == "GH39" ~ "β-xylosidase-like",
      Parent_family == "GH43" ~ "xylan / arabinose-glycan-active",
      Parent_family == "GH53" ~ "endo-β-1,4-galactanase-like",
      Parent_family == "GH72" ~ "β-glucan remodeling-like",
      Parent_family == "GH78" ~ "α-L-rhamnosidase-like",
      Parent_family == "GH88" ~ "unsaturated glucuronyl hydrolase-like",
      Parent_family == "GH95" ~ "α-fucosidase-like",
      Parent_family == "GH105" ~ "unsaturated rhamnogalacturonyl hydrolase-like",
      Parent_family == "CE2" ~ "acetyl xylan esterase-like",
      Parent_family == "CE4" ~ "carbohydrate deacetylase-like",
      Parent_family == "CE12" ~ "pectin/acetyl esterase-like",
      Parent_family == "CBM34" ~ "starch-binding module",
      Parent_family == "CBM50" ~ "LysM / peptidoglycan-binding",
      Parent_family == "CBM68" ~ "carbohydrate-binding module",
      Parent_family == "CBM83" ~ "plant glycan-binding module",
      Parent_family == "GT4" ~ "glycosyltransferase family 4",
      Parent_family == "GT28" ~ "glycosyltransferase family 28",
      Parent_family == "GT51" ~ "peptidoglycan glycosyltransferase-like",
      Parent_family == "GT101" ~ "glycan biosynthesis/modification",
      Parent_family == "GT111" ~ "β-1,3-galactofuranosyltransferase-like",
      Parent_family == "PL12" ~ "GAG lyase-like",
      Parent_family == "AA4" ~ "vanillyl-alcohol oxidase-like",
      TRUE ~ paste0(Parent_family, " family-related")
    ),
    Feature_label = paste0(Feature, "\n", Functional_label)
  )


# ---------- 4. Substrate / functional-axis annotation ----------
# Categories with no selected feature will be removed later.
# GT families are assigned to glycan biosynthesis/modification, not to glycosaminoglycans.

cazy_subfamily_result <- cazy_subfamily_result %>%
  dplyr::mutate(
    Dietary_fiber = dplyr::case_when(
      Parent_family %in% c(
        "GH1", "GH3", "GH5", "GH8", "GH10", "GH11", "GH16",
        "GH26", "GH28", "GH30", "GH39", "GH43", "GH51", "GH53",
        "GH67", "GH72", "GH74", "GH78", "GH105", "GH106",
        "PL1", "PL4", "PL7", "PL8", "PL9", "PL10", "PL11",
        "PL13", "PL15", "PL17", "PL26", "PL30", "PL32", "PL35",
        "CE8", "CE12",
        "CBM3", "CBM6", "CBM13", "CBM22", "CBM23",
        "CBM27", "CBM32", "CBM44", "CBM65", "CBM77",
        "CBM79", "CBM83", "CBM86"
      ) ~ "Yes",
      TRUE ~ "No"
    ),
    Mucin = dplyr::case_when(
      Parent_family %in% c(
        "GH2", "GH20", "GH29", "GH33", "GH35",
        "GH84", "GH89", "GH92", "GH95", "GH99",
        "GH101", "GH109", "GH110", "GH112", "GH123",
        "CE4", "CE14",
        "CBM40", "CBM51", "CBM62"
      ) ~ "Yes",
      TRUE ~ "No"
    ),
    Peptidoglycan = dplyr::case_when(
      Parent_family %in% c(
        "GH18", "GH19", "GH23", "GH24", "GH25",
        "GH73", "GH102", "GH103", "GH104", "GH108",
        "GH170", "CBM50", "GT51"
      ) ~ "Yes",
      TRUE ~ "No"
    ),
    Glycogen_starch_alpha_glucan = dplyr::case_when(
      Parent_family %in% c(
        "GH13", "GH14", "GH15", "GH31", "GH57",
        "GH65", "GH77", "GH97",
        "CBM20", "CBM25", "CBM26", "CBM34", "CBM48", "CBM58"
      ) ~ "Yes",
      TRUE ~ "No"
    ),
    Fructan_inulin = dplyr::case_when(
      Parent_family %in% c("GH32", "GH68", "CBM38") ~ "Yes",
      TRUE ~ "No"
    ),
    Acetylated_plant_polysaccharides = dplyr::case_when(
      Parent_family %in% c("CE2", "CE8", "CE12") ~ "Yes",
      TRUE ~ "No"
    ),
    Glycan_biosynthesis_modification = dplyr::case_when(
      stringr::str_detect(Parent_family, "^GT") ~ "Yes",
      TRUE ~ "No"
    ),
    Auxiliary_redox_activity = dplyr::case_when(
      stringr::str_detect(Parent_family, "^AA") ~ "Yes",
      TRUE ~ "No"
    )
  )


# ---------- 5. Select features using both p-value and effect size ----------
# Selection logic:
#   1) remove sparse and low-abundance features
#   2) require nominal significance and minimum effect size
#   3) rank features by both p-value and absolute effect size
#   4) select up to top 10 pCR-enriched and top 10 non-pCR-enriched features
#
# Effect size:
#   median_pairwise_log2_diff_pCR_vs_non_pCR
#   positive = higher in pCR
#   negative = higher in non-pCR

subfamily_prevalence_cutoff <- 0.30
subfamily_mean_tpm_cutoff <- 1
subfamily_p_cutoff <- 0.05
subfamily_effect_cutoff <- 0.5
top_n_each_direction <- 10

cazy_subfamily_selected <- dplyr::bind_rows(
  cazy_subfamily_result %>%
    dplyr::filter(
      !is.na(p_value),
      p_value < subfamily_p_cutoff,
      prevalence >= subfamily_prevalence_cutoff,
      mean_TPM_all >= subfamily_mean_tpm_cutoff,
      abs(median_pairwise_log2_diff_pCR_vs_non_pCR) >= subfamily_effect_cutoff,
      Association == "pCR-enriched"
    ) %>%
    dplyr::mutate(
      p_rank = dplyr::dense_rank(p_value),
      effect_rank = dplyr::dense_rank(
        dplyr::desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR))
      ),
      selection_rank = p_rank + effect_rank,
      selection_score = -log10(pmax(p_value, .Machine$double.xmin)) *
        abs(median_pairwise_log2_diff_pCR_vs_non_pCR)
    ) %>%
    dplyr::arrange(
      selection_rank,
      p_value,
      dplyr::desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR)),
      dplyr::desc(prevalence),
      dplyr::desc(mean_TPM_all),
      Feature
    ) %>%
    dplyr::slice_head(n = top_n_each_direction),
  
  cazy_subfamily_result %>%
    dplyr::filter(
      !is.na(p_value),
      p_value < subfamily_p_cutoff,
      prevalence >= subfamily_prevalence_cutoff,
      mean_TPM_all >= subfamily_mean_tpm_cutoff,
      abs(median_pairwise_log2_diff_pCR_vs_non_pCR) >= subfamily_effect_cutoff,
      Association == "non-pCR-enriched"
    ) %>%
    dplyr::mutate(
      p_rank = dplyr::dense_rank(p_value),
      effect_rank = dplyr::dense_rank(
        dplyr::desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR))
      ),
      selection_rank = p_rank + effect_rank,
      selection_score = -log10(pmax(p_value, .Machine$double.xmin)) *
        abs(median_pairwise_log2_diff_pCR_vs_non_pCR)
    ) %>%
    dplyr::arrange(
      selection_rank,
      p_value,
      dplyr::desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR)),
      dplyr::desc(prevalence),
      dplyr::desc(mean_TPM_all),
      Feature
    ) %>%
    dplyr::slice_head(n = top_n_each_direction)
) %>%
  dplyr::distinct(Feature, .keep_all = TRUE) %>%
  dplyr::arrange(
    Association,
    selection_rank,
    p_value,
    dplyr::desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR)),
    Feature
  ) %>%
  dplyr::mutate(
    Feature_plot = factor(Feature, levels = Feature),
    Feature_label = factor(Feature_label, levels = Feature_label)
  )

cazy_subfamily_selected %>%
  dplyr::count(Association)

cazy_subfamily_selected %>%
  dplyr::select(
    Feature,
    Functional_label,
    Association,
    prevalence,
    mean_TPM_all,
    median_pairwise_log2_diff_pCR_vs_non_pCR,
    p_value,
    q_value,
    p_rank,
    effect_rank,
    selection_rank,
    selection_score
  )

# write_tsv(
#   cazy_subfamily_selected,
#   "260224 final Input file/CAZyme/CAZy_subfamily_selected_pvalue_effect_rank.tsv"
# )


# ---------- 6. Keep only substrate categories observed in selected features ----------
# Empty annotation rows are removed automatically.
# For example, Fructan/inulin or Glycosaminoglycans will appear only if at least one selected feature is assigned to them.

cazy_substrate_columns <- c(
  "Dietary_fiber",
  "Mucin",
  "Peptidoglycan",
  "Glycogen_starch_alpha_glucan",
  "Fructan_inulin",
  "Acetylated_plant_polysaccharides",
  "Glycan_biosynthesis_modification",
  "Auxiliary_redox_activity"
)

cazy_substrate_columns <- cazy_substrate_columns[
  colSums(cazy_subfamily_selected[, cazy_substrate_columns] == "Yes") > 0
]

cazy_substrate_columns


# ---------- 7. Statistical-significance bar plot ----------
# Bar height represents -log10(P) from the Wilcoxon rank-sum test.
# The dashed line indicates nominal P = 0.05, and colors indicate the
# direction of the median pairwise abundance difference.

p_cazy_subfamily_signif <- cazy_subfamily_selected %>%
  ggplot(
    aes(
      x = Feature_plot,
      y = neg_log10_p,
      fill = Association
    )
  ) +
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.25
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.45
  ) +
  scale_fill_manual(
    values = c(
      "pCR-enriched" = group_colors[["pCR"]],
      "non-pCR-enriched" = group_colors[["non_pCR"]]
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = NULL,
    y = expression(-log[10](italic(P))),
    fill = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 9),
    legend.position = "top",
    panel.grid = element_blank(),
    plot.margin = margin(5.5, 5.5, 0, 5.5)
  )

p_cazy_subfamily_signif


# ---------- 8. Substrate / functional-axis tile ----------
# Only categories represented among selected features are shown.

cazy_subfamily_substrate_tile <- cazy_subfamily_selected %>%
  dplyr::select(
    Feature_plot,
    Feature_label,
    dplyr::all_of(cazy_substrate_columns)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(cazy_substrate_columns),
    names_to = "Substrate",
    values_to = "Target"
  ) %>%
  dplyr::mutate(
    Substrate = factor(
      Substrate,
      levels = cazy_substrate_columns,
      labels = dplyr::case_when(
        cazy_substrate_columns == "Dietary_fiber" ~ "Dietary\nfiber",
        cazy_substrate_columns == "Mucin" ~ "Mucin",
        cazy_substrate_columns == "Peptidoglycan" ~ "Peptido-\nglycan",
        cazy_substrate_columns == "Glycogen_starch_alpha_glucan" ~ "Glycogen/\nstarch/\nα-glucan",
        cazy_substrate_columns == "Fructan_inulin" ~ "Fructan /\ninulin",
        cazy_substrate_columns == "Acetylated_plant_polysaccharides" ~ "Acetylated\nplant polysac.",
        cazy_substrate_columns == "Glycan_biosynthesis_modification" ~ "Glycan\nbiosyn./\nmodif.",
        cazy_substrate_columns == "Auxiliary_redox_activity" ~ "Auxiliary\nredox",
        TRUE ~ cazy_substrate_columns
      )
    ),
    Feature_plot = factor(
      Feature_plot,
      levels = levels(cazy_subfamily_selected$Feature_plot)
    ),
    Target = factor(Target, levels = c("No", "Yes"))
  )

p_cazy_subfamily_substrate_tile <- cazy_subfamily_substrate_tile %>%
  ggplot(aes(x = Feature_plot, y = Substrate, fill = Target)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_manual(
    values = c(
      "No" = "#F6F6F6",
      "Yes" = "#2F5D7E"
    )
  ) +
  scale_x_discrete(
    labels = setNames(
      as.character(cazy_subfamily_selected$Feature_label),
      as.character(cazy_subfamily_selected$Feature_plot)
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "right",
    axis.title = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 8,
      lineheight = 0.9
    ),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  ) +
  labs(fill = "Functional\naxis")

p_cazy_subfamily_substrate_tile


# ---------- 9. Combine panels ----------
# Top: -log10(p) bar plot
# Bottom: functional-axis tile
# Effect-size tile is intentionally removed.

p_cazy_subfamily_signif_compact <- ggpubr::ggarrange(
  p_cazy_subfamily_signif,
  p_cazy_subfamily_substrate_tile,
  ncol = 1,
  heights = c(2.5, 1.9),
  align = "v"
)

p_cazy_subfamily_signif_compact

# ggsave(
#   "figures/CAZy_subfamily_significance_functional_axis_before.svg",
#   p_cazy_subfamily_signif_compact,
#   width = 10.5,
#   height = 6.8,
#   device = "svg"
# )






# ---------- 1. Feature selection ----------
# Selection strategy:
# 1) remove sparse / very low-abundance features
# 2) require nominal significance and minimum effect size
# 3) rank within each direction using both p-value and effect size
# 4) select top features separately for pCR-enriched and non-pCR-enriched groups
#
# Positive effect size = enriched in pCR
# Negative effect size = enriched in non-pCR

subfamily_prevalence_cutoff <- 0.30
subfamily_mean_tpm_cutoff <- 1
subfamily_p_cutoff <- 0.10
subfamily_effect_cutoff <- 0.30
top_n_each_direction <- 8

cazy_subfamily_selected <- dplyr::bind_rows(
  cazy_subfamily_result %>%
    dplyr::filter(
      !is.na(p_value),
      p_value < subfamily_p_cutoff,
      prevalence >= subfamily_prevalence_cutoff,
      mean_TPM_all >= subfamily_mean_tpm_cutoff,
      abs(median_pairwise_log2_diff_pCR_vs_non_pCR) >= subfamily_effect_cutoff,
      Association == "pCR-enriched"
    ) %>%
    dplyr::mutate(
      p_rank = dplyr::dense_rank(p_value),
      effect_rank = dplyr::dense_rank(
        dplyr::desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR))
      ),
      selection_rank = p_rank + effect_rank
    ) %>%
    dplyr::arrange(
      selection_rank,
      p_value,
      dplyr::desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR)),
      dplyr::desc(prevalence),
      dplyr::desc(mean_TPM_all),
      Feature
    ) %>%
    dplyr::slice_head(n = top_n_each_direction),
  
  cazy_subfamily_result %>%
    dplyr::filter(
      !is.na(p_value),
      p_value < subfamily_p_cutoff,
      prevalence >= subfamily_prevalence_cutoff,
      mean_TPM_all >= subfamily_mean_tpm_cutoff,
      abs(median_pairwise_log2_diff_pCR_vs_non_pCR) >= subfamily_effect_cutoff,
      Association == "non-pCR-enriched"
    ) %>%
    dplyr::mutate(
      p_rank = dplyr::dense_rank(p_value),
      effect_rank = dplyr::dense_rank(
        dplyr::desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR))
      ),
      selection_rank = p_rank + effect_rank
    ) %>%
    dplyr::arrange(
      selection_rank,
      p_value,
      dplyr::desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR)),
      dplyr::desc(prevalence),
      dplyr::desc(mean_TPM_all),
      Feature
    ) %>%
    dplyr::slice_head(n = top_n_each_direction)
) %>%
  dplyr::distinct(Feature, .keep_all = TRUE) %>%
  dplyr::mutate(
    Association = factor(
      Association,
      levels = c("pCR-enriched", "non-pCR-enriched")
    )
  ) %>%
  dplyr::arrange(
    Association,
    selection_rank,
    p_value,
    dplyr::desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR)),
    Feature
  ) %>%
  dplyr::mutate(
    Feature_plot = factor(Feature, levels = Feature),
    Feature_label = factor(
      paste0(Feature, "\n", Functional_label),
      levels = paste0(Feature, "\n", Functional_label)
    )
  )

cazy_subfamily_selected %>%
  dplyr::count(Association)

write_tsv(
  cazy_subfamily_selected,
  "260224 final Input file/CAZyme/CAZy_subfamily_selected_tilepanel.tsv"
)

# ---------- 2. Annotation rows actually represented in selected features ----------

cazy_substrate_columns <- c(
  "Dietary_fiber",
  "Mucin",
  "Peptidoglycan",
  "Glycogen_starch_alpha_glucan",
  "Fructan_inulin",
  "Acetylated_plant_polysaccharides",
  "Glycan_biosynthesis_modification",
  "Auxiliary_redox_activity"
)

cazy_substrate_columns <- cazy_substrate_columns[cazy_substrate_columns %in% colnames(cazy_subfamily_selected)]

cazy_substrate_columns <- cazy_substrate_columns[
  colSums(cazy_subfamily_selected[, cazy_substrate_columns] == "Yes", na.rm = TRUE) > 0
]

cazy_substrate_columns

# ---------- 3. -log10(p) tile ----------

p_cazy_subfamily_p_tile <- cazy_subfamily_selected %>%
  ggplot(aes(x = Feature_plot, y = "-log10(p)", fill = neg_log10_p)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient(
    low = "#F7F7F7",
    high = "#3E5C76"
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "right",
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank()
  ) +
  labs(fill = expression(-log[10](italic(p))))

p_cazy_subfamily_p_tile

# ---------- 4. Effect-size tile ----------
# Positive values indicate higher abundance in pCR.
# Negative values indicate higher abundance in non-pCR.

p_cazy_subfamily_effect_tile <- cazy_subfamily_selected %>%
  ggplot(
    aes(
      x = Feature_plot,
      y = "Median pairwise\nlog2 difference",
      fill = median_pairwise_log2_diff_pCR_vs_non_pCR
    )
  ) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient2(
    low = "#DE7872",
    mid = "#FAFAFA",
    high = "#4FAE9A",
    midpoint = 0
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "right",
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank()
  ) +
  labs(fill = "Median pairwise\nlog2 difference")

p_cazy_subfamily_effect_tile


# ---------- 6. Combine panels ----------

p_cazy_subfamily_tilepanel <- ggpubr::ggarrange(
  p_cazy_subfamily_p_tile,
  p_cazy_subfamily_effect_tile,
  p_cazy_subfamily_substrate_tile,
  ncol = 1,
  heights = c(0.55, 0.55, 2.15),
  align = "v"
)

p_cazy_subfamily_tilepanel

ggsave(
  "figures/CAZy_subfamily_tilepanel_pvalue_effect_annotation.svg",
  p_cazy_subfamily_tilepanel,
  width = 9.5,
  height = 5.8,
  device = "svg"
)
Sys.setFileTime(
  "figures/CAZy_subfamily_tilepanel_pvalue_effect_annotation.svg",
  Sys.time()
)




# ---------- CAZyme chord diagram ----------
# This panel summarizes how many differential CAZyme subfamilies
# belong to each CAZy class in pCR-enriched and non-pCR-enriched groups.

library(circlize)

cazy_chord <- cazy_subfamily_result %>%
  dplyr::filter(
    !is.na(p_value),
    prevalence >= 0.2,
    mean_TPM_all >= 1,
    p_value < 0.05,
    abs(median_pairwise_log2_diff_pCR_vs_non_pCR) >= 0.30
  ) %>%
  dplyr::mutate(
    Direction = dplyr::case_when(
      median_pairwise_log2_diff_pCR_vs_non_pCR > 0 ~ "pCR",
      median_pairwise_log2_diff_pCR_vs_non_pCR < 0 ~ "non-pCR"
    ),
    Parent_family = stringr::str_remove(Feature, "_.*$"),
    CAZy_class = stringr::str_extract(Parent_family, "^[A-Z]+")
  ) %>%
  dplyr::filter(CAZy_class %in% c("GH", "GT", "PL", "CE", "CBM", "AA")) %>%
  dplyr::count(Direction, CAZy_class, name = "n")


svg("figures/CAZy_chord_pCR_non_pCR.svg", width = 7, height = 5)

circos.clear()
circos.par(gap.degree = c(8, 2, 2, 2, 2, 2, 2, 8))

chordDiagram(
  x = cazy_chord,
  order = c("pCR", "non-pCR", "GH", "PL", "CBM", "GT", "CE", "AA"),
  grid.col = c(
    "pCR" = "#4FAE9A",
    "non-pCR" = "#DE7872",
    "GH" = "#5DBB8A",
    "PL" = "#32A4A6",
    "CBM" = "#3F78A8",
    "GT" = "#5967AE",
    "CE" = "#7B6AAE",
    "AA" = "#8C5FA8"
  ),
  transparency = 0.45,
  annotationTrack = "grid",
  preAllocateTracks = 1
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    circos.text(
      x = mean(xlim),
      y = ylim[1] + 0.2,
      labels = sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.9
    )
  },
  bg.border = NA
)

dev.off()

table(m_before$TRG_1)


# ---------- Check differential CAZyme feature counts used for chord plot ----------

cazy_chord_features <- cazy_subfamily_result %>%
  dplyr::filter(
    !is.na(p_value),
    prevalence >= 0.5,
    mean_TPM_all >= 1,
    p_value < 0.10,
    abs(median_pairwise_log2_diff_pCR_vs_non_pCR) >= 0.20
  ) %>%
  dplyr::mutate(
    Direction = dplyr::case_when(
      median_pairwise_log2_diff_pCR_vs_non_pCR > 0 ~ "pCR",
      median_pairwise_log2_diff_pCR_vs_non_pCR < 0 ~ "non-pCR"
    ),
    Parent_family = stringr::str_remove(Feature, "_.*$"),
    CAZy_class = stringr::str_extract(Parent_family, "^[A-Z]+")
  ) %>%
  dplyr::filter(
    !is.na(Direction),
    CAZy_class %in% c("GH", "GT", "PL", "CE", "CBM", "AA")
  )

# Total number of differential features by direction
cazy_chord_features %>%
  dplyr::count(Direction)

# Number of differential features by direction and CAZy class
cazy_chord_features %>%
  dplyr::count(Direction, CAZy_class) %>%
  tidyr::pivot_wider(
    names_from = CAZy_class,
    values_from = n,
    values_fill = 0
  )

# Inspect individual features
cazy_chord_features %>%
  dplyr::select(
    Feature,
    Functional_label,
    Direction,
    CAZy_class,
    prevalence,
    mean_TPM_all,
    median_pairwise_log2_diff_pCR_vs_non_pCR,
    p_value,
    q_value
  ) %>%
  dplyr::arrange(Direction, p_value)


# ---------- Global CAZyme burden check: subfamily level ----------

cazy_subfamily_global <- cazy_subfamily_before %>%
  dplyr::group_by(SampleID, Response) %>%
  dplyr::summarise(
    total_CAZy_TPM = sum(TPM, na.rm = TRUE),
    detected_CAZy_features = sum(TPM > 0, na.rm = TRUE),
    mean_log_TPM = mean(log_TPM, na.rm = TRUE),
    median_log_TPM = median(log_TPM, na.rm = TRUE),
    .groups = "drop"
  )

cazy_subfamily_global %>%
  dplyr::group_by(Response) %>%
  dplyr::summarise(
    n = dplyr::n(),
    median_total_CAZy_TPM = median(total_CAZy_TPM),
    median_detected_CAZy_features = median(detected_CAZy_features),
    median_mean_log_TPM = median(mean_log_TPM),
    .groups = "drop"
  )

wilcox.test(total_CAZy_TPM ~ Response, data = cazy_subfamily_global, exact = FALSE)
wilcox.test(detected_CAZy_features ~ Response, data = cazy_subfamily_global, exact = FALSE)
wilcox.test(mean_log_TPM ~ Response, data = cazy_subfamily_global, exact = FALSE)


# ---------- Directional skew check ----------

cazy_subfamily_result %>%
  dplyr::filter(
    !is.na(median_pairwise_log2_diff_pCR_vs_non_pCR),
    prevalence >= 0.30,
    mean_TPM_all >= 1
  ) %>%
  dplyr::summarise(
    n_features = dplyr::n(),
    n_pCR_higher = sum(median_pairwise_log2_diff_pCR_vs_non_pCR > 0),
    n_non_pCR_higher = sum(median_pairwise_log2_diff_pCR_vs_non_pCR < 0),
    median_effect = median(median_pairwise_log2_diff_pCR_vs_non_pCR, na.rm = TRUE),
    mean_effect = mean(median_pairwise_log2_diff_pCR_vs_non_pCR, na.rm = TRUE)
  )

p_cazy_effect_distribution <- cazy_subfamily_result %>%
  dplyr::filter(
    !is.na(median_pairwise_log2_diff_pCR_vs_non_pCR),
    prevalence >= 0.30,
    mean_TPM_all >= 1
  ) %>%
  ggplot(aes(x = median_pairwise_log2_diff_pCR_vs_non_pCR)) +
  geom_histogram(bins = 40, color = "white") +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed") +
  labs(
    x = "Median pairwise log2 difference (pCR - non-pCR)",
    y = "Number of CAZyme subfamily features"
  ) +
  theme_classic(base_size = 12)

p_cazy_effect_distribution

# ggsave(
#   "figures/CAZy_subfamily_effect_size_distribution.svg",
#   p_cazy_effect_distribution,
#   width = 5.5,
#   height = 4.2,
#   device = "svg"
# )

# ---------- Chord feature count under different cutoffs ----------

cazy_chord_cutoff_check <- cazy_subfamily_result %>%
  dplyr::filter(
    !is.na(p_value),
    prevalence >= 0.30,
    mean_TPM_all >= 1
  ) %>%
  dplyr::mutate(
    Direction = dplyr::case_when(
      median_pairwise_log2_diff_pCR_vs_non_pCR > 0 ~ "pCR",
      median_pairwise_log2_diff_pCR_vs_non_pCR < 0 ~ "non-pCR"
    ),
    abs_effect = abs(median_pairwise_log2_diff_pCR_vs_non_pCR),
    cutoff_p_1 = p_value <= 1,
    cutoff_p_0_10 = p_value < 0.10,
    cutoff_p_0_05 = p_value < 0.05,
    cutoff_p_0_10_effect_0_30 = p_value < 0.10 & abs_effect >= 0.30,
    cutoff_p_0_05_effect_0_50 = p_value < 0.05 & abs_effect >= 0.50
  ) %>%
  dplyr::select(
    Feature,
    Direction,
    starts_with("cutoff_")
  ) %>%
  tidyr::pivot_longer(
    cols = starts_with("cutoff_"),
    names_to = "Cutoff",
    values_to = "Pass"
  ) %>%
  dplyr::filter(Pass) %>%
  dplyr::count(Cutoff, Direction) %>%
  tidyr::pivot_wider(
    names_from = Direction,
    values_from = n,
    values_fill = 0
  )

cazy_chord_cutoff_check

setdiff(m_before$SampleID, colnames(cazy_subfamily))
setdiff(colnames(cazy_subfamily)[-1], m_before$SampleID)

m_before %>%
  dplyr::count(TRG_1, Response)

stopifnot(
  all(
    c(
      "Feature",
      "median_TPM_pCR",
      "median_TPM_non_pCR",
      "median_pairwise_log2_diff_pCR_vs_non_pCR",
      "Association"
    ) %in% colnames(cazy_subfamily_result)
  )
)

cazy_subfamily_result %>%
  dplyr::select(
    Feature,
    median_TPM_pCR,
    median_TPM_non_pCR,
    median_pairwise_log2_diff_pCR_vs_non_pCR,
    Association
  ) %>%
  dplyr::arrange(desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR))) %>%
  head(20)


# ---------- Global CAZyme abundance check ----------

cazy_family_global <- cazy_family %>%
  dplyr::select(Feature, dplyr::all_of(m_before$SampleID)) %>%
  tidyr::pivot_longer(
    cols = -Feature,
    names_to = "SampleID",
    values_to = "TPM"
  ) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(
    total_CAZy_family_TPM = sum(TPM, na.rm = TRUE),
    detected_CAZy_families = sum(TPM > 0, na.rm = TRUE),
    mean_log_family_TPM = mean(log10(TPM + 1), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::inner_join(
    m_before %>%
      dplyr::select(SampleID, Response, TRG_1, TRG_score),
    by = "SampleID"
  )

cazy_family_global %>%
  dplyr::group_by(Response) %>%
  dplyr::summarise(
    n = dplyr::n(),
    median_total_CAZy_family_TPM = median(total_CAZy_family_TPM),
    median_detected_CAZy_families = median(detected_CAZy_families),
    median_mean_log_family_TPM = median(mean_log_family_TPM),
    .groups = "drop"
  )

wilcox.test(total_CAZy_family_TPM ~ Response, data = cazy_family_global, exact = FALSE)
wilcox.test(detected_CAZy_families ~ Response, data = cazy_family_global, exact = FALSE)
wilcox.test(mean_log_family_TPM ~ Response, data = cazy_family_global, exact = FALSE)

cazy_subfamily_global <- cazy_subfamily %>%
  dplyr::select(Feature, dplyr::all_of(m_before$SampleID)) %>%
  tidyr::pivot_longer(
    cols = -Feature,
    names_to = "SampleID",
    values_to = "TPM"
  ) %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(
    total_CAZy_subfamily_TPM = sum(TPM, na.rm = TRUE),
    detected_CAZy_subfamilies = sum(TPM > 0, na.rm = TRUE),
    mean_log_subfamily_TPM = mean(log10(TPM + 1), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::inner_join(
    m_before %>%
      dplyr::select(SampleID, Response, TRG_1, TRG_score),
    by = "SampleID"
  )

cazy_subfamily_global %>%
  dplyr::group_by(Response) %>%
  dplyr::summarise(
    n = dplyr::n(),
    median_total_CAZy_subfamily_TPM = median(total_CAZy_subfamily_TPM),
    median_detected_CAZy_subfamilies = median(detected_CAZy_subfamilies),
    median_mean_log_subfamily_TPM = median(mean_log_subfamily_TPM),
    .groups = "drop"
  )

wilcox.test(total_CAZy_subfamily_TPM ~ Response, data = cazy_subfamily_global, exact = FALSE)
wilcox.test(detected_CAZy_subfamilies ~ Response, data = cazy_subfamily_global, exact = FALSE)
wilcox.test(mean_log_subfamily_TPM ~ Response, data = cazy_subfamily_global, exact = FALSE)

# ---------- Directional skew of feature-wise effects ----------

cazy_subfamily_direction_check <- cazy_subfamily_result %>%
  dplyr::filter(
    !is.na(median_pairwise_log2_diff_pCR_vs_non_pCR),
    prevalence >= 0.30,
    mean_TPM_all >= 1
  ) %>%
  dplyr::mutate(
    Direction = dplyr::case_when(
      median_pairwise_log2_diff_pCR_vs_non_pCR > 0 ~ "pCR-higher",
      median_pairwise_log2_diff_pCR_vs_non_pCR < 0 ~ "non-pCR-higher",
      TRUE ~ "No difference"
    )
  )

cazy_subfamily_direction_check %>%
  dplyr::count(Direction)

cazy_subfamily_direction_check %>%
  dplyr::summarise(
    n_features = dplyr::n(),
    n_pCR_higher = sum(Direction == "pCR-higher"),
    n_non_pCR_higher = sum(Direction == "non-pCR-higher"),
    proportion_pCR_higher = n_pCR_higher / (n_pCR_higher + n_non_pCR_higher),
    median_effect = median(median_pairwise_log2_diff_pCR_vs_non_pCR, na.rm = TRUE),
    mean_effect = mean(median_pairwise_log2_diff_pCR_vs_non_pCR, na.rm = TRUE)
  )

binom.test(
  sum(cazy_subfamily_direction_check$Direction == "pCR-higher"),
  sum(cazy_subfamily_direction_check$Direction %in% c("pCR-higher", "non-pCR-higher")),
  p = 0.5
)

p_cazy_subfamily_effect_distribution <- cazy_subfamily_direction_check %>%
  ggplot(
    aes(x = median_pairwise_log2_diff_pCR_vs_non_pCR)
  ) +
  geom_histogram(
    bins = 40,
    color = "white"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  labs(
    x = "Median pairwise log2 difference (pCR - non-pCR)",
    y = "Number of CAZyme subfamily features"
  ) +
  theme_classic(base_size = 12)

p_cazy_subfamily_effect_distribution

ggsave(
  "figures/CAZy_subfamily_effect_direction_distribution.svg",
  p_cazy_subfamily_effect_distribution,
  width = 5.5,
  height = 4.2,
  device = "svg"
)
Sys.setFileTime(
  "figures/CAZy_subfamily_effect_direction_distribution.svg",
  Sys.time()
)



# ---------- CAZyme subfamily effect-size distribution colored by CAZy class ----------

cazy_subfamily_effect_class <- cazy_subfamily_result %>%
  dplyr::filter(
    !is.na(median_pairwise_log2_diff_pCR_vs_non_pCR),
    prevalence >= 0.30,
    mean_TPM_all >= 1
  ) %>%
  dplyr::mutate(
    Parent_family = stringr::str_remove(Feature, "_.*$"),
    CAZy_class = stringr::str_extract(Parent_family, "^[A-Z]+"),
    CAZy_class = factor(CAZy_class, levels = c("GH", "GT", "PL", "CE", "CBM", "AA"))
  ) %>%
  dplyr::filter(!is.na(CAZy_class))

cazy_class_colors <- c(
  "GH" = "#5DBB8A",
  "GT" = "#5E6FB1",
  "PL" = "#2FA7A0",
  "CE" = "#7A68B3",
  "CBM" = "#3F7FB3",
  "AA" = "#8A57A8"
)

p_cazy_subfamily_effect_hist_class <- cazy_subfamily_effect_class %>%
  ggplot(
    aes(
      x = median_pairwise_log2_diff_pCR_vs_non_pCR,
      fill = CAZy_class
    )
  ) +
  geom_histogram(
    bins = 45,
    color = "white",
    linewidth = 0.2,
    position = "stack"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  scale_fill_manual(values = cazy_class_colors) +
  labs(
    x = "Median pairwise log2 difference (pCR - non-pCR)",
    y = "Number of CAZyme subfamily features",
    fill = "CAZy class"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "right",
    panel.grid = element_blank()
  )

p_cazy_subfamily_effect_hist_class

ggsave(
  "figures/CAZy_subfamily_effect_size_distribution_by_class.svg",
  p_cazy_subfamily_effect_hist_class,
  width = 6.4,
  height = 4.6,
  device = "svg"
)
Sys.setFileTime(
  "figures/CAZy_subfamily_effect_size_distribution_by_class.svg",
  Sys.time()
)

p_cazy_subfamily_effect_hist_facet <- cazy_subfamily_effect_class %>%
  ggplot(
    aes(
      x = median_pairwise_log2_diff_pCR_vs_non_pCR,
      fill = CAZy_class
    )
  ) +
  geom_histogram(
    bins = 30,
    color = "white",
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  scale_fill_manual(values = cazy_class_colors) +
  facet_wrap(~ CAZy_class, ncol = 3, scales = "free_y") +
  labs(
    x = "Median pairwise log2 difference (pCR - non-pCR)",
    y = "Number of CAZyme subfamily features"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    legend.position = "none",
    panel.grid = element_blank()
  )

p_cazy_subfamily_effect_hist_facet

# ggsave(
#   "figures/CAZy_subfamily_effect_size_distribution_by_class_facet.svg",
#   p_cazy_subfamily_effect_hist_facet,
#   width = 7.2,
#   height = 5.6,
#   device = "svg"
# )

p_cazy_subfamily_effect_density_class <- cazy_subfamily_effect_class %>%
  ggplot(
    aes(
      x = median_pairwise_log2_diff_pCR_vs_non_pCR,
      color = CAZy_class
    )
  ) +
  geom_density(
    linewidth = 0.9,
    adjust = 1.1
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  scale_color_manual(values = cazy_class_colors) +
  labs(
    x = "Median pairwise log2 difference (pCR - non-pCR)",
    y = "Density",
    color = "CAZy class"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "right",
    panel.grid = element_blank()
  )

p_cazy_subfamily_effect_density_class

# ggsave(
#   "figures/CAZy_subfamily_effect_density_by_class.svg",
#   p_cazy_subfamily_effect_density_class,
#   width = 6.2,
#   height = 4.6,
#   device = "svg"
# )


# ==================== CAZyme subfamily effect-size histogram by CAZy class ====================

# This plot is descriptive.
# Each feature is a CAZyme subfamily.
# Positive effect size means higher abundance in pCR than in non-pCR.
# The annotation summarizes the overall directional skew across features.

# ---------- 1. Data preparation ----------

cazy_subfamily_effect_class <- cazy_subfamily_result %>%
  dplyr::filter(
    !is.na(median_pairwise_log2_diff_pCR_vs_non_pCR),
    prevalence >= 0.30,
    mean_TPM_all >= 1
  ) %>%
  dplyr::mutate(
    Parent_family = stringr::str_remove(Feature, "_.*$"),
    CAZy_class = stringr::str_extract(Parent_family, "^[A-Z]+")
  ) %>%
  dplyr::filter(
    CAZy_class %in% c("GH", "GT", "PL", "CE", "CBM", "AA")
  ) %>%
  dplyr::mutate(
    CAZy_class = factor(
      CAZy_class,
      levels = c("GH", "GT", "PL", "CE", "CBM", "AA"),
      labels = c(
        "GH: Glycoside hydrolases",
        "GT: Glycosyltransferases",
        "PL: Polysaccharide lyases",
        "CE: Carbohydrate esterases",
        "CBM: Carbohydrate-binding modules",
        "AA: Auxiliary activities"
      )
    )
  )

# ---------- 2. Overall directional-skew summary ----------
# This is a descriptive sign-based summary across CAZyme subfamily features.
# Because CAZyme features are not fully independent, the sign-test p-value should be interpreted descriptively.

cazy_subfamily_effect_overall <- cazy_subfamily_effect_class %>%
  dplyr::summarise(
    n_features = dplyr::n(),
    n_pCR_higher = sum(median_pairwise_log2_diff_pCR_vs_non_pCR > 0, na.rm = TRUE),
    n_non_pCR_higher = sum(median_pairwise_log2_diff_pCR_vs_non_pCR < 0, na.rm = TRUE),
    proportion_pCR_higher = n_pCR_higher / (n_pCR_higher + n_non_pCR_higher),
    median_effect = median(median_pairwise_log2_diff_pCR_vs_non_pCR, na.rm = TRUE),
    mean_effect = mean(median_pairwise_log2_diff_pCR_vs_non_pCR, na.rm = TRUE)
  )

cazy_subfamily_effect_sign_test <- binom.test(
  x = cazy_subfamily_effect_overall$n_pCR_higher,
  n = cazy_subfamily_effect_overall$n_pCR_higher + cazy_subfamily_effect_overall$n_non_pCR_higher,
  p = 0.5
)

cazy_subfamily_effect_annotation <- paste0(
  "Filtered subfamilies = ", cazy_subfamily_effect_overall$n_features, "\n",
  "Higher in pCR = ", cazy_subfamily_effect_overall$n_pCR_higher, "\n",
  "Higher in non-pCR = ", cazy_subfamily_effect_overall$n_non_pCR_higher, "\n",
  "Proportion higher in pCR = ", round(cazy_subfamily_effect_overall$proportion_pCR_higher, 3), "\n",
  "Median effect = ", round(cazy_subfamily_effect_overall$median_effect, 3), "\n",
  "Sign-test p = ", formatC(cazy_subfamily_effect_sign_test$p.value, format = "e", digits = 2)
)

# ---------- 3. Color palette ----------
# Colors are distinct but kept in a relatively soft publication-friendly tone.

cazy_class_colors <- c(
  "GH: Glycoside hydrolases" = "#5DBB8A",
  "GT: Glycosyltransferases" = "#5E6FB1",
  "PL: Polysaccharide lyases" = "#2FA7A0",
  "CE: Carbohydrate esterases" = "#7B6AAE",
  "CBM: Carbohydrate-binding modules" = "#3F7FB3",
  "AA: Auxiliary activities" = "#A16FB1"
)

# ---------- 4. Plot ----------
# Stacked histogram colored by CAZy class.
# Dashed vertical line indicates zero effect.
# Top-right annotation summarizes the overall directional skew.

p_cazy_effect_hist_class <- cazy_subfamily_effect_class %>%
  ggplot(
    aes(
      x = median_pairwise_log2_diff_pCR_vs_non_pCR,
      fill = CAZy_class
    )
  ) +
  geom_histogram(
    bins = 45,
    color = "white",
    linewidth = 0.2,
    position = "stack"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  annotate(
    "label",
    x = Inf,
    y = Inf,
    label = cazy_subfamily_effect_annotation,
    hjust = 1.02,
    vjust = 1.02,
    size = 3.4,
    label.size = 0.2,
    fill = "white",
    color = "black"
  ) +
  scale_fill_manual(values = cazy_class_colors) +
  labs(
    x = "Median pairwise log2 difference (pCR - non-pCR)",
    y = "Number of CAZyme subfamily features",
    fill = "CAZy class"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid = element_blank()
  )

p_cazy_effect_hist_class

# ggsave(
#   "figures/CAZy_subfamily_effect_size_distribution_by_class.svg",
#   p_cazy_effect_hist_class,
#   width = 7.6,
#   height = 5.0,
#   device = "svg"
# )

# ==================== CAZy class-wise median effect summary plot ====================

# This plot summarizes the feature-level effect-size distribution within each CAZy class.
# The point is the class-wise median effect size.
# Error bars show bootstrap 95% confidence intervals of the class-wise median.

# ---------- 1. Class-wise summary ----------
# Bootstrap is done within each CAZy class across subfamily features.
# Again, this is descriptive because subfamily features are not fully independent.

set.seed(123)

cazy_class_effect_summary <- cazy_subfamily_effect_class %>%
  dplyr::group_by(CAZy_class) %>%
  dplyr::summarise(
    n_features = dplyr::n(),
    n_pCR_higher = sum(median_pairwise_log2_diff_pCR_vs_non_pCR > 0, na.rm = TRUE),
    n_non_pCR_higher = sum(median_pairwise_log2_diff_pCR_vs_non_pCR < 0, na.rm = TRUE),
    proportion_pCR_higher = n_pCR_higher / (n_pCR_higher + n_non_pCR_higher),
    median_effect = median(median_pairwise_log2_diff_pCR_vs_non_pCR, na.rm = TRUE),
    boot_medians = list(
      replicate(
        2000,
        median(
          sample(
            median_pairwise_log2_diff_pCR_vs_non_pCR,
            size = dplyr::n(),
            replace = TRUE
          ),
          na.rm = TRUE
        )
      )
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    ci_low = purrr::map_dbl(boot_medians, ~ quantile(.x, 0.025, na.rm = TRUE)),
    ci_high = purrr::map_dbl(boot_medians, ~ quantile(.x, 0.975, na.rm = TRUE))
  ) %>%
  dplyr::select(-boot_medians) %>%
  dplyr::arrange(median_effect) %>%
  dplyr::mutate(
    CAZy_class = factor(CAZy_class, levels = CAZy_class)
  )

cazy_class_effect_summary

write_tsv(
  cazy_class_effect_summary,
  "260224 final Input file/CAZyme/CAZy_class_effect_summary.tsv"
)

# ---------- 2. Plot ----------
# Horizontal effect-size summary plot with point + bootstrap 95% CI.

# ---------- Fancy class-wise median effect summary plot ----------

effect_xlim <- max(
  abs(c(cazy_class_effect_summary$ci_low, cazy_class_effect_summary$ci_high)),
  na.rm = TRUE
) * 1.25

p_cazy_class_effect_summary <- cazy_class_effect_summary %>%
  ggplot(
    aes(
      y = CAZy_class,
      x = median_effect,
      xmin = ci_low,
      xmax = ci_high,
      color = CAZy_class
    )
  ) +
  annotate(
    "rect",
    xmin = -0.22,
    xmax = 0,
    ymin = -Inf,
    ymax = Inf,
    fill = "#DE7872",
    alpha = 0.045
  ) +
  annotate(
    "rect",
    xmin = 0,
    xmax = 0.22,
    ymin = -Inf,
    ymax = Inf,
    fill = "#4FAE9A",
    alpha = 0.045
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey30"
  ) +
  geom_segment(
    aes(
      x = ci_low,
      xend = ci_high,
      y = CAZy_class,
      yend = CAZy_class
    ),
    linewidth = 1.05,
    lineend = "round"
  ) +
  geom_point(size = 3.6) +
  geom_text(
    aes(
      x = ci_high,
      label = paste0("n = ", n_features)
    ),
    hjust = -0.18,
    size = 3.2,
    color = "black"
  ) +
  annotate(
    "text",
    x = -0.15,
    y = Inf,
    label = "Higher in non-pCR",
    vjust = 1.5,
    size = 3.3,
    color = "grey35"
  ) +
  annotate(
    "text",
    x = 0.15,
    y = Inf,
    label = "Higher in pCR",
    vjust = 1.5,
    size = 3.3,
    color = "grey35"
  ) +
  scale_color_manual(values = cazy_class_colors) +
  scale_y_discrete(
    labels = c(
      "GH: Glycoside hydrolases" = "GH",
      "GT: Glycosyltransferases" = "GT",
      "PL: Polysaccharide lyases" = "PL",
      "CE: Carbohydrate esterases" = "CE",
      "CBM: Carbohydrate-binding modules" = "CBM",
      "AA: Auxiliary activities" = "AA"
    )
  ) +
  scale_x_continuous(
    breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
    labels = c("-0.2", "-0.1", "0", "0.1", "0.2")
  ) +
  labs(
    x = "Median subfamily shift (pCR − non-pCR)",
    y = NULL
  ) +
  coord_cartesian(
    xlim = c(-0.22, 0.22),
    clip = "off"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 10.5, color = "black"),
    axis.title.x = element_text(size = 12, margin = margin(t = 8)),
    legend.position = "none",
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(8, 55, 24, 8)
  )
p_cazy_class_effect_summary

# ggsave(
#   "figures/CAZy_class_median_subfamily_shift_summary.svg",
#   p_cazy_class_effect_summary,
#   width = 4.5,
#   height = 4,
#   device = "svg"
# )



p_cazy_effect_hist_class <- cazy_subfamily_effect_class %>%
  ggplot(
    aes(
      x = median_pairwise_log2_diff_pCR_vs_non_pCR,
      fill = CAZy_class
    )
  ) +
  geom_histogram(
    binwidth = 0.2,
    boundary = 0,
    color = "white",
    linewidth = 0.2,
    position = "stack"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  annotate(
    "label",
    x = Inf,
    y = Inf,
    label = cazy_subfamily_effect_annotation,
    hjust = 1.02,
    vjust = 1.02,
    size = 3.4,
    label.size = 0.2,
    fill = "white",
    color = "black"
  ) +
  scale_fill_manual(values = cazy_class_colors) +
  scale_x_continuous(
    breaks = seq(-2, 2, by = 0.5),
    limits = c(-2, 2)
  ) +
  labs(
    x = "Median pairwise log2 difference (pCR - non-pCR)",
    y = "Number of CAZyme subfamily features",
    fill = "CAZy class"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid = element_blank()
  )

p_cazy_effect_hist_class

p_cazy_effect_hist_class <- cazy_subfamily_effect_class %>%
  ggplot(
    aes(
      x = median_pairwise_log2_diff_pCR_vs_non_pCR,
      fill = CAZy_class
    )
  ) +
  geom_histogram(
    binwidth = 0.2,
    boundary = 0,
    color = "white",
    linewidth = 0.2,
    position = "stack"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  annotate(
    "label",
    x = Inf,
    y = Inf,
    label = cazy_subfamily_effect_annotation,
    hjust = 1.02,
    vjust = 1.02,
    size = 3.4,
    label.size = 0.2,
    fill = "white",
    color = "black"
  ) +
  scale_fill_manual(values = cazy_class_colors) +
  scale_x_continuous(
    breaks = seq(-2, 2, by = 0.5)
  ) +
  coord_cartesian(xlim = c(-2, 2)) +
  labs(
    x = "Median pairwise log2 difference (pCR - non-pCR)",
    y = "Number of CAZyme subfamily features",
    fill = "CAZy class"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid = element_blank()
  )

p_cazy_effect_hist_class

# ggsave(
#   "figures/CAZy_subfamily_effect_size_distribution_by_class.svg",
#   p_cazy_effect_hist_class,
#   width = 7.6,
#   height = 5.0,
#   device = "svg"
# )




# ==================== CAZyme subfamily volcano plot ====================

library(ggrepel)

volcano_p_cutoff <- 0.05
volcano_top_label_n <- 8

cazy_subfamily_volcano <- cazy_subfamily_result %>%
  dplyr::filter(
    !is.na(p_value),
    !is.na(median_pairwise_log2_diff_pCR_vs_non_pCR),
    prevalence >= 0.30,
    mean_TPM_all >= 1
  ) %>%
  dplyr::mutate(
    log2_FC = median_pairwise_log2_diff_pCR_vs_non_pCR,
    neg_log10_p = -log10(pmax(p_value, .Machine$double.xmin)),
    Volcano_group = dplyr::case_when(
      p_value < volcano_p_cutoff & log2_FC > 0 ~ "CR-enriched",
      p_value < volcano_p_cutoff & log2_FC < 0 ~ "nonCR-enriched",
      TRUE ~ "Not significant"
    ),
    Volcano_group = factor(
      Volcano_group,
      levels = c("nonCR-enriched", "Not significant", "CR-enriched")
    ),
    label_score = neg_log10_p * abs(log2_FC)
  )

cazy_subfamily_volcano_label <- cazy_subfamily_volcano %>%
  dplyr::filter(Volcano_group != "Not significant") %>%
  dplyr::group_by(Volcano_group) %>%
  dplyr::arrange(
    dplyr::desc(label_score),
    p_value,
    dplyr::desc(abs(log2_FC)),
    Feature
  ) %>%
  dplyr::slice_head(n = volcano_top_label_n) %>%
  dplyr::ungroup()

cazy_subfamily_volcano %>%
  dplyr::count(Volcano_group)

cazy_subfamily_volcano_label %>%
  dplyr::select(
    Feature,
    Functional_label,
    Volcano_group,
    log2_FC,
    p_value,
    neg_log10_p,
    label_score
  )

write_tsv(
  cazy_subfamily_volcano,
  "260224 final Input file/CAZyme/CAZy_subfamily_volcano_result.tsv"
)

write_tsv(
  cazy_subfamily_volcano_label,
  "260224 final Input file/CAZyme/CAZy_subfamily_volcano_labeled_features.tsv"
)

# ---------- Volcano plot ----------

p_cazy_subfamily_volcano <- ggplot(
  cazy_subfamily_volcano,
  aes(x = log2_FC, y = neg_log10_p)
) +
  geom_point(
    data = cazy_subfamily_volcano %>%
      dplyr::filter(Volcano_group == "Not significant"),
    color = "grey78",
    size = 1.4,
    alpha = 0.55
  ) +
  geom_point(
    data = cazy_subfamily_volcano %>%
      dplyr::filter(Volcano_group != "Not significant"),
    aes(color = Volcano_group),
    size = 2.2,
    alpha = 0.88
  ) +
  geom_hline(
    yintercept = -log10(volcano_p_cutoff),
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.45,
    color = "grey35"
  ) +
  ggrepel::geom_text_repel(
    data = cazy_subfamily_volcano_label,
    aes(label = Feature, color = Volcano_group),
    size = 3.0,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.color = "grey45",
    segment.linewidth = 0.25,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "CR-enriched" = "#4FAE9A",
      "nonCR-enriched" = "#DE7872"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.06, 0.08))
  ) +
  labs(
    x = "Median pairwise log2 difference (CR − nonCR)",
    y = expression(-log[10](italic(p))),
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 10.5),
    legend.key.width = unit(0.8, "lines"),
    panel.grid = element_blank(),
    plot.margin = margin(8, 12, 8, 8)
  )

p_cazy_subfamily_volcano

# ggsave(
#   "figures/CAZy_subfamily_volcano_CR_nonCR.svg",
#   p_cazy_subfamily_volcano,
#   width = 3.8,
#   height = 5.0,
#   device = "svg"
# )


# ==================== Panel F. Representative CAZyme features ====================

library(ggbeeswarm)

group_colors <- c("pCR" = "#4FAE9A", "non_pCR" = "#DE7872")
group_labels <- c("pCR" = "pCR", "non_pCR" = "non-pCR")

dir.create("figures", showWarnings = FALSE)

# ---------- 1. Before sample metadata ----------

m_before <- m %>%
  dplyr::filter(
    TNT == "Before",
    TRG_1 %in% c("CR", "nonCR")
  ) %>%
  dplyr::mutate(
    Response = dplyr::case_when(
      TRG_1 == "CR" ~ "pCR",
      TRG_1 == "nonCR" ~ "non_pCR"
    ),
    Response = factor(Response, levels = c("pCR", "non_pCR"))
  ) %>%
  dplyr::arrange(Response, SampleID)


# ---------- 2. Check object formats ----------

if (!"Feature" %in% colnames(cazy_subfamily)) {
  cazy_subfamily <- cazy_subfamily %>%
    dplyr::rename(Feature = 1)
}

if (!"SampleID" %in% colnames(cazy_module)) {
  if ("sample_id" %in% colnames(cazy_module)) {
    cazy_module <- cazy_module %>%
      dplyr::rename(SampleID = sample_id)
  } else {
    cazy_module <- cazy_module %>%
      dplyr::rename(SampleID = 1)
  }
}


# ---------- 3. Select representative subfamilies ----------
# Selection score = -log10(p) × |effect size|
# pCR side: top 3 significant features
# non-pCR side: top 2 significant features
# If fewer features pass p < 0.05 in one direction, the strongest directional features are used.

cazy_subfamily_selected_pcr <- cazy_subfamily_result %>%
  dplyr::filter(
    !is.na(p_value),
    !is.na(median_pairwise_log2_diff_pCR_vs_non_pCR),
    prevalence >= 0.30,
    mean_TPM_all >= 1,
    p_value < 0.05,
    median_pairwise_log2_diff_pCR_vs_non_pCR > 0
  ) %>%
  dplyr::mutate(
    selection_score = -log10(pmax(p_value, .Machine$double.xmin)) *
      abs(median_pairwise_log2_diff_pCR_vs_non_pCR)
  ) %>%
  dplyr::arrange(
    dplyr::desc(selection_score),
    p_value,
    dplyr::desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR)),
    Feature
  ) %>%
  dplyr::slice_head(n = 3)

if (nrow(cazy_subfamily_selected_pcr) < 3) {
  cazy_subfamily_selected_pcr <- cazy_subfamily_result %>%
    dplyr::filter(
      !is.na(p_value),
      !is.na(median_pairwise_log2_diff_pCR_vs_non_pCR),
      prevalence >= 0.30,
      mean_TPM_all >= 1,
      median_pairwise_log2_diff_pCR_vs_non_pCR > 0
    ) %>%
    dplyr::mutate(
      selection_score = -log10(pmax(p_value, .Machine$double.xmin)) *
        abs(median_pairwise_log2_diff_pCR_vs_non_pCR)
    ) %>%
    dplyr::arrange(
      p_value,
      dplyr::desc(selection_score),
      dplyr::desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR)),
      Feature
    ) %>%
    dplyr::slice_head(n = 3)
}

cazy_subfamily_selected_nonpcr <- cazy_subfamily_result %>%
  dplyr::filter(
    !is.na(p_value),
    !is.na(median_pairwise_log2_diff_pCR_vs_non_pCR),
    prevalence >= 0.30,
    mean_TPM_all >= 1,
    p_value < 0.05,
    median_pairwise_log2_diff_pCR_vs_non_pCR < 0
  ) %>%
  dplyr::mutate(
    selection_score = -log10(pmax(p_value, .Machine$double.xmin)) *
      abs(median_pairwise_log2_diff_pCR_vs_non_pCR)
  ) %>%
  dplyr::arrange(
    dplyr::desc(selection_score),
    p_value,
    dplyr::desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR)),
    Feature
  ) %>%
  dplyr::slice_head(n = 2)

if (nrow(cazy_subfamily_selected_nonpcr) < 2) {
  cazy_subfamily_selected_nonpcr <- cazy_subfamily_result %>%
    dplyr::filter(
      !is.na(p_value),
      !is.na(median_pairwise_log2_diff_pCR_vs_non_pCR),
      prevalence >= 0.30,
      mean_TPM_all >= 1,
      median_pairwise_log2_diff_pCR_vs_non_pCR < 0
    ) %>%
    dplyr::mutate(
      selection_score = -log10(pmax(p_value, .Machine$double.xmin)) *
        abs(median_pairwise_log2_diff_pCR_vs_non_pCR)
    ) %>%
    dplyr::arrange(
      p_value,
      dplyr::desc(selection_score),
      dplyr::desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR)),
      Feature
    ) %>%
    dplyr::slice_head(n = 2)
}

cazy_subfamily_selected_panelF <- dplyr::bind_rows(
  cazy_subfamily_selected_pcr,
  cazy_subfamily_selected_nonpcr
) %>%
  dplyr::distinct(Feature, .keep_all = TRUE) %>%
  dplyr::mutate(
    Direction = dplyr::case_when(
      median_pairwise_log2_diff_pCR_vs_non_pCR > 0 ~ "pCR-enriched",
      median_pairwise_log2_diff_pCR_vs_non_pCR < 0 ~ "non-pCR-enriched"
    ),
    Panel = Feature
  )

cazy_subfamily_selected_panelF %>%
  dplyr::select(
    Feature,
    Direction,
    prevalence,
    mean_TPM_all,
    median_pairwise_log2_diff_pCR_vs_non_pCR,
    p_value,
    selection_score
  )

write_tsv(
  cazy_subfamily_selected_panelF,
  "260224 final Input file/CAZyme/CAZy_subfamily_selected_for_panelF.tsv"
)


# ---------- 4. Module panel ----------
# cazy_module has samples as rows and module scores as columns.

cazy_panelF_module <- cazy_module %>%
  dplyr::select(
    SampleID,
    Value = sucrose_starch_score
  ) %>%
  dplyr::inner_join(
    m_before %>%
      dplyr::select(SampleID, Response),
    by = "SampleID"
  ) %>%
  dplyr::mutate(
    Panel = "Starch/sucrose-associated\nCAZyme capacity",
    Data_type = "Module score"
  ) %>%
  dplyr::select(SampleID, Response, Panel, Data_type, Value)


# ---------- 5. Subfamily panels ----------
# Subfamily abundance is shown as log10(TPM + 1).

cazy_panelF_subfamily <- cazy_subfamily %>%
  dplyr::filter(Feature %in% cazy_subfamily_selected_panelF$Feature) %>%
  dplyr::select(Feature, dplyr::all_of(m_before$SampleID)) %>%
  tidyr::pivot_longer(
    cols = -Feature,
    names_to = "SampleID",
    values_to = "TPM"
  ) %>%
  dplyr::inner_join(
    cazy_subfamily_selected_panelF %>%
      dplyr::select(Feature, Panel),
    by = "Feature"
  ) %>%
  dplyr::inner_join(
    m_before %>%
      dplyr::select(SampleID, Response),
    by = "SampleID"
  ) %>%
  dplyr::mutate(
    Value = log10(TPM + 1),
    Data_type = "Subfamily abundance, log10(TPM + 1)"
  ) %>%
  dplyr::select(SampleID, Response, Panel, Data_type, Value)


# ==================== Panel F. Representative CAZyme features ====================

library(ggbeeswarm)

group_colors <- c("pCR" = "#4FAE9A", "non_pCR" = "#DE7872")
group_labels <- c("pCR" = "pCR", "non_pCR" = "non-pCR")

dir.create("figures", showWarnings = FALSE)

# ---------- 1. Before sample metadata ----------

m_before <- m %>%
  dplyr::filter(
    TNT == "Before",
    TRG_1 %in% c("CR", "nonCR")
  ) %>%
  dplyr::mutate(
    Response = dplyr::case_when(
      TRG_1 == "CR" ~ "pCR",
      TRG_1 == "nonCR" ~ "non_pCR"
    ),
    Response = factor(Response, levels = c("pCR", "non_pCR"))
  ) %>%
  dplyr::arrange(Response, SampleID)


# ---------- 2. Check object formats ----------

if (!"Feature" %in% colnames(cazy_subfamily)) {
  cazy_subfamily <- cazy_subfamily %>%
    dplyr::rename(Feature = 1)
}

if (!"SampleID" %in% colnames(cazy_module)) {
  if ("sample_id" %in% colnames(cazy_module)) {
    cazy_module <- cazy_module %>%
      dplyr::rename(SampleID = sample_id)
  } else {
    cazy_module <- cazy_module %>%
      dplyr::rename(SampleID = 1)
  }
}


# ---------- 3. Select representative subfamilies ----------
# Selection score = -log10(p) × |effect size|
# pCR side: top 3 significant features
# non-pCR side: top 2 significant features
# If fewer features pass p < 0.05 in one direction, the strongest directional features are used.

cazy_subfamily_selected_pcr <- cazy_subfamily_result %>%
  dplyr::filter(
    !is.na(p_value),
    !is.na(median_pairwise_log2_diff_pCR_vs_non_pCR),
    prevalence >= 0.30,
    mean_TPM_all >= 1,
    p_value < 0.05,
    median_pairwise_log2_diff_pCR_vs_non_pCR > 0
  ) %>%
  dplyr::mutate(
    selection_score = -log10(pmax(p_value, .Machine$double.xmin)) *
      abs(median_pairwise_log2_diff_pCR_vs_non_pCR)
  ) %>%
  dplyr::arrange(
    dplyr::desc(selection_score),
    p_value,
    dplyr::desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR)),
    Feature
  ) %>%
  dplyr::slice_head(n = 3)

if (nrow(cazy_subfamily_selected_pcr) < 3) {
  cazy_subfamily_selected_pcr <- cazy_subfamily_result %>%
    dplyr::filter(
      !is.na(p_value),
      !is.na(median_pairwise_log2_diff_pCR_vs_non_pCR),
      prevalence >= 0.30,
      mean_TPM_all >= 1,
      median_pairwise_log2_diff_pCR_vs_non_pCR > 0
    ) %>%
    dplyr::mutate(
      selection_score = -log10(pmax(p_value, .Machine$double.xmin)) *
        abs(median_pairwise_log2_diff_pCR_vs_non_pCR)
    ) %>%
    dplyr::arrange(
      p_value,
      dplyr::desc(selection_score),
      dplyr::desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR)),
      Feature
    ) %>%
    dplyr::slice_head(n = 3)
}

cazy_subfamily_selected_nonpcr <- cazy_subfamily_result %>%
  dplyr::filter(
    !is.na(p_value),
    !is.na(median_pairwise_log2_diff_pCR_vs_non_pCR),
    prevalence >= 0.30,
    mean_TPM_all >= 1,
    p_value < 0.05,
    median_pairwise_log2_diff_pCR_vs_non_pCR < 0
  ) %>%
  dplyr::mutate(
    selection_score = -log10(pmax(p_value, .Machine$double.xmin)) *
      abs(median_pairwise_log2_diff_pCR_vs_non_pCR)
  ) %>%
  dplyr::arrange(
    dplyr::desc(selection_score),
    p_value,
    dplyr::desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR)),
    Feature
  ) %>%
  dplyr::slice_head(n = 2)

if (nrow(cazy_subfamily_selected_nonpcr) < 2) {
  cazy_subfamily_selected_nonpcr <- cazy_subfamily_result %>%
    dplyr::filter(
      !is.na(p_value),
      !is.na(median_pairwise_log2_diff_pCR_vs_non_pCR),
      prevalence >= 0.30,
      mean_TPM_all >= 1,
      median_pairwise_log2_diff_pCR_vs_non_pCR < 0
    ) %>%
    dplyr::mutate(
      selection_score = -log10(pmax(p_value, .Machine$double.xmin)) *
        abs(median_pairwise_log2_diff_pCR_vs_non_pCR)
    ) %>%
    dplyr::arrange(
      p_value,
      dplyr::desc(selection_score),
      dplyr::desc(abs(median_pairwise_log2_diff_pCR_vs_non_pCR)),
      Feature
    ) %>%
    dplyr::slice_head(n = 2)
}

cazy_subfamily_selected_panelF <- dplyr::bind_rows(
  cazy_subfamily_selected_pcr,
  cazy_subfamily_selected_nonpcr
) %>%
  dplyr::distinct(Feature, .keep_all = TRUE) %>%
  dplyr::mutate(
    Direction = dplyr::case_when(
      median_pairwise_log2_diff_pCR_vs_non_pCR > 0 ~ "pCR-enriched",
      median_pairwise_log2_diff_pCR_vs_non_pCR < 0 ~ "non-pCR-enriched"
    ),
    Panel = Feature
  )

cazy_subfamily_selected_panelF %>%
  dplyr::select(
    Feature,
    Direction,
    prevalence,
    mean_TPM_all,
    median_pairwise_log2_diff_pCR_vs_non_pCR,
    p_value,
    selection_score
  )

write_tsv(
  cazy_subfamily_selected_panelF,
  "260224 final Input file/CAZyme/CAZy_subfamily_selected_for_panelF.tsv"
)


# ---------- 4. Module panel ----------
# cazy_module has samples as rows and module scores as columns.

cazy_panelF_module <- cazy_module %>%
  dplyr::select(
    SampleID,
    Value = sucrose_starch_score
  ) %>%
  dplyr::inner_join(
    m_before %>%
      dplyr::select(SampleID, Response),
    by = "SampleID"
  ) %>%
  dplyr::mutate(
    Panel = "Starch/sucrose-associated\nCAZyme capacity",
    Data_type = "Module score"
  ) %>%
  dplyr::select(SampleID, Response, Panel, Data_type, Value)


# ---------- 5. Subfamily panels ----------
# Subfamily abundance is shown as log10(TPM + 1).

cazy_panelF_subfamily <- cazy_subfamily %>%
  dplyr::filter(Feature %in% cazy_subfamily_selected_panelF$Feature) %>%
  dplyr::select(Feature, dplyr::all_of(m_before$SampleID)) %>%
  tidyr::pivot_longer(
    cols = -Feature,
    names_to = "SampleID",
    values_to = "TPM"
  ) %>%
  dplyr::inner_join(
    cazy_subfamily_selected_panelF %>%
      dplyr::select(Feature, Panel),
    by = "Feature"
  ) %>%
  dplyr::inner_join(
    m_before %>%
      dplyr::select(SampleID, Response),
    by = "SampleID"
  ) %>%
  dplyr::mutate(
    Value = log10(TPM + 1),
    Data_type = "Subfamily abundance, log10(TPM + 1)"
  ) %>%
  dplyr::select(SampleID, Response, Panel, Data_type, Value)


# ---------- 6. Combine panels ----------
# Keep feature name and y-axis unit as separate columns.
# Each small panel will have its own y-axis label.

cazy_panelF <- dplyr::bind_rows(
  cazy_panelF_module,
  cazy_panelF_subfamily
) %>%
  dplyr::mutate(
    Response = factor(Response, levels = c("pCR", "non_pCR")),
    Panel = factor(
      Panel,
      levels = c(
        "Starch/sucrose-associated\nCAZyme capacity",
        cazy_subfamily_selected_panelF$Panel
      )
    ),
    Y_label = dplyr::case_when(
      Data_type == "Module score" ~ "Module score",
      TRUE ~ "log10(TPM + 1)"
    )
  )

cazy_panelF %>%
  dplyr::count(Panel, Data_type, Y_label)

cazy_panelF %>%
  dplyr::filter(is.na(Panel))



# ---------- 7. Plot ----------
# Six independent small panels are generated and arranged as 2 rows × 3 columns.
# This allows each panel to carry its own small y-axis unit.

cazy_panelF_split <- split(
  cazy_panelF,
  cazy_panelF$Panel
)

p_cazy_panelF_list <- purrr::map(
  cazy_panelF_split,
  ~ ggplot(
    .x,
    aes(x = Response, y = Value, color = Response)
  ) +
    geom_crossbar(
      stat = "summary",
      fun.min = function(x) quantile(x, 0.25, na.rm = TRUE),
      fun = median,
      fun.max = function(x) quantile(x, 0.75, na.rm = TRUE),
      width = 0.34,
      fill = NA,
      linewidth = 0.62,
      fatten = 1.0
    ) +
    ggbeeswarm::geom_quasirandom(
      width = 0.14,
      size = 1.65,
      alpha = 0.9
    ) +
    ggpubr::stat_compare_means(
      method = "wilcox.test",
      label = "p.format",
      label.y.npc = 0.95,
      size = 2.65
    ) +
    scale_color_manual(values = group_colors) +
    scale_x_discrete(labels = group_labels) +
    labs(
      x = NULL,
      y = unique(.x$Y_label),
      title = unique(as.character(.x$Panel))
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(
        size = 8.8,
        hjust = 0.5,
        face = "plain",
        lineheight = 0.86
      ),
      axis.text.x = element_text(size = 8.0, color = "black"),
      axis.text.y = element_text(size = 6.2, color = "black"),
      axis.title.y = element_text(size = 7.2, color = "black", margin = margin(r = 2)),
      legend.position = "none",
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.36),
      axis.ticks = element_line(color = "black", linewidth = 0.28),
      plot.margin = margin(4, 4, 4, 4)
    )
)

p_cazy_panelF <- ggpubr::ggarrange(
  plotlist = p_cazy_panelF_list,
  nrow = 2,
  ncol = 3,
  align = "hv"
)

p_cazy_panelF

ggsave(
  "figures/F_CAZyme_representative_features.svg",
  p_cazy_panelF,
  width = 4,
  height = 4,
  device = "svg"
)
Sys.setFileTime(
  "figures/F_CAZyme_representative_features.svg",
  Sys.time()
)


# ==================== Panel F with GH32-excluded starch/alpha-glucan score ====================

# The five representative subfamily panels are unchanged. Only the first panel
# is replaced with the starch/alpha-glucan-associated capacity score calculated
# in Section 2-3 from GH13, GH65, GH77, GH97, CBM20, CBM25, CBM26, CBM34,
# and CBM48. GH32 is not included in this score.

cazy_panelF_without_gh32 <- dplyr::bind_rows(
  cazy_capacity_scores %>%
    dplyr::filter(Capacity == "Starch/α-glucan") %>%
    dplyr::transmute(
      SampleID,
      Response,
      Panel = "Starch/α-glucan-associated\nCAZyme capacity",
      Data_type = "GH32-excluded capacity score",
      Value = CAZyme_capacity_score
    ),
  cazy_panelF_subfamily
) %>%
  dplyr::mutate(
    Response = factor(Response, levels = c("pCR", "non_pCR")),
    Panel = factor(
      Panel,
      levels = c(
        "Starch/α-glucan-associated\nCAZyme capacity",
        cazy_subfamily_selected_panelF$Panel
      )
    ),
    Y_label = dplyr::case_when(
      Data_type == "GH32-excluded capacity score" ~ "Capacity score",
      TRUE ~ "log10(TPM + 1)"
    )
  )

cazy_panelF_without_gh32 %>%
  dplyr::count(Panel, Data_type, Y_label)

p_cazy_panelF_without_gh32_list <- cazy_panelF_without_gh32 %>%
  split(.$Panel) %>%
  purrr::map(
    ~ ggplot(
      .x,
      aes(x = Response, y = Value, color = Response)
    ) +
      geom_crossbar(
        stat = "summary",
        fun.min = function(x) quantile(x, 0.25, na.rm = TRUE),
        fun = median,
        fun.max = function(x) quantile(x, 0.75, na.rm = TRUE),
        width = 0.34,
        fill = NA,
        linewidth = 0.62,
        fatten = 1.0
      ) +
      ggbeeswarm::geom_quasirandom(
        width = 0.14,
        size = 1.65,
        alpha = 0.9
      ) +
      ggpubr::stat_compare_means(
        method = "wilcox.test",
        label = "p.format",
        label.y.npc = 0.95,
        size = 2.65
      ) +
      scale_color_manual(values = group_colors) +
      scale_x_discrete(labels = group_labels) +
      labs(
        x = NULL,
        y = unique(.x$Y_label),
        title = unique(as.character(.x$Panel))
      ) +
      theme_classic(base_size = 10) +
      theme(
        plot.title = element_text(
          size = 8.8,
          hjust = 0.5,
          face = "plain",
          lineheight = 0.86
        ),
        axis.text.x = element_text(size = 8.0, color = "black"),
        axis.text.y = element_text(size = 6.2, color = "black"),
        axis.title.y = element_text(
          size = 7.2,
          color = "black",
          margin = margin(r = 2)
        ),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.36),
        axis.ticks = element_line(color = "black", linewidth = 0.28),
        plot.margin = margin(4, 4, 4, 4)
      )
  )

p_cazy_panelF_without_gh32 <- ggpubr::ggarrange(
  plotlist = p_cazy_panelF_without_gh32_list,
  nrow = 2,
  ncol = 3,
  align = "hv"
)

p_cazy_panelF_without_gh32

ggsave(
  "figures/F_CAZyme_representative_features_without_GH32.svg",
  p_cazy_panelF_without_gh32,
  width = 4,
  height = 4,
  device = "svg"
)
Sys.setFileTime(
  "figures/F_CAZyme_representative_features_without_GH32.svg",
  Sys.time()
)
