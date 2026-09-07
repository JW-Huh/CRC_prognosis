#-----------------------------------------------------------------#
#                                                                 #
#     Metabolite figure: forest plot, ternary plot, beeswarm      #
#                                                                 #
#-----------------------------------------------------------------#

rm(list = ls())
options(stringsAsFactors = FALSE)
setwd("D:/2-연구/2-CRC metagenomics/")

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(patchwork)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

set.seed(123)

group_cols <- c(
  pCR = "#4FAE9A",
  non_pCR = "#DE7872",
  CR = "#4FAE9A",
  nonCR = "#DE7872"
)

B <- 2000


#-----------------------------------------------------------------#
# 1. Load and merge metabolite data
#-----------------------------------------------------------------#

metabolite_raw <- read.csv("input/metabolites.csv")

ID_matching <- read.csv("input/metadata_matching.csv")
ID_matching$SNU_ID <- sub(
  ".*?(AL\\d{3}-\\d{3})_S_(\\d{2})$",
  "\\1_\\2",
  ID_matching$SNU_ID
)

metabolite <- metabolite_raw %>%
  dplyr::rename(DNA_ID = SampleID) %>%
  dplyr::left_join(ID_matching, by = "DNA_ID") %>%
  dplyr::select(DNA_ID, SNU_ID, Sample_weight_g, DW_ul, Time, dplyr::everything())

m <- readxl::read_xlsx("260224 final Input file/metadata_최종.xlsx") %>%
  dplyr::mutate(
    SampleID = paste0("Sample_", SampleID),
    TRG = ifelse(is.na(TRG), "nearCR", TRG),
    Pre_Op_Tstage_bin = ifelse(Pre_Op_Tstage >= 3, "1", "0"),
    Pre_Op_Nstage_bin = ifelse(Pre_Op_Nstage >= 1, "1", "0"),
    Age_bin = ifelse(Age >= 60, "1", "0"),
    CEA = ifelse(CEA == "<0.5", 0.5, CEA),
    CEA = as.numeric(CEA),
    CEA_bin = ifelse(CEA < 5.0, "0", "1"),
    BMI_bin = ifelse(BMI > 25.0, "1", "0")
  ) %>%
  dplyr::filter(!is.na(TRG_score))

metabolite <- metabolite %>%
  merge(
    m %>%
      dplyr::mutate(
        SNU_ID = ifelse(
          TNT == "Before",
          paste0(SNU_ID, "_01"),
          paste0(SNU_ID, "_02")
        )
      ),
    by = "SNU_ID"
  )


#-----------------------------------------------------------------#
# 2. Define metabolite columns and baseline table
#-----------------------------------------------------------------#

metabolite_cols <- c(
  "Acetate", "Propionate", "Isobutyrate", "Butyrate",
  "Isovalerate", "Valerate",
  "Acetaminophen", "Aminophenol", "Anthranilic_acid", "Dopamine",
  "gamma_Aminobutyric_acid", "Glutamic_acid", "Histamine",
  "Indole", "Indole_acetamide",
  "Indole_acetic_acid", "Indole_lactic_acid", "Indolepropionic_acid",
  "Kynurenic_acid", "Nicotinic_acid", "Picolinic_acid",
  "Tryptamine", "Tryptophan", "Xanthurenic_acid",
  "Cholic_acid", "Chenodeoxycholic_acid",
  "Glycocholic_acid", "Glycochenodeoxycholic_acid",
  "Taurocholic_acid", "Taurochenodeoxycholic_acid",
  "Deoxycholic_acid", "Hyodeoxycholic_acid",
  "Lithocholi_acid", "Ursodeoxycholic_acid",
  "Glycodeoxycholic_acid", "Glycolithocholic_acid",
  "Glycoursodeoxycholic_acid",
  "Taurodeoxycholic_acid", "Taurolithocholic_acid",
  "Tauroursodeoxycholic_acid_Taurohyodeoxycholic_acid"
)

metabolite_before <- metabolite %>%
  dplyr::filter(TNT == "Before", TRG_1 %in% c("CR", "nonCR")) %>%
  dplyr::select(TRG_1, SampleID, dplyr::all_of(metabolite_cols))

prevalence <- colSums(metabolite_before[, -c(1:2)] > 0, na.rm = TRUE)
valid_metabolite_cols <- names(prevalence[prevalence >= ceiling(0.20 * nrow(metabolite_before))])

metabolite_before_valid <- metabolite_before %>%
  dplyr::select(TRG_1, SampleID, dplyr::all_of(valid_metabolite_cols))


#-----------------------------------------------------------------#
# 3. Forest plot: Hedges' g for selected metabolites
#-----------------------------------------------------------------#

screen_tbl <- data.frame()

for (met in valid_metabolite_cols) {
  
  x <- as.numeric(metabolite_before_valid[metabolite_before_valid$TRG_1 == "CR", met])
  y <- as.numeric(metabolite_before_valid[metabolite_before_valid$TRG_1 == "nonCR", met])
  
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  p <- tryCatch(
    stats::wilcox.test(x, y, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  
  screen_tbl <- rbind(
    screen_tbl,
    data.frame(
      metabolite = met,
      n_pCR = length(x),
      n_non_pCR = length(y),
      wilcox_p = p,
      median_pCR = stats::median(x, na.rm = TRUE),
      median_non_pCR = stats::median(y, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  )
}

selected_mets <- screen_tbl %>%
  dplyr::filter(is.finite(wilcox_p), wilcox_p < 0.10) %>%
  dplyr::arrange(wilcox_p) %>%
  dplyr::slice_head(n = 7) %>%
  dplyr::pull(metabolite)

effect_tbl <- data.frame()

for (met in selected_mets) {
  
  x <- as.numeric(metabolite_before_valid[metabolite_before_valid$TRG_1 == "CR", met])
  y <- as.numeric(metabolite_before_valid[metabolite_before_valid$TRG_1 == "nonCR", met])
  
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  p <- tryCatch(
    stats::wilcox.test(x, y, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  
  n1 <- length(x)
  n2 <- length(y)
  df <- n1 + n2 - 2
  
  pooled_sd <- sqrt(
    ((n1 - 1) * stats::var(x, na.rm = TRUE) +
       (n2 - 1) * stats::var(y, na.rm = TRUE)) / df
  )
  
  hedges_g <- ifelse(
    is.finite(pooled_sd) && pooled_sd > 0,
    ((mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)) / pooled_sd) *
      (1 - 3 / (4 * df - 1)),
    NA_real_
  )
  
  boot_g <- replicate(B, {
    
    xb <- sample(x, size = n1, replace = TRUE)
    yb <- sample(y, size = n2, replace = TRUE)
    
    boot_df <- n1 + n2 - 2
    
    boot_sd <- sqrt(
      ((n1 - 1) * stats::var(xb, na.rm = TRUE) +
         (n2 - 1) * stats::var(yb, na.rm = TRUE)) / boot_df
    )
    
    ifelse(
      is.finite(boot_sd) && boot_sd > 0,
      ((mean(xb, na.rm = TRUE) - mean(yb, na.rm = TRUE)) / boot_sd) *
        (1 - 3 / (4 * boot_df - 1)),
      NA_real_
    )
  })
  
  effect_tbl <- rbind(
    effect_tbl,
    data.frame(
      metabolite = met,
      metabolite_label = gsub("_", " ", met),
      n_pCR = n1,
      n_non_pCR = n2,
      median_pCR = stats::median(x, na.rm = TRUE),
      median_non_pCR = stats::median(y, na.rm = TRUE),
      mean_pCR = mean(x, na.rm = TRUE),
      mean_non_pCR = mean(y, na.rm = TRUE),
      wilcox_p = p,
      hedges_g = hedges_g,
      hedges_g_low = unname(stats::quantile(boot_g, 0.025, na.rm = TRUE)),
      hedges_g_high = unname(stats::quantile(boot_g, 0.975, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  )
}

effect_tbl$metabolite_label <- gsub(
  "Lithocholi acid",
  "Lithocholic acid",
  effect_tbl$metabolite_label
)

effect_tbl$direction <- ifelse(
  effect_tbl$hedges_g >= 0,
  "Higher in pCR",
  "Higher in non_pCR"
)


effect_tbl$p_lab <- ifelse(
  is.na(effect_tbl$wilcox_p),
  "NA",
  ifelse(
    effect_tbl$wilcox_p < 0.001,
    "<0.001",
    formatC(effect_tbl$wilcox_p, format = "f", digits = 3)
  )
)

effect_tbl$pCR_lab <- format(signif(effect_tbl$median_pCR, 3), trim = TRUE)
effect_tbl$non_pCR_lab <- format(signif(effect_tbl$median_non_pCR, 3), trim = TRUE)

effect_tbl$metabolite_label <- factor(
  effect_tbl$metabolite_label,
  levels = effect_tbl$metabolite_label[order(effect_tbl$hedges_g)]
)

x_min <- min(effect_tbl$hedges_g_low, na.rm = TRUE)
x_max <- max(effect_tbl$hedges_g_high, na.rm = TRUE)
x_span <- x_max - x_min

x_pCR <- x_min - 0.98 * x_span
x_non_pCR <- x_min - 0.58 * x_span
x_pval <- x_max + 0.33 * x_span
y_top <- nrow(effect_tbl) + 0.75

p_forest <-
  ggplot2::ggplot(effect_tbl, ggplot2::aes(x = hedges_g, y = metabolite_label)) +
  
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = 2,
    color = "grey55",
    linewidth = 0.5
  ) +
  
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin = hedges_g_low,
      xmax = hedges_g_high,
      color = direction
    ),
    height = 0.14,
    linewidth = 0.72
  ) +
  
  ggplot2::geom_point(
    ggplot2::aes(color = direction),
    size = 2.7
  ) +
  
  ggplot2::geom_text(
    ggplot2::aes(x = x_pCR, label = pCR_lab),
    hjust = 1,
    size = 3.0,
    color = group_cols["pCR"]
  ) +
  
  ggplot2::geom_text(
    ggplot2::aes(x = x_non_pCR, label = non_pCR_lab),
    hjust = 1,
    size = 3.0,
    color = group_cols["non_pCR"]
  ) +
  
  ggplot2::geom_text(
    ggplot2::aes(x = x_pval, label = p_lab),
    hjust = 0,
    size = 3.0,
    color = "black"
  ) +
  
  ggplot2::annotate(
    "text",
    x = x_pCR,
    y = y_top,
    label = "pCR\nmedian",
    hjust = 1,
    vjust = 0,
    size = 3.8,
    fontface = "bold",
    color = group_cols[["pCR"]]
  ) +
  
  ggplot2::annotate(
    "text",
    x = x_non_pCR,
    y = y_top,
    label = "non-pCR\nmedian",
    hjust = 1,
    vjust = 0,
    size = 3.8,
    fontface = "bold",
    color = group_cols[["non_pCR"]]
  ) +
  
  ggplot2::annotate(
    "text",
    x = 0,
    y = y_top,
    label = "Hedges' g (95% CI)",
    hjust = 0.5,
    vjust = 0,
    size = 4.3,
    fontface = "bold"
  ) +
  
  ggplot2::annotate(
    "text",
    x = x_pval,
    y = y_top,
    label = "Wilcoxon P",
    hjust = 0,
    vjust = 0,
    size = 3.8,
    fontface = "bold"
  ) +
  
  ggplot2::scale_color_manual(
    values = c(
      "Higher in pCR" = group_cols[["pCR"]],
      "Higher in non_pCR" = group_cols[["non_pCR"]]
    )
  ) +
  
  ggplot2::scale_y_discrete(
    expand = ggplot2::expansion(mult = c(0.02, 0.14))
  ) +
  
  ggplot2::labs(
    x = "Standardized effect size",
    y = NULL
  ) +
  
  ggplot2::coord_cartesian(
    xlim = c(x_min - 1.15 * x_span, x_max + 0.58 * x_span),
    clip = "off"
  ) +
  
  ggplot2::theme_classic(base_size = 11.5) +
  ggplot2::theme(
    legend.position = "none",
    axis.title.x = ggplot2::element_text(size = 12.5),
    axis.text.x = ggplot2::element_text(size = 8.8, color = "black"),
    axis.text.y = ggplot2::element_text(size = 10.2, color = "black", lineheight = 0.92),
    axis.line.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(14, 55, 10, 70)
  )

p_forest

# ggplot2::ggsave(
#   "figures/metabolite_forest_hedges_g_selected.svg",
#   p_forest,
#   width = 7.4,
#   height = 3.5,
#   device = "svg"
# )


#-----------------------------------------------------------------#
# 4. Ternary plot: ILA, IPA, IAA with four boxplot insets
#-----------------------------------------------------------------#

ternary_df <- metabolite_before_valid %>%
  dplyr::select(
    TRG_1,
    SampleID,
    ILA = Indole_lactic_acid,
    IPA = Indolepropionic_acid,
    IAA = Indole_acetic_acid
  ) %>%
  dplyr::filter(TRG_1 %in% c("CR", "nonCR")) %>%
  dplyr::filter(is.finite(ILA), is.finite(IPA), is.finite(IAA)) %>%
  dplyr::mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR")))

z_mat <- scale(as.matrix(ternary_df[, c("ILA", "IPA", "IAA")]))
colnames(z_mat) <- c("ILA_z", "IPA_z", "IAA_z")

z_exp <- exp(z_mat - apply(z_mat, 1, max, na.rm = TRUE))
w_mat <- z_exp / rowSums(z_exp, na.rm = TRUE)
colnames(w_mat) <- c("ILA_w", "IPA_w", "IAA_w")

ternary_df <- cbind(
  ternary_df,
  as.data.frame(z_mat),
  as.data.frame(w_mat)
)

ternary_df$x <- 0.5 * ternary_df$ILA_w + ternary_df$IAA_w
ternary_df$y <- sqrt(3) / 2 * ternary_df$ILA_w

triangle_df <- data.frame(
  x = c(0, 0.5, 1, 0),
  y = c(0, sqrt(3) / 2, 0, 0)
)

grid_df <- data.frame()

for (a in seq(0.2, 0.8, by = 0.2)) {
  
  grid_df <- rbind(
    grid_df,
    data.frame(
      line_id = paste0("ILA_", a),
      ILA_w = c(a, a),
      IPA_w = c(1 - a, 0),
      IAA_w = c(0, 1 - a)
    ),
    data.frame(
      line_id = paste0("IPA_", a),
      ILA_w = c(1 - a, 0),
      IPA_w = c(a, a),
      IAA_w = c(0, 1 - a)
    ),
    data.frame(
      line_id = paste0("IAA_", a),
      ILA_w = c(1 - a, 0),
      IPA_w = c(0, 1 - a),
      IAA_w = c(a, a)
    )
  )
}

grid_df$x <- 0.5 * grid_df$ILA_w + grid_df$IAA_w
grid_df$y <- sqrt(3) / 2 * grid_df$ILA_w

hull_df <- data.frame()

for (grp in c("CR", "nonCR")) {
  
  tmp <- ternary_df[ternary_df$TRG_1 == grp, , drop = FALSE]
  
  if (nrow(tmp) >= 3 && length(unique(tmp$x)) >= 2 && length(unique(tmp$y)) >= 2) {
    hull_df <- rbind(hull_df, tmp[chull(tmp$x, tmp$y), , drop = FALSE])
  }
}

centroid_df <- ternary_df %>%
  dplyr::group_by(TRG_1) %>%
  dplyr::summarise(
    x = mean(x, na.rm = TRUE),
    y = mean(y, na.rm = TRUE),
    .groups = "drop"
  )

pseudo <- min(
  c(
    ternary_df$IPA[ternary_df$IPA > 0],
    ternary_df$ILA[ternary_df$ILA > 0],
    ternary_df$IAA[ternary_df$IAA > 0]
  ),
  na.rm = TRUE
) / 2


#-----------------------------------------------------------------#
# Ternary plot inset: boxplot style update
#-----------------------------------------------------------------#

inset_df <- ternary_df %>%
  dplyr::mutate(
    TRG_plot = factor(
      TRG_1,
      levels = c("CR", "nonCR"),
      labels = c("pCR", "non-pCR")
    ),
    Ratio = log((IPA + ILA + pseudo) / (IAA + pseudo))
  )

p_inset_IPA <-
  ggplot2::ggplot(
    inset_df,
    ggplot2::aes(x = TRG_plot, y = IPA, fill = TRG_plot, color = TRG_plot)
  ) +
  ggplot2::geom_boxplot(
    width = 0.52,
    outlier.shape = NA,
    linewidth = 0.75,
    alpha = 0.18
  ) +
  ggpubr::stat_compare_means(
    comparisons = list(c("pCR", "non-pCR")),
    method = "wilcox.test",
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.8,
    size = 7.2
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "pCR" = scales::alpha(group_cols[["pCR"]], 0.18),
      "non-pCR" = scales::alpha(group_cols[["non_pCR"]], 0.18)
    )
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "pCR" = group_cols[["pCR"]],
      "non-pCR" = group_cols[["non_pCR"]]
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.10, 0.38))
  ) +
  ggplot2::labs(title = "IPA", x = NULL, y = NULL) +
  ggplot2::theme_classic(base_size = 7) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
    legend.position = "none",
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.8),
    plot.margin = ggplot2::margin(2, 2, 2, 2)
  )

p_inset_ILA <-
  ggplot2::ggplot(
    inset_df,
    ggplot2::aes(x = TRG_plot, y = ILA, fill = TRG_plot, color = TRG_plot)
  ) +
  ggplot2::geom_boxplot(
    width = 0.52,
    outlier.shape = NA,
    linewidth = 0.75,
    alpha = 0.18
  ) +
  ggpubr::stat_compare_means(
    comparisons = list(c("pCR", "non-pCR")),
    method = "wilcox.test",
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.8,
    size = 7.2
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "pCR" = scales::alpha(group_cols[["pCR"]], 0.18),
      "non-pCR" = scales::alpha(group_cols[["non_pCR"]], 0.18)
    )
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "pCR" = group_cols[["pCR"]],
      "non-pCR" = group_cols[["non_pCR"]]
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.10, 0.38))
  ) +
  ggplot2::labs(title = "ILA", x = NULL, y = NULL) +
  ggplot2::theme_classic(base_size = 7) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
    legend.position = "none",
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.8),
    plot.margin = ggplot2::margin(2, 2, 2, 2)
  )

p_inset_IAA <-
  ggplot2::ggplot(
    inset_df,
    ggplot2::aes(x = TRG_plot, y = IAA, fill = TRG_plot, color = TRG_plot)
  ) +
  ggplot2::geom_boxplot(
    width = 0.52,
    outlier.shape = NA,
    linewidth = 0.75,
    alpha = 0.18
  ) +
  ggpubr::stat_compare_means(
    comparisons = list(c("pCR", "non-pCR")),
    method = "wilcox.test",
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.8,
    size = 7.2
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "pCR" = scales::alpha(group_cols[["pCR"]], 0.18),
      "non-pCR" = scales::alpha(group_cols[["non_pCR"]], 0.18)
    )
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "pCR" = group_cols[["pCR"]],
      "non-pCR" = group_cols[["non_pCR"]]
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.10, 0.38))
  ) +
  ggplot2::labs(title = "IAA", x = NULL, y = NULL) +
  ggplot2::theme_classic(base_size = 7) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
    legend.position = "none",
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.8),
    plot.margin = ggplot2::margin(2, 2, 2, 2)
  )

p_inset_ratio <-
  ggplot2::ggplot(
    inset_df,
    ggplot2::aes(x = TRG_plot, y = Ratio, fill = TRG_plot, color = TRG_plot)
  ) +
  ggplot2::geom_boxplot(
    width = 0.52,
    outlier.shape = NA,
    linewidth = 0.75,
    alpha = 0.18
  ) +
  ggpubr::stat_compare_means(
    comparisons = list(c("pCR", "non-pCR")),
    method = "wilcox.test",
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.8,
    size = 7.2
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "pCR" = scales::alpha(group_cols[["pCR"]], 0.18),
      "non-pCR" = scales::alpha(group_cols[["non_pCR"]], 0.18)
    )
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "pCR" = group_cols[["pCR"]],
      "non-pCR" = group_cols[["non_pCR"]]
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.10, 0.38))
  ) +
  ggplot2::labs(title = "Log((IPA+ILA)/IAA)", x = NULL, y = NULL) +
  ggplot2::theme_classic(base_size = 7) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 8.5, face = "bold", hjust = 0.5),
    legend.position = "none",
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.8),
    plot.margin = ggplot2::margin(2, 2, 2, 2)
  )

p_ternary_inset <-
  (p_inset_IPA | p_inset_ILA) /
  (p_inset_IAA | p_inset_ratio) +
  patchwork::plot_layout(widths = c(1, 1), heights = c(1, 1))




#-----------------------------------------------------------------#
# Ternary plot with right-side inset
#-----------------------------------------------------------------#

p_ternary <-
  ggplot2::ggplot() +
  
  ggplot2::geom_path(
    data = grid_df,
    ggplot2::aes(x = x, y = y, group = line_id),
    color = "grey80",
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  
  ggplot2::geom_path(
    data = triangle_df,
    ggplot2::aes(x = x, y = y),
    color = "black",
    linewidth = 1.15,
    lineend = "round"
  ) +
  
  ggplot2::geom_polygon(
    data = hull_df,
    ggplot2::aes(x = x, y = y, group = TRG_1, fill = TRG_1, color = TRG_1),
    alpha = 0.11,
    linewidth = 0.8
  ) +
  
  ggplot2::geom_point(
    data = ternary_df,
    ggplot2::aes(x = x, y = y, fill = TRG_1),
    shape = 21,
    color = "white",
    stroke = 0.7,
    size = 3.4,
    alpha = 0.98
  ) +
  
  ggplot2::geom_point(
    data = centroid_df,
    ggplot2::aes(x = x, y = y, fill = TRG_1),
    shape = 23,
    color = "black",
    stroke = 0.9,
    size = 5.2
  ) +
  
  ggplot2::annotate(
    "text",
    x = 0.5,
    y = sqrt(3) / 2 + 0.055,
    label = "ILA",
    size = 11,
    fontface = "bold"
  ) +
  
  ggplot2::annotate(
    "text",
    x = -0.035,
    y = -0.03,
    label = "IPA",
    size = 11,
    fontface = "bold",
    hjust = 1
  ) +
  
  ggplot2::annotate(
    "text",
    x = 1.035,
    y = -0.03,
    label = "IAA",
    size = 11,
    fontface = "bold",
    hjust = 0
  ) +
  
  ggplot2::annotate(
    "text",
    x = 0.21,
    y = 0.63,
    label = "CR-enriched\nILA/IPA-biased",
    color = group_cols[["CR"]],
    size = 5.2,
    lineheight = 0.95
  ) +
  
  ggplot2::annotate(
    "text",
    x = 0.86,
    y = 0.21,
    label = "nonCR-enriched\nIAA-biased",
    color = group_cols[["nonCR"]],
    size = 5.2,
    lineheight = 0.95
  ) +
  
  ggplot2::scale_fill_manual(
    values = c(
      "CR" = group_cols[["CR"]],
      "nonCR" = group_cols[["nonCR"]]
    ),
    breaks = c("CR", "nonCR"),
    labels = c("pCR", "non-pCR")
  ) +
  
  ggplot2::scale_color_manual(
    values = c(
      "CR" = group_cols[["CR"]],
      "nonCR" = group_cols[["nonCR"]]
    ),
    breaks = c("CR", "nonCR")
  ) +
  
  ggplot2::guides(
    color = "none",
    fill = ggplot2::guide_legend(
      title = NULL,
      override.aes = list(
        shape = 21,
        color = "white",
        size = 6
      )
    )
  ) +
  
  ggplot2::coord_equal(
    xlim = c(-0.08, 1.45),
    ylim = c(-0.16, 0.96),
    clip = "off"
  ) +
  
  ggplot2::theme_void(base_size = 12) +
  ggplot2::theme(
    legend.position = c(0.045, 0.95),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.text = ggplot2::element_text(size = 16),
    plot.margin = ggplot2::margin(18, 28, 32, 28)
  )

p_ternary_with_inset <-
  p_ternary +
  patchwork::inset_element(
    p_ternary_inset,
    left = 0.63,
    bottom = 0.18,
    right = 0.98,
    top = 0.72,
    align_to = "plot"
  )

p_ternary_with_inset

ggplot2::ggsave(
  "figures/ternary_ILA_IPA_IAA_before_CR_nonCR_with_2x2_inset.svg",
  plot = p_ternary_with_inset,
  width = 7.4,
  height = 6.4,
  device = "svg",
  bg = "white"
)



#-----------------------------------------------------------------#
# 5. Selected feature beeswarm plots with median line
#-----------------------------------------------------------------#

feature_plot_df <- metabolite_before_valid %>%
  dplyr::select(TRG_1, SampleID, dplyr::all_of(selected_mets)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(selected_mets),
    names_to = "Metabolite",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    TRG_plot = dplyr::recode(TRG_1, "CR" = "pCR", "nonCR" = "non_pCR"),
    TRG_plot = factor(TRG_plot, levels = c("pCR", "non_pCR")),
    MetaboliteLabel = gsub("_", " ", Metabolite),
    MetaboliteLabel = gsub("Lithocholi acid", "Lithocholic acid", MetaboliteLabel)
  )

feature_plot_df$MetaboliteLabel <- factor(
  feature_plot_df$MetaboliteLabel,
  levels = effect_tbl$metabolite_label[order(effect_tbl$wilcox_p)]
)

p_feature_beeswarm <-
  ggplot2::ggplot(
    feature_plot_df,
    ggplot2::aes(x = TRG_plot, y = Value, color = TRG_plot)
  ) +
  
  ggbeeswarm::geom_beeswarm(
    cex = 2.1,
    alpha = 0.88,
    priority = "density"
  ) +
  
  ggplot2::stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.45,
    linewidth = 0.45,
    color = "black"
  ) +
  
  ggplot2::facet_wrap(
    ~ MetaboliteLabel,
    scales = "free_y",
    nrow = 1
  ) +
  
  ggplot2::scale_color_manual(
    values = c("pCR" = group_cols["pCR"], "non_pCR" = group_cols["non_pCR"])
  ) +
  
  ggplot2::scale_x_discrete(
    labels = c("pCR" = "pCR", "non_pCR" = "non-pCR")
  ) +
  
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.08, 0.18))
  ) +
  
  ggplot2::labs(x = NULL, y = NULL) +
  
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 9.5, face = "bold"),
    axis.text.x = ggplot2::element_text(size = 9.5, color = "black", face = "bold"),
    axis.text.y = ggplot2::element_text(size = 7, color = "black"),
    axis.ticks.y = ggplot2::element_line(linewidth = 0.25, color = "black"),
    axis.ticks.length.y = grid::unit(1.5, "pt"),
    legend.position = "none",
    aspect.ratio = 1.45
  )

p_feature_beeswarm

ggplot2::ggsave(
  "figures/metabolite_selected_features_beeswarm_median.svg",
  plot = p_feature_beeswarm,
  width = 7.8,
  height = 2.4,
  device = "svg"
)


#-----------------------------------------------------------------#
# 6. Store final objects
#-----------------------------------------------------------------#

metabolite_figure_results <- list(
  metabolite = metabolite,
  metabolite_before_valid = metabolite_before_valid,
  screen_tbl = screen_tbl,
  effect_tbl = effect_tbl,
  ternary_df = ternary_df,
  inset_df = inset_df,
  feature_plot_df = feature_plot_df,
  p_forest = p_forest,
  p_ternary_with_inset = p_ternary_with_inset,
  p_feature_beeswarm = p_feature_beeswarm
)

rm(
  list = setdiff(
    ls(),
    c(
      "metabolite_figure_results",
      "p_forest",
      "p_ternary_with_inset",
      "p_feature_beeswarm"
    )
  )
)

gc()