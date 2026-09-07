
########## Perturbation analysis

setwd("D:/2-연구/2-CRC metagenomics/")
options(java.parameters = "-Xmx64g", stringsAsFactors = F)
set.seed(123)

library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggforce)
library(purrr)
library(broom)
library(ggrepel)



ko_selected <- readRDS(file = "input/260619 KO_category_v3.RDS")

colnames(ko_selected)
# [1] "KO"                "Wilcoxon_p"        "CR_mean"           "nonCR_mean"       
# [5] "Mean_abun"         "Log2FC"            "name"              "upper_category"   
# [9] "assigned_category" "confidence"        "brite_category"    "primary_text"     
# [13] "search_text"       "-log10(p)"         "index" 

head(ko_selected)
table(ko_selected$upper_category)
table(ko_selected$assigned_category)

#-----------------------------------------------------------------#
#                                                                 #
#     KO category perturbation analysis using Higher Criticism     #
#                                                                 #
#-----------------------------------------------------------------#

## Input object:
##   ko_selected
##
## Required columns:
##   KO
##   Wilcoxon_p
##   Log2FC
##   assigned_category
##   upper_category
##
## Concept:
##   For each KO category, test whether its p-value distribution is shifted
##   toward smaller p-values compared with random KO sets of the same size.
##
## This is a feature-level competitive enrichment analysis.
## It does not re-permute sample labels. Instead, it preserves the observed
## KO-level p-value background and tests category-level perturbation.

set.seed(123)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

category_col <- "assigned_category"
p_col <- "Wilcoxon_p"
direction_col <- "Log2FC"

n_perm <- 10000
min_set_size <- 10

threshold_grid <- c(
  0.001, 0.003, 0.005,
  0.01, 0.03, 0.05,
  0.10, 0.20, 0.50
)

max_curve_categories <- 8


#-----------------------------------------------------------------#
# 1. Prepare input
#-----------------------------------------------------------------#

ko_perturb <- as.data.frame(ko_selected, stringsAsFactors = FALSE)

ko_perturb <- ko_perturb[
  is.finite(ko_perturb[[p_col]]) &
    !is.na(ko_perturb[[category_col]]) &
    ko_perturb[[category_col]] != "",
  ,
  drop = FALSE
]

ko_perturb$p_value <- as.numeric(ko_perturb[[p_col]])

ko_perturb$p_value <- pmin(
  pmax(ko_perturb$p_value, 1e-300),
  1
)

ko_perturb$log_p <- -log10(ko_perturb$p_value)

ko_perturb$category <- as.character(ko_perturb[[category_col]])

if (direction_col %in% colnames(ko_perturb)) {
  ko_perturb$direction_value <- as.numeric(ko_perturb[[direction_col]])
} else {
  ko_perturb$direction_value <- NA_real_
}

p_all <- ko_perturb$p_value
logp_all <- ko_perturb$log_p

categories <- sort(unique(ko_perturb$category))

ko_category_perturbation <- data.frame()
ko_category_curve <- data.frame()


#-----------------------------------------------------------------#
# 2. Category-wise Higher Criticism perturbation test
#-----------------------------------------------------------------#

for (cat_name in categories) {
  
  cat_idx <- ko_perturb$category == cat_name
  p_obs <- sort(ko_perturb$p_value[cat_idx])
  n_cat <- length(p_obs)
  
  if (n_cat < min_set_size) {
    next
  }
  
  rank_obs <- seq_len(n_cat) / n_cat
  
  ## Higher Criticism statistic.
  ## This captures the maximal standardized excess of small p-values.
  hc_valid <- p_obs > 0 & p_obs < 0.5 & rank_obs > p_obs
  
  if (any(hc_valid)) {
    hc_obs <- max(
      sqrt(n_cat) *
        (rank_obs[hc_valid] - p_obs[hc_valid]) /
        sqrt(p_obs[hc_valid] * (1 - p_obs[hc_valid])),
      na.rm = TRUE
    )
  } else {
    hc_obs <- 0
  }
  
  ## One-sided KS-like statistic.
  ## This captures the maximal cumulative left-shift of p-values.
  ks_obs <- max(rank_obs - p_obs, na.rm = TRUE)
  
  ## Mean -log10(P), useful as an intuitive perturbation burden.
  mean_logp_obs <- mean(-log10(p_obs), na.rm = TRUE)
  
  prop_005_obs <- mean(p_obs <= 0.05, na.rm = TRUE)
  prop_020_obs <- mean(p_obs <= 0.20, na.rm = TRUE)
  
  median_log2fc_obs <- stats::median(
    ko_perturb$direction_value[cat_idx],
    na.rm = TRUE
  )
  
  mean_abs_log2fc_obs <- mean(
    abs(ko_perturb$direction_value[cat_idx]),
    na.rm = TRUE
  )
  
  null_hc <- numeric(n_perm)
  null_ks <- numeric(n_perm)
  null_mean_logp <- numeric(n_perm)
  null_curve <- matrix(
    NA_real_,
    nrow = n_perm,
    ncol = length(threshold_grid)
  )
  
  for (b in seq_len(n_perm)) {
    
    p_null <- sort(sample(p_all, n_cat, replace = FALSE))
    rank_null <- seq_len(n_cat) / n_cat
    
    hc_null_valid <- p_null > 0 & p_null < 0.5 & rank_null > p_null
    
    if (any(hc_null_valid)) {
      null_hc[b] <- max(
        sqrt(n_cat) *
          (rank_null[hc_null_valid] - p_null[hc_null_valid]) /
          sqrt(p_null[hc_null_valid] * (1 - p_null[hc_null_valid])),
        na.rm = TRUE
      )
    } else {
      null_hc[b] <- 0
    }
    
    null_ks[b] <- max(rank_null - p_null, na.rm = TRUE)
    null_mean_logp[b] <- mean(-log10(p_null), na.rm = TRUE)
    
    null_curve[b, ] <- findInterval(threshold_grid, p_null) / n_cat
  }
  
  hc_perm_p <- (sum(null_hc >= hc_obs, na.rm = TRUE) + 1) / (n_perm + 1)
  ks_perm_p <- (sum(null_ks >= ks_obs, na.rm = TRUE) + 1) / (n_perm + 1)
  mean_logp_perm_p <- (
    sum(null_mean_logp >= mean_logp_obs, na.rm = TRUE) + 1
  ) / (n_perm + 1)
  
  hc_z <- ifelse(
    stats::sd(null_hc, na.rm = TRUE) > 0,
    (hc_obs - mean(null_hc, na.rm = TRUE)) / stats::sd(null_hc, na.rm = TRUE),
    NA_real_
  )
  
  ks_z <- ifelse(
    stats::sd(null_ks, na.rm = TRUE) > 0,
    (ks_obs - mean(null_ks, na.rm = TRUE)) / stats::sd(null_ks, na.rm = TRUE),
    NA_real_
  )
  
  mean_logp_z <- ifelse(
    stats::sd(null_mean_logp, na.rm = TRUE) > 0,
    (
      mean_logp_obs - mean(null_mean_logp, na.rm = TRUE)
    ) / stats::sd(null_mean_logp, na.rm = TRUE),
    NA_real_
  )
  
  ko_category_perturbation <- rbind(
    ko_category_perturbation,
    data.frame(
      category = cat_name,
      n_KO = n_cat,
      
      HC_score = hc_obs,
      HC_null_mean = mean(null_hc, na.rm = TRUE),
      HC_z = hc_z,
      HC_perm_p = hc_perm_p,
      
      KS_score = ks_obs,
      KS_null_mean = mean(null_ks, na.rm = TRUE),
      KS_z = ks_z,
      KS_perm_p = ks_perm_p,
      
      mean_log10P = mean_logp_obs,
      mean_log10P_null_mean = mean(null_mean_logp, na.rm = TRUE),
      mean_log10P_z = mean_logp_z,
      mean_log10P_perm_p = mean_logp_perm_p,
      
      prop_P005 = prop_005_obs,
      prop_P020 = prop_020_obs,
      
      median_Log2FC = median_log2fc_obs,
      mean_abs_Log2FC = mean_abs_log2fc_obs,
      
      stringsAsFactors = FALSE
    )
  )
  
  ko_category_curve <- rbind(
    ko_category_curve,
    data.frame(
      category = cat_name,
      n_KO = n_cat,
      threshold = threshold_grid,
      neg_log10_threshold = -log10(threshold_grid),
      observed_fraction = findInterval(threshold_grid, p_obs) / n_cat,
      null_mean = apply(null_curve, 2, mean, na.rm = TRUE),
      null_low = apply(null_curve, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
      null_high = apply(null_curve, 2, stats::quantile, probs = 0.975, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  )
}

ko_category_perturbation$HC_q <- stats::p.adjust(
  ko_category_perturbation$HC_perm_p,
  method = "BH"
)

ko_category_perturbation$KS_q <- stats::p.adjust(
  ko_category_perturbation$KS_perm_p,
  method = "BH"
)

ko_category_perturbation$mean_log10P_q <- stats::p.adjust(
  ko_category_perturbation$mean_log10P_perm_p,
  method = "BH"
)

ko_category_perturbation <- ko_category_perturbation[
  order(
    ko_category_perturbation$HC_perm_p,
    -ko_category_perturbation$HC_z,
    -ko_category_perturbation$mean_log10P_z
  ),
  ,
  drop = FALSE
]

ko_category_perturbation$category_label <- sapply(
  ko_category_perturbation$category,
  function(x) {
    paste(strwrap(x, width = 34), collapse = "\n")
  }
)

ko_category_perturbation$minus_log10_HC_p <- -log10(
  pmax(ko_category_perturbation$HC_perm_p, 1 / (n_perm + 1))
)


#-----------------------------------------------------------------#
# 3. Summary perturbation plot
#-----------------------------------------------------------------#

ko_category_perturbation$category_label <- factor(
  ko_category_perturbation$category_label,
  levels = rev(ko_category_perturbation$category_label)
)

p_ko_perturbation_summary <-
  ggplot2::ggplot(
    ko_category_perturbation,
    ggplot2::aes(
      x = category_label,
      y = HC_z
    )
  ) +
  
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey60"
  ) +
  
  ggplot2::geom_point(
    ggplot2::aes(
      size = n_KO,
      fill = minus_log10_HC_p
    ),
    shape = 21,
    color = "black",
    stroke = 0.35,
    alpha = 0.92
  ) +
  
  ggplot2::coord_flip() +
  
  ggplot2::scale_fill_gradient(
    low = "white",
    high = "#B64E4E",
    name = "-log10\nperm. P"
  ) +
  
  ggplot2::scale_size_continuous(
    range = c(2.5, 8.5),
    name = "Number\nof KOs"
  ) +
  
  ggplot2::labs(
    x = NULL,
    y = "Higher Criticism perturbation score, z",
    title = "KO category-level perturbation sensitivity"
  ) +
  
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 8.5, color = "black"),
    axis.text.x = ggplot2::element_text(size = 9, color = "black"),
    axis.title.x = ggplot2::element_text(size = 10),
    plot.title = ggplot2::element_text(size = 11, face = "bold"),
    legend.position = "right",
    plot.margin = ggplot2::margin(4, 6, 4, 4)
  )

p_ko_perturbation_summary


#-----------------------------------------------------------------#
# 4. Enrichment curve plot for top categories
#-----------------------------------------------------------------#

top_curve_categories <- ko_category_perturbation$category[
  seq_len(min(max_curve_categories, nrow(ko_category_perturbation)))
]

ko_category_curve_plot <- ko_category_curve[
  ko_category_curve$category %in% top_curve_categories,
  ,
  drop = FALSE
]

ko_category_curve_plot$category_label <- ko_category_perturbation$category_label[
  match(
    ko_category_curve_plot$category,
    ko_category_perturbation$category
  )
]

ko_category_curve_plot$category_label <- factor(
  ko_category_curve_plot$category_label,
  levels = ko_category_perturbation$category_label[
    ko_category_perturbation$category %in% top_curve_categories
  ]
)

p_ko_perturbation_curve <-
  ggplot2::ggplot(
    ko_category_curve_plot,
    ggplot2::aes(
      x = neg_log10_threshold
    )
  ) +
  
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = null_low,
      ymax = null_high
    ),
    fill = "grey85",
    alpha = 0.85
  ) +
  
  ggplot2::geom_line(
    ggplot2::aes(y = null_mean),
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey40"
  ) +
  
  ggplot2::geom_line(
    ggplot2::aes(y = observed_fraction),
    linewidth = 0.75,
    color = "black"
  ) +
  
  ggplot2::geom_point(
    ggplot2::aes(y = observed_fraction),
    size = 1.8,
    color = "black"
  ) +
  
  ggplot2::facet_wrap(
    ~ category_label,
    scales = "free_y",
    ncol = 2
  ) +
  
  ggplot2::scale_x_continuous(
    breaks = -log10(c(0.20, 0.10, 0.05, 0.01, 0.001)),
    labels = c("0.20", "0.10", "0.05", "0.01", "0.001")
  ) +
  
  ggplot2::labs(
    x = "P-value threshold",
    y = "Fraction of KOs with P below threshold",
    title = "Category-level cumulative perturbation curves"
  ) +
  
  ggplot2::theme_classic(base_size = 10.5) +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 8.5, face = "bold"),
    axis.text.x = ggplot2::element_text(size = 8, color = "black"),
    axis.text.y = ggplot2::element_text(size = 8, color = "black"),
    axis.title = ggplot2::element_text(size = 9),
    plot.title = ggplot2::element_text(size = 11, face = "bold"),
    panel.spacing = grid::unit(0.7, "lines"),
    plot.margin = ggplot2::margin(4, 6, 4, 4)
  )

p_ko_perturbation_curve


#-----------------------------------------------------------------#
# 5. Save outputs
#-----------------------------------------------------------------#

ggplot2::ggsave(
  filename = "figures/KO_category_perturbation_HC_summary.svg",
  plot = p_ko_perturbation_summary,
  device = svglite::svglite,
  width = 6.6,
  height = 5.8,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  filename = "figures/KO_category_perturbation_cumulative_curves.svg",
  plot = p_ko_perturbation_curve,
  device = svglite::svglite,
  width = 7.0,
  height = 6.4,
  units = "in",
  bg = "white"
)

write.csv(
  ko_category_perturbation,
  file = "figures/KO_category_perturbation_HC_summary.csv",
  row.names = FALSE
)

write.csv(
  ko_category_curve,
  file = "figures/KO_category_perturbation_cumulative_curve_data.csv",
  row.names = FALSE
)


#-----------------------------------------------------------------#
# 6. Keep final objects
#-----------------------------------------------------------------#

ko_perturbation_results <- list(
  input_table = ko_perturb,
  summary = ko_category_perturbation,
  curve = ko_category_curve,
  summary_plot = p_ko_perturbation_summary,
  curve_plot = p_ko_perturbation_curve,
  settings = list(
    category_col = category_col,
    p_col = p_col,
    direction_col = direction_col,
    n_perm = n_perm,
    min_set_size = min_set_size,
    threshold_grid = threshold_grid
  )
)

ko_perturbation_results$summary













#----------------------------#
#-----------------------------------------------------------------#
#                                                                 #
#   Robust KO category perturbation analysis                      #
#   using size-matched permutation and trimmed statistics          #
#                                                                 #
#-----------------------------------------------------------------#

## Input:
##   ko_selected
##
## Required columns:
##   KO
##   Wilcoxon_p
##   Log2FC
##   assigned_category
##
## Main idea:
##   Test whether each KO category shows a global shift toward smaller
##   CR/non-CR P-values compared with random KO sets of the same size.
##
## Primary statistic:
##   AUC of cumulative low-P enrichment curve.
##
## Robustness statistics:
##   1) truncated Higher Criticism
##   2) trimmed mean -log10(P)
##   3) top-50% mean -log10(P), after dropping the single most extreme KO
##
## This avoids relying on a single binary P-value threshold.

set.seed(123)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

category_col <- "assigned_category"
p_col <- "Wilcoxon_p"
direction_col <- "Log2FC"

n_perm <- 10000
min_set_size <- 15

## Thresholds used to build the cumulative perturbation curve.
## Avoid extremely tiny thresholds because they are often driven by one KO.
threshold_grid <- c(
  0.005,
  0.01,
  0.02,
  0.05,
  0.10,
  0.20,
  0.30,
  0.50
)

## For truncated Higher Criticism.
## P-values smaller than hc_min_p are not used for HC calculation,
## reducing the influence of a single extremely low P-value.
hc_min_p <- 0.01
hc_max_p <- 0.50

## For trimmed mean -log10(P).
trim_fraction <- 0.10

## For top-half statistic.
## If a category has enough KOs, remove the single most extreme KO first.
drop_top_extreme_if_n_at_least <- 10

max_curve_categories <- 8


#-----------------------------------------------------------------#
# 1. Prepare input table
#-----------------------------------------------------------------#

ko_perturb <- as.data.frame(ko_selected, stringsAsFactors = FALSE)

ko_perturb <- ko_perturb[
  is.finite(ko_perturb[[p_col]]) &
    !is.na(ko_perturb[[category_col]]) &
    ko_perturb[[category_col]] != "",
  ,
  drop = FALSE
]

ko_perturb$p_value <- as.numeric(ko_perturb[[p_col]])

ko_perturb$p_value <- pmin(
  pmax(ko_perturb$p_value, 1e-300),
  1
)

ko_perturb$log_p <- -log10(ko_perturb$p_value)
ko_perturb$category <- as.character(ko_perturb[[category_col]])

if (direction_col %in% colnames(ko_perturb)) {
  ko_perturb$direction_value <- as.numeric(ko_perturb[[direction_col]])
} else {
  ko_perturb$direction_value <- NA_real_
}

p_all <- ko_perturb$p_value
logp_all <- ko_perturb$log_p

categories <- sort(unique(ko_perturb$category))

ko_category_perturbation <- data.frame()
ko_category_curve <- data.frame()


#-----------------------------------------------------------------#
# 2. Category-wise robust perturbation test
#-----------------------------------------------------------------#

for (cat_name in categories) {
  
  cat_idx <- ko_perturb$category == cat_name
  
  p_obs <- sort(ko_perturb$p_value[cat_idx])
  logp_obs <- -log10(p_obs)
  n_cat <- length(p_obs)
  
  if (n_cat < min_set_size) {
    next
  }
  
  rank_obs <- seq_len(n_cat) / n_cat
  
  #-----------------------------#
  # 2-1. Truncated Higher Criticism
  #-----------------------------#
  
  hc_valid <- p_obs >= hc_min_p &
    p_obs <= hc_max_p &
    rank_obs > p_obs
  
  if (any(hc_valid)) {
    hc_trunc_obs <- max(
      sqrt(n_cat) *
        (rank_obs[hc_valid] - p_obs[hc_valid]) /
        sqrt(p_obs[hc_valid] * (1 - p_obs[hc_valid])),
      na.rm = TRUE
    )
  } else {
    hc_trunc_obs <- 0
  }
  
  #-----------------------------#
  # 2-2. Cumulative AUC score
  #-----------------------------#
  
  observed_curve <- findInterval(threshold_grid, p_obs) / n_cat
  
  ## This AUC is a discrete average over multiple P-value thresholds.
  ## It is less sensitive to a single extreme P-value than classic HC.
  auc_obs <- mean(observed_curve, na.rm = TRUE)
  
  #-----------------------------#
  # 2-3. Trimmed mean -log10(P)
  #-----------------------------#
  
  trimmed_mean_logp_obs <- mean(
    logp_obs,
    trim = trim_fraction,
    na.rm = TRUE
  )
  
  #-----------------------------#
  # 2-4. Top 50% mean -log10(P), robust version
  #-----------------------------#
  
  logp_obs_sorted <- sort(logp_obs, decreasing = TRUE)
  
  drop_top_n <- ifelse(
    n_cat >= drop_top_extreme_if_n_at_least,
    1,
    0
  )
  
  top50_n <- max(1, floor(n_cat * 0.50))
  
  top50_index <- seq(
    from = drop_top_n + 1,
    length.out = min(top50_n, n_cat - drop_top_n)
  )
  
  top50_mean_logp_obs <- mean(
    logp_obs_sorted[top50_index],
    na.rm = TRUE
  )
  
  #-----------------------------#
  # 2-5. Directional summaries
  #-----------------------------#
  
  median_log2fc_obs <- stats::median(
    ko_perturb$direction_value[cat_idx],
    na.rm = TRUE
  )
  
  mean_abs_log2fc_obs <- mean(
    abs(ko_perturb$direction_value[cat_idx]),
    na.rm = TRUE
  )
  
  prop_005_obs <- mean(p_obs <= 0.05, na.rm = TRUE)
  prop_020_obs <- mean(p_obs <= 0.20, na.rm = TRUE)
  
  n_p005_obs <- sum(p_obs <= 0.05, na.rm = TRUE)
  n_p020_obs <- sum(p_obs <= 0.20, na.rm = TRUE)
  
  #-----------------------------#
  # 2-6. Size-matched permutation
  #-----------------------------#
  
  null_hc_trunc <- numeric(n_perm)
  null_auc <- numeric(n_perm)
  null_trimmed_mean_logp <- numeric(n_perm)
  null_top50_mean_logp <- numeric(n_perm)
  
  null_curve <- matrix(
    NA_real_,
    nrow = n_perm,
    ncol = length(threshold_grid)
  )
  
  for (b in seq_len(n_perm)) {
    
    p_null <- sort(sample(p_all, n_cat, replace = FALSE))
    logp_null <- -log10(p_null)
    
    rank_null <- seq_len(n_cat) / n_cat
    
    hc_null_valid <- p_null >= hc_min_p &
      p_null <= hc_max_p &
      rank_null > p_null
    
    if (any(hc_null_valid)) {
      null_hc_trunc[b] <- max(
        sqrt(n_cat) *
          (rank_null[hc_null_valid] - p_null[hc_null_valid]) /
          sqrt(p_null[hc_null_valid] * (1 - p_null[hc_null_valid])),
        na.rm = TRUE
      )
    } else {
      null_hc_trunc[b] <- 0
    }
    
    null_curve[b, ] <- findInterval(threshold_grid, p_null) / n_cat
    
    null_auc[b] <- mean(null_curve[b, ], na.rm = TRUE)
    
    null_trimmed_mean_logp[b] <- mean(
      logp_null,
      trim = trim_fraction,
      na.rm = TRUE
    )
    
    logp_null_sorted <- sort(logp_null, decreasing = TRUE)
    
    null_top50_mean_logp[b] <- mean(
      logp_null_sorted[top50_index],
      na.rm = TRUE
    )
  }
  
  #-----------------------------#
  # 2-7. Empirical P-values
  #-----------------------------#
  
  hc_trunc_perm_p <- (
    sum(null_hc_trunc >= hc_trunc_obs, na.rm = TRUE) + 1
  ) / (n_perm + 1)
  
  auc_perm_p <- (
    sum(null_auc >= auc_obs, na.rm = TRUE) + 1
  ) / (n_perm + 1)
  
  trimmed_mean_logp_perm_p <- (
    sum(null_trimmed_mean_logp >= trimmed_mean_logp_obs, na.rm = TRUE) + 1
  ) / (n_perm + 1)
  
  top50_mean_logp_perm_p <- (
    sum(null_top50_mean_logp >= top50_mean_logp_obs, na.rm = TRUE) + 1
  ) / (n_perm + 1)
  
  hc_trunc_z <- ifelse(
    stats::sd(null_hc_trunc, na.rm = TRUE) > 0,
    (
      hc_trunc_obs - mean(null_hc_trunc, na.rm = TRUE)
    ) / stats::sd(null_hc_trunc, na.rm = TRUE),
    NA_real_
  )
  
  auc_z <- ifelse(
    stats::sd(null_auc, na.rm = TRUE) > 0,
    (
      auc_obs - mean(null_auc, na.rm = TRUE)
    ) / stats::sd(null_auc, na.rm = TRUE),
    NA_real_
  )
  
  trimmed_mean_logp_z <- ifelse(
    stats::sd(null_trimmed_mean_logp, na.rm = TRUE) > 0,
    (
      trimmed_mean_logp_obs -
        mean(null_trimmed_mean_logp, na.rm = TRUE)
    ) / stats::sd(null_trimmed_mean_logp, na.rm = TRUE),
    NA_real_
  )
  
  top50_mean_logp_z <- ifelse(
    stats::sd(null_top50_mean_logp, na.rm = TRUE) > 0,
    (
      top50_mean_logp_obs -
        mean(null_top50_mean_logp, na.rm = TRUE)
    ) / stats::sd(null_top50_mean_logp, na.rm = TRUE),
    NA_real_
  )
  
  #-----------------------------#
  # 2-8. Store summary
  #-----------------------------#
  
  ko_category_perturbation <- rbind(
    ko_category_perturbation,
    data.frame(
      category = cat_name,
      n_KO = n_cat,
      
      AUC_score = auc_obs,
      AUC_null_mean = mean(null_auc, na.rm = TRUE),
      AUC_z = auc_z,
      AUC_perm_p = auc_perm_p,
      
      HC_trunc_score = hc_trunc_obs,
      HC_trunc_null_mean = mean(null_hc_trunc, na.rm = TRUE),
      HC_trunc_z = hc_trunc_z,
      HC_trunc_perm_p = hc_trunc_perm_p,
      
      trimmed_mean_log10P = trimmed_mean_logp_obs,
      trimmed_mean_log10P_null_mean = mean(null_trimmed_mean_logp, na.rm = TRUE),
      trimmed_mean_log10P_z = trimmed_mean_logp_z,
      trimmed_mean_log10P_perm_p = trimmed_mean_logp_perm_p,
      
      top50_mean_log10P = top50_mean_logp_obs,
      top50_mean_log10P_null_mean = mean(null_top50_mean_logp, na.rm = TRUE),
      top50_mean_log10P_z = top50_mean_logp_z,
      top50_mean_log10P_perm_p = top50_mean_logp_perm_p,
      
      prop_P005 = prop_005_obs,
      prop_P020 = prop_020_obs,
      n_P005 = n_p005_obs,
      n_P020 = n_p020_obs,
      
      median_Log2FC = median_log2fc_obs,
      mean_abs_Log2FC = mean_abs_log2fc_obs,
      
      low_count_flag = n_cat < 20,
      single_hit_flag = n_p020_obs < 2,
      
      stringsAsFactors = FALSE
    )
  )
  
  #-----------------------------#
  # 2-9. Store curve data
  #-----------------------------#
  
  ko_category_curve <- rbind(
    ko_category_curve,
    data.frame(
      category = cat_name,
      n_KO = n_cat,
      threshold = threshold_grid,
      neg_log10_threshold = -log10(threshold_grid),
      observed_fraction = observed_curve,
      null_mean = apply(null_curve, 2, mean, na.rm = TRUE),
      null_low = apply(null_curve, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
      null_high = apply(null_curve, 2, stats::quantile, probs = 0.975, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  )
}


#-----------------------------------------------------------------#
# 3. Multiple-testing correction and ranking
#-----------------------------------------------------------------#

ko_category_perturbation$AUC_q <- stats::p.adjust(
  ko_category_perturbation$AUC_perm_p,
  method = "BH"
)

ko_category_perturbation$HC_trunc_q <- stats::p.adjust(
  ko_category_perturbation$HC_trunc_perm_p,
  method = "BH"
)

ko_category_perturbation$trimmed_mean_log10P_q <- stats::p.adjust(
  ko_category_perturbation$trimmed_mean_log10P_perm_p,
  method = "BH"
)

ko_category_perturbation$top50_mean_log10P_q <- stats::p.adjust(
  ko_category_perturbation$top50_mean_log10P_perm_p,
  method = "BH"
)

ko_category_perturbation <- ko_category_perturbation[
  order(
    -ko_category_perturbation$AUC_z,
    ko_category_perturbation$AUC_perm_p,
    -ko_category_perturbation$top50_mean_log10P_z,
    -ko_category_perturbation$HC_trunc_z
  ),
  ,
  drop = FALSE
]

ko_category_perturbation$category_label <- sapply(
  ko_category_perturbation$category,
  function(x) {
    paste(strwrap(x, width = 34), collapse = "\n")
  }
)

ko_category_perturbation$category_label <- factor(
  ko_category_perturbation$category_label,
  levels = rev(ko_category_perturbation$category_label)
)

ko_category_perturbation$minus_log10_AUC_p <- -log10(
  pmax(ko_category_perturbation$AUC_perm_p, 1 / (n_perm + 1))
)

ko_category_perturbation$evidence_class <- ifelse(
  ko_category_perturbation$single_hit_flag,
  "single-hit sensitive",
  ifelse(
    ko_category_perturbation$low_count_flag,
    "low-count category",
    "distributed signal"
  )
)


#-----------------------------------------------------------------#
# 4. Summary plot
#-----------------------------------------------------------------#

ko_category_perturbation$category_label <- factor(
  ko_category_perturbation$category_label,
  levels = rev(ko_category_perturbation$category_label)
)

p_ko_perturbation_summary <-
  ggplot2::ggplot(
    ko_category_perturbation,
    ggplot2::aes(
      x = category_label,
      y = AUC_z
    )
  ) +
  
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey60"
  ) +
  
  ggplot2::geom_point(
    ggplot2::aes(
      size = n_KO,
      fill = minus_log10_AUC_p,
      shape = evidence_class
    ),
    color = "black",
    stroke = 0.35,
    alpha = 0.92
  ) +
  
  ggplot2::coord_flip() +
  
  ggplot2::scale_fill_gradient(
    low = "white",
    high = "#B64E4E",
    name = "-log10\nperm. P"
  ) +
  
  ggplot2::scale_size_continuous(
    range = c(2.4, 8.2),
    name = "Number\nof KOs"
  ) +
  
  ggplot2::scale_shape_manual(
    values = c(
      "distributed signal" = 21,
      "low-count category" = 24,
      "single-hit sensitive" = 22
    ),
    name = "Evidence type"
  ) +
  
  ggplot2::labs(
    x = NULL,
    y = "AUC-based perturbation score, z",
    title = "KO category-level perturbation sensitivity"
  ) +
  
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    ## category labels after coord_flip()
    axis.text.y = ggplot2::element_text(
      size = 9.8,
      color = "black"
    ),
    
    ## numeric x-axis labels after coord_flip()
    axis.text.x = ggplot2::element_text(
      size = 8.2,
      color = "black"
    ),
    
    axis.title.x = ggplot2::element_text(
      size = 9.2
    ),
    
    plot.title = ggplot2::element_text(
      size = 11,
      face = "bold"
    ),
    
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 8.7),
    legend.text = ggplot2::element_text(size = 7.6),
    
    plot.margin = ggplot2::margin(4, 6, 4, 4)
  )



p_ko_perturbation_summary


#-----------------------------------------------------------------#
# 5. Cumulative perturbation curve plot
#-----------------------------------------------------------------#

top_curve_categories <- ko_category_perturbation$category[
  seq_len(min(max_curve_categories, nrow(ko_category_perturbation)))
]

ko_category_curve_plot <- ko_category_curve[
  ko_category_curve$category %in% top_curve_categories,
  ,
  drop = FALSE
]

ko_category_curve_plot$category_label <- ko_category_perturbation$category_label[
  match(
    ko_category_curve_plot$category,
    ko_category_perturbation$category
  )
]

ko_category_curve_plot$category_label <- factor(
  ko_category_curve_plot$category_label,
  levels = ko_category_perturbation$category_label[
    ko_category_perturbation$category %in% top_curve_categories
  ]
)

p_ko_perturbation_curve <-
  ggplot2::ggplot(
    ko_category_curve_plot,
    ggplot2::aes(
      x = neg_log10_threshold
    )
  ) +
  
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = null_low,
      ymax = null_high
    ),
    fill = "grey85",
    alpha = 0.85
  ) +
  
  ggplot2::geom_line(
    ggplot2::aes(y = null_mean),
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey40"
  ) +
  
  ggplot2::geom_line(
    ggplot2::aes(y = observed_fraction),
    linewidth = 0.75,
    color = "black"
  ) +
  
  ggplot2::geom_point(
    ggplot2::aes(y = observed_fraction),
    size = 1.8,
    color = "black"
  ) +
  
  ggplot2::facet_wrap(
    ~ category_label,
    scales = "free_y",
    ncol = 2
  ) +
  
  ggplot2::scale_x_continuous(
    breaks = -log10(c(0.50, 0.20, 0.10, 0.05, 0.01, 0.005)),
    labels = c("0.50", "0.20", "0.10", "0.05", "0.01", "0.005")
  ) +
  
  ggplot2::labs(
    x = "P-value threshold",
    y = "Fraction of KOs with P below threshold",
    title = "Cumulative KO perturbation curves"
  ) +
  
  ggplot2::theme_classic(base_size = 10.5) +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 8.5, face = "bold"),
    axis.text.x = ggplot2::element_text(size = 8, color = "black"),
    axis.text.y = ggplot2::element_text(size = 8, color = "black"),
    axis.title = ggplot2::element_text(size = 9),
    plot.title = ggplot2::element_text(size = 11, face = "bold"),
    panel.spacing = grid::unit(0.7, "lines"),
    plot.margin = ggplot2::margin(4, 6, 4, 4)
  )

p_ko_perturbation_curve


#-----------------------------------------------------------------#
# 6. Save outputs
#-----------------------------------------------------------------#

ggplot2::ggsave(
  filename = "figures/KO_category_robust_perturbation_summary.svg",
  plot = p_ko_perturbation_summary,
  device = svglite::svglite,
  width = 5.1,
  height = 5.8,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  filename = "figures/KO_category_robust_perturbation_curves.svg",
  plot = p_ko_perturbation_curve,
  device = svglite::svglite,
  width = 7.0,
  height = 6.4,
  units = "in",
  bg = "white"
)

write.csv(
  ko_category_perturbation,
  file = "figures/KO_category_robust_perturbation_summary.csv",
  row.names = FALSE
)

write.csv(
  ko_category_curve,
  file = "figures/KO_category_robust_perturbation_curve_data.csv",
  row.names = FALSE
)


#-----------------------------------------------------------------#
# 7. Store final object
#-----------------------------------------------------------------#

ko_perturbation_results <- list(
  input_table = ko_perturb,
  summary = ko_category_perturbation,
  curve = ko_category_curve,
  summary_plot = p_ko_perturbation_summary,
  curve_plot = p_ko_perturbation_curve,
  settings = list(
    category_col = category_col,
    p_col = p_col,
    direction_col = direction_col,
    n_perm = n_perm,
    min_set_size = min_set_size,
    threshold_grid = threshold_grid,
    hc_min_p = hc_min_p,
    hc_max_p = hc_max_p,
    trim_fraction = trim_fraction
  )
)

ko_perturbation_results$summary