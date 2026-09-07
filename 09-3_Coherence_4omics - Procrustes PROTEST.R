#-----------------------------------------------------------------#
#
# Four-omics pairwise Procrustes / PROTEST analysis
#
# Input:
#   input/coherence_species_ko_metabolite_host_data_metabolite50_host50pct_no_scaling.RData
#
# Main biological-chain comparisons:
#   Species vs KEGG ortholog
#   KEGG ortholog vs Metabolite
#   Metabolite vs Host RNA-seq
#
# Outputs:
#   1. Four-by-four matrix: upper Procrustes scatter and lower summary
#   2. Upper-triangle-only Procrustes matrix with in-panel statistics
#   3. Compact forest-style statistic summary
#   4. Four-block GPA consensus with three biological-chain pairwise panels
#   5. Sample-level distance-to-consensus violin/box/beeswarm plot
#   6. Modality agreement with SubjectID-cluster bootstrap 95% CI
#
# Ordination preprocessing:
#   - This script uses the pair-specific PCoA coordinates stored by the
#     original input-processing workflow.
#   - Species / KO: Bray-Curtis distance on relative abundance.
#   - Metabolite: Euclidean distance on the transformed and feature-scaled
#     concentration table.
#   - Host RNA-seq: genes selected by cumulative VST variance (primary 50%);
#     no gene-wise z-scaling; Euclidean distance followed by PCoA.
#   - Procrustes: symmetric = TRUE, which equalizes configuration scale.
#
# Sample matching:
#   - Every plot uses all SubjectID-Timepoint samples containing both data
#     types in the relevant pair.
#   - Baseline-After RT completeness is not required.
#   - Baseline and After RT are not separated in the scatter plots.
#
# Interpretation:
#   - Reversing the order of a pair does not provide an independent test
#     under symmetric Procrustes. The nine-panel ordered display is therefore
#     a visualization supplement, not nine independent comparisons.
#
#-----------------------------------------------------------------#

rm(list = ls())
options(stringsAsFactors = FALSE)
setwd("D:/2-연구/2-CRC metagenomics/")

source(
  "7-11.0 Coherence_4omics - common utilities.R"
)

required_packages <- c(
  "vegan",
  "ape",
  "ggplot2",
  "patchwork",
  "svglite",
  "scales",
  "ggbeeswarm"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages first: ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    call. = FALSE
  )
}

dir.create(
  "figures/coherence",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/coherence",
  recursive = TRUE,
  showWarnings = FALSE
)

coherence_data <- coherence_load_data(
  path = paste0(
    "input/coherence_species_ko_metabolite_host_data_",
    "metabolite50_host50pct_no_scaling.RData"
  )
)

if (
  is.null(
    coherence_data$host_sensitivity
  ) ||
  is.null(
    coherence_data$host_sensitivity$host_feature_scaling
  ) ||
  !identical(
    as.character(
      coherence_data$host_sensitivity$host_feature_scaling
    ),
    "none"
  )
) {
  stop(
    paste0(
      "The loaded input does not document host_feature_scaling = 'none'. ",
      "Run input processing v6 before this Procrustes script."
    ),
    call. = FALSE
  )
}

if (
  is.null(
    coherence_data$analysis_settings$host_selection_mode
  ) ||
  !identical(
    as.character(
      coherence_data$analysis_settings$host_selection_mode
    ),
    "cumulative_variance"
  )
) {
  stop(
    paste0(
      "The loaded input was not generated with cumulative-variance ",
      "Host RNA-seq feature selection."
    ),
    call. = FALSE
  )
}

if (
  is.null(
    coherence_data$host_sensitivity$primary_cumulative_fraction
  ) ||
  !isTRUE(
    all.equal(
      as.numeric(
        coherence_data$host_sensitivity$primary_cumulative_fraction
      ),
      0.50,
      tolerance = 1e-8
    )
  )
) {
  stop(
    paste0(
      "The loaded input is not the primary 50% cumulative VST-variance ",
      "Host RNA-seq object."
    ),
    call. = FALSE
  )
}

required_pair_objects <- c(
  "species_ko",
  "species_metabolite",
  "species_host",
  "ko_metabolite",
  "ko_host",
  "metabolite_host"
)

missing_pair_objects <- setdiff(
  required_pair_objects,
  names(coherence_data$pairs)
)

if (length(missing_pair_objects) > 0) {
  stop(
    paste0(
      "The input RData lacks the following pair objects: ",
      paste(
        missing_pair_objects,
        collapse = ", "
      ),
      ". Re-run the updated all-pairs input-processing script."
    ),
    call. = FALSE
  )
}

set.seed(20260718)

protest_permutations <- 9999
response_comparison_permutations <- 999
host_gene_sensitivity_permutations <- 999
protest_max_axes <- 5
gpa_max_axes <- 5
gpa_permutations <- 9999
gpa_agreement_bootstraps <- 999
gpa_tolerance <- 1e-10
gpa_max_iterations <- 500

layer_labels <- c(
  species = "Species",
  ko = "KEGG ortholog",
  metabolite = "Metabolite",
  host = "Host RNA-seq"
)

# Muted, response-neutral pastel colors. Color and shape both encode data type.
layer_colors <- c(
  Species = "#79B3A3",
  `KEGG ortholog` = "#D8A46F",
  Metabolite = "#82A8C7",
  `Host RNA-seq` = "#B29AC6"
)

layer_shapes <- c(
  Species = 21,
  `KEGG ortholog` = 24,
  Metabolite = 22,
  `Host RNA-seq` = 23
)


#=================================================================#
# Original pair-specific ordination
#=================================================================#
#
# The input-processing workflow already generated pair-specific PCoA scores
# using the intended modality-specific distances. This script does not rebuild
# or rescale the feature tables before ordination.
#
# Species / KO:
#   relative abundance -> Bray-Curtis distance -> PCoA
#
# Metabolite:
#   transformed, feature-scaled concentration table
#   -> Euclidean distance -> PCoA
#
# Host RNA-seq:
#   DESeq2 blind VST
#   -> genes ranked by sample-wise VST variance
#   -> smallest gene set reaching 50% cumulative gene-wise variance
#   -> no gene-wise z-scaling
#   -> Euclidean distance -> PCoA
#
# Euclidean PCoA of the unscaled VST matrix has the same sample geometry as
# centered PCA, apart from arbitrary axis sign and rotation.
#
# The same pair-specific common samples are used for both configurations.
# Symmetric Procrustes subsequently standardizes the total configuration scale,
# but it does not alter the internal geometry established by each distance.
#
#-----------------------------------------------------------------#

original_score_matrix <- function(
    pair_object,
    block_name
) {
  score_matrix <- coherence_score_matrix(
    pair_object,
    block_name
  )
  
  score_matrix <- as.matrix(
    score_matrix
  )
  
  storage.mode(score_matrix) <- "numeric"
  
  complete_axes <- colSums(
    is.finite(score_matrix)
  ) == nrow(score_matrix)
  
  score_matrix <- score_matrix[
    ,
    complete_axes,
    drop = FALSE
  ]
  
  score_matrix
}

original_ordination_preprocessing <- data.frame(
  block = c(
    "species",
    "ko",
    "metabolite",
    "host"
  ),
  transformation = c(
    "Relative abundance",
    "Relative abundance",
    "Transformed and feature-scaled concentration table",
    paste0(
      "DESeq2 blind VST; cumulative-variance gene selection; ",
      "no gene-wise z-scaling"
    )
  ),
  distance = c(
    "Bray-Curtis",
    "Bray-Curtis",
    "Euclidean",
    "Euclidean"
  ),
  n_samples = vapply(
    coherence_data$modalities[
      c(
        "species",
        "ko",
        "metabolite",
        "host"
      )
    ],
    function(x) {
      nrow(x$sample_meta)
    },
    integer(1)
  ),
  n_features = vapply(
    coherence_data$modalities[
      c(
        "species",
        "ko",
        "metabolite",
        "host"
      )
    ],
    function(x) {
      as.integer(x$n_features)
    },
    integer(1)
  ),
  stringsAsFactors = FALSE
)

original_ordination_preprocessing$host_primary_gene_count <-
  as.integer(
    coherence_data$host_sensitivity$primary_top_n
  )

original_ordination_preprocessing$host_primary_cumulative_variance <-
  as.numeric(
    coherence_data$host_sensitivity$primary_cumulative_fraction
  )

original_ordination_preprocessing$host_gene_wise_scaling <-
  as.character(
    coherence_data$host_sensitivity$host_feature_scaling
  )

primary_pair_specs <- data.frame(
  pair_name = c(
    "species_ko",
    "ko_metabolite",
    "metabolite_host"
  ),
  block_x = c(
    "species",
    "ko",
    "metabolite"
  ),
  block_y = c(
    "ko",
    "metabolite",
    "host"
  ),
  pair_label = c(
    "Species vs KEGG ortholog",
    "KEGG ortholog vs Metabolite",
    "Metabolite vs Host RNA-seq"
  ),
  stringsAsFactors = FALSE
)

all_pair_specs <- data.frame(
  pair_name = c(
    "species_ko",
    "species_metabolite",
    "species_host",
    "ko_metabolite",
    "ko_host",
    "metabolite_host"
  ),
  block_x = c(
    "species",
    "species",
    "species",
    "ko",
    "ko",
    "metabolite"
  ),
  block_y = c(
    "ko",
    "metabolite",
    "host",
    "metabolite",
    "host",
    "host"
  ),
  pair_label = c(
    "Species vs KEGG ortholog",
    "Species vs Metabolite",
    "Species vs Host RNA-seq",
    "KEGG ortholog vs Metabolite",
    "KEGG ortholog vs Host RNA-seq",
    "Metabolite vs Host RNA-seq"
  ),
  stringsAsFactors = FALSE
)

ordered_display_specs <- data.frame(
  pair_name = c(
    "species_ko",
    "species_metabolite",
    "species_host",
    "species_ko",
    "ko_metabolite",
    "ko_host",
    "species_metabolite",
    "ko_metabolite",
    "metabolite_host"
  ),
  block_x = c(
    "species",
    "species",
    "species",
    "ko",
    "ko",
    "ko",
    "metabolite",
    "metabolite",
    "metabolite"
  ),
  block_y = c(
    "ko",
    "metabolite",
    "host",
    "species",
    "metabolite",
    "host",
    "species",
    "ko",
    "host"
  ),
  plot_title = c(
    "Species vs KEGG ortholog",
    "Species vs Metabolite",
    "Species vs Host RNA-seq",
    "KEGG ortholog vs Species",
    "KEGG ortholog vs Metabolite",
    "KEGG ortholog vs Host RNA-seq",
    "Metabolite vs Species",
    "Metabolite vs KEGG ortholog",
    "Metabolite vs Host RNA-seq"
  ),
  row_group = factor(
    rep(
      c(
        "Species as first block",
        "KEGG ortholog as first block",
        "Metabolite as first block"
      ),
      each = 3
    ),
    levels = c(
      "Species as first block",
      "KEGG ortholog as first block",
      "Metabolite as first block"
    )
  ),
  stringsAsFactors = FALSE
)

subset_levels <- c(
  "All",
  "Baseline",
  "After RT",
  "pCR",
  "non-pCR"
)

all_pair_summary <- data.frame()
all_pair_results <- list()

for (i in seq_len(nrow(all_pair_specs))) {
  pair_object <- coherence_get_pair(
    coherence_data,
    all_pair_specs$pair_name[i]
  )
  
  X_all <- original_score_matrix(
    pair_object,
    all_pair_specs$block_x[i]
  )
  
  Y_all <- original_score_matrix(
    pair_object,
    all_pair_specs$block_y[i]
  )
  
  common_all <- pair_object$sample_meta$SampleID[
    pair_object$sample_meta$SampleID %in%
      intersect(
        rownames(X_all),
        rownames(Y_all)
      )
  ]
  
  X_all <- X_all[
    common_all,
    ,
    drop = FALSE
  ]
  
  Y_all <- Y_all[
    common_all,
    ,
    drop = FALSE
  ]
  
  all_pair_results[[all_pair_specs$pair_name[i]]] <- list()
  
  for (subset_name in subset_levels) {
    subset_ids <- coherence_subset_ids(
      pair_object,
      subset_name
    )
    
    common_samples <- pair_object$sample_meta$SampleID[
      pair_object$sample_meta$SampleID %in%
        intersect(
          common_all,
          subset_ids
        )
    ]
    
    axes_used <- min(
      protest_max_axes,
      ncol(X_all),
      ncol(Y_all),
      length(common_samples) - 1
    )
    
    if (
      length(common_samples) < 4 ||
      axes_used < 2
    ) {
      all_pair_summary <- rbind(
        all_pair_summary,
        data.frame(
          pair_name = all_pair_specs$pair_name[i],
          pair = all_pair_specs$pair_label[i],
          subset = subset_name,
          n = length(common_samples),
          axes_used = max(
            0,
            axes_used
          ),
          procrustes_ss = NA_real_,
          protest_r = NA_real_,
          protest_p = NA_real_,
          stringsAsFactors = FALSE
        )
      )
      
      next
    }
    
    X <- X_all[
      common_samples,
      seq_len(axes_used),
      drop = FALSE
    ]
    
    Y <- Y_all[
      common_samples,
      seq_len(axes_used),
      drop = FALSE
    ]
    
    proc <- vegan::procrustes(
      X,
      Y,
      symmetric = TRUE
    )
    
    set.seed(
      20260718 +
        i * 100 +
        match(
          subset_name,
          subset_levels
        )
    )
    
    prot <- vegan::protest(
      X,
      Y,
      permutations = protest_permutations,
      symmetric = TRUE
    )
    
    all_pair_results[[all_pair_specs$pair_name[i]]][[subset_name]] <- list(
      procrustes = proc,
      protest = prot,
      samples = common_samples,
      axes_used = axes_used
    )
    
    all_pair_summary <- rbind(
      all_pair_summary,
      data.frame(
        pair_name = all_pair_specs$pair_name[i],
        pair = all_pair_specs$pair_label[i],
        subset = subset_name,
        n = length(common_samples),
        axes_used = axes_used,
        procrustes_ss = unname(
          proc$ss
        ),
        protest_r = unname(
          prot$t0
        ),
        protest_p = unname(
          prot$signif
        ),
        stringsAsFactors = FALSE
      )
    )
  }
}

rownames(all_pair_summary) <- NULL

all_pair_summary$subset <- factor(
  all_pair_summary$subset,
  levels = subset_levels
)

all_pair_summary$pair <- factor(
  all_pair_summary$pair,
  levels = all_pair_specs$pair_label
)


#=================================================================#
# 1. Pair-level summaries for the pooled "All" comparison
#=================================================================#

all_summary_all <- all_pair_summary[
  all_pair_summary$subset == "All",
  ,
  drop = FALSE
]

all_summary_all$pair <- factor(
  as.character(all_summary_all$pair),
  levels = all_pair_specs$pair_label
)

all_summary_all$protest_q <- stats::p.adjust(
  all_summary_all$protest_p,
  method = "BH"
)

pair_group_stats <- data.frame()

for (i in seq_len(nrow(all_pair_specs))) {
  pair_object <- coherence_get_pair(
    coherence_data,
    all_pair_specs$pair_name[i]
  )
  
  all_result <- all_pair_results[[all_pair_specs$pair_name[i]]][["All"]]
  
  if (is.null(all_result)) {
    pair_group_stats <- rbind(
      pair_group_stats,
      data.frame(
        pair_name = all_pair_specs$pair_name[i],
        block_x = all_pair_specs$block_x[i],
        block_y = all_pair_specs$block_y[i],
        pair = all_pair_specs$pair_label[i],
        n_total = 0,
        n_subjects = 0,
        n_pCR = 0,
        n_non_pCR = 0,
        residual_median = NA_real_,
        residual_mean = NA_real_,
        midpoint_centroid_distance = NA_real_,
        midpoint_centroid_distance_p = NA_real_,
        stringsAsFactors = FALSE
      )
    )
    next
  }
  
  proc <- all_result$procrustes
  sample_ids <- all_result$samples
  
  meta_plot <- pair_object$sample_meta[
    match(
      sample_ids,
      pair_object$sample_meta$SampleID
    ),
    ,
    drop = FALSE
  ]
  
  meta_plot$TRG_plot <- factor(
    as.character(meta_plot$TRG_plot),
    levels = c(
      "pCR",
      "non_pCR"
    )
  )
  
  residual_length <- sqrt(
    rowSums(
      (as.matrix(proc$X) - as.matrix(proc$Yrot))^2
    )
  )
  
  midpoint_df <- data.frame(
    SampleID = sample_ids,
    SubjectID = as.character(meta_plot$SubjectID),
    Axis1 = (proc$X[, 1] + proc$Yrot[, 1]) / 2,
    Axis2 = (proc$X[, 2] + proc$Yrot[, 2]) / 2,
    TRG_plot = as.character(meta_plot$TRG_plot),
    stringsAsFactors = FALSE
  )
  
  subject_midpoint_df <- stats::aggregate(
    cbind(
      Axis1,
      Axis2
    ) ~ SubjectID + TRG_plot,
    data = midpoint_df,
    FUN = mean
  )
  
  subject_midpoint_df$TRG_plot <- factor(
    subject_midpoint_df$TRG_plot,
    levels = c(
      "pCR",
      "non_pCR"
    )
  )
  
  group_p <- NA_real_
  group_distance <- NA_real_
  
  if (
    sum(subject_midpoint_df$TRG_plot == "pCR") >= 2 &&
    sum(subject_midpoint_df$TRG_plot == "non_pCR") >= 2
  ) {
    centroid_pcr <- c(
      mean(
        subject_midpoint_df$Axis1[subject_midpoint_df$TRG_plot == "pCR"],
        na.rm = TRUE
      ),
      mean(
        subject_midpoint_df$Axis2[subject_midpoint_df$TRG_plot == "pCR"],
        na.rm = TRUE
      )
    )
    
    centroid_non <- c(
      mean(
        subject_midpoint_df$Axis1[subject_midpoint_df$TRG_plot == "non_pCR"],
        na.rm = TRUE
      ),
      mean(
        subject_midpoint_df$Axis2[subject_midpoint_df$TRG_plot == "non_pCR"],
        na.rm = TRUE
      )
    )
    
    group_distance <- sqrt(
      sum(
        (centroid_pcr - centroid_non)^2
      )
    )
    
    set.seed(
      20262718 + i
    )
    
    perm_distance <- rep(
      NA_real_,
      response_comparison_permutations
    )
    
    for (perm_index in seq_len(response_comparison_permutations)) {
      perm_group <- sample(
        subject_midpoint_df$TRG_plot
      )
      
      perm_centroid_pcr <- c(
        mean(
          subject_midpoint_df$Axis1[perm_group == "pCR"],
          na.rm = TRUE
        ),
        mean(
          subject_midpoint_df$Axis2[perm_group == "pCR"],
          na.rm = TRUE
        )
      )
      
      perm_centroid_non <- c(
        mean(
          subject_midpoint_df$Axis1[perm_group == "non_pCR"],
          na.rm = TRUE
        ),
        mean(
          subject_midpoint_df$Axis2[perm_group == "non_pCR"],
          na.rm = TRUE
        )
      )
      
      perm_distance[perm_index] <- sqrt(
        sum(
          (perm_centroid_pcr - perm_centroid_non)^2
        )
      )
    }
    
    group_p <- (
      1 +
        sum(
          perm_distance >= group_distance,
          na.rm = TRUE
        )
    ) / (
      1 +
        sum(
          is.finite(perm_distance)
        )
    )
  }
  
  pair_group_stats <- rbind(
    pair_group_stats,
    data.frame(
      pair_name = all_pair_specs$pair_name[i],
      block_x = all_pair_specs$block_x[i],
      block_y = all_pair_specs$block_y[i],
      pair = all_pair_specs$pair_label[i],
      n_total = length(sample_ids),
      n_subjects = nrow(subject_midpoint_df),
      n_pCR = sum(subject_midpoint_df$TRG_plot == "pCR"),
      n_non_pCR = sum(subject_midpoint_df$TRG_plot == "non_pCR"),
      residual_median = stats::median(
        residual_length,
        na.rm = TRUE
      ),
      residual_mean = mean(
        residual_length,
        na.rm = TRUE
      ),
      midpoint_centroid_distance = group_distance,
      midpoint_centroid_distance_p = group_p,
      stringsAsFactors = FALSE
    )
  )
}

pair_group_stats$pair <- factor(
  as.character(pair_group_stats$pair),
  levels = all_pair_specs$pair_label
)

pair_display_summary <- merge(
  all_summary_all,
  pair_group_stats,
  by = c(
    "pair_name",
    "pair"
  ),
  all.x = TRUE,
  sort = FALSE
)

pair_display_summary$pair <- factor(
  as.character(pair_display_summary$pair),
  levels = all_pair_specs$pair_label
)

matrix_summary_bootstraps <- 999
pair_r_bootstrap_summary <- data.frame()

for (i in seq_len(nrow(all_pair_specs))) {
  pair_name_current <-
    all_pair_specs$pair_name[i]
  
  pair_object <- coherence_get_pair(
    coherence_data,
    pair_name_current
  )
  
  X_all <- original_score_matrix(
    pair_object,
    all_pair_specs$block_x[i]
  )
  
  Y_all <- original_score_matrix(
    pair_object,
    all_pair_specs$block_y[i]
  )
  
  common_samples <- pair_object$sample_meta$SampleID[
    pair_object$sample_meta$SampleID %in%
      intersect(
        rownames(X_all),
        rownames(Y_all)
      )
  ]
  
  axes_used <- min(
    protest_max_axes,
    ncol(X_all),
    ncol(Y_all),
    length(common_samples) - 1
  )
  
  if (
    length(common_samples) < 4 ||
    axes_used < 2
  ) {
    pair_r_bootstrap_summary <- rbind(
      pair_r_bootstrap_summary,
      data.frame(
        pair_name = pair_name_current,
        observed_r = NA_real_,
        r_ci_low = NA_real_,
        r_ci_high = NA_real_,
        bootstrap_se = NA_real_,
        bootstrap_bias = NA_real_,
        ci_method =
          "Observed-r-centered normal interval using SubjectID-cluster bootstrap SE",
        bootstrap_replicates =
          matrix_summary_bootstraps,
        bootstrap_unit =
          "SubjectID cluster",
        bootstrap_seed =
          20267718 + i,
        stringsAsFactors = FALSE
      )
    )
    
    next
  }
  
  X <- X_all[
    common_samples,
    seq_len(axes_used),
    drop = FALSE
  ]
  
  Y <- Y_all[
    common_samples,
    seq_len(axes_used),
    drop = FALSE
  ]
  
  meta_bootstrap <- pair_object$sample_meta[
    match(
      common_samples,
      pair_object$sample_meta$SampleID
    ),
    ,
    drop = FALSE
  ]
  
  bootstrap_clusters <- if (
    "SubjectID" %in%
    colnames(meta_bootstrap)
  ) {
    as.character(
      meta_bootstrap$SubjectID
    )
  } else {
    as.character(
      meta_bootstrap$SampleID
    )
  }
  
  unique_clusters <- unique(
    bootstrap_clusters
  )
  
  set.seed(
    20267718 + i
  )
  
  bootstrap_r <- rep(
    NA_real_,
    matrix_summary_bootstraps
  )
  
  for (
    bootstrap_index in
    seq_len(
      matrix_summary_bootstraps
    )
  ) {
    sampled_clusters <- sample(
      unique_clusters,
      size = length(unique_clusters),
      replace = TRUE
    )
    
    bootstrap_sample_index <- unlist(
      lapply(
        sampled_clusters,
        function(cluster_id) {
          which(
            bootstrap_clusters ==
              cluster_id
          )
        }
      ),
      use.names = FALSE
    )
    
    bootstrap_proc <- tryCatch(
      vegan::procrustes(
        X[
          bootstrap_sample_index,
          ,
          drop = FALSE
        ],
        Y[
          bootstrap_sample_index,
          ,
          drop = FALSE
        ],
        symmetric = TRUE
      ),
      error = function(e) {
        NULL
      }
    )
    
    if (is.null(bootstrap_proc)) {
      next
    }
    
    bootstrap_r[
      bootstrap_index
    ] <- sqrt(
      max(
        0,
        1 -
          as.numeric(
            bootstrap_proc$ss
          )
      )
    )
  }
  
  observed_r <- all_summary_all$protest_r[
    match(
      pair_name_current,
      all_summary_all$pair_name
    )
  ]
  
  finite_bootstrap_r <- bootstrap_r[
    is.finite(
      bootstrap_r
    )
  ]
  
  if (
    length(
      finite_bootstrap_r
    ) >= 50 &&
    is.finite(
      observed_r
    )
  ) {
    bootstrap_se <- stats::sd(
      finite_bootstrap_r,
      na.rm = TRUE
    )
    
    r_ci <- c(
      max(
        0,
        observed_r -
          stats::qnorm(
            0.975
          ) *
          bootstrap_se
      ),
      min(
        1,
        observed_r +
          stats::qnorm(
            0.975
          ) *
          bootstrap_se
      )
    )
    
    bootstrap_bias <- mean(
      finite_bootstrap_r,
      na.rm = TRUE
    ) -
      observed_r
  } else {
    bootstrap_se <- NA_real_
    bootstrap_bias <- NA_real_
    r_ci <- c(
      NA_real_,
      NA_real_
    )
  }
  
  pair_r_bootstrap_summary <- rbind(
    pair_r_bootstrap_summary,
    data.frame(
      pair_name =
        pair_name_current,
      observed_r =
        observed_r,
      r_ci_low =
        r_ci[1],
      r_ci_high =
        r_ci[2],
      bootstrap_se =
        bootstrap_se,
      bootstrap_bias =
        bootstrap_bias,
      ci_method =
        "Observed-r-centered normal interval using SubjectID-cluster bootstrap SE",
      bootstrap_replicates =
        matrix_summary_bootstraps,
      bootstrap_unit =
        "SubjectID cluster",
      bootstrap_seed =
        20267718 + i,
      stringsAsFactors = FALSE
    )
  )
}

pair_display_summary <- merge(
  pair_display_summary,
  pair_r_bootstrap_summary,
  by = "pair_name",
  all.x = TRUE,
  sort = FALSE
)

pair_display_summary$pair <- factor(
  as.character(
    pair_display_summary$pair
  ),
  levels =
    all_pair_specs$pair_label
)


#=================================================================#
# 2. Procrustes matrix visualizations
#=================================================================#
#
# Figure A:
#   - Upper triangle: Procrustes scatter.
#   - Diagonal: data-type labels.
#   - Lower triangle: compact PROTEST-r profile with bootstrap 95% CI,
#     permutation P value, and sample count.
#
# Figure B:
#   - Same diagonal and upper-triangle scatter panels.
#   - Lower triangle removed.
#   - r, 95% CI, P value, and n are printed inside each scatter panel.
#
# pCR/non-pCR is not mapped to any aesthetic. Both fill color and shape
# distinguish the four data types.
#
#-----------------------------------------------------------------#

modality_order <- c(
  "species",
  "ko",
  "metabolite",
  "host"
)

modality_labels <- unname(
  layer_labels[modality_order]
)

format_matrix_p <- function(p_value) {
  vapply(
    p_value,
    function(x) {
      if (
        length(x) == 0 ||
        !is.finite(x)
      ) {
        return("P = NA")
      }
      
      paste0(
        "P ",
        ifelse(
          x < 0.001,
          "< 0.001",
          paste0(
            "= ",
            format.pval(
              x,
              digits = 2,
              eps = 0.001
            )
          )
        )
      )
    },
    character(1)
  )
}

make_diagonal_cell <- function(display_label) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "rect",
      xmin = 0,
      xmax = 1,
      ymin = 0,
      ymax = 1,
      fill = "white",
      color = "grey68",
      linewidth = 0.55
    ) +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = display_label,
      size = 4.8,
      fontface = "plain"
    ) +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.margin = ggplot2::margin(
        1.5,
        1.5,
        1.5,
        1.5
      )
    )
}

make_blank_cell <- function() {
  ggplot2::ggplot() +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.margin = ggplot2::margin(
        1.5,
        1.5,
        1.5,
        1.5
      )
    )
}

make_summary_profile <- function(pair_row) {
  if (nrow(pair_row) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "rect",
          xmin = 0,
          xmax = 1,
          ymin = 0,
          ymax = 1,
          fill = "grey96",
          color = NA
        ) +
        ggplot2::annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "Not available",
          size = 3.6
        ) +
        ggplot2::coord_fixed(
          ratio = 1,
          xlim = c(0, 1),
          ylim = c(0, 1),
          expand = FALSE
        ) +
        ggplot2::theme_void()
    )
  }
  
  r_value <- as.numeric(
    pair_row$protest_r[1]
  )
  r_ci_low <- as.numeric(
    pair_row$r_ci_low[1]
  )
  r_ci_high <- as.numeric(
    pair_row$r_ci_high[1]
  )
  p_value <- as.numeric(
    pair_row$protest_p[1]
  )
  
  evidence_supported <-
    is.finite(p_value) &&
    p_value < 0.05
  
  accent_color <- if (evidence_supported) {
    "#365F7D"
  } else {
    "grey58"
  }
  
  point_fill <- if (evidence_supported) {
    "#365F7D"
  } else {
    "white"
  }
  
  profile_start <- 0.16
  profile_end <- 0.84
  profile_width <- profile_end - profile_start
  
  scale_r <- function(value) {
    if (!is.finite(value)) {
      return(profile_start)
    }
    
    profile_start +
      profile_width *
      min(
        max(value, 0),
        1
      )
  }
  
  r_point_x <- scale_r(r_value)
  r_low_x <- scale_r(
    ifelse(
      is.finite(r_ci_low),
      r_ci_low,
      r_value
    )
  )
  r_high_x <- scale_r(
    ifelse(
      is.finite(r_ci_high),
      r_ci_high,
      r_value
    )
  )
  
  ci_label <- if (
    is.finite(r_ci_low) &&
    is.finite(r_ci_high)
  ) {
    paste0(
      "95% CI ",
      sprintf("%.2f", r_ci_low),
      "–",
      sprintf("%.2f", r_ci_high)
    )
  } else {
    "95% CI unavailable"
  }
  
  support_label <- paste0(
    format_matrix_p(
      pair_row$protest_p[1]
    ),
    ", n = ",
    pair_row$n[1]
  )
  
  ggplot2::ggplot() +
    ggplot2::annotate(
      "rect",
      xmin = 0,
      xmax = 1,
      ymin = 0,
      ymax = 1,
      fill = "grey96",
      color = NA
    ) +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.70,
      label = paste0(
        "r = ",
        sprintf("%.2f", r_value)
      ),
      size = 4.25,
      fontface = "bold"
    ) +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.56,
      label = ci_label,
      size = 3.25,
      color = "grey30"
    ) +
    ggplot2::annotate(
      "segment",
      x = profile_start,
      xend = profile_end,
      y = 0.36,
      yend = 0.36,
      linewidth = 0.55,
      color = "grey82"
    ) +
    ggplot2::annotate(
      "segment",
      x = r_low_x,
      xend = r_high_x,
      y = 0.36,
      yend = 0.36,
      linewidth = 1.65,
      color = accent_color
    ) +
    ggplot2::annotate(
      "segment",
      x = r_low_x,
      xend = r_low_x,
      y = 0.325,
      yend = 0.395,
      linewidth = 0.62,
      color = accent_color
    ) +
    ggplot2::annotate(
      "segment",
      x = r_high_x,
      xend = r_high_x,
      y = 0.325,
      yend = 0.395,
      linewidth = 0.62,
      color = accent_color
    ) +
    ggplot2::annotate(
      "point",
      x = r_point_x,
      y = 0.36,
      shape = 21,
      size = 3.35,
      fill = point_fill,
      color = "grey18",
      stroke = 0.65
    ) +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.15,
      label = support_label,
      size = 3.05,
      fontface = ifelse(
        evidence_supported,
        "bold",
        "plain"
      ),
      color = "grey20"
    ) +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.margin = ggplot2::margin(
        1.5,
        1.5,
        1.5,
        1.5
      )
    )
}

format_scatter_statistics <- function(pair_row) {
  if (nrow(pair_row) == 0) {
    return("'Statistics unavailable'")
  }
  
  ci_text <- if (
    is.finite(pair_row$r_ci_low[1]) &&
    is.finite(pair_row$r_ci_high[1])
  ) {
    sprintf(
      "r = %.2f [95%% CI %.2f-%.2f]",
      pair_row$protest_r[1],
      pair_row$r_ci_low[1],
      pair_row$r_ci_high[1]
    )
  } else {
    sprintf(
      "r = %.2f [95%% CI unavailable]",
      pair_row$protest_r[1]
    )
  }
  
  p_value <- as.numeric(pair_row$protest_p[1])
  n_value <- as.integer(pair_row$n[1])
  
  p_expression <- if (!is.finite(p_value)) {
    "italic(p) == NA"
  } else if (p_value < 0.001) {
    "italic(p) < 0.001"
  } else {
    paste0(
      "italic(p) == ",
      format(
        round(p_value, 3),
        nsmall = ifelse(p_value < 0.01, 3, 2),
        trim = TRUE,
        scientific = FALSE
      )
    )
  }
  
  paste0(
    "atop('",
    ci_text,
    "', ",
    p_expression,
    "*', n = ",
    n_value,
    "')"
  )
}

matrix_panels_full <- list()
matrix_panels_upper <- list()

for (row_index in seq_along(modality_order)) {
  for (col_index in seq_along(modality_order)) {
    row_block <- modality_order[row_index]
    col_block <- modality_order[col_index]
    panel_key <- paste0(
      row_block,
      "__",
      col_block
    )
    
    if (row_index == col_index) {
      diagonal_plot <- make_diagonal_cell(
        modality_labels[row_index]
      )
      
      matrix_panels_full[[panel_key]] <- diagonal_plot
      matrix_panels_upper[[panel_key]] <- diagonal_plot
      next
    }
    
    pair_match <- all_pair_specs[
      (all_pair_specs$block_x == row_block &
         all_pair_specs$block_y == col_block) |
        (all_pair_specs$block_x == col_block &
           all_pair_specs$block_y == row_block),
      ,
      drop = FALSE
    ]
    
    if (nrow(pair_match) != 1) {
      matrix_panels_full[[panel_key]] <- make_blank_cell()
      matrix_panels_upper[[panel_key]] <- make_blank_cell()
      next
    }
    
    pair_name_current <- pair_match$pair_name[1]
    pair_row <- pair_display_summary[
      pair_display_summary$pair_name == pair_name_current,
      ,
      drop = FALSE
    ]
    
    if (col_index > row_index) {
      pair_object <- coherence_get_pair(
        coherence_data,
        pair_name_current
      )
      
      X_all <- original_score_matrix(
        pair_object,
        row_block
      )
      
      Y_all <- original_score_matrix(
        pair_object,
        col_block
      )
      
      common_samples <- pair_object$sample_meta$SampleID[
        pair_object$sample_meta$SampleID %in%
          intersect(
            rownames(X_all),
            rownames(Y_all)
          )
      ]
      
      axes_used <- min(
        protest_max_axes,
        ncol(X_all),
        ncol(Y_all),
        length(common_samples) - 1
      )
      
      if (
        length(common_samples) < 4 ||
        axes_used < 2
      ) {
        na_plot <- ggplot2::ggplot() +
          ggplot2::annotate(
            "rect",
            xmin = 0,
            xmax = 1,
            ymin = 0,
            ymax = 1,
            fill = "grey97",
            color = "grey80",
            linewidth = 0.55
          ) +
          ggplot2::annotate(
            "text",
            x = 0.5,
            y = 0.5,
            label = "Not enough\nsamples",
            size = 3.7
          ) +
          ggplot2::coord_fixed(
            ratio = 1,
            xlim = c(0, 1),
            ylim = c(0, 1),
            expand = FALSE
          ) +
          ggplot2::theme_void() +
          ggplot2::theme(
            aspect.ratio = 1
          )
        
        matrix_panels_full[[panel_key]] <- na_plot
        matrix_panels_upper[[panel_key]] <- na_plot
        next
      }
      
      X <- X_all[
        common_samples,
        seq_len(axes_used),
        drop = FALSE
      ]
      Y <- Y_all[
        common_samples,
        seq_len(axes_used),
        drop = FALSE
      ]
      
      proc <- vegan::procrustes(
        X,
        Y,
        symmetric = TRUE
      )
      
      segment_data <- data.frame(
        x = proc$X[, 1],
        y = proc$X[, 2],
        xend = proc$Yrot[, 1],
        yend = proc$Yrot[, 2],
        stringsAsFactors = FALSE
      )
      
      point_data <- rbind(
        data.frame(
          Axis1 = proc$X[, 1],
          Axis2 = proc$X[, 2],
          Layer = unname(layer_labels[row_block]),
          stringsAsFactors = FALSE
        ),
        data.frame(
          Axis1 = proc$Yrot[, 1],
          Axis2 = proc$Yrot[, 2],
          Layer = unname(layer_labels[col_block]),
          stringsAsFactors = FALSE
        )
      )
      
      point_data$Layer <- factor(
        point_data$Layer,
        levels = unname(layer_labels)
      )
      
      x_range <- range(
        c(
          segment_data$x,
          segment_data$xend
        ),
        na.rm = TRUE
      )
      
      y_range <- range(
        c(
          segment_data$y,
          segment_data$yend
        ),
        na.rm = TRUE
      )
      
      x_center <- mean(x_range)
      y_center <- mean(y_range)
      half_span <- 0.5 * max(
        diff(x_range),
        diff(y_range)
      )
      
      if (
        !is.finite(half_span) ||
        half_span <= 0
      ) {
        half_span <- 1
      }
      
      half_span <- half_span * 1.08
      x_limits <- x_center + c(
        -half_span,
        half_span
      )
      y_limits <- y_center + c(
        -half_span,
        half_span
      )
      
      scatter_basic <- ggplot2::ggplot() +
        ggplot2::geom_segment(
          data = segment_data,
          ggplot2::aes(
            x = x,
            y = y,
            xend = xend,
            yend = yend
          ),
          color = "grey79",
          linewidth = 0.32
        ) +
        ggplot2::geom_point(
          data = point_data,
          ggplot2::aes(
            x = Axis1,
            y = Axis2,
            shape = Layer,
            fill = Layer
          ),
          color = "grey25",
          size = 2.55,
          stroke = 0.60
        ) +
        ggplot2::scale_shape_manual(
          values = layer_shapes,
          breaks = unname(layer_labels),
          drop = FALSE,
          name = "Data type"
        ) +
        ggplot2::scale_fill_manual(
          values = layer_colors,
          breaks = unname(layer_labels),
          drop = FALSE,
          name = "Data type"
        ) +
        ggplot2::coord_fixed(
          ratio = 1,
          xlim = x_limits,
          ylim = y_limits,
          expand = FALSE,
          clip = "off"
        ) +
        ggplot2::theme_classic(
          base_size = 9.8
        ) +
        ggplot2::theme(
          aspect.ratio = 1,
          axis.title = ggplot2::element_blank(),
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          legend.position = "none",
          plot.margin = ggplot2::margin(
            1.5,
            1.5,
            1.5,
            1.5
          )
        )
      
      scatter_annotated <- scatter_basic +
        ggplot2::annotate(
          "label",
          x = x_limits[1] + 0.045 * diff(x_limits),
          y = y_limits[2] - 0.045 * diff(y_limits),
          label = format_scatter_statistics(pair_row),
          parse = TRUE,
          hjust = 0,
          vjust = 1,
          size = 2.35,
          lineheight = 0.92,
          color = "grey15",
          fill = scales::alpha("white", 0.84),
          linewidth = 0.20,
          label.padding = grid::unit(0.10, "lines"),
          label.r = grid::unit(0.08, "lines")
        )
      
      matrix_panels_full[[panel_key]] <- scatter_basic
      matrix_panels_upper[[panel_key]] <- scatter_annotated
    } else {
      matrix_panels_full[[panel_key]] <- make_summary_profile(
        pair_row
      )
      matrix_panels_upper[[panel_key]] <- make_blank_cell()
    }
  }
}

matrix_panel_order <- c()
for (row_block in modality_order) {
  for (col_block in modality_order) {
    matrix_panel_order <- c(
      matrix_panel_order,
      paste0(
        row_block,
        "__",
        col_block
      )
    )
  }
}

p_procrustes_legend <- ggplot2::ggplot(
  data.frame(
    Axis1 = 0,
    Axis2 = 0,
    Layer = factor(
      unname(layer_labels),
      levels = unname(layer_labels)
    )
  ),
  ggplot2::aes(
    x = Axis1,
    y = Axis2,
    shape = Layer,
    fill = Layer
  )
) +
  ggplot2::geom_point(
    size = 3.3,
    color = "grey25",
    stroke = 0.60
  ) +
  ggplot2::scale_shape_manual(
    values = layer_shapes,
    breaks = unname(layer_labels),
    drop = FALSE,
    name = "Data type"
  ) +
  ggplot2::scale_fill_manual(
    values = layer_colors,
    breaks = unname(layer_labels),
    drop = FALSE,
    name = "Data type"
  ) +
  ggplot2::guides(
    shape = ggplot2::guide_legend(
      override.aes = list(
        size = 3.6,
        color = "grey25",
        stroke = 0.60
      )
    )
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position = "right",
    legend.title = ggplot2::element_text(
      size = 10.2,
      face = "bold"
    ),
    legend.text = ggplot2::element_text(
      size = 9.4
    ),
    legend.key.height = grid::unit(
      12,
      "pt"
    )
  )

p_procrustes_matrix_full_core <- patchwork::wrap_plots(
  matrix_panels_full[
    matrix_panel_order
  ],
  ncol = 4,
  nrow = 4,
  widths = rep(1, 4),
  heights = rep(1, 4)
)

p_procrustes_matrix_full <- patchwork::wrap_plots(
  patchwork::wrap_elements(
    full = p_procrustes_matrix_full_core
  ),
  patchwork::wrap_elements(
    full = p_procrustes_legend
  ),
  nrow = 1,
  widths = c(
    4,
    0.62
  )
)

ggplot2::ggsave(
  filename = paste0(
    "figures/coherence/",
    "Procrustes_4x4_data_type_colors_",
    "upper_scatter_lower_statistics_v17.svg"
  ),
  plot = p_procrustes_matrix_full,
  width = 11.4,
  height = 9.8,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

p_procrustes_matrix_upper_core <- patchwork::wrap_plots(
  matrix_panels_upper[
    matrix_panel_order
  ],
  ncol = 4,
  nrow = 4,
  widths = rep(1, 4),
  heights = rep(1, 4)
)

p_procrustes_matrix_upper <- patchwork::wrap_plots(
  patchwork::wrap_elements(
    full = p_procrustes_matrix_upper_core
  ),
  patchwork::wrap_elements(
    full = p_procrustes_legend
  ),
  nrow = 1,
  widths = c(
    4,
    0.62
  )
)

ggplot2::ggsave(
  filename = paste0(
    "figures/coherence/",
    "Procrustes_4x4_data_type_colors_",
    "upper_triangle_inpanel_statistics_v17.svg"
  ),
  plot = p_procrustes_matrix_upper,
  width = 11.4,
  height = 9.8,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)


#=================================================================#
# 3. Compact forest-style PROTEST-r summary
#=================================================================#

forest_summary <- pair_display_summary[
  ,
  c(
    "pair_name",
    "pair",
    "n",
    "protest_r",
    "r_ci_low",
    "r_ci_high",
    "protest_p"
  ),
  drop = FALSE
]

forest_summary$pair <- factor(
  as.character(
    forest_summary$pair
  ),
  levels = rev(
    all_pair_specs$pair_label
  )
)

forest_summary$evidence <- ifelse(
  is.finite(
    forest_summary$protest_p
  ) &
    forest_summary$protest_p < 0.05,
  "P < 0.05",
  "P ≥ 0.05"
)

forest_summary$annotation <- paste0(
  format_matrix_p(
    forest_summary$protest_p
  ),
  ", n = ",
  forest_summary$n
)

p_procrustes_forest <- ggplot2::ggplot(
  forest_summary,
  ggplot2::aes(
    x = protest_r,
    y = pair,
    color = evidence,
    fill = evidence
  )
) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin = r_ci_low,
      xmax = r_ci_high
    ),
    height = 0.12,
    linewidth = 0.72,
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    shape = 21,
    size = 3.0,
    stroke = 0.70,
    na.rm = TRUE
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      x = 1.02,
      label = annotation
    ),
    hjust = 0,
    size = 3.15,
    color = "grey20",
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(
    values = c(
      `P < 0.05` = "#365F7D",
      `P ≥ 0.05` = "grey58"
    ),
    guide = "none"
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      `P < 0.05` = "#365F7D",
      `P ≥ 0.05` = "white"
    ),
    guide = "none"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      1.17
    ),
    breaks = seq(
      0,
      1,
      by = 0.2
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0)
    ),
    name = "PROTEST correlation-like statistic (r)"
  ) +
  ggplot2::scale_y_discrete(
    expand = ggplot2::expansion(
      add = c(0.22, 0.22)
    )
  ) +
  ggplot2::labs(
    y = NULL,
    caption = paste0(
      "Intervals: observed-r-centered 95% CIs from SubjectID-cluster bootstrap SEs; ",
      "P values: 9,999 PROTEST permutations."
    )
  ) +
  ggplot2::theme_classic(
    base_size = 10.5
  ) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(
      face = "bold",
      size = 9.2
    ),
    axis.text.x = ggplot2::element_text(
      size = 8.8
    ),
    axis.title.x = ggplot2::element_text(
      size = 9.8,
      margin = ggplot2::margin(
        t = 4
      )
    ),
    axis.ticks.y = ggplot2::element_blank(),
    plot.caption = ggplot2::element_text(
      hjust = 0,
      size = 7.6,
      margin = ggplot2::margin(
        t = 4
      )
    ),
    plot.margin = ggplot2::margin(
      3,
      30,
      3,
      3
    )
  ) +
  ggplot2::coord_cartesian(
    clip = "off"
  )

ggplot2::ggsave(
  filename =
    "figures/coherence/Procrustes_PROTEST_r_compact_forest_v17.svg",
  plot = p_procrustes_forest,
  width = 7.5,
  height = 3.55,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)



#=================================================================#
# 4. Four-block generalized Procrustes analysis
#=================================================================#
#
# This section extends the pairwise Procrustes analysis to all four omics
# configurations simultaneously. It is a generalized Procrustes analysis
# (GPA) of modality-specific PCoA configurations, not a direct concatenation
# of the original feature tables.
#
# All four blocks use the same four-way common samples and the same number of
# ordination axes. Each configuration is centered and normalized, iteratively
# rotated/reflected toward a common consensus, and finally referred to the
# principal axes of that consensus.
#
# The primary visualization is a consensus-spoke plot:
#   - small grey cross: cross-omics consensus position of one sample
#   - colored point: modality-specific position after GPA alignment
#   - spoke length: modality-specific disagreement from the consensus
#
# Global inference:
#   - statistic: total GPA residual sum of squares; smaller values indicate
#     stronger four-way correspondence
#   - permutation: Species is fixed; the sample labels of the other three
#     blocks are independently permuted
#   - when SubjectID and Timepoint are available, whole subjects are permuted
#     only among subjects with the same observed Timepoint pattern
#
# The current input stores pair-specific PCoA scores. For each modality, this
# code automatically selects the available pair object containing all four-way
# samples and the largest number of retained axes. A future input-processing
# version may instead store a dedicated four-way-common-sample PCoA for each
# modality; that would remove this source-selection step.
#
#-----------------------------------------------------------------#

fit_gpa_configurations <- function(
    configurations,
    tolerance = 1e-10,
    max_iterations = 500,
    calculate_pairwise = TRUE
) {
  if (length(configurations) < 3) {
    stop(
      "Generalized Procrustes analysis requires at least three configurations.",
      call. = FALSE
    )
  }
  
  if (
    length(
      unique(
        vapply(
          configurations,
          nrow,
          integer(1)
        )
      )
    ) != 1 ||
    length(
      unique(
        vapply(
          configurations,
          ncol,
          integer(1)
        )
      )
    ) != 1
  ) {
    stop(
      "All GPA configurations must have identical row and column dimensions.",
      call. = FALSE
    )
  }
  
  standardized <- lapply(
    configurations,
    function(x) {
      x <- as.matrix(x)
      storage.mode(x) <- "numeric"
      
      x <- sweep(
        x,
        2,
        colMeans(
          x,
          na.rm = TRUE
        ),
        FUN = "-"
      )
      
      configuration_norm <- sqrt(
        sum(
          x^2,
          na.rm = TRUE
        )
      )
      
      if (
        !is.finite(configuration_norm) ||
        configuration_norm <= 0
      ) {
        stop(
          "A GPA configuration has zero or non-finite total variation.",
          call. = FALSE
        )
      }
      
      x / configuration_norm
    }
  )
  
  consensus <- Reduce(
    "+",
    standardized
  ) / length(
    standardized
  )
  
  consensus <- sweep(
    consensus,
    2,
    colMeans(
      consensus,
      na.rm = TRUE
    ),
    FUN = "-"
  )
  
  consensus_norm <- sqrt(
    sum(
      consensus^2,
      na.rm = TRUE
    )
  )
  
  if (
    !is.finite(consensus_norm) ||
    consensus_norm <= 0
  ) {
    stop(
      "The initial GPA consensus has zero or non-finite total variation.",
      call. = FALSE
    )
  }
  
  consensus <- consensus / consensus_norm
  previous_loss <- Inf
  converged <- FALSE
  iteration_used <- max_iterations
  
  for (iteration_index in seq_len(max_iterations)) {
    aligned <- lapply(
      standardized,
      function(x) {
        rotation_svd <- svd(
          crossprod(
            x,
            consensus
          )
        )
        
        x %*%
          rotation_svd$u %*%
          t(
            rotation_svd$v
          )
      }
    )
    
    updated_consensus <- Reduce(
      "+",
      aligned
    ) / length(
      aligned
    )
    
    updated_consensus <- sweep(
      updated_consensus,
      2,
      colMeans(
        updated_consensus,
        na.rm = TRUE
      ),
      FUN = "-"
    )
    
    updated_norm <- sqrt(
      sum(
        updated_consensus^2,
        na.rm = TRUE
      )
    )
    
    if (
      !is.finite(updated_norm) ||
      updated_norm <= 0
    ) {
      stop(
        "The GPA consensus became degenerate during iteration.",
        call. = FALSE
      )
    }
    
    updated_consensus <-
      updated_consensus /
      updated_norm
    
    current_loss <- sum(
      vapply(
        aligned,
        function(x) {
          sum(
            (
              x -
                updated_consensus
            )^2,
            na.rm = TRUE
          )
        },
        numeric(1)
      )
    )
    
    if (
      is.finite(previous_loss) &&
      abs(
        previous_loss -
        current_loss
      ) < tolerance
    ) {
      consensus <- updated_consensus
      converged <- TRUE
      iteration_used <- iteration_index
      break
    }
    
    consensus <- updated_consensus
    previous_loss <- current_loss
  }
  
  aligned <- lapply(
    standardized,
    function(x) {
      rotation_svd <- svd(
        crossprod(
          x,
          consensus
        )
      )
      
      x %*%
        rotation_svd$u %*%
        t(
          rotation_svd$v
        )
    }
  )
  
  principal_rotation <- svd(
    consensus,
    nu = 0,
    nv = ncol(
      consensus
    )
  )$v
  
  consensus <-
    consensus %*%
    principal_rotation
  
  aligned <- lapply(
    aligned,
    function(x) {
      x %*%
        principal_rotation
    }
  )
  
  residual_ss <- vapply(
    aligned,
    function(x) {
      sum(
        (
          x -
            consensus
        )^2,
        na.rm = TRUE
      )
    },
    numeric(1)
  )
  
  consensus_similarity <- vapply(
    aligned,
    function(x) {
      sum(
        x *
          consensus,
        na.rm = TRUE
      ) /
        sqrt(
          sum(
            x^2,
            na.rm = TRUE
          ) *
            sum(
              consensus^2,
              na.rm = TRUE
            )
        )
    },
    numeric(1)
  )
  
  consensus_similarity <- pmax(
    -1,
    pmin(
      1,
      consensus_similarity
    )
  )
  
  pairwise_r <- matrix(
    NA_real_,
    nrow = length(
      aligned
    ),
    ncol = length(
      aligned
    ),
    dimnames = list(
      names(
        aligned
      ),
      names(
        aligned
      )
    )
  )
  
  diag(
    pairwise_r
  ) <- 1
  
  if (calculate_pairwise) {
    for (
      first_index in
      seq_len(
        length(
          aligned
        ) - 1
      )
    ) {
      for (
        second_index in
        seq.int(
          first_index + 1,
          length(
            aligned
          )
        )
      ) {
        pairwise_proc <- vegan::procrustes(
          aligned[[first_index]],
          aligned[[second_index]],
          symmetric = TRUE
        )
        
        pairwise_r[
          first_index,
          second_index
        ] <- sqrt(
          max(
            0,
            1 -
              as.numeric(
                pairwise_proc$ss
              )
          )
        )
        
        pairwise_r[
          second_index,
          first_index
        ] <- pairwise_r[
          first_index,
          second_index
        ]
      }
    }
  }
  
  consensus_axis_fraction <- colSums(
    consensus^2,
    na.rm = TRUE
  )
  
  consensus_axis_fraction <-
    consensus_axis_fraction /
    sum(
      consensus_axis_fraction,
      na.rm = TRUE
    )
  
  list(
    aligned = aligned,
    consensus = consensus,
    residual_ss = residual_ss,
    total_residual_ss = sum(
      residual_ss,
      na.rm = TRUE
    ),
    consensus_similarity =
      consensus_similarity,
    overall_agreement = mean(
      consensus_similarity,
      na.rm = TRUE
    ),
    pairwise_r = pairwise_r,
    mean_pairwise_r = if (
      calculate_pairwise
    ) {
      mean(
        pairwise_r[
          upper.tri(
            pairwise_r
          )
        ],
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    consensus_axis_fraction =
      consensus_axis_fraction,
    converged = converged,
    iterations = iteration_used
  )
}


make_gpa_permutation_index <- function(
    sample_meta,
    use_subject_pattern
) {
  if (use_subject_pattern) {
    subject_id <- as.character(
      sample_meta$SubjectID
    )
    timepoint <- as.character(
      sample_meta$Timepoint
    )
    
    subject_pattern <- vapply(
      split(
        timepoint,
        subject_id
      ),
      function(x) {
        paste(
          sort(
            unique(
              x
            )
          ),
          collapse = "|"
        )
      },
      character(1)
    )
    
    subject_mapping <- setNames(
      names(
        subject_pattern
      ),
      names(
        subject_pattern
      )
    )
    
    for (
      pattern_name in
      unique(
        subject_pattern
      )
    ) {
      pattern_subjects <- names(
        subject_pattern
      )[
        subject_pattern ==
          pattern_name
      ]
      
      subject_mapping[
        pattern_subjects
      ] <- sample(
        pattern_subjects,
        length(
          pattern_subjects
        ),
        replace = FALSE
      )
    }
    
    source_key <- paste(
      subject_id,
      timepoint,
      sep = "|||"
    )
    
    target_key <- paste(
      unname(
        subject_mapping[
          subject_id
        ]
      ),
      timepoint,
      sep = "|||"
    )
    
    permutation_index <- match(
      target_key,
      source_key
    )
    
    if (
      all(
        is.finite(
          permutation_index
        )
      ) &&
      length(
        unique(
          permutation_index
        )
      ) ==
      nrow(
        sample_meta
      )
    ) {
      return(
        permutation_index
      )
    }
  }
  
  sample.int(
    nrow(
      sample_meta
    )
  )
}


gpa_pair_candidates <- list(
  species = c(
    "species_ko",
    "species_metabolite",
    "species_host"
  ),
  ko = c(
    "species_ko",
    "ko_metabolite",
    "ko_host"
  ),
  metabolite = c(
    "species_metabolite",
    "ko_metabolite",
    "metabolite_host"
  ),
  host = c(
    "species_host",
    "ko_host",
    "metabolite_host"
  )
)

gpa_common_samples <- Reduce(
  intersect,
  lapply(
    coherence_data$modalities[
      modality_order
    ],
    function(x) {
      as.character(
        x$sample_meta$SampleID
      )
    }
  )
)

gpa_common_samples <- as.character(
  coherence_data$modalities$species$sample_meta$SampleID[
    coherence_data$modalities$species$sample_meta$SampleID %in%
      gpa_common_samples
  ]
)

if (length(gpa_common_samples) < 4) {
  stop(
    "Fewer than four samples contain all four omics data types.",
    call. = FALSE
  )
}

gpa_configuration_sources <- data.frame()
gpa_score_candidates <- list()

for (block_name in modality_order) {
  candidate_summary <- data.frame()
  candidate_scores <- list()
  
  for (
    pair_name_current in
    gpa_pair_candidates[[block_name]]
  ) {
    pair_object <- coherence_get_pair(
      coherence_data,
      pair_name_current
    )
    
    score_current <- tryCatch(
      original_score_matrix(
        pair_object,
        block_name
      ),
      error = function(e) {
        NULL
      }
    )
    
    if (
      is.null(
        score_current
      ) ||
      !all(
        gpa_common_samples %in%
        rownames(
          score_current
        )
      )
    ) {
      next
    }
    
    candidate_scores[[pair_name_current]] <-
      score_current
    
    candidate_summary <- rbind(
      candidate_summary,
      data.frame(
        block = block_name,
        pair_name = pair_name_current,
        n_score_samples = nrow(
          score_current
        ),
        n_score_axes = ncol(
          score_current
        ),
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (nrow(candidate_summary) == 0) {
    stop(
      paste0(
        "No pair-specific score matrix for ",
        block_name,
        " contains every four-way common sample."
      ),
      call. = FALSE
    )
  }
  
  selected_index <- order(
    -candidate_summary$n_score_axes,
    -candidate_summary$n_score_samples,
    candidate_summary$pair_name
  )[1]
  
  selected_pair <-
    candidate_summary$pair_name[
      selected_index
    ]
  
  gpa_score_candidates[[block_name]] <-
    candidate_scores[[selected_pair]]
  
  gpa_configuration_sources <- rbind(
    gpa_configuration_sources,
    candidate_summary[
      selected_index,
      ,
      drop = FALSE
    ]
  )
}

gpa_axes_used <- min(
  gpa_max_axes,
  length(
    gpa_common_samples
  ) - 1,
  vapply(
    gpa_score_candidates,
    ncol,
    integer(1)
  )
)

if (gpa_axes_used < 2) {
  stop(
    "Fewer than two common ordination axes are available for four-block GPA.",
    call. = FALSE
  )
}

gpa_configurations <- lapply(
  gpa_score_candidates,
  function(x) {
    x[
      gpa_common_samples,
      seq_len(
        gpa_axes_used
      ),
      drop = FALSE
    ]
  }
)

names(
  gpa_configurations
) <- modality_order

gpa_sample_meta <- coherence_data$modalities$species$sample_meta[
  match(
    gpa_common_samples,
    coherence_data$modalities$species$sample_meta$SampleID
  ),
  ,
  drop = FALSE
]

gpa_use_subject_pattern_permutation <-
  all(
    c(
      "SubjectID",
      "Timepoint"
    ) %in%
      colnames(
        gpa_sample_meta
      )
  ) &&
  all(
    stats::complete.cases(
      gpa_sample_meta[
        ,
        c(
          "SubjectID",
          "Timepoint"
        ),
        drop = FALSE
      ]
    )
  ) &&
  !anyDuplicated(
    paste(
      gpa_sample_meta$SubjectID,
      gpa_sample_meta$Timepoint,
      sep = "|||"
    )
  )

gpa_permutation_method <- if (
  gpa_use_subject_pattern_permutation
) {
  paste0(
    "Independent SubjectID permutations within identical ",
    "Timepoint-availability patterns"
  )
} else {
  "Independent unrestricted sample-label permutations"
}

gpa_4omics_result <- fit_gpa_configurations(
  gpa_configurations,
  tolerance = gpa_tolerance,
  max_iterations = gpa_max_iterations,
  calculate_pairwise = TRUE
)

set.seed(
  20260728
)

gpa_permuted_residual_ss <- rep(
  NA_real_,
  gpa_permutations
)

for (
  permutation_index in
  seq_len(
    gpa_permutations
  )
) {
  permuted_configurations <-
    gpa_configurations
  
  for (
    block_index in
    2:length(
      permuted_configurations
    )
  ) {
    permuted_configurations[[block_index]] <- permuted_configurations[[block_index]][
      make_gpa_permutation_index(
        gpa_sample_meta,
        gpa_use_subject_pattern_permutation
      ),
      ,
      drop = FALSE
    ]
  }
  
  gpa_permuted_residual_ss[
    permutation_index
  ] <- fit_gpa_configurations(
    permuted_configurations,
    tolerance = gpa_tolerance,
    max_iterations = gpa_max_iterations,
    calculate_pairwise = FALSE
  )$total_residual_ss
}

gpa_global_p <- (
  1 +
    sum(
      gpa_permuted_residual_ss <=
        gpa_4omics_result$total_residual_ss,
      na.rm = TRUE
    )
) / (
  1 +
    sum(
      is.finite(
        gpa_permuted_residual_ss
      )
    )
)

gpa_global_summary <- data.frame(
  n_samples = length(
    gpa_common_samples
  ),
  axes_used = gpa_axes_used,
  n_blocks = length(
    gpa_configurations
  ),
  total_residual_ss =
    gpa_4omics_result$total_residual_ss,
  overall_agreement =
    gpa_4omics_result$overall_agreement,
  mean_pairwise_r =
    gpa_4omics_result$mean_pairwise_r,
  permutation_p =
    gpa_global_p,
  permutations =
    gpa_permutations,
  permutation_method =
    gpa_permutation_method,
  converged =
    gpa_4omics_result$converged,
  iterations =
    gpa_4omics_result$iterations,
  stringsAsFactors = FALSE
)

gpa_block_summary <- data.frame(
  block = modality_order,
  data_type = unname(
    layer_labels[
      modality_order
    ]
  ),
  residual_ss = unname(
    gpa_4omics_result$residual_ss[
      modality_order
    ]
  ),
  consensus_similarity = unname(
    gpa_4omics_result$consensus_similarity[
      modality_order
    ]
  ),
  score_source_pair = gpa_configuration_sources$pair_name[
    match(
      modality_order,
      gpa_configuration_sources$block
    )
  ],
  stringsAsFactors = FALSE
)

gpa_block_summary$data_type <- factor(
  gpa_block_summary$data_type,
  levels = rev(
    unname(
      layer_labels[
        modality_order
      ]
    )
  )
)

gpa_pairwise_summary <- data.frame()

for (
  first_index in
  seq_len(
    length(
      modality_order
    ) - 1
  )
) {
  for (
    second_index in
    seq.int(
      first_index + 1,
      length(
        modality_order
      )
    )
  ) {
    gpa_pairwise_summary <- rbind(
      gpa_pairwise_summary,
      data.frame(
        block_1 =
          modality_order[
            first_index
          ],
        block_2 =
          modality_order[
            second_index
          ],
        data_type_1 =
          unname(
            layer_labels[
              modality_order[
                first_index
              ]
            ]
          ),
        data_type_2 =
          unname(
            layer_labels[
              modality_order[
                second_index
              ]
            ]
          ),
        procrustes_r =
          gpa_4omics_result$pairwise_r[
            first_index,
            second_index
          ],
        stringsAsFactors = FALSE
      )
    )
  }
}

gpa_consensus_plot_data <- data.frame(
  SampleID = gpa_common_samples,
  Consensus1 =
    gpa_4omics_result$consensus[
      ,
      1
    ],
  Consensus2 =
    gpa_4omics_result$consensus[
      ,
      2
    ],
  stringsAsFactors = FALSE
)

gpa_point_plot_data <- do.call(
  rbind,
  lapply(
    modality_order,
    function(block_name) {
      data.frame(
        SampleID = gpa_common_samples,
        Axis1 =
          gpa_4omics_result$aligned[[block_name]][
            ,
            1
          ],
        Axis2 =
          gpa_4omics_result$aligned[[block_name]][
            ,
            2
          ],
        Layer =
          unname(
            layer_labels[
              block_name
            ]
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)

gpa_point_plot_data$Layer <- factor(
  gpa_point_plot_data$Layer,
  levels = unname(
    layer_labels[
      modality_order
    ]
  )
)

gpa_segment_plot_data <- merge(
  gpa_point_plot_data,
  gpa_consensus_plot_data,
  by = "SampleID",
  all.x = TRUE,
  sort = FALSE
)

gpa_residual_distance_data <- do.call(
  rbind,
  lapply(
    modality_order,
    function(block_name) {
      data.frame(
        SampleID = gpa_common_samples,
        SubjectID = if (
          "SubjectID" %in% colnames(gpa_sample_meta)
        ) {
          as.character(gpa_sample_meta$SubjectID)
        } else {
          NA_character_
        },
        Timepoint = if (
          "Timepoint" %in% colnames(gpa_sample_meta)
        ) {
          as.character(gpa_sample_meta$Timepoint)
        } else {
          NA_character_
        },
        block = block_name,
        data_type = unname(layer_labels[block_name]),
        residual_distance = sqrt(
          rowSums(
            (
              gpa_4omics_result$aligned[[block_name]] -
                gpa_4omics_result$consensus
            )^2,
            na.rm = TRUE
          )
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

gpa_residual_distance_data$data_type <- factor(
  gpa_residual_distance_data$data_type,
  levels = unname(
    layer_labels[
      modality_order
    ]
  )
)

gpa_residual_box_summary <- do.call(
  rbind,
  lapply(
    levels(
      gpa_residual_distance_data$data_type
    ),
    function(data_type_current) {
      residual_current <- gpa_residual_distance_data$residual_distance[
        gpa_residual_distance_data$data_type == data_type_current
      ]
      
      data.frame(
        data_type = data_type_current,
        q1 = unname(
          stats::quantile(
            residual_current,
            0.25,
            na.rm = TRUE
          )
        ),
        median = stats::median(
          residual_current,
          na.rm = TRUE
        ),
        q3 = unname(
          stats::quantile(
            residual_current,
            0.75,
            na.rm = TRUE
          )
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

gpa_residual_box_summary$data_type <- factor(
  gpa_residual_box_summary$data_type,
  levels = levels(
    gpa_residual_distance_data$data_type
  )
)

# Bootstrap uncertainty for modality-to-consensus agreement. The consensus is
# re-estimated in every bootstrap replicate. SubjectID is used as the cluster
# whenever available so that repeated sample-timepoint observations remain
# together; otherwise, SampleID is used as the resampling unit.
gpa_bootstrap_uses_subject <-
  "SubjectID" %in%
  colnames(gpa_sample_meta) &&
  all(
    !is.na(
      gpa_sample_meta$SubjectID
    )
  )

gpa_bootstrap_cluster <- if (
  gpa_bootstrap_uses_subject
) {
  as.character(
    gpa_sample_meta$SubjectID
  )
} else {
  as.character(
    gpa_sample_meta$SampleID
  )
}

gpa_bootstrap_unit <- if (
  gpa_bootstrap_uses_subject
) {
  "SubjectID cluster"
} else {
  "SampleID"
}

gpa_unique_bootstrap_clusters <- unique(
  gpa_bootstrap_cluster
)

gpa_agreement_bootstrap_values <- matrix(
  NA_real_,
  nrow = gpa_agreement_bootstraps,
  ncol = length(modality_order),
  dimnames = list(
    NULL,
    modality_order
  )
)

set.seed(
  20260729
)

for (
  bootstrap_index in
  seq_len(
    gpa_agreement_bootstraps
  )
) {
  sampled_clusters <- sample(
    gpa_unique_bootstrap_clusters,
    size = length(
      gpa_unique_bootstrap_clusters
    ),
    replace = TRUE
  )
  
  bootstrap_sample_index <- unlist(
    lapply(
      sampled_clusters,
      function(cluster_id) {
        which(
          gpa_bootstrap_cluster ==
            cluster_id
        )
      }
    ),
    use.names = FALSE
  )
  
  bootstrap_fit <- tryCatch(
    fit_gpa_configurations(
      lapply(
        gpa_configurations,
        function(x) {
          x[
            bootstrap_sample_index,
            ,
            drop = FALSE
          ]
        }
      ),
      tolerance = gpa_tolerance,
      max_iterations = gpa_max_iterations,
      calculate_pairwise = FALSE
    ),
    error = function(e) {
      NULL
    }
  )
  
  if (is.null(bootstrap_fit)) {
    next
  }
  
  gpa_agreement_bootstrap_values[
    bootstrap_index,
    modality_order
  ] <- bootstrap_fit$consensus_similarity[
    modality_order
  ]
}

gpa_agreement_bootstrap_summary <- do.call(
  rbind,
  lapply(
    modality_order,
    function(block_name) {
      bootstrap_values <- gpa_agreement_bootstrap_values[
        ,
        block_name
      ]
      
      bootstrap_values <- bootstrap_values[
        is.finite(
          bootstrap_values
        )
      ]
      
      data.frame(
        block = block_name,
        data_type = unname(
          layer_labels[
            block_name
          ]
        ),
        agreement = unname(
          gpa_4omics_result$consensus_similarity[
            block_name
          ]
        ),
        ci_low = if (
          length(bootstrap_values) >= 50
        ) {
          unname(
            stats::quantile(
              bootstrap_values,
              0.025,
              na.rm = TRUE
            )
          )
        } else {
          NA_real_
        },
        ci_high = if (
          length(bootstrap_values) >= 50
        ) {
          unname(
            stats::quantile(
              bootstrap_values,
              0.975,
              na.rm = TRUE
            )
          )
        } else {
          NA_real_
        },
        bootstrap_replicates = length(
          bootstrap_values
        ),
        bootstrap_unit = gpa_bootstrap_unit,
        stringsAsFactors = FALSE
      )
    }
  )
)

gpa_agreement_bootstrap_summary$data_type <- factor(
  gpa_agreement_bootstrap_summary$data_type,
  levels = rev(
    unname(
      layer_labels[
        modality_order
      ]
    )
  )
)

gpa_axis_1_fraction <-
  gpa_4omics_result$consensus_axis_fraction[
    1
  ]

gpa_axis_2_fraction <-
  gpa_4omics_result$consensus_axis_fraction[
    2
  ]

gpa_p_label <- if (
  is.finite(
    gpa_global_p
  ) &&
  gpa_global_p < 0.001
) {
  "P < 0.001"
} else {
  paste0(
    "P = ",
    format.pval(
      gpa_global_p,
      digits = 2,
      eps = 0.001
    )
  )
}

gpa_global_summary_text <- paste0(
  "Global four-omics fit\n",
  "n = ",
  length(
    gpa_common_samples
  ),
  ", axes = ",
  gpa_axes_used,
  "\nAgreement = ",
  sprintf(
    "%.2f",
    gpa_4omics_result$overall_agreement
  ),
  "\nMean all-pair r = ",
  sprintf(
    "%.2f",
    gpa_4omics_result$mean_pairwise_r
  ),
  "\nTotal residual SS = ",
  sprintf(
    "%.3f",
    gpa_4omics_result$total_residual_ss
  ),
  "\nPermutation ",
  gpa_p_label
)

gpa_x_range <- range(
  c(
    gpa_segment_plot_data$Axis1,
    gpa_segment_plot_data$Consensus1
  ),
  na.rm = TRUE
)

gpa_y_range <- range(
  c(
    gpa_segment_plot_data$Axis2,
    gpa_segment_plot_data$Consensus2
  ),
  na.rm = TRUE
)

gpa_x_span <- diff(
  gpa_x_range
)

gpa_y_span <- diff(
  gpa_y_range
)

if (
  !is.finite(
    gpa_x_span
  ) ||
  gpa_x_span <= 0
) {
  gpa_x_span <- 1
}

if (
  !is.finite(
    gpa_y_span
  ) ||
  gpa_y_span <= 0
) {
  gpa_y_span <- 1
}

p_gpa_4omics_consensus <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = gpa_segment_plot_data,
    ggplot2::aes(
      x = Consensus1,
      y = Consensus2,
      xend = Axis1,
      yend = Axis2,
      color = Layer
    ),
    linewidth = 0.34,
    alpha = 0.48
  ) +
  ggplot2::geom_point(
    data = gpa_consensus_plot_data,
    ggplot2::aes(
      x = Consensus1,
      y = Consensus2
    ),
    shape = 4,
    size = 1.35,
    stroke = 0.45,
    color = "grey40"
  ) +
  ggplot2::geom_point(
    data = gpa_point_plot_data,
    ggplot2::aes(
      x = Axis1,
      y = Axis2,
      shape = Layer,
      fill = Layer
    ),
    size = 2.75,
    color = "grey25",
    stroke = 0.62
  ) +
  ggplot2::annotate(
    "label",
    x = gpa_x_range[2] - 0.03 * gpa_x_span,
    y = gpa_y_range[1] + 0.07 * gpa_y_span,
    label = gpa_global_summary_text,
    hjust = 1,
    vjust = 0,
    size = 2.95,
    lineheight = 1.02,
    fill = scales::alpha(
      "white",
      0.88
    ),
    color = "grey15",
    linewidth = 0.24,
    label.padding = grid::unit(
      0.13,
      "lines"
    ),
    label.r = grid::unit(
      0.08,
      "lines"
    )
  ) +
  ggplot2::scale_shape_manual(
    values = layer_shapes,
    breaks = unname(
      layer_labels[
        modality_order
      ]
    ),
    drop = FALSE,
    name = "Data type"
  ) +
  ggplot2::scale_fill_manual(
    values = layer_colors,
    breaks = unname(
      layer_labels[
        modality_order
      ]
    ),
    drop = FALSE,
    name = "Data type"
  ) +
  ggplot2::scale_color_manual(
    values = layer_colors,
    breaks = unname(
      layer_labels[
        modality_order
      ]
    ),
    drop = FALSE,
    guide = "none"
  ) +
  ggplot2::labs(
    title =
      "Four-omics generalized Procrustes consensus",
    subtitle = paste0(
      "Each spoke links one data-type-specific sample position ",
      "to its cross-omics consensus."
    ),
    x = paste0(
      "Consensus axis 1 (",
      scales::percent(
        gpa_axis_1_fraction,
        accuracy = 0.1
      ),
      ")"
    ),
    y = paste0(
      "Consensus axis 2 (",
      scales::percent(
        gpa_axis_2_fraction,
        accuracy = 0.1
      ),
      ")"
    )
  ) +
  ggplot2::coord_fixed(
    ratio = 1,
    expand = TRUE,
    clip = "off"
  ) +
  ggplot2::theme_classic(
    base_size = 10.5
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 11.8
    ),
    plot.subtitle = ggplot2::element_text(
      size = 9.3,
      color = "grey35"
    ),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(
      face = "bold"
    ),
    legend.box = "horizontal",
    plot.margin = ggplot2::margin(
      5,
      7,
      4,
      5
    )
  )

p_gpa_pair_species_ko <-
  matrix_panels_upper[[
    "species__ko"
  ]] +
  ggplot2::labs(
    title = "Species-KEGG ortholog"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "plain",
      hjust = 0.5,
      size = 9.1,
      margin = ggplot2::margin(
        b = 2
      )
    )
  )

p_gpa_pair_ko_metabolite <-
  matrix_panels_upper[[
    "ko__metabolite"
  ]] +
  ggplot2::labs(
    title = "KEGG ortholog-Metabolite"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "plain",
      hjust = 0.5,
      size = 9.1,
      margin = ggplot2::margin(
        b = 2
      )
    )
  )

p_gpa_pair_metabolite_host <-
  matrix_panels_upper[[
    "metabolite__host"
  ]] +
  ggplot2::labs(
    title = "Metabolite-Host RNA-seq"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "plain",
      hjust = 0.5,
      size = 9.1,
      margin = ggplot2::margin(
        b = 2
      )
    )
  )

p_gpa_global_summary_box <- ggplot2::ggplot() +
  ggplot2::annotate(
    "rect",
    xmin = 0,
    xmax = 1,
    ymin = 0,
    ymax = 1,
    fill = "grey96",
    color = "grey78",
    linewidth = 0.55
  ) +
  ggplot2::annotate(
    "text",
    x = 0.06,
    y = 0.90,
    label = gpa_global_summary_text,
    hjust = 0,
    vjust = 1,
    size = 3.25,
    lineheight = 1.03,
    color = "grey15"
  ) +
  ggplot2::coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE,
    clip = "off"
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    plot.margin = ggplot2::margin(
      2,
      2,
      2,
      2
    )
  )

p_gpa_pairwise_chain_column <- patchwork::wrap_plots(
  p_gpa_pair_species_ko,
  p_gpa_pair_ko_metabolite,
  p_gpa_pair_metabolite_host,
  ncol = 1,
  heights = c(
    1,
    1,
    1
  )
)

p_gpa_4omics_chain_overview <- patchwork::wrap_plots(
  p_gpa_4omics_consensus,
  p_gpa_pairwise_chain_column,
  nrow = 1,
  widths = c(
    1.95,
    1
  )
)

p_gpa_residual_distribution <- ggplot2::ggplot(
  gpa_residual_distance_data,
  ggplot2::aes(
    x = data_type,
    y = residual_distance
  )
) +
  ggplot2::geom_violin(
    ggplot2::aes(
      fill = data_type
    ),
    trim = FALSE,
    scale = "width",
    width = 0.62,
    alpha = 0.24,
    color = NA
  ) +
  ggplot2::geom_crossbar(
    data = gpa_residual_box_summary,
    ggplot2::aes(
      x = data_type,
      y = median,
      ymin = q1,
      ymax = q3,
      fill = data_type,
      color = data_type
    ),
    width = 0.15,
    linewidth = 0.72,
    alpha = 0.30,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  ggbeeswarm::geom_quasirandom(
    ggplot2::aes(
      fill = data_type,
      color = data_type
    ),
    shape = 21,
    width = 0.09,
    size = 1.95,
    alpha = 0.88,
    stroke = 0.42,
    show.legend = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = layer_colors,
    breaks = unname(
      layer_labels[
        modality_order
      ]
    ),
    drop = FALSE,
    guide = "none"
  ) +
  ggplot2::scale_color_manual(
    values = layer_colors,
    breaks = unname(
      layer_labels[
        modality_order
      ]
    ),
    drop = FALSE,
    guide = "none"
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(
        0.03,
        0.08
      )
    )
  ) +
  ggplot2::labs(
    title = "Sample-level disagreement from the four-omics consensus",
    subtitle = paste0(
      "Euclidean distance across all ",
      gpa_axes_used,
      " retained GPA axes; lower values indicate closer agreement."
    ),
    x = NULL,
    y = "Distance to GPA consensus"
  ) +
  ggplot2::theme_classic(
    base_size = 10.5
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 11.3
    ),
    plot.subtitle = ggplot2::element_text(
      size = 8.9,
      color = "grey35"
    ),
    axis.text.x = ggplot2::element_text(
      face = "bold",
      size = 9.0
    ),
    axis.ticks.x = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(
      5,
      6,
      4,
      5
    )
  )

finite_agreement_ci <- c(
  gpa_agreement_bootstrap_summary$ci_low,
  gpa_agreement_bootstrap_summary$ci_high,
  gpa_agreement_bootstrap_summary$agreement
)

finite_agreement_ci <- finite_agreement_ci[
  is.finite(
    finite_agreement_ci
  )
]

gpa_agreement_x_limits <- if (
  length(finite_agreement_ci) > 0
) {
  c(
    max(
      0,
      min(finite_agreement_ci) - 0.06
    ),
    min(
      1.05,
      max(finite_agreement_ci) + 0.07
    )
  )
} else {
  c(0, 1.05)
}

gpa_agreement_bootstrap_summary$label_x <- pmin(
  gpa_agreement_x_limits[2] - 0.01,
  ifelse(
    is.finite(
      gpa_agreement_bootstrap_summary$ci_high
    ),
    gpa_agreement_bootstrap_summary$ci_high +
      0.014,
    gpa_agreement_bootstrap_summary$agreement +
      0.014
  )
)

p_gpa_agreement_bootstrap_ci <- ggplot2::ggplot(
  gpa_agreement_bootstrap_summary,
  ggplot2::aes(
    y = data_type
  )
) +
  ggplot2::geom_vline(
    xintercept = gpa_4omics_result$overall_agreement,
    linetype = "dashed",
    linewidth = 0.55,
    color = "grey55"
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = ci_low,
      xend = ci_high,
      yend = data_type,
      color = data_type
    ),
    linewidth = 1.00,
    lineend = "round",
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = agreement,
      fill = data_type
    ),
    shape = 21,
    size = 3.4,
    color = "grey20",
    stroke = 0.60
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      x = label_x,
      label = sprintf(
        "%.2f",
        agreement
      ),
      color = data_type
    ),
    hjust = 0,
    size = 3.0,
    show.legend = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = layer_colors,
    guide = "none"
  ) +
  ggplot2::scale_color_manual(
    values = layer_colors,
    guide = "none"
  ) +
  ggplot2::scale_x_continuous(
    limits = gpa_agreement_x_limits,
    breaks = scales::pretty_breaks(
      n = 4
    ),
    expand = ggplot2::expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  ggplot2::scale_y_discrete(
    expand = ggplot2::expansion(
      add = c(
        0.12,
        0.12
      )
    )
  ) +
  ggplot2::labs(
    title = "Agreement with the four-omics consensus",
    subtitle = paste0(
      "Points are observed values; horizontal intervals are percentile 95% CIs from ",
      scales::comma(
        gpa_agreement_bootstraps
      ),
      " ",
      gpa_bootstrap_unit,
      " bootstrap resamples."
    ),
    x = "Procrustes cosine similarity",
    y = NULL
  ) +
  ggplot2::coord_cartesian(
    clip = "off"
  ) +
  ggplot2::theme_classic(
    base_size = 10.2
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 11.1
    ),
    plot.subtitle = ggplot2::element_text(
      size = 8.6,
      color = "grey35"
    ),
    axis.text.y = ggplot2::element_text(
      face = "bold",
      size = 9.0
    ),
    axis.ticks.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(
      5,
      12,
      4,
      5
    )
  )

# Standalone export dimensions are intentionally specified directly below.
# To change only panel 2 or panel 3 width later, edit the corresponding
# ggsave(width = ..., height = ...) values; the ggplot objects have no fixed
# aspect ratio and will reflow without code changes.
ggplot2::ggsave(
  filename = paste0(
    "figures/coherence/",
    "Generalized_Procrustes_4omics_",
    "consensus_v20.svg"
  ),
  plot = p_gpa_4omics_consensus,
  width = 7.2,
  height = 6.25,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

ggplot2::ggsave(
  filename = paste0(
    "figures/coherence/",
    "Generalized_Procrustes_4omics_",
    "consensus_pairwise_chain_v20.svg"
  ),
  plot = p_gpa_4omics_chain_overview,
  width = 11.3,
  height = 7.45,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

ggplot2::ggsave(
  filename = paste0(
    "figures/coherence/",
    "Generalized_Procrustes_4omics_",
    "sample_residual_distribution_v20.svg"
  ),
  plot = p_gpa_residual_distribution,
  width = 4.1,
  height = 4.25,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

ggplot2::ggsave(
  filename = paste0(
    "figures/coherence/",
    "Generalized_Procrustes_4omics_",
    "agreement_bootstrap_CI_v20.svg"
  ),
  plot = p_gpa_agreement_bootstrap_ci,
  width = 4.8,
  height = 3.05,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)


#=================================================================#
# 5. Host gene-count sensitivity for Metabolite vs Host RNA-seq
#=================================================================#
#
# The primary analysis uses the host_top_n selected in input processing.
# This sensitivity analysis repeats only the host ordination and the
# Metabolite–Host Procrustes fit across multiple unsupervised VST-variance
# thresholds. Metabolite preprocessing, pairwise samples, and the first five
# positive PCoA axes remain fixed.
#
#-----------------------------------------------------------------#

make_positive_pcoa_scores <- function(
    distance_object,
    max_axes
) {
  pcoa_object <- ape::pcoa(
    distance_object
  )
  
  positive_axes <- which(
    pcoa_object$values$Eigenvalues > 1e-8
  )
  
  if (length(positive_axes) == 0) {
    return(
      NULL
    )
  }
  
  positive_axes <- positive_axes[
    seq_len(
      min(
        max_axes,
        length(
          positive_axes
        )
      )
    )
  ]
  
  score_matrix <- as.matrix(
    pcoa_object$vectors[
      ,
      positive_axes,
      drop = FALSE
    ]
  )
  
  colnames(score_matrix) <- paste0(
    "Axis",
    seq_len(
      ncol(
        score_matrix
      )
    )
  )
  
  score_matrix
}

host_gene_sensitivity_summary <- data.frame()
host_gene_sensitivity_plot_data <- list()
host_gene_sensitivity_plots <- list()
p_host_gene_sensitivity_grid <- NULL
p_host_gene_sensitivity_r <- NULL

if (
  !is.null(
    coherence_data$host_sensitivity
  ) &&
  !is.null(
    coherence_data$host_sensitivity$vst_qc_all
  )
) {
  host_vst_qc_all <- as.matrix(
    coherence_data$host_sensitivity$vst_qc_all
  )
  
  host_gene_variance_all <-
    coherence_data$host_sensitivity$gene_variance
  
  host_gene_counts <- unique(
    as.integer(
      coherence_data$host_sensitivity$candidate_top_n
    )
  )
  
  host_gene_counts <- sort(
    host_gene_counts[
      is.finite(
        host_gene_counts
      ) &
        host_gene_counts >= 2 &
        host_gene_counts <= ncol(
          host_vst_qc_all
        )
    ]
  )
  
  metabolite_host_pair <- coherence_get_pair(
    coherence_data,
    "metabolite_host"
  )
  
  metabolite_scores_all <- original_score_matrix(
    metabolite_host_pair,
    "metabolite"
  )
  
  sensitivity_common_samples <-
    metabolite_host_pair$sample_meta$SampleID[
      metabolite_host_pair$sample_meta$SampleID %in%
        intersect(
          rownames(
            metabolite_scores_all
          ),
          rownames(
            host_vst_qc_all
          )
        )
    ]
  
  sensitivity_meta <- metabolite_host_pair$sample_meta[
    match(
      sensitivity_common_samples,
      metabolite_host_pair$sample_meta$SampleID
    ),
    ,
    drop = FALSE
  ]
  
  sensitivity_meta$TRG_plot <- factor(
    as.character(
      sensitivity_meta$TRG_plot
    ),
    levels = c(
      "pCR",
      "non_pCR"
    )
  )
  
  for (
    sensitivity_index in
    seq_along(
      host_gene_counts
    )
  ) {
    gene_count <- host_gene_counts[
      sensitivity_index
    ]
    
    selected_genes <- utils::head(
      as.character(
        host_gene_variance_all$Gene
      ),
      gene_count
    )
    
    host_analysis_all <- as.matrix(
      host_vst_qc_all[
        ,
        selected_genes,
        drop = FALSE
      ]
    )
    
    complete_features <- colSums(
      is.finite(
        host_analysis_all
      )
    ) == nrow(
      host_analysis_all
    )
    
    host_analysis_all <- host_analysis_all[
      ,
      complete_features,
      drop = FALSE
    ]
    
    if (ncol(host_analysis_all) < 2) {
      next
    }
    
    # DESeq2 VST values are used without gene-wise z-scaling. Euclidean
    # PCoA on these values is equivalent to PCA of the centered sample-by-
    # gene matrix, up to axis sign/rotation.
    host_distance <- stats::dist(
      host_analysis_all[
        sensitivity_common_samples,
        ,
        drop = FALSE
      ],
      method = "euclidean"
    )
    
    host_scores <- make_positive_pcoa_scores(
      host_distance,
      protest_max_axes
    )
    
    if (is.null(host_scores)) {
      next
    }
    
    axes_used <- min(
      protest_max_axes,
      ncol(
        metabolite_scores_all
      ),
      ncol(
        host_scores
      ),
      length(
        sensitivity_common_samples
      ) - 1
    )
    
    if (axes_used < 2) {
      next
    }
    
    metabolite_scores <- metabolite_scores_all[
      sensitivity_common_samples,
      seq_len(
        axes_used
      ),
      drop = FALSE
    ]
    
    host_scores <- host_scores[
      sensitivity_common_samples,
      seq_len(
        axes_used
      ),
      drop = FALSE
    ]
    
    sensitivity_proc <- vegan::procrustes(
      metabolite_scores,
      host_scores,
      symmetric = TRUE
    )
    
    set.seed(
      20268718 +
        sensitivity_index
    )
    
    sensitivity_prot <- vegan::protest(
      metabolite_scores,
      host_scores,
      permutations =
        host_gene_sensitivity_permutations,
      symmetric = TRUE
    )
    
    gene_count_label <- if (
      gene_count == ncol(
        host_vst_qc_all
      )
    ) {
      paste0(
        "All QC genes\n(n = ",
        scales::comma(
          gene_count
        ),
        ")"
      )
    } else {
      paste0(
        "Top ",
        scales::comma(
          gene_count
        ),
        " genes\n(",
        scales::percent(
          host_gene_variance_all$cumulative_variance_fraction[
            gene_count
          ],
          accuracy = 1
        ),
        " cumulative variance)"
      )
    }
    
    host_gene_sensitivity_summary <- rbind(
      host_gene_sensitivity_summary,
      data.frame(
        gene_count = gene_count,
        gene_count_label =
          gene_count_label,
        n_samples = length(
          sensitivity_common_samples
        ),
        axes_used = axes_used,
        protest_r = unname(
          sensitivity_prot$t0
        ),
        protest_p = unname(
          sensitivity_prot$signif
        ),
        permutations =
          host_gene_sensitivity_permutations,
        primary_setting =
          gene_count ==
          coherence_data$host_sensitivity$primary_top_n,
        variance_elbow_setting =
          gene_count ==
          coherence_data$host_sensitivity$variance_elbow_rank,
        cumulative_variance_fraction =
          host_gene_variance_all$cumulative_variance_fraction[
            gene_count
          ],
        stringsAsFactors = FALSE
      )
    )
    
    sensitivity_segment <- data.frame(
      x = sensitivity_proc$X[
        ,
        1
      ],
      y = sensitivity_proc$X[
        ,
        2
      ],
      xend = sensitivity_proc$Yrot[
        ,
        1
      ],
      yend = sensitivity_proc$Yrot[
        ,
        2
      ],
      stringsAsFactors = FALSE
    )
    
    sensitivity_points <- rbind(
      data.frame(
        Axis1 = sensitivity_proc$X[
          ,
          1
        ],
        Axis2 = sensitivity_proc$X[
          ,
          2
        ],
        Layer = "Metabolite",
        TRG_plot = sensitivity_meta$TRG_plot,
        stringsAsFactors = FALSE
      ),
      data.frame(
        Axis1 = sensitivity_proc$Yrot[
          ,
          1
        ],
        Axis2 = sensitivity_proc$Yrot[
          ,
          2
        ],
        Layer = "Host RNA-seq",
        TRG_plot = sensitivity_meta$TRG_plot,
        stringsAsFactors = FALSE
      )
    )
    
    sensitivity_points$Layer <- factor(
      sensitivity_points$Layer,
      levels = unname(
        layer_labels
      )
    )
    
    host_gene_sensitivity_plot_data[[as.character(
      gene_count
    )]
    ] <- list(
      segment = sensitivity_segment,
      points = sensitivity_points,
      label = gene_count_label,
      r = unname(
        sensitivity_prot$t0
      ),
      p = unname(
        sensitivity_prot$signif
      )
    )
  }
  
  if (nrow(host_gene_sensitivity_summary) > 0) {
    sensitivity_limit <- max(
      unlist(
        lapply(
          host_gene_sensitivity_plot_data,
          function(x) {
            abs(
              c(
                x$segment$x,
                x$segment$y,
                x$segment$xend,
                x$segment$yend
              )
            )
          }
        )
      ),
      na.rm = TRUE
    )
    
    if (
      !is.finite(
        sensitivity_limit
      ) ||
      sensitivity_limit <= 0
    ) {
      sensitivity_limit <- 1
    }
    
    sensitivity_limit <-
      sensitivity_limit * 1.08
    
    for (
      sensitivity_index in
      seq_len(
        nrow(
          host_gene_sensitivity_summary
        )
      )
    ) {
      gene_count <-
        host_gene_sensitivity_summary$gene_count[
          sensitivity_index
        ]
      
      plot_data <- host_gene_sensitivity_plot_data[[as.character(
        gene_count
      )]
      ]
      
      panel_annotation <- paste0(
        "r = ",
        sprintf(
          "%.2f",
          plot_data$r
        ),
        "\nP = ",
        format.pval(
          plot_data$p,
          digits = 2,
          eps = 0.001
        )
      )
      
      host_gene_sensitivity_plots[[as.character(
        gene_count
      )]
      ] <- ggplot2::ggplot() +
        ggplot2::geom_segment(
          data = plot_data$segment,
          ggplot2::aes(
            x = x,
            y = y,
            xend = xend,
            yend = yend
          ),
          color = "grey80",
          linewidth = 0.32
        ) +
        ggplot2::geom_point(
          data = plot_data$points,
          ggplot2::aes(
            x = Axis1,
            y = Axis2,
            shape = Layer,
            fill = Layer
          ),
          color = "grey25",
          size = 2.65,
          stroke = 0.60
        ) +
        ggplot2::annotate(
          "text",
          x = -sensitivity_limit * 0.96,
          y = sensitivity_limit * 0.96,
          label = panel_annotation,
          hjust = 0,
          vjust = 1,
          size = 3.0,
          fontface = "bold"
        ) +
        ggplot2::scale_shape_manual(
          values = layer_shapes,
          breaks = c(
            "Metabolite",
            "Host RNA-seq"
          ),
          drop = TRUE,
          name = "Data type"
        ) +
        ggplot2::scale_fill_manual(
          values = layer_colors,
          breaks = c(
            "Metabolite",
            "Host RNA-seq"
          ),
          drop = TRUE,
          name = "Data type"
        ) +
        ggplot2::labs(
          title = plot_data$label,
          x = NULL,
          y = NULL
        ) +
        ggplot2::coord_fixed(
          ratio = 1,
          xlim = c(
            -sensitivity_limit,
            sensitivity_limit
          ),
          ylim = c(
            -sensitivity_limit,
            sensitivity_limit
          ),
          expand = FALSE
        ) +
        ggplot2::theme_classic(
          base_size = 9.2
        ) +
        ggplot2::theme(
          aspect.ratio = 1,
          plot.title = ggplot2::element_text(
            face = "bold",
            hjust = 0.5,
            size = 9.4
          ),
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          legend.position = "right"
        )
    }
    
    p_host_gene_sensitivity_grid <-
      patchwork::wrap_plots(
        host_gene_sensitivity_plots[
          as.character(
            host_gene_sensitivity_summary$gene_count
          )
        ],
        ncol = 3,
        guides = "collect"
      ) &
      ggplot2::theme(
        legend.position = "right"
      )
    
    ggplot2::ggsave(
      paste0(
        "figures/coherence/",
        "Procrustes_metabolite_host_",
        "host_gene_count_scatter_grid.svg"
      ),
      p_host_gene_sensitivity_grid,
      width = 11.2,
      height = 7.0,
      units = "in",
      device = svglite::svglite,
      bg = "white"
    )
    
    host_gene_sensitivity_summary$setting <- ifelse(
      host_gene_sensitivity_summary$primary_setting,
      "Primary",
      "Sensitivity"
    )
    
    p_host_gene_sensitivity_r <- ggplot2::ggplot(
      host_gene_sensitivity_summary,
      ggplot2::aes(
        x = gene_count,
        y = protest_r,
        group = 1
      )
    ) +
      ggplot2::geom_line(
        linewidth = 0.75,
        color = "grey45"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(
          fill = setting
        ),
        shape = 21,
        size = 3.5,
        color = "grey20",
        stroke = 0.70
      ) +
      ggplot2::geom_text(
        ggplot2::aes(
          label = paste0(
            "P=",
            format.pval(
              protest_p,
              digits = 2,
              eps = 0.001
            )
          )
        ),
        nudge_y = 0.025,
        size = 2.9,
        check_overlap = TRUE
      ) +
      ggplot2::geom_vline(
        xintercept =
          coherence_data$host_sensitivity$primary_top_n,
        linewidth = 0.65,
        color = "grey20"
      ) +
      ggplot2::geom_vline(
        xintercept =
          coherence_data$host_sensitivity$variance_elbow_rank,
        linewidth = 0.55,
        linetype = "longdash",
        color = "grey55",
        na.rm = TRUE
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          Primary = "#365F7D",
          Sensitivity = "white"
        ),
        name = NULL
      ) +
      ggplot2::scale_x_continuous(
        labels = scales::label_comma(),
        breaks = host_gene_sensitivity_summary$gene_count
      ) +
      ggplot2::scale_y_continuous(
        limits = c(
          0,
          1
        ),
        breaks = seq(
          0,
          1,
          by = 0.2
        )
      ) +
      ggplot2::labs(
        x = "Number of highest-VST-variance host genes retained",
        y = "Metabolite–Host PROTEST r",
        caption = paste0(
          "Solid vertical line: primary ",
          scales::percent(
            coherence_data$host_sensitivity$primary_cumulative_fraction,
            accuracy = 1
          ),
          " cumulative marginal VST-variance setting (n = ",
          scales::comma(
            coherence_data$host_sensitivity$primary_top_n
          ),
          "); long-dashed line: variance-rank elbow. ",
          "Host VST genes are not gene-wise z-scaled. ",
          "P values use 999 permutations for sensitivity analysis."
        )
      ) +
      ggplot2::theme_classic(
        base_size = 10.5
      ) +
      ggplot2::theme(
        legend.position = "top",
        axis.text.x = ggplot2::element_text(
          angle = 35,
          hjust = 1
        ),
        plot.caption = ggplot2::element_text(
          hjust = 0,
          size = 8.2
        )
      )
    
    ggplot2::ggsave(
      paste0(
        "figures/coherence/",
        "Procrustes_metabolite_host_",
        "host_gene_count_r_sensitivity.svg"
      ),
      p_host_gene_sensitivity_r,
      width = 7.6,
      height = 4.7,
      units = "in",
      device = svglite::svglite,
      bg = "white"
    )
  }
} else {
  warning(
    "host_sensitivity data were not found in coherence_data; ",
    "the host gene-count sensitivity analysis was skipped."
  )
}

#=================================================================#
# 6. Save results
#=================================================================#

write.csv(
  original_ordination_preprocessing,
  "results/coherence/PROTEST_original_ordination_preprocessing.csv",
  row.names = FALSE
)

write.csv(
  all_pair_summary,
  "results/coherence/PROTEST_original_ordination_all_six_pairs_summary.csv",
  row.names = FALSE
)

write.csv(
  pair_display_summary,
  "results/coherence/PROTEST_original_ordination_pair_display_summary.csv",
  row.names = FALSE
)

write.csv(
  pair_group_stats,
  "results/coherence/PROTEST_original_ordination_group_separation_summary.csv",
  row.names = FALSE
)

write.csv(
  forest_summary,
  "results/coherence/PROTEST_original_ordination_statistics_forest_summary.csv",
  row.names = FALSE
)

write.csv(
  pair_r_bootstrap_summary,
  "results/coherence/PROTEST_pair_r_subject_cluster_bootstrap_CI.csv",
  row.names = FALSE
)


write.csv(
  gpa_configuration_sources,
  "results/coherence/GPA_4omics_configuration_sources.csv",
  row.names = FALSE
)

write.csv(
  gpa_global_summary,
  "results/coherence/GPA_4omics_global_summary.csv",
  row.names = FALSE
)

write.csv(
  gpa_block_summary,
  "results/coherence/GPA_4omics_block_summary.csv",
  row.names = FALSE
)

write.csv(
  gpa_pairwise_summary,
  "results/coherence/GPA_4omics_pairwise_r_summary.csv",
  row.names = FALSE
)

write.csv(
  gpa_residual_distance_data,
  "results/coherence/GPA_4omics_sample_residual_distances.csv",
  row.names = FALSE
)

write.csv(
  gpa_agreement_bootstrap_summary,
  "results/coherence/GPA_4omics_agreement_bootstrap_CI.csv",
  row.names = FALSE
)

if (nrow(host_gene_sensitivity_summary) > 0) {
  write.csv(
    host_gene_sensitivity_summary,
    paste0(
      "results/coherence/",
      "PROTEST_metabolite_host_gene_count_sensitivity.csv"
    ),
    row.names = FALSE
  )
}

save(
  original_ordination_preprocessing,
  all_pair_summary,
  all_pair_results,
  pair_display_summary,
  pair_r_bootstrap_summary,
  pair_group_stats,
  forest_summary,
  host_gene_sensitivity_summary,
  host_gene_sensitivity_plot_data,
  host_gene_sensitivity_plots,
  p_host_gene_sensitivity_grid,
  p_host_gene_sensitivity_r,
  p_procrustes_matrix_full_core,
  p_procrustes_matrix_upper_core,
  p_procrustes_legend,
  p_procrustes_matrix_full,
  p_procrustes_matrix_upper,
  p_procrustes_forest,
  gpa_configuration_sources,
  gpa_configurations,
  gpa_4omics_result,
  gpa_permuted_residual_ss,
  gpa_global_summary,
  gpa_block_summary,
  gpa_pairwise_summary,
  gpa_consensus_plot_data,
  gpa_point_plot_data,
  gpa_segment_plot_data,
  gpa_residual_distance_data,
  gpa_residual_box_summary,
  gpa_agreement_bootstrap_values,
  gpa_agreement_bootstrap_summary,
  p_gpa_4omics_consensus,
  p_gpa_pair_species_ko,
  p_gpa_pair_ko_metabolite,
  p_gpa_pair_metabolite_host,
  p_gpa_global_summary_box,
  p_gpa_pairwise_chain_column,
  p_gpa_4omics_chain_overview,
  p_gpa_residual_distribution,
  p_gpa_agreement_bootstrap_ci,
  file =
    "results/coherence/Procrustes_PROTEST_4omics_v20_GPA_panels.RData"
)

rm(
  list = setdiff(
    ls(),
    c(
      "coherence_data",
      "original_ordination_preprocessing",
      "all_pair_summary",
      "all_pair_results",
      "pair_display_summary",
      "pair_r_bootstrap_summary",
      "pair_group_stats",
      "forest_summary",
      "host_gene_sensitivity_summary",
      "host_gene_sensitivity_plot_data",
      "host_gene_sensitivity_plots",
      "p_host_gene_sensitivity_grid",
      "p_host_gene_sensitivity_r",
      "p_procrustes_matrix_full_core",
      "p_procrustes_matrix_upper_core",
      "p_procrustes_legend",
      "p_procrustes_matrix_full",
      "p_procrustes_matrix_upper",
      "p_procrustes_forest",
      "gpa_configuration_sources",
      "gpa_configurations",
      "gpa_4omics_result",
      "gpa_permuted_residual_ss",
      "gpa_global_summary",
      "gpa_block_summary",
      "gpa_pairwise_summary",
      "gpa_consensus_plot_data",
      "gpa_point_plot_data",
      "gpa_segment_plot_data",
      "gpa_residual_distance_data",
      "gpa_residual_box_summary",
      "gpa_agreement_bootstrap_values",
      "gpa_agreement_bootstrap_summary",
      "p_gpa_4omics_consensus",
      "p_gpa_pair_species_ko",
      "p_gpa_pair_ko_metabolite",
      "p_gpa_pair_metabolite_host",
      "p_gpa_global_summary_box",
      "p_gpa_pairwise_chain_column",
      "p_gpa_4omics_chain_overview",
      "p_gpa_residual_distribution",
      "p_gpa_agreement_bootstrap_ci"
    )
  )
)

gc()
