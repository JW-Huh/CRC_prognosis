#-----------------------------------------------------------------#
#
# Four-layer coherence input processing
#
# Biological sequence:
#   Species -> Microbial function (KO) -> Metabolite -> Host RNA-seq
#
# Primary biological-chain pairwise objects:
#   species_ko
#   ko_metabolite
#   metabolite_host
#
# Additional sensitivity pairwise objects:
#   species_metabolite
#   species_host
#   ko_host
#
# Matching rules:
#   1. Pair sample-level:
#      all SubjectID + Timepoint overlaps for the relevant pair.
#      No Before-After completeness requirement.
#   2. Pair longitudinal:
#      pair-specific subjects with one Before and one Ongoing overlap.
#   3. Four-block sample-level:
#      all SubjectID + Timepoint overlaps across all four layers.
#   4. Four-block longitudinal:
#      subjects with one Before and one Ongoing four-block overlap.
#
# KO and species originate from the same metagenomic library. Therefore,
# KO does not normally reduce the sample count relative to species, but it is
# retained as a distinct biological view of microbial functional potential.
#
# Main output:
#   input/coherence_species_ko_metabolite_host_data_metabolite50_host50pct_no_scaling.RData
#
# The finalized UpSet plot is not regenerated.
#
#-----------------------------------------------------------------#

rm(list = ls())
options(stringsAsFactors = FALSE)
setwd("D:/2-연구/2-CRC metagenomics/")

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ape)
  library(DESeq2)
  library(patchwork)
  library(svglite)
})

dir.create(
  "results/coherence_input",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "input",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures/coherence",
  recursive = TRUE,
  showWarnings = FALSE
)


#=================================================================#
# 0. Analysis settings
#=================================================================#

time_levels <- c(
  "Before",
  "Ongoing"
)

species_prev <- 0.20
species_top_n <- NULL

ko_prev <- 0.20
ko_top_n <- NULL

metabolite_min_finite_fraction <- 0.50
metabolite_min_detected_fraction <- 0.50
metabolite_floor_min_repeats <- 2
metabolite_floor_relative_tolerance <- 1e-8
metabolite_floor_absolute_tolerance <- 1e-12
metabolite_transform <- "abs_z"
# metabolite_transform <- "log10_z"

host_min_count <- 10
host_min_sample_fraction <- 0.20

# Host RNA-seq feature selection is unsupervised and is based only on
# gene-wise variance after DESeq2 VST. Differential-expression status and
# treatment-response labels are not used.
#
# The primary set retains the smallest number of highest-variance genes whose
# marginal gene-wise VST variances account for 50% of the total across-gene
# variance sum. This is an explicit data-adaptive filtering rule, not PCA
# explained variance. Sensitivity analyses include fixed and cumulative-
# variance thresholds.
host_selection_mode <- "cumulative_variance"
host_primary_cumulative_fraction <- 0.50
host_sensitivity_cumulative_fractions <- c(
  0.50,
  0.75,
  0.80,
  0.90,
  0.95
)
host_sensitivity_fixed_top_n <- c(
  5000,
  12000
)
host_include_all_qc_genes_in_sensitivity <- TRUE
host_top_n <- NA_integer_

pcoa_max_axes <- 5
min_pairwise_samples <- 4
min_paired_subjects <- 3
min_four_block_samples <- 4

metabolite_start_name <- "Acetate"
metabolite_end_name <-
  "Tauroursodeoxycholic_acid_Taurohyodeoxycholic_acid"

ko_file <-
  "260224 final script/Input/merged_genefamilies_KO_named.tsv"


#=================================================================#
# 1. Reused input-processing functions
#=================================================================#

standardize_timepoint <- function(x) {
  x <- stringr::str_to_lower(
    stringr::str_trim(as.character(x))
  )
  
  dplyr::case_when(
    x %in% c(
      "before",
      "pre",
      "baseline",
      "pre-rt",
      "pre_rt"
    ) ~ "Before",
    x %in% c(
      "ongoing",
      "post",
      "after",
      "post-rt",
      "post_rt",
      "post rt",
      "during rt"
    ) ~ "Ongoing",
    TRUE ~ NA_character_
  )
}


make_trg_plot <- function(x) {
  dplyr::case_when(
    as.character(x) == "CR" ~ "pCR",
    as.character(x) == "nonCR" ~ "non_pCR",
    TRUE ~ NA_character_
  )
}


check_subject_timepoint_duplicates <- function(
    meta,
    omic_name
) {
  meta %>%
    dplyr::count(
      SubjectID,
      Timepoint,
      name = "n_samples"
    ) %>%
    dplyr::filter(n_samples > 1) %>%
    dplyr::mutate(
      Omic = omic_name,
      .before = 1
    )
}


make_pcoa_result <- function(
    distance_object,
    max_axes = 5
) {
  pcoa_object <- ape::pcoa(distance_object)
  
  positive_axes <- which(
    pcoa_object$values$Eigenvalues > 1e-8
  )
  
  if (length(positive_axes) == 0) {
    return(
      list(
        object = pcoa_object,
        scores = NULL
      )
    )
  }
  
  positive_axes <- positive_axes[
    seq_len(
      min(
        max_axes,
        length(positive_axes)
      )
    )
  ]
  
  pcoa_scores <- as.data.frame(
    pcoa_object$vectors[
      ,
      positive_axes,
      drop = FALSE
    ],
    check.names = FALSE
  )
  
  colnames(pcoa_scores) <- paste0(
    "Axis",
    seq_len(ncol(pcoa_scores))
  )
  
  pcoa_scores$SampleID <- rownames(
    pcoa_scores
  )
  
  list(
    object = pcoa_object,
    scores = pcoa_scores
  )
}


calculate_rank_curve_elbow <- function(sorted_variance) {
  sorted_variance <- as.numeric(
    sorted_variance
  )
  
  keep_index <- is.finite(
    sorted_variance
  ) &
    sorted_variance > 0
  
  sorted_variance <- sorted_variance[
    keep_index
  ]
  
  if (length(sorted_variance) < 20) {
    return(
      NA_integer_
    )
  }
  
  x <- seq_along(
    sorted_variance
  )
  y <- log10(
    sorted_variance
  )
  
  x_scaled <- (
    x - min(x)
  ) / (
    max(x) - min(x)
  )
  
  y_scaled <- (
    y - min(y)
  ) / (
    max(y) - min(y)
  )
  
  # Perpendicular distance from the line connecting the first and last
  # points of the ranked log-variance curve. Extreme 2.5% ranks are excluded
  # so the selected elbow is not forced to an endpoint.
  line_distance <- abs(
    x_scaled +
      y_scaled -
      1
  ) / sqrt(2)
  
  eligible_index <- x >= ceiling(
    0.025 * length(x)
  ) &
    x <= floor(
      0.975 * length(x)
    )
  
  if (!any(eligible_index)) {
    return(
      NA_integer_
    )
  }
  
  eligible_rank <- x[
    eligible_index
  ]
  eligible_rank[
    which.max(
      line_distance[
        eligible_index
      ]
    )
  ]
}


prepare_compositional_block <- function(
    data_matrix,
    block_name,
    prevalence_cutoff
) {
  data_matrix <- as.matrix(data_matrix)
  storage.mode(data_matrix) <- "numeric"
  
  filtered_matrix <- data_matrix[
    ,
    colMeans(
      data_matrix > 0,
      na.rm = TRUE
    ) >= prevalence_cutoff &
      colSums(
        data_matrix,
        na.rm = TRUE
      ) > 0,
    drop = FALSE
  ]
  
  if (ncol(filtered_matrix) == 0) {
    stop(
      "No features remained in the ",
      block_name,
      " block after prevalence filtering."
    )
  }
  
  sample_total <- rowSums(
    filtered_matrix,
    na.rm = TRUE
  )
  
  if (
    any(
      !is.finite(sample_total) |
      sample_total <= 0
    )
  ) {
    stop(
      "At least one ",
      block_name,
      " sample has zero total abundance."
    )
  }
  
  relative_matrix <- sweep(
    filtered_matrix,
    1,
    sample_total,
    "/"
  )
  
  transformed_matrix <- sqrt(
    relative_matrix
  )
  
  feature_variance <- apply(
    transformed_matrix,
    2,
    stats::var,
    na.rm = TRUE
  )
  
  keep_features <- names(
    feature_variance
  )[
    is.finite(feature_variance) &
      feature_variance > 0
  ]
  
  filtered_matrix <- filtered_matrix[
    ,
    keep_features,
    drop = FALSE
  ]
  
  relative_matrix <- relative_matrix[
    ,
    keep_features,
    drop = FALSE
  ]
  
  transformed_matrix <- transformed_matrix[
    ,
    keep_features,
    drop = FALSE
  ]
  
  scaled_matrix <- as.matrix(
    base::scale(transformed_matrix)
  )
  
  storage.mode(scaled_matrix) <- "numeric"
  
  block_dist <- vegan::vegdist(
    relative_matrix,
    method = "bray"
  )
  
  block_pcoa <- make_pcoa_result(
    block_dist,
    pcoa_max_axes
  )
  
  list(
    raw = data_matrix,
    filtered = filtered_matrix,
    relative = relative_matrix,
    transformed = transformed_matrix,
    table = scaled_matrix,
    dist = block_dist,
    pcoa = block_pcoa$object,
    pcoa_scores = block_pcoa$scores,
    n_features = ncol(scaled_matrix)
  )
}



qc_metabolite_detection_floor <- function(
    data_matrix,
    min_finite_fraction,
    min_detected_fraction,
    floor_min_repeats,
    floor_relative_tolerance,
    floor_absolute_tolerance
) {
  data_matrix <- as.matrix(data_matrix)
  storage.mode(data_matrix) <- "numeric"
  data_matrix[!is.finite(data_matrix)] <- NA_real_
  
  if (
    any(
      data_matrix < 0,
      na.rm = TRUE
    )
  ) {
    stop(
      "Negative metabolite concentrations were found. ",
      "The detection-floor QC assumes non-negative concentration data.",
      call. = FALSE
    )
  }
  
  qc_summary <- data.frame()
  qc_matrix <- data_matrix
  
  for (feature_name in colnames(data_matrix)) {
    x <- data_matrix[, feature_name]
    finite_index <- is.finite(x)
    finite_x <- x[finite_index]
    
    finite_fraction <- mean(finite_index)
    
    if (length(finite_x) == 0) {
      qc_summary <- rbind(
        qc_summary,
        data.frame(
          metabolite = feature_name,
          n_samples = length(x),
          n_finite = 0,
          finite_fraction = 0,
          floor_value = NA_real_,
          zero_floor = NA,
          floor_repeated = FALSE,
          floor_count = 0,
          floor_fraction = NA_real_,
          detected_count = 0,
          detected_fraction = 0,
          n_unique_finite = 0,
          n_unique_detected = 0,
          variance_after_floor_replacement = NA_real_,
          floor_replacement_value = NA_real_,
          keep = FALSE,
          exclusion_reason = "no finite values",
          stringsAsFactors = FALSE
        )
      )
      
      next
    }
    
    floor_value <- min(
      finite_x,
      na.rm = TRUE
    )
    
    floor_tolerance <- max(
      abs(floor_value) *
        floor_relative_tolerance,
      floor_absolute_tolerance
    )
    
    floor_index <-
      finite_index &
      abs(
        x -
          floor_value
      ) <=
      floor_tolerance
    
    floor_count <- sum(
      floor_index
    )
    
    zero_floor <-
      abs(
        floor_value
      ) <=
      floor_absolute_tolerance
    
    floor_repeated <-
      zero_floor ||
      floor_count >=
      floor_min_repeats
    
    if (floor_repeated) {
      detected_index <-
        finite_index &
        !floor_index
    } else {
      detected_index <- finite_index
    }
    
    detected_x <- x[
      detected_index
    ]
    
    detected_fraction <- mean(
      detected_index
    )
    
    floor_fraction <- if (floor_repeated) {
      mean(
        floor_index
      )
    } else {
      0
    }
    
    replacement_value <- NA_real_
    
    if (floor_repeated) {
      if (
        is.finite(
          floor_value
        ) &&
        floor_value > 0
      ) {
        replacement_value <-
          floor_value / 2
      } else {
        positive_detected <- detected_x[
          is.finite(
            detected_x
          ) &
            detected_x > 0
        ]
        
        replacement_value <- if (
          length(
            positive_detected
          ) > 0
        ) {
          min(
            positive_detected,
            na.rm = TRUE
          ) / 2
        } else {
          0
        }
      }
      
      qc_matrix[
        floor_index,
        feature_name
      ] <- replacement_value
    }
    
    variance_after_replacement <- stats::var(
      qc_matrix[
        ,
        feature_name
      ],
      na.rm = TRUE
    )
    
    keep_finite <-
      finite_fraction >=
      min_finite_fraction
    
    keep_detected <-
      detected_fraction >=
      min_detected_fraction
    
    keep_unique <-
      length(
        unique(
          detected_x[
            is.finite(
              detected_x
            )
          ]
        )
      ) >= 2
    
    keep_variance <-
      is.finite(
        variance_after_replacement
      ) &&
      variance_after_replacement > 0
    
    keep_feature <-
      keep_finite &&
      keep_detected &&
      keep_unique &&
      keep_variance
    
    exclusion_reason <- if (keep_feature) {
      "retained"
    } else {
      paste(
        c(
          if (!keep_finite) {
            paste0(
              "finite fraction < ",
              min_finite_fraction
            )
          },
          if (!keep_detected) {
            paste0(
              "detected fraction above repeated floor < ",
              min_detected_fraction
            )
          },
          if (!keep_unique) {
            "fewer than two unique detected values"
          },
          if (!keep_variance) {
            "zero or non-finite variance"
          }
        ),
        collapse = "; "
      )
    }
    
    qc_summary <- rbind(
      qc_summary,
      data.frame(
        metabolite = feature_name,
        n_samples = length(x),
        n_finite = sum(
          finite_index
        ),
        finite_fraction =
          finite_fraction,
        floor_value =
          floor_value,
        zero_floor =
          zero_floor,
        floor_repeated =
          floor_repeated,
        floor_count =
          floor_count,
        floor_fraction =
          floor_fraction,
        detected_count = sum(
          detected_index
        ),
        detected_fraction =
          detected_fraction,
        n_unique_finite = length(
          unique(
            finite_x
          )
        ),
        n_unique_detected = length(
          unique(
            detected_x[
              is.finite(
                detected_x
              )
            ]
          )
        ),
        variance_after_floor_replacement =
          variance_after_replacement,
        floor_replacement_value =
          replacement_value,
        keep =
          keep_feature,
        exclusion_reason =
          exclusion_reason,
        stringsAsFactors = FALSE
      )
    )
  }
  
  retained_features <- qc_summary$metabolite[
    qc_summary$keep
  ]
  
  if (length(retained_features) == 0) {
    stop(
      "No metabolite passed the detection-floor prevalence and variance QC.",
      call. = FALSE
    )
  }
  
  list(
    matrix = qc_matrix[
      ,
      retained_features,
      drop = FALSE
    ],
    summary = qc_summary,
    retained_features =
      retained_features
  )
}


prepare_continuous_block <- function(
    data_matrix,
    block_name,
    transform = c(
      "none",
      "log10"
    ),
    scale_features = TRUE
) {
  transform <- match.arg(transform)
  
  data_matrix <- as.matrix(data_matrix)
  storage.mode(data_matrix) <- "numeric"
  data_matrix[!is.finite(data_matrix)] <- NA_real_
  
  if (transform == "log10") {
    positive_values <- data_matrix[
      is.finite(data_matrix) &
        data_matrix > 0
    ]
    
    pseudocount <- ifelse(
      length(positive_values) > 0,
      min(positive_values) / 2,
      1e-8
    )
    
    transformed_matrix <- log10(
      data_matrix + pseudocount
    )
  } else {
    pseudocount <- NA_real_
    transformed_matrix <- data_matrix
  }
  
  for (j in seq_len(ncol(transformed_matrix))) {
    if (anyNA(transformed_matrix[, j])) {
      feature_median <- stats::median(
        transformed_matrix[, j],
        na.rm = TRUE
      )
      
      if (!is.finite(feature_median)) {
        feature_median <- 0
      }
      
      transformed_matrix[
        is.na(transformed_matrix[, j]),
        j
      ] <- feature_median
    }
  }
  
  feature_variance <- apply(
    transformed_matrix,
    2,
    stats::var,
    na.rm = TRUE
  )
  
  transformed_matrix <- transformed_matrix[
    ,
    is.finite(feature_variance) &
      feature_variance > 0,
    drop = FALSE
  ]
  
  if (ncol(transformed_matrix) == 0) {
    stop(
      "No non-constant features remained in the ",
      block_name,
      " block."
    )
  }
  
  analysis_matrix <- if (scale_features) {
    as.matrix(
      base::scale(
        transformed_matrix
      )
    )
  } else {
    transformed_matrix
  }
  
  storage.mode(analysis_matrix) <- "numeric"
  
  block_dist <- stats::dist(
    analysis_matrix,
    method = "euclidean"
  )
  
  block_pcoa <- make_pcoa_result(
    block_dist,
    pcoa_max_axes
  )
  
  list(
    raw = data_matrix,
    transformed = transformed_matrix,
    table = analysis_matrix,
    dist = block_dist,
    pcoa = block_pcoa$object,
    pcoa_scores = block_pcoa$scores,
    pseudocount = pseudocount,
    feature_scaling = ifelse(
      scale_features,
      "gene/feature-wise z-score",
      "none"
    ),
    n_features = ncol(analysis_matrix)
  )
}


prepare_named_block <- function(
    block_name,
    data_matrix
) {
  if (block_name == "species") {
    return(
      prepare_compositional_block(
        data_matrix,
        "species",
        species_prev
      )
    )
  }
  
  if (block_name == "ko") {
    return(
      prepare_compositional_block(
        data_matrix,
        "KO",
        ko_prev
      )
    )
  }
  
  if (block_name == "metabolite") {
    return(
      prepare_continuous_block(
        data_matrix,
        "metabolite",
        ifelse(
          metabolite_transform == "log10_z",
          "log10",
          "none"
        ),
        scale_features = TRUE
      )
    )
  }
  
  if (block_name == "host") {
    return(
      prepare_continuous_block(
        data_matrix,
        "host RNA-seq",
        "none",
        scale_features = FALSE
      )
    )
  }
  
  stop(
    "Unknown block: ",
    block_name
  )
}


find_complete_longitudinal_subjects <- function(meta) {
  meta %>%
    dplyr::count(
      SubjectID,
      Timepoint,
      name = "n"
    ) %>%
    dplyr::filter(n == 1) %>%
    dplyr::count(
      SubjectID,
      name = "n_timepoints"
    ) %>%
    dplyr::filter(
      n_timepoints ==
        length(time_levels)
    ) %>%
    dplyr::pull(SubjectID)
}


make_delta_result <- function(
    analysis_table,
    longitudinal_meta
) {
  if (nrow(longitudinal_meta) == 0) {
    return(
      list(
        table = matrix(
          numeric(0),
          nrow = 0,
          ncol = ncol(analysis_table),
          dimnames = list(
            character(0),
            colnames(analysis_table)
          )
        ),
        dist = NULL,
        pcoa = NULL,
        pcoa_scores = NULL
      )
    )
  }
  
  before_meta <- longitudinal_meta %>%
    dplyr::filter(
      Timepoint == "Before"
    ) %>%
    dplyr::arrange(SubjectID)
  
  after_meta <- longitudinal_meta %>%
    dplyr::filter(
      Timepoint == "Ongoing"
    ) %>%
    dplyr::arrange(SubjectID)
  
  if (
    !identical(
      before_meta$SubjectID,
      after_meta$SubjectID
    )
  ) {
    stop(
      "Before and Ongoing subjects are not aligned."
    )
  }
  
  delta_table <-
    analysis_table[
      after_meta$SampleID,
      ,
      drop = FALSE
    ] -
    analysis_table[
      before_meta$SampleID,
      ,
      drop = FALSE
    ]
  
  rownames(delta_table) <-
    before_meta$SubjectID
  
  if (nrow(delta_table) >= 2) {
    delta_dist <- stats::dist(
      delta_table,
      method = "euclidean"
    )
    
    delta_pcoa <- make_pcoa_result(
      delta_dist,
      pcoa_max_axes
    )
  } else {
    delta_dist <- NULL
    delta_pcoa <- list(
      object = NULL,
      scores = NULL
    )
  }
  
  list(
    table = delta_table,
    dist = delta_dist,
    pcoa = delta_pcoa$object,
    pcoa_scores = delta_pcoa$scores
  )
}


make_sample_set_index <- function(meta) {
  list(
    all_samples = as.character(
      meta$SampleID
    ),
    baseline_samples = as.character(
      meta$SampleID[
        meta$Timepoint == "Before"
      ]
    ),
    after_rt_samples = as.character(
      meta$SampleID[
        meta$Timepoint == "Ongoing"
      ]
    ),
    pCR_samples = as.character(
      meta$SampleID[
        meta$TRG_plot == "pCR"
      ]
    ),
    non_pCR_samples = as.character(
      meta$SampleID[
        meta$TRG_plot == "non_pCR"
      ]
    )
  )
}


build_modality_object <- function(
    block_name,
    meta,
    data_matrix
) {
  sample_meta <- meta %>%
    dplyr::arrange(
      SubjectID,
      match(
        Timepoint,
        time_levels
      )
    ) %>%
    dplyr::mutate(
      SampleID = paste(
        SubjectID,
        Timepoint,
        sep = "__"
      )
    )
  
  if (
    anyDuplicated(
      sample_meta$SampleID
    ) > 0
  ) {
    stop(
      block_name,
      ": duplicated standardized SampleID."
    )
  }
  
  matched_matrix <- data_matrix[
    sample_meta$OmicSampleID,
    ,
    drop = FALSE
  ]
  
  rownames(matched_matrix) <-
    sample_meta$SampleID
  
  block <- prepare_named_block(
    block_name,
    matched_matrix
  )
  
  longitudinal_subjects <-
    find_complete_longitudinal_subjects(
      sample_meta
    )
  
  longitudinal_meta <- sample_meta %>%
    dplyr::filter(
      SubjectID %in%
        longitudinal_subjects
    ) %>%
    dplyr::arrange(
      SubjectID,
      match(
        Timepoint,
        time_levels
      )
    )
  
  delta <- make_delta_result(
    block$table,
    longitudinal_meta
  )
  
  modality_object <- list(
    block_name = block_name,
    sample_meta = sample_meta,
    sample_sets = make_sample_set_index(
      sample_meta
    ),
    longitudinal = list(
      meta = longitudinal_meta,
      subjects = longitudinal_subjects,
      sample_ids = longitudinal_meta$SampleID,
      delta_table = delta$table,
      delta_pcoa_scores =
        delta$pcoa_scores,
      status = ifelse(
        length(longitudinal_subjects) >=
          min_paired_subjects,
        "ready",
        "exploratory"
      )
    ),
    raw = block$raw,
    transformed = block$transformed,
    table = block$table,
    dist = block$dist,
    pcoa = block$pcoa,
    pcoa_scores = block$pcoa_scores,
    n_features = block$n_features
  )
  
  if (
    block_name %in% c(
      "species",
      "ko"
    )
  ) {
    modality_object$filtered <-
      block$filtered
    modality_object$relative <-
      block$relative
  }
  
  modality_object
}


build_pair_object <- function(
    pair_name,
    block_x_name,
    block_y_name,
    modality_x,
    modality_y
) {
  meta_x <- modality_x$sample_meta %>%
    dplyr::transmute(
      SubjectID,
      Timepoint,
      sample_id_x = SampleID,
      trg_x = as.character(TRG_plot)
    )
  
  meta_y <- modality_y$sample_meta %>%
    dplyr::transmute(
      SubjectID,
      Timepoint,
      sample_id_y = SampleID,
      trg_y = as.character(TRG_plot)
    )
  
  sample_meta <- dplyr::inner_join(
    meta_x,
    meta_y,
    by = c(
      "SubjectID",
      "Timepoint"
    )
  ) %>%
    dplyr::mutate(
      SampleID = paste(
        SubjectID,
        Timepoint,
        sep = "__"
      ),
      TRG_plot = dplyr::coalesce(
        trg_x,
        trg_y
      ),
      trg_discordant =
        !is.na(trg_x) &
        !is.na(trg_y) &
        trg_x != trg_y
    )
  
  if (any(sample_meta$trg_discordant)) {
    stop(
      pair_name,
      ": discordant TRG labels."
    )
  }
  
  sample_meta <- sample_meta %>%
    dplyr::arrange(
      SubjectID,
      match(
        Timepoint,
        time_levels
      )
    ) %>%
    dplyr::mutate(
      TRG_plot = factor(
        TRG_plot,
        levels = c(
          "non_pCR",
          "pCR"
        )
      )
    ) %>%
    dplyr::select(
      SampleID,
      SubjectID,
      Timepoint,
      TRG_plot,
      sample_id_x,
      sample_id_y
    )
  
  longitudinal_subjects <-
    find_complete_longitudinal_subjects(
      sample_meta
    )
  
  longitudinal_meta <- sample_meta %>%
    dplyr::filter(
      SubjectID %in%
        longitudinal_subjects
    ) %>%
    dplyr::arrange(
      SubjectID,
      match(
        Timepoint,
        time_levels
      )
    )
  
  pair_summary <- tibble::tibble(
    pair = pair_name,
    block_x = block_x_name,
    block_y = block_y_name,
    n_overlap_samples =
      nrow(sample_meta),
    n_overlap_baseline_samples =
      sum(
        sample_meta$Timepoint ==
          "Before"
      ),
    n_overlap_after_rt_samples =
      sum(
        sample_meta$Timepoint ==
          "Ongoing"
      ),
    n_overlap_subjects =
      dplyr::n_distinct(
        sample_meta$SubjectID
      ),
    n_longitudinal_subjects =
      length(
        longitudinal_subjects
      ),
    n_longitudinal_samples =
      nrow(longitudinal_meta)
  )
  
  if (nrow(sample_meta) < 2) {
    return(
      list(
        pair_name = pair_name,
        block_names = c(
          block_x_name,
          block_y_name
        ),
        sample_meta = sample_meta,
        meta = sample_meta,
        longitudinal = list(
          meta = longitudinal_meta,
          subjects =
            longitudinal_subjects,
          status = "insufficient"
        ),
        summary = pair_summary,
        status =
          "insufficient sample-level overlap"
      )
    )
  }
  
  sample_ids <- sample_meta$SampleID
  
  x_table <- modality_x$table[
    sample_ids,
    ,
    drop = FALSE
  ]
  
  y_table <- modality_y$table[
    sample_ids,
    ,
    drop = FALSE
  ]
  
  if (
    block_x_name %in% c(
      "species",
      "ko"
    )
  ) {
    x_dist <- vegan::vegdist(
      modality_x$relative[
        sample_ids,
        ,
        drop = FALSE
      ],
      method = "bray"
    )
  } else {
    x_dist <- stats::dist(
      x_table,
      method = "euclidean"
    )
  }
  
  if (
    block_y_name %in% c(
      "species",
      "ko"
    )
  ) {
    y_dist <- vegan::vegdist(
      modality_y$relative[
        sample_ids,
        ,
        drop = FALSE
      ],
      method = "bray"
    )
  } else {
    y_dist <- stats::dist(
      y_table,
      method = "euclidean"
    )
  }
  
  x_pcoa <- make_pcoa_result(
    x_dist,
    pcoa_max_axes
  )
  
  y_pcoa <- make_pcoa_result(
    y_dist,
    pcoa_max_axes
  )
  
  delta_x <- make_delta_result(
    x_table,
    longitudinal_meta
  )
  
  delta_y <- make_delta_result(
    y_table,
    longitudinal_meta
  )
  
  pair_summary$n_features_x <-
    ncol(x_table)
  
  pair_summary$n_features_y <-
    ncol(y_table)
  
  pair_object <- list(
    pair_name = pair_name,
    block_names = c(
      block_x_name,
      block_y_name
    ),
    sample_meta = sample_meta,
    meta = sample_meta,
    sample_sets = make_sample_set_index(
      sample_meta
    ),
    longitudinal = list(
      meta = longitudinal_meta,
      subjects =
        longitudinal_subjects,
      sample_ids =
        longitudinal_meta$SampleID,
      status = ifelse(
        length(
          longitudinal_subjects
        ) >= min_paired_subjects,
        "ready",
        "exploratory"
      )
    ),
    summary = pair_summary,
    status = ifelse(
      nrow(sample_meta) >=
        min_pairwise_samples,
      "ready",
      "insufficient sample-level overlap"
    )
  )
  
  pair_object[[paste0(
    block_x_name,
    "_table"
  )]] <- x_table
  
  pair_object[[paste0(
    block_y_name,
    "_table"
  )]] <- y_table
  
  pair_object[[paste0(
    block_x_name,
    "_dist"
  )]] <- x_dist
  
  pair_object[[paste0(
    block_y_name,
    "_dist"
  )]] <- y_dist
  
  pair_object[[paste0(
    block_x_name,
    "_pcoa"
  )]] <- x_pcoa$object
  
  pair_object[[paste0(
    block_y_name,
    "_pcoa"
  )]] <- y_pcoa$object
  
  pair_object[[paste0(
    block_x_name,
    "_pcoa_scores"
  )]] <- x_pcoa$scores
  
  pair_object[[paste0(
    block_y_name,
    "_pcoa_scores"
  )]] <- y_pcoa$scores
  
  pair_object[[paste0(
    block_x_name,
    "_delta_table"
  )]] <- delta_x$table
  
  pair_object[[paste0(
    block_y_name,
    "_delta_table"
  )]] <- delta_y$table
  
  pair_object[[paste0(
    block_x_name,
    "_delta_pcoa_scores"
  )]] <- delta_x$pcoa_scores
  
  pair_object[[paste0(
    block_y_name,
    "_delta_pcoa_scores"
  )]] <- delta_y$pcoa_scores
  
  pair_object
}


#=================================================================#
# 2. Load metadata, taxonomy, metabolites, KO, and RNA-seq inputs
#=================================================================#

load(
  "input/metabolite_metadata.RData"
)

required_objects <- c(
  "m",
  "s",
  "metabolite"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1)
  )
]

if (length(missing_objects) > 0) {
  stop(
    "Missing objects in metabolite_metadata.RData: ",
    paste(
      missing_objects,
      collapse = ", "
    )
  )
}

if (!file.exists(ko_file)) {
  stop(
    "KO abundance file was not found: ",
    ko_file
  )
}

rna <- read.csv(
  "host_RNAseq/TNT_Expression_Profile.GRCh38.gene.csv",
  check.names = FALSE
)

rna_meta <- read.csv(
  "host_RNAseq/meta_RNA.csv",
  check.names = FALSE
)

rna_map <- read.csv(
  "host_RNAseq/chart_numb_matching.csv",
  check.names = FALSE
)


#=================================================================#
# 3. Species abundance matrix
#=================================================================#

species_sample_cols <- grep(
  "^Sample_",
  colnames(s),
  value = TRUE
)

if (
  !"Species" %in% colnames(s) ||
  length(species_sample_cols) == 0
) {
  stop(
    "Species names or sample columns were not found."
  )
}

species_raw_all <- s %>%
  dplyr::select(
    Species,
    dplyr::all_of(
      species_sample_cols
    )
  ) %>%
  dplyr::filter(
    !is.na(Species),
    Species != ""
  ) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(
        species_sample_cols
      ),
      ~ suppressWarnings(
        as.numeric(
          as.character(.x)
        )
      )
    )
  ) %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(
        species_sample_cols
      ),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  tibble::column_to_rownames(
    "Species"
  ) %>%
  as.matrix() %>%
  base::t()

storage.mode(species_raw_all) <-
  "numeric"

meta_species <- m %>%
  dplyr::transmute(
    OmicSampleID =
      as.character(SampleID),
    SubjectID =
      as.character(SubjectID),
    Timepoint =
      standardize_timepoint(TNT),
    TRG_1 =
      as.character(TRG_1),
    TRG_plot =
      make_trg_plot(TRG_1)
  ) %>%
  dplyr::filter(
    OmicSampleID %in%
      rownames(species_raw_all),
    Timepoint %in% time_levels,
    !is.na(TRG_plot)
  ) %>%
  dplyr::distinct(
    OmicSampleID,
    .keep_all = TRUE
  )

meta_species <- meta_species %>%
  dplyr::filter(
    OmicSampleID %in%
      rownames(species_raw_all)[
        rowSums(
          species_raw_all,
          na.rm = TRUE
        ) > 0
      ]
  )

species_raw_all <- species_raw_all[
  meta_species$OmicSampleID,
  ,
  drop = FALSE
]

species_prevalence_all <- colMeans(
  species_raw_all > 0,
  na.rm = TRUE
)

species_keep_global <- names(
  species_prevalence_all
)[
  species_prevalence_all >=
    species_prev
]

if (length(species_keep_global) == 0) {
  stop("No species feature passed the global prevalence filter.")
}

species_relative_for_rank <-
  species_raw_all[
    ,
    species_keep_global,
    drop = FALSE
  ]

species_relative_for_rank <- sweep(
  species_relative_for_rank,
  1,
  rowSums(
    species_relative_for_rank,
    na.rm = TRUE
  ),
  "/"
)

if (
  !is.null(species_top_n) &&
  ncol(
    species_relative_for_rank
  ) > species_top_n
) {
  species_variance_all <- apply(
    species_relative_for_rank,
    2,
    stats::var,
    na.rm = TRUE
  )
  
  species_keep_global <- names(
    sort(
      species_variance_all,
      decreasing = TRUE
    )
  )[
    seq_len(species_top_n)
  ]
}

species_raw_all <- species_raw_all[
  ,
  species_keep_global,
  drop = FALSE
]

rm(
  species_relative_for_rank
)


#=================================================================#
# 4. KO abundance matrix and annotation
#=================================================================#

ko <- readr::read_tsv(
  ko_file,
  show_col_types = FALSE
) %>%
  dplyr::rename(
    GeneFamily_raw =
      `# Gene Family`
  ) %>%
  tidyr::separate(
    GeneFamily_raw,
    into = c(
      "KO",
      "Taxon"
    ),
    sep = "\\|",
    fill = "right",
    extra = "merge",
    remove = FALSE
  ) %>%
  dplyr::mutate(
    Taxon = dplyr::na_if(
      Taxon,
      ""
    ),
    Feature_level = ifelse(
      is.na(Taxon),
      "KO_total",
      "KO_stratified"
    )
  )

ko_sample_cols <- grep(
  "^Sample_",
  colnames(ko),
  value = TRUE
)

if (length(ko_sample_cols) == 0) {
  stop(
    "No KO sample columns starting with Sample_ were found."
  )
}

ko_total_clean <- ko %>%
  dplyr::filter(
    Feature_level == "KO_total"
  ) %>%
  dplyr::mutate(
    KO_clean = KO %>%
      stringr::str_trim() %>%
      stringr::str_remove(
        "^['\"]"
      ) %>%
      stringr::str_remove(
        "['\"]$"
      ) %>%
      stringr::str_squish(),
    KO_number =
      stringr::str_extract(
        KO_clean,
        "^K\\d{5}"
      ),
    EC_number =
      stringr::str_match(
        KO_clean,
        "\\[EC:([^\\]]+)\\]"
      )[, 2],
    Protein_name = KO_clean %>%
      stringr::str_remove(
        "^K\\d{5}:\\s*"
      ) %>%
      stringr::str_remove(
        "\\s*\\[EC:[^\\]]+\\]\\s*$"
      ) %>%
      stringr::str_squish(),
    KO_text =
      stringr::str_to_lower(
        KO_clean
      ),
    exclude_reason =
      dplyr::case_when(
        stringr::str_detect(
          KO_clean,
          "^UNGROUPED|^NO_NAME|^UNMAPPED"
        ) ~ "non-informative",
        stringr::str_detect(
          KO_text,
          paste0(
            "large subunit ribosomal protein|",
            "small subunit ribosomal protein"
          )
        ) ~ "ribosomal protein",
        TRUE ~ NA_character_
      )
  ) %>%
  dplyr::filter(
    is.na(exclude_reason),
    !is.na(KO_number)
  )

ko_annotation <- ko_total_clean %>%
  dplyr::select(
    KO_number,
    Protein_name,
    EC_number,
    KO_clean
  ) %>%
  dplyr::distinct(
    KO_number,
    .keep_all = TRUE
  )

ko_raw_all <- ko_total_clean %>%
  dplyr::select(
    KO_number,
    dplyr::all_of(
      ko_sample_cols
    )
  ) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(
        ko_sample_cols
      ),
      ~ suppressWarnings(
        as.numeric(
          as.character(.x)
        )
      )
    )
  ) %>%
  dplyr::group_by(
    KO_number
  ) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(
        ko_sample_cols
      ),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  tibble::column_to_rownames(
    "KO_number"
  ) %>%
  as.matrix() %>%
  base::t()

storage.mode(ko_raw_all) <-
  "numeric"

meta_ko <- m %>%
  dplyr::transmute(
    OmicSampleID =
      as.character(SampleID),
    SubjectID =
      as.character(SubjectID),
    Timepoint =
      standardize_timepoint(TNT),
    TRG_1 =
      as.character(TRG_1),
    TRG_plot =
      make_trg_plot(TRG_1)
  ) %>%
  dplyr::filter(
    OmicSampleID %in%
      rownames(ko_raw_all),
    Timepoint %in% time_levels,
    !is.na(TRG_plot)
  ) %>%
  dplyr::distinct(
    OmicSampleID,
    .keep_all = TRUE
  )

meta_ko <- meta_ko %>%
  dplyr::filter(
    OmicSampleID %in%
      rownames(ko_raw_all)[
        rowSums(
          ko_raw_all,
          na.rm = TRUE
        ) > 0
      ]
  )

ko_raw_all <- ko_raw_all[
  meta_ko$OmicSampleID,
  ,
  drop = FALSE
]

ko_prevalence_all <- colMeans(
  ko_raw_all > 0,
  na.rm = TRUE
)

ko_keep_global <- names(
  ko_prevalence_all
)[
  ko_prevalence_all >= ko_prev
]

if (length(ko_keep_global) == 0) {
  stop("No KO feature passed the global prevalence filter.")
}

ko_relative_for_rank <- ko_raw_all[
  ,
  ko_keep_global,
  drop = FALSE
]

ko_relative_for_rank <- sweep(
  ko_relative_for_rank,
  1,
  rowSums(
    ko_relative_for_rank,
    na.rm = TRUE
  ),
  "/"
)

if (
  !is.null(ko_top_n) &&
  ncol(
    ko_relative_for_rank
  ) > ko_top_n
) {
  ko_variance_all <- apply(
    ko_relative_for_rank,
    2,
    stats::var,
    na.rm = TRUE
  )
  
  ko_keep_global <- names(
    sort(
      ko_variance_all,
      decreasing = TRUE
    )
  )[
    seq_len(ko_top_n)
  ]
}

ko_raw_all <- ko_raw_all[
  ,
  ko_keep_global,
  drop = FALSE
]

ko_annotation <- ko_annotation %>%
  dplyr::filter(
    KO_number %in%
      ko_keep_global
  )

write.csv(
  ko_annotation,
  "results/coherence_input/KO_feature_annotation.csv",
  row.names = FALSE
)

rm(
  ko_relative_for_rank
)


#=================================================================#
# 5. Metabolite abundance matrix
#=================================================================#

metabolite_start_idx <- match(
  metabolite_start_name,
  colnames(metabolite)
)

metabolite_end_idx <- match(
  metabolite_end_name,
  colnames(metabolite)
)

if (
  is.na(metabolite_start_idx) ||
  is.na(metabolite_end_idx) ||
  metabolite_start_idx >
  metabolite_end_idx
) {
  stop(
    "Could not identify the metabolite column range."
  )
}

metabolite_cols <- colnames(
  metabolite
)[
  metabolite_start_idx:
    metabolite_end_idx
]

if ("Time" %in% colnames(metabolite)) {
  metabolite_time_col <- "Time"
} else if (
  "TNT" %in% colnames(metabolite)
) {
  metabolite_time_col <- "TNT"
} else {
  stop(
    "Neither Time nor TNT was found in metabolite."
  )
}

metabolite_raw_all <- metabolite %>%
  dplyr::select(
    SampleID,
    dplyr::all_of(
      metabolite_cols
    )
  ) %>%
  dplyr::filter(
    !is.na(SampleID),
    SampleID != ""
  ) %>%
  dplyr::distinct(
    SampleID,
    .keep_all = TRUE
  ) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(
        metabolite_cols
      ),
      ~ suppressWarnings(
        as.numeric(
          as.character(.x)
        )
      )
    )
  ) %>%
  tibble::column_to_rownames(
    "SampleID"
  ) %>%
  as.matrix()

storage.mode(metabolite_raw_all) <-
  "numeric"

meta_metabolite <- metabolite %>%
  dplyr::transmute(
    OmicSampleID =
      as.character(SampleID),
    SubjectID =
      as.character(SubjectID),
    Timepoint =
      standardize_timepoint(
        .data[[metabolite_time_col]]
      ),
    TRG_1 =
      as.character(TRG_1),
    TRG_plot =
      make_trg_plot(TRG_1)
  ) %>%
  dplyr::filter(
    OmicSampleID %in%
      rownames(
        metabolite_raw_all
      ),
    Timepoint %in% time_levels,
    !is.na(TRG_plot)
  ) %>%
  dplyr::distinct(
    OmicSampleID,
    .keep_all = TRUE
  )

metabolite_valid_samples <- rownames(
  metabolite_raw_all
)[
  rowMeans(
    is.finite(
      metabolite_raw_all
    )
  ) >=
    metabolite_min_finite_fraction
]

meta_metabolite <- meta_metabolite %>%
  dplyr::filter(
    OmicSampleID %in%
      metabolite_valid_samples
  )

metabolite_raw_all <-
  metabolite_raw_all[
    meta_metabolite$OmicSampleID,
    ,
    drop = FALSE
  ]

metabolite_qc_result <-
  qc_metabolite_detection_floor(
    metabolite_raw_all,
    min_finite_fraction =
      metabolite_min_finite_fraction,
    min_detected_fraction =
      metabolite_min_detected_fraction,
    floor_min_repeats =
      metabolite_floor_min_repeats,
    floor_relative_tolerance =
      metabolite_floor_relative_tolerance,
    floor_absolute_tolerance =
      metabolite_floor_absolute_tolerance
  )

metabolite_raw_all <-
  metabolite_qc_result$matrix

metabolite_qc_summary <-
  metabolite_qc_result$summary

metabolite_keep_global <-
  metabolite_qc_result$retained_features

write.csv(
  metabolite_qc_summary,
  "results/coherence_input/metabolite_detection_floor_QC.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    metabolite =
      metabolite_keep_global,
    stringsAsFactors = FALSE
  ),
  "results/coherence_input/metabolite_features_retained_after_QC.csv",
  row.names = FALSE
)

message(
  "Metabolite QC: retained ",
  length(
    metabolite_keep_global
  ),
  " of ",
  nrow(
    metabolite_qc_summary
  ),
  " metabolites."
)


#=================================================================#
# 6. Host RNA-seq VST expression matrix
#=================================================================#

m_rna_all <- m %>%
  dplyr::mutate(
    SNU_AL_ID =
      as.character(SNU_ID),
    Timepoint =
      standardize_timepoint(TNT)
  ) %>%
  dplyr::inner_join(
    rna_map %>%
      dplyr::select(
        -dplyr::any_of(
          c(
            "sex",
            "age"
          )
        )
      ),
    by = "SNU_AL_ID"
  ) %>%
  dplyr::mutate(
    Chart_matching = ifelse(
      Timepoint == "Before",
      paste0(
        Chart_numb,
        "_1"
      ),
      paste0(
        Chart_numb,
        "_2"
      )
    )
  ) %>%
  dplyr::inner_join(
    rna_meta %>%
      dplyr::select(
        RNA_sample_id,
        Chart_numb,
        timepoint
      ) %>%
      dplyr::mutate(
        RNA_Timepoint =
          standardize_timepoint(
            timepoint
          ),
        Chart_matching = ifelse(
          RNA_Timepoint ==
            "Before",
          paste0(
            Chart_numb,
            "_1"
          ),
          paste0(
            Chart_numb,
            "_2"
          )
        )
      ) %>%
      dplyr::select(
        -Chart_numb
      ),
    by = "Chart_matching"
  ) %>%
  dplyr::mutate(
    SubjectID =
      as.character(SubjectID),
    RNA_Timepoint =
      standardize_timepoint(
        timepoint
      ),
    TRG_plot =
      make_trg_plot(TRG_1)
  ) %>%
  dplyr::filter(
    Timepoint %in% time_levels,
    RNA_Timepoint %in% time_levels,
    Timepoint == RNA_Timepoint,
    !is.na(TRG_plot)
  ) %>%
  dplyr::distinct(
    RNA_sample_id,
    .keep_all = TRUE
  )

rna_cols_all <- m_rna_all %>%
  dplyr::distinct(
    RNA_sample_id
  ) %>%
  dplyr::mutate(
    count_col = paste0(
      RNA_sample_id,
      "_Read_Count"
    )
  ) %>%
  dplyr::filter(
    count_col %in%
      colnames(rna)
  )

if (nrow(rna_cols_all) == 0) {
  stop(
    "No mapped host RNA-seq sample was found."
  )
}

cnt_rna_all <- rna %>%
  dplyr::mutate(
    Gene = dplyr::case_when(
      !is.na(Gene_Symbol) &
        Gene_Symbol != "" ~
        as.character(
          Gene_Symbol
        ),
      !is.na(Gene_ID) &
        Gene_ID != "" ~
        as.character(
          Gene_ID
        ),
      TRUE ~
        as.character(
          Transcript_ID
        )
    )
  ) %>%
  dplyr::filter(
    !is.na(Gene),
    Gene != ""
  ) %>%
  dplyr::select(
    Gene,
    dplyr::all_of(
      rna_cols_all$count_col
    )
  ) %>%
  dplyr::group_by(
    Gene
  ) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(
        rna_cols_all$count_col
      ),
      ~ sum(
        suppressWarnings(
          as.numeric(.x)
        ),
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  ) %>%
  tibble::column_to_rownames(
    "Gene"
  ) %>%
  as.matrix()

cnt_rna_all <- round(
  cnt_rna_all
)

storage.mode(cnt_rna_all) <-
  "integer"

colnames(cnt_rna_all) <-
  rna_cols_all$RNA_sample_id[
    match(
      colnames(cnt_rna_all),
      rna_cols_all$count_col
    )
  ]

m_rna_all <- m_rna_all %>%
  dplyr::filter(
    RNA_sample_id %in%
      colnames(cnt_rna_all)
  ) %>%
  dplyr::arrange(
    match(
      RNA_sample_id,
      colnames(cnt_rna_all)
    )
  )

cnt_rna_all <- cnt_rna_all[
  ,
  m_rna_all$RNA_sample_id,
  drop = FALSE
]

host_count_keep <- rowSums(
  cnt_rna_all >= host_min_count,
  na.rm = TRUE
) >=
  ceiling(
    host_min_sample_fraction *
      ncol(cnt_rna_all)
  )

if (sum(host_count_keep) == 0) {
  stop("No host gene passed the count-based expression filter.")
}

host_col_data <- data.frame(
  intercept = rep(
    1L,
    ncol(cnt_rna_all)
  ),
  row.names =
    colnames(cnt_rna_all),
  check.names = FALSE
)

host_dds <-
  DESeq2::DESeqDataSetFromMatrix(
    countData =
      cnt_rna_all[
        host_count_keep,
        ,
        drop = FALSE
      ],
    colData =
      host_col_data,
    design = ~ 1
  )

host_vst <-
  DESeq2::varianceStabilizingTransformation(
    host_dds,
    blind = TRUE
  )

host_vst_all <- base::t(
  SummarizedExperiment::assay(
    host_vst
  )
)

storage.mode(host_vst_all) <-
  "numeric"

host_gene_mean <- colMeans(
  host_vst_all,
  na.rm = TRUE
)

host_gene_variance <- apply(
  host_vst_all,
  2,
  stats::var,
  na.rm = TRUE
)

host_gene_variance_all <- data.frame(
  Gene = names(
    host_gene_variance
  ),
  mean_vst = as.numeric(
    host_gene_mean[
      names(
        host_gene_variance
      )
    ]
  ),
  variance = as.numeric(
    host_gene_variance
  ),
  stringsAsFactors = FALSE
)

host_gene_variance_all <-
  host_gene_variance_all[
    is.finite(
      host_gene_variance_all$mean_vst
    ) &
      is.finite(
        host_gene_variance_all$variance
      ) &
      host_gene_variance_all$variance > 0,
    ,
    drop = FALSE
  ]

host_gene_variance_all <-
  host_gene_variance_all[
    order(
      -host_gene_variance_all$variance
    ),
    ,
    drop = FALSE
  ]

host_gene_variance_all$variance_rank <-
  seq_len(
    nrow(
      host_gene_variance_all
    )
  )

host_gene_variance_all$log10_variance <- log10(
  host_gene_variance_all$variance
)

host_gene_variance_all$cumulative_variance_fraction <- cumsum(
  host_gene_variance_all$variance
) / sum(
  host_gene_variance_all$variance
)

host_variance_elbow_rank <- calculate_rank_curve_elbow(
  host_gene_variance_all$variance
)

host_cumulative_rank_50 <- which(
  host_gene_variance_all$cumulative_variance_fraction >= 0.50
)[1]

host_cumulative_rank_75 <- which(
  host_gene_variance_all$cumulative_variance_fraction >= 0.75
)[1]

host_cumulative_rank_80 <- which(
  host_gene_variance_all$cumulative_variance_fraction >= 0.80
)[1]

host_cumulative_rank_90 <- which(
  host_gene_variance_all$cumulative_variance_fraction >= 0.90
)[1]

host_cumulative_rank_95 <- which(
  host_gene_variance_all$cumulative_variance_fraction >= 0.95
)[1]

if (
  identical(
    host_selection_mode,
    "cumulative_variance"
  )
) {
  host_top_n <- which(
    host_gene_variance_all$cumulative_variance_fraction >=
      host_primary_cumulative_fraction
  )[1]
  
  if (
    !is.finite(
      host_top_n
    )
  ) {
    stop(
      "The host cumulative-variance threshold could not be calculated.",
      call. = FALSE
    )
  }
} else if (
  identical(
    host_selection_mode,
    "variance_elbow"
  )
) {
  host_top_n <- host_variance_elbow_rank
} else {
  stop(
    "Unknown host_selection_mode: ",
    host_selection_mode,
    call. = FALSE
  )
}

host_sensitivity_cumulative_ranks <- vapply(
  host_sensitivity_cumulative_fractions,
  function(fraction_value) {
    which(
      host_gene_variance_all$cumulative_variance_fraction >=
        fraction_value
    )[1]
  },
  integer(1)
)

host_sensitivity_top_n_actual <- unique(
  c(
    host_sensitivity_fixed_top_n,
    host_sensitivity_cumulative_ranks,
    host_top_n
  )
)



if (host_include_all_qc_genes_in_sensitivity) {
  host_sensitivity_top_n_actual <- unique(
    c(
      host_sensitivity_top_n_actual,
      nrow(
        host_gene_variance_all
      )
    )
  )
}

host_sensitivity_top_n_actual <- sort(
  host_sensitivity_top_n_actual[
    host_sensitivity_top_n_actual >= 2
  ]
)

host_variance_threshold_summary <- data.frame(
  selection_mode = host_selection_mode,
  variance_definition =
    "Sample-wise variance of blind DESeq2 VST expression",
  primary_rule =
    "Smallest highest-variance gene set reaching 50% cumulative gene-wise variance",
  gene_wise_scaling_after_selection =
    "none",
  selected_top_n = host_top_n,
  primary_cumulative_fraction =
    host_primary_cumulative_fraction,
  selected_cumulative_variance_fraction =
    host_gene_variance_all$cumulative_variance_fraction[
      host_top_n
    ],
  variance_elbow_rank = host_variance_elbow_rank,
  cumulative_variance_50_rank = host_cumulative_rank_50,
  cumulative_variance_75_rank = host_cumulative_rank_75,
  cumulative_variance_80_rank = host_cumulative_rank_80,
  cumulative_variance_90_rank = host_cumulative_rank_90,
  cumulative_variance_95_rank = host_cumulative_rank_95,
  n_genes_after_count_and_variance_QC = nrow(
    host_gene_variance_all
  ),
  selected_variance_threshold = host_gene_variance_all$variance[
    host_top_n
  ],
  stringsAsFactors = FALSE
)

write.csv(
  host_gene_variance_all,
  "results/coherence_input/host_gene_VST_variance_all_QC_genes.csv",
  row.names = FALSE
)

write.csv(
  host_variance_threshold_summary,
  "results/coherence_input/host_gene_variance_threshold_summary.csv",
  row.names = FALSE
)

host_variance_marker_df <- data.frame(
  marker = c(
    "Primary",
    "Elbow",
    "50% cumulative",
    "75% cumulative",
    "80% cumulative",
    "90% cumulative",
    "95% cumulative"
  ),
  rank = c(
    host_top_n,
    host_variance_elbow_rank,
    host_cumulative_rank_50,
    host_cumulative_rank_75,
    host_cumulative_rank_80,
    host_cumulative_rank_90,
    host_cumulative_rank_95
  ),
  stringsAsFactors = FALSE
)

host_variance_marker_df <- host_variance_marker_df[
  is.finite(
    host_variance_marker_df$rank
  ) &
    host_variance_marker_df$rank >= 1 &
    host_variance_marker_df$rank <= nrow(
      host_gene_variance_all
    ),
  ,
  drop = FALSE
]

host_variance_marker_df$variance <-
  host_gene_variance_all$variance[
    host_variance_marker_df$rank
  ]

p_host_variance_rank <- ggplot2::ggplot(
  host_gene_variance_all,
  ggplot2::aes(
    x = variance_rank,
    y = variance
  )
) +
  ggplot2::geom_line(
    linewidth = 0.55,
    color = "grey30"
  ) +
  ggplot2::geom_vline(
    xintercept = host_sensitivity_top_n_actual,
    linewidth = 0.35,
    linetype = "dotted",
    color = "grey72"
  ) +
  ggplot2::geom_vline(
    data = host_variance_marker_df[
      host_variance_marker_df$marker %in% c(
        "Primary",
        "Elbow"
      ),
      ,
      drop = FALSE
    ],
    ggplot2::aes(
      xintercept = rank,
      linetype = marker
    ),
    linewidth = 0.75,
    color = "grey15"
  ) +
  ggplot2::scale_y_log10() +
  ggplot2::scale_linetype_manual(
    values = c(
      Primary = "solid",
      Elbow = "longdash"
    ),
    name = NULL
  ) +
  ggplot2::labs(
    x = "Gene rank by VST variance",
    y = "Gene-wise VST variance",
    caption = paste0(
      "Primary 50% cumulative-variance rank = ",
      host_top_n,
      "; elbow rank = ",
      ifelse(
        is.finite(
          host_variance_elbow_rank
        ),
        host_variance_elbow_rank,
        "NA"
      ),
      ". Dotted lines indicate gene-count sensitivity settings."
    )
  ) +
  ggplot2::theme_classic(
    base_size = 10.5
  ) +
  ggplot2::theme(
    legend.position = "top",
    plot.caption = ggplot2::element_text(
      hjust = 0,
      size = 8.2
    )
  )

ggplot2::ggsave(
  "figures/coherence/Host_RNAseq_VST_gene_variance_rank.svg",
  p_host_variance_rank,
  width = 6.8,
  height = 4.6,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

p_host_variance_cumulative <- ggplot2::ggplot(
  host_gene_variance_all,
  ggplot2::aes(
    x = variance_rank,
    y = cumulative_variance_fraction
  )
) +
  ggplot2::geom_line(
    linewidth = 0.65,
    color = "grey30"
  ) +
  ggplot2::geom_hline(
    yintercept = c(
      0.50,
      0.75,
      0.80,
      0.90,
      0.95
    ),
    linewidth = 0.35,
    linetype = "dotted",
    color = "grey70"
  ) +
  ggplot2::geom_vline(
    xintercept = host_top_n,
    linewidth = 0.75,
    color = "grey15"
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_percent(
      accuracy = 1
    ),
    limits = c(
      0,
      1
    )
  ) +
  ggplot2::labs(
    x = "Number of highest-variance genes retained",
    y = "Cumulative fraction of total VST variance",
    caption = paste0(
      "Ranks reaching 50%, 75%, 80%, 90%, and 95% of total marginal VST variance: ",
      host_cumulative_rank_50,
      ", ",
      host_cumulative_rank_75,
      ", ",
      host_cumulative_rank_80,
      ", ",
      host_cumulative_rank_90,
      ", and ",
      host_cumulative_rank_95,
      "."
    )
  ) +
  ggplot2::theme_classic(
    base_size = 10.5
  ) +
  ggplot2::theme(
    plot.caption = ggplot2::element_text(
      hjust = 0,
      size = 8.2
    )
  )

ggplot2::ggsave(
  "figures/coherence/Host_RNAseq_VST_gene_variance_cumulative.svg",
  p_host_variance_cumulative,
  width = 6.8,
  height = 4.6,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

host_mean_variance_bin <- host_gene_variance_all %>%
  dplyr::mutate(
    mean_bin = dplyr::ntile(
      mean_vst,
      50
    )
  ) %>%
  dplyr::group_by(
    mean_bin
  ) %>%
  dplyr::summarise(
    mean_vst = stats::median(
      mean_vst,
      na.rm = TRUE
    ),
    variance_median = stats::median(
      variance,
      na.rm = TRUE
    ),
    variance_q25 = stats::quantile(
      variance,
      0.25,
      na.rm = TRUE
    ),
    variance_q75 = stats::quantile(
      variance,
      0.75,
      na.rm = TRUE
    ),
    n_genes = dplyr::n(),
    .groups = "drop"
  )

p_host_mean_variance <- ggplot2::ggplot(
  host_mean_variance_bin,
  ggplot2::aes(
    x = mean_vst,
    y = variance_median
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = variance_q25,
      ymax = variance_q75
    ),
    fill = "grey85",
    color = NA
  ) +
  ggplot2::geom_line(
    linewidth = 0.7,
    color = "grey25"
  ) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    x = "Mean DESeq2 VST expression",
    y = "VST variance across samples",
    caption = paste(
      "Points are summarized into 50 equal-count mean-expression bins;",
      "the line and ribbon show the median and interquartile range."
    )
  ) +
  ggplot2::theme_classic(
    base_size = 10.5
  ) +
  ggplot2::theme(
    plot.caption = ggplot2::element_text(
      hjust = 0,
      size = 8.2
    )
  )

ggplot2::ggsave(
  "figures/coherence/Host_RNAseq_VST_mean_variance.svg",
  p_host_mean_variance,
  width = 6.8,
  height = 4.6,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

p_host_variance_diagnostics <- patchwork::wrap_plots(
  p_host_variance_rank,
  p_host_variance_cumulative,
  p_host_mean_variance,
  nrow = 1,
  widths = c(
    1,
    1,
    1
  )
)

ggplot2::ggsave(
  "figures/coherence/Host_RNAseq_VST_variance_diagnostics.svg",
  p_host_variance_diagnostics,
  width = 16.8,
  height = 4.7,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

host_vst_qc_all <- host_vst_all[
  ,
  host_gene_variance_all$Gene,
  drop = FALSE
]

host_gene_selection <- host_gene_variance_all[
  seq_len(
    host_top_n
  ),
  ,
  drop = FALSE
]

host_gene_selection$selected_primary <- TRUE

host_vst_all <- host_vst_qc_all[
  ,
  host_gene_selection$Gene,
  drop = FALSE
]

meta_host <- m_rna_all %>%
  dplyr::transmute(
    OmicSampleID =
      as.character(
        RNA_sample_id
      ),
    SubjectID =
      as.character(
        SubjectID
      ),
    Timepoint,
    TRG_1 =
      as.character(
        TRG_1
      ),
    TRG_plot =
      as.character(
        TRG_plot
      )
  ) %>%
  dplyr::filter(
    OmicSampleID %in%
      rownames(host_vst_all),
    Timepoint %in% time_levels,
    !is.na(TRG_plot)
  ) %>%
  dplyr::distinct(
    OmicSampleID,
    .keep_all = TRUE
  )

host_vst_all <- host_vst_all[
  meta_host$OmicSampleID,
  ,
  drop = FALSE
]

host_vst_qc_all <- host_vst_qc_all[
  meta_host$OmicSampleID,
  ,
  drop = FALSE
]

host_sensitivity_sample_meta <- meta_host %>%
  dplyr::mutate(
    SampleID = paste(
      SubjectID,
      Timepoint,
      sep = "__"
    )
  ) %>%
  dplyr::select(
    SampleID,
    OmicSampleID,
    SubjectID,
    Timepoint,
    TRG_1,
    TRG_plot
  )

rownames(host_vst_qc_all) <-
  host_sensitivity_sample_meta$SampleID

write.csv(
  host_gene_selection,
  "results/coherence_input/host_gene_variance_selection.csv",
  row.names = FALSE
)


#=================================================================#
# 7. Duplicate checks and four modality objects
#=================================================================#

duplicate_subject_timepoint <- dplyr::bind_rows(
  check_subject_timepoint_duplicates(
    meta_species,
    "Species"
  ),
  check_subject_timepoint_duplicates(
    meta_ko,
    "KO"
  ),
  check_subject_timepoint_duplicates(
    meta_metabolite,
    "Metabolite"
  ),
  check_subject_timepoint_duplicates(
    meta_host,
    "Host_RNAseq"
  )
)

write.csv(
  duplicate_subject_timepoint,
  "results/coherence_input/coherence_input_duplicate_subject_timepoint.csv",
  row.names = FALSE
)

if (nrow(duplicate_subject_timepoint) > 0) {
  stop(
    "Duplicated SubjectID-Timepoint records were found. ",
    "Inspect the duplicate CSV."
  )
}

coherence_modalities <- list(
  species = build_modality_object(
    "species",
    meta_species,
    species_raw_all
  ),
  ko = build_modality_object(
    "ko",
    meta_ko,
    ko_raw_all
  ),
  metabolite = build_modality_object(
    "metabolite",
    meta_metabolite,
    metabolite_raw_all
  ),
  host = build_modality_object(
    "host",
    meta_host,
    host_vst_all
  )
)


#=================================================================#
# 8. All six non-redundant pairwise objects
#=================================================================#
#
# The first three pairs constitute the primary biological chain.
# The remaining three are retained for sensitivity and supplementary plots.
#
# Pair objects always contain all SubjectID-Timepoint overlaps for the
# corresponding two data types. Before-After completeness is required only
# for the nested longitudinal object.

coherence_pairs <- list(
  species_ko = build_pair_object(
    "species_ko",
    "species",
    "ko",
    coherence_modalities$species,
    coherence_modalities$ko
  ),
  species_metabolite = build_pair_object(
    "species_metabolite",
    "species",
    "metabolite",
    coherence_modalities$species,
    coherence_modalities$metabolite
  ),
  species_host = build_pair_object(
    "species_host",
    "species",
    "host",
    coherence_modalities$species,
    coherence_modalities$host
  ),
  ko_metabolite = build_pair_object(
    "ko_metabolite",
    "ko",
    "metabolite",
    coherence_modalities$ko,
    coherence_modalities$metabolite
  ),
  ko_host = build_pair_object(
    "ko_host",
    "ko",
    "host",
    coherence_modalities$ko,
    coherence_modalities$host
  ),
  metabolite_host = build_pair_object(
    "metabolite_host",
    "metabolite",
    "host",
    coherence_modalities$metabolite,
    coherence_modalities$host
  )
)

coherence_pair_summary <-
  dplyr::bind_rows(
    lapply(
      coherence_pairs,
      function(x) x$summary
    )
  )

write.csv(
  coherence_pair_summary,
  "results/coherence_input/coherence_all_pair_sample_counts.csv",
  row.names = FALSE
)


#=================================================================#
# 9. Four-block sample-level and longitudinal object
#=================================================================#

four_block_meta <-
  coherence_modalities$species$sample_meta %>%
  dplyr::transmute(
    SubjectID,
    Timepoint,
    SampleID,
    trg_species =
      as.character(TRG_plot)
  ) %>%
  dplyr::inner_join(
    coherence_modalities$ko$sample_meta %>%
      dplyr::transmute(
        SubjectID,
        Timepoint,
        trg_ko =
          as.character(TRG_plot)
      ),
    by = c(
      "SubjectID",
      "Timepoint"
    )
  ) %>%
  dplyr::inner_join(
    coherence_modalities$metabolite$sample_meta %>%
      dplyr::transmute(
        SubjectID,
        Timepoint,
        trg_metabolite =
          as.character(TRG_plot)
      ),
    by = c(
      "SubjectID",
      "Timepoint"
    )
  ) %>%
  dplyr::inner_join(
    coherence_modalities$host$sample_meta %>%
      dplyr::transmute(
        SubjectID,
        Timepoint,
        trg_host =
          as.character(TRG_plot)
      ),
    by = c(
      "SubjectID",
      "Timepoint"
    )
  ) %>%
  dplyr::mutate(
    TRG_plot = dplyr::coalesce(
      trg_species,
      trg_ko,
      trg_metabolite,
      trg_host
    ),
    trg_discordant =
      (
        !is.na(trg_species) &
          !is.na(trg_ko) &
          trg_species != trg_ko
      ) |
      (
        !is.na(trg_species) &
          !is.na(trg_metabolite) &
          trg_species != trg_metabolite
      ) |
      (
        !is.na(trg_species) &
          !is.na(trg_host) &
          trg_species != trg_host
      )
  )

if (any(four_block_meta$trg_discordant)) {
  stop(
    "Discordant TRG labels were found in the four-block overlap."
  )
}

four_block_meta <- four_block_meta %>%
  dplyr::arrange(
    SubjectID,
    match(
      Timepoint,
      time_levels
    )
  ) %>%
  dplyr::mutate(
    TRG_plot = factor(
      TRG_plot,
      levels = c(
        "non_pCR",
        "pCR"
      )
    )
  ) %>%
  dplyr::select(
    SampleID,
    SubjectID,
    Timepoint,
    TRG_plot
  )

four_block_longitudinal_subjects <-
  find_complete_longitudinal_subjects(
    four_block_meta
  )

four_block_longitudinal_meta <-
  four_block_meta %>%
  dplyr::filter(
    SubjectID %in%
      four_block_longitudinal_subjects
  ) %>%
  dplyr::arrange(
    SubjectID,
    match(
      Timepoint,
      time_levels
    )
  )

four_block <- list(
  block_names = c(
    "species",
    "ko",
    "metabolite",
    "host"
  ),
  sample_meta = four_block_meta,
  sample_sets = make_sample_set_index(
    four_block_meta
  ),
  longitudinal = list(
    meta =
      four_block_longitudinal_meta,
    subjects =
      four_block_longitudinal_subjects,
    status = ifelse(
      length(
        four_block_longitudinal_subjects
      ) >= min_paired_subjects,
      "ready",
      "exploratory"
    )
  ),
  species_table =
    coherence_modalities$species$table[
      four_block_meta$SampleID,
      ,
      drop = FALSE
    ],
  ko_table =
    coherence_modalities$ko$table[
      four_block_meta$SampleID,
      ,
      drop = FALSE
    ],
  metabolite_table =
    coherence_modalities$metabolite$table[
      four_block_meta$SampleID,
      ,
      drop = FALSE
    ],
  host_table =
    coherence_modalities$host$table[
      four_block_meta$SampleID,
      ,
      drop = FALSE
    ],
  status = ifelse(
    nrow(four_block_meta) >=
      min_four_block_samples,
    "ready",
    "insufficient four-block overlap"
  )
)

four_block_summary <- data.frame(
  analysis_set =
    "species_ko_metabolite_host",
  n_overlap_samples =
    nrow(four_block_meta),
  n_overlap_baseline_samples =
    sum(
      four_block_meta$Timepoint ==
        "Before"
    ),
  n_overlap_after_rt_samples =
    sum(
      four_block_meta$Timepoint ==
        "Ongoing"
    ),
  n_overlap_subjects =
    dplyr::n_distinct(
      four_block_meta$SubjectID
    ),
  n_longitudinal_subjects =
    length(
      four_block_longitudinal_subjects
    ),
  status = four_block$status,
  stringsAsFactors = FALSE
)

write.csv(
  four_block_summary,
  "results/coherence_input/coherence_four_block_sample_counts.csv",
  row.names = FALSE
)


#=================================================================#
# 10. Save one essential object
#=================================================================#

coherence_data <- list(
  schema_version = "3.5",
  biological_sequence = c(
    "species",
    "ko",
    "metabolite",
    "host"
  ),
  interpretation = paste(
    "Layered concordance of taxonomy, microbial functional potential,",
    "fecal metabolites, and host transcription; not a causal model."
  ),
  matching_principles = list(
    pair_sample_level =
      "All SubjectID-Timepoint overlaps for the relevant pair.",
    pair_longitudinal =
      "Pair-specific subjects with both Before and Ongoing.",
    four_block_sample_level =
      "All SubjectID-Timepoint overlaps across four layers.",
    four_block_longitudinal =
      "Subjects with both timepoints across all four layers."
  ),
  analysis_settings = list(
    time_levels = time_levels,
    species_prev = species_prev,
    species_top_n = species_top_n,
    ko_prev = ko_prev,
    ko_top_n = ko_top_n,
    metabolite_transform =
      metabolite_transform,
    metabolite_min_finite_fraction =
      metabolite_min_finite_fraction,
    metabolite_min_detected_fraction =
      metabolite_min_detected_fraction,
    metabolite_floor_min_repeats =
      metabolite_floor_min_repeats,
    metabolite_zero_is_detection_floor = TRUE,
    host_selection_mode =
      host_selection_mode,
    host_primary_cumulative_fraction =
      host_primary_cumulative_fraction,
    host_primary_top_n =
      host_top_n,
    host_top_n = host_top_n,
    host_variance_elbow_rank =
      host_variance_elbow_rank,
    host_sensitivity_top_n =
      host_sensitivity_top_n_actual,
    pcoa_max_axes = pcoa_max_axes
  ),
  feature_annotation = list(
    ko = ko_annotation,
    metabolite_qc =
      metabolite_qc_summary,
    host = host_gene_selection,
    host_variance_all =
      host_gene_variance_all,
    host_variance_threshold =
      host_variance_threshold_summary
  ),
  host_sensitivity = list(
    vst_qc_all =
      host_vst_qc_all,
    sample_meta =
      host_sensitivity_sample_meta,
    gene_variance =
      host_gene_variance_all,
    candidate_top_n =
      host_sensitivity_top_n_actual,
    primary_top_n =
      host_top_n,
    primary_cumulative_fraction =
      host_primary_cumulative_fraction,
    variance_elbow_rank =
      host_variance_elbow_rank,
    cumulative_rank_50 =
      host_cumulative_rank_50,
    cumulative_rank_75 =
      host_cumulative_rank_75,
    cumulative_rank_80 =
      host_cumulative_rank_80,
    cumulative_rank_90 =
      host_cumulative_rank_90,
    cumulative_rank_95 =
      host_cumulative_rank_95,
    host_feature_scaling = "none"
  ),
  modalities = coherence_modalities,
  pairs = coherence_pairs,
  pair_summary = coherence_pair_summary,
  four_block = four_block,
  four_block_summary = four_block_summary
)

feature_count_after_qc <- data.frame(
  block = c(
    "species",
    "ko",
    "metabolite",
    "host"
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
      nrow(
        x$sample_meta
      )
    },
    integer(1)
  ),
  n_features_used_for_distance = vapply(
    coherence_data$modalities[
      c(
        "species",
        "ko",
        "metabolite",
        "host"
      )
    ],
    function(x) {
      as.integer(
        x$n_features
      )
    },
    integer(1)
  ),
  feature_selection_cap = c(
    ifelse(
      is.null(
        species_top_n
      ),
      "none",
      as.character(
        species_top_n
      )
    ),
    ifelse(
      is.null(
        ko_top_n
      ),
      "none",
      as.character(
        ko_top_n
      )
    ),
    "not applicable",
    ifelse(
      is.null(
        host_top_n
      ),
      "none",
      as.character(
        host_top_n
      )
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  feature_count_after_qc,
  "results/coherence_input/coherence_feature_counts_after_QC.csv",
  row.names = FALSE
)


save(
  coherence_data,
  file =
    "input/coherence_species_ko_metabolite_host_data_metabolite50_host50pct_no_scaling.RData"
)

message(
  "\nSaved: ",
  "input/coherence_species_ko_metabolite_host_data_metabolite50_host50pct_no_scaling.RData"
)

print(
  coherence_data$pair_summary
)

print(
  coherence_data$four_block_summary
)

rm(
  list = setdiff(
    ls(),
    "coherence_data"
  )
)

gc()
