#-----------------------------------------------------------------#
#
# Reused utilities for four-layer coherence analyses
#
# Layers:
#   Species -> Microbial function (KO) -> Metabolite -> Host RNA-seq
#
# Expected input:
#   input/coherence_species_ko_metabolite_host_data.RData
#
# The input RData must contain:
#   coherence_data
#
# Only functions reused by at least two downstream scripts are retained here.
#
#-----------------------------------------------------------------#

options(stringsAsFactors = FALSE)

coherence_load_data <- function(
  path = "input/coherence_species_ko_metabolite_host_data.RData"
) {
  if (!file.exists(path)) {
    stop(
      "Input file was not found: ",
      path,
      call. = FALSE
    )
  }

  input_env <- new.env(parent = emptyenv())
  load(path, envir = input_env)

  if (!exists(
    "coherence_data",
    envir = input_env,
    inherits = FALSE
  )) {
    stop(
      "`coherence_data` was not found in: ",
      path,
      call. = FALSE
    )
  }

  coherence_data <- input_env$coherence_data

  required_pairs <- c(
    "species_ko",
    "ko_metabolite",
    "metabolite_host"
  )

  if (
    is.null(coherence_data$pairs) ||
    !all(required_pairs %in% names(coherence_data$pairs))
  ) {
    stop(
      "The input does not contain the three primary chain pairs: ",
      paste(required_pairs, collapse = ", "),
      call. = FALSE
    )
  }

  coherence_data
}


coherence_get_pair <- function(coherence_data, pair_name) {
  if (!pair_name %in% names(coherence_data$pairs)) {
    stop(
      "Unknown pair name: ",
      pair_name,
      call. = FALSE
    )
  }

  pair_object <- coherence_data$pairs[[pair_name]]

  if (!identical(pair_object$status, "ready")) {
    stop(
      "Pair object is not ready: ",
      pair_name,
      " [status: ",
      pair_object$status,
      "]",
      call. = FALSE
    )
  }

  pair_object
}


coherence_score_matrix <- function(
  pair_object,
  block_name,
  delta = FALSE
) {
  object_name <- paste0(
    block_name,
    if (delta) "_delta_pcoa_scores" else "_pcoa_scores"
  )

  if (
    !object_name %in% names(pair_object) ||
    is.null(pair_object[[object_name]])
  ) {
    stop(
      "Score object was not found: ",
      pair_object$pair_name,
      "$",
      object_name,
      call. = FALSE
    )
  }

  score_df <- as.data.frame(
    pair_object[[object_name]],
    check.names = FALSE
  )

  if ("SampleID" %in% colnames(score_df)) {
    rownames(score_df) <- as.character(
      score_df$SampleID
    )
  }

  axis_cols <- grep(
    "^Axis",
    colnames(score_df),
    value = TRUE
  )

  if (length(axis_cols) < 2) {
    stop(
      object_name,
      " contains fewer than two ordination axes.",
      call. = FALSE
    )
  }

  score_matrix <- as.matrix(
    score_df[, axis_cols, drop = FALSE]
  )
  storage.mode(score_matrix) <- "numeric"

  score_matrix <- score_matrix[
    ,
    colSums(is.finite(score_matrix)) ==
      nrow(score_matrix),
    drop = FALSE
  ]

  if (ncol(score_matrix) < 2) {
    stop(
      object_name,
      " contains fewer than two complete axes.",
      call. = FALSE
    )
  }

  score_matrix
}


coherence_table_matrix <- function(
  pair_object,
  block_name,
  delta = FALSE
) {
  object_name <- paste0(
    block_name,
    if (delta) "_delta_table" else "_table"
  )

  if (
    !object_name %in% names(pair_object) ||
    is.null(pair_object[[object_name]])
  ) {
    stop(
      "Analysis table was not found: ",
      pair_object$pair_name,
      "$",
      object_name,
      call. = FALSE
    )
  }

  x <- as.matrix(pair_object[[object_name]])
  storage.mode(x) <- "numeric"

  x <- x[
    ,
    apply(x, 2, stats::var, na.rm = TRUE) > 0,
    drop = FALSE
  ]

  if (ncol(x) < 2) {
    stop(
      object_name,
      " contains fewer than two non-constant features.",
      call. = FALSE
    )
  }

  x
}


coherence_subset_ids <- function(pair_object, subset_name) {
  meta <- pair_object$sample_meta

  if (subset_name == "All") {
    return(as.character(meta$SampleID))
  }

  if (subset_name == "Baseline") {
    return(
      as.character(
        meta$SampleID[meta$Timepoint == "Before"]
      )
    )
  }

  if (subset_name == "After RT") {
    return(
      as.character(
        meta$SampleID[meta$Timepoint == "Ongoing"]
      )
    )
  }

  if (subset_name == "pCR") {
    return(
      as.character(
        meta$SampleID[meta$TRG_plot == "pCR"]
      )
    )
  }

  if (subset_name == "non-pCR") {
    return(
      as.character(
        meta$SampleID[meta$TRG_plot == "non_pCR"]
      )
    )
  }

  stop(
    "Unknown subset: ",
    subset_name,
    call. = FALSE
  )
}


coherence_p_to_star <- function(p) {
  ifelse(
    is.na(p),
    "",
    ifelse(
      p < 0.001,
      "***",
      ifelse(
        p < 0.01,
        "**",
        ifelse(p < 0.05, "*", "")
      )
    )
  )
}


coherence_primary_pairs <- data.frame(
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
    "Species vs Microbial function",
    "Microbial function vs Metabolite",
    "Metabolite vs Host RNA-seq"
  ),
  stringsAsFactors = FALSE
)

coherence_layer_labels <- c(
  species = "Species",
  ko = "Microbial function",
  metabolite = "Metabolite",
  host = "Host RNA-seq"
)

coherence_layer_shapes <- c(
  Species = 21,
  `Microbial function` = 24,
  Metabolite = 22,
  `Host RNA-seq` = 23
)

coherence_group_colors <- c(
  pCR = "#56A99B",
  non_pCR = "#D47472"
)

coherence_group_colors_faint <- c(
  pCR = "#B9D9D2",
  non_pCR = "#E8C0BE"
)

coherence_block_colors <- c(
  Species = "#AFC5D8",
  `Microbial function` = "#C6B8D7",
  Metabolite = "#BFD6C5",
  `Host RNA-seq` = "#D9C1B3"
)
