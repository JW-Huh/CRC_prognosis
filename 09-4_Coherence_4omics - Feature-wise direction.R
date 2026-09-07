#-----------------------------------------------------------------#
#
# Four-layer feature-wise directional effect analysis
#
# Layers:
#   Species
#   Microbial function (KO)
#   Metabolite
#   Host RNA-seq
#
# For each feature and timepoint:
#   Hedges' g = standardized pCR - non-pCR difference
#
# Display candidates:
#   - same effect direction at Baseline and After RT
#   - ranked by the smaller absolute Hedges' g across the two timepoints
#
# Circle size:
#   |Hedges' g|
#
# Circle fill:
#   signed -log10(P); green indicates higher in pCR and red indicates
#   higher in non-pCR.
#
#-----------------------------------------------------------------#

rm(list = ls())
options(stringsAsFactors = FALSE)
setwd("D:/2-연구/2-CRC metagenomics/")

source(
  "7-11.0 Coherence_4omics - common utilities.R"
)

required_packages <- c(
  "ggplot2",
  "svglite",
  "scales"
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
    "Install: ",
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

coherence_data <- coherence_load_data()

feature_display_n <- c(
  species = 8,
  ko = 8,
  metabolite = 6,
  host = 10
)

minimum_abs_hedges_g <- 0.50

hedges_g <- function(
  x_pcr,
  x_non
) {
  x_pcr <- x_pcr[
    is.finite(x_pcr)
  ]

  x_non <- x_non[
    is.finite(x_non)
  ]

  n_pcr <- length(x_pcr)
  n_non <- length(x_non)

  if (
    n_pcr < 2 ||
    n_non < 2
  ) {
    return(NA_real_)
  }

  pooled_variance <- (
    (n_pcr - 1) *
      stats::var(x_pcr) +
      (n_non - 1) *
        stats::var(x_non)
  ) / (
    n_pcr +
      n_non -
      2
  )

  if (
    !is.finite(
      pooled_variance
    ) ||
    pooled_variance <= 0
  ) {
    return(NA_real_)
  }

  d <- (
    mean(x_pcr) -
      mean(x_non)
  ) /
    sqrt(
      pooled_variance
    )

  correction <- 1 -
    3 /
      (
        4 *
          (
            n_pcr +
              n_non
          ) -
          9
      )

  correction * d
}


safe_wilcox <- function(
  x_pcr,
  x_non
) {
  x_pcr <- x_pcr[
    is.finite(x_pcr)
  ]

  x_non <- x_non[
    is.finite(x_non)
  ]

  if (
    length(x_pcr) < 2 ||
    length(x_non) < 2 ||
    length(
      unique(
        c(
          x_pcr,
          x_non
        )
      )
    ) < 2
  ) {
    return(NA_real_)
  }

  tryCatch(
    stats::wilcox.test(
      x_pcr,
      x_non,
      exact = FALSE
    )$p.value,
    error = function(e) {
      NA_real_
    }
  )
}


feature_results <- data.frame()

for (
  block_name in c(
    "species",
    "ko",
    "metabolite",
    "host"
  )
) {
  modality_object <-
    coherence_data$modalities[[block_name]]

  x <- as.matrix(
    modality_object$table
  )

  meta <- modality_object$sample_meta[
    match(
      rownames(x),
      modality_object$sample_meta$SampleID
    ),
    ,
    drop = FALSE
  ]

  for (
    timepoint_value in c(
      "Before",
      "Ongoing"
    )
  ) {
    idx <- meta$Timepoint ==
      timepoint_value

    x_time <- x[
      idx,
      ,
      drop = FALSE
    ]

    group_time <-
      as.character(
        meta$TRG_plot[idx]
      )

    g_value <- vapply(
      seq_len(ncol(x_time)),
      function(j) {
        hedges_g(
          x_time[
            group_time == "pCR",
            j
          ],
          x_time[
            group_time == "non_pCR",
            j
          ]
        )
      },
      numeric(1)
    )

    p_value <- vapply(
      seq_len(ncol(x_time)),
      function(j) {
        safe_wilcox(
          x_time[
            group_time == "pCR",
            j
          ],
          x_time[
            group_time == "non_pCR",
            j
          ]
        )
      },
      numeric(1)
    )

    feature_results <- rbind(
      feature_results,
      data.frame(
        block = block_name,
        feature =
          colnames(x_time),
        timepoint =
          ifelse(
            timepoint_value ==
              "Before",
            "Baseline",
            "After RT"
          ),
        n_pCR =
          sum(
            group_time ==
              "pCR"
          ),
        n_non_pCR =
          sum(
            group_time ==
              "non_pCR"
          ),
        hedges_g = g_value,
        p = p_value,
        stringsAsFactors = FALSE
      )
    )
  }
}

feature_results$FDR <-
  ave(
    feature_results$p,
    interaction(
      feature_results$block,
      feature_results$timepoint
    ),
    FUN = function(x) {
      stats::p.adjust(
        x,
        method = "BH"
      )
    }
  )

baseline_results <- feature_results[
  feature_results$timepoint ==
    "Baseline",
  ,
  drop = FALSE
]

after_results <- feature_results[
  feature_results$timepoint ==
    "After RT",
  ,
  drop = FALSE
]

names(
  baseline_results
)[
  names(
    baseline_results
  ) %in%
    c(
      "n_pCR",
      "n_non_pCR",
      "hedges_g",
      "p",
      "FDR"
    )
] <- paste0(
  names(
    baseline_results
  )[
    names(
      baseline_results
    ) %in%
      c(
        "n_pCR",
        "n_non_pCR",
        "hedges_g",
        "p",
        "FDR"
      )
  ],
  "_baseline"
)

names(
  after_results
)[
  names(
    after_results
  ) %in%
    c(
      "n_pCR",
      "n_non_pCR",
      "hedges_g",
      "p",
      "FDR"
    )
] <- paste0(
  names(
    after_results
  )[
    names(
      after_results
    ) %in%
      c(
        "n_pCR",
        "n_non_pCR",
        "hedges_g",
        "p",
        "FDR"
      )
  ],
  "_after"
)

feature_consistency <- merge(
  baseline_results[
    ,
    setdiff(
      names(
        baseline_results
      ),
      "timepoint"
    ),
    drop = FALSE
  ],
  after_results[
    ,
    setdiff(
      names(
        after_results
      ),
      "timepoint"
    ),
    drop = FALSE
  ],
  by = c(
    "block",
    "feature"
  ),
  all = TRUE
)

feature_consistency$same_direction <-
  is.finite(
    feature_consistency$
      hedges_g_baseline
  ) &
  is.finite(
    feature_consistency$
      hedges_g_after
  ) &
  sign(
    feature_consistency$
      hedges_g_baseline
  ) ==
    sign(
      feature_consistency$
        hedges_g_after
    )

feature_consistency$minimum_abs_g <-
  pmin(
    abs(
      feature_consistency$
        hedges_g_baseline
    ),
    abs(
      feature_consistency$
        hedges_g_after
    ),
    na.rm = TRUE
  )

feature_consistency$mean_abs_g <-
  rowMeans(
    cbind(
      abs(
        feature_consistency$
          hedges_g_baseline
      ),
      abs(
        feature_consistency$
          hedges_g_after
      )
    ),
    na.rm = TRUE
  )

selected_features <- data.frame()

for (
  block_name in names(
    feature_display_n
  )
) {
  block_candidates <-
    feature_consistency[
      feature_consistency$block ==
        block_name &
        feature_consistency$
          same_direction &
        feature_consistency$
          minimum_abs_g >=
            minimum_abs_hedges_g,
      ,
      drop = FALSE
    ]

  block_candidates <-
    block_candidates[
      order(
        -block_candidates$
          minimum_abs_g,
        -block_candidates$
          mean_abs_g
      ),
      ,
      drop = FALSE
    ]

  if (nrow(block_candidates) == 0) {
    block_candidates <-
      feature_consistency[
        feature_consistency$block ==
          block_name &
          feature_consistency$
            same_direction,
        ,
        drop = FALSE
      ]

    block_candidates <-
      block_candidates[
        order(
          -block_candidates$
            minimum_abs_g,
          -block_candidates$
            mean_abs_g
        ),
        ,
        drop = FALSE
      ]
  }

  selected_features <- rbind(
    selected_features,
    head(
      block_candidates,
      feature_display_n[
        block_name
      ]
    )
  )
}

selected_long <- merge(
  feature_results,
  selected_features[
    ,
    c(
      "block",
      "feature",
      "minimum_abs_g",
      "mean_abs_g"
    ),
    drop = FALSE
  ],
  by = c(
    "block",
    "feature"
  ),
  all = FALSE
)

selected_long$feature_label <-
  selected_long$feature

ko_annotation <-
  coherence_data$feature_annotation$ko

ko_label <- ko_annotation$Protein_name[
  match(
    selected_long$feature,
    ko_annotation$KO_number
  )
]

selected_long$feature_label[
  selected_long$block == "ko" &
    !is.na(ko_label) &
    ko_label != ""
] <- ko_label[
  selected_long$block == "ko" &
    !is.na(ko_label) &
    ko_label != ""
]

selected_long$block_label <-
  unname(
    coherence_layer_labels[
      selected_long$block
    ]
  )

selected_long$timepoint <- factor(
  selected_long$timepoint,
  levels = c(
    "After RT",
    "Baseline"
  )
)

selected_long$signed_log10_p <-
  sign(
    selected_long$hedges_g
  ) *
  -log10(
    pmax(
      selected_long$p,
      1e-4
    )
  )

selected_long$significance <- coherence_p_to_star(
  selected_long$p
)

selected_long$feature_label <- factor(
  selected_long$feature_label,
  levels = unique(
    selected_long$feature_label[
      order(
        selected_long$block_label,
        -selected_long$
          minimum_abs_g
      )
    ]
  )
)

p_feature_direction <-
  ggplot2::ggplot(
    selected_long,
    ggplot2::aes(
      x = feature_label,
      y = timepoint
    )
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      size = abs(hedges_g),
      fill = signed_log10_p
    ),
    shape = 21,
    color = "grey25",
    stroke = 0.65
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = significance
    ),
    size = 2.6
  ) +
  ggplot2::facet_grid(
    ~ block_label,
    scales = "free_x",
    space = "free_x"
  ) +
  ggplot2::scale_fill_gradient2(
    low =
      coherence_group_colors[
        "non_pCR"
      ],
    mid = "white",
    high =
      coherence_group_colors[
        "pCR"
      ],
    midpoint = 0,
    breaks = c(
      -2,
      -1,
      0,
      1,
      2
    ),
    labels = c(
      "0.01",
      "0.1",
      "0",
      "0.1",
      "0.01"
    ),
    limits = c(
      -2.5,
      2.5
    ),
    oob = scales::squish,
    name =
      "Directional P-value"
  ) +
  ggplot2::scale_size_continuous(
    range = c(
      2.5,
      8
    ),
    breaks = c(
      0.5,
      1.0,
      1.5
    ),
    limits = c(
      0,
      max(
        1.5,
        abs(
          selected_long$hedges_g
        ),
        na.rm = TRUE
      )
    ),
    name = "|Hedges' g|"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    caption = paste(
      "Features were selected by concordant Hedges' g direction at",
      "Baseline and After RT, ranked by the smaller absolute effect."
    )
  ) +
  ggplot2::theme_classic(
    base_size = 10
  ) +
  ggplot2::theme(
    strip.background =
      ggplot2::element_blank(),
    strip.text =
      ggplot2::element_text(
        face = "bold",
        size = 10.5
      ),
    axis.text.x =
      ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 7.5
      ),
    axis.text.y =
      ggplot2::element_text(
        face = "bold"
      ),
    axis.ticks =
      ggplot2::element_blank(),
    legend.position = "bottom",
    plot.caption =
      ggplot2::element_text(
        hjust = 0,
        size = 8
      )
  )

ggplot2::ggsave(
  "figures/coherence/Feature_direction_4omics_D.svg",
  p_feature_direction,
  width = 11.0,
  height = 4.5,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

write.csv(
  feature_results,
  "results/coherence/Feature_direction_all_results.csv",
  row.names = FALSE
)

write.csv(
  feature_consistency,
  "results/coherence/Feature_direction_consistency_results.csv",
  row.names = FALSE
)

write.csv(
  selected_long,
  "results/coherence/Feature_direction_selected_features.csv",
  row.names = FALSE
)

save(
  feature_results,
  feature_consistency,
  selected_features,
  selected_long,
  p_feature_direction,
  file =
    "results/coherence/Feature_direction_results.RData"
)

rm(
  list = setdiff(
    ls(),
    c(
      "coherence_data",
      "feature_results",
      "feature_consistency",
      "selected_features",
      "selected_long",
      "p_feature_direction"
    )
  )
)

gc()
