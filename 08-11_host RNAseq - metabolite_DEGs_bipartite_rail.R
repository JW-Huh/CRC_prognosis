rm(list = ls())
options(stringsAsFactors = FALSE)

#=================================================================#
# Fig. 5D 9.8.4. De novo metabolite–host DEG association map
#=================================================================#
#
# Purpose:
#   Recalculate metabolite–host associations using all host DEGs rather
#   than restricting the analysis to the Figure 5C gene panel.
#
# Statistical workflow:
#   1. Select a generous host-gene set using DESeq2 P < 0.10 or
#      Wilcoxon P < 0.10, consistent with the upstream axis screening.
#   2. Calculate Spearman correlations for all metabolite–DEG pairs.
#   3. Use nominal Spearman P < 0.10 for exploratory edge eligibility,
#      while retaining P < 0.05, P < 0.01, and FDR < 0.10 as stronger tiers.
#   4. Use ordinary Spearman correlation as the primary association.
#   5. Calculate a secondary TRG-adjusted residual-rank correlation by:
#        - rank-transforming each metabolite and gene,
#        - regressing both rank vectors on TRG_plot across all samples,
#        - correlating the resulting residuals.
#      This is not an average of separate pCR and non-pCR correlations.
#   6. Require leave-one-out sign stability before showing adjusted support.
#   7. Use Panel C genes only as a saved annotation, not as a visual encoding.
#   8. Exclude readthrough/antisense/uncharacterized genes from the main plot.
#   9. Prioritize CRC, treatment-response, barrier, immune-trafficking,
#      and gut immunometabolism genes using literature-guided rules.
#  10. Prioritize differential metabolites and preserve metabolite-class diversity.
#
# Main visual encoding:
#   - Left rail: pCR-enriched host genes
#   - Center rail: selected metabolites
#   - Right rail: non-pCR-enriched host genes
#   - Edge color: ordinary Spearman rho
#   - Edge width: capped -log10 nominal Spearman P
#   - All edges are solid because the primary analysis is Spearman
#   - Diamond at edge midpoint: stable TRG-adjusted residual-rank support
#   - Metabolite node fill: metabolite log2FC, pCR versus non-pCR
#   - Host-gene label color: functional annotation category
#
# Important:
#   The positive DEG direction is assumed to mean higher expression in pCR.
#   Confirm that this matches the contrast used in the upstream DEG table.
#
# Output:
#   figures/Fig5D_metabolite_host_DEG_rail_<analysis_set>.svg
#
#   host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/
#     Fig5D_8.4_all_metabolite_DEG_correlations_<analysis_set>.csv
#     Fig5D_8.4_significant_metabolite_DEG_edges_<analysis_set>.csv
#     Fig5D_8.4_displayed_metabolite_DEG_edges_<analysis_set>.csv
#     Fig5D_8.4_selected_metabolites_<analysis_set>.csv
#     Fig5D_8.4_selected_host_DEGs_<analysis_set>.csv
#     Fig5D_8.4_diagnostics_<analysis_set>.csv
#
#   Fig5D_9.8.4_metabolite_host_DEG_rail_<analysis_set>.RData
#
#=================================================================#

setwd("D:/2-연구/2-CRC metagenomics/")

for (pkg in c(
  "dplyr", "tidyr", "tibble", "stringr",
  "ggplot2", "svglite", "ggnewscale"
)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg, call. = FALSE)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(svglite)
  library(ggnewscale)
})

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create(
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

analysis_set <- "before"
# analysis_set <- "all_available"

# The current screening RData stores the DEG table as host_deg_results.
# Set the expression object explicitly if automatic detection is ambiguous.
host_expression_object_name <- NULL
host_deg_object_name <- "host_deg_results"

# Optional additional RData containing the DEG result table.
# Example:
# host_deg_file <- "host_RNAseq/results/Figure5C_host_DEG_results.RData"
host_deg_file <- NULL

# When a separate DEG table is unavailable, use unique genes from
# met_gene_stats as the upstream host-DEG universe. This fallback is valid
# only when met_gene_stats was generated from an already DEG-filtered set.
allow_met_gene_stats_gene_fallback <- TRUE

# Generous host-gene eligibility for exploratory axis screening.
# A gene is retained when DESeq2 P < 0.10 or Wilcoxon P < 0.10.
# No FDR or absolute log2FC threshold is imposed at this stage.
gene_p_cutoff <- 0.10

# Pairwise association filtering.
# This is an exploratory network. Eligible edges use nominal Spearman
# P < 0.10, while P < 0.05, P < 0.01, and FDR < 0.10 are retained as
# progressively stronger evidence tiers.
spearman_p_cutoff <- 0.10
spearman_moderate_p_cutoff <- 0.05
spearman_strong_p_cutoff <- 0.01
spearman_relaxed_fdr_cutoff <- 0.10
spearman_rho_cutoff <- 0.30

# The adjusted result is a secondary sensitivity analysis.
# It is shown only when the residual-rank correlation retains the same
# direction, has nominal P < 0.05, absolute rho >= 0.25, and at least
# 80% of leave-one-out correlations retain the same sign.
adjusted_p_cutoff <- 0.05
adjusted_rho_cutoff <- 0.25
adjusted_loo_sign_fraction_cutoff <- 0.80

# Covariates used for residualized partial Spearman correlation.
# TRG_plot is retained because metabolites and DEGs were both screened
# against treatment-response status.
adjustment_covariates <- c("TRG_plot")

min_complete_samples <- 8L
maximum_missing_fraction <- 0.20

max_selected_metabolites <- 8L
max_metabolites_per_class <- 2L
max_selected_genes <- 44L
max_edges_per_metabolite <- 8L
max_edges_per_gene <- 2L

# Constant edge width in points. ggplot2 converts this to mm below.
# Change only this value when a thinner or thicker Illustrator appearance
# is required after the full figure is resized.
edge_linewidth_pt <- 0.45

positive_deg_label <- "Higher in pCR"
negative_deg_label <- "Higher in non-pCR"

# Butyrate is prioritized when it has at least one eligible biological edge.
# Lithocholic acid receives a smaller literature-based priority but is not
# forced when the data do not support an eligible association.
manual_priority_metabolites <- c("butyrate")
literature_priority_metabolites <- c("lithocholic acid")

# Explicit main-figure exclusions. These metabolites remain in the complete
# correlation table but cannot enter the displayed network.
manual_metabolite_exclude <- c("kynurenic acid")

# Genes specifically requested or previously retained in the triad review.
# These genes are prioritized only when they have an eligible metabolite
# correlation in the present de novo analysis.
manual_priority_genes <- c(
  "CCR1", "TJP1", "S100A9", "CXCL1", "MARCO",
  "SIPA1", "BATF3", "LGR5", "GADD45B", "BOK",
  "RUNX3", "NSD3", "GPR137B", "SART3",
  paste0("IGFBP", 1:7)
)

# Remove readthrough, antisense, uncharacterized, and manually rejected genes
# from the main figure. They remain in the complete correlation table.
manual_gene_exclude <- c(
  "SENP3-EIF4A1", "DGCR11", "SP2-AS1", "C9ORF16",
  "UCA1", "MALAT1", "NEAT1", "XIST"
)

# Panel B category colors, with an additional treatment-response category.
category_palette <- c(
  "Tumor remodeling/CRC" = "#DF8A82",
  "Treatment response" = "#F4A261",
  "Metabolite signaling" = "#A8C98B",
  "Barrier/mucus/AMP" = "#D9B44A",
  "Innate immunity" = "#9DC7DD",
  "Adaptive immunity" = "#B7A3CE",
  "Gut homing" = "#77C1B5"
)

# Exact category assignments take precedence over regex rules.
manual_gene_category <- tribble(
  ~Gene, ~Functional_category, ~biological_priority,
  "CCR1",   "Gut homing",             4,
  "TJP1",   "Barrier/mucus/AMP",       4,
  "S100A9", "Innate immunity",         4,
  "CXCL1",  "Innate immunity",         4,
  "MARCO",  "Innate immunity",         4,
  "SIPA1",  "Adaptive immunity",       4,
  "BATF3",  "Adaptive immunity",       4,
  "LGR5",   "Tumor remodeling/CRC",    4,
  "GADD45B","Treatment response",      4,
  "BOK",    "Treatment response",      4,
  "RUNX3",  "Adaptive immunity",       4,
  "NSD3",   "Tumor remodeling/CRC",    3,
  "GPR137B","Metabolite signaling",    3,
  "SART3",  "Treatment response",      3,
  "IGFBP1", "Treatment response",      4,
  "IGFBP2", "Treatment response",      4,
  "IGFBP3", "Treatment response",      4,
  "IGFBP4", "Treatment response",      4,
  "IGFBP5", "Treatment response",      4,
  "IGFBP6", "Treatment response",      4,
  "IGFBP7", "Treatment response",      4
)

# Literature-guided rules adapted from the previous metabolite–host
# screening script. These are ranking annotations, not statistical evidence.
host_function_rules <- tribble(
  ~Gene_regex, ~Functional_category, ~biological_priority,
  
  "^(LGR5|AXIN2|NOTUM|NKD1|MYC|MKI67|TOP2A|CDK1|CCNB1|CCNB2|BUB1|BUB1B|AURKA|AURKB|CDC20|PLK1)$",
  "Tumor remodeling/CRC", 3,
  
  "^(SNAI1|SNAI2|ZEB1|ZEB2|TWIST1|VIM|FN1|MMP1|MMP2|MMP3|MMP7|MMP9|MMP14|COL1A1|COL1A2|COL3A1|COL22A1|COMP|PDLIM4|CRK|MSX1|NSD3)$",
  "Tumor remodeling/CRC", 3,
  
  "^(BAX|BAK1|BCL2|BCL2L1|BOK|CASP3|CASP7|CASP8|CASP9|FAS|FASLG|TNFAIP1|GADD45A|GADD45B|GADD45G|IGFBP[1-7])$",
  "Treatment response", 3,
  
  "^(SOD2|HMOX1|NQO1|TXNRD1|GPX2|NOX1|DUOX2|ATM|ATR|BRCA1|BRCA2|RAD51|PARP1)$",
  "Treatment response", 3,
  
  "^(AHR|ARNT|CYP1A1|CYP1B1|NR1H4|GPBAR1|NR0B2|FGF19|SLC10A2|CYP7A1|S1PR2|S1PR3|HCAR2|HCAR3|NIACR1|FFAR2|FFAR3|GPR41|GPR43|GPR137B)$",
  "Metabolite signaling", 3,
  
  "^(SLC2A1|HK2|LDHA|PDK1|CA9|VEGFA|ENO1|ALDOA|PKM|SLC16A1|SLC16A3|SLC16A4)$",
  "Metabolite signaling", 2,
  
  "^(TJP1|TJP2|TJP3|OCLN|EPCAM|SLC9A3|MUC1|MUC2|MUC5AC|MUC5B|REG1A|REG1B|REG3A|REG3G|DEFA[1-9]|DEFB[0-9]+|CLDN[1-9][0-9]*)$",
  "Barrier/mucus/AMP", 3,
  
  "^(S100A8|S100A9|MARCO|TLR[1-9][0-9]*|NOD1|NOD2|MYD88|IRAK[1-4]|NLRP3|IL1B|TNF|CXCL1|CXCL2|CXCL3|CXCL5|CXCL8|CCL2|CCL3|CCL4|CCL5|C3|C5AR1|FCGR[1-3][A-Z]*)$",
  "Innate immunity", 3,
  
  "^(BATF3|SIPA1|RUNX3|CD3D|CD3E|CD3G|CD4|CD8A|CD8B|TRAC|TRBC1|TRBC2|KLRK1|NKG7|GNLY|GZMB|PRF1|TBX21|GATA3|FOXP3|IL7R)$",
  "Adaptive immunity", 3,
  
  "^(CCR1|CCR2|CCR5|CCR6|CCR7|CCR9|CXCR3|CXCR4|CX3CR1|ITGA4|ITGB7|ICAM1|VCAM1|SELE|SELL|SELPLG|CCL25|CX3CL1)$",
  "Gut homing", 3
)


#=================================================================#
# 1. Reusable functions
#=================================================================#

extract_analysis_set <- function(x, analysis_set) {
  if (
    is.list(x) &&
    !is.data.frame(x) &&
    !is.matrix(x) &&
    analysis_set %in% names(x)
  ) {
    return(x[[analysis_set]])
  }
  
  x
}

normalize_gene_id <- function(x) {
  str_to_upper(
    str_trim(
      str_remove(
        as.character(x),
        "\\.[0-9]+$"
      )
    )
  )
}

normalize_metabolite_id <- function(x) {
  x <- str_to_lower(
    str_squish(
      str_replace_all(
        as.character(x),
        "_",
        " "
      )
    )
  )
  
  dplyr::recode(
    x,
    "lithocholi" = "lithocholic acid",
    "lithocholic" = "lithocholic acid",
    "lithocholicacid" = "lithocholic acid"
  )
}

pick_column <- function(candidates, available_columns) {
  matched_position <- match(
    str_to_lower(candidates),
    str_to_lower(available_columns)
  )
  
  matched_position <- matched_position[
    !is.na(matched_position)
  ]
  
  if (length(matched_position) == 0) {
    return(character())
  }
  
  available_columns[matched_position[1]]
}

rank_columns <- function(x) {
  x <- as.matrix(x)
  
  ranked <- apply(
    x,
    2,
    function(z) {
      rank(
        z,
        ties.method = "average",
        na.last = "keep"
      )
    }
  )
  
  if (is.null(dim(ranked))) {
    ranked <- matrix(
      ranked,
      ncol = 1,
      dimnames = list(
        rownames(x),
        colnames(x)
      )
    )
  } else {
    rownames(ranked) <- rownames(x)
    colnames(ranked) <- colnames(x)
  }
  
  ranked
}

calculate_cross_correlations <- function(x, y, min_complete_samples) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  rho <- suppressWarnings(
    cor(
      x,
      y,
      use = "pairwise.complete.obs",
      method = "pearson"
    )
  )
  
  n_complete <- crossprod(
    1L * is.finite(x),
    1L * is.finite(y)
  )
  
  test_statistic <- rho * sqrt(
    pmax(n_complete - 2, 0) /
      pmax(1 - rho^2, .Machine$double.eps)
  )
  
  p_value <- 2 * pt(
    -abs(test_statistic),
    df = pmax(n_complete - 2, 1)
  )
  
  invalid <- (
    n_complete < min_complete_samples |
      !is.finite(rho)
  )
  
  rho[invalid] <- NA_real_
  p_value[invalid] <- NA_real_
  
  list(
    rho = rho,
    p = p_value,
    n = n_complete
  )
}

residualize_rank_matrix <- function(x, design, min_complete_samples) {
  x <- as.matrix(x)
  
  residual_matrix <- matrix(
    NA_real_,
    nrow = nrow(x),
    ncol = ncol(x),
    dimnames = dimnames(x)
  )
  
  for (j in seq_len(ncol(x))) {
    keep <- (
      is.finite(x[, j]) &
        complete.cases(design)
    )
    
    if (
      sum(keep) < min_complete_samples ||
      length(unique(x[keep, j])) <= 1
    ) {
      next
    }
    
    design_rank <- qr(
      design[keep, , drop = FALSE]
    )$rank
    
    if (sum(keep) < design_rank + 3L) {
      next
    }
    
    residual_matrix[keep, j] <- lm.fit(
      x = design[keep, , drop = FALSE],
      y = rank(
        x[keep, j],
        ties.method = "average"
      )
    )$residuals
  }
  
  residual_matrix
}


#=================================================================#
# 2. Load matched metabolite, metadata, host-expression, and DEG objects
#=================================================================#

loaded_object_names <- load(
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5D_microbe_metabolite_host_axis_screening_results.RData"
)

if (!is.null(host_deg_file)) {
  if (!file.exists(host_deg_file)) {
    stop(
      "The specified host_deg_file does not exist: ",
      host_deg_file,
      call. = FALSE
    )
  }
  
  loaded_object_names <- unique(
    c(
      loaded_object_names,
      load(host_deg_file)
    )
  )
}

if (
  !exists("sample_metadata_results") ||
  !exists("met_use_raw_results")
) {
  stop(
    "The 9.8.1 RData must contain sample_metadata_results and ",
    "met_use_raw_results.",
    call. = FALSE
  )
}

if (
  is.null(sample_metadata_results[[analysis_set]]) ||
  is.null(met_use_raw_results[[analysis_set]])
) {
  stop(
    "The selected analysis_set is absent from sample_metadata_results ",
    "or met_use_raw_results.",
    call. = FALSE
  )
}

metadata <- as.data.frame(
  sample_metadata_results[[analysis_set]],
  check.names = FALSE
)

metabolite_raw <- as.data.frame(
  met_use_raw_results[[analysis_set]],
  check.names = FALSE
)

if (
  is.null(rownames(metadata)) ||
  is.null(rownames(metabolite_raw))
) {
  stop(
    "Metadata and metabolite matrices must use sample IDs as row names.",
    call. = FALSE
  )
}

if (!"TRG_plot" %in% colnames(metadata)) {
  stop(
    "TRG_plot is required for differential metabolite statistics and ",
    "the default adjusted analysis.",
    call. = FALSE
  )
}

if (
  !"pCR" %in% as.character(metadata$TRG_plot) ||
  !"non_pCR" %in% as.character(metadata$TRG_plot)
) {
  stop(
    "TRG_plot must contain both 'pCR' and 'non_pCR'. Current values: ",
    paste(
      sort(
        unique(
          as.character(metadata$TRG_plot)
        )
      ),
      collapse = ", "
    ),
    call. = FALSE
  )
}

#-----------------------------------------------------------------#
# 2.1. Detect the host-expression object
#-----------------------------------------------------------------#

if (!is.null(host_expression_object_name)) {
  if (!exists(host_expression_object_name)) {
    stop(
      "The specified host_expression_object_name was not loaded: ",
      host_expression_object_name,
      call. = FALSE
    )
  }
  
  host_expression_object <- extract_analysis_set(
    get(host_expression_object_name),
    analysis_set
  )
} else {
  expression_candidate_names <- loaded_object_names[
    !loaded_object_names %in% c(
      "sample_metadata_results",
      "met_use_raw_results",
      "met_gene_stats"
    )
  ]
  
  expression_candidates <- bind_rows(
    lapply(
      expression_candidate_names,
      function(object_name) {
        object_value <- extract_analysis_set(
          get(object_name),
          analysis_set
        )
        
        if (
          !is.matrix(object_value) &&
          !is.data.frame(object_value)
        ) {
          return(NULL)
        }
        
        object_value <- try(
          as.data.frame(
            object_value,
            check.names = FALSE
          ),
          silent = TRUE
        )
        
        if (
          inherits(object_value, "try-error") ||
          is.null(rownames(object_value)) ||
          is.null(colnames(object_value))
        ) {
          return(NULL)
        }
        
        row_overlap <- sum(
          rownames(object_value) %in%
            rownames(metadata)
        )
        
        column_overlap <- sum(
          colnames(object_value) %in%
            rownames(metadata)
        )
        
        if (
          max(
            row_overlap,
            column_overlap
          ) < min_complete_samples
        ) {
          return(NULL)
        }
        
        feature_count <- if (
          row_overlap >= column_overlap
        ) {
          ncol(object_value)
        } else {
          nrow(object_value)
        }
        
        tibble(
          object_name = object_name,
          sample_overlap = max(
            row_overlap,
            column_overlap
          ),
          feature_count = feature_count
        )
      }
    )
  )
  
  if (nrow(expression_candidates) == 0) {
    stop(
      "A host-expression matrix could not be detected automatically. ",
      "Set host_expression_object_name near the top of the script. ",
      "Loaded objects: ",
      paste(
        loaded_object_names,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  expression_candidates <- expression_candidates %>%
    arrange(
      desc(feature_count),
      desc(sample_overlap)
    )
  
  host_expression_object_name <-
    expression_candidates$object_name[1]
  
  host_expression_object <- extract_analysis_set(
    get(host_expression_object_name),
    analysis_set
  )
}

host_expression <- as.data.frame(
  host_expression_object,
  check.names = FALSE
)

row_sample_overlap <- sum(
  rownames(host_expression) %in%
    rownames(metadata)
)

column_sample_overlap <- sum(
  colnames(host_expression) %in%
    rownames(metadata)
)

if (
  row_sample_overlap < min_complete_samples &&
  column_sample_overlap < min_complete_samples
) {
  stop(
    "The selected host-expression object does not contain sufficient ",
    "sample IDs matching the metadata: ",
    host_expression_object_name,
    call. = FALSE
  )
}

if (column_sample_overlap > row_sample_overlap) {
  host_expression <- as.data.frame(
    t(as.matrix(host_expression)),
    check.names = FALSE
  )
}

message(
  "Host-expression object: ",
  host_expression_object_name,
  " [",
  nrow(host_expression),
  " samples × ",
  ncol(host_expression),
  " features]"
)

#-----------------------------------------------------------------#
# 2.2. Detect the host-DEG result table
#-----------------------------------------------------------------#

host_deg_source <- NULL

if (!is.null(host_deg_object_name)) {
  if (!exists(host_deg_object_name)) {
    stop(
      "The specified host_deg_object_name was not loaded: ",
      host_deg_object_name,
      call. = FALSE
    )
  }
  
  host_deg_object <- extract_analysis_set(
    get(host_deg_object_name),
    analysis_set
  )
  
  host_deg_source <- host_deg_object_name
} else {
  deg_candidate_names <- loaded_object_names[
    !loaded_object_names %in% c(
      "sample_metadata_results",
      "met_use_raw_results",
      "met_gene_stats",
      host_expression_object_name
    )
  ]
  
  deg_candidates <- bind_rows(
    lapply(
      deg_candidate_names,
      function(object_name) {
        object_value <- extract_analysis_set(
          get(object_name),
          analysis_set
        )
        
        object_value <- try(
          as.data.frame(
            object_value,
            check.names = FALSE
          ),
          silent = TRUE
        )
        
        if (
          inherits(object_value, "try-error") ||
          nrow(object_value) == 0 ||
          ncol(object_value) == 0
        ) {
          return(NULL)
        }
        
        columns_lower <- str_to_lower(
          colnames(object_value)
        )
        
        has_gene_column <- any(
          columns_lower %in% str_to_lower(
            c(
              "Gene", "gene", "gene_symbol", "Gene_symbol",
              "symbol", "SYMBOL", "external_gene_name",
              "gene_name", "Gene_name"
            )
          )
        ) || (
          !is.null(rownames(object_value)) &&
            !all(
              rownames(object_value) ==
                as.character(seq_len(nrow(object_value)))
            )
        )
        
        has_log2fc <- any(
          columns_lower %in% str_to_lower(
            c(
              "log2FoldChange", "log2FC", "logFC",
              "gene_log2FC", "host_logFC", "host_log2FC"
            )
          )
        )
        
        has_fdr <- any(
          columns_lower %in% str_to_lower(
            c(
              "padj", "FDR", "adj.P.Val",
              "qvalue", "gene_FDR", "host_DESeq2_FDR"
            )
          )
        )
        
        if (
          !has_gene_column ||
          !has_log2fc ||
          !has_fdr
        ) {
          return(NULL)
        }
        
        tibble(
          object_name = object_name,
          n_rows = nrow(object_value)
        )
      }
    )
  )
  
  if (nrow(deg_candidates) > 0) {
    deg_candidates <- deg_candidates %>%
      arrange(
        desc(n_rows)
      )
    
    host_deg_object_name <-
      deg_candidates$object_name[1]
    
    host_deg_object <- extract_analysis_set(
      get(host_deg_object_name),
      analysis_set
    )
    
    host_deg_source <- host_deg_object_name
  } else if (
    allow_met_gene_stats_gene_fallback
  ) {
    fallback_correlation_object_name <- c(
      "met_gene_all_results",
      "met_gene_results",
      "met_gene_stats"
    )[
      vapply(
        c(
          "met_gene_all_results",
          "met_gene_results",
          "met_gene_stats"
        ),
        function(object_name) {
          exists(object_name) &&
            is.data.frame(get(object_name)) &&
            "Gene" %in% colnames(get(object_name))
        },
        logical(1)
      )
    ][1]
    
    if (
      length(fallback_correlation_object_name) == 1 &&
      !is.na(fallback_correlation_object_name)
    ) {
      host_deg_object <- as.data.frame(
        get(fallback_correlation_object_name),
        check.names = FALSE
      ) %>%
        distinct(
          Gene,
          .keep_all = TRUE
        )
      
      host_deg_object_name <- paste0(
        fallback_correlation_object_name,
        "_gene_universe"
      )
      
      host_deg_source <- paste0(
        "Unique genes from ",
        fallback_correlation_object_name,
        "; assumes the upstream table was already restricted to host DEGs"
      )
      
      warning(
        "No separate DEG result table was detected. Unique genes from ",
        fallback_correlation_object_name,
        " will be used as the host-DEG universe. Gene-level log2FC/FDR ",
        "filtering will be applied only when those columns are available."
      )
    } else {
      stop(
        "A host-DEG result table could not be detected, and no usable ",
        "metabolite-gene correlation table was available for fallback.",
        call. = FALSE
      )
    }
  } else {
    loaded_object_columns <- vapply(
      loaded_object_names,
      function(object_name) {
        object_value <- try(
          extract_analysis_set(
            get(object_name),
            analysis_set
          ),
          silent = TRUE
        )
        
        if (
          inherits(object_value, "try-error") ||
          (
            !is.data.frame(object_value) &&
            !is.matrix(object_value)
          )
        ) {
          return(
            paste0(
              object_name,
              ": <not a table>"
            )
          )
        }
        
        paste0(
          object_name,
          ": ",
          paste(
            head(
              colnames(object_value),
              12
            ),
            collapse = ", "
          )
        )
      },
      character(1)
    )
    
    stop(
      "A host-DEG result table could not be detected. Set host_deg_file ",
      "and host_deg_object_name, or enable the met_gene_stats fallback. ",
      "Loaded objects and their first columns: ",
      paste(
        loaded_object_columns,
        collapse = " | "
      ),
      call. = FALSE
    )
  }
}

host_deg <- as.data.frame(
  host_deg_object,
  check.names = FALSE
)

message(
  "Host-DEG source: ",
  host_deg_source,
  " [",
  nrow(host_deg),
  " rows before standardization]"
)


#=================================================================#
# 3. Standardize DEG identifiers and select all host DEGs
#=================================================================#

gene_column <- pick_column(
  c(
    "Gene", "gene", "gene_symbol", "Gene_symbol",
    "symbol", "SYMBOL", "external_gene_name",
    "gene_name", "Gene_name"
  ),
  colnames(host_deg)
)

if (length(gene_column) == 0) {
  if (
    !is.null(rownames(host_deg)) &&
    !all(
      rownames(host_deg) ==
      as.character(seq_len(nrow(host_deg)))
    )
  ) {
    host_deg$Gene_from_rownames <- rownames(host_deg)
    gene_column <- "Gene_from_rownames"
  } else {
    stop(
      "The DEG source has no recognizable gene column and no usable ",
      "gene row names.",
      call. = FALSE
    )
  }
}

log2fc_column <- pick_column(
  c(
    "log2FoldChange", "log2FC", "logFC",
    "gene_log2FC", "host_logFC", "host_log2FC"
  ),
  colnames(host_deg)
)

fdr_column <- pick_column(
  c(
    "padj", "FDR", "adj.P.Val",
    "qvalue", "gene_FDR", "host_DESeq2_FDR"
  ),
  colnames(host_deg)
)

pvalue_column <- pick_column(
  c(
    "host_DESeq2_p", "pvalue", "P.Value",
    "p_value", "p", "gene_p"
  ),
  colnames(host_deg)
)

wilcox_p_column <- pick_column(
  c(
    "host_wilcox_p", "wilcox_p", "gene_wilcox_p"
  ),
  colnames(host_deg)
)

gene_statistics_available <- (
  length(log2fc_column) > 0 &&
    (
      length(pvalue_column) > 0 ||
        length(wilcox_p_column) > 0
    )
)

host_deg <- host_deg %>%
  transmute(
    Gene_original = as.character(
      .data[[gene_column[1]]]
    ),
    Gene = normalize_gene_id(
      .data[[gene_column[1]]]
    ),
    gene_log2FC = if (
      length(log2fc_column) > 0
    ) {
      suppressWarnings(
        as.numeric(
          .data[[log2fc_column[1]]]
        )
      )
    } else {
      NA_real_
    },
    gene_p = if (
      length(pvalue_column) > 0
    ) {
      suppressWarnings(
        as.numeric(
          .data[[pvalue_column[1]]]
        )
      )
    } else {
      NA_real_
    },
    gene_wilcox_p = if (
      length(wilcox_p_column) > 0
    ) {
      suppressWarnings(
        as.numeric(
          .data[[wilcox_p_column[1]]]
        )
      )
    } else {
      NA_real_
    },
    gene_FDR = if (
      length(fdr_column) > 0
    ) {
      suppressWarnings(
        as.numeric(
          .data[[fdr_column[1]]]
        )
      )
    } else {
      NA_real_
    }
  ) %>%
  mutate(
    gene_min_p = pmin(
      gene_p,
      gene_wilcox_p,
      na.rm = TRUE
    ),
    gene_min_p = if_else(
      is.infinite(gene_min_p),
      NA_real_,
      gene_min_p
    )
  ) %>%
  filter(
    !is.na(Gene),
    Gene != ""
  ) %>%
  arrange(
    gene_min_p,
    desc(
      abs(gene_log2FC)
    )
  ) %>%
  distinct(
    Gene,
    .keep_all = TRUE
  )

if (gene_statistics_available) {
  host_deg <- host_deg %>%
    filter(
      is.finite(gene_log2FC),
      (
        is.finite(gene_p) &
          gene_p < gene_p_cutoff
      ) |
        (
          is.finite(gene_wilcox_p) &
            gene_wilcox_p < gene_p_cutoff
        )
    )
  
  if (nrow(host_deg) == 0) {
    stop(
      "No host gene met DESeq2 P < ",
      gene_p_cutoff,
      " or Wilcoxon P < ",
      gene_p_cutoff,
      ".",
      call. = FALSE
    )
  }
} else {
  warning(
    "Gene-level log2FC and nominal P-value columns were not available. ",
    "All unique genes in the selected source will be treated as the ",
    "upstream host-gene universe."
  )
}

# Normalize expression-matrix gene identifiers.
colnames(host_expression) <- normalize_gene_id(
  colnames(host_expression)
)

# Collapse duplicated expression columns after identifier normalization.
if (anyDuplicated(colnames(host_expression)) > 0) {
  expression_sum <- rowsum(
    t(as.matrix(host_expression)),
    group = colnames(host_expression),
    reorder = FALSE,
    na.rm = TRUE
  )
  
  expression_n <- rowsum(
    t(1L * is.finite(as.matrix(host_expression))),
    group = colnames(host_expression),
    reorder = FALSE
  )
  
  expression_sum[expression_n == 0] <- NA_real_
  
  host_expression <- as.data.frame(
    t(expression_sum / expression_n),
    check.names = FALSE
  )
}

matched_deg_genes <- intersect(
  host_deg$Gene,
  colnames(host_expression)
)

if (length(matched_deg_genes) == 0) {
  stop(
    "None of the selected DEG identifiers matched the host-expression ",
    "matrix columns. Check whether one object uses Ensembl IDs and the ",
    "other uses gene symbols.",
    call. = FALSE
  )
}

host_deg <- host_deg %>%
  filter(
    Gene %in% matched_deg_genes
  )

host_expression <- host_expression[
  ,
  host_deg$Gene,
  drop = FALSE
]


#=================================================================#
# 4. Match samples and perform feature-level quality control
#=================================================================#

common_samples <- Reduce(
  intersect,
  list(
    rownames(metadata),
    rownames(metabolite_raw),
    rownames(host_expression)
  )
)

if (length(common_samples) < min_complete_samples) {
  stop(
    "Only ",
    length(common_samples),
    " samples are shared by metadata, metabolite, and host-expression ",
    "objects.",
    call. = FALSE
  )
}

metadata <- metadata[
  common_samples,
  ,
  drop = FALSE
]

metabolite_raw <- metabolite_raw[
  common_samples,
  ,
  drop = FALSE
]

host_expression <- host_expression[
  common_samples,
  ,
  drop = FALSE
]

if (
  any(
    !vapply(
      metabolite_raw,
      is.numeric,
      logical(1)
    )
  )
) {
  stop(
    "The metabolite matrix contains non-numeric columns: ",
    paste(
      colnames(metabolite_raw)[
        !vapply(
          metabolite_raw,
          is.numeric,
          logical(1)
        )
      ],
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (
  any(
    !vapply(
      host_expression,
      is.numeric,
      logical(1)
    )
  )
) {
  stop(
    "The host-expression matrix contains non-numeric columns.",
    call. = FALSE
  )
}

colnames(metabolite_raw) <- normalize_metabolite_id(
  colnames(metabolite_raw)
)

if (anyDuplicated(colnames(metabolite_raw)) > 0) {
  stop(
    "Metabolite-name normalization produced duplicated column names: ",
    paste(
      unique(
        colnames(metabolite_raw)[
          duplicated(colnames(metabolite_raw)) |
            duplicated(
              colnames(metabolite_raw),
              fromLast = TRUE
            )
        ]
      ),
      collapse = ", "
    ),
    call. = FALSE
  )
}

# Respect the upstream metabolite-QC decision when met_qc_results is available.
if (exists("met_qc_results")) {
  met_qc_results <- as.data.frame(
    extract_analysis_set(
      met_qc_results,
      analysis_set
    ),
    check.names = FALSE
  )
  
  if (
    all(
      c(
        "Metabolite",
        "keep_metabolite"
      ) %in% colnames(met_qc_results)
    )
  ) {
    metabolite_raw <- metabolite_raw[
      ,
      colnames(metabolite_raw) %in%
        normalize_metabolite_id(
          met_qc_results$Metabolite[
            met_qc_results$keep_metabolite %in% TRUE
          ]
        ),
      drop = FALSE
    ]
  }
}

# Treat a positive assay-specific minimum as the analytical zero.
# Subtracting a constant does not change Spearman ranks, but it prevents
# an assay floor from distorting median ratios and metabolite log2FC.
metabolite_floor_stats <- tibble(
  Metabolite_key = colnames(metabolite_raw),
  assay_floor = vapply(
    metabolite_raw,
    function(x) {
      if (any(is.finite(x))) {
        min(
          x[is.finite(x)]
        )
      } else {
        NA_real_
      }
    },
    numeric(1)
  ),
  n_at_assay_floor = vapply(
    metabolite_raw,
    function(x) {
      if (any(is.finite(x))) {
        sum(
          x[is.finite(x)] ==
            min(
              x[is.finite(x)]
            )
        )
      } else {
        0L
      }
    },
    integer(1)
  )
)

for (i in seq_len(ncol(metabolite_raw))) {
  if (
    is.finite(
      metabolite_floor_stats$assay_floor[i]
    ) &&
    metabolite_floor_stats$assay_floor[i] > 0
  ) {
    finite_values <- is.finite(
      metabolite_raw[
        ,
        i
      ]
    )
    
    metabolite_raw[
      finite_values,
      i
    ] <- pmax(
      metabolite_raw[
        finite_values,
        i
      ] -
        metabolite_floor_stats$assay_floor[i],
      0
    )
  }
}

metabolite_keep <- vapply(
  metabolite_raw,
  function(x) {
    mean(!is.finite(x)) <= maximum_missing_fraction &&
      sum(is.finite(x)) >= min_complete_samples &&
      length(unique(x[is.finite(x)])) > 1
  },
  logical(1)
)

gene_keep <- vapply(
  host_expression,
  function(x) {
    mean(!is.finite(x)) <= maximum_missing_fraction &&
      sum(is.finite(x)) >= min_complete_samples &&
      length(unique(x[is.finite(x)])) > 1
  },
  logical(1)
)

metabolite_raw <- metabolite_raw[
  ,
  metabolite_keep,
  drop = FALSE
]

metabolite_floor_stats <- metabolite_floor_stats %>%
  mutate(
    floor_corrected =
      is.finite(assay_floor) &
      assay_floor > 0,
    retained_after_variance_filter =
      Metabolite_key %in%
      colnames(metabolite_raw)
  )

host_expression <- host_expression[
  ,
  gene_keep,
  drop = FALSE
]

host_deg <- host_deg %>%
  filter(
    Gene %in% colnames(host_expression)
  )

if (
  ncol(metabolite_raw) == 0 ||
  ncol(host_expression) == 0
) {
  stop(
    "No metabolite or DEG feature remained after missingness and ",
    "variance filtering.",
    call. = FALSE
  )
}


#=================================================================#
# 5. Metabolite differential statistics
#=================================================================#

metabolite_node_stats <- tibble(
  Metabolite_key = colnames(metabolite_raw),
  Metabolite_label = str_to_sentence(
    colnames(metabolite_raw)
  ),
  metabolite_log2FC = NA_real_,
  metabolite_p = NA_real_
)

for (i in seq_len(nrow(metabolite_node_stats))) {
  pcr_values <- metabolite_raw[
    as.character(metadata$TRG_plot) == "pCR",
    metabolite_node_stats$Metabolite_key[i],
    drop = TRUE
  ]
  
  non_pcr_values <- metabolite_raw[
    as.character(metadata$TRG_plot) == "non_pCR",
    metabolite_node_stats$Metabolite_key[i],
    drop = TRUE
  ]
  
  pcr_values <- pcr_values[
    is.finite(pcr_values)
  ]
  
  non_pcr_values <- non_pcr_values[
    is.finite(non_pcr_values)
  ]
  
  if (
    length(pcr_values) > 0 &&
    length(non_pcr_values) > 0
  ) {
    pcr_median <- median(
      pcr_values
    )
    
    non_pcr_median <- median(
      non_pcr_values
    )
    
    if (
      pcr_median >= 0 &&
      non_pcr_median >= 0
    ) {
      metabolite_node_stats$metabolite_log2FC[i] <- log2(
        (pcr_median + 1e-06) /
          (non_pcr_median + 1e-06)
      )
    }
  }
  
  if (
    length(pcr_values) >= 2 &&
    length(non_pcr_values) >= 2
  ) {
    metabolite_node_stats$metabolite_p[i] <- suppressWarnings(
      wilcox.test(
        pcr_values,
        non_pcr_values,
        exact = FALSE
      )$p.value
    )
  }
}

metabolite_node_stats <- metabolite_node_stats %>%
  mutate(
    metabolite_FDR = p.adjust(
      metabolite_p,
      method = "BH"
    ),
    Metabolite_class = case_when(
      Metabolite_key %in% c(
        "acetate", "propionate", "butyrate",
        "isobutyrate", "isovalerate", "valerate",
        "gamma aminobutyric acid"
      ) ~ "SCFA / related",
      
      str_detect(
        Metabolite_key,
        "cholic|deoxycholic|ursodeoxycholic|lithocholic|tauro|glyco"
      ) ~ "Bile acid",
      
      Metabolite_key %in% c(
        "indolepropionic acid", "indole lactic acid",
        "indole acetic acid", "indole", "tryptamine",
        "tryptophan", "kynurenic acid", "xanthurenic acid",
        "nicotinic acid"
      ) ~ "Tryptophan / indole",
      
      TRUE ~ "Other metabolite"
    )
  )


#=================================================================#
# 6. Calculate all unadjusted metabolite–DEG Spearman correlations
#=================================================================#

metabolite_rank <- rank_columns(
  metabolite_raw
)

gene_rank <- rank_columns(
  host_expression
)

spearman_result <- calculate_cross_correlations(
  metabolite_rank,
  gene_rank,
  min_complete_samples
)

valid_pairs <- which(
  is.finite(spearman_result$rho),
  arr.ind = TRUE
)

metabolite_gene_correlations <- tibble(
  Metabolite_key = rownames(
    spearman_result$rho
  )[
    valid_pairs[, 1]
  ],
  Gene = colnames(
    spearman_result$rho
  )[
    valid_pairs[, 2]
  ],
  n_complete = spearman_result$n[
    valid_pairs
  ],
  spearman_rho = spearman_result$rho[
    valid_pairs
  ],
  spearman_p = spearman_result$p[
    valid_pairs
  ]
) %>%
  mutate(
    spearman_FDR = p.adjust(
      spearman_p,
      method = "BH"
    )
  )


#=================================================================#
# 7. Covariate-adjusted partial Spearman correlations
#=================================================================#

adjustment_covariates_used <- adjustment_covariates[
  adjustment_covariates %in%
    colnames(metadata)
]

adjustment_covariates_used <- adjustment_covariates_used[
  vapply(
    adjustment_covariates_used,
    function(variable_name) {
      variable_values <- metadata[
        ,
        variable_name
      ]
      
      sum(
        !is.na(variable_values)
      ) >= min_complete_samples &&
        dplyr::n_distinct(
          variable_values[
            !is.na(variable_values)
          ]
        ) > 1
    },
    logical(1)
  )
]

metabolite_gene_correlations$adjusted_n <- NA_integer_
metabolite_gene_correlations$adjusted_rho <- NA_real_
metabolite_gene_correlations$adjusted_p <- NA_real_
metabolite_gene_correlations$adjusted_FDR <- NA_real_

if (length(adjustment_covariates_used) > 0) {
  adjustment_complete <- complete.cases(
    metadata[
      ,
      adjustment_covariates_used,
      drop = FALSE
    ]
  )
  
  if (sum(adjustment_complete) < min_complete_samples) {
    stop(
      "Too few samples have complete adjustment covariates: ",
      paste(
        adjustment_covariates_used,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  adjustment_metadata <- droplevels(
    metadata[
      adjustment_complete,
      adjustment_covariates_used,
      drop = FALSE
    ]
  )
  
  adjustment_design <- model.matrix(
    ~ .,
    data = adjustment_metadata
  )
  
  metabolite_residuals <- residualize_rank_matrix(
    metabolite_raw[
      adjustment_complete,
      ,
      drop = FALSE
    ],
    adjustment_design,
    min_complete_samples
  )
  
  gene_residuals <- residualize_rank_matrix(
    host_expression[
      adjustment_complete,
      ,
      drop = FALSE
    ],
    adjustment_design,
    min_complete_samples
  )
  
  adjusted_result <- calculate_cross_correlations(
    metabolite_residuals,
    gene_residuals,
    min_complete_samples
  )
  
  # Partial-correlation inference must account for the covariate-model
  # degrees of freedom. The ordinary n - 2 correlation test is too liberal.
  adjusted_df <- adjusted_result$n -
    qr(adjustment_design)$rank -
    1
  
  adjusted_result$p <- 2 * pt(
    -abs(
      adjusted_result$rho * sqrt(
        pmax(
          adjusted_df,
          0
        ) /
          pmax(
            1 - adjusted_result$rho^2,
            .Machine$double.eps
          )
      )
    ),
    df = pmax(
      adjusted_df,
      1
    )
  )
  
  adjusted_result$p[
    adjusted_df <= 0 |
      !is.finite(
        adjusted_result$rho
      )
  ] <- NA_real_
  
  adjusted_row <- match(
    metabolite_gene_correlations$Metabolite_key,
    rownames(adjusted_result$rho)
  )
  
  adjusted_column <- match(
    metabolite_gene_correlations$Gene,
    colnames(adjusted_result$rho)
  )
  
  adjusted_index <- cbind(
    adjusted_row,
    adjusted_column
  )
  
  metabolite_gene_correlations$adjusted_n <-
    adjusted_result$n[
      adjusted_index
    ]
  
  metabolite_gene_correlations$adjusted_rho <-
    adjusted_result$rho[
      adjusted_index
    ]
  
  metabolite_gene_correlations$adjusted_p <-
    adjusted_result$p[
      adjusted_index
    ]
  
  metabolite_gene_correlations$adjusted_FDR <- p.adjust(
    metabolite_gene_correlations$adjusted_p,
    method = "BH"
  )
} else {
  warning(
    "None of the requested adjustment covariates was available and ",
    "informative. Adjusted correlations will remain NA."
  )
}


#=================================================================#
# 8. Add DEG, metabolite, and Figure 5C annotations
#=================================================================#

panel_c_genes <- c(
  "MSX1", "NOTUM", "UCA1", "NKD1", "HTRA1", "PLOD1",
  "LIG4", "CLPTM1", "METTL14", "ZBTB7B",
  "PIK3CB", "CRK", "KLRK1", "RAB27A",
  "NINJ1", "CX3CL1", "ADAM8", "ICAM1"
)

metabolite_gene_correlations <- metabolite_gene_correlations %>%
  left_join(
    host_deg,
    by = "Gene"
  ) %>%
  left_join(
    metabolite_node_stats,
    by = "Metabolite_key"
  ) %>%
  mutate(
    Panel_C_gene = Gene %in% panel_c_genes,
    
    adjusted_direction_consistent =
      is.finite(adjusted_rho) &
      sign(adjusted_rho) == sign(spearman_rho),
    
    spearman_high_confidence =
      spearman_p < spearman_strong_p_cutoff |
      spearman_FDR < spearman_relaxed_fdr_cutoff,
    
    spearman_evidence_tier = case_when(
      spearman_FDR < spearman_relaxed_fdr_cutoff ~ "FDR < 0.10",
      spearman_p < spearman_strong_p_cutoff ~ "Nominal P < 0.01",
      spearman_p < spearman_moderate_p_cutoff ~ "Nominal P < 0.05",
      spearman_p < spearman_p_cutoff ~ "Nominal P < 0.10",
      TRUE ~ "Not selected"
    ),
    
    adjusted_support =
      adjusted_direction_consistent &
      abs(adjusted_rho) >= adjusted_rho_cutoff &
      adjusted_p < adjusted_p_cutoff,
    
    exclude_from_main =
      Gene %in% manual_gene_exclude |
      str_detect(
        Gene,
        regex(
          paste0(
            "-AS[0-9]*$|^LINC[0-9]+$|^LOC[0-9]+$|",
            "^C[0-9XY]+ORF[0-9]+$|^AC[0-9]+|^AL[0-9]+|^AP[0-9]+|",
            "^RP[0-9]+-|^MIR[0-9]+|^SNOR[AD][0-9]+|^RNU[0-9]+|",
            "^MT-T|^MT-R|^HCG[0-9]+|",
            "^RPL[0-9A-Z]*P[0-9]+$|^RPS[0-9A-Z]*P[0-9]+$|",
            "^KRT[0-9A-Z]*P[0-9]+$|^EEF1A1P[0-9]+$|",
            "^GAPDHP[0-9]+$|^HNRNP[A-Z0-9]*P[0-9]+$"
          ),
          ignore_case = TRUE
        )
      ),
    
    manual_priority_gene =
      Gene %in% manual_priority_genes,
    
    Functional_category = NA_character_,
    biological_priority = 0
  )

for (i in seq_len(nrow(manual_gene_category))) {
  hit <- metabolite_gene_correlations$Gene ==
    normalize_gene_id(manual_gene_category$Gene[i])
  
  metabolite_gene_correlations$Functional_category[hit] <-
    manual_gene_category$Functional_category[i]
  
  metabolite_gene_correlations$biological_priority[hit] <-
    manual_gene_category$biological_priority[i]
}

for (i in seq_len(nrow(host_function_rules))) {
  hit <- is.na(
    metabolite_gene_correlations$Functional_category
  ) & str_detect(
    metabolite_gene_correlations$Gene,
    regex(
      host_function_rules$Gene_regex[i],
      ignore_case = TRUE
    )
  )
  
  metabolite_gene_correlations$Functional_category[hit] <-
    host_function_rules$Functional_category[i]
  
  metabolite_gene_correlations$biological_priority[hit] <-
    host_function_rules$biological_priority[i]
}

metabolite_gene_correlations <- metabolite_gene_correlations %>%
  mutate(
    biologically_annotated =
      !is.na(Functional_category) &
      biological_priority > 0,
    
    Functional_category = factor(
      Functional_category,
      levels = names(category_palette)
    )
  )


#=================================================================#
# 10. Save the complete association table
#=================================================================#

write.csv(
  metabolite_gene_correlations %>%
    arrange(
      spearman_FDR,
      desc(abs(spearman_rho))
    ),
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_8.4_all_metabolite_DEG_correlations_",
    analysis_set,
    ".csv"
  ),
  row.names = FALSE
)


#=================================================================#
# 11. Identify exploratory, biologically interpretable main-figure edges
#=================================================================#

significant_edges <- metabolite_gene_correlations %>%
  filter(
    is.finite(spearman_rho),
    is.finite(spearman_p),
    abs(spearman_rho) >= spearman_rho_cutoff,
    spearman_p < spearman_p_cutoff
  )

if (nrow(significant_edges) == 0) {
  stop(
    "No metabolite–DEG pair met nominal Spearman P < ",
    spearman_p_cutoff,
    " and absolute rho >= ",
    spearman_rho_cutoff,
    ". The complete correlation table was saved.",
    call. = FALSE
  )
}

biological_edges <- significant_edges %>%
  filter(
    !exclude_from_main,
    biologically_annotated
  )

if (nrow(biological_edges) == 0) {
  stop(
    "Exploratory correlations were detected, but none matched the ",
    "biological annotation rules. Expand host_function_rules rather than ",
    "plotting uncharacterized genes.",
    call. = FALSE
  )
}

write.csv(
  significant_edges %>%
    arrange(
      desc(biologically_annotated),
      desc(manual_priority_gene),
      desc(biological_priority),
      desc(spearman_high_confidence),
      spearman_p,
      spearman_FDR,
      desc(abs(spearman_rho))
    ),
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_8.4_significant_metabolite_DEG_edges_",
    analysis_set,
    ".csv"
  ),
  row.names = FALSE
)

# Diagnose why specifically requested genes are present or absent.
priority_gene_diagnostics <- tibble(
  Gene = manual_priority_genes
) %>%
  distinct() %>%
  left_join(
    host_deg %>%
      select(
        Gene,
        gene_log2FC,
        gene_p,
        gene_wilcox_p,
        gene_min_p,
        gene_FDR
      ),
    by = "Gene",
    relationship = "one-to-one"
  ) %>%
  left_join(
    metabolite_gene_correlations %>%
      group_by(Gene) %>%
      summarise(
        tested_pairs = n(),
        best_spearman_p = min(
          spearman_p,
          na.rm = TRUE
        ),
        best_abs_spearman_rho = max(
          abs(spearman_rho),
          na.rm = TRUE
        ),
        n_p_lt_0_10 = sum(
          spearman_p < 0.10 &
            abs(spearman_rho) >= spearman_rho_cutoff,
          na.rm = TRUE
        ),
        n_p_lt_0_05 = sum(
          spearman_p < 0.05 &
            abs(spearman_rho) >= spearman_rho_cutoff,
          na.rm = TRUE
        ),
        n_adjusted_support = sum(
          adjusted_support,
          na.rm = TRUE
        ),
        .groups = "drop"
      ),
    by = "Gene",
    relationship = "one-to-one"
  ) %>%
  mutate(
    in_DEG_universe = !is.na(gene_log2FC),
    has_eligible_edge = dplyr::coalesce(
      n_p_lt_0_10,
      0L
    ) > 0
  )

#=================================================================#
# 12. Select metabolites using differential signal and class diversity
#
# Metabolites with fewer than three eligible biological gene edges are
# excluded from the ordinary candidate pool. After edge selection, the same
# threshold is applied again to the final displayed gene-pair count.
# Butyrate is retained as the only explicit final-display exception.
#=================================================================#

metabolite_selection <- biological_edges %>%
  group_by(
    Metabolite_key,
    Metabolite_label,
    Metabolite_class,
    metabolite_log2FC,
    metabolite_p,
    metabolite_FDR
  ) %>%
  summarise(
    n_biological_genes = n_distinct(Gene),
    n_manual_priority_genes = n_distinct(
      Gene[manual_priority_gene]
    ),
    summed_biological_priority = sum(
      biological_priority,
      na.rm = TRUE
    ),
    n_high_confidence_edges = sum(
      spearman_high_confidence,
      na.rm = TRUE
    ),
    best_spearman_p = min(
      spearman_p,
      na.rm = TRUE
    ),
    mean_abs_spearman_rho = mean(
      abs(spearman_rho),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    manual_metabolite_priority =
      Metabolite_key %in%
      manual_priority_metabolites,
    
    literature_metabolite_priority =
      Metabolite_key %in%
      literature_priority_metabolites,
    
    differential_tier = case_when(
      is.finite(metabolite_FDR) &
        metabolite_FDR < 0.10 ~ 4L,
      is.finite(metabolite_p) &
        metabolite_p < 0.05 ~ 3L,
      is.finite(metabolite_p) &
        metabolite_p < 0.10 ~ 2L,
      is.finite(metabolite_log2FC) &
        abs(metabolite_log2FC) >= 0.50 ~ 1L,
      TRUE ~ 0L
    )
  )

metabolite_selection <- metabolite_selection %>%
  group_by(
    Metabolite_class
  ) %>%
  arrange(
    desc(differential_tier),
    desc(
      abs(
        metabolite_log2FC
      )
    ),
    desc(n_manual_priority_genes),
    desc(summed_biological_priority),
    desc(n_high_confidence_edges),
    best_spearman_p,
    .by_group = TRUE
  ) %>%
  mutate(
    metabolite_class_rank = row_number()
  ) %>%
  ungroup()

selected_metabolites <- metabolite_selection %>%
  filter(
    !Metabolite_key %in%
      manual_metabolite_exclude,
    (
      n_biological_genes >= 3 &
        metabolite_class_rank <=
        max_metabolites_per_class
    ) |
      manual_metabolite_priority |
      literature_metabolite_priority
  ) %>%
  arrange(
    desc(manual_metabolite_priority),
    desc(differential_tier),
    desc(
      abs(
        metabolite_log2FC
      )
    ),
    desc(literature_metabolite_priority),
    desc(n_manual_priority_genes),
    desc(summed_biological_priority),
    desc(n_high_confidence_edges),
    best_spearman_p,
    desc(mean_abs_spearman_rho)
  ) %>%
  slice_head(
    n = max_selected_metabolites
  ) %>%
  mutate(
    selection_order = row_number()
  )

#=================================================================#
# 13. Select host DEGs primarily by Spearman evidence
#=================================================================#

selected_metabolite_edge_pool <- biological_edges %>%
  filter(
    Metabolite_key %in%
      selected_metabolites$Metabolite_key
  )

gene_selection <- selected_metabolite_edge_pool %>%
  group_by(
    Gene
  ) %>%
  summarise(
    Gene_original = first(
      Gene_original
    ),
    gene_log2FC = first(
      gene_log2FC
    ),
    gene_p = first(
      gene_p
    ),
    gene_wilcox_p = first(
      gene_wilcox_p
    ),
    gene_min_p = first(
      gene_min_p
    ),
    gene_FDR = first(
      gene_FDR
    ),
    Panel_C_gene = first(
      Panel_C_gene
    ),
    Functional_category = first(
      Functional_category
    ),
    biological_priority = max(
      biological_priority,
      na.rm = TRUE
    ),
    manual_priority_gene = any(
      manual_priority_gene
    ),
    n_significant_metabolites = n_distinct(
      Metabolite_key
    ),
    n_high_confidence_edges = sum(
      spearman_high_confidence,
      na.rm = TRUE
    ),
    best_spearman_p = min(
      spearman_p,
      na.rm = TRUE
    ),
    maximum_abs_spearman_rho = max(
      abs(spearman_rho),
      na.rm = TRUE
    ),
    mean_abs_spearman_rho = mean(
      abs(spearman_rho),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    DEG_side = case_when(
      is.finite(gene_log2FC) &
        gene_log2FC > 0 ~ positive_deg_label,
      TRUE ~ negative_deg_label
    )
  ) %>%
  arrange(
    best_spearman_p,
    desc(maximum_abs_spearman_rho),
    desc(n_high_confidence_edges),
    desc(n_significant_metabolites),
    desc(biological_priority),
    gene_FDR,
    desc(manual_priority_gene),
    Gene
  )

# Use a modestly larger pCR quota. Manual priority genes are not forced;
# they are used only after the correlation and biological-ranking terms.
selected_genes <- bind_rows(
  gene_selection %>%
    filter(
      DEG_side == positive_deg_label
    ) %>%
    slice_head(n = 24),
  
  gene_selection %>%
    filter(
      DEG_side == negative_deg_label
    ) %>%
    slice_head(n = 20),
  
  gene_selection
) %>%
  distinct(
    Gene,
    .keep_all = TRUE
  ) %>%
  slice_head(
    n = max_selected_genes
  )

#=================================================================#
# 14. Select strongest displayed edges for each selected gene
#=================================================================#

edge_pool <- selected_metabolite_edge_pool %>%
  filter(
    Gene %in%
      selected_genes$Gene
  ) %>%
  mutate(
    edge_id = paste(
      Metabolite_key,
      Gene,
      sep = "|||"
    )
  )

# Retain the strongest eligible metabolite edge for every selected gene.
# This prevents biologically prioritized genes from disappearing simply
# because another gene has a slightly smaller correlation P value.
displayed_edges <- edge_pool %>%
  group_by(
    Gene
  ) %>%
  arrange(
    spearman_p,
    desc(abs(spearman_rho)),
    desc(spearman_high_confidence),
    desc(biological_priority),
    desc(manual_priority_gene),
    .by_group = TRUE
  ) %>%
  slice_head(n = 1) %>%
  ungroup()

maximum_display_edges <- min(
  nrow(edge_pool),
  nrow(selected_metabolites) *
    max_edges_per_metabolite,
  nrow(selected_genes) *
    max_edges_per_gene
)

while (
  nrow(displayed_edges) <
  maximum_display_edges
) {
  metabolite_degree <- displayed_edges %>%
    count(
      Metabolite_key,
      name = "metabolite_degree"
    )
  
  gene_degree <- displayed_edges %>%
    count(
      Gene,
      name = "gene_degree"
    )
  
  remaining_edges <- edge_pool %>%
    anti_join(
      displayed_edges %>%
        select(edge_id),
      by = "edge_id"
    ) %>%
    left_join(
      metabolite_degree,
      by = "Metabolite_key"
    ) %>%
    left_join(
      gene_degree,
      by = "Gene"
    ) %>%
    mutate(
      metabolite_degree = coalesce(
        metabolite_degree,
        0L
      ),
      gene_degree = coalesce(
        gene_degree,
        0L
      )
    ) %>%
    filter(
      metabolite_degree <
        max_edges_per_metabolite,
      gene_degree <
        max_edges_per_gene
    )
  
  if (nrow(remaining_edges) == 0) {
    break
  }
  
  next_edge <- remaining_edges %>%
    arrange(
      spearman_p,
      desc(abs(spearman_rho)),
      desc(spearman_high_confidence),
      desc(biological_priority),
      gene_FDR,
      desc(manual_priority_gene)
    ) %>%
    slice_head(n = 1) %>%
    select(
      -metabolite_degree,
      -gene_degree
    )
  
  displayed_edges <- bind_rows(
    displayed_edges,
    next_edge
  )
}

displayed_edges <- displayed_edges %>%
  distinct(
    edge_id,
    .keep_all = TRUE
  ) %>%
  filter(
    !Metabolite_key %in%
      manual_metabolite_exclude
  ) %>%
  group_by(
    Metabolite_key
  ) %>%
  mutate(
    n_displayed_gene_pairs = n_distinct(
      Gene
    )
  ) %>%
  ungroup() %>%
  filter(
    n_displayed_gene_pairs >= 3 |
      Metabolite_key %in%
      manual_priority_metabolites
  )

# Do not refill empty metabolite slots after this final filter. This keeps
# the network limited to metabolites with at least three displayed gene
# pairs, with butyrate retained as the only explicit exception.

selected_metabolites <- selected_metabolites %>%
  filter(
    Metabolite_key %in%
      displayed_edges$Metabolite_key
  ) %>%
  distinct(
    Metabolite_key,
    .keep_all = TRUE
  ) %>%
  arrange(
    desc(differential_tier),
    desc(
      abs(
        metabolite_log2FC
      )
    ),
    selection_order
  ) %>%
  mutate(
    metabolite_display_order = as.numeric(
      row_number()
    ),
    metabolite_display_order = case_when(
      Metabolite_key == "glutamic acid" &
        "indole acetic acid" %in%
        Metabolite_key ~
        match(
          "indole acetic acid",
          Metabolite_key
        ) - 0.5,
      TRUE ~ metabolite_display_order
    )
  ) %>%
  arrange(
    metabolite_display_order
  )

selected_genes <- selected_genes %>%
  filter(
    Gene %in%
      displayed_edges$Gene
  ) %>%
  distinct(
    Gene,
    .keep_all = TRUE
  ) %>%
  left_join(
    displayed_edges %>%
      count(
        Gene,
        name = "n_displayed_edges"
      ),
    by = "Gene",
    relationship = "many-to-one"
  ) %>%
  left_join(
    displayed_edges %>%
      arrange(
        Gene,
        spearman_p,
        desc(abs(spearman_rho))
      ) %>%
      group_by(
        Gene
      ) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      transmute(
        Gene,
        anchor_metabolite = Metabolite_key,
        anchor_metabolite_order = match(
          Metabolite_key,
          selected_metabolites$Metabolite_key
        ),
        anchor_spearman_p = spearman_p,
        anchor_abs_spearman_rho = abs(spearman_rho)
      ),
    by = "Gene",
    relationship = "many-to-one"
  ) %>%
  arrange(
    DEG_side,
    anchor_metabolite_order,
    Functional_category,
    anchor_spearman_p,
    desc(anchor_abs_spearman_rho),
    best_spearman_p,
    Gene
  )

# The main network is selected exclusively by ordinary Spearman evidence
# and biological annotation. Adjusted support is added only after selection.
displayed_edges$adjusted_loo_sign_fraction <- NA_real_

if (
  exists("metabolite_residuals") &&
  exists("gene_residuals")
) {
  for (i in seq_len(nrow(displayed_edges))) {
    residual_metabolite <- metabolite_residuals[
      ,
      displayed_edges$Metabolite_key[i]
    ]
    
    residual_gene <- gene_residuals[
      ,
      displayed_edges$Gene[i]
    ]
    
    complete_residuals <- which(
      is.finite(residual_metabolite) &
        is.finite(residual_gene)
    )
    
    if (
      length(complete_residuals) >=
      min_complete_samples + 1L &&
      is.finite(
        displayed_edges$adjusted_rho[i]
      )
    ) {
      leave_one_out_rho <- vapply(
        complete_residuals,
        function(removed_sample) {
          retained_samples <- setdiff(
            complete_residuals,
            removed_sample
          )
          
          suppressWarnings(
            cor(
              residual_metabolite[
                retained_samples
              ],
              residual_gene[
                retained_samples
              ],
              method = "pearson"
            )
          )
        },
        numeric(1)
      )
      
      displayed_edges$adjusted_loo_sign_fraction[i] <- mean(
        sign(leave_one_out_rho) ==
          sign(
            displayed_edges$adjusted_rho[i]
          ),
        na.rm = TRUE
      )
    }
  }
}

displayed_edges <- displayed_edges %>%
  mutate(
    adjusted_plot_support =
      adjusted_support &
      is.finite(
        adjusted_loo_sign_fraction
      ) &
      adjusted_loo_sign_fraction >=
      adjusted_loo_sign_fraction_cutoff
  )

priority_gene_diagnostics <- priority_gene_diagnostics %>%
  mutate(
    selected_for_figure =
      Gene %in% selected_genes$Gene
  )

write.csv(
  priority_gene_diagnostics,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_8.4_priority_gene_diagnostics_",
    analysis_set,
    ".csv"
  ),
  row.names = FALSE
)

#=================================================================#
# 15. Rail coordinates
#=================================================================#

if (anyDuplicated(selected_metabolites$Metabolite_key) > 0) {
  stop(
    "selected_metabolites contains duplicated Metabolite_key values.",
    call. = FALSE
  )
}

if (anyDuplicated(selected_genes$Gene) > 0) {
  stop(
    "selected_genes contains duplicated Gene values.",
    call. = FALSE
  )
}

x_left <- -1.55
x_center <- 0
x_right <- 1.55

selected_genes <- selected_genes %>%
  mutate(
    gene_x = case_when(
      DEG_side == positive_deg_label ~ x_left,
      TRUE ~ x_right
    )
  )

left_genes <- selected_genes %>%
  filter(
    DEG_side == positive_deg_label
  ) %>%
  arrange(
    anchor_metabolite_order,
    Functional_category,
    anchor_spearman_p,
    desc(anchor_abs_spearman_rho),
    best_spearman_p,
    Gene
  )

right_genes <- selected_genes %>%
  filter(
    DEG_side == negative_deg_label
  ) %>%
  arrange(
    anchor_metabolite_order,
    Functional_category,
    anchor_spearman_p,
    desc(anchor_abs_spearman_rho),
    best_spearman_p,
    Gene
  )

y_top <- 1 + (
  max(
    nrow(left_genes),
    nrow(right_genes),
    nrow(selected_metabolites),
    1L
  ) - 1
) * 0.74

if (nrow(left_genes) > 0) {
  left_genes$y <- seq(
    from = y_top,
    to = 1,
    length.out = nrow(left_genes)
  )
}

if (nrow(right_genes) > 0) {
  right_genes$y <- seq(
    from = y_top,
    to = 1,
    length.out = nrow(right_genes)
  )
}

selected_genes <- bind_rows(
  left_genes,
  right_genes
) %>%
  arrange(
    gene_x,
    desc(y)
  )

selected_metabolites <- selected_metabolites %>%
  mutate(
    y = seq(
      from = y_top,
      to = 1,
      length.out = nrow(selected_metabolites)
    ),
    metabolite_x = x_center
  )

displayed_edges <- displayed_edges %>%
  left_join(
    selected_metabolites %>%
      select(
        Metabolite_key,
        metabolite_y = y,
        metabolite_x
      ),
    by = "Metabolite_key",
    relationship = "many-to-one"
  ) %>%
  left_join(
    selected_genes %>%
      select(
        Gene,
        gene_y = y,
        gene_x,
        DEG_side,
        Functional_category
      ),
    by = "Gene",
    relationship = "many-to-one"
  ) %>%
  arrange(
    match(
      Metabolite_key,
      selected_metabolites$Metabolite_key
    ),
    DEG_side,
    desc(abs(spearman_rho))
  )

if (
  any(
    !is.finite(displayed_edges$metabolite_y) |
    !is.finite(displayed_edges$gene_y)
  )
) {
  stop(
    "One or more displayed edges has a missing rail coordinate. ",
    "Check metabolite and gene identifier matching before plotting.",
    call. = FALSE
  )
}

edge_curve_data <- displayed_edges %>%
  mutate(
    edge_plot_order = row_number(),
    control1_x = metabolite_x + 0.35 * (gene_x - metabolite_x),
    control2_x = metabolite_x + 0.75 * (gene_x - metabolite_x),
    control1_y = metabolite_y,
    control2_y = gene_y
  ) %>%
  tidyr::uncount(
    weights = 51,
    .id = "curve_point"
  ) %>%
  mutate(
    curve_t = (
      curve_point - 1
    ) / 50,
    
    curve_x =
      (1 - curve_t)^3 * metabolite_x +
      3 * (1 - curve_t)^2 * curve_t * control1_x +
      3 * (1 - curve_t) * curve_t^2 * control2_x +
      curve_t^3 * gene_x,
    
    curve_y =
      (1 - curve_t)^3 * metabolite_y +
      3 * (1 - curve_t)^2 * curve_t * control1_y +
      3 * (1 - curve_t) * curve_t^2 * control2_y +
      curve_t^3 * gene_y
  ) %>%
  arrange(
    edge_plot_order,
    curve_point
  )

top_label_y <- y_top + 0.75

#=================================================================#
# 16. Draw the structured bipartite rail network
#=================================================================#

if (
  all(
    !is.finite(
      selected_metabolites$metabolite_log2FC
    )
  )
) {
  metabolite_fc_limit <- 1
} else {
  metabolite_fc_limit <- max(
    1,
    as.numeric(
      quantile(
        abs(
          selected_metabolites$metabolite_log2FC[
            is.finite(
              selected_metabolites$metabolite_log2FC
            )
          ]
        ),
        0.90,
        na.rm = TRUE
      )
    )
  )
}

selected_metabolites <- selected_metabolites %>%
  mutate(
    metabolite_log2FC_plot = pmax(
      pmin(
        metabolite_log2FC,
        metabolite_fc_limit
      ),
      -metabolite_fc_limit
    )
  )

p_metabolite_host_rail <- ggplot() +
  annotate(
    "segment",
    x = x_left,
    xend = x_left,
    y = 0.65,
    yend = y_top + 0.35,
    color = "grey80",
    linewidth = 0.25
  ) +
  annotate(
    "segment",
    x = x_center,
    xend = x_center,
    y = 0.65,
    yend = y_top + 0.35,
    color = "grey78",
    linewidth = 0.28
  ) +
  annotate(
    "segment",
    x = x_right,
    xend = x_right,
    y = 0.65,
    yend = y_top + 0.35,
    color = "grey80",
    linewidth = 0.25
  ) +
  geom_path(
    data = edge_curve_data,
    aes(
      x = curve_x,
      y = curve_y,
      group = edge_id,
      color = spearman_rho
    ),
    linewidth = edge_linewidth_pt * 0.3527778,
    lineend = "round",
    linejoin = "round",
    alpha = 0.92,
    na.rm = TRUE
  ) +
  geom_point(
    data = edge_curve_data %>%
      filter(
        curve_point == 26,
        adjusted_plot_support
      ),
    aes(
      x = curve_x,
      y = curve_y,
      shape = "TRG-adjusted residual-rank support"
    ),
    inherit.aes = FALSE,
    size = 1.65,
    fill = "white",
    color = "grey15",
    stroke = 0.40
  ) +
  geom_point(
    data = selected_metabolites,
    aes(
      x = metabolite_x,
      y = y,
      fill = metabolite_log2FC_plot
    ),
    shape = 21,
    size = 4.6,
    color = "grey30",
    stroke = 0.55
  ) +
  geom_text(
    data = selected_metabolites,
    aes(
      x = metabolite_x,
      y = y,
      label = Metabolite_label
    ),
    nudge_y = 0.38,
    hjust = 0.5,
    size = 2.95,
    color = "grey15"
  ) +
  annotate(
    "text",
    x = x_left,
    y = top_label_y,
    label = "pCR-enriched DEGs",
    fontface = "bold",
    size = 3.3
  ) +
  annotate(
    "text",
    x = x_center,
    y = top_label_y,
    label = "Metabolites",
    fontface = "bold",
    size = 3.3
  ) +
  annotate(
    "text",
    x = x_right,
    y = top_label_y,
    label = "non-pCR-enriched DEGs",
    fontface = "bold",
    size = 3.3
  ) +
  scale_color_gradient2(
    low = "#5A83AD",
    mid = "#EDEAE6",
    high = "#B86B68",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman rho",
    guide = guide_colorbar(
      order = 3,
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(30, "mm")
    )
  ) +
  scale_fill_gradient2(
    low = "#D97A6C",
    mid = "white",
    high = "#5AB4AC",
    midpoint = 0,
    limits = c(
      -metabolite_fc_limit,
      metabolite_fc_limit
    ),
    name = "Metabolite log2FC
pCR vs non-pCR",
    guide = guide_colorbar(
      order = 4,
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(4, "mm"),
      barheight = grid::unit(30, "mm")
    )
  ) +
  scale_shape_manual(
    values = c(
      "TRG-adjusted residual-rank support" = 23
    ),
    name = "Sensitivity analysis",
    guide = guide_legend(
      order = 1,
      override.aes = list(
        fill = "white",
        color = "grey15",
        size = 2.2,
        stroke = 0.45
      )
    )
  ) +
  ggnewscale::new_scale_fill() +
  geom_point(
    data = selected_genes,
    aes(
      x = gene_x,
      y = y,
      fill = Functional_category
    ),
    shape = 21,
    size = 2.45,
    color = "grey35",
    stroke = 0.35
  ) +
  geom_text(
    data = selected_genes %>%
      filter(
        gene_x < 0
      ),
    aes(
      x = gene_x - 0.08,
      y = y,
      label = Gene
    ),
    hjust = 1,
    size = 3.0,
    color = "grey20"
  ) +
  geom_text(
    data = selected_genes %>%
      filter(
        gene_x > 0
      ),
    aes(
      x = gene_x + 0.08,
      y = y,
      label = Gene
    ),
    hjust = 0,
    size = 3.0,
    color = "grey20"
  ) +
  scale_fill_manual(
    values = category_palette,
    breaks = names(category_palette),
    drop = FALSE,
    name = "Category",
    guide = guide_legend(
      order = 2,
      ncol = 1,
      byrow = TRUE,
      override.aes = list(
        shape = 21,
        size = 4,
        color = "grey35",
        stroke = 0.35
      )
    )
  ) +
  coord_cartesian(
    xlim = c(x_left - 1.10, x_right + 1.20),
    ylim = c(
      0.45,
      top_label_y + 0.35
    ),
    clip = "off"
  ) +
  theme_void(
    base_size = 10
  ) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.justification = "center",
    legend.title = element_text(
      size = 8.5
    ),
    legend.text = element_text(
      size = 8
    ),
    legend.spacing.y = grid::unit(
      1.3,
      "mm"
    ),
    plot.margin = margin(
      10,
      140,
      10,
      140
    )
  )

p_metabolite_host_rail

ggsave(
  paste0(
    "figures/Fig5D_metabolite_host_DEG_rail_",
    analysis_set,
    ".svg"
  ),
  p_metabolite_host_rail,
  width = 12.0,
  height = max(
    6.1,
    0.22 * y_top + 2.45
  ),
  device = "svg"
)

#=================================================================#
# 17. Save displayed data, diagnostics, and essential R objects
#=================================================================#

write.csv(
  metabolite_floor_stats,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_8.4_metabolite_assay_floor_",
    analysis_set,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  displayed_edges %>%
    select(
      -edge_id
    ),
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_8.4_displayed_metabolite_DEG_edges_",
    analysis_set,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  selected_metabolites,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_8.4_selected_metabolites_",
    analysis_set,
    ".csv"
  ),
  row.names = FALSE
)

write.csv(
  selected_genes,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_8.4_selected_host_DEGs_",
    analysis_set,
    ".csv"
  ),
  row.names = FALSE
)

diagnostics <- tibble(
  Metric = c(
    "Analysis set",
    "Host-expression object",
    "Host-DEG source",
    "Gene-level nominal DEG statistics available",
    "Matched samples",
    "Metabolites tested",
    paste0(
      "Host genes with DESeq2 or Wilcoxon P < ",
      gene_p_cutoff
    ),
    "Host DEGs tested",
    "All tested metabolite-DEG pairs",
    "Nominal Spearman P < 0.10",
    paste0(
      "Nominal Spearman P < ",
      spearman_p_cutoff,
      " and absolute rho >= ",
      spearman_rho_cutoff
    ),
    paste0(
      "Nominal Spearman P < ",
      spearman_strong_p_cutoff,
      " and absolute rho >= ",
      spearman_rho_cutoff
    ),
    paste0(
      "Spearman FDR < ",
      spearman_relaxed_fdr_cutoff,
      " and absolute rho >= ",
      spearman_rho_cutoff
    ),
    "High-confidence eligible edges",
    "Preliminary TRG-adjusted residual-rank edges",
    "Displayed LOO-stable adjusted-support markers",
    "Butyrate selected",
    "Lithocholic acid selected",
    "Adjustment covariates",
    "Biologically annotated eligible edges",
    "Priority genes with eligible edges",
    "Priority genes selected for figure",
    "Selected metabolites",
    "Selected host DEGs",
    "Selected pCR-enriched DEGs",
    "Selected non-pCR-enriched DEGs",
    "Displayed edges",
    "Displayed Panel C genes (annotation only)"
  ),
  Value = as.character(
    c(
      analysis_set,
      if (
        is.null(host_expression_object_name) ||
        length(host_expression_object_name) == 0
      ) {
        "Not specified"
      } else {
        host_expression_object_name
      },
      if (
        is.null(host_deg_source) ||
        length(host_deg_source) == 0
      ) {
        "Not specified"
      } else {
        host_deg_source
      },
      gene_statistics_available,
      length(common_samples),
      ncol(metabolite_raw),
      nrow(host_deg),
      ncol(host_expression),
      nrow(metabolite_gene_correlations),
      sum(
        metabolite_gene_correlations$spearman_p < 0.10,
        na.rm = TRUE
      ),
      nrow(significant_edges),
      sum(
        metabolite_gene_correlations$spearman_p <
          spearman_strong_p_cutoff &
          abs(metabolite_gene_correlations$spearman_rho) >=
          spearman_rho_cutoff,
        na.rm = TRUE
      ),
      sum(
        metabolite_gene_correlations$spearman_FDR <
          spearman_relaxed_fdr_cutoff &
          abs(metabolite_gene_correlations$spearman_rho) >=
          spearman_rho_cutoff,
        na.rm = TRUE
      ),
      sum(
        significant_edges$spearman_high_confidence,
        na.rm = TRUE
      ),
      sum(
        significant_edges$adjusted_support,
        na.rm = TRUE
      ),
      sum(
        displayed_edges$adjusted_plot_support,
        na.rm = TRUE
      ),
      "butyrate" %in%
        selected_metabolites$Metabolite_key,
      "lithocholic acid" %in%
        selected_metabolites$Metabolite_key,
      if (
        length(adjustment_covariates_used) > 0
      ) {
        paste(
          adjustment_covariates_used,
          collapse = ", "
        )
      } else {
        "None"
      },
      nrow(biological_edges),
      sum(
        priority_gene_diagnostics$has_eligible_edge,
        na.rm = TRUE
      ),
      sum(
        priority_gene_diagnostics$selected_for_figure,
        na.rm = TRUE
      ),
      nrow(selected_metabolites),
      nrow(selected_genes),
      sum(
        selected_genes$DEG_side == positive_deg_label,
        na.rm = TRUE
      ),
      sum(
        selected_genes$DEG_side == negative_deg_label,
        na.rm = TRUE
      ),
      nrow(displayed_edges),
      sum(
        selected_genes$Panel_C_gene,
        na.rm = TRUE
      )
    )
  )
)

write.csv(
  diagnostics,
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/tables/",
    "Fig5D_8.4_diagnostics_",
    analysis_set,
    ".csv"
  ),
  row.names = FALSE
)

print(
  diagnostics,
  n = Inf
)

save(
  analysis_set,
  host_expression_object_name,
  host_deg_object_name,
  host_deg_source,
  gene_statistics_available,
  adjustment_covariates_used,
  diagnostics,
  host_deg,
  metabolite_floor_stats,
  metabolite_node_stats,
  metabolite_gene_correlations,
  significant_edges,
  biological_edges,
  priority_gene_diagnostics,
  metabolite_selection,
  selected_metabolites,
  selected_genes,
  displayed_edges,
  edge_curve_data,
  category_palette,
  manual_gene_category,
  host_function_rules,
  manual_metabolite_exclude,
  p_metabolite_host_rail,
  file = paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
    "Fig5D_9.8.4_metabolite_host_DEG_rail_",
    analysis_set,
    ".RData"
  )
)

message(
  "De novo metabolite–host DEG analysis completed: ",
  ncol(metabolite_raw),
  " metabolites × ",
  ncol(host_expression),
  " DEGs; ",
  nrow(significant_edges),
  " nominal-P eligible edges; ",
  sum(
    displayed_edges$adjusted_plot_support,
    na.rm = TRUE
  ),
  " LOO-stable adjusted-support markers; ",
  nrow(displayed_edges),
  " displayed edges."
)
