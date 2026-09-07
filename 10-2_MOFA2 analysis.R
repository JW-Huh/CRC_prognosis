#-----------------------------------------------------------------#
# Four-omics MOFA2 analysis (publication atlas, full-feature preprocessing, clustered heatmaps, and factor-pair PERMANOVA)
#
# Biological sequence:
#   Species -> KEGG ortholog -> Metabolite -> Host RNA-seq
#
# Input:
#   input/coherence_species_ko_metabolite_host_data_
#   metabolite50_host50pct_no_scaling.RData
#
# Analysis objective:
#   - infer unsupervised latent factors shared across multiple omics
#   - distinguish broadly shared, pair-shared, and view-specific factors
#   - retain samples with partially missing assays by encoding missing views as NA
#   - export factor scores, variance explained, and feature-loading rankings
#
# Important interpretation:
#   - pCR/non-pCR and Timepoint are not used to train the model
#   - these variables are attached only after model fitting for interpretation
#   - SubjectID-Timepoint observations are treated as samples by standard MOFA;
#     repeated-measure dependence should be considered in downstream testing
#   - Species and KO originate from the same metagenomic library, so factors
#     shared only by these two views may partly reflect common measurement origin
#
# Python environment:
#   Run 7-12.0_Setup_MOFA2_conda_windows_R43.R once before this script.
#   This analysis deliberately disables the legacy basilisk environment and
#   uses a pinned conda-forge environment to avoid h5py/HDF5 DLL conflicts.
#
#-----------------------------------------------------------------#
#
# v50 standalone end-to-end workflow
#   - starts only from the coherence_data RData written by the input-processing script
#   - constructs MOFA matrices and validates/reuses the established v30 HDF5 model
#     because v50 changes interpretation, statistical screening, and publication figures
#     rather than the fitted latent model
#   - adds active-view-consistent feature displays, response- and feature-informed
#     network selection, exhaustive Factor7 pair scanning, and three-factor PERMANOVA
#   - never loads v16-v29 analysis/checkpoint RData files
#   - writes the focused publication_v50 figure and result set
#

rm(list = ls())
options(stringsAsFactors = FALSE)
setwd("D:/2-연구/2-CRC metagenomics/")

# The Python interpreter must be selected before the MOFA2 namespace is loaded.
# Do not use requireNamespace("MOFA2") before this configuration block:
# MOFA2 1.12.1 can initialize its bundled basilisk Python during namespace load.
mofa_conda_env <- "mofa2_py310_070"

Sys.setenv(
  RETICULATE_AUTOCONFIGURE = "FALSE",
  RETICULATE_USE_MANAGED_VENV = "no"
)

required_packages <- c(
  "MOFA2",
  "reticulate",
  "BiocManager",
  "dplyr",
  "tidyr",
  "stringr",
  "readr",
  "tidyselect",
  "ggplot2",
  "patchwork",
  "svglite",
  "scales",
  "nlme",
  "MASS",
  "vegan",
  "permute",
  "ggbeeswarm"
)

# Revalidate the per-session temporary directory. On Windows/RStudio the
# session temp directory can occasionally be removed while R is still open;
# functions that cache metadata in tempdir() then fail with a libloc_*.rds
# connection error. tempdir(check = TRUE) recreates it when necessary.
mofa_session_tempdir <- tempdir(check = TRUE)

if (!dir.exists(mofa_session_tempdir)) {
  dir.create(
    mofa_session_tempdir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

if (!dir.exists(mofa_session_tempdir)) {
  stop(
    "R's session temporary directory is unavailable: ",
    mofa_session_tempdir,
    ". Restart R/RStudio and check TEMP/TMP permissions.",
    call. = FALSE
  )
}

message(
  "Validated R session temporary directory: ",
  normalizePath(
    mofa_session_tempdir,
    winslash = "/",
    mustWork = TRUE
  )
)

# Check only the packages required by this workflow. system.file(package = ...)
# locates installed packages through find.package() without loading their
# namespaces, so MOFA2's basilisk environment is not initialized here. It also
# avoids installed.packages()' metadata cache in tempdir().
package_is_installed <- vapply(
  required_packages,
  function(package_name) {
    nzchar(
      system.file(
        package = package_name
      )
    )
  },
  logical(1)
)

missing_packages <- required_packages[
  !package_is_installed
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

# The environment was created by a conda executable embedded under the
# basilisk cache. In a later R session, reticulate may discover a different
# conda installation, so condaenv_exists(<name>) can return a false negative.
# Use the exact Python path recorded by the one-time setup script instead.
mofa_python_override <- NULL
mofa_setup_record <-
  "results/mofa/MOFA_conda_environment_setup.txt"

mofa_python_candidates <- character(0)

if (
  !is.null(mofa_python_override) &&
  nzchar(mofa_python_override)
) {
  mofa_python_candidates <- c(
    mofa_python_candidates,
    mofa_python_override
  )
}

if (file.exists(mofa_setup_record)) {
  mofa_setup_lines <- readLines(
    mofa_setup_record,
    warn = FALSE
  )

  mofa_recorded_python <- sub(
    "^python=",
    "",
    mofa_setup_lines[
      grepl(
        "^python=",
        mofa_setup_lines
      )
    ]
  )

  mofa_python_candidates <- c(
    mofa_python_candidates,
    mofa_recorded_python
  )
}

# Fallback discovery for the basilisk conda location used by the setup script.
if (.Platform$OS.type == "windows") {
  mofa_basilisk_cache <- file.path(
    Sys.getenv("LOCALAPPDATA"),
    "R",
    "cache",
    "R",
    "basilisk"
  )

  if (dir.exists(mofa_basilisk_cache)) {
    mofa_python_candidates <- c(
      mofa_python_candidates,
      Sys.glob(
        file.path(
          mofa_basilisk_cache,
          "*",
          "0",
          "envs",
          mofa_conda_env,
          "python.exe"
        )
      )
    )
  }
}

mofa_python_candidates <- unique(
  mofa_python_candidates[
    nzchar(mofa_python_candidates) &
    file.exists(mofa_python_candidates)
  ]
)

if (length(mofa_python_candidates) == 0) {
  stop(
    paste0(
      "The previously created MOFA Python executable could not be found. ",
      "The setup script should only be run once. First check whether ",
      mofa_setup_record,
      " exists and contains a valid 'python=' path. If the environment was ",
      "moved or the R cache was deleted, set mofa_python_override near the ",
      "top of this script to the existing python.exe path."
    ),
    call. = FALSE
  )
}

mofa_python <- normalizePath(
  mofa_python_candidates[1],
  winslash = "/",
  mustWork = TRUE
)

message(
  "Using persisted MOFA Python: ",
  mofa_python
)

mofa_conda_prefix <- dirname(
  mofa_python
)

# Windows does not automatically reproduce all DLL-search changes that occur
# during `conda activate` when Python is embedded inside R through reticulate.
# Register every relevant conda DLL directory before NumPy, SciPy, or h5py is
# imported. The Python-side handles are retained for the whole R session.
mofa_dll_directories <- unique(
  c(
    mofa_conda_prefix,
    file.path(
      mofa_conda_prefix,
      "DLLs"
    ),
    file.path(
      mofa_conda_prefix,
      "Library",
      "mingw-w64",
      "bin"
    ),
    file.path(
      mofa_conda_prefix,
      "Library",
      "usr",
      "bin"
    ),
    file.path(
      mofa_conda_prefix,
      "Library",
      "bin"
    ),
    file.path(
      mofa_conda_prefix,
      "Scripts"
    ),
    file.path(
      mofa_conda_prefix,
      "bin"
    )
  )
)

mofa_dll_directories <- mofa_dll_directories[
  dir.exists(
    mofa_dll_directories
  )
]

if (.Platform$OS.type == "windows") {
  Sys.setenv(
    CONDA_PREFIX = mofa_conda_prefix,
    CONDA_DEFAULT_ENV = mofa_conda_env,
    CONDA_DLL_SEARCH_MODIFICATION_ENABLE = "1",
    PATH = paste(
      c(
        mofa_dll_directories,
        Sys.getenv("PATH")
      ),
      collapse = .Platform$path.sep
    )
  )
}

# RETICULATE_PYTHON is prescriptive and overrides later environment discovery.
Sys.setenv(
  RETICULATE_PYTHON = mofa_python
)

# A Python binary path is being supplied, so use_python() is clearer than
# use_condaenv(). This call requests the interpreter but does not import MOFA2.
reticulate::use_python(
  mofa_python,
  required = TRUE
)

# Python 3.8+ exposes os.add_dll_directory() specifically for dependencies of
# imported extension modules. This is more reliable than PATH alone for an
# embedded conda Python on Windows. Keep returned handles in a Python global;
# closing or garbage-collecting a handle removes the corresponding directory.
if (
  .Platform$OS.type == "windows" &&
  length(
    mofa_dll_directories
  ) > 0
) {
  mofa_dll_python_literal <- paste0(
    "[",
    paste(
      sprintf(
        "r'%s'",
        gsub(
          "'",
          "\\'",
          normalizePath(
            mofa_dll_directories,
            winslash = "/",
            mustWork = TRUE
          ),
          fixed = TRUE
        )
      ),
      collapse = ", "
    ),
    "]"
  )

  reticulate::py_run_string(
    paste0(
      "import os\n",
      "_mofa_dll_directory_handles = []\n",
      "for _mofa_dll_dir in ",
      mofa_dll_python_literal,
      ":\n",
      "    if os.path.isdir(_mofa_dll_dir):\n",
      "        _mofa_dll_directory_handles.append(",
      "os.add_dll_directory(_mofa_dll_dir))\n"
    ),
    local = FALSE,
    convert = FALSE
  )
}

# Fail before loading MOFA2 or processing data if NumPy, h5py, or the mofapy2
# entry point cannot be imported from the selected environment.
mofa_python_preflight <- tryCatch(
  {
    reticulate::import("numpy", convert = FALSE)
    reticulate::import("h5py", convert = FALSE)
    reticulate::import(
      "mofapy2.run.entry_point",
      convert = FALSE
    )
    TRUE
  },
  error = function(e) {
    stop(
      paste0(
        "MOFA Python preflight failed in environment '",
        mofa_conda_env,
        "'. The selected interpreter exists, but a compiled Python dependency ",
        "could not load inside the R process. Run this script in a newly ",
        "restarted R session. Original error: ",
        conditionMessage(e)
      ),
      call. = FALSE
    )
  }
)

mofa_active_python <- normalizePath(
  reticulate::py_config()$python,
  winslash = "/",
  mustWork = TRUE
)

if (!identical(
  tolower(mofa_active_python),
  tolower(mofa_python)
)) {
  stop(
    paste0(
      "reticulate initialized an unexpected Python interpreter. Expected: ",
      mofa_python,
      "; active: ",
      mofa_active_python,
      ". Close all RStudio windows and run this script in a clean session."
    ),
    call. = FALSE
  )
}

# Load MOFA2 only after reticulate has been bound to the external conda Python.
suppressPackageStartupMessages({
  library(MOFA2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(nlme)
})

dir.create(
  "results/mofa",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures/mofa",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures/mofa/internal_v16",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures/mofa/publication_v16",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/mofa/publication_v16",
  recursive = TRUE,
  showWarnings = FALSE
)

mofa_input_rdata <- paste0(
  "input/coherence_species_ko_metabolite_host_data_",
  "metabolite50_host50pct_no_scaling.RData"
)

if (!file.exists(mofa_input_rdata)) {
  stop(
    "MOFA input-processing output was not found: ",
    mofa_input_rdata,
    call. = FALSE
  )
}

load(mofa_input_rdata)
message("Loaded standalone MOFA input: ", mofa_input_rdata)

if (!exists("coherence_data")) {
  stop(
    "The input RData does not contain coherence_data.",
    call. = FALSE
  )
}


if (
  !is.null(coherence_data$schema_version) &&
  !identical(as.character(coherence_data$schema_version), "3.5")
) {
  warning(
    "The input coherence_data schema version is ",
    as.character(coherence_data$schema_version),
    "; this script was reviewed against schema version 3.5.",
    call. = FALSE
  )
}

required_views <- c(
  "species",
  "ko",
  "metabolite",
  "host"
)

if (
  is.null(coherence_data$modalities) ||
  !all(required_views %in% names(coherence_data$modalities))
) {
  stop(
    "coherence_data$modalities does not contain all four required views.",
    call. = FALSE
  )
}

mofa_input_preflight <- dplyr::bind_rows(
  lapply(
    required_views,
    function(view_name) {
      modality <- coherence_data$modalities[[view_name]]
      required_modality_fields <- c(
        "sample_meta",
        "transformed",
        "table"
      )
      missing_modality_fields <- setdiff(
        required_modality_fields,
        names(modality)
      )

      metadata_columns <- if (
        is.null(modality$sample_meta)
      ) {
        character(0)
      } else {
        colnames(modality$sample_meta)
      }

      missing_metadata_columns <- setdiff(
        c(
          "SampleID",
          "SubjectID",
          "Timepoint",
          "TRG_plot"
        ),
        metadata_columns
      )

      matrix_for_check <- if (
        !is.null(modality$table)
      ) {
        as.matrix(modality$table)
      } else if (
        !is.null(modality$transformed)
      ) {
        as.matrix(modality$transformed)
      } else {
        matrix(
          numeric(0),
          nrow = 0,
          ncol = 0
        )
      }

      data.frame(
        view = view_name,
        missing_fields = paste(
          missing_modality_fields,
          collapse = ";"
        ),
        missing_metadata_columns = paste(
          missing_metadata_columns,
          collapse = ";"
        ),
        n_samples_metadata = if (
          is.null(modality$sample_meta)
        ) {
          0L
        } else {
          nrow(modality$sample_meta)
        },
        n_samples_matrix = nrow(matrix_for_check),
        n_features_matrix = ncol(matrix_for_check),
        sample_names_unique =
          !is.null(rownames(matrix_for_check)) &&
          anyDuplicated(rownames(matrix_for_check)) == 0,
        feature_names_unique =
          !is.null(colnames(matrix_for_check)) &&
          anyDuplicated(colnames(matrix_for_check)) == 0,
        passed =
          length(missing_modality_fields) == 0 &&
          length(missing_metadata_columns) == 0 &&
          nrow(matrix_for_check) > 0 &&
          ncol(matrix_for_check) > 0 &&
          !is.null(rownames(matrix_for_check)) &&
          !is.null(colnames(matrix_for_check)) &&
          anyDuplicated(rownames(matrix_for_check)) == 0 &&
          anyDuplicated(colnames(matrix_for_check)) == 0,
        stringsAsFactors = FALSE
      )
    }
  )
)

if (!all(mofa_input_preflight$passed)) {
  write.csv(
    mofa_input_preflight,
    "results/mofa/MOFA_v30_input_preflight.csv",
    row.names = FALSE
  )

  stop(
    paste0(
      "The coherence_data input failed the v30 preflight. Inspect ",
      "results/mofa/MOFA_v30_input_preflight.csv."
    ),
    call. = FALSE
  )
}


#=================================================================#
# 1. Analysis settings
#=================================================================#

set.seed(20260801)

view_labels <- c(
  species = "Species",
  ko = "KEGG ortholog",
  metabolite = "Metabolite",
  host = "Host RNA-seq"
)

view_colors <- c(
  Species = "#79B3A3",
  `KEGG ortholog` = "#D8A46F",
  Metabolite = "#82A8C7",
  `Host RNA-seq` = "#B29AC6"
)

mofa_response_colors <- c(
  pCR = "#4FAE9A",
  non_pCR = "#DE7872",
  CR = "#4FAE9A",
  nonCR = "#DE7872"
)

timepoint_display_labels <- c(
  Before = "Baseline",
  Ongoing = "After RT"
)

timepoint_display_colors <- c(
  Baseline = "#667A8A",
  `After RT` = "#C59A5B"
)

# A sample is retained when at least this many assays are measured.
# Missing assays are represented by NA and are handled by MOFA2.
mofa_min_views_per_sample <- 2

# Primary model scope. The all-time model describes global molecular structure;
# response testing is performed on baseline samples and paired changes below.
# This avoids interpreting duplicated SubjectID-Timepoint observations as
# independent evidence for pCR/non-pCR separation.
mofa_primary_training_scope <- "all_time"

# Primary preprocessing is deliberately aligned with the feature matrices used
# to construct the Procrustes input object. This avoids introducing a second,
# unvalidated transformation solely for MOFA and makes discrepancies between
# Procrustes and MOFA attributable to the integration model rather than to a
# different preprocessing pipeline.
#
# Primary profile:
#   Species / KO : prevalence-QC relative abundance -> square-root (Hellinger)
#                  -> feature-wise z-score, exactly as stored in $table.
#   Metabolite   : detection-floor-QC concentration -> feature-wise z-score
#                  (the upstream abs_z profile), exactly as stored in $table.
#   Host RNA-seq : blind DESeq2 VST, cumulative-variance primary gene set,
#                  no gene-wise z-scaling, exactly as stored in $table.
#
# Alternative transformations remain available as explicit sensitivity models;
# they are not selected adaptively feature by feature in the primary model.
mofa_primary_preprocessing_profile <- "procrustes_aligned"
mofa_microbiome_transform <- "procrustes_hellinger_z"
# Alternatives: "simple_multiplicative_clr", "bayesian_count_clr"
mofa_microbiome_zero_replacement_fraction <- 0.65
mofa_metabolite_transform <- "procrustes_abs_z"
# Alternatives: "log10_z", "pareto"
# Retained only for the legacy diagnostic helper; no adaptive transformation is
# applied in the primary model.
mofa_metabolite_skewness_threshold <- 1.00
mofa_host_feature_scaling <- "none"
# Alternative: "z_score"

# Alternative preprocessing profiles are exposed above as explicit settings.
# Each alternative must use a distinct HDF5 path and be fitted as a sensitivity
# model; the primary script does not switch transformations adaptively.

# Primary feature rule for v16: within each upstream Procrustes feature universe,
# rank non-constant features by pre-MOFA variance. Large views retain the smallest
# feature set reaching the balanced cumulative-variance target; views with 100 or fewer variable
# features are retained in full. Large views use a balanced cumulative-variance
# rule with explicit feature caps. This reduces unequal likelihood contribution
# from KO/host while preserving the small metabolite view. Remaining imbalance is
# handled by view scaling, ARD, spike-and-slab loading shrinkage, and QC audits.
mofa_feature_selection_mode <- "balanced_cumulative_variance"
mofa_host_use_all_qc_genes <- FALSE

# The curves below remain diagnostic. The primary target is 95%; a view-specific
# maximum can stop earlier when retaining the full 95% tail would recreate a
# severe dimensional imbalance.
mofa_feature_sensitivity_fractions <- c(
  0.50,
  0.75,
  0.90,
  0.95,
  0.99,
  1.00
)

# Optional alternative model settings are exported for transparent sensitivity
# analyses. They do not affect the primary model unless
# mofa_feature_selection_mode is changed manually.
mofa_hvf_variance_quantile <- c(
  species = 0.50,
  ko = 0.50,
  metabolite = 0.00,
  host = 0.50
)
mofa_hvf_min_features <- c(
  species = 50,
  ko = 300,
  metabolite = 1,
  host = 500
)
mofa_hvf_max_features <- c(
  species = 500,
  ko = 1500,
  metabolite = Inf,
  host = 2000
)
mofa_keep_all_if_n_features_at_most <- 100

# Primary feature retention:
#   - views with >100 variable features: retain the smallest variance-ranked
#     feature set reaching 95% cumulative variance, subject to the min/max caps;
#   - views with <=100 variable features: retain all non-constant features.
mofa_cumulative_variance_target <- 0.95

mofa_min_observed_samples <- 10
mofa_initial_factors <- 12
mofa_drop_factor_threshold <- 0.005
mofa_use_spikeslab_weights <- TRUE
# Publication-level active-view calls require R2 >=2% in a view. A 3% threshold
# is retained as a stricter audit. These thresholds affect interpretation and
# display only, not model fitting.
mofa_active_view_r2 <- 0.02
mofa_strict_active_view_r2 <- 0.03
mofa_publication_min_active_views <- 2
mofa_top_features_per_direction <- 5

# Longitudinal and response-analysis settings.
mofa_factor_permutations <- 1999
mofa_pair_permutations <- 999
mofa_min_paired_subjects <- 6
mofa_min_paired_per_response <- 3
mofa_subject_specificity_icc_threshold <- 0.50
mofa_selected_factor_override <- NULL
mofa_selected_pair_override <- NULL
mofa_heatmap_response_order <- c("pCR", "non_pCR")
mofa_heatmap_timepoint_order <- c("Before", "Ongoing")
mofa_heatmap_effect_metric <- "hedges_g"
mofa_pair_permanova_permutations <- 999
mofa_pair_maps_to_show <- 2
mofa_publication_pair_count <- 2
mofa_network_factor_count <- 5
mofa_network_features_per_view_direction <- 2
mofa_network_shared_loading_threshold <- 0.35

# Optional sensitivity models. These use the same selected features but train
# separate MOFA models on baseline, ongoing, subject-mean, and paired-delta
# matrices. They are one-seed sensitivity fits and are reused after creation.
# FALSE is the practical first-run default because the full host feature universe
# makes these additional fits computationally expensive. Set TRUE after the
# primary v16 model and post-hoc figures have been checked.
mofa_run_scope_sensitivity_models <- FALSE
# A four-view complete-case model directly tests whether missing-view structure
# fragmented shared factors in the all-available model. It is fitted with one
# seed and re-used after creation.
mofa_run_four_view_complete_model <- FALSE
mofa_scope_sensitivity_seed <- 20260831
mofa_scope_min_samples <- 10
mofa_scope_initial_factors <- 10

# Main-display filters affect labels/figures only, never model fitting or the
# complete exported loading table. "Uncharacterized" is filtered only when an
# explicit annotation description contains that term; symbols such as C10orf99
# are not removed solely because their name contains "orf".
mofa_exclude_pseudogene_from_display <- TRUE
mofa_exclude_readthrough_from_display <- TRUE
mofa_exclude_uncharacterized_from_display <- TRUE

mofa_host_annotation_file <-
  "host_RNAseq/TNT_Expression_Profile.GRCh38.gene.csv"

mofa_host_response_rdata_candidates <- c(
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
    "Fig5_host_RNAseq_inputs.RData"
  ),
  paste0(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/",
    "Fig5B_C_downstream_inputs.RData"
  )
)

# Multiple initializations reduce dependence on a single variational solution.
mofa_seeds <- c(20260801)

mofa_use_basilisk <- FALSE
mofa_save_training_data <- TRUE

# Retain the validated core calculations, but write only the focused publication_v30 figures.
mofa_save_legacy_figures <- FALSE


# v30 uses a dedicated model path. Older v14-v29 HDF5/RData outputs are never
# treated as inputs because they may encode different samples, features, or
# preprocessing settings.
mofa_best_model_path <-
  "results/mofa/MOFA_4omics_v30_balanced_primary.hdf5"
mofa_seed_summary_path <-
  "results/mofa/MOFA_v30_balanced_primary_seed_ELBO_summary.csv"
mofa_reuse_existing_best_model <- FALSE




mofa_literature_consensus <- data.frame(
  study = c(
    "Argelaguet et al., Molecular Systems Biology 2018 / MOFA2 guidance",
    "Mohr et al., Nature Communications 2024",
    "van den Berg et al., BMC Genomics 2006",
    "Martin-Fernandez et al., Statistical Modelling 2015",
    "Current CRC study: Procrustes-aligned primary MOFA"
  ),
  data_context = c(
    "General multi-omics factor-analysis framework",
    "Longitudinal 16S microbiome plus plasma metabolome; approximately 40 participants",
    "Metabolomics preprocessing benchmark",
    "Compositional count-zero replacement methodology",
    "Shotgun species, KO, targeted fecal metabolites, host bulk RNA-seq"
  ),
  preprocessing = c(
    paste(
      "Remove library-size and technical effects; variance-stabilize count data;",
      "filter highly variable features; consider stronger filtering for much larger views"
    ),
    paste(
      "ASVs filtered, aggregated to genus, CLR-scaled; metabolites QC-filtered,",
      "log10-transformed and Pareto-scaled; 53 taxa and 138 metabolites entered MOFA"
    ),
    paste(
      "Compared autoscaling, Pareto, range, VAST and transformations;",
      "no universally optimal scaling and autoscaling/range performed well in their data"
    ),
    paste(
      "Bayesian-multiplicative replacement for count zeros;",
      "non-zero components are multiplicatively adjusted to preserve closure"
    ),
    paste(
      "Use the exact upstream Procrustes feature tables as primary input;",
      "evaluate CLR, metabolite log/Pareto and host z-scaling only as sensitivities"
    )
  ),
  implication_for_current_analysis = c(
    "Do not interpret view-specific factors as failure; quantify feature imbalance and sharedness explicitly",
    "Genus aggregation was a dimension/sparsity choice for 16S ASVs, not a universal requirement for shotgun species",
    "Keep z-scaling as primary for consistency; use Pareto only as an explicit metabolite sensitivity",
    "Avoid half-minimum replacement as the only analysis; Bayesian count replacement requires integer-like counts",
    "Differences from Procrustes now reflect model estimands rather than incompatible preprocessing"
  ),
  stringsAsFactors = FALSE
)


#=================================================================#
# 2. Sample metadata and assay availability
#=================================================================#

mofa_metadata_long <- dplyr::bind_rows(
  lapply(
    required_views,
    function(view_name) {
      coherence_data$modalities[[view_name]]$sample_meta %>%
        dplyr::transmute(
          SampleID = as.character(SampleID),
          SubjectID = as.character(SubjectID),
          Timepoint = as.character(Timepoint),
          TRG_plot = as.character(TRG_plot),
          view = view_name
        )
    }
  )
)

mofa_metadata_conflicts <- mofa_metadata_long %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(
    n_subject = dplyr::n_distinct(
      SubjectID[!is.na(SubjectID)]
    ),
    n_timepoint = dplyr::n_distinct(
      Timepoint[!is.na(Timepoint)]
    ),
    n_trg = dplyr::n_distinct(
      TRG_plot[!is.na(TRG_plot)]
    ),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    n_subject > 1 |
      n_timepoint > 1 |
      n_trg > 1
  )

if (nrow(mofa_metadata_conflicts) > 0) {
  write.csv(
    mofa_metadata_conflicts,
    "results/mofa/MOFA_sample_metadata_conflicts.csv",
    row.names = FALSE
  )

  stop(
    "Conflicting metadata were found for standardized SampleID values. ",
    "Inspect results/mofa/MOFA_sample_metadata_conflicts.csv.",
    call. = FALSE
  )
}

mofa_sample_metadata <- mofa_metadata_long %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(
    SubjectID = dplyr::first(
      SubjectID[!is.na(SubjectID)]
    ),
    Timepoint = dplyr::first(
      Timepoint[!is.na(Timepoint)]
    ),
    TRG_plot = dplyr::first(
      TRG_plot[!is.na(TRG_plot)]
    ),
    .groups = "drop"
  )

mofa_availability_long <- tidyr::expand_grid(
  SampleID = mofa_sample_metadata$SampleID,
  view = required_views
) %>%
  dplyr::left_join(
    mofa_metadata_long %>%
      dplyr::distinct(
        SampleID,
        view
      ) %>%
      dplyr::mutate(
        available = 1L
      ),
    by = c(
      "SampleID",
      "view"
    )
  ) %>%
  dplyr::mutate(
    available = dplyr::coalesce(
      available,
      0L
    )
  )

mofa_availability_wide <- mofa_availability_long %>%
  tidyr::pivot_wider(
    names_from = view,
    values_from = available
  ) %>%
  dplyr::mutate(
    n_views = rowSums(
      dplyr::across(
        dplyr::all_of(required_views)
      )
    )
  )

mofa_sample_ids <- mofa_availability_wide$SampleID[
  mofa_availability_wide$n_views >=
    mofa_min_views_per_sample
]

if (length(mofa_sample_ids) < 15) {
  stop(
    "Fewer than 15 samples have at least ",
    mofa_min_views_per_sample,
    " measured views. MOFA inference would be poorly supported.",
    call. = FALSE
  )
}

mofa_sample_metadata <- mofa_sample_metadata[
  match(
    mofa_sample_ids,
    mofa_sample_metadata$SampleID
  ),
  ,
  drop = FALSE
]

mofa_sample_metadata$Timepoint <- factor(
  mofa_sample_metadata$Timepoint,
  levels = c(
    "Before",
    "Ongoing"
  )
)

mofa_sample_metadata$TRG_plot <- factor(
  mofa_sample_metadata$TRG_plot,
  levels = c(
    "non_pCR",
    "pCR"
  )
)

mofa_availability_long <- mofa_availability_long %>%
  dplyr::filter(
    SampleID %in% mofa_sample_ids
  ) %>%
  dplyr::mutate(
    view_label = factor(
      unname(
        view_labels[view]
      ),
      levels = rev(
        unname(view_labels)
      )
    ),
    SampleID = factor(
      SampleID,
      levels = mofa_sample_metadata$SampleID
    )
  )

write.csv(
  mofa_sample_metadata,
  "results/mofa/MOFA_v30_sample_metadata.csv",
  row.names = FALSE
)


#=================================================================#
# 3. MOFA-specific feature selection and matrix construction
#=================================================================#
#
# The coherence input stores samples in rows and features in columns.
# This section uses the transformed, non-ordination matrices:
#   Species / KO: square-root relative abundance
#   Metabolite: transformed post-QC concentration matrix
#   Host: blind DESeq2 VST expression
#
# Features are selected by variance within each view, then standardized to
# mean 0 and SD 1 across measured samples. View-level missingness is retained
# as NA after expanding each matrix to the union of retained samples.

calculate_sample_skewness <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) < 4) {
    return(NA_real_)
  }

  x_sd <- stats::sd(x)

  if (!is.finite(x_sd) || x_sd <= 0) {
    return(0)
  }

  mean(
    (
      (x - mean(x)) /
        x_sd
    ) ^ 3
  )
}

multiplicative_replace_zeros <- function(
    data_matrix,
    replacement_fraction = 0.50
) {
  data_matrix <- as.matrix(data_matrix)
  storage.mode(data_matrix) <- "numeric"

  replaced_matrix <- data_matrix
  replacement_values <- rep(
    NA_real_,
    nrow(data_matrix)
  )

  for (sample_index in seq_len(nrow(data_matrix))) {
    x <- data_matrix[sample_index, ]
    x[!is.finite(x)] <- 0

    sample_total <- sum(x)

    if (!is.finite(sample_total) || sample_total <= 0) {
      next
    }

    x <- x / sample_total
    zero_index <- x <= 0
    positive_x <- x[!zero_index]

    if (!any(zero_index)) {
      replaced_matrix[sample_index, ] <- x
      replacement_values[sample_index] <- 0
      next
    }

    if (length(positive_x) == 0) {
      next
    }

    zero_count <- sum(zero_index)
    replacement_value <- min(
      positive_x,
      na.rm = TRUE
    ) * replacement_fraction

    replacement_value <- min(
      replacement_value,
      0.50 / zero_count
    )

    remaining_mass <- 1 -
      zero_count * replacement_value

    if (!is.finite(remaining_mass) || remaining_mass <= 0) {
      replacement_value <- 0.25 / zero_count
      remaining_mass <- 0.75
    }

    x[zero_index] <- replacement_value
    x[!zero_index] <-
      x[!zero_index] /
      sum(x[!zero_index]) *
      remaining_mass

    replaced_matrix[sample_index, ] <- x
    replacement_values[sample_index] <- replacement_value
  }

  list(
    matrix = replaced_matrix,
    replacement_values = replacement_values
  )
}

adaptive_transform_metabolites <- function(data_matrix) {
  data_matrix <- as.matrix(data_matrix)
  storage.mode(data_matrix) <- "numeric"
  transformed_matrix <- data_matrix

  diagnostics <- data.frame(
    feature = colnames(data_matrix),
    skewness_before = NA_real_,
    skewness_after = NA_real_,
    transformation = "none",
    pseudocount = NA_real_,
    stringsAsFactors = FALSE
  )

  for (feature_index in seq_len(ncol(data_matrix))) {
    x <- data_matrix[, feature_index]
    finite_x <- x[is.finite(x)]

    diagnostics$skewness_before[feature_index] <-
      calculate_sample_skewness(finite_x)

    use_log <-
      identical(
        mofa_metabolite_transform,
        "adaptive_log"
      ) &&
      length(finite_x) >= 4 &&
      all(finite_x >= 0) &&
      is.finite(
        diagnostics$skewness_before[feature_index]
      ) &&
      diagnostics$skewness_before[feature_index] >=
        mofa_metabolite_skewness_threshold

    if (use_log) {
      positive_x <- finite_x[finite_x > 0]
      pseudocount <- if (length(positive_x) > 0) {
        min(positive_x) / 2
      } else {
        1e-08
      }

      x <- log10(x + pseudocount)
      diagnostics$transformation[feature_index] <- "log10"
      diagnostics$pseudocount[feature_index] <- pseudocount
    }

    if (anyNA(x)) {
      feature_median <- stats::median(
        x,
        na.rm = TRUE
      )

      if (!is.finite(feature_median)) {
        feature_median <- 0
      }

      x[is.na(x)] <- feature_median
    }

    transformed_matrix[, feature_index] <- x
    diagnostics$skewness_after[feature_index] <-
      calculate_sample_skewness(x)
  }

  list(
    matrix = transformed_matrix,
    diagnostics = diagnostics
  )
}


z_score_columns <- function(data_matrix) {
  data_matrix <- as.matrix(data_matrix)
  storage.mode(data_matrix) <- "numeric"

  feature_mean <- colMeans(
    data_matrix,
    na.rm = TRUE
  )
  feature_sd <- apply(
    data_matrix,
    2,
    stats::sd,
    na.rm = TRUE
  )

  keep <- is.finite(feature_sd) &
    feature_sd > 0
  data_matrix <- data_matrix[, keep, drop = FALSE]
  feature_mean <- feature_mean[keep]
  feature_sd <- feature_sd[keep]

  data_matrix <- sweep(
    data_matrix,
    2,
    feature_mean,
    FUN = "-"
  )
  data_matrix <- sweep(
    data_matrix,
    2,
    feature_sd,
    FUN = "/"
  )

  data_matrix
}

pareto_scale_columns <- function(data_matrix) {
  data_matrix <- as.matrix(data_matrix)
  storage.mode(data_matrix) <- "numeric"

  feature_mean <- colMeans(
    data_matrix,
    na.rm = TRUE
  )
  feature_sd <- apply(
    data_matrix,
    2,
    stats::sd,
    na.rm = TRUE
  )

  keep <- is.finite(feature_sd) &
    feature_sd > 0
  data_matrix <- data_matrix[, keep, drop = FALSE]
  feature_mean <- feature_mean[keep]
  feature_sd <- feature_sd[keep]

  data_matrix <- sweep(
    data_matrix,
    2,
    feature_mean,
    FUN = "-"
  )
  data_matrix <- sweep(
    data_matrix,
    2,
    sqrt(feature_sd),
    FUN = "/"
  )

  data_matrix
}

log10_z_columns <- function(data_matrix) {
  data_matrix <- as.matrix(data_matrix)
  storage.mode(data_matrix) <- "numeric"
  pseudocount <- rep(
    NA_real_,
    ncol(data_matrix)
  )

  for (feature_index in seq_len(ncol(data_matrix))) {
    x <- data_matrix[, feature_index]
    positive_x <- x[
      is.finite(x) &
        x > 0
    ]

    pseudocount[feature_index] <- if (
      length(positive_x) > 0
    ) {
      min(positive_x) / 2
    } else {
      1e-08
    }

    data_matrix[, feature_index] <- log10(
      x + pseudocount[feature_index]
    )
  }

  list(
    matrix = z_score_columns(data_matrix),
    pseudocount = stats::setNames(
      pseudocount,
      colnames(data_matrix)
    )
  )
}

bayesian_count_clr <- function(data_matrix) {
  if (!requireNamespace(
    "zCompositions",
    quietly = TRUE
  )) {
    stop(
      "Install zCompositions before using mofa_microbiome_transform = ",
      "'bayesian_count_clr'.",
      call. = FALSE
    )
  }

  data_matrix <- as.matrix(data_matrix)
  storage.mode(data_matrix) <- "numeric"

  finite_values <- data_matrix[
    is.finite(data_matrix)
  ]
  integer_like <- length(finite_values) > 0 &&
    max(
      abs(
        finite_values -
          round(finite_values)
      ),
      na.rm = TRUE
    ) < 1e-08

  if (!integer_like) {
    stop(
      paste0(
        "Bayesian GBM zero replacement is a count-data model, but the selected ",
        "microbiome matrix is not integer-like. Use the Procrustes-aligned ",
        "Hellinger profile or simple_multiplicative_clr instead."
      ),
      call. = FALSE
    )
  }

  imputed_proportions <- as.matrix(
    zCompositions::cmultRepl(
      round(data_matrix),
      label = 0,
      method = "GBM",
      output = "prop",
      suppress.print = TRUE
    )
  )

  log_matrix <- log(imputed_proportions)
  clr_matrix <- sweep(
    log_matrix,
    1,
    rowMeans(
      log_matrix,
      na.rm = TRUE
    ),
    FUN = "-"
  )

  z_score_columns(clr_matrix)
}

prepare_mofa_view <- function(view_name) {
  transformation_diagnostics <- data.frame()
  clr_replacement_summary <- data.frame()
  feature_scaling_description <- NA_character_

  if (view_name %in% c("species", "ko")) {
    variance_reference_matrix <- as.matrix(
      coherence_data$modalities[[view_name]]$transformed
    )

    if (identical(
      mofa_microbiome_transform,
      "procrustes_hellinger_z"
    )) {
      view_matrix <- as.matrix(
        coherence_data$modalities[[view_name]]$table
      )
      transformation_description <- paste0(
        "Procrustes-aligned: prevalence-QC relative abundance; ",
        "square-root (Hellinger); feature-wise z-score"
      )
      feature_scaling_description <- "upstream z-score"
    } else if (identical(
      mofa_microbiome_transform,
      "simple_multiplicative_clr"
    )) {
      relative_matrix <- as.matrix(
        coherence_data$modalities[[view_name]]$relative
      )
      replacement_result <- multiplicative_replace_zeros(
        relative_matrix,
        replacement_fraction =
          mofa_microbiome_zero_replacement_fraction
      )
      log_matrix <- log(
        replacement_result$matrix
      )
      clr_matrix <- sweep(
        log_matrix,
        1,
        rowMeans(
          log_matrix,
          na.rm = TRUE
        ),
        FUN = "-"
      )
      view_matrix <- z_score_columns(
        clr_matrix
      )
      variance_reference_matrix <- clr_matrix
      clr_replacement_summary <- data.frame(
        view = view_name,
        SampleID = rownames(view_matrix),
        zero_replacement =
          replacement_result$replacement_values,
        method = "simple multiplicative",
        stringsAsFactors = FALSE
      )
      transformation_description <- paste0(
        "Sensitivity: simple multiplicative replacement (fraction = ",
        mofa_microbiome_zero_replacement_fraction,
        "); CLR; feature-wise z-score"
      )
      feature_scaling_description <- "z-score after CLR"
    } else if (identical(
      mofa_microbiome_transform,
      "bayesian_count_clr"
    )) {
      count_matrix <- as.matrix(
        coherence_data$modalities[[view_name]]$filtered
      )
      view_matrix <- bayesian_count_clr(
        count_matrix
      )
      variance_reference_matrix <- view_matrix
      transformation_description <- paste0(
        "Sensitivity: zCompositions GBM Bayesian-multiplicative count-zero ",
        "replacement; CLR; feature-wise z-score"
      )
      feature_scaling_description <- "z-score after Bayesian CLR"
    } else {
      stop(
        "Unknown microbiome transformation: ",
        mofa_microbiome_transform,
        call. = FALSE
      )
    }
  } else if (view_name == "metabolite") {
    variance_reference_matrix <- as.matrix(
      coherence_data$modalities$metabolite$transformed
    )

    if (identical(
      mofa_metabolite_transform,
      "procrustes_abs_z"
    )) {
      view_matrix <- as.matrix(
        coherence_data$modalities$metabolite$table
      )
      transformation_description <- paste0(
        "Procrustes-aligned: detection-floor-QC concentration; ",
        "no log transformation; feature-wise z-score"
      )
      feature_scaling_description <- "upstream z-score"
    } else if (identical(
      mofa_metabolite_transform,
      "log10_z"
    )) {
      log_result <- log10_z_columns(
        variance_reference_matrix
      )
      view_matrix <- log_result$matrix
      transformation_description <- paste0(
        "Sensitivity: feature-specific half-minimum pseudocount; ",
        "log10 transformation; feature-wise z-score"
      )
      feature_scaling_description <- "log10 plus z-score"
    } else if (identical(
      mofa_metabolite_transform,
      "pareto"
    )) {
      view_matrix <- pareto_scale_columns(
        variance_reference_matrix
      )
      transformation_description <- paste0(
        "Sensitivity: mean-centering and Pareto scaling; no log transformation"
      )
      feature_scaling_description <- "Pareto scaling"
    } else {
      stop(
        "Unknown metabolite transformation: ",
        mofa_metabolite_transform,
        call. = FALSE
      )
    }

    transformation_diagnostics <- data.frame(
      view = view_name,
      feature = colnames(
        variance_reference_matrix
      ),
      skewness_input = vapply(
        seq_len(
          ncol(variance_reference_matrix)
        ),
        function(feature_index) {
          calculate_sample_skewness(
            variance_reference_matrix[, feature_index]
          )
        },
        numeric(1)
      ),
      primary_transformation =
        mofa_metabolite_transform,
      stringsAsFactors = FALSE
    )
  } else if (view_name == "host") {
    if (
      isTRUE(mofa_host_use_all_qc_genes) &&
      !is.null(
        coherence_data$host_sensitivity$vst_qc_all
      )
    ) {
      variance_reference_matrix <- as.matrix(
        coherence_data$host_sensitivity$vst_qc_all
      )
    } else {
      variance_reference_matrix <- as.matrix(
        coherence_data$modalities$host$transformed
      )
    }

    if (identical(
      mofa_host_feature_scaling,
      "none"
    )) {
      view_matrix <- variance_reference_matrix
      transformation_description <- paste0(
        "Procrustes-aligned: blind DESeq2 VST; cumulative-variance primary ",
        "host feature set; no gene-wise z-scaling"
      )
      feature_scaling_description <- "none; MOFA centers groups and scales views"
    } else if (identical(
      mofa_host_feature_scaling,
      "z_score"
    )) {
      view_matrix <- z_score_columns(
        variance_reference_matrix
      )
      transformation_description <- paste0(
        "Sensitivity: blind DESeq2 VST followed by gene-wise z-score"
      )
      feature_scaling_description <- "gene-wise z-score"
    } else {
      stop(
        "Unknown host feature scaling: ",
        mofa_host_feature_scaling,
        call. = FALSE
      )
    }
  } else {
    stop(
      "Unknown view: ",
      view_name,
      call. = FALSE
    )
  }

  storage.mode(view_matrix) <- "numeric"
  storage.mode(variance_reference_matrix) <- "numeric"

  measured_samples <- Reduce(
    intersect,
    list(
      mofa_sample_ids,
      rownames(view_matrix),
      rownames(variance_reference_matrix)
    )
  )

  if (length(measured_samples) < mofa_min_observed_samples) {
    stop(
      view_name,
      ": fewer than ",
      mofa_min_observed_samples,
      " measured samples remain.",
      call. = FALSE
    )
  }

  common_features <- intersect(
    colnames(view_matrix),
    colnames(variance_reference_matrix)
  )

  view_matrix <- view_matrix[
    measured_samples,
    common_features,
    drop = FALSE
  ]
  variance_reference_matrix <- variance_reference_matrix[
    measured_samples,
    common_features,
    drop = FALSE
  ]

  observed_n <- colSums(
    is.finite(view_matrix)
  )
  prepared_variance <- apply(
    view_matrix,
    2,
    stats::var,
    na.rm = TRUE
  )
  feature_variance <- apply(
    variance_reference_matrix,
    2,
    stats::var,
    na.rm = TRUE
  )

  keep_features <- names(feature_variance)[
    observed_n >= mofa_min_observed_samples &
      is.finite(prepared_variance) &
      prepared_variance > 0 &
      is.finite(feature_variance) &
      feature_variance > 0
  ]

  if (length(keep_features) == 0) {
    stop(
      view_name,
      ": no variable features remain for MOFA.",
      call. = FALSE
    )
  }

  keep_features <- keep_features[
    order(
      feature_variance[keep_features],
      decreasing = TRUE
    )
  ]

  variance_rank_table <- data.frame(
    view = view_name,
    feature = keep_features,
    variance_before_mofa_scaling = as.numeric(
      feature_variance[keep_features]
    ),
    variance_rank = seq_along(keep_features),
    stringsAsFactors = FALSE
  )
  variance_rank_table$cumulative_variance_fraction <-
    cumsum(
      variance_rank_table$variance_before_mofa_scaling
    ) /
    sum(
      variance_rank_table$variance_before_mofa_scaling
    )

  variance_cutoff <- NA_real_

  minimum_feature_count <- min(
    length(keep_features),
    as.integer(
      mofa_hvf_min_features[view_name]
    )
  )

  maximum_feature_count <-
    mofa_hvf_max_features[view_name]

  if (identical(
    mofa_feature_selection_mode,
    "all_qc_variable_features"
  )) {
    selected_features <- keep_features
  } else if (identical(
    mofa_feature_selection_mode,
    "balanced_cumulative_variance"
  )) {
    if (
      length(keep_features) <=
        mofa_keep_all_if_n_features_at_most
    ) {
      selected_feature_count <- length(keep_features)
    } else {
      selected_feature_count <- which(
        variance_rank_table$cumulative_variance_fraction >=
          mofa_cumulative_variance_target
      )[1]

      if (!is.finite(selected_feature_count)) {
        selected_feature_count <- length(keep_features)
      }

      selected_feature_count <- max(
        selected_feature_count,
        minimum_feature_count
      )

      if (is.finite(maximum_feature_count)) {
        selected_feature_count <- min(
          selected_feature_count,
          as.integer(maximum_feature_count)
        )
      }
    }

    selected_features <- keep_features[
      seq_len(selected_feature_count)
    ]

    variance_cutoff <- min(
      feature_variance[selected_features],
      na.rm = TRUE
    )
  } else {
    variance_cutoff <- if (
      length(keep_features) <=
        mofa_keep_all_if_n_features_at_most
    ) {
      -Inf
    } else {
      stats::quantile(
        feature_variance[keep_features],
        probs = mofa_hvf_variance_quantile[view_name],
        na.rm = TRUE,
        names = FALSE
      )
    }

    selected_features <- keep_features[
      feature_variance[keep_features] >=
        variance_cutoff
    ]

    if (length(selected_features) < minimum_feature_count) {
      selected_features <- keep_features[
        seq_len(minimum_feature_count)
      ]
    }

    if (
      is.finite(maximum_feature_count) &&
      length(selected_features) >
        maximum_feature_count
    ) {
      selected_features <- selected_features[
        seq_len(
          as.integer(maximum_feature_count)
        )
      ]
    }
  }

  variance_rank_table$selected <-
    variance_rank_table$feature %in%
    selected_features

  sensitivity_fraction_table <- dplyr::bind_rows(
    lapply(
      mofa_feature_sensitivity_fractions,
      function(fraction_value) {
        feature_count <- which(
          variance_rank_table$cumulative_variance_fraction >=
            fraction_value
        )[1]

        if (!is.finite(feature_count)) {
          feature_count <- nrow(
            variance_rank_table
          )
        }

        data.frame(
          view = view_name,
          cumulative_variance_target =
            fraction_value,
          n_features = feature_count,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  view_matrix <- view_matrix[
    ,
    selected_features,
    drop = FALSE
  ]

  expanded_matrix <- matrix(
    NA_real_,
    nrow = length(mofa_sample_ids),
    ncol = ncol(view_matrix),
    dimnames = list(
      mofa_sample_ids,
      colnames(view_matrix)
    )
  )
  expanded_matrix[
    rownames(view_matrix),
    ] <- view_matrix

  list(
    matrix = base::t(expanded_matrix),
    summary = data.frame(
      view = view_name,
      view_label = unname(
        view_labels[view_name]
      ),
      n_samples_measured =
        length(measured_samples),
      n_features_input =
        ncol(variance_reference_matrix),
      n_features_variable =
        length(keep_features),
      n_features_selected =
        ncol(view_matrix),
      selection_mode =
        mofa_feature_selection_mode,
      target_cumulative_variance_fraction = if (
        identical(
          mofa_feature_selection_mode,
          "balanced_cumulative_variance"
        ) &&
        length(keep_features) >
          mofa_keep_all_if_n_features_at_most
      ) {
        mofa_cumulative_variance_target
      } else {
        1
      },
      variance_quantile =
        mofa_hvf_variance_quantile[view_name],
      variance_cutoff = variance_cutoff,
      minimum_features =
        mofa_hvf_min_features[view_name],
      maximum_features =
        mofa_hvf_max_features[view_name],
      achieved_cumulative_variance_fraction =
        max(
          variance_rank_table$cumulative_variance_fraction[
            variance_rank_table$selected
          ],
          na.rm = TRUE
        ),
      transformation =
        transformation_description,
      feature_scaling =
        feature_scaling_description,
      stringsAsFactors = FALSE
    ),
    selected_feature_variance =
      variance_rank_table,
    sensitivity_fraction_table =
      sensitivity_fraction_table,
    transformation_diagnostics =
      transformation_diagnostics,
    clr_replacement_summary =
      clr_replacement_summary
  )
}

mofa_view_preparation <- lapply(
  required_views,
  prepare_mofa_view
)

names(mofa_view_preparation) <- required_views


# v16 feature-count audit. This records where the metabolite view becomes
# 25 features, rather than assuming that the number was created by the
# cumulative-variance rule.
count_matrix_features <- function(x) {
  if (is.null(x)) {
    return(NA_integer_)
  }

  x <- as.matrix(x)

  if (length(dim(x)) != 2) {
    return(NA_integer_)
  }

  ncol(x)
}

mofa_feature_count_audit <- dplyr::bind_rows(
  lapply(
    required_views,
    function(view_name) {
      modality <- coherence_data$modalities[[view_name]]

      data.frame(
        view = view_name,
        n_features_raw =
          count_matrix_features(modality$raw),
        n_features_filtered =
          count_matrix_features(modality$filtered),
        n_features_relative =
          count_matrix_features(modality$relative),
        n_features_transformed =
          count_matrix_features(modality$transformed),
        n_features_table =
          count_matrix_features(modality$table),
        n_features_variable =
          mofa_view_preparation[[view_name]]$summary$n_features_variable,
        n_features_selected =
          mofa_view_preparation[[view_name]]$summary$n_features_selected,
        selection_mode =
          mofa_view_preparation[[view_name]]$summary$selection_mode,
        target_cumulative_variance_fraction =
          mofa_view_preparation[[view_name]]$summary$
            target_cumulative_variance_fraction,
        achieved_cumulative_variance_fraction =
          mofa_view_preparation[[view_name]]$summary$
            achieved_cumulative_variance_fraction,
        stringsAsFactors = FALSE
      )
    }
  )
)

mofa_data <- lapply(
  mofa_view_preparation,
  function(x) x$matrix
)

mofa_feature_selection_summary <- dplyr::bind_rows(
  lapply(
    mofa_view_preparation,
    function(x) x$summary
  )
)

mofa_selected_feature_variance <- dplyr::bind_rows(
  lapply(
    mofa_view_preparation,
    function(x) x$selected_feature_variance
  )
)

mofa_feature_selection_fraction_table <-
  dplyr::bind_rows(
    lapply(
      mofa_view_preparation,
      function(x) x$sensitivity_fraction_table
    )
  )

mofa_metabolite_distribution_diagnostics <-
  dplyr::bind_rows(
    lapply(
      mofa_view_preparation,
      function(x) x$transformation_diagnostics
    )
  )

mofa_clr_replacement_summary <-
  dplyr::bind_rows(
    lapply(
      mofa_view_preparation,
      function(x) x$clr_replacement_summary
    )
  )

if (
  !all(
    vapply(
      mofa_data,
      function(x) {
        identical(
          colnames(x),
          mofa_sample_ids
        )
      },
      logical(1)
    )
  )
) {
  stop(
    "MOFA view matrices do not share the same ordered sample columns.",
    call. = FALSE
  )
}


#=================================================================#
# 4. MOFA2 1.12.x / mofapy2 0.7.0 compatibility runner
#=================================================================#
#
# MOFA2 1.12.1 calls a dotted Python attribute named
# "mofapy2$run.entry_point". With newer reticulate builds this may be
# interpreted as one literal attribute and fail even though the submodule
# mofapy2.run.entry_point is installed correctly. The runner below imports the
# Python submodule explicitly and otherwise reproduces the MOFA2 training call.
# It can be removed after upgrading to a recent MOFA2/Bioconductor release.
#
#-----------------------------------------------------------------#

run_mofa_reticulate_compat <- function(
    object,
    outfile,
    save_data
) {
  if (!methods::is(object, "MOFA")) {
    stop(
      "'object' has to be an instance of MOFA.",
      call. = FALSE
    )
  }

  entry_point_module <- reticulate::import(
    "mofapy2.run.entry_point",
    convert = FALSE
  )

  mofa_entrypoint <- entry_point_module$entry_point()

  mofa_entrypoint$set_data_options(
    scale_views = object@data_options$scale_views,
    scale_groups = object@data_options$scale_groups
  )

  if (
    "samples_metadata" %in%
      methods::slotNames(object)
  ) {
    mofa_entrypoint$data_opts$samples_metadata <-
      reticulate::r_to_py(
        lapply(
          object@data_options$groups,
          function(group_name) {
            object@samples_metadata[
              object@samples_metadata$group ==
                group_name,
              ,
              drop = FALSE
            ]
          }
        )
      )
  }

  if (
    "features_metadata" %in%
      methods::slotNames(object)
  ) {
    mofa_entrypoint$data_opts$features_metadata <-
      reticulate::r_to_py(
        unname(
          lapply(
            object@data_options$views,
            function(view_name) {
              object@features_metadata[
                object@features_metadata$view ==
                  view_name,
                ,
                drop = FALSE
              ]
            }
          )
        )
      )
  }

  mofa_entrypoint$set_data_matrix(
    data = reticulate::r_to_py(
      unname(
        lapply(
          object@data,
          function(view_data) {
            unname(
              lapply(
                view_data,
                function(group_matrix) {
                  reticulate::r_to_py(
                    t(group_matrix)
                  )
                }
              )
            )
          }
        )
      )
    ),
    likelihoods = unname(
      object@model_options$likelihoods
    ),
    views_names = reticulate::r_to_py(
      as.list(
        object@data_options$views
      )
    ),
    groups_names = reticulate::r_to_py(
      as.list(
        object@data_options$groups
      )
    ),
    samples_names = reticulate::r_to_py(
      unname(
        lapply(
          object@data[[1]],
          colnames
        )
      )
    ),
    features_names = reticulate::r_to_py(
      unname(
        lapply(
          object@data,
          function(view_data) {
            rownames(
              view_data[[1]]
            )
          }
        )
      )
    )
  )

  mofa_entrypoint$set_model_options(
    factors = object@model_options$num_factors,
    spikeslab_factors =
      object@model_options$spikeslab_factors,
    spikeslab_weights =
      object@model_options$spikeslab_weights,
    ard_factors =
      object@model_options$ard_factors,
    ard_weights =
      object@model_options$ard_weights
  )

  mofa_entrypoint$set_train_options(
    iter = object@training_options$maxiter,
    convergence_mode =
      object@training_options$convergence_mode,
    dropR2 =
      object@training_options$drop_factor_threshold,
    startELBO = object@training_options$startELBO,
    freqELBO = object@training_options$freqELBO,
    seed = object@training_options$seed,
    gpu_mode = object@training_options$gpu_mode,
    verbose = object@training_options$verbose,
    outfile = object@training_options$outfile,
    save_interrupted =
      object@training_options$save_interrupted
  )

  if (isTRUE(object@training_options$stochastic)) {
    mofa_entrypoint$set_stochastic_options(
      learning_rate =
        object@stochastic_options$learning_rate,
      forgetting_rate =
        object@stochastic_options$forgetting_rate,
      batch_size =
        object@stochastic_options$batch_size,
      start_stochastic =
        object@stochastic_options$start_stochastic
    )
  }

  mofa_entrypoint$build()
  mofa_entrypoint$run()
  mofa_entrypoint$save(
    outfile,
    save_data = save_data
  )

  invisible(NULL)
}


run_mofa_compat <- function(
    object,
    outfile,
    save_data = TRUE,
    use_basilisk = FALSE
) {
  if (!methods::is(object, "MOFA")) {
    stop(
      "'object' has to be an instance of MOFA.",
      call. = FALSE
    )
  }

  if (object@status == "trained") {
    stop(
      "The model is already trained.",
      call. = FALSE
    )
  }

  if (
    length(outfile) != 1 ||
    is.na(outfile) ||
    !nzchar(outfile)
  ) {
    stop(
      "A non-empty .hdf5 outfile must be supplied.",
      call. = FALSE
    )
  }

  outfile <- normalizePath(
    outfile,
    winslash = "/",
    mustWork = FALSE
  )

  if (file.exists(outfile)) {
    file.remove(outfile)
  }

  if (use_basilisk) {
    mofa_basilisk_environment <- get(
      "mofa_env",
      envir = asNamespace("MOFA2")
    )

    basilisk_process <- basilisk::basiliskStart(
      mofa_basilisk_environment
    )

    on.exit(
      basilisk::basiliskStop(
        basilisk_process
      ),
      add = TRUE
    )

    basilisk::basiliskRun(
      basilisk_process,
      function(
          object,
          outfile,
          save_data,
          runner
      ) {
        runner(
          object = object,
          outfile = outfile,
          save_data = save_data
        )
      },
      object = object,
      outfile = outfile,
      save_data = save_data,
      runner = run_mofa_reticulate_compat
    )
  } else {
    run_mofa_reticulate_compat(
      object = object,
      outfile = outfile,
      save_data = save_data
    )
  }

  if (!file.exists(outfile)) {
    stop(
      "MOFA training ended without creating the expected HDF5 file: ",
      outfile,
      call. = FALSE
    )
  }

  MOFA2::load_model(outfile)
}


# Record the R and Python versions used for training.
python_metadata <- reticulate::import(
  "importlib.metadata",
  convert = TRUE
)

python_package_version <- function(package_name) {
  tryCatch(
    as.character(
      python_metadata$version(package_name)
    ),
    error = function(e) NA_character_
  )
}

mofa_software_versions <- data.frame(
  component = c(
    "R",
    "Bioconductor",
    "MOFA2",
    "reticulate",
    "Python",
    "numpy",
    "scipy",
    "pandas",
    "h5py",
    "scikit-learn",
    "anndata",
    "mofapy2"
  ),
  version = c(
    paste(
      R.version$major,
      R.version$minor,
      sep = "."
    ),
    as.character(
      BiocManager::version()
    ),
    as.character(
      utils::packageVersion("MOFA2")
    ),
    as.character(
      utils::packageVersion("reticulate")
    ),
    as.character(
      reticulate::py_config()$version
    ),
    python_package_version("numpy"),
    python_package_version("scipy"),
    python_package_version("pandas"),
    python_package_version("h5py"),
    python_package_version("scikit-learn"),
    python_package_version("anndata"),
    python_package_version("mofapy2")
  ),
  environment = c(
    rep(NA_character_, 4),
    rep(mofa_conda_env, 8)
  ),
  stringsAsFactors = FALSE
)


#=================================================================#
# 5. Train MOFA or validate/reuse the dedicated v30 model
#=================================================================#

# Reuse is permitted only when the HDF5 contains exactly the current views,
# sample identifiers, and selected feature identifiers. A stale or incompatible
# model is ignored and overwritten after a new fit.
mofa_existing_model <- NULL
mofa_model_validation <- data.frame(
  check = c(
    "model_file_exists",
    "model_loadable",
    "views_match",
    "samples_match",
    "features_match"
  ),
  passed = c(
    file.exists(mofa_best_model_path),
    FALSE,
    FALSE,
    FALSE,
    FALSE
  ),
  detail = NA_character_,
  stringsAsFactors = FALSE
)

if (file.exists(mofa_best_model_path)) {
  mofa_existing_model <- tryCatch(
    MOFA2::load_model(mofa_best_model_path),
    error = function(e) {
      mofa_model_validation$detail[
        mofa_model_validation$check == "model_loadable"
      ] <- conditionMessage(e)
      NULL
    }
  )

  if (!is.null(mofa_existing_model)) {
    mofa_model_validation$passed[
      mofa_model_validation$check == "model_loadable"
    ] <- TRUE

    existing_views <- as.character(
      MOFA2::views_names(mofa_existing_model)
    )
    existing_samples_by_group <- MOFA2::samples_names(
      mofa_existing_model
    )
    existing_samples <- as.character(
      unlist(existing_samples_by_group, use.names = FALSE)
    )
    existing_features <- MOFA2::features_names(
      mofa_existing_model
    )

    if (is.null(names(existing_features))) {
      names(existing_features) <- existing_views[
        seq_along(existing_features)
      ]
    }

    views_match <- setequal(
      existing_views,
      required_views
    )
    samples_match <- setequal(
      existing_samples,
      mofa_sample_ids
    )
    features_match <- views_match &&
      all(
        vapply(
          required_views,
          function(view_name) {
            view_name %in% names(existing_features) &&
              setequal(
                as.character(existing_features[[view_name]]),
                rownames(mofa_data[[view_name]])
              )
          },
          logical(1)
        )
      )

    mofa_model_validation$passed[
      mofa_model_validation$check == "views_match"
    ] <- views_match
    mofa_model_validation$passed[
      mofa_model_validation$check == "samples_match"
    ] <- samples_match
    mofa_model_validation$passed[
      mofa_model_validation$check == "features_match"
    ] <- features_match

    mofa_model_validation$detail[
      mofa_model_validation$check == "views_match"
    ] <- paste(existing_views, collapse = ";")
    mofa_model_validation$detail[
      mofa_model_validation$check == "samples_match"
    ] <- paste0(
      "model=", length(existing_samples),
      "; current=", length(mofa_sample_ids)
    )
    mofa_model_validation$detail[
      mofa_model_validation$check == "features_match"
    ] <- paste(
      vapply(
        required_views,
        function(view_name) {
          paste0(
            view_name,
            ": model=",
            if (view_name %in% names(existing_features)) {
              length(existing_features[[view_name]])
            } else {
              0
            },
            ", current=",
            nrow(mofa_data[[view_name]])
          )
        },
        character(1)
      ),
      collapse = "; "
    )

    mofa_reuse_existing_best_model <-
      isTRUE(views_match) &&
      isTRUE(samples_match) &&
      isTRUE(features_match)
  }
}

if (isTRUE(mofa_reuse_existing_best_model)) {
  message(
    "Reusing validated v30 MOFA model: ",
    mofa_best_model_path
  )

  mofa_model <- mofa_existing_model

  if (file.exists(mofa_seed_summary_path)) {
    mofa_seed_summary <- read.csv(
      mofa_seed_summary_path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    mofa_seed_summary <- data.frame(
      seed = NA_integer_,
      final_elbo = NA_real_,
      model_path = mofa_best_model_path,
      selected = TRUE,
      stringsAsFactors = FALSE
    )
  }
} else {
  if (file.exists(mofa_best_model_path)) {
    message(
      "Existing v30 HDF5 is incompatible with the current input and will be replaced."
    )
    file.remove(mofa_best_model_path)
  }

  mofa_seed_summary <- data.frame()

  for (seed_value in mofa_seeds) {
    mofa_object_seed <- MOFA2::create_mofa(
      mofa_data
    )

    mofa_data_options <- MOFA2::get_default_data_options(
      mofa_object_seed
    )
    mofa_data_options$scale_views <- TRUE
    mofa_data_options$scale_groups <- FALSE
    mofa_data_options$center_groups <- TRUE
    mofa_data_options$use_float32 <- TRUE

    mofa_model_options <- MOFA2::get_default_model_options(
      mofa_object_seed
    )
    mofa_model_options$num_factors <- mofa_initial_factors
    mofa_model_options$likelihoods[] <- "gaussian"
    mofa_model_options$ard_weights <- TRUE
    mofa_model_options$ard_factors <- FALSE
    mofa_model_options$spikeslab_weights <-
      mofa_use_spikeslab_weights
    mofa_model_options$spikeslab_factors <- FALSE

    mofa_training_options <- MOFA2::get_default_training_options(
      mofa_object_seed
    )
    mofa_training_options$maxiter <- 5000
    mofa_training_options$convergence_mode <- "slow"
    mofa_training_options$drop_factor_threshold <-
      mofa_drop_factor_threshold
    mofa_training_options$startELBO <- 1
    mofa_training_options$freqELBO <- 5
    mofa_training_options$verbose <- TRUE
    mofa_training_options$seed <- seed_value
    mofa_training_options$stochastic <- FALSE
    mofa_training_options$gpu_mode <- FALSE

    mofa_object_seed <- MOFA2::prepare_mofa(
      object = mofa_object_seed,
      data_options = mofa_data_options,
      model_options = mofa_model_options,
      training_options = mofa_training_options
    )

    seed_model_path <- paste0(
      "results/mofa/MOFA_4omics_v30_seed_",
      seed_value,
      ".hdf5"
    )

    if (file.exists(seed_model_path)) {
      file.remove(seed_model_path)
    }

    set.seed(seed_value)

    mofa_model_seed <- run_mofa_compat(
      object = mofa_object_seed,
      outfile = seed_model_path,
      save_data = mofa_save_training_data,
      use_basilisk = mofa_use_basilisk
    )

    seed_elbo <- as.numeric(
      unlist(
        MOFA2::get_elbo(mofa_model_seed)
      )
    )
    seed_elbo <- seed_elbo[is.finite(seed_elbo)]

    mofa_seed_summary <- rbind(
      mofa_seed_summary,
      data.frame(
        seed = seed_value,
        final_elbo = if (length(seed_elbo) > 0) {
          utils::tail(seed_elbo, 1)
        } else {
          NA_real_
        },
        model_path = seed_model_path,
        stringsAsFactors = FALSE
      )
    )

    rm(
      mofa_object_seed,
      mofa_model_seed
    )
    gc()
  }

  if (!any(is.finite(mofa_seed_summary$final_elbo))) {
    stop(
      "No v30 MOFA initialization returned a finite final ELBO.",
      call. = FALSE
    )
  }

  best_seed_row <- which.max(
    mofa_seed_summary$final_elbo
  )
  mofa_seed_summary$selected <- FALSE
  mofa_seed_summary$selected[best_seed_row] <- TRUE

  write.csv(
    mofa_seed_summary,
    mofa_seed_summary_path,
    row.names = FALSE
  )

  if (!file.copy(
    from = mofa_seed_summary$model_path[best_seed_row],
    to = mofa_best_model_path,
    overwrite = TRUE
  )) {
    stop(
      "The selected v30 seed model could not be copied to the primary HDF5 path.",
      call. = FALSE
    )
  }

  mofa_model <- MOFA2::load_model(
    mofa_best_model_path
  )
  mofa_reuse_existing_best_model <- FALSE
}

rm(mofa_existing_model)

mofa_v30_model_settings_audit <- data.frame(
  setting = c(
    "input_rdata",
    "feature_selection_mode",
    "cumulative_variance_target",
    "species_min_features",
    "species_max_features",
    "ko_min_features",
    "ko_max_features",
    "metabolite_min_features",
    "metabolite_max_features",
    "host_min_features",
    "host_max_features",
    "initial_factors",
    "scale_views",
    "ard_weights",
    "spikeslab_weights",
    "maxiter",
    "convergence_mode",
    "active_view_r2",
    "strict_active_view_r2",
    "model_hdf5",
    "reused_validated_model"
  ),
  value = c(
    mofa_input_rdata,
    mofa_feature_selection_mode,
    as.character(mofa_cumulative_variance_target),
    as.character(mofa_hvf_min_features["species"]),
    as.character(mofa_hvf_max_features["species"]),
    as.character(mofa_hvf_min_features["ko"]),
    as.character(mofa_hvf_max_features["ko"]),
    as.character(mofa_hvf_min_features["metabolite"]),
    as.character(mofa_hvf_max_features["metabolite"]),
    as.character(mofa_hvf_min_features["host"]),
    as.character(mofa_hvf_max_features["host"]),
    as.character(mofa_initial_factors),
    "TRUE",
    "TRUE",
    as.character(mofa_use_spikeslab_weights),
    "5000",
    "slow",
    as.character(mofa_active_view_r2),
    as.character(mofa_strict_active_view_r2),
    mofa_best_model_path,
    as.character(mofa_reuse_existing_best_model)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  mofa_v30_model_settings_audit,
  "results/mofa/MOFA_v30_model_settings_audit.csv",
  row.names = FALSE
)

# samples_metadata() requires the metadata group labels to match the group
# names stored inside the trained HDF5 model. Do not assume the default group
# name, because matrix-list input may use "group1" rather than "single_group".
model_group_names <- MOFA2::groups_names(
  mofa_model
)

model_sample_names_by_group <- MOFA2::samples_names(
  mofa_model
)

if (
  length(model_group_names) !=
    length(model_sample_names_by_group)
) {
  stop(
    "The number of MOFA group names does not match the sample-name groups.",
    call. = FALSE
  )
}

model_sample_index <- do.call(
  rbind,
  lapply(
    seq_along(model_sample_names_by_group),
    function(group_index) {
      data.frame(
        sample = as.character(
          model_sample_names_by_group[[group_index]]
        ),
        group = rep(
          as.character(
            model_group_names[group_index]
          ),
          length(
            model_sample_names_by_group[[group_index]]
          )
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(model_sample_index) <- NULL
model_sample_names <- model_sample_index$sample

if (anyDuplicated(model_sample_names) > 0) {
  stop(
    "Duplicated sample names were found across MOFA groups.",
    call. = FALSE
  )
}

mofa_model_metadata <- mofa_sample_metadata[
  match(
    model_sample_names,
    mofa_sample_metadata$SampleID
  ),
  ,
  drop = FALSE
]

if (anyNA(mofa_model_metadata$SampleID)) {
  stop(
    "At least one sample stored in the MOFA model is absent from the metadata.",
    call. = FALSE
  )
}

mofa_model_metadata$sample <- model_sample_index$sample
mofa_model_metadata$group <- model_sample_index$group

mofa_model_metadata <- mofa_model_metadata[
  ,
  c(
    "sample",
    "group",
    setdiff(
      colnames(mofa_model_metadata),
      c(
        "sample",
        "group"
      )
    )
  ),
  drop = FALSE
]

samples_metadata(mofa_model) <- mofa_model_metadata


#=================================================================#
# 6. Extract factor scores, variance explained, and feature weights
#=================================================================#

mofa_factor_list <- MOFA2::get_factors(
  mofa_model,
  groups = "all",
  factors = "all",
  scale = FALSE,
  as.data.frame = FALSE
)

if (length(mofa_factor_list) != length(model_group_names)) {
  stop(
    "The number of extracted factor-score matrices does not match the model groups.",
    call. = FALSE
  )
}

mofa_factor_matrix <- do.call(
  rbind,
  lapply(
    seq_along(mofa_factor_list),
    function(group_index) {
      factor_matrix <- as.matrix(
        mofa_factor_list[[group_index]]
      )

      expected_samples <-
        model_sample_names_by_group[[group_index]]

      if (
        !all(
          expected_samples %in%
            rownames(factor_matrix)
        ) &&
        all(
          expected_samples %in%
            colnames(factor_matrix)
        )
      ) {
        factor_matrix <- base::t(
          factor_matrix
        )
      }

      factor_matrix[
        expected_samples,
        ,
        drop = FALSE
      ]
    }
  )
)

mofa_factor_matrix <- mofa_factor_matrix[
  model_sample_names,
  ,
  drop = FALSE
]

mofa_factor_scores <- data.frame(
  SampleID = rownames(mofa_factor_matrix),
  mofa_factor_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    mofa_sample_metadata,
    by = "SampleID"
  )

# Reproduce the trained-model quality-control diagnostic. MOFA2 issues the
# warning when a factor has |correlation| > 0.75 with the per-sample sum of
# values in at least one view. This table identifies the responsible view(s).
mofa_factor_total_level_correlations <- do.call(
  rbind,
  lapply(
    required_views,
    function(view_name) {
      view_total_level <- colSums(
        mofa_data[[view_name]],
        na.rm = TRUE
      )

      data.frame(
        view = view_name,
        view_label = unname(
          view_labels[view_name]
        ),
        factor = colnames(mofa_factor_matrix),
        correlation = vapply(
          colnames(mofa_factor_matrix),
          function(factor_name) {
            common_samples <- intersect(
              names(view_total_level),
              rownames(mofa_factor_matrix)
            )

            keep <- is.finite(
              view_total_level[common_samples]
            ) &
              is.finite(
                mofa_factor_matrix[
                  common_samples,
                  factor_name
                ]
              )

            if (sum(keep) < 4) {
              return(NA_real_)
            }

            stats::cor(
              view_total_level[common_samples][keep],
              mofa_factor_matrix[
                common_samples,
                factor_name
              ][keep],
              method = "pearson"
            )
          },
          numeric(1)
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

mofa_factor_total_level_correlations$flag_abs_r_gt_0_75 <-
  abs(
    mofa_factor_total_level_correlations$correlation
  ) > 0.75

write.csv(
  mofa_factor_total_level_correlations,
  paste0(
    "results/mofa/",
    "MOFA_v30_factor_total_level_correlations.csv"
  ),
  row.names = FALSE
)

mofa_variance_explained_raw <- MOFA2::get_variance_explained(
  mofa_model,
  groups = "all",
  views = "all",
  factors = "all",
  as.data.frame = TRUE
)

# MOFA2 1.12.x returns a named list when as.data.frame = TRUE:
#   $r2_per_factor: group, view, factor, value
#   $r2_total:      group, view, value
# Newer or development versions may return a data.frame directly. Standardise
# both structures here and use proportions internally (0-1), while retaining
# percentage columns for transparent export.
if (
  is.list(mofa_variance_explained_raw) &&
  "r2_per_factor" %in% names(mofa_variance_explained_raw)
) {
  mofa_variance_explained <- as.data.frame(
    mofa_variance_explained_raw$r2_per_factor,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  mofa_variance_total <- if (
    "r2_total" %in% names(mofa_variance_explained_raw)
  ) {
    as.data.frame(
      mofa_variance_explained_raw$r2_total,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    NULL
  }
} else if (is.data.frame(mofa_variance_explained_raw)) {
  mofa_variance_explained <- as.data.frame(
    mofa_variance_explained_raw,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  mofa_variance_total <- NULL
} else {
  stop(
    "Unsupported MOFA variance-explained return structure: ",
    paste(
      class(mofa_variance_explained_raw),
      collapse = ", "
    ),
    call. = FALSE
  )
}

r2_value_column <- intersect(
  c(
    "r2",
    "R2",
    "value",
    "variance_explained"
  ),
  colnames(mofa_variance_explained)
)

if (
  !all(
    c(
      "view",
      "factor"
    ) %in% colnames(mofa_variance_explained)
  ) ||
  length(r2_value_column) == 0
) {
  stop(
    paste0(
      "Unexpected columns in MOFA factor-wise variance output: ",
      paste(
        colnames(mofa_variance_explained),
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

mofa_variance_explained$r2_raw <- as.numeric(
  mofa_variance_explained[[
    r2_value_column[1]
  ]]
)

finite_r2_values <- mofa_variance_explained$r2_raw[
  is.finite(
    mofa_variance_explained$r2_raw
  )
]

mofa_r2_returned_as_percent <-
  length(finite_r2_values) > 0 &&
  max(
    finite_r2_values,
    na.rm = TRUE
  ) > 1

mofa_variance_explained$r2 <- if (
  mofa_r2_returned_as_percent
) {
  mofa_variance_explained$r2_raw / 100
} else {
  mofa_variance_explained$r2_raw
}

mofa_variance_explained$r2_percent <-
  100 * mofa_variance_explained$r2

mofa_variance_explained$view <- as.character(
  mofa_variance_explained$view
)

mofa_variance_explained$factor <- as.character(
  mofa_variance_explained$factor
)

if (!"group" %in% colnames(mofa_variance_explained)) {
  mofa_variance_explained$group <- model_group_names[1]
}

if (is.null(mofa_variance_total)) {
  mofa_variance_total <- mofa_variance_explained %>%
    dplyr::group_by(
      group,
      view
    ) %>%
    dplyr::summarise(
      total_r2 = sum(
        r2,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      total_r2_raw = total_r2,
      total_r2_percent = 100 * total_r2,
      total_r2_source = "sum_of_factorwise_r2"
    )
} else {
  total_r2_value_column <- intersect(
    c(
      "r2",
      "R2",
      "value",
      "variance_explained"
    ),
    colnames(mofa_variance_total)
  )

  if (
    !"view" %in% colnames(mofa_variance_total) ||
    length(total_r2_value_column) == 0
  ) {
    stop(
      paste0(
        "Unexpected columns in MOFA total variance output: ",
        paste(
          colnames(mofa_variance_total),
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }

  mofa_variance_total$total_r2_raw <- as.numeric(
    mofa_variance_total[[
      total_r2_value_column[1]
    ]]
  )

  mofa_variance_total$total_r2 <- if (
    mofa_r2_returned_as_percent
  ) {
    mofa_variance_total$total_r2_raw / 100
  } else {
    mofa_variance_total$total_r2_raw
  }

  mofa_variance_total$total_r2_percent <-
    100 * mofa_variance_total$total_r2

  mofa_variance_total$total_r2_source <-
    "MOFA_get_variance_explained_r2_total"

  mofa_variance_total$view <- as.character(
    mofa_variance_total$view
  )

  if (!"group" %in% colnames(mofa_variance_total)) {
    mofa_variance_total$group <- model_group_names[1]
  }
}

mofa_variance_explained$view_label <- factor(
  unname(
    view_labels[
      as.character(
        mofa_variance_explained$view
      )
    ]
  ),
  levels = unname(view_labels)
)


mofa_variance_total$view_label <- factor(
  unname(
    view_labels[
      as.character(
        mofa_variance_total$view
      )
    ]
  ),
  levels = unname(view_labels)
)

mofa_factor_summary <- mofa_variance_explained %>%
  dplyr::group_by(factor) %>%
  dplyr::summarise(
    total_r2 = sum(
      r2,
      na.rm = TRUE
    ),
    mean_r2 = mean(
      r2,
      na.rm = TRUE
    ),
    active_views = sum(
      r2 >= mofa_active_view_r2,
      na.rm = TRUE
    ),
    strongest_view = as.character(
      view[
        which.max(r2)
      ]
    ),
    strongest_view_r2 = max(
      r2,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    factor_type = dplyr::case_when(
      active_views >= 3 ~ "Broadly shared",
      active_views == 2 ~ "Pair-shared",
      active_views == 1 ~ "View-specific",
      TRUE ~ "Low-activity"
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(total_r2)
  )

mofa_factor_order <- as.character(
  mofa_factor_summary$factor
)

mofa_variance_explained$factor <- factor(
  as.character(
    mofa_variance_explained$factor
  ),
  levels = rev(mofa_factor_order)
)

mofa_weights_list <- MOFA2::get_weights(
  mofa_model,
  views = "all",
  factors = "all",
  abs = FALSE,
  scale = FALSE,
  as.data.frame = FALSE
)

mofa_feature_weights <- dplyr::bind_rows(
  lapply(
    names(mofa_weights_list),
    function(view_name) {
      weight_matrix <- as.matrix(
        mofa_weights_list[[view_name]]
      )

      expected_features <- rownames(
        mofa_data[[view_name]]
      )

      if (
        !all(
          expected_features %in%
            rownames(weight_matrix)
        ) &&
        all(
          expected_features %in%
            colnames(weight_matrix)
        )
      ) {
        weight_matrix <- base::t(
          weight_matrix
        )
      }

      weight_matrix <- weight_matrix[
        expected_features,
        ,
        drop = FALSE
      ]

      data.frame(
        view = view_name,
        feature = rep(
          rownames(weight_matrix),
          times = ncol(weight_matrix)
        ),
        factor = rep(
          colnames(weight_matrix),
          each = nrow(weight_matrix)
        ),
        weight = as.numeric(weight_matrix),
        stringsAsFactors = FALSE
      )
    }
  )
) %>%
  dplyr::mutate(
    view_label = unname(
      view_labels[view]
    ),
    abs_weight = abs(weight),
    direction = ifelse(
      weight >= 0,
      "Positive",
      "Negative"
    )
  ) %>%
  dplyr::group_by(
    view,
    factor
  ) %>%
  dplyr::arrange(
    dplyr::desc(abs_weight),
    .by_group = TRUE
  ) %>%
  dplyr::mutate(
    absolute_rank = dplyr::row_number()
  ) %>%
  dplyr::ungroup()

mofa_feature_weights$feature_label <- mofa_feature_weights$feature

if (
  !is.null(coherence_data$feature_annotation$ko) &&
  all(
    c(
      "KO_number",
      "Protein_name"
    ) %in%
      colnames(
        coherence_data$feature_annotation$ko
      )
  )
) {
  ko_label <- coherence_data$feature_annotation$ko %>%
    dplyr::transmute(
      feature = as.character(KO_number),
      annotated_label = ifelse(
        is.na(Protein_name) |
          Protein_name == "",
        feature,
        paste0(
          feature,
          ": ",
          Protein_name
        )
      )
    )

  mofa_feature_weights <- mofa_feature_weights %>%
    dplyr::left_join(
      ko_label,
      by = "feature"
    ) %>%
    dplyr::mutate(
      feature_label = ifelse(
        view == "ko" &
          !is.na(annotated_label),
        annotated_label,
        feature_label
      )
    ) %>%
    dplyr::select(
      -annotated_label
    )
}

mofa_top_features <- mofa_feature_weights %>%
  dplyr::group_by(
    view,
    factor,
    direction
  ) %>%
  dplyr::slice_max(
    order_by = abs_weight,
    n = mofa_top_features_per_direction,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

write.csv(
  mofa_factor_scores,
  "results/mofa/MOFA_v30_factor_scores.csv",
  row.names = FALSE
)

write.csv(
  mofa_variance_explained,
  "results/mofa/MOFA_v30_variance_explained.csv",
  row.names = FALSE
)


write.csv(
  mofa_variance_total,
  "results/mofa/MOFA_v30_total_variance_explained.csv",
  row.names = FALSE
)


#=================================================================#
# 7. Feature-selection diagnostics and factor stability
#=================================================================#

mofa_selected_feature_variance$view_label <- factor(
  unname(
    view_labels[
      mofa_selected_feature_variance$view
    ]
  ),
  levels = unname(view_labels)
)

mofa_feature_selection_fraction_table <- dplyr::bind_rows(
  lapply(
    split(
      mofa_selected_feature_variance,
      mofa_selected_feature_variance$view
    ),
    function(view_variance_table) {
      dplyr::bind_rows(
        lapply(
          mofa_feature_sensitivity_fractions,
          function(target_fraction) {
            selected_rank <- if (target_fraction >= 1) {
              nrow(view_variance_table)
            } else {
              which(
                view_variance_table$cumulative_variance_fraction >=
                  target_fraction
              )[1]
            }

            data.frame(
              view = view_variance_table$view[1],
              target_cumulative_variance_fraction = target_fraction,
              n_features_required = selected_rank,
              achieved_cumulative_variance_fraction =
                view_variance_table$cumulative_variance_fraction[selected_rank],
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
) %>%
  dplyr::mutate(
    view_label = factor(
      unname(view_labels[view]),
      levels = unname(view_labels)
    ),
    target_label = scales::percent(
      target_cumulative_variance_fraction,
      accuracy = 1
    )
  )

mofa_feature_endpoint_data <-
  mofa_selected_feature_variance %>%
  dplyr::group_by(view, view_label) %>%
  dplyr::slice_max(
    variance_rank,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

p_mofa_feature_selection <- ggplot2::ggplot(
  mofa_selected_feature_variance,
  ggplot2::aes(
    x = variance_rank,
    y = cumulative_variance_fraction,
    color = view_label
  )
) +
  ggplot2::geom_hline(
    yintercept = c(0.50, 0.75, 0.90, 0.95, 0.99),
    linetype = "dotted",
    linewidth = 0.25,
    color = "grey78"
  ) +
  ggplot2::geom_line(
    linewidth = 0.78
  ) +
  ggplot2::geom_point(
    data = mofa_feature_selection_fraction_table %>%
      dplyr::filter(
        target_cumulative_variance_fraction %in%
          c(0.50, 0.75, 0.90, 0.95, 0.99)
      ),
    ggplot2::aes(
      x = n_features_required,
      y = achieved_cumulative_variance_fraction
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 1.75,
    fill = "white",
    color = "grey25",
    stroke = 0.45
  ) +
  ggplot2::geom_point(
    data = mofa_feature_endpoint_data,
    size = 2.35,
    color = "grey15"
  ) +
  ggplot2::geom_text(
    data = mofa_feature_endpoint_data,
    ggplot2::aes(
      label = scales::comma(variance_rank)
    ),
    hjust = 1.08,
    vjust = -0.55,
    size = 2.55,
    color = "grey20",
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(
    ~ view_label,
    scales = "free_x",
    nrow = 1
  ) +
  ggplot2::scale_color_manual(
    values = view_colors,
    guide = "none"
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, 1.02),
    breaks = c(0, 0.50, 0.75, 0.90, 1.00),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::scale_x_continuous(
    labels = scales::comma
  ) +
  ggplot2::labs(
    x = "Feature rank by variance before MOFA scaling",
    y = "Cumulative across-feature variance"
  ) +
  ggplot2::theme_classic(base_size = 9.8) +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "plain", size = 9.1),
    axis.title = ggplot2::element_text(size = 9.0),
    axis.text = ggplot2::element_text(size = 8.0),
    plot.margin = ggplot2::margin(4, 6, 4, 4)
  )

mofa_factor_correlation <- stats::cor(
  mofa_factor_matrix,
  use = "pairwise.complete.obs"
)

mofa_factor_correlation_plot_data <- as.data.frame(
  as.table(
    mofa_factor_correlation
  ),
  stringsAsFactors = FALSE
)

colnames(
  mofa_factor_correlation_plot_data
) <- c(
  "factor_x",
  "factor_y",
  "correlation"
)

p_mofa_factor_correlation <- ggplot2::ggplot(
  mofa_factor_correlation_plot_data,
  ggplot2::aes(
    x = factor_x,
    y = factor_y,
    fill = correlation
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.30
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf(
        "%.2f",
        correlation
      )
    ),
    size = 2.6
  ) +
  ggplot2::scale_fill_gradient2(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "#3B75AF",
    mid = "white",
    high = "#C44E52",
    midpoint = 0,
    limits = c(
      -1,
      1
    ),
    name = "Pearson r"
  ) +
  ggplot2::coord_fixed() +
  ggplot2::labs(
    title = "Correlation among MOFA factor scores",
    subtitle = "High off-diagonal correlations indicate partially redundant sample gradients.",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 10.0
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 8.6,
      color = "grey35"
    ),
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      face = "plain"
    ),
    axis.text.y = ggplot2::element_text(
      face = "plain"
    ),
    axis.ticks = ggplot2::element_blank()
  )

mofa_seed_factor_stability <- data.frame()

if (
  nrow(mofa_seed_summary) > 1 &&
  any(
    file.exists(
      mofa_seed_summary$model_path
    )
  )
) {
  for (seed_index in seq_len(nrow(mofa_seed_summary))) {
    if (
      !file.exists(
        mofa_seed_summary$model_path[seed_index]
      ) ||
      isTRUE(
        mofa_seed_summary$selected[seed_index]
      )
    ) {
      next
    }

    seed_model_for_stability <- tryCatch(
      MOFA2::load_model(
        mofa_seed_summary$model_path[seed_index]
      ),
      error = function(e) {
        NULL
      }
    )

    if (is.null(seed_model_for_stability)) {
      next
    }

    seed_factor_list <- MOFA2::get_factors(
      seed_model_for_stability,
      groups = "all",
      factors = "all",
      scale = FALSE,
      as.data.frame = FALSE
    )

    seed_factor_matrix <- as.matrix(
      seed_factor_list[[1]]
    )

    if (
      !all(
        model_sample_names %in%
          rownames(seed_factor_matrix)
      ) &&
      all(
        model_sample_names %in%
          colnames(seed_factor_matrix)
      )
    ) {
      seed_factor_matrix <- base::t(
        seed_factor_matrix
      )
    }

    common_seed_samples <- intersect(
      model_sample_names,
      rownames(seed_factor_matrix)
    )

    if (length(common_seed_samples) < 10) {
      next
    }

    factor_correlation_matrix <- abs(
      stats::cor(
        mofa_factor_matrix[
          common_seed_samples,
          ,
          drop = FALSE
        ],
        seed_factor_matrix[
          common_seed_samples,
          ,
          drop = FALSE
        ],
        use = "pairwise.complete.obs"
      )
    )

    mofa_seed_factor_stability <- dplyr::bind_rows(
      mofa_seed_factor_stability,
      data.frame(
        reference_factor = rownames(
          factor_correlation_matrix
        ),
        comparison_seed =
          mofa_seed_summary$seed[seed_index],
        best_matching_factor = apply(
          factor_correlation_matrix,
          1,
          function(x) {
            colnames(
              factor_correlation_matrix
            )[which.max(x)]
          }
        ),
        absolute_score_correlation = apply(
          factor_correlation_matrix,
          1,
          max,
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    )
  }
}


#=================================================================#
# 8. Factor-score association with Timepoint and treatment response
#=================================================================#

mofa_factor_scores <- mofa_factor_scores %>%
  dplyr::left_join(
    mofa_availability_wide %>%
      dplyr::select(
        SampleID,
        dplyr::all_of(required_views),
        n_views
      ),
    by = "SampleID"
  )

mofa_factor_scores$Timepoint <- stats::relevel(
  factor(
    mofa_factor_scores$Timepoint,
    levels = c(
      "Before",
      "Ongoing"
    )
  ),
  ref = "Before"
)

mofa_factor_scores$TRG_plot <- stats::relevel(
  factor(
    mofa_factor_scores$TRG_plot,
    levels = c(
      "non_pCR",
      "pCR"
    )
  ),
  ref = "non_pCR"
)

extract_factor_model_terms <- function(
    model_fit,
    factor_name,
    analysis_set,
    n_samples,
    n_subjects
) {
  coefficient_table <- as.data.frame(
    summary(model_fit)$tTable,
    check.names = FALSE
  )

  coefficient_table$term <- rownames(
    coefficient_table
  )

  rownames(coefficient_table) <- NULL

  coefficient_table %>%
    dplyr::filter(
      term %in% c(
        "TimepointOngoing",
        "TRG_plotpCR",
        "TimepointOngoing:TRG_plotpCR"
      )
    ) %>%
    dplyr::transmute(
      factor = factor_name,
      analysis_set = analysis_set,
      term = term,
      estimate = Value,
      std_error = `Std.Error`,
      degrees_freedom = DF,
      statistic = `t-value`,
      model_p = `p-value`,
      ci_low = estimate -
        stats::qt(
          0.975,
          degrees_freedom
        ) * std_error,
      ci_high = estimate +
        stats::qt(
          0.975,
          degrees_freedom
        ) * std_error,
      n_samples = n_samples,
      n_subjects = n_subjects
    )
}

fit_factor_model <- function(
    factor_name,
    analysis_set = "primary_all_samples"
) {
  factor_model_data <- mofa_factor_scores %>%
    dplyr::transmute(
      SampleID,
      SubjectID,
      Timepoint,
      TRG_plot,
      n_views,
      species,
      ko,
      metabolite,
      host,
      factor_score = .data[[factor_name]]
    ) %>%
    dplyr::filter(
      is.finite(factor_score),
      !is.na(SubjectID),
      !is.na(Timepoint),
      !is.na(TRG_plot)
    )

  if (analysis_set == "strongest_view_measured") {
    strongest_view_current <- mofa_factor_summary$strongest_view[
      match(
        factor_name,
        mofa_factor_summary$factor
      )
    ]

    factor_model_data <- factor_model_data %>%
      dplyr::filter(
        .data[[strongest_view_current]] == 1
      )
  }

  if (analysis_set == "four_view_complete") {
    factor_model_data <- factor_model_data %>%
      dplyr::filter(
        n_views == length(required_views)
      )
  }

  if (
    dplyr::n_distinct(
      factor_model_data$SubjectID
    ) < 8 ||
    min(
      table(
        factor_model_data %>%
          dplyr::distinct(
            SubjectID,
            TRG_plot
          ) %>%
          dplyr::pull(TRG_plot)
      )
    ) < 3
  ) {
    return(
      list(
        fit = NULL,
        data = factor_model_data,
        coefficients = data.frame()
      )
    )
  }

  factor_model_data$n_views_centered <-
    factor_model_data$n_views -
    mean(
      factor_model_data$n_views,
      na.rm = TRUE
    )

  model_formula <- if (
    analysis_set == "primary_plus_n_views"
  ) {
    factor_score ~
      Timepoint * TRG_plot +
      n_views_centered
  } else {
    factor_score ~
      Timepoint * TRG_plot
  }

  model_fit <- tryCatch(
    nlme::lme(
      fixed = model_formula,
      random = ~ 1 | SubjectID,
      data = factor_model_data,
      method = "REML",
      na.action = stats::na.omit,
      control = nlme::lmeControl(
        returnObject = TRUE,
        maxIter = 100,
        msMaxIter = 100
      )
    ),
    error = function(e) {
      NULL
    }
  )

  if (is.null(model_fit)) {
    return(
      list(
        fit = NULL,
        data = factor_model_data,
        coefficients = data.frame()
      )
    )
  }

  list(
    fit = model_fit,
    data = factor_model_data,
    coefficients = extract_factor_model_terms(
      model_fit,
      factor_name,
      analysis_set,
      nrow(factor_model_data),
      dplyr::n_distinct(
        factor_model_data$SubjectID
      )
    )
  )
}

mofa_factor_model_results <- list()
mofa_factor_clinical_associations <- data.frame()

for (factor_name in colnames(mofa_factor_matrix)) {
  for (analysis_set in c(
    "primary_all_samples",
    "primary_plus_n_views",
    "strongest_view_measured",
    "four_view_complete"
  )) {
    factor_model_result <- fit_factor_model(
      factor_name,
      analysis_set
    )

    mofa_factor_model_results[[
      paste(
        factor_name,
        analysis_set,
        sep = "__"
      )
    ]] <- factor_model_result

    mofa_factor_clinical_associations <- dplyr::bind_rows(
      mofa_factor_clinical_associations,
      factor_model_result$coefficients
    )
  }
}

set.seed(20260809)

mofa_factor_permutation_results <- data.frame()

for (factor_name in colnames(mofa_factor_matrix)) {
  primary_result <- mofa_factor_model_results[[
    paste(
      factor_name,
      "primary_all_samples",
      sep = "__"
    )
  ]]

  if (
    is.null(primary_result$fit) ||
    nrow(primary_result$coefficients) == 0
  ) {
    next
  }

  observed_terms <- primary_result$coefficients$term
  observed_estimates <- stats::setNames(
    primary_result$coefficients$estimate,
    observed_terms
  )

  subject_response <- primary_result$data %>%
    dplyr::distinct(
      SubjectID,
      TRG_plot
    ) %>%
    dplyr::arrange(SubjectID)

  permutation_estimates <- matrix(
    NA_real_,
    nrow = mofa_factor_permutations,
    ncol = length(observed_terms),
    dimnames = list(
      NULL,
      observed_terms
    )
  )

  for (permutation_index in seq_len(mofa_factor_permutations)) {
    permuted_response <- subject_response
    permuted_response$TRG_plot <- sample(
      subject_response$TRG_plot,
      size = nrow(subject_response),
      replace = FALSE
    )

    permutation_data <- primary_result$data %>%
      dplyr::select(
        -TRG_plot
      ) %>%
      dplyr::left_join(
        permuted_response,
        by = "SubjectID"
      )

    permutation_data$TRG_plot <- stats::relevel(
      factor(
        permutation_data$TRG_plot,
        levels = c(
          "non_pCR",
          "pCR"
        )
      ),
      ref = "non_pCR"
    )

    permutation_fit <- tryCatch(
      nlme::lme(
        fixed = factor_score ~
          Timepoint * TRG_plot,
        random = ~ 1 | SubjectID,
        data = permutation_data,
        method = "REML",
        na.action = stats::na.omit,
        control = nlme::lmeControl(
          returnObject = TRUE,
          maxIter = 60,
          msMaxIter = 60
        )
      ),
      error = function(e) {
        NULL
      }
    )

    if (is.null(permutation_fit)) {
      next
    }

    permutation_table <- as.data.frame(
      summary(permutation_fit)$tTable,
      check.names = FALSE
    )

    available_terms <- intersect(
      observed_terms,
      rownames(permutation_table)
    )

    permutation_estimates[
      permutation_index,
      available_terms
    ] <- permutation_table[
      available_terms,
      "Value"
    ]
  }

  for (term_name in observed_terms) {
    valid_permutations <- permutation_estimates[
      is.finite(
        permutation_estimates[, term_name]
      ),
      term_name
    ]

    mofa_factor_permutation_results <- dplyr::bind_rows(
      mofa_factor_permutation_results,
      data.frame(
        factor = factor_name,
        term = term_name,
        observed_estimate = observed_estimates[term_name],
        permutation_p = if (
          length(valid_permutations) > 0
        ) {
          (
            1 +
              sum(
                abs(valid_permutations) >=
                  abs(observed_estimates[term_name])
              )
          ) /
          (
            1 +
              length(valid_permutations)
          )
        } else {
          NA_real_
        },
        valid_permutations = length(valid_permutations),
        stringsAsFactors = FALSE
      )
    )
  }
}

mofa_factor_clinical_associations <-
  mofa_factor_clinical_associations %>%
  dplyr::left_join(
    mofa_factor_permutation_results,
    by = c(
      "factor",
      "term"
    )
  ) %>%
  dplyr::mutate(
    permutation_p = ifelse(
      analysis_set == "primary_all_samples",
      permutation_p,
      NA_real_
    ),
    valid_permutations = ifelse(
      analysis_set == "primary_all_samples",
      valid_permutations,
      NA_integer_
    )
  ) %>%
  dplyr::group_by(
    analysis_set,
    term
  ) %>%
  dplyr::mutate(
    model_fdr = stats::p.adjust(
      model_p,
      method = "BH"
    ),
    permutation_fdr = ifelse(
      analysis_set == "primary_all_samples",
      stats::p.adjust(
        permutation_p,
        method = "BH"
      ),
      NA_real_
    )
  ) %>%
  dplyr::ungroup()


#=================================================================#
# 9. Baseline pCR versus non-pCR feature-level statistics
#=================================================================#

calculate_baseline_response_stats <- function(view_name) {
  selected_features <- rownames(
    mofa_data[[view_name]]
  )

  view_meta <- coherence_data$modalities[[view_name]]$sample_meta %>%
    dplyr::filter(
      Timepoint == "Before",
      TRG_plot %in% c(
        "pCR",
        "non_pCR"
      ),
      SampleID %in% colnames(
        mofa_data[[view_name]]
      )
    ) %>%
    dplyr::distinct(
      SampleID,
      .keep_all = TRUE
    )

  if (
    sum(view_meta$TRG_plot == "pCR") < 2 ||
    sum(view_meta$TRG_plot == "non_pCR") < 2
  ) {
    return(
      data.frame(
        view = view_name,
        feature = selected_features,
        response_effect = NA_real_,
        effect_metric = NA_character_,
        p_value = NA_real_,
        fdr = NA_real_,
        n_pCR = sum(
          view_meta$TRG_plot == "pCR"
        ),
        n_non_pCR = sum(
          view_meta$TRG_plot == "non_pCR"
        ),
        stringsAsFactors = FALSE
      )
    )
  }

  model_scale_matrix <- base::t(
    mofa_data[[view_name]][
      selected_features,
      view_meta$SampleID,
      drop = FALSE
    ]
  )

  if (view_name %in% c("species", "ko")) {
    response_scale_matrix <-
      coherence_data$modalities[[view_name]]$relative[
        view_meta$SampleID,
        selected_features,
        drop = FALSE
      ]
    effect_metric <- "log2 ratio of mean relative abundance"
  } else if (view_name == "metabolite") {
    response_scale_matrix <-
      coherence_data$modalities[[view_name]]$raw[
        view_meta$SampleID,
        selected_features,
        drop = FALSE
      ]
    effect_metric <- "log2 ratio of mean concentration"
  } else {
    response_scale_matrix <-
      coherence_data$modalities[[view_name]]$transformed[
        view_meta$SampleID,
        selected_features,
        drop = FALSE
      ]
    effect_metric <- "mean VST difference"
  }

  response_scale_matrix <- as.matrix(
    response_scale_matrix
  )

  response_effect <- vapply(
    selected_features,
    function(feature_name) {
      feature_values <- response_scale_matrix[
        ,
        feature_name
      ]

      if (view_name == "host") {
        return(
          mean(
            feature_values[
              view_meta$TRG_plot == "pCR"
            ],
            na.rm = TRUE
          ) -
          mean(
            feature_values[
              view_meta$TRG_plot == "non_pCR"
            ],
            na.rm = TRUE
          )
        )
      }

      positive_values <- feature_values[
        is.finite(feature_values) &
          feature_values > 0
      ]

      pseudocount <- if (
        length(positive_values) > 0
      ) {
        min(
          positive_values,
          na.rm = TRUE
        ) / 2
      } else {
        1e-08
      }

      log2(
        (
          mean(
            feature_values[
              view_meta$TRG_plot == "pCR"
            ],
            na.rm = TRUE
          ) + pseudocount
        ) /
        (
          mean(
            feature_values[
              view_meta$TRG_plot == "non_pCR"
            ],
            na.rm = TRUE
          ) + pseudocount
        )
      )
    },
    numeric(1)
  )

  p_value <- vapply(
    selected_features,
    function(feature_name) {
      feature_values <- model_scale_matrix[
        ,
        feature_name
      ]

      suppressWarnings(
        tryCatch(
          stats::wilcox.test(
            feature_values[
              view_meta$TRG_plot == "pCR"
            ],
            feature_values[
              view_meta$TRG_plot == "non_pCR"
            ],
            exact = FALSE
          )$p.value,
          error = function(e) {
            NA_real_
          }
        )
      )
    },
    numeric(1)
  )

  data.frame(
    view = view_name,
    feature = selected_features,
    response_effect = response_effect,
    effect_metric = effect_metric,
    p_value = p_value,
    fdr = stats::p.adjust(
      p_value,
      method = "BH"
    ),
    n_pCR = sum(
      view_meta$TRG_plot == "pCR"
    ),
    n_non_pCR = sum(
      view_meta$TRG_plot == "non_pCR"
    ),
    stringsAsFactors = FALSE
  )
}

mofa_feature_response_stats <- dplyr::bind_rows(
  lapply(
    required_views,
    calculate_baseline_response_stats
  )
)

host_response_environment <- new.env(
  parent = emptyenv()
)

host_response_source <- NA_character_

for (host_response_file in mofa_host_response_rdata_candidates) {
  if (!file.exists(host_response_file)) {
    next
  }

  load(
    host_response_file,
    envir = host_response_environment
  )

  if (
    exists(
      "deg_rna",
      envir = host_response_environment,
      inherits = FALSE
    )
  ) {
    host_response_source <- host_response_file
    break
  }
}

host_log2fc_orientation <- data.frame(
  source = host_response_source,
  orientation_correlation = NA_real_,
  log2fc_flipped = FALSE,
  stringsAsFactors = FALSE
)

if (!is.na(host_response_source)) {
  host_deg_rna <- get(
    "deg_rna",
    envir = host_response_environment,
    inherits = FALSE
  )

  if (
    all(
      c(
        "Gene",
        "log2FoldChange",
        "pval_DESeq2"
      ) %in%
      colnames(host_deg_rna)
    )
  ) {
    host_deg_stats <- host_deg_rna %>%
      dplyr::transmute(
        feature = as.character(Gene),
        host_log2fc_original = as.numeric(
          log2FoldChange
        ),
        host_p_value = as.numeric(
          pval_DESeq2
        ),
        host_fdr = if (
          "FDR_DESeq2" %in%
          colnames(host_deg_rna)
        ) {
          as.numeric(FDR_DESeq2)
        } else {
          stats::p.adjust(
            host_p_value,
            method = "BH"
          )
        }
      ) %>%
      dplyr::distinct(
        feature,
        .keep_all = TRUE
      )

    host_vst_direction <- mofa_feature_response_stats %>%
      dplyr::filter(
        view == "host"
      ) %>%
      dplyr::select(
        feature,
        vst_mean_difference = response_effect
      )

    host_orientation_check <- host_deg_stats %>%
      dplyr::inner_join(
        host_vst_direction,
        by = "feature"
      ) %>%
      dplyr::filter(
        is.finite(host_log2fc_original),
        is.finite(vst_mean_difference)
      )

    if (nrow(host_orientation_check) >= 20) {
      host_log2fc_orientation$orientation_correlation <-
        stats::cor(
          host_orientation_check$host_log2fc_original,
          host_orientation_check$vst_mean_difference,
          method = "spearman"
        )
    }

    host_log2fc_orientation$log2fc_flipped <-
      is.finite(
        host_log2fc_orientation$orientation_correlation
      ) &&
      host_log2fc_orientation$orientation_correlation < 0

    if (host_log2fc_orientation$log2fc_flipped) {
      host_deg_stats$host_log2fc_original <-
        -host_deg_stats$host_log2fc_original
    }

    mofa_feature_response_stats <-
      mofa_feature_response_stats %>%
      dplyr::left_join(
        host_deg_stats,
        by = "feature"
      ) %>%
      dplyr::mutate(
        response_effect = ifelse(
          view == "host" &
            is.finite(host_log2fc_original),
          host_log2fc_original,
          response_effect
        ),
        effect_metric = ifelse(
          view == "host" &
            is.finite(host_log2fc_original),
          "DESeq2 log2FC oriented as pCR/non-pCR",
          effect_metric
        ),
        p_value = ifelse(
          view == "host" &
            is.finite(host_p_value),
          host_p_value,
          p_value
        ),
        fdr = ifelse(
          view == "host" &
            is.finite(host_fdr),
          host_fdr,
          fdr
        )
      ) %>%
      dplyr::select(
        -dplyr::any_of(
          c(
            "host_log2fc_original",
            "host_p_value",
            "host_fdr"
          )
        )
      )
  }
}


#=================================================================#
# 10. Host-gene annotation and publication-display filtering
#=================================================================#

mofa_host_annotation <- data.frame(
  feature = character(0),
  annotation_biotype = character(0),
  annotation_description = character(0),
  stringsAsFactors = FALSE
)

if (
  file.exists(mofa_host_annotation_file) &&
  requireNamespace(
    "readr",
    quietly = TRUE
  ) &&
  requireNamespace(
    "tidyselect",
    quietly = TRUE
  )
) {
  host_annotation_header <- colnames(
    readr::read_csv(
      mofa_host_annotation_file,
      n_max = 0,
      show_col_types = FALSE,
      progress = FALSE
    )
  )

  host_annotation_columns <- intersect(
    c(
      "Gene_Symbol",
      "Gene_ID",
      "Transcript_ID",
      "Gene_Type",
      "Gene_Biotype",
      "gene_biotype",
      "Biotype",
      "Description",
      "Gene_Description",
      "Product",
      "Gene_Name"
    ),
    host_annotation_header
  )

  if (length(host_annotation_columns) > 0) {
    host_annotation_raw <- readr::read_csv(
      mofa_host_annotation_file,
      col_select = tidyselect::all_of(
        host_annotation_columns
      ),
      show_col_types = FALSE,
      progress = FALSE
    ) %>%
      as.data.frame(
        check.names = FALSE
      )

    first_nonempty <- function(x) {
      x <- as.character(x)
      x <- x[
        !is.na(x) &
          nzchar(x)
      ]
      if (length(x) == 0) {
        return(NA_character_)
      }
      x[1]
    }

    host_annotation_raw$feature <- dplyr::coalesce(
      if (
        "Gene_Symbol" %in%
        colnames(host_annotation_raw)
      ) {
        as.character(
          host_annotation_raw$Gene_Symbol
        )
      } else {
        rep(
          NA_character_,
          nrow(host_annotation_raw)
        )
      },
      if (
        "Gene_ID" %in%
        colnames(host_annotation_raw)
      ) {
        as.character(
          host_annotation_raw$Gene_ID
        )
      } else {
        rep(
          NA_character_,
          nrow(host_annotation_raw)
        )
      },
      if (
        "Transcript_ID" %in%
        colnames(host_annotation_raw)
      ) {
        as.character(
          host_annotation_raw$Transcript_ID
        )
      } else {
        rep(
          NA_character_,
          nrow(host_annotation_raw)
        )
      }
    )

    biotype_columns <- intersect(
      c(
        "Gene_Type",
        "Gene_Biotype",
        "gene_biotype",
        "Biotype"
      ),
      colnames(host_annotation_raw)
    )

    description_columns <- intersect(
      c(
        "Description",
        "Gene_Description",
        "Product",
        "Gene_Name"
      ),
      colnames(host_annotation_raw)
    )

    host_annotation_raw$annotation_biotype <- if (
      length(biotype_columns) > 0
    ) {
      apply(
        host_annotation_raw[
          ,
          biotype_columns,
          drop = FALSE
        ],
        1,
        first_nonempty
      )
    } else {
      NA_character_
    }

    host_annotation_raw$annotation_description <- if (
      length(description_columns) > 0
    ) {
      apply(
        host_annotation_raw[
          ,
          description_columns,
          drop = FALSE
        ],
        1,
        first_nonempty
      )
    } else {
      NA_character_
    }

    mofa_host_annotation <- host_annotation_raw %>%
      dplyr::filter(
        !is.na(feature),
        nzchar(feature)
      ) %>%
      dplyr::group_by(feature) %>%
      dplyr::summarise(
        annotation_biotype = first_nonempty(
          annotation_biotype
        ),
        annotation_description = first_nonempty(
          annotation_description
        ),
        .groups = "drop"
      )
  }
}

mofa_feature_weights <- mofa_feature_weights %>%
  dplyr::left_join(
    mofa_host_annotation,
    by = "feature"
  ) %>%
  dplyr::mutate(
    annotation_text = stringr::str_to_lower(
      paste(
        dplyr::coalesce(feature, ""),
        dplyr::coalesce(feature_label, ""),
        dplyr::coalesce(annotation_biotype, ""),
        dplyr::coalesce(annotation_description, "")
      )
    ),
    display_exclusion_reason = dplyr::case_when(
      view == "host" &
        mofa_exclude_pseudogene_from_display &
        stringr::str_detect(
          annotation_text,
          "pseudogene"
        ) ~ "pseudogene",
      view == "host" &
        mofa_exclude_readthrough_from_display &
        stringr::str_detect(
          annotation_text,
          "readthrough"
        ) ~ "readthrough transcript",
      view == "host" &
        mofa_exclude_uncharacterized_from_display &
        stringr::str_detect(
          annotation_text,
          "uncharacterized|unknown function|poorly characterized|hypothetical"
        ) ~ "explicitly uncharacterized annotation",
      TRUE ~ NA_character_
    ),
    display_eligible = is.na(
      display_exclusion_reason
    )
  ) %>%
  dplyr::group_by(
    view,
    factor
  ) %>%
  dplyr::mutate(
    weight_within_view = weight /
      max(
        abs(weight),
        na.rm = TRUE
      )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(
    mofa_feature_response_stats,
    by = c(
      "view",
      "feature"
    )
  ) %>%
  dplyr::mutate(
    neg_log10_p = -log10(
      pmax(
        p_value,
        .Machine$double.xmin
      )
    )
  )

mofa_top_features <- mofa_feature_weights %>%
  dplyr::filter(
    display_eligible
  ) %>%
  dplyr::group_by(
    view,
    factor,
    direction
  ) %>%
  dplyr::slice_max(
    order_by = abs_weight,
    n = mofa_top_features_per_direction,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

write.csv(
  mofa_feature_weights,
  "results/mofa/MOFA_v30_feature_weights_annotated.csv",
  row.names = FALSE
)


#=================================================================#
# 11. Factor ranking and integrated visualization data
#=================================================================#

mofa_factor_r2_wide <- mofa_variance_explained %>%
  dplyr::mutate(
    factor = as.character(factor)
  ) %>%
  dplyr::select(
    factor,
    view,
    r2
  ) %>%
  tidyr::pivot_wider(
    names_from = view,
    values_from = r2,
    values_fill = 0
  )

for (view_name in required_views) {
  if (!view_name %in% colnames(mofa_factor_r2_wide)) {
    mofa_factor_r2_wide[[view_name]] <- 0
  }
}

mofa_factor_balance <- mofa_variance_explained %>%
  dplyr::mutate(
    factor = as.character(factor)
  ) %>%
  dplyr::group_by(factor) %>%
  dplyr::summarise(
    r2_sum = sum(
      r2,
      na.rm = TRUE
    ),
    strongest_view_fraction = ifelse(
      sum(r2, na.rm = TRUE) > 0,
      max(r2, na.rm = TRUE) /
        sum(r2, na.rm = TRUE),
      NA_real_
    ),
    view_entropy = {
      p <- r2 /
        sum(
          r2,
          na.rm = TRUE
        )
      p <- p[
        is.finite(p) &
          p > 0
      ]
      if (length(p) <= 1) {
        0
      } else {
        -sum(
          p * log(p)
        ) /
          log(
            length(required_views)
          )
      }
    },
    active_views_primary = sum(
      r2 >= mofa_active_view_r2,
      na.rm = TRUE
    ),
    active_views_strict = sum(
      r2 >= mofa_strict_active_view_r2,
      na.rm = TRUE
    ),
    minimum_active_r2 = if (
      any(
        r2 >= mofa_active_view_r2,
        na.rm = TRUE
      )
    ) {
      min(
        r2[
          r2 >= mofa_active_view_r2
        ],
        na.rm = TRUE
      )
    } else {
      0
    },
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    mofa_factor_r2_wide,
    by = "factor"
  ) %>%
  dplyr::mutate(
    effective_views = exp(
      view_entropy *
        log(
          length(required_views)
        )
    ),
    balanced_shared_score =
      r2_sum *
      view_entropy *
      (1 - strongest_view_fraction),
    microbial_chain_geomean = (
      pmax(species, 1e-06) *
        pmax(ko, 1e-06) *
        pmax(metabolite, 1e-06)
    ) ^ (1 / 3),
    metabolite_host_bridge = sqrt(
      pmax(metabolite, 0) *
        pmax(host, 0)
    ),
    pan_omic_geomean = (
      pmax(species, 1e-06) *
        pmax(ko, 1e-06) *
        pmax(metabolite, 1e-06) *
        pmax(host, 1e-06)
    ) ^ (1 / 4),
    shared_factor_priority =
      balanced_shared_score +
      microbial_chain_geomean +
      metabolite_host_bridge
  )

mofa_factor_qc_summary <-
  mofa_factor_total_level_correlations %>%
  dplyr::group_by(factor) %>%
  dplyr::summarise(
    total_level_qc_flag = any(
      flag_abs_r_gt_0_75,
      na.rm = TRUE
    ),
    total_level_qc_views = paste(
      view[
        flag_abs_r_gt_0_75 %in% TRUE
      ],
      collapse = ";"
    ),
    .groups = "drop"
  )

mofa_primary_clinical_priority <-
  mofa_factor_clinical_associations %>%
  dplyr::filter(
    analysis_set == "primary_all_samples",
    term %in% c(
      "TRG_plotpCR",
      "TimepointOngoing:TRG_plotpCR"
    )
  ) %>%
  dplyr::group_by(factor) %>%
  dplyr::summarise(
    minimum_response_permutation_p = min(
      permutation_p,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    minimum_response_permutation_p = ifelse(
      is.infinite(
        minimum_response_permutation_p
      ),
      NA_real_,
      minimum_response_permutation_p
    )
  )

mofa_factor_summary <- mofa_factor_summary %>%
  dplyr::left_join(
    mofa_factor_balance,
    by = "factor"
  ) %>%
  dplyr::left_join(
    mofa_factor_qc_summary,
    by = "factor"
  ) %>%
  dplyr::left_join(
    mofa_primary_clinical_priority,
    by = "factor"
  ) %>%
  dplyr::arrange(
    dplyr::desc(total_r2)
  )

write.csv(
  mofa_factor_summary,
  "results/mofa/MOFA_v30_factor_summary.csv",
  row.names = FALSE
)

p_mofa_shared_factor_priority <- ggplot2::ggplot(
  mofa_factor_summary,
  ggplot2::aes(
    x = strongest_view_fraction,
    y = effective_views,
    size = total_r2,
    fill = microbial_chain_geomean,
    label = factor
  )
) +
  ggplot2::geom_point(
    shape = 21,
    color = "grey25",
    stroke = 0.45,
    alpha = 0.90
  ) +
  ggplot2::geom_text(
    size = 2.8,
    vjust = -0.85,
    check_overlap = TRUE
  ) +
  ggplot2::scale_size_continuous(
    range = c(2.5, 7.0),
    labels = scales::label_percent(
      accuracy = 1
    ),
    name = "Sum of view R2"
  ) +
  ggplot2::scale_fill_gradient(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "white",
    high = "grey25",
    labels = scales::label_percent(
      accuracy = 0.1
    ),
    name = paste0(
      "Species-KO-metabolite\n",
      "geometric mean R2"
    )
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 1),
    labels = scales::label_percent(
      accuracy = 1
    )
  ) +
  ggplot2::scale_y_continuous(
    limits = c(1, length(required_views)),
    breaks = seq_len(
      length(required_views)
    )
  ) +
  ggplot2::labs(
    x = "Fraction of factor R2 carried by the strongest view",
    y = "Effective number of contributing views"
  ) +
  ggplot2::theme_classic(
    base_size = 10.5
  ) +
  ggplot2::theme(
    legend.position = "right",
    plot.margin = ggplot2::margin(
      5, 8, 5, 5
    )
  )

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename = "figures/mofa/internal_v16/MOFA_v16_shared_factor_priority.svg",
  plot = p_mofa_shared_factor_priority,
  width = 7.2,
  height = 4.6,
  device = "svg"
)
}

if (
  !is.null(mofa_selected_factor_override) &&
  mofa_selected_factor_override %in%
    mofa_factor_summary$factor
) {
  mofa_selected_factor <-
    mofa_selected_factor_override
  mofa_selected_factor_reason <-
    "User-specified override"
} else {
  clinically_prioritized_factors <-
    mofa_factor_summary %>%
    dplyr::filter(
      active_views >= 2,
      !dplyr::coalesce(
        total_level_qc_flag,
        FALSE
      ),
      is.finite(
        minimum_response_permutation_p
      ),
      minimum_response_permutation_p < 0.10
    ) %>%
    dplyr::arrange(
      minimum_response_permutation_p,
      dplyr::desc(
        shared_factor_priority
      )
    )

  if (nrow(clinically_prioritized_factors) > 0) {
    mofa_selected_factor <-
      clinically_prioritized_factors$factor[1]
    mofa_selected_factor_reason <- paste0(
      "Exploratory clinical prioritization: minimum subject-permutation p = ",
      format.pval(
        clinically_prioritized_factors$minimum_response_permutation_p[1],
        digits = 2,
        eps = 0.001
      )
    )
  } else {
    mofa_selected_factor <-
      mofa_factor_summary %>%
      dplyr::filter(
        active_views >= 2,
        !dplyr::coalesce(
          total_level_qc_flag,
          FALSE
        )
      ) %>%
      dplyr::arrange(
        dplyr::desc(
          balanced_shared_score
        )
      ) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::pull(factor)

    if (length(mofa_selected_factor) == 0) {
      mofa_selected_factor <-
        mofa_factor_summary$factor[1]
    }

    mofa_selected_factor_reason <-
      "Highest balanced shared-variance score among factors without a total-level QC flag"
  }
}

mofa_selected_factor <- as.character(
  mofa_selected_factor
)

mofa_factor_selection_record <- data.frame(
  selected_factor = mofa_selected_factor,
  selection_reason = mofa_selected_factor_reason,
  stringsAsFactors = FALSE
)


#=================================================================#
# 12. Publication-oriented figures
#=================================================================#

mofa_response_colors <- c(
  pCR = "#4FAE9A",
  non_pCR = "#DE7872",
  CR = "#4FAE9A",
  nonCR = "#DE7872"
)

response_colors <- mofa_response_colors

timepoint_colors <- c(
  Before = "#667A8A",
  Ongoing = "#C59A5B"
)

p_mofa_data_availability <- ggplot2::ggplot(
  mofa_availability_long,
  ggplot2::aes(
    x = SampleID,
    y = view_label,
    fill = factor(available)
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.22
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      `0` = "grey92",
      `1` = "grey30"
    ),
    labels = c(
      `0` = "Missing",
      `1` = "Measured"
    ),
    name = NULL
  ) +
  ggplot2::labs(
    title = "Four-omics data availability used for MOFA",
    subtitle = paste0(
      "Samples measured in at least ",
      mofa_min_views_per_sample,
      " views are retained; missing views are encoded as NA."
    ),
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 10.5
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 9.0,
      color = "grey35"
    ),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(
      face = "plain"
    ),
    legend.position = "top"
  )

mofa_variance_explained$factor <- factor(
  as.character(
    mofa_variance_explained$factor
  ),
  levels = rev(
    mofa_factor_summary$factor
  )
)

p_mofa_variance_heatmap <- ggplot2::ggplot(
  mofa_variance_explained,
  ggplot2::aes(
    x = view_label,
    y = factor,
    fill = r2
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.35
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = ifelse(
        r2 >= 0.005,
        scales::percent(
          r2,
          accuracy = 0.1
        ),
        ""
      )
    ),
    size = 2.6
  ) +
  ggplot2::scale_fill_gradient(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "white",
    high = "grey20",
    labels = scales::label_percent(
      accuracy = 1
    ),
    name = "Factor R²"
  ) +
  ggplot2::labs(
    title = "Variance explained by latent factor and omics view",
    subtitle = paste0(
      "Active view threshold: R² >= ",
      scales::percent(
        mofa_active_view_r2,
        accuracy = 1
      ),
      "."
    ),
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 10.5
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 8.8,
      color = "grey35"
    ),
    axis.text.x = ggplot2::element_text(
      face = "plain",
      angle = 25,
      hjust = 1
    ),
    axis.text.y = ggplot2::element_text(
      face = "plain"
    ),
    axis.ticks = ggplot2::element_blank(),
    legend.position = "right"
  )

mofa_view_total_r2 <- mofa_variance_total %>%
  dplyr::group_by(
    view,
    view_label
  ) %>%
  dplyr::summarise(
    total_r2 = mean(
      total_r2,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

p_mofa_view_total_r2 <- ggplot2::ggplot(
  mofa_view_total_r2,
  ggplot2::aes(
    x = view_label,
    y = total_r2,
    fill = view_label
  )
) +
  ggplot2::geom_col(
    width = 0.62,
    color = "grey25",
    linewidth = 0.42
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = scales::percent(
        total_r2,
        accuracy = 0.1
      )
    ),
    vjust = -0.25,
    size = 3.0
  ) +
  ggplot2::scale_fill_manual(
    values = view_colors,
    guide = "none"
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_percent(
      accuracy = 1
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.12)
    )
  ) +
  ggplot2::labs(
    title = "Total model variance explained",
    x = NULL,
    y = "Total R²"
  ) +
  ggplot2::theme_classic(
    base_size = 9.8
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 10.5
    ),
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(
      4,
      45,
      0,
      42
    )
  )

p_mofa_variance_overview <-
  p_mofa_view_total_r2 /
  p_mofa_variance_heatmap +
  patchwork::plot_layout(
    heights = c(
      0.72,
      3.0
    )
  )

mofa_primary_clinical_plot_data <-
  mofa_factor_clinical_associations %>%
  dplyr::filter(
    analysis_set == "primary_all_samples"
  ) %>%
  dplyr::mutate(
    factor = factor(
      factor,
      levels = rev(
        mofa_factor_summary$factor
      )
    ),
    term_label = factor(
      dplyr::recode(
        term,
        TimepointOngoing = "Timepoint: Ongoing vs Before",
        TRG_plotpCR = "Response: pCR vs non-pCR at Before",
        `TimepointOngoing:TRG_plotpCR` =
          "Interaction: differential treatment change"
      ),
      levels = c(
        "Timepoint: Ongoing vs Before",
        "Response: pCR vs non-pCR at Before",
        "Interaction: differential treatment change"
      )
    ),
    significance_label = dplyr::case_when(
      is.finite(permutation_p) &
        permutation_p < 0.001 ~ "***",
      is.finite(permutation_p) &
        permutation_p < 0.01 ~ "**",
      is.finite(permutation_p) &
        permutation_p < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

p_mofa_factor_clinical_forest <- ggplot2::ggplot(
  mofa_primary_clinical_plot_data,
  ggplot2::aes(
    x = estimate,
    y = factor
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    color = "grey70",
    linewidth = 0.45
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin = ci_low,
      xmax = ci_high
    ),
    height = 0.14,
    linewidth = 0.65
  ) +
  ggplot2::geom_point(
    shape = 21,
    size = 2.8,
    fill = "white",
    stroke = 0.65
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = significance_label
    ),
    nudge_y = 0.28,
    size = 3.0
  ) +
  ggplot2::facet_wrap(
    ~ term_label,
    scales = "free_x",
    nrow = 1
  ) +
  ggplot2::labs(
    title = "Post-hoc clinical associations of MOFA factor scores",
    subtitle = paste0(
      "Linear mixed models include SubjectID random intercepts; stars use ",
      scales::comma(
        mofa_factor_permutations
      ),
      " SubjectID-level response permutations."
    ),
    x = "Factor-score coefficient (95% CI)",
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 10.2
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 8.7,
      color = "grey35"
    ),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(
      face = "plain",
      size = 9.0
    ),
    axis.text.y = ggplot2::element_text(
      face = "plain"
    )
  )

mofa_clinical_tile_data <-
  mofa_primary_clinical_plot_data %>%
  dplyr::mutate(
    standardized_effect = estimate /
      std_error,
    neg_log10_permutation_p = -log10(
      pmax(
        permutation_p,
        .Machine$double.xmin
      )
    )
  )

clinical_effect_limit <- max(
  abs(
    mofa_clinical_tile_data$standardized_effect
  ),
  na.rm = TRUE
)

if (
  !is.finite(clinical_effect_limit) ||
  clinical_effect_limit <= 0
) {
  clinical_effect_limit <- 1
}

p_mofa_clinical_effect_tiles <- ggplot2::ggplot(
  mofa_clinical_tile_data,
  ggplot2::aes(
    x = term_label,
    y = factor,
    fill = standardized_effect
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.35
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      size = neg_log10_permutation_p
    ),
    shape = 21,
    fill = "white",
    color = "grey25",
    stroke = 0.40
  ) +
  ggplot2::scale_fill_gradient2(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "#3B75AF",
    mid = "white",
    high = "#C44E52",
    midpoint = 0,
    limits = c(
      -clinical_effect_limit,
      clinical_effect_limit
    ),
    name = "Coefficient / SE"
  ) +
  ggplot2::scale_size_continuous(
    range = c(
      0.5,
      4.0
    ),
    name = "-log10 permutation p"
  ) +
  ggplot2::labs(
    title = "Clinical-effect fingerprint",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 9.8
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    axis.text.x = ggplot2::element_text(
      angle = 30,
      hjust = 1,
      size = 8.3
    ),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )

p_mofa_factor_landscape <-
  p_mofa_variance_heatmap +
  p_mofa_clinical_effect_tiles +
  patchwork::plot_layout(
    widths = c(
      1.55,
      1
    )
  )

mofa_score_long <- mofa_factor_scores %>%
  dplyr::select(
    SampleID,
    SubjectID,
    Timepoint,
    TRG_plot,
    n_views,
    dplyr::all_of(
      colnames(mofa_factor_matrix)
    )
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      colnames(mofa_factor_matrix)
    ),
    names_to = "factor",
    values_to = "factor_score"
  ) %>%
  dplyr::group_by(factor) %>%
  dplyr::mutate(
    factor_score_z = as.numeric(
      base::scale(factor_score)
    )
  ) %>%
  dplyr::ungroup()

mofa_factor_names <- colnames(mofa_factor_matrix)
mofa_factor_names <- mofa_factor_names[
  !is.na(mofa_factor_names) &
    nzchar(mofa_factor_names) &
    grepl(
      "^Factor[0-9]+$",
      mofa_factor_names
    )
]

order_heatmap_group <- function(group_data, factor_names) {
  if (nrow(group_data) <= 2) {
    return(as.character(group_data$SampleID))
  }

  group_matrix <- as.matrix(
    group_data[
      ,
      factor_names,
      drop = FALSE
    ]
  )

  group_matrix <- scale(group_matrix)
  group_matrix[!is.finite(group_matrix)] <- 0

  cluster_result <- stats::hclust(
    stats::dist(group_matrix),
    method = "ward.D2"
  )

  as.character(group_data$SampleID[cluster_result$order])
}

sample_order_for_heatmap <- unlist(
  lapply(
    mofa_heatmap_response_order,
    function(response_name) {
      unlist(
        lapply(
          mofa_heatmap_timepoint_order,
          function(timepoint_name) {
            group_data <- mofa_factor_scores %>%
              dplyr::filter(
                as.character(TRG_plot) == response_name,
                as.character(Timepoint) == timepoint_name
              )

            order_heatmap_group(
              group_data,
              mofa_factor_names
            )
          }
        ),
        use.names = FALSE
      )
    }
  ),
  use.names = FALSE
)

sample_order_for_heatmap <- unique(
  c(
    sample_order_for_heatmap,
    setdiff(
      as.character(mofa_factor_scores$SampleID),
      sample_order_for_heatmap
    )
  )
)

mofa_heatmap_group_sizes <- mofa_factor_scores %>%
  dplyr::mutate(
    response_order = match(
      as.character(TRG_plot),
      mofa_heatmap_response_order
    ),
    timepoint_order = match(
      as.character(Timepoint),
      mofa_heatmap_timepoint_order
    )
  ) %>%
  dplyr::count(
    response_order,
    timepoint_order,
    .drop = FALSE,
    name = "n_samples"
  ) %>%
  dplyr::filter(n_samples > 0) %>%
  dplyr::arrange(response_order, timepoint_order)

mofa_heatmap_group_boundaries <- head(
  cumsum(mofa_heatmap_group_sizes$n_samples) + 0.5,
  -1
)

mofa_score_long$SampleID <- factor(
  mofa_score_long$SampleID,
  levels = sample_order_for_heatmap
)

mofa_score_long$factor <- factor(
  mofa_score_long$factor,
  levels = rev(mofa_factor_summary$factor)
)

calculate_hedges_g <- function(x_group1, x_group0) {
  x_group1 <- x_group1[is.finite(x_group1)]
  x_group0 <- x_group0[is.finite(x_group0)]

  if (length(x_group1) < 2 || length(x_group0) < 2) {
    return(NA_real_)
  }

  pooled_sd <- sqrt(
    (
      (length(x_group1) - 1) * stats::var(x_group1) +
        (length(x_group0) - 1) * stats::var(x_group0)
    ) /
      (length(x_group1) + length(x_group0) - 2)
  )

  if (!is.finite(pooled_sd) || pooled_sd <= 0) {
    return(NA_real_)
  }

  cohens_d <- (mean(x_group1) - mean(x_group0)) / pooled_sd
  correction <- 1 - 3 /
    (4 * (length(x_group1) + length(x_group0)) - 9)

  correction * cohens_d
}

calculate_paired_dz <- function(before_value, ongoing_value) {
  paired_delta <- ongoing_value - before_value
  paired_delta <- paired_delta[is.finite(paired_delta)]

  if (length(paired_delta) < 3) {
    return(NA_real_)
  }

  delta_sd <- stats::sd(paired_delta)

  if (!is.finite(delta_sd) || delta_sd <= 0) {
    return(NA_real_)
  }

  mean(paired_delta) / delta_sd
}

calculate_shifted_log2_ratio <- function(group1, group0) {
  finite_values <- c(group1, group0)
  finite_values <- finite_values[is.finite(finite_values)]

  if (length(finite_values) == 0) {
    return(NA_real_)
  }

  shift_value <- 1 - min(finite_values)
  numerator <- mean(group1 + shift_value, na.rm = TRUE)
  denominator <- mean(group0 + shift_value, na.rm = TRUE)

  if (
    !is.finite(numerator) ||
    !is.finite(denominator) ||
    numerator <= 0 ||
    denominator <= 0
  ) {
    return(NA_real_)
  }

  log2(numerator / denominator)
}

mofa_factor_display_effects <- dplyr::bind_rows(
  lapply(
    mofa_factor_names,
    function(factor_name) {
      baseline_data <- mofa_factor_scores %>%
        dplyr::filter(
          as.character(Timepoint) == "Before",
          !is.na(TRG_plot)
        ) %>%
        dplyr::distinct(SubjectID, .keep_all = TRUE)

      pcr_values <- baseline_data[
        as.character(baseline_data$TRG_plot) == "pCR",
        factor_name,
        drop = TRUE
      ]
      non_pcr_values <- baseline_data[
        as.character(baseline_data$TRG_plot) == "non_pCR",
        factor_name,
        drop = TRUE
      ]

      paired_wide <- mofa_factor_scores %>%
        dplyr::filter(
          as.character(Timepoint) %in% c("Before", "Ongoing")
        ) %>%
        dplyr::select(
          SubjectID,
          Timepoint,
          dplyr::all_of(factor_name)
        ) %>%
        tidyr::pivot_wider(
          names_from = Timepoint,
          values_from = dplyr::all_of(factor_name)
        ) %>%
        dplyr::filter(
          is.finite(Before),
          is.finite(Ongoing)
        )

      data.frame(
        factor = factor_name,
        baseline_hedges_g = calculate_hedges_g(
          pcr_values,
          non_pcr_values
        ),
        baseline_shifted_log2_ratio = calculate_shifted_log2_ratio(
          pcr_values,
          non_pcr_values
        ),
        paired_dz = calculate_paired_dz(
          paired_wide$Before,
          paired_wide$Ongoing
        ),
        paired_shifted_log2_ratio = calculate_shifted_log2_ratio(
          paired_wide$Ongoing,
          paired_wide$Before
        ),
        n_baseline_pCR = sum(is.finite(pcr_values)),
        n_baseline_non_pCR = sum(is.finite(non_pcr_values)),
        n_paired = nrow(paired_wide),
        stringsAsFactors = FALSE
      )
    }
  )
)

if (identical(mofa_heatmap_effect_metric, "shifted_log2_ratio")) {
  mofa_factor_display_effects$response_display_effect <-
    mofa_factor_display_effects$baseline_shifted_log2_ratio
  mofa_factor_display_effects$timepoint_display_effect <-
    mofa_factor_display_effects$paired_shifted_log2_ratio
  response_effect_label <- "pCR/non-pCR shifted log2 ratio"
  timepoint_effect_label <- "Ongoing/Before shifted log2 ratio"
} else {
  mofa_factor_display_effects$response_display_effect <-
    mofa_factor_display_effects$baseline_hedges_g
  mofa_factor_display_effects$timepoint_display_effect <-
    mofa_factor_display_effects$paired_dz
  response_effect_label <- "pCR vs non-pCR at Before (Hedges g)"
  timepoint_effect_label <- "Ongoing vs Before (paired dz)"
}

p_mofa_factor_score_heatmap_core <- ggplot2::ggplot(
  mofa_score_long,
  ggplot2::aes(
    x = SampleID,
    y = factor,
    fill = factor_score_z
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.14
  ) +
  ggplot2::geom_vline(
    xintercept = mofa_heatmap_group_boundaries,
    linewidth = 0.42,
    color = "grey40"
  ) +
  ggplot2::scale_fill_gradient2(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "#3B75AF",
    mid = "white",
    high = "#C44E52",
    midpoint = 0,
    name = "Factor score
(z-score)"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 9.6) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(face = "plain", size = 8.2),
    legend.position = "right",
    plot.margin = ggplot2::margin(0, 3, 4, 4)
  )

mofa_heatmap_metadata <- mofa_factor_scores %>%
  dplyr::mutate(
    SampleID = factor(SampleID, levels = sample_order_for_heatmap),
    Response = factor(
      as.character(TRG_plot),
      levels = mofa_heatmap_response_order
    ),
    `Time point` = factor(
      as.character(Timepoint),
      levels = mofa_heatmap_timepoint_order
    )
  )

p_mofa_heatmap_timepoint_annotation <- ggplot2::ggplot(
  mofa_heatmap_metadata,
  ggplot2::aes(
    x = SampleID,
    y = "Time point",
    fill = `Time point`
  )
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.14) +
  ggplot2::geom_vline(
    xintercept = mofa_heatmap_group_boundaries,
    linewidth = 0.42,
    color = "grey40"
  ) +
  ggplot2::scale_fill_manual(
    values = timepoint_colors,
    drop = FALSE,
    name = "Time point"
  ) +
  ggplot2::theme_void(base_size = 8.8) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(
      face = "plain",
      color = "grey15"
    ),
    legend.position = "top",
    legend.direction = "horizontal",
    plot.margin = ggplot2::margin(0, 3, 0, 4)
  )

p_mofa_heatmap_response_annotation <- ggplot2::ggplot(
  mofa_heatmap_metadata,
  ggplot2::aes(
    x = SampleID,
    y = "Response",
    fill = Response
  )
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.14) +
  ggplot2::geom_vline(
    xintercept = mofa_heatmap_group_boundaries,
    linewidth = 0.42,
    color = "grey40"
  ) +
  ggplot2::scale_fill_manual(
    values = mofa_response_colors,
    drop = FALSE,
    name = "Response"
  ) +
  ggplot2::theme_void(base_size = 8.8) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(
      face = "plain",
      color = "grey15"
    ),
    legend.position = "top",
    legend.direction = "horizontal",
    plot.margin = ggplot2::margin(0, 3, 0, 4)
  )

mofa_heatmap_availability <- mofa_availability_long %>%
  dplyr::filter(SampleID %in% sample_order_for_heatmap) %>%
  dplyr::mutate(
    SampleID = factor(SampleID, levels = sample_order_for_heatmap),
    view_label = factor(
      unname(view_labels[view]),
      levels = rev(unname(view_labels))
    ),
    availability_fill = ifelse(
      available == 1,
      as.character(view_label),
      "Missing"
    )
  )

p_mofa_heatmap_availability_annotation <- ggplot2::ggplot(
  mofa_heatmap_availability,
  ggplot2::aes(
    x = SampleID,
    y = view_label,
    fill = availability_fill
  )
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.14) +
  ggplot2::geom_vline(
    xintercept = mofa_heatmap_group_boundaries,
    linewidth = 0.42,
    color = "grey40"
  ) +
  ggplot2::scale_fill_manual(
    values = c(view_colors, Missing = "grey94"),
    guide = "none"
  ) +
  ggplot2::labs(x = NULL, y = "Measured views") +
  ggplot2::theme_void(base_size = 8.5) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(
      face = "plain",
      color = "grey15",
      size = 7.7
    ),
    axis.title.y = ggplot2::element_text(
      face = "bold",
      size = 8.3,
      angle = 90
    ),
    plot.margin = ggplot2::margin(0, 3, 1, 4)
  )

p_mofa_factor_score_annotation <-
  p_mofa_heatmap_timepoint_annotation /
  p_mofa_heatmap_response_annotation /
  p_mofa_heatmap_availability_annotation +
  patchwork::plot_layout(heights = c(0.42, 0.42, 1.15), guides = "collect") &
  ggplot2::theme(legend.position = "top")

mofa_factor_effect_tile_data <- mofa_factor_display_effects %>%
  dplyr::select(
    factor,
    response_display_effect,
    timepoint_display_effect
  ) %>%
  tidyr::pivot_longer(
    cols = c(response_display_effect, timepoint_display_effect),
    names_to = "contrast",
    values_to = "effect"
  ) %>%
  dplyr::mutate(
    factor = factor(
      factor,
      levels = rev(mofa_factor_summary$factor)
    ),
    contrast = factor(
      contrast,
      levels = c(
        "response_display_effect",
        "timepoint_display_effect"
      ),
      labels = c(
        response_effect_label,
        timepoint_effect_label
      )
    )
  )

factor_effect_limit <- max(
  abs(mofa_factor_effect_tile_data$effect),
  na.rm = TRUE
)

if (!is.finite(factor_effect_limit) || factor_effect_limit <= 0) {
  factor_effect_limit <- 1
}

p_mofa_factor_score_effect_tiles <- ggplot2::ggplot(
  mofa_factor_effect_tile_data,
  ggplot2::aes(
    x = contrast,
    y = factor,
    fill = effect
  )
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.30) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = ifelse(is.finite(effect), sprintf("%.2f", effect), "")
    ),
    size = 2.45
  ) +
  ggplot2::scale_fill_gradient2(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "#3B75AF",
    mid = "white",
    high = "#C44E52",
    midpoint = 0,
    limits = c(-factor_effect_limit, factor_effect_limit),
    oob = scales::squish,
    name = "Standardized
effect"
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_classic(base_size = 8.3) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      size = 7.0
    ),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    legend.position = "right",
    plot.margin = ggplot2::margin(0, 4, 4, 2)
  )

p_mofa_factor_score_heatmap <-
  (
    p_mofa_factor_score_annotation |
      patchwork::plot_spacer()
  ) /
  (
    p_mofa_factor_score_heatmap_core |
      p_mofa_factor_score_effect_tiles
  ) +
  patchwork::plot_layout(
    heights = c(0.98, 3.2),
    widths = c(4.5, 1.25)
  )

p_mofa_paired_factor_change <- ggplot2::ggplot(
  mofa_score_long,
  ggplot2::aes(
    x = Timepoint,
    y = factor_score,
    group = SubjectID,
    color = TRG_plot
  )
) +
  ggplot2::geom_line(
    alpha = 0.30,
    linewidth = 0.45
  ) +
  ggplot2::geom_point(
    size = 1.7,
    alpha = 0.82
  ) +
  ggplot2::facet_wrap(
    ~ factor,
    scales = "free_y",
    ncol = 4
  ) +
  ggplot2::scale_color_manual(
    values = response_colors,
    name = "Response"
  ) +
  ggplot2::labs(
    title = "Paired treatment trajectories of MOFA factor scores",
    subtitle = "Lines connect Before and Ongoing observations from the same subject.",
    x = NULL,
    y = "MOFA factor score"
  ) +
  ggplot2::theme_classic(
    base_size = 9.8
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 8.7,
      color = "grey35"
    ),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(
      face = "plain"
    ),
    legend.position = "top"
  )

mofa_factor_pair_groups <- split(
  mofa_factor_summary$factor,
  ceiling(
    seq_along(
      mofa_factor_summary$factor
    ) / 2
  )
)

p_mofa_factor_pair_maps <- list()

for (pair_index in seq_along(mofa_factor_pair_groups)) {
  pair_factors <- mofa_factor_pair_groups[[pair_index]]

  if (length(pair_factors) < 2) {
    next
  }

  x_factor_summary <- mofa_factor_summary[
    match(
      pair_factors[1],
      mofa_factor_summary$factor
    ),
    ,
    drop = FALSE
  ]

  y_factor_summary <- mofa_factor_summary[
    match(
      pair_factors[2],
      mofa_factor_summary$factor
    ),
    ,
    drop = FALSE
  ]

  p_mofa_factor_pair_maps[[pair_index]] <-
    ggplot2::ggplot(
      mofa_factor_scores,
      ggplot2::aes(
        x = .data[[pair_factors[1]]],
        y = .data[[pair_factors[2]]],
        fill = TRG_plot,
        shape = Timepoint
      )
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.35,
      color = "grey85"
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.35,
      color = "grey85"
    ) +
    ggplot2::geom_point(
      size = 2.65,
      color = "grey25",
      stroke = 0.55,
      alpha = 0.88
    ) +
    ggplot2::scale_fill_manual(
      values = response_colors,
      name = "Response"
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        Before = 21,
        Ongoing = 24
      ),
      name = "Timepoint"
    ) +
    ggplot2::labs(
      title = paste(
        pair_factors,
        collapse = " versus "
      ),
      subtitle = paste0(
        pair_factors[1],
        ": ",
        view_labels[
          x_factor_summary$strongest_view
        ],
        " dominant (",
        scales::percent(
          x_factor_summary$strongest_view_r2,
          accuracy = 0.1
        ),
        "); ",
        pair_factors[2],
        ": ",
        view_labels[
          y_factor_summary$strongest_view
        ],
        " dominant (",
        scales::percent(
          y_factor_summary$strongest_view_r2,
          accuracy = 0.1
        ),
        ")."
      ),
      x = pair_factors[1],
      y = pair_factors[2]
    ) +
    ggplot2::theme_classic(
      base_size = 10.0
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(
        size = 7.9,
        color = "grey35"
      ),
      legend.position = "bottom"
    )

  if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
    filename = paste0(
      "figures/mofa/internal_v16/MOFA_v16_factor_pair_",
      pair_factors[1],
      "_",
      pair_factors[2],
      ".svg"
    ),
    plot = p_mofa_factor_pair_maps[[pair_index]],
    width = 5.4,
    height = 4.8,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}
}

p_mofa_factor_pair_maps_all <- patchwork::wrap_plots(
  p_mofa_factor_pair_maps,
  ncol = 2,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "MOFA factor-pair sample maps",
    subtitle = "Pairs follow total factor-level variance explained; clinical labels are overlaid only after unsupervised fitting."
  ) &
  ggplot2::theme(
    legend.position = "bottom"
  )

build_factor_feature_plot <- function(
    factor_name,
    top_features_per_direction = mofa_top_features_per_direction,
    add_annotation = TRUE
) {
  active_views <- mofa_variance_explained %>%
    dplyr::filter(
      as.character(factor) == factor_name,
      r2 >= mofa_active_view_r2
    ) %>%
    dplyr::arrange(
      match(
        as.character(view),
        required_views
      )
    ) %>%
    dplyr::pull(view) %>%
    as.character()

  if (length(active_views) == 0) {
    active_views <- mofa_factor_summary$strongest_view[
      match(
        factor_name,
        mofa_factor_summary$factor
      )
    ]
  }

  factor_weight_data <- mofa_feature_weights %>%
    dplyr::filter(
      factor == factor_name,
      view %in% active_views,
      display_eligible
    ) %>%
    dplyr::group_by(
      view,
      direction
    ) %>%
    dplyr::slice_max(
      abs_weight,
      n = top_features_per_direction,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      view_label = factor(
        view_label,
        levels = unname(
          view_labels[active_views]
        )
      ),
      feature_plot = paste(
        view,
        feature_label,
        sep = "|||"
      )
    ) %>%
    dplyr::arrange(
      match(
        view,
        active_views
      ),
      weight_within_view
    )

  factor_weight_data$feature_plot <- factor(
    factor_weight_data$feature_plot,
    levels = unique(
      factor_weight_data$feature_plot
    )
  )

  response_effect_limit <- stats::quantile(
    abs(
      factor_weight_data$response_effect
    ),
    probs = 0.95,
    na.rm = TRUE,
    names = FALSE
  )

  if (
    !is.finite(response_effect_limit) ||
    response_effect_limit <= 0
  ) {
    response_effect_limit <- 1
  }

  p_weight <- ggplot2::ggplot(
    factor_weight_data,
    ggplot2::aes(
      x = weight_within_view,
      y = feature_plot,
      fill = view_label
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.42,
      color = "grey60"
    ) +
    ggplot2::geom_col(
      width = 0.72,
      color = "grey25",
      linewidth = 0.25
    ) +
    ggplot2::facet_grid(
      view_label ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_y_discrete(
      labels = function(x) {
        sub(
          "^[^|]+[|][|][|]",
          "",
          x
        )
      }
    ) +
    ggplot2::scale_x_continuous(
      limits = c(
        -1.05,
        1.05
      ),
      breaks = c(
        -1,
        -0.5,
        0,
        0.5,
        1
      )
    ) +
    ggplot2::scale_fill_manual(
      values = view_colors,
      guide = "none"
    ) +
    ggplot2::labs(
      x = "Within-view normalized MOFA weight",
      y = NULL
    ) +
    ggplot2::theme_classic(
      base_size = 9.6
    ) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_text(
        face = "plain",
        angle = 0,
        hjust = 0
      ),
      axis.text.y = ggplot2::element_text(
        size = 7.6
      ),
      panel.spacing.y = grid::unit(
        0.65,
        "lines"
      )
    )

  p_response_effect <- ggplot2::ggplot(
    factor_weight_data,
    ggplot2::aes(
      x = "Effect",
      y = feature_plot,
      fill = response_effect
    )
  ) +
    ggplot2::geom_tile(
      color = "white",
      linewidth = 0.30
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = ifelse(
          is.finite(response_effect),
          sprintf(
            "%.2f",
            response_effect
          ),
          ""
        )
      ),
      size = 2.35
    ) +
    ggplot2::facet_grid(
      view_label ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_fill_gradient2(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
      low = "#3B75AF",
      mid = "white",
      high = "#C44E52",
      midpoint = 0,
      limits = c(
        -response_effect_limit,
        response_effect_limit
      ),
      oob = scales::squish,
      name = "pCR/non-pCR
effect"
    ) +
    ggplot2::labs(
      x = "pCR vs non-pCR
log2FC/effect",
      y = NULL
    ) +
    ggplot2::theme_void(
      base_size = 9.0
    ) +
    ggplot2::theme(
      strip.text = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(
        4,
        2,
        22,
        2
      ),
      axis.title.x = ggplot2::element_text(
        size = 8.0
      )
    )

  p_response_p <- ggplot2::ggplot(
    factor_weight_data,
    ggplot2::aes(
      x = "P",
      y = feature_plot,
      fill = neg_log10_p
    )
  ) +
    ggplot2::geom_tile(
      color = "white",
      linewidth = 0.30
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = ifelse(
          is.finite(neg_log10_p),
          sprintf(
            "%.1f",
            neg_log10_p
          ),
          ""
        )
      ),
      size = 2.35
    ) +
    ggplot2::facet_grid(
      view_label ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_fill_gradient(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
      low = "white",
      high = "grey20",
      name = "-log10(p)"
    ) +
    ggplot2::labs(
      x = "-log10(p)",
      y = NULL
    ) +
    ggplot2::theme_void(
      base_size = 9.0
    ) +
    ggplot2::theme(
      strip.text = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(
        4,
        2,
        22,
        2
      ),
      axis.title.x = ggplot2::element_text(
        size = 8.0
      )
    )

  combined_plot <-
    p_weight +
    p_response_effect +
    p_response_p +
    patchwork::plot_layout(
      widths = c(
        6.2,
        1.15,
        1.05
      ),
      guides = "collect"
    ) &
    ggplot2::theme(
      legend.position = "bottom"
    )

  if (isTRUE(add_annotation)) {
    combined_plot <- combined_plot +
      patchwork::plot_annotation(
        title = paste0(
          factor_name,
          ": active-view feature program"
        ),
        subtitle = paste0(
          "Only views with R² >= ",
          scales::percent(
            mofa_active_view_r2,
            accuracy = 1
          ),
          " are displayed. Weight bars are normalized within each view; response tiles use baseline pCR versus non-pCR contrasts."
        )
      )
  }

  list(
    plot = combined_plot,
    data = factor_weight_data,
    active_views = active_views
  )
}

p_mofa_top_weights <- list()
mofa_factor_feature_plot_data <- list()

for (factor_name in mofa_factor_summary$factor) {
  factor_plot_result <- build_factor_feature_plot(
    factor_name
  )

  p_mofa_top_weights[[factor_name]] <-
    factor_plot_result$plot
  mofa_factor_feature_plot_data[[factor_name]] <-
    factor_plot_result$data

  figure_height <- min(
    13.5,
    max(
      5.2,
      2.3 +
        0.235 *
        nrow(
          factor_plot_result$data
        )
    )
  )

  if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
    filename = paste0(
      "figures/mofa/internal_v16/MOFA_v16_",
      factor_name,
      "_weights_response_tiles.svg"
    ),
    plot = factor_plot_result$plot,
    width = 12.4,
    height = figure_height,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}
}

selected_factor_feature_plot <- build_factor_feature_plot(
  mofa_selected_factor,
  top_features_per_direction = 3,
  add_annotation = FALSE
)

selected_factor_r2_data <-
  mofa_variance_explained %>%
  dplyr::filter(
    as.character(factor) ==
      mofa_selected_factor
  )

p_selected_factor_r2 <- ggplot2::ggplot(
  selected_factor_r2_data,
  ggplot2::aes(
    x = view_label,
    y = r2,
    fill = view_label
  )
) +
  ggplot2::geom_col(
    width = 0.62,
    color = "grey25",
    linewidth = 0.35
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = scales::percent(
        r2,
        accuracy = 0.1
      )
    ),
    vjust = -0.28,
    size = 2.8
  ) +
  ggplot2::scale_fill_manual(
    values = view_colors,
    guide = "none"
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_percent(
      accuracy = 1
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.12)
    )
  ) +
  ggplot2::labs(
    title = paste0(
      mofa_selected_factor,
      " omics fingerprint"
    ),
    subtitle = mofa_selected_factor_reason,
    x = NULL,
    y = "Variance explained"
  ) +
  ggplot2::theme_classic(
    base_size = 9.5
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 7.8,
      color = "grey35"
    ),
    axis.text.x = ggplot2::element_text(
      angle = 25,
      hjust = 1,
      face = "plain"
    )
  )

selected_factor_clinical_data <-
  mofa_primary_clinical_plot_data %>%
  dplyr::filter(
    as.character(factor) ==
      mofa_selected_factor
  )

p_selected_factor_clinical <- ggplot2::ggplot(
  selected_factor_clinical_data,
  ggplot2::aes(
    x = estimate,
    y = term_label
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    color = "grey70",
    linewidth = 0.42
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin = ci_low,
      xmax = ci_high
    ),
    height = 0.12,
    linewidth = 0.62
  ) +
  ggplot2::geom_point(
    shape = 21,
    size = 2.7,
    fill = "white"
  ) +
  ggplot2::labs(
    title = "Clinical association",
    x = "Factor-score coefficient (95% CI)",
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 9.3
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    axis.text.y = ggplot2::element_text(
      size = 8.0
    )
  )

p_mofa_selected_factor_summary <-
  (
    (
      p_selected_factor_r2 /
      p_selected_factor_clinical +
      patchwork::plot_layout(
        heights = c(
          1,
          1.15
        )
      )
    ) |
    selected_factor_feature_plot$plot
  ) +
  patchwork::plot_layout(
    widths = c(
      1.0,
      2.65
    )
  ) +
  patchwork::plot_annotation(
    title = paste0(
      "Selected factor biological summary: ",
      mofa_selected_factor
    ),
    subtitle = "Selection is exploratory and should be reported together with the prespecified mixed-model and permutation results."
  )

#-----------------------------------------------------------------#
# Figure exports
#-----------------------------------------------------------------#

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename = "figures/mofa/internal_v16/MOFA_v16_feature_selection_diagnostics.svg",
  plot = p_mofa_feature_selection,
  width = 11.5,
  height = 3.8,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename = "figures/mofa/internal_v16/MOFA_v16_factor_score_correlations.svg",
  plot = p_mofa_factor_correlation,
  width = 5.8,
  height = 5.2,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename = "figures/mofa/internal_v16/MOFA_v16_data_availability.svg",
  plot = p_mofa_data_availability,
  width = 9.2,
  height = 3.5,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename = "figures/mofa/internal_v16/MOFA_v16_variance_overview.svg",
  plot = p_mofa_variance_overview,
  width = 7.7,
  height = 6.2,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename = "figures/mofa/internal_v16/MOFA_v16_factor_clinical_forest.svg",
  plot = p_mofa_factor_clinical_forest,
  width = 11.5,
  height = 5.0,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename = "figures/mofa/internal_v16/MOFA_v16_factor_landscape.svg",
  plot = p_mofa_factor_landscape,
  width = 12.0,
  height = 5.5,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename = "figures/mofa/internal_v16/MOFA_v16_factor_score_heatmap.svg",
  plot = p_mofa_factor_score_heatmap,
  width = 11.5,
  height = 5.4,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename = "figures/mofa/internal_v16/MOFA_v16_paired_factor_trajectories.svg",
  plot = p_mofa_paired_factor_change,
  width = 12.0,
  height = 6.4,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename = "figures/mofa/internal_v16/MOFA_v16_factor_pair_maps_all.svg",
  plot = p_mofa_factor_pair_maps_all,
  width = 11.2,
  height = 9.0,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename = paste0(
    "figures/mofa/internal_v16/MOFA_v16_selected_factor_summary_",
    mofa_selected_factor,
    ".svg"
  ),
  plot = p_mofa_selected_factor_summary,
  width = 16.0,
  height = 8.0,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}


#=================================================================#
# 13. Longitudinal robustness and subject-specificity diagnostics
#=================================================================#
#
# The primary MOFA fit is intentionally unsupervised and uses all available
# SubjectID-Timepoint observations. This section does not alter that fit. It
# quantifies whether individual factors are dominated by stable inter-individual
# differences and then separates three distinct clinical questions:
#   1. baseline pCR/non-pCR difference;
#   2. overall paired Before-to-Ongoing change;
#   3. differential paired change between pCR and non-pCR subjects.
#
# A random SubjectID intercept corrects dependence among repeated observations;
# it is not a feature supplied to MOFA and therefore cannot make the latent
# model learn SubjectID labels.
#-----------------------------------------------------------------#

mofa_factor_names <- colnames(
  mofa_factor_matrix
)
mofa_factor_names <- mofa_factor_names[
  !is.na(mofa_factor_names) &
    nzchar(mofa_factor_names) &
    grepl(
      "^Factor[0-9]+$",
      mofa_factor_names
    )
]

calculate_factor_subject_specificity <- function(
    factor_name
) {
  factor_data <- mofa_factor_scores %>%
    dplyr::transmute(
      SampleID,
      SubjectID,
      Timepoint = as.character(Timepoint),
      factor_score = .data[[factor_name]]
    ) %>%
    dplyr::filter(
      is.finite(factor_score),
      !is.na(SubjectID),
      Timepoint %in% c(
        "Before",
        "Ongoing"
      )
    )

  repeated_subjects <- factor_data %>%
    dplyr::count(
      SubjectID,
      name = "n_observations"
    ) %>%
    dplyr::filter(
      n_observations >= 2
    ) %>%
    dplyr::pull(SubjectID)

  repeated_data <- factor_data %>%
    dplyr::filter(
      SubjectID %in% repeated_subjects
    )

  random_intercept_fit <- if (
    length(repeated_subjects) >= 5
  ) {
    tryCatch(
      nlme::lme(
        fixed = factor_score ~ 1,
        random = ~ 1 | SubjectID,
        data = repeated_data,
        method = "REML",
        na.action = stats::na.omit,
        control = nlme::lmeControl(
          returnObject = TRUE,
          maxIter = 100,
          msMaxIter = 100
        )
      ),
      error = function(e) {
        NULL
      }
    )
  } else {
    NULL
  }

  subject_variance <- NA_real_
  residual_variance <- NA_real_
  icc <- NA_real_

  if (!is.null(random_intercept_fit)) {
    variance_values <- suppressWarnings(
      as.numeric(
        nlme::VarCorr(
          random_intercept_fit
        )[, "Variance"]
      )
    )

    variance_values <- variance_values[
      is.finite(variance_values)
    ]

    if (length(variance_values) >= 2) {
      subject_variance <- variance_values[1]
      residual_variance <- utils::tail(
        variance_values,
        1
      )
      icc <- subject_variance /
        (
          subject_variance +
            residual_variance
        )
    }
  }

  paired_data <- repeated_data %>%
    dplyr::select(
      SubjectID,
      Timepoint,
      factor_score
    ) %>%
    dplyr::distinct(
      SubjectID,
      Timepoint,
      .keep_all = TRUE
    ) %>%
    tidyr::pivot_wider(
      names_from = Timepoint,
      values_from = factor_score
    ) %>%
    dplyr::filter(
      is.finite(Before),
      is.finite(Ongoing)
    )

  paired_correlation <- if (
    nrow(paired_data) >= 5
  ) {
    suppressWarnings(
      stats::cor(
        paired_data$Before,
        paired_data$Ongoing,
        method = "spearman"
      )
    )
  } else {
    NA_real_
  }

  mean_absolute_change <- if (
    nrow(paired_data) > 0
  ) {
    mean(
      abs(
        paired_data$Ongoing -
          paired_data$Before
      ),
      na.rm = TRUE
    )
  } else {
    NA_real_
  }

  nearest_same_subject_rate <- NA_real_

  if (nrow(repeated_data) >= 6) {
    nearest_same_subject <- logical(
      nrow(repeated_data)
    )

    for (sample_index in seq_len(nrow(repeated_data))) {
      score_distance <- abs(
        repeated_data$factor_score -
          repeated_data$factor_score[sample_index]
      )
      score_distance[sample_index] <- Inf
      nearest_index <- which.min(
        score_distance
      )
      nearest_same_subject[sample_index] <-
        repeated_data$SubjectID[nearest_index] ==
        repeated_data$SubjectID[sample_index]
    }

    nearest_same_subject_rate <- mean(
      nearest_same_subject,
      na.rm = TRUE
    )
  }

  data.frame(
    factor = factor_name,
    n_samples = nrow(factor_data),
    n_repeated_subjects = length(
      repeated_subjects
    ),
    subject_variance = subject_variance,
    residual_variance = residual_variance,
    icc = icc,
    paired_spearman = paired_correlation,
    mean_absolute_paired_change =
      mean_absolute_change,
    nearest_same_subject_rate =
      nearest_same_subject_rate,
    subject_dominated =
      is.finite(icc) &&
      icc >=
        mofa_subject_specificity_icc_threshold,
    stringsAsFactors = FALSE
  )
}

mofa_factor_subject_specificity <-
  dplyr::bind_rows(
    lapply(
      mofa_factor_names,
      calculate_factor_subject_specificity
    )
  )

p_mofa_subject_specificity <- ggplot2::ggplot(
  mofa_factor_subject_specificity,
  ggplot2::aes(
    x = factor,
    y = icc,
    fill = subject_dominated
  )
) +
  ggplot2::geom_hline(
    yintercept =
      mofa_subject_specificity_icc_threshold,
    linetype = "dashed",
    linewidth = 0.55,
    color = "grey45"
  ) +
  ggplot2::geom_col(
    width = 0.68,
    color = "grey25",
    linewidth = 0.35
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = ifelse(
        is.finite(icc),
        sprintf(
          "%.2f",
          icc
        ),
        "NA"
      )
    ),
    vjust = -0.35,
    size = 3.0
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      `FALSE` = "grey70",
      `TRUE` = "grey25"
    ),
    guide = "none"
  ) +
  ggplot2::scale_y_continuous(
    limits = c(
      0,
      1.05
    ),
    breaks = seq(
      0,
      1,
      by = 0.25
    ),
    expand = ggplot2::expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  ggplot2::labs(
    title = "Subject specificity of MOFA factors",
    subtitle = paste0(
      "ICC is estimated from a random-intercept model among repeated subjects; ",
      "ICC >= ",
      mofa_subject_specificity_icc_threshold,
      " is flagged as subject-dominated."
    ),
    x = NULL,
    y = "Subject-level intraclass correlation"
  ) +
  ggplot2::theme_classic(
    base_size = 10.5
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 8.8,
      color = "grey35"
    ),
    axis.text.x = ggplot2::element_text(
      face = "plain"
    )
  )


#=================================================================#
# 14. Hierarchical clinical analysis of factor scores
#=================================================================#
#
# Primary response question:
#   baseline factor score ~ pCR status
#
# Primary longitudinal questions:
#   paired delta ~ 1                         (overall treatment change)
#   paired delta ~ pCR status                (differential treatment change)
#
# Supportive model when paired sampling is sufficient:
#   factor score ~ Timepoint * pCR + (1 | SubjectID)
#
# The interaction term tests whether Before-to-Ongoing change differs between
# response groups. It is not needed for a baseline-response comparison and is
# fitted only when both response groups contain enough paired subjects.
#-----------------------------------------------------------------#

extract_standard_model_row <- function(
    coefficient_table,
    coefficient_name,
    factor_name,
    model_type,
    interpretation,
    n_samples,
    n_subjects,
    n_paired_subjects = NA_integer_
) {
  if (
    is.null(coefficient_table) ||
    !coefficient_name %in%
      rownames(coefficient_table)
  ) {
    return(
      data.frame()
    )
  }

  coefficient_row <- coefficient_table[
    coefficient_name,
    ,
    drop = FALSE
  ]

  estimate <- as.numeric(
    coefficient_row[1, "Estimate"]
  )

  standard_error <- as.numeric(
    coefficient_row[1, "Std. Error"]
  )

  statistic_name <- intersect(
    c(
      "t value",
      "z value"
    ),
    colnames(coefficient_row)
  )

  p_name <- intersect(
    c(
      "Pr(>|t|)",
      "Pr(>|z|)"
    ),
    colnames(coefficient_row)
  )

  statistic <- if (length(statistic_name) > 0) {
    as.numeric(
      coefficient_row[1, statistic_name[1]]
    )
  } else {
    estimate / standard_error
  }

  model_p <- if (length(p_name) > 0) {
    as.numeric(
      coefficient_row[1, p_name[1]]
    )
  } else {
    2 * stats::pnorm(
      -abs(statistic)
    )
  }

  data.frame(
    factor = factor_name,
    model_type = model_type,
    term = coefficient_name,
    interpretation = interpretation,
    estimate = estimate,
    std_error = standard_error,
    statistic = statistic,
    model_p = model_p,
    ci_low = estimate -
      stats::qnorm(0.975) *
      standard_error,
    ci_high = estimate +
      stats::qnorm(0.975) *
      standard_error,
    n_samples = n_samples,
    n_subjects = n_subjects,
    n_paired_subjects =
      n_paired_subjects,
    stringsAsFactors = FALSE
  )
}

extract_lme_model_row <- function(
    model_fit,
    coefficient_name,
    factor_name,
    model_type,
    interpretation,
    n_samples,
    n_subjects,
    n_paired_subjects
) {
  if (is.null(model_fit)) {
    return(
      data.frame()
    )
  }

  coefficient_table <- as.data.frame(
    summary(model_fit)$tTable,
    check.names = FALSE
  )

  if (!coefficient_name %in%
      rownames(coefficient_table)) {
    return(
      data.frame()
    )
  }

  coefficient_row <- coefficient_table[
    coefficient_name,
    ,
    drop = FALSE
  ]

  degrees_freedom <- as.numeric(
    coefficient_row[1, "DF"]
  )

  estimate <- as.numeric(
    coefficient_row[1, "Value"]
  )
  standard_error <- as.numeric(
    coefficient_row[1, "Std.Error"]
  )

  data.frame(
    factor = factor_name,
    model_type = model_type,
    term = coefficient_name,
    interpretation = interpretation,
    estimate = estimate,
    std_error = standard_error,
    statistic = as.numeric(
      coefficient_row[1, "t-value"]
    ),
    model_p = as.numeric(
      coefficient_row[1, "p-value"]
    ),
    ci_low = estimate -
      stats::qt(
        0.975,
        degrees_freedom
      ) *
      standard_error,
    ci_high = estimate +
      stats::qt(
        0.975,
        degrees_freedom
      ) *
      standard_error,
    n_samples = n_samples,
    n_subjects = n_subjects,
    n_paired_subjects =
      n_paired_subjects,
    stringsAsFactors = FALSE
  )
}

permute_group_effect <- function(
    model_data,
    outcome_name,
    include_n_views = TRUE,
    permutations = mofa_factor_permutations
) {
  if (
    nrow(model_data) < 6 ||
    length(
      unique(
        model_data$TRG_plot
      )
    ) < 2
  ) {
    return(
      list(
        observed = NA_real_,
        permutation_p = NA_real_,
        valid_permutations = 0L
      )
    )
  }

  permutation_data <- model_data
  permutation_data$outcome <-
    permutation_data[[outcome_name]]
  permutation_data$TRG_plot <- stats::relevel(
    factor(
      permutation_data$TRG_plot,
      levels = c(
        "non_pCR",
        "pCR"
      )
    ),
    ref = "non_pCR"
  )

  if (include_n_views) {
    permutation_data$n_views_centered <-
      permutation_data$n_views -
      mean(
        permutation_data$n_views,
        na.rm = TRUE
      )

    observed_fit <- tryCatch(
      stats::lm(
        outcome ~
          TRG_plot +
          n_views_centered,
        data = permutation_data
      ),
      error = function(e) {
        NULL
      }
    )
  } else {
    observed_fit <- tryCatch(
      stats::lm(
        outcome ~ TRG_plot,
        data = permutation_data
      ),
      error = function(e) {
        NULL
      }
    )
  }

  if (is.null(observed_fit)) {
    return(
      list(
        observed = NA_real_,
        permutation_p = NA_real_,
        valid_permutations = 0L
      )
    )
  }

  observed_estimate <- stats::coef(
    observed_fit
  )[
    "TRG_plotpCR"
  ]

  permuted_estimates <- rep(
    NA_real_,
    permutations
  )

  for (permutation_index in seq_len(permutations)) {
    permuted_data <- permutation_data
    permuted_data$TRG_plot <- sample(
      permutation_data$TRG_plot,
      size = nrow(permutation_data),
      replace = FALSE
    )

    permuted_fit <- tryCatch(
      if (include_n_views) {
        stats::lm(
          outcome ~
            TRG_plot +
            n_views_centered,
          data = permuted_data
        )
      } else {
        stats::lm(
          outcome ~ TRG_plot,
          data = permuted_data
        )
      },
      error = function(e) {
        NULL
      }
    )

    if (!is.null(permuted_fit)) {
      permuted_estimates[permutation_index] <-
        stats::coef(
          permuted_fit
        )[
          "TRG_plotpCR"
        ]
    }
  }

  permuted_estimates <- permuted_estimates[
    is.finite(permuted_estimates)
  ]

  list(
    observed = observed_estimate,
    permutation_p = if (
      length(permuted_estimates) > 0
    ) {
      (
        1 +
          sum(
            abs(permuted_estimates) >=
              abs(observed_estimate)
          )
      ) /
      (
        1 +
          length(permuted_estimates)
      )
    } else {
      NA_real_
    },
    valid_permutations =
      length(permuted_estimates)
  )
}

sign_flip_time_effect <- function(
    paired_delta,
    permutations = mofa_factor_permutations
) {
  paired_delta <- paired_delta[
    is.finite(paired_delta)
  ]

  if (length(paired_delta) < 4) {
    return(
      list(
        observed = NA_real_,
        permutation_p = NA_real_,
        valid_permutations = 0L
      )
    )
  }

  observed_mean <- mean(
    paired_delta
  )

  permuted_means <- replicate(
    permutations,
    mean(
      paired_delta *
        sample(
          c(
            -1,
            1
          ),
          size = length(paired_delta),
          replace = TRUE
        )
    )
  )

  list(
    observed = observed_mean,
    permutation_p =
      (
        1 +
          sum(
            abs(permuted_means) >=
              abs(observed_mean)
          )
      ) /
      (
        1 +
          length(permuted_means)
      ),
    valid_permutations =
      length(permuted_means)
  )
}

mofa_factor_hierarchical_models <- data.frame()
mofa_factor_hierarchical_permutations <- data.frame()
mofa_factor_paired_data <- list()

set.seed(20260901)

for (factor_name in mofa_factor_names) {
  factor_data <- mofa_factor_scores %>%
    dplyr::transmute(
      SampleID,
      SubjectID,
      Timepoint = factor(
        as.character(Timepoint),
        levels = c(
          "Before",
          "Ongoing"
        )
      ),
      TRG_plot = factor(
        as.character(TRG_plot),
        levels = c(
          "non_pCR",
          "pCR"
        )
      ),
      n_views,
      factor_score = .data[[factor_name]]
    ) %>%
    dplyr::filter(
      is.finite(factor_score),
      !is.na(SubjectID),
      !is.na(Timepoint),
      !is.na(TRG_plot)
    )

  baseline_data <- factor_data %>%
    dplyr::filter(
      Timepoint == "Before"
    ) %>%
    dplyr::distinct(
      SubjectID,
      .keep_all = TRUE
    ) %>%
    dplyr::mutate(
      n_views_centered = n_views -
        mean(
          n_views,
          na.rm = TRUE
        )
    )

  ongoing_data <- factor_data %>%
    dplyr::filter(
      Timepoint == "Ongoing"
    ) %>%
    dplyr::distinct(
      SubjectID,
      .keep_all = TRUE
    ) %>%
    dplyr::mutate(
      n_views_centered = n_views -
        mean(
          n_views,
          na.rm = TRUE
        )
    )

  before_for_pair <- factor_data %>%
    dplyr::filter(
      Timepoint == "Before"
    ) %>%
    dplyr::transmute(
      SubjectID,
      TRG_plot,
      factor_before = factor_score,
      n_views_before = n_views
    ) %>%
    dplyr::distinct(
      SubjectID,
      .keep_all = TRUE
    )

  ongoing_for_pair <- factor_data %>%
    dplyr::filter(
      Timepoint == "Ongoing"
    ) %>%
    dplyr::transmute(
      SubjectID,
      factor_ongoing = factor_score,
      n_views_ongoing = n_views
    ) %>%
    dplyr::distinct(
      SubjectID,
      .keep_all = TRUE
    )

  paired_data <- dplyr::inner_join(
    before_for_pair,
    ongoing_for_pair,
    by = "SubjectID"
  ) %>%
    dplyr::mutate(
      factor_delta = factor_ongoing -
        factor_before,
      n_views = pmin(
        n_views_before,
        n_views_ongoing
      )
    )

  mofa_factor_paired_data[[factor_name]] <-
    paired_data

  if (
    nrow(baseline_data) >= 8 &&
    min(
      table(
        baseline_data$TRG_plot
      )
    ) >= 3
  ) {
    baseline_fit <- stats::lm(
      factor_score ~
        TRG_plot +
        n_views_centered,
      data = baseline_data
    )

    baseline_row <- extract_standard_model_row(
      summary(baseline_fit)$coefficients,
      "TRG_plotpCR",
      factor_name,
      "baseline_response",
      "pCR minus non-pCR at Before",
      nrow(baseline_data),
      dplyr::n_distinct(
        baseline_data$SubjectID
      )
    )

    mofa_factor_hierarchical_models <-
      dplyr::bind_rows(
        mofa_factor_hierarchical_models,
        baseline_row
      )

    baseline_permutation <-
      permute_group_effect(
        baseline_data,
        outcome_name = "factor_score",
        include_n_views = TRUE
      )

    mofa_factor_hierarchical_permutations <-
      dplyr::bind_rows(
        mofa_factor_hierarchical_permutations,
        data.frame(
          factor = factor_name,
          model_type = "baseline_response",
          term = "TRG_plotpCR",
          observed_estimate =
            baseline_permutation$observed,
          permutation_p =
            baseline_permutation$permutation_p,
          valid_permutations =
            baseline_permutation$valid_permutations,
          stringsAsFactors = FALSE
        )
      )
  }

  if (
    nrow(ongoing_data) >= 8 &&
    min(
      table(
        ongoing_data$TRG_plot
      )
    ) >= 3
  ) {
    ongoing_fit <- stats::lm(
      factor_score ~
        TRG_plot +
        n_views_centered,
      data = ongoing_data
    )

    ongoing_row <- extract_standard_model_row(
      summary(ongoing_fit)$coefficients,
      "TRG_plotpCR",
      factor_name,
      "ongoing_response",
      "pCR minus non-pCR at Ongoing",
      nrow(ongoing_data),
      dplyr::n_distinct(
        ongoing_data$SubjectID
      )
    )

    mofa_factor_hierarchical_models <-
      dplyr::bind_rows(
        mofa_factor_hierarchical_models,
        ongoing_row
      )

    ongoing_permutation <-
      permute_group_effect(
        ongoing_data,
        outcome_name = "factor_score",
        include_n_views = TRUE
      )

    mofa_factor_hierarchical_permutations <-
      dplyr::bind_rows(
        mofa_factor_hierarchical_permutations,
        data.frame(
          factor = factor_name,
          model_type = "ongoing_response",
          term = "TRG_plotpCR",
          observed_estimate =
            ongoing_permutation$observed,
          permutation_p =
            ongoing_permutation$permutation_p,
          valid_permutations =
            ongoing_permutation$valid_permutations,
          stringsAsFactors = FALSE
        )
      )
  }

  if (
    nrow(paired_data) >=
      mofa_min_paired_subjects
  ) {
    paired_time_fit <- stats::lm(
      factor_delta ~ 1,
      data = paired_data
    )

    paired_time_row <- extract_standard_model_row(
      summary(paired_time_fit)$coefficients,
      "(Intercept)",
      factor_name,
      "paired_overall_change",
      "Mean Ongoing minus Before change",
      nrow(paired_data) * 2,
      nrow(paired_data),
      nrow(paired_data)
    )

    mofa_factor_hierarchical_models <-
      dplyr::bind_rows(
        mofa_factor_hierarchical_models,
        paired_time_row
      )

    time_permutation <- sign_flip_time_effect(
      paired_data$factor_delta
    )

    mofa_factor_hierarchical_permutations <-
      dplyr::bind_rows(
        mofa_factor_hierarchical_permutations,
        data.frame(
          factor = factor_name,
          model_type = "paired_overall_change",
          term = "(Intercept)",
          observed_estimate =
            time_permutation$observed,
          permutation_p =
            time_permutation$permutation_p,
          valid_permutations =
            time_permutation$valid_permutations,
          stringsAsFactors = FALSE
        )
      )
  }

  paired_response_counts <- table(
    paired_data$TRG_plot
  )

  if (
    nrow(paired_data) >=
      mofa_min_paired_subjects &&
    length(paired_response_counts) == 2 &&
    min(paired_response_counts) >=
      mofa_min_paired_per_response
  ) {
    paired_response_fit <- stats::lm(
      factor_delta ~ TRG_plot,
      data = paired_data
    )

    paired_response_row <-
      extract_standard_model_row(
        summary(
          paired_response_fit
        )$coefficients,
        "TRG_plotpCR",
        factor_name,
        "paired_differential_change",
        paste0(
          "Difference in Ongoing-Before change: ",
          "pCR minus non-pCR"
        ),
        nrow(paired_data) * 2,
        nrow(paired_data),
        nrow(paired_data)
      )

    mofa_factor_hierarchical_models <-
      dplyr::bind_rows(
        mofa_factor_hierarchical_models,
        paired_response_row
      )

    paired_permutation <-
      permute_group_effect(
        paired_data,
        outcome_name = "factor_delta",
        include_n_views = FALSE
      )

    mofa_factor_hierarchical_permutations <-
      dplyr::bind_rows(
        mofa_factor_hierarchical_permutations,
        data.frame(
          factor = factor_name,
          model_type =
            "paired_differential_change",
          term = "TRG_plotpCR",
          observed_estimate =
            paired_permutation$observed,
          permutation_p =
            paired_permutation$permutation_p,
          valid_permutations =
            paired_permutation$valid_permutations,
          stringsAsFactors = FALSE
        )
      )
  }

  n_paired_by_response <- paired_data %>%
    dplyr::count(
      TRG_plot,
      name = "n"
    )

  interaction_supported <-
    nrow(paired_data) >=
      mofa_min_paired_subjects &&
    nrow(n_paired_by_response) == 2 &&
    min(n_paired_by_response$n) >=
      mofa_min_paired_per_response

  if (
    dplyr::n_distinct(
      factor_data$SubjectID
    ) >= 8
  ) {
    factor_data$n_views_centered <-
      factor_data$n_views -
      mean(
        factor_data$n_views,
        na.rm = TRUE
      )

    additive_fit <- tryCatch(
      nlme::lme(
        fixed = factor_score ~
          Timepoint +
          TRG_plot +
          n_views_centered,
        random = ~ 1 | SubjectID,
        data = factor_data,
        method = "REML",
        na.action = stats::na.omit,
        control = nlme::lmeControl(
          returnObject = TRUE,
          maxIter = 100,
          msMaxIter = 100
        )
      ),
      error = function(e) {
        NULL
      }
    )

    mofa_factor_hierarchical_models <-
      dplyr::bind_rows(
        mofa_factor_hierarchical_models,
        extract_lme_model_row(
          additive_fit,
          "TimepointOngoing",
          factor_name,
          "mixed_additive",
          paste0(
            "Adjusted Ongoing-Before difference ",
            "without interaction"
          ),
          nrow(factor_data),
          dplyr::n_distinct(
            factor_data$SubjectID
          ),
          nrow(paired_data)
        ),
        extract_lme_model_row(
          additive_fit,
          "TRG_plotpCR",
          factor_name,
          "mixed_additive",
          paste0(
            "Adjusted pCR-non-pCR difference ",
            "without interaction"
          ),
          nrow(factor_data),
          dplyr::n_distinct(
            factor_data$SubjectID
          ),
          nrow(paired_data)
        )
      )

    if (interaction_supported) {
      interaction_fit <- tryCatch(
        nlme::lme(
          fixed = factor_score ~
            Timepoint * TRG_plot +
            n_views_centered,
          random = ~ 1 | SubjectID,
          data = factor_data,
          method = "REML",
          na.action = stats::na.omit,
          control = nlme::lmeControl(
            returnObject = TRUE,
            maxIter = 100,
            msMaxIter = 100
          )
        ),
        error = function(e) {
          NULL
        }
      )

      mofa_factor_hierarchical_models <-
        dplyr::bind_rows(
          mofa_factor_hierarchical_models,
          extract_lme_model_row(
            interaction_fit,
            "TimepointOngoing:TRG_plotpCR",
            factor_name,
            "mixed_interaction",
            paste0(
              "Differential time change: ",
              "pCR minus non-pCR"
            ),
            nrow(factor_data),
            dplyr::n_distinct(
              factor_data$SubjectID
            ),
            nrow(paired_data)
          )
        )
    }
  }
}

mofa_factor_hierarchical_models <-
  mofa_factor_hierarchical_models %>%
  dplyr::left_join(
    mofa_factor_hierarchical_permutations,
    by = c(
      "factor",
      "model_type",
      "term"
    )
  ) %>%
  dplyr::group_by(
    model_type
  ) %>%
  dplyr::mutate(
    model_fdr = stats::p.adjust(
      model_p,
      method = "BH"
    ),
    permutation_fdr = stats::p.adjust(
      permutation_p,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup()

safe_wilcoxon_rank_sum <- function(
    x,
    group
) {
  keep <- is.finite(x) &
    !is.na(group)
  x <- x[keep]
  group <- droplevels(
    factor(
      as.character(group[keep]),
      levels = c(
        "non_pCR",
        "pCR"
      )
    )
  )

  if (
    length(x) < 6 ||
    nlevels(group) != 2 ||
    min(table(group)) < 3
  ) {
    return(NA_real_)
  }

  suppressWarnings(
    tryCatch(
      stats::wilcox.test(
        x ~ group,
        exact = FALSE,
        correct = FALSE
      )$p.value,
      error = function(e) {
        NA_real_
      }
    )
  )
}

safe_wilcoxon_paired <- function(
    before,
    after
) {
  keep <- is.finite(before) &
    is.finite(after)
  before <- before[keep]
  after <- after[keep]

  if (length(before) < 4) {
    return(NA_real_)
  }

  suppressWarnings(
    tryCatch(
      stats::wilcox.test(
        after,
        before,
        paired = TRUE,
        exact = FALSE,
        correct = FALSE
      )$p.value,
      error = function(e) {
        NA_real_
      }
    )
  )
}

mofa_factor_wilcoxon_tests <- dplyr::bind_rows(
  lapply(
    mofa_factor_names,
    function(factor_name) {
      factor_data <- mofa_factor_scores %>%
        dplyr::transmute(
          SubjectID,
          Timepoint = as.character(Timepoint),
          TRG_plot = factor(
            as.character(TRG_plot),
            levels = c(
              "non_pCR",
              "pCR"
            )
          ),
          factor_score = .data[[factor_name]]
        ) %>%
        dplyr::filter(
          is.finite(factor_score),
          !is.na(SubjectID),
          !is.na(TRG_plot),
          Timepoint %in% c(
            "Before",
            "Ongoing"
          )
        )

      baseline_data <- factor_data %>%
        dplyr::filter(
          Timepoint == "Before"
        ) %>%
        dplyr::distinct(
          SubjectID,
          .keep_all = TRUE
        )

      ongoing_data <- factor_data %>%
        dplyr::filter(
          Timepoint == "Ongoing"
        ) %>%
        dplyr::distinct(
          SubjectID,
          .keep_all = TRUE
        )

      subject_mean_data <- factor_data %>%
        dplyr::group_by(
          SubjectID,
          TRG_plot
        ) %>%
        dplyr::summarise(
          factor_score = mean(
            factor_score,
            na.rm = TRUE
          ),
          n_timepoints = dplyr::n(),
          .groups = "drop"
        )

      paired_data <- dplyr::inner_join(
        baseline_data %>%
          dplyr::transmute(
            SubjectID,
            TRG_plot,
            factor_before = factor_score
          ),
        ongoing_data %>%
          dplyr::transmute(
            SubjectID,
            factor_after = factor_score
          ),
        by = "SubjectID"
      ) %>%
        dplyr::mutate(
          factor_delta = factor_after -
            factor_before
        )

      data.frame(
        factor = factor_name,
        comparison = c(
          "overall_response_subject_mean",
          "baseline_response",
          "ongoing_response",
          "paired_overall_change",
          "paired_differential_change"
        ),
        comparison_label = c(
          "pCR vs non-pCR; subject mean across available timepoints",
          "pCR vs non-pCR at Baseline",
          "pCR vs non-pCR after RT",
          "Within-subject After RT minus Baseline change",
          "pCR vs non-pCR difference in paired change"
        ),
        wilcoxon_effect = c(
          median(
            subject_mean_data$factor_score[
              subject_mean_data$TRG_plot == "pCR"
            ],
            na.rm = TRUE
          ) -
            median(
              subject_mean_data$factor_score[
                subject_mean_data$TRG_plot == "non_pCR"
              ],
              na.rm = TRUE
            ),
          median(
            baseline_data$factor_score[
              baseline_data$TRG_plot == "pCR"
            ],
            na.rm = TRUE
          ) -
            median(
              baseline_data$factor_score[
                baseline_data$TRG_plot == "non_pCR"
              ],
              na.rm = TRUE
            ),
          median(
            ongoing_data$factor_score[
              ongoing_data$TRG_plot == "pCR"
            ],
            na.rm = TRUE
          ) -
            median(
              ongoing_data$factor_score[
                ongoing_data$TRG_plot == "non_pCR"
              ],
              na.rm = TRUE
            ),
          median(
            paired_data$factor_delta,
            na.rm = TRUE
          ),
          median(
            paired_data$factor_delta[
              paired_data$TRG_plot == "pCR"
            ],
            na.rm = TRUE
          ) -
            median(
              paired_data$factor_delta[
                paired_data$TRG_plot == "non_pCR"
              ],
              na.rm = TRUE
            )
        ),
        wilcoxon_p = c(
          safe_wilcoxon_rank_sum(
            subject_mean_data$factor_score,
            subject_mean_data$TRG_plot
          ),
          safe_wilcoxon_rank_sum(
            baseline_data$factor_score,
            baseline_data$TRG_plot
          ),
          safe_wilcoxon_rank_sum(
            ongoing_data$factor_score,
            ongoing_data$TRG_plot
          ),
          safe_wilcoxon_paired(
            paired_data$factor_before,
            paired_data$factor_after
          ),
          safe_wilcoxon_rank_sum(
            paired_data$factor_delta,
            paired_data$TRG_plot
          )
        ),
        n_subjects = c(
          dplyr::n_distinct(
            subject_mean_data$SubjectID
          ),
          dplyr::n_distinct(
            baseline_data$SubjectID
          ),
          dplyr::n_distinct(
            ongoing_data$SubjectID
          ),
          nrow(paired_data),
          nrow(paired_data)
        ),
        n_pCR = c(
          sum(
            subject_mean_data$TRG_plot == "pCR"
          ),
          sum(
            baseline_data$TRG_plot == "pCR"
          ),
          sum(
            ongoing_data$TRG_plot == "pCR"
          ),
          sum(
            paired_data$TRG_plot == "pCR"
          ),
          sum(
            paired_data$TRG_plot == "pCR"
          )
        ),
        n_non_pCR = c(
          sum(
            subject_mean_data$TRG_plot == "non_pCR"
          ),
          sum(
            baseline_data$TRG_plot == "non_pCR"
          ),
          sum(
            ongoing_data$TRG_plot == "non_pCR"
          ),
          sum(
            paired_data$TRG_plot == "non_pCR"
          ),
          sum(
            paired_data$TRG_plot == "non_pCR"
          )
        ),
        test = c(
          "Wilcoxon rank-sum; one subject-level mean per subject",
          "Wilcoxon rank-sum",
          "Wilcoxon rank-sum",
          "Wilcoxon signed-rank; paired by SubjectID",
          "Wilcoxon rank-sum on SubjectID-level deltas"
        ),
        stringsAsFactors = FALSE
      )
    }
  )
) %>%
  dplyr::group_by(
    comparison
  ) %>%
  dplyr::mutate(
    wilcoxon_fdr = stats::p.adjust(
      wilcoxon_p,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup()

mofa_factor_hierarchical_models <-
  mofa_factor_hierarchical_models %>%
  dplyr::left_join(
    mofa_factor_wilcoxon_tests %>%
      dplyr::select(
        factor,
        comparison,
        wilcoxon_effect,
        wilcoxon_p,
        wilcoxon_fdr,
        wilcoxon_test = test
      ),
    by = c(
      "factor" = "factor",
      "model_type" = "comparison"
    )
  )

write.csv(
  mofa_factor_wilcoxon_tests,
  "results/mofa/MOFA_v30_factor_wilcoxon_tests.csv",
  row.names = FALSE
)

write.csv(
  mofa_factor_hierarchical_models,
  "results/mofa/MOFA_v30_factor_hierarchical_models.csv",
  row.names = FALSE
)

hierarchical_plot_data <-
  mofa_factor_hierarchical_models %>%
  dplyr::filter(
    model_type %in% c(
      "baseline_response",
      "paired_overall_change",
      "paired_differential_change"
    )
  ) %>%
  dplyr::mutate(
    model_label = factor(
      dplyr::recode(
        model_type,
        baseline_response =
          "Baseline: pCR - non-pCR",
        paired_overall_change =
          "Paired change: Ongoing - Before",
        paired_differential_change =
          paste0(
            "Differential paired change: ",
            "pCR - non-pCR"
          )
      ),
      levels = c(
        "Baseline: pCR - non-pCR",
        "Paired change: Ongoing - Before",
        paste0(
          "Differential paired change: ",
          "pCR - non-pCR"
        )
      )
    ),
    factor = factor(
      factor,
      levels = rev(
        as.character(
          mofa_factor_summary$factor
        )
      )
    )
  )

p_mofa_hierarchical_forest <- ggplot2::ggplot(
  hierarchical_plot_data,
  ggplot2::aes(
    x = estimate,
    y = factor
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.48,
    color = "grey55"
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin = ci_low,
      xmax = ci_high
    ),
    height = 0.12,
    linewidth = 0.65
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      size = -log10(
        pmax(
          permutation_p,
          .Machine$double.xmin
        )
      )
    ),
    shape = 21,
    fill = "white",
    color = "grey15",
    stroke = 0.55
  ) +
  ggplot2::facet_grid(
    . ~ model_label,
    scales = "free_x"
  ) +
  ggplot2::scale_size_continuous(
    range = c(
      2.2,
      5.0
    ),
    name = "-log10 permutation p"
  ) +
  ggplot2::labs(
    title = "MOFA factor associations separated by biological question",
    subtitle = paste0(
      "Baseline response, paired overall change, and differential paired change ",
      "are not interchangeable estimands."
    ),
    x = "Factor-score effect",
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 10.3
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 8.8,
      color = "grey35"
    ),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(
      face = "plain"
    ),
    legend.position = "bottom"
  )


#=================================================================#
# 15. Factor-pair atlas and response-oriented pair selection
#=================================================================#
#
# All factor combinations are evaluated on one baseline sample per subject.
# Pair selection is exploratory and is reported with a max-statistic permutation
# p-value that accounts for searching across all factor combinations.
#-----------------------------------------------------------------#

calculate_auc_rank <- function(
    observed_group,
    predicted_score
) {
  keep_index <-
    !is.na(observed_group) &
    is.finite(predicted_score)

  observed_group <- factor(
    observed_group[keep_index],
    levels = c(
      "non_pCR",
      "pCR"
    )
  )
  predicted_score <- predicted_score[
    keep_index
  ]

  n_positive <- sum(
    observed_group == "pCR"
  )
  n_negative <- sum(
    observed_group == "non_pCR"
  )

  if (
    n_positive == 0 ||
    n_negative == 0
  ) {
    return(NA_real_)
  }

  score_rank <- rank(
    predicted_score,
    ties.method = "average"
  )

  auc <- (
    sum(
      score_rank[
        observed_group == "pCR"
      ]
    ) -
      n_positive *
      (
        n_positive + 1
      ) /
      2
  ) /
  (
    n_positive *
      n_negative
  )

  as.numeric(auc)
}

loocv_pair_auc <- function(
    pair_data,
    factor_x,
    factor_y
) {
  pair_data <- pair_data %>%
    dplyr::filter(
      is.finite(.data[[factor_x]]),
      is.finite(.data[[factor_y]]),
      !is.na(TRG_plot)
    )

  if (
    nrow(pair_data) < 8 ||
    min(
      table(
        pair_data$TRG_plot
      )
    ) < 3
  ) {
    return(NA_real_)
  }

  prediction <- rep(
    NA_real_,
    nrow(pair_data)
  )

  for (test_index in seq_len(nrow(pair_data))) {
    training_data <- pair_data[
      -test_index,
      ,
      drop = FALSE
    ]
    test_data <- pair_data[
      test_index,
      ,
      drop = FALSE
    ]

    training_matrix <- as.matrix(
      training_data[
        ,
        c(
          factor_x,
          factor_y
        ),
        drop = FALSE
      ]
    )

    test_matrix <- as.matrix(
      test_data[
        ,
        c(
          factor_x,
          factor_y
        ),
        drop = FALSE
      ]
    )

    training_mean <- colMeans(
      training_matrix
    )
    training_sd <- apply(
      training_matrix,
      2,
      stats::sd
    )
    training_sd[
      !is.finite(training_sd) |
      training_sd <= 0
    ] <- 1

    training_matrix <- sweep(
      sweep(
        training_matrix,
        2,
        training_mean,
        "-"
      ),
      2,
      training_sd,
      "/"
    )

    test_matrix <- sweep(
      sweep(
        test_matrix,
        2,
        training_mean,
        "-"
      ),
      2,
      training_sd,
      "/"
    )

    lda_fit <- tryCatch(
      MASS::lda(
        x = training_matrix,
        grouping = training_data$TRG_plot
      ),
      error = function(e) {
        NULL
      }
    )

    if (!is.null(lda_fit)) {
      lda_prediction <- MASS:::predict.lda(
        lda_fit,
        newdata = test_matrix
      )$posterior

      if (
        "pCR" %in%
        colnames(lda_prediction)
      ) {
        prediction[test_index] <-
          lda_prediction[1, "pCR"]
      }
    }
  }

  calculate_auc_rank(
    pair_data$TRG_plot,
    prediction
  )
}

standardize_matrix_columns <- function(data_matrix) {
  data_matrix <- as.matrix(data_matrix)
  column_mean <- colMeans(data_matrix, na.rm = TRUE)
  column_sd <- apply(data_matrix, 2, stats::sd, na.rm = TRUE)
  column_sd[!is.finite(column_sd) | column_sd <= 0] <- 1

  sweep(
    sweep(data_matrix, 2, column_mean, FUN = "-"),
    2,
    column_sd,
    FUN = "/"
  )
}

baseline_factor_data <- mofa_factor_scores %>%
  dplyr::filter(
    as.character(Timepoint) == "Before",
    !is.na(TRG_plot)
  ) %>%
  dplyr::distinct(
    SubjectID,
    .keep_all = TRUE
  )

factor_pair_matrix <- utils::combn(
  mofa_factor_names,
  2
)

factor_r2_matrix <- mofa_variance_explained %>%
  dplyr::mutate(
    factor = as.character(factor),
    view = as.character(view)
  ) %>%
  dplyr::select(
    factor,
    view,
    r2
  ) %>%
  tidyr::pivot_wider(
    names_from = view,
    values_from = r2,
    values_fill = 0
  ) %>%
  as.data.frame(
    check.names = FALSE
  )

rownames(factor_r2_matrix) <-
  factor_r2_matrix$factor
factor_r2_matrix$factor <- NULL
factor_r2_matrix <- as.matrix(
  factor_r2_matrix
)

factor_pair_rows <- vector(
  "list",
  ncol(factor_pair_matrix)
)

for (pair_index in seq_len(ncol(factor_pair_matrix))) {
  factor_x <- factor_pair_matrix[1, pair_index]
  factor_y <- factor_pair_matrix[2, pair_index]

  pair_data <- baseline_factor_data %>%
    dplyr::select(
      TRG_plot,
      dplyr::all_of(
        c(
          factor_x,
          factor_y
        )
      )
    ) %>%
    dplyr::filter(
      is.finite(.data[[factor_x]]),
      is.finite(.data[[factor_y]])
    )

  standardized_scores <- standardize_matrix_columns(
    pair_data[
      ,
      c(
        factor_x,
        factor_y
      ),
      drop = FALSE
    ]
  )

  group_difference <- if (
    nrow(pair_data) >= 6 &&
    length(
      unique(
        pair_data$TRG_plot
      )
    ) == 2
  ) {
    colMeans(
      standardized_scores[
        pair_data$TRG_plot == "pCR",
        ,
        drop = FALSE
      ]
    ) -
      colMeans(
        standardized_scores[
          pair_data$TRG_plot == "non_pCR",
          ,
          drop = FALSE
        ]
      )
  } else {
    c(
      NA_real_,
      NA_real_
    )
  }

  r2_x <- factor_r2_matrix[
    factor_x,
    required_views,
    drop = TRUE
  ]
  r2_y <- factor_r2_matrix[
    factor_y,
    required_views,
    drop = TRUE
  ]

  fingerprint_cosine <- if (
    sqrt(sum(r2_x ^ 2)) > 0 &&
    sqrt(sum(r2_y ^ 2)) > 0
  ) {
    sum(r2_x * r2_y) /
      (
        sqrt(sum(r2_x ^ 2)) *
          sqrt(sum(r2_y ^ 2))
      )
  } else {
    NA_real_
  }

  factor_pair_rows[[pair_index]] <-
    data.frame(
      pair_id = paste(
        factor_x,
        factor_y,
        sep = "__"
      ),
      factor_x = factor_x,
      factor_y = factor_y,
      n_baseline = nrow(pair_data),
      score_correlation = suppressWarnings(
        stats::cor(
          pair_data[[factor_x]],
          pair_data[[factor_y]],
          method = "spearman"
        )
      ),
      view_fingerprint_cosine =
        fingerprint_cosine,
      centroid_distance = sqrt(
        sum(
          group_difference ^ 2
        )
      ),
      response_vector_x =
        group_difference[1],
      response_vector_y =
        group_difference[2],
      loocv_lda_auc = loocv_pair_auc(
        baseline_factor_data,
        factor_x,
        factor_y
      ),
      stringsAsFactors = FALSE
    )
}

mofa_factor_pair_summary <- dplyr::bind_rows(
  factor_pair_rows
)

set.seed(20260902)

pair_observed_distance <- stats::setNames(
  mofa_factor_pair_summary$centroid_distance,
  mofa_factor_pair_summary$pair_id
)

pair_permuted_distance <- matrix(
  NA_real_,
  nrow = mofa_pair_permutations,
  ncol = nrow(mofa_factor_pair_summary),
  dimnames = list(
    NULL,
    mofa_factor_pair_summary$pair_id
  )
)

baseline_score_matrix <- as.matrix(
  baseline_factor_data[
    ,
    mofa_factor_names,
    drop = FALSE
  ]
)

baseline_score_matrix <- standardize_matrix_columns(
  baseline_score_matrix
)

for (permutation_index in seq_len(
  mofa_pair_permutations
)) {
  permuted_group <- sample(
    baseline_factor_data$TRG_plot,
    size = nrow(baseline_factor_data),
    replace = FALSE
  )

  permuted_difference <- colMeans(
    baseline_score_matrix[
      permuted_group == "pCR",
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  ) -
    colMeans(
      baseline_score_matrix[
        permuted_group == "non_pCR",
        ,
        drop = FALSE
      ],
      na.rm = TRUE
    )

  for (pair_index in seq_len(
    nrow(mofa_factor_pair_summary)
  )) {
    pair_factors <- c(
      mofa_factor_pair_summary$factor_x[pair_index],
      mofa_factor_pair_summary$factor_y[pair_index]
    )

    pair_permuted_distance[
      permutation_index,
      pair_index
    ] <- sqrt(
      sum(
        permuted_difference[
          pair_factors
        ] ^ 2
      )
    )
  }
}

max_permuted_distance <- apply(
  pair_permuted_distance,
  1,
  max,
  na.rm = TRUE
)

mofa_factor_pair_summary$pair_permutation_p <-
  vapply(
    seq_len(
      nrow(mofa_factor_pair_summary)
    ),
    function(pair_index) {
      permuted_values <- pair_permuted_distance[
        ,
        pair_index
      ]
      permuted_values <- permuted_values[
        is.finite(permuted_values)
      ]

      (
        1 +
          sum(
            permuted_values >=
              mofa_factor_pair_summary$centroid_distance[
                pair_index
              ]
          )
      ) /
      (
        1 +
          length(permuted_values)
      )
    },
    numeric(1)
  )

mofa_factor_pair_summary$max_statistic_p <-
  vapply(
    mofa_factor_pair_summary$centroid_distance,
    function(observed_distance) {
      (
        1 +
          sum(
            max_permuted_distance >=
              observed_distance
          )
      ) /
      (
        1 +
          length(max_permuted_distance)
      )
    },
    numeric(1)
  )

mofa_factor_pair_summary <-
  mofa_factor_pair_summary %>%
  dplyr::left_join(
    mofa_factor_subject_specificity %>%
      dplyr::select(
        factor,
        icc_x = icc,
        subject_dominated_x =
          subject_dominated
      ),
    by = c(
      "factor_x" = "factor"
    )
  ) %>%
  dplyr::left_join(
    mofa_factor_subject_specificity %>%
      dplyr::select(
        factor,
        icc_y = icc,
        subject_dominated_y =
          subject_dominated
      ),
    by = c(
      "factor_y" = "factor"
    )
  ) %>%
  dplyr::mutate(
    both_subject_dominated =
      dplyr::coalesce(
        subject_dominated_x,
        FALSE
      ) &
      dplyr::coalesce(
        subject_dominated_y,
        FALSE
      )
  ) %>%
  dplyr::mutate(
    response_rank = dplyr::min_rank(
      dplyr::desc(loocv_lda_auc)
    ),
    view_similarity_rank = dplyr::min_rank(
      dplyr::desc(view_fingerprint_cosine)
    ),
    view_complementarity_rank = dplyr::min_rank(
      view_fingerprint_cosine
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(loocv_lda_auc),
    max_statistic_p,
    dplyr::desc(centroid_distance)
  )

mofa_factor_pair_strategy <- dplyr::bind_rows(
  mofa_factor_pair_summary %>%
    dplyr::slice_min(
      response_rank,
      n = 5,
      with_ties = FALSE
    ) %>%
    dplyr::mutate(strategy = "Response-oriented"),
  mofa_factor_pair_summary %>%
    dplyr::filter(abs(score_correlation) < 0.80) %>%
    dplyr::slice_min(
      view_similarity_rank,
      n = 5,
      with_ties = FALSE
    ) %>%
    dplyr::mutate(strategy = "Similar omics fingerprint"),
  mofa_factor_pair_summary %>%
    dplyr::slice_min(
      view_complementarity_rank,
      n = 5,
      with_ties = FALSE
    ) %>%
    dplyr::mutate(strategy = "Complementary omics fingerprint")
) %>%
  dplyr::distinct(strategy, pair_id, .keep_all = TRUE)

make_factor_axis_label <- function(factor_name) {
  factor_row <- mofa_factor_summary %>%
    dplyr::filter(
      as.character(factor) == factor_name
    )

  if (nrow(factor_row) == 0) {
    return(factor_name)
  }

  paste0(
    factor_name,
    " (",
    unname(
      view_labels[
        factor_row$strongest_view[1]
      ]
    ),
    " ",
    scales::percent(
      factor_row$strongest_view_r2[1],
      accuracy = 0.1
    ),
    ")"
  )
}

if (
  !is.null(mofa_selected_pair_override) &&
  length(mofa_selected_pair_override) == 2 &&
  all(
    mofa_selected_pair_override %in%
      mofa_factor_names
  )
) {
  mofa_selected_factor_pair <-
    as.character(
      mofa_selected_pair_override
    )
} else {
  pair_candidates <- mofa_factor_pair_summary %>%
    dplyr::filter(
      !both_subject_dominated
    )

  if (nrow(pair_candidates) == 0) {
    pair_candidates <-
      mofa_factor_pair_summary
  }

  mofa_selected_factor_pair <- c(
    pair_candidates$factor_x[1],
    pair_candidates$factor_y[1]
  )
}

selected_pair_row <- mofa_factor_pair_summary %>%
  dplyr::filter(
    factor_x ==
      mofa_selected_factor_pair[1],
    factor_y ==
      mofa_selected_factor_pair[2]
  )

if (nrow(selected_pair_row) == 0) {
  selected_pair_row <- mofa_factor_pair_summary %>%
    dplyr::filter(
      factor_x ==
        mofa_selected_factor_pair[2],
      factor_y ==
        mofa_selected_factor_pair[1]
    )
}

selected_pair_data <- baseline_factor_data %>%
  dplyr::select(
    SampleID,
    SubjectID,
    TRG_plot,
    dplyr::all_of(
      mofa_selected_factor_pair
    )
  )

selected_pair_x <-
  mofa_selected_factor_pair[1]
selected_pair_y <-
  mofa_selected_factor_pair[2]

p_selected_pair_scatter <- ggplot2::ggplot(
  selected_pair_data,
  ggplot2::aes(
    x = .data[[selected_pair_x]],
    y = .data[[selected_pair_y]],
    fill = TRG_plot
  )
) +
  ggplot2::geom_point(
    shape = 21,
    size = 3.4,
    color = "grey20",
    stroke = 0.65,
    alpha = 0.88
  ) +
  ggplot2::stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 5.0,
    color = "black",
    stroke = 0.9
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      non_pCR = "grey75",
      pCR = "white"
    ),
    name = "Response"
  ) +
  ggplot2::labs(
    title = paste0(
      "Baseline response-oriented factor pair: ",
      selected_pair_x,
      " and ",
      selected_pair_y
    ),
    subtitle = paste0(
      "LOOCV LDA AUC = ",
      sprintf(
        "%.2f",
        selected_pair_row$loocv_lda_auc[1]
      ),
      "; pair permutation p = ",
      format.pval(
        selected_pair_row$pair_permutation_p[1],
        digits = 2,
        eps = 0.001
      ),
      "; search-adjusted max-statistic p = ",
      format.pval(
        selected_pair_row$max_statistic_p[1],
        digits = 2,
        eps = 0.001
      )
    ),
    x = make_factor_axis_label(
      selected_pair_x
    ),
    y = make_factor_axis_label(
      selected_pair_y
    )
  ) +
  ggplot2::theme_classic(
    base_size = 10.5
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 8.5,
      color = "grey35"
    ),
    legend.position = "bottom"
  )

p_selected_pair_x <- ggplot2::ggplot(
  selected_pair_data,
  ggplot2::aes(
    x = TRG_plot,
    y = .data[[selected_pair_x]],
    fill = TRG_plot
  )
) +
  ggplot2::geom_boxplot(
    width = 0.54,
    outlier.shape = NA,
    alpha = 0.55
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(
      width = 0.08,
      height = 0
    ),
    shape = 21,
    size = 2.0
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      non_pCR = "grey75",
      pCR = "white"
    ),
    guide = "none"
  ) +
  ggplot2::labs(
    x = NULL,
    y = selected_pair_x
  ) +
  ggplot2::theme_classic(
    base_size = 9.2
  )

p_selected_pair_y <- ggplot2::ggplot(
  selected_pair_data,
  ggplot2::aes(
    x = TRG_plot,
    y = .data[[selected_pair_y]],
    fill = TRG_plot
  )
) +
  ggplot2::geom_boxplot(
    width = 0.54,
    outlier.shape = NA,
    alpha = 0.55
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(
      width = 0.08,
      height = 0
    ),
    shape = 21,
    size = 2.0
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      non_pCR = "grey75",
      pCR = "white"
    ),
    guide = "none"
  ) +
  ggplot2::labs(
    x = NULL,
    y = selected_pair_y
  ) +
  ggplot2::theme_classic(
    base_size = 9.2
  )

p_mofa_selected_pair_map <-
  p_selected_pair_scatter |
  (
    p_selected_pair_x /
      p_selected_pair_y
  ) +
  patchwork::plot_layout(
    widths = c(
      2.2,
      1
    )
  )

ranked_pair_rows <- utils::head(
  mofa_factor_pair_summary,
  6
)

ranked_pair_plot_list <- vector(
  "list",
  nrow(ranked_pair_rows)
)

for (pair_index in seq_len(
  nrow(ranked_pair_rows)
)) {
  factor_x <- ranked_pair_rows$factor_x[
    pair_index
  ]
  factor_y <- ranked_pair_rows$factor_y[
    pair_index
  ]

  ranked_pair_plot_list[[pair_index]] <-
    ggplot2::ggplot(
      baseline_factor_data,
      ggplot2::aes(
        x = .data[[factor_x]],
        y = .data[[factor_y]],
        fill = TRG_plot
      )
    ) +
    ggplot2::geom_point(
      shape = 21,
      size = 2.7,
      color = "grey20",
      stroke = 0.55,
      alpha = 0.86
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        non_pCR = "grey75",
        pCR = "white"
      ),
      name = "Response"
    ) +
    ggplot2::labs(
      title = paste(
        factor_x,
        "vs",
        factor_y
      ),
      subtitle = paste0(
        "AUC ",
        sprintf(
          "%.2f",
          ranked_pair_rows$loocv_lda_auc[
            pair_index
          ]
        ),
        "; max-p ",
        format.pval(
          ranked_pair_rows$max_statistic_p[
            pair_index
          ],
          digits = 2,
          eps = 0.001
        )
      ),
      x = make_factor_axis_label(
        factor_x
      ),
      y = make_factor_axis_label(
        factor_y
      )
    ) +
    ggplot2::theme_classic(
      base_size = 9.1
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "plain",
        size = 9.5
      ),
      plot.subtitle = ggplot2::element_text(
        size = 7.7,
        color = "grey35"
      )
    )
}

p_mofa_ranked_pair_maps <- patchwork::wrap_plots(
  ranked_pair_plot_list,
  ncol = 3,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "Top baseline factor pairs for pCR/non-pCR separation",
    subtitle = paste0(
      "The display is exploratory; max-statistic permutation p-values account ",
      "for screening all factor pairs."
    )
  ) &
  ggplot2::theme(
    legend.position = "bottom"
  )


pair_landscape_label_data <-
  mofa_factor_pair_summary %>%
  dplyr::filter(
    response_rank <= 5 |
      view_similarity_rank <= 3 |
      view_complementarity_rank <= 3
  )

p_mofa_factor_pair_landscape <- ggplot2::ggplot(
  mofa_factor_pair_summary,
  ggplot2::aes(
    x = view_fingerprint_cosine,
    y = loocv_lda_auc,
    size = centroid_distance,
    fill = -log10(
      pmax(max_statistic_p, .Machine$double.xmin)
    )
  )
) +
  ggplot2::geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey60"
  ) +
  ggplot2::geom_point(
    shape = 21,
    color = "grey20",
    stroke = 0.45,
    alpha = 0.85
  ) +
  ggplot2::geom_text(
    data = pair_landscape_label_data,
    ggplot2::aes(
      label = paste(factor_x, factor_y, sep = "-")
    ),
    size = 2.65,
    vjust = -0.75,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  ggplot2::scale_fill_gradient(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "white",
    high = "grey20",
    name = "-log10 search-adjusted p"
  ) +
  ggplot2::scale_size_continuous(
    range = c(2.0, 5.5),
    name = "Baseline centroid distance"
  ) +
  ggplot2::labs(
    title = "Factor-pair landscape",
    subtitle = paste0(
      "Pairs can be prioritized for response separation, similar omics fingerprints, ",
      "or complementary omics fingerprints."
    ),
    x = "Cosine similarity of factor-by-view variance fingerprints",
    y = "Baseline pCR/non-pCR LOOCV LDA AUC"
  ) +
  ggplot2::theme_classic(base_size = 10.2) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 8.6, color = "grey35"),
    legend.position = "bottom"
  )


#=================================================================#
# 16. Paired feature changes and feature-centered latent axes
#=================================================================#

calculate_feature_paired_change <- function(
    view_name
) {
  view_matrix <- base::t(
    mofa_data[[view_name]]
  )

  measured_view_samples <- rownames(view_matrix)[
    rowSums(
      is.finite(view_matrix)
    ) > 0
  ]

  view_metadata <- mofa_sample_metadata %>%
    dplyr::filter(
      SampleID %in%
        measured_view_samples,
      Timepoint %in% c(
        "Before",
        "Ongoing"
      ),
      !is.na(TRG_plot)
    )

  before_meta <- view_metadata %>%
    dplyr::filter(
      Timepoint == "Before"
    ) %>%
    dplyr::select(
      SubjectID,
      TRG_plot,
      BeforeSampleID = SampleID
    ) %>%
    dplyr::distinct(
      SubjectID,
      .keep_all = TRUE
    )

  ongoing_meta <- view_metadata %>%
    dplyr::filter(
      Timepoint == "Ongoing"
    ) %>%
    dplyr::select(
      SubjectID,
      OngoingSampleID = SampleID
    ) %>%
    dplyr::distinct(
      SubjectID,
      .keep_all = TRUE
    )

  paired_meta <- dplyr::inner_join(
    before_meta,
    ongoing_meta,
    by = "SubjectID"
  )

  if (nrow(paired_meta) < 3) {
    return(
      data.frame()
    )
  }

  delta_matrix <-
    view_matrix[
      paired_meta$OngoingSampleID,
      ,
      drop = FALSE
    ] -
    view_matrix[
      paired_meta$BeforeSampleID,
      ,
      drop = FALSE
    ]

  rownames(delta_matrix) <-
    paired_meta$SubjectID

  feature_rows <- vector(
    "list",
    ncol(delta_matrix)
  )

  for (feature_index in seq_len(
    ncol(delta_matrix)
  )) {
    feature_name <- colnames(delta_matrix)[
      feature_index
    ]
    feature_delta <- delta_matrix[
      ,
      feature_index
    ]

    time_p <- suppressWarnings(
      tryCatch(
        stats::wilcox.test(
          feature_delta,
          mu = 0,
          exact = FALSE
        )$p.value,
        error = function(e) {
          NA_real_
        }
      )
    )

    response_p <- if (
      sum(
        paired_meta$TRG_plot == "pCR"
      ) >= 2 &&
      sum(
        paired_meta$TRG_plot == "non_pCR"
      ) >= 2
    ) {
      suppressWarnings(
        tryCatch(
          stats::wilcox.test(
            feature_delta[
              paired_meta$TRG_plot == "pCR"
            ],
            feature_delta[
              paired_meta$TRG_plot == "non_pCR"
            ],
            exact = FALSE
          )$p.value,
          error = function(e) {
            NA_real_
          }
        )
      )
    } else {
      NA_real_
    }

    feature_rows[[feature_index]] <-
      data.frame(
        view = view_name,
        feature = feature_name,
        n_paired = nrow(paired_meta),
        n_pCR_paired = sum(
          paired_meta$TRG_plot == "pCR"
        ),
        n_non_pCR_paired = sum(
          paired_meta$TRG_plot == "non_pCR"
        ),
        mean_delta = mean(
          feature_delta,
          na.rm = TRUE
        ),
        time_p = time_p,
        differential_delta =
          mean(
            feature_delta[
              paired_meta$TRG_plot == "pCR"
            ],
            na.rm = TRUE
          ) -
          mean(
            feature_delta[
              paired_meta$TRG_plot == "non_pCR"
            ],
            na.rm = TRUE
          ),
        differential_p = response_p,
        stringsAsFactors = FALSE
      )
  }

  dplyr::bind_rows(
    feature_rows
  ) %>%
    dplyr::mutate(
      time_fdr = stats::p.adjust(
        time_p,
        method = "BH"
      ),
      differential_fdr = stats::p.adjust(
        differential_p,
        method = "BH"
      )
    )
}

mofa_feature_paired_change_stats <-
  dplyr::bind_rows(
    lapply(
      required_views,
      calculate_feature_paired_change
    )
  )

baseline_factor_effect <-
  mofa_factor_hierarchical_models %>%
  dplyr::filter(
    model_type == "baseline_response"
  ) %>%
  dplyr::select(
    factor,
    baseline_factor_effect = estimate,
    baseline_factor_permutation_p =
      permutation_p
  )

paired_factor_effect <-
  mofa_factor_hierarchical_models %>%
  dplyr::filter(
    model_type ==
      "paired_differential_change"
  ) %>%
  dplyr::select(
    factor,
    paired_factor_effect = estimate,
    paired_factor_permutation_p =
      permutation_p
  )

baseline_factor_effect$baseline_effect_z <-
  as.numeric(
    scale(
      baseline_factor_effect$baseline_factor_effect
    )
  )

paired_factor_effect$paired_effect_z <-
  if (
    nrow(paired_factor_effect) >= 2
  ) {
    as.numeric(
      scale(
        paired_factor_effect$paired_factor_effect
      )
    )
  } else {
    NA_real_
  }

mofa_feature_factor_long <-
  mofa_feature_weights %>%
  dplyr::mutate(
    factor = as.character(factor)
  ) %>%
  dplyr::left_join(
    mofa_variance_explained %>%
      dplyr::transmute(
        view = as.character(view),
        factor = as.character(factor),
        factor_view_r2 = r2
      ),
    by = c(
      "view",
      "factor"
    )
  ) %>%
  dplyr::group_by(
    view,
    factor
  ) %>%
  dplyr::mutate(
    absolute_weight_percentile =
      dplyr::percent_rank(
        abs_weight
      ),
    recurrent_high_weight =
      absolute_weight_percentile >= 0.95 &
      factor_view_r2 >=
        mofa_active_view_r2
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(
    baseline_factor_effect,
    by = "factor"
  ) %>%
  dplyr::left_join(
    paired_factor_effect,
    by = "factor"
  ) %>%
  dplyr::mutate(
    baseline_axis_component =
      weight_within_view *
      baseline_effect_z,
    paired_axis_component =
      weight_within_view *
      paired_effect_z
  )

mofa_feature_centered_summary <-
  mofa_feature_factor_long %>%
  dplyr::group_by(
    view,
    view_label,
    feature,
    feature_label,
    display_eligible
  ) %>%
  dplyr::summarise(
    recurrent_factor_count = sum(
      recurrent_high_weight,
      na.rm = TRUE
    ),
    recurrent_strength = sum(
      abs(weight_within_view) *
        recurrent_high_weight,
      na.rm = TRUE
    ),
    recurrent_sign_consistency = if (
      sum(
        recurrent_high_weight,
        na.rm = TRUE
      ) > 0
    ) {
      abs(
        mean(
          sign(
            weight_within_view[
              recurrent_high_weight
            ]
          ),
          na.rm = TRUE
        )
      )
    } else {
      NA_real_
    },
    baseline_latent_response_axis = if (
      any(is.finite(baseline_axis_component))
    ) {
      sum(
        baseline_axis_component,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    paired_latent_response_axis = if (
      any(is.finite(paired_axis_component))
    ) {
      sum(
        paired_axis_component,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    maximum_absolute_weight = max(
      abs_weight,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    mofa_feature_response_stats %>%
      dplyr::select(
        view,
        feature,
        baseline_direct_effect =
          response_effect,
        baseline_direct_p = p_value,
        baseline_direct_fdr = fdr
      ),
    by = c(
      "view",
      "feature"
    )
  ) %>%
  dplyr::left_join(
    mofa_feature_paired_change_stats %>%
      dplyr::select(
        view,
        feature,
        mean_paired_change = mean_delta,
        paired_time_p = time_p,
        paired_differential_effect =
          differential_delta,
        paired_differential_p =
          differential_p
      ),
    by = c(
      "view",
      "feature"
    )
  ) %>%
  dplyr::group_by(
    view
  ) %>%
  dplyr::mutate(
    baseline_axis_rank = rank(
      -abs(
        baseline_latent_response_axis
      ),
      ties.method = "min"
    ),
    paired_axis_rank = rank(
      -abs(
        paired_latent_response_axis
      ),
      ties.method = "min"
    )
  ) %>%
  dplyr::ungroup()

recurrent_feature_selection <-
  mofa_feature_centered_summary %>%
  dplyr::filter(
    display_eligible
  ) %>%
  dplyr::group_by(
    view
  ) %>%
  dplyr::arrange(
    dplyr::desc(recurrent_factor_count),
    dplyr::desc(recurrent_strength),
    dplyr::desc(
      abs(
        baseline_latent_response_axis
      )
    ),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(
    n = 10
  ) %>%
  dplyr::ungroup()

recurrent_heatmap_data <-
  mofa_feature_factor_long %>%
  dplyr::semi_join(
    recurrent_feature_selection,
    by = c(
      "view",
      "feature"
    )
  ) %>%
  dplyr::mutate(
    factor = factor(
      factor,
      levels = as.character(
        mofa_factor_summary$factor
      )
    ),
    feature_plot = paste(
      view,
      feature_label,
      sep = "|||"
    )
  )

feature_plot_levels <- recurrent_feature_selection %>%
  dplyr::arrange(
    match(
      view,
      required_views
    ),
    recurrent_factor_count,
    recurrent_strength
  ) %>%
  dplyr::mutate(
    feature_plot = paste(
      view,
      feature_label,
      sep = "|||"
    )
  ) %>%
  dplyr::pull(feature_plot) %>%
  unique()

recurrent_heatmap_data$feature_plot <- factor(
  recurrent_heatmap_data$feature_plot,
  levels = feature_plot_levels
)

p_recurrent_weight_heatmap <- ggplot2::ggplot(
  recurrent_heatmap_data,
  ggplot2::aes(
    x = factor,
    y = feature_plot,
    fill = weight_within_view
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.25
  ) +
  ggplot2::facet_grid(
    view_label ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::scale_y_discrete(
    labels = function(x) {
      sub(
        "^[^|]+[|][|][|]",
        "",
        x
      )
    }
  ) +
  ggplot2::scale_fill_gradient2(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "grey20",
    mid = "white",
    high = "grey55",
    midpoint = 0,
    limits = c(
      -1,
      1
    ),
    name = "Normalized
weight"
  ) +
  ggplot2::labs(
    title = "Recurrent high-impact features across MOFA factors",
    subtitle = paste0(
      "Features are prioritized by recurrence among the top 5% of weights in ",
      "active factor-view combinations."
    ),
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 9.6
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 8.4,
      color = "grey35"
    ),
    strip.background = ggplot2::element_blank(),
    strip.text.y = ggplot2::element_text(
      angle = 0,
      hjust = 0,
      face = "plain"
    ),
    axis.text.x = ggplot2::element_text(
      face = "plain"
    ),
    axis.text.y = ggplot2::element_text(
      size = 7.2
    )
  )

recurrent_axis_tile_data <- recurrent_feature_selection %>%
  dplyr::mutate(
    feature_plot = factor(
      paste(
        view,
        feature_label,
        sep = "|||"
      ),
      levels = feature_plot_levels
    )
  ) %>%
  dplyr::select(
    view_label,
    feature_plot,
    baseline_latent_response_axis,
    paired_latent_response_axis
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      baseline_latent_response_axis,
      paired_latent_response_axis
    ),
    names_to = "axis_type",
    values_to = "axis_score"
  ) %>%
  dplyr::mutate(
    axis_type = dplyr::recode(
      axis_type,
      baseline_latent_response_axis =
        "Baseline response axis",
      paired_latent_response_axis =
        "Differential-change axis"
    )
  )

axis_limit <- stats::quantile(
  abs(
    recurrent_axis_tile_data$axis_score
  ),
  probs = 0.95,
  na.rm = TRUE,
  names = FALSE
)

if (!is.finite(axis_limit) || axis_limit <= 0) {
  axis_limit <- 1
}

p_recurrent_axis_tiles <- ggplot2::ggplot(
  recurrent_axis_tile_data,
  ggplot2::aes(
    x = axis_type,
    y = feature_plot,
    fill = axis_score
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.25
  ) +
  ggplot2::facet_grid(
    view_label ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::scale_fill_gradient2(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "#3B75AF",
    mid = "white",
    high = "#C44E52",
    midpoint = 0,
    limits = c(
      -axis_limit,
      axis_limit
    ),
    oob = scales::squish,
    name = "Latent-axis
score"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_void(
    base_size = 8.7
  ) +
  ggplot2::theme(
    strip.text = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      angle = 25,
      hjust = 1
    )
  )

p_mofa_feature_recurrence_heatmap <-
  p_recurrent_weight_heatmap |
  p_recurrent_axis_tiles +
  patchwork::plot_layout(
    widths = c(
      5.0,
      1.5
    ),
    guides = "collect"
  ) &
  ggplot2::theme(
    legend.position = "bottom"
  )

selected_pair_response_vector <- c(
  selected_pair_row$response_vector_x[1],
  selected_pair_row$response_vector_y[1]
)
names(selected_pair_response_vector) <- c(
  selected_pair_row$factor_x[1],
  selected_pair_row$factor_y[1]
)
selected_pair_response_vector <-
  selected_pair_response_vector[
    mofa_selected_factor_pair
  ]

selected_pair_feature_metadata <-
  mofa_feature_weights %>%
  dplyr::filter(
    factor %in%
      mofa_selected_factor_pair,
    display_eligible
  ) %>%
  dplyr::select(
    view,
    view_label,
    feature,
    feature_label,
    response_effect,
    p_value,
    neg_log10_p
  ) %>%
  dplyr::distinct(
    view,
    feature,
    .keep_all = TRUE
  )

selected_pair_weight_wide <-
  mofa_feature_weights %>%
  dplyr::filter(
    factor %in%
      mofa_selected_factor_pair,
    display_eligible
  ) %>%
  dplyr::select(
    view,
    feature,
    factor,
    weight_within_view
  ) %>%
  tidyr::pivot_wider(
    names_from = factor,
    values_from = weight_within_view,
    values_fill = 0,
    values_fn = dplyr::first
  ) %>%
  dplyr::left_join(
    selected_pair_feature_metadata,
    by = c(
      "view",
      "feature"
    )
  )

mofa_selected_pair_feature_axis <-
  selected_pair_weight_wide %>%
  dplyr::mutate(
    pair_axis_score =
      .data[[mofa_selected_factor_pair[1]]] *
      selected_pair_response_vector[
        mofa_selected_factor_pair[1]
      ] +
      .data[[mofa_selected_factor_pair[2]]] *
      selected_pair_response_vector[
        mofa_selected_factor_pair[2]
      ]
  ) %>%
  dplyr::group_by(
    view
  ) %>%
  dplyr::mutate(
    pair_axis_denominator = max(
      abs(pair_axis_score),
      na.rm = TRUE
    ),
    pair_axis_score_normalized = ifelse(
      is.finite(pair_axis_denominator) &
        pair_axis_denominator > 0,
      pair_axis_score /
        pair_axis_denominator,
      0
    )
  ) %>%
  dplyr::ungroup()

selected_pair_axis_plot_data <-
  mofa_selected_pair_feature_axis %>%
  dplyr::group_by(
    view,
    view_label
  ) %>%
  dplyr::slice_max(
    order_by = abs(
      pair_axis_score_normalized
    ),
    n = 10,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    feature_plot = paste(
      view,
      feature_label,
      sep = "|||"
    )
  ) %>%
  dplyr::arrange(
    match(
      view,
      required_views
    ),
    pair_axis_score_normalized
  )

selected_pair_axis_plot_data$feature_plot <-
  factor(
    selected_pair_axis_plot_data$feature_plot,
    levels = unique(
      selected_pair_axis_plot_data$feature_plot
    )
  )

p_pair_axis_weight <- ggplot2::ggplot(
  selected_pair_axis_plot_data,
  ggplot2::aes(
    x = pair_axis_score_normalized,
    y = feature_plot,
    fill = view_label
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    color = "grey55",
    linewidth = 0.45
  ) +
  ggplot2::geom_col(
    width = 0.70,
    color = "grey25",
    linewidth = 0.25
  ) +
  ggplot2::facet_grid(
    view_label ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::scale_y_discrete(
    labels = function(x) {
      sub(
        "^[^|]+[|][|][|]",
        "",
        x
      )
    }
  ) +
  ggplot2::scale_fill_manual(
    values = view_colors,
    guide = "none"
  ) +
  ggplot2::labs(
    title = paste0(
      "Features aligned with the ",
      selected_pair_x,
      "-",
      selected_pair_y,
      " baseline response axis"
    ),
    subtitle = paste0(
      "The axis projects the pCR/non-pCR centroid difference through MOFA ",
      "weights; it is exploratory rather than an independent feature test."
    ),
    x = "Within-view normalized pair-axis contribution",
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 9.5
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 8.2,
      color = "grey35"
    ),
    strip.background = ggplot2::element_blank(),
    strip.text.y = ggplot2::element_text(
      angle = 0,
      hjust = 0,
      face = "plain"
    ),
    axis.text.y = ggplot2::element_text(
      size = 7.2
    )
  )

p_pair_axis_direct_effect <- ggplot2::ggplot(
  selected_pair_axis_plot_data,
  ggplot2::aes(
    x = "Direct",
    y = feature_plot,
    fill = response_effect
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.25
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = ifelse(
        is.finite(response_effect),
        sprintf(
          "%.2f",
          response_effect
        ),
        ""
      )
    ),
    size = 2.25
  ) +
  ggplot2::facet_grid(
    view_label ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::scale_fill_gradient2(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "#3B75AF",
    mid = "white",
    high = "#C44E52",
    midpoint = 0,
    name = "Direct
baseline effect"
  ) +
  ggplot2::labs(
    x = "pCR vs non-pCR",
    y = NULL
  ) +
  ggplot2::theme_void(
    base_size = 8.5
  ) +
  ggplot2::theme(
    strip.text = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_text(
      size = 7.8
    )
  )

p_pair_axis_p <- ggplot2::ggplot(
  selected_pair_axis_plot_data,
  ggplot2::aes(
    x = "P",
    y = feature_plot,
    fill = neg_log10_p
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.25
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = ifelse(
        is.finite(neg_log10_p),
        sprintf(
          "%.1f",
          neg_log10_p
        ),
        ""
      )
    ),
    size = 2.25
  ) +
  ggplot2::facet_grid(
    view_label ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::scale_fill_gradient(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "white",
    high = "grey20",
    name = "-log10(p)"
  ) +
  ggplot2::labs(
    x = "-log10(p)",
    y = NULL
  ) +
  ggplot2::theme_void(
    base_size = 8.5
  ) +
  ggplot2::theme(
    strip.text = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_text(
      size = 7.8
    )
  )

p_mofa_selected_pair_feature_axis <-
  p_pair_axis_weight |
  p_pair_axis_direct_effect |
  p_pair_axis_p +
  patchwork::plot_layout(
    widths = c(
      6.0,
      1.2,
      1.0
    ),
    guides = "collect"
  ) &
  ggplot2::theme(
    legend.position = "bottom"
  )


#=================================================================#
# 17. Factor-feature constellation: a compact multi-factor view
#=================================================================#

constellation_edges <-
  mofa_feature_factor_long %>%
  dplyr::filter(
    factor_view_r2 >=
      mofa_active_view_r2,
    display_eligible
  ) %>%
  dplyr::group_by(
    view,
    factor
  ) %>%
  dplyr::slice_max(
    order_by = abs_weight,
    n = 2,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    factor_y = match(
      factor,
      rev(
        as.character(
          mofa_factor_summary$factor
        )
      )
    ),
    direction = factor(
      direction,
      levels = c(
        "Negative",
        "Positive"
      )
    )
  ) %>%
  dplyr::group_by(
    view
  ) %>%
  dplyr::arrange(
    factor_y,
    dplyr::desc(abs_weight),
    .by_group = TRUE
  ) %>%
  dplyr::mutate(
    feature_y = seq(
      min(factor_y, na.rm = TRUE),
      max(factor_y, na.rm = TRUE),
      length.out = dplyr::n()
    )
  ) %>%
  dplyr::ungroup()

constellation_factor_nodes <-
  constellation_edges %>%
  dplyr::distinct(
    view,
    view_label,
    factor,
    factor_y
  )

constellation_feature_nodes <-
  constellation_edges %>%
  dplyr::distinct(
    view,
    view_label,
    feature_label,
    feature_y,
    direction
  )

p_mofa_factor_feature_constellation <-
  ggplot2::ggplot() +
  ggplot2::geom_curve(
    data = constellation_edges,
    ggplot2::aes(
      x = 0,
      y = factor_y,
      xend = 1,
      yend = feature_y,
      color = direction,
      linewidth = abs(
        weight_within_view
      )
    ),
    curvature = 0.12,
    alpha = 0.58
  ) +
  ggplot2::geom_point(
    data = constellation_factor_nodes,
    ggplot2::aes(
      x = 0,
      y = factor_y
    ),
    shape = 21,
    size = 2.8,
    fill = "white",
    color = "grey20"
  ) +
  ggplot2::geom_text(
    data = constellation_factor_nodes,
    ggplot2::aes(
      x = -0.04,
      y = factor_y,
      label = factor
    ),
    hjust = 1,
    size = 2.8
  ) +
  ggplot2::geom_point(
    data = constellation_feature_nodes,
    ggplot2::aes(
      x = 1,
      y = feature_y,
      color = direction
    ),
    size = 1.8
  ) +
  ggplot2::geom_text(
    data = constellation_feature_nodes,
    ggplot2::aes(
      x = 1.04,
      y = feature_y,
      label = feature_label
    ),
    hjust = 0,
    size = 2.35
  ) +
  ggplot2::facet_wrap(
    ~ view_label,
    nrow = 1,
    scales = "free_y"
  ) +
  ggplot2::scale_color_manual(
    values = c(
      Negative = "grey25",
      Positive = "grey60"
    ),
    name = "Weight sign"
  ) +
  ggplot2::scale_linewidth_continuous(
    range = c(
      0.25,
      1.35
    ),
    guide = "none"
  ) +
  ggplot2::coord_cartesian(
    xlim = c(
      -0.32,
      1.75
    ),
    clip = "off"
  ) +
  ggplot2::labs(
    title = "Factor-feature constellation",
    subtitle = paste0(
      "Each active factor-view combination contributes its two highest-weight ",
      "publication-eligible features."
    ),
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_void(
    base_size = 9.3
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 8.3,
      color = "grey35"
    ),
    strip.text = ggplot2::element_text(
      face = "plain"
    ),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(
      5,
      105,
      5,
      38
    )
  )


#=================================================================#
# 18. Separate-scope MOFA sensitivity models
#=================================================================#
#
# These models answer whether the global all-time factors are reproducible when
# Before and Ongoing samples are separated, repeated measurements are averaged
# to one stable subject profile, or paired differences are analyzed directly.
# Paired-delta views with insufficient paired subjects, commonly Host RNA-seq in
# this cohort, are automatically excluded rather than imputed.
#-----------------------------------------------------------------#

make_scope_matrix_list <- function(
    scope_name
) {
  if (scope_name %in% c(
    "baseline",
    "ongoing"
  )) {
    target_timepoint <- if (
      scope_name == "baseline"
    ) {
      "Before"
    } else {
      "Ongoing"
    }

    scope_ids <- mofa_sample_metadata$SampleID[
      as.character(
        mofa_sample_metadata$Timepoint
      ) == target_timepoint
    ]

    matrix_list <- lapply(
      mofa_data,
      function(view_matrix) {
        view_matrix[
          ,
          intersect(
            scope_ids,
            colnames(view_matrix)
          ),
          drop = FALSE
        ]
      }
    )
  } else if (scope_name == "four_view_complete") {
    four_view_ids <- mofa_availability_wide %>%
      dplyr::filter(
        dplyr::if_all(
          dplyr::all_of(required_views),
          ~ .x == 1L
        )
      ) %>%
      dplyr::pull(SampleID)

    matrix_list <- lapply(
      mofa_data,
      function(view_matrix) {
        view_matrix[
          ,
          intersect(
            four_view_ids,
            colnames(view_matrix)
          ),
          drop = FALSE
        ]
      }
    )
  } else if (scope_name == "subject_mean") {
    scope_subjects <- unique(
      mofa_sample_metadata$SubjectID
    )

    matrix_list <- lapply(
      mofa_data,
      function(view_matrix) {
        subject_matrix <- matrix(
          NA_real_,
          nrow = nrow(view_matrix),
          ncol = length(scope_subjects),
          dimnames = list(
            rownames(view_matrix),
            scope_subjects
          )
        )

        for (subject_id in scope_subjects) {
          subject_sample_ids <-
            mofa_sample_metadata$SampleID[
              mofa_sample_metadata$SubjectID ==
                subject_id
            ]

          subject_sample_ids <- intersect(
            subject_sample_ids,
            colnames(view_matrix)
          )

          if (length(subject_sample_ids) == 0) {
            next
          }

          subject_values <- view_matrix[
            ,
            subject_sample_ids,
            drop = FALSE
          ]

          subject_mean <- rowMeans(
            subject_values,
            na.rm = TRUE
          )
          subject_mean[
            rowSums(
              is.finite(subject_values)
            ) == 0
          ] <- NA_real_

          subject_matrix[, subject_id] <-
            subject_mean
        }

        subject_matrix
      }
    )
  } else if (scope_name == "within_subject_centered") {
    matrix_list <- lapply(
      mofa_data,
      function(view_matrix) {
        centered_matrix <- matrix(
          NA_real_,
          nrow = nrow(view_matrix),
          ncol = ncol(view_matrix),
          dimnames = dimnames(view_matrix)
        )

        for (subject_id in unique(
          mofa_sample_metadata$SubjectID
        )) {
          subject_sample_ids <-
            mofa_sample_metadata$SampleID[
              mofa_sample_metadata$SubjectID ==
                subject_id
            ]
          subject_sample_ids <- intersect(
            subject_sample_ids,
            colnames(view_matrix)
          )

          measured_subject_samples <-
            subject_sample_ids[
              colSums(
                is.finite(
                  view_matrix[
                    ,
                    subject_sample_ids,
                    drop = FALSE
                  ]
                )
              ) > 0
            ]

          if (length(measured_subject_samples) < 2) {
            next
          }

          subject_values <- view_matrix[
            ,
            measured_subject_samples,
            drop = FALSE
          ]
          subject_mean <- rowMeans(
            subject_values,
            na.rm = TRUE
          )
          subject_mean[
            rowSums(
              is.finite(subject_values)
            ) == 0
          ] <- NA_real_

          centered_matrix[
            ,
            measured_subject_samples
          ] <- sweep(
            subject_values,
            1,
            subject_mean,
            FUN = "-"
          )
        }

        centered_matrix
      }
    )
  } else if (scope_name == "paired_delta") {
    scope_subjects <- unique(
      mofa_sample_metadata$SubjectID
    )

    matrix_list <- lapply(
      mofa_data,
      function(view_matrix) {
        delta_matrix <- matrix(
          NA_real_,
          nrow = nrow(view_matrix),
          ncol = length(scope_subjects),
          dimnames = list(
            rownames(view_matrix),
            scope_subjects
          )
        )

        for (subject_id in scope_subjects) {
          subject_meta <- mofa_sample_metadata %>%
            dplyr::filter(
              SubjectID == subject_id,
              SampleID %in%
                colnames(view_matrix)
            )

          before_id <- subject_meta$SampleID[
            as.character(
              subject_meta$Timepoint
            ) == "Before"
          ]
          ongoing_id <- subject_meta$SampleID[
            as.character(
              subject_meta$Timepoint
            ) == "Ongoing"
          ]

          if (
            length(before_id) != 1 ||
            length(ongoing_id) != 1
          ) {
            next
          }

          before_values <- view_matrix[
            ,
            before_id,
            drop = TRUE
          ]
          ongoing_values <- view_matrix[
            ,
            ongoing_id,
            drop = TRUE
          ]

          valid_index <-
            is.finite(before_values) &
            is.finite(ongoing_values)

          delta_matrix[
            valid_index,
            subject_id
          ] <-
            ongoing_values[valid_index] -
            before_values[valid_index]
        }

        delta_matrix
      }
    )
  } else {
    stop(
      "Unknown scope: ",
      scope_name,
      call. = FALSE
    )
  }

  names(matrix_list) <- names(mofa_data)

  view_measurement_counts <- vapply(
    matrix_list,
    function(view_matrix) {
      sum(
        colSums(
          is.finite(view_matrix)
        ) > 0
      )
    },
    integer(1)
  )

  matrix_list <- matrix_list[
    view_measurement_counts >=
      mofa_scope_min_samples
  ]

  if (length(matrix_list) < 2) {
    return(NULL)
  }

  common_column_order <- unique(
    unlist(
      lapply(
        matrix_list,
        colnames
      ),
      use.names = FALSE
    )
  )

  for (view_name in names(matrix_list)) {
    missing_columns <- setdiff(
      common_column_order,
      colnames(
        matrix_list[[view_name]]
      )
    )

    if (length(missing_columns) > 0) {
      missing_matrix <- matrix(
        NA_real_,
        nrow = nrow(
          matrix_list[[view_name]]
        ),
        ncol = length(missing_columns),
        dimnames = list(
          rownames(
            matrix_list[[view_name]]
          ),
          missing_columns
        )
      )

      matrix_list[[view_name]] <- cbind(
        matrix_list[[view_name]],
        missing_matrix
      )
    }

    matrix_list[[view_name]] <-
      matrix_list[[view_name]][
        ,
        common_column_order,
        drop = FALSE
      ]
  }

  sample_view_count <- Reduce(
    "+",
    lapply(
      matrix_list,
      function(view_matrix) {
        as.integer(
          colSums(
            is.finite(view_matrix)
          ) > 0
        )
      }
    )
  )

  keep_columns <- common_column_order[
    sample_view_count >= 2
  ]

  if (length(keep_columns) < mofa_scope_min_samples) {
    return(NULL)
  }

  lapply(
    matrix_list,
    function(view_matrix) {
      view_matrix[
        ,
        keep_columns,
        drop = FALSE
      ]
    }
  )
}

fit_scope_mofa_model <- function(
    scope_name,
    scope_data
) {
  scope_model_path <- paste0(
    "results/mofa/MOFA_v16_scope_",
    scope_name,
    ".hdf5"
  )

  if (file.exists(scope_model_path)) {
    return(
      MOFA2::load_model(
        scope_model_path
      )
    )
  }

  scope_object <- MOFA2::create_mofa(
    scope_data
  )

  scope_data_options <-
    MOFA2::get_default_data_options(
      scope_object
    )
  scope_data_options$scale_views <- TRUE
  scope_data_options$scale_groups <- FALSE
  scope_data_options$center_groups <- TRUE
  scope_data_options$use_float32 <- TRUE

  scope_model_options <-
    MOFA2::get_default_model_options(
      scope_object
    )
  scope_model_options$num_factors <- min(
    mofa_scope_initial_factors,
    ncol(scope_data[[1]]) - 1
  )
  scope_model_options$likelihoods[] <-
    "gaussian"
  scope_model_options$ard_weights <- TRUE
  scope_model_options$ard_factors <- FALSE
  scope_model_options$spikeslab_weights <-
    FALSE
  scope_model_options$spikeslab_factors <-
    FALSE

  scope_training_options <-
    MOFA2::get_default_training_options(
      scope_object
    )
  scope_training_options$maxiter <- 2500
  scope_training_options$convergence_mode <-
    "medium"
  scope_training_options$drop_factor_threshold <-
    mofa_drop_factor_threshold
  scope_training_options$startELBO <- 1
  scope_training_options$freqELBO <- 5
  scope_training_options$verbose <- TRUE
  scope_training_options$seed <-
    mofa_scope_sensitivity_seed
  scope_training_options$stochastic <- FALSE
  scope_training_options$gpu_mode <- FALSE

  scope_object <- MOFA2::prepare_mofa(
    object = scope_object,
    data_options = scope_data_options,
    model_options = scope_model_options,
    training_options = scope_training_options
  )

  run_mofa_compat(
    object = scope_object,
    outfile = scope_model_path,
    save_data = FALSE,
    use_basilisk = mofa_use_basilisk
  )
}

extract_scope_weights <- function(
    scope_model,
    scope_name
) {
  scope_weight_list <- MOFA2::get_weights(
    scope_model,
    views = "all",
    factors = "all",
    abs = FALSE,
    scale = FALSE,
    as.data.frame = FALSE
  )

  dplyr::bind_rows(
    lapply(
      names(scope_weight_list),
      function(view_name) {
        weight_matrix <- as.matrix(
          scope_weight_list[[view_name]]
        )

        expected_scope_features <-
          MOFA2::features_names(
            scope_model
          )[[view_name]]

        if (
          !all(
            expected_scope_features %in%
              rownames(weight_matrix)
          ) &&
          all(
            expected_scope_features %in%
              colnames(weight_matrix)
          )
        ) {
          weight_matrix <- base::t(
            weight_matrix
          )
        }

        weight_matrix <- weight_matrix[
          intersect(
            expected_scope_features,
            rownames(weight_matrix)
          ),
          ,
          drop = FALSE
        ]

        data.frame(
          scope = scope_name,
          view = view_name,
          feature = rep(
            rownames(weight_matrix),
            times = ncol(weight_matrix)
          ),
          scope_factor = rep(
            colnames(weight_matrix),
            each = nrow(weight_matrix)
          ),
          scope_weight = as.numeric(
            weight_matrix
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )
}

extract_scope_variance <- function(
    scope_model,
    scope_name
) {
  variance_raw <- MOFA2::get_variance_explained(
    scope_model,
    groups = "all",
    views = "all",
    factors = "all",
    as.data.frame = TRUE
  )

  variance_data <- if (
    is.list(variance_raw) &&
    "r2_per_factor" %in%
      names(variance_raw)
  ) {
    variance_raw$r2_per_factor
  } else {
    variance_raw
  }

  r2_column <- intersect(
    c(
      "r2",
      "value"
    ),
    colnames(variance_data)
  )[1]

  variance_data %>%
    dplyr::transmute(
      scope = scope_name,
      view = as.character(view),
      scope_factor = as.character(factor),
      r2 = ifelse(
        .data[[r2_column]] > 1,
        .data[[r2_column]] / 100,
        .data[[r2_column]]
      )
    )
}

mofa_scope_model_summary <- data.frame()
mofa_scope_factor_similarity <- data.frame()
mofa_scope_models <- list()

mofa_scope_names_to_run <- c(
  if (isTRUE(mofa_run_four_view_complete_model)) {
    "four_view_complete"
  },
  if (isTRUE(mofa_run_scope_sensitivity_models)) {
    c(
      "baseline",
      "ongoing",
      "subject_mean",
      "within_subject_centered",
      "paired_delta"
    )
  }
)

if (length(mofa_scope_names_to_run) > 0) {
  for (scope_name in mofa_scope_names_to_run) {
    message(
      "Preparing scope sensitivity model: ",
      scope_name
    )

    scope_data <- make_scope_matrix_list(
      scope_name
    )

    if (is.null(scope_data)) {
      mofa_scope_model_summary <-
        dplyr::bind_rows(
          mofa_scope_model_summary,
          data.frame(
            scope = scope_name,
            status = "skipped: insufficient samples/views",
            n_samples = NA_integer_,
            n_views = NA_integer_,
            views = NA_character_,
            stringsAsFactors = FALSE
          )
        )
      next
    }

    scope_model <- tryCatch(
      fit_scope_mofa_model(
        scope_name,
        scope_data
      ),
      error = function(e) {
        message(
          "Scope model failed for ",
          scope_name,
          ": ",
          conditionMessage(e)
        )
        NULL
      }
    )

    if (is.null(scope_model)) {
      mofa_scope_model_summary <-
        dplyr::bind_rows(
          mofa_scope_model_summary,
          data.frame(
            scope = scope_name,
            status = "failed",
            n_samples = ncol(
              scope_data[[1]]
            ),
            n_views = length(scope_data),
            views = paste(
              names(scope_data),
              collapse = ";"
            ),
            stringsAsFactors = FALSE
          )
        )
      next
    }

    mofa_scope_models[[scope_name]] <-
      scope_model

    scope_weights <- extract_scope_weights(
      scope_model,
      scope_name
    )
    scope_variance <- extract_scope_variance(
      scope_model,
      scope_name
    )

    mofa_scope_model_summary <-
      dplyr::bind_rows(
        mofa_scope_model_summary,
        data.frame(
          scope = scope_name,
          status = "completed",
          n_samples = ncol(
            scope_data[[1]]
          ),
          n_views = length(scope_data),
          views = paste(
            names(scope_data),
            collapse = ";"
          ),
          n_factors = dplyr::n_distinct(
            scope_variance$scope_factor
          ),
          stringsAsFactors = FALSE
        )
      )

    for (primary_factor in mofa_factor_names) {
      for (scope_factor_name in unique(
        scope_weights$scope_factor
      )) {
        similarity_by_view <-
          mofa_feature_weights %>%
          dplyr::filter(
            factor == primary_factor,
            view %in%
              names(scope_data)
          ) %>%
          dplyr::select(
            view,
            feature,
            primary_weight = weight
          ) %>%
          dplyr::inner_join(
            scope_weights %>%
              dplyr::filter(
                .data$scope_factor ==
                  scope_factor_name
              ) %>%
              dplyr::select(
                view,
                feature,
                scope_weight
              ),
            by = c(
              "view",
              "feature"
            )
          ) %>%
          dplyr::group_by(
            view
          ) %>%
          dplyr::summarise(
            n_common_features = dplyr::n(),
            cosine = if (
              dplyr::n() >= 10 &&
              sqrt(
                sum(
                  primary_weight ^ 2
                )
              ) > 0 &&
              sqrt(
                sum(
                  scope_weight ^ 2
                )
              ) > 0
            ) {
              sum(
                primary_weight *
                  scope_weight
              ) /
              (
                sqrt(
                  sum(
                    primary_weight ^ 2
                  )
                ) *
                sqrt(
                  sum(
                    scope_weight ^ 2
                  )
                )
              )
            } else {
              NA_real_
            },
            .groups = "drop"
          ) %>%
          dplyr::filter(
            is.finite(cosine)
          )

        mofa_scope_factor_similarity <-
          dplyr::bind_rows(
            mofa_scope_factor_similarity,
            data.frame(
              scope = scope_name,
              primary_factor =
                primary_factor,
              scope_factor = scope_factor_name,
              mean_absolute_cosine = if (
                nrow(similarity_by_view) > 0
              ) {
                mean(
                  abs(
                    similarity_by_view$cosine
                  ),
                  na.rm = TRUE
                )
              } else {
                NA_real_
              },
              n_compared_views =
                nrow(similarity_by_view),
              stringsAsFactors = FALSE
            )
          )
      }
    }
  }
}

if (nrow(mofa_scope_factor_similarity) > 0) {
  mofa_scope_factor_similarity <-
    mofa_scope_factor_similarity %>%
    dplyr::filter(
      is.finite(mean_absolute_cosine)
    ) %>%
    dplyr::group_by(
      scope,
      primary_factor
    ) %>%
    dplyr::mutate(
      best_match =
        mean_absolute_cosine ==
        max(
          mean_absolute_cosine,
          na.rm = TRUE
        )
    ) %>%
    dplyr::ungroup()
} else {
  # Preserve a stable schema so downstream filter(best_match) calls are valid
  # even when all optional scope models are skipped or fail.
  mofa_scope_factor_similarity <- data.frame(
    scope = character(0),
    primary_factor = character(0),
    scope_factor = character(0),
    mean_absolute_cosine = numeric(0),
    n_compared_views = integer(0),
    best_match = logical(0),
    stringsAsFactors = FALSE
  )
}

scope_stability_plot_data <-
  mofa_scope_factor_similarity %>%
  dplyr::filter(
    best_match
  ) %>%
  dplyr::mutate(
    primary_factor = factor(
      primary_factor,
      levels = rev(
        as.character(
          mofa_factor_summary$factor
        )
      )
    ),
    scope = factor(
      scope,
      levels = c(
        "four_view_complete",
        "baseline",
        "ongoing",
        "subject_mean",
        "within_subject_centered",
        "paired_delta"
      )
    )
  )

if (nrow(scope_stability_plot_data) > 0) {
  p_mofa_scope_factor_stability <-
    ggplot2::ggplot(
      scope_stability_plot_data,
      ggplot2::aes(
        x = scope,
        y = primary_factor,
        fill = mean_absolute_cosine
      )
    ) +
    ggplot2::geom_tile(
      color = "white",
      linewidth = 0.35
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0(
          scope_factor,
          "\n",
          sprintf(
            "%.2f",
            mean_absolute_cosine
          )
        )
      ),
      size = 2.55
    ) +
    ggplot2::scale_fill_gradient(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
      low = "white",
      high = "grey25",
      limits = c(
        0,
        1
      ),
      name = "Absolute
loading cosine"
    ) +
    ggplot2::labs(
      title = "Recovery of all-time factors in separate-scope MOFA models",
      subtitle = paste0(
        "Each cell shows the best-matching scope factor and the mean absolute ",
        "loading cosine across shared views."
      ),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_classic(
      base_size = 10.0
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(
        size = 8.5,
        color = "grey35"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 25,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(
        face = "plain"
      ),
      axis.ticks = ggplot2::element_blank()
    )
} else {
  p_mofa_scope_factor_stability <-
    ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      label = "Scope sensitivity models were not available."
    ) +
    ggplot2::theme_void()
}


#=================================================================#
# 19. Export longitudinal and feature-atlas figures
#=================================================================#

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename =
    "figures/mofa/internal_v16/MOFA_v16_factor_subject_specificity.svg",
  plot = p_mofa_subject_specificity,
  width = 7.2,
  height = 4.2,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename =
    "figures/mofa/internal_v16/MOFA_v16_factor_hierarchical_forest.svg",
  plot = p_mofa_hierarchical_forest,
  width = 12.5,
  height = 5.8,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename =
    "figures/mofa/internal_v16/MOFA_v16_selected_factor_pair_map.svg",
  plot = p_mofa_selected_pair_map,
  width = 9.2,
  height = 5.2,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename =
    "figures/mofa/internal_v16/MOFA_v16_ranked_factor_pair_maps.svg",
  plot = p_mofa_ranked_pair_maps,
  width = 12.0,
  height = 7.8,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}


if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename =
    "figures/mofa/internal_v16/MOFA_v16_factor_pair_landscape.svg",
  plot = p_mofa_factor_pair_landscape,
  width = 7.8,
  height = 5.4,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename =
    "figures/mofa/internal_v16/MOFA_v16_feature_recurrence_heatmap.svg",
  plot = p_mofa_feature_recurrence_heatmap,
  width = 13.5,
  height = 10.5,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename =
    "figures/mofa/internal_v16/MOFA_v16_selected_pair_feature_axis.svg",
  plot = p_mofa_selected_pair_feature_axis,
  width = 14.0,
  height = 11.5,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename =
    "figures/mofa/internal_v16/MOFA_v16_factor_feature_constellation.svg",
  plot = p_mofa_factor_feature_constellation,
  width = 17.0,
  height = 8.5,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename =
    "figures/mofa/internal_v16/MOFA_v16_scope_factor_stability.svg",
  plot = p_mofa_scope_factor_stability,
  width = 7.8,
  height = 5.4,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}



#=================================================================#
# 20. Publication response plots, pairwise PERMANOVA, and all-factor summaries
#=================================================================#

# Factor scores are centered latent coordinates and can be negative. A literal
# log2(mean group 1 / mean group 0) is therefore not invariant and may be
# undefined. Publication figures use Hedges g for baseline response and paired
# dz for Ongoing-Before change. Shifted log2 ratios are exported only as a
# descriptive sensitivity column in MOFA_v16_factor_display_effects.csv.

mofa_factor_response_baseline_long <- mofa_factor_scores %>%
  dplyr::filter(
    as.character(Timepoint) == "Before",
    !is.na(TRG_plot)
  ) %>%
  dplyr::distinct(SubjectID, .keep_all = TRUE) %>%
  dplyr::select(
    SampleID,
    SubjectID,
    TRG_plot,
    dplyr::all_of(mofa_factor_names)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(mofa_factor_names),
    names_to = "factor",
    values_to = "factor_score"
  ) %>%
  dplyr::mutate(
    factor = factor(
      factor,
      levels = mofa_factor_summary$factor
    ),
    TRG_plot = factor(
      as.character(TRG_plot),
      levels = mofa_heatmap_response_order
    )
  )

mofa_factor_response_all_long <- mofa_factor_scores %>%
  dplyr::filter(!is.na(TRG_plot)) %>%
  dplyr::select(
    SampleID,
    SubjectID,
    Timepoint,
    TRG_plot,
    dplyr::all_of(mofa_factor_names)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(mofa_factor_names),
    names_to = "factor",
    values_to = "factor_score"
  ) %>%
  dplyr::mutate(
    factor = factor(
      factor,
      levels = mofa_factor_summary$factor
    ),
    TRG_plot = factor(
      as.character(TRG_plot),
      levels = mofa_heatmap_response_order
    ),
    Timepoint = factor(
      as.character(Timepoint),
      levels = mofa_heatmap_timepoint_order
    )
  )

p_mofa_factor_response_baseline <- ggplot2::ggplot(
  mofa_factor_response_baseline_long,
  ggplot2::aes(
    x = TRG_plot,
    y = factor_score,
    fill = TRG_plot,
    color = TRG_plot
  )
) +
  ggplot2::geom_violin(
    width = 0.72,
    trim = FALSE,
    alpha = 0.18,
    linewidth = 0.35
  ) +
  ggplot2::geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.42,
    linewidth = 0.42
  ) +
  ggbeeswarm::geom_quasirandom(
    width = 0.11,
    shape = 21,
    size = 1.75,
    stroke = 0.40,
    alpha = 0.90
  ) +
  ggplot2::facet_wrap(
    ~ factor,
    scales = "free_y",
    ncol = 3
  ) +
  ggplot2::scale_fill_manual(
    values = mofa_response_colors,
    drop = FALSE,
    guide = "none"
  ) +
  ggplot2::scale_color_manual(
    values = mofa_response_colors,
    drop = FALSE,
    guide = "none"
  ) +
  ggplot2::labs(
    x = NULL,
    y = "MOFA factor score"
  ) +
  ggplot2::theme_classic(base_size = 9.2) +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "plain", size = 8.7),
    axis.text.x = ggplot2::element_text(face = "plain", size = 8.0),
    plot.margin = ggplot2::margin(4, 5, 4, 4)
  )

p_mofa_factor_response_all <- ggplot2::ggplot(
  mofa_factor_response_all_long,
  ggplot2::aes(
    x = TRG_plot,
    y = factor_score,
    fill = TRG_plot,
    color = TRG_plot,
    shape = Timepoint
  )
) +
  ggplot2::geom_violin(
    width = 0.72,
    trim = FALSE,
    alpha = 0.15,
    linewidth = 0.35
  ) +
  ggplot2::geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.34,
    linewidth = 0.42
  ) +
  ggbeeswarm::geom_quasirandom(
    width = 0.11,
    size = 1.75,
    stroke = 0.42,
    alpha = 0.88
  ) +
  ggplot2::facet_wrap(
    ~ factor,
    scales = "free_y",
    ncol = 3
  ) +
  ggplot2::scale_fill_manual(
    values = mofa_response_colors,
    drop = FALSE,
    name = "Response"
  ) +
  ggplot2::scale_color_manual(
    values = mofa_response_colors,
    drop = FALSE,
    guide = "none"
  ) +
  ggplot2::scale_shape_manual(
    values = c(Before = 21, Ongoing = 24),
    drop = FALSE,
    name = "Time point"
  ) +
  ggplot2::labs(
    x = NULL,
    y = "MOFA factor score"
  ) +
  ggplot2::theme_classic(base_size = 9.2) +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "plain", size = 8.7),
    axis.text.x = ggplot2::element_text(face = "plain", size = 8.0),
    legend.position = "top",
    plot.margin = ggplot2::margin(4, 5, 4, 4)
  )

run_factor_pair_permanova <- function(
    pair_data,
    factor_x,
    factor_y,
    group_column,
    block_column = NULL,
    permutations = mofa_pair_permanova_permutations
) {
  pair_data <- pair_data %>%
    dplyr::filter(
      is.finite(.data[[factor_x]]),
      is.finite(.data[[factor_y]]),
      !is.na(.data[[group_column]])
    )

  if (
    nrow(pair_data) < 6 ||
    dplyr::n_distinct(pair_data[[group_column]]) < 2
  ) {
    return(
      data.frame(
        n_samples = nrow(pair_data),
        pseudo_f = NA_real_,
        r2 = NA_real_,
        permutation_p = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }

  score_matrix <- as.matrix(
    pair_data[
      ,
      c(factor_x, factor_y),
      drop = FALSE
    ]
  )
  score_matrix <- scale(score_matrix)
  score_matrix[!is.finite(score_matrix)] <- 0
  distance_object <- stats::dist(score_matrix)

  model_data <- data.frame(
    pair_group = factor(pair_data[[group_column]])
  )

  permutation_scheme <- permutations

  if (
    !is.null(block_column) &&
    block_column %in% colnames(pair_data)
  ) {
    permutation_scheme <- permute::how(
      nperm = permutations,
      blocks = factor(pair_data[[block_column]])
    )
  }

  permanova_fit <- tryCatch(
    vegan::adonis2(
      distance_object ~ pair_group,
      data = model_data,
      permutations = permutation_scheme,
      by = NULL
    ),
    error = function(e) {
      NULL
    }
  )

  if (is.null(permanova_fit)) {
    return(
      data.frame(
        n_samples = nrow(pair_data),
        pseudo_f = NA_real_,
        r2 = NA_real_,
        permutation_p = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }

  data.frame(
    n_samples = nrow(pair_data),
    pseudo_f = as.numeric(permanova_fit$F[1]),
    r2 = as.numeric(permanova_fit$R2[1]),
    permutation_p = as.numeric(permanova_fit$`Pr(>F)`[1]),
    stringsAsFactors = FALSE
  )
}

mofa_paired_subjects_for_permanova <- mofa_factor_scores %>%
  dplyr::filter(
    as.character(Timepoint) %in% c("Before", "Ongoing")
  ) %>%
  dplyr::count(SubjectID, Timepoint, name = "n") %>%
  dplyr::filter(n == 1) %>%
  dplyr::count(SubjectID, name = "n_timepoints") %>%
  dplyr::filter(n_timepoints == 2) %>%
  dplyr::pull(SubjectID)

mofa_paired_timepoint_data <- mofa_factor_scores %>%
  dplyr::filter(
    SubjectID %in% mofa_paired_subjects_for_permanova,
    as.character(Timepoint) %in% c("Before", "Ongoing")
  )

mofa_pair_permanova_rows <- vector(
  "list",
  ncol(factor_pair_matrix)
)

for (pair_index in seq_len(ncol(factor_pair_matrix))) {
  factor_x <- factor_pair_matrix[1, pair_index]
  factor_y <- factor_pair_matrix[2, pair_index]

  response_result <- run_factor_pair_permanova(
    baseline_factor_data,
    factor_x,
    factor_y,
    group_column = "TRG_plot",
    permutations = mofa_pair_permanova_permutations
  )

  timepoint_result <- run_factor_pair_permanova(
    mofa_paired_timepoint_data,
    factor_x,
    factor_y,
    group_column = "Timepoint",
    block_column = "SubjectID",
    permutations = mofa_pair_permanova_permutations
  )

  mofa_pair_permanova_rows[[pair_index]] <- data.frame(
    pair_id = paste(factor_x, factor_y, sep = "__"),
    factor_x = factor_x,
    factor_y = factor_y,
    response_n = response_result$n_samples,
    response_permanova_f = response_result$pseudo_f,
    response_permanova_r2 = response_result$r2,
    response_permanova_p = response_result$permutation_p,
    timepoint_n = timepoint_result$n_samples,
    timepoint_permanova_f = timepoint_result$pseudo_f,
    timepoint_permanova_r2 = timepoint_result$r2,
    timepoint_permanova_p = timepoint_result$permutation_p,
    stringsAsFactors = FALSE
  )
}

mofa_factor_pair_permanova <- dplyr::bind_rows(
  mofa_pair_permanova_rows
) %>%
  dplyr::mutate(
    response_permanova_fdr = stats::p.adjust(
      response_permanova_p,
      method = "BH"
    ),
    timepoint_permanova_fdr = stats::p.adjust(
      timepoint_permanova_p,
      method = "BH"
    )
  ) %>%
  dplyr::left_join(
    mofa_factor_pair_summary %>%
      dplyr::select(
        pair_id,
        score_correlation,
        view_fingerprint_cosine,
        centroid_distance,
        loocv_lda_auc,
        pair_permutation_p,
        max_statistic_p
      ),
    by = "pair_id"
  ) %>%
  dplyr::left_join(
    mofa_factor_summary %>%
      dplyr::select(
        factor_x = factor,
        factor_x_total_r2 = total_r2
      ),
    by = "factor_x"
  ) %>%
  dplyr::left_join(
    mofa_factor_summary %>%
      dplyr::select(
        factor_y = factor,
        factor_y_total_r2 = total_r2
      ),
    by = "factor_y"
  )

safe_standardize_vector <- function(x) {
  x_mean <- mean(x, na.rm = TRUE)
  x_sd <- stats::sd(x, na.rm = TRUE)

  if (!is.finite(x_sd) || x_sd <= 0) {
    return(rep(0, length(x)))
  }

  output <- (x - x_mean) / x_sd
  output[!is.finite(output)] <- 0
  output
}

mofa_factor_pair_permanova <- mofa_factor_pair_permanova %>%
  dplyr::mutate(
    structural_r2_sum = factor_x_total_r2 + factor_y_total_r2,
    nonredundancy = 1 - abs(score_correlation),
    overall_pair_score =
      safe_standardize_vector(response_permanova_r2) +
      safe_standardize_vector(timepoint_permanova_r2) +
      safe_standardize_vector(structural_r2_sum) +
      safe_standardize_vector(nonredundancy)
  ) %>%
  dplyr::arrange(
    response_permanova_fdr,
    dplyr::desc(response_permanova_r2)
  )

select_best_pair_row <- function(
    result_table,
    p_column,
    r2_column
) {
  candidate_table <- result_table %>%
    dplyr::filter(
      is.finite(.data[[r2_column]])
    ) %>%
    dplyr::arrange(
      .data[[p_column]],
      dplyr::desc(.data[[r2_column]])
    )

  if (nrow(candidate_table) == 0) {
    candidate_table <- result_table %>%
      dplyr::arrange(
        dplyr::desc(.data[[r2_column]])
      )
  }

  candidate_table %>%
    dplyr::slice_head(n = 1)
}

mofa_best_response_pair <- select_best_pair_row(
  mofa_factor_pair_permanova,
  "response_permanova_fdr",
  "response_permanova_r2"
)

mofa_best_timepoint_pair <- select_best_pair_row(
  mofa_factor_pair_permanova,
  "timepoint_permanova_fdr",
  "timepoint_permanova_r2"
)

mofa_best_overall_pair <- mofa_factor_pair_permanova %>%
  dplyr::arrange(dplyr::desc(overall_pair_score)) %>%
  dplyr::slice_head(n = 1)

make_pair_permanova_plot <- function(
    pair_row,
    analysis_type = c("response", "timepoint", "overall")
) {
  analysis_type <- match.arg(analysis_type)
  factor_x <- pair_row$factor_x[1]
  factor_y <- pair_row$factor_y[1]

  if (analysis_type == "response") {
    plot_data <- baseline_factor_data
    fill_variable <- "TRG_plot"
    fill_values <- mofa_response_colors
    legend_name <- "Response"
    plot_subtitle <- paste0(
      "Baseline PERMANOVA: R2 = ",
      sprintf("%.3f", pair_row$response_permanova_r2[1]),
      "; p = ",
      format.pval(pair_row$response_permanova_p[1], digits = 2, eps = 0.001),
      "; BH p = ",
      format.pval(pair_row$response_permanova_fdr[1], digits = 2, eps = 0.001)
    )
  } else if (analysis_type == "timepoint") {
    plot_data <- mofa_paired_timepoint_data
    fill_variable <- "Timepoint"
    fill_values <- timepoint_colors
    legend_name <- "Time point"
    plot_subtitle <- paste0(
      "Subject-blocked PERMANOVA: R2 = ",
      sprintf("%.3f", pair_row$timepoint_permanova_r2[1]),
      "; p = ",
      format.pval(pair_row$timepoint_permanova_p[1], digits = 2, eps = 0.001),
      "; BH p = ",
      format.pval(pair_row$timepoint_permanova_fdr[1], digits = 2, eps = 0.001)
    )
  } else {
    plot_data <- mofa_factor_scores
    fill_variable <- "TRG_plot"
    fill_values <- mofa_response_colors
    legend_name <- "Response"
    plot_subtitle <- paste0(
      "Overall pair score = ",
      sprintf("%.2f", pair_row$overall_pair_score[1]),
      "; score correlation = ",
      sprintf("%.2f", pair_row$score_correlation[1])
    )
  }

  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data[[factor_x]],
      y = .data[[factor_y]],
      fill = .data[[fill_variable]]
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.36,
      linetype = "dashed",
      color = "grey68"
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.36,
      linetype = "dashed",
      color = "grey68"
    )

  if (analysis_type == "timepoint") {
    plot_object <- plot_object +
      ggplot2::geom_path(
        ggplot2::aes(group = SubjectID),
        color = "grey65",
        linewidth = 0.38,
        alpha = 0.62
      ) +
      ggplot2::geom_point(
        ggplot2::aes(shape = Timepoint),
        size = 2.75,
        color = "grey20",
        stroke = 0.58,
        alpha = 0.90
      ) +
      ggplot2::scale_shape_manual(
        values = c(Before = 21, Ongoing = 24),
        guide = "none"
      )
  } else if (analysis_type == "overall") {
    plot_object <- plot_object +
      ggplot2::geom_point(
        ggplot2::aes(shape = Timepoint),
        size = 2.85,
        color = "grey20",
        stroke = 0.58,
        alpha = 0.90
      ) +
      ggplot2::scale_shape_manual(
        values = c(Before = 21, Ongoing = 24),
        name = "Time point"
      )
  } else {
    plot_object <- plot_object +
      ggplot2::geom_point(
        shape = 21,
        size = 2.85,
        color = "grey20",
        stroke = 0.58,
        alpha = 0.90
      )
  }

  plot_object +
    ggplot2::scale_fill_manual(
      values = fill_values,
      drop = FALSE,
      name = legend_name
    ) +
    ggplot2::labs(
      title = paste(factor_x, "vs", factor_y),
      subtitle = plot_subtitle,
      x = make_factor_axis_label(factor_x),
      y = make_factor_axis_label(factor_y)
    ) +
    ggplot2::theme_classic(base_size = 9.4) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "plain", size = 9.6),
      plot.subtitle = ggplot2::element_text(size = 7.7, color = "grey35"),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(4, 5, 4, 4)
    )
}

p_mofa_best_response_pair <- make_pair_permanova_plot(
  mofa_best_response_pair,
  "response"
)

p_mofa_best_timepoint_pair <- make_pair_permanova_plot(
  mofa_best_timepoint_pair,
  "timepoint"
)

p_mofa_best_overall_pair <- make_pair_permanova_plot(
  mofa_best_overall_pair,
  "overall"
)

p_mofa_pair_permanova_summary <-
  p_mofa_best_response_pair |
  p_mofa_best_timepoint_pair |
  p_mofa_best_overall_pair

build_factor_biological_summary_publication <- function(factor_name) {
  factor_r2_data <- mofa_variance_explained %>%
    dplyr::filter(
      as.character(factor) == factor_name
    )

  p_factor_r2 <- ggplot2::ggplot(
    factor_r2_data,
    ggplot2::aes(
      x = view_label,
      y = r2,
      fill = view_label
    )
  ) +
    ggplot2::geom_col(
      width = 0.64,
      color = "grey25",
      linewidth = 0.32
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = scales::percent(r2, accuracy = 0.1)
      ),
      vjust = -0.28,
      size = 2.65
    ) +
    ggplot2::scale_fill_manual(
      values = view_colors,
      guide = "none"
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.13))
    ) +
    ggplot2::labs(
      title = factor_name,
      x = NULL,
      y = "Variance explained"
    ) +
    ggplot2::theme_classic(base_size = 8.8) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 10.0),
      axis.text.x = ggplot2::element_text(
        angle = 28,
        hjust = 1,
        size = 7.2
      ),
      axis.text.y = ggplot2::element_text(size = 7.2),
      axis.title.y = ggplot2::element_text(size = 8.0),
      plot.margin = ggplot2::margin(3, 4, 2, 4)
    )

  factor_clinical_data <- mofa_primary_clinical_plot_data %>%
    dplyr::filter(
      as.character(factor) == factor_name
    )

  p_factor_clinical <- ggplot2::ggplot(
    factor_clinical_data,
    ggplot2::aes(
      x = estimate,
      y = term_label
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      color = "grey70",
      linewidth = 0.40
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = ci_low,
        xmax = ci_high
      ),
      height = 0.12,
      linewidth = 0.55
    ) +
    ggplot2::geom_point(
      shape = 21,
      size = 2.45,
      fill = "white",
      stroke = 0.55
    ) +
    ggplot2::labs(
      x = "Coefficient (95% CI)",
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 8.5) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 7.0),
      axis.text.x = ggplot2::element_text(size = 7.0),
      axis.title.x = ggplot2::element_text(size = 7.8),
      plot.margin = ggplot2::margin(2, 4, 3, 4)
    )

  factor_feature_plot <- build_factor_feature_plot(
    factor_name,
    top_features_per_direction = 3,
    add_annotation = FALSE
  )

  (
    (
      p_factor_r2 /
        p_factor_clinical +
        patchwork::plot_layout(heights = c(0.95, 1.10))
    ) |
      factor_feature_plot$plot
  ) +
    patchwork::plot_layout(widths = c(0.95, 2.75))
}

p_mofa_all_factor_biological_summaries <- list()

for (factor_name in mofa_factor_summary$factor) {
  p_mofa_all_factor_biological_summaries[[factor_name]] <-
    build_factor_biological_summary_publication(factor_name)

  factor_summary_height <- max(
    5.0,
    min(
      8.5,
      4.7 +
        0.08 *
        nrow(mofa_factor_feature_plot_data[[factor_name]])
    )
  )

  if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
    filename = paste0(
      "figures/mofa/internal_v16/MOFA_v16_",
      factor_name,
      "_biological_summary.svg"
    ),
    plot = p_mofa_all_factor_biological_summaries[[factor_name]],
    width = 13.4,
    height = factor_summary_height,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}
}

# Replace the earlier selected-factor composite with the publication-layout
# version while retaining the transparent automatic selection record.
p_mofa_selected_factor_summary <-
  p_mofa_all_factor_biological_summaries[[mofa_selected_factor]]

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename = paste0(
    "figures/mofa/internal_v16/MOFA_v16_selected_factor_summary_",
    mofa_selected_factor,
    ".svg"
  ),
  plot = p_mofa_selected_factor_summary,
  width = 13.4,
  height = 5.8,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename = "figures/mofa/internal_v16/MOFA_v16_factor_score_heatmap_clustered.svg",
  plot = p_mofa_factor_score_heatmap,
  width = 13.2,
  height = 6.2,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename = "figures/mofa/internal_v16/MOFA_v16_factor_response_baseline.svg",
  plot = p_mofa_factor_response_baseline,
  width = 8.8,
  height = 8.2,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename = "figures/mofa/internal_v16/MOFA_v16_factor_response_all_samples.svg",
  plot = p_mofa_factor_response_all,
  width = 8.8,
  height = 8.2,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename = "figures/mofa/internal_v16/MOFA_v16_best_factor_pairs_PERMANOVA.svg",
  plot = p_mofa_pair_permanova_summary,
  width = 15.0,
  height = 5.2,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}



#=================================================================#
# 21. v16 compact publication figures and inferential displays
#=================================================================#
#
# The broad atlas above is retained under figures/mofa/internal_v16.
# This section creates the compact figures intended for manuscript assembly.
#
# Statistical principles:
#   - pCR versus non-pCR pairwise displays: two-sided Wilcoxon rank-sum tests;
#   - Baseline versus After RT pairwise displays: paired Wilcoxon signed-rank tests;
#   - all-time response displays use one mean factor score per SubjectID so repeated
#     observations are not treated as independent;
#   - mixed models and subject-level permutations remain exported as sensitivity analyses;
#   - factor-pair maps use PERMANOVA because the response is two-dimensional;
#   - compact figures include factors with >=2 views at R2 >=2%; 3% is a strict audit.
#-----------------------------------------------------------------#

factor_numeric_order <- function(x) {
  suppressWarnings(
    as.integer(
      gsub(
        "[^0-9]",
        "",
        as.character(x)
      )
    )
  )
}

mofa_publication_factor_audit <-
  mofa_variance_explained %>%
  dplyr::mutate(
    factor = as.character(
      factor
    )
  ) %>%
  dplyr::group_by(
    factor
  ) %>%
  dplyr::summarise(
    n_views_ge_primary = sum(
      r2 >= mofa_active_view_r2,
      na.rm = TRUE
    ),
    n_views_ge_strict = sum(
      r2 >= mofa_strict_active_view_r2,
      na.rm = TRUE
    ),
    total_r2 = sum(
      r2,
      na.rm = TRUE
    ),
    strongest_view_r2 = max(
      r2,
      na.rm = TRUE
    ),
    strongest_view = as.character(
      view[
        which.max(
          r2
        )
      ][1]
    ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    mofa_factor_hierarchical_models %>%
      dplyr::filter(
        model_type ==
          "baseline_response"
      ) %>%
      dplyr::transmute(
        factor = as.character(
          factor
        ),
        baseline_response_effect =
          estimate /
          std_error,
        baseline_response_permutation_p =
          permutation_p
      ),
    by = "factor"
  ) %>%
  dplyr::mutate(
    include_publication =
      n_views_ge_primary >=
      mofa_publication_min_active_views,
    exclusion_reason =
      dplyr::case_when(
        include_publication ~
          NA_character_,
        n_views_ge_primary == 1 ~
          "Only one view reached R2 >= 2%",
        TRUE ~
          "No view pair reached R2 >= 2%"
      )
  )

mofa_publication_factor_order <-
  mofa_publication_factor_audit %>%
  dplyr::filter(
    include_publication
  ) %>%
  dplyr::arrange(
    factor_numeric_order(
      factor
    )
  ) %>%
  dplyr::pull(
    factor
  )

if (
  length(
    mofa_publication_factor_order
  ) < 2
) {
  warning(
    paste0(
      "Fewer than two factors had at least two views with R2 >= 2%. ",
      "The publication display falls back to the four factors with the largest total R2."
    )
  )

  mofa_publication_factor_order <-
    mofa_publication_factor_audit %>%
    dplyr::arrange(
      dplyr::desc(
        n_views_ge_primary
      ),
      dplyr::desc(
        total_r2
      )
    ) %>%
    dplyr::slice_head(
      n = min(
        4,
        dplyr::n()
      )
    ) %>%
    dplyr::arrange(
      factor_numeric_order(
        factor
      )
    ) %>%
    dplyr::pull(
      factor
    )
}

mofa_publication_factor_order <- unique(
  as.character(
    mofa_publication_factor_order
  )
)
mofa_publication_factor_order <-
  mofa_publication_factor_order[
    !is.na(mofa_publication_factor_order) &
      nzchar(mofa_publication_factor_order) &
      mofa_publication_factor_order %in%
        mofa_factor_names
  ]

view_labels_compact <- c(
  Species = "Species",
  `KEGG ortholog` = "KEGG\northolog",
  Metabolite = "Metabolite",
  `Host RNA-seq` = "Host\nRNA-seq"
)

pvalue_class <- function(p) {
  dplyr::case_when(
    !is.finite(p) ~ "Not tested",
    p < 0.001 ~ "p < 0.001",
    p < 0.01 ~ "p < 0.01",
    p < 0.05 ~ "p < 0.05",
    p < 0.10 ~ "0.05 <= p < 0.10",
    TRUE ~ "p >= 0.10"
  )
}

pvalue_symbol <- function(p) {
  dplyr::case_when(
    !is.finite(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

format_permutation_p <- function(p) {
  ifelse(
    is.finite(p),
    format.pval(
      p,
      digits = 2,
      eps = 0.001
    ),
    "NA"
  )
}

first_numeric_or_na <- function(x) {
  if (
    length(x) == 0 ||
    !is.finite(x[1])
  ) {
    return(NA_real_)
  }

  as.numeric(x[1])
}


#-----------------------------------------------------------------#
# 21.1 Feature-retention diagnostic
#-----------------------------------------------------------------#

mofa_v16_feature_target_points <-
  mofa_feature_selection_fraction_table %>%
  dplyr::filter(
    target_cumulative_variance_fraction %in%
      c(
        0.50,
        0.75,
        0.90,
        0.95,
        0.99,
        1.00
      )
  ) %>%
  dplyr::mutate(
    view_label = factor(
      unname(
        view_labels[view]
      ),
      levels = unname(
        view_labels
      )
    )
  )

mofa_v16_selected_feature_points <-
  mofa_feature_selection_summary %>%
  dplyr::transmute(
    view,
    view_label = factor(
      as.character(
        view_label
      ),
      levels = unname(
        view_labels
      )
    ),
    variance_rank =
      as.integer(
        n_features_selected
      )
  ) %>%
  dplyr::left_join(
    mofa_selected_feature_variance %>%
      dplyr::select(
        view,
        variance_rank,
        cumulative_variance_fraction
      ),
    by = c(
      "view",
      "variance_rank"
    )
  )

p_mofa_v16_feature_selection <-
  ggplot2::ggplot(
    mofa_selected_feature_variance,
    ggplot2::aes(
      x = variance_rank,
      y = cumulative_variance_fraction,
      color = view_label
    )
  ) +
  ggplot2::geom_hline(
    yintercept = c(
      0.50,
      0.75,
      0.90,
      0.95,
      0.99
    ),
    linetype = "dotted",
    linewidth = 0.25,
    color = "grey80"
  ) +
  ggplot2::geom_line(
    linewidth = 0.78,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(
    data =
      mofa_v16_feature_target_points,
    ggplot2::aes(
      x = n_features_required,
      y =
        achieved_cumulative_variance_fraction
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 1.8,
    fill = "white",
    color = "grey30",
    stroke = 0.45
  ) +
  ggplot2::geom_point(
    data =
      mofa_v16_selected_feature_points,
    ggplot2::aes(
      x = variance_rank,
      y = cumulative_variance_fraction
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 3.4,
    fill = "#F2C94C",
    color = "grey20",
    stroke = 0.55
  ) +
  ggplot2::geom_text(
    data =
      mofa_v16_selected_feature_points,
    ggplot2::aes(
      x = variance_rank,
      y = cumulative_variance_fraction,
      label = paste0(
        "n = ",
        scales::comma(
          variance_rank
        )
      )
    ),
    inherit.aes = FALSE,
    hjust = -0.08,
    vjust = -0.55,
    size = 2.75
  ) +
  ggplot2::facet_wrap(
    ~ view_label,
    scales = "free_x",
    nrow = 1
  ) +
  ggplot2::scale_color_manual(
    values = view_colors,
    guide = "none"
  ) +
  ggplot2::scale_y_continuous(
    labels =
      scales::label_percent(
        accuracy = 1
      ),
    limits = c(
      0,
      1.02
    ),
    breaks = c(
      0,
      0.50,
      0.75,
      0.90,
      0.95,
      0.99,
      1.00
    ),
    expand =
      ggplot2::expansion(
        mult = c(
          0,
          0
        )
      )
  ) +
  ggplot2::scale_x_continuous(
    labels = scales::comma
  ) +
  ggplot2::labs(
    title =
      "Variance-ranked feature retention before MOFA scaling",
    subtitle = paste0(
      "Open circles: 50%, 75%, 90%, 95%, 99%, and 100%; ",
      "the fifth open circle is 99%. Gold: retained feature count. ",
      "Views with <= ",
      mofa_keep_all_if_n_features_at_most,
      " variable features are retained in full."
    ),
    x =
      "Feature rank by variance",
    y =
      "Cumulative across-feature variance"
  ) +
  ggplot2::theme_classic(
    base_size = 9.8
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    plot.subtitle =
      ggplot2::element_text(
        size = 8.1,
        color = "grey35"
      ),
    strip.background =
      ggplot2::element_blank(),
    strip.text =
      ggplot2::element_text(
        face = "plain",
        size = 9.0
      ),
    axis.text =
      ggplot2::element_text(
        size = 7.8
      )
  )

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename =
    "figures/mofa/publication_v16/MOFA_v16_feature_selection.svg",
  plot =
    p_mofa_v16_feature_selection,
  width = 11.2,
  height = 3.6,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}


#-----------------------------------------------------------------#
# 21.2 Four-omics availability: horizontal and vertical
#-----------------------------------------------------------------#

mofa_v16_availability_data <-
  mofa_availability_long %>%
  dplyr::left_join(
    mofa_sample_metadata %>%
      dplyr::select(
        SampleID,
        SubjectID,
        Timepoint,
        TRG_plot
      ),
    by = "SampleID"
  ) %>%
  dplyr::arrange(
    Timepoint,
    TRG_plot,
    SubjectID,
    SampleID
  ) %>%
  dplyr::mutate(
    SampleID = factor(
      as.character(
        SampleID
      ),
      levels = unique(
        as.character(
          SampleID
        )
      )
    ),
    view_label = factor(
      as.character(
        view_label
      ),
      levels = rev(
        unname(
          view_labels
        )
      )
    )
  )

p_mofa_v16_availability_horizontal <-
  ggplot2::ggplot(
    mofa_v16_availability_data,
    ggplot2::aes(
      x = SampleID,
      y = view_label,
      fill = factor(
        available
      )
    )
  ) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.20
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      `0` = "grey92",
      `1` = "grey25"
    ),
    labels = c(
      `0` = "Missing",
      `1` = "Measured"
    ),
    name = NULL
  ) +
  ggplot2::labs(
    title =
      "Four-omics data availability used for MOFA",
    subtitle = paste0(
      "Samples measured in at least ",
      mofa_min_views_per_sample,
      " views are retained."
    ),
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 10.2
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    plot.subtitle =
      ggplot2::element_text(
        size = 8.4,
        color = "grey35"
      ),
    axis.text.x =
      ggplot2::element_blank(),
    axis.ticks.x =
      ggplot2::element_blank(),
    axis.text.y =
      ggplot2::element_text(
        face = "plain",
        size = 8.8
      ),
    legend.position = "top",
    plot.margin =
      ggplot2::margin(
        3,
        4,
        3,
        4
      )
  )

p_mofa_v16_availability_vertical <-
  ggplot2::ggplot(
    mofa_v16_availability_data,
    ggplot2::aes(
      x = view_label,
      y = SampleID,
      fill = factor(
        available
      )
    )
  ) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.20
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      `0` = "grey92",
      `1` = "grey25"
    ),
    labels = c(
      `0` = "Missing",
      `1` = "Measured"
    ),
    name = NULL
  ) +
  ggplot2::scale_x_discrete(
    labels =
      view_labels_compact
  ) +
  ggplot2::labs(
    title =
      "Four-omics data availability",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 10.0
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    axis.text.x =
      ggplot2::element_text(
        angle = 0,
        hjust = 0.5,
        face = "plain",
        size = 8.5
      ),
    axis.text.y =
      ggplot2::element_blank(),
    axis.ticks.y =
      ggplot2::element_blank(),
    legend.position = "top"
  )

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename =
    "figures/mofa/publication_v16/MOFA_v16_data_availability_horizontal.svg",
  plot =
    p_mofa_v16_availability_horizontal,
  width = 8.0,
  height = 2.25,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename =
    "figures/mofa/publication_v16/MOFA_v16_data_availability_vertical.svg",
  plot =
    p_mofa_v16_availability_vertical,
  width = 3.8,
  height = 6.0,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}


#-----------------------------------------------------------------#
# 21.3 Compact variance-explained overview
#-----------------------------------------------------------------#

mofa_v16_variance_data <-
  mofa_variance_explained %>%
  dplyr::filter(
    as.character(
      factor
    ) %in%
      mofa_publication_factor_order
  ) %>%
  dplyr::mutate(
    factor = factor(
      as.character(
        factor
      ),
      levels = rev(
        mofa_publication_factor_order
      )
    ),
    view_label = factor(
      as.character(
        view_label
      ),
      levels = unname(
        view_labels
      )
    ),
    cell_label = dplyr::case_when(
      !is.finite(r2) ~ "",
      r2 > 0 & r2 < 0.0005 ~ "<0.1%",
      TRUE ~ scales::percent(
        r2,
        accuracy = 0.1
      )
    )
  )

mofa_v16_total_r2 <-
  mofa_variance_total %>%
  dplyr::group_by(
    view,
    view_label
  ) %>%
  dplyr::summarise(
    total_r2 = mean(
      total_r2,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    view_label = factor(
      as.character(
        view_label
      ),
      levels = unname(
        view_labels
      )
    )
  )

p_mofa_v16_total_r2 <-
  ggplot2::ggplot(
    mofa_v16_total_r2,
    ggplot2::aes(
      x = view_label,
      y = total_r2,
      fill = view_label
    )
  ) +
  ggplot2::geom_col(
    width = 0.66,
    color = "grey25",
    linewidth = 0.35
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label =
        scales::percent(
          total_r2,
          accuracy = 0.1
        )
    ),
    vjust = -0.30,
    size = 2.8
  ) +
  ggplot2::scale_fill_manual(
    values = view_colors,
    guide = "none"
  ) +
  ggplot2::scale_x_discrete(
    labels =
      view_labels_compact
  ) +
  ggplot2::scale_y_continuous(
    labels =
      scales::label_percent(
        accuracy = 1
      ),
    expand =
      ggplot2::expansion(
        mult = c(
          0,
          0.15
        )
      )
  ) +
  ggplot2::labs(
    title =
      "Total model variance explained",
    x = NULL,
    y = "Total R²"
  ) +
  ggplot2::theme_classic(
    base_size = 9.2
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 9.8
      ),
    axis.text.x =
      ggplot2::element_text(
        angle = 0,
        hjust = 0.5,
        size = 8.0
      )
  )

p_mofa_v16_variance_heatmap <-
  ggplot2::ggplot(
    mofa_v16_variance_data,
    ggplot2::aes(
      x = view_label,
      y = factor,
      fill = r2
    )
  ) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.32
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = cell_label
    ),
    size = 2.45
  ) +
  ggplot2::scale_fill_gradient(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "white",
    high = "grey20",
    labels =
      scales::label_percent(
        accuracy = 1
      ),
    name = "Factor R²"
  ) +
  ggplot2::scale_x_discrete(
    labels =
      view_labels_compact
  ) +
  ggplot2::labs(
    title =
      "Variance explained by latent factor and omics view",
    subtitle = paste0(
      "Active-view display threshold: R² >= ",
      scales::percent(
        mofa_active_view_r2,
        accuracy = 1
      ),
      ". View dominance is descriptive and does not exclude a factor."
    ),
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 9.2
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 10.0
      ),
    plot.subtitle =
      ggplot2::element_text(
        size = 7.7,
        color = "grey35"
      ),
    axis.text.x =
      ggplot2::element_text(
        angle = 0,
        hjust = 0.5,
        face = "plain",
        size = 8.0
      ),
    axis.text.y =
      ggplot2::element_text(
        face = "plain",
        size = 7.8
      ),
    axis.ticks =
      ggplot2::element_blank()
  )

p_mofa_v16_variance_overview <-
  p_mofa_v16_total_r2 /
  p_mofa_v16_variance_heatmap +
  patchwork::plot_layout(
    heights = c(
      0.34,
      1
    )
  )

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename =
    "figures/mofa/publication_v16/MOFA_v16_variance_overview_compact.svg",
  plot =
    p_mofa_v16_variance_overview,
  width = 5.3,
  height = max(
    4.8,
    2.5 +
      0.28 *
      length(
        mofa_publication_factor_order
      )
  ),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}


#-----------------------------------------------------------------#
# 21.4 Biological-question-specific factor associations
#-----------------------------------------------------------------#
# Tile colour retains the prespecified model coefficient / SE so the direction
# and standardized model signal remain visible. Point symbols use the primary
# non-parametric tests: rank-sum for response contrasts and signed-rank for the
# paired timepoint contrast.

mofa_v16_association_data <-
  mofa_factor_hierarchical_models %>%
  dplyr::filter(
    factor %in%
      mofa_publication_factor_order,
    model_type %in%
      c(
        "baseline_response",
        "ongoing_response",
        "paired_overall_change"
      )
  ) %>%
  dplyr::mutate(
    standardized_effect = dplyr::if_else(
      is.finite(std_error) &
        std_error > 0,
      estimate / std_error,
      NA_real_
    ),
    p_class = pvalue_class(
      wilcoxon_p
    ),
    p_symbol = pvalue_symbol(
      wilcoxon_p
    ),
    question = dplyr::recode(
      model_type,
      baseline_response =
        "pCR vs non-pCR\nat Baseline",
      ongoing_response =
        "pCR vs non-pCR\nafter RT",
      paired_overall_change =
        "After RT vs Baseline\npaired change"
    ),
    question = factor(
      question,
      levels = c(
        "pCR vs non-pCR\nat Baseline",
        "pCR vs non-pCR\nafter RT",
        "After RT vs Baseline\npaired change"
      )
    ),
    factor = factor(
      factor,
      levels = rev(
        mofa_publication_factor_order
      )
    )
  )

mofa_v16_effect_limit <- max(
  abs(
    mofa_v16_association_data$standardized_effect
  ),
  na.rm = TRUE
)

if (
  !is.finite(mofa_v16_effect_limit) ||
  mofa_v16_effect_limit <= 0
) {
  mofa_v16_effect_limit <- 1
}

p_mofa_v16_factor_associations <-
  ggplot2::ggplot(
    mofa_v16_association_data,
    ggplot2::aes(
      x = question,
      y = factor,
      fill = standardized_effect
    )
  ) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.32
  ) +
  ggplot2::geom_point(
    data = mofa_v16_association_data %>%
      dplyr::filter(
        !is.finite(wilcoxon_p) |
          wilcoxon_p >= 0.10
      ),
    shape = 21,
    size = 1.85,
    fill = "white",
    color = "grey30",
    stroke = 0.42
  ) +
  ggplot2::geom_point(
    data = mofa_v16_association_data %>%
      dplyr::filter(
        is.finite(wilcoxon_p),
        wilcoxon_p >= 0.05,
        wilcoxon_p < 0.10
      ),
    shape = 21,
    size = 1.95,
    fill = "grey72",
    color = "grey35",
    stroke = 0.42
  ) +
  ggplot2::geom_point(
    data = mofa_v16_association_data %>%
      dplyr::filter(
        is.finite(wilcoxon_p),
        wilcoxon_p < 0.05
      ),
    shape = 21,
    size = 2.55,
    fill = "white",
    color = "black",
    stroke = 0.62
  ) +
  ggplot2::geom_text(
    data = mofa_v16_association_data %>%
      dplyr::filter(
        is.finite(wilcoxon_p),
        wilcoxon_p < 0.05
      ),
    ggplot2::aes(
      label = p_symbol
    ),
    color = "black",
    size = 2.0,
    vjust = 0.63
  ) +
  ggplot2::scale_fill_gradient2(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "#3B75AF",
    mid = "white",
    high = "#C44E52",
    midpoint = 0,
    limits = c(
      -mofa_v16_effect_limit,
      mofa_v16_effect_limit
    ),
    oob = scales::squish,
    name = "Model\ncoefficient / SE"
  ) +
  ggplot2::labs(
    title =
      "MOFA factor associations by biological question",
    subtitle = paste0(
      "Open circle: Wilcoxon p >= 0.10; grey circle: 0.05 <= p < 0.10; ",
      "open circle: p < 0.05. Response contrasts use rank-sum tests; ",
      "paired change uses the signed-rank test."
    ),
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 8.8
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 9.8
    ),
    plot.subtitle = ggplot2::element_text(
      size = 7.1,
      color = "grey35"
    ),
    axis.text.x = ggplot2::element_text(
      hjust = 0.5,
      size = 7.2
    ),
    axis.text.y = ggplot2::element_text(
      face = "plain",
      size = 7.5
    ),
    legend.key.height = grid::unit(
      1.2,
      "lines"
    ),
    plot.margin = ggplot2::margin(
      4,
      4,
      4,
      4
    )
  )

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename =
    "figures/mofa/publication_v16/MOFA_v16_factor_associations_main.svg",
  plot = p_mofa_v16_factor_associations,
  width = 4.45,
  height = max(
    4.2,
    2.0 +
      0.27 *
      length(
        mofa_publication_factor_order
      )
  ),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}


#-----------------------------------------------------------------#
# 21.5 All-factor pairwise distribution plots
#-----------------------------------------------------------------#

add_factor_annotation_position <- function(
    plot_data
) {
  plot_data %>%
    dplyr::group_by(
      factor
    ) %>%
    dplyr::summarise(
      y_max = max(
        factor_score,
        na.rm = TRUE
      ),
      y_min = min(
        factor_score,
        na.rm = TRUE
      ),
      y_position = y_max +
        0.13 *
        pmax(
          y_max - y_min,
          1
        ),
      .groups = "drop"
    )
}

mofa_v16_overall_response_long <-
  mofa_factor_scores %>%
  dplyr::filter(
    !is.na(TRG_plot)
  ) %>%
  dplyr::group_by(
    SubjectID,
    TRG_plot
  ) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(
        mofa_factor_names
      ),
      ~ mean(
        .x,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      mofa_factor_names
    ),
    names_to = "factor",
    values_to = "factor_score"
  ) %>%
  dplyr::filter(
    factor %in%
      mofa_publication_factor_order
  ) %>%
  dplyr::mutate(
    factor = factor(
      factor,
      levels = mofa_publication_factor_order
    ),
    TRG_plot = factor(
      as.character(TRG_plot),
      levels = c(
        "pCR",
        "non_pCR"
      )
    )
  )

mofa_v16_baseline_long <-
  mofa_factor_scores %>%
  dplyr::filter(
    as.character(Timepoint) == "Before",
    !is.na(TRG_plot)
  ) %>%
  dplyr::distinct(
    SubjectID,
    .keep_all = TRUE
  ) %>%
  dplyr::select(
    SubjectID,
    TRG_plot,
    dplyr::all_of(
      mofa_factor_names
    )
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      mofa_factor_names
    ),
    names_to = "factor",
    values_to = "factor_score"
  ) %>%
  dplyr::filter(
    factor %in%
      mofa_publication_factor_order
  ) %>%
  dplyr::mutate(
    factor = factor(
      factor,
      levels = mofa_publication_factor_order
    ),
    TRG_plot = factor(
      as.character(TRG_plot),
      levels = c(
        "pCR",
        "non_pCR"
      )
    )
  )

mofa_v16_after_rt_long <-
  mofa_factor_scores %>%
  dplyr::filter(
    as.character(Timepoint) == "Ongoing",
    !is.na(TRG_plot)
  ) %>%
  dplyr::distinct(
    SubjectID,
    .keep_all = TRUE
  ) %>%
  dplyr::select(
    SubjectID,
    TRG_plot,
    dplyr::all_of(
      mofa_factor_names
    )
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      mofa_factor_names
    ),
    names_to = "factor",
    values_to = "factor_score"
  ) %>%
  dplyr::filter(
    factor %in%
      mofa_publication_factor_order
  ) %>%
  dplyr::mutate(
    factor = factor(
      factor,
      levels = mofa_publication_factor_order
    ),
    TRG_plot = factor(
      as.character(TRG_plot),
      levels = c(
        "pCR",
        "non_pCR"
      )
    )
  )

mofa_v16_paired_subject_ids <-
  mofa_factor_scores %>%
  dplyr::filter(
    as.character(Timepoint) %in%
      c(
        "Before",
        "Ongoing"
      )
  ) %>%
  dplyr::distinct(
    SubjectID,
    Timepoint
  ) %>%
  dplyr::count(
    SubjectID,
    name = "n_timepoints"
  ) %>%
  dplyr::filter(
    n_timepoints == 2
  ) %>%
  dplyr::pull(
    SubjectID
  )

mofa_v16_paired_long <-
  mofa_factor_scores %>%
  dplyr::filter(
    SubjectID %in%
      mofa_v16_paired_subject_ids,
    as.character(Timepoint) %in%
      c(
        "Before",
        "Ongoing"
      )
  ) %>%
  dplyr::select(
    SubjectID,
    Timepoint,
    TRG_plot,
    dplyr::all_of(
      mofa_factor_names
    )
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      mofa_factor_names
    ),
    names_to = "factor",
    values_to = "factor_score"
  ) %>%
  dplyr::filter(
    factor %in%
      mofa_publication_factor_order
  ) %>%
  dplyr::mutate(
    factor = factor(
      factor,
      levels = mofa_publication_factor_order
    ),
    Timepoint_display = factor(
      dplyr::recode(
        as.character(Timepoint),
        Before = "Baseline",
        Ongoing = "After RT"
      ),
      levels = c(
        "Baseline",
        "After RT"
      )
    )
  )

make_wilcoxon_annotation <- function(
    plot_data,
    comparison_name,
    prefix = "Wilcoxon p = "
) {
  add_factor_annotation_position(
    plot_data
  ) %>%
    dplyr::left_join(
      mofa_factor_wilcoxon_tests %>%
        dplyr::filter(
          comparison ==
            comparison_name
        ) %>%
        dplyr::select(
          factor,
          wilcoxon_p
        ),
      by = "factor"
    ) %>%
    dplyr::mutate(
      label = paste0(
        prefix,
        format_permutation_p(
          wilcoxon_p
        )
      )
    )
}

mofa_v16_overall_annotation <-
  make_wilcoxon_annotation(
    mofa_v16_overall_response_long,
    "overall_response_subject_mean"
  )

mofa_v16_baseline_annotation <-
  make_wilcoxon_annotation(
    mofa_v16_baseline_long,
    "baseline_response"
  )

mofa_v16_after_rt_annotation <-
  make_wilcoxon_annotation(
    mofa_v16_after_rt_long,
    "ongoing_response"
  )

mofa_v16_paired_annotation <-
  make_wilcoxon_annotation(
    mofa_v16_paired_long,
    "paired_overall_change",
    prefix = "Paired Wilcoxon p = "
  )

build_response_violin <- function(
    plot_data,
    annotation_data,
    title_text
) {
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = TRG_plot,
      y = factor_score,
      fill = TRG_plot
    )
  ) +
    ggplot2::geom_violin(
      width = 0.42,
      trim = TRUE,
      alpha = 0.25,
      color = "grey35",
      linewidth = 0.30
    ) +
    ggplot2::geom_boxplot(
      width = 0.12,
      outlier.shape = NA,
      alpha = 0.42,
      linewidth = 0.32
    ) +
    ggbeeswarm::geom_quasirandom(
      width = 0.07,
      shape = 21,
      size = 1.35,
      color = "grey25",
      stroke = 0.30,
      alpha = 0.90
    ) +
    ggplot2::geom_text(
      data = annotation_data,
      ggplot2::aes(
        x = 1.5,
        y = y_position,
        label = label
      ),
      inherit.aes = FALSE,
      size = 2.15
    ) +
    ggplot2::facet_wrap(
      ~ factor,
      scales = "free_y",
      ncol = 4
    ) +
    ggplot2::scale_fill_manual(
      values = mofa_response_colors[
        c(
          "pCR",
          "non_pCR"
        )
      ],
      guide = "none"
    ) +
    ggplot2::coord_cartesian(
      clip = "off"
    ) +
    ggplot2::labs(
      title = title_text,
      x = NULL,
      y = "MOFA factor score"
    ) +
    ggplot2::theme_classic(
      base_size = 8.4
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 9.3
      ),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "plain",
        size = 7.8
      ),
      axis.text.x = ggplot2::element_text(
        face = "plain",
        size = 7.0
      ),
      panel.spacing = grid::unit(
        0.48,
        "lines"
      ),
      plot.margin = ggplot2::margin(
        4,
        4,
        4,
        4
      )
    )
}

p_mofa_v16_overall_response <-
  build_response_violin(
    mofa_v16_overall_response_long,
    mofa_v16_overall_annotation,
    "pCR versus non-pCR: subject mean across available timepoints"
  )

p_mofa_v16_baseline_response <-
  build_response_violin(
    mofa_v16_baseline_long,
    mofa_v16_baseline_annotation,
    "pCR versus non-pCR at Baseline"
  )

p_mofa_v16_after_rt_response <-
  build_response_violin(
    mofa_v16_after_rt_long,
    mofa_v16_after_rt_annotation,
    "pCR versus non-pCR after RT"
  )

p_mofa_v16_paired_change <-
  ggplot2::ggplot(
    mofa_v16_paired_long,
    ggplot2::aes(
      x = Timepoint_display,
      y = factor_score,
      group = SubjectID
    )
  ) +
  ggplot2::geom_violin(
    ggplot2::aes(
      fill = Timepoint_display,
      group = Timepoint_display
    ),
    width = 0.42,
    trim = TRUE,
    alpha = 0.18,
    color = "grey40",
    linewidth = 0.30
  ) +
  ggplot2::geom_line(
    ggplot2::aes(
      color = TRG_plot
    ),
    linewidth = 0.31,
    alpha = 0.52
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      fill = TRG_plot
    ),
    shape = 21,
    size = 1.30,
    color = "grey25",
    stroke = 0.28
  ) +
  ggplot2::geom_text(
    data = mofa_v16_paired_annotation,
    ggplot2::aes(
      x = 1.5,
      y = y_position,
      label = label
    ),
    inherit.aes = FALSE,
    size = 2.12
  ) +
  ggplot2::facet_wrap(
    ~ factor,
    scales = "free_y",
    ncol = 4
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      timepoint_display_colors,
      mofa_response_colors[
        c(
          "pCR",
          "non_pCR"
        )
      ]
    ),
    guide = "none"
  ) +
  ggplot2::scale_color_manual(
    values = mofa_response_colors[
      c(
        "pCR",
        "non_pCR"
      )
    ],
    guide = "none"
  ) +
  ggplot2::coord_cartesian(
    clip = "off"
  ) +
  ggplot2::labs(
    title = "Within-subject treatment change",
    subtitle = paste0(
      "Only subjects measured at both Baseline and After RT are shown; ",
      "the signed-rank test evaluates whether the paired After RT - Baseline change is centred at zero."
    ),
    x = NULL,
    y = "MOFA factor score"
  ) +
  ggplot2::theme_classic(
    base_size = 8.4
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 9.3
    ),
    plot.subtitle = ggplot2::element_text(
      size = 7.2,
      color = "grey35"
    ),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(
      face = "plain",
      size = 7.8
    ),
    axis.text.x = ggplot2::element_text(
      size = 6.9
    ),
    panel.spacing = grid::unit(
      0.48,
      "lines"
    )
  )

mofa_v16_distribution_height <- max(
  5.2,
  1.55 +
    1.28 *
    ceiling(
      length(
        mofa_publication_factor_order
      ) /
        4
    )
)

for (
  plot_spec in list(
    list(
      filename = "MOFA_v16_factor_response_overall_subject_mean.svg",
      plot = p_mofa_v16_overall_response
    ),
    list(
      filename = "MOFA_v16_factor_response_baseline.svg",
      plot = p_mofa_v16_baseline_response
    ),
    list(
      filename = "MOFA_v16_factor_response_after_RT.svg",
      plot = p_mofa_v16_after_rt_response
    ),
    list(
      filename = "MOFA_v16_factor_paired_change.svg",
      plot = p_mofa_v16_paired_change
    )
  )
) {
  if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
    filename = paste0(
      "figures/mofa/publication_v16/",
      plot_spec$filename
    ),
    plot = plot_spec$plot,
    width = 6.35,
    height = mofa_v16_distribution_height,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}
}


#-----------------------------------------------------------------#
# 21.6 Subject repeatability, not an automatic exclusion criterion
#-----------------------------------------------------------------#

mofa_v16_subject_repeatability <-
  mofa_factor_subject_specificity %>%
  dplyr::mutate(
    factor = factor(
      factor,
      levels =
        mofa_publication_factor_order
    ),
    repeatability_class =
      dplyr::case_when(
        !is.finite(
          icc
        ) ~
          "Not estimable",
        icc >= 0.75 ~
          "High repeatability",
        icc >= 0.50 ~
          "Moderate repeatability",
        TRUE ~
          "Low repeatability"
      )
  )

p_mofa_v16_subject_repeatability <-
  ggplot2::ggplot(
    mofa_v16_subject_repeatability,
    ggplot2::aes(
      x = factor,
      y = icc,
      fill = repeatability_class
    )
  ) +
  ggplot2::geom_hline(
    yintercept = 0.50,
    linetype = "dashed",
    linewidth = 0.48,
    color = "grey45"
  ) +
  ggplot2::geom_col(
    width = 0.66,
    color = "grey25",
    linewidth = 0.32
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = ifelse(
        is.finite(
          icc
        ),
        sprintf(
          "%.2f",
          icc
        ),
        "NA"
      )
    ),
    vjust = -0.32,
    size = 2.65
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      `Low repeatability` =
        "grey82",
      `Moderate repeatability` =
        "grey58",
      `High repeatability` =
        "grey25",
      `Not estimable` =
        "white"
    ),
    name = NULL
  ) +
  ggplot2::scale_y_continuous(
    limits = c(
      0,
      1.05
    ),
    breaks = seq(
      0,
      1,
      by = 0.25
    ),
    expand =
      ggplot2::expansion(
        mult = c(
          0,
          0
        )
      )
  ) +
  ggplot2::labs(
    title =
      "Repeatability of MOFA factor scores across time",
    subtitle = paste0(
      "ICC >= 0.50 indicates a stable between-subject component; ",
      "it does not demonstrate that one subject drives the factor and is not an exclusion rule."
    ),
    x = NULL,
    y = "Subject-level ICC"
  ) +
  ggplot2::theme_classic(
    base_size = 9.4
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold"
      ),
    plot.subtitle =
      ggplot2::element_text(
        size = 7.9,
        color = "grey35"
      ),
    axis.text.x =
      ggplot2::element_text(
        face = "plain",
        size = 7.5
      ),
    legend.position = "top"
  )

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename =
    "figures/mofa/publication_v16/MOFA_v16_factor_subject_repeatability.svg",
  plot =
    p_mofa_v16_subject_repeatability,
  width = 7.2,
  height = 3.8,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}


#-----------------------------------------------------------------#
# 21.7 Response-selected factor-pair maps
#-----------------------------------------------------------------#
# The primary response map uses one subject-level mean coordinate per subject,
# thereby incorporating every available timepoint without treating repeated
# observations as independent. A separate Baseline-only map is retained.

square_pair_limits <- function(
    x,
    y
) {
  limit_value <- max(
    abs(
      c(
        x,
        y
      )
    ),
    na.rm = TRUE
  )

  if (!is.finite(limit_value) || limit_value <= 0) {
    limit_value <- 1
  }

  c(
    -1.08 * limit_value,
    1.08 * limit_value
  )
}

convex_hull_data <- function(
    data,
    x_name,
    y_name,
    group_name
) {
  data %>%
    dplyr::filter(
      is.finite(.data[[x_name]]),
      is.finite(.data[[y_name]]),
      !is.na(.data[[group_name]])
    ) %>%
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(group_name)
      )
    ) %>%
    dplyr::group_modify(
      function(.x, .y) {
        if (nrow(.x) < 3) {
          return(
            .x[0, , drop = FALSE]
          )
        }

        hull_index <- grDevices::chull(
          .x[[x_name]],
          .x[[y_name]]
        )

        .x[
          c(
            hull_index,
            hull_index[1]
          ),
          ,
          drop = FALSE
        ]
      }
    ) %>%
    dplyr::ungroup()
}

factor_axis_label_v16 <- function(
    factor_name
) {
  factor_view <- mofa_variance_explained %>%
    dplyr::filter(
      as.character(factor) ==
        factor_name
    ) %>%
    dplyr::slice_max(
      order_by = r2,
      n = 1,
      with_ties = FALSE
    )

  if (nrow(factor_view) == 0) {
    return(factor_name)
  }

  paste0(
    factor_name,
    " (",
    as.character(factor_view$view_label[1]),
    " ",
    scales::percent(
      factor_view$r2[1],
      accuracy = 0.1
    ),
    ")"
  )
}

mofa_v16_subject_mean_factor_data <-
  mofa_factor_scores %>%
  dplyr::filter(
    !is.na(TRG_plot)
  ) %>%
  dplyr::group_by(
    SubjectID,
    TRG_plot
  ) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(
        mofa_factor_names
      ),
      ~ mean(
        .x,
        na.rm = TRUE
      )
    ),
    n_timepoints = dplyr::n(),
    .groups = "drop"
  )

mofa_v16_baseline_factor_data <-
  mofa_factor_scores %>%
  dplyr::filter(
    as.character(Timepoint) == "Before",
    !is.na(TRG_plot)
  ) %>%
  dplyr::distinct(
    SubjectID,
    .keep_all = TRUE
  )

calculate_factor_pair_permanova <- function(
    analysis_data,
    scope_label
) {
  pair_list <- utils::combn(
    mofa_publication_factor_order,
    2,
    simplify = FALSE
  )

  dplyr::bind_rows(
    lapply(
      pair_list,
      function(pair_factors) {
        pair_data <- analysis_data %>%
          dplyr::select(
            SubjectID,
            TRG_plot,
            dplyr::all_of(
              pair_factors
            )
          ) %>%
          dplyr::filter(
            is.finite(.data[[pair_factors[1]]]),
            is.finite(.data[[pair_factors[2]]]),
            !is.na(TRG_plot)
          )

        if (
          nrow(pair_data) < 8 ||
          dplyr::n_distinct(pair_data$TRG_plot) != 2 ||
          min(table(pair_data$TRG_plot)) < 3
        ) {
          return(
            data.frame(
              scope = scope_label,
              pair_id = paste(
                pair_factors,
                collapse = "__"
              ),
              factor_x = pair_factors[1],
              factor_y = pair_factors[2],
              n_subjects = nrow(pair_data),
              permanova_f = NA_real_,
              permanova_r2 = NA_real_,
              permanova_p = NA_real_,
              stringsAsFactors = FALSE
            )
          )
        }

        pair_matrix <- scale(
          as.matrix(
            pair_data[
              ,
              pair_factors,
              drop = FALSE
            ]
          )
        )
        pair_matrix[!is.finite(pair_matrix)] <- 0

        set.seed(
          20261001 +
            match(
              paste(pair_factors, collapse = "__"),
              vapply(
                pair_list,
                paste,
                collapse = "__",
                FUN.VALUE = character(1)
              )
            ) +
            ifelse(
              scope_label == "Overall subject mean",
              0,
              1000
            )
        )

        permanova_result <- tryCatch(
          vegan::adonis2(
            stats::dist(pair_matrix) ~ TRG_plot,
            data = pair_data,
            permutations =
              mofa_pair_permanova_permutations
          ),
          error = function(e) {
            NULL
          }
        )

        data.frame(
          scope = scope_label,
          pair_id = paste(
            pair_factors,
            collapse = "__"
          ),
          factor_x = pair_factors[1],
          factor_y = pair_factors[2],
          n_subjects = nrow(pair_data),
          permanova_f = if (
            is.null(permanova_result)
          ) {
            NA_real_
          } else {
            as.numeric(permanova_result$F[1])
          },
          permanova_r2 = if (
            is.null(permanova_result)
          ) {
            NA_real_
          } else {
            as.numeric(permanova_result$R2[1])
          },
          permanova_p = if (
            is.null(permanova_result)
          ) {
            NA_real_
          } else {
            as.numeric(permanova_result$`Pr(>F)`[1])
          },
          stringsAsFactors = FALSE
        )
      }
    )
  ) %>%
    dplyr::group_by(
      scope
    ) %>%
    dplyr::mutate(
      permanova_fdr = stats::p.adjust(
        permanova_p,
        method = "BH"
      )
    ) %>%
    dplyr::ungroup()
}

mofa_v16_pair_permanova_atlas <- dplyr::bind_rows(
  calculate_factor_pair_permanova(
    mofa_v16_subject_mean_factor_data,
    "Overall subject mean"
  ),
  calculate_factor_pair_permanova(
    mofa_v16_baseline_factor_data,
    "Baseline only"
  )
)

mofa_v30_pair_permanova_atlas <-
  mofa_v16_pair_permanova_atlas

mofa_v16_selected_pair_rows <-
  mofa_v16_pair_permanova_atlas %>%
  dplyr::filter(
    is.finite(permanova_p)
  ) %>%
  dplyr::group_by(
    scope
  ) %>%
  dplyr::arrange(
    permanova_p,
    dplyr::desc(permanova_r2),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(
    n = 1
  ) %>%
  dplyr::ungroup()

build_response_pair_map_v16 <- function(
    pair_row
) {
  factor_x <- as.character(
    pair_row$factor_x[1]
  )
  factor_y <- as.character(
    pair_row$factor_y[1]
  )
  scope_label <- as.character(
    pair_row$scope[1]
  )

  source_data <- if (
    scope_label == "Overall subject mean"
  ) {
    mofa_v16_subject_mean_factor_data
  } else {
    mofa_v16_baseline_factor_data
  }

  pair_data <- source_data %>%
    dplyr::select(
      SubjectID,
      TRG_plot,
      dplyr::all_of(
        c(
          factor_x,
          factor_y
        )
      )
    ) %>%
    dplyr::filter(
      is.finite(.data[[factor_x]]),
      is.finite(.data[[factor_y]]),
      !is.na(TRG_plot)
    ) %>%
    dplyr::mutate(
      factor_x_z = as.numeric(
        scale(.data[[factor_x]])
      ),
      factor_y_z = as.numeric(
        scale(.data[[factor_y]])
      ),
      TRG_plot = factor(
        as.character(TRG_plot),
        levels = c(
          "pCR",
          "non_pCR"
        )
      )
    )

  pair_hull <- convex_hull_data(
    pair_data,
    "factor_x_z",
    "factor_y_z",
    "TRG_plot"
  )

  pair_limits <- square_pair_limits(
    pair_data$factor_x_z,
    pair_data$factor_y_z
  )

  ggplot2::ggplot(
    pair_data,
    ggplot2::aes(
      x = factor_x_z,
      y = factor_y_z
    )
  ) +
  ggplot2::geom_polygon(
    data = pair_hull,
    ggplot2::aes(
      group = TRG_plot,
      fill = TRG_plot,
      color = TRG_plot
    ),
    alpha = 0.11,
    linewidth = 0.58,
    show.legend = FALSE
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.32,
    linetype = "dashed",
    color = "grey72"
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.32,
    linetype = "dashed",
    color = "grey72"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      fill = TRG_plot
    ),
    shape = 21,
    size = 2.45,
    color = "grey20",
    stroke = 0.50,
    alpha = 0.92
  ) +
  ggplot2::stat_summary(
    ggplot2::aes(
      fill = TRG_plot
    ),
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3.9,
    color = "black",
    stroke = 0.72
  ) +
  ggplot2::scale_fill_manual(
    values = mofa_response_colors[
      c(
        "pCR",
        "non_pCR"
      )
    ],
    name = "Response"
  ) +
  ggplot2::scale_color_manual(
    values = mofa_response_colors[
      c(
        "pCR",
        "non_pCR"
      )
    ],
    guide = "none"
  ) +
  ggplot2::scale_x_continuous(
    limits = pair_limits
  ) +
  ggplot2::scale_y_continuous(
    limits = pair_limits
  ) +
  ggplot2::coord_fixed(
    ratio = 1
  ) +
  ggplot2::labs(
    title = paste0(
      factor_x,
      " vs ",
      factor_y
    ),
    subtitle = paste0(
      scope_label,
      "; PERMANOVA R² = ",
      sprintf(
        "%.2f",
        pair_row$permanova_r2[1]
      ),
      ", p = ",
      format_permutation_p(
        pair_row$permanova_p[1]
      ),
      "; n = ",
      pair_row$n_subjects[1],
      " subjects"
    ),
    x = paste0(
      factor_axis_label_v16(factor_x),
      " (z-score)"
    ),
    y = paste0(
      factor_axis_label_v16(factor_y),
      " (z-score)"
    )
  ) +
  ggplot2::theme_classic(
    base_size = 9.0
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 9.8
    ),
    plot.subtitle = ggplot2::element_text(
      size = 7.3,
      color = "grey35"
    ),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(
      4,
      4,
      4,
      4
    )
  )
}

p_mofa_v16_response_pair_maps <- lapply(
  seq_len(
    nrow(mofa_v16_selected_pair_rows)
  ),
  function(pair_index) {
    build_response_pair_map_v16(
      mofa_v16_selected_pair_rows[
        pair_index,
        ,
        drop = FALSE
      ]
    )
  }
)

if (length(p_mofa_v16_response_pair_maps) > 0) {
  p_mofa_v16_selected_response_pairs <-
    patchwork::wrap_plots(
      p_mofa_v16_response_pair_maps,
      ncol = min(
        2,
        length(p_mofa_v16_response_pair_maps)
      ),
      guides = "collect"
    ) +
    patchwork::plot_annotation(
      title = "Response-selected MOFA factor-pair maps",
      subtitle = paste0(
        "The overall map uses one mean coordinate per subject across available timepoints; ",
        "the second map uses Baseline only. Pair screening remains exploratory."
      )
    ) &
    ggplot2::theme(
      legend.position = "bottom"
    )

  if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
    filename = paste0(
      "figures/mofa/publication_v16/",
      "MOFA_v16_selected_response_pair_maps.svg"
    ),
    plot = p_mofa_v16_selected_response_pairs,
    width = 4.35 *
      min(
        2,
        length(p_mofa_v16_response_pair_maps)
      ),
    height = 4.45,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}
}


#-----------------------------------------------------------------#
# 21.8 Compact selected-factor summary
#-----------------------------------------------------------------#

mofa_v16_selected_factor_candidates <-
  mofa_factor_wilcoxon_tests %>%
  dplyr::filter(
    comparison ==
      "baseline_response",
    factor %in%
      mofa_publication_factor_order
  ) %>%
  dplyr::arrange(
    wilcoxon_p,
    dplyr::desc(
      abs(
        wilcoxon_effect
      )
    )
  )

mofa_v16_selected_factor <- if (
  mofa_selected_factor %in%
    mofa_publication_factor_order
) {
  as.character(
    mofa_selected_factor
  )
} else if (
  nrow(
    mofa_v16_selected_factor_candidates
  ) > 0
) {
  as.character(
    mofa_v16_selected_factor_candidates$factor[1]
  )
} else {
  as.character(
    mofa_publication_factor_order[1]
  )
}

mofa_v16_selected_factor_r2 <-
  mofa_variance_explained %>%
  dplyr::filter(
    as.character(
      factor
    ) ==
      mofa_v16_selected_factor
  ) %>%
  dplyr::mutate(
    view_label = factor(
      as.character(
        view_label
      ),
      levels = rev(
        unname(
          view_labels
        )
      )
    )
  )

p_mofa_v16_selected_factor_r2 <-
  ggplot2::ggplot(
    mofa_v16_selected_factor_r2,
    ggplot2::aes(
      x = r2,
      y = view_label,
      fill = view_label
    )
  ) +
  ggplot2::geom_col(
    width = 0.60,
    color = "grey25",
    linewidth = 0.32
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label =
        scales::percent(
          r2,
          accuracy = 0.1
        )
    ),
    hjust = -0.12,
    size = 2.55
  ) +
  ggplot2::scale_fill_manual(
    values = view_colors,
    guide = "none"
  ) +
  ggplot2::scale_x_continuous(
    labels =
      scales::label_percent(
        accuracy = 1
      ),
    expand =
      ggplot2::expansion(
        mult = c(
          0,
          0.18
        )
      )
  ) +
  ggplot2::labs(
    title =
      mofa_v16_selected_factor,
    subtitle =
      "Variance explained by view",
    x = "R²",
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 8.4
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 10.5
      ),
    plot.subtitle =
      ggplot2::element_text(
        size = 7.6,
        color = "grey35"
      )
  )

mofa_v16_selected_baseline <-
  mofa_v16_baseline_long %>%
  dplyr::filter(
    as.character(
      factor
    ) ==
      mofa_v16_selected_factor
  )

mofa_v16_selected_paired <-
  mofa_v16_paired_long %>%
  dplyr::filter(
    as.character(
      factor
    ) ==
      mofa_v16_selected_factor
  )

mofa_v16_selected_baseline_p <-
  first_numeric_or_na(
    mofa_factor_wilcoxon_tests %>%
      dplyr::filter(
        factor ==
          mofa_v16_selected_factor,
        comparison ==
          "baseline_response"
      ) %>%
      dplyr::slice_head(
        n = 1
      ) %>%
      dplyr::pull(
        wilcoxon_p
      )
  )

mofa_v16_selected_paired_p <-
  first_numeric_or_na(
    mofa_factor_wilcoxon_tests %>%
      dplyr::filter(
        factor ==
          mofa_v16_selected_factor,
        comparison ==
          "paired_overall_change"
      ) %>%
      dplyr::slice_head(
        n = 1
      ) %>%
      dplyr::pull(
        wilcoxon_p
      )
  )

p_mofa_v16_selected_baseline <-
  ggplot2::ggplot(
    mofa_v16_selected_baseline,
    ggplot2::aes(
      x = TRG_plot,
      y = factor_score,
      fill = TRG_plot
    )
  ) +
  ggplot2::geom_violin(
    width = 0.42,
    trim = TRUE,
    alpha = 0.28,
    color = "grey35",
    linewidth = 0.32
  ) +
  ggplot2::geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    alpha = 0.40,
    linewidth = 0.34
  ) +
  ggbeeswarm::geom_quasirandom(
    width = 0.07,
    shape = 21,
    size = 1.38,
    color = "grey25",
    stroke = 0.30
  ) +
  ggplot2::scale_fill_manual(
    values =
      mofa_response_colors[
        c(
          "pCR",
          "non_pCR"
        )
      ],
    guide = "none"
  ) +
  ggplot2::labs(
    title =
      "Response at Baseline",
    subtitle = paste0(
      "Wilcoxon p = ",
      format_permutation_p(
        mofa_v16_selected_baseline_p
      )
    ),
    x = NULL,
    y = "Factor score"
  ) +
  ggplot2::theme_classic(
    base_size = 8.2
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 8.8
      ),
    plot.subtitle =
      ggplot2::element_text(
        size = 7.3,
        color = "grey35"
      )
  )

p_mofa_v16_selected_paired <-
  ggplot2::ggplot(
    mofa_v16_selected_paired,
    ggplot2::aes(
      x = Timepoint_display,
      y = factor_score,
      group = SubjectID
    )
  ) +
  ggplot2::geom_violin(
    ggplot2::aes(
      fill =
        Timepoint_display,
      group =
        Timepoint_display
    ),
    width = 0.42,
    trim = TRUE,
    alpha = 0.20,
    color = "grey40",
    linewidth = 0.32
  ) +
  ggplot2::geom_line(
    ggplot2::aes(
      color = TRG_plot
    ),
    linewidth = 0.34,
    alpha = 0.55
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      fill = TRG_plot
    ),
    shape = 21,
    size = 1.45,
    color = "grey25",
    stroke = 0.30
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      timepoint_display_colors,
      mofa_response_colors[
        c(
          "pCR",
          "non_pCR"
        )
      ]
    ),
    guide = "none"
  ) +
  ggplot2::scale_color_manual(
    values =
      mofa_response_colors[
        c(
          "pCR",
          "non_pCR"
        )
      ],
    guide = "none"
  ) +
  ggplot2::labs(
    title =
      "Paired treatment change",
    subtitle = paste0(
      "Wilcoxon p = ",
      format_permutation_p(
        mofa_v16_selected_paired_p
      )
    ),
    x = NULL,
    y = "Factor score"
  ) +
  ggplot2::theme_classic(
    base_size = 8.2
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 8.8
      ),
    plot.subtitle =
      ggplot2::element_text(
        size = 7.3,
        color = "grey35"
      ),
    axis.text.x =
      ggplot2::element_text(
        angle = 0,
        hjust = 0.5,
        size = 7.4
      )
  )

mofa_v16_selected_features <-
  mofa_feature_weights %>%
  dplyr::filter(
    factor ==
      mofa_v16_selected_factor,
    display_eligible
  ) %>%
  dplyr::group_by(
    view,
    direction
  ) %>%
  dplyr::slice_max(
    order_by = abs_weight,
    n = 3,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    feature_label_short =
      stringr::str_trunc(
        feature_label,
        width = 42
      ),
    feature_plot = factor(
      paste(
        view,
        feature_label_short,
        sep = "|||"
      ),
      levels = rev(
        unique(
          paste(
            view,
            feature_label_short,
            sep = "|||"
          )
        )
      )
    ),
    view_label = factor(
      as.character(
        view_label
      ),
      levels = unname(
        view_labels
      )
    )
  )

p_mofa_v16_selected_features <-
  ggplot2::ggplot(
    mofa_v16_selected_features,
    ggplot2::aes(
      x = weight_within_view,
      y = feature_plot,
      fill = view_label
    )
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    color = "grey55"
  ) +
  ggplot2::geom_col(
    width = 0.64,
    color = "grey25",
    linewidth = 0.24
  ) +
  ggplot2::facet_grid(
    view_label ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::scale_y_discrete(
    labels = function(x) {
      sub(
        "^[^|]+[|][|][|]",
        "",
        x
      )
    }
  ) +
  ggplot2::scale_fill_manual(
    values = view_colors,
    guide = "none"
  ) +
  ggplot2::labs(
    title =
      "Top within-view feature loadings",
    x =
      "Normalized MOFA loading",
    y = NULL
  ) +
  ggplot2::theme_classic(
    base_size = 8.0
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 9.0
      ),
    strip.background =
      ggplot2::element_blank(),
    strip.text.y =
      ggplot2::element_text(
        angle = 0,
        face = "bold",
        size = 7.0
      ),
    axis.text.y =
      ggplot2::element_text(
        size = 6.7
      ),
    panel.spacing.y =
      grid::unit(
        0.45,
        "lines"
      )
  )

p_mofa_v16_selected_factor_summary <-
  (
    p_mofa_v16_selected_factor_r2 /
    (
      p_mofa_v16_selected_baseline +
      p_mofa_v16_selected_paired
    )
  ) |
  p_mofa_v16_selected_features +
  patchwork::plot_layout(
    widths = c(
      0.95,
      1.55
    ),
    heights = c(
      0.70,
      1
    )
  ) +
  patchwork::plot_annotation(
    title = paste0(
      "Compact biological summary: ",
      mofa_v16_selected_factor
    ),
    subtitle = paste0(
      "Factor selection is exploratory. A dominant view does not invalidate smaller ",
      "cross-view contributions that pass the ",
      scales::percent(
        mofa_active_view_r2,
        accuracy = 1
      ),
      " threshold."
    )
  ) &
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold"
      )
  )

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename = paste0(
    "figures/mofa/publication_v16/",
    "MOFA_v16_selected_factor_summary_",
    mofa_v16_selected_factor,
    ".svg"
  ),
  plot =
    p_mofa_v16_selected_factor_summary,
  width = 8.0,
  height = 4.7,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}


#-----------------------------------------------------------------#
# 21.9 Clustered multi-factor feature-loading network
#-----------------------------------------------------------------#
# Feature nodes are clustered within each omics view according to their loading
# profiles across the selected factors. The resulting feature modules occupy the
# centre of the plot, while factor nodes surround them on the outside.

mofa_v16_response_factor_rank <-
  mofa_factor_wilcoxon_tests %>%
  dplyr::filter(
    comparison ==
      "baseline_response",
    factor %in%
      mofa_publication_factor_order
  ) %>%
  dplyr::arrange(
    wilcoxon_p,
    dplyr::desc(
      abs(
        wilcoxon_effect
      )
    )
  )

mofa_v16_pair_network_factors <- if (
  exists("mofa_v16_selected_pair_rows") &&
  nrow(mofa_v16_selected_pair_rows) > 0
) {
  unique(
    c(
      mofa_v16_selected_pair_rows$factor_x,
      mofa_v16_selected_pair_rows$factor_y
    )
  )
} else {
  character(0)
}

mofa_v16_shared_factor_rank <-
  mofa_factor_balance %>%
  dplyr::filter(
    factor %in%
      mofa_publication_factor_order
  ) %>%
  dplyr::arrange(
    dplyr::desc(shared_factor_priority),
    dplyr::desc(r2_sum)
  )

mofa_v16_network_factors <- unique(
  c(
    mofa_v16_pair_network_factors,
    utils::head(
      mofa_v16_response_factor_rank$factor,
      2
    ),
    if (
      mofa_v16_selected_factor %in%
        mofa_publication_factor_order
    ) {
      mofa_v16_selected_factor
    } else {
      character(0)
    },
    mofa_v16_shared_factor_rank$factor
  )
)

mofa_v16_network_factors <- utils::head(
  mofa_v16_network_factors[
    mofa_v16_network_factors %in%
      mofa_publication_factor_order
  ],
  mofa_network_factor_count
)

if (length(mofa_v16_network_factors) < 2) {
  mofa_v16_network_factors <- utils::head(
    mofa_publication_factor_order,
    min(
      mofa_network_factor_count,
      length(mofa_publication_factor_order)
    )
  )
}

mofa_v16_active_factor_views <-
  mofa_variance_explained %>%
  dplyr::mutate(
    factor = as.character(factor),
    view = as.character(view)
  ) %>%
  dplyr::filter(
    factor %in%
      mofa_v16_network_factors,
    r2 >= mofa_active_view_r2
  ) %>%
  dplyr::select(
    factor,
    view,
    view_r2 = r2
  )

mofa_v16_network_seed_edges <-
  mofa_feature_weights %>%
  dplyr::filter(
    factor %in%
      mofa_v16_network_factors,
    display_eligible
  ) %>%
  dplyr::inner_join(
    mofa_v16_active_factor_views,
    by = c(
      "factor",
      "view"
    )
  ) %>%
  dplyr::group_by(
    factor,
    view,
    direction
  ) %>%
  dplyr::slice_max(
    order_by = abs_weight,
    n = mofa_network_features_per_view_direction,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    feature_node_id = paste(
      view,
      feature,
      sep = "::"
    ),
    edge_id = paste(
      factor,
      feature_node_id,
      sep = "__"
    )
  )

mofa_v16_network_edges <-
  mofa_feature_weights %>%
  dplyr::filter(
    factor %in%
      mofa_v16_network_factors,
    display_eligible
  ) %>%
  dplyr::inner_join(
    mofa_v16_active_factor_views,
    by = c(
      "factor",
      "view"
    )
  ) %>%
  dplyr::mutate(
    feature_node_id = paste(
      view,
      feature,
      sep = "::"
    ),
    edge_id = paste(
      factor,
      feature_node_id,
      sep = "__"
    )
  ) %>%
  dplyr::filter(
    feature_node_id %in%
      unique(
        mofa_v16_network_seed_edges$feature_node_id
      )
  ) %>%
  dplyr::filter(
    abs(weight_within_view) >=
      mofa_network_shared_loading_threshold |
      edge_id %in%
        mofa_v16_network_seed_edges$edge_id
  ) %>%
  dplyr::mutate(
    loading_sign = ifelse(
      weight >= 0,
      "Positive",
      "Negative"
    ),
    loading_strength = abs(
      weight_within_view
    )
  )

cluster_network_features <- function(
    view_edges
) {
  feature_profile <- view_edges %>%
    dplyr::select(
      feature_node_id,
      factor,
      weight_within_view
    ) %>%
    tidyr::pivot_wider(
      names_from = factor,
      values_from = weight_within_view,
      values_fill = 0
    )

  for (
    factor_name in setdiff(
      mofa_v16_network_factors,
      colnames(feature_profile)
    )
  ) {
    feature_profile[[factor_name]] <- 0
  }

  profile_matrix <- as.matrix(
    feature_profile[
      ,
      mofa_v16_network_factors,
      drop = FALSE
    ]
  )
  rownames(profile_matrix) <-
    feature_profile$feature_node_id

  n_modules <- min(
    3,
    max(
      1,
      ceiling(
        nrow(profile_matrix) / 7
      )
    )
  )

  module <- if (
    nrow(profile_matrix) <= 2 ||
    n_modules == 1
  ) {
    rep(1L, nrow(profile_matrix))
  } else {
    stats::cutree(
      stats::hclust(
        stats::dist(profile_matrix),
        method = "ward.D2"
      ),
      k = n_modules
    )
  }

  data.frame(
    feature_node_id = rownames(profile_matrix),
    module = as.integer(module),
    stringsAsFactors = FALSE
  )
}

mofa_v16_network_feature_nodes <-
  mofa_v16_network_edges %>%
  dplyr::group_by(
    feature_node_id,
    view,
    view_label,
    feature,
    feature_label
  ) %>%
  dplyr::summarise(
    n_connected_factors =
      dplyr::n_distinct(factor),
    maximum_loading = max(
      loading_strength,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::group_by(view) %>%
  dplyr::group_modify(
    function(.x, .y) {
      module_table <- cluster_network_features(
        mofa_v16_network_edges %>%
          dplyr::filter(
            view == .y$view[[1]]
          )
      )

      dplyr::left_join(
        .x,
        module_table,
        by = "feature_node_id"
      )
    }
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    view_order = match(
      view,
      required_views
    ),
    view_center_x = c(
      species = -0.46,
      ko = -0.15,
      metabolite = 0.16,
      host = 0.47
    )[view],
    shared_feature =
      n_connected_factors >= 2,
    feature_label_short =
      stringr::str_trunc(
        feature_label,
        width = 30
      ),
    view_label = factor(
      as.character(view_label),
      levels = unname(view_labels)
    )
  ) %>%
  dplyr::group_by(
    view
  ) %>%
  dplyr::mutate(
    module_count = max(
      module,
      na.rm = TRUE
    ),
    module_center_y = if (
      module_count[1] == 1
    ) {
      0
    } else {
      seq(
        0.62,
        -0.62,
        length.out = module_count[1]
      )[module]
    }
  ) %>%
  dplyr::group_by(
    view,
    module
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_connected_factors),
    dplyr::desc(maximum_loading),
    .by_group = TRUE
  ) %>%
  dplyr::mutate(
    node_angle = if (
      dplyr::n() == 1
    ) {
      0
    } else {
      seq(
        0,
        2 * pi,
        length.out = dplyr::n() + 1
      )[seq_len(dplyr::n())]
    },
    x = view_center_x +
      0.065 * cos(node_angle),
    y = module_center_y +
      0.090 * sin(node_angle)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    label_x = x + ifelse(
      x < 0,
      -0.035,
      0.035
    ),
    label_hjust = ifelse(
      x < 0,
      1,
      0
    )
  )

mofa_v16_network_modules <-
  mofa_v16_network_feature_nodes %>%
  dplyr::group_by(
    view,
    view_label,
    module
  ) %>%
  dplyr::summarise(
    xmin = min(x) - 0.09,
    xmax = max(x) + 0.09,
    ymin = min(y) - 0.12,
    ymax = max(y) + 0.12,
    module_label = paste0(
      as.character(view_label[1]),
      " M",
      module[1]
    ),
    .groups = "drop"
  )

factor_angles <- seq(
  pi * 0.78,
  pi * 0.78 +
    2 * pi,
  length.out =
    length(mofa_v16_network_factors) + 1
)[seq_len(length(mofa_v16_network_factors))]

mofa_v16_network_factor_nodes <-
  data.frame(
    factor = mofa_v16_network_factors,
    x = 1.18 * cos(factor_angles),
    y = 0.88 * sin(factor_angles),
    stringsAsFactors = FALSE
  ) %>%
  dplyr::left_join(
    mofa_v16_response_factor_rank %>%
      dplyr::select(
        factor,
        wilcoxon_p
      ) %>%
      dplyr::distinct(
        factor,
        .keep_all = TRUE
      ),
    by = "factor"
  ) %>%
  dplyr::mutate(
    factor_label = paste0(
      factor,
      ifelse(
        is.finite(wilcoxon_p),
        paste0(
          "\nBaseline Wilcoxon p = ",
          format_permutation_p(wilcoxon_p)
        ),
        ""
      )
    )
  )

mofa_v16_network_plot_edges <-
  mofa_v16_network_edges %>%
  dplyr::left_join(
    mofa_v16_network_factor_nodes %>%
      dplyr::select(
        factor,
        x_factor = x,
        y_factor = y
      ),
    by = "factor"
  ) %>%
  dplyr::left_join(
    mofa_v16_network_feature_nodes %>%
      dplyr::select(
        feature_node_id,
        x_feature = x,
        y_feature = y
      ),
    by = "feature_node_id"
  )

p_mofa_v16_multifactor_feature_network <-
  ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = mofa_v16_network_modules,
    ggplot2::aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = view_label
    ),
    color = NA,
    alpha = 0.10
  ) +
  ggplot2::geom_curve(
    data = mofa_v16_network_plot_edges %>%
      dplyr::filter(
        x_factor < 0
      ),
    ggplot2::aes(
      x = x_factor,
      y = y_factor,
      xend = x_feature,
      yend = y_feature,
      color = loading_sign,
      linewidth = loading_strength
    ),
    curvature = 0.12,
    alpha = 0.42,
    lineend = "round"
  ) +
  ggplot2::geom_curve(
    data = mofa_v16_network_plot_edges %>%
      dplyr::filter(
        x_factor >= 0
      ),
    ggplot2::aes(
      x = x_factor,
      y = y_factor,
      xend = x_feature,
      yend = y_feature,
      color = loading_sign,
      linewidth = loading_strength
    ),
    curvature = -0.12,
    alpha = 0.42,
    lineend = "round"
  ) +
  ggplot2::geom_label(
    data = mofa_v16_network_factor_nodes,
    ggplot2::aes(
      x = x,
      y = y,
      label = factor_label
    ),
    fill = "#FFF2B3",
    color = "grey10",
    label.size = 0.32,
    label.padding = grid::unit(
      0.12,
      "lines"
    ),
    size = 2.45,
    fontface = "bold"
  ) +
  ggplot2::geom_point(
    data = mofa_v16_network_feature_nodes,
    ggplot2::aes(
      x = x,
      y = y,
      shape = view_label,
      fill = view_label,
      size = maximum_loading
    ),
    color = "grey20",
    stroke = 0.48
  ) +
  ggplot2::geom_point(
    data = mofa_v16_network_feature_nodes %>%
      dplyr::filter(
        shared_feature
      ),
    ggplot2::aes(
      x = x,
      y = y
    ),
    shape = 21,
    size = 4.7,
    fill = NA,
    color = "#7A3E9D",
    stroke = 0.85
  ) +
  ggplot2::geom_text(
    data = mofa_v16_network_feature_nodes,
    ggplot2::aes(
      x = label_x,
      y = y,
      label = feature_label_short,
      hjust = label_hjust
    ),
    size = 2.05
  ) +
  ggplot2::geom_text(
    data = mofa_v16_network_modules,
    ggplot2::aes(
      x = (xmin + xmax) / 2,
      y = ymax + 0.035,
      label = module_label,
      color = view_label
    ),
    fontface = "bold",
    size = 2.15,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(
    values = c(
      Positive = "#D55E00",
      Negative = "#0072B2",
      view_colors
    ),
    name = "Loading sign"
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      Species = 21,
      `KEGG ortholog` = 22,
      Metabolite = 23,
      `Host RNA-seq` = 24
    ),
    name = "Omics view"
  ) +
  ggplot2::scale_fill_manual(
    values = view_colors,
    name = "Omics view"
  ) +
  ggplot2::scale_size_continuous(
    range = c(
      2.2,
      4.0
    ),
    guide = "none"
  ) +
  ggplot2::scale_linewidth_continuous(
    range = c(
      0.28,
      1.05
    ),
    guide = "none"
  ) +
  ggplot2::coord_cartesian(
    xlim = c(
      -1.45,
      1.45
    ),
    ylim = c(
      -1.08,
      1.08
    ),
    clip = "off"
  ) +
  ggplot2::labs(
    title = "Clustered multi-factor feature-loading network",
    subtitle = paste0(
      "Central modules group features with similar loading profiles across the selected factors; ",
      "edges are MOFA loadings and do not represent direct feature-feature interactions."
    ),
    caption = paste0(
      "Factor-label 'Baseline Wilcoxon p' is the two-sided Wilcoxon rank-sum p-value ",
      "for pCR versus non-pCR factor scores at Baseline. Purple rings mark features connected to multiple factors."
    )
  ) +
  ggplot2::theme_void(
    base_size = 8.8
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 10.0
    ),
    plot.subtitle = ggplot2::element_text(
      size = 7.5,
      color = "grey35"
    ),
    plot.caption = ggplot2::element_text(
      size = 6.8,
      color = "grey35",
      hjust = 0
    ),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(
      8,
      34,
      8,
      34
    )
  )

if (isTRUE(mofa_save_legacy_figures)) {
ggplot2::ggsave(
  filename = paste0(
    "figures/mofa/publication_v16/",
    "MOFA_v16_multifactor_feature_loading_network.svg"
  ),
  plot = p_mofa_v16_multifactor_feature_network,
  width = 11.2,
  height = 7.2,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
}


#-----------------------------------------------------------------#
# 21.10 Factor interpretation audit
#-----------------------------------------------------------------#

mofa_v16_factor_interpretation <-
  mofa_factor_balance %>%
  dplyr::mutate(
    active_views_primary =
      rowSums(
        dplyr::across(
          dplyr::all_of(
            required_views
          ),
          ~ .x >=
            mofa_active_view_r2
        ),
        na.rm = TRUE
      ),
    dominance_class =
      dplyr::case_when(
        strongest_view_fraction >= 0.80 &
          active_views_primary >= 2 ~
          "Dominant view with smaller cross-view contributions",
        strongest_view_fraction >= 0.80 ~
          "Predominantly view-specific",
        active_views_primary >= 3 ~
          "Broadly shared",
        active_views_primary == 2 ~
          "Pair-shared",
        TRUE ~
          "Low cross-view activity"
      ),
    interpretation_note =
      dplyr::case_when(
        strongest_view_fraction >= 0.80 &
          active_views_primary >= 2 ~
          paste0(
            "Do not infer causality from unequal R2. The smaller views co-vary with ",
            "the dominant program but directionality requires temporal or mechanistic evidence."
          ),
        strongest_view_fraction >= 0.80 ~
          paste0(
            "Retain as a valid view-dominant latent factor unless a separate QC or stability ",
            "diagnostic indicates artefact."
          ),
        TRUE ~
          "Interpret as a shared covariance pattern; causal order is not identified by MOFA."
      )
  ) %>%
  dplyr::left_join(
    mofa_factor_subject_specificity %>%
      dplyr::select(
        factor,
        icc,
        paired_spearman,
        nearest_same_subject_rate
      ),
    by = "factor"
  )

mofa_v16_publication_object_names <- intersect(
  c(
    "mofa_feature_count_audit",
    "mofa_v16_factor_interpretation",
    "mofa_factor_wilcoxon_tests",
    "mofa_v16_association_data",
    "mofa_v16_pair_permanova_atlas",
    "mofa_v16_subject_mean_factor_data",
    "mofa_publication_factor_audit",
    "mofa_v16_selected_pair_rows",
    "mofa_v16_network_edges",
    "mofa_v16_network_feature_nodes",
    "mofa_v16_selected_factor",
    "p_mofa_v16_feature_selection",
    "p_mofa_v16_availability_horizontal",
    "p_mofa_v16_availability_vertical",
    "p_mofa_v16_variance_overview",
    "p_mofa_v16_factor_associations",
    "p_mofa_v16_overall_response",
    "p_mofa_v16_baseline_response",
    "p_mofa_v16_after_rt_response",
    "p_mofa_v16_paired_change",
    "p_mofa_v16_subject_repeatability",
    "p_mofa_v16_selected_factor_summary",
    "p_mofa_v16_selected_response_pairs",
    "p_mofa_v16_multifactor_feature_network"
  ),
  ls()
)

save(
  list =
    mofa_v16_publication_object_names,
  file = paste0(
    "results/mofa/publication_v16/",
    "MOFA_v16_publication_objects.RData"
  )
)



#-----------------------------------------------------------------#
# Save completed core analysis before publication rendering.
# This file is an output/recovery aid only; the script never loads it.
#-----------------------------------------------------------------#

save(
  mofa_sample_metadata,
  mofa_availability_wide,
  mofa_feature_selection_summary,
  mofa_selected_feature_variance,
  mofa_model_metadata,
  mofa_factor_scores,
  mofa_factor_total_level_correlations,
  mofa_variance_explained,
  mofa_variance_total,
  mofa_factor_summary,
  mofa_factor_balance,
  mofa_factor_wilcoxon_tests,
  mofa_factor_paired_data,
  mofa_feature_weights,
  mofa_v30_pair_permanova_atlas,
  file = "results/mofa/MOFA_4omics_analysis_v30_core.RData"
)


#=================================================================#
# 22. v50 publication analysis and response-integrated visualization
#=================================================================#

mofa_v50_figure_dir <- "figures/mofa/publication_v50"
mofa_v50_result_dir <- "results/mofa/publication_v50"
dir.create(mofa_v50_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mofa_v50_result_dir, recursive = TRUE, showWarnings = FALSE)

# The pairwise atlas is calculated upstream from the validated latent scores.
# Alias it once so every publication-layer object and output follows v50 naming.
mofa_v50_pair_permanova_atlas <- get(
  "mofa_v30_pair_permanova_atlas",
  envir = environment(),
  inherits = FALSE
)

# Screening settings used only for the v50 interpretation layer.
# Feature display reuses mofa_active_view_r2, network size reuses
# mofa_network_factor_count, and final tests reuse the 999-permutation setting.
mofa_v50_triple_screen_permutations <- 999L
mofa_v50_outlier_sensitivity_permutations <- 999L

required_views <- c("species", "ko", "metabolite", "host")
view_labels <- c(
  species = "Species",
  ko = "KEGG ortholog",
  metabolite = "Metabolite",
  host = "Host RNA-seq"
)
view_colors <- c(
  Species = "#79B3A3",
  `KEGG ortholog` = "#D8A46F",
  Metabolite = "#82A8C7",
  `Host RNA-seq` = "#B29AC6"
)
mofa_group_colors <- c(
  pCR = "#4FAE9A",
  non_pCR = "#DE7872",
  CR = "#4FAE9A",
  nonCR = "#DE7872"
)
timepoint_display_colors <- c(
  Baseline = "#667A8A",
  `After RT` = "#C59A5B"
)

mofa_v50_source_revision <- "v33_response_integrated with robust 3D annotation and point rendering"

mofa_v50_required_objects <- c(
  "mofa_factor_scores",
  "mofa_variance_explained",
  "mofa_factor_wilcoxon_tests",
  "mofa_factor_paired_data",
  "mofa_feature_weights",
  "mofa_factor_balance",
  "mofa_v50_pair_permanova_atlas",
  "mofa_factor_total_level_correlations"
)
# Do not use vapply(required_names, exists, inherits = FALSE) here.
# In that pattern, exists() is evaluated from vapply's function frame and can
# falsely report every object as missing even when the objects exist in the
# script environment. Compare names directly with the current environment.
mofa_v50_execution_environment <- environment()
mofa_v50_missing_objects <- setdiff(
  mofa_v50_required_objects,
  ls(
    envir = mofa_v50_execution_environment,
    all.names = TRUE
  )
)
if (length(mofa_v50_missing_objects) > 0) {
  stop(
    paste0(
      "The upstream v30 analysis stage did not create required objects: ",
      paste(mofa_v50_missing_objects, collapse = ", "),
      ". This script does not load any prior analysis RData; inspect the first ",
      "earlier error in the same source() run."
    ),
    call. = FALSE
  )
}

mofa_v50_sort_factors <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & grepl("^Factor[0-9]+$", x)]
  x[order(as.integer(sub("^Factor", "", x)))]
}

mofa_v50_format_p <- function(p) {
  ifelse(
    !is.finite(p),
    "NA",
    ifelse(
      p < 0.001,
      "<0.001",
      formatC(p, format = "f", digits = 3)
    )
  )
}

# Convert p values without truncating values above 1 on the -log10 scale.
# The upper legend limit is fixed at 1 only when every observed value is below 1;
# otherwise it follows the observed maximum.
mofa_v50_neglog10_p <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  positive_p <- p[is.finite(p) & p > 0]
  zero_floor <- if (length(positive_p) > 0) {
    min(positive_p, na.rm = TRUE) / 10
  } else {
    1e-06
  }
  p_safe <- ifelse(is.finite(p) & p > 0, p, ifelse(is.finite(p) & p == 0, zero_floor, NA_real_))
  -log10(p_safe)
}

mofa_v50_neglog10_upper <- function(x, minimum = 1) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x) & x >= 0]
  if (length(x) == 0) {
    return(minimum)
  }
  max(minimum, max(x, na.rm = TRUE))
}

mofa_v50_neglog10_breaks <- function(upper) {
  upper <- max(1, as.numeric(upper))
  breaks <- pretty(c(0, upper), n = 4)
  unique(c(0, breaks[breaks > 0 & breaks < upper], upper))
}

# Reused by the feature-weight bar plot and the multi-factor network.
# Species bin identifiers (for example, GGB9634_SGB15093) remain upright;
# conventional binomial species names are italicized.
mofa_v50_plotmath_feature_label <- function(view, label) {
  view <- as.character(view)
  label <- as.character(label)

  vapply(
    seq_along(label),
    function(index) {
      is_species_name <-
        view[index] == "species" &&
        !grepl("^(GGB|SGB)[0-9]+", label[index])

      display_label <- if (is_species_name) {
        gsub("_", " ", label[index], fixed = TRUE)
      } else {
        label[index]
      }

      quoted_label <- encodeString(display_label, quote = '"')

      if (is_species_name) {
        paste0("italic(", quoted_label, ")")
      } else {
        quoted_label
      }
    },
    character(1)
  )
}

mofa_v50_factor_order <- mofa_v50_sort_factors(
  grep("^Factor[0-9]+$", colnames(mofa_factor_scores), value = TRUE)
)
if (length(mofa_v50_factor_order) == 0) {
  stop("No Factor1, Factor2, ... columns were found in mofa_factor_scores.", call. = FALSE)
}


mofa_v50_factor_view_display_audit <- mofa_variance_explained %>%
  dplyr::mutate(
    factor = as.character(factor),
    view = as.character(view),
    view_label = as.character(view_label),
    feature_display_active = is.finite(r2) & r2 >= mofa_active_view_r2
  ) %>%
  dplyr::arrange(
    match(factor, mofa_v50_factor_order),
    match(view, required_views)
  )

write.csv(
  mofa_v50_factor_view_display_audit,
  file.path(mofa_v50_result_dir, "MOFA_v50_factor_view_feature_display_audit.csv"),
  row.names = FALSE
)

# MOFA does not provide a penalty that forces every factor to be equally shared
# across views. The v16-compatible model uses scale_views and weight ARD; ARD is
# allowed to retain view-specific factors. v30 therefore keeps the validated model
# and adds explicit global-level QC plus a sharedness-penalized interpretation rank.
# Separate pCR and non-pCR MOFA fits are not used as the primary test because factor
# signs, orders, and rotations are not directly identifiable across independent fits.
# Response-specific biology is tested on subject-level After RT minus Baseline deltas.
mofa_v50_view_balance_strategy <- data.frame(
  component = c(
    "View scaling",
    "Weight ARD",
    "Spike-slab weights",
    "Global-level QC",
    "Sharedness-penalized ranking",
    "Response-specific inference"
  ),
  setting = c(
    "scale_views = TRUE in the v16-compatible fit",
    "ard_weights = TRUE; view-specific factors remain permissible",
    "FALSE in the reused v16 model",
    "Flag |cor(sample total, factor)| >= 0.75",
    "R2 sum penalized by strongest-view fraction, effective views, and QC flag",
    "Wilcoxon rank-sum on subject-level paired deltas"
  ),
  interpretation = c(
    "Equalizes overall view variance but does not guarantee equal factor sharing",
    "Regularizes factor-view loadings and can preserve biologically valid view specificity",
    "No additional point-mass sparsity prior was imposed in the validated model",
    "Separates global abundance or detection gradients from specific programs",
    "Prevents one-view/global-level factors from automatically dominating publication selection",
    "Directly tests whether treatment-associated factor direction differs by response group"
  ),
  stringsAsFactors = FALSE
)
#-----------------------------------------------------------------#
# 22.1 View-dominance and global-level QC
#-----------------------------------------------------------------#

mofa_v50_global_level_diagnostics <- mofa_factor_total_level_correlations %>%
  dplyr::transmute(
    view = as.character(view),
    view_label = as.character(view_label),
    factor = as.character(factor),
    metric = "sample_sum",
    correlation = as.numeric(correlation),
    qc_flag = is.finite(correlation) & abs(correlation) >= 0.75
  )

# Add mean, detected-feature, and zero-fraction diagnostics from the exact
# matrices supplied to the primary MOFA model.
if (exists("mofa_data", inherits = FALSE)) {
  mofa_v50_extra_global_diagnostics <- dplyr::bind_rows(
    lapply(required_views, function(view_name) {
      view_matrix <- as.matrix(mofa_data[[view_name]])
      common_samples <- intersect(colnames(view_matrix), mofa_factor_scores$SampleID)
      if (length(common_samples) < 4) {
        return(NULL)
      }
      factor_matrix <- as.matrix(
        mofa_factor_scores[
          match(common_samples, mofa_factor_scores$SampleID),
          mofa_v50_factor_order,
          drop = FALSE
        ]
      )
      rownames(factor_matrix) <- common_samples
      metric_values <- list(
        sample_mean = colMeans(view_matrix[, common_samples, drop = FALSE], na.rm = TRUE),
        detected_features = colSums(
          is.finite(view_matrix[, common_samples, drop = FALSE]) &
            view_matrix[, common_samples, drop = FALSE] != 0
        ),
        zero_fraction = colMeans(
          is.finite(view_matrix[, common_samples, drop = FALSE]) &
            view_matrix[, common_samples, drop = FALSE] == 0
        )
      )
      dplyr::bind_rows(lapply(names(metric_values), function(metric_name) {
        values <- metric_values[[metric_name]]
        data.frame(
          view = view_name,
          view_label = unname(view_labels[view_name]),
          factor = mofa_v50_factor_order,
          metric = metric_name,
          correlation = vapply(mofa_v50_factor_order, function(factor_name) {
            keep <- is.finite(values) & is.finite(factor_matrix[, factor_name])
            if (sum(keep) < 4) {
              return(NA_real_)
            }
            stats::cor(values[keep], factor_matrix[keep, factor_name])
          }, numeric(1)),
          stringsAsFactors = FALSE
        )
      }))
    })
  ) %>%
    dplyr::mutate(qc_flag = is.finite(correlation) & abs(correlation) >= 0.75)

  mofa_v50_global_level_diagnostics <- dplyr::bind_rows(
    mofa_v50_global_level_diagnostics,
    mofa_v50_extra_global_diagnostics
  )
}

mofa_v50_global_factor_flag <- mofa_v50_global_level_diagnostics %>%
  dplyr::group_by(factor) %>%
  dplyr::summarise(
    global_level_flag = any(qc_flag, na.rm = TRUE),
    maximum_global_abs_r = if (any(is.finite(correlation))) {
      max(abs(correlation), na.rm = TRUE)
    } else {
      NA_real_
    },
    .groups = "drop"
  )

mofa_v50_view_dominance_audit <- mofa_factor_balance %>%
  dplyr::mutate(factor = as.character(factor)) %>%
  dplyr::left_join(mofa_v50_global_factor_flag, by = "factor") %>%
  dplyr::mutate(
    global_level_flag = dplyr::coalesce(global_level_flag, FALSE),
    strongest_view_fraction = as.numeric(strongest_view_fraction),
    effective_views = as.numeric(effective_views),
    sharedness_penalized_r2 =
      r2_sum *
      pmax(0, 1 - strongest_view_fraction) *
      pmin(effective_views, length(required_views)) / length(required_views) *
      ifelse(global_level_flag, 0.50, 1.00),
    factor = factor(
      factor,
      levels = rev(mofa_v50_factor_order)
    )
  )

mofa_v50_qc_plot_data <- mofa_v50_global_level_diagnostics %>%
  dplyr::filter(metric == "sample_sum") %>%
  dplyr::mutate(
    factor = factor(factor, levels = rev(mofa_v50_factor_order)),
    view_label = factor(view_label, levels = unname(view_labels)),
    abs_correlation = abs(correlation)
  )

p_mofa_v50_global_level_qc <- ggplot2::ggplot(
  mofa_v50_qc_plot_data,
  ggplot2::aes(x = view_label, y = factor, fill = abs_correlation)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.30) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(is.finite(correlation), sprintf("%+.2f", correlation), "")),
    size = 2.4
  ) +
  ggplot2::geom_point(
    data = mofa_v50_qc_plot_data %>% dplyr::filter(qc_flag),
    shape = 21,
    size = 3.0,
    fill = NA,
    color = "black",
    stroke = 0.70
  ) +
  ggplot2::scale_fill_gradient(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "white",
    high = "#2166AC",
    limits = c(0, 1),
    oob = scales::squish,
    name = "|Correlation|"
  ) +
  ggplot2::labs(
    title = "MOFA factor correlation with sample-level total signal",
    subtitle = "Open circles mark |r| >= 0.75; these factors require technical-versus-biological interpretation.",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 9.2) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.7, color = "grey35"),
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, size = 7.8),
    axis.text.y = ggplot2::element_text(face = "plain", size = 8.0)
  )

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_global_level_QC.svg"),
  plot = p_mofa_v50_global_level_qc,
  width = 6.0,
  height = max(4.4, 1.8 + 0.28 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

p_mofa_v50_view_dominance <- ggplot2::ggplot(
  mofa_v50_view_dominance_audit,
  ggplot2::aes(
    x = strongest_view_fraction,
    y = factor,
    size = r2_sum,
    fill = effective_views,
    shape = global_level_flag
  )
) +
  ggplot2::geom_vline(xintercept = 0.80, linetype = "dashed", color = "grey55", linewidth = 0.45) +
  ggplot2::geom_point(color = "grey20", stroke = 0.55) +
  ggplot2::scale_fill_gradient(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "white",
    high = "#2166AC",
    limits = c(1, 4),
    oob = scales::squish,
    name = "Effective\nviews"
  ) +
  ggplot2::scale_shape_manual(
    values = c(`FALSE` = 21, `TRUE` = 23),
    labels = c(`FALSE` = "No global-level flag", `TRUE` = "Global-level flag"),
    name = NULL
  ) +
  ggplot2::scale_x_continuous(limits = c(0.25, 1), breaks = c(0.25, 0.50, 0.75, 1.00)) +
  ggplot2::scale_size_continuous(range = c(2.4, 5.2), guide = "none") +
  ggplot2::labs(
    title = "Factor sharedness and view dominance",
    subtitle = "Values near 1 indicate that one omics view explains most factor-associated variance.",
    x = "Strongest-view fraction of factor R²",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 9.2) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.7, color = "grey35"),
    axis.text.y = ggplot2::element_text(face = "plain", size = 7.4),
    legend.position = "right"
  )

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_view_dominance.svg"),
  plot = p_mofa_v50_view_dominance,
  width = 6.3,
  height = max(4.5, 1.8 + 0.27 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

#-----------------------------------------------------------------#
# 22.2 Numeric-order variance and p-value association maps
#-----------------------------------------------------------------#

mofa_v50_variance_plot_data <- mofa_variance_explained %>%
  dplyr::mutate(
    factor = factor(as.character(factor), levels = rev(mofa_v50_factor_order)),
    view_label = factor(as.character(view_label), levels = unname(view_labels))
  )

p_mofa_v50_variance_heatmap <- ggplot2::ggplot(
  mofa_v50_variance_plot_data,
  ggplot2::aes(x = view_label, y = factor, fill = r2)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.30) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(r2 >= 0.005, sprintf("%.1f", 100 * r2), "")),
    size = 2.65
  ) +
  ggplot2::scale_fill_gradient(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "white",
    high = "#2166AC",
    name = "Variance\nexplained"
  ) +
  ggplot2::labs(
    title = "Variance explained by factor and omics view",
    subtitle = "Cell labels are percentages; factors are ordered numerically.",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 9.2) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.7, color = "grey35"),
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, size = 7.8),
    axis.text.y = ggplot2::element_text(face = "plain", size = 8.0)
  )

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_variance_explained_numeric_order.svg"),
  plot = p_mofa_v50_variance_heatmap,
  width = 5.1,
  height = max(4.4, 1.8 + 0.28 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

mofa_v50_association_data <- mofa_factor_wilcoxon_tests %>%
  dplyr::filter(
    factor %in% mofa_v50_factor_order,
    comparison %in% c(
      "overall_response_subject_mean",
      "baseline_response",
      "ongoing_response",
      "paired_overall_change",
      "paired_differential_change"
    )
  ) %>%
  dplyr::mutate(
    comparison_label = dplyr::recode(
      comparison,
      overall_response_subject_mean = "Overall response\nsubject mean",
      baseline_response = "Response at\nBaseline",
      ongoing_response = "Response after\nRT",
      paired_overall_change = "Paired\nchange",
      paired_differential_change = "Differential\npaired change"
    ),
    comparison_label = factor(
      comparison_label,
      levels = c(
        "Overall response\nsubject mean",
        "Response at\nBaseline",
        "Response after\nRT",
        "Paired\nchange",
        "Differential\npaired change"
      )
    ),
    factor = factor(factor, levels = rev(mofa_v50_factor_order)),
    minus_log10_p = mofa_v50_neglog10_p(wilcoxon_p),
    effect_label = ifelse(is.finite(wilcoxon_effect), sprintf("%+.2f", wilcoxon_effect), "")
  )

mofa_v50_association_p_upper <- mofa_v50_neglog10_upper(
  mofa_v50_association_data$minus_log10_p
)
mofa_v50_association_p_breaks <- mofa_v50_neglog10_breaks(
  mofa_v50_association_p_upper
)

p_mofa_v50_association_pmap <- ggplot2::ggplot(
  mofa_v50_association_data,
  ggplot2::aes(x = comparison_label, y = factor, fill = minus_log10_p)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.30) +
  ggplot2::geom_text(ggplot2::aes(label = effect_label), size = 2.15) +
  ggplot2::scale_fill_gradient(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "white",
    high = "#2166AC",
    limits = c(0, mofa_v50_association_p_upper),
    breaks = mofa_v50_association_p_breaks,
    oob = scales::squish,
    name = expression(-log[10](italic(p)))
  ) +
  ggplot2::labs(
    title = "Non-parametric factor associations",
    subtitle = paste0(
      "Text is the signed median contrast; the -log10(p) scale has a minimum upper bound of 1 ",
      "and expands to ", formatC(mofa_v50_association_p_upper, format = "f", digits = 1),
      " when stronger evidence is present."
    ),
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 9.2) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.7, color = "grey35"),
    axis.text.x = ggplot2::element_text(size = 9.0, face = "plain"),
    axis.text.y = ggplot2::element_text(face = "plain", size = 7.4)
  )

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_association_pvalue_map.svg"),
  plot = p_mofa_v50_association_pmap,
  width = 6.3,
  height = max(4.5, 1.9 + 0.28 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

#-----------------------------------------------------------------#
# 22.3 Response-specific factor trajectories
#-----------------------------------------------------------------#

mofa_v50_factor_long <- mofa_factor_scores %>%
  dplyr::select(
    SampleID,
    SubjectID,
    Timepoint,
    TRG_plot,
    dplyr::all_of(mofa_v50_factor_order)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(mofa_v50_factor_order),
    names_to = "factor",
    values_to = "factor_score"
  ) %>%
  dplyr::mutate(
    factor = factor(factor, levels = mofa_v50_factor_order),
    Timepoint = as.character(Timepoint),
    TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR"))
  ) %>%
  dplyr::filter(
    is.finite(factor_score),
    !is.na(SubjectID),
    !is.na(TRG_plot)
  )

mofa_v50_overall_response_long <- mofa_v50_factor_long %>%
  dplyr::group_by(SubjectID, TRG_plot, factor) %>%
  dplyr::summarise(
    factor_score = mean(factor_score, na.rm = TRUE),
    .groups = "drop"
  )

mofa_v50_baseline_response_long <- mofa_v50_factor_long %>%
  dplyr::filter(Timepoint == "Before") %>%
  dplyr::group_by(SubjectID, TRG_plot, factor) %>%
  dplyr::summarise(factor_score = dplyr::first(factor_score), .groups = "drop")

mofa_v50_after_response_long <- mofa_v50_factor_long %>%
  dplyr::filter(Timepoint == "Ongoing") %>%
  dplyr::group_by(SubjectID, TRG_plot, factor) %>%
  dplyr::summarise(factor_score = dplyr::first(factor_score), .groups = "drop")

# Standardized pCR-versus-non-pCR enrichment by factor. Positive values indicate
# higher factor scores in pCR; negative values indicate higher scores in non-pCR.
mofa_v50_response_enrichment_data <- dplyr::bind_rows(
  mofa_v50_baseline_response_long %>% dplyr::mutate(Timepoint_display = "Baseline"),
  mofa_v50_after_response_long %>% dplyr::mutate(Timepoint_display = "After RT")
) %>%
  dplyr::group_by(factor, Timepoint_display) %>%
  dplyr::summarise(
    standardized_effect = calculate_hedges_g(
      factor_score[TRG_plot == "pCR"],
      factor_score[TRG_plot == "non_pCR"]
    ),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    mofa_factor_wilcoxon_tests %>%
      dplyr::filter(comparison %in% c("baseline_response", "ongoing_response")) %>%
      dplyr::transmute(
        factor = factor(as.character(factor), levels = mofa_v50_factor_order),
        Timepoint_display = ifelse(comparison == "baseline_response", "Baseline", "After RT"),
        wilcoxon_p,
        wilcoxon_fdr
      ),
    by = c("factor", "Timepoint_display")
  ) %>%
  dplyr::mutate(
    factor = factor(as.character(factor), levels = rev(mofa_v50_factor_order)),
    Timepoint_display = factor(Timepoint_display, levels = c("Baseline", "After RT")),
    minus_log10_p = mofa_v50_neglog10_p(wilcoxon_p),
    response_direction = factor(
      ifelse(standardized_effect >= 0, "pCR", "non_pCR"),
      levels = c("pCR", "non_pCR")
    ),
    significance_label = dplyr::case_when(
      is.finite(wilcoxon_p) & wilcoxon_p < 0.001 ~ "***",
      is.finite(wilcoxon_p) & wilcoxon_p < 0.01 ~ "**",
      is.finite(wilcoxon_p) & wilcoxon_p < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

mofa_v50_response_enrichment_p_upper <- mofa_v50_neglog10_upper(
  mofa_v50_response_enrichment_data$minus_log10_p
)

p_mofa_v50_response_enrichment <- ggplot2::ggplot(
  mofa_v50_response_enrichment_data,
  ggplot2::aes(
    x = standardized_effect,
    y = factor,
    color = response_direction,
    fill = response_direction,
    size = minus_log10_p
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.42
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend = standardized_effect,
      yend = factor,
      color = response_direction
    ),
    linewidth = 0.62,
    alpha = 0.72,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(
    shape = 21,
    color = "grey20",
    stroke = 0.55
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = significance_label),
    nudge_y = 0.25,
    size = 3.0,
    color = "black",
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~ Timepoint_display, nrow = 1) +
  ggplot2::scale_color_manual(
    values = mofa_group_colors[c("pCR", "non_pCR")],
    labels = c(pCR = "pCR-enriched", non_pCR = "non-pCR-enriched"),
    name = "Direction"
  ) +
  ggplot2::scale_fill_manual(
    values = mofa_group_colors[c("pCR", "non_pCR")],
    guide = "none"
  ) +
  ggplot2::scale_size_continuous(
    range = c(2.5, 5.2),
    limits = c(0, mofa_v50_response_enrichment_p_upper),
    breaks = mofa_v50_neglog10_breaks(mofa_v50_response_enrichment_p_upper),
    name = expression(-log[10](italic(p)))
  ) +
  ggplot2::labs(
    title = "pCR versus non-pCR factor enrichment",
    subtitle = paste0(
      "Hedges' g is oriented as pCR minus non-pCR. Direction uses the prespecified response colors; ",
      "point size shows the two-sided Wilcoxon rank-sum p value."
    ),
    x = "Standardized response enrichment (Hedges' g)",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 9.2) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.7, color = "grey35"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "plain"),
    axis.text.y = ggplot2::element_text(face = "plain", size = 7.6),
    legend.position = "right"
  )

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_enrichment.svg"),
  plot = p_mofa_v50_response_enrichment,
  width = 7.0,
  height = max(4.6, 1.9 + 0.28 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

mofa_v50_paired_factor_data <- dplyr::inner_join(
  mofa_v50_baseline_response_long %>%
    dplyr::transmute(SubjectID, TRG_plot, factor, Baseline = factor_score),
  mofa_v50_after_response_long %>%
    dplyr::transmute(SubjectID, factor, `After RT` = factor_score),
  by = c("SubjectID", "factor")
) %>%
  dplyr::mutate(
    factor_delta = `After RT` - Baseline,
    factor = factor(as.character(factor), levels = mofa_v50_factor_order)
  )

mofa_v50_paired_long <- mofa_v50_paired_factor_data %>%
  tidyr::pivot_longer(
    cols = c("Baseline", "After RT"),
    names_to = "Timepoint_display",
    values_to = "factor_score"
  ) %>%
  dplyr::mutate(
    Timepoint_display = factor(Timepoint_display, levels = c("Baseline", "After RT"))
  )

mofa_v50_make_annotation <- function(plot_data, comparison_name) {
  plot_data %>%
    dplyr::group_by(factor) %>%
    dplyr::summarise(
      y_min = min(factor_score, na.rm = TRUE),
      y_max = max(factor_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      y_position = y_max + 0.12 * pmax(y_max - y_min, 1)
    ) %>%
    dplyr::left_join(
      mofa_factor_wilcoxon_tests %>%
        dplyr::filter(comparison == comparison_name) %>%
        dplyr::transmute(
          factor = factor(as.character(factor), levels = mofa_v50_factor_order),
          wilcoxon_effect,
          wilcoxon_p
        ),
      by = "factor"
    ) %>%
    dplyr::mutate(
      label = paste0(
        "pCR - non-pCR = ", sprintf("%+.2f", wilcoxon_effect),
        "\nWilcoxon ",
        ifelse(
          !is.finite(wilcoxon_p),
          "p = NA",
          ifelse(wilcoxon_p < 0.001, "p < 0.001", paste0("p = ", mofa_v50_format_p(wilcoxon_p)))
        )
      )
    )
}

mofa_v50_build_response_plot <- function(plot_data, annotation_data, title_text) {
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = TRG_plot, y = factor_score, fill = TRG_plot)
  ) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.68, color = "grey30", linewidth = 0.38) +
    ggplot2::geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      fill = NA,
      coef = 0,
      staplewidth = 0,
      linewidth = 0.42
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.08, height = 0),
      shape = 21,
      size = 1.45,
      color = "grey20",
      stroke = 0.32,
      alpha = 0.85
    ) +
    ggplot2::geom_segment(
      data = annotation_data,
      ggplot2::aes(x = 1, xend = 2, y = y_position, yend = y_position),
      inherit.aes = FALSE,
      linewidth = 0.36,
      color = "grey25"
    ) +
    ggplot2::geom_segment(
      data = annotation_data,
      ggplot2::aes(
        x = 1,
        xend = 1,
        y = y_position,
        yend = y_position - 0.035 * pmax(y_max - y_min, 1)
      ),
      inherit.aes = FALSE,
      linewidth = 0.36,
      color = "grey25"
    ) +
    ggplot2::geom_segment(
      data = annotation_data,
      ggplot2::aes(
        x = 2,
        xend = 2,
        y = y_position,
        yend = y_position - 0.035 * pmax(y_max - y_min, 1)
      ),
      inherit.aes = FALSE,
      linewidth = 0.36,
      color = "grey25"
    ) +
    ggplot2::geom_text(
      data = annotation_data,
      ggplot2::aes(
        x = 1.5,
        y = y_position + 0.05 * pmax(y_max - y_min, 1),
        label = label
      ),
      inherit.aes = FALSE,
      size = 2.25,
      lineheight = 0.90,
      color = "grey20"
    ) +
    ggplot2::scale_fill_manual(
      values = mofa_group_colors[c("pCR", "non_pCR")],
      guide = "none",
      drop = FALSE
    ) +
    ggplot2::facet_wrap(~ factor, nrow = 2, scales = "free_y") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.30))) +
    ggplot2::labs(title = title_text, x = NULL, y = "MOFA factor score") +
    ggplot2::theme_classic(base_size = 8.8) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "plain", size = 8.0),
      axis.text.x = ggplot2::element_text(size = 9.0, face = "plain"),
      panel.spacing = grid::unit(0.75, "lines")
    )
}

mofa_v50_save_violin_variants <- function(
  plot_object,
  base_filename,
  width_current,
  height_current
) {
  for (width_multiplier in c(0.80, 0.70, 0.60)) {
    ggplot2::ggsave(
      filename = file.path(
        mofa_v50_figure_dir,
        paste0(
          tools::file_path_sans_ext(base_filename),
          "_w",
          sprintf("%02d", round(100 * width_multiplier)),
          ".svg"
        )
      ),
      plot = plot_object,
      width = width_current * width_multiplier,
      height = height_current,
      units = "in",
      device = svglite::svglite,
      bg = "white"
    )
  }
}

mofa_v50_overall_annotation <- mofa_v50_make_annotation(
  mofa_v50_overall_response_long,
  "overall_response_subject_mean"
)
mofa_v50_baseline_annotation <- mofa_v50_make_annotation(
  mofa_v50_baseline_response_long,
  "baseline_response"
)
mofa_v50_after_annotation <- mofa_v50_make_annotation(
  mofa_v50_after_response_long,
  "ongoing_response"
)

p_mofa_v50_overall_response <- mofa_v50_build_response_plot(
  mofa_v50_overall_response_long,
  mofa_v50_overall_annotation,
  "pCR versus non-pCR: subject mean across available timepoints"
)

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_overall_subject_mean.svg"),
  plot = p_mofa_v50_overall_response,
  width = 7.15,
  height = max(5.4, 2.25 * ceiling(length(mofa_v50_factor_order) / 5)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

mofa_v50_save_violin_variants(
  p_mofa_v50_overall_response,
  "MOFA_v50_factor_response_overall_subject_mean.svg",
  7.15,
  max(5.4, 2.25 * ceiling(length(mofa_v50_factor_order) / 5))
)

p_mofa_v50_baseline_response <- mofa_v50_build_response_plot(
  mofa_v50_baseline_response_long,
  mofa_v50_baseline_annotation,
  "pCR versus non-pCR at Baseline"
)

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_baseline.svg"),
  plot = p_mofa_v50_baseline_response,
  width = 7.15,
  height = max(5.4, 2.25 * ceiling(length(mofa_v50_factor_order) / 5)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

mofa_v50_save_violin_variants(
  p_mofa_v50_baseline_response,
  "MOFA_v50_factor_response_baseline.svg",
  7.15,
  max(5.4, 2.25 * ceiling(length(mofa_v50_factor_order) / 5))
)

p_mofa_v50_after_response <- mofa_v50_build_response_plot(
  mofa_v50_after_response_long,
  mofa_v50_after_annotation,
  "pCR versus non-pCR after radiotherapy"
)

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_after_RT.svg"),
  plot = p_mofa_v50_after_response,
  width = 7.15,
  height = max(5.4, 2.25 * ceiling(length(mofa_v50_factor_order) / 5)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

mofa_v50_save_violin_variants(
  p_mofa_v50_after_response,
  "MOFA_v50_factor_response_after_RT.svg",
  7.15,
  max(5.4, 2.25 * ceiling(length(mofa_v50_factor_order) / 5))
)

p_mofa_v50_paired_change <- ggplot2::ggplot(
  mofa_v50_paired_long,
  ggplot2::aes(x = Timepoint_display, y = factor_score, group = SubjectID)
) +
  ggplot2::geom_line(color = "grey72", linewidth = 0.42, alpha = 0.75) +
  ggplot2::geom_point(
    ggplot2::aes(fill = Timepoint_display),
    shape = 21,
    size = 1.75,
    color = "grey20",
    stroke = 0.38
  ) +
  ggplot2::scale_fill_manual(values = timepoint_display_colors, guide = "none") +
  ggplot2::facet_wrap(~ factor, ncol = 4, scales = "free_y") +
  ggplot2::labs(
    title = "Within-subject paired factor change",
    subtitle = "Lines and points are intentionally not coloured by response group.",
    x = NULL,
    y = "MOFA factor score"
  ) +
  ggplot2::theme_classic(base_size = 8.8) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.7, color = "grey35"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "plain", size = 8.0),
    panel.spacing = grid::unit(0.75, "lines")
  )

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_paired_change.svg"),
  plot = p_mofa_v50_paired_change,
  width = 10.2,
  height = max(5.4, 2.25 * ceiling(length(mofa_v50_factor_order) / 4)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

mofa_v50_delta_annotation <- mofa_v50_paired_factor_data %>%
  dplyr::group_by(factor) %>%
  dplyr::summarise(
    y_min = min(factor_delta, na.rm = TRUE),
    y_max = max(factor_delta, na.rm = TRUE),
    delta_difference =
      median(factor_delta[TRG_plot == "pCR"], na.rm = TRUE) -
      median(factor_delta[TRG_plot == "non_pCR"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    mofa_factor_wilcoxon_tests %>%
      dplyr::filter(comparison == "paired_differential_change") %>%
      dplyr::transmute(
        factor = factor(as.character(factor), levels = mofa_v50_factor_order),
        wilcoxon_p
      ),
    by = "factor"
  ) %>%
  dplyr::mutate(
    y_position = y_max + 0.14 * pmax(y_max - y_min, 1),
    label = paste0(
      "Delta difference = ", sprintf("%+.2f", delta_difference),
      "\nWilcoxon ",
      ifelse(
        !is.finite(wilcoxon_p),
        "p = NA",
        ifelse(wilcoxon_p < 0.001, "p < 0.001", paste0("p = ", mofa_v50_format_p(wilcoxon_p)))
      )
    )
  )

p_mofa_v50_paired_change_by_response <- ggplot2::ggplot(
  mofa_v50_paired_factor_data,
  ggplot2::aes(x = TRG_plot, y = factor_delta, fill = TRG_plot)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.40) +
  ggplot2::geom_violin(trim = FALSE, alpha = 0.68, color = "grey30", linewidth = 0.38) +
  ggplot2::geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    fill = NA,
    coef = 0,
    staplewidth = 0,
    linewidth = 0.42
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.08, height = 0),
    shape = 21,
    size = 1.55,
    color = "grey20",
    stroke = 0.32
  ) +
  ggplot2::geom_segment(
    data = mofa_v50_delta_annotation,
    ggplot2::aes(x = 1, xend = 2, y = y_position, yend = y_position),
    inherit.aes = FALSE,
    linewidth = 0.36,
    color = "grey25"
  ) +
  ggplot2::geom_text(
    data = mofa_v50_delta_annotation,
    ggplot2::aes(
      x = 1.5,
      y = y_position + 0.05 * pmax(y_max - y_min, 1),
      label = label
    ),
    inherit.aes = FALSE,
    size = 2.10,
    lineheight = 0.90
  ) +
  ggplot2::scale_fill_manual(
    values = mofa_group_colors[c("pCR", "non_pCR")],
    guide = "none",
    drop = FALSE
  ) +
  ggplot2::facet_wrap(~ factor, nrow = 2, scales = "free_y") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.30))) +
  ggplot2::labs(
    title = "Response-specific paired treatment change",
    subtitle = "Delta is After RT minus Baseline for each subject; positive effects indicate a larger increase in pCR.",
    x = NULL,
    y = "Within-subject factor delta"
  ) +
  ggplot2::theme_classic(base_size = 8.8) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.7, color = "grey35"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "plain", size = 8.0),
    axis.text.x = ggplot2::element_text(size = 9.0, face = "plain"),
    panel.spacing = grid::unit(0.75, "lines")
  )

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_paired_change_by_response.svg"),
  plot = p_mofa_v50_paired_change_by_response,
  width = 7.15,
  height = max(5.4, 2.30 * ceiling(length(mofa_v50_factor_order) / 5)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

mofa_v50_save_violin_variants(
  p_mofa_v50_paired_change_by_response,
  "MOFA_v50_factor_paired_change_by_response.svg",
  7.15,
  max(5.4, 2.30 * ceiling(length(mofa_v50_factor_order) / 5))
)

mofa_v50_paired_effect_data <- mofa_factor_wilcoxon_tests %>%
  dplyr::filter(
    comparison == "paired_differential_change",
    factor %in% mofa_v50_factor_order
  ) %>%
  dplyr::mutate(
    factor = factor(factor, levels = rev(mofa_v50_factor_order)),
    minus_log10_p = mofa_v50_neglog10_p(wilcoxon_p)
  )

mofa_v50_paired_effect_p_upper <- mofa_v50_neglog10_upper(
  mofa_v50_paired_effect_data$minus_log10_p
)
mofa_v50_paired_effect_p_breaks <- mofa_v50_neglog10_breaks(
  mofa_v50_paired_effect_p_upper
)

p_mofa_v50_paired_change_effects <- ggplot2::ggplot(
  mofa_v50_paired_effect_data,
  ggplot2::aes(x = wilcoxon_effect, y = factor, fill = minus_log10_p)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.45) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = wilcoxon_effect, yend = factor),
    color = "grey65",
    linewidth = 0.55
  ) +
  ggplot2::geom_point(shape = 21, size = 3.4, color = "grey20", stroke = 0.55) +
  ggplot2::scale_fill_gradient(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "white",
    high = "#2166AC",
    limits = c(0, mofa_v50_paired_effect_p_upper),
    breaks = mofa_v50_paired_effect_p_breaks,
    oob = scales::squish,
    name = expression(-log[10](italic(p)))
  ) +
  ggplot2::labs(
    title = "Differential paired-change effects",
    subtitle = "Effect = median delta in pCR minus median delta in non-pCR.",
    x = "Delta difference: pCR minus non-pCR",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 9.2) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.7, color = "grey35"),
    axis.text.y = ggplot2::element_text(face = "plain", size = 7.5)
  )

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_paired_change_effects.svg"),
  plot = p_mofa_v50_paired_change_effects,
  width = 6.0,
  height = max(4.5, 1.8 + 0.27 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

#-----------------------------------------------------------------#
# 22.4 Non-duplicated response-selected factor-pair maps
#-----------------------------------------------------------------#

mofa_v50_pair_atlas <- mofa_v50_pair_permanova_atlas %>%
  dplyr::mutate(
    factor_x = as.character(factor_x),
    factor_y = as.character(factor_y),
    pair_id = paste(pmin(factor_x, factor_y), pmax(factor_x, factor_y), sep = "__")
  )

mofa_v50_baseline_pair <- mofa_v50_pair_atlas %>%
  dplyr::filter(scope == "Baseline only", is.finite(permanova_p)) %>%
  dplyr::arrange(permanova_p, dplyr::desc(permanova_r2)) %>%
  dplyr::slice_head(n = 1)

mofa_v50_overall_pair_candidates <- mofa_v50_pair_atlas %>%
  dplyr::filter(scope == "Overall subject mean", is.finite(permanova_p)) %>%
  dplyr::arrange(permanova_p, dplyr::desc(permanova_r2))

if (nrow(mofa_v50_baseline_pair) > 0) {
  mofa_v50_overall_pair_nonduplicate <- mofa_v50_overall_pair_candidates %>%
    dplyr::filter(pair_id != mofa_v50_baseline_pair$pair_id[1])
  mofa_v50_overall_pair <- if (nrow(mofa_v50_overall_pair_nonduplicate) > 0) {
    mofa_v50_overall_pair_nonduplicate %>% dplyr::slice_head(n = 1)
  } else {
    mofa_v50_overall_pair_candidates %>% dplyr::slice_head(n = 1)
  }
} else {
  mofa_v50_overall_pair <- mofa_v50_overall_pair_candidates %>% dplyr::slice_head(n = 1)
}

mofa_v50_selected_pair_rows <- dplyr::bind_rows(
  mofa_v50_baseline_pair,
  mofa_v50_overall_pair
)

write.csv(
  mofa_v50_pair_atlas,
  file.path(mofa_v50_result_dir, "MOFA_v50_factor_pair_PERMANOVA_atlas.csv"),
  row.names = FALSE
)
mofa_v50_subject_mean_factor_data <- mofa_factor_scores %>%
  dplyr::group_by(SubjectID, TRG_plot) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(mofa_v50_factor_order),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

mofa_v50_baseline_factor_data <- mofa_factor_scores %>%
  dplyr::filter(as.character(Timepoint) == "Before") %>%
  dplyr::group_by(SubjectID, TRG_plot) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(mofa_v50_factor_order),
      ~ dplyr::first(.x)
    ),
    .groups = "drop"
  )

mofa_v50_build_pair_map <- function(pair_row) {
  factor_x <- as.character(pair_row$factor_x[1])
  factor_y <- as.character(pair_row$factor_y[1])
  scope_label <- as.character(pair_row$scope[1])

  pair_data <- (
    if (scope_label == "Baseline only") {
      mofa_v50_baseline_factor_data
    } else {
      mofa_v50_subject_mean_factor_data
    }
  ) %>%
    dplyr::transmute(
      SubjectID,
      TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR")),
      factor_x_value = .data[[factor_x]],
      factor_y_value = .data[[factor_y]]
    ) %>%
    dplyr::filter(
      is.finite(factor_x_value),
      is.finite(factor_y_value),
      !is.na(TRG_plot)
    ) %>%
    dplyr::mutate(
      factor_x_z = as.numeric(scale(factor_x_value)),
      factor_y_z = as.numeric(scale(factor_y_value))
    )

  hull_data <- pair_data %>%
    dplyr::group_by(TRG_plot) %>%
    dplyr::filter(dplyr::n() >= 3) %>%
    dplyr::slice(chull(factor_x_z, factor_y_z)) %>%
    dplyr::ungroup()

  permanova_label <- if (
    is.finite(pair_row$permanova_p[1]) &&
    pair_row$permanova_p[1] < 0.001
  ) {
    paste0(
      "PERMANOVA~R^2==", sprintf("%.2f", pair_row$permanova_r2[1]),
      "*','~~italic(p)<0.001*','~~italic(n)==", pair_row$n_subjects[1]
    )
  } else {
    paste0(
      "PERMANOVA~R^2==", sprintf("%.2f", pair_row$permanova_r2[1]),
      "*','~~italic(p)==", sprintf("%.3f", pair_row$permanova_p[1]),
      "*','~~italic(n)==", pair_row$n_subjects[1]
    )
  }

  ggplot2::ggplot(pair_data, ggplot2::aes(x = factor_x_z, y = factor_y_z)) +
    ggplot2::geom_polygon(
      data = hull_data,
      ggplot2::aes(group = TRG_plot, fill = TRG_plot, color = TRG_plot),
      alpha = 0.11,
      linewidth = 0.50,
      show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey72", linewidth = 0.32) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey72", linewidth = 0.32) +
    ggplot2::geom_point(
      ggplot2::aes(fill = TRG_plot),
      shape = 21,
      size = 2.05,
      color = "grey20",
      stroke = 0.45,
      alpha = 0.92
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = TRG_plot),
      fun = mean,
      geom = "point",
      shape = 21,
      size = 2.85,
      color = "black",
      stroke = 0.62
    ) +
    ggplot2::annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = permanova_label,
      parse = TRUE,
      hjust = -0.04,
      vjust = 1.20,
      size = 3.0,
      color = "black"
    ) +
    ggplot2::scale_fill_manual(
      values = mofa_group_colors[c("pCR", "non_pCR")],
      name = "Response"
    ) +
    ggplot2::scale_color_manual(
      values = mofa_group_colors[c("pCR", "non_pCR")],
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.08))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.10))) +
    ggplot2::labs(
      title = paste0(scope_label, ": ", factor_x, " versus ", factor_y),
      x = paste0(factor_x, " (z-score within ", scope_label, ")"),
      y = paste0(factor_y, " (z-score within ", scope_label, ")")
    ) +
    ggplot2::theme_classic(base_size = 9.0) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 9.5),
      legend.position = "bottom",
      aspect.ratio = 1
    )
}

p_mofa_v50_baseline_pair_map <- if (nrow(mofa_v50_baseline_pair) > 0) {
  mofa_v50_build_pair_map(mofa_v50_baseline_pair)
} else {
  NULL
}
p_mofa_v50_overall_pair_map <- if (nrow(mofa_v50_overall_pair) > 0) {
  mofa_v50_build_pair_map(mofa_v50_overall_pair)
} else {
  NULL
}

if (!is.null(p_mofa_v50_baseline_pair_map)) {
  ggplot2::ggsave(
    filename = file.path(mofa_v50_figure_dir, "MOFA_v50_selected_pair_map_baseline.svg"),
    plot = p_mofa_v50_baseline_pair_map,
    width = 4.6,
    height = 4.4,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

if (!is.null(p_mofa_v50_overall_pair_map)) {
  ggplot2::ggsave(
    filename = file.path(mofa_v50_figure_dir, "MOFA_v50_selected_pair_map_overall_subject_mean.svg"),
    plot = p_mofa_v50_overall_pair_map,
    width = 4.6,
    height = 4.4,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

mofa_v50_pair_plots <- Filter(Negate(is.null), list(
  p_mofa_v50_baseline_pair_map,
  p_mofa_v50_overall_pair_map
))
if (length(mofa_v50_pair_plots) > 0) {
  p_mofa_v50_selected_pair_maps <- patchwork::wrap_plots(
    mofa_v50_pair_plots,
    ncol = length(mofa_v50_pair_plots),
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = "Response-selected MOFA factor-pair maps",
      subtitle = "Baseline and overall maps use different subject-level estimands and non-duplicated factor pairs when available."
    ) &
    ggplot2::theme(legend.position = "bottom")

  ggplot2::ggsave(
    filename = file.path(mofa_v50_figure_dir, "MOFA_v50_selected_response_pair_maps.svg"),
    plot = p_mofa_v50_selected_pair_maps,
    width = 4.7 * length(mofa_v50_pair_plots),
    height = 4.7,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

#-----------------------------------------------------------------#
# 22.5 Selected-factor panels with full feature labels
#-----------------------------------------------------------------#

mofa_v50_selected_factor_candidates <- mofa_factor_wilcoxon_tests %>%
  dplyr::filter(
    comparison == "baseline_response",
    factor %in% mofa_v50_factor_order
  ) %>%
  dplyr::left_join(
    mofa_v50_view_dominance_audit %>%
      dplyr::transmute(
        factor = as.character(factor),
        global_level_flag,
        strongest_view_fraction,
        sharedness_penalized_r2
      ),
    by = "factor"
  ) %>%
  dplyr::arrange(
    global_level_flag,
    strongest_view_fraction > 0.80,
    wilcoxon_p,
    dplyr::desc(abs(wilcoxon_effect)),
    dplyr::desc(sharedness_penalized_r2)
  )

mofa_v50_selected_factor <- if (
  exists("mofa_selected_factor", inherits = FALSE) &&
  as.character(mofa_selected_factor) %in% mofa_v50_factor_order
) {
  as.character(mofa_selected_factor)
} else if (nrow(mofa_v50_selected_factor_candidates) > 0) {
  as.character(mofa_v50_selected_factor_candidates$factor[1])
} else {
  mofa_v50_factor_order[1]
}

mofa_v50_selected_factor_r2 <- mofa_variance_explained %>%
  dplyr::filter(as.character(factor) == mofa_v50_selected_factor) %>%
  dplyr::mutate(
    view_label = factor(as.character(view_label), levels = rev(unname(view_labels)))
  )

p_mofa_v50_selected_factor_r2 <- ggplot2::ggplot(
  mofa_v50_selected_factor_r2,
  ggplot2::aes(x = r2, y = view_label, fill = view_label)
) +
  ggplot2::geom_col(width = 0.62, color = "grey25", linewidth = 0.35) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.1f%%", 100 * r2)),
    hjust = -0.08,
    size = 2.7
  ) +
  ggplot2::scale_fill_manual(values = view_colors, guide = "none") +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
  ggplot2::labs(
    title = paste0(mofa_v50_selected_factor, ": variance explained"),
    x = "Variance explained",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 9.0) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_selected_factor_R2.svg"),
  plot = p_mofa_v50_selected_factor_r2,
  width = 4.8,
  height = 2.8,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

mofa_v50_selected_baseline_data <- mofa_v50_baseline_response_long %>%
  dplyr::filter(as.character(factor) == mofa_v50_selected_factor)

mofa_v50_selected_baseline_annotation <- mofa_v50_baseline_annotation %>%
  dplyr::filter(as.character(factor) == mofa_v50_selected_factor)

p_mofa_v50_selected_factor_baseline <- ggplot2::ggplot(
  mofa_v50_selected_baseline_data,
  ggplot2::aes(x = TRG_plot, y = factor_score, fill = TRG_plot)
) +
  ggplot2::geom_violin(trim = FALSE, alpha = 0.68, color = "grey30") +
  ggplot2::geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    fill = NA,
    coef = 0,
    staplewidth = 0,
    linewidth = 0.42
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.08),
    shape = 21,
    size = 1.8,
    color = "grey20"
  ) +
  ggplot2::geom_segment(
    data = mofa_v50_selected_baseline_annotation,
    ggplot2::aes(x = 1, xend = 2, y = y_position, yend = y_position),
    inherit.aes = FALSE,
    linewidth = 0.38,
    color = "grey25"
  ) +
  ggplot2::geom_text(
    data = mofa_v50_selected_baseline_annotation,
    ggplot2::aes(
      x = 1.5,
      y = y_position + 0.06 * pmax(y_max - y_min, 1),
      label = label
    ),
    inherit.aes = FALSE,
    size = 2.45,
    lineheight = 0.90
  ) +
  ggplot2::scale_fill_manual(values = mofa_group_colors[c("pCR", "non_pCR")], guide = "none") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.30))) +
  ggplot2::labs(
    title = paste0(mofa_v50_selected_factor, ": Baseline response contrast"),
    x = NULL,
    y = "Factor score"
  ) +
  ggplot2::theme_classic(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(size = 10.0, face = "plain")
  )

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_selected_factor_baseline.svg"),
  plot = p_mofa_v50_selected_factor_baseline,
  width = 3.0,
  height = 3.5,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

mofa_v50_save_violin_variants(
  p_mofa_v50_selected_factor_baseline,
  "MOFA_v50_selected_factor_baseline.svg",
  3.0,
  3.5
)

mofa_v50_selected_paired_data <- mofa_v50_paired_long %>%
  dplyr::filter(as.character(factor) == mofa_v50_selected_factor)

p_mofa_v50_selected_factor_paired <- ggplot2::ggplot(
  mofa_v50_selected_paired_data,
  ggplot2::aes(x = Timepoint_display, y = factor_score, group = SubjectID)
) +
  ggplot2::geom_line(color = "grey70", linewidth = 0.48) +
  ggplot2::geom_point(
    ggplot2::aes(fill = Timepoint_display),
    shape = 21,
    size = 2.0,
    color = "grey20"
  ) +
  ggplot2::scale_fill_manual(values = timepoint_display_colors, guide = "none") +
  ggplot2::labs(
    title = paste0(mofa_v50_selected_factor, ": paired change"),
    x = NULL,
    y = "Factor score"
  ) +
  ggplot2::theme_classic(base_size = 9.0) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_selected_factor_paired.svg"),
  plot = p_mofa_v50_selected_factor_paired,
  width = 3.5,
  height = 3.2,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

mofa_v50_selected_loading_data <- mofa_feature_weights %>%
  dplyr::filter(
    factor == mofa_v50_selected_factor,
    display_eligible
  ) %>%
  dplyr::inner_join(
    mofa_v50_factor_view_display_audit %>%
      dplyr::filter(
        factor == mofa_v50_selected_factor,
        feature_display_active
      ) %>%
      dplyr::select(factor, view, view_r2 = r2),
    by = c("factor", "view")
  ) %>%
  dplyr::group_by(view, direction) %>%
  dplyr::slice_max(order_by = abs(weight_within_view), n = 5, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    view_label = factor(view_label, levels = unname(view_labels)),
    feature_plot_id = paste(view, feature, sep = "::"),
    feature_label_plotmath = mofa_v50_plotmath_feature_label(view, feature_label)
  ) %>%
  dplyr::arrange(
    match(view, required_views),
    dplyr::desc(weight_within_view)
  )

mofa_v50_selected_loading_data$feature_plot_id <- factor(
  mofa_v50_selected_loading_data$feature_plot_id,
  levels = rev(unique(mofa_v50_selected_loading_data$feature_plot_id))
)

p_mofa_v50_selected_factor_loadings <- ggplot2::ggplot(
  mofa_v50_selected_loading_data,
  ggplot2::aes(x = weight_within_view, y = feature_plot_id, fill = view_label)
) +
  ggplot2::geom_vline(xintercept = 0, color = "grey65", linewidth = 0.40) +
  ggplot2::geom_col(width = 0.66, color = "grey25", linewidth = 0.25) +
  ggplot2::facet_grid(view_label ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_y_discrete(
    labels = function(x) {
      parse(
        text = stats::setNames(
          mofa_v50_selected_loading_data$feature_label_plotmath,
          as.character(mofa_v50_selected_loading_data$feature_plot_id)
        )[x]
      )
    }
  ) +
  ggplot2::scale_x_continuous(
    limits = c(-1.05, 1.05),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::scale_fill_manual(values = view_colors, guide = "none") +
  ggplot2::labs(
    title = paste0(mofa_v50_selected_factor, ": strongest feature weights"),
    subtitle = paste0(
      "Absolute weight indicates the strength of feature-factor association; ",
      "the sign indicates direction. Features are displayed only for views with R² >= ",
      scales::percent(mofa_active_view_r2, accuracy = 1), "."
    ),
    x = "Scaled MOFA feature weight",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.7, color = "grey35"),
    strip.background = ggplot2::element_blank(),
    strip.text.y = ggplot2::element_text(face = "bold", angle = 0),
    axis.text.y = ggplot2::element_text(size = 7.5),
    panel.spacing.y = grid::unit(0.65, "lines")
  )

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_selected_factor_feature_loadings.svg"),
  plot = p_mofa_v50_selected_factor_loadings,
  width = 9.2,
  height = max(7.0, 0.21 * nrow(mofa_v50_selected_loading_data) + 3.0),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

p_mofa_v50_selected_factor_summary <-
  (p_mofa_v50_selected_factor_r2 | p_mofa_v50_selected_factor_baseline | p_mofa_v50_selected_factor_paired) /
  p_mofa_v50_selected_factor_loadings +
  patchwork::plot_layout(heights = c(0.42, 1.0)) +
  patchwork::plot_annotation(
    title = paste0("Selected MOFA factor summary: ", mofa_v50_selected_factor)
  )

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_selected_factor_summary.svg"),
  plot = p_mofa_v50_selected_factor_summary,
  width = 10.8,
  height = max(10.2, 0.21 * nrow(mofa_v50_selected_loading_data) + 6.5),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

#-----------------------------------------------------------------#
# 22.6 Compact radial factor-feature network
#-----------------------------------------------------------------#

# Publication selection rule for network factors:
#   1) at least one active view (R² >= 2%); multi-view factors are preferred;
#   2) no sample-total/global-level QC flag;
#   3) evidence from the factor-level pCR/non-pCR contrast, an active-view
#      representative feature, or both.
# Five factors are retained so the network is informative without becoming a
# complete all-factor hairball. Representative response features are forced into
# the network even when they are not the largest absolute loading.
mofa_v50_network_response_rank <- mofa_factor_wilcoxon_tests %>%
  dplyr::filter(
    factor %in% mofa_v50_factor_order,
    comparison %in% c(
      "overall_response_subject_mean",
      "baseline_response",
      "ongoing_response",
      "paired_differential_change"
    )
  ) %>%
  dplyr::group_by(factor) %>%
  dplyr::arrange(wilcoxon_p, dplyr::desc(abs(wilcoxon_effect)), .by_group = TRUE) %>%
  dplyr::summarise(
    minimum_response_p = dplyr::first(wilcoxon_p[is.finite(wilcoxon_p)], default = NA_real_),
    minimum_response_fdr = dplyr::first(wilcoxon_fdr[is.finite(wilcoxon_fdr)], default = NA_real_),
    strongest_response_effect = dplyr::first(wilcoxon_effect[is.finite(wilcoxon_p)], default = NA_real_),
    strongest_response_comparison = dplyr::first(comparison[is.finite(wilcoxon_p)], default = NA_character_),
    .groups = "drop"
  )

mofa_v50_network_feature_rank <- mofa_feature_weights %>%
  dplyr::filter(
    factor %in% mofa_v50_factor_order,
    display_eligible,
    is.finite(p_value)
  ) %>%
  dplyr::inner_join(
    mofa_v50_factor_view_display_audit %>%
      dplyr::filter(feature_display_active) %>%
      dplyr::select(factor, view, view_r2 = r2),
    by = c("factor", "view")
  ) %>%
  dplyr::mutate(
    feature_response_score = abs(weight_within_view) *
      pmin(mofa_v50_neglog10_p(p_value), 5)
  ) %>%
  dplyr::group_by(factor) %>%
  dplyr::arrange(
    dplyr::desc(is.finite(fdr) & fdr < 0.10),
    dplyr::desc(feature_response_score),
    p_value,
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    factor,
    best_response_feature_view = view,
    best_response_feature = feature,
    best_response_feature_label = feature_label,
    best_response_feature_weight = weight_within_view,
    best_response_feature_effect = response_effect,
    best_response_feature_p = p_value,
    best_response_feature_fdr = fdr,
    best_response_feature_score = feature_response_score
  )

mofa_v50_network_selection_audit <- mofa_v50_view_dominance_audit %>%
  dplyr::mutate(factor = as.character(factor)) %>%
  dplyr::filter(factor %in% mofa_v50_factor_order) %>%
  dplyr::left_join(mofa_v50_network_response_rank, by = "factor") %>%
  dplyr::left_join(mofa_v50_network_feature_rank, by = "factor") %>%
  dplyr::mutate(
    eligible_active_view = active_views_primary >= 1,
    eligible_multiview = active_views_primary >= 2,
    eligible_qc = !dplyr::coalesce(global_level_flag, FALSE),
    factor_response_supported = is.finite(minimum_response_p) & minimum_response_p < 0.10,
    feature_response_supported =
      (is.finite(best_response_feature_fdr) & best_response_feature_fdr < 0.10) |
      (
        is.finite(best_response_feature_p) & best_response_feature_p < 0.05 &
          is.finite(best_response_feature_weight) & abs(best_response_feature_weight) >= 0.35
      ),
    selection_tier = dplyr::case_when(
      eligible_multiview & eligible_qc & factor_response_supported & feature_response_supported ~ 1L,
      eligible_multiview & eligible_qc & factor_response_supported ~ 2L,
      eligible_multiview & eligible_qc & feature_response_supported ~ 3L,
      eligible_active_view & eligible_qc & (factor_response_supported | feature_response_supported) ~ 4L,
      eligible_multiview & eligible_qc ~ 5L,
      eligible_active_view & eligible_qc ~ 6L,
      TRUE ~ 7L
    )
  ) %>%
  dplyr::arrange(
    selection_tier,
    minimum_response_p,
    best_response_feature_p,
    dplyr::desc(sharedness_penalized_r2),
    match(factor, mofa_v50_factor_order)
  )

mofa_v50_network_factors <- utils::head(
  mofa_v50_network_selection_audit$factor[
    mofa_v50_network_selection_audit$eligible_active_view &
      mofa_v50_network_selection_audit$eligible_qc
  ],
  min(mofa_network_factor_count, length(mofa_v50_factor_order))
)

if (length(mofa_v50_network_factors) < 2) {
  stop("Fewer than two QC-passing active factors are available for the multi-factor network.", call. = FALSE)
}

mofa_v50_network_selection_audit <- mofa_v50_network_selection_audit %>%
  dplyr::mutate(
    selected_for_network = factor %in% mofa_v50_network_factors,
    selection_reason = dplyr::case_when(
      selected_for_network & selection_tier == 1L ~ "Factor- and feature-level response evidence; multi-view",
      selected_for_network & selection_tier == 2L ~ "Factor-level response evidence; multi-view",
      selected_for_network & selection_tier == 3L ~ "Representative-feature response evidence; multi-view",
      selected_for_network & selection_tier == 4L ~ "Factor- or feature-level response evidence; active view",
      selected_for_network & selection_tier == 5L ~ "Fallback: QC-passing multi-view factor",
      selected_for_network ~ "Fallback: QC-passing active factor",
      TRUE ~ "Not selected"
    )
  )

write.csv(
  mofa_v50_network_selection_audit,
  file.path(mofa_v50_result_dir, "MOFA_v50_network_factor_selection_audit.csv"),
  row.names = FALSE
)

mofa_v50_active_factor_views <- mofa_variance_explained %>%
  dplyr::mutate(factor = as.character(factor), view = as.character(view)) %>%
  dplyr::filter(
    factor %in% mofa_v50_network_factors,
    r2 >= mofa_active_view_r2
  ) %>%
  dplyr::select(factor, view, view_r2 = r2)


mofa_v50_network_seed_edges <- dplyr::bind_rows(
  mofa_feature_weights %>%
    dplyr::filter(
      factor %in% mofa_v50_network_factors,
      display_eligible
    ) %>%
    dplyr::inner_join(mofa_v50_active_factor_views, by = c("factor", "view")) %>%
    dplyr::group_by(factor, view, direction) %>%
    dplyr::slice_max(
      order_by = abs(weight_within_view),
      n = mofa_network_features_per_view_direction,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup(),
  mofa_feature_weights %>%
    dplyr::inner_join(
      mofa_v50_network_selection_audit %>%
        dplyr::filter(
          selected_for_network,
          is.finite(best_response_feature_p)
        ) %>%
        dplyr::select(
          factor,
          view = best_response_feature_view,
          feature = best_response_feature
        ),
      by = c("factor", "view", "feature")
    ) %>%
    dplyr::inner_join(mofa_v50_active_factor_views, by = c("factor", "view"))
) %>%
  dplyr::distinct(factor, view, feature, .keep_all = TRUE) %>%
  dplyr::mutate(
    feature_node_id = paste(view, feature, sep = "::"),
    edge_id = paste(factor, feature_node_id, sep = "__")
  )

mofa_v50_network_edges <- mofa_feature_weights %>%
  dplyr::filter(
    factor %in% mofa_v50_network_factors,
    display_eligible
  ) %>%
  dplyr::inner_join(mofa_v50_active_factor_views, by = c("factor", "view")) %>%
  dplyr::mutate(
    feature_node_id = paste(view, feature, sep = "::"),
    edge_id = paste(factor, feature_node_id, sep = "__")
  ) %>%
  dplyr::filter(
    feature_node_id %in% unique(mofa_v50_network_seed_edges$feature_node_id),
    abs(weight_within_view) >= mofa_network_shared_loading_threshold |
      edge_id %in% mofa_v50_network_seed_edges$edge_id
  ) %>%
  dplyr::mutate(
    loading_sign = ifelse(weight >= 0, "Positive", "Negative"),
    loading_strength = abs(weight_within_view)
  )

if (nrow(mofa_v50_network_edges) == 0) {
  stop("No feature edges passed the multi-factor network display rule.", call. = FALSE)
}

# Modules are data-derived clusters of signed feature-weight profiles within each view.
mofa_v50_module_membership <- dplyr::bind_rows(
  lapply(
    split(mofa_v50_network_edges, mofa_v50_network_edges$view),
    function(view_edge_data) {
      profile <- view_edge_data %>%
        dplyr::select(feature_node_id, factor, weight_within_view) %>%
        tidyr::pivot_wider(
          names_from = factor,
          values_from = weight_within_view,
          values_fill = 0
        )

      for (factor_name in setdiff(mofa_v50_network_factors, colnames(profile))) {
        profile[[factor_name]] <- 0
      }

      profile_matrix <- as.matrix(
        profile[, mofa_v50_network_factors, drop = FALSE]
      )
      rownames(profile_matrix) <- profile$feature_node_id

      module <- if (nrow(profile_matrix) <= 2) {
        rep(1L, nrow(profile_matrix))
      } else {
        stats::cutree(
          stats::hclust(stats::dist(profile_matrix), method = "ward.D2"),
          k = min(2L, max(1L, ceiling(nrow(profile_matrix) / 7)))
        )
      }

      data.frame(
        feature_node_id = rownames(profile_matrix),
        module = as.integer(module),
        stringsAsFactors = FALSE
      )
    }
  )
)

mofa_v50_network_feature_nodes <- mofa_v50_network_edges %>%
  dplyr::group_by(feature_node_id, view, view_label, feature, feature_label) %>%
  dplyr::summarise(
    n_connected_factors = dplyr::n_distinct(factor),
    maximum_loading = max(loading_strength, na.rm = TRUE),
    minimum_response_p = if (any(is.finite(p_value))) min(p_value, na.rm = TRUE) else NA_real_,
    minimum_response_fdr = if (any(is.finite(fdr))) min(fdr, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  dplyr::left_join(mofa_v50_module_membership, by = "feature_node_id") %>%
  dplyr::mutate(
    view_order = match(view, required_views),
    module_label = paste0(as.character(view_label), " M", module),
    shared_feature = n_connected_factors >= 2,
    response_feature = is.finite(minimum_response_p) & minimum_response_p < 0.05,
    feature_label_plotmath = mofa_v50_plotmath_feature_label(view, feature_label)
  ) %>%
  dplyr::arrange(
    view_order,
    module,
    dplyr::desc(n_connected_factors),
    dplyr::desc(maximum_loading)
  )

mofa_v50_network_modules <- mofa_v50_network_feature_nodes %>%
  dplyr::distinct(view, view_label, view_order, module, module_label) %>%
  dplyr::arrange(view_order, module) %>%
  dplyr::mutate(
    module_index = dplyr::row_number(),
    module_angle = pi / 2 - 2 * pi * (module_index - 1) / dplyr::n(),
    module_x = 2.45 * cos(module_angle),
    module_y = 2.05 * sin(module_angle)
  )

mofa_v50_network_feature_nodes <- mofa_v50_network_feature_nodes %>%
  dplyr::left_join(
    mofa_v50_network_modules %>%
      dplyr::select(module_label, module_x, module_y),
    by = "module_label"
  ) %>%
  dplyr::group_by(module_label) %>%
  dplyr::mutate(
    feature_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.38,
    x = module_x,
    y = module_y + feature_offset,
    label_x = x + ifelse(module_x < -0.25, -0.52, 0.52),
    label_hjust = ifelse(module_x < -0.25, 1, 0)
  ) %>%
  dplyr::ungroup()

mofa_v50_network_modules <- mofa_v50_network_modules %>%
  dplyr::left_join(
    mofa_v50_network_feature_nodes %>%
      dplyr::group_by(module_label) %>%
      dplyr::summarise(
        xmin = min(x) - 0.34,
        xmax = max(x) + 0.34,
        ymin = min(y) - 0.23,
        ymax = max(y) + 0.23,
        label_y = max(y) + 0.42,
        .groups = "drop"
      ),
    by = "module_label"
  )

mofa_v50_network_factor_nodes <- data.frame(
  factor = mofa_v50_network_factors,
  factor_angle = pi / 2 -
    2 * pi * (seq_along(mofa_v50_network_factors) - 1) /
    length(mofa_v50_network_factors),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    mofa_v50_network_selection_audit %>%
      dplyr::select(
        factor,
        strongest_response_effect,
        minimum_response_p,
        strongest_response_comparison,
        best_response_feature_label,
        best_response_feature_p
      ),
    by = "factor"
  ) %>%
  dplyr::mutate(
    x = 5.35 * cos(factor_angle),
    y = 4.10 * sin(factor_angle),
    factor_label = paste0(
      factor,
      ifelse(
        is.finite(minimum_response_p),
        paste0(
          "\npCR - non-pCR = ", sprintf("%+.2f", strongest_response_effect),
          "\nWilcoxon ",
          ifelse(
            !is.finite(minimum_response_p),
            "p = NA",
            ifelse(minimum_response_p < 0.001, "p < 0.001", paste0("p = ", mofa_v50_format_p(minimum_response_p)))
          )
        ),
        ""
      )
    )
  )

mofa_v50_network_plot_edges <- mofa_v50_network_edges %>%
  dplyr::left_join(
    mofa_v50_network_factor_nodes %>%
      dplyr::select(factor, x_factor = x, y_factor = y),
    by = "factor"
  ) %>%
  dplyr::left_join(
    mofa_v50_network_feature_nodes %>%
      dplyr::select(feature_node_id, x_feature = x, y_feature = y),
    by = "feature_node_id"
  )

p_mofa_v50_multifactor_feature_network <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = mofa_v50_network_modules,
    ggplot2::aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = view_label
    ),
    alpha = 0.12,
    color = "grey72",
    linewidth = 0.28
  ) +
  ggplot2::geom_curve(
    data = mofa_v50_network_plot_edges,
    ggplot2::aes(
      x = x_factor,
      y = y_factor,
      xend = x_feature,
      yend = y_feature,
      color = loading_sign,
      linewidth = loading_strength
    ),
    curvature = 0.05,
    alpha = 0.43,
    lineend = "round"
  ) +
  ggplot2::geom_label(
    data = mofa_v50_network_factor_nodes,
    ggplot2::aes(x = x, y = y, label = factor_label),
    fill = "#FFF2B3",
    color = "grey10",
    label.size = 0.28,
    size = 2.35,
    lineheight = 0.88,
    fontface = "bold"
  ) +
  ggplot2::geom_point(
    data = mofa_v50_network_feature_nodes,
    ggplot2::aes(x = x, y = y, fill = view_label, size = maximum_loading),
    shape = 21,
    color = "grey20",
    stroke = 0.42
  ) +
  ggplot2::geom_point(
    data = mofa_v50_network_feature_nodes %>% dplyr::filter(shared_feature),
    ggplot2::aes(x = x, y = y),
    shape = 21,
    size = 3.8,
    fill = NA,
    color = "#7A3E9D",
    stroke = 0.72
  ) +
  ggplot2::geom_text(
    data = mofa_v50_network_feature_nodes %>% dplyr::filter(response_feature),
    ggplot2::aes(x = x, y = y, label = "*"),
    nudge_y = 0.24,
    size = 3.2,
    fontface = "bold",
    color = "black"
  ) +
  ggplot2::geom_text(
    data = mofa_v50_network_feature_nodes,
    ggplot2::aes(
      x = label_x,
      y = y,
      label = feature_label_plotmath,
      hjust = label_hjust
    ),
    parse = TRUE,
    size = 2.10
  ) +
  ggplot2::geom_text(
    data = mofa_v50_network_modules,
    ggplot2::aes(x = module_x, y = label_y, label = module_label),
    fontface = "bold",
    size = 2.35,
    color = "grey20"
  ) +
  ggplot2::scale_color_manual(
    values = c(Positive = "#D55E00", Negative = "#0072B2"),
    name = "Feature-weight sign"
  ) +
  ggplot2::scale_fill_manual(values = view_colors, name = "Omics view") +
  ggplot2::scale_linewidth_continuous(range = c(0.28, 1.00), guide = "none") +
  ggplot2::scale_size_continuous(range = c(2.2, 3.7), guide = "none") +
  ggplot2::coord_equal(
    xlim = c(-7.2, 7.2),
    ylim = c(-5.4, 5.4),
    clip = "off"
  ) +
  ggplot2::labs(
    title = "Multi-factor feature-weight network",
    subtitle = paste0(
      "Five QC-passing factors are ranked by factor-level and active-view feature-level response evidence; ",
      "multi-view factors are preferred and no feature is shown from a view with R² < 2%."
    ),
    caption = paste0(
      "Modules are Ward.D2 clusters of signed feature-weight profiles within each omics view, not predefined pathways. ",
      "Purple outlines denote features connected to multiple factors; asterisks mark nominal p < 0.05 pCR/non-pCR feature contrasts."
    ),
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_void(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.6, color = "grey35"),
    plot.caption = ggplot2::element_text(size = 7.1, color = "grey35", hjust = 0),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(5, 170, 5, 170)
  )

if (FALSE) ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_multifactor_feature_loading_network.svg"),
  plot = p_mofa_v50_multifactor_feature_network,
  width = 12.5,
  height = 9.2,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)


#-----------------------------------------------------------------#
# 22.7 Integrated factor interpretation atlas
#-----------------------------------------------------------------#

# The main atlas emphasizes response enrichment, while retaining the omics
# variance fingerprint and a compact non-parametric association panel.
p_mofa_v50_factor_interpretation_atlas <-
  (
    p_mofa_v50_response_enrichment |
      p_mofa_v50_variance_heatmap
  ) /
  p_mofa_v50_association_pmap +
  patchwork::plot_layout(
    widths = c(1.75, 1.0),
    heights = c(1.65, 1.0),
    guides = "collect"
  ) +
  patchwork::plot_annotation(
    title = "MOFA factor interpretation atlas",
    subtitle = paste0(
      "Response enrichment is the primary clinical summary. Variance explained defines the active omics context; ",
      "the lower panel provides complementary Wilcoxon contrasts."
    )
  ) &
  ggplot2::theme(legend.position = "right")

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_variance_association_atlas.svg"),
  plot = p_mofa_v50_factor_interpretation_atlas,
  width = 12.0,
  height = max(8.4, 3.3 + 0.45 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

#-----------------------------------------------------------------#
# 22.8 Factor7-anchored exhaustive two-factor response scan
#-----------------------------------------------------------------#

mofa_v50_anchor_factor <- if ("Factor7" %in% mofa_v50_factor_order) {
  "Factor7"
} else {
  mofa_factor_wilcoxon_tests %>%
    dplyr::filter(
      comparison == "overall_response_subject_mean",
      factor %in% mofa_v50_factor_order,
      is.finite(wilcoxon_p)
    ) %>%
    dplyr::arrange(wilcoxon_p, dplyr::desc(abs(wilcoxon_effect))) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::pull(factor)
}

mofa_v50_anchor_pair_diagnostics <- dplyr::bind_rows(
  lapply(
    setdiff(mofa_v50_factor_order, mofa_v50_anchor_factor),
    function(partner_factor) {
      pair_factors <- c(mofa_v50_anchor_factor, partner_factor)
      pair_data <- mofa_v50_subject_mean_factor_data %>%
        dplyr::select(SubjectID, TRG_plot, dplyr::all_of(pair_factors)) %>%
        dplyr::filter(
          is.finite(.data[[pair_factors[1]]]),
          is.finite(.data[[pair_factors[2]]]),
          !is.na(TRG_plot)
        ) %>%
        dplyr::mutate(
          TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR"))
        )

      score_matrix <- scale(as.matrix(pair_data[, pair_factors, drop = FALSE]))
      score_matrix[!is.finite(score_matrix)] <- 0
      distance_object <- stats::dist(score_matrix)

      set.seed(20261101 + as.integer(sub("^Factor", "", partner_factor)))
      permanova_fit <- tryCatch(
        vegan::adonis2(
          distance_object ~ TRG_plot,
          data = pair_data,
          permutations = mofa_pair_permanova_permutations
        ),
        error = function(e) NULL
      )

      dispersion_fit <- tryCatch(
        vegan::betadisper(distance_object, pair_data$TRG_plot),
        error = function(e) NULL
      )
      dispersion_test <- if (is.null(dispersion_fit)) {
        NULL
      } else {
        tryCatch(
          vegan::permutest(
            dispersion_fit,
            permutations = mofa_v50_outlier_sensitivity_permutations
          ),
          error = function(e) NULL
        )
      }

      outlier_rows <- as.integer(
        unlist(
          lapply(
            levels(pair_data$TRG_plot),
            function(group_name) {
              group_rows <- which(pair_data$TRG_plot == group_name)
              group_matrix <- score_matrix[group_rows, , drop = FALSE]
              group_distance <- sqrt(
                rowSums(
                  (group_matrix - rep(colMeans(group_matrix), each = nrow(group_matrix)))^2
                )
              )
              group_rows[which.max(group_distance)]
            }
          )
        )
      )

      sensitivity_rows <- dplyr::bind_rows(
        lapply(
          outlier_rows,
          function(excluded_row) {
            reduced_data <- pair_data[-excluded_row, , drop = FALSE]
            reduced_matrix <- scale(as.matrix(reduced_data[, pair_factors, drop = FALSE]))
            reduced_matrix[!is.finite(reduced_matrix)] <- 0
            set.seed(20261201 + excluded_row + as.integer(sub("^Factor", "", partner_factor)))
            reduced_fit <- tryCatch(
              vegan::adonis2(
                stats::dist(reduced_matrix) ~ TRG_plot,
                data = reduced_data,
                permutations = mofa_v50_outlier_sensitivity_permutations
              ),
              error = function(e) NULL
            )
            data.frame(
              excluded_subject = pair_data$SubjectID[excluded_row],
              permanova_r2 = if (is.null(reduced_fit)) NA_real_ else as.numeric(reduced_fit$R2[1]),
              permanova_p = if (is.null(reduced_fit)) NA_real_ else as.numeric(reduced_fit$`Pr(>F)`[1]),
              stringsAsFactors = FALSE
            )
          }
        )
      )

      data.frame(
        anchor_factor = mofa_v50_anchor_factor,
        partner_factor = partner_factor,
        n_subjects = nrow(pair_data),
        permanova_f = if (is.null(permanova_fit)) NA_real_ else as.numeric(permanova_fit$F[1]),
        permanova_r2 = if (is.null(permanova_fit)) NA_real_ else as.numeric(permanova_fit$R2[1]),
        permanova_p = if (is.null(permanova_fit)) NA_real_ else as.numeric(permanova_fit$`Pr(>F)`[1]),
        dispersion_p = if (is.null(dispersion_test)) NA_real_ else as.numeric(dispersion_test$tab$`Pr(>F)`[1]),
        outlier_candidates = paste(pair_data$SubjectID[outlier_rows], collapse = ";"),
        sensitivity_max_p = if (any(is.finite(sensitivity_rows$permanova_p))) max(sensitivity_rows$permanova_p, na.rm = TRUE) else NA_real_,
        sensitivity_min_r2 = if (any(is.finite(sensitivity_rows$permanova_r2))) min(sensitivity_rows$permanova_r2, na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  )
) %>%
  dplyr::mutate(
    permanova_fdr = stats::p.adjust(permanova_p, method = "BH"),
    robust_to_candidate_outliers =
      is.finite(sensitivity_max_p) & sensitivity_max_p < 0.10,
    partner_number = as.integer(sub("^Factor", "", partner_factor))
  ) %>%
  dplyr::arrange(permanova_p, dplyr::desc(permanova_r2))

write.csv(
  mofa_v50_anchor_pair_diagnostics,
  file.path(mofa_v50_result_dir, "MOFA_v50_Factor7_pair_scan_PERMANOVA.csv"),
  row.names = FALSE
)

mofa_v50_build_anchor_pair_plot <- function(partner_factor) {
  diagnostic_row <- mofa_v50_anchor_pair_diagnostics %>%
    dplyr::filter(partner_factor == .env$partner_factor)

  pair_data <- mofa_v50_subject_mean_factor_data %>%
    dplyr::transmute(
      SubjectID,
      TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR")),
      anchor_value = .data[[mofa_v50_anchor_factor]],
      partner_value = .data[[partner_factor]]
    ) %>%
    dplyr::filter(
      is.finite(anchor_value),
      is.finite(partner_value),
      !is.na(TRG_plot)
    ) %>%
    dplyr::mutate(
      anchor_z = as.numeric(scale(anchor_value)),
      partner_z = as.numeric(scale(partner_value))
    )

  hull_data <- pair_data %>%
    dplyr::group_by(TRG_plot) %>%
    dplyr::filter(dplyr::n() >= 3) %>%
    dplyr::slice(chull(anchor_z, partner_z)) %>%
    dplyr::ungroup()

  annotation_label <- paste0(
    "PERMANOVA R² = ", sprintf("%.2f", diagnostic_row$permanova_r2),
    "\n", ifelse(diagnostic_row$permanova_p < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", diagnostic_row$permanova_p))),
    "\nDispersion p = ", mofa_v50_format_p(diagnostic_row$dispersion_p),
    "\nOutlier-check max p = ", mofa_v50_format_p(diagnostic_row$sensitivity_max_p)
  )

  ggplot2::ggplot(pair_data, ggplot2::aes(x = anchor_z, y = partner_z)) +
    ggplot2::geom_polygon(
      data = hull_data,
      ggplot2::aes(group = TRG_plot, fill = TRG_plot, color = TRG_plot),
      alpha = 0.10,
      linewidth = 0.42,
      show.legend = FALSE
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey75", linewidth = 0.30) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey75", linewidth = 0.30) +
    ggplot2::geom_point(
      ggplot2::aes(fill = TRG_plot),
      shape = 21,
      size = 1.85,
      color = "grey20",
      stroke = 0.40,
      alpha = 0.90
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = TRG_plot),
      fun = mean,
      geom = "point",
      shape = 21,
      size = 2.55,
      color = "black",
      stroke = 0.58
    ) +
    ggplot2::annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = annotation_label,
      hjust = -0.03,
      vjust = 1.08,
      size = 2.35,
      lineheight = 0.90,
      color = "black"
    ) +
    ggplot2::scale_fill_manual(values = mofa_group_colors[c("pCR", "non_pCR")], name = "Response") +
    ggplot2::scale_color_manual(values = mofa_group_colors[c("pCR", "non_pCR")], guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.08))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.16))) +
    ggplot2::labs(
      title = paste0(mofa_v50_anchor_factor, " + ", partner_factor),
      x = paste0(mofa_v50_anchor_factor, " subject-mean score (z)"),
      y = paste0(partner_factor, " subject-mean score (z)")
    ) +
    ggplot2::theme_classic(base_size = 8.3) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 8.8),
      legend.position = "bottom",
      aspect.ratio = 1
    )
}

mofa_v50_anchor_pair_plots <- lapply(
  setdiff(mofa_v50_factor_order, mofa_v50_anchor_factor),
  mofa_v50_build_anchor_pair_plot
)

p_mofa_v50_anchor_pair_scan <- patchwork::wrap_plots(
  mofa_v50_anchor_pair_plots,
  ncol = 3,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = paste0(mofa_v50_anchor_factor, "-anchored overall-response factor-pair scan"),
    subtitle = paste0(
      "All remaining factors are evaluated on one subject-level mean per subject. ",
      "Dispersion and candidate-outlier sensitivity are shown with each PERMANOVA result."
    )
  ) &
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_Factor7_pair_scan.svg"),
  plot = p_mofa_v50_anchor_pair_scan,
  width = 12.0,
  height = 11.0,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

#-----------------------------------------------------------------#
# 22.9 Exhaustive three-factor PERMANOVA and static 3D projections
#-----------------------------------------------------------------#

mofa_v50_factor_triples <- utils::combn(
  mofa_v50_factor_order,
  3,
  simplify = FALSE
)

mofa_v50_triple_permanova_atlas <- dplyr::bind_rows(
  lapply(
    seq_along(mofa_v50_factor_triples),
    function(triple_index) {
      triple_factors <- mofa_v50_factor_triples[[triple_index]]
      triple_data <- mofa_v50_subject_mean_factor_data %>%
        dplyr::select(SubjectID, TRG_plot, dplyr::all_of(triple_factors)) %>%
        dplyr::filter(
          dplyr::if_all(dplyr::all_of(triple_factors), is.finite),
          !is.na(TRG_plot)
        ) %>%
        dplyr::mutate(
          TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR"))
        )

      score_matrix <- scale(as.matrix(triple_data[, triple_factors, drop = FALSE]))
      score_matrix[!is.finite(score_matrix)] <- 0
      distance_object <- stats::dist(score_matrix)

      set.seed(20261301 + triple_index)
      permanova_fit <- tryCatch(
        vegan::adonis2(
          distance_object ~ TRG_plot,
          data = triple_data,
          permutations = mofa_v50_triple_screen_permutations
        ),
        error = function(e) NULL
      )
      dispersion_fit <- tryCatch(
        vegan::betadisper(distance_object, triple_data$TRG_plot),
        error = function(e) NULL
      )
      dispersion_test <- if (is.null(dispersion_fit)) {
        NULL
      } else {
        tryCatch(
          vegan::permutest(
            dispersion_fit,
            permutations = mofa_v50_triple_screen_permutations
          ),
          error = function(e) NULL
        )
      }

      data.frame(
        triple_id = paste(triple_factors, collapse = "__"),
        factor_x = triple_factors[1],
        factor_y = triple_factors[2],
        factor_z = triple_factors[3],
        n_subjects = nrow(triple_data),
        permanova_f = if (is.null(permanova_fit)) NA_real_ else as.numeric(permanova_fit$F[1]),
        permanova_r2 = if (is.null(permanova_fit)) NA_real_ else as.numeric(permanova_fit$R2[1]),
        permanova_p = if (is.null(permanova_fit)) NA_real_ else as.numeric(permanova_fit$`Pr(>F)`[1]),
        dispersion_p = if (is.null(dispersion_test)) NA_real_ else as.numeric(dispersion_test$tab$`Pr(>F)`[1]),
        stringsAsFactors = FALSE
      )
    }
  )
) %>%
  dplyr::mutate(
    permanova_fdr = stats::p.adjust(permanova_p, method = "BH"),
    dispersion_acceptable = !is.finite(dispersion_p) | dispersion_p >= 0.05
  ) %>%
  dplyr::arrange(
    dplyr::desc(dispersion_acceptable),
    permanova_p,
    dplyr::desc(permanova_r2)
  )

write.csv(
  mofa_v50_triple_permanova_atlas,
  file.path(mofa_v50_result_dir, "MOFA_v50_three_factor_PERMANOVA_atlas.csv"),
  row.names = FALSE
)

mofa_v50_prespecified_triple <- c("Factor2", "Factor4", "Factor7")
mofa_v50_prespecified_triple <- mofa_v50_prespecified_triple[
  mofa_v50_prespecified_triple %in% mofa_v50_factor_order
]

mofa_v50_best_triple <- mofa_v50_triple_permanova_atlas %>%
  dplyr::filter(
    is.finite(permanova_p),
    dispersion_acceptable
  ) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::select(factor_x, factor_y, factor_z) %>%
  unlist(use.names = FALSE) %>%
  as.character()

if (length(mofa_v50_best_triple) != 3) {
  mofa_v50_best_triple <- mofa_v50_triple_permanova_atlas %>%
    dplyr::filter(is.finite(permanova_p)) %>%
    dplyr::arrange(permanova_p, dplyr::desc(permanova_r2)) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::select(factor_x, factor_y, factor_z) %>%
    unlist(use.names = FALSE) %>%
    as.character()
}

mofa_v50_selected_triples <- unique(
  list(
    if (length(mofa_v50_prespecified_triple) == 3) mofa_v50_prespecified_triple else mofa_v50_best_triple,
    mofa_v50_best_triple
  )
)

mofa_v50_selected_triple_diagnostics <- dplyr::bind_rows(
  lapply(
    seq_along(mofa_v50_selected_triples),
    function(selected_index) {
      triple_factors <- mofa_v50_selected_triples[[selected_index]]
      triple_data <- mofa_v50_subject_mean_factor_data %>%
        dplyr::select(SubjectID, TRG_plot, dplyr::all_of(triple_factors)) %>%
        dplyr::filter(
          dplyr::if_all(dplyr::all_of(triple_factors), is.finite),
          !is.na(TRG_plot)
        ) %>%
        dplyr::mutate(
          TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR"))
        )

      score_matrix <- scale(as.matrix(triple_data[, triple_factors, drop = FALSE]))
      score_matrix[!is.finite(score_matrix)] <- 0
      distance_object <- stats::dist(score_matrix)

      set.seed(20261501 + selected_index)
      final_fit <- tryCatch(
        vegan::adonis2(
          distance_object ~ TRG_plot,
          data = triple_data,
          permutations = mofa_pair_permanova_permutations
        ),
        error = function(e) NULL
      )
      dispersion_fit <- tryCatch(
        vegan::betadisper(distance_object, triple_data$TRG_plot),
        error = function(e) NULL
      )
      dispersion_test <- if (is.null(dispersion_fit)) {
        NULL
      } else {
        tryCatch(
          vegan::permutest(
            dispersion_fit,
            permutations = mofa_pair_permanova_permutations
          ),
          error = function(e) NULL
        )
      }

      outlier_rows <- as.integer(
        unlist(
          lapply(
            levels(triple_data$TRG_plot),
            function(group_name) {
              group_rows <- which(triple_data$TRG_plot == group_name)
              group_matrix <- score_matrix[group_rows, , drop = FALSE]
              group_distance <- sqrt(
                rowSums(
                  (group_matrix - rep(colMeans(group_matrix), each = nrow(group_matrix)))^2
                )
              )
              group_rows[which.max(group_distance)]
            }
          )
        )
      )

      sensitivity_rows <- dplyr::bind_rows(
        lapply(
          outlier_rows,
          function(excluded_row) {
            reduced_data <- triple_data[-excluded_row, , drop = FALSE]
            reduced_matrix <- scale(as.matrix(reduced_data[, triple_factors, drop = FALSE]))
            reduced_matrix[!is.finite(reduced_matrix)] <- 0
            set.seed(20261601 + selected_index + excluded_row)
            reduced_fit <- tryCatch(
              vegan::adonis2(
                stats::dist(reduced_matrix) ~ TRG_plot,
                data = reduced_data,
                permutations = mofa_v50_outlier_sensitivity_permutations
              ),
              error = function(e) NULL
            )
            data.frame(
              excluded_subject = triple_data$SubjectID[excluded_row],
              permanova_r2 = if (is.null(reduced_fit)) NA_real_ else as.numeric(reduced_fit$R2[1]),
              permanova_p = if (is.null(reduced_fit)) NA_real_ else as.numeric(reduced_fit$`Pr(>F)`[1]),
              stringsAsFactors = FALSE
            )
          }
        )
      )

      data.frame(
        selection = if (
          length(mofa_v50_prespecified_triple) == 3 &&
          identical(triple_factors, mofa_v50_prespecified_triple)
        ) {
          "Prespecified Factor2-Factor4-Factor7"
        } else {
          "Best dispersion-acceptable triple"
        },
        triple_id = paste(triple_factors, collapse = "__"),
        factor_x = triple_factors[1],
        factor_y = triple_factors[2],
        factor_z = triple_factors[3],
        n_subjects = nrow(triple_data),
        permanova_f = if (is.null(final_fit)) NA_real_ else as.numeric(final_fit$F[1]),
        permanova_r2 = if (is.null(final_fit)) NA_real_ else as.numeric(final_fit$R2[1]),
        permanova_p = if (is.null(final_fit)) NA_real_ else as.numeric(final_fit$`Pr(>F)`[1]),
        dispersion_p = if (is.null(dispersion_test)) NA_real_ else as.numeric(dispersion_test$tab$`Pr(>F)`[1]),
        outlier_candidates = paste(triple_data$SubjectID[outlier_rows], collapse = ";"),
        sensitivity_max_p = if (any(is.finite(sensitivity_rows$permanova_p))) max(sensitivity_rows$permanova_p, na.rm = TRUE) else NA_real_,
        sensitivity_min_r2 = if (any(is.finite(sensitivity_rows$permanova_r2))) min(sensitivity_rows$permanova_r2, na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  mofa_v50_selected_triple_diagnostics,
  file.path(mofa_v50_result_dir, "MOFA_v50_selected_three_factor_diagnostics.csv"),
  row.names = FALSE
)

mofa_v50_build_3d_projection <- function(triple_factors, diagnostic_row) {
  plot_data <- mofa_v50_subject_mean_factor_data %>%
    dplyr::select(SubjectID, TRG_plot, dplyr::all_of(triple_factors)) %>%
    dplyr::filter(
      dplyr::if_all(dplyr::all_of(triple_factors), is.finite),
      !is.na(TRG_plot)
    ) %>%
    dplyr::mutate(
      TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR"))
    )

  score_matrix <- scale(as.matrix(plot_data[, triple_factors, drop = FALSE]))
  score_matrix[!is.finite(score_matrix)] <- 0
  azimuth <- 40 * pi / 180
  elevation <- 24 * pi / 180

  plot_data$x3 <- score_matrix[, 1]
  plot_data$y3 <- score_matrix[, 2]
  plot_data$z3 <- score_matrix[, 3]
  plot_data$x_projected <- cos(azimuth) * plot_data$x3 - sin(azimuth) * plot_data$y3
  plot_data$y_projected <-
    sin(elevation) * (sin(azimuth) * plot_data$x3 + cos(azimuth) * plot_data$y3) +
    cos(elevation) * plot_data$z3

  axis_length <- max(abs(score_matrix), na.rm = TRUE) * 0.78
  axis_matrix <- rbind(
    c(0, 0, 0),
    c(axis_length, 0, 0),
    c(0, axis_length, 0),
    c(0, 0, axis_length)
  )
  axis_data <- data.frame(
    axis = c("origin", triple_factors),
    x_projected = cos(azimuth) * axis_matrix[, 1] - sin(azimuth) * axis_matrix[, 2],
    y_projected =
      sin(elevation) * (sin(azimuth) * axis_matrix[, 1] + cos(azimuth) * axis_matrix[, 2]) +
      cos(elevation) * axis_matrix[, 3],
    stringsAsFactors = FALSE
  )

  hull_data <- plot_data %>%
    dplyr::group_by(TRG_plot) %>%
    dplyr::filter(dplyr::n() >= 3) %>%
    dplyr::slice(chull(x_projected, y_projected)) %>%
    dplyr::ungroup()

  annotation_label <- paste0(
    "3D PERMANOVA R² = ", sprintf("%.2f", diagnostic_row$permanova_r2),
    "\n", ifelse(diagnostic_row$permanova_p < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", diagnostic_row$permanova_p))),
    "\nDispersion p = ", mofa_v50_format_p(diagnostic_row$dispersion_p),
    "\nOutlier-check max p = ", mofa_v50_format_p(diagnostic_row$sensitivity_max_p)
  )

  ggplot2::ggplot(plot_data, ggplot2::aes(x = x_projected, y = y_projected)) +
    ggplot2::geom_polygon(
      data = hull_data,
      ggplot2::aes(group = TRG_plot, fill = TRG_plot, color = TRG_plot),
      alpha = 0.10,
      linewidth = 0.48,
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = axis_data[-1, ],
      ggplot2::aes(
        x = axis_data$x_projected[1],
        y = axis_data$y_projected[1],
        xend = x_projected,
        yend = y_projected
      ),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = grid::unit(0.10, "inches")),
      linewidth = 0.45,
      color = "grey35"
    ) +
    ggplot2::geom_text(
      data = axis_data[-1, ],
      ggplot2::aes(x = x_projected, y = y_projected, label = axis),
      inherit.aes = FALSE,
      nudge_y = 0.10,
      size = 3.0,
      fontface = "bold"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = TRG_plot),
      shape = 21,
      size = 2.5,
      color = "grey20",
      stroke = 0.48,
      alpha = 0.92
    ) +
    ggplot2::stat_summary(
      ggplot2::aes(fill = TRG_plot),
      fun = mean,
      geom = "point",
      shape = 21,
      size = 3.4,
      color = "black",
      stroke = 0.66
    ) +
    ggplot2::annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = annotation_label,
      hjust = -0.03,
      vjust = 1.05,
      size = 2.8,
      lineheight = 0.92,
      color = "black"
    ) +
    ggplot2::scale_fill_manual(values = mofa_group_colors[c("pCR", "non_pCR")], name = "Response") +
    ggplot2::scale_color_manual(values = mofa_group_colors[c("pCR", "non_pCR")], guide = "none") +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(
      title = diagnostic_row$selection,
      subtitle = paste(triple_factors, collapse = " + "),
      x = NULL,
      y = NULL,
      caption = "PERMANOVA is calculated in the original standardized three-factor space; the panel is a fixed isometric projection."
    ) +
    ggplot2::theme_classic(base_size = 9.0) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(face = "bold", size = 8.5),
      plot.caption = ggplot2::element_text(size = 7.0, color = "grey35", hjust = 0),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(5, 12, 5, 5)
    )
}

mofa_v50_3d_plots <- lapply(
  seq_len(nrow(mofa_v50_selected_triple_diagnostics)),
  function(selected_index) {
    diagnostic_row <- mofa_v50_selected_triple_diagnostics[selected_index, , drop = FALSE]
    mofa_v50_build_3d_projection(
      c(diagnostic_row$factor_x, diagnostic_row$factor_y, diagnostic_row$factor_z),
      diagnostic_row
    )
  }
)

p_mofa_v50_three_factor_response_maps <- patchwork::wrap_plots(
  mofa_v50_3d_plots,
  ncol = length(mofa_v50_3d_plots),
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "Three-factor response separation in overall subject-mean MOFA space",
    subtitle = "The prespecified Factor2-Factor4-Factor7 combination is compared with the best dispersion-acceptable triple from all combinations."
  ) &
  ggplot2::theme(legend.position = "bottom")

if (FALSE) ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_three_factor_response_3D_maps.svg"),
  plot = p_mofa_v50_three_factor_response_maps,
  width = 6.0 * length(mofa_v50_3d_plots),
  height = 5.7,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)


#-----------------------------------------------------------------#
# 22.7b Atlas refinement and 3D/network updates for v50
#-----------------------------------------------------------------#

mofa_v50_coarse_neglog10_breaks <- function(upper_bound) {
  breaks_out <- c(1, 2)
  breaks_out <- breaks_out[breaks_out <= upper_bound + 1e-8]
  if (length(breaks_out) == 0) {
    breaks_out <- upper_bound
  }
  unique(breaks_out)
}

# Refined variance heatmap for the atlas
p_mofa_v50_variance_heatmap_atlas <- ggplot2::ggplot(
  mofa_v50_variance_plot_data,
  ggplot2::aes(x = view_label, y = factor, fill = r2)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.32) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(r2 >= 0.005, sprintf("%.1f", 100 * r2), "")),
    size = 3.05
  ) +
  ggplot2::scale_fill_gradient(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "white",
    high = "#2166AC",
    name = "Variance\nexplained"
  ) +
  ggplot2::labs(
    title = "Variance explained by factor and omics view",
    subtitle = "Cell labels are percentages.",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 9.4) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.7, color = "grey35"),
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, size = 8.3),
    axis.text.y = ggplot2::element_text(face = "plain", size = 8.4),
    legend.position = "right"
  )

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_variance_explained_numeric_order.svg"),
  plot = p_mofa_v50_variance_heatmap_atlas,
  width = 4.7,
  height = max(4.4, 1.8 + 0.28 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

p_mofa_v50_response_enrichment_atlas <- ggplot2::ggplot(
  mofa_v50_response_enrichment_data,
  ggplot2::aes(
    x = standardized_effect,
    y = factor,
    color = response_direction,
    fill = response_direction,
    size = minus_log10_p
  )
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.42) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = standardized_effect, yend = factor),
    linewidth = 0.62,
    alpha = 0.72,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(shape = 21, color = "grey20", stroke = 0.55) +
  ggplot2::geom_text(
    ggplot2::aes(label = significance_label),
    nudge_y = 0.25,
    size = 2.9,
    color = "black",
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~ Timepoint_display, nrow = 1) +
  ggplot2::scale_color_manual(
    values = mofa_response_colors[c("pCR", "non_pCR")],
    labels = c(pCR = "pCR-enriched", non_pCR = "non-pCR-enriched"),
    name = "Direction"
  ) +
  ggplot2::scale_fill_manual(
    values = mofa_response_colors[c("pCR", "non_pCR")],
    guide = "none"
  ) +
  ggplot2::scale_size_continuous(
    range = c(2.5, 5.2),
    limits = c(0, mofa_v50_response_enrichment_p_upper),
    breaks = mofa_v50_coarse_neglog10_breaks(mofa_v50_response_enrichment_p_upper),
    name = expression(-log[10](italic(p)))
  ) +
  ggplot2::labs(
    title = "pCR versus non-pCR factor enrichment",
    subtitle = "Hedges' g is oriented as pCR minus non-pCR.",
    x = "Standardized response enrichment (Hedges' g)",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 9.2) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.6, color = "grey35"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "plain"),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    legend.position = "right"
  )

mofa_v50_association_delta_data <- mofa_v50_association_data %>%
  dplyr::filter(
    comparison_label %in% c(
      "Paired\nchange",
      "Differential\npaired change"
    )
  ) %>%
  dplyr::mutate(
    comparison_label = factor(
      comparison_label,
      levels = c("Paired\nchange", "Differential\npaired change")
    )
  )

p_mofa_v50_association_delta_heatmap <- ggplot2::ggplot(
  mofa_v50_association_delta_data,
  ggplot2::aes(x = comparison_label, y = factor, fill = minus_log10_p)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.30) +
  ggplot2::geom_text(ggplot2::aes(label = effect_label), size = 2.35) +
  ggplot2::scale_fill_gradient(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "white",
    high = "black",
    limits = c(0, mofa_v50_association_p_upper),
    breaks = mofa_v50_coarse_neglog10_breaks(mofa_v50_association_p_upper),
    oob = scales::squish,
    name = expression(-log[10](italic(p)))
  ) +
  ggplot2::labs(
    title = "Paired-change associations",
    subtitle = "Text is the signed median contrast.",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 9.2) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.6, color = "grey35"),
    axis.text.x = ggplot2::element_text(size = 8.6, face = "plain"),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    legend.position = "right"
  )

p_mofa_v50_factor_interpretation_atlas <-
  (p_mofa_v50_variance_heatmap_atlas | p_mofa_v50_response_enrichment_atlas | p_mofa_v50_association_delta_heatmap) +
  patchwork::plot_layout(widths = c(0.92, 1.05, 0.62), guides = "collect") +
  patchwork::plot_annotation(
    title = "MOFA factor interpretation atlas",
    subtitle = paste0(
      "Variance explained provides the omics context, response enrichment is the primary clinical summary, ",
      "and the right heatmap highlights paired-change associations."
    )
  ) &
  ggplot2::theme(legend.position = "right")

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_variance_association_atlas.svg"),
  plot = p_mofa_v50_factor_interpretation_atlas,
  width = 10.8,
  height = max(6.2, 3.1 + 0.36 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

# Restrict three-factor screening to Factor7 plus two factors from 4/8/2/9/10.
mofa_v50_anchor_factor <- if ("Factor7" %in% mofa_v50_factor_order) {
  "Factor7"
} else {
  mofa_v50_anchor_factor
}

mofa_v50_triple_candidate_pool <- intersect(
  c("Factor4", "Factor8", "Factor2", "Factor9", "Factor10"),
  setdiff(mofa_v50_factor_order, mofa_v50_anchor_factor)
)

mofa_v50_factor_triples <- if (length(mofa_v50_triple_candidate_pool) >= 2) {
  lapply(
    utils::combn(mofa_v50_triple_candidate_pool, 2, simplify = FALSE),
    function(partners) c(mofa_v50_anchor_factor, partners)
  )
} else {
  list(c(mofa_v50_factor_order[1:3]))
}

mofa_v50_triple_permanova_atlas <- dplyr::bind_rows(
  lapply(
    seq_along(mofa_v50_factor_triples),
    function(triple_index) {
      triple_factors <- mofa_v50_factor_triples[[triple_index]]
      triple_data <- mofa_v50_subject_mean_factor_data %>%
        dplyr::select(SubjectID, TRG_plot, dplyr::all_of(triple_factors)) %>%
        dplyr::filter(
          dplyr::if_all(dplyr::all_of(triple_factors), is.finite),
          !is.na(TRG_plot)
        ) %>%
        dplyr::mutate(
          TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR"))
        )

      score_matrix <- scale(as.matrix(triple_data[, triple_factors, drop = FALSE]))
      score_matrix[!is.finite(score_matrix)] <- 0
      distance_object <- stats::dist(score_matrix)

      set.seed(20261701 + triple_index)
      permanova_fit <- tryCatch(
        vegan::adonis2(
          distance_object ~ TRG_plot,
          data = triple_data,
          permutations = mofa_v50_outlier_sensitivity_permutations
        ),
        error = function(e) NULL
      )
      dispersion_fit <- tryCatch(
        vegan::betadisper(distance_object, triple_data$TRG_plot),
        error = function(e) NULL
      )
      dispersion_test <- if (is.null(dispersion_fit)) {
        NULL
      } else {
        tryCatch(
          vegan::permutest(
            dispersion_fit,
            permutations = mofa_v50_outlier_sensitivity_permutations
          ),
          error = function(e) NULL
        )
      }

      outlier_rows <- as.integer(
        unlist(
          lapply(
            levels(triple_data$TRG_plot),
            function(group_name) {
              group_rows <- which(triple_data$TRG_plot == group_name)
              group_matrix <- score_matrix[group_rows, , drop = FALSE]
              group_distance <- sqrt(
                rowSums((group_matrix - rep(colMeans(group_matrix), each = nrow(group_matrix)))^2)
              )
              group_rows[which.max(group_distance)]
            }
          )
        )
      )

      sensitivity_rows <- dplyr::bind_rows(
        lapply(
          outlier_rows,
          function(excluded_row) {
            reduced_data <- triple_data[-excluded_row, , drop = FALSE]
            reduced_matrix <- scale(as.matrix(reduced_data[, triple_factors, drop = FALSE]))
            reduced_matrix[!is.finite(reduced_matrix)] <- 0
            set.seed(20261801 + triple_index + excluded_row)
            reduced_fit <- tryCatch(
              vegan::adonis2(
                stats::dist(reduced_matrix) ~ TRG_plot,
                data = reduced_data,
                permutations = mofa_v50_outlier_sensitivity_permutations
              ),
              error = function(e) NULL
            )
            data.frame(
              excluded_subject = triple_data$SubjectID[excluded_row],
              permanova_r2 = if (is.null(reduced_fit)) NA_real_ else as.numeric(reduced_fit$R2[1]),
              permanova_p = if (is.null(reduced_fit)) NA_real_ else as.numeric(reduced_fit$`Pr(>F)`[1]),
              stringsAsFactors = FALSE
            )
          }
        )
      )

      data.frame(
        triple_id = paste(triple_factors, collapse = "__"),
        factor_x = triple_factors[1],
        factor_y = triple_factors[2],
        factor_z = triple_factors[3],
        n_subjects = nrow(triple_data),
        permanova_f = if (is.null(permanova_fit)) NA_real_ else as.numeric(permanova_fit$F[1]),
        permanova_r2 = if (is.null(permanova_fit)) NA_real_ else as.numeric(permanova_fit$R2[1]),
        permanova_p = if (is.null(permanova_fit)) NA_real_ else as.numeric(permanova_fit$`Pr(>F)`[1]),
        dispersion_p = if (is.null(dispersion_test)) NA_real_ else as.numeric(dispersion_test$tab$`Pr(>F)`[1]),
        outlier_candidates = paste(triple_data$SubjectID[outlier_rows], collapse = ";"),
        sensitivity_max_p = if (any(is.finite(sensitivity_rows$permanova_p))) max(sensitivity_rows$permanova_p, na.rm = TRUE) else NA_real_,
        sensitivity_min_r2 = if (any(is.finite(sensitivity_rows$permanova_r2))) min(sensitivity_rows$permanova_r2, na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  )
) %>%
  dplyr::mutate(
    permanova_fdr = stats::p.adjust(permanova_p, method = "BH"),
    dispersion_ok = is.na(dispersion_p) | dispersion_p >= 0.05,
    outlier_robust = is.na(sensitivity_max_p) | sensitivity_max_p < 0.10
  ) %>%
  dplyr::arrange(dplyr::desc(dispersion_ok), permanova_p, dplyr::desc(permanova_r2), dplyr::desc(outlier_robust))

write.csv(
  mofa_v50_triple_permanova_atlas,
  file.path(mofa_v50_result_dir, "MOFA_v50_three_factor_PERMANOVA_atlas.csv"),
  row.names = FALSE
)

mofa_v50_selected_triple_diagnostics <- mofa_v50_triple_permanova_atlas %>%
  dplyr::slice_head(n = 1)

write.csv(
  mofa_v50_selected_triple_diagnostics,
  file.path(mofa_v50_result_dir, "MOFA_v50_selected_three_factor_diagnostics.csv"),
  row.names = FALSE
)

mofa_v50_rotate_project_3d <- function(x, y, z, azimuth_deg = 46, elevation_deg = 24, distance = 8.5) {
  azimuth <- azimuth_deg * pi / 180
  elevation <- elevation_deg * pi / 180
  x1 <- cos(azimuth) * x - sin(azimuth) * y
  y1 <- sin(azimuth) * x + cos(azimuth) * y
  z1 <- z
  x2 <- x1
  y2 <- cos(elevation) * y1 - sin(elevation) * z1
  z2 <- sin(elevation) * y1 + cos(elevation) * z1
  scale_factor <- distance / (distance - z2 + 1e-06)
  data.frame(
    x = x2 * scale_factor,
    y = y2 * scale_factor,
    depth = z2,
    scale_factor = scale_factor,
    stringsAsFactors = FALSE
  )
}

mofa_v50_box_edges <- function(x_range, y_range, z_range) {
  corners <- expand.grid(
    x = x_range,
    y = y_range,
    z = z_range,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  edges <- list(
    c(1, 2), c(1, 3), c(1, 5),
    c(2, 4), c(2, 6), c(3, 4),
    c(3, 7), c(4, 8), c(5, 6),
    c(5, 7), c(6, 8), c(7, 8)
  )
  dplyr::bind_rows(
    lapply(
      seq_along(edges),
      function(edge_index) {
        edge_rows <- edges[[edge_index]]
        data.frame(
          edge_id = edge_index,
          x = corners$x[edge_rows],
          y = corners$y[edge_rows],
          z = corners$z[edge_rows],
          point_order = c(1, 2),
          stringsAsFactors = FALSE
        )
      }
    )
  )
}

mofa_v50_make_ellipsoid_paths <- function(score_matrix, group_vector, level = 0.68) {
  circle_parameter <- seq(0, 2 * pi, length.out = 120)
  basis_list <- list(
    cbind(cos(circle_parameter), sin(circle_parameter), 0),
    cbind(cos(circle_parameter), 0, sin(circle_parameter)),
    cbind(0, cos(circle_parameter), sin(circle_parameter))
  )
  dplyr::bind_rows(
    lapply(
      levels(group_vector),
      function(group_name) {
        group_matrix <- score_matrix[group_vector == group_name, , drop = FALSE]
        if (nrow(group_matrix) < 4) {
          return(NULL)
        }
        center <- colMeans(group_matrix)
        covariance <- stats::cov(group_matrix)
        eig <- eigen(covariance, symmetric = TRUE)
        radii <- sqrt(pmax(eig$values, 1e-06)) * sqrt(stats::qchisq(level, df = 3))
        transform_matrix <- eig$vectors %*% diag(radii, nrow = 3)
        dplyr::bind_rows(
          lapply(
            seq_along(basis_list),
            function(loop_index) {
              loop_xyz <- t(center + t(transform_matrix %*% t(basis_list[[loop_index]])))
              data.frame(
                TRG_plot = group_name,
                loop = paste(group_name, loop_index, sep = "__"),
                x3 = loop_xyz[, 1],
                y3 = loop_xyz[, 2],
                z3 = loop_xyz[, 3],
                stringsAsFactors = FALSE
              )
            }
          )
        )
      }
    )
  )
}

mofa_v50_build_3d_scatter <- function(triple_factors, diagnostic_row) {
  plot_data <- mofa_v50_subject_mean_factor_data %>%
    dplyr::select(SubjectID, TRG_plot, dplyr::all_of(triple_factors)) %>%
    dplyr::filter(
      dplyr::if_all(dplyr::all_of(triple_factors), is.finite),
      !is.na(TRG_plot)
    ) %>%
    dplyr::mutate(TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR")))

  score_matrix <- scale(as.matrix(plot_data[, triple_factors, drop = FALSE]))
  score_matrix[!is.finite(score_matrix)] <- 0
  plot_data$x3 <- score_matrix[, 1]
  plot_data$y3 <- score_matrix[, 2]
  plot_data$z3 <- score_matrix[, 3]

  span <- apply(score_matrix, 2, range)
  pad <- apply(score_matrix, 2, function(v) 0.18 * diff(range(v)))
  x_range <- c(span[1, 1] - pad[1], span[2, 1] + pad[1])
  y_range <- c(span[1, 2] - pad[2], span[2, 2] + pad[2])
  z_range <- c(span[1, 3] - pad[3], span[2, 3] + pad[3])
  floor_z <- z_range[1]

  projected_points <- mofa_v50_rotate_project_3d(plot_data$x3, plot_data$y3, plot_data$z3)
  projected_shadow <- mofa_v50_rotate_project_3d(plot_data$x3, plot_data$y3, rep(floor_z, nrow(plot_data)))
  plot_data$x_proj <- projected_points$x
  plot_data$y_proj <- projected_points$y
  plot_data$shadow_x <- projected_shadow$x
  plot_data$shadow_y <- projected_shadow$y
  finite_depth <- is.finite(projected_points$depth)
  if (sum(finite_depth) >= 2 && diff(range(projected_points$depth[finite_depth])) > 0) {
    plot_data$depth_scale <- scales::rescale(
      projected_points$depth,
      to = c(2.0, 3.2),
      from = range(projected_points$depth[finite_depth])
    )
  } else {
    plot_data$depth_scale <- rep(2.6, nrow(plot_data))
  }
  plot_data$depth_scale[!is.finite(plot_data$depth_scale)] <- 2.6
  plot_data <- plot_data %>%
    dplyr::filter(
      is.finite(x_proj),
      is.finite(y_proj),
      is.finite(shadow_x),
      is.finite(shadow_y),
      is.finite(depth_scale)
    )

  box_data <- mofa_v50_box_edges(x_range, y_range, z_range)
  box_projected <- mofa_v50_rotate_project_3d(box_data$x, box_data$y, box_data$z)
  box_data$x_proj <- box_projected$x
  box_data$y_proj <- box_projected$y

  axis_data <- data.frame(
    axis = triple_factors,
    x3 = c(x_range[2], x_range[1], x_range[1]),
    y3 = c(y_range[1], y_range[2], y_range[1]),
    z3 = c(floor_z, floor_z, z_range[2]),
    stringsAsFactors = FALSE
  )
  axis_projected <- mofa_v50_rotate_project_3d(axis_data$x3, axis_data$y3, axis_data$z3)
  axis_data$x_proj <- axis_projected$x
  axis_data$y_proj <- axis_projected$y

  mean_points <- plot_data %>%
    dplyr::group_by(TRG_plot) %>%
    dplyr::summarise(
      x3 = mean(x3),
      y3 = mean(y3),
      z3 = mean(z3),
      .groups = "drop"
    )
  mean_projected <- mofa_v50_rotate_project_3d(mean_points$x3, mean_points$y3, mean_points$z3)
  mean_points$x_proj <- mean_projected$x
  mean_points$y_proj <- mean_projected$y

  ellipsoid_paths <- mofa_v50_make_ellipsoid_paths(score_matrix, plot_data$TRG_plot)
  if (nrow(ellipsoid_paths) > 0) {
    ellipsoid_projected <- mofa_v50_rotate_project_3d(ellipsoid_paths$x3, ellipsoid_paths$y3, ellipsoid_paths$z3)
    ellipsoid_paths$x_proj <- ellipsoid_projected$x
    ellipsoid_paths$y_proj <- ellipsoid_projected$y
  }

  annotation_label <- paste0(
    paste(triple_factors, collapse = " + "),
    "\nPERMANOVA R² = ", sprintf("%.2f", diagnostic_row$permanova_r2),
    "\nPERMANOVA ",
    ifelse(
      is.na(diagnostic_row$permanova_p),
      "p = NA",
      ifelse(
        diagnostic_row$permanova_p < 0.001,
        "p < 0.001",
        paste0("p = ", sprintf("%.3f", diagnostic_row$permanova_p))
      )
    ),
    "\nDispersion ",
    ifelse(
      is.na(diagnostic_row$dispersion_p),
      "p = NA",
      paste0("p = ", mofa_v50_format_p(diagnostic_row$dispersion_p))
    ),
    "\nOutlier-check max ",
    ifelse(
      is.na(diagnostic_row$sensitivity_max_p),
      "p = NA",
      paste0("p = ", mofa_v50_format_p(diagnostic_row$sensitivity_max_p))
    )
  )

  x_limits <- range(c(box_data$x_proj, plot_data$x_proj, plot_data$shadow_x), na.rm = TRUE)
  y_limits <- range(c(box_data$y_proj, plot_data$y_proj, plot_data$shadow_y), na.rm = TRUE)
  x_margin <- 0.08 * diff(x_limits)
  y_margin <- 0.12 * diff(y_limits)

  ggplot2::ggplot() +
    ggplot2::geom_path(
      data = box_data,
      ggplot2::aes(x = x_proj, y = y_proj, group = edge_id),
      color = "grey74",
      linewidth = 0.38
    ) +
    ggplot2::geom_segment(
      data = plot_data,
      ggplot2::aes(x = shadow_x, y = shadow_y, xend = x_proj, yend = y_proj, color = TRG_plot),
      linewidth = 0.34,
      alpha = 0.28,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(x = shadow_x, y = shadow_y, fill = TRG_plot),
      shape = 21,
      size = 1.95,
      color = NA,
      alpha = 0.18,
      show.legend = FALSE
    ) +
    ggplot2::geom_path(
      data = ellipsoid_paths,
      ggplot2::aes(x = x_proj, y = y_proj, group = loop, color = TRG_plot),
      linewidth = 0.46,
      alpha = 0.70,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(x = x_proj, y = y_proj, fill = TRG_plot),
      size = plot_data$depth_scale,
      shape = 21,
      color = "grey20",
      stroke = 0.44,
      alpha = 0.93,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = mean_points,
      ggplot2::aes(x = x_proj, y = y_proj, fill = TRG_plot),
      shape = 21,
      size = 3.5,
      color = "black",
      stroke = 0.65
    ) +
    ggplot2::geom_text(
      data = axis_data,
      ggplot2::aes(x = x_proj, y = y_proj, label = axis),
      fontface = "bold",
      size = 3.1,
      nudge_y = 0.12
    ) +
    ggplot2::annotate(
      "text",
      x = x_limits[1],
      y = y_limits[2],
      label = annotation_label,
      parse = FALSE,
      hjust = 0,
      vjust = 1,
      lineheight = 0.92,
      size = 2.75,
      color = "black"
    ) +
    ggplot2::scale_fill_manual(values = mofa_response_colors[c("pCR", "non_pCR")], name = "Response") +
    ggplot2::scale_color_manual(values = mofa_response_colors[c("pCR", "non_pCR")], guide = "none") +
    ggplot2::coord_equal(
      xlim = c(x_limits[1] - x_margin, x_limits[2] + x_margin),
      ylim = c(y_limits[1] - y_margin, y_limits[2] + y_margin),
      clip = "off"
    ) +
    ggplot2::labs(
      title = "Overall-response 3D MOFA separation",
      subtitle = "Perspective scatter with floor shadows and group covariance wireframes.",
      x = NULL,
      y = NULL,
      caption = "Factor7 is fixed; the other two factors are chosen from Factor4, Factor8, Factor2, Factor9, and Factor10."
    ) +
    ggplot2::theme_void(base_size = 9.2) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 7.8, color = "grey35"),
      plot.caption = ggplot2::element_text(size = 7.1, color = "grey35", hjust = 0),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(5, 8, 5, 5)
    )
}

p_mofa_v50_three_factor_response_maps <- mofa_v50_build_3d_scatter(
  c(
    mofa_v50_selected_triple_diagnostics$factor_x[1],
    mofa_v50_selected_triple_diagnostics$factor_y[1],
    mofa_v50_selected_triple_diagnostics$factor_z[1]
  ),
  mofa_v50_selected_triple_diagnostics[1, , drop = FALSE]
)

if (FALSE) ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_three_factor_response_3D_maps.svg"),
  plot = p_mofa_v50_three_factor_response_maps,
  width = 7.1,
  height = 6.1,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

# Triple-specific network with factors outside and feature modules inside.
mofa_v50_triple_network_factors <- c(
  mofa_v50_selected_triple_diagnostics$factor_x[1],
  mofa_v50_selected_triple_diagnostics$factor_y[1],
  mofa_v50_selected_triple_diagnostics$factor_z[1]
)

mofa_v50_triple_active_factor_views <- mofa_variance_explained %>%
  dplyr::mutate(factor = as.character(factor), view = as.character(view)) %>%
  dplyr::filter(
    factor %in% mofa_v50_triple_network_factors,
    r2 >= mofa_active_view_r2
  ) %>%
  dplyr::select(factor, view, view_r2 = r2)

mofa_v50_triple_network_seed_edges <- dplyr::bind_rows(
  mofa_feature_weights %>%
    dplyr::filter(
      factor %in% mofa_v50_triple_network_factors,
      display_eligible
    ) %>%
    dplyr::inner_join(mofa_v50_triple_active_factor_views, by = c("factor", "view")) %>%
    dplyr::group_by(factor, view, direction) %>%
    dplyr::slice_max(
      order_by = abs(weight_within_view),
      n = mofa_network_features_per_view_direction,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup(),
  mofa_feature_weights %>%
    dplyr::filter(
      factor %in% mofa_v50_triple_network_factors,
      display_eligible,
      is.finite(p_value),
      p_value < 0.05
    ) %>%
    dplyr::inner_join(mofa_v50_triple_active_factor_views, by = c("factor", "view")) %>%
    dplyr::group_by(factor) %>%
    dplyr::arrange(p_value, dplyr::desc(abs(weight_within_view)), .by_group = TRUE) %>%
    dplyr::slice_head(n = 2) %>%
    dplyr::ungroup()
) %>%
  dplyr::distinct(factor, view, feature, .keep_all = TRUE) %>%
  dplyr::mutate(
    feature_node_id = paste(view, feature, sep = "::"),
    edge_id = paste(factor, feature_node_id, sep = "__")
  )

mofa_v50_triple_network_edges <- mofa_feature_weights %>%
  dplyr::filter(
    factor %in% mofa_v50_triple_network_factors,
    display_eligible
  ) %>%
  dplyr::inner_join(mofa_v50_triple_active_factor_views, by = c("factor", "view")) %>%
  dplyr::mutate(
    feature_node_id = paste(view, feature, sep = "::"),
    edge_id = paste(factor, feature_node_id, sep = "__")
  ) %>%
  dplyr::filter(
    feature_node_id %in% unique(mofa_v50_triple_network_seed_edges$feature_node_id),
    abs(weight_within_view) >= mofa_network_shared_loading_threshold |
      edge_id %in% mofa_v50_triple_network_seed_edges$edge_id
  ) %>%
  dplyr::mutate(
    loading_sign = ifelse(weight >= 0, "Positive", "Negative"),
    loading_strength = abs(weight_within_view)
  )

mofa_v50_triple_module_membership <- dplyr::bind_rows(
  lapply(
    split(mofa_v50_triple_network_edges, mofa_v50_triple_network_edges$view),
    function(view_edge_data) {
      profile <- view_edge_data %>%
        dplyr::select(feature_node_id, factor, weight_within_view) %>%
        tidyr::pivot_wider(
          names_from = factor,
          values_from = weight_within_view,
          values_fill = 0
        )
      for (factor_name in setdiff(mofa_v50_triple_network_factors, colnames(profile))) {
        profile[[factor_name]] <- 0
      }
      profile_matrix <- as.matrix(profile[, mofa_v50_triple_network_factors, drop = FALSE])
      rownames(profile_matrix) <- profile$feature_node_id
      module <- if (nrow(profile_matrix) <= 2) {
        rep(1L, nrow(profile_matrix))
      } else {
        stats::cutree(
          stats::hclust(stats::dist(profile_matrix), method = "ward.D2"),
          k = min(2L, max(1L, ceiling(nrow(profile_matrix) / 7)))
        )
      }
      data.frame(
        feature_node_id = rownames(profile_matrix),
        module = as.integer(module),
        stringsAsFactors = FALSE
      )
    }
  )
)

mofa_v50_triple_network_feature_nodes <- mofa_v50_triple_network_edges %>%
  dplyr::group_by(feature_node_id, view, view_label, feature, feature_label) %>%
  dplyr::summarise(
    n_connected_factors = dplyr::n_distinct(factor),
    maximum_loading = max(loading_strength, na.rm = TRUE),
    minimum_response_p = if (any(is.finite(p_value))) min(p_value, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  dplyr::left_join(mofa_v50_triple_module_membership, by = "feature_node_id") %>%
  dplyr::mutate(
    view_order = match(view, required_views),
    module_label = paste0(as.character(view_label), " M", module),
    shared_feature = n_connected_factors >= 2,
    response_feature = is.finite(minimum_response_p) & minimum_response_p < 0.05,
    feature_label_plotmath = mofa_v50_plotmath_feature_label(view, feature_label)
  ) %>%
  dplyr::arrange(view_order, module, dplyr::desc(n_connected_factors), dplyr::desc(maximum_loading))

mofa_v50_triple_network_modules <- mofa_v50_triple_network_feature_nodes %>%
  dplyr::distinct(view, view_label, view_order, module, module_label) %>%
  dplyr::arrange(view_order, module) %>%
  dplyr::mutate(
    module_index = dplyr::row_number(),
    module_angle = pi / 2 - 2 * pi * (module_index - 1) / dplyr::n(),
    module_x = 2.25 * cos(module_angle),
    module_y = 1.90 * sin(module_angle)
  )

mofa_v50_triple_network_feature_nodes <- mofa_v50_triple_network_feature_nodes %>%
  dplyr::left_join(
    mofa_v50_triple_network_modules %>%
      dplyr::select(module_label, module_x, module_y),
    by = "module_label"
  ) %>%
  dplyr::group_by(module_label) %>%
  dplyr::mutate(
    feature_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.34,
    x = module_x,
    y = module_y + feature_offset,
    label_x = x + ifelse(module_x < -0.25, -0.46, 0.46),
    label_hjust = ifelse(module_x < -0.25, 1, 0)
  ) %>%
  dplyr::ungroup()

mofa_v50_triple_network_modules <- mofa_v50_triple_network_modules %>%
  dplyr::left_join(
    mofa_v50_triple_network_feature_nodes %>%
      dplyr::group_by(module_label) %>%
      dplyr::summarise(
        xmin = min(x) - 0.30,
        xmax = max(x) + 0.30,
        ymin = min(y) - 0.22,
        ymax = max(y) + 0.22,
        label_y = max(y) + 0.38,
        .groups = "drop"
      ),
    by = "module_label"
  )

mofa_v50_triple_network_factor_nodes <- data.frame(
  factor = mofa_v50_triple_network_factors,
  factor_angle = pi / 2 - 2 * pi * (seq_along(mofa_v50_triple_network_factors) - 1) /
    length(mofa_v50_triple_network_factors),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    mofa_factor_wilcoxon_tests %>%
      dplyr::filter(
        comparison == "overall_response_subject_mean",
        factor %in% mofa_v50_triple_network_factors
      ) %>%
      dplyr::select(factor, strongest_response_effect = wilcoxon_effect, minimum_response_p = wilcoxon_p),
    by = "factor"
  ) %>%
  dplyr::mutate(
    x = 5.20 * cos(factor_angle),
    y = 4.10 * sin(factor_angle),
    factor_label = paste0(
      factor,
      ifelse(
        is.finite(minimum_response_p),
        paste0(
          "\npCR - non-pCR = ", sprintf("%+.2f", strongest_response_effect),
          "\nWilcoxon ",
          ifelse(minimum_response_p < 0.001, "p < 0.001", paste0("p = ", mofa_v50_format_p(minimum_response_p)))
        ),
        ""
      )
    )
  )

mofa_v50_triple_network_plot_edges <- mofa_v50_triple_network_edges %>%
  dplyr::left_join(
    mofa_v50_triple_network_factor_nodes %>% dplyr::select(factor, x_factor = x, y_factor = y),
    by = "factor"
  ) %>%
  dplyr::left_join(
    mofa_v50_triple_network_feature_nodes %>% dplyr::select(feature_node_id, x_feature = x, y_feature = y),
    by = "feature_node_id"
  )

p_mofa_v50_selected_triple_feature_network <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = mofa_v50_triple_network_modules,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = view_label),
    alpha = 0.12,
    color = "grey72",
    linewidth = 0.28
  ) +
  ggplot2::geom_curve(
    data = mofa_v50_triple_network_plot_edges,
    ggplot2::aes(x = x_factor, y = y_factor, xend = x_feature, yend = y_feature, color = loading_sign, linewidth = loading_strength),
    curvature = 0.08,
    alpha = 0.43,
    lineend = "round"
  ) +
  ggplot2::geom_label(
    data = mofa_v50_triple_network_factor_nodes,
    ggplot2::aes(x = x, y = y, label = factor_label),
    fill = "#FFF2B3",
    color = "grey10",
    label.size = 0.28,
    size = 2.55,
    lineheight = 0.88,
    fontface = "bold"
  ) +
  ggplot2::geom_point(
    data = mofa_v50_triple_network_feature_nodes,
    ggplot2::aes(x = x, y = y, fill = view_label, size = maximum_loading),
    shape = 21,
    color = "grey20",
    stroke = 0.42
  ) +
  ggplot2::geom_point(
    data = mofa_v50_triple_network_feature_nodes %>% dplyr::filter(shared_feature),
    ggplot2::aes(x = x, y = y),
    shape = 21,
    size = 3.6,
    fill = NA,
    color = "#7A3E9D",
    stroke = 0.72
  ) +
  ggplot2::geom_text(
    data = mofa_v50_triple_network_feature_nodes %>% dplyr::filter(response_feature),
    ggplot2::aes(x = x, y = y, label = "*"),
    nudge_y = 0.24,
    size = 3.1,
    fontface = "bold",
    color = "black"
  ) +
  ggplot2::geom_text(
    data = mofa_v50_triple_network_feature_nodes,
    ggplot2::aes(x = label_x, y = y, label = feature_label_plotmath, hjust = label_hjust),
    parse = TRUE,
    size = 2.20
  ) +
  ggplot2::geom_text(
    data = mofa_v50_triple_network_modules,
    ggplot2::aes(x = module_x, y = label_y, label = module_label),
    fontface = "bold",
    size = 2.35,
    color = "grey20"
  ) +
  ggplot2::scale_color_manual(values = c(Positive = "#D55E00", Negative = "#0072B2"), name = "Feature-weight sign") +
  ggplot2::scale_fill_manual(values = view_colors, name = "Omics view") +
  ggplot2::scale_linewidth_continuous(range = c(0.28, 1.00), guide = "none") +
  ggplot2::scale_size_continuous(range = c(2.3, 3.8), guide = "none") +
  ggplot2::coord_equal(xlim = c(-7.0, 7.0), ylim = c(-5.2, 5.2), clip = "off") +
  ggplot2::labs(
    title = "Selected three-factor feature-weight network",
    subtitle = "Feature modules are placed inside; factor labels are placed outside. Only active views with R² >= 2% contribute features.",
    caption = "Purple outlines denote features shared by multiple factors; asterisks mark nominal p < 0.05 pCR/non-pCR feature contrasts.",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_void(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.6, color = "grey35"),
    plot.caption = ggplot2::element_text(size = 7.1, color = "grey35", hjust = 0),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(5, 165, 5, 165)
  )

if (FALSE) ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_selected_triple_multifactor_feature_loading_network.svg"),
  plot = p_mofa_v50_selected_triple_feature_network,
  width = 12.0,
  height = 8.8,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)


#-----------------------------------------------------------------#
# 22.9 v50 display refinements for violin, 3D and network figures
#-----------------------------------------------------------------#

mofa_v50_make_annotation <- function(plot_data, comparison_name) {
  plot_data %>%
    dplyr::group_by(factor) %>%
    dplyr::summarise(
      y_min = min(factor_score, na.rm = TRUE),
      y_max = max(factor_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      y_position = y_max + 0.20 * pmax(y_max - y_min, 1)
    ) %>%
    dplyr::left_join(
      mofa_factor_wilcoxon_tests %>%
        dplyr::filter(comparison == comparison_name) %>%
        dplyr::transmute(
          factor = factor(as.character(factor), levels = mofa_v50_factor_order),
          wilcoxon_effect,
          wilcoxon_p
        ),
      by = "factor"
    ) %>%
    dplyr::mutate(
      label = paste0(
        "pCR - non-pCR = ", sprintf("%+.2f", wilcoxon_effect),
        "\np = ",
        ifelse(
          !is.finite(wilcoxon_p),
          "NA",
          ifelse(wilcoxon_p < 0.001, "< 0.001", mofa_v50_format_p(wilcoxon_p))
        )
      )
    )
}

mofa_v50_build_response_plot <- function(plot_data, annotation_data, title_text) {
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = TRG_plot, y = factor_score, fill = TRG_plot)
  ) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.68, color = "grey30", linewidth = 0.38) +
    ggplot2::geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      fill = NA,
      coef = 0,
      staplewidth = 0,
      linewidth = 0.42
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.08, height = 0),
      shape = 21,
      size = 1.45,
      color = "grey20",
      stroke = 0.32,
      alpha = 0.85
    ) +
    ggplot2::geom_segment(
      data = annotation_data,
      ggplot2::aes(x = 1, xend = 2, y = y_position, yend = y_position),
      inherit.aes = FALSE,
      linewidth = 0.36,
      color = "grey25"
    ) +
    ggplot2::geom_segment(
      data = annotation_data,
      ggplot2::aes(
        x = 1,
        xend = 1,
        y = y_position,
        yend = y_position - 0.035 * pmax(y_max - y_min, 1)
      ),
      inherit.aes = FALSE,
      linewidth = 0.36,
      color = "grey25"
    ) +
    ggplot2::geom_segment(
      data = annotation_data,
      ggplot2::aes(
        x = 2,
        xend = 2,
        y = y_position,
        yend = y_position - 0.035 * pmax(y_max - y_min, 1)
      ),
      inherit.aes = FALSE,
      linewidth = 0.36,
      color = "grey25"
    ) +
    ggplot2::geom_text(
      data = annotation_data,
      ggplot2::aes(
        x = 1.5,
        y = y_position + 0.08 * pmax(y_max - y_min, 1),
        label = label
      ),
      inherit.aes = FALSE,
      size = 2.25,
      lineheight = 0.92,
      color = "grey20"
    ) +
    ggplot2::scale_fill_manual(
      values = mofa_group_colors[c("pCR", "non_pCR")],
      guide = "none",
      drop = FALSE
    ) +
    ggplot2::facet_wrap(~ factor, nrow = 2, scales = "free_y") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.42))) +
    ggplot2::labs(title = title_text, x = NULL, y = "MOFA factor score") +
    ggplot2::theme_classic(base_size = 8.8) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "plain", size = 8.0),
      axis.text.x = ggplot2::element_text(size = 9.0, face = "plain"),
      axis.text.y = ggplot2::element_text(face = "plain"),
      panel.spacing = grid::unit(0.75, "lines")
    )
}

mofa_v50_save_violin_variants <- function(
  plot_object,
  base_filename,
  width_current,
  height_current
) {
  stem <- tools::file_path_sans_ext(base_filename)
  unlink(file.path(mofa_v50_figure_dir, paste0(stem, c("_w60.svg", "_w70.svg", "_w80.svg"))))
  ggplot2::ggsave(
    filename = file.path(mofa_v50_figure_dir, paste0(stem, "_w90.svg")),
    plot = plot_object,
    width = width_current * 0.90,
    height = height_current,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

mofa_v50_overall_annotation <- mofa_v50_make_annotation(
  mofa_v50_overall_response_long,
  "overall_response_subject_mean"
)
mofa_v50_baseline_annotation <- mofa_v50_make_annotation(
  mofa_v50_baseline_response_long,
  "baseline_response"
)
mofa_v50_after_annotation <- mofa_v50_make_annotation(
  mofa_v50_after_response_long,
  "ongoing_response"
)

p_mofa_v50_overall_response <- mofa_v50_build_response_plot(
  mofa_v50_overall_response_long,
  mofa_v50_overall_annotation,
  "pCR versus non-pCR: subject mean across available timepoints"
)
ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_overall_subject_mean.svg"),
  plot = p_mofa_v50_overall_response,
  width = 7.15,
  height = max(5.4, 2.25 * ceiling(length(mofa_v50_factor_order) / 5)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
mofa_v50_save_violin_variants(
  p_mofa_v50_overall_response,
  "MOFA_v50_factor_response_overall_subject_mean.svg",
  7.15,
  max(5.4, 2.25 * ceiling(length(mofa_v50_factor_order) / 5))
)

p_mofa_v50_baseline_response <- mofa_v50_build_response_plot(
  mofa_v50_baseline_response_long,
  mofa_v50_baseline_annotation,
  "pCR versus non-pCR at Baseline"
)
ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_baseline.svg"),
  plot = p_mofa_v50_baseline_response,
  width = 7.15,
  height = max(5.4, 2.25 * ceiling(length(mofa_v50_factor_order) / 5)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
mofa_v50_save_violin_variants(
  p_mofa_v50_baseline_response,
  "MOFA_v50_factor_response_baseline.svg",
  7.15,
  max(5.4, 2.25 * ceiling(length(mofa_v50_factor_order) / 5))
)

p_mofa_v50_after_response <- mofa_v50_build_response_plot(
  mofa_v50_after_response_long,
  mofa_v50_after_annotation,
  "pCR versus non-pCR after radiotherapy"
)
ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_after_RT.svg"),
  plot = p_mofa_v50_after_response,
  width = 7.15,
  height = max(5.4, 2.25 * ceiling(length(mofa_v50_factor_order) / 5)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
mofa_v50_save_violin_variants(
  p_mofa_v50_after_response,
  "MOFA_v50_factor_response_after_RT.svg",
  7.15,
  max(5.4, 2.25 * ceiling(length(mofa_v50_factor_order) / 5))
)

mofa_v50_delta_annotation <- mofa_v50_paired_factor_data %>%
  dplyr::group_by(factor) %>%
  dplyr::summarise(
    y_min = min(factor_delta, na.rm = TRUE),
    y_max = max(factor_delta, na.rm = TRUE),
    delta_difference =
      median(factor_delta[TRG_plot == "pCR"], na.rm = TRUE) -
      median(factor_delta[TRG_plot == "non_pCR"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    mofa_factor_wilcoxon_tests %>%
      dplyr::filter(comparison == "paired_differential_change") %>%
      dplyr::transmute(
        factor = factor(as.character(factor), levels = mofa_v50_factor_order),
        wilcoxon_p
      ),
    by = "factor"
  ) %>%
  dplyr::mutate(
    y_position = y_max + 0.20 * pmax(y_max - y_min, 1),
    label = paste0(
      "Delta difference = ", sprintf("%+.2f", delta_difference),
      "\np = ",
      ifelse(
        !is.finite(wilcoxon_p),
        "NA",
        ifelse(wilcoxon_p < 0.001, "< 0.001", mofa_v50_format_p(wilcoxon_p))
      )
    )
  )

p_mofa_v50_paired_change_by_response <- ggplot2::ggplot(
  mofa_v50_paired_factor_data,
  ggplot2::aes(x = TRG_plot, y = factor_delta, fill = TRG_plot)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.40) +
  ggplot2::geom_violin(trim = FALSE, alpha = 0.68, color = "grey30", linewidth = 0.38) +
  ggplot2::geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    fill = NA,
    coef = 0,
    staplewidth = 0,
    linewidth = 0.42
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.08, height = 0),
    shape = 21,
    size = 1.55,
    color = "grey20",
    stroke = 0.32
  ) +
  ggplot2::geom_segment(
    data = mofa_v50_delta_annotation,
    ggplot2::aes(x = 1, xend = 2, y = y_position, yend = y_position),
    inherit.aes = FALSE,
    linewidth = 0.36,
    color = "grey25"
  ) +
  ggplot2::geom_text(
    data = mofa_v50_delta_annotation,
    ggplot2::aes(
      x = 1.5,
      y = y_position + 0.08 * pmax(y_max - y_min, 1),
      label = label
    ),
    inherit.aes = FALSE,
    size = 2.10,
    lineheight = 0.92
  ) +
  ggplot2::scale_fill_manual(
    values = mofa_group_colors[c("pCR", "non_pCR")],
    guide = "none",
    drop = FALSE
  ) +
  ggplot2::facet_wrap(~ factor, nrow = 2, scales = "free_y") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.42))) +
  ggplot2::labs(
    title = "Response-specific paired treatment change",
    subtitle = "Delta is After RT minus Baseline for each subject; positive effects indicate a larger increase in pCR.",
    x = NULL,
    y = "Within-subject factor delta"
  ) +
  ggplot2::theme_classic(base_size = 8.8) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.7, color = "grey35"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "plain", size = 8.0),
    axis.text.x = ggplot2::element_text(size = 9.0, face = "plain"),
    axis.text.y = ggplot2::element_text(face = "plain"),
    panel.spacing = grid::unit(0.75, "lines")
  )
ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_paired_change_by_response.svg"),
  plot = p_mofa_v50_paired_change_by_response,
  width = 7.15,
  height = max(5.4, 2.30 * ceiling(length(mofa_v50_factor_order) / 5)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
mofa_v50_save_violin_variants(
  p_mofa_v50_paired_change_by_response,
  "MOFA_v50_factor_paired_change_by_response.svg",
  7.15,
  max(5.4, 2.30 * ceiling(length(mofa_v50_factor_order) / 5))
)

mofa_v50_selected_baseline_data <- mofa_v50_baseline_response_long %>%
  dplyr::filter(as.character(factor) == mofa_v50_selected_factor)

mofa_v50_selected_baseline_annotation <- mofa_v50_baseline_annotation %>%
  dplyr::filter(as.character(factor) == mofa_v50_selected_factor)

p_mofa_v50_selected_factor_baseline <- ggplot2::ggplot(
  mofa_v50_selected_baseline_data,
  ggplot2::aes(x = TRG_plot, y = factor_score, fill = TRG_plot)
) +
  ggplot2::geom_violin(trim = FALSE, alpha = 0.68, color = "grey30") +
  ggplot2::geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    fill = NA,
    coef = 0,
    staplewidth = 0,
    linewidth = 0.42
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.08),
    shape = 21,
    size = 1.8,
    color = "grey20"
  ) +
  ggplot2::geom_segment(
    data = mofa_v50_selected_baseline_annotation,
    ggplot2::aes(x = 1, xend = 2, y = y_position, yend = y_position),
    inherit.aes = FALSE,
    linewidth = 0.38,
    color = "grey25"
  ) +
  ggplot2::geom_text(
    data = mofa_v50_selected_baseline_annotation,
    ggplot2::aes(
      x = 1.5,
      y = y_position + 0.08 * pmax(y_max - y_min, 1),
      label = label
    ),
    inherit.aes = FALSE,
    size = 2.45,
    lineheight = 0.92
  ) +
  ggplot2::scale_fill_manual(values = mofa_group_colors[c("pCR", "non_pCR")], guide = "none") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.42))) +
  ggplot2::labs(
    title = paste0(mofa_v50_selected_factor, ": Baseline response contrast"),
    x = NULL,
    y = "Factor score"
  ) +
  ggplot2::theme_classic(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(size = 10.0, face = "plain"),
    axis.text.y = ggplot2::element_text(face = "plain")
  )
ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_selected_factor_baseline.svg"),
  plot = p_mofa_v50_selected_factor_baseline,
  width = 3.0,
  height = 3.5,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)
mofa_v50_save_violin_variants(
  p_mofa_v50_selected_factor_baseline,
  "MOFA_v50_selected_factor_baseline.svg",
  3.0,
  3.5
)

# Make module labels more readable in the five-factor network.
mofa_v50_network_modules <- mofa_v50_network_modules %>%
  dplyr::arrange(view_order, module) %>%
  dplyr::mutate(
    module_index = dplyr::row_number(),
    module_angle = pi / 2 - 2 * pi * (module_index - 1) / dplyr::n(),
    module_x = 2.80 * cos(module_angle),
    module_y = 2.30 * sin(module_angle)
  )

mofa_v50_network_feature_nodes <- mofa_v50_network_feature_nodes %>%
  dplyr::select(-dplyr::any_of(c("module_x", "module_y", "x", "y", "label_x", "label_hjust", "feature_offset"))) %>%
  dplyr::left_join(
    mofa_v50_network_modules %>% dplyr::select(module_label, module_x, module_y),
    by = "module_label"
  ) %>%
  dplyr::group_by(module_label) %>%
  dplyr::mutate(
    feature_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.34,
    x = module_x,
    y = module_y + feature_offset,
    label_x = x + ifelse(module_x < -0.25, -0.62, 0.62),
    label_hjust = ifelse(module_x < -0.25, 1, 0)
  ) %>%
  dplyr::ungroup()

mofa_v50_network_modules <- mofa_v50_network_modules %>%
  dplyr::select(-dplyr::any_of(c("xmin", "xmax", "ymin", "ymax", "label_y"))) %>%
  dplyr::left_join(
    mofa_v50_network_feature_nodes %>%
      dplyr::group_by(module_label) %>%
      dplyr::summarise(
        xmin = min(x) - 0.36,
        xmax = max(x) + 0.36,
        ymin = min(y) - 0.24,
        ymax = max(y) + 0.24,
        label_y = max(y) + 0.54,
        .groups = "drop"
      ),
    by = "module_label"
  )

mofa_v50_network_factor_nodes <- mofa_v50_network_factor_nodes %>%
  dplyr::mutate(
    x = 6.20 * cos(factor_angle),
    y = 4.90 * sin(factor_angle),
    factor_label = paste0(
      factor,
      ifelse(
        is.finite(minimum_response_p),
        paste0(
          "\npCR - non-pCR = ", sprintf("%+.2f", strongest_response_effect),
          "\np = ",
          ifelse(
            !is.finite(minimum_response_p),
            "NA",
            ifelse(minimum_response_p < 0.001, "< 0.001", mofa_v50_format_p(minimum_response_p))
          )
        ),
        ""
      )
    )
  )

mofa_v50_network_plot_edges <- mofa_v50_network_edges %>%
  dplyr::left_join(
    mofa_v50_network_factor_nodes %>% dplyr::select(factor, x_factor = x, y_factor = y),
    by = "factor"
  ) %>%
  dplyr::left_join(
    mofa_v50_network_feature_nodes %>% dplyr::select(feature_node_id, x_feature = x, y_feature = y),
    by = "feature_node_id"
  )

p_mofa_v50_multifactor_feature_network <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = mofa_v50_network_modules,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = view_label),
    alpha = 0.10,
    color = "grey72",
    linewidth = 0.28
  ) +
  ggplot2::geom_curve(
    data = mofa_v50_network_plot_edges,
    ggplot2::aes(x = x_factor, y = y_factor, xend = x_feature, yend = y_feature, color = loading_sign, linewidth = loading_strength),
    curvature = 0.04,
    alpha = 0.40,
    lineend = "round"
  ) +
  ggplot2::geom_label(
    data = mofa_v50_network_factor_nodes,
    ggplot2::aes(x = x, y = y, label = factor_label),
    fill = "#FFF2B3",
    color = "grey10",
    label.size = 0.28,
    size = 2.35,
    lineheight = 0.90,
    fontface = "plain"
  ) +
  ggplot2::geom_point(
    data = mofa_v50_network_feature_nodes,
    ggplot2::aes(x = x, y = y, fill = view_label, size = maximum_loading),
    shape = 21,
    color = "grey20",
    stroke = 0.42
  ) +
  ggplot2::geom_point(
    data = mofa_v50_network_feature_nodes %>% dplyr::filter(shared_feature),
    ggplot2::aes(x = x, y = y),
    shape = 21,
    size = 3.8,
    fill = NA,
    color = "#7A3E9D",
    stroke = 0.72
  ) +
  ggplot2::geom_text(
    data = mofa_v50_network_feature_nodes %>% dplyr::filter(response_feature),
    ggplot2::aes(x = x, y = y, label = "*"),
    nudge_y = 0.24,
    size = 3.0,
    fontface = "plain",
    color = "black"
  ) +
  ggplot2::geom_text(
    data = mofa_v50_network_feature_nodes,
    ggplot2::aes(x = label_x, y = y, label = feature_label_plotmath, hjust = label_hjust),
    parse = TRUE,
    size = 2.05
  ) +
  ggplot2::geom_label(
    data = mofa_v50_network_modules,
    ggplot2::aes(x = module_x, y = label_y, label = module_label),
    fill = "white",
    color = "grey15",
    label.size = 0.18,
    size = 2.55,
    fontface = "plain",
    label.padding = grid::unit(0.10, "lines")
  ) +
  ggplot2::scale_color_manual(values = c(Positive = "#D55E00", Negative = "#0072B2"), name = "Feature-weight sign") +
  ggplot2::scale_fill_manual(values = view_colors, name = "Omics view") +
  ggplot2::scale_linewidth_continuous(range = c(0.28, 1.00), guide = "none") +
  ggplot2::scale_size_continuous(range = c(2.1, 3.5), guide = "none") +
  ggplot2::coord_equal(xlim = c(-8.0, 8.0), ylim = c(-6.0, 6.0), clip = "off") +
  ggplot2::labs(
    title = "Multi-factor feature-weight network",
    subtitle = paste0(
      "Five QC-passing factors are ranked primarily by response p-value evidence, while active-view feature p-values and multi-view support are used to break ties. ",
      "No feature is shown from a view with R² < 2%."
    ),
    caption = paste0(
      "Modules are Ward.D2 clusters of signed feature-weight profiles within each omics view. ",
      "Purple outlines denote features connected to multiple factors; asterisks mark nominal p < 0.05 pCR/non-pCR feature contrasts."
    ),
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_void(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.6, color = "grey35"),
    plot.caption = ggplot2::element_text(size = 7.1, color = "grey35", hjust = 0),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(5, 200, 5, 200)
  )
if (FALSE) ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_multifactor_feature_loading_network.svg"),
  plot = p_mofa_v50_multifactor_feature_network,
  width = 13.0,
  height = 9.5,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

# Signed median contrast note for the paired-change association heatmap / atlas.
p_mofa_v50_association_delta_heatmap <- p_mofa_v50_association_delta_heatmap +
  ggplot2::labs(
    subtitle = paste0(
      "Text is the signed median contrast, i.e. the median of the first group minus the median of the second group; ",
      "for paired-change panels it summarizes the median within-subject delta contrast."
    ),
    caption = "Signed median contrast = median(group A) - median(group B); positive values indicate a larger value in the first-named group."
  )

p_mofa_v50_factor_interpretation_atlas <-
  (p_mofa_v50_variance_heatmap_atlas | p_mofa_v50_response_enrichment_atlas | p_mofa_v50_association_delta_heatmap) +
  patchwork::plot_layout(widths = c(0.92, 1.05, 0.62), guides = "collect") +
  patchwork::plot_annotation(
    title = "MOFA factor interpretation atlas",
    subtitle = paste0(
      "Variance explained provides the omics context, response enrichment is the primary clinical summary, ",
      "and the right heatmap highlights paired-change associations."
    )
  ) &
  ggplot2::theme(legend.position = "right")

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_variance_association_atlas.svg"),
  plot = p_mofa_v50_factor_interpretation_atlas,
  width = 10.8,
  height = max(6.2, 3.1 + 0.36 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

# Re-select the three-factor combinations using the top five overall-response Wilcoxon factors, with Factor7 fixed.
mofa_v50_overall_response_rank <- mofa_factor_wilcoxon_tests %>%
  dplyr::filter(
    comparison == "overall_response_subject_mean",
    factor %in% mofa_v50_factor_order,
    is.finite(wilcoxon_p)
  ) %>%
  dplyr::arrange(wilcoxon_p, dplyr::desc(abs(wilcoxon_effect))) %>%
  dplyr::mutate(factor = as.character(factor))

mofa_v50_top5_response_factors <- unique(c(
  if ("Factor7" %in% mofa_v50_factor_order) "Factor7" else character(0),
  mofa_v50_overall_response_rank$factor
))
mofa_v50_top5_response_factors <- mofa_v50_top5_response_factors[seq_len(min(5, length(mofa_v50_top5_response_factors)))]
if (!("Factor7" %in% mofa_v50_top5_response_factors) && "Factor7" %in% mofa_v50_factor_order) {
  mofa_v50_top5_response_factors <- c("Factor7", setdiff(mofa_v50_top5_response_factors, "Factor7"))
  mofa_v50_top5_response_factors <- unique(mofa_v50_top5_response_factors)[seq_len(min(5, length(unique(mofa_v50_top5_response_factors))))]
}
mofa_v50_anchor_factor <- if ("Factor7" %in% mofa_v50_top5_response_factors) "Factor7" else mofa_v50_top5_response_factors[1]
mofa_v50_triple_candidate_pool <- setdiff(mofa_v50_top5_response_factors, mofa_v50_anchor_factor)
if (length(mofa_v50_triple_candidate_pool) < 2) {
  mofa_v50_triple_candidate_pool <- setdiff(mofa_v50_factor_order, mofa_v50_anchor_factor)[seq_len(min(4, length(setdiff(mofa_v50_factor_order, mofa_v50_anchor_factor))))]
}
mofa_v50_factor_triples <- lapply(
  utils::combn(mofa_v50_triple_candidate_pool, 2, simplify = FALSE),
  function(partners) c(mofa_v50_anchor_factor, partners)
)

mofa_v50_variance_explained_long <- mofa_variance_explained %>%
  dplyr::mutate(factor = as.character(factor), view = as.character(view))
if (!"view_label" %in% colnames(mofa_v50_variance_explained_long)) {
  if (exists("view_labels")) {
    mofa_v50_variance_explained_long <- mofa_v50_variance_explained_long %>%
      dplyr::left_join(
        data.frame(view = names(view_labels), view_label = unname(view_labels), stringsAsFactors = FALSE),
        by = "view"
      )
  } else {
    mofa_v50_variance_explained_long$view_label <- mofa_v50_variance_explained_long$view
  }
}
mofa_v50_variance_explained_long$view_label[is.na(mofa_v50_variance_explained_long$view_label)] <- mofa_v50_variance_explained_long$view[is.na(mofa_v50_variance_explained_long$view_label)]

mofa_v50_host_active_factors <- unique(
  mofa_v50_variance_explained_long %>%
    dplyr::filter(
      r2 >= mofa_active_view_r2,
      view_label == "Host RNA-seq" | grepl("host", view, ignore.case = TRUE)
    ) %>%
    dplyr::pull(factor)
)

mofa_v50_triple_permanova_atlas <- dplyr::bind_rows(
  lapply(
    seq_along(mofa_v50_factor_triples),
    function(triple_index) {
      triple_factors <- mofa_v50_factor_triples[[triple_index]]
      triple_data <- mofa_v50_subject_mean_factor_data %>%
        dplyr::select(SubjectID, TRG_plot, dplyr::all_of(triple_factors)) %>%
        dplyr::filter(
          dplyr::if_all(dplyr::all_of(triple_factors), is.finite),
          !is.na(TRG_plot)
        ) %>%
        dplyr::mutate(TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR")))

      score_matrix <- scale(as.matrix(triple_data[, triple_factors, drop = FALSE]))
      score_matrix[!is.finite(score_matrix)] <- 0
      distance_object <- stats::dist(score_matrix)

      set.seed(20261701 + triple_index)
      permanova_fit <- tryCatch(
        vegan::adonis2(
          distance_object ~ TRG_plot,
          data = triple_data,
          permutations = mofa_v50_outlier_sensitivity_permutations
        ),
        error = function(e) NULL
      )
      dispersion_fit <- tryCatch(
        vegan::betadisper(distance_object, triple_data$TRG_plot),
        error = function(e) NULL
      )
      dispersion_test <- if (is.null(dispersion_fit)) NULL else tryCatch(
        vegan::permutest(dispersion_fit, permutations = mofa_v50_outlier_sensitivity_permutations),
        error = function(e) NULL
      )

      outlier_rows <- as.integer(
        unlist(
          lapply(
            levels(triple_data$TRG_plot),
            function(group_name) {
              group_rows <- which(triple_data$TRG_plot == group_name)
              group_matrix <- score_matrix[group_rows, , drop = FALSE]
              group_distance <- sqrt(
                rowSums((group_matrix - rep(colMeans(group_matrix), each = nrow(group_matrix)))^2)
              )
              group_rows[which.max(group_distance)]
            }
          )
        )
      )

      sensitivity_rows <- dplyr::bind_rows(
        lapply(
          outlier_rows,
          function(excluded_row) {
            reduced_data <- triple_data[-excluded_row, , drop = FALSE]
            reduced_matrix <- scale(as.matrix(reduced_data[, triple_factors, drop = FALSE]))
            reduced_matrix[!is.finite(reduced_matrix)] <- 0
            set.seed(20261801 + triple_index + excluded_row)
            reduced_fit <- tryCatch(
              vegan::adonis2(
                stats::dist(reduced_matrix) ~ TRG_plot,
                data = reduced_data,
                permutations = mofa_v50_outlier_sensitivity_permutations
              ),
              error = function(e) NULL
            )
            data.frame(
              excluded_subject = triple_data$SubjectID[excluded_row],
              permanova_r2 = if (is.null(reduced_fit)) NA_real_ else as.numeric(reduced_fit$R2[1]),
              permanova_p = if (is.null(reduced_fit)) NA_real_ else as.numeric(reduced_fit$`Pr(>F)`[1]),
              stringsAsFactors = FALSE
            )
          }
        )
      )

      data.frame(
        triple_id = paste(triple_factors, collapse = "__"),
        factor_x = triple_factors[1],
        factor_y = triple_factors[2],
        factor_z = triple_factors[3],
        n_subjects = nrow(triple_data),
        permanova_f = if (is.null(permanova_fit)) NA_real_ else as.numeric(permanova_fit$F[1]),
        permanova_r2 = if (is.null(permanova_fit)) NA_real_ else as.numeric(permanova_fit$R2[1]),
        permanova_p = if (is.null(permanova_fit)) NA_real_ else as.numeric(permanova_fit$`Pr(>F)`[1]),
        dispersion_p = if (is.null(dispersion_test)) NA_real_ else as.numeric(dispersion_test$tab$`Pr(>F)`[1]),
        outlier_candidates = paste(triple_data$SubjectID[outlier_rows], collapse = ";"),
        sensitivity_max_p = if (any(is.finite(sensitivity_rows$permanova_p))) max(sensitivity_rows$permanova_p, na.rm = TRUE) else NA_real_,
        sensitivity_min_r2 = if (any(is.finite(sensitivity_rows$permanova_r2))) min(sensitivity_rows$permanova_r2, na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  )
) %>%
  dplyr::mutate(
    partner_has_host_view = factor_y %in% mofa_v50_host_active_factors | factor_z %in% mofa_v50_host_active_factors,
    permanova_fdr = stats::p.adjust(permanova_p, method = "BH"),
    dispersion_ok = is.na(dispersion_p) | dispersion_p >= 0.05,
    outlier_robust = is.na(sensitivity_max_p) | sensitivity_max_p < 0.10
  ) %>%
  dplyr::arrange(
    dplyr::desc(partner_has_host_view),
    dplyr::desc(dispersion_ok),
    permanova_p,
    dplyr::desc(permanova_r2),
    dplyr::desc(outlier_robust)
  )

write.csv(
  mofa_v50_triple_permanova_atlas,
  file.path(mofa_v50_result_dir, "MOFA_v50_three_factor_PERMANOVA_atlas.csv"),
  row.names = FALSE
)

mofa_v50_selected_triple_diagnostics <- mofa_v50_triple_permanova_atlas %>%
  dplyr::slice_head(n = min(3, nrow(mofa_v50_triple_permanova_atlas))) %>%
  dplyr::mutate(selection_rank = dplyr::row_number())

write.csv(
  mofa_v50_selected_triple_diagnostics,
  file.path(mofa_v50_result_dir, "MOFA_v50_selected_three_factor_diagnostics.csv"),
  row.names = FALSE
)

mofa_v50_build_3d_scatter <- function(triple_factors, diagnostic_row) {
  plot_data <- mofa_v50_subject_mean_factor_data %>%
    dplyr::select(SubjectID, TRG_plot, dplyr::all_of(triple_factors)) %>%
    dplyr::filter(
      dplyr::if_all(dplyr::all_of(triple_factors), is.finite),
      !is.na(TRG_plot)
    ) %>%
    dplyr::mutate(TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR")))

  score_matrix <- scale(as.matrix(plot_data[, triple_factors, drop = FALSE]))
  score_matrix[!is.finite(score_matrix)] <- 0
  plot_data$x3 <- score_matrix[, 1]
  plot_data$y3 <- score_matrix[, 2]
  plot_data$z3 <- score_matrix[, 3]

  span <- apply(score_matrix, 2, range)
  pad <- apply(score_matrix, 2, function(v) max(0.20 * diff(range(v)), 0.5))
  x_range <- c(span[1, 1] - pad[1], span[2, 1] + pad[1])
  y_range <- c(span[1, 2] - pad[2], span[2, 2] + pad[2])
  z_range <- c(span[1, 3] - pad[3], span[2, 3] + pad[3])
  floor_z <- z_range[1]
  side_x <- x_range[1]
  back_y <- y_range[1]

  projected_points <- mofa_v50_rotate_project_3d(plot_data$x3, plot_data$y3, plot_data$z3, azimuth_deg = 52, elevation_deg = 26, distance = 10.5)
  projected_floor <- mofa_v50_rotate_project_3d(plot_data$x3, plot_data$y3, rep(floor_z, nrow(plot_data)), azimuth_deg = 52, elevation_deg = 26, distance = 10.5)
  projected_side <- mofa_v50_rotate_project_3d(rep(side_x, nrow(plot_data)), plot_data$y3, plot_data$z3, azimuth_deg = 52, elevation_deg = 26, distance = 10.5)
  projected_back <- mofa_v50_rotate_project_3d(plot_data$x3, rep(back_y, nrow(plot_data)), plot_data$z3, azimuth_deg = 52, elevation_deg = 26, distance = 10.5)

  plot_data$x_proj <- projected_points$x
  plot_data$y_proj <- projected_points$y
  plot_data$floor_x <- projected_floor$x
  plot_data$floor_y <- projected_floor$y
  plot_data$side_x_proj <- projected_side$x
  plot_data$side_y_proj <- projected_side$y
  plot_data$back_x_proj <- projected_back$x
  plot_data$back_y_proj <- projected_back$y

  finite_depth <- is.finite(projected_points$depth)
  if (sum(finite_depth) >= 2 && diff(range(projected_points$depth[finite_depth])) > 0) {
    plot_data$depth_scale <- scales::rescale(projected_points$depth, to = c(2.0, 3.3), from = range(projected_points$depth[finite_depth]))
  } else {
    plot_data$depth_scale <- rep(2.6, nrow(plot_data))
  }
  plot_data$depth_scale[!is.finite(plot_data$depth_scale)] <- 2.6
  plot_data <- plot_data %>%
    dplyr::filter(
      is.finite(x_proj), is.finite(y_proj),
      is.finite(floor_x), is.finite(floor_y),
      is.finite(side_x_proj), is.finite(side_y_proj),
      is.finite(back_x_proj), is.finite(back_y_proj),
      is.finite(depth_scale)
    ) %>%
    dplyr::arrange(x3 + y3 + z3)

  box_data <- mofa_v50_box_edges(x_range, y_range, z_range)
  box_projected <- mofa_v50_rotate_project_3d(box_data$x, box_data$y, box_data$z, azimuth_deg = 52, elevation_deg = 26, distance = 10.5)
  box_data$x_proj <- box_projected$x
  box_data$y_proj <- box_projected$y

  axis_data <- data.frame(
    axis = triple_factors,
    x3 = c(x_range[2], side_x, side_x),
    y3 = c(back_y, y_range[2], back_y),
    z3 = c(floor_z, floor_z, z_range[2]),
    stringsAsFactors = FALSE
  )
  axis_projected <- mofa_v50_rotate_project_3d(axis_data$x3, axis_data$y3, axis_data$z3, azimuth_deg = 52, elevation_deg = 26, distance = 10.5)
  axis_data$x_proj <- axis_projected$x
  axis_data$y_proj <- axis_projected$y

  mean_points <- plot_data %>%
    dplyr::group_by(TRG_plot) %>%
    dplyr::summarise(
      x3 = mean(x3), y3 = mean(y3), z3 = mean(z3), .groups = "drop"
    )
  mean_projected <- mofa_v50_rotate_project_3d(mean_points$x3, mean_points$y3, mean_points$z3, azimuth_deg = 52, elevation_deg = 26, distance = 10.5)
  mean_points$x_proj <- mean_projected$x
  mean_points$y_proj <- mean_projected$y

  annotation_label <- paste0(
    "Rank ", diagnostic_row$selection_rank[1], ": ", paste(triple_factors, collapse = " + "),
    "\nPERMANOVA R² = ", sprintf("%.2f", diagnostic_row$permanova_r2[1]),
    "\nPERMANOVA p = ", ifelse(is.na(diagnostic_row$permanova_p[1]), "NA", ifelse(diagnostic_row$permanova_p[1] < 0.001, "< 0.001", sprintf("%.3f", diagnostic_row$permanova_p[1]))),
    "\nDispersion p = ", ifelse(is.na(diagnostic_row$dispersion_p[1]), "NA", mofa_v50_format_p(diagnostic_row$dispersion_p[1])),
    "\nHost RNA-seq partner = ", ifelse(diagnostic_row$partner_has_host_view[1], "yes", "no")
  )

  x_limits <- range(c(box_data$x_proj, plot_data$x_proj, plot_data$floor_x, plot_data$side_x_proj, plot_data$back_x_proj), na.rm = TRUE)
  y_limits <- range(c(box_data$y_proj, plot_data$y_proj, plot_data$floor_y, plot_data$side_y_proj, plot_data$back_y_proj), na.rm = TRUE)
  x_margin <- 0.08 * diff(x_limits)
  y_margin <- 0.10 * diff(y_limits)

  ggplot2::ggplot() +
    ggplot2::geom_path(
      data = box_data,
      ggplot2::aes(x = x_proj, y = y_proj, group = edge_id),
      color = "grey72",
      linewidth = 0.40
    ) +
    ggplot2::geom_segment(
      data = plot_data,
      ggplot2::aes(x = floor_x, y = floor_y, xend = x_proj, yend = y_proj, color = TRG_plot),
      linewidth = 0.26,
      alpha = 0.18,
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = plot_data,
      ggplot2::aes(x = side_x_proj, y = side_y_proj, xend = x_proj, yend = y_proj, color = TRG_plot),
      linewidth = 0.18,
      alpha = 0.12,
      show.legend = FALSE
    ) +
    ggplot2::geom_segment(
      data = plot_data,
      ggplot2::aes(x = back_x_proj, y = back_y_proj, xend = x_proj, yend = y_proj, color = TRG_plot),
      linewidth = 0.18,
      alpha = 0.12,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(x = floor_x, y = floor_y, fill = TRG_plot),
      shape = 21, size = 1.55, color = NA, alpha = 0.15, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(x = side_x_proj, y = side_y_proj, fill = TRG_plot),
      shape = 21, size = 1.40, color = NA, alpha = 0.10, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(x = back_x_proj, y = back_y_proj, fill = TRG_plot),
      shape = 21, size = 1.40, color = NA, alpha = 0.10, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(x = x_proj, y = y_proj, fill = TRG_plot),
      size = plot_data$depth_scale,
      shape = 21,
      color = "grey20",
      stroke = 0.44,
      alpha = 0.95,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = mean_points,
      ggplot2::aes(x = x_proj, y = y_proj, fill = TRG_plot),
      shape = 21,
      size = 3.8,
      color = "black",
      stroke = 0.68
    ) +
    ggplot2::geom_label(
      data = axis_data,
      ggplot2::aes(x = x_proj, y = y_proj, label = axis),
      size = 2.75,
      label.size = 0.15,
      fill = "white",
      fontface = "plain"
    ) +
    ggplot2::annotate(
      "text",
      x = x_limits[1],
      y = y_limits[2],
      label = annotation_label,
      hjust = 0,
      vjust = 1,
      lineheight = 0.92,
      size = 2.68,
      color = "black"
    ) +
    ggplot2::scale_fill_manual(values = mofa_response_colors[c("pCR", "non_pCR")], name = "Response") +
    ggplot2::scale_color_manual(values = mofa_response_colors[c("pCR", "non_pCR")], guide = "none") +
    ggplot2::coord_equal(
      xlim = c(x_limits[1] - x_margin, x_limits[2] + x_margin),
      ylim = c(y_limits[1] - y_margin, y_limits[2] + y_margin),
      clip = "off"
    ) +
    ggplot2::labs(
      title = "Overall-response three-factor separation",
      subtitle = "Cabinet-projection scatter with floor and wall projections; lighter wall/floor points are orthogonal shadows.",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_void(base_size = 9.2) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 9.5),
      plot.subtitle = ggplot2::element_text(size = 7.5, color = "grey35"),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(5, 8, 5, 5)
    )
}

mofa_v50_three_factor_plot_list <- lapply(
  seq_len(nrow(mofa_v50_selected_triple_diagnostics)),
  function(i) {
    triple_factors <- c(
      mofa_v50_selected_triple_diagnostics$factor_x[i],
      mofa_v50_selected_triple_diagnostics$factor_y[i],
      mofa_v50_selected_triple_diagnostics$factor_z[i]
    )
    mofa_v50_build_3d_scatter(triple_factors, mofa_v50_selected_triple_diagnostics[i, , drop = FALSE])
  }
)

p_mofa_v50_three_factor_response_maps <- patchwork::wrap_plots(
  mofa_v50_three_factor_plot_list,
  nrow = 1,
  guides = "collect"
)
if (FALSE) ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_three_factor_response_3D_maps.svg"),
  plot = p_mofa_v50_three_factor_response_maps,
  width = max(6.2, 5.8 * length(mofa_v50_three_factor_plot_list)),
  height = 6.0,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

# Rebuild the selected-triple network using the top-ranked triple and clearer module labels.
mofa_v50_triple_network_factors <- c(
  mofa_v50_selected_triple_diagnostics$factor_x[1],
  mofa_v50_selected_triple_diagnostics$factor_y[1],
  mofa_v50_selected_triple_diagnostics$factor_z[1]
)

mofa_v50_triple_active_factor_views <- mofa_variance_explained %>%
  dplyr::mutate(factor = as.character(factor), view = as.character(view)) %>%
  dplyr::filter(
    factor %in% mofa_v50_triple_network_factors,
    r2 >= mofa_active_view_r2
  ) %>%
  dplyr::select(factor, view, view_r2 = r2)

mofa_v50_triple_network_seed_edges <- dplyr::bind_rows(
  mofa_feature_weights %>%
    dplyr::filter(factor %in% mofa_v50_triple_network_factors, display_eligible) %>%
    dplyr::inner_join(mofa_v50_triple_active_factor_views, by = c("factor", "view")) %>%
    dplyr::group_by(factor, view, direction) %>%
    dplyr::slice_max(order_by = abs(weight_within_view), n = mofa_network_features_per_view_direction, with_ties = FALSE) %>%
    dplyr::ungroup(),
  mofa_feature_weights %>%
    dplyr::filter(
      factor %in% mofa_v50_triple_network_factors,
      display_eligible,
      is.finite(p_value),
      p_value < 0.05
    ) %>%
    dplyr::inner_join(mofa_v50_triple_active_factor_views, by = c("factor", "view")) %>%
    dplyr::group_by(factor) %>%
    dplyr::arrange(p_value, dplyr::desc(abs(weight_within_view)), .by_group = TRUE) %>%
    dplyr::slice_head(n = 2) %>%
    dplyr::ungroup()
) %>%
  dplyr::distinct(factor, view, feature, .keep_all = TRUE) %>%
  dplyr::mutate(
    feature_node_id = paste(view, feature, sep = "::"),
    edge_id = paste(factor, feature_node_id, sep = "__")
  )

mofa_v50_triple_network_edges <- mofa_feature_weights %>%
  dplyr::filter(factor %in% mofa_v50_triple_network_factors, display_eligible) %>%
  dplyr::inner_join(mofa_v50_triple_active_factor_views, by = c("factor", "view")) %>%
  dplyr::mutate(
    feature_node_id = paste(view, feature, sep = "::"),
    edge_id = paste(factor, feature_node_id, sep = "__")
  ) %>%
  dplyr::filter(
    feature_node_id %in% unique(mofa_v50_triple_network_seed_edges$feature_node_id),
    abs(weight_within_view) >= mofa_network_shared_loading_threshold |
      edge_id %in% mofa_v50_triple_network_seed_edges$edge_id
  ) %>%
  dplyr::mutate(
    loading_sign = ifelse(weight >= 0, "Positive", "Negative"),
    loading_strength = abs(weight_within_view)
  )

mofa_v50_triple_module_membership <- dplyr::bind_rows(
  lapply(
    split(mofa_v50_triple_network_edges, mofa_v50_triple_network_edges$view),
    function(view_edge_data) {
      profile <- view_edge_data %>%
        dplyr::select(feature_node_id, factor, weight_within_view) %>%
        tidyr::pivot_wider(names_from = factor, values_from = weight_within_view, values_fill = 0)
      for (factor_name in setdiff(mofa_v50_triple_network_factors, colnames(profile))) profile[[factor_name]] <- 0
      profile_matrix <- as.matrix(profile[, mofa_v50_triple_network_factors, drop = FALSE])
      rownames(profile_matrix) <- profile$feature_node_id
      module <- if (nrow(profile_matrix) <= 2) rep(1L, nrow(profile_matrix)) else stats::cutree(
        stats::hclust(stats::dist(profile_matrix), method = "ward.D2"),
        k = min(2L, max(1L, ceiling(nrow(profile_matrix) / 7)))
      )
      data.frame(feature_node_id = rownames(profile_matrix), module = as.integer(module), stringsAsFactors = FALSE)
    }
  )
)

mofa_v50_triple_network_feature_nodes <- mofa_v50_triple_network_edges %>%
  dplyr::group_by(feature_node_id, view, view_label, feature, feature_label) %>%
  dplyr::summarise(
    n_connected_factors = dplyr::n_distinct(factor),
    maximum_loading = max(loading_strength, na.rm = TRUE),
    minimum_response_p = if (any(is.finite(p_value))) min(p_value, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  dplyr::left_join(mofa_v50_triple_module_membership, by = "feature_node_id") %>%
  dplyr::mutate(
    view_order = match(view, required_views),
    module_label = paste0(as.character(view_label), " M", module),
    shared_feature = n_connected_factors >= 2,
    response_feature = is.finite(minimum_response_p) & minimum_response_p < 0.05,
    feature_label_plotmath = mofa_v50_plotmath_feature_label(view, feature_label)
  ) %>%
  dplyr::arrange(view_order, module, dplyr::desc(n_connected_factors), dplyr::desc(maximum_loading))

mofa_v50_triple_network_modules <- mofa_v50_triple_network_feature_nodes %>%
  dplyr::distinct(view, view_label, view_order, module, module_label) %>%
  dplyr::arrange(view_order, module) %>%
  dplyr::mutate(
    module_index = dplyr::row_number(),
    module_angle = pi / 2 - 2 * pi * (module_index - 1) / dplyr::n(),
    module_x = 2.70 * cos(module_angle),
    module_y = 2.25 * sin(module_angle)
  )

mofa_v50_triple_network_feature_nodes <- mofa_v50_triple_network_feature_nodes %>%
  dplyr::left_join(mofa_v50_triple_network_modules %>% dplyr::select(module_label, module_x, module_y), by = "module_label") %>%
  dplyr::group_by(module_label) %>%
  dplyr::mutate(
    feature_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.32,
    x = module_x,
    y = module_y + feature_offset,
    label_x = x + ifelse(module_x < -0.25, -0.60, 0.60),
    label_hjust = ifelse(module_x < -0.25, 1, 0)
  ) %>%
  dplyr::ungroup()

mofa_v50_triple_network_modules <- mofa_v50_triple_network_modules %>%
  dplyr::left_join(
    mofa_v50_triple_network_feature_nodes %>%
      dplyr::group_by(module_label) %>%
      dplyr::summarise(
        xmin = min(x) - 0.34,
        xmax = max(x) + 0.34,
        ymin = min(y) - 0.22,
        ymax = max(y) + 0.22,
        label_y = max(y) + 0.52,
        .groups = "drop"
      ),
    by = "module_label"
  )

mofa_v50_triple_network_factor_nodes <- data.frame(
  factor = mofa_v50_triple_network_factors,
  factor_angle = pi / 2 - 2 * pi * (seq_along(mofa_v50_triple_network_factors) - 1) / length(mofa_v50_triple_network_factors),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    mofa_factor_wilcoxon_tests %>%
      dplyr::filter(comparison == "overall_response_subject_mean", factor %in% mofa_v50_triple_network_factors) %>%
      dplyr::select(factor, strongest_response_effect = wilcoxon_effect, minimum_response_p = wilcoxon_p),
    by = "factor"
  ) %>%
  dplyr::mutate(
    x = 6.10 * cos(factor_angle),
    y = 4.80 * sin(factor_angle),
    factor_label = paste0(
      factor,
      ifelse(
        is.finite(minimum_response_p),
        paste0(
          "\npCR - non-pCR = ", sprintf("%+.2f", strongest_response_effect),
          "\np = ",
          ifelse(minimum_response_p < 0.001, "< 0.001", mofa_v50_format_p(minimum_response_p))
        ),
        ""
      )
    )
  )

mofa_v50_triple_network_plot_edges <- mofa_v50_triple_network_edges %>%
  dplyr::left_join(mofa_v50_triple_network_factor_nodes %>% dplyr::select(factor, x_factor = x, y_factor = y), by = "factor") %>%
  dplyr::left_join(mofa_v50_triple_network_feature_nodes %>% dplyr::select(feature_node_id, x_feature = x, y_feature = y), by = "feature_node_id")

p_mofa_v50_selected_triple_feature_network <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = mofa_v50_triple_network_modules,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = view_label),
    alpha = 0.10,
    color = "grey72",
    linewidth = 0.28
  ) +
  ggplot2::geom_curve(
    data = mofa_v50_triple_network_plot_edges,
    ggplot2::aes(x = x_factor, y = y_factor, xend = x_feature, yend = y_feature, color = loading_sign, linewidth = loading_strength),
    curvature = 0.04,
    alpha = 0.40,
    lineend = "round"
  ) +
  ggplot2::geom_label(
    data = mofa_v50_triple_network_factor_nodes,
    ggplot2::aes(x = x, y = y, label = factor_label),
    fill = "#FFF2B3",
    color = "grey10",
    label.size = 0.28,
    size = 2.55,
    lineheight = 0.90,
    fontface = "plain"
  ) +
  ggplot2::geom_point(
    data = mofa_v50_triple_network_feature_nodes,
    ggplot2::aes(x = x, y = y, fill = view_label, size = maximum_loading),
    shape = 21,
    color = "grey20",
    stroke = 0.42
  ) +
  ggplot2::geom_point(
    data = mofa_v50_triple_network_feature_nodes %>% dplyr::filter(shared_feature),
    ggplot2::aes(x = x, y = y),
    shape = 21,
    size = 3.6,
    fill = NA,
    color = "#7A3E9D",
    stroke = 0.72
  ) +
  ggplot2::geom_text(
    data = mofa_v50_triple_network_feature_nodes %>% dplyr::filter(response_feature),
    ggplot2::aes(x = x, y = y, label = "*"),
    nudge_y = 0.24,
    size = 3.0,
    fontface = "plain",
    color = "black"
  ) +
  ggplot2::geom_text(
    data = mofa_v50_triple_network_feature_nodes,
    ggplot2::aes(x = label_x, y = y, label = feature_label_plotmath, hjust = label_hjust),
    parse = TRUE,
    size = 2.08
  ) +
  ggplot2::geom_label(
    data = mofa_v50_triple_network_modules,
    ggplot2::aes(x = module_x, y = label_y, label = module_label),
    fill = "white",
    color = "grey15",
    label.size = 0.18,
    size = 2.55,
    fontface = "plain",
    label.padding = grid::unit(0.10, "lines")
  ) +
  ggplot2::scale_color_manual(values = c(Positive = "#D55E00", Negative = "#0072B2"), name = "Feature-weight sign") +
  ggplot2::scale_fill_manual(values = view_colors, name = "Omics view") +
  ggplot2::scale_linewidth_continuous(range = c(0.28, 1.00), guide = "none") +
  ggplot2::scale_size_continuous(range = c(2.2, 3.7), guide = "none") +
  ggplot2::coord_equal(xlim = c(-8.0, 8.0), ylim = c(-6.0, 6.0), clip = "off") +
  ggplot2::labs(
    title = "Selected three-factor feature-weight network",
    subtitle = "The highest-ranked three-factor combination is shown; feature modules are placed inside and factor labels outside. Module labels are drawn explicitly to remain visible.",
    caption = "Purple outlines denote features shared by multiple factors; asterisks mark nominal p < 0.05 pCR/non-pCR feature contrasts.",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_void(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.6, color = "grey35"),
    plot.caption = ggplot2::element_text(size = 7.1, color = "grey35", hjust = 0),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(5, 190, 5, 190)
  )
if (FALSE) ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_selected_triple_multifactor_feature_loading_network.svg"),
  plot = p_mofa_v50_selected_triple_feature_network,
  width = 12.8,
  height = 9.1,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)


#-----------------------------------------------------------------#
# 22.9 v50 three-dimensional and network refinements
#-----------------------------------------------------------------#

mofa_v50_energy_distance_test <- function(score_matrix, group_vector, permutations = 999) {
  group_vector <- as.factor(group_vector)
  if (nlevels(group_vector) != 2) {
    return(data.frame(energy_statistic = NA_real_, energy_p = NA_real_, stringsAsFactors = FALSE))
  }
  distance_matrix <- as.matrix(stats::dist(score_matrix))
  compute_energy <- function(group_labels) {
    level_a <- levels(group_labels)[1]
    level_b <- levels(group_labels)[2]
    index_a <- which(group_labels == level_a)
    index_b <- which(group_labels == level_b)
    n_a <- length(index_a)
    n_b <- length(index_b)
    if (n_a < 2 || n_b < 2) return(NA_real_)
    d_ab <- mean(distance_matrix[index_a, index_b, drop = FALSE])
    d_aa <- mean(distance_matrix[index_a, index_a, drop = FALSE])
    d_bb <- mean(distance_matrix[index_b, index_b, drop = FALSE])
    2 * d_ab - d_aa - d_bb
  }
  observed <- compute_energy(group_vector)
  if (!is.finite(observed)) {
    return(data.frame(energy_statistic = NA_real_, energy_p = NA_real_, stringsAsFactors = FALSE))
  }
  permuted <- replicate(
    permutations,
    compute_energy(factor(sample(as.character(group_vector)), levels = levels(group_vector)))
  )
  p_value <- (1 + sum(permuted >= observed, na.rm = TRUE)) / (1 + sum(is.finite(permuted)))
  data.frame(energy_statistic = observed, energy_p = p_value, stringsAsFactors = FALSE)
}

mofa_v50_draw_s3d_panel <- function(triple_factors, diagnostic_row, point_cex = 0.95) {
  if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
    plot.new()
    text(0.5, 0.6, "Package 'scatterplot3d' is required for v50 3D plots.")
    text(0.5, 0.48, paste(triple_factors, collapse = " + "), cex = 0.9)
    return(invisible(NULL))
  }
  plot_data <- mofa_v50_subject_mean_factor_data %>%
    dplyr::select(SubjectID, TRG_plot, dplyr::all_of(triple_factors)) %>%
    dplyr::filter(
      dplyr::if_all(dplyr::all_of(triple_factors), is.finite),
      !is.na(TRG_plot)
    ) %>%
    dplyr::mutate(TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR")))

  score_matrix <- scale(as.matrix(plot_data[, triple_factors, drop = FALSE]))
  score_matrix[!is.finite(score_matrix)] <- 0
  plot_data$x <- score_matrix[, 1]
  plot_data$y <- score_matrix[, 2]
  plot_data$z <- score_matrix[, 3]
  plot_data$color <- unname(mofa_response_colors[as.character(plot_data$TRG_plot)])
  plot_data$line_color <- grDevices::adjustcolor(plot_data$color, alpha.f = 0.20)

  s3d <- scatterplot3d::scatterplot3d(
    x = plot_data$x,
    y = plot_data$y,
    z = plot_data$z,
    pch = 16,
    color = plot_data$color,
    cex.symbols = point_cex,
    type = "h",
    lty.hplot = 1,
    mar = c(2.2, 2.4, 2.4, 1.6),
    main = paste0("Rank ", diagnostic_row$selection_rank[1], ": ", paste(triple_factors, collapse = " + ")),
    xlab = triple_factors[1],
    ylab = triple_factors[2],
    zlab = triple_factors[3],
    angle = 52,
    scale.y = 1.0,
    box = TRUE,
    grid = TRUE
  )

  # Re-draw drop lines and points for clearer group coloring.
  for (i in seq_len(nrow(plot_data))) {
    xy_point <- s3d$xyz.convert(plot_data$x[i], plot_data$y[i], plot_data$z[i])
    xy_floor <- s3d$xyz.convert(plot_data$x[i], plot_data$y[i], min(plot_data$z, na.rm = TRUE))
    graphics::segments(xy_floor$x, xy_floor$y, xy_point$x, xy_point$y, col = plot_data$line_color[i], lwd = 0.6)
  }
  point_xy <- s3d$xyz.convert(plot_data$x, plot_data$y, plot_data$z)
  graphics::points(point_xy$x, point_xy$y, pch = 21, bg = plot_data$color, col = "grey20", cex = point_cex)

  center_data <- plot_data %>%
    dplyr::group_by(TRG_plot) %>%
    dplyr::summarise(x = mean(x), y = mean(y), z = mean(z), .groups = "drop")
  center_xy <- s3d$xyz.convert(center_data$x, center_data$y, center_data$z)
  graphics::points(center_xy$x, center_xy$y, pch = 21, bg = unname(mofa_response_colors[as.character(center_data$TRG_plot)]), col = "black", cex = 1.4)

  stats_label <- paste0(
    "PERMANOVA p = ", ifelse(is.na(diagnostic_row$permanova_p[1]), "NA", ifelse(diagnostic_row$permanova_p[1] < 0.001, "<0.001", sprintf("%.3f", diagnostic_row$permanova_p[1]))),
    " | R2 = ", ifelse(is.na(diagnostic_row$permanova_r2[1]), "NA", sprintf("%.2f", diagnostic_row$permanova_r2[1])),
    "\nEnergy p = ", ifelse(is.na(diagnostic_row$energy_p[1]), "NA", ifelse(diagnostic_row$energy_p[1] < 0.001, "<0.001", sprintf("%.3f", diagnostic_row$energy_p[1]))),
    " | E = ", ifelse(is.na(diagnostic_row$energy_statistic[1]), "NA", sprintf("%.2f", diagnostic_row$energy_statistic[1])),
    "\nHost RNA-seq partner count = ", diagnostic_row$host_partner_count[1]
  )
  graphics::mtext(stats_label, side = 3, line = -1.0, adj = 0.02, cex = 0.72)
  graphics::legend(
    "bottomleft",
    legend = c("pCR", "non-pCR"),
    pt.bg = unname(mofa_response_colors[c("pCR", "non_pCR")]),
    pch = 21,
    col = "grey20",
    pt.cex = 1.0,
    cex = 0.75,
    bty = "n"
  )
  invisible(NULL)
}

# Host RNA-seq active factors for selection preferences.
mofa_v50_variance_explained_long <- mofa_variance_explained %>%
  dplyr::mutate(factor = as.character(factor), view = as.character(view))
if (!"view_label" %in% colnames(mofa_v50_variance_explained_long)) {
  if (exists("view_labels")) {
    mofa_v50_variance_explained_long <- mofa_v50_variance_explained_long %>%
      dplyr::left_join(
        data.frame(view = names(view_labels), view_label = unname(view_labels), stringsAsFactors = FALSE),
        by = "view"
      )
  } else {
    mofa_v50_variance_explained_long$view_label <- mofa_v50_variance_explained_long$view
  }
}
mofa_v50_variance_explained_long$view_label[is.na(mofa_v50_variance_explained_long$view_label)] <- mofa_v50_variance_explained_long$view[is.na(mofa_v50_variance_explained_long$view_label)]

mofa_v50_host_active_factors <- unique(
  mofa_v50_variance_explained_long %>%
    dplyr::filter(
      r2 >= mofa_active_view_r2,
      view_label == "Host RNA-seq" | grepl("host", view, ignore.case = TRUE)
    ) %>%
    dplyr::pull(factor)
)

# Re-rank five-factor network selection to prefer host RNA-seq-active factors.
mofa_v50_network_selection_audit <- mofa_v50_network_selection_audit %>%
  dplyr::mutate(
    host_rnaseq_active = factor %in% mofa_v50_host_active_factors
  ) %>%
  dplyr::arrange(
    selection_tier,
    dplyr::desc(host_rnaseq_active),
    minimum_response_p,
    best_response_feature_p,
    dplyr::desc(sharedness_penalized_r2),
    match(factor, mofa_v50_factor_order)
  )

mofa_v50_network_factors <- utils::head(
  mofa_v50_network_selection_audit$factor[
    mofa_v50_network_selection_audit$eligible_active_view &
      mofa_v50_network_selection_audit$eligible_qc
  ],
  min(mofa_network_factor_count, length(mofa_v50_factor_order))
)

mofa_v50_network_selection_audit <- mofa_v50_network_selection_audit %>%
  dplyr::mutate(selected_for_network = factor %in% mofa_v50_network_factors)

write.csv(
  mofa_v50_network_selection_audit,
  file.path(mofa_v50_result_dir, "MOFA_v50_network_factor_selection_audit.csv"),
  row.names = FALSE
)

# Rebuild the 5-factor network with more space and orthogonal pipe-like edges.
mofa_v50_active_factor_views <- mofa_variance_explained %>%
  dplyr::mutate(factor = as.character(factor), view = as.character(view)) %>%
  dplyr::filter(
    factor %in% mofa_v50_network_factors,
    r2 >= mofa_active_view_r2
  ) %>%
  dplyr::select(factor, view, view_r2 = r2)

mofa_v50_network_seed_edges <- dplyr::bind_rows(
  mofa_feature_weights %>%
    dplyr::filter(factor %in% mofa_v50_network_factors, display_eligible) %>%
    dplyr::inner_join(mofa_v50_active_factor_views, by = c("factor", "view")) %>%
    dplyr::group_by(factor, view, direction) %>%
    dplyr::slice_max(order_by = abs(weight_within_view), n = mofa_network_features_per_view_direction, with_ties = FALSE) %>%
    dplyr::ungroup(),
  mofa_feature_weights %>%
    dplyr::filter(
      factor %in% mofa_v50_network_factors,
      display_eligible,
      is.finite(p_value),
      p_value < 0.05
    ) %>%
    dplyr::inner_join(mofa_v50_active_factor_views, by = c("factor", "view")) %>%
    dplyr::group_by(factor) %>%
    dplyr::arrange(p_value, dplyr::desc(abs(weight_within_view)), .by_group = TRUE) %>%
    dplyr::slice_head(n = 2) %>%
    dplyr::ungroup()
) %>%
  dplyr::distinct(factor, view, feature, .keep_all = TRUE) %>%
  dplyr::mutate(feature_node_id = paste(view, feature, sep = "::"), edge_id = paste(factor, feature_node_id, sep = "__"))

mofa_v50_network_edges <- mofa_feature_weights %>%
  dplyr::filter(factor %in% mofa_v50_network_factors, display_eligible) %>%
  dplyr::inner_join(mofa_v50_active_factor_views, by = c("factor", "view")) %>%
  dplyr::mutate(feature_node_id = paste(view, feature, sep = "::"), edge_id = paste(factor, feature_node_id, sep = "__")) %>%
  dplyr::filter(
    feature_node_id %in% unique(mofa_v50_network_seed_edges$feature_node_id),
    abs(weight_within_view) >= mofa_network_shared_loading_threshold | edge_id %in% mofa_v50_network_seed_edges$edge_id
  ) %>%
  dplyr::mutate(loading_sign = ifelse(weight >= 0, "Positive", "Negative"), loading_strength = abs(weight_within_view))

mofa_v50_module_membership <- dplyr::bind_rows(
  lapply(
    split(mofa_v50_network_edges, mofa_v50_network_edges$view),
    function(view_edge_data) {
      profile <- view_edge_data %>%
        dplyr::select(feature_node_id, factor, weight_within_view) %>%
        tidyr::pivot_wider(names_from = factor, values_from = weight_within_view, values_fill = 0)
      for (factor_name in setdiff(mofa_v50_network_factors, colnames(profile))) profile[[factor_name]] <- 0
      profile_matrix <- as.matrix(profile[, mofa_v50_network_factors, drop = FALSE])
      rownames(profile_matrix) <- profile$feature_node_id
      module <- if (nrow(profile_matrix) <= 2) rep(1L, nrow(profile_matrix)) else stats::cutree(
        stats::hclust(stats::dist(profile_matrix), method = "ward.D2"),
        k = min(2L, max(1L, ceiling(nrow(profile_matrix) / 7)))
      )
      data.frame(feature_node_id = rownames(profile_matrix), module = as.integer(module), stringsAsFactors = FALSE)
    }
  )
)

mofa_v50_network_feature_nodes <- mofa_v50_network_edges %>%
  dplyr::group_by(feature_node_id, view, view_label, feature, feature_label) %>%
  dplyr::summarise(
    n_connected_factors = dplyr::n_distinct(factor),
    maximum_loading = max(loading_strength, na.rm = TRUE),
    minimum_response_p = if (any(is.finite(p_value))) min(p_value, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  dplyr::left_join(mofa_v50_module_membership, by = "feature_node_id") %>%
  dplyr::mutate(
    view_order = match(view, required_views),
    module_label = paste0(as.character(view_label), " M", module),
    shared_feature = n_connected_factors >= 2,
    response_feature = is.finite(minimum_response_p) & minimum_response_p < 0.05,
    feature_label_plotmath = mofa_v50_plotmath_feature_label(view, feature_label)
  ) %>%
  dplyr::arrange(view_order, module, dplyr::desc(n_connected_factors), dplyr::desc(maximum_loading))

mofa_v50_network_modules <- mofa_v50_network_feature_nodes %>%
  dplyr::distinct(view, view_label, view_order, module, module_label) %>%
  dplyr::arrange(view_order, module) %>%
  dplyr::mutate(
    module_index = dplyr::row_number(),
    module_angle = pi / 2 - 2 * pi * (module_index - 1) / dplyr::n(),
    module_x = 3.10 * cos(module_angle),
    module_y = 2.60 * sin(module_angle)
  )

mofa_v50_network_feature_nodes <- mofa_v50_network_feature_nodes %>%
  dplyr::left_join(mofa_v50_network_modules %>% dplyr::select(module_label, module_x, module_y), by = "module_label") %>%
  dplyr::group_by(module_label) %>%
  dplyr::mutate(
    feature_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.34,
    x = module_x,
    y = module_y + feature_offset,
    label_x = x + ifelse(module_x < -0.25, -0.72, 0.72),
    label_hjust = ifelse(module_x < -0.25, 1, 0)
  ) %>%
  dplyr::ungroup()

mofa_v50_network_modules <- mofa_v50_network_modules %>%
  dplyr::left_join(
    mofa_v50_network_feature_nodes %>%
      dplyr::group_by(module_label) %>%
      dplyr::summarise(
        xmin = min(x) - 0.38,
        xmax = max(x) + 0.38,
        ymin = min(y) - 0.26,
        ymax = max(y) + 0.26,
        label_y = max(y) + 0.58,
        .groups = "drop"
      ),
    by = "module_label"
  )

mofa_v50_network_factor_nodes <- data.frame(
  factor = mofa_v50_network_factors,
  factor_angle = pi / 2 - 2 * pi * (seq_along(mofa_v50_network_factors) - 1) / length(mofa_v50_network_factors),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    mofa_v50_network_selection_audit %>%
      dplyr::select(factor, strongest_response_effect, minimum_response_p, host_rnaseq_active),
    by = "factor"
  ) %>%
  dplyr::mutate(
    x = 7.00 * cos(factor_angle),
    y = 5.55 * sin(factor_angle),
    elbow_x = ifelse(x > 0, 4.65, -4.65),
    factor_label = paste0(
      factor,
      ifelse(host_rnaseq_active, "\nHost RNA-seq active", ""),
      ifelse(
        is.finite(minimum_response_p),
        paste0(
          "\npCR - non-pCR = ", sprintf("%+.2f", strongest_response_effect),
          "\np = ", ifelse(minimum_response_p < 0.001, "< 0.001", mofa_v50_format_p(minimum_response_p))
        ),
        ""
      )
    )
  )

mofa_v50_network_plot_edges <- mofa_v50_network_edges %>%
  dplyr::left_join(mofa_v50_network_factor_nodes %>% dplyr::select(factor, x_factor = x, y_factor = y, elbow_x), by = "factor") %>%
  dplyr::left_join(mofa_v50_network_feature_nodes %>% dplyr::select(feature_node_id, x_feature = x, y_feature = y), by = "feature_node_id") %>%
  dplyr::mutate(mid_x = elbow_x)

p_mofa_v50_multifactor_feature_network <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = mofa_v50_network_modules,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = view_label),
    alpha = 0.10,
    color = "grey72",
    linewidth = 0.28
  ) +
  ggplot2::geom_segment(
    data = mofa_v50_network_plot_edges,
    ggplot2::aes(x = x_factor, y = y_factor, xend = mid_x, yend = y_factor, color = loading_sign, linewidth = loading_strength),
    alpha = 0.42, lineend = "round"
  ) +
  ggplot2::geom_segment(
    data = mofa_v50_network_plot_edges,
    ggplot2::aes(x = mid_x, y = y_factor, xend = mid_x, yend = y_feature, color = loading_sign, linewidth = loading_strength),
    alpha = 0.42, lineend = "round"
  ) +
  ggplot2::geom_segment(
    data = mofa_v50_network_plot_edges,
    ggplot2::aes(x = mid_x, y = y_feature, xend = x_feature, yend = y_feature, color = loading_sign, linewidth = loading_strength),
    alpha = 0.42, lineend = "round"
  ) +
  ggplot2::geom_label(
    data = mofa_v50_network_factor_nodes,
    ggplot2::aes(x = x, y = y, label = factor_label),
    fill = "#FFF2B3", color = "grey10", label.size = 0.28,
    size = 2.45, lineheight = 0.90, fontface = "plain"
  ) +
  ggplot2::geom_point(
    data = mofa_v50_network_feature_nodes,
    ggplot2::aes(x = x, y = y, fill = view_label, size = maximum_loading),
    shape = 21, color = "grey20", stroke = 0.42
  ) +
  ggplot2::geom_point(
    data = mofa_v50_network_feature_nodes %>% dplyr::filter(shared_feature),
    ggplot2::aes(x = x, y = y),
    shape = 21, size = 3.9, fill = NA, color = "#7A3E9D", stroke = 0.76
  ) +
  ggplot2::geom_text(
    data = mofa_v50_network_feature_nodes %>% dplyr::filter(response_feature),
    ggplot2::aes(x = x, y = y, label = "*"),
    nudge_y = 0.25, size = 3.0, color = "black"
  ) +
  ggplot2::geom_text(
    data = mofa_v50_network_feature_nodes,
    ggplot2::aes(x = label_x, y = y, label = feature_label_plotmath, hjust = label_hjust),
    parse = TRUE, size = 2.06
  ) +
  ggplot2::geom_label(
    data = mofa_v50_network_modules,
    ggplot2::aes(x = module_x, y = label_y, label = module_label),
    fill = "white", color = "grey15", label.size = 0.18, size = 2.60,
    label.padding = grid::unit(0.10, "lines"), fontface = "plain"
  ) +
  ggplot2::scale_color_manual(values = c(Positive = "#D55E00", Negative = "#0072B2"), name = "Feature-weight sign") +
  ggplot2::scale_fill_manual(values = view_colors, name = "Omics view") +
  ggplot2::scale_linewidth_continuous(range = c(0.26, 0.95), guide = "none") +
  ggplot2::scale_size_continuous(range = c(2.1, 3.5), guide = "none") +
  ggplot2::coord_equal(xlim = c(-9.0, 9.0), ylim = c(-6.8, 6.8), clip = "off") +
  ggplot2::labs(
    title = "Multi-factor feature-weight network",
    subtitle = paste0(
      "Five QC-passing factors are selected with response evidence as the primary criterion, ",
      "while Host RNA-seq-active factors are preferentially retained when scores are comparable."
    ),
    caption = paste0(
      "Modules are Ward.D2 clusters of signed feature-weight profiles within each omics view. ",
      "Pipe-like orthogonal connectors are used for readability; purple outlines denote shared features."
    ),
    x = NULL, y = NULL
  ) +
  ggplot2::theme_void(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.7, color = "grey35"),
    plot.caption = ggplot2::element_text(size = 7.1, color = "grey35", hjust = 0),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(5, 220, 5, 220)
  )

if (FALSE) ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_multifactor_feature_loading_network.svg"),
  plot = p_mofa_v50_multifactor_feature_network,
  width = 13.8, height = 9.9, units = "in", device = svglite::svglite, bg = "white"
)

# All 15 Factor7-fixed triples using the requested candidates.
mofa_v50_anchor_factor <- if ("Factor7" %in% mofa_v50_factor_order) "Factor7" else mofa_v50_factor_order[1]
mofa_v50_triple_candidate_pool <- intersect(
  c("Factor1", "Factor4", "Factor2", "Factor8", "Factor9", "Factor10"),
  setdiff(mofa_v50_factor_order, mofa_v50_anchor_factor)
)
mofa_v50_factor_triples <- lapply(
  utils::combn(mofa_v50_triple_candidate_pool, 2, simplify = FALSE),
  function(partners) c(mofa_v50_anchor_factor, partners)
)

mofa_v50_triple_permanova_atlas <- dplyr::bind_rows(
  lapply(
    seq_along(mofa_v50_factor_triples),
    function(triple_index) {
      triple_factors <- mofa_v50_factor_triples[[triple_index]]
      triple_data <- mofa_v50_subject_mean_factor_data %>%
        dplyr::select(SubjectID, TRG_plot, dplyr::all_of(triple_factors)) %>%
        dplyr::filter(dplyr::if_all(dplyr::all_of(triple_factors), is.finite), !is.na(TRG_plot)) %>%
        dplyr::mutate(TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR")))
      score_matrix <- scale(as.matrix(triple_data[, triple_factors, drop = FALSE]))
      score_matrix[!is.finite(score_matrix)] <- 0
      distance_object <- stats::dist(score_matrix)

      set.seed(20261901 + triple_index)
      permanova_fit <- tryCatch(
        vegan::adonis2(distance_object ~ TRG_plot, data = triple_data, permutations = mofa_v50_outlier_sensitivity_permutations),
        error = function(e) NULL
      )
      dispersion_fit <- tryCatch(vegan::betadisper(distance_object, triple_data$TRG_plot), error = function(e) NULL)
      dispersion_test <- if (is.null(dispersion_fit)) NULL else tryCatch(
        vegan::permutest(dispersion_fit, permutations = mofa_v50_outlier_sensitivity_permutations),
        error = function(e) NULL
      )
      energy_test <- mofa_v50_energy_distance_test(score_matrix, triple_data$TRG_plot, permutations = 999L)

      host_partner_count <- sum(triple_factors[-1] %in% mofa_v50_host_active_factors)
      data.frame(
        triple_id = paste(triple_factors, collapse = "__"),
        factor_x = triple_factors[1], factor_y = triple_factors[2], factor_z = triple_factors[3],
        n_subjects = nrow(triple_data),
        permanova_f = if (is.null(permanova_fit)) NA_real_ else as.numeric(permanova_fit$F[1]),
        permanova_r2 = if (is.null(permanova_fit)) NA_real_ else as.numeric(permanova_fit$R2[1]),
        permanova_p = if (is.null(permanova_fit)) NA_real_ else as.numeric(permanova_fit$`Pr(>F)`[1]),
        dispersion_p = if (is.null(dispersion_test)) NA_real_ else as.numeric(dispersion_test$tab$`Pr(>F)`[1]),
        energy_statistic = energy_test$energy_statistic[1],
        energy_p = energy_test$energy_p[1],
        host_partner_count = host_partner_count,
        partner_has_host_view = host_partner_count >= 1,
        stringsAsFactors = FALSE
      )
    }
  )
) %>%
  dplyr::mutate(
    dispersion_ok = is.na(dispersion_p) | dispersion_p >= 0.05,
    selection_rank = dplyr::row_number()
  ) %>%
  dplyr::arrange(
    dplyr::desc(host_partner_count),
    dplyr::desc(dispersion_ok),
    energy_p,
    permanova_p,
    dplyr::desc(energy_statistic),
    dplyr::desc(permanova_r2)
  ) %>%
  dplyr::mutate(selection_rank = dplyr::row_number())

write.csv(
  mofa_v50_triple_permanova_atlas,
  file.path(mofa_v50_result_dir, "MOFA_v50_three_factor_PERMANOVA_atlas.csv"),
  row.names = FALSE
)

mofa_v50_selected_triple_diagnostics <- mofa_v50_triple_permanova_atlas %>%
  dplyr::slice_head(n = min(3, nrow(mofa_v50_triple_permanova_atlas)))

write.csv(
  mofa_v50_selected_triple_diagnostics,
  file.path(mofa_v50_result_dir, "MOFA_v50_selected_three_factor_diagnostics.csv"),
  row.names = FALSE
)

# Save all combinations in a multi-panel SVG and the top 3 separately.
if (FALSE) {
svglite::svglite(
  file.path(mofa_v50_figure_dir, "MOFA_v50_three_factor_response_3D_all_combinations.svg"),
  width = 14.0,
  height = 20.0,
  bg = "white"
)
graphics::par(mfrow = c(5, 3), mar = c(1.6, 1.8, 2.2, 0.6), oma = c(0.2, 0.2, 0.8, 0.2), xpd = NA)
for (i in seq_len(nrow(mofa_v50_triple_permanova_atlas))) {
  mofa_v50_draw_s3d_panel(
    c(
      mofa_v50_triple_permanova_atlas$factor_x[i],
      mofa_v50_triple_permanova_atlas$factor_y[i],
      mofa_v50_triple_permanova_atlas$factor_z[i]
    ),
    mofa_v50_triple_permanova_atlas[i, , drop = FALSE],
    point_cex = 0.80
  )
}
graphics::mtext("Factor7-fixed three-factor response maps for all requested combinations", side = 3, outer = TRUE, line = -0.2, font = 2)
grDevices::dev.off()
}


if (FALSE) {
svglite::svglite(
  file.path(mofa_v50_figure_dir, "MOFA_v50_three_factor_response_3D_top3.svg"),
  width = 15.0,
  height = 5.2,
  bg = "white"
)
graphics::par(mfrow = c(1, max(1, nrow(mofa_v50_selected_triple_diagnostics))), mar = c(1.6, 1.8, 2.3, 0.8), xpd = NA)
for (i in seq_len(nrow(mofa_v50_selected_triple_diagnostics))) {
  mofa_v50_draw_s3d_panel(
    c(
      mofa_v50_selected_triple_diagnostics$factor_x[i],
      mofa_v50_selected_triple_diagnostics$factor_y[i],
      mofa_v50_selected_triple_diagnostics$factor_z[i]
    ),
    mofa_v50_selected_triple_diagnostics[i, , drop = FALSE],
    point_cex = 0.95
  )
}
grDevices::dev.off()
}


# Rebuild the selected-triple network using the top-ranked triple, with extra space and pipe-like edges.
mofa_v50_triple_network_factors <- c(
  mofa_v50_selected_triple_diagnostics$factor_x[1],
  mofa_v50_selected_triple_diagnostics$factor_y[1],
  mofa_v50_selected_triple_diagnostics$factor_z[1]
)

mofa_v50_triple_active_factor_views <- mofa_variance_explained %>%
  dplyr::mutate(factor = as.character(factor), view = as.character(view)) %>%
  dplyr::filter(factor %in% mofa_v50_triple_network_factors, r2 >= mofa_active_view_r2) %>%
  dplyr::select(factor, view, view_r2 = r2)

mofa_v50_triple_network_seed_edges <- dplyr::bind_rows(
  mofa_feature_weights %>%
    dplyr::filter(factor %in% mofa_v50_triple_network_factors, display_eligible) %>%
    dplyr::inner_join(mofa_v50_triple_active_factor_views, by = c("factor", "view")) %>%
    dplyr::group_by(factor, view, direction) %>%
    dplyr::slice_max(order_by = abs(weight_within_view), n = mofa_network_features_per_view_direction, with_ties = FALSE) %>%
    dplyr::ungroup(),
  mofa_feature_weights %>%
    dplyr::filter(
      factor %in% mofa_v50_triple_network_factors,
      display_eligible,
      is.finite(p_value),
      p_value < 0.05
    ) %>%
    dplyr::inner_join(mofa_v50_triple_active_factor_views, by = c("factor", "view")) %>%
    dplyr::group_by(factor) %>%
    dplyr::arrange(p_value, dplyr::desc(abs(weight_within_view)), .by_group = TRUE) %>%
    dplyr::slice_head(n = 2) %>%
    dplyr::ungroup()
) %>%
  dplyr::distinct(factor, view, feature, .keep_all = TRUE) %>%
  dplyr::mutate(feature_node_id = paste(view, feature, sep = "::"), edge_id = paste(factor, feature_node_id, sep = "__"))

mofa_v50_triple_network_edges <- mofa_feature_weights %>%
  dplyr::filter(factor %in% mofa_v50_triple_network_factors, display_eligible) %>%
  dplyr::inner_join(mofa_v50_triple_active_factor_views, by = c("factor", "view")) %>%
  dplyr::mutate(feature_node_id = paste(view, feature, sep = "::"), edge_id = paste(factor, feature_node_id, sep = "__")) %>%
  dplyr::filter(
    feature_node_id %in% unique(mofa_v50_triple_network_seed_edges$feature_node_id),
    abs(weight_within_view) >= mofa_network_shared_loading_threshold | edge_id %in% mofa_v50_triple_network_seed_edges$edge_id
  ) %>%
  dplyr::mutate(loading_sign = ifelse(weight >= 0, "Positive", "Negative"), loading_strength = abs(weight_within_view))

mofa_v50_triple_module_membership <- dplyr::bind_rows(
  lapply(
    split(mofa_v50_triple_network_edges, mofa_v50_triple_network_edges$view),
    function(view_edge_data) {
      profile <- view_edge_data %>%
        dplyr::select(feature_node_id, factor, weight_within_view) %>%
        tidyr::pivot_wider(names_from = factor, values_from = weight_within_view, values_fill = 0)
      for (factor_name in setdiff(mofa_v50_triple_network_factors, colnames(profile))) profile[[factor_name]] <- 0
      profile_matrix <- as.matrix(profile[, mofa_v50_triple_network_factors, drop = FALSE])
      rownames(profile_matrix) <- profile$feature_node_id
      module <- if (nrow(profile_matrix) <= 2) rep(1L, nrow(profile_matrix)) else stats::cutree(
        stats::hclust(stats::dist(profile_matrix), method = "ward.D2"),
        k = min(2L, max(1L, ceiling(nrow(profile_matrix) / 7)))
      )
      data.frame(feature_node_id = rownames(profile_matrix), module = as.integer(module), stringsAsFactors = FALSE)
    }
  )
)

mofa_v50_triple_network_feature_nodes <- mofa_v50_triple_network_edges %>%
  dplyr::group_by(feature_node_id, view, view_label, feature, feature_label) %>%
  dplyr::summarise(
    n_connected_factors = dplyr::n_distinct(factor),
    maximum_loading = max(loading_strength, na.rm = TRUE),
    minimum_response_p = if (any(is.finite(p_value))) min(p_value, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  dplyr::left_join(mofa_v50_triple_module_membership, by = "feature_node_id") %>%
  dplyr::mutate(
    view_order = match(view, required_views),
    module_label = paste0(as.character(view_label), " M", module),
    shared_feature = n_connected_factors >= 2,
    response_feature = is.finite(minimum_response_p) & minimum_response_p < 0.05,
    feature_label_plotmath = mofa_v50_plotmath_feature_label(view, feature_label)
  ) %>%
  dplyr::arrange(view_order, module, dplyr::desc(n_connected_factors), dplyr::desc(maximum_loading))

mofa_v50_triple_network_modules <- mofa_v50_triple_network_feature_nodes %>%
  dplyr::distinct(view, view_label, view_order, module, module_label) %>%
  dplyr::arrange(view_order, module) %>%
  dplyr::mutate(
    module_index = dplyr::row_number(),
    module_angle = pi / 2 - 2 * pi * (module_index - 1) / dplyr::n(),
    module_x = 3.05 * cos(module_angle),
    module_y = 2.55 * sin(module_angle)
  )

mofa_v50_triple_network_feature_nodes <- mofa_v50_triple_network_feature_nodes %>%
  dplyr::left_join(mofa_v50_triple_network_modules %>% dplyr::select(module_label, module_x, module_y), by = "module_label") %>%
  dplyr::group_by(module_label) %>%
  dplyr::mutate(
    feature_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.33,
    x = module_x,
    y = module_y + feature_offset,
    label_x = x + ifelse(module_x < -0.25, -0.70, 0.70),
    label_hjust = ifelse(module_x < -0.25, 1, 0)
  ) %>%
  dplyr::ungroup()

mofa_v50_triple_network_modules <- mofa_v50_triple_network_modules %>%
  dplyr::left_join(
    mofa_v50_triple_network_feature_nodes %>%
      dplyr::group_by(module_label) %>%
      dplyr::summarise(
        xmin = min(x) - 0.36,
        xmax = max(x) + 0.36,
        ymin = min(y) - 0.24,
        ymax = max(y) + 0.24,
        label_y = max(y) + 0.56,
        .groups = "drop"
      ),
    by = "module_label"
  )

mofa_v50_triple_network_factor_nodes <- data.frame(
  factor = mofa_v50_triple_network_factors,
  factor_angle = pi / 2 - 2 * pi * (seq_along(mofa_v50_triple_network_factors) - 1) / length(mofa_v50_triple_network_factors),
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    mofa_factor_wilcoxon_tests %>%
      dplyr::filter(comparison == "overall_response_subject_mean", factor %in% mofa_v50_triple_network_factors) %>%
      dplyr::select(factor, strongest_response_effect = wilcoxon_effect, minimum_response_p = wilcoxon_p),
    by = "factor"
  ) %>%
  dplyr::mutate(
    host_rnaseq_active = factor %in% mofa_v50_host_active_factors,
    x = 6.85 * cos(factor_angle),
    y = 5.40 * sin(factor_angle),
    elbow_x = ifelse(x > 0, 4.50, -4.50),
    factor_label = paste0(
      factor,
      ifelse(host_rnaseq_active, "\nHost RNA-seq active", ""),
      ifelse(
        is.finite(minimum_response_p),
        paste0(
          "\npCR - non-pCR = ", sprintf("%+.2f", strongest_response_effect),
          "\np = ", ifelse(minimum_response_p < 0.001, "< 0.001", mofa_v50_format_p(minimum_response_p))
        ),
        ""
      )
    )
  )

mofa_v50_triple_network_plot_edges <- mofa_v50_triple_network_edges %>%
  dplyr::left_join(mofa_v50_triple_network_factor_nodes %>% dplyr::select(factor, x_factor = x, y_factor = y, elbow_x), by = "factor") %>%
  dplyr::left_join(mofa_v50_triple_network_feature_nodes %>% dplyr::select(feature_node_id, x_feature = x, y_feature = y), by = "feature_node_id") %>%
  dplyr::mutate(mid_x = elbow_x)

p_mofa_v50_selected_triple_feature_network <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = mofa_v50_triple_network_modules,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = view_label),
    alpha = 0.10, color = "grey72", linewidth = 0.28
  ) +
  ggplot2::geom_segment(
    data = mofa_v50_triple_network_plot_edges,
    ggplot2::aes(x = x_factor, y = y_factor, xend = mid_x, yend = y_factor, color = loading_sign, linewidth = loading_strength),
    alpha = 0.40, lineend = "round"
  ) +
  ggplot2::geom_segment(
    data = mofa_v50_triple_network_plot_edges,
    ggplot2::aes(x = mid_x, y = y_factor, xend = mid_x, yend = y_feature, color = loading_sign, linewidth = loading_strength),
    alpha = 0.40, lineend = "round"
  ) +
  ggplot2::geom_segment(
    data = mofa_v50_triple_network_plot_edges,
    ggplot2::aes(x = mid_x, y = y_feature, xend = x_feature, yend = y_feature, color = loading_sign, linewidth = loading_strength),
    alpha = 0.40, lineend = "round"
  ) +
  ggplot2::geom_label(
    data = mofa_v50_triple_network_factor_nodes,
    ggplot2::aes(x = x, y = y, label = factor_label),
    fill = "#FFF2B3", color = "grey10", label.size = 0.28,
    size = 2.55, lineheight = 0.90, fontface = "plain"
  ) +
  ggplot2::geom_point(
    data = mofa_v50_triple_network_feature_nodes,
    ggplot2::aes(x = x, y = y, fill = view_label, size = maximum_loading),
    shape = 21, color = "grey20", stroke = 0.42
  ) +
  ggplot2::geom_point(
    data = mofa_v50_triple_network_feature_nodes %>% dplyr::filter(shared_feature),
    ggplot2::aes(x = x, y = y),
    shape = 21, size = 3.7, fill = NA, color = "#7A3E9D", stroke = 0.72
  ) +
  ggplot2::geom_text(
    data = mofa_v50_triple_network_feature_nodes %>% dplyr::filter(response_feature),
    ggplot2::aes(x = x, y = y, label = "*"),
    nudge_y = 0.24, size = 3.0, color = "black"
  ) +
  ggplot2::geom_text(
    data = mofa_v50_triple_network_feature_nodes,
    ggplot2::aes(x = label_x, y = y, label = feature_label_plotmath, hjust = label_hjust),
    parse = TRUE, size = 2.10
  ) +
  ggplot2::geom_label(
    data = mofa_v50_triple_network_modules,
    ggplot2::aes(x = module_x, y = label_y, label = module_label),
    fill = "white", color = "grey15", label.size = 0.18, size = 2.60,
    label.padding = grid::unit(0.10, "lines"), fontface = "plain"
  ) +
  ggplot2::scale_color_manual(values = c(Positive = "#D55E00", Negative = "#0072B2"), name = "Feature-weight sign") +
  ggplot2::scale_fill_manual(values = view_colors, name = "Omics view") +
  ggplot2::scale_linewidth_continuous(range = c(0.26, 0.95), guide = "none") +
  ggplot2::scale_size_continuous(range = c(2.2, 3.7), guide = "none") +
  ggplot2::coord_equal(xlim = c(-8.8, 8.8), ylim = c(-6.6, 6.6), clip = "off") +
  ggplot2::labs(
    title = "Selected three-factor feature-weight network",
    subtitle = "The top-ranked Factor7-fixed combination is shown. Host RNA-seq-active factors are explicitly flagged when present.",
    caption = "Pipe-like orthogonal connectors are used for readability; purple outlines denote features shared by multiple factors.",
    x = NULL, y = NULL
  ) +
  ggplot2::theme_void(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.6, color = "grey35"),
    plot.caption = ggplot2::element_text(size = 7.1, color = "grey35", hjust = 0),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(5, 210, 5, 210)
  )

if (FALSE) ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_selected_triple_multifactor_feature_loading_network.svg"),
  plot = p_mofa_v50_selected_triple_feature_network,
  width = 13.2, height = 9.4, units = "in", device = svglite::svglite, bg = "white"
)


#-----------------------------------------------------------------#
# 22.9 v50 final display refinement: violin, network, and 3D plots
#-----------------------------------------------------------------#

# Remove all violin-width variant outputs from previous versions.
unlink(Sys.glob(file.path(mofa_v50_figure_dir, "*_w60.svg")))
unlink(Sys.glob(file.path(mofa_v50_figure_dir, "*_w70.svg")))
unlink(Sys.glob(file.path(mofa_v50_figure_dir, "*_w80.svg")))
unlink(Sys.glob(file.path(mofa_v50_figure_dir, "*_w90.svg")))

mofa_v50_save_violin_variants <- function(...) invisible(NULL)

mofa_v50_make_p_only_annotation <- function(plot_data, comparison_name, value_column = "factor_score") {
  y_value <- rlang::sym(value_column)
  plot_data %>%
    dplyr::group_by(factor) %>%
    dplyr::summarise(
      y_min = min(!!y_value, na.rm = TRUE),
      y_max = max(!!y_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      y_position = y_max + 0.18 * pmax(y_max - y_min, 1)
    ) %>%
    dplyr::left_join(
      mofa_factor_wilcoxon_tests %>%
        dplyr::filter(comparison == comparison_name) %>%
        dplyr::transmute(
          factor = factor(as.character(factor), levels = mofa_v50_factor_order),
          wilcoxon_p
        ),
      by = "factor"
    ) %>%
    dplyr::mutate(
      label = paste0(
        "p = ",
        ifelse(
          !is.finite(wilcoxon_p),
          "NA",
          ifelse(wilcoxon_p < 0.001, "< 0.001", mofa_v50_format_p(wilcoxon_p))
        )
      )
    )
}

mofa_v50_build_response_plot <- function(plot_data, annotation_data, title_text) {
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = TRG_plot, y = factor_score, fill = TRG_plot)
  ) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.68, color = "grey30", linewidth = 0.38) +
    ggplot2::geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      fill = NA,
      coef = 0,
      staplewidth = 0,
      linewidth = 0.42
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.08, height = 0),
      shape = 21,
      size = 1.45,
      color = "grey20",
      stroke = 0.32,
      alpha = 0.85
    ) +
    ggplot2::geom_segment(
      data = annotation_data,
      ggplot2::aes(x = 1, xend = 2, y = y_position, yend = y_position),
      inherit.aes = FALSE,
      linewidth = 0.36,
      color = "grey25"
    ) +
    ggplot2::geom_segment(
      data = annotation_data,
      ggplot2::aes(x = 1, xend = 1, y = y_position, yend = y_position - 0.03 * pmax(y_max - y_min, 1)),
      inherit.aes = FALSE,
      linewidth = 0.36,
      color = "grey25"
    ) +
    ggplot2::geom_segment(
      data = annotation_data,
      ggplot2::aes(x = 2, xend = 2, y = y_position, yend = y_position - 0.03 * pmax(y_max - y_min, 1)),
      inherit.aes = FALSE,
      linewidth = 0.36,
      color = "grey25"
    ) +
    ggplot2::geom_text(
      data = annotation_data,
      ggplot2::aes(x = 1.5, y = y_position + 0.06 * pmax(y_max - y_min, 1), label = label),
      inherit.aes = FALSE,
      size = 2.25,
      lineheight = 0.90,
      color = "grey20"
    ) +
    ggplot2::scale_fill_manual(values = mofa_group_colors[c("pCR", "non_pCR")], guide = "none", drop = FALSE) +
    ggplot2::facet_wrap(~ factor, nrow = 2, scales = "free_y") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.36))) +
    ggplot2::labs(title = title_text, x = NULL, y = "MOFA factor score") +
    ggplot2::theme_classic(base_size = 8.8) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "plain", size = 8.0),
      axis.text.x = ggplot2::element_text(size = 9.0, face = "plain"),
      axis.text.y = ggplot2::element_text(face = "plain"),
      panel.spacing = grid::unit(0.75, "lines")
    )
}

mofa_v50_overall_annotation <- mofa_v50_make_p_only_annotation(
  mofa_v50_overall_response_long,
  "overall_response_subject_mean"
)
mofa_v50_baseline_annotation <- mofa_v50_make_p_only_annotation(
  mofa_v50_baseline_response_long,
  "baseline_response"
)
mofa_v50_after_annotation <- mofa_v50_make_p_only_annotation(
  mofa_v50_after_response_long,
  "ongoing_response"
)

p_mofa_v50_overall_response <- mofa_v50_build_response_plot(
  mofa_v50_overall_response_long,
  mofa_v50_overall_annotation,
  "pCR versus non-pCR: subject mean across available timepoints"
)
ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_overall_subject_mean.svg"),
  plot = p_mofa_v50_overall_response,
  width = 7.15,
  height = max(5.4, 2.25 * ceiling(length(mofa_v50_factor_order) / 5)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

p_mofa_v50_baseline_response <- mofa_v50_build_response_plot(
  mofa_v50_baseline_response_long,
  mofa_v50_baseline_annotation,
  "pCR versus non-pCR at Baseline"
)
ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_baseline.svg"),
  plot = p_mofa_v50_baseline_response,
  width = 7.15,
  height = max(5.4, 2.25 * ceiling(length(mofa_v50_factor_order) / 5)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

p_mofa_v50_after_response <- mofa_v50_build_response_plot(
  mofa_v50_after_response_long,
  mofa_v50_after_annotation,
  "pCR versus non-pCR after radiotherapy"
)
ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_after_RT.svg"),
  plot = p_mofa_v50_after_response,
  width = 7.15,
  height = max(5.4, 2.25 * ceiling(length(mofa_v50_factor_order) / 5)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

mofa_v50_delta_annotation <- mofa_v50_make_p_only_annotation(
  mofa_v50_paired_factor_data,
  "paired_differential_change",
  value_column = "factor_delta"
)

p_mofa_v50_paired_change_by_response <- ggplot2::ggplot(
  mofa_v50_paired_factor_data,
  ggplot2::aes(x = TRG_plot, y = factor_delta, fill = TRG_plot)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.40) +
  ggplot2::geom_violin(trim = FALSE, alpha = 0.68, color = "grey30", linewidth = 0.38) +
  ggplot2::geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    fill = NA,
    coef = 0,
    staplewidth = 0,
    linewidth = 0.42
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.08, height = 0),
    shape = 21,
    size = 1.55,
    color = "grey20",
    stroke = 0.32
  ) +
  ggplot2::geom_segment(
    data = mofa_v50_delta_annotation,
    ggplot2::aes(x = 1, xend = 2, y = y_position, yend = y_position),
    inherit.aes = FALSE,
    linewidth = 0.36,
    color = "grey25"
  ) +
  ggplot2::geom_text(
    data = mofa_v50_delta_annotation,
    ggplot2::aes(x = 1.5, y = y_position + 0.06 * pmax(y_max - y_min, 1), label = label),
    inherit.aes = FALSE,
    size = 2.10,
    lineheight = 0.90
  ) +
  ggplot2::scale_fill_manual(values = mofa_group_colors[c("pCR", "non_pCR")], guide = "none", drop = FALSE) +
  ggplot2::facet_wrap(~ factor, nrow = 2, scales = "free_y") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.36))) +
  ggplot2::labs(
    title = "Response-specific paired treatment change",
    subtitle = "Delta is After RT minus Baseline for each subject.",
    x = NULL,
    y = "Within-subject factor delta"
  ) +
  ggplot2::theme_classic(base_size = 8.8) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.7, color = "grey35"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "plain", size = 8.0),
    axis.text.x = ggplot2::element_text(size = 9.0, face = "plain"),
    axis.text.y = ggplot2::element_text(face = "plain"),
    panel.spacing = grid::unit(0.75, "lines")
  )

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_paired_change_by_response.svg"),
  plot = p_mofa_v50_paired_change_by_response,
  width = 7.15,
  height = max(5.4, 2.30 * ceiling(length(mofa_v50_factor_order) / 5)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

mofa_v50_selected_baseline_annotation <- mofa_v50_baseline_annotation %>%
  dplyr::filter(as.character(factor) == mofa_v50_selected_factor)

p_mofa_v50_selected_factor_baseline <- ggplot2::ggplot(
  mofa_v50_baseline_response_long %>% dplyr::filter(as.character(factor) == mofa_v50_selected_factor),
  ggplot2::aes(x = TRG_plot, y = factor_score, fill = TRG_plot)
) +
  ggplot2::geom_violin(trim = FALSE, alpha = 0.68, color = "grey30") +
  ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, fill = NA, coef = 0, staplewidth = 0, linewidth = 0.42) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08), shape = 21, size = 1.8, color = "grey20") +
  ggplot2::geom_segment(
    data = mofa_v50_selected_baseline_annotation,
    ggplot2::aes(x = 1, xend = 2, y = y_position, yend = y_position),
    inherit.aes = FALSE,
    linewidth = 0.38,
    color = "grey25"
  ) +
  ggplot2::geom_text(
    data = mofa_v50_selected_baseline_annotation,
    ggplot2::aes(x = 1.5, y = y_position + 0.06 * pmax(y_max - y_min, 1), label = label),
    inherit.aes = FALSE,
    size = 2.45,
    lineheight = 0.90
  ) +
  ggplot2::scale_fill_manual(values = mofa_group_colors[c("pCR", "non_pCR")], guide = "none") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.36))) +
  ggplot2::labs(title = paste0(mofa_v50_selected_factor, ": Baseline response contrast"), x = NULL, y = "Factor score") +
  ggplot2::theme_classic(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(size = 10.0, face = "plain"),
    axis.text.y = ggplot2::element_text(face = "plain")
  )

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_selected_factor_baseline.svg"),
  plot = p_mofa_v50_selected_factor_baseline,
  width = 3.0,
  height = 3.5,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

# Rebuild 5-factor network with simpler module labels, less overlap, and no statistics in factor labels.
mofa_v50_network_modules <- mofa_v50_network_feature_nodes %>%
  dplyr::distinct(view, view_label, view_order, module) %>%
  dplyr::arrange(view_order, module) %>%
  dplyr::group_by(view) %>%
  dplyr::mutate(module_short = paste0("M", dplyr::row_number())) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    module_index = dplyr::row_number(),
    module_angle = pi / 2 - 2 * pi * (module_index - 1) / dplyr::n(),
    module_x = 3.35 * cos(module_angle),
    module_y = 2.85 * sin(module_angle)
  )

mofa_v50_network_feature_nodes <- mofa_v50_network_feature_nodes %>%
  dplyr::select(-dplyr::any_of(c("module_x", "module_y", "x", "y", "label_x", "label_hjust", "feature_offset"))) %>%
  dplyr::left_join(mofa_v50_network_modules %>% dplyr::select(view, module, module_short, module_x, module_y), by = c("view", "module")) %>%
  dplyr::group_by(view, module) %>%
  dplyr::mutate(
    feature_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.40,
    x = module_x,
    y = module_y + feature_offset,
    label_x = x + ifelse(module_x < -0.25, -0.82, 0.82),
    label_hjust = ifelse(module_x < -0.25, 1, 0)
  ) %>%
  dplyr::ungroup()

mofa_v50_network_modules <- mofa_v50_network_modules %>%
  dplyr::left_join(
    mofa_v50_network_feature_nodes %>%
      dplyr::group_by(view, module, module_short) %>%
      dplyr::summarise(
        xmin = min(x) - 0.40,
        xmax = max(x) + 0.40,
        ymin = min(y) - 0.28,
        ymax = max(y) + 0.28,
        label_y = max(y) + 0.60,
        .groups = "drop"
      ),
    by = c("view", "module", "module_short")
  )

mofa_v50_network_factor_nodes <- data.frame(
  factor = mofa_v50_network_factors,
  factor_angle = pi / 2 - 2 * pi * (seq_along(mofa_v50_network_factors) - 1) / length(mofa_v50_network_factors),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    factor_rank = dplyr::row_number(),
    x = 7.30 * cos(factor_angle),
    y = 5.90 * sin(factor_angle),
    elbow_x = ifelse(x > 0, 5.30 + 0.18 * (factor_rank - 3), -5.30 - 0.18 * (factor_rank - 3)),
    factor_label = factor
  )

mofa_v50_network_plot_edges <- mofa_v50_network_edges %>%
  dplyr::left_join(mofa_v50_network_factor_nodes %>% dplyr::select(factor, x_factor = x, y_factor = y, elbow_x, factor_rank), by = "factor") %>%
  dplyr::left_join(mofa_v50_network_feature_nodes %>% dplyr::select(feature_node_id, x_feature = x, y_feature = y), by = "feature_node_id") %>%
  dplyr::group_by(factor) %>%
  dplyr::mutate(edge_rank = dplyr::row_number(), edge_offset = (edge_rank - (dplyr::n() + 1) / 2) * 0.03, mid_x = elbow_x + edge_offset) %>%
  dplyr::ungroup()

p_mofa_v50_multifactor_feature_network <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = mofa_v50_network_modules,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = view_label),
    alpha = 0.10,
    color = "grey72",
    linewidth = 0.28
  ) +
  ggplot2::geom_segment(
    data = mofa_v50_network_plot_edges,
    ggplot2::aes(x = x_factor, y = y_factor, xend = mid_x, yend = y_factor, color = loading_sign, linewidth = loading_strength),
    alpha = 0.42,
    lineend = "round"
  ) +
  ggplot2::geom_segment(
    data = mofa_v50_network_plot_edges,
    ggplot2::aes(x = mid_x, y = y_factor, xend = mid_x, yend = y_feature, color = loading_sign, linewidth = loading_strength),
    alpha = 0.42,
    lineend = "round"
  ) +
  ggplot2::geom_segment(
    data = mofa_v50_network_plot_edges,
    ggplot2::aes(x = mid_x, y = y_feature, xend = x_feature, yend = y_feature, color = loading_sign, linewidth = loading_strength),
    alpha = 0.42,
    lineend = "round"
  ) +
  ggplot2::geom_label(
    data = mofa_v50_network_factor_nodes,
    ggplot2::aes(x = x, y = y, label = factor_label),
    fill = "#FFF2B3",
    color = "grey10",
    label.size = 0.28,
    size = 2.65,
    lineheight = 0.90,
    fontface = "plain"
  ) +
  ggplot2::geom_point(
    data = mofa_v50_network_feature_nodes,
    ggplot2::aes(x = x, y = y, fill = view_label, size = maximum_loading),
    shape = 21,
    color = "grey20",
    stroke = 0.42
  ) +
  ggplot2::geom_point(
    data = mofa_v50_network_feature_nodes %>% dplyr::filter(shared_feature),
    ggplot2::aes(x = x, y = y),
    shape = 21,
    size = 3.9,
    fill = NA,
    color = "#7A3E9D",
    stroke = 0.76
  ) +
  ggplot2::geom_text(
    data = mofa_v50_network_feature_nodes %>% dplyr::filter(response_feature),
    ggplot2::aes(x = x, y = y, label = "*"),
    nudge_y = 0.26,
    size = 3.0,
    color = "black"
  ) +
  ggplot2::geom_text(
    data = mofa_v50_network_feature_nodes,
    ggplot2::aes(x = label_x, y = y, label = feature_label_plotmath, hjust = label_hjust),
    parse = TRUE,
    size = 1.95
  ) +
  ggplot2::geom_label(
    data = mofa_v50_network_modules,
    ggplot2::aes(x = module_x, y = label_y, label = module_short, fill = view_label),
    color = "grey15",
    label.size = 0.16,
    size = 2.45,
    label.padding = grid::unit(0.09, "lines"),
    fontface = "plain",
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(values = c(Positive = "#D55E00", Negative = "#0072B2"), name = "Feature-weight sign") +
  ggplot2::scale_fill_manual(
    values = view_colors,
    name = "Feature group",
    guide = ggplot2::guide_legend(override.aes = list(shape = 21, size = 4, color = "grey20", alpha = 1))
  ) +
  ggplot2::scale_linewidth_continuous(range = c(0.24, 0.88), guide = "none") +
  ggplot2::scale_size_continuous(range = c(2.0, 3.4), guide = "none") +
  ggplot2::coord_equal(xlim = c(-9.5, 9.5), ylim = c(-7.2, 7.2), clip = "off") +
  ggplot2::labs(
    title = "Multi-factor feature-weight network",
    subtitle = "Modules are shown as M1/M2 labels only; feature groups are identified in the legend.",
    caption = "Pipe-like orthogonal connectors are slightly offset to reduce overlap; purple outlines denote shared features.",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_void(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.7, color = "grey35"),
    plot.caption = ggplot2::element_text(size = 7.1, color = "grey35", hjust = 0),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(5, 235, 5, 235)
  )

if (FALSE) ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_multifactor_feature_loading_network.svg"),
  plot = p_mofa_v50_multifactor_feature_network,
  width = 14.4,
  height = 10.3,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

# Rebuild top-ranked 3-factor network similarly.
mofa_v50_triple_network_factors <- c(
  mofa_v50_selected_triple_diagnostics$factor_x[1],
  mofa_v50_selected_triple_diagnostics$factor_y[1],
  mofa_v50_selected_triple_diagnostics$factor_z[1]
)

mofa_v50_triple_network_modules <- mofa_v50_triple_network_feature_nodes %>%
  dplyr::distinct(view, view_label, view_order, module) %>%
  dplyr::arrange(view_order, module) %>%
  dplyr::group_by(view) %>%
  dplyr::mutate(module_short = paste0("M", dplyr::row_number())) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    module_index = dplyr::row_number(),
    module_angle = pi / 2 - 2 * pi * (module_index - 1) / dplyr::n(),
    module_x = 3.15 * cos(module_angle),
    module_y = 2.70 * sin(module_angle)
  )

mofa_v50_triple_network_feature_nodes <- mofa_v50_triple_network_feature_nodes %>%
  dplyr::select(-dplyr::any_of(c("module_x", "module_y", "x", "y", "label_x", "label_hjust", "feature_offset"))) %>%
  dplyr::left_join(mofa_v50_triple_network_modules %>% dplyr::select(view, module, module_short, module_x, module_y), by = c("view", "module")) %>%
  dplyr::group_by(view, module) %>%
  dplyr::mutate(
    feature_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.38,
    x = module_x,
    y = module_y + feature_offset,
    label_x = x + ifelse(module_x < -0.25, -0.78, 0.78),
    label_hjust = ifelse(module_x < -0.25, 1, 0)
  ) %>%
  dplyr::ungroup()

mofa_v50_triple_network_modules <- mofa_v50_triple_network_modules %>%
  dplyr::left_join(
    mofa_v50_triple_network_feature_nodes %>%
      dplyr::group_by(view, module, module_short) %>%
      dplyr::summarise(
        xmin = min(x) - 0.38,
        xmax = max(x) + 0.38,
        ymin = min(y) - 0.26,
        ymax = max(y) + 0.26,
        label_y = max(y) + 0.58,
        .groups = "drop"
      ),
    by = c("view", "module", "module_short")
  )

mofa_v50_triple_network_factor_nodes <- data.frame(
  factor = mofa_v50_triple_network_factors,
  factor_angle = pi / 2 - 2 * pi * (seq_along(mofa_v50_triple_network_factors) - 1) / length(mofa_v50_triple_network_factors),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    factor_rank = dplyr::row_number(),
    x = 7.00 * cos(factor_angle),
    y = 5.60 * sin(factor_angle),
    elbow_x = ifelse(x > 0, 5.10 + 0.16 * (factor_rank - 2), -5.10 - 0.16 * (factor_rank - 2)),
    factor_label = factor
  )

mofa_v50_triple_network_plot_edges <- mofa_v50_triple_network_edges %>%
  dplyr::left_join(mofa_v50_triple_network_factor_nodes %>% dplyr::select(factor, x_factor = x, y_factor = y, elbow_x), by = "factor") %>%
  dplyr::left_join(mofa_v50_triple_network_feature_nodes %>% dplyr::select(feature_node_id, x_feature = x, y_feature = y), by = "feature_node_id") %>%
  dplyr::group_by(factor) %>%
  dplyr::mutate(edge_rank = dplyr::row_number(), edge_offset = (edge_rank - (dplyr::n() + 1) / 2) * 0.035, mid_x = elbow_x + edge_offset) %>%
  dplyr::ungroup()

p_mofa_v50_selected_triple_feature_network <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = mofa_v50_triple_network_modules,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = view_label),
    alpha = 0.10, color = "grey72", linewidth = 0.28
  ) +
  ggplot2::geom_segment(
    data = mofa_v50_triple_network_plot_edges,
    ggplot2::aes(x = x_factor, y = y_factor, xend = mid_x, yend = y_factor, color = loading_sign, linewidth = loading_strength),
    alpha = 0.40, lineend = "round"
  ) +
  ggplot2::geom_segment(
    data = mofa_v50_triple_network_plot_edges,
    ggplot2::aes(x = mid_x, y = y_factor, xend = mid_x, yend = y_feature, color = loading_sign, linewidth = loading_strength),
    alpha = 0.40, lineend = "round"
  ) +
  ggplot2::geom_segment(
    data = mofa_v50_triple_network_plot_edges,
    ggplot2::aes(x = mid_x, y = y_feature, xend = x_feature, yend = y_feature, color = loading_sign, linewidth = loading_strength),
    alpha = 0.40, lineend = "round"
  ) +
  ggplot2::geom_label(
    data = mofa_v50_triple_network_factor_nodes,
    ggplot2::aes(x = x, y = y, label = factor_label),
    fill = "#FFF2B3", color = "grey10", label.size = 0.28,
    size = 2.70, lineheight = 0.90, fontface = "plain"
  ) +
  ggplot2::geom_point(
    data = mofa_v50_triple_network_feature_nodes,
    ggplot2::aes(x = x, y = y, fill = view_label, size = maximum_loading),
    shape = 21, color = "grey20", stroke = 0.42
  ) +
  ggplot2::geom_point(
    data = mofa_v50_triple_network_feature_nodes %>% dplyr::filter(shared_feature),
    ggplot2::aes(x = x, y = y),
    shape = 21, size = 3.7, fill = NA, color = "#7A3E9D", stroke = 0.72
  ) +
  ggplot2::geom_text(
    data = mofa_v50_triple_network_feature_nodes %>% dplyr::filter(response_feature),
    ggplot2::aes(x = x, y = y, label = "*"),
    nudge_y = 0.24, size = 3.0, color = "black"
  ) +
  ggplot2::geom_text(
    data = mofa_v50_triple_network_feature_nodes,
    ggplot2::aes(x = label_x, y = y, label = feature_label_plotmath, hjust = label_hjust),
    parse = TRUE, size = 2.00
  ) +
  ggplot2::geom_label(
    data = mofa_v50_triple_network_modules,
    ggplot2::aes(x = module_x, y = label_y, label = module_short, fill = view_label),
    color = "grey15", label.size = 0.16, size = 2.45,
    label.padding = grid::unit(0.09, "lines"), fontface = "plain", show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(values = c(Positive = "#D55E00", Negative = "#0072B2"), name = "Feature-weight sign") +
  ggplot2::scale_fill_manual(
    values = view_colors,
    name = "Feature group",
    guide = ggplot2::guide_legend(override.aes = list(shape = 21, size = 4, color = "grey20", alpha = 1))
  ) +
  ggplot2::scale_linewidth_continuous(range = c(0.24, 0.88), guide = "none") +
  ggplot2::scale_size_continuous(range = c(2.1, 3.5), guide = "none") +
  ggplot2::coord_equal(xlim = c(-9.2, 9.2), ylim = c(-7.0, 7.0), clip = "off") +
  ggplot2::labs(
    title = "Selected three-factor feature-weight network",
    subtitle = "The top-ranked Factor7-fixed combination is shown with simplified module labels and offset pipe-like connectors.",
    caption = "Feature groups are identified by the legend only; purple outlines denote features shared by multiple factors.",
    x = NULL, y = NULL
  ) +
  ggplot2::theme_void(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.6, color = "grey35"),
    plot.caption = ggplot2::element_text(size = 7.1, color = "grey35", hjust = 0),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(5, 220, 5, 220)
  )

if (FALSE) ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_selected_triple_multifactor_feature_loading_network.svg"),
  plot = p_mofa_v50_selected_triple_feature_network,
  width = 13.8,
  height = 9.8,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

# Refine the 3D panels: larger fixed point size and non-overlapping legend placement.
mofa_v50_draw_s3d_panel <- function(triple_factors, diagnostic_row, point_cex = 1.18) {
  if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
    plot.new()
    text(0.5, 0.6, "Package 'scatterplot3d' is required for v50 3D plots.")
    text(0.5, 0.48, paste(triple_factors, collapse = " + "), cex = 0.9)
    return(invisible(NULL))
  }
  plot_data <- mofa_v50_subject_mean_factor_data %>%
    dplyr::select(SubjectID, TRG_plot, dplyr::all_of(triple_factors)) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(triple_factors), is.finite), !is.na(TRG_plot)) %>%
    dplyr::mutate(TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR")))

  score_matrix <- scale(as.matrix(plot_data[, triple_factors, drop = FALSE]))
  score_matrix[!is.finite(score_matrix)] <- 0
  plot_data$x <- score_matrix[, 1]
  plot_data$y <- score_matrix[, 2]
  plot_data$z <- score_matrix[, 3]
  plot_data$color <- unname(mofa_response_colors[as.character(plot_data$TRG_plot)])
  plot_data$line_color <- grDevices::adjustcolor(plot_data$color, alpha.f = 0.22)

  s3d <- scatterplot3d::scatterplot3d(
    x = plot_data$x,
    y = plot_data$y,
    z = plot_data$z,
    pch = 16,
    color = plot_data$color,
    cex.symbols = point_cex,
    type = "h",
    lty.hplot = 1,
    mar = c(2.0, 2.0, 3.0, 4.8),
    main = paste0("Rank ", diagnostic_row$selection_rank[1], ": ", paste(triple_factors, collapse = " + ")),
    xlab = triple_factors[1],
    ylab = triple_factors[2],
    zlab = triple_factors[3],
    angle = 52,
    scale.y = 1.0,
    box = TRUE,
    grid = TRUE
  )

  z_floor <- min(plot_data$z, na.rm = TRUE)
  for (i in seq_len(nrow(plot_data))) {
    xy_point <- s3d$xyz.convert(plot_data$x[i], plot_data$y[i], plot_data$z[i])
    xy_floor <- s3d$xyz.convert(plot_data$x[i], plot_data$y[i], z_floor)
    graphics::segments(xy_floor$x, xy_floor$y, xy_point$x, xy_point$y, col = plot_data$line_color[i], lwd = 0.7)
  }
  point_xy <- s3d$xyz.convert(plot_data$x, plot_data$y, plot_data$z)
  graphics::points(point_xy$x, point_xy$y, pch = 21, bg = plot_data$color, col = "grey20", cex = point_cex)

  stats_label <- paste0(
    "PERMANOVA p = ", ifelse(is.na(diagnostic_row$permanova_p[1]), "NA", ifelse(diagnostic_row$permanova_p[1] < 0.001, "<0.001", sprintf("%.3f", diagnostic_row$permanova_p[1]))),
    " | R2 = ", ifelse(is.na(diagnostic_row$permanova_r2[1]), "NA", sprintf("%.2f", diagnostic_row$permanova_r2[1])),
    "\nEnergy p = ", ifelse(is.na(diagnostic_row$energy_p[1]), "NA", ifelse(diagnostic_row$energy_p[1] < 0.001, "<0.001", sprintf("%.3f", diagnostic_row$energy_p[1]))),
    " | E = ", ifelse(is.na(diagnostic_row$energy_statistic[1]), "NA", sprintf("%.2f", diagnostic_row$energy_statistic[1]))
  )
  graphics::mtext(stats_label, side = 3, line = 0.2, adj = 0.02, cex = 0.72)
  graphics::legend(
    x = grconvertX(1.02, from = "npc", to = "user"),
    y = grconvertY(0.95, from = "npc", to = "user"),
    legend = c("pCR", "non-pCR"),
    pt.bg = unname(mofa_response_colors[c("pCR", "non_pCR")]),
    pch = 21,
    col = "grey20",
    pt.cex = 1.2,
    cex = 0.78,
    bty = "n",
    xpd = NA,
    xjust = 0,
    yjust = 1
  )
  invisible(NULL)
}

svglite::svglite(
  file.path(mofa_v50_figure_dir, "MOFA_v50_three_factor_response_3D_all_combinations.svg"),
  width = 15.0,
  height = 20.4,
  bg = "white"
)
graphics::par(mfrow = c(5, 3), mar = c(1.8, 1.8, 2.9, 3.6), oma = c(0.2, 0.2, 0.8, 0.2), xpd = NA)
for (i in seq_len(nrow(mofa_v50_triple_permanova_atlas))) {
  mofa_v50_draw_s3d_panel(
    c(mofa_v50_triple_permanova_atlas$factor_x[i], mofa_v50_triple_permanova_atlas$factor_y[i], mofa_v50_triple_permanova_atlas$factor_z[i]),
    mofa_v50_triple_permanova_atlas[i, , drop = FALSE],
    point_cex = 1.05
  )
}
graphics::mtext("Factor7-fixed three-factor response maps for all requested combinations", side = 3, outer = TRUE, line = -0.2, font = 2)
grDevices::dev.off()

if (FALSE) {
svglite::svglite(
  file.path(mofa_v50_figure_dir, "MOFA_v50_three_factor_response_3D_top3.svg"),
  width = 16.2,
  height = 5.6,
  bg = "white"
)
graphics::par(mfrow = c(1, max(1, nrow(mofa_v50_selected_triple_diagnostics))), mar = c(1.8, 1.8, 2.9, 4.0), xpd = NA)
for (i in seq_len(nrow(mofa_v50_selected_triple_diagnostics))) {
  mofa_v50_draw_s3d_panel(
    c(mofa_v50_selected_triple_diagnostics$factor_x[i], mofa_v50_selected_triple_diagnostics$factor_y[i], mofa_v50_selected_triple_diagnostics$factor_z[i]),
    mofa_v50_selected_triple_diagnostics[i, , drop = FALSE],
    point_cex = 1.28
  )
}
grDevices::dev.off()
}



#-----------------------------------------------------------------#
# 22.9 v50 network expansion and de-overlap refinement
#-----------------------------------------------------------------#

if (FALSE) {
  # Obsolete inherited network-rendering block disabled in v50.

mofa_v50_build_network_plot <- function(
  factor_set,
  plot_title,
  plot_subtitle = NULL,
  figure_width = 14.0,
  figure_height = 10.0
) {
  factor_set <- unique(factor_set)
  factor_set <- factor_set[factor_set %in% mofa_v50_factor_order]
  if (length(factor_set) < 2) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Not enough factors to build the network.")
    )
  }

  active_factor_views <- mofa_variance_explained %>%
    dplyr::mutate(factor = as.character(factor), view = as.character(view)) %>%
    dplyr::filter(factor %in% factor_set, r2 >= mofa_active_view_r2) %>%
    dplyr::select(factor, view, view_r2 = r2)

  seed_edges <- dplyr::bind_rows(
    mofa_feature_weights %>%
      dplyr::filter(factor %in% factor_set, display_eligible) %>%
      dplyr::inner_join(active_factor_views, by = c("factor", "view")) %>%
      dplyr::group_by(factor, view, direction) %>%
      dplyr::slice_max(order_by = abs(weight_within_view), n = mofa_network_features_per_view_direction, with_ties = FALSE) %>%
      dplyr::ungroup(),
    mofa_feature_weights %>%
      dplyr::filter(factor %in% factor_set, display_eligible, is.finite(p_value), p_value < 0.05) %>%
      dplyr::inner_join(active_factor_views, by = c("factor", "view")) %>%
      dplyr::group_by(factor) %>%
      dplyr::arrange(p_value, dplyr::desc(abs(weight_within_view)), .by_group = TRUE) %>%
      dplyr::slice_head(n = 2) %>%
      dplyr::ungroup()
  ) %>%
    dplyr::distinct(factor, view, feature, .keep_all = TRUE) %>%
    dplyr::mutate(
      feature_node_id = paste(view, feature, sep = "::"),
      edge_id = paste(factor, feature_node_id, sep = "__")
    )

  network_edges <- mofa_feature_weights %>%
    dplyr::filter(factor %in% factor_set, display_eligible) %>%
    dplyr::inner_join(active_factor_views, by = c("factor", "view")) %>%
    dplyr::mutate(
      feature_node_id = paste(view, feature, sep = "::"),
      edge_id = paste(factor, feature_node_id, sep = "__")
    ) %>%
    dplyr::filter(
      feature_node_id %in% unique(seed_edges$feature_node_id),
      abs(weight_within_view) >= mofa_network_shared_loading_threshold | edge_id %in% seed_edges$edge_id
    ) %>%
    dplyr::mutate(
      loading_sign = ifelse(weight >= 0, "Positive", "Negative"),
      loading_strength = abs(weight_within_view)
    )

  if (nrow(network_edges) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = paste(plot_title, "\n(no eligible edges)"))
    )
  }

  module_membership <- dplyr::bind_rows(
    lapply(
      split(network_edges, network_edges$view),
      function(view_edge_data) {
        profile <- view_edge_data %>%
          dplyr::select(feature_node_id, factor, weight_within_view) %>%
          tidyr::pivot_wider(names_from = factor, values_from = weight_within_view, values_fill = 0)
        for (factor_name in setdiff(factor_set, colnames(profile))) profile[[factor_name]] <- 0
        profile_matrix <- as.matrix(profile[, factor_set, drop = FALSE])
        rownames(profile_matrix) <- profile$feature_node_id
        module <- if (nrow(profile_matrix) <= 2) {
          rep(1L, nrow(profile_matrix))
        } else {
          stats::cutree(
            stats::hclust(stats::dist(profile_matrix), method = "ward.D2"),
            k = min(2L, max(1L, ceiling(nrow(profile_matrix) / 7)))
          )
        }
        data.frame(feature_node_id = rownames(profile_matrix), module = as.integer(module), stringsAsFactors = FALSE)
      }
    )
  )

  feature_nodes <- network_edges %>%
    dplyr::group_by(feature_node_id, view, view_label, feature, feature_label) %>%
    dplyr::summarise(
      n_connected_factors = dplyr::n_distinct(factor),
      maximum_loading = max(loading_strength, na.rm = TRUE),
      minimum_response_p = if (any(is.finite(p_value))) min(p_value, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::left_join(module_membership, by = "feature_node_id") %>%
    dplyr::mutate(
      view_order = match(view, required_views),
      shared_feature = n_connected_factors >= 2,
      response_feature = is.finite(minimum_response_p) & minimum_response_p < 0.05,
      feature_label_plotmath = mofa_v50_plotmath_feature_label(view, feature_label)
    ) %>%
    dplyr::arrange(view_order, module, dplyr::desc(n_connected_factors), dplyr::desc(maximum_loading))

  modules <- feature_nodes %>%
    dplyr::distinct(view, view_label, view_order, module) %>%
    dplyr::arrange(view_order, module) %>%
    dplyr::group_by(view) %>%
    dplyr::mutate(module_short = paste0("M", dplyr::row_number())) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      module_index = dplyr::row_number(),
      module_angle = pi / 2 - 2 * pi * (module_index - 1) / dplyr::n(),
      module_x = 3.55 * cos(module_angle),
      module_y = 3.05 * sin(module_angle)
    )

  feature_nodes <- feature_nodes %>%
    dplyr::left_join(modules %>% dplyr::select(view, module, module_short, module_x, module_y), by = c("view", "module")) %>%
    dplyr::group_by(view, module) %>%
    dplyr::mutate(
      feature_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.48,
      x = module_x,
      y = module_y + feature_offset,
      label_x = x + ifelse(module_x < -0.25, -0.95, 0.95),
      label_hjust = ifelse(module_x < -0.25, 1, 0)
    ) %>%
    dplyr::ungroup()

  modules <- modules %>%
    dplyr::left_join(
      feature_nodes %>%
        dplyr::group_by(view, module, module_short) %>%
        dplyr::summarise(
          xmin = min(x) - 0.44,
          xmax = max(x) + 0.44,
          ymin = min(y) - 0.30,
          ymax = max(y) + 0.30,
          label_y = max(y) + 0.64,
          .groups = "drop"
        ),
      by = c("view", "module", "module_short")
    )

  factor_nodes <- data.frame(
    factor = factor_set,
    factor_angle = pi / 2 - 2 * pi * (seq_along(factor_set) - 1) / length(factor_set),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      factor_rank = dplyr::row_number(),
      x = 7.95 * cos(factor_angle),
      y = 6.40 * sin(factor_angle),
      anchor_x = ifelse(x > 0, x - 0.55, x + 0.55),
      factor_label = factor
    )

  plot_edges <- network_edges %>%
    dplyr::left_join(factor_nodes %>% dplyr::select(factor, x_factor = x, y_factor = y, anchor_x, factor_rank), by = "factor") %>%
    dplyr::left_join(feature_nodes %>% dplyr::select(feature_node_id, x_feature = x, y_feature = y), by = "feature_node_id") %>%
    dplyr::group_by(factor) %>%
    dplyr::arrange(y_feature, .by_group = TRUE) %>%
    dplyr::mutate(
      factor_edge_rank = dplyr::row_number(),
      factor_edge_offset = (factor_edge_rank - (dplyr::n() + 1) / 2) * 0.11
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(feature_node_id) %>%
    dplyr::arrange(y_factor, .by_group = TRUE) %>%
    dplyr::mutate(
      feature_edge_rank = dplyr::row_number(),
      feature_edge_offset = (feature_edge_rank - (dplyr::n() + 1) / 2) * 0.09
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      x_start = anchor_x,
      y_start = y_factor + factor_edge_offset,
      x_end = x_feature,
      y_end = y_feature + feature_edge_offset,
      side = ifelse(x_factor > 0, "right", "left")
    )

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = modules,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = view_label),
      alpha = 0.10,
      color = "grey72",
      linewidth = 0.28
    ) +
    ggplot2::geom_segment(
      data = plot_edges,
      ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = loading_sign, linewidth = loading_strength),
      alpha = 0.30,
      lineend = "round"
    ) +
    ggplot2::geom_label(
      data = factor_nodes,
      ggplot2::aes(x = x, y = y, label = factor_label),
      fill = "#FFF2B3",
      color = "grey10",
      label.size = 0.28,
      size = 2.75,
      fontface = "plain"
    ) +
    ggplot2::geom_point(
      data = feature_nodes,
      ggplot2::aes(x = x, y = y, fill = view_label, size = maximum_loading),
      shape = 21,
      color = "grey20",
      stroke = 0.42
    ) +
    ggplot2::geom_point(
      data = feature_nodes %>% dplyr::filter(shared_feature),
      ggplot2::aes(x = x, y = y),
      shape = 21,
      size = 3.9,
      fill = NA,
      color = "#7A3E9D",
      stroke = 0.76
    ) +
    ggplot2::geom_text(
      data = feature_nodes %>% dplyr::filter(response_feature),
      ggplot2::aes(x = x, y = y, label = "*"),
      nudge_y = 0.28,
      size = 3.0,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = feature_nodes,
      ggplot2::aes(x = label_x, y = y, label = feature_label_plotmath, hjust = label_hjust),
      parse = TRUE,
      size = 1.90
    ) +
    ggplot2::geom_label(
      data = modules,
      ggplot2::aes(x = module_x, y = label_y, label = module_short, fill = view_label),
      color = "grey15",
      label.size = 0.16,
      size = 2.40,
      label.padding = grid::unit(0.09, "lines"),
      fontface = "plain",
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = c(Positive = "#D55E00", Negative = "#0072B2"), name = "Feature-weight sign") +
    ggplot2::scale_fill_manual(
      values = view_colors,
      name = "Feature group",
      guide = ggplot2::guide_legend(override.aes = list(shape = 21, size = 4, color = "grey20", alpha = 1))
    ) +
    ggplot2::scale_linewidth_continuous(range = c(0.22, 0.78), guide = "none") +
    ggplot2::scale_size_continuous(range = c(2.0, 3.5), guide = "none") +
    ggplot2::coord_equal(xlim = c(-10.4, 10.4), ylim = c(-7.8, 7.8), clip = "off") +
    ggplot2::labs(
      title = plot_title,
      subtitle = plot_subtitle,
      caption = "Direct factor-to-feature connectors use stronger start/end offsets to reduce overlap. Purple outlines denote features shared by multiple factors.",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_void(base_size = 9.0) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 7.6, color = "grey35"),
      plot.caption = ggplot2::element_text(size = 7.1, color = "grey35", hjust = 0),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(5, 255, 5, 255)
    )

  attr(p, "network_feature_nodes") <- feature_nodes
  attr(p, "network_edges") <- network_edges
  attr(p, "network_modules") <- modules
  attr(p, "network_factor_nodes") <- factor_nodes
  p
}

# Re-save the main 5-factor network with stronger edge de-overlap and direct connectors.
p_mofa_v50_multifactor_feature_network <- mofa_v50_build_network_plot(
  factor_set = mofa_v50_network_factors,
  plot_title = "Multi-factor feature-weight network",
  plot_subtitle = "Five selected factors with direct de-overlapped connectors; feature groups are identified in the legend."
)
if (FALSE) ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_multifactor_feature_loading_network.svg"),
  plot = p_mofa_v50_multifactor_feature_network,
  width = 14.8,
  height = 10.6,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

# Re-save the top-ranked selected-triple network.
mofa_v50_triple_network_factors <- c(
  mofa_v50_selected_triple_diagnostics$factor_x[1],
  mofa_v50_selected_triple_diagnostics$factor_y[1],
  mofa_v50_selected_triple_diagnostics$factor_z[1]
)
p_mofa_v50_selected_triple_feature_network <- mofa_v50_build_network_plot(
  factor_set = mofa_v50_triple_network_factors,
  plot_title = "Selected three-factor feature-weight network",
  plot_subtitle = "Top-ranked Factor7-fixed combination with direct de-overlapped connectors."
)
if (FALSE) ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_selected_triple_multifactor_feature_loading_network.svg"),
  plot = p_mofa_v50_selected_triple_feature_network,
  width = 14.2,
  height = 10.1,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

# Draw the five user-requested three-factor networks.
mofa_v50_requested_network_triplets <- list(
  c("Factor1", "Factor4", "Factor7"),
  c("Factor1", "Factor7", "Factor10"),
  c("Factor1", "Factor2", "Factor7"),
  c("Factor4", "Factor7", "Factor9"),
  c("Factor2", "Factor4", "Factor7")
)

mofa_v50_requested_network_labels <- vapply(
  mofa_v50_requested_network_triplets,
  function(x) paste(x, collapse = " + "),
  character(1)
)

mofa_v50_requested_network_plot_list <- lapply(
  seq_along(mofa_v50_requested_network_triplets),
  function(i) {
    factor_set <- mofa_v50_requested_network_triplets[[i]]
    factor_label <- paste(factor_set, collapse = " + ")
    plot_object <- mofa_v50_build_network_plot(
      factor_set = factor_set,
      plot_title = paste0("Feature-weight network: ", factor_label),
      plot_subtitle = "Direct factor-to-feature connectors with stronger offsets are used to reduce label and edge overlap."
    )
    file_stub <- paste(factor_set, collapse = "_")
    ggplot2::ggsave(
      filename = file.path(mofa_v50_figure_dir, paste0("MOFA_v50_network_", file_stub, ".svg")),
      plot = plot_object,
      width = 13.8,
      height = 9.8,
      units = "in",
      device = svglite::svglite,
      bg = "white"
    )
    plot_object
  }
)

p_mofa_v50_requested_networks_panel <- patchwork::wrap_plots(
  mofa_v50_requested_network_plot_list,
  ncol = 2,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "Requested three-factor feature-weight networks",
    subtitle = "Prespecified Factor7-containing combinations requested for visual comparison."
  )

# Combined requested-network patchwork disabled in v50; individual SVGs are saved below.



}

#-----------------------------------------------------------------#
# 22.9 v50 focused refinement: network readability, 3D annotation, paired-change focus, and compact atlas
#-----------------------------------------------------------------#

mofa_v50_factor_axis_label <- function(factor_name, max_views = 2) {
  view_info <- mofa_variance_explained %>%
    dplyr::mutate(factor = as.character(factor), view = as.character(view)) %>%
    dplyr::filter(factor == factor_name, is.finite(r2), r2 > 0) %>%
    dplyr::arrange(dplyr::desc(r2))
  if (nrow(view_info) == 0) return(factor_name)
  if (exists("view_labels")) {
    view_info$view_label <- dplyr::recode(view_info$view, !!!view_labels, .default = view_info$view)
  } else {
    view_info$view_label <- view_info$view
  }
  view_text <- paste0(view_info$view_label[seq_len(min(max_views, nrow(view_info)))], " ", sprintf("%.1f%%", 100 * view_info$r2[seq_len(min(max_views, nrow(view_info)))]) )
  paste0(factor_name, "\n", paste(view_text, collapse = "; "))
}

mofa_v50_draw_s3d_panel <- function(triple_factors, diagnostic_row, point_cex = 1.45, draw_all_plane_lines = TRUE) {
  if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
    plot.new()
    text(0.5, 0.6, "Package 'scatterplot3d' is required for v50 3D plots.")
    text(0.5, 0.48, paste(triple_factors, collapse = " + "), cex = 0.9)
    return(invisible(NULL))
  }
  plot_data <- mofa_v50_subject_mean_factor_data %>%
    dplyr::select(SubjectID, TRG_plot, dplyr::all_of(triple_factors)) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(triple_factors), is.finite), !is.na(TRG_plot)) %>%
    dplyr::mutate(TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR")))

  score_matrix <- scale(as.matrix(plot_data[, triple_factors, drop = FALSE]))
  score_matrix[!is.finite(score_matrix)] <- 0
  plot_data$x <- score_matrix[, 1]
  plot_data$y <- score_matrix[, 2]
  plot_data$z <- score_matrix[, 3]
  plot_data$color <- unname(mofa_response_colors[as.character(plot_data$TRG_plot)])
  plot_data$line_color <- grDevices::adjustcolor(plot_data$color, alpha.f = 0.12)
  plot_data$drop_color <- grDevices::adjustcolor(plot_data$color, alpha.f = 0.22)

  xlab_text <- mofa_v50_factor_axis_label(triple_factors[1])
  ylab_text <- mofa_v50_factor_axis_label(triple_factors[2])
  zlab_text <- mofa_v50_factor_axis_label(triple_factors[3])

  s3d <- scatterplot3d::scatterplot3d(
    x = plot_data$x,
    y = plot_data$y,
    z = plot_data$z,
    pch = 16,
    color = plot_data$color,
    cex.symbols = point_cex,
    type = "h",
    lty.hplot = 1,
    mar = c(2.2, 2.3, 3.3, 5.2),
    main = paste0("Rank ", diagnostic_row$selection_rank[1], ": ", paste(triple_factors, collapse = " + ")),
    xlab = xlab_text,
    ylab = ylab_text,
    zlab = zlab_text,
    angle = 52,
    scale.y = 1.0,
    box = TRUE,
    grid = TRUE
  )

  x_wall <- min(plot_data$x, na.rm = TRUE)
  y_wall <- max(plot_data$y, na.rm = TRUE)
  z_floor <- min(plot_data$z, na.rm = TRUE)

  for (i in seq_len(nrow(plot_data))) {
    xy_point <- s3d$xyz.convert(plot_data$x[i], plot_data$y[i], plot_data$z[i])
    xy_floor <- s3d$xyz.convert(plot_data$x[i], plot_data$y[i], z_floor)
    graphics::segments(xy_floor$x, xy_floor$y, xy_point$x, xy_point$y, col = plot_data$drop_color[i], lwd = 0.75)
    if (draw_all_plane_lines) {
      xy_xwall <- s3d$xyz.convert(x_wall, plot_data$y[i], plot_data$z[i])
      xy_ywall <- s3d$xyz.convert(plot_data$x[i], y_wall, plot_data$z[i])
      graphics::segments(xy_xwall$x, xy_xwall$y, xy_point$x, xy_point$y, col = plot_data$line_color[i], lwd = 0.55)
      graphics::segments(xy_ywall$x, xy_ywall$y, xy_point$x, xy_point$y, col = plot_data$line_color[i], lwd = 0.55)
    }
  }
  point_xy <- s3d$xyz.convert(plot_data$x, plot_data$y, plot_data$z)
  graphics::points(point_xy$x, point_xy$y, pch = 21, bg = plot_data$color, col = "grey20", cex = point_cex)

  stats_label <- paste0(
    "PERMANOVA p = ", ifelse(is.na(diagnostic_row$permanova_p[1]), "NA", ifelse(diagnostic_row$permanova_p[1] < 0.001, "<0.001", sprintf("%.3f", diagnostic_row$permanova_p[1]))),
    " | R2 = ", ifelse(is.na(diagnostic_row$permanova_r2[1]), "NA", sprintf("%.2f", diagnostic_row$permanova_r2[1])),
    "\nEnergy p = ", ifelse(is.na(diagnostic_row$energy_p[1]), "NA", ifelse(diagnostic_row$energy_p[1] < 0.001, "<0.001", sprintf("%.3f", diagnostic_row$energy_p[1]))),
    " | E = ", ifelse(is.na(diagnostic_row$energy_statistic[1]), "NA", sprintf("%.2f", diagnostic_row$energy_statistic[1]))
  )
  graphics::mtext(stats_label, side = 3, line = 0.35, adj = 0.02, cex = 0.72)
  graphics::legend(
    x = grconvertX(1.04, from = "npc", to = "user"),
    y = grconvertY(0.94, from = "npc", to = "user"),
    legend = c("pCR", "non-pCR"),
    pt.bg = unname(mofa_response_colors[c("pCR", "non_pCR")]),
    pch = 21,
    col = "grey20",
    pt.cex = 1.3,
    cex = 0.80,
    bty = "n",
    xpd = NA,
    xjust = 0,
    yjust = 1
  )
  invisible(NULL)
}

# add requested Factor1+Factor4+Factor7 panel alongside the top 3 combinations.
mofa_v50_manual_three_factor_rows <- dplyr::bind_rows(
  mofa_v50_selected_triple_diagnostics,
  mofa_v50_triple_permanova_atlas %>%
    dplyr::filter(
      factor_x == "Factor7",
      (factor_y == "Factor1" & factor_z == "Factor4") | (factor_y == "Factor4" & factor_z == "Factor1")
    ) %>%
    dplyr::slice_head(n = 1)
) %>% dplyr::distinct(triple_id, .keep_all = TRUE)

svglite::svglite(
  file.path(mofa_v50_figure_dir, "MOFA_v50_three_factor_response_3D_top3_plus_F1_F4_F7.svg"),
  width = 20.5,
  height = 5.8,
  bg = "white"
)
graphics::par(mfrow = c(1, max(1, nrow(mofa_v50_manual_three_factor_rows))), mar = c(2.0, 2.1, 3.2, 4.6), xpd = NA)
for (i in seq_len(nrow(mofa_v50_manual_three_factor_rows))) {
  mofa_v50_draw_s3d_panel(
    c(mofa_v50_manual_three_factor_rows$factor_x[i], mofa_v50_manual_three_factor_rows$factor_y[i], mofa_v50_manual_three_factor_rows$factor_z[i]),
    mofa_v50_manual_three_factor_rows[i, , drop = FALSE],
    point_cex = 1.55,
    draw_all_plane_lines = TRUE
  )
}
grDevices::dev.off()

# rebuild the top3 figure with larger points and multi-plane faint projections.
svglite::svglite(
  file.path(mofa_v50_figure_dir, "MOFA_v50_three_factor_response_3D_top3.svg"),
  width = 16.6,
  height = 5.8,
  bg = "white"
)
graphics::par(mfrow = c(1, max(1, nrow(mofa_v50_selected_triple_diagnostics))), mar = c(2.0, 2.1, 3.2, 4.6), xpd = NA)
for (i in seq_len(nrow(mofa_v50_selected_triple_diagnostics))) {
  mofa_v50_draw_s3d_panel(
    c(mofa_v50_selected_triple_diagnostics$factor_x[i], mofa_v50_selected_triple_diagnostics$factor_y[i], mofa_v50_selected_triple_diagnostics$factor_z[i]),
    mofa_v50_selected_triple_diagnostics[i, , drop = FALSE],
    point_cex = 1.55,
    draw_all_plane_lines = TRUE
  )
}
grDevices::dev.off()

mofa_v50_build_network_plot <- function(
  factor_set,
  plot_title,
  plot_subtitle = NULL
) {
  factor_set <- unique(factor_set)
  factor_set <- factor_set[factor_set %in% mofa_v50_factor_order]
  if (length(factor_set) < 2) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Not enough factors to build the network."))
  }

  active_factor_views <- mofa_variance_explained %>%
    dplyr::mutate(factor = as.character(factor), view = as.character(view)) %>%
    dplyr::filter(factor %in% factor_set, r2 >= mofa_active_view_r2) %>%
    dplyr::select(factor, view, view_r2 = r2)

  seed_edges <- dplyr::bind_rows(
    mofa_feature_weights %>%
      dplyr::filter(factor %in% factor_set, display_eligible) %>%
      dplyr::inner_join(active_factor_views, by = c("factor", "view")) %>%
      dplyr::group_by(factor, view, direction) %>%
      dplyr::slice_max(order_by = abs(weight_within_view), n = mofa_network_features_per_view_direction, with_ties = FALSE) %>%
      dplyr::ungroup(),
    mofa_feature_weights %>%
      dplyr::filter(factor %in% factor_set, display_eligible, is.finite(p_value), p_value < 0.05) %>%
      dplyr::inner_join(active_factor_views, by = c("factor", "view")) %>%
      dplyr::group_by(factor) %>%
      dplyr::arrange(p_value, dplyr::desc(abs(weight_within_view)), .by_group = TRUE) %>%
      dplyr::slice_head(n = 2) %>%
      dplyr::ungroup()
  ) %>%
    dplyr::distinct(factor, view, feature, .keep_all = TRUE) %>%
    dplyr::mutate(feature_node_id = paste(view, feature, sep = "::"), edge_id = paste(factor, feature_node_id, sep = "__"))

  network_edges <- mofa_feature_weights %>%
    dplyr::filter(factor %in% factor_set, display_eligible) %>%
    dplyr::inner_join(active_factor_views, by = c("factor", "view")) %>%
    dplyr::mutate(feature_node_id = paste(view, feature, sep = "::"), edge_id = paste(factor, feature_node_id, sep = "__")) %>%
    dplyr::filter(feature_node_id %in% unique(seed_edges$feature_node_id), abs(weight_within_view) >= mofa_network_shared_loading_threshold | edge_id %in% seed_edges$edge_id) %>%
    dplyr::mutate(loading_sign = ifelse(weight >= 0, "Positive", "Negative"), loading_strength = abs(weight_within_view))

  module_membership <- dplyr::bind_rows(lapply(split(network_edges, network_edges$view), function(view_edge_data) {
    profile <- view_edge_data %>%
      dplyr::select(feature_node_id, factor, weight_within_view) %>%
      tidyr::pivot_wider(names_from = factor, values_from = weight_within_view, values_fill = 0)
    for (factor_name in setdiff(factor_set, colnames(profile))) profile[[factor_name]] <- 0
    profile_matrix <- as.matrix(profile[, factor_set, drop = FALSE])
    rownames(profile_matrix) <- profile$feature_node_id
    module <- if (nrow(profile_matrix) <= 2) rep(1L, nrow(profile_matrix)) else stats::cutree(stats::hclust(stats::dist(profile_matrix), method = "ward.D2"), k = min(2L, max(1L, ceiling(nrow(profile_matrix) / 7))))
    data.frame(feature_node_id = rownames(profile_matrix), module = as.integer(module), stringsAsFactors = FALSE)
  }))

  feature_nodes <- network_edges %>%
    dplyr::group_by(feature_node_id, view, view_label, feature, feature_label) %>%
    dplyr::summarise(
      n_connected_factors = dplyr::n_distinct(factor),
      maximum_loading = max(loading_strength, na.rm = TRUE),
      minimum_response_p = if (any(is.finite(p_value))) min(p_value, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::left_join(module_membership, by = "feature_node_id") %>%
    dplyr::mutate(
      view_order = match(view, required_views),
      shared_feature = n_connected_factors >= 2,
      response_feature = is.finite(minimum_response_p) & minimum_response_p < 0.05,
      feature_label_plotmath = mofa_v50_plotmath_feature_label(view, feature_label)
    ) %>%
    dplyr::arrange(view_order, module, dplyr::desc(n_connected_factors), dplyr::desc(maximum_loading))

  modules <- feature_nodes %>%
    dplyr::distinct(view, view_label, view_order, module) %>%
    dplyr::arrange(view_order, module) %>%
    dplyr::group_by(view) %>%
    dplyr::mutate(module_short = paste0("M", dplyr::row_number())) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(module_index = dplyr::row_number(), module_angle = pi / 2 - 2 * pi * (module_index - 1) / dplyr::n(), module_x = 3.75 * cos(module_angle), module_y = 3.20 * sin(module_angle))

  feature_nodes <- feature_nodes %>%
    dplyr::left_join(modules %>% dplyr::select(view, module, module_short, module_x, module_y), by = c("view", "module")) %>%
    dplyr::group_by(view, module) %>%
    dplyr::mutate(feature_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.52, x = module_x, y = module_y + feature_offset, label_x = x + ifelse(module_x < -0.25, -1.05, 1.05), label_hjust = ifelse(module_x < -0.25, 1, 0)) %>%
    dplyr::ungroup()

  modules <- modules %>%
    dplyr::left_join(feature_nodes %>% dplyr::group_by(view, module, module_short) %>% dplyr::summarise(xmin = min(x) - 0.46, xmax = max(x) + 0.46, ymin = min(y) - 0.32, ymax = max(y) + 0.32, label_y = max(y) + 0.74, .groups = "drop"), by = c("view", "module", "module_short"))

  factor_nodes <- data.frame(factor = factor_set, factor_angle = pi / 2 - 2 * pi * (seq_along(factor_set) - 1) / length(factor_set), stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      x = 8.40 * cos(factor_angle),
      y = 6.65 * sin(factor_angle),
      x_hub = 7.10 * cos(factor_angle),
      y_hub = 5.55 * sin(factor_angle),
      factor_label = factor
    )

  plot_edges <- network_edges %>%
    dplyr::left_join(factor_nodes %>% dplyr::select(factor, x_hub, y_hub), by = "factor") %>%
    dplyr::left_join(feature_nodes %>% dplyr::select(feature_node_id, x_feature = x, y_feature = y), by = "feature_node_id") %>%
    dplyr::group_by(factor) %>%
    dplyr::arrange(y_feature, .by_group = TRUE) %>%
    dplyr::mutate(start_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.14) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(feature_node_id) %>%
    dplyr::arrange(y_hub, .by_group = TRUE) %>%
    dplyr::mutate(end_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.11) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(x_start = x_hub, y_start = y_hub + start_offset, x_end = x_feature, y_end = y_feature + end_offset)

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = modules, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = view_label), alpha = 0.10, color = "grey72", linewidth = 0.28) +
    ggplot2::geom_segment(data = factor_nodes, ggplot2::aes(x = x_hub, y = y_hub, xend = x, yend = y), color = "grey70", linewidth = 0.36) +
    ggplot2::geom_segment(data = plot_edges, ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = loading_sign, linewidth = loading_strength), alpha = 0.24, lineend = "round") +
    ggplot2::geom_label(data = factor_nodes, ggplot2::aes(x = x, y = y, label = factor_label), fill = "#FFF2B3", color = "grey10", label.size = 0.28, size = 2.75, fontface = "plain") +
    ggplot2::geom_point(data = factor_nodes, ggplot2::aes(x = x_hub, y = y_hub), shape = 21, size = 1.9, stroke = 0.25, fill = "white", color = "grey35") +
    ggplot2::geom_point(data = feature_nodes, ggplot2::aes(x = x, y = y, fill = view_label, size = maximum_loading), shape = 21, color = "grey20", stroke = 0.42) +
    ggplot2::geom_point(data = feature_nodes %>% dplyr::filter(shared_feature), ggplot2::aes(x = x, y = y), shape = 21, size = 3.9, fill = NA, color = "#7A3E9D", stroke = 0.76) +
    ggplot2::geom_text(data = feature_nodes %>% dplyr::filter(response_feature), ggplot2::aes(x = x, y = y, label = "*"), nudge_y = 0.30, size = 3.0, color = "black") +
    ggplot2::geom_label(data = modules, ggplot2::aes(x = module_x, y = label_y, label = module_short, fill = view_label), color = "grey15", label.size = 0.16, size = 2.36, label.padding = grid::unit(0.08, "lines"), fontface = "plain", show.legend = FALSE)

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = feature_nodes,
      ggplot2::aes(x = label_x, y = y, label = feature_label_plotmath, hjust = label_hjust),
      parse = TRUE,
      size = 1.88,
      min.segment.length = 0,
      seed = 20261040,
      force = 0.8,
      max.overlaps = 120,
      max.time = 1.0,
      max.iter = 3000,
      direction = "y",
      box.padding = 0.10,
      point.padding = 0.08,
      segment.alpha = 0.18,
      segment.size = 0.18
    )
  } else {
    p <- p + ggplot2::geom_text(data = feature_nodes, ggplot2::aes(x = label_x, y = y, label = feature_label_plotmath, hjust = label_hjust), parse = TRUE, size = 1.88)
  }

  p +
    ggplot2::scale_color_manual(values = c(Positive = "#D55E00", Negative = "#0072B2"), name = "Feature-weight sign") +
    ggplot2::scale_fill_manual(values = view_colors, name = "Feature group", guide = ggplot2::guide_legend(override.aes = list(shape = 21, size = 4, color = "grey20", alpha = 1))) +
    ggplot2::scale_linewidth_continuous(range = c(0.18, 0.68), guide = "none") +
    ggplot2::scale_size_continuous(range = c(2.0, 3.4), guide = "none") +
    ggplot2::coord_equal(xlim = c(-11.2, 11.2), ylim = c(-8.2, 8.2), clip = "off") +
    ggplot2::labs(title = plot_title, subtitle = plot_subtitle, caption = "Factor hubs and repelled feature labels are used to minimize overlaps. Purple outlines denote features shared by multiple factors.", x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 9.0) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), plot.subtitle = ggplot2::element_text(size = 7.6, color = "grey35"), plot.caption = ggplot2::element_text(size = 7.0, color = "grey35", hjust = 0), legend.position = "bottom", plot.margin = ggplot2::margin(5, 270, 5, 270))
}

# rebuild requested network figures with stronger de-overlap logic.
p_mofa_v50_multifactor_feature_network <- mofa_v50_build_network_plot(
  factor_set = mofa_v50_network_factors,
  plot_title = "Multi-factor feature-weight network",
  plot_subtitle = "Expanded spacing and hub-linked direct connectors improve factor-feature alignment."
)
ggplot2::ggsave(file.path(mofa_v50_figure_dir, "MOFA_v50_multifactor_feature_loading_network.svg"), p_mofa_v50_multifactor_feature_network, width = 15.2, height = 10.8, units = "in", device = svglite::svglite, bg = "white")

p_mofa_v50_selected_triple_feature_network <- mofa_v50_build_network_plot(
  factor_set = c(mofa_v50_selected_triple_diagnostics$factor_x[1], mofa_v50_selected_triple_diagnostics$factor_y[1], mofa_v50_selected_triple_diagnostics$factor_z[1]),
  plot_title = "Selected three-factor feature-weight network",
  plot_subtitle = "Top-ranked Factor7-fixed combination with repelled labels and hub-linked direct connectors."
)
ggplot2::ggsave(file.path(mofa_v50_figure_dir, "MOFA_v50_selected_triple_multifactor_feature_loading_network.svg"), p_mofa_v50_selected_triple_feature_network, width = 14.6, height = 10.2, units = "in", device = svglite::svglite, bg = "white")

mofa_v50_requested_network_triplets <- list(
  c("Factor1", "Factor4", "Factor7"),
  c("Factor1", "Factor7", "Factor10"),
  c("Factor1", "Factor2", "Factor7"),
  c("Factor4", "Factor7", "Factor9"),
  c("Factor2", "Factor4", "Factor7")
)

mofa_v50_requested_network_plot_list <- vector(
  "list",
  length(mofa_v50_requested_network_triplets)
)

for (i in seq_along(mofa_v50_requested_network_triplets)) {
  factor_set <- mofa_v50_requested_network_triplets[[i]]
  factor_label <- paste(factor_set, collapse = " + ")
  output_file <- file.path(
    mofa_v50_figure_dir,
    paste0("MOFA_v50_network_", paste(factor_set, collapse = "_"), ".svg")
  )
  start_time <- Sys.time()
  message(
    sprintf(
      "[%d/%d] Rendering network: %s",
      i,
      length(mofa_v50_requested_network_triplets),
      factor_label
    )
  )
  plot_object <- mofa_v50_build_network_plot(
    factor_set = factor_set,
    plot_title = paste0("Feature-weight network: ", factor_label),
    plot_subtitle = "Repelled labels and hub-linked direct connectors are used to reduce overlap."
  )
  ggplot2::ggsave(
    filename = output_file,
    plot = plot_object,
    width = 14.4,
    height = 10.0,
    units = "in",
    device = svglite::svglite,
    bg = "white",
    limitsize = FALSE
  )
  message(
    sprintf(
      "Saved %s (%.1f min)",
      basename(output_file),
      as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    )
  )
  mofa_v50_requested_network_plot_list[[i]] <- output_file
  rm(plot_object)
  invisible(gc())
}

# A single five-network patchwork is intentionally not rendered in v50.
# Each network uses coord_equal(), many parsed labels, and repulsion. Combining
# all five forces patchwork to solve fixed-aspect layouts and re-render every
# network a second time, which is both slow and prone to the misleading
# "RStudio Plots window may be too small" error during ggsave().
p_mofa_v50_requested_networks_panel <- NULL
write.csv(
  data.frame(
    factor_combination = vapply(
      mofa_v50_requested_network_triplets,
      paste,
      collapse = " + ",
      FUN.VALUE = character(1)
    ),
    svg_file = unlist(mofa_v50_requested_network_plot_list),
    stringsAsFactors = FALSE
  ),
  file.path(mofa_v50_result_dir, "MOFA_v50_requested_network_files.csv"),
  row.names = FALSE
)
message(
  "Skipped the combined requested-network patchwork; five individual SVG files were saved instead."
)

# Focused paired-change plots: strongest positive/negative paired overall change and differential paired change.
mofa_v50_paired_overall_rank <- mofa_factor_wilcoxon_tests %>%
  dplyr::filter(comparison == "paired_overall_change", factor %in% mofa_v50_factor_order, is.finite(wilcoxon_p), is.finite(wilcoxon_effect))
mofa_v50_paired_diff_rank <- mofa_factor_wilcoxon_tests %>%
  dplyr::filter(comparison == "paired_differential_change", factor %in% mofa_v50_factor_order, is.finite(wilcoxon_p), is.finite(wilcoxon_effect))

mofa_v50_best_paired_positive <- mofa_v50_paired_overall_rank %>% dplyr::filter(wilcoxon_effect > 0) %>% dplyr::arrange(wilcoxon_p, dplyr::desc(wilcoxon_effect)) %>% dplyr::slice_head(n = 1) %>% dplyr::pull(factor)
mofa_v50_best_paired_negative <- mofa_v50_paired_overall_rank %>% dplyr::filter(wilcoxon_effect < 0) %>% dplyr::arrange(wilcoxon_p, wilcoxon_effect) %>% dplyr::slice_head(n = 1) %>% dplyr::pull(factor)
mofa_v50_best_diff_positive <- mofa_v50_paired_diff_rank %>% dplyr::filter(wilcoxon_effect > 0) %>% dplyr::arrange(wilcoxon_p, dplyr::desc(wilcoxon_effect)) %>% dplyr::slice_head(n = 1) %>% dplyr::pull(factor)
mofa_v50_best_diff_negative <- mofa_v50_paired_diff_rank %>% dplyr::filter(wilcoxon_effect < 0) %>% dplyr::arrange(wilcoxon_p, wilcoxon_effect) %>% dplyr::slice_head(n = 1) %>% dplyr::pull(factor)
if (length(mofa_v50_best_paired_positive) == 0) mofa_v50_best_paired_positive <- mofa_v50_paired_overall_rank$factor[1]
if (length(mofa_v50_best_paired_negative) == 0) mofa_v50_best_paired_negative <- mofa_v50_paired_overall_rank$factor[min(2, nrow(mofa_v50_paired_overall_rank))]
if (length(mofa_v50_best_diff_positive) == 0) mofa_v50_best_diff_positive <- mofa_v50_paired_diff_rank$factor[1]
if (length(mofa_v50_best_diff_negative) == 0) mofa_v50_best_diff_negative <- mofa_v50_paired_diff_rank$factor[min(2, nrow(mofa_v50_paired_diff_rank))]

mofa_v50_make_signed_rank_label <- function(data_frame) {
  if (nrow(data_frame) < 3) return("paired p = NA")
  wt <- tryCatch(stats::wilcox.test(data_frame$factor_score[data_frame$Timepoint_display == levels(data_frame$Timepoint_display)[2]], data_frame$factor_score[data_frame$Timepoint_display == levels(data_frame$Timepoint_display)[1]], paired = TRUE, exact = FALSE), error = function(e) NULL)
  if (is.null(wt)) return("paired p = NA")
  paste0("paired p = ", ifelse(wt$p.value < 0.001, "< 0.001", sprintf("%.3f", wt$p.value)))
}

mofa_v50_build_paired_factor_plot <- function(factor_name, split_by_response = FALSE, plot_title = NULL) {
  plot_data <- mofa_v50_paired_long %>% dplyr::filter(as.character(factor) == factor_name)
  if (!split_by_response) {
    label_text <- mofa_v50_make_signed_rank_label(plot_data)
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Timepoint_display, y = factor_score, group = SubjectID)) +
      ggplot2::geom_line(color = "grey72", linewidth = 0.42, alpha = 0.8) +
      ggplot2::geom_point(ggplot2::aes(fill = Timepoint_display), shape = 21, size = 1.9, color = "grey20", stroke = 0.28) +
      ggplot2::scale_fill_manual(values = timepoint_display_colors, guide = "none") +
      ggplot2::annotate("text", x = 1.5, y = max(plot_data$factor_score, na.rm = TRUE) + 0.12 * diff(range(plot_data$factor_score, na.rm = TRUE)), label = label_text, size = 2.6) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.18))) +
      ggplot2::labs(title = plot_title %||% factor_name, x = NULL, y = "Factor score") +
      ggplot2::theme_classic(base_size = 8.9) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), axis.text.x = ggplot2::element_text(face = "plain"))
    return(p)
  }
  label_df <- plot_data %>%
    dplyr::group_by(TRG_plot) %>%
    dplyr::group_modify(~ data.frame(label = mofa_v50_make_signed_rank_label(.x), y = max(.x$factor_score, na.rm = TRUE) + 0.12 * diff(range(.x$factor_score, na.rm = TRUE))))
  ggplot2::ggplot(plot_data, ggplot2::aes(x = Timepoint_display, y = factor_score, group = SubjectID)) +
    ggplot2::geom_line(color = "grey72", linewidth = 0.40, alpha = 0.8) +
    ggplot2::geom_point(ggplot2::aes(fill = Timepoint_display), shape = 21, size = 1.8, color = "grey20", stroke = 0.28) +
    ggplot2::geom_text(data = label_df, ggplot2::aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 2.5) +
    ggplot2::scale_fill_manual(values = timepoint_display_colors, guide = "none") +
    ggplot2::facet_wrap(~ TRG_plot, nrow = 1, labeller = ggplot2::as_labeller(c(pCR = "pCR", non_pCR = "non-pCR"))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.20))) +
    ggplot2::labs(title = plot_title %||% factor_name, x = NULL, y = "Factor score") +
    ggplot2::theme_classic(base_size = 8.8) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"), strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "plain"), axis.text.x = ggplot2::element_text(face = "plain"))
}

`%||%` <- function(x, y) if (!is.null(x)) x else y

mofa_v50_paired_focus_plots <- list(
  mofa_v50_build_paired_factor_plot(mofa_v50_best_paired_positive, FALSE, paste0(mofa_v50_best_paired_positive, ": strongest positive paired change")),
  mofa_v50_build_paired_factor_plot(mofa_v50_best_paired_negative, FALSE, paste0(mofa_v50_best_paired_negative, ": strongest negative paired change")),
  mofa_v50_build_paired_factor_plot(mofa_v50_best_diff_positive, TRUE, paste0(mofa_v50_best_diff_positive, ": positive differential paired change")),
  mofa_v50_build_paired_factor_plot(mofa_v50_best_diff_negative, TRUE, paste0(mofa_v50_best_diff_negative, ": negative differential paired change"))
)

p_mofa_v50_paired_change_focus_panel <- patchwork::wrap_plots(mofa_v50_paired_focus_plots, ncol = 2) + patchwork::plot_annotation(title = "Focused paired-change factor displays", subtitle = "Overall paired change is tested by the paired Wilcoxon signed-rank test; differential paired-change panels show within-group paired tests for pCR and non-pCR separately.")
ggplot2::ggsave(file.path(mofa_v50_figure_dir, "MOFA_v50_paired_change_focus_panel.svg"), p_mofa_v50_paired_change_focus_panel, width = 10.2, height = 8.0, units = "in", device = svglite::svglite, bg = "white")

# Compact the atlas.
p_mofa_v50_factor_interpretation_atlas <-
  (p_mofa_v50_variance_heatmap_atlas | p_mofa_v50_response_enrichment_atlas | p_mofa_v50_association_delta_heatmap) +
  patchwork::plot_layout(widths = c(0.80, 0.92, 0.54), guides = "collect") +
  patchwork::plot_annotation(
    title = "MOFA factor interpretation atlas",
    subtitle = paste0(
      "Variance explained provides the omics context, response enrichment is the primary clinical summary, ",
      "and the right heatmap highlights paired-change associations."
    )
  ) &
  ggplot2::theme(legend.position = "right")

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_variance_association_atlas.svg"),
  plot = p_mofa_v50_factor_interpretation_atlas,
  width = 9.6,
  height = max(5.6, 2.8 + 0.31 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)


#-----------------------------------------------------------------#
# 22.9 v50 focused layout refinement for enrichment, networks, and Factor1+4+7 displays
#-----------------------------------------------------------------#

#------------------------------#
# 22.9.1 Response enrichment: narrower width and coarser size legend
#------------------------------#

mofa_v50_response_size_bucket <- function(x) {
  x <- ifelse(is.finite(x), x, 0)
  as.character(pmin(2, pmax(0, round(x))))
}

mofa_v50_response_enrichment_data <- mofa_v50_response_enrichment_data %>%
  dplyr::mutate(
    minus_log10_p_bucket = factor(
      mofa_v50_response_size_bucket(minus_log10_p),
      levels = c("0", "1", "2")
    )
  )

p_mofa_v50_response_enrichment <- ggplot2::ggplot(
  mofa_v50_response_enrichment_data,
  ggplot2::aes(
    x = standardized_effect,
    y = factor,
    color = response_direction,
    fill = response_direction,
    size = minus_log10_p_bucket
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.40
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = standardized_effect, yend = factor, color = response_direction),
    linewidth = 0.58,
    alpha = 0.74,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(shape = 21, color = "grey20", stroke = 0.55) +
  ggplot2::geom_text(
    ggplot2::aes(label = significance_label),
    nudge_y = 0.22,
    size = 2.9,
    color = "black",
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~ Timepoint_display, nrow = 1) +
  ggplot2::scale_color_manual(
    values = mofa_group_colors[c("pCR", "non_pCR")],
    labels = c(pCR = "pCR-enriched", non_pCR = "non-pCR-enriched"),
    name = "Direction"
  ) +
  ggplot2::scale_fill_manual(values = mofa_group_colors[c("pCR", "non_pCR")], guide = "none") +
  ggplot2::scale_size_manual(
    values = c("0" = 2.1, "1" = 4.1, "2" = 6.0),
    name = expression(-log[10](italic(p))),
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "pCR versus non-pCR factor enrichment",
    subtitle = paste0(
      "Hedges' g is oriented as pCR minus non-pCR. Point size uses a coarse ",
      "-log10(p) scale of 0, 1, and 2."
    ),
    x = "Standardized response enrichment (Hedges' g)",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.6, color = "grey35"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "plain"),
    axis.text.y = ggplot2::element_text(face = "plain", size = 7.5),
    legend.position = "right"
  )

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_enrichment.svg"),
  plot = p_mofa_v50_response_enrichment,
  width = 4.8,
  height = max(4.5, 1.8 + 0.28 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

p_mofa_v50_response_enrichment_atlas <- ggplot2::ggplot(
  mofa_v50_response_enrichment_data,
  ggplot2::aes(
    x = standardized_effect,
    y = factor,
    color = response_direction,
    fill = response_direction,
    size = minus_log10_p_bucket
  )
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.38) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = standardized_effect, yend = factor, color = response_direction),
    linewidth = 0.56,
    alpha = 0.74,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(shape = 21, color = "grey20", stroke = 0.50) +
  ggplot2::facet_wrap(~ Timepoint_display, nrow = 1) +
  ggplot2::scale_color_manual(
    values = mofa_group_colors[c("pCR", "non_pCR")],
    labels = c(pCR = "pCR-enriched", non_pCR = "non-pCR-enriched"),
    name = "Direction"
  ) +
  ggplot2::scale_fill_manual(values = mofa_group_colors[c("pCR", "non_pCR")], guide = "none") +
  ggplot2::scale_size_manual(values = c("0" = 2.0, "1" = 3.8, "2" = 5.5), name = expression(-log[10](italic(p))), drop = FALSE) +
  ggplot2::labs(x = "Response enrichment (Hedges' g)", y = NULL) +
  ggplot2::theme_classic(base_size = 8.8) +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "plain"),
    axis.text.y = ggplot2::element_text(face = "plain", size = 7.3),
    legend.position = "right"
  )

#------------------------------#
# 22.9.2 Simpler all-multi-factor network layout: factors on one side, module-feature blocks on the other
#------------------------------#

mofa_v50_build_network_plot <- function(
  factor_set,
  plot_title,
  plot_subtitle = NULL,
  edge_spread = 1.0,
  factor_layout = c("left"),
  sort_focus = NULL
) {
  factor_layout <- match.arg(factor_layout)
  factor_set <- unique(factor_set)
  factor_set <- factor_set[factor_set %in% mofa_v50_factor_order]
  if (length(factor_set) < 2) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Not enough factors to build the network."))
  }

  active_factor_views <- mofa_variance_explained %>%
    dplyr::mutate(factor = as.character(factor), view = as.character(view)) %>%
    dplyr::filter(factor %in% factor_set, r2 >= mofa_active_view_r2) %>%
    dplyr::select(factor, view, view_r2 = r2)

  seed_edges <- dplyr::bind_rows(
    mofa_feature_weights %>%
      dplyr::filter(factor %in% factor_set, display_eligible) %>%
      dplyr::inner_join(active_factor_views, by = c("factor", "view")) %>%
      dplyr::group_by(factor, view, direction) %>%
      dplyr::slice_max(order_by = abs(weight_within_view), n = mofa_network_features_per_view_direction, with_ties = FALSE) %>%
      dplyr::ungroup(),
    mofa_feature_weights %>%
      dplyr::filter(factor %in% factor_set, display_eligible, is.finite(p_value), p_value < 0.05) %>%
      dplyr::inner_join(active_factor_views, by = c("factor", "view")) %>%
      dplyr::group_by(factor) %>%
      dplyr::arrange(p_value, dplyr::desc(abs(weight_within_view)), .by_group = TRUE) %>%
      dplyr::slice_head(n = 2) %>%
      dplyr::ungroup()
  ) %>%
    dplyr::distinct(factor, view, feature, .keep_all = TRUE) %>%
    dplyr::mutate(feature_node_id = paste(view, feature, sep = "::"), edge_id = paste(factor, feature_node_id, sep = "__"))

  network_edges <- mofa_feature_weights %>%
    dplyr::filter(factor %in% factor_set, display_eligible) %>%
    dplyr::inner_join(active_factor_views, by = c("factor", "view")) %>%
    dplyr::mutate(feature_node_id = paste(view, feature, sep = "::"), edge_id = paste(factor, feature_node_id, sep = "__")) %>%
    dplyr::filter(
      feature_node_id %in% unique(seed_edges$feature_node_id),
      abs(weight_within_view) >= mofa_network_shared_loading_threshold | edge_id %in% seed_edges$edge_id
    ) %>%
    dplyr::mutate(loading_sign = ifelse(weight >= 0, "Positive", "Negative"), loading_strength = abs(weight_within_view))

  if (nrow(network_edges) == 0) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No eligible edges."))
  }

  module_membership <- dplyr::bind_rows(lapply(split(network_edges, network_edges$view), function(view_edge_data) {
    profile <- view_edge_data %>%
      dplyr::select(feature_node_id, factor, weight_within_view) %>%
      tidyr::pivot_wider(names_from = factor, values_from = weight_within_view, values_fill = 0)
    for (factor_name in setdiff(factor_set, colnames(profile))) profile[[factor_name]] <- 0
    profile_matrix <- as.matrix(profile[, factor_set, drop = FALSE])
    rownames(profile_matrix) <- profile$feature_node_id
    module <- if (nrow(profile_matrix) <= 2) rep(1L, nrow(profile_matrix)) else stats::cutree(stats::hclust(stats::dist(profile_matrix), method = "ward.D2"), k = min(2L, max(1L, ceiling(nrow(profile_matrix) / 7))))
    data.frame(feature_node_id = rownames(profile_matrix), module = as.integer(module), stringsAsFactors = FALSE)
  }))

  feature_nodes <- network_edges %>%
    dplyr::group_by(feature_node_id, view, view_label, feature, feature_label) %>%
    dplyr::summarise(
      n_connected_factors = dplyr::n_distinct(factor),
      maximum_loading = max(loading_strength, na.rm = TRUE),
      minimum_response_p = if (any(is.finite(p_value))) min(p_value, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::left_join(module_membership, by = "feature_node_id") %>%
    dplyr::mutate(
      view_order = match(view, required_views),
      shared_feature = n_connected_factors >= 2,
      response_feature = is.finite(minimum_response_p) & minimum_response_p < 0.05,
      feature_label_plotmath = mofa_v50_plotmath_feature_label(view, feature_label)
    ) %>%
    dplyr::arrange(view_order, module, dplyr::desc(n_connected_factors), dplyr::desc(maximum_loading), feature_label)

  modules <- feature_nodes %>%
    dplyr::distinct(view, view_label, view_order, module) %>%
    dplyr::arrange(view_order, module) %>%
    dplyr::group_by(view) %>%
    dplyr::mutate(module_short = paste0("M", dplyr::row_number())) %>%
    dplyr::ungroup()

  feature_nodes <- feature_nodes %>%
    dplyr::left_join(modules, by = c("view", "view_label", "view_order", "module"))

  module_size_df <- feature_nodes %>%
    dplyr::group_by(view_order, view_label, module, module_short) %>%
    dplyr::summarise(n_features = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(view_order, module)

  gap_between_modules <- 1.2
  current_top <- 0
  module_rows <- vector("list", nrow(module_size_df))
  for (i in seq_len(nrow(module_size_df))) {
    n_i <- module_size_df$n_features[i]
    top_y <- -current_top
    bottom_y <- -(current_top + n_i - 1)
    module_rows[[i]] <- data.frame(
      view_order = module_size_df$view_order[i],
      view_label = module_size_df$view_label[i],
      module = module_size_df$module[i],
      module_short = module_size_df$module_short[i],
      n_features = n_i,
      y_top = top_y + 0.35,
      y_bottom = bottom_y - 0.35,
      header_y = top_y + 0.78,
      center_y = mean(c(top_y, bottom_y)),
      stringsAsFactors = FALSE
    )
    current_top <- current_top + n_i + gap_between_modules
  }
  module_layout <- dplyr::bind_rows(module_rows)

  feature_nodes <- feature_nodes %>%
    dplyr::left_join(module_layout, by = c("view_order", "view_label", "module", "module_short")) %>%
    dplyr::group_by(view_order, module) %>%
    dplyr::mutate(
      feature_row_index = dplyr::row_number(),
      y = y_top - 0.75 - (feature_row_index - 1),
      x_node = 2.45,
      x_label = 2.90
    ) %>%
    dplyr::ungroup()

  total_y_range <- range(feature_nodes$y)
  factor_nodes <- data.frame(
    factor = factor_set,
    factor_rank = seq_along(factor_set),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      y = seq(from = max(total_y_range) - 0.3, to = min(total_y_range) + 0.3, length.out = length(factor_set)),
      x = -5.20,
      factor_label = factor
    )

  plot_edges <- network_edges %>%
    dplyr::left_join(factor_nodes %>% dplyr::select(factor, x_factor = x, y_factor = y), by = "factor") %>%
    dplyr::left_join(feature_nodes %>% dplyr::select(feature_node_id, x_feature = x_node, y_feature = y), by = "feature_node_id") %>%
    dplyr::group_by(factor) %>%
    dplyr::arrange(y_feature, .by_group = TRUE) %>%
    dplyr::mutate(start_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.09 * edge_spread) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(feature_node_id) %>%
    dplyr::arrange(y_factor, .by_group = TRUE) %>%
    dplyr::mutate(end_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.08 * edge_spread) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(x_start = x_factor + 0.55, y_start = y_factor + start_offset, x_end = x_feature - 0.20, y_end = y_feature + end_offset)

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = module_layout,
      ggplot2::aes(xmin = 1.55, xmax = 2.70, ymin = y_bottom, ymax = y_top, fill = view_label),
      alpha = 0.10,
      color = "grey72",
      linewidth = 0.30
    ) +
    ggplot2::geom_segment(
      data = plot_edges,
      ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = loading_sign, linewidth = loading_strength),
      alpha = 0.30,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      data = factor_nodes,
      ggplot2::aes(x = x, y = y),
      shape = 21,
      size = 10.0,
      stroke = 0.70,
      fill = "white",
      color = "grey20"
    ) +
    ggplot2::geom_text(
      data = factor_nodes,
      ggplot2::aes(x = x, y = y, label = factor_label),
      size = 2.8,
      fontface = "plain"
    ) +
    ggplot2::geom_label(
      data = module_layout,
      ggplot2::aes(x = 1.72, y = header_y, label = module_short, fill = view_label),
      color = "grey15",
      label.size = 0.16,
      size = 2.35,
      label.padding = grid::unit(0.08, "lines"),
      fontface = "plain",
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = feature_nodes,
      ggplot2::aes(x = x_node, y = y, fill = view_label, size = maximum_loading),
      shape = 21,
      color = "grey20",
      stroke = 0.42
    ) +
    ggplot2::geom_point(
      data = feature_nodes %>% dplyr::filter(shared_feature),
      ggplot2::aes(x = x_node, y = y),
      shape = 21,
      size = 3.8,
      fill = NA,
      color = "#7A3E9D",
      stroke = 0.75
    ) +
    ggplot2::geom_text(
      data = feature_nodes %>% dplyr::filter(response_feature),
      ggplot2::aes(x = x_node, y = y, label = "*"),
      nudge_y = 0.24,
      size = 2.8,
      color = "black"
    )

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = feature_nodes,
      ggplot2::aes(x = x_label, y = y, label = feature_label_plotmath),
      parse = TRUE,
      size = 1.92,
      direction = "y",
      hjust = 0,
      nudge_x = 0.45,
      min.segment.length = 0,
      seed = 20261042,
      force = 1.0,
      max.overlaps = Inf,
      box.padding = 0.07,
      point.padding = 0.10,
      segment.alpha = 0.20,
      segment.size = 0.18
    )
  } else {
    p <- p + ggplot2::geom_text(
      data = feature_nodes,
      ggplot2::aes(x = x_label + 0.45, y = y, label = feature_label_plotmath),
      parse = TRUE,
      hjust = 0,
      size = 1.92
    )
  }

  p +
    ggplot2::scale_color_manual(values = c(Positive = "#D55E00", Negative = "#0072B2"), name = "Feature-weight sign") +
    ggplot2::scale_fill_manual(values = view_colors, name = "Feature group", guide = ggplot2::guide_legend(override.aes = list(shape = 21, size = 4.0, color = "grey20", alpha = 1))) +
    ggplot2::scale_linewidth_continuous(range = c(0.22, 0.76), guide = "none") +
    ggplot2::scale_size_continuous(range = c(2.0, 3.6), guide = "none") +
    ggplot2::coord_cartesian(xlim = c(-6.2, 7.4), ylim = c(min(feature_nodes$y) - 1.2, max(feature_nodes$y) + 1.2), clip = "off") +
    ggplot2::labs(
      title = plot_title,
      subtitle = plot_subtitle,
      caption = "Factors are shown as circles on the left; module-feature blocks are aligned on the right. Purple outlines denote features shared by multiple factors.",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_void(base_size = 9.0) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 7.6, color = "grey35"),
      plot.caption = ggplot2::element_text(size = 7.0, color = "grey35", hjust = 0),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(5, 180, 5, 5)
    )
}

# Main 5-factor network in the new layout.
p_mofa_v50_multifactor_feature_network <- mofa_v50_build_network_plot(
  factor_set = mofa_v50_network_factors,
  plot_title = "Multi-factor feature-weight network",
  plot_subtitle = "Factors are consolidated on the left and module-feature blocks are stacked on the right."
)
ggplot2::ggsave(
  file.path(mofa_v50_figure_dir, "MOFA_v50_multifactor_feature_loading_network.svg"),
  p_mofa_v50_multifactor_feature_network,
  width = 12.2,
  height = 10.2,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

# Selected top-ranked 3-factor network.
p_mofa_v50_selected_triple_feature_network <- mofa_v50_build_network_plot(
  factor_set = c(mofa_v50_selected_triple_diagnostics$factor_x[1], mofa_v50_selected_triple_diagnostics$factor_y[1], mofa_v50_selected_triple_diagnostics$factor_z[1]),
  plot_title = "Selected three-factor feature-weight network",
  plot_subtitle = "The same left-versus-right layout is used for the top-ranked Factor7-fixed combination."
)
ggplot2::ggsave(
  file.path(mofa_v50_figure_dir, "MOFA_v50_selected_triple_multifactor_feature_loading_network.svg"),
  p_mofa_v50_selected_triple_feature_network,
  width = 11.8,
  height = 9.6,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

# Requested fixed triplets in the same layout.
mofa_v50_requested_network_triplets <- list(
  c("Factor1", "Factor4", "Factor7"),
  c("Factor1", "Factor7", "Factor10"),
  c("Factor1", "Factor2", "Factor7"),
  c("Factor4", "Factor7", "Factor9"),
  c("Factor2", "Factor4", "Factor7")
)

mofa_v50_requested_network_plot_list <- vector("list", length(mofa_v50_requested_network_triplets))
for (i in seq_along(mofa_v50_requested_network_triplets)) {
  factor_set <- mofa_v50_requested_network_triplets[[i]]
  factor_label <- paste(factor_set, collapse = " + ")
  output_file <- file.path(mofa_v50_figure_dir, paste0("MOFA_v50_network_", paste(factor_set, collapse = "_"), ".svg"))
  message(sprintf("[v50 network %d/%d] %s", i, length(mofa_v50_requested_network_triplets), factor_label))
  plot_object <- mofa_v50_build_network_plot(
    factor_set = factor_set,
    plot_title = paste0("Feature-weight network: ", factor_label),
    plot_subtitle = "Left-stacked factor circles and right-aligned module-feature blocks are used to reduce clutter."
  )
  ggplot2::ggsave(output_file, plot_object, width = 11.8, height = 9.6, units = "in", device = svglite::svglite, bg = "white", limitsize = FALSE)
  mofa_v50_requested_network_plot_list[[i]] <- output_file
  rm(plot_object)
  invisible(gc())
}
write.csv(
  data.frame(
    factor_combination = vapply(mofa_v50_requested_network_triplets, paste, collapse = " + ", FUN.VALUE = character(1)),
    svg_file = unlist(mofa_v50_requested_network_plot_list),
    stringsAsFactors = FALSE
  ),
  file.path(mofa_v50_result_dir, "MOFA_v50_requested_network_files.csv"),
  row.names = FALSE
)

# Focused pretty versions for Factor1 + Factor4 + Factor7.
mofa_v50_f147_variants <- data.frame(
  variant_name = c("balanced", "wide_offset", "tight_offset"),
  edge_spread = c(1.0, 1.45, 0.75),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(mofa_v50_f147_variants))) {
  variant_plot <- mofa_v50_build_network_plot(
    factor_set = c("Factor1", "Factor4", "Factor7"),
    plot_title = paste0("Feature-weight network: Factor1 + Factor4 + Factor7 (", mofa_v50_f147_variants$variant_name[i], ")"),
    plot_subtitle = paste0("Variant with edge-spread multiplier = ", sprintf("%.2f", mofa_v50_f147_variants$edge_spread[i]), "."),
    edge_spread = mofa_v50_f147_variants$edge_spread[i]
  )
  ggplot2::ggsave(
    file.path(mofa_v50_figure_dir, paste0("MOFA_v50_network_Factor1_Factor4_Factor7_", mofa_v50_f147_variants$variant_name[i], ".svg")),
    variant_plot,
    width = 11.8,
    height = 9.6,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

#------------------------------#
# 22.9.3 Factor1+Factor4+Factor7 stand-alone 3D displays
#------------------------------#

mofa_v50_find_triple_row <- function(target_factors) {
  candidates <- mofa_v50_triple_permanova_atlas %>%
    dplyr::rowwise() %>%
    dplyr::mutate(match_target = setequal(c(factor_x, factor_y, factor_z), target_factors)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(match_target)
  if (nrow(candidates) == 0) return(mofa_v50_selected_triple_diagnostics[1, , drop = FALSE])
  candidates[1, , drop = FALSE]
}

mofa_v50_draw_three_factor_focus <- function(triple_factors, diagnostic_row, point_cex = 1.95, add_centroids = TRUE, draw_all_plane_lines = TRUE) {
  if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
    plot.new(); text(0.5, 0.5, "Package 'scatterplot3d' is required.")
    return(invisible(NULL))
  }
  plot_data <- mofa_v50_subject_mean_factor_data %>%
    dplyr::select(SubjectID, TRG_plot, dplyr::all_of(triple_factors)) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(triple_factors), is.finite), !is.na(TRG_plot)) %>%
    dplyr::mutate(TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR")))

  score_matrix <- scale(as.matrix(plot_data[, triple_factors, drop = FALSE]))
  score_matrix[!is.finite(score_matrix)] <- 0
  plot_data$x <- score_matrix[, 1]
  plot_data$y <- score_matrix[, 2]
  plot_data$z <- score_matrix[, 3]
  plot_data$color <- unname(mofa_response_colors[as.character(plot_data$TRG_plot)])
  plot_data$line_color <- grDevices::adjustcolor(plot_data$color, alpha.f = 0.10)
  plot_data$drop_color <- grDevices::adjustcolor(plot_data$color, alpha.f = 0.20)

  s3d <- scatterplot3d::scatterplot3d(
    x = plot_data$x,
    y = plot_data$y,
    z = plot_data$z,
    pch = 16,
    color = plot_data$color,
    cex.symbols = point_cex,
    type = "h",
    lty.hplot = 1,
    mar = c(2.4, 2.5, 3.6, 5.5),
    main = paste(triple_factors, collapse = " + "),
    xlab = mofa_v50_factor_axis_label(triple_factors[1]),
    ylab = mofa_v50_factor_axis_label(triple_factors[2]),
    zlab = mofa_v50_factor_axis_label(triple_factors[3]),
    angle = 52,
    scale.y = 1.0,
    box = TRUE,
    grid = TRUE
  )

  x_wall <- min(plot_data$x, na.rm = TRUE)
  y_wall <- max(plot_data$y, na.rm = TRUE)
  z_floor <- min(plot_data$z, na.rm = TRUE)
  for (i in seq_len(nrow(plot_data))) {
    xy_point <- s3d$xyz.convert(plot_data$x[i], plot_data$y[i], plot_data$z[i])
    xy_floor <- s3d$xyz.convert(plot_data$x[i], plot_data$y[i], z_floor)
    graphics::segments(xy_floor$x, xy_floor$y, xy_point$x, xy_point$y, col = plot_data$drop_color[i], lwd = 0.85)
    if (draw_all_plane_lines) {
      xy_xwall <- s3d$xyz.convert(x_wall, plot_data$y[i], plot_data$z[i])
      xy_ywall <- s3d$xyz.convert(plot_data$x[i], y_wall, plot_data$z[i])
      graphics::segments(xy_xwall$x, xy_xwall$y, xy_point$x, xy_point$y, col = plot_data$line_color[i], lwd = 0.60)
      graphics::segments(xy_ywall$x, xy_ywall$y, xy_point$x, xy_point$y, col = plot_data$line_color[i], lwd = 0.60)
    }
  }
  point_xy <- s3d$xyz.convert(plot_data$x, plot_data$y, plot_data$z)
  graphics::points(point_xy$x, point_xy$y, pch = 21, bg = plot_data$color, col = "grey20", cex = point_cex)

  if (add_centroids) {
    centroid_df <- plot_data %>%
      dplyr::group_by(TRG_plot) %>%
      dplyr::summarise(x = mean(x), y = mean(y), z = mean(z), .groups = "drop") %>%
      dplyr::mutate(color = unname(mofa_response_colors[as.character(TRG_plot)]))
    centroid_xy <- s3d$xyz.convert(centroid_df$x, centroid_df$y, centroid_df$z)
    graphics::points(centroid_xy$x, centroid_xy$y, pch = 22, bg = grDevices::adjustcolor(centroid_df$color, alpha.f = 0.65), col = "black", cex = 2.35)
  }

  stats_label <- paste0(
    "PERMANOVA p = ", ifelse(is.na(diagnostic_row$permanova_p[1]), "NA", ifelse(diagnostic_row$permanova_p[1] < 0.001, "<0.001", sprintf("%.3f", diagnostic_row$permanova_p[1]))),
    " | R2 = ", ifelse(is.na(diagnostic_row$permanova_r2[1]), "NA", sprintf("%.2f", diagnostic_row$permanova_r2[1])),
    "\nEnergy p = ", ifelse(is.na(diagnostic_row$energy_p[1]), "NA", ifelse(diagnostic_row$energy_p[1] < 0.001, "<0.001", sprintf("%.3f", diagnostic_row$energy_p[1]))),
    " | E = ", ifelse(is.na(diagnostic_row$energy_statistic[1]), "NA", sprintf("%.2f", diagnostic_row$energy_statistic[1]))
  )
  graphics::mtext(stats_label, side = 3, line = 0.35, adj = 0.02, cex = 0.74)
  graphics::legend(
    x = grconvertX(1.05, from = "npc", to = "user"),
    y = grconvertY(0.94, from = "npc", to = "user"),
    legend = c("pCR", "non-pCR", "group centroid"),
    pt.bg = c(unname(mofa_response_colors[c("pCR", "non_pCR")]), "grey75"),
    pch = c(21, 21, 22),
    col = c("grey20", "grey20", "black"),
    pt.cex = c(1.5, 1.5, 1.7),
    cex = 0.82,
    bty = "n",
    xpd = NA,
    xjust = 0,
    yjust = 1
  )
}

mofa_v50_f147_row <- mofa_v50_find_triple_row(c("Factor1", "Factor4", "Factor7"))
cat(
  sprintf(
    "Energy-distance test (Factor1 + Factor4 + Factor7; 999 permutations): E = %.6f, p = %.4f\n",
    mofa_v50_f147_row$energy_statistic[1],
    mofa_v50_f147_row$energy_p[1]
  )
)
svglite::svglite(file.path(mofa_v50_figure_dir, "MOFA_v50_three_factor_response_Factor1_Factor4_Factor7.svg"), width = 6.9, height = 6.2, bg = "white")
graphics::par(mar = c(2.1, 2.2, 3.3, 5.0), xpd = NA)
mofa_v50_draw_three_factor_focus(c("Factor1", "Factor4", "Factor7"), mofa_v50_f147_row, point_cex = 2.00, add_centroids = TRUE, draw_all_plane_lines = TRUE)
grDevices::dev.off()

# Factor1+Factor4+Factor7 with representative feature arrows.
mofa_v50_select_f147_feature_vectors <- function(n_features = 4) {
  raw_tbl <- mofa_feature_weights %>%
    dplyr::filter(factor %in% c("Factor1", "Factor4", "Factor7"), display_eligible, is.finite(weight_within_view))
  summary_tbl <- raw_tbl %>%
    dplyr::group_by(view, feature, feature_label) %>%
    dplyr::summarise(
      n_factors = dplyr::n_distinct(factor),
      best_p = if (any(is.finite(p_value))) min(p_value, na.rm = TRUE) else 1,
      max_abs_loading = max(abs(weight_within_view), na.rm = TRUE),
      combined_score = max_abs_loading + pmax(0, -log10(best_p + 1e-12)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(best_p, dplyr::desc(n_factors), dplyr::desc(combined_score))
  if (sum(summary_tbl$n_factors >= 2) >= n_features) {
    summary_tbl <- summary_tbl %>% dplyr::filter(n_factors >= 2)
  }
  selected_tbl <- summary_tbl %>% dplyr::slice_head(n = n_features)
  raw_tbl %>%
    dplyr::semi_join(selected_tbl, by = c("view", "feature", "feature_label")) %>%
    dplyr::select(view, feature, feature_label, factor, weight_within_view) %>%
    tidyr::pivot_wider(names_from = factor, values_from = weight_within_view, values_fill = 0) %>%
    dplyr::mutate(
      arrow_length = sqrt(Factor1^2 + Factor4^2 + Factor7^2),
      scale_factor = ifelse(max(arrow_length, na.rm = TRUE) > 0, 2.1 / max(arrow_length, na.rm = TRUE), 1),
      x = Factor1 * scale_factor,
      y = Factor4 * scale_factor,
      z = Factor7 * scale_factor
    )
}

mofa_v50_draw_f147_feature_arrow_plot <- function(feature_vector_tbl, point_cex = 1.9) {
  if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
    plot.new(); text(0.5, 0.5, "Package 'scatterplot3d' is required.")
    return(invisible(NULL))
  }
  plot_data <- mofa_v50_subject_mean_factor_data %>%
    dplyr::select(SubjectID, TRG_plot, Factor1, Factor4, Factor7) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(c("Factor1", "Factor4", "Factor7")), is.finite), !is.na(TRG_plot)) %>%
    dplyr::mutate(TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR")))
  score_matrix <- scale(as.matrix(plot_data[, c("Factor1", "Factor4", "Factor7"), drop = FALSE]))
  score_matrix[!is.finite(score_matrix)] <- 0
  plot_data$x <- score_matrix[, 1]; plot_data$y <- score_matrix[, 2]; plot_data$z <- score_matrix[, 3]
  plot_data$color <- unname(mofa_response_colors[as.character(plot_data$TRG_plot)])

  s3d <- scatterplot3d::scatterplot3d(
    x = plot_data$x, y = plot_data$y, z = plot_data$z,
    pch = 16, color = plot_data$color, cex.symbols = point_cex,
    type = "h", lty.hplot = 1,
    mar = c(2.2, 2.3, 3.4, 5.6),
    main = "Factor1 + Factor4 + Factor7 with representative feature arrows",
    xlab = mofa_v50_factor_axis_label("Factor1"),
    ylab = mofa_v50_factor_axis_label("Factor4"),
    zlab = mofa_v50_factor_axis_label("Factor7"),
    angle = 52, scale.y = 1.0, box = TRUE, grid = TRUE
  )
  point_xy <- s3d$xyz.convert(plot_data$x, plot_data$y, plot_data$z)
  graphics::points(point_xy$x, point_xy$y, pch = 21, bg = plot_data$color, col = "grey20", cex = point_cex)

  if (nrow(feature_vector_tbl) > 0) {
    feature_vector_tbl$arrow_col <- unname(view_colors[feature_vector_tbl$view])
    for (i in seq_len(nrow(feature_vector_tbl))) {
      p0 <- s3d$xyz.convert(0, 0, 0)
      p1 <- s3d$xyz.convert(feature_vector_tbl$x[i], feature_vector_tbl$y[i], feature_vector_tbl$z[i])
      graphics::arrows(p0$x, p0$y, p1$x, p1$y, length = 0.08, lwd = 1.2, col = feature_vector_tbl$arrow_col[i])
      graphics::text(p1$x, p1$y, labels = feature_vector_tbl$feature_label[i], pos = 4, cex = 0.62, col = "black")
    }
  }
  graphics::legend(
    x = grconvertX(1.05, from = "npc", to = "user"),
    y = grconvertY(0.94, from = "npc", to = "user"),
    legend = c("pCR", "non-pCR", unique(feature_vector_tbl$view)),
    pt.bg = c(unname(mofa_response_colors[c("pCR", "non_pCR")]), rep(NA, length(unique(feature_vector_tbl$view)))),
    pch = c(21, 21, rep(NA, length(unique(feature_vector_tbl$view)))),
    col = c("grey20", "grey20", unname(view_colors[unique(feature_vector_tbl$view)])),
    lty = c(NA, NA, rep(1, length(unique(feature_vector_tbl$view)))),
    lwd = c(NA, NA, rep(1.2, length(unique(feature_vector_tbl$view)))),
    pt.cex = 1.3,
    cex = 0.78,
    bty = "n",
    xpd = NA,
    xjust = 0,
    yjust = 1
  )
}

mofa_v50_f147_feature_vectors <- mofa_v50_select_f147_feature_vectors(4)
svglite::svglite(file.path(mofa_v50_figure_dir, "MOFA_v50_three_factor_response_Factor1_Factor4_Factor7_feature_arrows.svg"), width = 7.5, height = 6.4, bg = "white")
graphics::par(mar = c(2.2, 2.3, 3.3, 5.3), xpd = NA)
mofa_v50_draw_f147_feature_arrow_plot(mofa_v50_f147_feature_vectors, point_cex = 1.95)
grDevices::dev.off()

#------------------------------#
# 22.9.4 More compact atlas re-save using the refined enrichment panel
#------------------------------#

p_mofa_v50_factor_interpretation_atlas <-
  (p_mofa_v50_variance_heatmap_atlas | p_mofa_v50_response_enrichment_atlas | p_mofa_v50_association_delta_heatmap) +
  patchwork::plot_layout(widths = c(0.74, 0.56, 0.46), guides = "collect") +
  patchwork::plot_annotation(
    title = "MOFA factor interpretation atlas",
    subtitle = paste0(
      "Variance explained provides the omics context, response enrichment is the primary clinical summary, ",
      "and the right heatmap highlights paired-change associations."
    )
  ) &
  ggplot2::theme(legend.position = "right")

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_variance_association_atlas.svg"),
  plot = p_mofa_v50_factor_interpretation_atlas,
  width = 8.4,
  height = max(5.0, 2.5 + 0.28 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)


#-----------------------------------------------------------------#
# 22.9 v50 publication refinements: atlas sizing, cleaned networks, and main 3D panels
#-----------------------------------------------------------------#

#------------------------------#
# 22.9.1 Atlas refinement
#------------------------------#

mofa_v50_response_size_bucket <- function(x) {
  x <- ifelse(is.finite(x), x, 0)
  as.character(pmin(2, pmax(0, round(x))))
}

mofa_v50_response_enrichment_data <- mofa_v50_response_enrichment_data %>%
  dplyr::mutate(
    minus_log10_p_bucket = factor(
      mofa_v50_response_size_bucket(minus_log10_p),
      levels = c("0", "1", "2")
    )
  )

p_mofa_v50_variance_heatmap_atlas <- ggplot2::ggplot(
  mofa_v50_variance_plot_data,
  ggplot2::aes(x = view_label, y = factor, fill = r2)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.32) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(r2 >= 0.005, sprintf("%.1f", 100 * r2), "")),
    size = 2.65
  ) +
  ggplot2::scale_fill_gradient(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),low = "white", high = "#2166AC", name = "Variance\nexplained") +
  ggplot2::labs(title = "Variance explained by factor and omics view", subtitle = "Cell labels are percentages.", x = NULL, y = NULL) +
  ggplot2::theme_classic(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.1, color = "grey35"),
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, size = 7.7),
    axis.text.y = ggplot2::element_text(face = "plain", size = 7.8),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 7.1),
    legend.text = ggplot2::element_text(size = 6.8),
    legend.key.height = grid::unit(0.32, "in")
  )

p_mofa_v50_response_enrichment_atlas <- ggplot2::ggplot(
  mofa_v50_response_enrichment_data,
  ggplot2::aes(
    x = standardized_effect,
    y = factor,
    color = response_direction,
    fill = response_direction,
    size = minus_log10_p_bucket
  )
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.38) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = standardized_effect, yend = factor),
    linewidth = 0.58,
    alpha = 0.74,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(shape = 21, color = "grey20", stroke = 0.48) +
  ggplot2::geom_text(
    ggplot2::aes(label = significance_label),
    nudge_y = 0.20,
    size = 2.5,
    color = "black",
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~ Timepoint_display, nrow = 1) +
  ggplot2::scale_color_manual(
    values = mofa_response_colors[c("pCR", "non_pCR")],
    labels = c(pCR = "pCR-enriched", non_pCR = "non-pCR-enriched"),
    name = "Direction"
  ) +
  ggplot2::scale_fill_manual(values = mofa_response_colors[c("pCR", "non_pCR")], guide = "none") +
  ggplot2::scale_size_manual(
    values = c("0" = 1.7, "1" = 3.3, "2" = 4.9),
    name = expression(-log[10](italic(p))),
    drop = FALSE
  ) +
  ggplot2::labs(title = "pCR versus non-pCR factor enrichment", subtitle = "Hedges' g is oriented as pCR minus non-pCR.", x = "Response enrichment (Hedges' g)", y = NULL) +
  ggplot2::theme_classic(base_size = 8.9) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.0, color = "grey35"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "plain", size = 7.2),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 7.0),
    legend.text = ggplot2::element_text(size = 6.7),
    legend.key.height = grid::unit(0.26, "in"),
    plot.margin = ggplot2::margin(5, 8, 5, 1)
  )

p_mofa_v50_association_delta_heatmap <- ggplot2::ggplot(
  mofa_v50_association_delta_data,
  ggplot2::aes(x = comparison_label, y = factor, fill = minus_log10_p)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.30) +
  ggplot2::geom_text(ggplot2::aes(label = effect_label), size = 2.65) +
  ggplot2::scale_fill_gradient(
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128),
    low = "white",
    high = "black",
    limits = c(0, mofa_v50_association_p_upper),
    breaks = mofa_v50_coarse_neglog10_breaks(mofa_v50_association_p_upper),
    oob = scales::squish,
    name = expression(-log[10](italic(p)))
  ) +
  ggplot2::labs(title = "Paired-change associations", subtitle = "Text is the signed median contrast.", x = NULL, y = NULL) +
  ggplot2::theme_classic(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.1, color = "grey35"),
    axis.text.x = ggplot2::element_text(size = 7.7, face = "plain"),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 7.0),
    legend.text = ggplot2::element_text(size = 6.7),
    legend.key.height = grid::unit(0.32, "in")
  )

p_mofa_v50_factor_interpretation_atlas <-
  (p_mofa_v50_variance_heatmap_atlas | p_mofa_v50_response_enrichment_atlas | p_mofa_v50_association_delta_heatmap) +
  patchwork::plot_layout(widths = c(0.68, 0.52, 0.44), guides = "collect") +
  patchwork::plot_annotation(
    title = "MOFA factor interpretation atlas",
    subtitle = paste0(
      "Variance explained provides the omics context, response enrichment is the primary clinical summary, ",
      "and the right heatmap highlights paired-change associations."
    )
  ) &
  ggplot2::theme(legend.position = "right")

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_variance_association_atlas.svg"),
  plot = p_mofa_v50_factor_interpretation_atlas,
  width = 8.1,
  height = max(5.2, 2.6 + 0.30 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

p_mofa_v50_factor_response_variance_atlas <-
  (p_mofa_v50_variance_heatmap_atlas | p_mofa_v50_response_enrichment_atlas) +
  patchwork::plot_layout(widths = c(0.72, 0.58), guides = "collect") +
  patchwork::plot_annotation(
    title = "MOFA factor variance and response enrichment",
    subtitle = "Variance explained heatmap and the response-enrichment lollipop summary."
  ) &
  ggplot2::theme(legend.position = "right")

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_variance_atlas.svg"),
  plot = p_mofa_v50_factor_response_variance_atlas,
  width = 6.6,
  height = max(5.0, 2.5 + 0.30 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

#------------------------------#
# 22.9.2 Clean publication network layout
#------------------------------#

mofa_v50_build_network_plot <- function(
  factor_set,
  plot_title,
  plot_subtitle = NULL,
  edge_spread = 1.0,
  factor_fill = "#EFE7D3"
) {
  factor_set <- unique(factor_set)
  factor_set <- factor_set[factor_set %in% mofa_v50_factor_order]
  if (length(factor_set) < 2) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Not enough factors to build the network."))
  }

  active_factor_views <- mofa_variance_explained %>%
    dplyr::mutate(factor = as.character(factor), view = as.character(view)) %>%
    dplyr::filter(factor %in% factor_set, r2 >= mofa_active_view_r2) %>%
    dplyr::select(factor, view, view_r2 = r2)

  seed_edges <- dplyr::bind_rows(
    mofa_feature_weights %>%
      dplyr::filter(factor %in% factor_set, display_eligible) %>%
      dplyr::inner_join(active_factor_views, by = c("factor", "view")) %>%
      dplyr::group_by(factor, view, direction) %>%
      dplyr::slice_max(order_by = abs(weight_within_view), n = mofa_network_features_per_view_direction, with_ties = FALSE) %>%
      dplyr::ungroup(),
    mofa_feature_weights %>%
      dplyr::filter(factor %in% factor_set, display_eligible, is.finite(p_value), p_value < 0.05) %>%
      dplyr::inner_join(active_factor_views, by = c("factor", "view")) %>%
      dplyr::group_by(factor) %>%
      dplyr::arrange(p_value, dplyr::desc(abs(weight_within_view)), .by_group = TRUE) %>%
      dplyr::slice_head(n = 2) %>%
      dplyr::ungroup()
  ) %>%
    dplyr::distinct(factor, view, feature, .keep_all = TRUE) %>%
    dplyr::mutate(feature_node_id = paste(view, feature, sep = "::"), edge_id = paste(factor, feature_node_id, sep = "__"))

  network_edges <- mofa_feature_weights %>%
    dplyr::filter(factor %in% factor_set, display_eligible) %>%
    dplyr::inner_join(active_factor_views, by = c("factor", "view")) %>%
    dplyr::mutate(feature_node_id = paste(view, feature, sep = "::"), edge_id = paste(factor, feature_node_id, sep = "__")) %>%
    dplyr::filter(
      feature_node_id %in% unique(seed_edges$feature_node_id),
      abs(weight_within_view) >= mofa_network_shared_loading_threshold | edge_id %in% seed_edges$edge_id
    ) %>%
    dplyr::mutate(loading_sign = ifelse(weight >= 0, "Positive", "Negative"), loading_strength = abs(weight_within_view))

  if (nrow(network_edges) == 0) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No eligible edges."))
  }

  module_membership <- dplyr::bind_rows(lapply(split(network_edges, network_edges$view), function(view_edge_data) {
    profile <- view_edge_data %>%
      dplyr::select(feature_node_id, factor, weight_within_view) %>%
      tidyr::pivot_wider(names_from = factor, values_from = weight_within_view, values_fill = 0)
    for (factor_name in setdiff(factor_set, colnames(profile))) profile[[factor_name]] <- 0
    profile_matrix <- as.matrix(profile[, factor_set, drop = FALSE])
    rownames(profile_matrix) <- profile$feature_node_id
    module <- if (nrow(profile_matrix) <= 2) rep(1L, nrow(profile_matrix)) else stats::cutree(stats::hclust(stats::dist(profile_matrix), method = "ward.D2"), k = min(2L, max(1L, ceiling(nrow(profile_matrix) / 7))))
    data.frame(feature_node_id = rownames(profile_matrix), module = as.integer(module), stringsAsFactors = FALSE)
  }))

  feature_nodes <- network_edges %>%
    dplyr::group_by(feature_node_id, view, view_label, feature, feature_label) %>%
    dplyr::summarise(
      n_connected_factors = dplyr::n_distinct(factor),
      minimum_response_p = if (any(is.finite(p_value))) min(p_value, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::left_join(module_membership, by = "feature_node_id") %>%
    dplyr::mutate(
      view_order = match(view, required_views),
      shared_feature = n_connected_factors >= 2,
      response_feature = is.finite(minimum_response_p) & minimum_response_p < 0.05,
      feature_label_plotmath = mofa_v50_plotmath_feature_label(view, feature_label)
    ) %>%
    dplyr::arrange(view_order, module, dplyr::desc(n_connected_factors), feature_label)

  modules <- feature_nodes %>%
    dplyr::distinct(view, view_label, view_order, module) %>%
    dplyr::arrange(view_order, module) %>%
    dplyr::group_by(view) %>%
    dplyr::mutate(module_short = paste0("M", dplyr::row_number())) %>%
    dplyr::ungroup()

  feature_nodes <- feature_nodes %>%
    dplyr::left_join(modules, by = c("view", "view_label", "view_order", "module"))

  module_size_df <- feature_nodes %>%
    dplyr::group_by(view_order, view_label, module, module_short) %>%
    dplyr::summarise(n_features = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(view_order, module)

  gap_between_modules <- 0.55
  current_top <- 0
  module_rows <- vector("list", nrow(module_size_df))
  for (i in seq_len(nrow(module_size_df))) {
    n_i <- module_size_df$n_features[i]
    top_y <- -current_top
    bottom_y <- -(current_top + n_i - 1)
    module_rows[[i]] <- data.frame(
      view_order = module_size_df$view_order[i],
      view_label = module_size_df$view_label[i],
      module = module_size_df$module[i],
      module_short = module_size_df$module_short[i],
      n_features = n_i,
      y_top = top_y + 0.36,
      y_bottom = bottom_y - 0.36,
      header_y = top_y + 0.72,
      stringsAsFactors = FALSE
    )
    current_top <- current_top + n_i + gap_between_modules
  }
  module_layout <- dplyr::bind_rows(module_rows)

  feature_nodes <- feature_nodes %>%
    dplyr::left_join(module_layout, by = c("view_order", "view_label", "module", "module_short")) %>%
    dplyr::group_by(view_order, module) %>%
    dplyr::mutate(
      feature_row_index = dplyr::row_number(),
      y = y_top - 0.74 - (feature_row_index - 1),
      x_node = 2.25,
      x_label = 2.82
    ) %>%
    dplyr::ungroup()

  total_y_range <- range(feature_nodes$y)
  factor_nodes <- data.frame(factor = factor_set, stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      y = seq(from = max(total_y_range) - 0.2, to = min(total_y_range) + 0.2, length.out = length(factor_set)),
      x = -3.70,
      factor_label = factor
    )

  plot_edges <- network_edges %>%
    dplyr::left_join(factor_nodes %>% dplyr::select(factor, x_factor = x, y_factor = y), by = "factor") %>%
    dplyr::left_join(feature_nodes %>% dplyr::select(feature_node_id, x_feature = x_node, y_feature = y), by = "feature_node_id") %>%
    dplyr::group_by(factor) %>%
    dplyr::arrange(y_feature, .by_group = TRUE) %>%
    dplyr::mutate(start_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.08 * edge_spread) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(feature_node_id) %>%
    dplyr::arrange(y_factor, .by_group = TRUE) %>%
    dplyr::mutate(end_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.07 * edge_spread) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(x_start = x_factor + 0.48, y_start = y_factor + start_offset, x_end = x_feature - 0.10, y_end = y_feature + end_offset)

  ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = module_layout,
      ggplot2::aes(xmin = 1.98, xmax = 2.50, ymin = y_bottom, ymax = y_top, fill = view_label),
      alpha = 0.11,
      color = "grey72",
      linewidth = 0.30
    ) +
    ggplot2::geom_segment(
      data = plot_edges,
      ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = loading_sign, linewidth = loading_strength),
      alpha = 0.28,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      data = factor_nodes,
      ggplot2::aes(x = x, y = y),
      shape = 21,
      size = 11.5,
      stroke = 0.65,
      fill = factor_fill,
      color = "grey25"
    ) +
    ggplot2::geom_text(
      data = factor_nodes,
      ggplot2::aes(x = x, y = y, label = factor_label),
      size = 2.45,
      fontface = "plain"
    ) +
    ggplot2::geom_label(
      data = module_layout,
      ggplot2::aes(x = 2.06, y = header_y, label = module_short, fill = view_label),
      color = "grey15",
      label.size = 0.15,
      size = 2.25,
      label.padding = grid::unit(0.07, "lines"),
      fontface = "plain",
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = feature_nodes,
      ggplot2::aes(x = x_node, y = y, fill = view_label),
      shape = 21,
      size = 3.0,
      color = "grey20",
      stroke = 0.42
    ) +
    ggplot2::geom_point(
      data = feature_nodes %>% dplyr::filter(shared_feature),
      ggplot2::aes(x = x_node, y = y),
      shape = 21,
      size = 3.5,
      fill = NA,
      color = "#7A3E9D",
      stroke = 0.72
    ) +
    ggplot2::geom_text(
      data = feature_nodes %>% dplyr::filter(response_feature),
      ggplot2::aes(x = x_node, y = y, label = "*"),
      size = 2.1,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = feature_nodes,
      ggplot2::aes(x = x_label, y = y, label = feature_label_plotmath),
      parse = TRUE,
      hjust = 0,
      size = 1.92
    ) +
    ggplot2::scale_color_manual(values = c(Positive = "#D55E00", Negative = "#0072B2"), name = "Feature-weight sign") +
    ggplot2::scale_fill_manual(values = view_colors, name = "Feature group", guide = ggplot2::guide_legend(override.aes = list(shape = 21, size = 4.0, color = "grey20", alpha = 1))) +
    ggplot2::scale_linewidth_continuous(range = c(0.24, 0.74), guide = "none") +
    ggplot2::coord_cartesian(xlim = c(-4.8, 6.1), ylim = c(min(feature_nodes$y) - 0.8, max(feature_nodes$y) + 0.9), clip = "off") +
    ggplot2::labs(
      title = plot_title,
      subtitle = plot_subtitle,
      caption = "Factors are shown as labeled circles. Module boxes tightly enclose only the feature nodes; purple outlines denote features shared by multiple factors.",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_void(base_size = 8.9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 7.4, color = "grey35"),
      plot.caption = ggplot2::element_text(size = 6.8, color = "grey35", hjust = 0),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 7.0),
      legend.text = ggplot2::element_text(size = 6.8),
      plot.margin = ggplot2::margin(4, 150, 4, 4)
    )
}

# publication-ready cleaned networks
p_mofa_v50_network_f147 <- mofa_v50_build_network_plot(
  factor_set = c("Factor1", "Factor4", "Factor7"),
  plot_title = "Feature-weight network: Factor1 + Factor4 + Factor7",
  plot_subtitle = "Publication-style layout with tightly enclosed module boxes and direct factor-to-feature connectors."
)
ggplot2::ggsave(file.path(mofa_v50_figure_dir, "MOFA_v50_network_Factor1_Factor4_Factor7_publication.svg"), p_mofa_v50_network_f147, width = 10.6, height = 9.2, units = "in", device = svglite::svglite, bg = "white")

p_mofa_v50_network_f247 <- mofa_v50_build_network_plot(
  factor_set = c("Factor2", "Factor4", "Factor7"),
  plot_title = "Feature-weight network: Factor2 + Factor4 + Factor7",
  plot_subtitle = "Publication-style layout with tightly enclosed module boxes and direct factor-to-feature connectors."
)
ggplot2::ggsave(file.path(mofa_v50_figure_dir, "MOFA_v50_network_Factor2_Factor4_Factor7_publication.svg"), p_mofa_v50_network_f247, width = 10.6, height = 9.2, units = "in", device = svglite::svglite, bg = "white")

# also re-save the generic main network outputs using the cleaned layout
p_mofa_v50_multifactor_feature_network <- mofa_v50_build_network_plot(
  factor_set = mofa_v50_network_factors,
  plot_title = "Multi-factor feature-weight network",
  plot_subtitle = "Factors are shown on the left and module-feature stacks are aligned on the right to reduce clutter."
)
ggplot2::ggsave(file.path(mofa_v50_figure_dir, "MOFA_v50_multifactor_feature_loading_network.svg"), p_mofa_v50_multifactor_feature_network, width = 11.0, height = 10.2, units = "in", device = svglite::svglite, bg = "white")

p_mofa_v50_selected_triple_feature_network <- mofa_v50_build_network_plot(
  factor_set = c(mofa_v50_selected_triple_diagnostics$factor_x[1], mofa_v50_selected_triple_diagnostics$factor_y[1], mofa_v50_selected_triple_diagnostics$factor_z[1]),
  plot_title = "Selected three-factor feature-weight network",
  plot_subtitle = "Cleaned publication layout for the top-ranked Factor7-fixed combination."
)
ggplot2::ggsave(file.path(mofa_v50_figure_dir, "MOFA_v50_selected_triple_multifactor_feature_loading_network.svg"), p_mofa_v50_selected_triple_feature_network, width = 10.6, height = 9.1, units = "in", device = svglite::svglite, bg = "white")

#------------------------------#
# 22.9.3 Main 3D response plots for Factor1+4+7 and Factor2+4+7
#------------------------------#

mofa_v50_find_triple_row <- function(target_factors) {
  candidates <- mofa_v50_triple_permanova_atlas %>%
    dplyr::rowwise() %>%
    dplyr::mutate(match_target = setequal(c(factor_x, factor_y, factor_z), target_factors)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(match_target)
  if (nrow(candidates) == 0) return(mofa_v50_selected_triple_diagnostics[1, , drop = FALSE])
  candidates[1, , drop = FALSE]
}

mofa_v50_draw_cube <- function(s3d_obj, x, y, z, side = 0.22, fill_col = NA, border_col = "black", line_lwd = 1.0) {
  pts <- expand.grid(dx = c(-side, side), dy = c(-side, side), dz = c(-side, side))
  pts$x <- x + pts$dx; pts$y <- y + pts$dy; pts$z <- z + pts$dz
  xy <- s3d_obj$xyz.convert(pts$x, pts$y, pts$z)
  coords <- cbind(x = xy$x, y = xy$y)
  edge_idx <- list(
    c(1,2), c(1,3), c(1,5), c(2,4), c(2,6), c(3,4), c(3,7), c(4,8),
    c(5,6), c(5,7), c(6,8), c(7,8)
  )
  for (e in edge_idx) graphics::segments(coords[e[1],1], coords[e[1],2], coords[e[2],1], coords[e[2],2], col = border_col, lwd = line_lwd)
  center_xy <- s3d_obj$xyz.convert(x, y, z)
  graphics::points(center_xy$x, center_xy$y, pch = 22, bg = fill_col, col = border_col, cex = 1.55)
}

mofa_v50_draw_three_factor_publication <- function(triple_factors, diagnostic_row, point_cex = 1.95) {
  if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
    plot.new(); text(0.5, 0.5, "Package 'scatterplot3d' is required.")
    return(invisible(NULL))
  }
  plot_data <- mofa_v50_subject_mean_factor_data %>%
    dplyr::select(SubjectID, TRG_plot, dplyr::all_of(triple_factors)) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(triple_factors), is.finite), !is.na(TRG_plot)) %>%
    dplyr::mutate(TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR")))
  score_matrix <- scale(as.matrix(plot_data[, triple_factors, drop = FALSE]))
  score_matrix[!is.finite(score_matrix)] <- 0
  plot_data$x <- score_matrix[, 1]; plot_data$y <- score_matrix[, 2]; plot_data$z <- score_matrix[, 3]
  plot_data$color <- unname(mofa_response_colors[as.character(plot_data$TRG_plot)])

  s3d <- scatterplot3d::scatterplot3d(
    x = plot_data$x, y = plot_data$y, z = plot_data$z,
    pch = 16, color = plot_data$color, cex.symbols = point_cex,
    type = "p",
    mar = c(2.2, 2.3, 3.4, 5.3),
    main = paste(triple_factors, collapse = " + "),
    xlab = mofa_v50_factor_axis_label(triple_factors[1]),
    ylab = mofa_v50_factor_axis_label(triple_factors[2]),
    zlab = mofa_v50_factor_axis_label(triple_factors[3]),
    angle = 52, scale.y = 1.0, box = TRUE, grid = TRUE
  )
  point_xy <- s3d$xyz.convert(plot_data$x, plot_data$y, plot_data$z)
  graphics::points(point_xy$x, point_xy$y, pch = 21, bg = plot_data$color, col = "grey20", cex = point_cex)

  z_floor <- min(plot_data$z, na.rm = TRUE)
  centroid_df <- plot_data %>%
    dplyr::group_by(TRG_plot) %>%
    dplyr::summarise(x = mean(x), y = mean(y), z = mean(z), .groups = "drop") %>%
    dplyr::mutate(color = unname(mofa_response_colors[as.character(TRG_plot)]), drop_col = grDevices::adjustcolor(color, alpha.f = 0.42))

  for (i in seq_len(nrow(centroid_df))) {
    centroid_xy <- s3d$xyz.convert(centroid_df$x[i], centroid_df$y[i], centroid_df$z[i])
    floor_xy <- s3d$xyz.convert(centroid_df$x[i], centroid_df$y[i], z_floor)
    graphics::segments(floor_xy$x, floor_xy$y, centroid_xy$x, centroid_xy$y, col = centroid_df$drop_col[i], lwd = 1.25)
    mofa_v50_draw_cube(s3d, centroid_df$x[i], centroid_df$y[i], centroid_df$z[i], side = 0.18, fill_col = grDevices::adjustcolor(centroid_df$color[i], alpha.f = 0.70), border_col = "black", line_lwd = 1.0)
  }

  stats_label <- paste0(
    "PERMANOVA p = ", ifelse(is.na(diagnostic_row$permanova_p[1]), "NA", ifelse(diagnostic_row$permanova_p[1] < 0.001, "<0.001", sprintf("%.3f", diagnostic_row$permanova_p[1]))),
    " | R2 = ", ifelse(is.na(diagnostic_row$permanova_r2[1]), "NA", sprintf("%.2f", diagnostic_row$permanova_r2[1])),
    "\nEnergy p = ", ifelse(is.na(diagnostic_row$energy_p[1]), "NA", ifelse(diagnostic_row$energy_p[1] < 0.001, "<0.001", sprintf("%.3f", diagnostic_row$energy_p[1]))),
    " | E = ", ifelse(is.na(diagnostic_row$energy_statistic[1]), "NA", sprintf("%.2f", diagnostic_row$energy_statistic[1]))
  )
  graphics::mtext(stats_label, side = 3, line = 0.38, adj = 0.02, cex = 0.74)
  graphics::legend(
    x = grconvertX(1.05, from = "npc", to = "user"),
    y = grconvertY(0.94, from = "npc", to = "user"),
    legend = c("pCR", "non-pCR", "group centroid cube"),
    pt.bg = c(unname(mofa_response_colors[c("pCR", "non_pCR")]), "grey75"),
    pch = c(21, 21, 22),
    col = c("grey20", "grey20", "black"),
    pt.cex = c(1.5, 1.5, 1.6),
    cex = 0.80,
    bty = "n",
    xpd = NA,
    xjust = 0,
    yjust = 1
  )
}

mofa_v50_f147_row <- mofa_v50_find_triple_row(c("Factor1", "Factor4", "Factor7"))
svglite::svglite(file.path(mofa_v50_figure_dir, "MOFA_v50_three_factor_response_Factor1_Factor4_Factor7_publication.svg"), width = 6.9, height = 6.2, bg = "white")
graphics::par(mar = c(2.1, 2.2, 3.2, 5.1), xpd = NA)
mofa_v50_draw_three_factor_publication(c("Factor1", "Factor4", "Factor7"), mofa_v50_f147_row, point_cex = 2.05)
grDevices::dev.off()

mofa_v50_f247_row <- mofa_v50_find_triple_row(c("Factor2", "Factor4", "Factor7"))
svglite::svglite(file.path(mofa_v50_figure_dir, "MOFA_v50_three_factor_response_Factor2_Factor4_Factor7_publication.svg"), width = 6.9, height = 6.2, bg = "white")
graphics::par(mar = c(2.1, 2.2, 3.2, 5.1), xpd = NA)
mofa_v50_draw_three_factor_publication(c("Factor2", "Factor4", "Factor7"), mofa_v50_f247_row, point_cex = 2.05)
grDevices::dev.off()

p_mofa_v50_main_three_factor_publication <- NULL
svglite::svglite(file.path(mofa_v50_figure_dir, "MOFA_v50_three_factor_response_main_pair.svg"), width = 13.2, height = 6.3, bg = "white")
graphics::par(mfrow = c(1, 2), mar = c(2.1, 2.2, 3.2, 5.1), xpd = NA)
mofa_v50_draw_three_factor_publication(c("Factor1", "Factor4", "Factor7"), mofa_v50_f147_row, point_cex = 2.00)
mofa_v50_draw_three_factor_publication(c("Factor2", "Factor4", "Factor7"), mofa_v50_f247_row, point_cex = 2.00)
grDevices::dev.off()


#-----------------------------------------------------------------#
# 22.9 v50 publication refinements: atlas compatibility, main networks, and main 3D / feature-arrow panels
#-----------------------------------------------------------------#

#------------------------------#
# 22.9.1 Atlas refinement with Illustrator-friendly stepped legends
#------------------------------#

mofa_v50_response_size_bucket <- function(x) {
  x <- ifelse(is.finite(x), x, 0)
  factor(as.character(pmin(2, pmax(0, round(x)))), levels = c("0", "1", "2"))
}

mofa_v50_response_enrichment_data <- mofa_v50_response_enrichment_data %>%
  dplyr::mutate(minus_log10_p_bucket = mofa_v50_response_size_bucket(minus_log10_p))

mofa_v50_variance_upper <- max(mofa_v50_variance_plot_data$r2, na.rm = TRUE)
mofa_v50_variance_breaks <- unique(round(pretty(c(0, mofa_v50_variance_upper), n = 4), 3))
mofa_v50_variance_breaks <- mofa_v50_variance_breaks[mofa_v50_variance_breaks >= 0 & mofa_v50_variance_breaks <= mofa_v50_variance_upper + 1e-8]
if (length(mofa_v50_variance_breaks) < 3) mofa_v50_variance_breaks <- c(0, mofa_v50_variance_upper / 2, mofa_v50_variance_upper)

mofa_v50_assoc_breaks <- mofa_v50_coarse_neglog10_breaks(mofa_v50_association_p_upper)
if (length(mofa_v50_assoc_breaks) < 3) mofa_v50_assoc_breaks <- c(0, mofa_v50_association_p_upper / 2, mofa_v50_association_p_upper)

mofa_v50_enrich_xlim <- range(mofa_v50_response_enrichment_data$standardized_effect, finite = TRUE)
mofa_v50_enrich_xlim <- c(mofa_v50_enrich_xlim[1] - 0.18, mofa_v50_enrich_xlim[2] + 0.18)

p_mofa_v50_variance_heatmap_atlas <- ggplot2::ggplot(
  mofa_v50_variance_plot_data,
  ggplot2::aes(x = view_label, y = factor, fill = r2)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.32) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(r2 >= 0.005, sprintf("%.1f", 100 * r2), "")),
    size = 2.6
  ) +
  ggplot2::scale_fill_stepsn(
    colors = grDevices::colorRampPalette(c("white", "#2166AC"))(7),
    limits = c(0, mofa_v50_variance_upper),
    breaks = mofa_v50_variance_breaks,
    show.limits = TRUE,
    guide = ggplot2::guide_coloursteps(
      title = "Variance\nexplained",
      even.steps = FALSE,
      show.limits = TRUE,
      barheight = grid::unit(0.95, "in"),
      barwidth = grid::unit(0.20, "in")
    )
  ) +
  ggplot2::labs(
    title = "Variance explained by factor and omics view",
    subtitle = "Cell labels are percentages.",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.0, color = "grey35"),
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, size = 7.5),
    axis.text.y = ggplot2::element_text(face = "plain", size = 7.7),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 6.8),
    legend.text = ggplot2::element_text(size = 6.5),
    legend.key.height = grid::unit(0.22, "in"),
    legend.margin = ggplot2::margin(0, 0, 0, 0)
  )

p_mofa_v50_response_enrichment_atlas <- ggplot2::ggplot(
  mofa_v50_response_enrichment_data,
  ggplot2::aes(
    x = standardized_effect,
    y = factor,
    color = response_direction,
    fill = response_direction,
    size = minus_log10_p_bucket
  )
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.36) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = standardized_effect, yend = factor),
    linewidth = 0.56,
    alpha = 0.76,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(shape = 21, color = "grey20", stroke = 0.46) +
  ggplot2::geom_text(
    ggplot2::aes(label = significance_label),
    nudge_y = 0.18,
    size = 2.4,
    color = "black",
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~ Timepoint_display, nrow = 1) +
  ggplot2::scale_x_continuous(limits = mofa_v50_enrich_xlim, expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
  ggplot2::scale_color_manual(
    values = mofa_response_colors[c("pCR", "non_pCR")],
    labels = c(pCR = "pCR-enriched", non_pCR = "non-pCR-enriched"),
    name = "Direction"
  ) +
  ggplot2::scale_fill_manual(values = mofa_response_colors[c("pCR", "non_pCR")], guide = "none") +
  ggplot2::scale_size_manual(
    values = c("0" = 1.35, "1" = 2.45, "2" = 3.45),
    name = expression(-log[10](italic(p))),
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "pCR versus non-pCR factor enrichment",
    subtitle = "Hedges' g is oriented as pCR minus non-pCR.",
    x = "Response enrichment (Hedges' g)",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 8.8) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6.9, color = "grey35"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "plain", size = 7.0),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 6.8),
    legend.text = ggplot2::element_text(size = 6.4),
    legend.key.height = grid::unit(0.20, "in"),
    plot.margin = ggplot2::margin(5, 5, 5, 1)
  )

p_mofa_v50_association_delta_heatmap <- ggplot2::ggplot(
  mofa_v50_association_delta_data,
  ggplot2::aes(x = comparison_label, y = factor, fill = minus_log10_p)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.30) +
  ggplot2::geom_text(ggplot2::aes(label = effect_label), size = 2.6) +
  ggplot2::scale_fill_stepsn(
    colors = grDevices::colorRampPalette(c("white", "black"))(7),
    limits = c(0, mofa_v50_association_p_upper),
    breaks = mofa_v50_assoc_breaks,
    show.limits = TRUE,
    oob = scales::squish,
    guide = ggplot2::guide_coloursteps(
      title = expression(-log[10](italic(p))),
      even.steps = FALSE,
      show.limits = TRUE,
      barheight = grid::unit(0.95, "in"),
      barwidth = grid::unit(0.20, "in")
    )
  ) +
  ggplot2::labs(
    title = "Paired-change associations",
    subtitle = "Text is the signed median contrast.",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.0, color = "grey35"),
    axis.text.x = ggplot2::element_text(size = 7.5, face = "plain"),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 6.8),
    legend.text = ggplot2::element_text(size = 6.4),
    legend.key.height = grid::unit(0.22, "in")
  )

p_mofa_v50_factor_interpretation_atlas <-
  (p_mofa_v50_variance_heatmap_atlas | p_mofa_v50_response_enrichment_atlas | p_mofa_v50_association_delta_heatmap) +
  patchwork::plot_layout(widths = c(0.62, 0.42, 0.38), guides = "collect") +
  patchwork::plot_annotation(
    title = "MOFA factor interpretation atlas",
    subtitle = paste0(
      "Variance explained provides the omics context, response enrichment is the primary clinical summary, ",
      "and the right heatmap highlights paired-change associations."
    )
  ) &
  ggplot2::theme(legend.position = "right")

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_variance_association_atlas.svg"),
  plot = p_mofa_v50_factor_interpretation_atlas,
  width = 7.7,
  height = max(5.1, 2.6 + 0.30 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

p_mofa_v50_factor_response_variance_atlas <-
  (p_mofa_v50_variance_heatmap_atlas | p_mofa_v50_response_enrichment_atlas) +
  patchwork::plot_layout(widths = c(0.66, 0.44), guides = "collect") +
  patchwork::plot_annotation(
    title = "MOFA factor variance and response enrichment",
    subtitle = "Variance explained heatmap and the response-enrichment summary."
  ) &
  ggplot2::theme(legend.position = "right")

ggplot2::ggsave(
  filename = file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_variance_atlas.svg"),
  plot = p_mofa_v50_factor_response_variance_atlas,
  width = 6.1,
  height = max(4.9, 2.5 + 0.30 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

#------------------------------#
# 22.9.2 Factor color summary for networks
#------------------------------#

mofa_v50_factor_fill_summary <- mofa_v50_response_enrichment_data %>%
  dplyr::group_by(factor) %>%
  # Each factor has Baseline and After RT response-enrichment estimates.
  # Use the timepoint with the strongest finite Wilcoxon evidence to assign
  # the factor-circle color; break ties by the larger absolute Hedges' g.
  dplyr::arrange(
    dplyr::desc(is.finite(wilcoxon_p)),
    wilcoxon_p,
    dplyr::desc(abs(standardized_effect)),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    enrichment_label = dplyr::case_when(
      !is.finite(standardized_effect) ~ "neutral",
      standardized_effect >= 0 ~ "pCR-enriched",
      TRUE ~ "non-pCR-enriched"
    ),
    factor_fill = dplyr::case_when(
      enrichment_label == "pCR-enriched" ~ unname(mofa_response_colors[["pCR"]]),
      enrichment_label == "non-pCR-enriched" ~ unname(mofa_response_colors[["non_pCR"]]),
      TRUE ~ "#EFE7D3"
    )
  ) %>%
  dplyr::select(
    factor,
    Timepoint_display,
    standardized_effect,
    wilcoxon_p,
    enrichment_label,
    factor_fill
  )

mofa_v50_select_network_edges <- function(factor_set) {
  active_factor_views <- mofa_variance_explained %>%
    dplyr::mutate(factor = as.character(factor), view = as.character(view)) %>%
    dplyr::filter(factor %in% factor_set, r2 >= mofa_active_view_r2) %>%
    dplyr::select(factor, view, view_r2 = r2)

  seed_edges <- dplyr::bind_rows(
    mofa_feature_weights %>%
      dplyr::filter(factor %in% factor_set, display_eligible) %>%
      dplyr::inner_join(active_factor_views, by = c("factor", "view")) %>%
      dplyr::group_by(factor, view, direction) %>%
      dplyr::slice_max(order_by = abs(weight_within_view), n = mofa_network_features_per_view_direction, with_ties = FALSE) %>%
      dplyr::ungroup(),
    mofa_feature_weights %>%
      dplyr::filter(factor %in% factor_set, display_eligible, is.finite(p_value), p_value < 0.05) %>%
      dplyr::inner_join(active_factor_views, by = c("factor", "view")) %>%
      dplyr::group_by(factor) %>%
      dplyr::arrange(p_value, dplyr::desc(abs(weight_within_view)), .by_group = TRUE) %>%
      dplyr::slice_head(n = 2) %>%
      dplyr::ungroup()
  ) %>%
    dplyr::distinct(factor, view, feature, .keep_all = TRUE) %>%
    dplyr::mutate(feature_node_id = paste(view, feature, sep = "::"), edge_id = paste(factor, feature_node_id, sep = "__"))

  network_edges <- mofa_feature_weights %>%
    dplyr::filter(factor %in% factor_set, display_eligible) %>%
    dplyr::inner_join(active_factor_views, by = c("factor", "view")) %>%
    dplyr::mutate(feature_node_id = paste(view, feature, sep = "::"), edge_id = paste(factor, feature_node_id, sep = "__")) %>%
    dplyr::filter(
      feature_node_id %in% unique(seed_edges$feature_node_id),
      abs(weight_within_view) >= mofa_network_shared_loading_threshold | edge_id %in% seed_edges$edge_id
    ) %>%
    dplyr::mutate(loading_sign = ifelse(weight >= 0, "Positive", "Negative"), loading_strength = abs(weight_within_view))
  network_edges
}

mofa_v50_prepare_network_layout <- function(factor_set, orientation = c("vertical", "horizontal")) {
  orientation <- match.arg(orientation)
  factor_set <- unique(factor_set)
  factor_set <- factor_set[factor_set %in% mofa_v50_factor_order]
  network_edges <- mofa_v50_select_network_edges(factor_set)
  if (nrow(network_edges) == 0) return(NULL)

  module_membership <- dplyr::bind_rows(lapply(split(network_edges, network_edges$view), function(view_edge_data) {
    profile <- view_edge_data %>%
      dplyr::select(feature_node_id, factor, weight_within_view) %>%
      tidyr::pivot_wider(names_from = factor, values_from = weight_within_view, values_fill = 0)
    for (factor_name in setdiff(factor_set, colnames(profile))) profile[[factor_name]] <- 0
    profile_matrix <- as.matrix(profile[, factor_set, drop = FALSE])
    rownames(profile_matrix) <- profile$feature_node_id
    module <- if (nrow(profile_matrix) <= 2) rep(1L, nrow(profile_matrix)) else stats::cutree(stats::hclust(stats::dist(profile_matrix), method = "ward.D2"), k = min(2L, max(1L, ceiling(nrow(profile_matrix) / 7))))
    data.frame(feature_node_id = rownames(profile_matrix), module = as.integer(module), stringsAsFactors = FALSE)
  }))

  feature_nodes <- network_edges %>%
    dplyr::group_by(feature_node_id, view, view_label, feature, feature_label) %>%
    dplyr::summarise(
      n_connected_factors = dplyr::n_distinct(factor),
      minimum_response_p = if (any(is.finite(p_value))) min(p_value, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::left_join(module_membership, by = "feature_node_id") %>%
    dplyr::mutate(
      view_order = match(view, required_views),
      shared_feature = n_connected_factors >= 2,
      response_feature = is.finite(minimum_response_p) & minimum_response_p < 0.05,
      feature_label_plotmath = mofa_v50_plotmath_feature_label(view, feature_label)
    ) %>%
    dplyr::arrange(view_order, module, dplyr::desc(n_connected_factors), feature_label)

  modules <- feature_nodes %>%
    dplyr::distinct(view, view_label, view_order, module) %>%
    dplyr::arrange(view_order, module) %>%
    dplyr::group_by(view) %>%
    dplyr::mutate(module_short = paste0("M", dplyr::row_number())) %>%
    dplyr::ungroup()
  feature_nodes <- feature_nodes %>% dplyr::left_join(modules, by = c("view", "view_label", "view_order", "module"))

  module_size_df <- feature_nodes %>%
    dplyr::group_by(view_order, view_label, module, module_short) %>%
    dplyr::summarise(n_features = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(view_order, module)

  if (orientation == "vertical") {
    gap_between_modules <- 0.60
    current_top <- 0
    module_rows <- vector("list", nrow(module_size_df))
    for (i in seq_len(nrow(module_size_df))) {
      n_i <- module_size_df$n_features[i]
      y_values <- -(current_top + seq_len(n_i) - 1)
      module_rows[[i]] <- data.frame(
        view_order = module_size_df$view_order[i],
        view_label = module_size_df$view_label[i],
        module = module_size_df$module[i],
        module_short = module_size_df$module_short[i],
        x_node = 2.18,
        x_label = 2.72,
        stringsAsFactors = FALSE
      )
      feature_nodes$y[feature_nodes$view_order == module_size_df$view_order[i] & feature_nodes$module == module_size_df$module[i]] <- y_values
      feature_nodes$x_node[feature_nodes$view_order == module_size_df$view_order[i] & feature_nodes$module == module_size_df$module[i]] <- 2.18
      feature_nodes$x_label[feature_nodes$view_order == module_size_df$view_order[i] & feature_nodes$module == module_size_df$module[i]] <- 2.72
      current_top <- current_top + n_i + gap_between_modules
    }
    module_layout <- feature_nodes %>%
      dplyr::group_by(view_order, view_label, module, module_short) %>%
      dplyr::summarise(
        xmin = 1.94,
        xmax = 2.42,
        ymin = min(y) - 0.46,
        ymax = max(y) + 0.46,
        module_x = 2.03,
        module_y = max(y) + 0.62,
        .groups = "drop"
      )
    factor_nodes <- data.frame(factor = factor_set, stringsAsFactors = FALSE) %>%
      dplyr::mutate(
        y = seq(from = max(feature_nodes$y) - 0.1, to = min(feature_nodes$y) + 0.1, length.out = length(factor_set)),
        x = -3.55,
        x_anchor = -2.92,
        y_anchor = y
      )
    plot_edges <- network_edges %>%
      dplyr::left_join(factor_nodes %>% dplyr::select(factor, x_factor = x, x_anchor, y_factor = y, y_anchor), by = "factor") %>%
      dplyr::left_join(feature_nodes %>% dplyr::select(feature_node_id, x_feature = x_node, y_feature = y), by = "feature_node_id") %>%
      dplyr::group_by(factor) %>%
      dplyr::arrange(y_feature, .by_group = TRUE) %>%
      dplyr::mutate(start_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.075) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(feature_node_id) %>%
      dplyr::arrange(y_factor, .by_group = TRUE) %>%
      dplyr::mutate(end_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.060) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(x_start = x_anchor, y_start = y_anchor + start_offset, x_end = x_feature - 0.10, y_end = y_feature + end_offset)
    layout_spec <- list(
      factor_nodes = factor_nodes,
      feature_nodes = feature_nodes,
      module_layout = module_layout,
      plot_edges = plot_edges,
      xlim = c(-4.7, 6.0),
      ylim = c(min(feature_nodes$y) - 0.8, max(feature_nodes$y) + 1.0)
    )
  } else {
    gap_between_modules <- 0.90
    current_left <- 0
    for (i in seq_len(nrow(module_size_df))) {
      n_i <- module_size_df$n_features[i]
      y_values <- rev(seq_len(n_i))
      x_value <- current_left
      feature_nodes$y[feature_nodes$view_order == module_size_df$view_order[i] & feature_nodes$module == module_size_df$module[i]] <- y_values
      feature_nodes$x_node[feature_nodes$view_order == module_size_df$view_order[i] & feature_nodes$module == module_size_df$module[i]] <- x_value
      feature_nodes$x_label[feature_nodes$view_order == module_size_df$view_order[i] & feature_nodes$module == module_size_df$module[i]] <- x_value + 0.44
      current_left <- current_left + 1.45 + gap_between_modules
    }
    module_layout <- feature_nodes %>%
      dplyr::group_by(view_order, view_label, module, module_short) %>%
      dplyr::summarise(
        xmin = unique(x_node) - 0.26,
        xmax = unique(x_node) + 0.26,
        ymin = min(y) - 0.46,
        ymax = max(y) + 0.46,
        module_x = unique(x_node),
        module_y = max(y) + 0.68,
        .groups = "drop"
      )
    factor_nodes <- data.frame(factor = factor_set, stringsAsFactors = FALSE) %>%
      dplyr::mutate(
        x = seq(from = min(feature_nodes$x_node) - 0.6, to = max(feature_nodes$x_node) + 0.6, length.out = length(factor_set)),
        y = max(module_layout$module_y) + 1.45,
        x_anchor = x,
        y_anchor = y - 0.62
      )
    plot_edges <- network_edges %>%
      dplyr::left_join(factor_nodes %>% dplyr::select(factor, x_factor = x, x_anchor, y_factor = y, y_anchor), by = "factor") %>%
      dplyr::left_join(feature_nodes %>% dplyr::select(feature_node_id, x_feature = x_node, y_feature = y), by = "feature_node_id") %>%
      dplyr::group_by(factor) %>%
      dplyr::arrange(x_feature, .by_group = TRUE) %>%
      dplyr::mutate(start_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.090) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(feature_node_id) %>%
      dplyr::arrange(x_factor, .by_group = TRUE) %>%
      dplyr::mutate(end_offset = (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.060) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(x_start = x_anchor + start_offset, y_start = y_anchor, x_end = x_feature + end_offset, y_end = y_feature + 0.10)
    layout_spec <- list(
      factor_nodes = factor_nodes,
      feature_nodes = feature_nodes,
      module_layout = module_layout,
      plot_edges = plot_edges,
      xlim = c(min(module_layout$xmin) - 1.0, max(module_layout$xmax) + 4.4),
      ylim = c(min(feature_nodes$y) - 0.9, max(factor_nodes$y) + 0.8)
    )
  }

  factor_nodes <- layout_spec$factor_nodes %>%
    dplyr::left_join(mofa_v50_factor_fill_summary, by = "factor") %>%
    dplyr::mutate(
      enrichment_label = ifelse(is.na(enrichment_label), "neutral", enrichment_label),
      factor_fill = ifelse(is.na(factor_fill), "#EFE7D3", factor_fill)
    )
  layout_spec$factor_nodes <- factor_nodes
  layout_spec$network_edges <- network_edges
  layout_spec
}

mofa_v50_build_network_plot <- function(
  factor_set,
  plot_title,
  plot_subtitle = NULL,
  orientation = c("vertical", "horizontal"),
  module_label_angle = 0
) {
  orientation <- match.arg(orientation)
  layout_spec <- mofa_v50_prepare_network_layout(factor_set, orientation = orientation)
  if (is.null(layout_spec)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No eligible edges."))
  }

  factor_nodes <- layout_spec$factor_nodes
  feature_nodes <- layout_spec$feature_nodes
  module_layout <- layout_spec$module_layout
  plot_edges <- layout_spec$plot_edges
  xlim_values <- layout_spec$xlim
  ylim_values <- layout_spec$ylim

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = module_layout,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = view_label),
      alpha = 0.11,
      color = "grey72",
      linewidth = 0.30
    ) +
    ggplot2::geom_segment(
      data = plot_edges,
      ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = loading_sign, linewidth = loading_strength),
      alpha = 0.28,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      data = factor_nodes,
      ggplot2::aes(x = x, y = y, fill = enrichment_label),
      shape = 21,
      size = 11.2,
      stroke = 0.68,
      color = "grey25"
    ) +
    ggplot2::geom_text(
      data = factor_nodes,
      ggplot2::aes(x = x, y = y, label = factor),
      size = 2.35,
      fontface = "plain"
    ) +
    ggplot2::geom_label(
      data = module_layout,
      ggplot2::aes(x = module_x, y = module_y, label = module_short, fill = view_label),
      color = "grey15",
      label.size = 0.15,
      size = 2.20,
      label.padding = grid::unit(0.07, "lines"),
      fontface = "plain",
      angle = module_label_angle,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = feature_nodes,
      ggplot2::aes(x = x_node, y = y, fill = view_label),
      shape = 21,
      size = 3.0,
      color = "grey20",
      stroke = 0.42
    ) +
    ggplot2::geom_point(
      data = feature_nodes %>% dplyr::filter(shared_feature),
      ggplot2::aes(x = x_node, y = y),
      shape = 21,
      size = 3.45,
      fill = NA,
      color = "#7A3E9D",
      stroke = 0.72
    ) +
    ggplot2::geom_text(
      data = feature_nodes %>% dplyr::filter(response_feature),
      ggplot2::aes(x = x_node, y = y, label = "*"),
      size = 2.1,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = feature_nodes,
      ggplot2::aes(x = x_label, y = y, label = feature_label_plotmath),
      parse = TRUE,
      hjust = 0,
      size = 1.90
    ) +
    ggplot2::scale_color_manual(values = c(Positive = "#D55E00", Negative = "#0072B2"), name = "Feature-weight sign") +
    ggplot2::scale_fill_manual(
      values = c(view_colors, "pCR-enriched" = mofa_response_colors[["pCR"]], "non-pCR-enriched" = mofa_response_colors[["non_pCR"]], "neutral" = "#EFE7D3"),
      breaks = c("pCR-enriched", "non-pCR-enriched", names(view_colors)),
      labels = c("pCR-enriched factor", "non-pCR-enriched factor", names(view_colors)),
      name = NULL,
      guide = ggplot2::guide_legend(override.aes = list(shape = 21, size = 4.0, color = "grey25", alpha = 1))
    ) +
    ggplot2::scale_linewidth_continuous(range = c(0.24, 0.74), guide = "none") +
    ggplot2::coord_cartesian(xlim = xlim_values, ylim = ylim_values, clip = "off") +
    ggplot2::labs(
      title = plot_title,
      subtitle = plot_subtitle,
      caption = "Factor circles are colored by response-enrichment direction. Module boxes tightly enclose feature nodes, and purple outlines denote features shared by multiple factors.",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_void(base_size = 8.9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 7.3, color = "grey35"),
      plot.caption = ggplot2::element_text(size = 6.7, color = "grey35", hjust = 0),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 6.8),
      legend.text = ggplot2::element_text(size = 6.6),
      plot.margin = if (orientation == "vertical") ggplot2::margin(4, 145, 4, 4) else ggplot2::margin(4, 80, 4, 4)
    )
  p
}

mofa_v50_save_network_bundle <- function(factor_set, file_stub) {
  p_vertical <- mofa_v50_build_network_plot(
    factor_set = factor_set,
    plot_title = paste0("Feature-weight network: ", paste(factor_set, collapse = " + ")),
    plot_subtitle = "Vertical publication layout.",
    orientation = "vertical",
    module_label_angle = 0
  )
  ggplot2::ggsave(file.path(mofa_v50_figure_dir, paste0("MOFA_v50_network_", file_stub, "_vertical.svg")), p_vertical, width = 10.4, height = 9.1, units = "in", device = svglite::svglite, bg = "white")

  p_horizontal <- mofa_v50_build_network_plot(
    factor_set = factor_set,
    plot_title = paste0("Feature-weight network: ", paste(factor_set, collapse = " + ")),
    plot_subtitle = "Horizontal publication layout.",
    orientation = "horizontal",
    module_label_angle = 0
  )
  ggplot2::ggsave(file.path(mofa_v50_figure_dir, paste0("MOFA_v50_network_", file_stub, "_horizontal.svg")), p_horizontal, width = 14.8, height = 6.8, units = "in", device = svglite::svglite, bg = "white")

  p_horizontal_m90 <- mofa_v50_build_network_plot(
    factor_set = factor_set,
    plot_title = paste0("Feature-weight network: ", paste(factor_set, collapse = " + ")),
    plot_subtitle = "Horizontal layout with 90-degree module labels.",
    orientation = "horizontal",
    module_label_angle = 90
  )
  ggplot2::ggsave(file.path(mofa_v50_figure_dir, paste0("MOFA_v50_network_", file_stub, "_horizontal_module90.svg")), p_horizontal_m90, width = 14.8, height = 6.8, units = "in", device = svglite::svglite, bg = "white")
  invisible(list(vertical = p_vertical, horizontal = p_horizontal, horizontal_module90 = p_horizontal_m90))
}

# main network bundles
mofa_v50_network_bundle_f147 <- mofa_v50_save_network_bundle(c("Factor1", "Factor4", "Factor7"), "Factor1_Factor4_Factor7_publication")
mofa_v50_network_bundle_f247 <- mofa_v50_save_network_bundle(c("Factor2", "Factor4", "Factor7"), "Factor2_Factor4_Factor7_publication")

# representative generic outputs
p_mofa_v50_multifactor_feature_network <- mofa_v50_build_network_plot(
  factor_set = mofa_v50_network_factors,
  plot_title = "Multi-factor feature-weight network",
  plot_subtitle = "Vertical publication layout.",
  orientation = "vertical"
)
ggplot2::ggsave(file.path(mofa_v50_figure_dir, "MOFA_v50_multifactor_feature_loading_network.svg"), p_mofa_v50_multifactor_feature_network, width = 10.8, height = 9.9, units = "in", device = svglite::svglite, bg = "white")

#------------------------------#
# 22.9.3 Main 3D response plots and feature-arrow versions
#------------------------------#

mofa_v50_find_triple_row <- function(target_factors) {
  candidates <- mofa_v50_triple_permanova_atlas %>%
    dplyr::rowwise() %>%
    dplyr::mutate(match_target = setequal(c(factor_x, factor_y, factor_z), target_factors)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(match_target)
  if (nrow(candidates) == 0) return(mofa_v50_selected_triple_diagnostics[1, , drop = FALSE])
  candidates[1, , drop = FALSE]
}

mofa_v50_draw_filled_cube <- function(s3d_obj, x, y, z, side = 0.14, fill_col = "grey70") {
  v <- data.frame(
    x = x + c(-side, side, side, -side, -side, side, side, -side),
    y = y + c(-side, -side, side, side, -side, -side, side, side),
    z = z + c(-side, -side, -side, -side, side, side, side, side)
  )
  p <- s3d_obj$xyz.convert(v$x, v$y, v$z)
  coords <- cbind(p$x, p$y)
  face_fill <- grDevices::adjustcolor(fill_col, alpha.f = 0.72)
  face_fill_light <- grDevices::adjustcolor(fill_col, alpha.f = 0.54)
  face_fill_dark <- grDevices::adjustcolor(fill_col, alpha.f = 0.86)
  graphics::polygon(coords[c(1,2,6,5),1], coords[c(1,2,6,5),2], col = face_fill, border = "black", lwd = 0.7)
  graphics::polygon(coords[c(2,3,7,6),1], coords[c(2,3,7,6),2], col = face_fill_dark, border = "black", lwd = 0.7)
  graphics::polygon(coords[c(5,6,7,8),1], coords[c(5,6,7,8),2], col = face_fill_light, border = "black", lwd = 0.7)
}

mofa_v50_draw_three_factor_publication <- function(triple_factors, diagnostic_row, point_cex = 2.00) {
  if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
    plot.new(); text(0.5, 0.5, "Package 'scatterplot3d' is required.")
    return(invisible(NULL))
  }
  plot_data <- mofa_v50_subject_mean_factor_data %>%
    dplyr::select(SubjectID, TRG_plot, dplyr::all_of(triple_factors)) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(triple_factors), is.finite), !is.na(TRG_plot)) %>%
    dplyr::mutate(TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR")))
  score_matrix <- scale(as.matrix(plot_data[, triple_factors, drop = FALSE]))
  score_matrix[!is.finite(score_matrix)] <- 0
  plot_data$x <- score_matrix[, 1]; plot_data$y <- score_matrix[, 2]; plot_data$z <- score_matrix[, 3]
  plot_data$color <- unname(mofa_response_colors[as.character(plot_data$TRG_plot)])

  s3d <- scatterplot3d::scatterplot3d(
    x = plot_data$x, y = plot_data$y, z = plot_data$z,
    pch = 16, color = plot_data$color, cex.symbols = point_cex,
    type = "p",
    mar = c(2.2, 2.3, 3.3, 5.3),
    main = paste(triple_factors, collapse = " + "),
    xlab = mofa_v50_factor_axis_label(triple_factors[1]),
    ylab = mofa_v50_factor_axis_label(triple_factors[2]),
    zlab = mofa_v50_factor_axis_label(triple_factors[3]),
    angle = 52, scale.y = 1.0, box = TRUE, grid = TRUE
  )
  point_xy <- s3d$xyz.convert(plot_data$x, plot_data$y, plot_data$z)
  graphics::points(point_xy$x, point_xy$y, pch = 21, bg = plot_data$color, col = "grey20", cex = point_cex)

  z_floor <- min(plot_data$z, na.rm = TRUE)
  centroid_df <- plot_data %>%
    dplyr::group_by(TRG_plot) %>%
    dplyr::summarise(x = mean(x), y = mean(y), z = mean(z), .groups = "drop") %>%
    dplyr::mutate(color = unname(mofa_response_colors[as.character(TRG_plot)]), drop_col = grDevices::adjustcolor(color, alpha.f = 0.62))
  for (i in seq_len(nrow(centroid_df))) {
    centroid_xy <- s3d$xyz.convert(centroid_df$x[i], centroid_df$y[i], centroid_df$z[i])
    floor_xy <- s3d$xyz.convert(centroid_df$x[i], centroid_df$y[i], z_floor)
    graphics::segments(floor_xy$x, floor_xy$y, centroid_xy$x, centroid_xy$y, col = centroid_df$drop_col[i], lwd = 1.35)
    mofa_v50_draw_filled_cube(s3d, centroid_df$x[i], centroid_df$y[i], centroid_df$z[i], side = 0.11, fill_col = centroid_df$color[i])
  }

  stats_label <- paste0(
    "PERMANOVA p = ", ifelse(is.na(diagnostic_row$permanova_p[1]), "NA", ifelse(diagnostic_row$permanova_p[1] < 0.001, "<0.001", sprintf("%.3f", diagnostic_row$permanova_p[1]))),
    " | R2 = ", ifelse(is.na(diagnostic_row$permanova_r2[1]), "NA", sprintf("%.2f", diagnostic_row$permanova_r2[1])),
    "\nEnergy p = ", ifelse(is.na(diagnostic_row$energy_p[1]), "NA", ifelse(diagnostic_row$energy_p[1] < 0.001, "<0.001", sprintf("%.3f", diagnostic_row$energy_p[1]))),
    " | E = ", ifelse(is.na(diagnostic_row$energy_statistic[1]), "NA", sprintf("%.2f", diagnostic_row$energy_statistic[1]))
  )
  graphics::mtext(stats_label, side = 3, line = 0.38, adj = 0.02, cex = 0.74)
  graphics::legend(
    x = grconvertX(1.05, from = "npc", to = "user"),
    y = grconvertY(0.94, from = "npc", to = "user"),
    legend = c("pCR", "non-pCR", "group centroid cube"),
    pt.bg = c(unname(mofa_response_colors[c("pCR", "non_pCR")]), "grey75"),
    pch = c(21, 21, 22),
    col = c("grey20", "grey20", "black"),
    pt.cex = c(1.5, 1.5, 1.5),
    cex = 0.80,
    bty = "n",
    xpd = NA,
    xjust = 0,
    yjust = 1
  )
}

mofa_v50_select_feature_vectors <- function(target_factors, n_features = 4) {
  raw_tbl <- mofa_feature_weights %>%
    dplyr::filter(factor %in% target_factors, display_eligible, is.finite(weight_within_view))
  summary_tbl <- raw_tbl %>%
    dplyr::group_by(view, feature, feature_label) %>%
    dplyr::summarise(
      n_factors = dplyr::n_distinct(factor),
      best_p = if (any(is.finite(p_value))) min(p_value, na.rm = TRUE) else 1,
      max_abs_loading = max(abs(weight_within_view), na.rm = TRUE),
      combined_score = max_abs_loading + pmax(0, -log10(best_p + 1e-12)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(best_p, dplyr::desc(n_factors), dplyr::desc(combined_score))
  if (sum(summary_tbl$n_factors >= 2) >= min(2, n_features)) {
    summary_tbl <- summary_tbl %>% dplyr::filter(n_factors >= 2)
  }
  selected_tbl <- summary_tbl %>% dplyr::slice_head(n = n_features)
  wide_tbl <- raw_tbl %>%
    dplyr::semi_join(selected_tbl, by = c("view", "feature", "feature_label")) %>%
    dplyr::select(view, feature, feature_label, factor, weight_within_view) %>%
    tidyr::pivot_wider(names_from = factor, values_from = weight_within_view, values_fill = 0)
  for (factor_name in target_factors) if (!(factor_name %in% colnames(wide_tbl))) wide_tbl[[factor_name]] <- 0
  vec_tbl <- wide_tbl %>%
    dplyr::mutate(
      arrow_length = sqrt(.data[[target_factors[1]]] ^ 2 + .data[[target_factors[2]]] ^ 2 + .data[[target_factors[3]]] ^ 2),
      scale_factor = ifelse(max(arrow_length, na.rm = TRUE) > 0, 2.0 / max(arrow_length, na.rm = TRUE), 1),
      x = .data[[target_factors[1]]] * scale_factor,
      y = .data[[target_factors[2]]] * scale_factor,
      z = .data[[target_factors[3]]] * scale_factor
    )
  vec_tbl
}

mofa_v50_draw_feature_arrow_plot <- function(target_factors, feature_vector_tbl, point_cex = 1.95) {
  if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
    plot.new(); text(0.5, 0.5, "Package 'scatterplot3d' is required.")
    return(invisible(NULL))
  }
  plot_data <- mofa_v50_subject_mean_factor_data %>%
    dplyr::select(SubjectID, TRG_plot, dplyr::all_of(target_factors)) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(target_factors), is.finite), !is.na(TRG_plot)) %>%
    dplyr::mutate(TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR")))
  score_matrix <- scale(as.matrix(plot_data[, target_factors, drop = FALSE]))
  score_matrix[!is.finite(score_matrix)] <- 0
  plot_data$x <- score_matrix[, 1]; plot_data$y <- score_matrix[, 2]; plot_data$z <- score_matrix[, 3]
  plot_data$color <- unname(mofa_response_colors[as.character(plot_data$TRG_plot)])
  s3d <- scatterplot3d::scatterplot3d(
    x = plot_data$x, y = plot_data$y, z = plot_data$z,
    pch = 16, color = plot_data$color, cex.symbols = point_cex,
    type = "p", mar = c(2.2, 2.3, 3.3, 5.5),
    main = paste0(paste(target_factors, collapse = " + "), " with representative feature arrows"),
    xlab = mofa_v50_factor_axis_label(target_factors[1]),
    ylab = mofa_v50_factor_axis_label(target_factors[2]),
    zlab = mofa_v50_factor_axis_label(target_factors[3]),
    angle = 52, scale.y = 1.0, box = TRUE, grid = TRUE
  )
  point_xy <- s3d$xyz.convert(plot_data$x, plot_data$y, plot_data$z)
  graphics::points(point_xy$x, point_xy$y, pch = 21, bg = plot_data$color, col = "grey20", cex = point_cex)
  if (nrow(feature_vector_tbl) > 0) {
    feature_vector_tbl$arrow_col <- unname(view_colors[feature_vector_tbl$view])
    for (i in seq_len(nrow(feature_vector_tbl))) {
      p0 <- s3d$xyz.convert(0, 0, 0)
      p1 <- s3d$xyz.convert(feature_vector_tbl$x[i], feature_vector_tbl$y[i], feature_vector_tbl$z[i])
      graphics::arrows(p0$x, p0$y, p1$x, p1$y, length = 0.08, lwd = 1.2, col = feature_vector_tbl$arrow_col[i])
      graphics::text(p1$x, p1$y, labels = feature_vector_tbl$feature_label[i], pos = 4, cex = 0.62, col = "black")
    }
  }
  graphics::legend(
    x = grconvertX(1.05, from = "npc", to = "user"),
    y = grconvertY(0.94, from = "npc", to = "user"),
    legend = c("pCR", "non-pCR", unique(feature_vector_tbl$view)),
    pt.bg = c(unname(mofa_response_colors[c("pCR", "non_pCR")]), rep(NA, length(unique(feature_vector_tbl$view)))),
    pch = c(21, 21, rep(NA, length(unique(feature_vector_tbl$view)))),
    col = c("grey20", "grey20", unname(view_colors[unique(feature_vector_tbl$view)])),
    lty = c(NA, NA, rep(1, length(unique(feature_vector_tbl$view)))),
    lwd = c(NA, NA, rep(1.2, length(unique(feature_vector_tbl$view)))),
    pt.cex = 1.3,
    cex = 0.78,
    bty = "n",
    xpd = NA,
    xjust = 0,
    yjust = 1
  )
}

mofa_v50_save_3d_bundle <- function(target_factors, file_stub) {
  diag_row <- mofa_v50_find_triple_row(target_factors)
  svglite::svglite(file.path(mofa_v50_figure_dir, paste0("MOFA_v50_three_factor_response_", file_stub, "_publication.svg")), width = 6.9, height = 6.2, bg = "white")
  graphics::par(mar = c(2.1, 2.2, 3.2, 5.1), xpd = NA)
  mofa_v50_draw_three_factor_publication(target_factors, diag_row, point_cex = 2.05)
  grDevices::dev.off()

  vec_tbl <- mofa_v50_select_feature_vectors(target_factors, n_features = 4)
  svglite::svglite(file.path(mofa_v50_figure_dir, paste0("MOFA_v50_three_factor_response_", file_stub, "_feature_arrows_publication.svg")), width = 7.5, height = 6.4, bg = "white")
  graphics::par(mar = c(2.2, 2.3, 3.3, 5.3), xpd = NA)
  mofa_v50_draw_feature_arrow_plot(target_factors, vec_tbl, point_cex = 1.95)
  grDevices::dev.off()
  invisible(list(diagnostic = diag_row, features = vec_tbl))
}

mofa_v50_3d_bundle_f147 <- mofa_v50_save_3d_bundle(c("Factor1", "Factor4", "Factor7"), "Factor1_Factor4_Factor7")
mofa_v50_3d_bundle_f247 <- mofa_v50_save_3d_bundle(c("Factor2", "Factor4", "Factor7"), "Factor2_Factor4_Factor7")

mofa_v50_f147_row <- mofa_v50_find_triple_row(c("Factor1", "Factor4", "Factor7"))
mofa_v50_f247_row <- mofa_v50_find_triple_row(c("Factor2", "Factor4", "Factor7"))
svglite::svglite(file.path(mofa_v50_figure_dir, "MOFA_v50_three_factor_response_main_pair.svg"), width = 13.2, height = 6.3, bg = "white")
graphics::par(mfrow = c(1, 2), mar = c(2.1, 2.2, 3.2, 5.1), xpd = NA)
mofa_v50_draw_three_factor_publication(c("Factor1", "Factor4", "Factor7"), mofa_v50_f147_row, point_cex = 2.00)
mofa_v50_draw_three_factor_publication(c("Factor2", "Factor4", "Factor7"), mofa_v50_f247_row, point_cex = 2.00)
grDevices::dev.off()


#-----------------------------------------------------------------#
# 22.9 v50 curated publication outputs: continuous atlas scales, tangent-edge networks, and organized subdirectories
#-----------------------------------------------------------------#

#------------------------------#
# 22.9.1 Curated subdirectories
#------------------------------#

mofa_v50_factor_figure_dir <- file.path(mofa_v50_figure_dir, "factor")
mofa_v50_three_d_figure_dir <- file.path(mofa_v50_figure_dir, "3d_scatter")
mofa_v50_two_d_figure_dir <- file.path(mofa_v50_figure_dir, "2d_scatter")
mofa_v50_network_figure_dir <- file.path(mofa_v50_figure_dir, "network")
mofa_v50_network_publication_dir <- file.path(mofa_v50_network_figure_dir, "feature_network_publication")
mofa_v50_network_test_dir <- file.path(mofa_v50_network_figure_dir, "feature_network_test")

for (dir_path in c(
  mofa_v50_factor_figure_dir,
  mofa_v50_three_d_figure_dir,
  mofa_v50_two_d_figure_dir,
  mofa_v50_network_figure_dir,
  mofa_v50_network_publication_dir,
  mofa_v50_network_test_dir
)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

# Clean root-level duplicates for curated outputs so the final deliverables are easier to locate.
unlink(Sys.glob(file.path(mofa_v50_figure_dir, "MOFA_v50_network_*.svg")))
unlink(Sys.glob(file.path(mofa_v50_figure_dir, "MOFA_v50_three_factor_response_Factor1_Factor4_Factor7*.svg")))
unlink(Sys.glob(file.path(mofa_v50_figure_dir, "MOFA_v50_three_factor_response_Factor2_Factor4_Factor7*.svg")))
unlink(Sys.glob(file.path(mofa_v50_figure_dir, "MOFA_v50_three_factor_response_main_pair.svg")))
unlink(Sys.glob(file.path(mofa_v50_figure_dir, "MOFA_v50_multifactor_feature_loading_network.svg")))
unlink(Sys.glob(file.path(mofa_v50_figure_dir, "MOFA_v50_selected_triple_multifactor_feature_loading_network.svg")))
unlink(Sys.glob(file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_variance_association_atlas.svg")))
unlink(Sys.glob(file.path(mofa_v50_figure_dir, "MOFA_v50_factor_response_variance_atlas.svg")))

#------------------------------#
# 22.9.2 Factor atlas with continuous color scales
#------------------------------#

mofa_v50_response_size_bucket <- function(x) {
  x <- ifelse(is.finite(x), x, 0)
  factor(as.character(pmin(2, pmax(0, round(x)))), levels = c("0", "1", "2"))
}

mofa_v50_response_enrichment_data <- mofa_v50_response_enrichment_data %>%
  dplyr::mutate(minus_log10_p_bucket = mofa_v50_response_size_bucket(minus_log10_p))

mofa_v50_variance_upper <- max(mofa_v50_variance_plot_data$r2, na.rm = TRUE)
mofa_v50_variance_breaks <- pretty(c(0, mofa_v50_variance_upper), n = 4)
mofa_v50_variance_breaks <- unique(mofa_v50_variance_breaks[mofa_v50_variance_breaks >= 0 & mofa_v50_variance_breaks <= mofa_v50_variance_upper + 1e-8])
if (length(mofa_v50_variance_breaks) < 3) mofa_v50_variance_breaks <- c(0, mofa_v50_variance_upper / 2, mofa_v50_variance_upper)

mofa_v50_assoc_breaks <- pretty(c(0, mofa_v50_association_p_upper), n = 4)
mofa_v50_assoc_breaks <- unique(mofa_v50_assoc_breaks[mofa_v50_assoc_breaks >= 0 & mofa_v50_assoc_breaks <= mofa_v50_association_p_upper + 1e-8])
if (length(mofa_v50_assoc_breaks) < 3) mofa_v50_assoc_breaks <- c(0, mofa_v50_association_p_upper / 2, mofa_v50_association_p_upper)

mofa_v50_enrich_xlim <- range(mofa_v50_response_enrichment_data$standardized_effect, finite = TRUE)
mofa_v50_enrich_xlim <- c(mofa_v50_enrich_xlim[1] - 0.18, mofa_v50_enrich_xlim[2] + 0.18)

p_mofa_v50_variance_heatmap_atlas <- ggplot2::ggplot(
  mofa_v50_variance_plot_data,
  ggplot2::aes(x = view_label, y = factor, fill = r2)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.32) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(r2 >= 0.005, sprintf("%.1f", 100 * r2), "")),
    size = 2.6
  ) +
  ggplot2::scale_fill_gradientn(
    colors = c("white", "#DCE6F2", "#9FBFE3", "#2166AC"),
    limits = c(0, mofa_v50_variance_upper),
    breaks = mofa_v50_variance_breaks,
    name = "Variance\nexplained",
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128, 
      barheight = grid::unit(0.95, "in"),
      barwidth = grid::unit(0.20, "in")
    )
  ) +
  ggplot2::labs(
    title = "Variance explained by factor and omics view",
    subtitle = "Cell labels are percentages.",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.0, color = "grey35"),
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, size = 7.4),
    axis.text.y = ggplot2::element_text(face = "plain", size = 7.7),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 6.8),
    legend.text = ggplot2::element_text(size = 6.5)
  )

p_mofa_v50_response_enrichment_atlas <- ggplot2::ggplot(
  mofa_v50_response_enrichment_data,
  ggplot2::aes(
    x = standardized_effect,
    y = factor,
    color = response_direction,
    fill = response_direction,
    size = minus_log10_p_bucket
  )
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.36) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = standardized_effect, yend = factor),
    linewidth = 0.56,
    alpha = 0.76,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(shape = 21, color = "grey20", stroke = 0.46) +
  ggplot2::geom_text(
    ggplot2::aes(label = significance_label),
    nudge_y = 0.18,
    size = 2.4,
    color = "black",
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~ Timepoint_display, nrow = 1) +
  ggplot2::scale_x_continuous(limits = mofa_v50_enrich_xlim, expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
  ggplot2::scale_color_manual(
    values = mofa_response_colors[c("pCR", "non_pCR")],
    labels = c(pCR = "pCR-enriched", non_pCR = "non-pCR-enriched"),
    name = "Direction"
  ) +
  ggplot2::scale_fill_manual(values = mofa_response_colors[c("pCR", "non_pCR")], guide = "none") +
  ggplot2::scale_size_manual(
    values = c("0" = 1.25, "1" = 2.25, "2" = 3.15),
    name = expression(-log[10](italic(p))),
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "pCR versus non-pCR factor enrichment",
    subtitle = "Hedges' g is oriented as pCR minus non-pCR.",
    x = "Response enrichment (Hedges' g)",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 8.8) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 6.9, color = "grey35"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "plain", size = 7.0),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 6.8),
    legend.text = ggplot2::element_text(size = 6.4),
    plot.margin = ggplot2::margin(5, 5, 5, 1)
  )

p_mofa_v50_association_delta_heatmap <- ggplot2::ggplot(
  mofa_v50_association_delta_data,
  ggplot2::aes(x = comparison_label, y = factor, fill = minus_log10_p)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.30) +
  ggplot2::geom_text(ggplot2::aes(label = effect_label), size = 2.6) +
  ggplot2::scale_fill_gradientn(
    colors = c("white", "#B3B3B3", "#5E5E5E", "black"),
    limits = c(0, mofa_v50_association_p_upper),
    breaks = mofa_v50_assoc_breaks,
    oob = scales::squish,
    name = expression(-log[10](italic(p))),
    guide = ggplot2::guide_colourbar(raster = FALSE, nbin = 128, 
      barheight = grid::unit(0.95, "in"),
      barwidth = grid::unit(0.20, "in")
    )
  ) +
  ggplot2::labs(
    title = "Paired-change associations",
    subtitle = "Text is the signed median contrast.",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 9.0) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 7.0, color = "grey35"),
    axis.text.x = ggplot2::element_text(size = 7.5, face = "plain"),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 6.8),
    legend.text = ggplot2::element_text(size = 6.4)
  )

p_mofa_v50_factor_interpretation_atlas <-
  (p_mofa_v50_variance_heatmap_atlas | p_mofa_v50_response_enrichment_atlas | p_mofa_v50_association_delta_heatmap) +
  patchwork::plot_layout(widths = c(0.54, 0.42, 0.40), guides = "collect") +
  patchwork::plot_annotation(
    title = "MOFA factor interpretation atlas",
    subtitle = paste0(
      "Variance explained provides the omics context, response enrichment is the primary clinical summary, ",
      "and the right heatmap highlights paired-change associations."
    )
  ) &
  ggplot2::theme(legend.position = "right")

ggplot2::ggsave(
  filename = file.path(mofa_v50_factor_figure_dir, "MOFA_v50_factor_response_variance_association_atlas.svg"),
  plot = p_mofa_v50_factor_interpretation_atlas,
  width = 7.5,
  height = max(5.1, 2.6 + 0.30 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

p_mofa_v50_factor_response_variance_atlas <-
  (p_mofa_v50_variance_heatmap_atlas | p_mofa_v50_response_enrichment_atlas) +
  patchwork::plot_layout(widths = c(0.58, 0.42), guides = "collect") +
  patchwork::plot_annotation(
    title = "MOFA factor variance and response enrichment",
    subtitle = "Variance explained heatmap and the response-enrichment summary."
  ) &
  ggplot2::theme(legend.position = "right")

ggplot2::ggsave(
  filename = file.path(mofa_v50_factor_figure_dir, "MOFA_v50_factor_response_variance_atlas.svg"),
  plot = p_mofa_v50_factor_response_variance_atlas,
  width = 5.9,
  height = max(4.9, 2.5 + 0.30 * length(mofa_v50_factor_order)),
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

#------------------------------#
# 22.9.3 Publication network plots only
#------------------------------#

mofa_v50_factor_fill_summary <- mofa_v50_response_enrichment_data %>%
  dplyr::group_by(factor) %>%
  dplyr::arrange(
    dplyr::desc(is.finite(wilcoxon_p)),
    wilcoxon_p,
    dplyr::desc(abs(standardized_effect)),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    enrichment_label = dplyr::case_when(
      !is.finite(standardized_effect) ~ "neutral",
      standardized_effect >= 0 ~ "pCR-enriched",
      TRUE ~ "non-pCR-enriched"
    ),
    factor_fill = dplyr::case_when(
      enrichment_label == "pCR-enriched" ~ unname(mofa_response_colors[["pCR"]]),
      enrichment_label == "non-pCR-enriched" ~ unname(mofa_response_colors[["non_pCR"]]),
      TRUE ~ "#EFE7D3"
    )
  ) %>%
  dplyr::select(factor, enrichment_label, factor_fill)

mofa_v50_select_network_edges <- function(factor_set) {
  active_factor_views <- mofa_variance_explained %>%
    dplyr::mutate(factor = as.character(factor), view = as.character(view)) %>%
    dplyr::filter(factor %in% factor_set, r2 >= mofa_active_view_r2) %>%
    dplyr::select(factor, view, view_r2 = r2)

  seed_edges <- dplyr::bind_rows(
    mofa_feature_weights %>%
      dplyr::filter(factor %in% factor_set, display_eligible) %>%
      dplyr::inner_join(active_factor_views, by = c("factor", "view")) %>%
      dplyr::group_by(factor, view, direction) %>%
      dplyr::slice_max(order_by = abs(weight_within_view), n = mofa_network_features_per_view_direction, with_ties = FALSE) %>%
      dplyr::ungroup(),
    mofa_feature_weights %>%
      dplyr::filter(factor %in% factor_set, display_eligible, is.finite(p_value), p_value < 0.05) %>%
      dplyr::inner_join(active_factor_views, by = c("factor", "view")) %>%
      dplyr::group_by(factor) %>%
      dplyr::arrange(p_value, dplyr::desc(abs(weight_within_view)), .by_group = TRUE) %>%
      dplyr::slice_head(n = 2) %>%
      dplyr::ungroup()
  ) %>%
    dplyr::distinct(factor, view, feature, .keep_all = TRUE) %>%
    dplyr::mutate(feature_node_id = paste(view, feature, sep = "::"), edge_id = paste(factor, feature_node_id, sep = "__"))

  mofa_feature_weights %>%
    dplyr::filter(factor %in% factor_set, display_eligible) %>%
    dplyr::inner_join(active_factor_views, by = c("factor", "view")) %>%
    dplyr::mutate(feature_node_id = paste(view, feature, sep = "::"), edge_id = paste(factor, feature_node_id, sep = "__")) %>%
    dplyr::filter(
      feature_node_id %in% unique(seed_edges$feature_node_id),
      abs(weight_within_view) >= mofa_network_shared_loading_threshold | edge_id %in% seed_edges$edge_id
    ) %>%
    dplyr::mutate(loading_sign = ifelse(weight >= 0, "Positive", "Negative"), loading_strength = abs(weight_within_view))
}

mofa_v50_attach_line_to_circles <- function(x0, y0, x1, y1, r0, r1) {
  dx <- x1 - x0
  dy <- y1 - y0
  d <- sqrt(dx^2 + dy^2)
  d[!is.finite(d) | d == 0] <- 1
  data.frame(
    x_start = x0 + r0 * dx / d,
    y_start = y0 + r0 * dy / d,
    x_end = x1 - r1 * dx / d,
    y_end = y1 - r1 * dy / d
  )
}

mofa_v50_prepare_network_layout <- function(factor_set, orientation = c("vertical", "horizontal")) {
  orientation <- match.arg(orientation)
  factor_set <- unique(factor_set)
  factor_set <- factor_set[factor_set %in% mofa_v50_factor_order]
  network_edges <- mofa_v50_select_network_edges(factor_set)
  if (nrow(network_edges) == 0) return(NULL)

  module_membership <- dplyr::bind_rows(lapply(split(network_edges, network_edges$view), function(view_edge_data) {
    profile <- view_edge_data %>%
      dplyr::select(feature_node_id, factor, weight_within_view) %>%
      tidyr::pivot_wider(names_from = factor, values_from = weight_within_view, values_fill = 0)
    for (factor_name in setdiff(factor_set, colnames(profile))) profile[[factor_name]] <- 0
    profile_matrix <- as.matrix(profile[, factor_set, drop = FALSE])
    rownames(profile_matrix) <- profile$feature_node_id
    module <- if (nrow(profile_matrix) <= 2) rep(1L, nrow(profile_matrix)) else stats::cutree(stats::hclust(stats::dist(profile_matrix), method = "ward.D2"), k = min(2L, max(1L, ceiling(nrow(profile_matrix) / 7))))
    data.frame(feature_node_id = rownames(profile_matrix), module = as.integer(module), stringsAsFactors = FALSE)
  }))

  feature_nodes <- network_edges %>%
    dplyr::group_by(feature_node_id, view, view_label, feature, feature_label) %>%
    dplyr::summarise(
      n_connected_factors = dplyr::n_distinct(factor),
      minimum_response_p = if (any(is.finite(p_value))) min(p_value, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::left_join(module_membership, by = "feature_node_id") %>%
    dplyr::mutate(
      view_order = match(view, required_views),
      shared_feature = n_connected_factors >= 2,
      response_feature = is.finite(minimum_response_p) & minimum_response_p < 0.05,
      feature_label_plotmath = mofa_v50_plotmath_feature_label(view, feature_label)
    ) %>%
    dplyr::arrange(view_order, module, dplyr::desc(n_connected_factors), feature_label)

  modules <- feature_nodes %>%
    dplyr::distinct(view, view_label, view_order, module) %>%
    dplyr::arrange(view_order, module) %>%
    dplyr::group_by(view) %>%
    dplyr::mutate(module_short = paste0("M", dplyr::row_number())) %>%
    dplyr::ungroup()
  feature_nodes <- feature_nodes %>% dplyr::left_join(modules, by = c("view", "view_label", "view_order", "module"))

  module_size_df <- feature_nodes %>%
    dplyr::group_by(view_order, view_label, module, module_short) %>%
    dplyr::summarise(n_features = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(view_order, module)

  if (orientation == "vertical") {
    gap_between_modules <- 0.62
    current_top <- 0
    feature_nodes$x_node <- NA_real_
    feature_nodes$x_label <- NA_real_
    feature_nodes$y <- NA_real_
    for (i in seq_len(nrow(module_size_df))) {
      idx <- feature_nodes$view_order == module_size_df$view_order[i] & feature_nodes$module == module_size_df$module[i]
      n_i <- module_size_df$n_features[i]
      feature_nodes$y[idx] <- -(current_top + seq_len(n_i) - 1)
      feature_nodes$x_node[idx] <- 2.10
      feature_nodes$x_label[idx] <- 2.60
      current_top <- current_top + n_i + gap_between_modules
    }
    module_layout <- feature_nodes %>%
      dplyr::group_by(view_order, view_label, module, module_short) %>%
      dplyr::summarise(
        xmin = 1.88,
        xmax = 2.34,
        ymin = min(y) - 0.38,
        ymax = max(y) + 0.38,
        module_x = 1.96,
        module_y = max(y) + 0.54,
        .groups = "drop"
      )
    factor_nodes <- data.frame(factor = factor_set, stringsAsFactors = FALSE) %>%
      dplyr::left_join(mofa_v50_factor_fill_summary, by = "factor") %>%
      dplyr::mutate(
        enrichment_label = ifelse(is.na(enrichment_label), "neutral", enrichment_label),
        factor_fill = ifelse(is.na(factor_fill), "#EFE7D3", factor_fill),
        x = -3.55,
        y = seq(from = max(feature_nodes$y) - 0.1, to = min(feature_nodes$y) + 0.1, length.out = length(factor_set))
      )
    plot_edges <- network_edges %>%
      dplyr::left_join(factor_nodes %>% dplyr::select(factor, x_factor = x, y_factor = y), by = "factor") %>%
      dplyr::left_join(feature_nodes %>% dplyr::select(feature_node_id, x_feature = x_node, y_feature = y), by = "feature_node_id") %>%
      dplyr::group_by(feature_node_id) %>%
      dplyr::arrange(y_factor, .by_group = TRUE) %>%
      dplyr::mutate(feature_y_target = y_feature + (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.07) %>%
      dplyr::ungroup()
    attach <- mofa_v50_attach_line_to_circles(
      x0 = plot_edges$x_factor,
      y0 = plot_edges$y_factor,
      x1 = plot_edges$x_feature,
      y1 = plot_edges$feature_y_target,
      r0 = 0.31,
      r1 = 0.10
    )
    plot_edges <- dplyr::bind_cols(plot_edges, attach)
    list(
      factor_nodes = factor_nodes,
      feature_nodes = feature_nodes,
      module_layout = module_layout,
      plot_edges = plot_edges,
      xlim = c(-4.5, 5.8),
      ylim = c(min(feature_nodes$y) - 0.8, max(feature_nodes$y) + 0.95)
    )
  } else {
    gap_between_modules <- 1.05
    current_left <- 0
    feature_nodes$x_node <- NA_real_
    feature_nodes$x_label <- NA_real_
    feature_nodes$y <- NA_real_
    for (i in seq_len(nrow(module_size_df))) {
      idx <- feature_nodes$view_order == module_size_df$view_order[i] & feature_nodes$module == module_size_df$module[i]
      n_i <- module_size_df$n_features[i]
      feature_nodes$x_node[idx] <- current_left
      feature_nodes$x_label[idx] <- current_left + 0.42
      feature_nodes$y[idx] <- rev(seq_len(n_i))
      current_left <- current_left + 1.40 + gap_between_modules
    }
    module_layout <- feature_nodes %>%
      dplyr::group_by(view_order, view_label, module, module_short) %>%
      dplyr::summarise(
        xmin = unique(x_node) - 0.24,
        xmax = unique(x_node) + 0.24,
        ymin = min(y) - 0.38,
        ymax = max(y) + 0.38,
        module_x = unique(x_node),
        module_y = max(y) + 0.56,
        .groups = "drop"
      )
    factor_nodes <- data.frame(factor = factor_set, stringsAsFactors = FALSE) %>%
      dplyr::left_join(mofa_v50_factor_fill_summary, by = "factor") %>%
      dplyr::mutate(
        enrichment_label = ifelse(is.na(enrichment_label), "neutral", enrichment_label),
        factor_fill = ifelse(is.na(factor_fill), "#EFE7D3", factor_fill),
        x = seq(from = min(feature_nodes$x_node) - 0.4, to = max(feature_nodes$x_node) + 0.4, length.out = length(factor_set)),
        y = max(module_layout$module_y) + 1.32
      )
    plot_edges <- network_edges %>%
      dplyr::left_join(factor_nodes %>% dplyr::select(factor, x_factor = x, y_factor = y), by = "factor") %>%
      dplyr::left_join(feature_nodes %>% dplyr::select(feature_node_id, x_feature = x_node, y_feature = y), by = "feature_node_id") %>%
      dplyr::group_by(feature_node_id) %>%
      dplyr::arrange(x_factor, .by_group = TRUE) %>%
      dplyr::mutate(feature_x_target = x_feature + (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.07) %>%
      dplyr::ungroup()
    attach <- mofa_v50_attach_line_to_circles(
      x0 = plot_edges$x_factor,
      y0 = plot_edges$y_factor,
      x1 = plot_edges$feature_x_target,
      y1 = plot_edges$y_feature,
      r0 = 0.32,
      r1 = 0.10
    )
    plot_edges <- dplyr::bind_cols(plot_edges, attach)
    list(
      factor_nodes = factor_nodes,
      feature_nodes = feature_nodes,
      module_layout = module_layout,
      plot_edges = plot_edges,
      xlim = c(min(module_layout$xmin) - 0.8, max(module_layout$xmax) + 4.2),
      ylim = c(min(feature_nodes$y) - 0.8, max(factor_nodes$y) + 0.9)
    )
  }
}

mofa_v50_build_network_plot <- function(factor_set, plot_title, plot_subtitle = NULL, orientation = c("vertical", "horizontal"), module_label_angle = 0) {
  orientation <- match.arg(orientation)
  layout_spec <- mofa_v50_prepare_network_layout(factor_set, orientation = orientation)
  if (is.null(layout_spec)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No eligible edges."))
  }
  factor_nodes <- layout_spec$factor_nodes
  feature_nodes <- layout_spec$feature_nodes
  module_layout <- layout_spec$module_layout
  plot_edges <- layout_spec$plot_edges

  ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = module_layout,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = view_label),
      alpha = 0.11,
      color = "grey72",
      linewidth = 0.30
    ) +
    ggplot2::geom_segment(
      data = plot_edges,
      ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = loading_sign, linewidth = loading_strength),
      alpha = 0.28,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      data = factor_nodes,
      ggplot2::aes(x = x, y = y, fill = enrichment_label),
      shape = 21,
      size = 11.2,
      stroke = 0.68,
      color = "grey25"
    ) +
    ggplot2::geom_text(
      data = factor_nodes,
      ggplot2::aes(x = x, y = y, label = factor),
      size = 2.35,
      fontface = "plain"
    ) +
    ggplot2::geom_label(
      data = module_layout,
      ggplot2::aes(x = module_x, y = module_y, label = module_short, fill = view_label),
      color = "grey15",
      label.size = 0.15,
      size = 2.20,
      label.padding = grid::unit(0.07, "lines"),
      fontface = "plain",
      angle = module_label_angle,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = feature_nodes,
      ggplot2::aes(x = x_node, y = y, fill = view_label),
      shape = 21,
      size = 3.0,
      color = "grey20",
      stroke = 0.42
    ) +
    ggplot2::geom_point(
      data = feature_nodes %>% dplyr::filter(shared_feature),
      ggplot2::aes(x = x_node, y = y),
      shape = 21,
      size = 3.45,
      fill = NA,
      color = "#7A3E9D",
      stroke = 0.72
    ) +
    ggplot2::geom_text(
      data = feature_nodes %>% dplyr::filter(response_feature),
      ggplot2::aes(x = x_node, y = y, label = "*"),
      size = 2.1,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = feature_nodes,
      ggplot2::aes(x = x_label, y = y, label = feature_label_plotmath),
      parse = TRUE,
      hjust = 0,
      size = 1.90
    ) +
    ggplot2::scale_color_manual(values = c(Positive = "#D55E00", Negative = "#0072B2"), name = "Feature-weight sign") +
    ggplot2::scale_fill_manual(
      values = c(view_colors, "pCR-enriched" = mofa_response_colors[["pCR"]], "non-pCR-enriched" = mofa_response_colors[["non_pCR"]], "neutral" = "#EFE7D3"),
      breaks = c("pCR-enriched", "non-pCR-enriched", names(view_colors)),
      labels = c("pCR-enriched factor", "non-pCR-enriched factor", names(view_colors)),
      name = NULL,
      guide = ggplot2::guide_legend(override.aes = list(shape = 21, size = 4.0, color = "grey25", alpha = 1))
    ) +
    ggplot2::scale_linewidth_continuous(range = c(0.24, 0.74), guide = "none") +
    ggplot2::coord_cartesian(xlim = layout_spec$xlim, ylim = layout_spec$ylim, clip = "off") +
    ggplot2::labs(
      title = plot_title,
      subtitle = plot_subtitle,
      caption = "Factor circles are colored by response-enrichment direction. Module boxes tightly enclose feature nodes, and purple outlines denote features shared by multiple factors.",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_void(base_size = 8.9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 7.3, color = "grey35"),
      plot.caption = ggplot2::element_text(size = 6.7, color = "grey35", hjust = 0),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 6.8),
      legend.text = ggplot2::element_text(size = 6.6),
      plot.margin = if (orientation == "vertical") ggplot2::margin(4, 145, 4, 4) else ggplot2::margin(4, 80, 4, 4)
    )
}

mofa_v50_save_publication_network_set <- function(factor_set, file_stub) {
  vertical_plot <- mofa_v50_build_network_plot(
    factor_set = factor_set,
    plot_title = paste0("Feature-weight network: ", paste(factor_set, collapse = " + ")),
    plot_subtitle = "Vertical publication layout.",
    orientation = "vertical",
    module_label_angle = 0
  )
  ggplot2::ggsave(
    file.path(mofa_v50_network_publication_dir, paste0("MOFA_v50_network_", file_stub, "_vertical.svg")),
    vertical_plot,
    width = 10.2,
    height = 9.0,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )

  horizontal_plot <- mofa_v50_build_network_plot(
    factor_set = factor_set,
    plot_title = paste0("Feature-weight network: ", paste(factor_set, collapse = " + ")),
    plot_subtitle = "Horizontal publication layout.",
    orientation = "horizontal",
    module_label_angle = 0
  )
  ggplot2::ggsave(
    file.path(mofa_v50_network_publication_dir, paste0("MOFA_v50_network_", file_stub, "_horizontal.svg")),
    horizontal_plot,
    width = 14.2,
    height = 6.6,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )

  horizontal_rot_plot <- mofa_v50_build_network_plot(
    factor_set = factor_set,
    plot_title = paste0("Feature-weight network: ", paste(factor_set, collapse = " + ")),
    plot_subtitle = "Horizontal publication layout with rotated module labels.",
    orientation = "horizontal",
    module_label_angle = 90
  )
  ggplot2::ggsave(
    file.path(mofa_v50_network_publication_dir, paste0("MOFA_v50_network_", file_stub, "_horizontal_module90.svg")),
    horizontal_rot_plot,
    width = 14.2,
    height = 6.6,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

mofa_v50_save_publication_network_set(c("Factor1", "Factor4", "Factor7"), "Factor1_Factor4_Factor7_publication")
mofa_v50_save_publication_network_set(c("Factor2", "Factor4", "Factor7"), "Factor2_Factor4_Factor7_publication")

#------------------------------#
# 22.9.4 Publication 3D plots with centroid cube guides in x, y, and z
#------------------------------#

mofa_v50_find_triple_row <- function(target_factors) {
  candidates <- mofa_v50_triple_permanova_atlas %>%
    dplyr::rowwise() %>%
    dplyr::mutate(match_target = setequal(c(factor_x, factor_y, factor_z), target_factors)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(match_target)
  if (nrow(candidates) == 0) return(mofa_v50_selected_triple_diagnostics[1, , drop = FALSE])
  candidates[1, , drop = FALSE]
}

mofa_v50_draw_filled_cube <- function(s3d_obj, x, y, z, side = 0.11, fill_col = "grey70") {
  v <- data.frame(
    x = x + c(-side, side, side, -side, -side, side, side, -side),
    y = y + c(-side, -side, side, side, -side, -side, side, side),
    z = z + c(-side, -side, -side, -side, side, side, side, side)
  )
  p <- s3d_obj$xyz.convert(v$x, v$y, v$z)
  coords <- cbind(p$x, p$y)
  face_fill <- grDevices::adjustcolor(fill_col, alpha.f = 0.72)
  face_fill_light <- grDevices::adjustcolor(fill_col, alpha.f = 0.54)
  face_fill_dark <- grDevices::adjustcolor(fill_col, alpha.f = 0.86)
  graphics::polygon(coords[c(1,2,6,5),1], coords[c(1,2,6,5),2], col = face_fill, border = "black", lwd = 0.7)
  graphics::polygon(coords[c(2,3,7,6),1], coords[c(2,3,7,6),2], col = face_fill_dark, border = "black", lwd = 0.7)
  graphics::polygon(coords[c(5,6,7,8),1], coords[c(5,6,7,8),2], col = face_fill_light, border = "black", lwd = 0.7)
}

mofa_v50_draw_three_factor_publication <- function(triple_factors, diagnostic_row, point_cex = 2.00) {
  if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
    plot.new(); text(0.5, 0.5, "Package 'scatterplot3d' is required.")
    return(invisible(NULL))
  }
  plot_data <- mofa_v50_subject_mean_factor_data %>%
    dplyr::select(SubjectID, TRG_plot, dplyr::all_of(triple_factors)) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(triple_factors), is.finite), !is.na(TRG_plot)) %>%
    dplyr::mutate(TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR")))
  score_matrix <- scale(as.matrix(plot_data[, triple_factors, drop = FALSE]))
  score_matrix[!is.finite(score_matrix)] <- 0
  plot_data$x <- score_matrix[, 1]
  plot_data$y <- score_matrix[, 2]
  plot_data$z <- score_matrix[, 3]
  plot_data$color <- unname(mofa_response_colors[as.character(plot_data$TRG_plot)])

  s3d <- scatterplot3d::scatterplot3d(
    x = plot_data$x, y = plot_data$y, z = plot_data$z,
    pch = 16, color = plot_data$color, cex.symbols = point_cex,
    type = "p",
    mar = c(2.2, 2.3, 3.3, 5.3),
    main = paste(triple_factors, collapse = " + "),
    xlab = mofa_v50_factor_axis_label(triple_factors[1]),
    ylab = mofa_v50_factor_axis_label(triple_factors[2]),
    zlab = mofa_v50_factor_axis_label(triple_factors[3]),
    angle = 52, scale.y = 1.0, box = TRUE, grid = TRUE
  )
  point_xy <- s3d$xyz.convert(plot_data$x, plot_data$y, plot_data$z)
  graphics::points(point_xy$x, point_xy$y, pch = 21, bg = plot_data$color, col = "grey20", cex = point_cex)

  z_floor <- min(plot_data$z, na.rm = TRUE)
  x_wall <- min(plot_data$x, na.rm = TRUE)
  y_wall <- max(plot_data$y, na.rm = TRUE)
  centroid_df <- plot_data %>%
    dplyr::group_by(TRG_plot) %>%
    dplyr::summarise(x = mean(x), y = mean(y), z = mean(z), .groups = "drop") %>%
    dplyr::mutate(color = unname(mofa_response_colors[as.character(TRG_plot)]), guide_col = grDevices::adjustcolor(color, alpha.f = 0.62))
  for (i in seq_len(nrow(centroid_df))) {
    centroid_xy <- s3d$xyz.convert(centroid_df$x[i], centroid_df$y[i], centroid_df$z[i])
    floor_xy <- s3d$xyz.convert(centroid_df$x[i], centroid_df$y[i], z_floor)
    xwall_xy <- s3d$xyz.convert(x_wall, centroid_df$y[i], centroid_df$z[i])
    ywall_xy <- s3d$xyz.convert(centroid_df$x[i], y_wall, centroid_df$z[i])
    graphics::segments(floor_xy$x, floor_xy$y, centroid_xy$x, centroid_xy$y, col = centroid_df$guide_col[i], lwd = 1.30)
    graphics::segments(xwall_xy$x, xwall_xy$y, centroid_xy$x, centroid_xy$y, col = centroid_df$guide_col[i], lwd = 1.15)
    graphics::segments(ywall_xy$x, ywall_xy$y, centroid_xy$x, centroid_xy$y, col = centroid_df$guide_col[i], lwd = 1.15)
    mofa_v50_draw_filled_cube(s3d, centroid_df$x[i], centroid_df$y[i], centroid_df$z[i], side = 0.10, fill_col = centroid_df$color[i])
  }

  stats_label <- paste0(
    "PERMANOVA p = ", ifelse(is.na(diagnostic_row$permanova_p[1]), "NA", ifelse(diagnostic_row$permanova_p[1] < 0.001, "<0.001", sprintf("%.3f", diagnostic_row$permanova_p[1]))),
    " | R2 = ", ifelse(is.na(diagnostic_row$permanova_r2[1]), "NA", sprintf("%.2f", diagnostic_row$permanova_r2[1])),
    "\nEnergy p = ", ifelse(is.na(diagnostic_row$energy_p[1]), "NA", ifelse(diagnostic_row$energy_p[1] < 0.001, "<0.001", sprintf("%.3f", diagnostic_row$energy_p[1]))),
    " | E = ", ifelse(is.na(diagnostic_row$energy_statistic[1]), "NA", sprintf("%.2f", diagnostic_row$energy_statistic[1]))
  )
  graphics::mtext(stats_label, side = 3, line = 0.38, adj = 0.02, cex = 0.74)
  graphics::legend(
    x = grconvertX(1.05, from = "npc", to = "user"),
    y = grconvertY(0.94, from = "npc", to = "user"),
    legend = c("pCR", "non-pCR", "group centroid cube"),
    pt.bg = c(unname(mofa_response_colors[c("pCR", "non_pCR")]), "grey75"),
    pch = c(21, 21, 22),
    col = c("grey20", "grey20", "black"),
    pt.cex = c(1.5, 1.5, 1.5),
    cex = 0.80,
    bty = "n",
    xpd = NA,
    xjust = 0,
    yjust = 1
  )
}

mofa_v50_select_feature_vectors <- function(target_factors, n_features = 4) {
  raw_tbl <- mofa_feature_weights %>%
    dplyr::filter(factor %in% target_factors, display_eligible, is.finite(weight_within_view))
  summary_tbl <- raw_tbl %>%
    dplyr::group_by(view, feature, feature_label) %>%
    dplyr::summarise(
      n_factors = dplyr::n_distinct(factor),
      best_p = if (any(is.finite(p_value))) min(p_value, na.rm = TRUE) else 1,
      max_abs_loading = max(abs(weight_within_view), na.rm = TRUE),
      combined_score = max_abs_loading + pmax(0, -log10(best_p + 1e-12)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(best_p, dplyr::desc(n_factors), dplyr::desc(combined_score))
  if (sum(summary_tbl$n_factors >= 2) >= min(2, n_features)) summary_tbl <- summary_tbl %>% dplyr::filter(n_factors >= 2)
  selected_tbl <- summary_tbl %>% dplyr::slice_head(n = n_features)
  wide_tbl <- raw_tbl %>%
    dplyr::semi_join(selected_tbl, by = c("view", "feature", "feature_label")) %>%
    dplyr::select(view, feature, feature_label, factor, weight_within_view) %>%
    tidyr::pivot_wider(names_from = factor, values_from = weight_within_view, values_fill = 0)
  for (factor_name in target_factors) if (!(factor_name %in% colnames(wide_tbl))) wide_tbl[[factor_name]] <- 0
  wide_tbl %>%
    dplyr::mutate(
      arrow_length = sqrt(.data[[target_factors[1]]] ^ 2 + .data[[target_factors[2]]] ^ 2 + .data[[target_factors[3]]] ^ 2),
      scale_factor = ifelse(max(arrow_length, na.rm = TRUE) > 0, 2.0 / max(arrow_length, na.rm = TRUE), 1),
      x = .data[[target_factors[1]]] * scale_factor,
      y = .data[[target_factors[2]]] * scale_factor,
      z = .data[[target_factors[3]]] * scale_factor
    )
}

mofa_v50_draw_feature_arrow_plot <- function(target_factors, feature_vector_tbl, point_cex = 1.95) {
  if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
    plot.new(); text(0.5, 0.5, "Package 'scatterplot3d' is required.")
    return(invisible(NULL))
  }
  plot_data <- mofa_v50_subject_mean_factor_data %>%
    dplyr::select(SubjectID, TRG_plot, dplyr::all_of(target_factors)) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(target_factors), is.finite), !is.na(TRG_plot)) %>%
    dplyr::mutate(TRG_plot = factor(as.character(TRG_plot), levels = c("pCR", "non_pCR")))
  score_matrix <- scale(as.matrix(plot_data[, target_factors, drop = FALSE]))
  score_matrix[!is.finite(score_matrix)] <- 0
  plot_data$x <- score_matrix[, 1]; plot_data$y <- score_matrix[, 2]; plot_data$z <- score_matrix[, 3]
  plot_data$color <- unname(mofa_response_colors[as.character(plot_data$TRG_plot)])
  s3d <- scatterplot3d::scatterplot3d(
    x = plot_data$x, y = plot_data$y, z = plot_data$z,
    pch = 16, color = plot_data$color, cex.symbols = point_cex,
    type = "p", mar = c(2.2, 2.3, 3.3, 5.5),
    main = paste0(paste(target_factors, collapse = " + "), " with representative feature arrows"),
    xlab = mofa_v50_factor_axis_label(target_factors[1]),
    ylab = mofa_v50_factor_axis_label(target_factors[2]),
    zlab = mofa_v50_factor_axis_label(target_factors[3]),
    angle = 52, scale.y = 1.0, box = TRUE, grid = TRUE
  )
  point_xy <- s3d$xyz.convert(plot_data$x, plot_data$y, plot_data$z)
  graphics::points(point_xy$x, point_xy$y, pch = 21, bg = plot_data$color, col = "grey20", cex = point_cex)
  if (nrow(feature_vector_tbl) > 0) {
    feature_vector_tbl$arrow_col <- unname(view_colors[feature_vector_tbl$view])
    for (i in seq_len(nrow(feature_vector_tbl))) {
      p0 <- s3d$xyz.convert(0, 0, 0)
      p1 <- s3d$xyz.convert(feature_vector_tbl$x[i], feature_vector_tbl$y[i], feature_vector_tbl$z[i])
      graphics::arrows(p0$x, p0$y, p1$x, p1$y, length = 0.08, lwd = 1.2, col = feature_vector_tbl$arrow_col[i])
      graphics::text(p1$x, p1$y, labels = feature_vector_tbl$feature_label[i], pos = 4, cex = 0.62, col = "black")
    }
  }
  graphics::legend(
    x = grconvertX(1.05, from = "npc", to = "user"),
    y = grconvertY(0.94, from = "npc", to = "user"),
    legend = c("pCR", "non-pCR", unique(feature_vector_tbl$view)),
    pt.bg = c(unname(mofa_response_colors[c("pCR", "non_pCR")]), rep(NA, length(unique(feature_vector_tbl$view)))),
    pch = c(21, 21, rep(NA, length(unique(feature_vector_tbl$view)))),
    col = c("grey20", "grey20", unname(view_colors[unique(feature_vector_tbl$view)])),
    lty = c(NA, NA, rep(1, length(unique(feature_vector_tbl$view)))),
    lwd = c(NA, NA, rep(1.2, length(unique(feature_vector_tbl$view)))),
    pt.cex = 1.3,
    cex = 0.78,
    bty = "n",
    xpd = NA,
    xjust = 0,
    yjust = 1
  )
}

mofa_v50_save_3d_bundle <- function(target_factors, file_stub) {
  diag_row <- mofa_v50_find_triple_row(target_factors)
  svglite::svglite(file.path(mofa_v50_three_d_figure_dir, paste0("MOFA_v50_three_factor_response_", file_stub, "_publication.svg")), width = 6.9, height = 6.2, bg = "white")
  graphics::par(mar = c(2.1, 2.2, 3.2, 5.1), xpd = NA)
  mofa_v50_draw_three_factor_publication(target_factors, diag_row, point_cex = 2.05)
  grDevices::dev.off()

  vec_tbl <- mofa_v50_select_feature_vectors(target_factors, n_features = 4)
  svglite::svglite(file.path(mofa_v50_three_d_figure_dir, paste0("MOFA_v50_three_factor_response_", file_stub, "_feature_arrows_publication.svg")), width = 7.5, height = 6.4, bg = "white")
  graphics::par(mar = c(2.2, 2.3, 3.3, 5.3), xpd = NA)
  mofa_v50_draw_feature_arrow_plot(target_factors, vec_tbl, point_cex = 1.95)
  grDevices::dev.off()
}

mofa_v50_save_3d_bundle(c("Factor1", "Factor4", "Factor7"), "Factor1_Factor4_Factor7")
mofa_v50_save_3d_bundle(c("Factor2", "Factor4", "Factor7"), "Factor2_Factor4_Factor7")

mofa_v50_f147_row <- mofa_v50_find_triple_row(c("Factor1", "Factor4", "Factor7"))
mofa_v50_f247_row <- mofa_v50_find_triple_row(c("Factor2", "Factor4", "Factor7"))
svglite::svglite(file.path(mofa_v50_three_d_figure_dir, "MOFA_v50_three_factor_response_main_pair.svg"), width = 13.2, height = 6.3, bg = "white")
graphics::par(mfrow = c(1, 2), mar = c(2.1, 2.2, 3.2, 5.1), xpd = NA)
mofa_v50_draw_three_factor_publication(c("Factor1", "Factor4", "Factor7"), mofa_v50_f147_row, point_cex = 2.00)
mofa_v50_draw_three_factor_publication(c("Factor2", "Factor4", "Factor7"), mofa_v50_f247_row, point_cex = 2.00)
grDevices::dev.off()

#-----------------------------------------------------------------#
# 22.10 Input/model provenance audit
#-----------------------------------------------------------------#

write.csv(
  data.frame(
    item = c(
      "input_rdata",
      "input_schema_version",
      "primary_model_hdf5",
      "reused_validated_v30_model",
      "publication_revision",
      "legacy_figures_saved"
    ),
    value = c(
      mofa_input_rdata,
      ifelse(is.null(coherence_data$schema_version), NA_character_, as.character(coherence_data$schema_version)),
      mofa_best_model_path,
      as.character(isTRUE(mofa_reuse_existing_best_model) && file.exists(mofa_best_model_path)),
      mofa_v50_source_revision,
      as.character(isTRUE(mofa_save_legacy_figures))
    ),
    stringsAsFactors = FALSE
  ),
  file.path(mofa_v50_result_dir, "MOFA_v50_input_model_provenance.csv"),
  row.names = FALSE
)

#-----------------------------------------------------------------#
# 22.11 Save v50 publication objects without deleting the workspace
#-----------------------------------------------------------------#

mofa_v50_publication_objects <- intersect(
  c(
    "mofa_v50_factor_order",
    "mofa_v50_view_balance_strategy",
    "mofa_v50_global_level_diagnostics",
    "mofa_v50_view_dominance_audit",
    "mofa_v50_association_data",
    "mofa_v50_response_enrichment_data",
    "mofa_v50_paired_factor_data",
    "mofa_v50_pair_atlas",
    "mofa_v50_selected_pair_rows",
    "mofa_v50_selected_factor",
    "mofa_v50_network_selection_audit",
    "mofa_v50_network_edges",
    "mofa_v50_network_feature_nodes",
    "mofa_v50_network_modules",
    "p_mofa_v50_global_level_qc",
    "p_mofa_v50_view_dominance",
    "p_mofa_v50_variance_heatmap",
    "p_mofa_v50_association_pmap",
    "p_mofa_v50_response_enrichment",
    "p_mofa_v50_overall_response",
    "p_mofa_v50_baseline_response",
    "p_mofa_v50_after_response",
    "p_mofa_v50_paired_change",
    "p_mofa_v50_paired_change_by_response",
    "p_mofa_v50_paired_change_effects",
    "p_mofa_v50_selected_factor_r2",
    "p_mofa_v50_selected_factor_baseline",
    "p_mofa_v50_selected_factor_paired",
    "p_mofa_v50_selected_factor_loadings",
    "p_mofa_v50_selected_factor_summary",
    "p_mofa_v50_multifactor_feature_network",
    "mofa_v50_source_revision",
    "mofa_v50_factor_view_display_audit",
    "mofa_v50_network_feature_rank",
    "p_mofa_v50_factor_interpretation_atlas",
    "p_mofa_v50_variance_heatmap_atlas",
    "p_mofa_v50_response_enrichment_atlas",
    "p_mofa_v50_association_delta_heatmap",
    "mofa_v50_association_delta_data",
    "mofa_v50_anchor_factor",
    "mofa_v50_anchor_pair_diagnostics",
    "p_mofa_v50_anchor_pair_scan",
    "mofa_v50_triple_permanova_atlas",
    "mofa_v50_selected_triple_diagnostics",
    "mofa_v50_triple_network_edges",
    "mofa_v50_triple_network_feature_nodes",
    "mofa_v50_triple_network_modules",
    "p_mofa_v50_three_factor_response_maps",
    "p_mofa_v50_selected_triple_feature_network"
  ),
  ls()
)

save(
  list = mofa_v50_publication_objects,
  file = file.path(mofa_v50_result_dir, "MOFA_v50_publication_objects.RData")
)

mofa_v50_essential_objects <- intersect(
  c(
    "mofa_sample_metadata",
    "mofa_availability_wide",
    "mofa_factor_scores",
    "mofa_factor_total_level_correlations",
    "mofa_variance_explained",
    "mofa_variance_total",
    "mofa_factor_summary",
    "mofa_factor_balance",
    "mofa_factor_wilcoxon_tests",
    "mofa_factor_paired_data",
    "mofa_feature_weights",
    "mofa_v50_pair_permanova_atlas",
    mofa_v50_publication_objects
  ),
  ls()
)

save(
  list = mofa_v50_essential_objects,
  file = "results/mofa/MOFA_4omics_analysis_v50_refined.RData"
)

message(
  "v50 completed using the validated v30 latent model. Publication figures: ",
  normalizePath(mofa_v50_figure_dir, winslash = "/", mustWork = FALSE)
)




#-----------------------------------------------------------------#
# 22.11 v50 horizontal network refinement
#-----------------------------------------------------------------#

mofa_v50_prepare_network_layout <- function(factor_set, orientation = c("vertical", "horizontal")) {
  orientation <- match.arg(orientation)
  factor_set <- unique(factor_set)
  factor_set <- factor_set[factor_set %in% mofa_v50_factor_order]
  network_edges <- mofa_v50_select_network_edges(factor_set)
  if (nrow(network_edges) == 0) return(NULL)

  module_membership <- dplyr::bind_rows(lapply(split(network_edges, network_edges$view), function(view_edge_data) {
    profile <- view_edge_data %>%
      dplyr::select(feature_node_id, factor, weight_within_view) %>%
      tidyr::pivot_wider(names_from = factor, values_from = weight_within_view, values_fill = 0)
    for (factor_name in setdiff(factor_set, colnames(profile))) profile[[factor_name]] <- 0
    profile_matrix <- as.matrix(profile[, factor_set, drop = FALSE])
    rownames(profile_matrix) <- profile$feature_node_id
    module <- if (nrow(profile_matrix) <= 2) rep(1L, nrow(profile_matrix)) else stats::cutree(stats::hclust(stats::dist(profile_matrix), method = "ward.D2"), k = min(2L, max(1L, ceiling(nrow(profile_matrix) / 7))))
    data.frame(feature_node_id = rownames(profile_matrix), module = as.integer(module), stringsAsFactors = FALSE)
  }))

  feature_nodes <- network_edges %>%
    dplyr::group_by(feature_node_id, view, view_label, feature, feature_label) %>%
    dplyr::summarise(
      n_connected_factors = dplyr::n_distinct(factor),
      minimum_response_p = if (any(is.finite(p_value))) min(p_value, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::left_join(module_membership, by = "feature_node_id") %>%
    dplyr::mutate(
      view_order = match(view, required_views),
      shared_feature = n_connected_factors >= 2,
      response_feature = is.finite(minimum_response_p) & minimum_response_p < 0.05,
      feature_label_plotmath = mofa_v50_plotmath_feature_label(view, feature_label)
    ) %>%
    dplyr::arrange(view_order, module, dplyr::desc(n_connected_factors), feature_label)

  modules <- feature_nodes %>%
    dplyr::distinct(view, view_label, view_order, module) %>%
    dplyr::arrange(view_order, module) %>%
    dplyr::group_by(view) %>%
    dplyr::mutate(module_short = paste0("M", dplyr::row_number())) %>%
    dplyr::ungroup()
  feature_nodes <- feature_nodes %>%
    dplyr::left_join(modules, by = c("view", "view_label", "view_order", "module"))

  module_size_df <- feature_nodes %>%
    dplyr::group_by(view_order, view_label, module, module_short) %>%
    dplyr::summarise(n_features = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(view_order, module)

  if (orientation == "vertical") {
    gap_between_modules <- 0.62
    current_top <- 0
    feature_nodes$x_node <- NA_real_
    feature_nodes$x_label <- NA_real_
    feature_nodes$y <- NA_real_
    feature_nodes$label_angle <- 0
    feature_nodes$label_hjust <- 0
    feature_nodes$label_vjust <- 0.5
    for (i in seq_len(nrow(module_size_df))) {
      idx <- feature_nodes$view_order == module_size_df$view_order[i] & feature_nodes$module == module_size_df$module[i]
      n_i <- module_size_df$n_features[i]
      feature_nodes$y[idx] <- -(current_top + seq_len(n_i) - 1)
      feature_nodes$x_node[idx] <- 2.05
      feature_nodes$x_label[idx] <- 2.52
      current_top <- current_top + n_i + gap_between_modules
    }
    module_layout <- feature_nodes %>%
      dplyr::group_by(view_order, view_label, module, module_short) %>%
      dplyr::summarise(
        xmin = 1.82,
        xmax = 2.30,
        ymin = min(y) - 0.38,
        ymax = max(y) + 0.38,
        module_x = 1.90,
        module_y = max(y) + 0.54,
        .groups = "drop"
      )
    factor_nodes <- data.frame(factor = factor_set, stringsAsFactors = FALSE) %>%
      dplyr::left_join(mofa_v50_factor_fill_summary, by = "factor") %>%
      dplyr::mutate(
        enrichment_label = ifelse(is.na(enrichment_label), "neutral", enrichment_label),
        factor_fill = ifelse(is.na(factor_fill), "#EFE7D3", factor_fill),
        x = -3.05,
        y = seq(from = max(feature_nodes$y) - 0.1, to = min(feature_nodes$y) + 0.1, length.out = length(factor_set))
      )
    plot_edges <- network_edges %>%
      dplyr::left_join(factor_nodes %>% dplyr::select(factor, x_factor = x, y_factor = y), by = "factor") %>%
      dplyr::left_join(feature_nodes %>% dplyr::select(feature_node_id, x_feature = x_node, y_feature = y), by = "feature_node_id") %>%
      dplyr::group_by(feature_node_id) %>%
      dplyr::arrange(y_factor, .by_group = TRUE) %>%
      dplyr::mutate(feature_y_target = y_feature + (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.07) %>%
      dplyr::ungroup()
    attach <- mofa_v50_attach_line_to_circles(
      x0 = plot_edges$x_factor, y0 = plot_edges$y_factor,
      x1 = plot_edges$x_feature, y1 = plot_edges$feature_y_target,
      r0 = 0.31, r1 = 0.10
    )
    plot_edges <- dplyr::bind_cols(plot_edges, attach)
    list(
      factor_nodes = factor_nodes,
      feature_nodes = feature_nodes,
      module_layout = module_layout,
      plot_edges = plot_edges,
      xlim = c(-4.0, 5.7),
      ylim = c(min(feature_nodes$y) - 0.8, max(feature_nodes$y) + 0.95),
      label_mode = "vertical"
    )
  } else {
    module_gap <- 0.85
    node_step <- 0.48
    current_left <- 0
    base_y <- 0
    feature_nodes$x_node <- NA_real_
    feature_nodes$x_label <- NA_real_
    feature_nodes$y <- NA_real_
    feature_nodes$label_angle <- 90
    feature_nodes$label_hjust <- 1
    feature_nodes$label_vjust <- 0.50
    for (i in seq_len(nrow(module_size_df))) {
      idx <- feature_nodes$view_order == module_size_df$view_order[i] & feature_nodes$module == module_size_df$module[i]
      n_i <- module_size_df$n_features[i]
      x_positions <- current_left + seq(0, by = node_step, length.out = n_i)
      feature_nodes$x_node[idx] <- x_positions
      feature_nodes$x_label[idx] <- x_positions
      feature_nodes$y[idx] <- base_y
      current_left <- max(x_positions) + module_gap
    }
    module_layout <- feature_nodes %>%
      dplyr::group_by(view_order, view_label, module, module_short) %>%
      dplyr::summarise(
        xmin = min(x_node) - 0.18,
        xmax = max(x_node) + 0.18,
        ymin = min(y) - 0.22,
        ymax = max(y) + 0.22,
        module_x = (min(x_node) + max(x_node)) / 2,
        module_y = max(y) + 0.30,
        .groups = "drop"
      )
    factor_nodes <- data.frame(factor = factor_set, stringsAsFactors = FALSE) %>%
      dplyr::left_join(mofa_v50_factor_fill_summary, by = "factor") %>%
      dplyr::mutate(
        enrichment_label = ifelse(is.na(enrichment_label), "neutral", enrichment_label),
        factor_fill = ifelse(is.na(factor_fill), "#EFE7D3", factor_fill),
        x = seq(from = min(module_layout$xmin) + 0.2, to = max(module_layout$xmax) - 0.2, length.out = length(factor_set)),
        y = max(module_layout$ymax) + 2.20
      )
    plot_edges <- network_edges %>%
      dplyr::left_join(factor_nodes %>% dplyr::select(factor, x_factor = x, y_factor = y), by = "factor") %>%
      dplyr::left_join(feature_nodes %>% dplyr::select(feature_node_id, x_feature = x_node, y_feature = y), by = "feature_node_id") %>%
      dplyr::group_by(feature_node_id) %>%
      dplyr::arrange(x_factor, .by_group = TRUE) %>%
      dplyr::mutate(feature_x_target = x_feature + (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.05) %>%
      dplyr::ungroup()
    attach <- mofa_v50_attach_line_to_circles(
      x0 = plot_edges$x_factor, y0 = plot_edges$y_factor,
      x1 = plot_edges$feature_x_target, y1 = plot_edges$y_feature,
      r0 = 0.31, r1 = 0.10
    )
    plot_edges <- dplyr::bind_cols(plot_edges, attach)
    list(
      factor_nodes = factor_nodes,
      feature_nodes = feature_nodes,
      module_layout = module_layout,
      plot_edges = plot_edges,
      xlim = c(min(module_layout$xmin) - 0.45, max(module_layout$xmax) + 0.45),
      ylim = c(-3.55, max(factor_nodes$y) + 0.55),
      label_mode = "horizontal"
    )
  }
}

mofa_v50_build_network_plot <- function(factor_set, plot_title, plot_subtitle = NULL, orientation = c("vertical", "horizontal"), module_label_angle = 0) {
  orientation <- match.arg(orientation)
  layout_spec <- mofa_v50_prepare_network_layout(factor_set, orientation = orientation)
  if (is.null(layout_spec)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No eligible edges."))
  }
  factor_nodes <- layout_spec$factor_nodes
  feature_nodes <- layout_spec$feature_nodes
  module_layout <- layout_spec$module_layout
  plot_edges <- layout_spec$plot_edges

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = module_layout,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = view_label),
      alpha = 0.11,
      color = "grey72",
      linewidth = 0.30
    ) +
    ggplot2::geom_segment(
      data = plot_edges,
      ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = loading_sign, linewidth = loading_strength),
      alpha = 0.30,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      data = factor_nodes,
      ggplot2::aes(x = x, y = y, fill = enrichment_label),
      shape = 21,
      size = 11.0,
      stroke = 0.68,
      color = "grey25"
    ) +
    ggplot2::geom_text(
      data = factor_nodes,
      ggplot2::aes(x = x, y = y, label = factor),
      size = 2.35,
      fontface = "plain",
      angle = 0
    ) +
    ggplot2::geom_label(
      data = module_layout,
      ggplot2::aes(x = module_x, y = module_y, label = module_short, fill = view_label),
      color = "grey15",
      label.size = 0.15,
      size = 2.10,
      label.padding = grid::unit(0.07, "lines"),
      fontface = "plain",
      angle = module_label_angle,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = feature_nodes,
      ggplot2::aes(x = x_node, y = y, fill = view_label),
      shape = 21,
      size = 3.0,
      color = "grey20",
      stroke = 0.42
    ) +
    ggplot2::geom_point(
      data = feature_nodes %>% dplyr::filter(shared_feature),
      ggplot2::aes(x = x_node, y = y),
      shape = 21,
      size = 3.45,
      fill = NA,
      color = "#7A3E9D",
      stroke = 0.72
    ) +
    ggplot2::geom_text(
      data = feature_nodes %>% dplyr::filter(response_feature),
      ggplot2::aes(x = x_node, y = y, label = "*"),
      size = 2.1,
      color = "black"
    ) +
    ggplot2::scale_color_manual(values = c(Positive = "#D55E00", Negative = "#0072B2"), name = "Feature-weight sign") +
    ggplot2::scale_fill_manual(
      values = c(view_colors, "pCR-enriched" = mofa_response_colors[["pCR"]], "non-pCR-enriched" = mofa_response_colors[["non_pCR"]], "neutral" = "#EFE7D3"),
      breaks = c("pCR-enriched", "non-pCR-enriched", names(view_colors)),
      labels = c("pCR-enriched factor", "non-pCR-enriched factor", names(view_colors)),
      name = NULL,
      guide = ggplot2::guide_legend(override.aes = list(shape = 21, size = 4.0, color = "grey25", alpha = 1))
    ) +
    ggplot2::scale_linewidth_continuous(range = c(0.24, 0.74), guide = "none") +
    ggplot2::coord_cartesian(xlim = layout_spec$xlim, ylim = layout_spec$ylim, clip = "off") +
    ggplot2::labs(
      title = plot_title,
      subtitle = plot_subtitle,
      caption = "Factor circles are colored by response-enrichment direction. Module boxes tightly enclose feature nodes, and purple outlines denote features shared by multiple factors.",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_void(base_size = 8.9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 7.3, color = "grey35"),
      plot.caption = ggplot2::element_text(size = 6.7, color = "grey35", hjust = 0),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 6.8),
      legend.text = ggplot2::element_text(size = 6.6),
      plot.margin = if (orientation == "vertical") ggplot2::margin(4, 145, 4, 4) else ggplot2::margin(4, 20, 18, 4)
    )

  if (orientation == "horizontal") {
    p <- p +
      ggplot2::geom_text(
        data = feature_nodes,
        ggplot2::aes(x = x_label, y = y - 0.30, label = feature_label_plotmath),
        parse = TRUE,
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 1.85
      )
  } else {
    p <- p +
      ggplot2::geom_text(
        data = feature_nodes,
        ggplot2::aes(x = x_label, y = y, label = feature_label_plotmath),
        parse = TRUE,
        hjust = 0,
        size = 1.90
      )
  }
  p
}

mofa_v50_save_publication_network_set <- function(factor_set, file_stub) {
  vertical_plot <- mofa_v50_build_network_plot(
    factor_set = factor_set,
    plot_title = paste0("Feature-weight network: ", paste(factor_set, collapse = " + ")),
    plot_subtitle = "Vertical publication layout.",
    orientation = "vertical",
    module_label_angle = 0
  )
  ggplot2::ggsave(
    file.path(mofa_v50_network_publication_dir, paste0("MOFA_v50_network_", file_stub, "_vertical.svg")),
    vertical_plot,
    width = 10.2,
    height = 9.0,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )

  horizontal_plot <- mofa_v50_build_network_plot(
    factor_set = factor_set,
    plot_title = paste0("Feature-weight network: ", paste(factor_set, collapse = " + ")),
    plot_subtitle = "Horizontal publication layout.",
    orientation = "horizontal",
    module_label_angle = 0
  )
  ggplot2::ggsave(
    file.path(mofa_v50_network_publication_dir, paste0("MOFA_v50_network_", file_stub, "_horizontal.svg")),
    horizontal_plot,
    width = 16.8,
    height = 7.2,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )

  horizontal_rot_plot <- mofa_v50_build_network_plot(
    factor_set = factor_set,
    plot_title = paste0("Feature-weight network: ", paste(factor_set, collapse = " + ")),
    plot_subtitle = "Horizontal publication layout with rotated module labels.",
    orientation = "horizontal",
    module_label_angle = 90
  )
  ggplot2::ggsave(
    file.path(mofa_v50_network_publication_dir, paste0("MOFA_v50_network_", file_stub, "_horizontal_module90.svg")),
    horizontal_rot_plot,
    width = 16.8,
    height = 7.2,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

# overwrite the main publication network outputs with the refined horizontal style
mofa_v50_save_publication_network_set(c("Factor1", "Factor4", "Factor7"), "Factor1_Factor4_Factor7_publication")
mofa_v50_save_publication_network_set(c("Factor2", "Factor4", "Factor7"), "Factor2_Factor4_Factor7_publication")


#-----------------------------------------------------------------#
# 22.12 v50 horizontal network width-series refinement
#-----------------------------------------------------------------#

mofa_v50_prepare_network_layout <- function(
  factor_set,
  orientation = c("vertical", "horizontal"),
  horizontal_width_scale = 1.0,
  horizontal_factor_gap = 2.65
) {
  orientation <- match.arg(orientation)
  factor_set <- unique(factor_set)
  factor_set <- factor_set[factor_set %in% mofa_v50_factor_order]
  network_edges <- mofa_v50_select_network_edges(factor_set)
  if (nrow(network_edges) == 0) return(NULL)

  module_membership <- dplyr::bind_rows(lapply(split(network_edges, network_edges$view), function(view_edge_data) {
    profile <- view_edge_data %>%
      dplyr::select(feature_node_id, factor, weight_within_view) %>%
      tidyr::pivot_wider(names_from = factor, values_from = weight_within_view, values_fill = 0)
    for (factor_name in setdiff(factor_set, colnames(profile))) profile[[factor_name]] <- 0
    profile_matrix <- as.matrix(profile[, factor_set, drop = FALSE])
    rownames(profile_matrix) <- profile$feature_node_id
    module <- if (nrow(profile_matrix) <= 2) rep(1L, nrow(profile_matrix)) else stats::cutree(stats::hclust(stats::dist(profile_matrix), method = "ward.D2"), k = min(2L, max(1L, ceiling(nrow(profile_matrix) / 7))))
    data.frame(feature_node_id = rownames(profile_matrix), module = as.integer(module), stringsAsFactors = FALSE)
  }))

  feature_nodes <- network_edges %>%
    dplyr::group_by(feature_node_id, view, view_label, feature, feature_label) %>%
    dplyr::summarise(
      n_connected_factors = dplyr::n_distinct(factor),
      minimum_response_p = if (any(is.finite(p_value))) min(p_value, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::left_join(module_membership, by = "feature_node_id") %>%
    dplyr::mutate(
      view_order = match(view, required_views),
      shared_feature = n_connected_factors >= 2,
      response_feature = is.finite(minimum_response_p) & minimum_response_p < 0.05,
      feature_label_plotmath = mofa_v50_plotmath_feature_label(view, feature_label)
    ) %>%
    dplyr::arrange(view_order, module, dplyr::desc(n_connected_factors), feature_label)

  modules <- feature_nodes %>%
    dplyr::distinct(view, view_label, view_order, module) %>%
    dplyr::arrange(view_order, module) %>%
    dplyr::group_by(view) %>%
    dplyr::mutate(module_short = paste0("M", dplyr::row_number())) %>%
    dplyr::ungroup()
  feature_nodes <- feature_nodes %>%
    dplyr::left_join(modules, by = c("view", "view_label", "view_order", "module"))

  module_size_df <- feature_nodes %>%
    dplyr::group_by(view_order, view_label, module, module_short) %>%
    dplyr::summarise(n_features = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(view_order, module)

  if (orientation == "vertical") {
    gap_between_modules <- 0.62
    current_top <- 0
    feature_nodes$x_node <- NA_real_
    feature_nodes$x_label <- NA_real_
    feature_nodes$y <- NA_real_
    for (i in seq_len(nrow(module_size_df))) {
      idx <- feature_nodes$view_order == module_size_df$view_order[i] & feature_nodes$module == module_size_df$module[i]
      n_i <- module_size_df$n_features[i]
      feature_nodes$y[idx] <- -(current_top + seq_len(n_i) - 1)
      feature_nodes$x_node[idx] <- 2.05
      feature_nodes$x_label[idx] <- 2.52
      current_top <- current_top + n_i + gap_between_modules
    }
    module_layout <- feature_nodes %>%
      dplyr::group_by(view_order, view_label, module, module_short) %>%
      dplyr::summarise(
        xmin = 1.82,
        xmax = 2.30,
        ymin = min(y) - 0.38,
        ymax = max(y) + 0.38,
        module_x = 1.90,
        module_y = max(y) + 0.54,
        .groups = "drop"
      )
    factor_nodes <- data.frame(factor = factor_set, stringsAsFactors = FALSE) %>%
      dplyr::left_join(mofa_v50_factor_fill_summary, by = "factor") %>%
      dplyr::mutate(
        enrichment_label = ifelse(is.na(enrichment_label), "neutral", enrichment_label),
        factor_fill = ifelse(is.na(factor_fill), "#EFE7D3", factor_fill),
        x = -3.05,
        y = seq(from = max(feature_nodes$y) - 0.1, to = min(feature_nodes$y) + 0.1, length.out = length(factor_set))
      )
    plot_edges <- network_edges %>%
      dplyr::left_join(factor_nodes %>% dplyr::select(factor, x_factor = x, y_factor = y), by = "factor") %>%
      dplyr::left_join(feature_nodes %>% dplyr::select(feature_node_id, x_feature = x_node, y_feature = y), by = "feature_node_id") %>%
      dplyr::group_by(feature_node_id) %>%
      dplyr::arrange(y_factor, .by_group = TRUE) %>%
      dplyr::mutate(feature_y_target = y_feature + (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.07) %>%
      dplyr::ungroup()
    attach <- mofa_v50_attach_line_to_circles(
      x0 = plot_edges$x_factor, y0 = plot_edges$y_factor,
      x1 = plot_edges$x_feature, y1 = plot_edges$feature_y_target,
      r0 = 0.31, r1 = 0.10
    )
    plot_edges <- dplyr::bind_cols(plot_edges, attach)
    return(list(
      factor_nodes = factor_nodes,
      feature_nodes = feature_nodes,
      module_layout = module_layout,
      plot_edges = plot_edges,
      xlim = c(-4.0, 5.7),
      ylim = c(min(feature_nodes$y) - 0.8, max(feature_nodes$y) + 0.95)
    ))
  }

  # horizontal layout: factor nodes on the top, compact long-format modules on the bottom
  module_gap <- 0.85
  node_step <- 0.48
  current_left <- 0
  base_y <- 0
  feature_nodes$x_node <- NA_real_
  feature_nodes$x_label <- NA_real_
  feature_nodes$y <- NA_real_
  for (i in seq_len(nrow(module_size_df))) {
    idx <- feature_nodes$view_order == module_size_df$view_order[i] & feature_nodes$module == module_size_df$module[i]
    n_i <- module_size_df$n_features[i]
    x_positions <- current_left + seq(0, by = node_step, length.out = n_i)
    feature_nodes$x_node[idx] <- x_positions
    feature_nodes$x_label[idx] <- x_positions
    feature_nodes$y[idx] <- base_y
    current_left <- max(x_positions) + module_gap
  }

  # horizontally compress around the center for 100/90/80/70% variants
  x_center <- mean(range(feature_nodes$x_node, na.rm = TRUE))
  feature_nodes$x_node <- x_center + (feature_nodes$x_node - x_center) * horizontal_width_scale
  feature_nodes$x_label <- x_center + (feature_nodes$x_label - x_center) * horizontal_width_scale

  module_layout <- feature_nodes %>%
    dplyr::group_by(view_order, view_label, module, module_short) %>%
    dplyr::summarise(
      xmin = min(x_node) - 0.18,
      xmax = max(x_node) + 0.18,
      ymin = min(y) - 0.22,
      ymax = max(y) + 0.22,
      module_x = (min(x_node) + max(x_node)) / 2,
      module_y = max(y) + 0.30,
      .groups = "drop"
    )

  factor_nodes <- data.frame(factor = factor_set, stringsAsFactors = FALSE) %>%
    dplyr::left_join(mofa_v50_factor_fill_summary, by = "factor") %>%
    dplyr::mutate(
      enrichment_label = ifelse(is.na(enrichment_label), "neutral", enrichment_label),
      factor_fill = ifelse(is.na(factor_fill), "#EFE7D3", factor_fill),
      x = seq(from = min(module_layout$xmin) + 0.25, to = max(module_layout$xmax) - 0.25, length.out = length(factor_set)),
      y = max(module_layout$ymax) + horizontal_factor_gap
    )

  plot_edges <- network_edges %>%
    dplyr::left_join(factor_nodes %>% dplyr::select(factor, x_factor = x, y_factor = y), by = "factor") %>%
    dplyr::left_join(feature_nodes %>% dplyr::select(feature_node_id, x_feature = x_node, y_feature = y), by = "feature_node_id") %>%
    dplyr::group_by(feature_node_id) %>%
    dplyr::arrange(x_factor, .by_group = TRUE) %>%
    dplyr::mutate(feature_x_target = x_feature + (dplyr::row_number() - (dplyr::n() + 1) / 2) * 0.05) %>%
    dplyr::ungroup()
  attach <- mofa_v50_attach_line_to_circles(
    x0 = plot_edges$x_factor, y0 = plot_edges$y_factor,
    x1 = plot_edges$feature_x_target, y1 = plot_edges$y_feature,
    r0 = 0.31, r1 = 0.10
  )
  plot_edges <- dplyr::bind_cols(plot_edges, attach)

  list(
    factor_nodes = factor_nodes,
    feature_nodes = feature_nodes,
    module_layout = module_layout,
    plot_edges = plot_edges,
    xlim = c(min(module_layout$xmin) - 0.35, max(module_layout$xmax) + 0.35),
    ylim = c(-3.40, max(factor_nodes$y) + 0.45)
  )
}

mofa_v50_build_network_plot <- function(
  factor_set,
  plot_title,
  plot_subtitle = NULL,
  orientation = c("vertical", "horizontal"),
  module_label_angle = 0,
  horizontal_width_scale = 1.0,
  horizontal_factor_gap = 2.65
) {
  orientation <- match.arg(orientation)
  layout_spec <- mofa_v50_prepare_network_layout(
    factor_set = factor_set,
    orientation = orientation,
    horizontal_width_scale = horizontal_width_scale,
    horizontal_factor_gap = horizontal_factor_gap
  )
  if (is.null(layout_spec)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No eligible edges."))
  }
  factor_nodes <- layout_spec$factor_nodes
  feature_nodes <- layout_spec$feature_nodes
  module_layout <- layout_spec$module_layout
  plot_edges <- layout_spec$plot_edges

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = module_layout,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = view_label),
      alpha = 0.11,
      color = "grey72",
      linewidth = 0.30
    ) +
    ggplot2::geom_segment(
      data = plot_edges,
      ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = loading_sign, linewidth = loading_strength),
      alpha = 0.30,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      data = factor_nodes,
      ggplot2::aes(x = x, y = y, fill = enrichment_label),
      shape = 21,
      size = 11.0,
      stroke = 0.68,
      color = "grey25"
    ) +
    ggplot2::geom_text(
      data = factor_nodes,
      ggplot2::aes(x = x, y = y, label = factor),
      size = 2.35,
      fontface = "plain"
    ) +
    ggplot2::geom_label(
      data = module_layout,
      ggplot2::aes(x = module_x, y = module_y, label = module_short, fill = view_label),
      color = "grey15",
      label.size = 0.15,
      size = 2.10,
      label.padding = grid::unit(0.07, "lines"),
      fontface = "plain",
      angle = module_label_angle,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = feature_nodes,
      ggplot2::aes(x = x_node, y = y, fill = view_label),
      shape = 21,
      size = 3.0,
      color = "grey20",
      stroke = 0.42
    ) +
    ggplot2::geom_point(
      data = feature_nodes %>% dplyr::filter(shared_feature),
      ggplot2::aes(x = x_node, y = y),
      shape = 21,
      size = 3.45,
      fill = NA,
      color = "#7A3E9D",
      stroke = 0.72
    ) +
    ggplot2::geom_text(
      data = feature_nodes %>% dplyr::filter(response_feature),
      ggplot2::aes(x = x_node, y = y, label = "*"),
      size = 2.1,
      color = "black"
    ) +
    ggplot2::scale_color_manual(values = c(Positive = "#D55E00", Negative = "#0072B2"), name = "Feature-weight sign") +
    ggplot2::scale_fill_manual(
      values = c(view_colors, "pCR-enriched" = mofa_response_colors[["pCR"]], "non-pCR-enriched" = mofa_response_colors[["non_pCR"]], "neutral" = "#EFE7D3"),
      breaks = c("pCR-enriched", "non-pCR-enriched", names(view_colors)),
      labels = c("pCR-enriched factor", "non-pCR-enriched factor", names(view_colors)),
      name = NULL,
      guide = ggplot2::guide_legend(override.aes = list(shape = 21, size = 4.0, color = "grey25", alpha = 1))
    ) +
    ggplot2::scale_linewidth_continuous(range = c(0.24, 0.74), guide = "none") +
    ggplot2::coord_cartesian(xlim = layout_spec$xlim, ylim = layout_spec$ylim, clip = "off") +
    ggplot2::labs(
      title = plot_title,
      subtitle = plot_subtitle,
      caption = "Factor circles are colored by response-enrichment direction. Module boxes tightly enclose feature nodes, and purple outlines denote features shared by multiple factors.",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_void(base_size = 8.9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 7.3, color = "grey35"),
      plot.caption = ggplot2::element_text(size = 6.7, color = "grey35", hjust = 0),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 6.8),
      legend.text = ggplot2::element_text(size = 6.6),
      plot.margin = if (orientation == "vertical") ggplot2::margin(4, 145, 4, 4) else ggplot2::margin(4, 18, 16, 4)
    )

  if (orientation == "horizontal") {
    p <- p +
      ggplot2::geom_text(
        data = feature_nodes,
        ggplot2::aes(x = x_label, y = y - 0.24, label = feature_label_plotmath),
        parse = TRUE,
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 1.82
      )
  } else {
    p <- p +
      ggplot2::geom_text(
        data = feature_nodes,
        ggplot2::aes(x = x_label, y = y, label = feature_label_plotmath),
        parse = TRUE,
        hjust = 0,
        size = 1.90
      )
  }
  p
}

mofa_v50_save_publication_network_set <- function(factor_set, file_stub) {
  vertical_plot <- mofa_v50_build_network_plot(
    factor_set = factor_set,
    plot_title = paste0("Feature-weight network: ", paste(factor_set, collapse = " + ")),
    plot_subtitle = "Vertical publication layout.",
    orientation = "vertical"
  )
  ggplot2::ggsave(
    file.path(mofa_v50_network_publication_dir, paste0("MOFA_v50_network_", file_stub, "_vertical.svg")),
    vertical_plot,
    width = 10.2,
    height = 9.0,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )

  horizontal_plot <- mofa_v50_build_network_plot(
    factor_set = factor_set,
    plot_title = paste0("Feature-weight network: ", paste(factor_set, collapse = " + ")),
    plot_subtitle = "Horizontal publication layout.",
    orientation = "horizontal",
    horizontal_width_scale = 1.0,
    horizontal_factor_gap = 2.65
  )
  ggplot2::ggsave(
    file.path(mofa_v50_network_publication_dir, paste0("MOFA_v50_network_", file_stub, "_horizontal.svg")),
    horizontal_plot,
    width = 13.2,
    height = 6.9,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )

  horizontal_rot_plot <- mofa_v50_build_network_plot(
    factor_set = factor_set,
    plot_title = paste0("Feature-weight network: ", paste(factor_set, collapse = " + ")),
    plot_subtitle = "Horizontal publication layout with rotated module labels.",
    orientation = "horizontal",
    module_label_angle = 90,
    horizontal_width_scale = 1.0,
    horizontal_factor_gap = 2.65
  )
  ggplot2::ggsave(
    file.path(mofa_v50_network_publication_dir, paste0("MOFA_v50_network_", file_stub, "_horizontal_module90.svg")),
    horizontal_rot_plot,
    width = 13.2,
    height = 6.9,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

mofa_v50_save_horizontal_width_series <- function(factor_set, file_stub, base_width = 13.2, base_height = 6.9) {
  width_specs <- data.frame(scale = c(1.0, 0.9, 0.8, 0.7), suffix = c("w100", "w090", "w080", "w070"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(width_specs))) {
    plot_obj <- mofa_v50_build_network_plot(
      factor_set = factor_set,
      plot_title = paste0("Feature-weight network: ", paste(factor_set, collapse = " + ")),
      plot_subtitle = paste0("Horizontal publication layout (", width_specs$suffix[i], ")."),
      orientation = "horizontal",
      horizontal_width_scale = width_specs$scale[i],
      horizontal_factor_gap = 2.65
    )
    ggplot2::ggsave(
      file.path(mofa_v50_network_publication_dir, paste0("MOFA_v50_network_", file_stub, "_horizontal_", width_specs$suffix[i], ".svg")),
      plot_obj,
      width = base_width * width_specs$scale[i],
      height = base_height,
      units = "in",
      device = svglite::svglite,
      bg = "white"
    )
  }
}

# overwrite the main publication network outputs with the refined layout
mofa_v50_save_publication_network_set(c("Factor1", "Factor4", "Factor7"), "Factor1_Factor4_Factor7_publication")
mofa_v50_save_publication_network_set(c("Factor2", "Factor4", "Factor7"), "Factor2_Factor4_Factor7_publication")

# width series requested for Factor1 + Factor4 + Factor7 horizontal network
mofa_v50_save_horizontal_width_series(
  factor_set = c("Factor1", "Factor4", "Factor7"),
  file_stub = "Factor1_Factor4_Factor7_publication"
)
