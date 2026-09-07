#-----------------------------------------------------------------#
# One-time MOFA2 Python environment setup for Windows / R 4.3
#
# Purpose:
#   - avoid the binary h5py/HDF5 DLL failure observed in the old
#     MOFA2 1.12.1 basilisk virtual environment
#   - install compiled numerical libraries consistently from conda-forge
#   - install only the pure-Python mofapy2 package with pip --no-deps
#
# Run this script once in a fresh R session. After successful completion,
# run 7-12.1_Coherence_4omics_MOFA_analysis_v3_conda.R.
#-----------------------------------------------------------------#

rm(list = ls())
options(stringsAsFactors = FALSE)
setwd("D:/2-연구/2-CRC metagenomics/")

if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate")
}

mofa_conda_env <- "mofa2_py310_070"
rebuild_environment <- FALSE
install_miniconda_if_missing <- TRUE

Sys.setenv(RETICULATE_AUTOCONFIGURE = "FALSE")

conda_binary <- tryCatch(
  reticulate::conda_binary(),
  error = function(e) NA_character_
)

if (!is.character(conda_binary) || length(conda_binary) != 1 ||
    is.na(conda_binary) || !file.exists(conda_binary)) {
  if (!install_miniconda_if_missing) {
    stop(
      "No usable conda installation was found. Install Miniconda/Miniforge, ",
      "or set install_miniconda_if_missing <- TRUE.",
      call. = FALSE
    )
  }

  reticulate::install_miniconda()
  conda_binary <- reticulate::conda_binary()
}

existing_envs <- reticulate::conda_list(conda = conda_binary)

if (
  rebuild_environment &&
  mofa_conda_env %in% existing_envs$name
) {
  reticulate::conda_remove(
    envname = mofa_conda_env,
    conda = conda_binary
  )

  existing_envs <- reticulate::conda_list(
    conda = conda_binary
  )
}

if (!mofa_conda_env %in% existing_envs$name) {
  # These versions reproduce the numerical stack bundled with the
  # MOFA2 1.12/1.20-era mofapy2 0.7.0 environment, but are installed
  # consistently from conda-forge rather than mixed binary sources.
  reticulate::conda_create(
    envname = mofa_conda_env,
    packages = c(
      "python==3.10",
      "pip",
      "numpy==1.23.1",
      "scipy==1.8.1",
      "pandas==1.4.3",
      "h5py==3.6.0",
      "scikit-learn==1.1.1",
      "dtw-python==1.2.2",
      "anndata==0.8.0",
      "natsort==8.4.0",
      "packaging"
    ),
    forge = TRUE,
    conda = conda_binary
  )

  # mofapy2 itself is pure Python. --no-deps prevents pip from replacing
  # conda-forge's compiled numpy/scipy/h5py binaries with incompatible wheels.
  reticulate::conda_install(
    envname = mofa_conda_env,
    packages = "mofapy2==0.7.0",
    pip = TRUE,
    pip_options = "--no-deps",
    conda = conda_binary
  )
}

mofa_python <- reticulate::conda_python(
  envname = mofa_conda_env,
  conda = conda_binary
)

if (length(mofa_python) != 1 || !file.exists(mofa_python)) {
  stop(
    "The conda environment was created, but its Python executable was not found.",
    call. = FALSE
  )
}

python_check_code <- paste(
  "import sys",
  "import numpy, scipy, pandas, h5py, sklearn, anndata",
  "from importlib.metadata import version as package_version",
  "from mofapy2.run import entry_point",
  "print('python=' + sys.version.replace(chr(10), ' '))",
  "print('numpy=' + numpy.__version__)",
  "print('scipy=' + scipy.__version__)",
  "print('pandas=' + pandas.__version__)",
  "print('h5py=' + h5py.__version__)",
  "print('scikit-learn=' + sklearn.__version__)",
  "print('anndata=' + anndata.__version__)",
  "print('mofapy2=' + package_version('mofapy2'))",
  sep = "; "
)

python_check <- system2(
  command = mofa_python,
  args = c(
    "-c",
    shQuote(python_check_code)
  ),
  stdout = TRUE,
  stderr = TRUE
)

python_check_status <- attr(
  python_check,
  "status"
)

if (
  !is.null(python_check_status) &&
  python_check_status != 0
) {
  stop(
    paste(
      c(
        "The new environment was created, but its import check failed:",
        python_check
      ),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

dir.create(
  "results/mofa",
  recursive = TRUE,
  showWarnings = FALSE
)

writeLines(
  c(
    paste0("environment=", mofa_conda_env),
    paste0("python=", normalizePath(mofa_python, winslash = "/")),
    python_check
  ),
  "results/mofa/MOFA_conda_environment_setup.txt"
)

message(
  paste(
    c(
      "MOFA Python environment is ready.",
      paste0("Environment: ", mofa_conda_env),
      paste0("Python: ", mofa_python),
      "Restart R, then run:",
      "7-12.1_Coherence_4omics_MOFA_analysis_v3_conda.R"
    ),
    collapse = "\n"
  )
)
