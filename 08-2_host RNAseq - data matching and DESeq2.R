
#-----------------------------------------------------------------#
#
#
#
# Host RNA-seq x targeted metabolite integration
#
# Baseline tumor samples, pCR vs non-pCR
#
#
#
#-----------------------------------------------------------------#

rm(list = ls())
options(stringsAsFactors = FALSE)
setwd("D:/2-연구/2-CRC metagenomics/")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

dir.create("host_RNAseq/results_clean_metabolite_host", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

#-----------------------------------------------------------------#

# 1. Load minimal microbiome data and external input files

#-----------------------------------------------------------------#

load("260224 final Input file/RData/m_s_t_only.RData")

if (!exists("m")) {
  stop("Object 'm' was not found after loading m_s_t_only.RData.")
}

if (!exists("s")) {
  warning("Object 's' was not found after loading m_s_t_only.RData.")
}

if (!exists("t")) {
  warning("Object 't' was not found after loading m_s_t_only.RData.")
}

met_raw <- read.csv("input/metabolites.csv", check.names = FALSE)
id_map <- read.csv("input/metadata_matching.csv", check.names = FALSE)

rna <- read.csv("host_RNAseq/TNT_Expression_Profile.GRCh38.gene.csv", check.names = FALSE)
rna_meta <- read.csv("host_RNAseq/meta_RNA.csv", check.names = FALSE)
rna_map <- read.csv("host_RNAseq/chart_numb_matching.csv", check.names = FALSE)

met_vars <- setdiff(colnames(met_raw), "SampleID")

#-----------------------------------------------------------------#

# 2. Build metabolite metadata table

#-----------------------------------------------------------------#

id_map$SNU_ID <- sub(
  ".*?(AL\\d{3}-\\d{3})_S_(\\d{2})$",
  "\\1_\\2",
  id_map$SNU_ID
)


m_met <- met_raw %>%
  dplyr::rename(DNA_ID = SampleID) %>%
  dplyr::inner_join(id_map, by = "DNA_ID") %>%
  dplyr::inner_join(
    m %>%
      dplyr::mutate(
        SNU_tp = ifelse(
          TNT == "Before",
          paste0(SNU_ID, "_01"),
          paste0(SNU_ID, "_02")
        )
      ),
    by = c("SNU_ID" = "SNU_tp")
  )

#-----------------------------------------------------------------#

# 3. Build RNA metadata table

#-----------------------------------------------------------------#

m_rna <- m %>%
  dplyr::mutate(SNU_AL_ID = SNU_ID) %>%
  dplyr::inner_join(
    rna_map %>%
      dplyr::select(-dplyr::any_of(c("sex", "age"))),
    by = "SNU_AL_ID"
  ) %>%
  dplyr::mutate(
    Chart_matching = ifelse(
      TNT == "Before",
      paste0(Chart_numb, "_1"),
      paste0(Chart_numb, "_2")
    )
  ) %>%
  dplyr::inner_join(
    rna_meta %>%
      dplyr::select(RNA_sample_id, Chart_numb, timepoint) %>%
      dplyr::mutate(
        Chart_matching = ifelse(
          tolower(timepoint) == "post",
          paste0(Chart_numb, "_2"),
          paste0(Chart_numb, "_1")
        )
      ) %>%
      dplyr::select(-Chart_numb),
    by = "Chart_matching"
  ) %>%
  dplyr::mutate(
    timepoint = tolower(as.character(timepoint)),
    TRG_plot = dplyr::case_when(
      TRG_1 == "CR" ~ "pCR",
      TRG_1 == "nonCR" ~ "non_pCR",
      TRUE ~ NA_character_
    ),
    TRG_plot = factor(TRG_plot, levels = c("non_pCR", "pCR")),
    SNU_tp = ifelse(
      TNT == "Before",
      paste0(SNU_ID, "_01"),
      paste0(SNU_ID, "_02")
    )
  ) %>%
  dplyr::filter(
    timepoint != "normal",
    TNT == "Before",
    !is.na(TRG_plot)
  ) %>%
  dplyr::distinct(RNA_sample_id, .keep_all = TRUE)


#-----------------------------------------------------------------#

# 4. Build RNA-metabolite matched metadata table

#-----------------------------------------------------------------#

m_int <- m_rna %>%
  dplyr::inner_join(
    m_met %>%
      dplyr::filter(TNT == "Before") %>%
      dplyr::select(
        SNU_tp = SNU_ID,
        DNA_ID,
        Sample_weight_g,
        DW_ul,
        Time,
        dplyr::all_of(met_vars)
      ),
    by = "SNU_tp"
  ) %>%
  dplyr::distinct(RNA_sample_id, .keep_all = TRUE)

rna_no_met <- m_rna %>%
  dplyr::anti_join(
    m_int %>%
      dplyr::distinct(RNA_sample_id),
    by = "RNA_sample_id"
  ) %>%
  dplyr::select(
    dplyr::any_of(c(
      "SampleID", "SNU_ID", "SNU_tp", "SubjectID",
      "Chart_numb", "RNA_sample_id", "TRG_plot"
    ))
  )

message("RNA metadata samples: ", dplyr::n_distinct(m_rna$RNA_sample_id))
message("RNA-metabolite matched samples: ", dplyr::n_distinct(m_int$RNA_sample_id))
message("RNA samples without matched metabolite data: ", nrow(rna_no_met))

print(
  m_rna %>%
    dplyr::count(TRG_plot, name = "n")
)

print(
  m_int %>%
    dplyr::count(TRG_plot, name = "n")
)

if (nrow(rna_no_met) > 0) {
  print(rna_no_met)
}

#-----------------------------------------------------------------#

# 5. Build RNA count matrix for RNA-only DEG analysis

# This keeps all RNA baseline tumor samples, not only metabolite-matched samples.

#-----------------------------------------------------------------#

rna_cols <- m_rna %>%
  dplyr::distinct(RNA_sample_id) %>%
  dplyr::mutate(count_col = paste0(RNA_sample_id, "_Read_Count")) %>%
  dplyr::filter(count_col %in% colnames(rna))

rna_no_count <- m_rna %>%
  dplyr::filter(!RNA_sample_id %in% rna_cols$RNA_sample_id) %>%
  dplyr::select(
    dplyr::any_of(c(
      "SampleID", "SNU_ID", "SNU_tp", "SubjectID",
      "Chart_numb", "RNA_sample_id", "TRG_plot"
    ))
  )

if (nrow(rna_cols) == 0) {
  stop("No matched RNA_sample_id found in RNA count columns.")
}

message("RNA metadata samples found in count matrix: ", nrow(rna_cols))

if (nrow(rna_no_count) > 0) {
  warning("Some RNA metadata samples were not found in the RNA count matrix.")
  print(rna_no_count)
}

cnt_rna <- rna %>%
  dplyr::mutate(
    Gene = dplyr::case_when(
      !is.na(Gene_Symbol) & Gene_Symbol != "" ~ as.character(Gene_Symbol),
      !is.na(Gene_ID) & Gene_ID != "" ~ as.character(Gene_ID),
      TRUE ~ as.character(Transcript_ID)
    )
  ) %>%
  dplyr::filter(!is.na(Gene), Gene != "") %>%
  dplyr::select(Gene, dplyr::all_of(rna_cols$count_col)) %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(rna_cols$count_col),
      ~ sum(as.numeric(.x), na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()

cnt_rna <- round(cnt_rna)
mode(cnt_rna) <- "integer"

colnames(cnt_rna) <- rna_cols$RNA_sample_id[
  match(colnames(cnt_rna), rna_cols$count_col)
]

col_rna <- m_rna %>%
  dplyr::filter(RNA_sample_id %in% colnames(cnt_rna)) %>%
  dplyr::distinct(RNA_sample_id, .keep_all = TRUE) %>%
  dplyr::arrange(match(RNA_sample_id, colnames(cnt_rna))) %>%
  as.data.frame()

cnt_rna <- cnt_rna[, col_rna$RNA_sample_id, drop = FALSE]

col_rna$TRG_plot <- factor(col_rna$TRG_plot, levels = c("non_pCR", "pCR"))
rownames(col_rna) <- col_rna$RNA_sample_id

stopifnot(identical(colnames(cnt_rna), rownames(col_rna)))

message("RNA-only count matrix: ", nrow(cnt_rna), " genes x ", ncol(cnt_rna), " samples")

print(
  col_rna %>%
    dplyr::count(TRG_plot, name = "n")
)

#-----------------------------------------------------------------#

# 6. Build matched RNA and metabolite matrices for integration

#-----------------------------------------------------------------#

col_int <- m_int %>%
  dplyr::filter(RNA_sample_id %in% colnames(cnt_rna)) %>%
  dplyr::distinct(RNA_sample_id, .keep_all = TRUE) %>%
  dplyr::arrange(match(RNA_sample_id, colnames(cnt_rna))) %>%
  as.data.frame()

cnt_int <- cnt_rna[, col_int$RNA_sample_id, drop = FALSE]

col_int$TRG_plot <- factor(col_int$TRG_plot, levels = c("non_pCR", "pCR"))
rownames(col_int) <- col_int$RNA_sample_id

stopifnot(identical(colnames(cnt_int), rownames(col_int)))

met_missing <- setdiff(met_vars, colnames(col_int))

if (length(met_missing) > 0) {
  stop(
    "The following metabolite variables are missing from col_int: ",
    paste(met_missing, collapse = ", ")
  )
}

met_int <- col_int[, met_vars, drop = FALSE] %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      ~ suppressWarnings(as.numeric(.x))
    )
  ) %>%
  as.matrix()

rownames(met_int) <- col_int$RNA_sample_id

stopifnot(identical(rownames(met_int), colnames(cnt_int)))

message("Matched RNA count matrix: ", nrow(cnt_int), " genes x ", ncol(cnt_int), " samples")
message("Matched metabolite matrix: ", nrow(met_int), " samples x ", ncol(met_int), " metabolites")

print(
  col_int %>%
    dplyr::count(TRG_plot, name = "n")
)


#-----------------------------------------------------------------#

# 7. Save integrated input object

#-----------------------------------------------------------------#
save(
  
  m, # Core microbiome metadata # 42 samples
  s, # Species abundance table
  t, # Strain abundance table
  m_met, # Metabolite metadata table # 32 samples
  m_rna, # RNA-seq metadata table # 20 samples
  m_int, # RNA-metabolite matched metadata table # 17 samples
  met_vars, # Metabolite feature names
  rna_cols, # RNA sample-to-count-column mapping table
  cnt_rna,  # RNA count matrix for RNA-only DEG analysis
  col_rna,  # DESeq2-compatible colData for RNA-only DEG analysis
  cnt_int,  # RNA count matrix for RNA-metabolite matched analysis
  col_int,  # DESeq2-compatible colData for RNA-metabolite matched analysis
  met_int,  # Matched metabolite matrix
  
  file = "host_RNAseq/results_clean_metabolite_host/host_met_rna_matched_input.RData"
)


load("host_RNAseq/results_clean_metabolite_host/host_met_rna_matched_input.RData")

ls()


