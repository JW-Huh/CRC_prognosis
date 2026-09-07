#-----------------------------------------------------------------#
#
# Host RNA-seq pathway analysis linked to metabolite profiles
# Baseline tumor samples, pCR vs non-pCR
#
# Input:
#   host_RNAseq/results_clean_metabolite_host/host_rna_DESeq2_pCR_vs_non_pCR.RData
#
# Required objects in the input RData:
#   deg_rna, vst_rna, vst_int, col_rna, col_int, met_int
#
# Main outputs:
#   1) Official and targeted gene-set database table
#   2) Preranked GSEA using DESeq2 statistic
#   3) Over-representation analysis using DE genes
#   4) sample-level ssGSEA/GSVA scores
#   5) pCR vs non-pCR comparison of pathway scores
#   6) metabolite-pathway score correlation
#
#-----------------------------------------------------------------#

rm(list = ls())
options(stringsAsFactors = FALSE)
setwd("D:/2-연구/2-CRC metagenomics/")

#-----------------------------------------------------------------#
# 0. Install packages if needed
#-----------------------------------------------------------------#

# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# 
# cran_pkgs <- c(
#   "dplyr", "tidyr", "tibble", "stringr", "purrr", "readr",
#   "ggplot2", "ggrepel", "ggbeeswarm", "forcats", "msigdbr"
# )
# 
# bioc_pkgs <- c(
#   "AnnotationDbi", "org.Hs.eg.db", "clusterProfiler", "fgsea", "GSVA"
# )
# 
# for (pkg in cran_pkgs) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     install.packages(pkg)
#   }
# }
# 
# for (pkg in bioc_pkgs) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     BiocManager::install(pkg, ask = FALSE, update = FALSE)
#   }
# }


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(ggbeeswarm)
  library(forcats)
  library(msigdbr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(fgsea)
  library(GSVA)
})




#-----------------------------------------------------------------#
# 1. Output folders
#-----------------------------------------------------------------#

dir.create("host_RNAseq/results_clean_metabolite_host/pathway_metabolite", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

#-----------------------------------------------------------------#
# 2. Load previous RNA-seq and matched metabolite input
#-----------------------------------------------------------------#

load("host_RNAseq/results_clean_metabolite_host/host_rna_DESeq2_pCR_vs_non_pCR.RData")

stopifnot(exists("deg_rna"))
stopifnot(exists("vst_rna"))
stopifnot(exists("vst_int"))
stopifnot(exists("col_rna"))
stopifnot(exists("col_int"))
stopifnot(exists("met_int"))

col_rna <- as.data.frame(col_rna)
col_int <- as.data.frame(col_int)
met_int <- as.data.frame(met_int)

col_rna$TRG_plot <- factor(col_rna$TRG_plot, levels = c("non_pCR", "pCR"))
col_int$TRG_plot <- factor(col_int$TRG_plot, levels = c("non_pCR", "pCR"))

if (!identical(colnames(vst_rna), rownames(col_rna))) {
  rownames(col_rna) <- col_rna$RNA_sample_id
  vst_rna <- vst_rna[, rownames(col_rna), drop = FALSE]
}

if (!identical(colnames(vst_int), rownames(col_int))) {
  rownames(col_int) <- col_int$RNA_sample_id
  vst_int <- vst_int[, rownames(col_int), drop = FALSE]
}

stopifnot(identical(colnames(vst_rna), rownames(col_rna)))
stopifnot(identical(colnames(vst_int), rownames(col_int)))
stopifnot(identical(rownames(met_int), colnames(vst_int)))

met_int <- met_int[, vapply(met_int, is.numeric, logical(1)), drop = FALSE]
met_int <- met_int[, colSums(is.finite(as.matrix(met_int))) >= 4, drop = FALSE]

message("RNA VST matrix: ", nrow(vst_rna), " genes x ", ncol(vst_rna), " samples")
message("Matched RNA VST matrix: ", nrow(vst_int), " genes x ", ncol(vst_int), " samples")
message("Matched metabolite matrix: ", nrow(met_int), " samples x ", ncol(met_int), " metabolites")

#-----------------------------------------------------------------#
# 3. Helper functions
#-----------------------------------------------------------------#

clean_ensembl_version <- function(x) {
  gsub("\\..*$", "", x)
}

collapse_matrix_by_gene <- function(mat, gene_symbol) {
  keep <- !is.na(gene_symbol) & gene_symbol != ""
  mat <- mat[keep, , drop = FALSE]
  gene_symbol <- gene_symbol[keep]
  mat_sum <- rowsum(mat, group = gene_symbol, reorder = FALSE)
  gene_n <- as.numeric(table(factor(gene_symbol, levels = rownames(mat_sum))))
  mat_sum / gene_n
}

standardize_gene_symbols <- function(vst_rna, vst_int, deg_rna) {
  gene_id <- rownames(vst_rna)
  gene_id_clean <- clean_ensembl_version(gene_id)
  is_ensembl <- mean(grepl("^ENSG", gene_id_clean)) > 0.50
  
  if (is_ensembl) {
    gene_map <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unique(gene_id_clean),
      keytype = "ENSEMBL",
      columns = c("ENSEMBL", "SYMBOL", "ENTREZID")
    ) %>%
      dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
      dplyr::distinct(ENSEMBL, SYMBOL, .keep_all = TRUE)
    
    gene_symbol <- gene_map$SYMBOL[match(gene_id_clean, gene_map$ENSEMBL)]
    vst_rna_symbol <- collapse_matrix_by_gene(vst_rna, gene_symbol)
    vst_int_symbol <- collapse_matrix_by_gene(vst_int, gene_symbol[match(rownames(vst_int), rownames(vst_rna))])
    
    deg_rna_symbol <- deg_rna %>%
      dplyr::mutate(Gene_input = Gene, Gene_clean = clean_ensembl_version(Gene)) %>%
      dplyr::left_join(gene_map %>% dplyr::select(Gene_clean = ENSEMBL, Gene_symbol = SYMBOL), by = "Gene_clean") %>%
      dplyr::mutate(Gene = Gene_symbol) %>%
      dplyr::filter(!is.na(Gene), Gene != "") %>%
      dplyr::group_by(Gene) %>%
      dplyr::arrange(dplyr::desc(abs(stat_DESeq2)), .by_group = TRUE) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup()
  } else {
    gene_symbol <- gene_id
    vst_rna_symbol <- vst_rna
    vst_int_symbol <- vst_int
    rownames(vst_rna_symbol) <- gene_symbol
    rownames(vst_int_symbol) <- rownames(vst_int)
    
    if (any(duplicated(rownames(vst_rna_symbol)))) {
      vst_rna_symbol <- collapse_matrix_by_gene(vst_rna_symbol, rownames(vst_rna_symbol))
    }
    if (any(duplicated(rownames(vst_int_symbol)))) {
      vst_int_symbol <- collapse_matrix_by_gene(vst_int_symbol, rownames(vst_int_symbol))
    }
    
    deg_rna_symbol <- deg_rna %>%
      dplyr::filter(!is.na(Gene), Gene != "") %>%
      dplyr::group_by(Gene) %>%
      dplyr::arrange(dplyr::desc(abs(stat_DESeq2)), .by_group = TRUE) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup()
  }
  
  list(
    vst_rna = vst_rna_symbol,
    vst_int = vst_int_symbol,
    deg_rna = deg_rna_symbol
  )
}

run_ssgsea <- function(expr_mat, gene_sets) {
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"
  
  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    param <- GSVA::ssgseaParam(
      exprData = expr_mat,
      geneSets = gene_sets,
      normalize = TRUE
    )
    score_mat <- GSVA::gsva(param, verbose = FALSE)
  } else {
    score_mat <- GSVA::gsva(
      expr = expr_mat,
      gset.idx.list = gene_sets,
      method = "ssgsea",
      kcdf = "Gaussian",
      abs.ranking = FALSE,
      ssgsea.norm = TRUE,
      verbose = FALSE
    )
  }
  score_mat
}

zscore_rows <- function(mat) {
  mat_z <- t(scale(t(mat)))
  mat_z[!is.finite(mat_z)] <- NA_real_
  mat_z
}

safe_wilcox <- function(x, g) {
  x1 <- x[g == "pCR"]
  x0 <- x[g == "non_pCR"]
  x1 <- x1[is.finite(x1)]
  x0 <- x0[is.finite(x0)]
  if (length(x1) < 2 || length(x0) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(x1, x0, exact = FALSE)$p.value)
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 4) {
    return(tibble::tibble(n = sum(keep), rho = NA_real_, pval = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  tibble::tibble(n = sum(keep), rho = unname(ct$estimate), pval = ct$p.value)
}

#-----------------------------------------------------------------#
# 4. Standardize gene identifiers to HGNC symbols
#-----------------------------------------------------------------#

std <- standardize_gene_symbols(vst_rna, vst_int, deg_rna)
vst_rna <- std$vst_rna
vst_int <- std$vst_int
deg_rna <- std$deg_rna
rm(std)

write.csv(
  tibble::tibble(Gene = rownames(vst_rna)),
  "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/expression_gene_symbols_used.csv",
  row.names = FALSE
)

#-----------------------------------------------------------------#
# 5. Download/load MSigDB gene sets via msigdbr
#-----------------------------------------------------------------#

msig <- msigdbr::msigdbr(species = "Homo sapiens")

if (!"gs_collection" %in% colnames(msig) && "gs_cat" %in% colnames(msig)) {
  msig <- msig %>% dplyr::rename(gs_collection = gs_cat)
}
if (!"gs_subcollection" %in% colnames(msig) && "gs_subcat" %in% colnames(msig)) {
  msig <- msig %>% dplyr::rename(gs_subcollection = gs_subcat)
}
if (!"db_version" %in% colnames(msig)) {
  msig$db_version <- NA_character_
}
if (!"gs_description" %in% colnames(msig)) {
  msig$gs_description <- NA_character_
}
if (!"gs_exact_source" %in% colnames(msig)) {
  msig$gs_exact_source <- NA_character_
}

msig <- msig %>%
  dplyr::select(
    gs_name,
    Gene = gene_symbol,
    gs_collection,
    gs_subcollection,
    gs_description,
    gs_exact_source,
    db_version
  ) %>%
  dplyr::filter(!is.na(Gene), Gene != "") %>%
  dplyr::distinct()

message("MSigDB version detected by msigdbr: ", paste(unique(msig$db_version), collapse = ", "))

write.csv(
  msig %>% dplyr::distinct(gs_collection, gs_subcollection, db_version) %>% dplyr::arrange(gs_collection, gs_subcollection),
  "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/MSigDB_collections_detected.csv",
  row.names = FALSE
)

#-----------------------------------------------------------------#
# 6. Preferred high-confidence gene sets
#-----------------------------------------------------------------#

preferred_gs <- tibble::tribble(
  ~Host_axis, ~gs_name,
  
  "Immune / inflammatory response", "HALLMARK_INFLAMMATORY_RESPONSE",
  "Immune / inflammatory response", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "Immune / inflammatory response", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "Immune / inflammatory response", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "Immune / inflammatory response", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "Immune / inflammatory response", "HALLMARK_COMPLEMENT",
  
  "Innate immune receptor downstream signaling", "REACTOME_TOLL_LIKE_RECEPTOR_CASCADES",
  "Innate immune receptor downstream signaling", "REACTOME_TOLL_LIKE_RECEPTOR_4_TLR4_CASCADE",
  "Innate immune receptor downstream signaling", "REACTOME_MYD88_DEPENDENT_TLR4_CASCADE",
  "Innate immune receptor downstream signaling", "REACTOME_TRIF_MEDIATED_TLR3_TLR4_SIGNALING",
  "Innate immune receptor downstream signaling", "REACTOME_NOD1_2_SIGNALING_PATHWAY",
  "Innate immune receptor downstream signaling", "REACTOME_NUCLEOTIDE_BINDING_DOMAIN_LEUCINE_RICH_REPEAT_CONTAINING_RECEPTOR_NLR_SIGNALING_PATHWAYS",
  "Innate immune receptor downstream signaling", "GOBP_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN",
  "Innate immune receptor downstream signaling", "GOBP_PATTERN_RECOGNITION_RECEPTOR_SIGNALING_PATHWAY",
  "Innate immune receptor downstream signaling", "GOBP_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "Innate immune receptor downstream signaling", "GOBP_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  "Innate immune receptor downstream signaling", "GOBP_INFLAMMASOME_COMPLEX_ASSEMBLY",
  
  "SCFA / butyrate / fatty-acid response", "GOBP_CELLULAR_RESPONSE_TO_BUTYRATE",
  "SCFA / butyrate / fatty-acid response", "GOBP_RESPONSE_TO_BUTYRATE",
  "SCFA / butyrate / fatty-acid response", "GOBP_CELLULAR_RESPONSE_TO_FATTY_ACID",
  "SCFA / butyrate / fatty-acid response", "GOBP_RESPONSE_TO_FATTY_ACID",
  "SCFA / butyrate / fatty-acid response", "GOBP_SHORT_CHAIN_FATTY_ACID_METABOLIC_PROCESS",
  "SCFA / butyrate / fatty-acid response", "HALLMARK_FATTY_ACID_METABOLISM",
  
  "Niacin / nicotinamide / NAD response", "GOMF_NICOTINIC_ACID_RECEPTOR_ACTIVITY",
  "Niacin / nicotinamide / NAD response", "GOBP_NICOTINAMIDE_NUCLEOTIDE_METABOLIC_PROCESS",
  "Niacin / nicotinamide / NAD response", "GOBP_NAD_METABOLIC_PROCESS",
  "Niacin / nicotinamide / NAD response", "REACTOME_METABOLISM_OF_WATER_SOLUBLE_VITAMINS_AND_COFACTORS",
  
  "Indole / AhR / xenobiotic response", "REACTOME_ARYL_HYDROCARBON_RECEPTOR_SIGNALLING",
  "Indole / AhR / xenobiotic response", "WP_ARYL_HYDROCARBON_RECEPTOR_PATHWAY_WP2586",
  "Indole / AhR / xenobiotic response", "HALLMARK_XENOBIOTIC_METABOLISM",
  "Indole / AhR / xenobiotic response", "GOBP_CELLULAR_RESPONSE_TO_XENOBIOTIC_STIMULUS",
  "Indole / AhR / xenobiotic response", "GOBP_RESPONSE_TO_XENOBIOTIC_STIMULUS",
  "Indole / AhR / xenobiotic response", "REACTOME_PHASE_I_FUNCTIONALIZATION_OF_COMPOUNDS",
  
  "Bile acid receptor / bile-acid response", "HALLMARK_BILE_ACID_METABOLISM",
  "Bile acid receptor / bile-acid response", "REACTOME_BILE_ACID_AND_BILE_SALT_METABOLISM",
  "Bile acid receptor / bile-acid response", "REACTOME_SYNTHESIS_OF_BILE_ACIDS_AND_BILE_SALTS",
  "Bile acid receptor / bile-acid response", "GOBP_CELLULAR_RESPONSE_TO_BILE_ACID",
  "Bile acid receptor / bile-acid response", "GOBP_RESPONSE_TO_BILE_ACID",
  "Bile acid receptor / bile-acid response", "GOBP_BILE_ACID_SECRETION",
  
  "Th1 / Th2 / Th17 / Treg lineage", "GOBP_T_HELPER_1_CELL_DIFFERENTIATION",
  "Th1 / Th2 / Th17 / Treg lineage", "GOBP_T_HELPER_2_CELL_DIFFERENTIATION",
  "Th1 / Th2 / Th17 / Treg lineage", "GOBP_T_HELPER_17_CELL_DIFFERENTIATION",
  "Th1 / Th2 / Th17 / Treg lineage", "GOBP_REGULATORY_T_CELL_DIFFERENTIATION",
  "Th1 / Th2 / Th17 / Treg lineage", "GOBP_ALPHA_BETA_T_CELL_DIFFERENTIATION",
  "Th1 / Th2 / Th17 / Treg lineage", "GOBP_T_CELL_DIFFERENTIATION",
  "Th1 / Th2 / Th17 / Treg lineage", "GOBP_T_CELL_ACTIVATION",
  
  "Cytokines / chemokines", "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
  "Cytokines / chemokines", "REACTOME_SIGNALING_BY_INTERLEUKINS",
  "Cytokines / chemokines", "REACTOME_INTERLEUKIN_1_SIGNALING",
  "Cytokines / chemokines", "REACTOME_INTERLEUKIN_6_SIGNALING",
  "Cytokines / chemokines", "REACTOME_INTERLEUKIN_10_SIGNALING",
  "Cytokines / chemokines", "REACTOME_INTERLEUKIN_17_SIGNALING",
  "Cytokines / chemokines", "REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES",
  "Cytokines / chemokines", "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY",
  "Cytokines / chemokines", "GOBP_CHEMOKINE_MEDIATED_SIGNALING_PATHWAY",
  
  "Gut homing / integrin / lymphocyte trafficking", "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
  "Gut homing / integrin / lymphocyte trafficking", "GOBP_INTEGRIN_MEDIATED_SIGNALING_PATHWAY",
  "Gut homing / integrin / lymphocyte trafficking", "GOBP_LEUKOCYTE_MIGRATION",
  "Gut homing / integrin / lymphocyte trafficking", "GOBP_LYMPHOCYTE_MIGRATION",
  "Gut homing / integrin / lymphocyte trafficking", "GOBP_T_CELL_MIGRATION",
  "Gut homing / integrin / lymphocyte trafficking", "GOBP_LEUKOCYTE_CELL_CELL_ADHESION",
  
  "Barrier / tight junction / epithelial integrity", "HALLMARK_APICAL_JUNCTION",
  "Barrier / tight junction / epithelial integrity", "HALLMARK_APICAL_SURFACE",
  "Barrier / tight junction / epithelial integrity", "GOBP_TIGHT_JUNCTION_ASSEMBLY",
  "Barrier / tight junction / epithelial integrity", "GOBP_CELL_CELL_JUNCTION_ASSEMBLY",
  "Barrier / tight junction / epithelial integrity", "GOBP_EPITHELIAL_CELL_DIFFERENTIATION",
  "Barrier / tight junction / epithelial integrity", "GOBP_INTESTINAL_EPITHELIAL_CELL_DIFFERENTIATION",
  
  "Mucus / goblet-cell / antimicrobial peptide program", "GOBP_REGULATION_OF_MUCUS_SECRETION",
  "Mucus / goblet-cell / antimicrobial peptide program", "GOBP_MUCUS_SECRETION",
  "Mucus / goblet-cell / antimicrobial peptide program", "GOBP_GOBLET_CELL_DIFFERENTIATION",
  "Mucus / goblet-cell / antimicrobial peptide program", "GOBP_ANTIMICROBIAL_HUMORAL_RESPONSE",
  "Mucus / goblet-cell / antimicrobial peptide program", "GOBP_DEFENSE_RESPONSE_TO_BACTERIUM",
  
  "Tumor remodeling / inflammation-linked CRC program", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "Tumor remodeling / inflammation-linked CRC program", "HALLMARK_TGF_BETA_SIGNALING",
  "Tumor remodeling / inflammation-linked CRC program", "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "Tumor remodeling / inflammation-linked CRC program", "HALLMARK_MYC_TARGETS_V1",
  "Tumor remodeling / inflammation-linked CRC program", "HALLMARK_MYC_TARGETS_V2",
  "Tumor remodeling / inflammation-linked CRC program", "HALLMARK_E2F_TARGETS",
  "Tumor remodeling / inflammation-linked CRC program", "HALLMARK_G2M_CHECKPOINT",
  "Tumor remodeling / inflammation-linked CRC program", "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "Tumor remodeling / inflammation-linked CRC program", "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
  "Tumor remodeling / inflammation-linked CRC program", "GOBP_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "Tumor remodeling / inflammation-linked CRC program", "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION",
  "Tumor remodeling / inflammation-linked CRC program", "GOBP_WNT_SIGNALING_PATHWAY",
  "Tumor remodeling / inflammation-linked CRC program", "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY",
  
  "Angiogenesis / hypoxia", "HALLMARK_ANGIOGENESIS",
  "Angiogenesis / hypoxia", "HALLMARK_HYPOXIA",
  "Angiogenesis / hypoxia", "GOBP_ANGIOGENESIS",
  "Angiogenesis / hypoxia", "GOBP_BLOOD_VESSEL_MORPHOGENESIS",
  "Angiogenesis / hypoxia", "GOBP_RESPONSE_TO_HYPOXIA",
  "Angiogenesis / hypoxia", "REACTOME_SIGNALING_BY_VEGF"
)

preferred_found <- msig %>%
  dplyr::inner_join(preferred_gs, by = "gs_name") %>%
  dplyr::mutate(Source_type = "Preferred_MSigDB")

preferred_missing <- preferred_gs %>%
  dplyr::anti_join(msig %>% dplyr::distinct(gs_name), by = "gs_name")

# write.csv(
#   preferred_missing,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/preferred_gene_sets_not_found_in_current_MSigDB.csv",
#   row.names = FALSE
# )

#-----------------------------------------------------------------#
# 7. Keyword-based rescue from MSigDB H/C2/C5/C7/C8
#-----------------------------------------------------------------#

axis_patterns <- tibble::tribble(
  ~Host_axis, ~pattern,
  "Immune / inflammatory response", "IMMUNE|INFLAMMATORY|INFLAMMATION|INTERFERON|COMPLEMENT|NF_KB|NFKB|TNFA|IL6|JAK_STAT",
  "Innate immune receptor downstream signaling", "TOLL_LIKE|TLR|NOD|NLR|INFLAMMASOME|PATTERN_RECOGNITION|BACTERIAL_ORIGIN|LIPOPOLYSACCHARIDE|MYD88|TICAM|TRIF",
  "SCFA / butyrate / fatty-acid response", "BUTYRATE|PROPIONATE|ACETATE|SHORT_CHAIN|FATTY_ACID|FATTY ACID|FFAR|HCAR|SLC5A8|MONOCARBOXYLATE",
  "Niacin / nicotinamide / NAD response", "NIACIN|NICOTINIC|NICOTINAMIDE|NAD|NADP|WATER_SOLUBLE_VITAMIN|VITAMIN_B3",
  "Indole / AhR / xenobiotic response", "ARYL_HYDROCARBON|AHR|XENOBIOTIC|TRYPTOPHAN|INDOLE|CYP1A1|CYP1B1|PXR|NR1I2",
  "Bile acid receptor / bile-acid response", "BILE_ACID|BILE ACID|BILE_SALT|BILE SALT|FXR|NR1H4|TGR5|GPBAR1|VDR|PXR|NR1I2",
  "Th1 / Th2 / Th17 / Treg lineage", "T_HELPER|TH1|TH2|TH17|REGULATORY_T_CELL|TREG|FOXP3|RORC|TBX21|GATA3|T_CELL_DIFFERENTIATION",
  "Cytokines / chemokines", "CYTOKINE|CHEMOKINE|INTERLEUKIN|IL1|IL6|IL10|IL17|IL22|TNF|IFN|CXCL|CCL",
  "Gut homing / integrin / lymphocyte trafficking", "INTEGRIN|LEUKOCYTE_MIGRATION|LYMPHOCYTE_MIGRATION|T_CELL_MIGRATION|CELL_CELL_ADHESION|TRAFFICKING|HOMING",
  "Barrier / tight junction / epithelial integrity", "TIGHT_JUNCTION|APICAL_JUNCTION|APICAL_SURFACE|CELL_CELL_JUNCTION|EPITHELIAL_BARRIER|EPITHELIAL_CELL_DIFFERENTIATION|ADHERENS",
  "Mucus / goblet-cell / antimicrobial peptide program", "MUCUS|MUCIN|GOBLET|ANTIMICROBIAL|DEFENSIN|DEFENSE_RESPONSE_TO_BACTERIUM|ANTIBACTERIAL",
  "Tumor remodeling / inflammation-linked CRC program", "EPITHELIAL_MESENCHYMAL|EXTRACELLULAR_MATRIX|MATRIX|COLLAGEN|MMP|WNT|BETA_CATENIN|MYC|E2F|G2M|KRAS|COLORECTAL|COLON_CANCER|TGF_BETA",
  "Angiogenesis / hypoxia", "ANGIOGENESIS|HYPOXIA|VEGF|BLOOD_VESSEL|VASCULATURE"
)

# Default keyword rescue is limited to Hallmark, C2 curated pathways, and C5 ontology.
# C7/C8 immune-cell signatures can be useful, but unrestricted keyword matching can
# return thousands of perturbation/cell-type signatures and slow down ssGSEA.
msig_for_keyword <- msig %>%
  dplyr::filter(
    gs_collection %in% c("H", "C2", "C5")
  ) %>%
  dplyr::mutate(
    text_for_search = paste(gs_name, gs_description, gs_exact_source, sep = " ")
  )

keyword_found <- purrr::map_dfr(seq_len(nrow(axis_patterns)), function(i) {
  msig_for_keyword %>%
    dplyr::filter(stringr::str_detect(text_for_search, stringr::regex(axis_patterns$pattern[i], ignore_case = TRUE))) %>%
    dplyr::mutate(
      Host_axis = axis_patterns$Host_axis[i],
      Source_type = "Keyword_MSigDB_rescue"
    )
}) %>%
  dplyr::select(-text_for_search) %>%
  dplyr::distinct()


#-----------------------------------------------------------------#
# 8. Direct GO-derived metabolite-response terms
#    These are added to prevent SCFA/butyrate or niacin-related terms from being lost
#    if the exact term is absent in the local MSigDB build.
#-----------------------------------------------------------------#

go_direct <- tibble::tribble(
  ~Host_axis, ~GO_ID, ~GO_name,
  "SCFA / butyrate / fatty-acid response", "GO:1903545", "cellular response to butyrate",
  "SCFA / butyrate / fatty-acid response", "GO:0071398", "cellular response to fatty acid",
  "SCFA / butyrate / fatty-acid response", "GO:0070542", "response to fatty acid",
  "Bile acid receptor / bile-acid response", "GO:1903413", "cellular response to bile acid",
  "Niacin / nicotinamide / NAD response", "GO:0070553", "nicotinic acid receptor activity",
  "Niacin / nicotinamide / NAD response", "GO:0046496", "nicotinamide nucleotide metabolic process",
  "Niacin / nicotinamide / NAD response", "GO:0019674", "NAD metabolic process"
)

go_direct_anno <- tryCatch(
  suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unique(go_direct$GO_ID),
      keytype = "GOALL",
      columns = c("SYMBOL", "ENTREZID", "GOALL", "EVIDENCEALL", "ONTOLOGYALL")
    )
  ),
  error = function(e) data.frame()
) %>%
  tibble::as_tibble() %>%
  dplyr::rename(
    GO_ID = GOALL,
    Evidence = EVIDENCEALL,
    Ontology = ONTOLOGYALL
  ) %>%
  dplyr::filter(!is.na(SYMBOL), !is.na(GO_ID)) %>%
  dplyr::distinct(GO_ID, SYMBOL, ENTREZID, Evidence, Ontology)

if (nrow(go_direct_anno) > 0) {
  go_direct_sets <- go_direct_anno %>%
    tibble::as_tibble() %>%
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
    dplyr::inner_join(
      go_direct %>%
        dplyr::distinct(GO_ID, Host_axis, GO_name),
      by = "GO_ID",
      relationship = "many-to-many"
    ) %>%
    dplyr::transmute(
      Host_axis,
      gs_name = paste0(
        "GO_DIRECT_",
        stringr::str_replace_all(toupper(GO_name), "[^A-Z0-9]+", "_")
      ),
      Gene = SYMBOL,
      gs_collection = "GO_DIRECT",
      gs_subcollection = Ontology,
      gs_description = GO_name,
      gs_exact_source = GO_ID,
      db_version = as.character(NA),
      Source_type = "Direct_GO_orgHs"
    ) %>%
    dplyr::distinct()
} else {
  go_direct_sets <- tibble::tibble(
    Host_axis = character(),
    gs_name = character(),
    Gene = character(),
    gs_collection = character(),
    gs_subcollection = character(),
    gs_description = character(),
    gs_exact_source = character(),
    db_version = character(),
    Source_type = character()
  )
}

#-----------------------------------------------------------------#
# 9. Optional minimal mechanistic marker panels
#    Keep these separate from official DB-derived sets.
#-----------------------------------------------------------------#

curated_marker_sets <- dplyr::bind_rows(
  tibble::tibble(
    Host_axis = "SCFA / butyrate / fatty-acid response",
    gs_name = "CURATED_SCFA_BUTYRATE_RECEPTOR_TRANSPORTER_MARKERS",
    Gene = c("HCAR2", "FFAR2", "FFAR3", "SLC5A8", "SLC16A1", "SLC16A3", "HDAC1", "HDAC2", "HDAC3"),
    gs_collection = "CURATED_MARKER",
    gs_subcollection = "Mechanism_informed",
    gs_description = "Minimal marker panel for SCFA/butyrate receptor, transporter, and HDAC-linked response",
    gs_exact_source = "Curated_marker_panel",
    db_version = as.character(NA),
    Source_type = "Curated_marker_sensitivity"
  ),
  tibble::tibble(
    Host_axis = "Niacin / nicotinamide / NAD response",
    gs_name = "CURATED_NIACIN_HCAR2_NAD_MARKERS",
    Gene = c("HCAR2", "NADSYN1", "NAMPT", "NAPRT", "NMNAT1", "NMNAT2", "NMNAT3", "NMRK1", "NMRK2", "SIRT1", "SIRT3"),
    gs_collection = "CURATED_MARKER",
    gs_subcollection = "Mechanism_informed",
    gs_description = "Minimal marker panel for nicotinic-acid receptor and NAD salvage/metabolism",
    gs_exact_source = "Curated_marker_panel",
    db_version = as.character(NA),
    Source_type = "Curated_marker_sensitivity"
  ),
  tibble::tibble(
    Host_axis = "Indole / AhR / xenobiotic response",
    gs_name = "CURATED_INDOLE_AHR_PXR_TARGET_MARKERS",
    Gene = c("AHR", "ARNT", "AHRR", "CYP1A1", "CYP1B1", "TIPARP", "NQO1", "ALDH1A1", "NR1I2", "CYP3A4", "ABCB1"),
    gs_collection = "CURATED_MARKER",
    gs_subcollection = "Mechanism_informed",
    gs_description = "Minimal marker panel for indole-linked AhR/PXR and xenobiotic response",
    gs_exact_source = "Curated_marker_panel",
    db_version = as.character(NA),
    Source_type = "Curated_marker_sensitivity"
  ),
  tibble::tibble(
    Host_axis = "Bile acid receptor / bile-acid response",
    gs_name = "CURATED_BILE_ACID_RECEPTOR_MARKERS",
    Gene = c("NR1H4", "GPBAR1", "NR0B2", "FGF19", "CYP3A4", "SLC10A2", "FABP6", "VDR", "NR1I2"),
    gs_collection = "CURATED_MARKER",
    gs_subcollection = "Mechanism_informed",
    gs_description = "Minimal marker panel for FXR, TGR5, VDR, PXR, and intestinal bile-acid response",
    gs_exact_source = "Curated_marker_panel",
    db_version = as.character(NA),
    Source_type = "Curated_marker_sensitivity"
  ),
  tibble::tibble(
    Host_axis = "Th1 / Th2 / Th17 / Treg lineage",
    gs_name = "CURATED_TH_MASTER_REGULATOR_MARKERS",
    Gene = c("TBX21", "GATA3", "RORC", "FOXP3", "STAT1", "STAT3", "STAT4", "STAT5A", "IL12RB2", "IL23R", "IL2RA"),
    gs_collection = "CURATED_MARKER",
    gs_subcollection = "Mechanism_informed",
    gs_description = "Minimal marker panel for Th1, Th2, Th17, and Treg master regulators",
    gs_exact_source = "Curated_marker_panel",
    db_version = as.character(NA),
    Source_type = "Curated_marker_sensitivity"
  ),
  tibble::tibble(
    Host_axis = "Gut homing / integrin / lymphocyte trafficking",
    gs_name = "CURATED_GUT_HOMING_INTEGRIN_MARKERS",
    Gene = c("ITGA4", "ITGB7", "ITGAE", "CCR6", "CCR9", "CXCR3", "CCL20", "MADCAM1"),
    gs_collection = "CURATED_MARKER",
    gs_subcollection = "Mechanism_informed",
    gs_description = "Minimal marker panel for gut-homing integrins and lymphocyte trafficking",
    gs_exact_source = "Curated_marker_panel",
    db_version = as.character(NA),
    Source_type = "Curated_marker_sensitivity"
  )
) %>%
  dplyr::distinct()

#-----------------------------------------------------------------#
# 10. Combine gene sets and filter by expression overlap
#-----------------------------------------------------------------#

all_gene_sets_long <- dplyr::bind_rows(
  preferred_found,
  keyword_found,
  go_direct_sets,
  curated_marker_sets
) %>%
  dplyr::distinct(Host_axis, gs_name, Gene, Source_type, .keep_all = TRUE)

all_gene_sets_long <- all_gene_sets_long %>%
  dplyr::mutate(
    Gene = dplyr::recode(Gene, GPR109A = "HCAR2"),
    pathway_id = paste(Host_axis, gs_name, sep = "__")
  ) %>%
  dplyr::distinct(Host_axis, gs_name, Gene, Source_type, .keep_all = TRUE)

gene_set_size_table <- all_gene_sets_long %>%
  dplyr::group_by(Host_axis, gs_name, Source_type, gs_collection, gs_subcollection, gs_description, gs_exact_source, db_version, pathway_id) %>%
  dplyr::summarise(
    n_total_genes = dplyr::n_distinct(Gene),
    n_genes_in_vst_rna = sum(unique(Gene) %in% rownames(vst_rna)),
    n_genes_in_vst_int = sum(unique(Gene) %in% rownames(vst_int)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Host_axis, Source_type, gs_name)

# write.csv(
#   all_gene_sets_long,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/all_targeted_gene_sets_long_unfiltered.csv",
#   row.names = FALSE
# )
# 
# write.csv(
#   gene_set_size_table,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/gene_set_size_and_expression_overlap.csv",
#   row.names = FALSE
# )

min_gs_size <- 5
max_gs_size <- 500

use_pathway_id <- gene_set_size_table %>%
  dplyr::filter(
    n_genes_in_vst_rna >= min_gs_size,
    n_genes_in_vst_rna <= max_gs_size
  ) %>%
  dplyr::pull(pathway_id)

gene_sets_use_long <- all_gene_sets_long %>%
  dplyr::filter(pathway_id %in% use_pathway_id, Gene %in% rownames(vst_rna)) %>%
  dplyr::distinct(pathway_id, Host_axis, gs_name, Gene, Source_type, .keep_all = TRUE)

gene_sets_use <- split(gene_sets_use_long$Gene, gene_sets_use_long$pathway_id)
gene_sets_use <- lapply(gene_sets_use, unique)

# write.csv(
#   gene_sets_use_long,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/gene_sets_use_long.csv",
#   row.names = FALSE
# )

message("Gene sets used after expression overlap filtering: ", length(gene_sets_use))

#-----------------------------------------------------------------#
# 11. Pre-ranked GSEA using DESeq2 statistic
#-----------------------------------------------------------------#

rank_tbl <- deg_rna %>%
  dplyr::filter(!is.na(Gene), Gene %in% rownames(vst_rna), !is.na(stat_DESeq2)) %>%
  dplyr::group_by(Gene) %>%
  dplyr::arrange(dplyr::desc(abs(stat_DESeq2)), .by_group = TRUE) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(dplyr::desc(stat_DESeq2))

rank_stat <- rank_tbl$stat_DESeq2
  names(rank_stat) <- rank_tbl$Gene
rank_stat <- sort(rank_stat, decreasing = TRUE)
rank_stat <- rank_stat[is.finite(rank_stat)]

fgsea_res <- fgsea::fgsea(
  pathways = gene_sets_use,
  stats = rank_stat,
  minSize = min_gs_size,
  maxSize = max_gs_size,
  eps = 0
  ) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    pathway_id = pathway,
    leadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  ) %>%
  dplyr::left_join(
    gene_set_size_table %>%
      dplyr::select(pathway_id, Host_axis, gs_name, Source_type, gs_collection, gs_subcollection, gs_description, gs_exact_source, n_total_genes, n_genes_in_vst_rna),
    by = "pathway_id"
  ) %>%
  dplyr::arrange(padj, pval, dplyr::desc(abs(NES)))

# write.csv(
#   fgsea_res,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/GSEA_preranked_DESeq2_stat_targeted_gene_sets.csv",
#   row.names = FALSE
# )

#-----------------------------------------------------------------#
# 12. Over-representation analysis using DE genes
#-----------------------------------------------------------------#

term2gene <- gene_sets_use_long %>%
  dplyr::select(term = pathway_id, gene = Gene) %>%
  dplyr::distinct()

term2name <- gene_sets_use_long %>%
  dplyr::distinct(term = pathway_id, name = gs_name)

bg_genes <- rank_tbl$Gene
up_genes <- deg_rna %>%
  dplyr::filter(Gene %in% bg_genes, !is.na(FDR_DESeq2), FDR_DESeq2 < 0.10, log2FoldChange > 0) %>%
  dplyr::pull(Gene) %>%
  unique()

down_genes <- deg_rna %>%
  dplyr::filter(Gene %in% bg_genes, !is.na(FDR_DESeq2), FDR_DESeq2 < 0.10, log2FoldChange < 0) %>%
  dplyr::pull(Gene) %>%
  unique()

if (length(up_genes) >= 5) {
  ora_up <- clusterProfiler::enricher(
    gene = up_genes,
    universe = bg_genes,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size
  ) %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(DEG_direction = "pCR_high")
} else {
  ora_up <- tibble::tibble()
}

if (length(down_genes) >= 5) {
  ora_down <- clusterProfiler::enricher(
    gene = down_genes,
    universe = bg_genes,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size
  ) %>%
    as.data.frame() %>%
    tibble::as_tibble() %>%
    dplyr::mutate(DEG_direction = "non_pCR_high")
} else {
  ora_down <- tibble::tibble()
}

# Collapse pathway-level annotations to one row per pathway before joining.
# This prevents unintended row expansion when the same pathway_id is mapped to
# multiple Host_axis or Source_type annotations.
gene_set_annot_table <- gene_set_size_table %>%
  dplyr::select(
    pathway_id,
    Host_axis,
    Source_type,
    gs_collection,
    gs_subcollection,
    gs_description,
    gs_exact_source
  ) %>%
  dplyr::group_by(pathway_id) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::everything(),
      ~ paste(unique(stats::na.omit(as.character(.x))), collapse = "; ")
    ),
    .groups = "drop"
  )

ora_res <- dplyr::bind_rows(ora_up, ora_down) %>%
  dplyr::rename(pathway_id = ID, gs_name = Description) %>%
  dplyr::left_join(
    gene_set_annot_table,
    by = "pathway_id",
    relationship = "many-to-one"
  ) %>%
  dplyr::arrange(DEG_direction, p.adjust, pvalue)



# write.csv(
#   ora_res,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/ORA_DEG_targeted_gene_sets.csv",
#   row.names = FALSE
# )

#-----------------------------------------------------------------#
# 13. Sample-level ssGSEA pathway scores
#-----------------------------------------------------------------#

score_rna <- run_ssgsea(vst_rna, gene_sets_use)
score_int <- run_ssgsea(vst_int, gene_sets_use)

score_rna_z <- zscore_rows(score_rna)
score_int_z <- zscore_rows(score_int)

# write.csv(
#   score_rna,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/ssGSEA_scores_all_RNA_samples.csv"
# )
# write.csv(
#   score_rna_z,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/ssGSEA_scores_all_RNA_samples_rowZ.csv"
# )
# write.csv(
#   score_int,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/ssGSEA_scores_matched_metabolite_samples.csv"
# )
# write.csv(
#   score_int_z,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/ssGSEA_scores_matched_metabolite_samples_rowZ.csv"
# )

#-----------------------------------------------------------------#
# 14. pCR vs non-pCR comparison of sample-level pathway scores
#-----------------------------------------------------------------#

score_group_df <- t(score_rna_z) %>%
  as.data.frame(check.names = FALSE) %>%
  tibble::rownames_to_column("RNA_sample_id") %>%
  dplyr::left_join(
    col_rna %>% tibble::rownames_to_column("RNA_sample_id_join") %>% dplyr::select(RNA_sample_id = RNA_sample_id_join, TRG_plot),
    by = "RNA_sample_id"
  ) %>%
  tidyr::pivot_longer(
    cols = -c(RNA_sample_id, TRG_plot),
    names_to = "pathway_id",
    values_to = "ssGSEA_z"
  )

score_group_stats <- score_group_df %>%
  dplyr::group_by(pathway_id) %>%
  dplyr::summarise(
    n_pCR = sum(TRG_plot == "pCR" & is.finite(ssGSEA_z)),
    n_non_pCR = sum(TRG_plot == "non_pCR" & is.finite(ssGSEA_z)),
    median_pCR = median(ssGSEA_z[TRG_plot == "pCR"], na.rm = TRUE),
    median_non_pCR = median(ssGSEA_z[TRG_plot == "non_pCR"], na.rm = TRUE),
    diff_median_pCR_minus_non_pCR = median_pCR - median_non_pCR,
    pval_wilcox = safe_wilcox(ssGSEA_z, TRG_plot),
    .groups = "drop"
  ) %>%
  dplyr::mutate(FDR_wilcox = p.adjust(pval_wilcox, method = "BH")) %>%
  dplyr::left_join(
    gene_set_size_table %>%
      dplyr::select(pathway_id, Host_axis, gs_name, Source_type, gs_collection, gs_subcollection, n_genes_in_vst_rna),
    by = "pathway_id"
  ) %>%
  dplyr::arrange(FDR_wilcox, pval_wilcox, dplyr::desc(abs(diff_median_pCR_minus_non_pCR)))

# write.csv(
#   score_group_stats,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/ssGSEA_pathway_score_pCR_vs_non_pCR_Wilcoxon.csv",
#   row.names = FALSE
# )

#-----------------------------------------------------------------#
# 15. Metabolite-pathway score Spearman correlation
#-----------------------------------------------------------------#

stopifnot(identical(colnames(score_int_z), rownames(met_int)))

score_int_sample <- t(score_int_z) %>%
  as.data.frame(check.names = FALSE) %>%
  tibble::rownames_to_column("RNA_sample_id")

met_sample <- met_int %>%
  as.data.frame(check.names = FALSE) %>%
  tibble::rownames_to_column("RNA_sample_id")

cor_input <- score_int_sample %>%
  dplyr::inner_join(met_sample, by = "RNA_sample_id")

pathway_cols <- colnames(score_int_sample)[colnames(score_int_sample) != "RNA_sample_id"]
met_cols <- colnames(met_sample)[colnames(met_sample) != "RNA_sample_id"]

# Collapse pathway-level annotations to one row per pathway before joining.
# This prevents unintended row expansion in metabolite-pathway correlation results.
pathway_cor_annot_table <- gene_set_size_table %>%
  dplyr::select(
    pathway_id,
    Host_axis,
    gs_name,
    Source_type,
    gs_collection,
    gs_subcollection,
    n_genes_in_vst_int
  ) %>%
  dplyr::group_by(pathway_id) %>%
  dplyr::summarise(
    Host_axis = paste(unique(stats::na.omit(as.character(Host_axis))), collapse = "; "),
    gs_name = paste(unique(stats::na.omit(as.character(gs_name))), collapse = "; "),
    Source_type = paste(unique(stats::na.omit(as.character(Source_type))), collapse = "; "),
    gs_collection = paste(unique(stats::na.omit(as.character(gs_collection))), collapse = "; "),
    gs_subcollection = paste(unique(stats::na.omit(as.character(gs_subcollection))), collapse = "; "),
    n_genes_in_vst_int = max(n_genes_in_vst_int, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    n_genes_in_vst_int = dplyr::if_else(
      is.infinite(n_genes_in_vst_int),
      NA_real_,
      as.numeric(n_genes_in_vst_int)
    )
  )

pathway_met_cor <- tidyr::expand_grid(
  pathway_id = pathway_cols,
  Metabolite = met_cols
) %>%
  dplyr::mutate(
    cor_tbl = purrr::map2(
      pathway_id,
      Metabolite,
      ~ safe_spearman(cor_input[[.x]], cor_input[[.y]])
    )
  ) %>%
  tidyr::unnest(cor_tbl) %>%
  dplyr::mutate(
    FDR_global = p.adjust(pval, method = "BH")
  ) %>%
  dplyr::left_join(
    pathway_cor_annot_table,
    by = "pathway_id",
    relationship = "many-to-one"
  ) %>%
  dplyr::group_by(Metabolite) %>%
  dplyr::mutate(FDR_within_metabolite = p.adjust(pval, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(pathway_id) %>%
  dplyr::mutate(FDR_within_pathway = p.adjust(pval, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(FDR_global, pval, dplyr::desc(abs(rho)))


# write.csv(
#   pathway_met_cor,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/metabolite_pathway_ssGSEA_Spearman_correlation.csv",
#   row.names = FALSE
# )

metabolite_focus_pattern <- "ACETATE|PROPIONATE|BUTYRATE|VALERATE|SCFA|INDOLE|IAA|ILA|IPA|TRYPTOPHAN|KYNURENINE|BILE|CHOLIC|DEOXYCHOLIC|LITHOCHOLIC|URSODEOXYCHOLIC|NIACIN|NICOTINIC|NICOTINAMIDE|NAD"

top_pathway_met_cor <- pathway_met_cor %>%
  dplyr::filter(
    stringr::str_detect(Metabolite, stringr::regex(metabolite_focus_pattern, ignore_case = TRUE)) |
      FDR_global < 0.25 |
      pval < 0.05
  ) %>%
  dplyr::arrange(FDR_global, pval, dplyr::desc(abs(rho)))

# write.csv(
#   top_pathway_met_cor,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/metabolite_pathway_ssGSEA_Spearman_correlation_focus.csv",
#   row.names = FALSE
# )

#-----------------------------------------------------------------#
# 16. Optional gene-level correlation for genes contained in targeted sets
#-----------------------------------------------------------------#

run_gene_metabolite_correlation <- TRUE

if (run_gene_metabolite_correlation) {
  targeted_genes <- unique(gene_sets_use_long$Gene)
  targeted_genes <- targeted_genes[targeted_genes %in% rownames(vst_int)]
  
  gene_expr_int <- t(vst_int[targeted_genes, , drop = FALSE]) %>%
    as.data.frame(check.names = FALSE) %>%
    tibble::rownames_to_column("RNA_sample_id")
  
  gene_cor_input <- gene_expr_int %>%
    dplyr::inner_join(met_sample, by = "RNA_sample_id")
  
  gene_met_cor <- tidyr::expand_grid(
    Gene = targeted_genes,
    Metabolite = met_cols
  ) %>%
    dplyr::mutate(
      cor_tbl = purrr::map2(
        Gene,
        Metabolite,
        ~ safe_spearman(gene_cor_input[[.x]], gene_cor_input[[.y]])
      )
    ) %>%
    tidyr::unnest(cor_tbl) %>%
    dplyr::mutate(FDR_global = p.adjust(pval, method = "BH")) %>%
    dplyr::left_join(
      gene_sets_use_long %>%
        dplyr::select(Host_axis, gs_name, Gene, Source_type) %>%
        dplyr::distinct(),
      by = "Gene"
    ) %>%
    dplyr::arrange(FDR_global, pval, dplyr::desc(abs(rho)))
  
  # write.csv(
  #   gene_met_cor,
  #   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/metabolite_targeted_gene_Spearman_correlation.csv",
  #   row.names = FALSE
  # )
}

gene_sets_use_long %>%
  dplyr::select(Host_axis, gs_name, Gene, Source_type) %>%
  dplyr::distinct() %>%
  dplyr::count(Gene, sort = TRUE) %>%
  dplyr::filter(n > 1)

#-----------------------------------------------------------------#
# 17. Figures
#-----------------------------------------------------------------#

group_cols <- c(
  pCR = "#4FAE9A",
  non_pCR = "#DE7872"
)

#-----------------------------------------------------------------#
# 17-1. Representative GSEA dot plot with curated family catalog
#-----------------------------------------------------------------#

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork", type = "binary")
}

#----------------------------#
# User-adjustable filters
#----------------------------#

fdr_cutoff <- 0.10
nes_cutoff <- 1.00

# Set TRUE only if you want to show key metabolite axes even when FDR is modest.
include_relaxed_metabolite <- FALSE
relaxed_metabolite_fdr_cutoff <- 0.25
relaxed_metabolite_nes_cutoff <- 0.80

#----------------------------#
# Label helper
#----------------------------#

format_gsea_label <- function(x) {
  x %>%
    stringr::str_remove("^HALLMARK_") %>%
    stringr::str_remove("^REACTOME_") %>%
    stringr::str_remove("^GOBP_") %>%
    stringr::str_remove("^GOCC_") %>%
    stringr::str_remove("^GOMF_") %>%
    stringr::str_remove("^WP_") %>%
    stringr::str_remove("^KEGG_") %>%
    stringr::str_remove("_WP[0-9]+$") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_lower() %>%
    stringr::str_to_sentence() %>%
    stringr::str_replace_all("\\bCrc\\b", "CRC") %>%
    stringr::str_replace_all("\\bScfa\\b", "SCFA") %>%
    stringr::str_replace_all("\\bAhr\\b", "AhR") %>%
    stringr::str_replace_all("\\bNfkb\\b", "NF-kB") %>%
    stringr::str_replace_all("\\bTnf\\b", "TNF") %>%
    stringr::str_replace_all("\\bIfn\\b", "IFN") %>%
    stringr::str_replace_all("\\bMyc\\b", "MYC") %>%
    stringr::str_replace_all("\\bEmt\\b", "EMT") %>%
    stringr::str_replace_all("\\bG2m\\b", "G2/M") %>%
    stringr::str_squish()
}

#----------------------------#
# Category colors
#----------------------------#

axis_cols <- c(
  "Tumor remodeling/CRC" = "#F08A84",   # coral/red
  "Metabolite signaling" = "#A8D39A",   # green
  "Innate immunity" = "#9FCBE3",        # blue
  "Adaptive immunity" = "#C3AFE8",      # purple
  "Barrier/mucus/AMP" = "#EBCB74",      # yellow
  "Gut homing" = "#83C9C1"              # teal
)

#----------------------------#
# Curated family catalog
#----------------------------#

term_catalog <- tibble::tribble(
  ~main_category, ~family_id, ~display_label, ~regex, ~preferred_terms, ~quota_order,
  
  # Metabolite signaling
  "Metabolite signaling", "MET_SCFA", "Butyrate/SCFA response",
  "BUTYRATE|SHORT_CHAIN|SCFA",
  "GOBP_CELLULAR_RESPONSE_TO_BUTYRATE;GOBP_SHORT_CHAIN_FATTY_ACID_METABOLIC_PROCESS",
  1,
  
  "Metabolite signaling", "MET_FA", "Fatty acid metabolism",
  "FATTY_ACID_METABOLISM|CELLULAR_RESPONSE_TO_FATTY_ACID",
  "HALLMARK_FATTY_ACID_METABOLISM;GOBP_CELLULAR_RESPONSE_TO_FATTY_ACID",
  2,
  
  "Metabolite signaling", "MET_AHR", "AhR/indole ligand response",
  "ARYL_HYDROCARBON|AHR|INDOLE|TRYPTOPHAN_DERIVATIVE|KYNURENINE",
  "REACTOME_ARYL_HYDROCARBON_RECEPTOR_SIGNALLING;WP_ARYL_HYDROCARBON_RECEPTOR_PATHWAY_WP2586",
  3,
  
  "Metabolite signaling", "MET_BILE", "Bile acid response/metabolism",
  "BILE_ACID|BILE_SALT|FXR|TGR5|GPBAR|NR1H4|GPBAR1",
  "GOBP_CELLULAR_RESPONSE_TO_BILE_ACID;HALLMARK_BILE_ACID_METABOLISM;REACTOME_BILE_ACID_AND_BILE_SALT_METABOLISM",
  4,
  
  "Metabolite signaling", "MET_NIACIN", "NAD/niacin metabolism",
  "NIACIN|NICOTINIC|NAD",
  "GOMF_NICOTINIC_ACID_RECEPTOR_ACTIVITY;GOBP_NAD_METABOLIC_PROCESS;GOBP_NAD_BIOSYNTHETIC_PROCESS",
  5,
  
  "Metabolite signaling", "MET_ROS", "Oxidative stress/ROS response",
  "REACTIVE_OXYGEN|OXIDATIVE_STRESS|RESPONSE_TO_OXIDATIVE|ROS|NRF2|NFE2L2|GLUTATHIONE|PEROXIDASE|SUPEROXIDE_DISMUTASE",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY;GOBP_RESPONSE_TO_OXIDATIVE_STRESS;GOBP_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS",
  6,
  
  "Metabolite signaling", "MET_CARB_DIGEST", "Carbohydrate digestion/absorption",
  "CARBOHYDRATE_DIGESTION|STARCH|SUCROSE|DISACCHARIDE|POLYSACCHARIDE|GLYCAN",
  "KEGG_CARBOHYDRATE_DIGESTION_AND_ABSORPTION;KEGG_STARCH_AND_SUCROSE_METABOLISM;REACTOME_DIGESTION_AND_ABSORPTION",
  7,
  
  "Metabolite signaling", "MET_HISTIDINE", "Histidine catabolism",
  "HISTIDINE|UROCANATE|HISTAMINE|IMIDAZOLE",
  "REACTOME_HISTIDINE_CATABOLISM;GOBP_HISTIDINE_CATABOLIC_PROCESS;GOBP_L_HISTIDINE_CATABOLIC_PROCESS;KEGG_HISTIDINE_METABOLISM",
  8,
  
  "Metabolite signaling", "MET_BCAA", "Branched-chain amino acid catabolism",
  "BRANCHED_CHAIN_AMINO_ACID|BRANCHED-CHAIN|BCAA|LEUCINE|ISOLEUCINE|VALINE",
  "REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM;GOBP_BRANCHED_CHAIN_AMINO_ACID_CATABOLIC_PROCESS;GOBP_LEUCINE_CATABOLIC_PROCESS;GOBP_VALINE_CATABOLIC_PROCESS;GOBP_ISOLEUCINE_CATABOLIC_PROCESS",
  9,
  
  "Metabolite signaling", "MET_PROTEIN_DIGEST", "Protein digestion/absorption",
  "PROTEIN_DIGESTION|DIGESTION_AND_ABSORPTION|PEPTIDE_TRANSPORT|AMINO_ACID_ABSORPTION",
  "KEGG_PROTEIN_DIGESTION_AND_ABSORPTION;REACTOME_DIGESTION_AND_ABSORPTION",
  10,
  
  "Metabolite signaling", "MET_AA_BROAD", "Amino acid transport/metabolism",
  "AMINO_ACID_TRANSPORT|AMINO_ACID_TRANSMEMBRANE|AMINO_ACID_METABOLIC_PROCESS|AMINO_ACID_CATABOLIC_PROCESS",
  "GOBP_AMINO_ACID_TRANSPORT;GOBP_AMINO_ACID_TRANSMEMBRANE_TRANSPORT;GOBP_CELLULAR_AMINO_ACID_METABOLIC_PROCESS;GOBP_AMINO_ACID_CATABOLIC_PROCESS",
  11,
  
  "Metabolite signaling", "MET_HEME", "Heme metabolism",
  "HEME_METABOLISM|HEME_METABOLIC_PROCESS|PORPHYRIN",
  "HALLMARK_HEME_METABOLISM;GOBP_HEME_METABOLIC_PROCESS;GOBP_PORPHYRIN_CONTAINING_COMPOUND_METABOLIC_PROCESS",
  12,
  
  
  # Innate immunity
  "Innate immunity", "INN_TLR_NOD", "TLR/NOD pattern-recognition signaling",
  "TOLL_LIKE|TLR|NOD1|NOD2",
  "REACTOME_TOLL_LIKE_RECEPTOR_CASCADES;REACTOME_NOD1_2_SIGNALING_PATHWAY",
  1,
  
  "Innate immunity", "INN_IFN", "Interferon response",
  "INTERFERON|IFN",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE;HALLMARK_INTERFERON_ALPHA_RESPONSE",
  2,
  
  "Innate immunity", "INN_TNF_IL6", "TNF-NF-kB / IL6-JAK-STAT signaling",
  "TNFA_SIGNALING_VIA_NFKB|IL6_JAK_STAT3_SIGNALING|TNFA|TNF|NF_KB|NFKB|IL6|JAK_STAT3",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB;HALLMARK_IL6_JAK_STAT3_SIGNALING",
  3,
  
  "Innate immunity", "INN_COMPLEMENT", "Complement",
  "COMPLEMENT",
  "HALLMARK_COMPLEMENT",
  4,
  
  # Adaptive immunity
  "Adaptive immunity", "ADP_TCELL_DIFF", "T-cell differentiation",
  "T_CELL_DIFFERENTIATION|T_HELPER|TH1|TH2|TH17|REGULATORY_T_CELL|TREG",
  "GOBP_T_CELL_DIFFERENTIATION;GOBP_T_HELPER_17_CELL_DIFFERENTIATION;GOBP_REGULATORY_T_CELL_DIFFERENTIATION",
  1,
  
  "Adaptive immunity", "ADP_LYMPH", "Lymphocyte-mediated immunity",
  "LYMPHOCYTE_MEDIATED_IMMUNITY",
  "GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY",
  2,
  
  "Adaptive immunity", "ADP_TACT", "T-cell activation",
  "ALPHA_BETA_T_CELL|T_CELL_ACTIVATION",
  "GOBP_ALPHA_BETA_T_CELL_ACTIVATION;GOBP_T_CELL_ACTIVATION",
  3,
  
  # Gut homing
  "Gut homing", "GUT_INTEGRIN", "Integrin-mediated gut homing",
  "INTEGRIN|ITGA4|ITGB7|CCR9|GUT_HOMING",
  "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
  1,
  
  "Gut homing", "GUT_MIGRATION", "Inflammatory leukocyte trafficking",
  "LEUKOCYTE_MIGRATION",
  "GOBP_LEUKOCYTE_MIGRATION_INVOLVED_IN_INFLAMMATORY_RESPONSE;GOBP_LEUKOCYTE_MIGRATION",
  2,
  
  # Barrier / mucus / AMP
  "Barrier/mucus/AMP", "BAR_JUNCTION", "Apical/tight junction",
  "APICAL_JUNCTION|APICAL_SURFACE|TIGHT_JUNCTION|CELL_CELL_JUNCTION",
  "HALLMARK_APICAL_JUNCTION;HALLMARK_APICAL_SURFACE;GOBP_TIGHT_JUNCTION_ORGANIZATION",
  1,
  
  "Barrier/mucus/AMP", "BAR_MUCUS", "Mucus/goblet-cell program",
  "MUCUS|MUCIN|GOBLET",
  "GOBP_REGULATION_OF_MUCUS_SECRETION;GOBP_MUCUS_SECRETION;GOBP_GOBLET_CELL_DIFFERENTIATION",
  2,
  
  "Barrier/mucus/AMP", "BAR_AMP", "Antimicrobial defense program",
  "ANTIMICROBIAL|DEFENSE_RESPONSE_TO_BACTERIUM|DEFENSIN",
  "GOBP_ANTIMICROBIAL_HUMORAL_RESPONSE;GOBP_DEFENSE_RESPONSE_TO_BACTERIUM",
  3,
  
  # Tumor remodeling / CRC
  "Tumor remodeling/CRC", "TUM_CRC_ADENOMA", "Colorectal adenoma signature",
  "ADENOMA",
  "SABATES_COLORECTAL_ADENOMA_UP",
  1,
  
  "Tumor remodeling/CRC", "TUM_EMT", "Epithelial-mesenchymal transition",
  "EPITHELIAL_MESENCHYMAL|EMT",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION;GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION",
  2,
  
  "Tumor remodeling/CRC", "TUM_ECM", "ECM/collagen remodeling",
  "EXTRACELLULAR_MATRIX|COLLAGEN",
  "GOBP_EXTRACELLULAR_MATRIX_ORGANIZATION;GOBP_COLLAGEN_FIBRIL_ORGANIZATION",
  3,
  
  "Tumor remodeling/CRC", "TUM_ANGIO", "Angiogenic remodeling",
  "ANGIOGENESIS",
  "HALLMARK_ANGIOGENESIS",
  4,
  
  "Tumor remodeling/CRC", "TUM_HYPOXIA", "Hypoxia response",
  "HYPOXIA",
  "HALLMARK_HYPOXIA",
  5,
  
  "Tumor remodeling/CRC", "TUM_WNT", "Wnt/beta-catenin signaling",
  "WNT|BETA_CATENIN",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  6,
  
  "Tumor remodeling/CRC", "TUM_MYC", "MYC targets",
  "MYC",
  "HALLMARK_MYC_TARGETS_V1",
  7,
  
  "Tumor remodeling/CRC", "TUM_G2M", "G2/M checkpoint",
  "G2M|G2_M|G2/M",
  "HALLMARK_G2M_CHECKPOINT",
  8,
  
  "Tumor remodeling/CRC", "TUM_REMODEL", "Tissue remodeling/migration",
  "WOUND_HEALING|TISSUE_REMODELING|CELL_MIGRATION|EPITHELIAL_CELL_MIGRATION",
  "GOBP_EPITHELIAL_CELL_MIGRATION;GOBP_CELL_MIGRATION;GOBP_WOUND_HEALING",
  9
)

category_quota <- tibble::tribble(
  ~main_category, ~category_quota,
  "Tumor remodeling/CRC", 9,
  "Metabolite signaling", 11,
  "Innate immunity", 4,
  "Adaptive immunity", 3,
  "Barrier/mucus/AMP", 3,
  "Gut homing", 2
)

#----------------------------#
# Pre-filter candidate terms
#----------------------------#

gsea_candidates <- fgsea_res %>%
  dplyr::filter(!is.na(padj), !is.na(NES)) %>%
  dplyr::mutate(
    gs_name_upper = toupper(gs_name),
    pathway_id_upper = toupper(pathway_id),
    gs_description_upper = toupper(dplyr::coalesce(gs_description, "")),
    Host_axis_upper = toupper(dplyr::coalesce(Host_axis, "")),
    term_text_upper = paste(gs_name_upper, pathway_id_upper, gs_description_upper, Host_axis_upper, sep = " | "),
    gs_collection_upper = toupper(dplyr::coalesce(gs_collection, "")),
    gs_subcollection_upper = toupper(dplyr::coalesce(gs_subcollection, "")),
    Source_type_upper = toupper(dplyr::coalesce(Source_type, "")),
    Direction = dplyr::case_when(
      NES > 0 ~ "pCR",
      NES < 0 ~ "non_pCR",
      TRUE ~ NA_character_
    ),
    neglog10_FDR = -log10(padj + 1e-300),
    neglog10_FDR_capped = pmin(neglog10_FDR, 10),
    FDR_symbol = dplyr::case_when(
      neglog10_FDR >= 10 ~ "+",
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ ""
    ),
    FDR_text_color = dplyr::case_when(
      neglog10_FDR_capped >= 5 ~ "white",
      TRUE ~ "grey20"
    ),
    is_allowed_source = dplyr::case_when(
      gs_collection_upper == "H" ~ TRUE,
      gs_collection_upper == "C5" ~ TRUE,
      gs_collection_upper == "C2" & gs_subcollection_upper != "CGP" ~ TRUE,
      Source_type_upper %in% c(
        "DIRECT_GO_ORGHS",
        "PREFERRED_MSIGDB",
        "KEYWORD_MSIGDB_RESCUE"
      ) ~ TRUE,
      TRUE ~ FALSE
    ),
    is_excluded_term =
      stringr::str_detect(term_text_upper, "^HP_|\\| HP_") |
      stringr::str_detect(term_text_upper, "OVARIAN|ENDOMETRIAL|THYMUS|BRAIN|NEURON|CARDIAC|MUSCLE|KIDNEY|LIVER|HEPATOCYTE|HEPATOCYTES|BRONCHIAL|LUNG|BREAST|PROSTATE|MELANOMA") |
      stringr::str_detect(term_text_upper, "NYSTAGMUS|HEAD_NODDING|IMMUNODEFICIENCY|POSTCOVID") |
      stringr::str_detect(term_text_upper, "CYTOKINESIS|CYTOKINETIC|MIDBODY|SPINDLE|CENTROSOME|CENTROMERE") |
      stringr::str_detect(term_text_upper, "NADPH_OXIDASE|ALCOHOL_DEHYDROGENASE|SUPEROXIDE_GENERATING_NADPH_OXIDASE") |
      stringr::str_detect(term_text_upper, "MANUMYCIN|CRPC|PID_|BILANGES") |
      stringr::str_detect(term_text_upper, "TRANSITION METAL|METAL ION TRANSPORT|TRANSMEMBRANE TRANSPORTER ACTIVITY") |
      stringr::str_detect(term_text_upper, "SPHINGOMYELIN|GLUCOSYLCERAMIDE") |
      stringr::str_detect(gs_name_upper, "_DN$") |
      stringr::str_detect(gs_name_upper, "MYC_TARGETS_V2"),
    pass_strict = padj < fdr_cutoff & abs(NES) >= nes_cutoff,
    pass_relaxed_metabolite = FALSE
  ) %>%
  dplyr::filter(
    !is.na(Direction),
    is_allowed_source,
    !is_excluded_term
  )

# Optional relaxed metabolite display
if (include_relaxed_metabolite) {
  gsea_candidates <- gsea_candidates %>%
    dplyr::mutate(
      pass_relaxed_metabolite =
        padj < relaxed_metabolite_fdr_cutoff &
        abs(NES) >= relaxed_metabolite_nes_cutoff &
        stringr::str_detect(
          term_text_upper,
          "BUTYRATE|SHORT_CHAIN|SCFA|FATTY_ACID|ARYL_HYDROCARBON|AHR|XENOBIOTIC|INDOLE|BILE_ACID|BILE_SALT|NIACIN|NICOTINIC|NAD"
        )
    )
}

#----------------------------#
# Match candidates to curated families
#----------------------------#

gsea_family_candidates <- term_catalog %>%
  dplyr::mutate(join_key = 1) %>%
  dplyr::inner_join(
    gsea_candidates %>% dplyr::mutate(join_key = 1),
    by = "join_key",
    relationship = "many-to-many"
  ) %>%
  dplyr::filter(
    stringr::str_detect(term_text_upper, regex) |
      purrr::map2_lgl(
        preferred_terms,
        gs_name,
        ~ .y %in% unlist(strsplit(.x, ";", fixed = TRUE))
      )
  ) %>%
  dplyr::mutate(
    preferred_vec = strsplit(preferred_terms, ";", fixed = TRUE),
    preferred_rank = purrr::map2_int(
      preferred_vec,
      gs_name,
      ~ {
        out <- match(.y, .x)
        ifelse(is.na(out), 999L, as.integer(out))
      }
    ),
    is_preferred = preferred_rank < 999,
    pass_display_filter = pass_strict |
      (main_category == "Metabolite signaling" & pass_relaxed_metabolite)
  ) %>%
  dplyr::filter(pass_display_filter)

# One representative term per family
gsea_family_representatives <- gsea_family_candidates %>%
  dplyr::group_by(family_id) %>%
  dplyr::arrange(
    dplyr::desc(is_preferred),
    preferred_rank,
    padj,
    dplyr::desc(abs(NES)),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup()

# Category quota
gsea_plot_df <- gsea_family_representatives %>%
  dplyr::left_join(
    category_quota %>%
      dplyr::rename(category_n = category_quota),
    by = "main_category",
    relationship = "many-to-one"
  ) %>%
  dplyr::group_by(main_category) %>%
  dplyr::arrange(
    quota_order,
    padj,
    dplyr::desc(abs(NES)),
    .by_group = TRUE
  ) %>%
  dplyr::mutate(
    category_rank = dplyr::row_number()
  ) %>%
  dplyr::filter(
    category_rank <= category_n
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(dplyr::desc(NES)) %>%
  dplyr::mutate(
    Direction = factor(Direction, levels = c("pCR", "non_pCR")),
    main_category = factor(
      main_category,
      levels = c(
        "Tumor remodeling/CRC",
        "Metabolite signaling",
        "Barrier/mucus/AMP",
        "Innate immunity",
        "Adaptive immunity",
        "Gut homing"
      )
    ),
    pathway_label = factor(display_label, levels = rev(display_label))
  )

message("Number of pathways plotted: ", nrow(gsea_plot_df))
print(gsea_plot_df %>% dplyr::count(main_category, Direction))

write.csv(
  gsea_plot_df %>%
    dplyr::select(
      main_category,
      family_id,
      display_label,
      pathway_id,
      gs_name,
      NES,
      pval,
      padj,
      neglog10_FDR,
      size,
      pass_strict,
      pass_relaxed_metabolite,
      Source_type,
      gs_collection,
      gs_subcollection,
      gs_description
    ),
  "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/GSEA_dotplot_representative_terms.csv",
  row.names = FALSE
)

#----------------------------#
# Diagnostic file: all metabolite candidates
#----------------------------#

metabolite_term_check <- gsea_family_candidates %>%
  dplyr::filter(main_category == "Metabolite signaling") %>%
  dplyr::arrange(family_id, padj, dplyr::desc(abs(NES))) %>%
  dplyr::select(
    family_id,
    display_label,
    gs_name,
    pathway_id,
    NES,
    pval,
    padj,
    neglog10_FDR,
    pass_strict,
    pass_relaxed_metabolite,
    size,
    gs_description
  )

# write.csv(
#   metabolite_term_check,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/GSEA_metabolite_term_check.csv",
#   row.names = FALSE
# )

#----------------------------#
# Plot panels
#----------------------------#

p_gsea_label <- gsea_plot_df %>%
  ggplot2::ggplot(
    ggplot2::aes(x = 1, y = pathway_label, label = pathway_label)
  ) +
  ggplot2::geom_text(hjust = 1, size = 3.1, color = "grey15") +
  ggplot2::scale_y_discrete(limits = levels(gsea_plot_df$pathway_label), drop = FALSE) +
  ggplot2::scale_x_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
  ggplot2::theme_void() +
  ggplot2::theme(
    plot.margin = ggplot2::margin(5.5, 1, 5.5, 5.5)
  )

p_axis_tile <- gsea_plot_df %>%
  ggplot2::ggplot(
    ggplot2::aes(x = 1, y = pathway_label, fill = main_category)
  ) +
  ggplot2::geom_tile(width = 1, height = 1, color = "white", linewidth = 0.15) +
  ggplot2::scale_y_discrete(limits = levels(gsea_plot_df$pathway_label), drop = FALSE, expand = c(0, 0)) +
  ggplot2::scale_x_continuous(breaks = 1, labels = "Axis", position = "top", expand = c(0, 0)) +
  ggplot2::scale_fill_manual(values = axis_cols, name = "Category") +
  ggplot2::theme_void() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 8, color = "grey20"),
    plot.margin = ggplot2::margin(5.5, 0, 5.5, 0)
  )

p_fdr_tile <- gsea_plot_df %>%
  ggplot2::ggplot(
    ggplot2::aes(x = 1, y = pathway_label, fill = neglog10_FDR_capped)
  ) +
  ggplot2::geom_tile(width = 1, height = 1, color = "white", linewidth = 0.15) +
  ggplot2::geom_text(
    ggplot2::aes(label = FDR_symbol, color = FDR_text_color),
    size = 2.6,
    fontface = "bold"
  ) +
  ggplot2::scale_color_identity() +
  ggplot2::scale_y_discrete(limits = levels(gsea_plot_df$pathway_label), drop = FALSE, expand = c(0, 0)) +
  ggplot2::scale_x_continuous(breaks = 1, labels = "FDR", position = "top", expand = c(0, 0)) +
  ggplot2::scale_fill_gradient(
    low = "white",
    high = "#08519C",
    limits = c(0, 10),
    oob = scales::squish,
    name = "-log10(FDR)\ncap = 10"
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 8, color = "grey20"),
    plot.margin = ggplot2::margin(5.5, 2, 5.5, 0)
  )

p_gsea_main <- gsea_plot_df %>%
  ggplot2::ggplot(
    ggplot2::aes(y = pathway_label, color = Direction)
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.3,
    color = "grey45"
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend = NES,
      yend = pathway_label
    ),
    linewidth = 0.55,
    alpha = 0.75
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = NES),
    size = 3.7,
    alpha = 0.95
  ) +
  ggplot2::scale_y_discrete(
    limits = levels(gsea_plot_df$pathway_label),
    drop = FALSE
  ) +
  ggplot2::scale_x_continuous(
    breaks = scales::pretty_breaks(n = 4),
    expand = ggplot2::expansion(mult = c(0.08, 0.08))
  ) +
  ggplot2::scale_color_manual(
    values = group_cols,
    breaks = c("pCR", "non_pCR"),
    labels = c("pCR-high", "non-pCR-high"),
    name = "Enriched in"
  ) +
  ggplot2::labs(
    x = "NES",
    y = NULL
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 8.5),
    axis.title.x = ggplot2::element_text(size = 9.5),
    legend.title = ggplot2::element_text(size = 9),
    legend.text = ggplot2::element_text(size = 8),
    plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 0)
  )

p_gsea_dot <- p_gsea_label +
  p_axis_tile +
  p_fdr_tile +
  p_gsea_main +
  patchwork::plot_layout(
    widths = c(4.8, 0.28, 0.32, 2.9),
    guides = "collect"
  ) &
  ggplot2::theme(legend.position = "right")

p_gsea_dot

ggplot2::ggsave(
  "figures/host_RNAseq_targeted_pathway_GSEA_lollipop_representative_catalog.svg",
  p_gsea_dot,
  width = 6,
  height = min(11, max(6.5, 0.30 * nrow(gsea_plot_df) + 1.2)-3),
  device = "svg"
)




# metabolic_candidate_check <- fgsea_res %>%
#   dplyr::mutate(
#     gs_name_upper = toupper(gs_name),
#     pathway_id_upper = toupper(pathway_id),
#     gs_description_upper = toupper(dplyr::coalesce(gs_description, "")),
#     Host_axis_upper = toupper(dplyr::coalesce(Host_axis, "")),
#     term_text_upper = paste(
#       gs_name_upper,
#       pathway_id_upper,
#       gs_description_upper,
#       Host_axis_upper,
#       sep = " | "
#     ),
#     candidate_axis = dplyr::case_when(
#       stringr::str_detect(term_text_upper, "HISTIDINE|UROCANATE|HISTAMINE|IMIDAZOLE") ~
#         "Histidine-related",
#       
#       stringr::str_detect(term_text_upper, "BRANCHED_CHAIN_AMINO_ACID|BRANCHED-CHAIN|BCAA|LEUCINE|ISOLEUCINE|VALINE") ~
#         "BCAA-related",
#       
#       stringr::str_detect(term_text_upper, "PROTEIN_DIGESTION|DIGESTION_AND_ABSORPTION|PEPTIDE_TRANSPORT|AMINO_ACID_ABSORPTION") ~
#         "Protein digestion/absorption",
#       
#       stringr::str_detect(term_text_upper, "AMINO_ACID_TRANSPORT|AMINO_ACID_TRANSMEMBRANE|AMINO_ACID_METABOLIC_PROCESS|AMINO_ACID_CATABOLIC_PROCESS") ~
#         "Broad amino acid metabolism",
#       
#       stringr::str_detect(term_text_upper, "HEME|PORPHYRIN|NITROSATIVE|NITRIC_OXIDE") ~
#         "Processed-meat proxy",
#       
#       stringr::str_detect(term_text_upper, "STARCH|SUCROSE|DISACCHARIDE|POLYSACCHARIDE|CARBOHYDRATE_DIGESTION|CARBOHYDRATE_ABSORPTION|GLUCAN") ~
#         "Carbohydrate-related",
#       
#       TRUE ~ NA_character_
#     )
#   ) %>%
#   dplyr::filter(!is.na(candidate_axis)) %>%
#   dplyr::select(
#     candidate_axis,
#     pathway_id,
#     gs_name,
#     NES,
#     pval,
#     padj,
#     size,
#     gs_collection,
#     gs_subcollection,
#     Source_type,
#     gs_description
#   ) %>%
#   dplyr::arrange(candidate_axis, padj, dplyr::desc(abs(NES)))
# 
# metabolic_candidate_check %>% 
#   as.data.frame() %>% 
#   pull(pathway_id) %>% 
#   unique()

# # Angiogenesis vs. hypoxia
# angi_hypoxia_genes <- gene_sets_use_long %>%
#   dplyr::filter(
#     gs_name %in% c(
#       "HALLMARK_ANGIOGENESIS",
#       "HALLMARK_HYPOXIA"
#     )
#   ) %>%
#   dplyr::select(gs_name, Gene) %>%
#   dplyr::distinct()
# 
# angi_genes <- angi_hypoxia_genes %>%
#   dplyr::filter(gs_name == "HALLMARK_ANGIOGENESIS") %>%
#   dplyr::pull(Gene) %>%
#   unique()
# 
# hypoxia_genes <- angi_hypoxia_genes %>%
#   dplyr::filter(gs_name == "HALLMARK_HYPOXIA") %>%
#   dplyr::pull(Gene) %>%
#   unique()
# 
# angi_hypoxia_overlap <- tibble::tibble(
#   set1 = "HALLMARK_ANGIOGENESIS",
#   set2 = "HALLMARK_HYPOXIA",
#   n_angiogenesis = length(angi_genes),
#   n_hypoxia = length(hypoxia_genes),
#   n_overlap = length(intersect(angi_genes, hypoxia_genes)),
#   jaccard = length(intersect(angi_genes, hypoxia_genes)) /
#     length(union(angi_genes, hypoxia_genes)),
#   overlap_genes = paste(intersect(angi_genes, hypoxia_genes), collapse = ";")
# )
# 
# angi_hypoxia_overlap
# # # A tibble: 1 × 7 --> 거의 겹치지 않음
# # set1              set2  n_angiogenesis n_hypoxia n_overlap jaccard overlap_genes
# # <chr>             <chr>          <int>     <int>     <int>   <dbl> <chr>        
# #   1 HALLMARK_ANGIOGE… HALL…             36       193         3  0.0133
# 
# fgsea_res %>%
#   dplyr::filter(
#     gs_name %in% c(
#       "HALLMARK_ANGIOGENESIS",
#       "HALLMARK_HYPOXIA"
#     )
#   ) %>%
#   dplyr::mutate(
#     leading_edge_genes = purrr::map_chr(
#       leadingEdge,
#       ~ paste(.x, collapse = ";")
#     )
#   ) %>%
#   dplyr::select(
#     gs_name,
#     NES,
#     pval,
#     padj,
#     leading_edge_genes
#   )

# # overlap / leading-edge check 
# msigdbr(species = "Homo sapiens") %>%
#   dplyr::filter(
#     gs_name %in% c(
#       "REACTOME_TOLL_LIKE_RECEPTOR_CASCADES",
#       "REACTOME_NOD1_2_SIGNALING_PATHWAY",
#       "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
#     )
#   ) %>%
#   dplyr::select(gs_name, gene_symbol) %>%
#   dplyr::distinct() -> t_nfkb_members
# 
# t_nfkb_members %>%
#   dplyr::count(gs_name)
# 
# # Jaccard overlap
# a <- t_nfkb_members %>%
#   dplyr::filter(gs_name == "REACTOME_TOLL_LIKE_RECEPTOR_CASCADES") %>%
#   dplyr::pull(gene_symbol) %>%
#   unique()
# 
# b <- t_nfkb_members %>%
#   dplyr::filter(gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB") %>%
#   dplyr::pull(gene_symbol) %>%
#   unique()
# 
# length(intersect(a, b)) / length(union(a, b))
# # [1] 0.04225352 # 낮음
# 
# # leading-edge check
# fgsea_res %>%
#   dplyr::filter(
#     gs_name %in% c(
#       "REACTOME_TOLL_LIKE_RECEPTOR_CASCADES",
#       "REACTOME_NOD1_2_SIGNALING_PATHWAY",
#       "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
#     )
#   ) %>%
#   dplyr::mutate(
#     leading_edge_genes = purrr::map_chr(leadingEdge, ~ paste(.x, collapse = ", "))
#   ) %>%
#   dplyr::select(gs_name, NES, padj, leading_edge_genes)

#---------------------------------------------------------------#
# Check overlap among innate immune gene-set families
#---------------------------------------------------------------#
# 
# innate_family_ids <- c(
#   "INN_TLR_NOD",
#   "INN_ACTIVATION",
#   "INN_REGULATION",
#   "INN_INFLAM",
#   "INN_IFN",
#   "INN_COMPLEMENT"
# )
# 
# innate_representatives <- gsea_family_representatives %>%
#   dplyr::filter(family_id %in% innate_family_ids) %>%
#   dplyr::distinct(
#     family_id,
#     display_label,
#     pathway_id,
#     gs_name,
#     NES,
#     padj,
#     .keep_all = TRUE
#   ) %>%
#   dplyr::arrange(family_id)
# 
# innate_representatives
# 
# innate_members <- gene_sets_use_long %>%
#   dplyr::semi_join(
#     innate_representatives %>%
#       dplyr::select(pathway_id, gs_name),
#     by = c("pathway_id", "gs_name")
#   ) %>%
#   dplyr::distinct(pathway_id, gs_name, Gene)
# 
# pairwise_jaccard <- function(df, id_col = "pathway_id", gene_col = "Gene") {
#   ids <- unique(df[[id_col]])
#   
#   purrr::map_dfr(
#     utils::combn(ids, 2, simplify = FALSE),
#     function(x) {
#       g1 <- df %>%
#         dplyr::filter(.data[[id_col]] == x[1]) %>%
#         dplyr::pull(.data[[gene_col]]) %>%
#         unique()
#       
#       g2 <- df %>%
#         dplyr::filter(.data[[id_col]] == x[2]) %>%
#         dplyr::pull(.data[[gene_col]]) %>%
#         unique()
#       
#       tibble::tibble(
#         pathway_id_1 = x[1],
#         pathway_id_2 = x[2],
#         n_1 = length(g1),
#         n_2 = length(g2),
#         n_overlap = length(intersect(g1, g2)),
#         jaccard = length(intersect(g1, g2)) / length(union(g1, g2)),
#         overlap_genes = paste(intersect(g1, g2), collapse = ";")
#       )
#     }
#   )
# }
# 
# innate_overlap <- pairwise_jaccard(innate_members) %>%
#   dplyr::left_join(
#     innate_representatives %>%
#       dplyr::select(
#         pathway_id_1 = pathway_id,
#         family_id_1 = family_id,
#         label_1 = display_label,
#         gs_name_1 = gs_name,
#         NES_1 = NES,
#         padj_1 = padj
#       ),
#     by = "pathway_id_1"
#   ) %>%
#   dplyr::left_join(
#     innate_representatives %>%
#       dplyr::select(
#         pathway_id_2 = pathway_id,
#         family_id_2 = family_id,
#         label_2 = display_label,
#         gs_name_2 = gs_name,
#         NES_2 = NES,
#         padj_2 = padj
#       ),
#     by = "pathway_id_2"
#   ) %>%
#   dplyr::select(
#     family_id_1,
#     label_1,
#     gs_name_1,
#     NES_1,
#     padj_1,
#     family_id_2,
#     label_2,
#     gs_name_2,
#     NES_2,
#     padj_2,
#     n_1,
#     n_2,
#     n_overlap,
#     jaccard,
#     overlap_genes
#   ) %>%
#   dplyr::arrange(dplyr::desc(jaccard), dplyr::desc(n_overlap))
# 
# innate_overlap %>% as.data.frame()
# 
# innate_leading_edge <- fgsea_res %>%
#   dplyr::semi_join(
#     innate_representatives %>%
#       dplyr::select(pathway_id, gs_name),
#     by = c("pathway_id", "gs_name")
#   ) %>%
#   dplyr::mutate(
#     leading_edge_genes = purrr::map_chr(
#       leadingEdge,
#       ~ paste(.x, collapse = ";")
#     )
#   ) %>%
#   dplyr::select(
#     pathway_id,
#     gs_name,
#     NES,
#     pval,
#     padj,
#     leading_edge_genes
#   ) %>%
#   dplyr::arrange(dplyr::desc(NES))
# 
# innate_leading_edge


#-----------------------------------------------------------------#
# 17-2. Pathway-score pCR vs non-pCR plot for curated significant pathways
#-----------------------------------------------------------------#

stopifnot(exists("gsea_plot_df"))

target_pathway_annot <- gsea_plot_df %>%
  dplyr::mutate(
    gsea_rank = dplyr::row_number()
  ) %>%
  dplyr::select(
    pathway_id,
    family_id,
    display_label,
    main_category,
    Direction,
    NES,
    padj,
    gsea_rank
  ) %>%
  dplyr::distinct(pathway_id, .keep_all = TRUE)

score_group_stats_one <- score_group_stats %>%
  dplyr::filter(!is.na(pval_wilcox)) %>%
  dplyr::group_by(pathway_id) %>%
  dplyr::arrange(
    pval_wilcox,
    dplyr::desc(abs(diff_median_pCR_minus_non_pCR)),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup()

pathway_box_n <- 24
wilcox_p_cutoff <- 0.05

plot_pathway_table <- target_pathway_annot %>%
  dplyr::left_join(
    score_group_stats_one %>%
      dplyr::select(
        pathway_id,
        pval_wilcox,
        diff_median_pCR_minus_non_pCR
      ),
    by = "pathway_id",
    relationship = "many-to-one"
  ) %>%
  dplyr::filter(
    !is.na(pval_wilcox),
    pval_wilcox < wilcox_p_cutoff
  ) %>%
  dplyr::mutate(
    pval_wilcox_FDR_curated = p.adjust(pval_wilcox, method = "BH"),
    plot_label = stringr::str_wrap(display_label, width = 24)
  ) %>%
  dplyr::arrange(
    pval_wilcox,
    dplyr::desc(abs(diff_median_pCR_minus_non_pCR)),
    gsea_rank
  ) %>%
  dplyr::slice_head(n = pathway_box_n) %>%
  dplyr::arrange(gsea_rank) %>%
  dplyr::mutate(
    plot_label = factor(plot_label, levels = plot_label)
  )

if (nrow(plot_pathway_table) == 0) {
  stop("No curated pathways passed the Wilcoxon P-value cutoff.")
}

plot_pathways <- plot_pathway_table %>%
  dplyr::pull(pathway_id)

pathway_box_df <- score_group_df %>%
  dplyr::filter(pathway_id %in% plot_pathways) %>%
  dplyr::inner_join(
    plot_pathway_table %>%
      dplyr::select(
        pathway_id,
        family_id,
        display_label,
        plot_label,
        main_category,
        Direction,
        NES,
        padj,
        pval_wilcox,
        pval_wilcox_FDR_curated,
        diff_median_pCR_minus_non_pCR
      ),
    by = "pathway_id",
    relationship = "many-to-one"
  ) %>%
  dplyr::mutate(
    TRG_plot = factor(TRG_plot, levels = c("non_pCR", "pCR"))
  )

pathway_p_label <- pathway_box_df %>%
  dplyr::group_by(plot_label) %>%
  dplyr::summarise(
    y_pos = max(ssGSEA_z, na.rm = TRUE) + 0.18 * diff(range(ssGSEA_z, na.rm = TRUE)),
    pval_wilcox = dplyr::first(pval_wilcox),
    pval_wilcox_FDR_curated = dplyr::first(pval_wilcox_FDR_curated),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_label = dplyr::case_when(
      pval_wilcox < 0.001 ~ "Wilcoxon P<0.001",
      TRUE ~ paste0("Wilcoxon P=", signif(pval_wilcox, 2))
    )
  )

p_pathway_box <- pathway_box_df %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = TRG_plot,
      y = ssGSEA_z,
      color = TRG_plot,
      fill = TRG_plot
    )
  ) +
  ggplot2::geom_boxplot(
    outlier.shape = NA,
    width = 0.52,
    alpha = 0.22,
    linewidth = 0.35
  ) +
  ggbeeswarm::geom_quasirandom(
    width = 0.12,
    size = 1.35,
    alpha = 0.85
  ) +
  ggplot2::geom_text(
    data = pathway_p_label,
    ggplot2::aes(
      x = 1.5,
      y = y_pos,
      label = p_label
    ),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  ggplot2::facet_wrap(
    ~ plot_label,
    scales = "free_y",
    ncol = 4
  ) +
  ggplot2::scale_color_manual(
    values = group_cols,
    breaks = c("non_pCR", "pCR"),
    labels = c("non-pCR", "pCR"),
    drop = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = group_cols,
    breaks = c("non_pCR", "pCR"),
    labels = c("non-pCR", "pCR"),
    drop = FALSE
  ) +
  ggplot2::labs(
    x = NULL,
    y = "ssGSEA score, row-wise z-score",
    color = "Response",
    fill = "Response"
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = ggplot2::element_text(size = 8),
    axis.title.y = ggplot2::element_text(size = 9),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 7.4, face = "bold"),
    legend.title = ggplot2::element_text(size = 9),
    legend.text = ggplot2::element_text(size = 8),
    plot.title = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(5.5, 10, 5.5, 5.5)
  )

p_pathway_box

ggplot2::ggsave(
  "figures/host_RNAseq_targeted_pathway_ssGSEA_boxplot_curated_wilcoxon_sig.svg",
  p_pathway_box,
  width = 11,
  height = max(6.5, 2.2 * ceiling(length(plot_pathways) / 4)),
  device = "svg"
)

# write.csv(
#   plot_pathway_table %>%
#     dplyr::select(
#       main_category,
#       family_id,
#       display_label,
#       pathway_id,
#       Direction,
#       NES,
#       padj,
#       pval_wilcox,
#       pval_wilcox_FDR_curated,
#       diff_median_pCR_minus_non_pCR
#     ),
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/ssGSEA_boxplot_curated_wilcoxon_sig_pathways.csv",
#   row.names = FALSE
# )

#-----------------------------------------------------------------#
# 17-2b. Hierarchically clustered heatmap for colorectal adenoma-up genes
#-----------------------------------------------------------------#

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", type = "binary")
}

stopifnot(exists("gsea_plot_df"))
stopifnot(exists("gene_sets_use_long"))
stopifnot(exists("fgsea_res"))
stopifnot(exists("vst_int"))
stopifnot(exists("col_int"))

gene_heatmap_n <- 40

adenoma_annot <- gsea_plot_df %>%
  dplyr::filter(
    stringr::str_detect(
      display_label,
      stringr::regex("adenoma", ignore_case = TRUE)
    )
  ) %>%
  dplyr::slice_head(n = 1)

if (nrow(adenoma_annot) == 0) {
  stop("No colorectal adenoma-related pathway was found in gsea_plot_df.")
}

adenoma_pathway_id <- adenoma_annot$pathway_id[1]
adenoma_gs_name <- adenoma_annot$gs_name[1]

adenoma_genes_all <- gene_sets_use_long %>%
  dplyr::filter(
    pathway_id == adenoma_pathway_id,
    gs_name == adenoma_gs_name
  ) %>%
  dplyr::pull(Gene) %>%
  unique()

adenoma_leading_edge <- fgsea_res %>%
  dplyr::filter(
    pathway_id == adenoma_pathway_id,
    gs_name == adenoma_gs_name
  ) %>%
  dplyr::pull(leadingEdge)

if (length(adenoma_leading_edge) > 0) {
  adenoma_leading_edge <- adenoma_leading_edge[[1]]
} else {
  adenoma_leading_edge <- character(0)
}

adenoma_genes_for_plot <- intersect(adenoma_leading_edge, rownames(vst_int))

if (length(adenoma_genes_for_plot) < 5) {
  adenoma_genes_for_plot <- intersect(adenoma_genes_all, rownames(vst_int))
}

if (length(adenoma_genes_for_plot) < 5) {
  stop("Too few adenoma signature genes were found in vst_int rownames.")
}

sample_annot <- col_int[colnames(vst_int), , drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::mutate(
    TRG_plot = factor(TRG_plot, levels = c("non_pCR", "pCR"))
  )

gene_diff_table <- tibble::tibble(
  Gene = adenoma_genes_for_plot,
  median_non_pCR = apply(
    vst_int[
      adenoma_genes_for_plot,
      sample_annot$Sample[sample_annot$TRG_plot == "non_pCR"],
      drop = FALSE
    ],
    1,
    median,
    na.rm = TRUE
  ),
  median_pCR = apply(
    vst_int[
      adenoma_genes_for_plot,
      sample_annot$Sample[sample_annot$TRG_plot == "pCR"],
      drop = FALSE
    ],
    1,
    median,
    na.rm = TRUE
  )
) %>%
  dplyr::mutate(
    diff_median_non_pCR_minus_pCR = median_non_pCR - median_pCR
  ) %>%
  dplyr::arrange(dplyr::desc(abs(diff_median_non_pCR_minus_pCR)))


#----------------------------#
# Gene-level statistics
#----------------------------#

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("ComplexHeatmap")
}

if (!requireNamespace("circlize", quietly = TRUE)) {
  install.packages("circlize", type = "binary")
}

gene_heatmap_n <- 35
gene_p_cutoff <- 0.1

gene_stat_tile_all <- purrr::map_dfr(
  adenoma_genes_for_plot,
  function(gene) {
    x <- tibble::tibble(
      Sample = colnames(vst_int),
      expr = as.numeric(vst_int[gene, colnames(vst_int)]),
      Response = sample_annot$TRG_plot[
        match(colnames(vst_int), sample_annot$Sample)
      ]
    ) %>%
      dplyr::filter(!is.na(expr), !is.na(Response))
    
    wt <- tryCatch(
      stats::wilcox.test(expr ~ Response, data = x, exact = FALSE),
      error = function(e) NULL
    )
    
    tibble::tibble(
      Gene = gene,
      median_non_pCR = median(x$expr[x$Response == "non_pCR"], na.rm = TRUE),
      median_pCR = median(x$expr[x$Response == "pCR"], na.rm = TRUE),
      logFC_non_pCR_minus_pCR =
        median(x$expr[x$Response == "non_pCR"], na.rm = TRUE) -
        median(x$expr[x$Response == "pCR"], na.rm = TRUE),
      pval_wilcox_gene = ifelse(is.null(wt), NA_real_, wt$p.value)
    )
  }
) %>%
  dplyr::mutate(
    FDR_wilcox_gene = p.adjust(pval_wilcox_gene, method = "BH"),
    neglog10_wilcox_p = -log10(pval_wilcox_gene + 1e-300),
    neglog10_wilcox_p_capped = pmin(neglog10_wilcox_p, 5),
    p_symbol = dplyr::case_when(
      pval_wilcox_gene < 0.001 ~ "***",
      pval_wilcox_gene < 0.01 ~ "**",
      pval_wilcox_gene < 0.05 ~ "*",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::arrange(
    pval_wilcox_gene,
    dplyr::desc(abs(logFC_non_pCR_minus_pCR))
  )

gene_stat_tile <- gene_stat_tile_all %>%
  dplyr::filter(
    !is.na(pval_wilcox_gene),
    pval_wilcox_gene < gene_p_cutoff
  ) %>%
  dplyr::slice_head(n = gene_heatmap_n)

if (nrow(gene_stat_tile) < 5) {
  stop("Fewer than 5 adenoma signature genes passed Wilcoxon P < 0.05.")
}

adenoma_plot_genes <- gene_stat_tile %>%
  dplyr::pull(Gene)

adenoma_expr_z <- vst_int[adenoma_plot_genes, , drop = FALSE]
adenoma_expr_z <- t(scale(t(adenoma_expr_z)))
adenoma_expr_z[adenoma_expr_z > 2.5] <- 2.5
adenoma_expr_z[adenoma_expr_z < -2.5] <- -2.5

gene_stat_tile_for_heatmap <- gene_stat_tile %>%
  dplyr::filter(Gene %in% rownames(adenoma_expr_z)) %>%
  dplyr::select(
    Gene,
    logFC_non_pCR_minus_pCR,
    neglog10_wilcox_p_capped,
    pval_wilcox_gene,
    FDR_wilcox_gene,
    p_symbol
  ) %>%
  tibble::column_to_rownames("Gene")

gene_stat_tile_for_heatmap <- gene_stat_tile_for_heatmap[
  rownames(adenoma_expr_z),
  ,
  drop = FALSE
]

# write.csv(
#   gene_stat_tile_all,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/colorectal_adenoma_signature_gene_logFC_wilcoxon_all.csv",
#   row.names = FALSE
# )
# 
# write.csv(
#   gene_stat_tile,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/colorectal_adenoma_signature_gene_logFC_wilcoxon_P005.csv",
#   row.names = FALSE
# )

#----------------------------#
# Sample annotation
#----------------------------#

adenoma_col_annot <- sample_annot %>%
  dplyr::select(Sample, Response = TRG_plot) %>%
  tibble::column_to_rownames("Sample")

adenoma_col_annot <- adenoma_col_annot[colnames(adenoma_expr_z), , drop = FALSE]

#----------------------------#
# Dendrograms
#----------------------------#

row_cor <- stats::cor(
  t(adenoma_expr_z),
  use = "pairwise.complete.obs"
)
row_cor[is.na(row_cor)] <- 0
diag(row_cor) <- 1

col_cor <- stats::cor(
  adenoma_expr_z,
  use = "pairwise.complete.obs"
)
col_cor[is.na(col_cor)] <- 0
diag(col_cor) <- 1

row_hc <- stats::hclust(
  stats::as.dist(1 - row_cor),
  method = "complete"
)

col_hc <- stats::hclust(
  stats::as.dist(1 - col_cor),
  method = "complete"
)

# Row branch rotation by non-pCR minus pCR logFC.
# This changes only left-right branch orientation, not clustering topology.
row_reorder_weight <- gene_stat_tile_for_heatmap[
  row_hc$labels,
  "logFC_non_pCR_minus_pCR"
]

row_dend <- stats::reorder(
  stats::as.dendrogram(row_hc),
  wts = row_reorder_weight,
  agglo.FUN = mean
)

# Column branch rotation: pCR left, non-pCR right.
# This reverses the previous directionality.
adenoma_sample_score <- colMeans(adenoma_expr_z, na.rm = TRUE)

col_reorder_weight <- adenoma_col_annot %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::mutate(
    Response_weight = dplyr::case_when(
      Response == "pCR" ~ 0,
      Response == "non_pCR" ~ 100,
      TRUE ~ 50
    ),
    Adenoma_score = adenoma_sample_score[Sample],
    reorder_weight =
      Response_weight +
      rank(Adenoma_score, ties.method = "average") / 1000
  ) %>%
  dplyr::select(Sample, reorder_weight) %>%
  tibble::deframe()

col_dend <- stats::reorder(
  stats::as.dendrogram(col_hc),
  wts = col_reorder_weight[col_hc$labels],
  agglo.FUN = mean
)

#----------------------------#
# Colors
#----------------------------#

expr_col_fun <- circlize::colorRamp2(
  c(-2.5, 0, 2.5),
  c("#3B6EA8", "white", "#B14A4A")
)

logFC_cap <- max(
  0.5,
  stats::quantile(
    abs(gene_stat_tile_for_heatmap$logFC_non_pCR_minus_pCR),
    probs = 0.95,
    na.rm = TRUE
  )
)

logFC_col_fun <- circlize::colorRamp2(
  c(-logFC_cap, 0, logFC_cap),
  c(group_cols[["pCR"]], "white", group_cols[["non_pCR"]])
)

p_col_fun <- circlize::colorRamp2(
  c(1.3, 2, 5),
  c("#E8EDF8", "#8EA6D9", "#08519C")
)

p_text_col <- dplyr::case_when(
  gene_stat_tile_for_heatmap$neglog10_wilcox_p_capped >= 2 ~ "white",
  TRUE ~ "grey20"
)

#----------------------------#
# ComplexHeatmap object
# Layout:
# Gene names | logFC tile | P tile with stars | expression heatmap | row dendrogram
#----------------------------#

ha_col <- ComplexHeatmap::HeatmapAnnotation(
  Response = adenoma_col_annot$Response,
  col = list(
    Response = c(
      "non_pCR" = group_cols[["non_pCR"]],
      "pCR" = group_cols[["pCR"]]
    )
  ),
  annotation_name_gp = grid::gpar(fontsize = 8),
  simple_anno_size = grid::unit(3.2, "mm")
)

ha_row_left <- ComplexHeatmap::rowAnnotation(
  Gene = ComplexHeatmap::anno_text(
    rownames(gene_stat_tile_for_heatmap),
    just = "right",
    location = grid::unit(1, "npc"),
    gp = grid::gpar(fontsize = 6.4)
  ),
  `logFC` = gene_stat_tile_for_heatmap$logFC_non_pCR_minus_pCR,
  `P` = ComplexHeatmap::anno_simple(
    gene_stat_tile_for_heatmap$neglog10_wilcox_p_capped,
    col = p_col_fun,
    pch = gene_stat_tile_for_heatmap$p_symbol,
    pt_gp = grid::gpar(
      col = p_text_col,
      fontsize = 6.2,
      fontface = "bold"
    ),
    pt_size = grid::unit(2.5, "mm")
  ),
  col = list(
    `logFC` = logFC_col_fun
  ),
  annotation_name_side = "top",
  annotation_name_rot = 0,
  annotation_name_gp = grid::gpar(fontsize = 7),
  simple_anno_size = grid::unit(4.2, "mm"),
  annotation_width = grid::unit.c(
    grid::unit(38, "mm"),
    grid::unit(4.2, "mm"),
    grid::unit(4.2, "mm")
  )
)

ht_adenoma <- ComplexHeatmap::Heatmap(
  adenoma_expr_z,
  name = "z",
  col = expr_col_fun,
  cluster_rows = row_dend,
  cluster_columns = col_dend,
  left_annotation = ha_row_left,
  top_annotation = ha_col,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_dend_side = "right",
  column_dend_side = "top",
  column_title = "Colorectal adenoma-up signature genes",
  column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Row z-score",
    at = c(-2, -1, 0, 1, 2)
  ),
  use_raster = FALSE
)

#----------------------------#
# Draw on screen first
#----------------------------#

grid::grid.newpage()

ComplexHeatmap::draw(
  ht_adenoma,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = FALSE
)

#----------------------------#
# Save as SVG
#----------------------------#

if (!requireNamespace("svglite", quietly = TRUE)) {
  install.packages("svglite", type = "binary")
}

svglite::svglite(
  "figures/host_RNAseq_colorectal_adenoma_signature_gene_heatmap_clustered_left_gene_logFC_wilcox_P005.svg",
  width = 7,
  height = max(5.2, 0.18 * nrow(adenoma_expr_z) + 1.2)
)

grid::grid.newpage()

ComplexHeatmap::draw(
  ht_adenoma,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = FALSE
)

grDevices::dev.off()

# while (grDevices::dev.cur() > 1) {
#   grDevices::dev.off()
# }



#-----------------------------------------------------------------#
# Save inputs for pathway gene heatmap batch script
#-----------------------------------------------------------------#

dir.create("host_RNAseq/results_clean_metabolite_host/pathway_gene_heatmaps",
           recursive = TRUE, showWarnings = FALSE)

save(
  gsea_plot_df,
  gene_sets_use_long,
  fgsea_res,
  vst_int,
  col_int,
  group_cols,
  file = "host_RNAseq/results_clean_metabolite_host/pathway_gene_heatmaps/pathway_gene_heatmap_inputs.RData"
)




#-----------------------------------------------------------------#
# 17-3. Colorectal adenoma-up genes × differential metabolites
#      Pearson/Spearman correlation and scatter plots
#      with metabolite prevalence filtering and literature-priority genes
#-----------------------------------------------------------------#

stopifnot(exists("adenoma_genes_all"))
stopifnot(exists("adenoma_plot_genes"))
stopifnot(exists("vst_int"))
stopifnot(exists("met_int"))
stopifnot(exists("col_int"))

metabolite_trend_p_cutoff <- 0.30
metabolite_prevalence_cutoff <- 0.50
min_detected_per_group <- 3

gene_met_scatter_n <- 12
gene_met_cor_abs_cutoff <- 0.45
gene_met_cor_p_cutoff <- 0.10

priority_gene_cor_abs_cutoff <- 0.30
priority_gene_cor_p_cutoff <- 0.25

# If raw metabolite values use 0 as non-detected, set this to 0.
# If NA represents missing/non-detected, keep NA_real_.
met_detect_min_value <- NA_real_

literature_priority_genes <- c(
  "NQO1",
  "IGFBP2",
  "CEMIP",
  "FOXQ1",
  "TEAD4",
  "GATA2-AS1",
  "ETV4",
  "NKD1"
)

genes_for_correlation <- unique(
  c(
    adenoma_plot_genes,
    intersect(
      literature_priority_genes,
      intersect(adenoma_genes_all, rownames(vst_int))
    )
  )
)

genes_for_correlation <- intersect(genes_for_correlation, rownames(vst_int))

if (length(genes_for_correlation) < 5) {
  stop("Too few colorectal adenoma-up genes available for correlation.")
}

#----------------------------#
# Prepare metabolite matrix
#----------------------------#

if (length(intersect(rownames(met_int), colnames(vst_int))) >= 5) {
  met_sample_mat <- met_int[
    intersect(rownames(met_int), colnames(vst_int)),
    ,
    drop = FALSE
  ]
} else if (length(intersect(colnames(met_int), colnames(vst_int))) >= 5) {
  met_sample_mat <- t(
    met_int[
      ,
      intersect(colnames(met_int), colnames(vst_int)),
      drop = FALSE
    ]
  )
} else {
  stop("Could not match samples between met_int and vst_int.")
}

met_sample_df <- as.data.frame(met_sample_mat, check.names = FALSE) %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::mutate(
    dplyr::across(
      -Sample,
      ~ suppressWarnings(as.numeric(.x))
    )
  )

metabolite_names <- met_sample_df %>%
  dplyr::select(-Sample) %>%
  dplyr::select(where(is.numeric)) %>%
  colnames()

met_sample_df <- met_sample_df %>%
  dplyr::select(Sample, dplyr::all_of(metabolite_names))

common_samples <- Reduce(
  intersect,
  list(
    colnames(vst_int),
    met_sample_df$Sample,
    rownames(col_int)
  )
)

if (length(common_samples) < 6) {
  stop("Too few matched samples for gene-metabolite correlation.")
}

sample_info_cor <- col_int[common_samples, , drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::mutate(
    TRG_plot = factor(TRG_plot, levels = c("non_pCR", "pCR"))
  )

met_sample_df <- met_sample_df %>%
  dplyr::filter(Sample %in% common_samples)

is_detected_value <- function(x) {
  if (is.na(met_detect_min_value)) {
    return(is.finite(x))
  }
  
  is.finite(x) & x > met_detect_min_value
}

#----------------------------#
# Metabolite prevalence and pCR vs non-pCR Wilcoxon test
#----------------------------#

metabolite_group_stats <- purrr::map_dfr(
  metabolite_names,
  function(met) {
    x <- met_sample_df %>%
      dplyr::select(Sample, value = dplyr::all_of(met)) %>%
      dplyr::inner_join(
        sample_info_cor %>% dplyr::select(Sample, TRG_plot),
        by = "Sample"
      ) %>%
      dplyr::filter(!is.na(TRG_plot))
    
    x <- x %>%
      dplyr::mutate(
        detected = is_detected_value(value)
      )
    
    x_test <- x %>%
      dplyr::filter(detected)
    
    n_non_pCR <- sum(x$TRG_plot == "non_pCR", na.rm = TRUE)
    n_pCR <- sum(x$TRG_plot == "pCR", na.rm = TRUE)
    
    n_detected_non_pCR <- sum(x$detected & x$TRG_plot == "non_pCR", na.rm = TRUE)
    n_detected_pCR <- sum(x$detected & x$TRG_plot == "pCR", na.rm = TRUE)
    
    prevalence_overall <- mean(x$detected, na.rm = TRUE)
    prevalence_non_pCR <- n_detected_non_pCR / n_non_pCR
    prevalence_pCR <- n_detected_pCR / n_pCR
    
    if (
      nrow(x_test) < 6 ||
      dplyr::n_distinct(x_test$TRG_plot) < 2 ||
      stats::sd(x_test$value, na.rm = TRUE) == 0
    ) {
      return(
        tibble::tibble(
          Metabolite = met,
          n_total = nrow(x),
          n_detected = sum(x$detected, na.rm = TRUE),
          n_detected_non_pCR = n_detected_non_pCR,
          n_detected_pCR = n_detected_pCR,
          prevalence_overall = prevalence_overall,
          prevalence_non_pCR = prevalence_non_pCR,
          prevalence_pCR = prevalence_pCR,
          median_non_pCR = NA_real_,
          median_pCR = NA_real_,
          diff_median_pCR_minus_non_pCR = NA_real_,
          pval_wilcox_met = NA_real_
        )
      )
    }
    
    wt <- tryCatch(
      stats::wilcox.test(value ~ TRG_plot, data = x_test, exact = FALSE),
      error = function(e) NULL
    )
    
    tibble::tibble(
      Metabolite = met,
      n_total = nrow(x),
      n_detected = sum(x$detected, na.rm = TRUE),
      n_detected_non_pCR = n_detected_non_pCR,
      n_detected_pCR = n_detected_pCR,
      prevalence_overall = prevalence_overall,
      prevalence_non_pCR = prevalence_non_pCR,
      prevalence_pCR = prevalence_pCR,
      median_non_pCR = median(x_test$value[x_test$TRG_plot == "non_pCR"], na.rm = TRUE),
      median_pCR = median(x_test$value[x_test$TRG_plot == "pCR"], na.rm = TRUE),
      diff_median_pCR_minus_non_pCR =
        median(x_test$value[x_test$TRG_plot == "pCR"], na.rm = TRUE) -
        median(x_test$value[x_test$TRG_plot == "non_pCR"], na.rm = TRUE),
      pval_wilcox_met = ifelse(is.null(wt), NA_real_, wt$p.value)
    )
  }
) %>%
  dplyr::mutate(
    FDR_wilcox_met = p.adjust(pval_wilcox_met, method = "BH"),
    pass_prevalence =
      prevalence_overall >= metabolite_prevalence_cutoff &
      n_detected_non_pCR >= min_detected_per_group &
      n_detected_pCR >= min_detected_per_group
  ) %>%
  dplyr::arrange(pval_wilcox_met)

# write.csv(
#   metabolite_group_stats,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/metabolite_pCR_vs_non_pCR_wilcoxon_prevalence.csv",
#   row.names = FALSE
# )

#----------------------------#
# Functional metabolite annotation
#----------------------------#

annotate_metabolite_axis <- function(x) {
  xu <- toupper(x)
  
  dplyr::case_when(
    stringr::str_detect(
      xu,
      "CHOL|CHOLATE|CHOLIC|DEOXYCHOL|LITHOCHOL|CHENODEOXY|URSODEOXY|GLYCOCHOL|TAUROCHOL|BILE"
    ) ~ "Bile acid",
    
    stringr::str_detect(
      xu,
      "ACETATE|PROPIONATE|BUTYRATE|VALERATE|ISOBUTYRATE|ISOVALERATE|SCFA"
    ) ~ "SCFA",
    
    stringr::str_detect(
      xu,
      "INDOLE|TRYPTOPHAN|KYNURENINE|IAA|IPA|ILA|SKATOLE"
    ) ~ "Indole/tryptophan",
    
    stringr::str_detect(
      xu,
      "HISTIDINE|HISTAMINE|UROCANATE|IMIDAZOLE"
    ) ~ "Histidine/imidazole",
    
    stringr::str_detect(
      xu,
      "LEUCINE|ISOLEUCINE|VALINE|BCAA|BRANCHED"
    ) ~ "BCAA",
    
    stringr::str_detect(
      xu,
      "ALANINE|GLYCINE|SERINE|THREONINE|METHIONINE|PHENYLALANINE|TYROSINE|ARGININE|ORNITHINE|CITRULLINE|PROLINE|ASPARTATE|GLUTAMATE|GLUTAMINE|LYSINE|AMINO"
    ) ~ "Amino acid",
    
    stringr::str_detect(
      xu,
      "SUCROSE|GLUCOSE|FRUCTOSE|MALTOSE|LACTOSE|STARCH|RIBOSE|XYLOSE|ARABINOSE|GALACTOSE|GLUCAN"
    ) ~ "Carbohydrate/sugar",
    
    stringr::str_detect(
      xu,
      "ACETAMINOPHEN|PARACETAMOL"
    ) ~ "Drug/exposure",
    
    stringr::str_detect(
      xu,
      "BENZOATE|PHENOL|CRESOL|SULFATE|GLUTATHIONE|CYSTEINE|TAURINE|HEME|PORPHYRIN|XENOBIOTIC"
    ) ~ "Xenobiotic/redox proxy",
    
    TRUE ~ "Other"
  )
}

#----------------------------#
# Select differential and sufficiently prevalent metabolites
#----------------------------#

metabolites_for_correlation <- metabolite_group_stats %>%
  dplyr::filter(
    pass_prevalence,
    !is.na(pval_wilcox_met),
    pval_wilcox_met < metabolite_trend_p_cutoff
  ) %>%
  dplyr::mutate(
    metabolite_axis = annotate_metabolite_axis(Metabolite),
    is_functionally_interpretable = metabolite_axis != "Other"
  ) %>%
  dplyr::arrange(
    dplyr::desc(is_functionally_interpretable),
    pval_wilcox_met,
    dplyr::desc(abs(diff_median_pCR_minus_non_pCR))
  ) %>%
  dplyr::pull(Metabolite)

if (length(metabolites_for_correlation) < 3) {
  stop("Fewer than 3 metabolites passed prevalence and pCR vs non-pCR trend filters.")
}

message("Metabolites used for gene-metabolite correlation: ", length(metabolites_for_correlation))

#----------------------------#
# Correlation helper
#----------------------------#

cor_one_pair <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  
  if (
    sum(keep) < 6 ||
    stats::sd(x[keep], na.rm = TRUE) == 0 ||
    stats::sd(y[keep], na.rm = TRUE) == 0
  ) {
    return(
      tibble::tibble(
        estimate = NA_real_,
        pval = NA_real_,
        n = sum(keep)
      )
    )
  }
  
  ct <- suppressWarnings(
    tryCatch(
      stats::cor.test(x[keep], y[keep], method = method, exact = FALSE),
      error = function(e) NULL
    )
  )
  
  if (is.null(ct)) {
    return(
      tibble::tibble(
        estimate = NA_real_,
        pval = NA_real_,
        n = sum(keep)
      )
    )
  }
  
  tibble::tibble(
    estimate = unname(ct$estimate),
    pval = ct$p.value,
    n = sum(keep)
  )
}

#----------------------------#
# Gene and metabolite z-score matrices
#----------------------------#

gene_expr_for_cor <- vst_int[genes_for_correlation, common_samples, drop = FALSE]
gene_expr_for_cor_z <- t(scale(t(gene_expr_for_cor)))
gene_expr_for_cor_z[gene_expr_for_cor_z > 2.5] <- 2.5
gene_expr_for_cor_z[gene_expr_for_cor_z < -2.5] <- -2.5

gene_z_df <- as.data.frame(
  t(gene_expr_for_cor_z),
  check.names = FALSE
) %>%
  tibble::rownames_to_column("Sample")

met_z_df <- met_sample_df %>%
  dplyr::filter(Sample %in% common_samples) %>%
  dplyr::select(Sample, dplyr::all_of(metabolites_for_correlation))

met_z_mat <- as.matrix(met_z_df[, metabolites_for_correlation, drop = FALSE])
met_z_mat <- scale(met_z_mat)

met_z_df <- as.data.frame(met_z_mat, check.names = FALSE) %>%
  dplyr::mutate(Sample = met_z_df$Sample) %>%
  dplyr::select(Sample, dplyr::everything())

#----------------------------#
# Pearson and Spearman correlations
#----------------------------#

gene_met_cor <- tidyr::expand_grid(
  Gene = genes_for_correlation,
  Metabolite = metabolites_for_correlation
) %>%
  dplyr::mutate(
    pearson = purrr::map2(
      Gene,
      Metabolite,
      ~ cor_one_pair(
        gene_z_df[[.x]],
        met_z_df[[.y]],
        method = "pearson"
      )
    ),
    spearman = purrr::map2(
      Gene,
      Metabolite,
      ~ cor_one_pair(
        gene_z_df[[.x]],
        met_z_df[[.y]],
        method = "spearman"
      )
    )
  ) %>%
  tidyr::unnest_wider(pearson, names_sep = "_") %>%
  tidyr::unnest_wider(spearman, names_sep = "_") %>%
  dplyr::rename(
    r_pearson = pearson_estimate,
    pval_pearson = pearson_pval,
    n_pearson = pearson_n,
    rho_spearman = spearman_estimate,
    pval_spearman = spearman_pval,
    n_spearman = spearman_n
  ) %>%
  dplyr::left_join(
    metabolite_group_stats %>%
      dplyr::mutate(
        metabolite_axis = annotate_metabolite_axis(Metabolite)
      ),
    by = "Metabolite",
    relationship = "many-to-one"
  ) %>%
  dplyr::mutate(
    FDR_pearson = p.adjust(pval_pearson, method = "BH"),
    FDR_spearman = p.adjust(pval_spearman, method = "BH"),
    concordant_direction = sign(r_pearson) == sign(rho_spearman),
    best_abs_coef = dplyr::case_when(
      is.na(r_pearson) & is.na(rho_spearman) ~ NA_real_,
      TRUE ~ pmax(abs(r_pearson), abs(rho_spearman), na.rm = TRUE)
    ),
    best_pval = dplyr::case_when(
      is.na(pval_pearson) & is.na(pval_spearman) ~ NA_real_,
      TRUE ~ pmin(pval_pearson, pval_spearman, na.rm = TRUE)
    ),
    is_literature_priority_gene = Gene %in% literature_priority_genes,
    functional_priority = dplyr::case_when(
      metabolite_axis %in% c(
        "Indole/tryptophan",
        "Histidine/imidazole",
        "BCAA",
        "Amino acid",
        "SCFA",
        "Carbohydrate/sugar",
        "Xenobiotic/redox proxy"
      ) ~ 2,
      metabolite_axis == "Bile acid" ~ 1,
      TRUE ~ 0
    ),
    ranking_score =
      functional_priority +
      as.numeric(is_literature_priority_gene) * 1.5 +
      best_abs_coef +
      (-log10(best_pval + 1e-300) / 5) +
      (-log10(pval_wilcox_met + 1e-300) / 10)
  ) %>%
  dplyr::arrange(
    dplyr::desc(is_literature_priority_gene),
    dplyr::desc(functional_priority),
    dplyr::desc(concordant_direction),
    best_pval,
    dplyr::desc(best_abs_coef),
    pval_wilcox_met
  )

# write.csv(
#   gene_met_cor,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/adenoma_gene_differential_metabolite_correlation_pearson_spearman_prevalence_filtered.csv",
#   row.names = FALSE
# )

#----------------------------#
# Select top pairs for scatter plots
#----------------------------#

scatter_candidates_standard <- gene_met_cor %>%
  dplyr::filter(
    !is.na(best_pval),
    !is.na(best_abs_coef),
    pass_prevalence,
    pval_wilcox_met < metabolite_trend_p_cutoff,
    concordant_direction,
    best_abs_coef >= gene_met_cor_abs_cutoff,
    best_pval < gene_met_cor_p_cutoff
  ) %>%
  dplyr::arrange(
    dplyr::desc(functional_priority),
    best_pval,
    dplyr::desc(best_abs_coef),
    pval_wilcox_met
  )

scatter_candidates_priority <- gene_met_cor %>%
  dplyr::filter(
    Gene %in% literature_priority_genes,
    !is.na(best_pval),
    !is.na(best_abs_coef),
    pass_prevalence,
    pval_wilcox_met < metabolite_trend_p_cutoff,
    concordant_direction,
    best_abs_coef >= priority_gene_cor_abs_cutoff,
    best_pval < priority_gene_cor_p_cutoff
  ) %>%
  dplyr::arrange(
    dplyr::desc(is_literature_priority_gene),
    best_pval,
    dplyr::desc(best_abs_coef),
    pval_wilcox_met
  ) %>%
  dplyr::slice_head(n = 6)

scatter_candidates <- dplyr::bind_rows(
  scatter_candidates_priority,
  scatter_candidates_standard
) %>%
  dplyr::distinct(Gene, Metabolite, .keep_all = TRUE) %>%
  dplyr::group_by(Metabolite) %>%
  dplyr::slice_head(n = 2) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Gene) %>%
  dplyr::slice_head(n = 2) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    dplyr::desc(is_literature_priority_gene),
    dplyr::desc(functional_priority),
    best_pval,
    dplyr::desc(best_abs_coef),
    pval_wilcox_met
  ) %>%
  dplyr::slice_head(n = gene_met_scatter_n) %>%
  dplyr::mutate(
    Metabolite_clean = stringr::str_replace_all(Metabolite, "_", " "),
    Metabolite_clean = stringr::str_trunc(Metabolite_clean, width = 30),
    pair_label = paste0(Gene, "\n", Metabolite_clean),
    stat_label = paste0(
      "Pearson r=", sprintf("%.2f", r_pearson),
      ", P=", signif(pval_pearson, 2),
      "\nSpearman rho=", sprintf("%.2f", rho_spearman),
      ", P=", signif(pval_spearman, 2),
      "\nMet P=", signif(pval_wilcox_met, 2),
      ", prev=", sprintf("%.0f", 100 * prevalence_overall), "%"
    )
  )

# write.csv(
#   scatter_candidates,
#   "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/adenoma_gene_differential_metabolite_scatter_candidates_prevalence_filtered.csv",
#   row.names = FALSE
# )

print(
  scatter_candidates %>%
    dplyr::select(
      Gene,
      Metabolite,
      metabolite_axis,
      prevalence_overall,
      n_detected_non_pCR,
      n_detected_pCR,
      r_pearson,
      pval_pearson,
      rho_spearman,
      pval_spearman,
      pval_wilcox_met,
      diff_median_pCR_minus_non_pCR,
      is_literature_priority_gene,
      best_abs_coef,
      best_pval
    )
)

#----------------------------#
# Scatter plot with linear regression
#----------------------------#

scatter_plot_df <- scatter_candidates %>%
  dplyr::select(
    Gene,
    Metabolite,
    pair_label
  ) %>%
  tidyr::expand_grid(Sample = common_samples) %>%
  dplyr::mutate(
    gene_z = purrr::map2_dbl(
      Gene,
      Sample,
      ~ gene_z_df[gene_z_df$Sample == .y, .x, drop = TRUE]
    ),
    metabolite_z = purrr::map2_dbl(
      Metabolite,
      Sample,
      ~ met_z_df[met_z_df$Sample == .y, .x, drop = TRUE]
    )
  ) %>%
  dplyr::left_join(
    sample_info_cor %>%
      dplyr::select(Sample, TRG_plot),
    by = "Sample"
  ) %>%
  dplyr::left_join(
    scatter_candidates %>%
      dplyr::select(pair_label, stat_label, metabolite_axis),
    by = "pair_label",
    relationship = "many-to-one"
  ) %>%
  dplyr::mutate(
    pair_label = factor(pair_label, levels = scatter_candidates$pair_label)
  )

scatter_stat_df <- scatter_candidates %>%
  dplyr::mutate(
    pair_label = factor(pair_label, levels = scatter_candidates$pair_label),
    x = -Inf,
    y = Inf
  )

p_adenoma_gene_met_scatter <- scatter_plot_df %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = metabolite_z,
      y = gene_z
    )
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = TRG_plot),
    size = 2.0,
    alpha = 0.85
  ) +
  ggplot2::geom_smooth(
    ggplot2::aes(group = 1),
    method = "lm",
    se = TRUE,
    linewidth = 0.45,
    color = "grey25"
  ) +
  ggplot2::geom_text(
    data = scatter_stat_df,
    ggplot2::aes(
      x = x,
      y = y,
      label = stat_label
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.05,
    size = 2.3
  ) +
  ggplot2::facet_wrap(
    ~ pair_label,
    scales = "free",
    ncol = 4
  ) +
  ggplot2::scale_color_manual(
    values = group_cols,
    breaks = c("non_pCR", "pCR"),
    labels = c("non-pCR", "pCR"),
    name = "Response"
  ) +
  ggplot2::labs(
    x = "Metabolite abundance, z-score",
    y = "Adenoma-up gene expression, row z-score"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 7.2, face = "bold"),
    axis.text = ggplot2::element_text(size = 8),
    axis.title = ggplot2::element_text(size = 9),
    legend.title = ggplot2::element_text(size = 9),
    legend.text = ggplot2::element_text(size = 8),
    plot.title = ggplot2::element_blank()
  )

p_adenoma_gene_met_scatter

ggplot2::ggsave(
  "figures/host_RNAseq_adenoma_signature_gene_differential_metabolite_scatter_top_prevalence_filtered.svg",
  p_adenoma_gene_met_scatter,
  width = 11,
  height = max(6, 2.8 * ceiling(nrow(scatter_candidates) / 4)),
  device = "svg"
)


#-----------------------------------------------------------------#
# 17-4. ILA-HK2/glycolysis axis check in human data
#-----------------------------------------------------------------#

stopifnot(exists("vst_int"))
stopifnot(exists("met_int"))
stopifnot(exists("col_int"))

#----------------------------#
# User-adjustable settings
#----------------------------#

# If zero means non-detected in metabolite table, keep 0.
# If NA is the only missing/non-detected marker, set this to NA_real_.
met_detect_min_value <- 0

use_log1p_metabolite <- FALSE

ila_gene_cor_p_cutoff <- 0.10
ila_gene_cor_abs_cutoff <- 0.30

#----------------------------#
# Target genes from the ILA paper figure
#----------------------------#

ila_glycolysis_genes_raw <- c(
  "HK2", "GPI", "PFKL", "ALDOA", "GAPDH",
  "PGK1", "ENO1", "PKM", "PKM2", "LDHA"
)

ila_oxphos_genes_raw <- c(
  "PDHA1", "PDHA", "PDHB", "CS", "ACO2",
  "IDH2", "DLST", "SDHB", "SUCLG2"
)

# Resolve a few figure-label aliases to HGNC symbols likely used in vst_int.
gene_alias_table <- tibble::tribble(
  ~input_gene, ~hgnc_gene,
  "PKM2", "PKM",
  "PDHA", "PDHA1"
)

resolve_gene_symbols <- function(x, row_ids) {
  x2 <- unique(c(
    x,
    gene_alias_table$hgnc_gene[match(x, gene_alias_table$input_gene)]
  ))
  x2 <- x2[!is.na(x2)]
  intersect(x2, row_ids)
}

ila_glycolysis_genes <- resolve_gene_symbols(
  ila_glycolysis_genes_raw,
  rownames(vst_int)
)

ila_oxphos_genes <- resolve_gene_symbols(
  ila_oxphos_genes_raw,
  rownames(vst_int)
)

ila_target_genes <- unique(c(
  ila_glycolysis_genes,
  ila_oxphos_genes
))

if (!"HK2" %in% rownames(vst_int)) {
  warning("HK2 was not found in rownames(vst_int). Check gene identifiers.")
}

if (length(ila_glycolysis_genes) < 3) {
  stop("Too few glycolysis genes were found in vst_int rownames.")
}

message("Glycolysis genes found: ", paste(ila_glycolysis_genes, collapse = ", "))
message("OXPHOS genes found: ", paste(ila_oxphos_genes, collapse = ", "))

#----------------------------#
# Prepare metabolite matrix and find ILA variable
#----------------------------#

if (length(intersect(rownames(met_int), colnames(vst_int))) >= 5) {
  met_sample_mat <- met_int[
    intersect(rownames(met_int), colnames(vst_int)),
    ,
    drop = FALSE
  ]
} else if (length(intersect(colnames(met_int), colnames(vst_int))) >= 5) {
  met_sample_mat <- t(
    met_int[
      ,
      intersect(colnames(met_int), colnames(vst_int)),
      drop = FALSE
    ]
  )
} else {
  stop("Could not match samples between met_int and vst_int.")
}

met_sample_df <- as.data.frame(met_sample_mat, check.names = FALSE) %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::mutate(
    dplyr::across(
      -Sample,
      ~ suppressWarnings(as.numeric(.x))
    )
  )

metabolite_names <- met_sample_df %>%
  dplyr::select(-Sample) %>%
  dplyr::select(where(is.numeric)) %>%
  colnames()

ila_candidates <- metabolite_names[
  stringr::str_detect(
    metabolite_names,
    stringr::regex("indole.*lactic|indole_lactic|indole-3-lactic|\\bILA\\b", ignore_case = TRUE)
  )
]

print(ila_candidates)

if ("Indole_lactic_acid" %in% ila_candidates) {
  ila_metabolite <- "Indole_lactic_acid"
} else if ("Indole-3-lactic_acid" %in% ila_candidates) {
  ila_metabolite <- "Indole-3-lactic_acid"
} else if (length(ila_candidates) >= 1) {
  ila_metabolite <- ila_candidates[1]
} else {
  stop("No ILA-like metabolite was found. Check metabolite names manually.")
}

message("Selected ILA metabolite: ", ila_metabolite)

#----------------------------#
# Match samples
#----------------------------#

common_samples <- Reduce(
  intersect,
  list(
    colnames(vst_int),
    met_sample_df$Sample,
    rownames(col_int)
  )
)

if (length(common_samples) < 6) {
  stop("Too few matched samples for ILA-gene analysis.")
}

sample_info_ila <- col_int[common_samples, , drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::mutate(
    TRG_plot = factor(TRG_plot, levels = c("non_pCR", "pCR"))
  )

ila_df <- met_sample_df %>%
  dplyr::filter(Sample %in% common_samples) %>%
  dplyr::select(Sample, ILA = dplyr::all_of(ila_metabolite)) %>%
  dplyr::inner_join(
    sample_info_ila %>% dplyr::select(Sample, TRG_plot),
    by = "Sample"
  ) %>%
  dplyr::mutate(
    ILA_for_cor = dplyr::case_when(
      use_log1p_metabolite ~ log1p(ILA),
      TRUE ~ ILA
    ),
    ILA_detected = dplyr::case_when(
      is.na(met_detect_min_value) ~ is.finite(ILA),
      TRUE ~ is.finite(ILA) & ILA > met_detect_min_value
    )
  )

ila_prevalence_table <- ila_df %>%
  dplyr::group_by(TRG_plot) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_detected = sum(ILA_detected, na.rm = TRUE),
    prevalence = mean(ILA_detected, na.rm = TRUE),
    median_ILA_all = median(ILA, na.rm = TRUE),
    median_ILA_detected = median(ILA[ILA_detected], na.rm = TRUE),
    .groups = "drop"
  )

print(ila_prevalence_table)

ila_wilcox <- tryCatch(
  stats::wilcox.test(ILA ~ TRG_plot, data = ila_df, exact = FALSE),
  error = function(e) NULL
)

ila_group_test <- tibble::tibble(
  Metabolite = ila_metabolite,
  n_total = nrow(ila_df),
  n_detected = sum(ila_df$ILA_detected, na.rm = TRUE),
  prevalence_overall = mean(ila_df$ILA_detected, na.rm = TRUE),
  median_pCR = median(ila_df$ILA[ila_df$TRG_plot == "pCR"], na.rm = TRUE),
  median_non_pCR = median(ila_df$ILA[ila_df$TRG_plot == "non_pCR"], na.rm = TRUE),
  diff_median_pCR_minus_non_pCR =
    median(ila_df$ILA[ila_df$TRG_plot == "pCR"], na.rm = TRUE) -
    median(ila_df$ILA[ila_df$TRG_plot == "non_pCR"], na.rm = TRUE),
  pval_wilcox_ILA = ifelse(is.null(ila_wilcox), NA_real_, ila_wilcox$p.value)
)

print(ila_group_test)

#----------------------------#
# Gene z-score and module scores
#----------------------------#

gene_expr <- vst_int[ila_target_genes, common_samples, drop = FALSE]
gene_expr_z <- t(scale(t(gene_expr)))
gene_expr_z[gene_expr_z > 2.5] <- 2.5
gene_expr_z[gene_expr_z < -2.5] <- -2.5

gene_z_df <- as.data.frame(t(gene_expr_z), check.names = FALSE) %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::inner_join(
    ila_df %>% dplyr::select(Sample, TRG_plot, ILA, ILA_for_cor, ILA_detected),
    by = "Sample"
  )

module_score_df <- gene_z_df %>%
  dplyr::mutate(
    Glycolysis_score = rowMeans(
      dplyr::across(dplyr::all_of(ila_glycolysis_genes)),
      na.rm = TRUE
    ),
    OXPHOS_score = dplyr::case_when(
      length(ila_oxphos_genes) >= 3 ~ rowMeans(
        dplyr::across(dplyr::all_of(ila_oxphos_genes)),
        na.rm = TRUE
      ),
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::select(
    Sample,
    TRG_plot,
    ILA,
    ILA_for_cor,
    ILA_detected,
    Glycolysis_score,
    OXPHOS_score
  )

#----------------------------#
# Correlation helper
#----------------------------#

cor_one_pair <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  
  if (
    sum(keep) < 6 ||
    stats::sd(x[keep], na.rm = TRUE) == 0 ||
    stats::sd(y[keep], na.rm = TRUE) == 0
  ) {
    return(
      tibble::tibble(
        estimate = NA_real_,
        pval = NA_real_,
        n = sum(keep)
      )
    )
  }
  
  ct <- suppressWarnings(
    tryCatch(
      stats::cor.test(x[keep], y[keep], method = method, exact = FALSE),
      error = function(e) NULL
    )
  )
  
  if (is.null(ct)) {
    return(
      tibble::tibble(
        estimate = NA_real_,
        pval = NA_real_,
        n = sum(keep)
      )
    )
  }
  
  tibble::tibble(
    estimate = unname(ct$estimate),
    pval = ct$p.value,
    n = sum(keep)
  )
}

#----------------------------#
# Gene-level correlation: all samples and detected-only
#----------------------------#

ila_gene_cor <- tidyr::expand_grid(
  Gene = ila_target_genes,
  analysis_set = c("all_samples", "ILA_detected_only")
) %>%
  dplyr::mutate(
    gene_class = dplyr::case_when(
      Gene %in% ila_glycolysis_genes ~ "Glycolysis",
      Gene %in% ila_oxphos_genes ~ "OXPHOS",
      TRUE ~ "Other"
    ),
    dat = purrr::map(
      analysis_set,
      ~ {
        if (.x == "ILA_detected_only") {
          gene_z_df %>% dplyr::filter(ILA_detected)
        } else {
          gene_z_df
        }
      }
    ),
    pearson = purrr::map2(
      Gene,
      dat,
      ~ cor_one_pair(
        .y$ILA_for_cor,
        .y[[.x]],
        method = "pearson"
      )
    ),
    spearman = purrr::map2(
      Gene,
      dat,
      ~ cor_one_pair(
        .y$ILA_for_cor,
        .y[[.x]],
        method = "spearman"
      )
    )
  ) %>%
  dplyr::select(-dat) %>%
  tidyr::unnest_wider(pearson, names_sep = "_") %>%
  tidyr::unnest_wider(spearman, names_sep = "_") %>%
  dplyr::rename(
    r_pearson = pearson_estimate,
    pval_pearson = pearson_pval,
    n_pearson = pearson_n,
    rho_spearman = spearman_estimate,
    pval_spearman = spearman_pval,
    n_spearman = spearman_n
  ) %>%
  dplyr::group_by(analysis_set) %>%
  dplyr::mutate(
    FDR_pearson = p.adjust(pval_pearson, method = "BH"),
    FDR_spearman = p.adjust(pval_spearman, method = "BH")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    ILA_expected_negative =
      gene_class == "Glycolysis" &
      rho_spearman < 0,
    best_abs_coef = pmax(abs(r_pearson), abs(rho_spearman), na.rm = TRUE),
    best_pval = pmin(pval_pearson, pval_spearman, na.rm = TRUE)
  ) %>%
  dplyr::arrange(
    analysis_set,
    gene_class,
    pval_spearman,
    dplyr::desc(abs(rho_spearman))
  )

#----------------------------#
# Gene-level pCR vs non-pCR Wilcoxon test
#----------------------------#

ila_gene_group_test <- purrr::map_dfr(
  ila_target_genes,
  function(gene) {
    x <- gene_z_df %>%
      dplyr::select(Sample, TRG_plot, expr_z = dplyr::all_of(gene)) %>%
      dplyr::filter(!is.na(expr_z), !is.na(TRG_plot))
    
    wt <- tryCatch(
      stats::wilcox.test(expr_z ~ TRG_plot, data = x, exact = FALSE),
      error = function(e) NULL
    )
    
    tibble::tibble(
      Gene = gene,
      gene_class = dplyr::case_when(
        gene %in% ila_glycolysis_genes ~ "Glycolysis",
        gene %in% ila_oxphos_genes ~ "OXPHOS",
        TRUE ~ "Other"
      ),
      median_pCR = median(x$expr_z[x$TRG_plot == "pCR"], na.rm = TRUE),
      median_non_pCR = median(x$expr_z[x$TRG_plot == "non_pCR"], na.rm = TRUE),
      diff_median_pCR_minus_non_pCR =
        median(x$expr_z[x$TRG_plot == "pCR"], na.rm = TRUE) -
        median(x$expr_z[x$TRG_plot == "non_pCR"], na.rm = TRUE),
      pval_wilcox_gene = ifelse(is.null(wt), NA_real_, wt$p.value)
    )
  }
) %>%
  dplyr::mutate(
    FDR_wilcox_gene = p.adjust(pval_wilcox_gene, method = "BH"),
    direction_consistent_with_ILA_suppression =
      gene_class == "Glycolysis" &
      ila_group_test$diff_median_pCR_minus_non_pCR > 0 &
      diff_median_pCR_minus_non_pCR < 0
  ) %>%
  dplyr::arrange(gene_class, pval_wilcox_gene)

#----------------------------#
# Module-level correlation and group test
#----------------------------#

module_cor <- tidyr::expand_grid(
  Module = c("Glycolysis_score", "OXPHOS_score"),
  analysis_set = c("all_samples", "ILA_detected_only")
) %>%
  dplyr::filter(
    Module %in% colnames(module_score_df)
  ) %>%
  dplyr::mutate(
    dat = purrr::map(
      analysis_set,
      ~ {
        if (.x == "ILA_detected_only") {
          module_score_df %>% dplyr::filter(ILA_detected)
        } else {
          module_score_df
        }
      }
    ),
    pearson = purrr::map2(
      Module,
      dat,
      ~ cor_one_pair(.y$ILA_for_cor, .y[[.x]], method = "pearson")
    ),
    spearman = purrr::map2(
      Module,
      dat,
      ~ cor_one_pair(.y$ILA_for_cor, .y[[.x]], method = "spearman")
    )
  ) %>%
  dplyr::select(-dat) %>%
  tidyr::unnest_wider(pearson, names_sep = "_") %>%
  tidyr::unnest_wider(spearman, names_sep = "_") %>%
  dplyr::rename(
    r_pearson = pearson_estimate,
    pval_pearson = pearson_pval,
    n_pearson = pearson_n,
    rho_spearman = spearman_estimate,
    pval_spearman = spearman_pval,
    n_spearman = spearman_n
  )

module_group_test <- purrr::map_dfr(
  c("Glycolysis_score", "OXPHOS_score"),
  function(mod) {
    if (!mod %in% colnames(module_score_df)) {
      return(NULL)
    }
    
    x <- module_score_df %>%
      dplyr::select(Sample, TRG_plot, score = dplyr::all_of(mod)) %>%
      dplyr::filter(!is.na(score), !is.na(TRG_plot))
    
    wt <- tryCatch(
      stats::wilcox.test(score ~ TRG_plot, data = x, exact = FALSE),
      error = function(e) NULL
    )
    
    tibble::tibble(
      Module = mod,
      median_pCR = median(x$score[x$TRG_plot == "pCR"], na.rm = TRUE),
      median_non_pCR = median(x$score[x$TRG_plot == "non_pCR"], na.rm = TRUE),
      diff_median_pCR_minus_non_pCR =
        median(x$score[x$TRG_plot == "pCR"], na.rm = TRUE) -
        median(x$score[x$TRG_plot == "non_pCR"], na.rm = TRUE),
      pval_wilcox_module = ifelse(is.null(wt), NA_real_, wt$p.value)
    )
  }
) %>%
  dplyr::mutate(
    FDR_wilcox_module = p.adjust(pval_wilcox_module, method = "BH")
  )

#----------------------------#
# Save result tables
#----------------------------#

ila_gene_result <- ila_gene_cor %>%
  dplyr::left_join(
    ila_gene_group_test,
    by = c("Gene", "gene_class"),
    relationship = "many-to-one"
  )

ila_gene_result %>% 
  as.data.frame() %>% 
  arrange(pval_spearman) %>% 
  head()

write.csv(
  ila_prevalence_table,
  "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/ILA_prevalence_by_response.csv",
  row.names = FALSE
)

write.csv(
  ila_group_test,
  "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/ILA_response_group_wilcoxon.csv",
  row.names = FALSE
)

write.csv(
  ila_gene_result,
  "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/ILA_HK2_glycolysis_gene_correlation_and_group_test.csv",
  row.names = FALSE
)

write.csv(
  module_cor,
  "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/ILA_glycolysis_module_correlation.csv",
  row.names = FALSE
)

write.csv(
  module_group_test,
  "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/ILA_glycolysis_module_response_group_test.csv",
  row.names = FALSE
)

print(ila_gene_result)
print(module_cor)
print(module_group_test)




#-----------------------------------------------------------------#
# 18. Save objects and session info
#-----------------------------------------------------------------#

save(
  all_gene_sets_long,
  gene_set_size_table,
  gene_sets_use_long,
  gene_sets_use,
  fgsea_res,
  ora_res,
  score_rna,
  score_rna_z,
  score_int,
  score_int_z,
  score_group_stats,
  pathway_met_cor,
  top_pathway_met_cor,
  file = "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/host_metabolite_pathway_analysis_results.RData"
)

save(
  deg_rna,
  vst_rna,
  vst_int,
  col_rna,
  col_int,
  msig,
  group_cols,
  file = "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/global_gsea_inputs.RData"
)

#-----------------------------------------------------------------#
# Save inputs for Figure 5 host RNA-seq script
#-----------------------------------------------------------------#
dir.create(
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq",
  recursive = TRUE,
  showWarnings = FALSE
)

save(
  deg_rna,
  vst_rna,
  vst_int,
  col_rna,
  col_int,
  met_int,
  gsea_plot_df,
  gene_sets_use_long,
  fgsea_res,
  group_cols,
  file = "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5_host_RNAseq_inputs.RData"
)



sink("host_RNAseq/results_clean_metabolite_host/pathway_metabolite/sessionInfo_pathway_metabolite.txt")
print(sessionInfo())
sink()

message("Analysis completed.")
message("Results saved in: host_RNAseq/results_clean_metabolite_host/pathway_metabolite")
message("Figures saved in: figures")

#-----------------------------------------------------------------#
# End
#-----------------------------------------------------------------#
