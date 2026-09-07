#-----------------------------------------------------------------#
#
# Figure 5. Host RNA-seq functional landscape
#
# Fig. 5A. Enrichment map / network
# Fig. 5B. Selected pathway gene heatmaps
# Fig. 5C. Microbe–metabolite–host gene scatter
#
#-----------------------------------------------------------------#

#-----------------------------------------------------------------#
#
# Fig. 5A. Host RNA-seq enrichment map
#
# Clean version:
#   - unified MSigDB GSEA only
#   - no bind with targeted fgsea_res
#   - no forced network_cluster x biological_axis splitting
#   - node color = pCR-high / non-pCR-high
#   - node size = -log10(FDR)
#   - edge = leading-edge gene overlap
#   - hull = Louvain community direction
#       the first plot retains the conservative 70% node/evidence rule
#       an additional plot uses a simple >= 80% majority rule
#       communities below 80% agreement are shown as mixed (gray)
#       these are descriptive visualization rules, not cluster-level tests
#   - text = Louvain cluster-level functional program
#
# Output:
#   figures/Fig5A_host_RNAseq_enrichment_map_mixed_direction.svg
#   figures/Fig5A_host_RNAseq_enrichment_map_majority80.svg
#   CSV audits of hypoxia/oxidative-stress pathways and overlapping genes
#
#-----------------------------------------------------------------#

rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("D:/2-연구/2-CRC metagenomics/")

#-----------------------------------------------------------------#
# 0. Packages
#-----------------------------------------------------------------#

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_pkgs <- c(
  "dplyr", "tidyr", "tibble", "stringr", "purrr", "readr",
  "ggplot2", "ggrepel", "forcats", "msigdbr",
  "igraph", "ggraph", "ggforce", "svglite"
)

bioc_pkgs <- c("fgsea")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, type = "binary")
  }
}

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(forcats)
  library(msigdbr)
  library(igraph)
  library(ggraph)
  library(ggforce)
  library(svglite)
  library(fgsea)
  library(grid)
})

#-----------------------------------------------------------------#
# 1. Output folders and input
#-----------------------------------------------------------------#

dir.create("figures", showWarnings = FALSE)
dir.create(
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq",
  recursive = TRUE,
  showWarnings = FALSE
)

load("host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5_host_RNAseq_inputs.RData")

stopifnot(exists("deg_rna"))
stopifnot(exists("vst_rna"))

if (!exists("group_cols")) {
  group_cols <- c(
    "non_pCR" = "#E07A73",
    "pCR" = "#5AB49B"
  )
}

if (!all(c("non_pCR", "pCR") %in% names(group_cols))) {
  group_cols <- c(
    "non_pCR" = "#E07A73",
    "pCR" = "#5AB49B"
  )
}

stopifnot("Gene" %in% colnames(deg_rna))
stopifnot("stat_DESeq2" %in% colnames(deg_rna))

#-----------------------------------------------------------------#
# 2. Helper functions
#-----------------------------------------------------------------#

format_pathway_label <- function(x) {
  x %>%
    stringr::str_remove("^HALLMARK_") %>%
    stringr::str_remove("^REACTOME_") %>%
    stringr::str_remove("^GOBP_") %>%
    stringr::str_remove("^GOCC_") %>%
    stringr::str_remove("^GOMF_") %>%
    stringr::str_remove("^WP_") %>%
    stringr::str_remove("^KEGG_") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_lower() %>%
    stringr::str_to_sentence() %>%
    stringr::str_replace_all("\\bNfkb\\b", "NF-kB") %>%
    stringr::str_replace_all("\\bTnf\\b", "TNF") %>%
    stringr::str_replace_all("\\bIfn\\b", "IFN") %>%
    stringr::str_replace_all("\\bIl6\\b", "IL-6") %>%
    stringr::str_replace_all("\\bJak\\b", "JAK") %>%
    stringr::str_replace_all("\\bStat\\b", "STAT") %>%
    stringr::str_replace_all("\\bTlr\\b", "TLR") %>%
    stringr::str_replace_all("\\bNod\\b", "NOD") %>%
    stringr::str_replace_all("\\bAhr\\b", "AhR") %>%
    stringr::str_replace_all("\\bNad\\b", "NAD") %>%
    stringr::str_replace_all("\\bEmt\\b", "EMT") %>%
    stringr::str_replace_all("\\bRos\\b", "ROS") %>%
    stringr::str_replace_all("\\bG2m\\b", "G2/M") %>%
    stringr::str_replace_all("\\bMyc\\b", "MYC") %>%
    stringr::str_squish()
}

assign_biological_axis <- function(x) {
  x_upper <- stringr::str_to_upper(x)
  
  dplyr::case_when(
    stringr::str_detect(
      x_upper,
      "COLORECTAL|ADENOMA|\\bCRC\\b|WNT|BETA.?CATENIN|\\bMYC(?:N|L)?\\b|ONCOGEN"
    ) ~ "Adenoma/CRC epithelial program",
    
    stringr::str_detect(
      x_upper,
      "G2.?M|E2F|CELL CYCLE|MITOTIC|DNA REPLICATION|CHROMOSOME SEGREGATION|CENTROSOME|SPINDLE|S PHASE|M PHASE"
    ) ~ "Cell-cycle checkpoint and proliferation",
    
    stringr::str_detect(
      x_upper,
      "EPITHELIAL.?MESENCHYMAL|\\bEMT\\b|ECM|COLLAGEN|INTEGRIN|EXTRACELLULAR MATRIX|MATRIX|MATRISOME|FIBROSIS|FIBROTIC|MIGRATION|INVASION|ADHESION"
    ) ~ "Epithelial–mesenchymal and matrix remodeling",
    
    stringr::str_detect(
      x_upper,
      "PROSTAGLANDIN|PROSTANOID|\\bPGE.?2\\b|\\bPTGS1\\b|\\bPTGS2\\b"
    ) ~ "Prostaglandin signaling and response",

    stringr::str_detect(
      x_upper,
      "OXIDATIVE|\\bROS\\b|HYPOXIA|REACTIVE OXYGEN|NRF2|STRESS RESPONSE"
    ) ~ "Hypoxia and oxidative-stress response",
    
    stringr::str_detect(
      x_upper,
      paste0(
        "T CELL|LYMPH|LEUKOCYTE|ANTIGEN|INTERFERON|\\bIFN\\b|",
        "IMMUNE|IMMUN|CYTOKINE|NF.?KB|\\bTNF[A-Z0-9]*\\b|IL.?6|JAK|",
        "\\bSTAT[0-9AB]*\\b|\\bTLR[0-9]*\\b|\\bNOD[0-9]*\\b|NOD.?LIKE"
      )
    ) ~ "Immune receptor and cytokine signaling",
    
    stringr::str_detect(
      x_upper,
      "MUCUS|GOBLET|MUCIN|BARRIER|TIGHT JUNCTION|APICAL|ANTIMICROBIAL|SECRETORY"
    ) ~ "Barrier and mucus program",
    
    stringr::str_detect(
      x_upper,
      paste0(
        "\\bAHR\\b|INDOLE|TRYPTOPHAN|\\bNAD\\b|\\bNADH\\b|",
        "\\bNADP\\b|\\bNADPH\\b|NIACIN|NICOTIN|LIPID|",
        "PHOSPHOLIPID|GLYCOPROTEIN|GLYCAN|METABOL"
      )
    ) ~ "Lipid/glycoprotein and metabolic remodeling",
    
    stringr::str_detect(
      x_upper,
      "RIBOSOM|TRANSLATION|TRANSLATIONAL|RNA PROCESSING|MRNA|PROTEIN SYNTHESIS"
    ) ~ "Ribosome and translation initiation",
    
    stringr::str_detect(
      x_upper,
      "GOLGI|VESICLE|ENDOSOME|MEMBRANE TRAFFICKING|COATED VESICLE|SECRETION"
    ) ~ "Golgi–endosomal vesicle trafficking",
    
    TRUE ~ "Other"
  )
}

normalize_leading_edge <- function(x) {
  if (is.list(x)) {
    return(
      lapply(x, function(z) {
        z <- as.character(unlist(z))
        z[!is.na(z) & z != ""]
      })
    )
  }
  
  if (is.character(x)) {
    return(
      lapply(x, function(z) {
        if (is.na(z) || z == "") {
          character(0)
        } else if (stringr::str_detect(z, ";")) {
          stringr::str_split(z, ";")[[1]] %>%
            stringr::str_trim() %>%
            .[. != ""]
        } else {
          z
        }
      })
    )
  }
  
  stop("Unsupported leadingEdge type.")
}

#-----------------------------------------------------------------#
# 3. MSigDB gene sets
#-----------------------------------------------------------------#

msig <- msigdbr::msigdbr(species = "Homo sapiens")

if (!"Gene" %in% colnames(msig)) {
  gene_col <- intersect(
    c("gene_symbol", "db_gene_symbol", "human_gene_symbol"),
    colnames(msig)
  )[1]
  
  if (is.na(gene_col)) {
    stop("No recognizable gene symbol column was found in msig.")
  }
  
  msig <- msig %>%
    dplyr::rename(Gene = dplyr::all_of(gene_col))
}

if (!"gs_collection" %in% colnames(msig) && "gs_cat" %in% colnames(msig)) {
  msig <- msig %>%
    dplyr::rename(gs_collection = gs_cat)
}

if (!"gs_subcollection" %in% colnames(msig) && "gs_subcat" %in% colnames(msig)) {
  msig <- msig %>%
    dplyr::rename(gs_subcollection = gs_subcat)
}

if (!"gs_description" %in% colnames(msig)) {
  msig$gs_description <- NA_character_
}

if (!"gs_exact_source" %in% colnames(msig)) {
  msig$gs_exact_source <- NA_character_
}

if (!"db_version" %in% colnames(msig)) {
  msig$db_version <- NA_character_
}

# Include:
#   H
#   C2:CP
#   C5:GO:BP
#   selected C2:CGP colorectal/adenoma signatures
msig_global_long <- msig %>%
  dplyr::select(
    gs_name,
    Gene,
    gs_collection,
    gs_subcollection,
    gs_description,
    gs_exact_source,
    db_version
  ) %>%
  dplyr::filter(
    !is.na(Gene),
    Gene != "",
    Gene %in% rownames(vst_rna)
  ) %>%
  dplyr::filter(
    gs_collection == "H" |
      (
        gs_collection == "C2" &
          stringr::str_detect(gs_subcollection, "^CP")
      ) |
      (
        gs_collection == "C2" &
          stringr::str_detect(gs_subcollection, "CGP") &
          stringr::str_detect(
            gs_name,
            stringr::regex("COLORECTAL|COLON|RECTAL|ADENOMA|\\bCRC\\b", ignore_case = TRUE)
          )
      ) |
      (
        gs_collection == "C5" &
          gs_subcollection == "GO:BP"
      )
  ) %>%
  dplyr::distinct()

gene_set_size_global <- msig_global_long %>%
  dplyr::group_by(
    gs_name,
    gs_collection,
    gs_subcollection,
    gs_description,
    gs_exact_source,
    db_version
  ) %>%
  dplyr::summarise(
    n_genes_total = dplyr::n_distinct(Gene),
    .groups = "drop"
  )

min_gs_size_global <- 10
max_gs_size_global <- 500

msig_global_use_long <- msig_global_long %>%
  dplyr::inner_join(
    gene_set_size_global %>%
      dplyr::filter(
        n_genes_total >= min_gs_size_global,
        n_genes_total <= max_gs_size_global
      ) %>%
      dplyr::select(gs_name),
    by = "gs_name"
  )

gene_sets_global <- split(
  msig_global_use_long$Gene,
  msig_global_use_long$gs_name
)

gene_sets_global <- lapply(gene_sets_global, unique)

message("Gene sets used for Fig. 5A: ", length(gene_sets_global))

#-----------------------------------------------------------------#
# 4. Preranked GSEA
#-----------------------------------------------------------------#

rank_tbl <- deg_rna %>%
  dplyr::filter(
    !is.na(Gene),
    Gene %in% rownames(vst_rna),
    !is.na(stat_DESeq2)
  ) %>%
  dplyr::group_by(Gene) %>%
  dplyr::arrange(dplyr::desc(abs(stat_DESeq2)), .by_group = TRUE) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(dplyr::desc(stat_DESeq2))

rank_stat <- rank_tbl$stat_DESeq2
names(rank_stat) <- rank_tbl$Gene
rank_stat <- sort(rank_stat, decreasing = TRUE)
rank_stat <- rank_stat[is.finite(rank_stat)]

if (length(rank_stat) < 1000) {
  warning("Fewer than 1000 ranked genes were available. Check gene identifiers.")
}

set.seed(20260703)

fgsea_global_raw <- fgsea::fgsea(
  pathways = gene_sets_global,
  stats = rank_stat,
  minSize = min_gs_size_global,
  maxSize = max_gs_size_global,
  eps = 0
)

fgsea_global_res <- fgsea_global_raw %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    gs_name = pathway,
    display_label = format_pathway_label(gs_name),
    biological_axis = assign_biological_axis(display_label),
    leadingEdge = normalize_leading_edge(leadingEdge),
    leadingEdge_string = vapply(
      leadingEdge,
      paste,
      collapse = ";",
      FUN.VALUE = character(1)
    )
  ) %>%
  dplyr::left_join(
    gene_set_size_global,
    by = "gs_name",
    relationship = "many-to-one"
  ) %>%
  dplyr::arrange(
    padj,
    pval,
    dplyr::desc(abs(NES))
  )

#-----------------------------------------------------------------#
# 5. Direction assignment
#-----------------------------------------------------------------#

# Direction is fixed based on DESeq2 contrast and sanity check:
#   NES > 0 = pCR-high
#   NES < 0 = non-pCR-high
#
# Critical check:
#   Sabates colorectal adenoma up should be non-pCR-high.
#   Sabates colorectal adenoma dn should be pCR-high.

fgsea_global_res <- fgsea_global_res %>%
  dplyr::mutate(
    Direction = dplyr::case_when(
      NES > 0 ~ "pCR_high",
      NES < 0 ~ "non_pCR_high",
      TRUE ~ NA_character_
    )
  )

direction_check <- fgsea_global_res %>%
  dplyr::filter(
    stringr::str_detect(
      display_label,
      stringr::regex(
        "Sabates colorectal adenoma up|Sabates colorectal adenoma dn|G2.?M|hypoxia|oxidative|\\bROS\\b|epithelial.*mesenchymal|EMT",
        ignore_case = TRUE
      )
    )
  ) %>%
  dplyr::select(
    display_label,
    NES,
    padj,
    Direction,
    biological_axis
  ) %>%
  dplyr::arrange(
    padj,
    dplyr::desc(abs(NES))
  )

print(direction_check, n = 100)

adenoma_up_check <- direction_check %>%
  dplyr::filter(
    stringr::str_detect(
      display_label,
      stringr::regex("Sabates colorectal adenoma up", ignore_case = TRUE)
    )
  )

if (nrow(adenoma_up_check) > 0) {
  if (!any(adenoma_up_check$Direction == "non_pCR_high")) {
    stop("Direction check failed: Sabates colorectal adenoma up is not assigned to non_pCR_high.")
  }
} else {
  warning("Sabates colorectal adenoma up was not found in direction_check.")
}

#-----------------------------------------------------------------#
# 6. Build map input from unified fgsea result only
#-----------------------------------------------------------------#

targeted_term_regex <- paste(
  c(
    "adenoma",
    "colorectal",
    "colon",
    "rectal",
    "\\bCRC\\b",
    "WNT",
    "beta.?catenin",
    "\\bMYC(?:N|L)?\\b",
    "G2.?M",
    "E2F",
    "cell cycle",
    "mitotic",
    "DNA replication",
    "epithelial.*mesenchymal",
    "\\bEMT\\b",
    "ECM",
    "collagen",
    "oxidative.*stress",
    "\\bROS\\b",
    "hypoxia",
    "TNF",
    "NF.?KB",
    "IL.?6",
    "JAK",
    "\\bSTAT[0-9AB]*\\b",
    "\\bTLR[0-9]*\\b",
    "\\bNOD[0-9]*\\b",
    "NOD.?LIKE",
    "cytokine",
    "chemokine"
  ),
  collapse = "|"
)

fgsea_map_res <- fgsea_global_res %>%
  dplyr::mutate(
    source_for_map = dplyr::case_when(
      stringr::str_detect(
        display_label,
        stringr::regex(targeted_term_regex, ignore_case = TRUE)
      ) ~ "targeted_keyword",
      TRUE ~ "global"
    ),
    node_label = display_label
  )

readr::write_csv(
  fgsea_map_res %>%
    dplyr::select(-leadingEdge),
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5A_unified_GSEA_results_with_direction.csv"
)

#-----------------------------------------------------------------#
# 7. Node selection
#-----------------------------------------------------------------#

gsea_network_fdr_cutoff <- 0.10
gsea_network_abs_nes_cutoff <- 1.30

targeted_fdr_cutoff <- 0.25
targeted_abs_nes_cutoff <- 1.20

gsea_network_top_n_per_direction <- 80

node_tbl <- fgsea_map_res %>%
  dplyr::filter(
    !is.na(padj),
    !is.na(Direction),
    biological_axis != "Other",
    (
      padj < gsea_network_fdr_cutoff &
        abs(NES) >= gsea_network_abs_nes_cutoff
    ) |
      (
        source_for_map == "targeted_keyword" &
          padj < targeted_fdr_cutoff &
          abs(NES) >= targeted_abs_nes_cutoff
      )
  ) %>%
  dplyr::group_by(Direction) %>%
  dplyr::arrange(
    padj,
    dplyr::desc(abs(NES)),
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = gsea_network_top_n_per_direction) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    node_id = gs_name,
    node_size = pmin(-log10(padj + 1e-300), 10),
    evidence_weight = abs(NES) * pmin(-log10(padj + 1e-300), 10),
    node_label = display_label
  )

node_tbl %>%
  dplyr::count(Direction, biological_axis, sort = TRUE) %>%
  print(n = 50)

node_tbl %>%
  dplyr::filter(
    stringr::str_detect(
      display_label,
      stringr::regex(
        "adenoma|colorectal|colon|rectal|\\bCRC\\b|G2|cell cycle|EMT|mesenchymal|\\bMYC(?:N|L)?\\b|WNT|hypoxia|oxidative|\\bROS\\b",
        ignore_case = TRUE
      )
    )
  ) %>%
  dplyr::select(
    display_label,
    NES,
    padj,
    Direction,
    biological_axis,
    source_for_map
  ) %>%
  dplyr::arrange(
    Direction,
    biological_axis,
    padj
  ) %>%
  print(n = 100)

# Audit broad substring rules that can produce biologically unrelated pathway
# annotations. TRUE values identify terms that the earlier rules could have
# captured without an exact biological keyword.
fgsea_global_res %>%
  dplyr::mutate(
    selected_in_map = gs_name %in% node_tbl$node_id,
    bare_ROS_false_positive =
      stringr::str_detect(display_label, stringr::regex("ros", ignore_case = TRUE)) &
      !stringr::str_detect(
        display_label,
        stringr::regex(
          "hypoxia|oxidative|\\bROS\\b|reactive oxygen|NRF2|stress response",
          ignore_case = TRUE
        )
      ),
    generic_TUMOR_rule_risk =
      stringr::str_detect(display_label, stringr::regex("tumor", ignore_case = TRUE)) &
      !stringr::str_detect(
        display_label,
        stringr::regex("colorectal|adenoma|\\bCRC\\b", ignore_case = TRUE)
      ),
    bare_STAT_false_positive =
      stringr::str_detect(display_label, stringr::regex("stat", ignore_case = TRUE)) &
      !stringr::str_detect(
        display_label,
        stringr::regex("\\bSTAT[0-9AB]*\\b", ignore_case = TRUE)
      ),
    bare_NAD_false_positive =
      stringr::str_detect(display_label, stringr::regex("nad", ignore_case = TRUE)) &
      !stringr::str_detect(
        display_label,
        stringr::regex(
          "\\bNAD\\b|\\bNADH\\b|\\bNADP\\b|\\bNADPH\\b",
          ignore_case = TRUE
        )
      ),
    bare_NOD_false_positive =
      stringr::str_detect(display_label, stringr::regex("nod", ignore_case = TRUE)) &
      !stringr::str_detect(
        display_label,
        stringr::regex(
          "\\bNOD[0-9]*\\b|NOD.?LIKE",
          ignore_case = TRUE
        )
      )
  ) %>%
  dplyr::filter(
    dplyr::if_any(
      dplyr::ends_with("_false_positive") | dplyr::ends_with("_rule_risk"),
      identity
    )
  ) %>%
  dplyr::select(
    gs_name,
    display_label,
    NES,
    pval,
    padj,
    Direction,
    biological_axis,
    selected_in_map,
    bare_ROS_false_positive,
    generic_TUMOR_rule_risk,
    bare_STAT_false_positive,
    bare_NAD_false_positive,
    bare_NOD_false_positive
  ) %>%
  dplyr::arrange(dplyr::desc(selected_in_map), padj) %>%
  readr::write_csv(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5A_biological_axis_keyword_audit.csv"
  )

if (nrow(node_tbl) < 5) {
  stop("Too few significant pathways for enrichment map. Relax FDR or NES cutoff.")
}

#-----------------------------------------------------------------#
# 8. Edge table: leading-edge gene overlap
#-----------------------------------------------------------------#

leading_edge_list <- fgsea_map_res$leadingEdge
names(leading_edge_list) <- fgsea_map_res$gs_name

edge_tbl <- utils::combn(node_tbl$node_id, 2, simplify = FALSE) %>%
  purrr::map_dfr(
    function(pair) {
      g1 <- pair[1]
      g2 <- pair[2]
      
      s1 <- unique(leading_edge_list[[g1]])
      s2 <- unique(leading_edge_list[[g2]])
      
      if (length(s1) == 0 || length(s2) == 0) {
        return(NULL)
      }
      
      overlap_n <- length(intersect(s1, s2))
      union_n <- length(union(s1, s2))
      
      if (union_n == 0) {
        return(NULL)
      }
      
      jaccard <- overlap_n / union_n
      overlap_coef <- overlap_n / min(length(s1), length(s2))
      
      tibble::tibble(
        from = g1,
        to = g2,
        overlap_n = overlap_n,
        jaccard = jaccard,
        overlap_coef = overlap_coef
      )
    }
  ) %>%
  dplyr::filter(
    overlap_n >= 3,
    jaccard >= 0.06
  )

if (nrow(edge_tbl) == 0) {
  stop("No edges remained. Lower overlap_n or jaccard cutoff.")
}

#-----------------------------------------------------------------#
# 9. Graph and Louvain clustering
#-----------------------------------------------------------------#

stopifnot(exists("node_tbl"))
stopifnot(exists("edge_tbl"))
stopifnot("node_id" %in% colnames(node_tbl))
stopifnot(all(c("from", "to") %in% colnames(edge_tbl)))

vertices_tbl <- node_tbl %>%
  dplyr::mutate(name = node_id) %>%
  dplyr::select(
    name,
    dplyr::everything(),
    -node_id
  )

missing_vertices <- setdiff(
  unique(c(edge_tbl$from, edge_tbl$to)),
  vertices_tbl$name
)

if (length(missing_vertices) > 0) {
  stop(
    "Some edge vertices are missing from vertices_tbl: ",
    paste(missing_vertices, collapse = ", ")
  )
}

g <- igraph::graph_from_data_frame(
  d = edge_tbl,
  vertices = vertices_tbl,
  directed = FALSE
)

cluster_membership <- igraph::cluster_louvain(
  g,
  weights = igraph::E(g)$jaccard
)$membership

cluster_tbl <- tibble::tibble(
  node_id = igraph::V(g)$name,
  network_cluster = as.integer(cluster_membership)
)

node_tbl <- node_tbl %>%
  dplyr::select(
    -dplyr::any_of(c(
      "network_cluster",
      "cluster_label",
      "dominant_direction",
      "cluster_direction",
      "cluster_direction_majority80",
      "n_nodes",
      "n_pCR",
      "n_non_pCR",
      "n_hypoxia_terms",
      "n_oxidative_terms",
      "pCR_prop",
      "pCR_weighted_prop",
      "median_fdr"
    ))
  ) %>%
  dplyr::left_join(
    cluster_tbl,
    by = "node_id"
  )

node_tbl %>%
  dplyr::count(network_cluster, Direction, sort = TRUE) %>%
  print(n = 100)

#-----------------------------------------------------------------#
# 10. Assign one label per Louvain cluster
#-----------------------------------------------------------------#
refine_cluster_label <- function(term_vec, axis_vec, direction_vec) {
  txt <- paste(term_vec, collapse = " | ")
  txt_upper <- stringr::str_to_upper(txt)
  hypoxia_term_n <- sum(
    stringr::str_detect(stringr::str_to_upper(term_vec), "HYPOXIA")
  )
  oxidative_term_n <- sum(
    stringr::str_detect(
      stringr::str_to_upper(term_vec),
      "OXIDATIVE|\\bROS\\b|REACTIVE OXYGEN"
    ) &
      !stringr::str_detect(stringr::str_to_upper(term_vec), "HYPOXIA")
  )
  prostaglandin_term_n <- sum(
    stringr::str_detect(
      stringr::str_to_upper(term_vec),
      "PROSTAGLANDIN|PROSTANOID|\\bPGE.?2\\b|\\bPTGS1\\b|\\bPTGS2\\b"
    )
  )
  
  #-----------------------------#
  # 1. MYC-related CRC signatures
  # Must be before generic CRC labels.
  #-----------------------------#
  
  has_myc_up <- stringr::str_detect(
    txt_upper,
    "\\bMYC\\b.*\\bUP\\b|\\bUP\\b.*\\bMYC\\b|MYC TARGET|HALLMARK MYC"
  )
  
  has_myc_dn <- stringr::str_detect(
    txt_upper,
    "\\bMYC\\b.*\\bDN\\b|\\bDN\\b.*\\bMYC\\b|MYC.*DOWN|DOWN.*MYC"
  )
  
  if (has_myc_up && !has_myc_dn) {
    return("MYC target\nprogram")
  }
  
  if (has_myc_dn && !has_myc_up) {
    return("MYC-down\n signature")
  }
  
  if (has_myc_up && has_myc_dn) {
    return("Mixed MYC\nCRC signatures")
  }
  
  #-----------------------------#
  # 2. EMT / ECM remodeling
  # Keep these before generic adenoma/CRC labels.
  #-----------------------------#
  
  if (stringr::str_detect(txt_upper, "EPITHELIAL.?MESENCHYMAL|\\bEMT\\b")) {
    return("Epithelial–mesenchymal\ntransition")
  }
  
  if (
    stringr::str_detect(
      txt_upper,
      "ECM|COLLAGEN|EXTRACELLULAR MATRIX|MATRIX|MATRISOME|FIBROSIS|FIBROTIC"
    )
  ) {
    return("ECM/collagen\nremodeling")
  }
  
  #-----------------------------#
  # 3. Specific CRC / epithelial programs
  #-----------------------------#
  
  if (stringr::str_detect(txt_upper, "ADENOMA.*UP|COLORECTAL ADENOMA.*UP|SABATES COLORECTAL ADENOMA UP")) {
    return("Colorectal adenoma\nsignature")
  }
  
  if (stringr::str_detect(txt_upper, "ADENOMA.*DN|COLORECTAL ADENOMA.*DN|SABATES COLORECTAL ADENOMA DN")) {
    return("Colorectal adenoma-down\nsignature")
  }
  
  if (stringr::str_detect(txt_upper, "ADENOMA")) {
    return("Colorectal adenoma\nsignature")
  }
  
  if (stringr::str_detect(txt_upper, "WNT|BETA.?CATENIN")) {
    return("Wnt/beta-catenin\nCRC program")
  }
  
  if (stringr::str_detect(txt_upper, "\\bMYC\\b")) {
    return("MYC-associated\nCRC signature")
  }
  
  #-----------------------------#
  # 4. Cell-cycle programs
  #-----------------------------#
  
  if (stringr::str_detect(txt_upper, "G2.?M")) {
    return("G2/M checkpoint")
  }
  
  if (stringr::str_detect(txt_upper, "E2F")) {
    return("E2F target\nprogram")
  }
  
  if (
    stringr::str_detect(
      txt_upper,
      "CELL CYCLE|MITOTIC|DNA REPLICATION|PROLIFERATION|CENTROSOME|SPINDLE"
    )
  ) {
    return("Cell-cycle\nprogression")
  }

  #-----------------------------#
  # 5. Hypoxia, oxidative-stress, and prostaglandin families
  # Label by the largest number of actual pathway terms, rather than by the
  # first keyword encountered. This prevents one prostaglandin term from
  # overriding a hypoxia-dominant community (and vice versa).
  #-----------------------------#

  response_family_counts <- c(
    hypoxia = hypoxia_term_n,
    oxidative = oxidative_term_n,
    prostaglandin = prostaglandin_term_n
  )

  if (max(response_family_counts) > 0) {
    dominant_response_families <- names(response_family_counts)[
      response_family_counts == max(response_family_counts)
    ]

    if (length(dominant_response_families) == 1) {
      return(
        dplyr::case_when(
          dominant_response_families == "hypoxia" ~ "Hypoxia\nresponse",
          dominant_response_families == "oxidative" ~ "Oxidative-stress /\nROS response",
          dominant_response_families == "prostaglandin" ~ "Prostaglandin\nresponse"
        )
      )
    }

    if (setequal(dominant_response_families, c("hypoxia", "oxidative"))) {
      return("Hypoxia / redox-\nrelated programs")
    }

    if (setequal(dominant_response_families, c("hypoxia", "prostaglandin"))) {
      return("Hypoxia / prostaglandin\nresponse")
    }

    if (setequal(dominant_response_families, c("oxidative", "prostaglandin"))) {
      return("Prostaglandin / redox\nresponse")
    }

    return("Mixed stress-response\nprograms")
  }
  
  #-----------------------------#
  # 6. Immune subprograms
  # More specific immune labels before generic immune labels.
  #-----------------------------#
  
  if (
    stringr::str_detect(
      txt_upper,
      "\\bTLR[0-9]*\\b|\\bNOD[0-9]*\\b|NOD.?LIKE|PATTERN RECOGNITION"
    )
  ) {
    return("TLR/NOD pattern-\nrecognition signaling")
  }
  
  if (
    stringr::str_detect(
      txt_upper,
      "\\bTNF[A-Z0-9]*\\b|NF.?KB|IL.?6|JAK|\\bSTAT[0-9AB]*\\b"
    )
  ) {
    return("TNF–NFkB /\nIL-6–JAK–STAT signaling")
  }
  
  if (stringr::str_detect(txt_upper, "CHEMOKINE")) {
    return("Chemokine\nsignaling")
  }
  
  if (stringr::str_detect(txt_upper, "CYTOKINE|INTERLEUKIN")) {
    return("Cytokine\nsignaling")
  }
  
  if (stringr::str_detect(txt_upper, "T CELL|T-CELL|LYMPHOCYTE|LYMPHOCYTE-MEDIATED|T-CELL ACTIVATION")) {
    return("T-cell / lymphocyte\nactivation")
  }
  
  if (stringr::str_detect(txt_upper, "INTERFERON|\\bIFN\\b")) {
    return("Interferon\nresponse")
  }
  
  if (stringr::str_detect(txt_upper, "IMMUNE|IMMUN|ANTIGEN")) {
    return("Immune receptor\nsignaling")
  }
  
  #-----------------------------#
  # 8. Barrier / mucus
  #-----------------------------#
  
  if (stringr::str_detect(txt_upper, "MUCUS|GOBLET|MUCIN")) {
    return("Mucus/goblet-cell\nprogram")
  }
  
  if (stringr::str_detect(txt_upper, "BARRIER|TIGHT JUNCTION|APICAL")) {
    return("Barrier / apical\njunction program")
  }
  
  #-----------------------------#
  # 9. Metabolism
  #-----------------------------#
  
  if (
    stringr::str_detect(
      txt_upper,
      "\\bNAD\\b|\\bNADH\\b|\\bNADP\\b|\\bNADPH\\b|NIACIN|NICOTIN"
    )
  ) {
    return("NAD/niacin\nmetabolism")
  }
  
  if (stringr::str_detect(txt_upper, "\\bAHR\\b|INDOLE")) {
    return("AhR/indole-ligand\nresponse")
  }
  
  if (stringr::str_detect(txt_upper, "LIPID|PHOSPHOLIPID|GLYCOPROTEIN|GLYCAN|METABOL")) {
    return("Lipid/glycoprotein\nmetabolic remodeling")
  }
  
  #-----------------------------#
  # 10. Translation / trafficking
  #-----------------------------#
  
  if (stringr::str_detect(txt_upper, "RIBOSOM")) {
    return("Ribosome-associated\ntranslation")
  }
  
  if (stringr::str_detect(txt_upper, "TRANSLATION|TRANSLATIONAL")) {
    return("Translation\ninitiation")
  }
  
  if (stringr::str_detect(txt_upper, "GOLGI|VESICLE|ENDOSOME|TRAFFICKING")) {
    return("Golgi–endosomal\nvesicle trafficking")
  }
  
  #-----------------------------#
  # 11. Generic CRC label
  # Keep this late.
  #-----------------------------#
  
  if (stringr::str_detect(txt_upper, "COLORECTAL|COLON|RECTAL|\\bCRC\\b")) {
    return("CRC epithelial\nprogram")
  }
  
  #-----------------------------#
  # 12. Fallback
  #-----------------------------#
  
  top_axis <- names(sort(table(axis_vec), decreasing = TRUE))[1]
  
  dplyr::case_when(
    top_axis == "Adenoma/CRC epithelial program" ~ "CRC epithelial\nprogram",
    top_axis == "Cell-cycle checkpoint and proliferation" ~ "Cell-cycle\nprogram",
    top_axis == "Epithelial–mesenchymal and matrix remodeling" ~ "Epithelial–mesenchymal\ntransition",
    top_axis == "Prostaglandin signaling and response" ~ "Prostaglandin\nresponse",
    top_axis == "Hypoxia and oxidative-stress response" ~ "Hypoxia / oxidative-\nstress programs",
    top_axis == "Immune receptor and cytokine signaling" ~ "Immune receptor\nsignaling",
    top_axis == "Barrier and mucus program" ~ "Barrier / mucus\nprogram",
    top_axis == "Lipid/glycoprotein and metabolic remodeling" ~ "Metabolic\nremodeling",
    top_axis == "Ribosome and translation initiation" ~ "Ribosome-associated\ntranslation",
    top_axis == "Golgi–endosomal vesicle trafficking" ~ "Golgi–endosomal\nvesicle trafficking",
    TRUE ~ top_axis
  )
}

#-----------------------------------------------------------------#
# 11. Re-assign refined cluster labels
#-----------------------------------------------------------------#

cluster_label_tbl <- node_tbl %>%
  dplyr::group_by(network_cluster) %>%
  dplyr::summarise(
    n_nodes = dplyr::n(),
    n_pCR = sum(Direction == "pCR_high"),
    n_non_pCR = sum(Direction == "non_pCR_high"),
    n_hypoxia_terms = sum(
      stringr::str_detect(
        stringr::str_to_upper(display_label),
        "HYPOXIA"
      )
    ),
    n_oxidative_terms = sum(
      stringr::str_detect(
        stringr::str_to_upper(display_label),
        "OXIDATIVE|\\bROS\\b|REACTIVE OXYGEN"
      ) &
        !stringr::str_detect(
          stringr::str_to_upper(display_label),
          "HYPOXIA"
        )
    ),
    n_prostaglandin_terms = sum(
      stringr::str_detect(
        stringr::str_to_upper(display_label),
        "PROSTAGLANDIN|PROSTANOID|\\bPGE.?2\\b|\\bPTGS1\\b|\\bPTGS2\\b"
      )
    ),
    pCR_prop = mean(Direction == "pCR_high"),
    pCR_weighted_prop = sum(
      evidence_weight[Direction == "pCR_high"]
    ) / sum(evidence_weight),
    median_fdr = median(padj, na.rm = TRUE),
    cluster_label = refine_cluster_label(
      term_vec = display_label,
      axis_vec = biological_axis,
      direction_vec = Direction
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    cluster_direction = dplyr::case_when(
      pCR_prop >= 0.70 & pCR_weighted_prop >= 0.70 ~ "pCR_high",
      pCR_prop <= 0.30 & pCR_weighted_prop <= 0.30 ~ "non_pCR_high",
      TRUE ~ "mixed"
    ),
    cluster_direction_majority80 = dplyr::case_when(
      pCR_prop >= 0.80 ~ "pCR_high",
      pCR_prop <= 0.20 ~ "non_pCR_high",
      TRUE ~ "mixed"
    )
  )

node_tbl <- node_tbl %>%
  dplyr::select(
    -dplyr::any_of(c(
      "cluster_label",
      "dominant_direction",
      "cluster_direction",
      "cluster_direction_majority80",
      "n_nodes",
      "n_pCR",
      "n_non_pCR",
      "n_hypoxia_terms",
      "n_oxidative_terms",
      "n_prostaglandin_terms",
      "pCR_prop",
      "pCR_weighted_prop",
      "median_fdr"
    ))
  ) %>%
  dplyr::left_join(
    cluster_label_tbl,
    by = "network_cluster"
  )

cluster_label_tbl %>%
  dplyr::arrange(network_cluster) %>%
  print(n = 100)


#-----------------------------------------------------------------#
# 12. Write node / edge / cluster label CV
#-----------------------------------------------------------------#

readr::write_csv(
  node_tbl,
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5A_enrichment_map_nodes_mixed_direction.csv"
)

readr::write_csv(
  edge_tbl,
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5A_enrichment_map_edges_mixed_direction.csv"
)

readr::write_csv(
  cluster_label_tbl,
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5A_enrichment_map_cluster_labels_mixed_direction.csv"
)

# Flat, list-column-free membership table for checking every pathway assigned
# to each Louvain cluster.
node_tbl %>%
  dplyr::select(
    network_cluster,
    cluster_label,
    cluster_direction,
    cluster_direction_majority80,
    n_nodes,
    n_pCR,
    n_non_pCR,
    n_hypoxia_terms,
    n_oxidative_terms,
    n_prostaglandin_terms,
    pCR_prop,
    pCR_weighted_prop,
    node_id,
    display_label,
    biological_axis,
    NES,
    pval,
    padj,
    Direction,
    n_genes_total,
    leadingEdge_string
  ) %>%
  dplyr::arrange(network_cluster, padj, dplyr::desc(abs(NES))) %>%
  readr::write_csv(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5A_cluster_membership_full.csv"
  )

# Diagnostic record of terms that the previous unbounded "ROS" expression
# would have captured by substring alone. These terms are not ROS pathways.
fgsea_global_res %>%
  dplyr::filter(
    stringr::str_detect(
      display_label,
      stringr::regex("ros", ignore_case = TRUE)
    ),
    !stringr::str_detect(
      display_label,
      stringr::regex(
        "hypoxia|oxidative|\\bROS\\b|reactive oxygen|NRF2|stress response",
        ignore_case = TRUE
      )
    )
  ) %>%
  dplyr::mutate(
    selected_after_fix = gs_name %in% node_tbl$node_id,
    audit_reason = paste0(
      "Previously matched by the substring ROS; now excluded unless another ",
      "explicit biological-axis rule applies."
    )
  ) %>%
  dplyr::select(
    gs_name,
    display_label,
    NES,
    pval,
    padj,
    Direction,
    selected_after_fix,
    audit_reason
  ) %>%
  dplyr::arrange(padj, dplyr::desc(abs(NES))) %>%
  readr::write_csv(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5A_previous_ROS_substring_false_positives.csv"
  )

# Direct audit of the prostaglandin-driven community that was previously
# mislabeled as hypoxia/oxidative stress.
fgsea_global_res %>%
  dplyr::filter(
    stringr::str_detect(
      display_label,
      stringr::regex(
        "prostaglandin|prostanoid|\\bPGE.?2\\b|\\bPTGS1\\b|\\bPTGS2\\b",
        ignore_case = TRUE
      )
    )
  ) %>%
  dplyr::mutate(
    selected_in_map = gs_name %in% node_tbl$node_id,
    leading_edge_n = lengths(leadingEdge)
  ) %>%
  dplyr::select(
    gs_name,
    display_label,
    gs_collection,
    gs_subcollection,
    NES,
    pval,
    padj,
    Direction,
    n_genes_total,
    leading_edge_n,
    selected_in_map,
    leadingEdge_string
  ) %>%
  dplyr::arrange(padj, dplyr::desc(abs(NES))) %>%
  readr::write_csv(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5A_all_prostaglandin_GSEA_results.csv"
  )

prostaglandin_cluster_ids <- node_tbl %>%
  dplyr::filter(
    stringr::str_detect(
      display_label,
      stringr::regex(
        "prostaglandin|prostanoid|\\bPGE.?2\\b|\\bPTGS1\\b|\\bPTGS2\\b",
        ignore_case = TRUE
      )
    )
  ) %>%
  dplyr::distinct(network_cluster) %>%
  dplyr::pull(network_cluster)

node_tbl %>%
  dplyr::filter(network_cluster %in% prostaglandin_cluster_ids) %>%
  dplyr::select(
    network_cluster,
    cluster_label,
    cluster_direction,
    cluster_direction_majority80,
    n_nodes,
    n_pCR,
    n_non_pCR,
    n_hypoxia_terms,
    n_oxidative_terms,
    n_prostaglandin_terms,
    pCR_prop,
    pCR_weighted_prop,
    node_id,
    display_label,
    biological_axis,
    NES,
    pval,
    padj,
    Direction,
    n_genes_total,
    leadingEdge_string
  ) %>%
  dplyr::arrange(network_cluster, padj, dplyr::desc(abs(NES))) %>%
  readr::write_csv(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5A_prostaglandin_cluster_pathways.csv"
  )

#-----------------------------------------------------------------#
# 12A. Audit all hypoxia / oxidative-stress GSEA results
#-----------------------------------------------------------------#

fgsea_global_res %>%
  dplyr::filter(
    stringr::str_detect(
      display_label,
      stringr::regex(
        "hypoxia|oxidative|\\bROS\\b|reactive oxygen",
        ignore_case = TRUE
      )
    )
  ) %>%
  dplyr::mutate(
    selected_in_map = gs_name %in% node_tbl$node_id
  ) %>%
  dplyr::left_join(
    node_tbl %>%
      dplyr::select(
        node_id,
        network_cluster,
        cluster_label,
        cluster_direction,
        cluster_direction_majority80
      ),
    by = c("gs_name" = "node_id")
  ) %>%
  dplyr::select(
    gs_name,
    display_label,
    NES,
    pval,
    padj,
    Direction,
    selected_in_map,
    network_cluster,
    cluster_label,
    cluster_direction,
    cluster_direction_majority80,
    n_genes_total,
    leadingEdge_string
  ) %>%
  dplyr::arrange(padj, dplyr::desc(abs(NES))) %>%
  readr::write_csv(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5A_all_hypoxia_oxidative_GSEA_results.csv"
  )

# Export the pre-existing objects that are subsequently passed to Fig. 5B.
# This permits an exact check of whether panel B used the same gene-set names
# and NES directions as the newly calculated fgsea_global_res used in panel A.
if (exists("fgsea_res")) {
  fig5b_fgsea_audit <- tibble::as_tibble(fgsea_res)
  fig5b_fgsea_text_cols <- names(fig5b_fgsea_audit)[
    vapply(fig5b_fgsea_audit, is.character, logical(1))
  ]

  if (length(fig5b_fgsea_text_cols) > 0) {
    fig5b_fgsea_keep <- Reduce(
      `|`,
      lapply(
        fig5b_fgsea_audit[fig5b_fgsea_text_cols],
        function(x) {
          stringr::str_detect(
            x,
            stringr::regex(
              "hypoxia|oxidative|\\bROS\\b|reactive oxygen",
              ignore_case = TRUE
            )
          )
        }
      )
    )

    fig5b_fgsea_keep <- tidyr::replace_na(fig5b_fgsea_keep, FALSE)
    fig5b_fgsea_audit <- fig5b_fgsea_audit[fig5b_fgsea_keep, , drop = FALSE]

    if ("NES" %in% colnames(fig5b_fgsea_audit)) {
      fig5b_fgsea_audit <- fig5b_fgsea_audit %>%
        dplyr::mutate(
          Direction_recomputed_from_NES = dplyr::case_when(
            NES > 0 ~ "pCR_high",
            NES < 0 ~ "non_pCR_high",
            TRUE ~ NA_character_
          )
        )
    }

    fig5b_fgsea_audit %>%
      dplyr::select(dplyr::where(~ !is.list(.x))) %>%
      readr::write_csv(
        "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5B_fgsea_res_hypoxia_oxidative_audit.csv"
      )
  }
}

if (exists("gsea_plot_df")) {
  fig5b_plot_audit <- tibble::as_tibble(gsea_plot_df)
  fig5b_plot_text_cols <- names(fig5b_plot_audit)[
    vapply(fig5b_plot_audit, is.character, logical(1))
  ]

  if (length(fig5b_plot_text_cols) > 0) {
    fig5b_plot_keep <- Reduce(
      `|`,
      lapply(
        fig5b_plot_audit[fig5b_plot_text_cols],
        function(x) {
          stringr::str_detect(
            x,
            stringr::regex(
              "hypoxia|oxidative|\\bROS\\b|reactive oxygen",
              ignore_case = TRUE
            )
          )
        }
      )
    )

    fig5b_plot_keep <- tidyr::replace_na(fig5b_plot_keep, FALSE)
    fig5b_plot_audit <- fig5b_plot_audit[fig5b_plot_keep, , drop = FALSE]

    if ("NES" %in% colnames(fig5b_plot_audit)) {
      fig5b_plot_audit <- fig5b_plot_audit %>%
        dplyr::mutate(
          Direction_recomputed_from_NES = dplyr::case_when(
            NES > 0 ~ "pCR_high",
            NES < 0 ~ "non_pCR_high",
            TRUE ~ NA_character_
          )
        )
    }

    fig5b_plot_audit %>%
      dplyr::select(dplyr::where(~ !is.list(.x))) %>%
      readr::write_csv(
        "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5B_gsea_plot_df_hypoxia_oxidative_audit.csv"
      )
  }
}

# Trace prostaglandin terms through the pre-existing objects used downstream
# for Panel B. This identifies whether the terms are absent from the tested
# collection (gene_sets_use_long), from fgsea_res, or only from the final
# plotting table (gsea_plot_df).
fig5b_prostaglandin_objects <- list()

if (exists("gene_sets_use_long")) {
  fig5b_prostaglandin_objects$gene_sets_use_long <-
    tibble::as_tibble(gene_sets_use_long)
}

if (exists("fgsea_res")) {
  fig5b_prostaglandin_objects$fgsea_res <- tibble::as_tibble(fgsea_res)
}

if (exists("gsea_plot_df")) {
  fig5b_prostaglandin_objects$gsea_plot_df <-
    tibble::as_tibble(gsea_plot_df)
}

fig5b_prostaglandin_stage_summary <- purrr::imap_dfr(
  fig5b_prostaglandin_objects,
  function(stage_df, stage_name) {
    text_cols <- names(stage_df)[
      vapply(stage_df, is.character, logical(1))
    ]

    row_hit <- rep(FALSE, nrow(stage_df))

    if (length(text_cols) > 0) {
      row_hit <- Reduce(
        `|`,
        lapply(
          stage_df[text_cols],
          function(x) {
            tidyr::replace_na(
              stringr::str_detect(
                x,
                stringr::regex(
                  "prostaglandin|prostanoid|\\bPGE.?2\\b|\\bPTGS1\\b|\\bPTGS2\\b",
                  ignore_case = TRUE
                )
              ),
              FALSE
            )
          }
        )
      )
    }

    tibble::tibble(
      stage = stage_name,
      n_total_rows = nrow(stage_df),
      n_prostaglandin_matching_rows = sum(row_hit)
    )
  }
)

fig5b_prostaglandin_stage_summary %>%
  readr::write_csv(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5B_prostaglandin_stage_summary.csv"
  )

fig5b_prostaglandin_match_details <- tibble::tibble(
  stage = character(),
  row_number = integer(),
  column = character(),
  matched_value = character()
)

for (stage_name in names(fig5b_prostaglandin_objects)) {
  stage_df <- fig5b_prostaglandin_objects[[stage_name]]
  text_cols <- names(stage_df)[
    vapply(stage_df, is.character, logical(1))
  ]

  for (column_name in text_cols) {
    hit <- tidyr::replace_na(
      stringr::str_detect(
        stage_df[[column_name]],
        stringr::regex(
          "prostaglandin|prostanoid|\\bPGE.?2\\b|\\bPTGS1\\b|\\bPTGS2\\b",
          ignore_case = TRUE
        )
      ),
      FALSE
    )

    if (any(hit)) {
      fig5b_prostaglandin_match_details <- dplyr::bind_rows(
        fig5b_prostaglandin_match_details,
        tibble::tibble(
          stage = stage_name,
          row_number = which(hit),
          column = column_name,
          matched_value = as.character(stage_df[[column_name]][hit])
        )
      )
    }
  }
}

fig5b_prostaglandin_match_details %>%
  readr::write_csv(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5B_prostaglandin_stage_matches.csv"
  )

#-----------------------------------------------------------------#
# 12B. Audit the Louvain community containing both term families
#-----------------------------------------------------------------#

hypoxia_oxidative_cluster_ids <- node_tbl %>%
  dplyr::group_by(network_cluster) %>%
  dplyr::summarise(
    has_hypoxia = any(
      stringr::str_detect(
        display_label,
        stringr::regex("hypoxia", ignore_case = TRUE)
      )
    ),
    has_oxidative = any(
      stringr::str_detect(
        display_label,
        stringr::regex(
          "oxidative|\\bROS\\b|reactive oxygen",
          ignore_case = TRUE
        )
      )
    ),
    .groups = "drop"
  ) %>%
  dplyr::filter(has_hypoxia, has_oxidative) %>%
  dplyr::pull(network_cluster)

tibble::tibble(
  audit_status = ifelse(
    length(hypoxia_oxidative_cluster_ids) == 0,
    "No Louvain community contains both hypoxia and oxidative-stress/ROS terms.",
    "At least one Louvain community contains both term families."
  ),
  n_matching_clusters = length(hypoxia_oxidative_cluster_ids),
  matching_cluster_ids = paste(hypoxia_oxidative_cluster_ids, collapse = ";")
) %>%
  readr::write_csv(
    "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5A_hypoxia_oxidative_cluster_audit_status.csv"
  )

if (length(hypoxia_oxidative_cluster_ids) == 0) {
  node_tbl %>%
    dplyr::slice(0) %>%
    dplyr::mutate(term_type = character()) %>%
    dplyr::select(
      network_cluster,
      cluster_label,
      cluster_direction,
      cluster_direction_majority80,
      n_nodes,
      n_pCR,
      n_non_pCR,
      n_hypoxia_terms,
      n_oxidative_terms,
      n_prostaglandin_terms,
      pCR_prop,
      pCR_weighted_prop,
      term_type,
      gs_name,
      display_label,
      gs_collection,
      gs_subcollection,
      NES,
      pval,
      padj,
      Direction,
      n_genes_total,
      leadingEdge_string
    ) %>%
    readr::write_csv(
      "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5A_hypoxia_oxidative_cluster_pathways.csv"
    )

  tibble::tibble(
    network_cluster = integer(),
    from = character(),
    to = character(),
    from_label = character(),
    to_label = character(),
    from_term_type = character(),
    to_term_type = character(),
    from_direction = character(),
    to_direction = character(),
    from_NES = double(),
    to_NES = double(),
    same_direction = logical(),
    database_overlap_n = integer(),
    database_jaccard = double(),
    database_overlap_genes = character(),
    analysis_set_overlap_n = integer(),
    analysis_set_jaccard = double(),
    analysis_set_overlap_genes = character(),
    leading_edge_overlap_n = integer(),
    leading_edge_jaccard = double(),
    leading_edge_overlap_genes = character(),
    edge_retained_in_map = logical()
  ) %>%
    readr::write_csv(
      "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5A_hypoxia_vs_oxidative_direct_overlap.csv"
    )

  warning("No Louvain community contained both hypoxia and oxidative-stress/ROS terms.")
} else {
  target_node_tbl <- node_tbl %>%
    dplyr::filter(network_cluster %in% hypoxia_oxidative_cluster_ids) %>%
    dplyr::mutate(
      term_type = dplyr::case_when(
        stringr::str_detect(
          display_label,
          stringr::regex("hypoxia", ignore_case = TRUE)
        ) &
          stringr::str_detect(
            display_label,
            stringr::regex(
              "oxidative|\\bROS\\b|reactive oxygen",
              ignore_case = TRUE
            )
          ) ~ "hypoxia_and_oxidative",
        stringr::str_detect(
          display_label,
          stringr::regex("hypoxia", ignore_case = TRUE)
        ) ~ "hypoxia",
        stringr::str_detect(
          display_label,
          stringr::regex(
            "oxidative|\\bROS\\b|reactive oxygen",
            ignore_case = TRUE
          )
        ) ~ "oxidative_stress_ROS",
        TRUE ~ "bridge_other"
      )
    )

  target_node_tbl %>%
    dplyr::select(
      network_cluster,
      cluster_label,
      cluster_direction,
      cluster_direction_majority80,
      n_nodes,
      n_pCR,
      n_non_pCR,
      n_hypoxia_terms,
      n_oxidative_terms,
      n_prostaglandin_terms,
      pCR_prop,
      pCR_weighted_prop,
      term_type,
      gs_name,
      display_label,
      gs_collection,
      gs_subcollection,
      NES,
      pval,
      padj,
      Direction,
      n_genes_total,
      leadingEdge_string
    ) %>%
    dplyr::arrange(network_cluster, term_type, padj) %>%
    readr::write_csv(
      "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5A_hypoxia_oxidative_cluster_pathways.csv"
    )

  target_leading_edge_long <- target_node_tbl %>%
    dplyr::select(
      gs_name,
      display_label,
      term_type,
      NES,
      padj,
      Direction,
      leadingEdge
    ) %>%
    tidyr::unnest_longer(leadingEdge, values_to = "Gene") %>%
    dplyr::filter(!is.na(Gene), Gene != "") %>%
    dplyr::distinct()

  target_leading_edge_long %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      n_cluster_pathways = dplyr::n_distinct(gs_name),
      pathway_names = paste(sort(unique(gs_name)), collapse = ";"),
      pathway_labels = paste(sort(unique(display_label)), collapse = ";"),
      term_types = paste(sort(unique(term_type)), collapse = ";"),
      directions = paste(sort(unique(Direction)), collapse = ";"),
      shared_by_hypoxia_and_oxidative =
        any(term_type %in% c("hypoxia", "hypoxia_and_oxidative")) &
        any(term_type %in% c("oxidative_stress_ROS", "hypoxia_and_oxidative")),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      rank_tbl %>%
        dplyr::select(Gene, stat_DESeq2),
      by = "Gene"
    ) %>%
    dplyr::arrange(
      dplyr::desc(shared_by_hypoxia_and_oxidative),
      dplyr::desc(n_cluster_pathways),
      dplyr::desc(abs(stat_DESeq2))
    ) %>%
    readr::write_csv(
      "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5A_hypoxia_oxidative_cluster_leading_edge_genes.csv"
    )

  msig %>%
    dplyr::filter(gs_name %in% target_node_tbl$gs_name) %>%
    dplyr::select(gs_name, Gene) %>%
    dplyr::distinct() %>%
    dplyr::left_join(
      target_node_tbl %>%
        dplyr::select(
          gs_name,
          display_label,
          term_type,
          NES,
          padj,
          Direction
        ),
      by = "gs_name"
    ) %>%
    dplyr::left_join(
      target_leading_edge_long %>%
        dplyr::select(gs_name, Gene) %>%
        dplyr::distinct() %>%
        dplyr::mutate(is_leading_edge = TRUE),
      by = c("gs_name", "Gene")
    ) %>%
    dplyr::mutate(
      is_leading_edge = tidyr::replace_na(is_leading_edge, FALSE)
    ) %>%
    dplyr::left_join(
      rank_tbl %>%
        dplyr::select(Gene, stat_DESeq2),
      by = "Gene"
    ) %>%
    dplyr::arrange(term_type, gs_name, dplyr::desc(is_leading_edge), Gene) %>%
    readr::write_csv(
      "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5A_hypoxia_oxidative_cluster_gene_membership.csv"
    )

  target_pair_tbl <- utils::combn(
    target_node_tbl$node_id,
    2,
    simplify = FALSE
  ) %>%
    purrr::map_dfr(
      function(pair) {
        g1 <- pair[1]
        g2 <- pair[2]

        le1 <- unique(leading_edge_list[[g1]])
        le2 <- unique(leading_edge_list[[g2]])
        analysis_set1 <- unique(gene_sets_global[[g1]])
        analysis_set2 <- unique(gene_sets_global[[g2]])
        database_set1 <- unique(msig$Gene[msig$gs_name == g1])
        database_set2 <- unique(msig$Gene[msig$gs_name == g2])

        le_overlap <- sort(intersect(le1, le2))
        analysis_overlap <- sort(intersect(analysis_set1, analysis_set2))
        database_overlap <- sort(intersect(database_set1, database_set2))

        row1 <- target_node_tbl %>% dplyr::filter(node_id == g1)
        row2 <- target_node_tbl %>% dplyr::filter(node_id == g2)

        tibble::tibble(
          network_cluster = row1$network_cluster[1],
          from = g1,
          to = g2,
          from_label = row1$display_label[1],
          to_label = row2$display_label[1],
          from_term_type = row1$term_type[1],
          to_term_type = row2$term_type[1],
          from_direction = row1$Direction[1],
          to_direction = row2$Direction[1],
          from_NES = row1$NES[1],
          to_NES = row2$NES[1],
          same_direction = row1$Direction[1] == row2$Direction[1],
          database_overlap_n = length(database_overlap),
          database_jaccard = length(database_overlap) /
            length(union(database_set1, database_set2)),
          database_overlap_genes = paste(database_overlap, collapse = ";"),
          analysis_set_overlap_n = length(analysis_overlap),
          analysis_set_jaccard = length(analysis_overlap) /
            length(union(analysis_set1, analysis_set2)),
          analysis_set_overlap_genes = paste(analysis_overlap, collapse = ";"),
          leading_edge_overlap_n = length(le_overlap),
          leading_edge_jaccard = length(le_overlap) /
            length(union(le1, le2)),
          leading_edge_overlap_genes = paste(le_overlap, collapse = ";"),
          edge_retained_in_map =
            length(le_overlap) >= 3 &
            length(le_overlap) / length(union(le1, le2)) >= 0.06
        )
      }
    )

  target_pair_tbl %>%
    dplyr::arrange(
      dplyr::desc(edge_retained_in_map),
      dplyr::desc(leading_edge_jaccard),
      dplyr::desc(leading_edge_overlap_n)
    ) %>%
    readr::write_csv(
      "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5A_hypoxia_oxidative_cluster_pairwise_overlap.csv"
    )

  target_pair_tbl %>%
    dplyr::filter(
      (
        from_term_type %in% c("hypoxia", "hypoxia_and_oxidative") &
          to_term_type %in% c("oxidative_stress_ROS", "hypoxia_and_oxidative")
      ) |
        (
          to_term_type %in% c("hypoxia", "hypoxia_and_oxidative") &
            from_term_type %in% c("oxidative_stress_ROS", "hypoxia_and_oxidative")
        )
    ) %>%
    dplyr::arrange(
      dplyr::desc(edge_retained_in_map),
      dplyr::desc(leading_edge_jaccard)
    ) %>%
    readr::write_csv(
      "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5A_hypoxia_vs_oxidative_direct_overlap.csv"
    )

  target_edge_tbl <- edge_tbl %>%
    dplyr::filter(
      from %in% target_node_tbl$node_id,
      to %in% target_node_tbl$node_id
    )

  target_g <- igraph::graph_from_data_frame(
    d = target_edge_tbl,
    vertices = target_node_tbl %>%
      dplyr::mutate(name = node_id) %>%
      dplyr::select(name, dplyr::everything(), -node_id),
    directed = FALSE
  )

  set.seed(20260703)

  p_hypoxia_oxidative_cluster_audit <- ggraph::ggraph(
    target_g,
    layout = "fr"
  ) +
    ggraph::geom_edge_link(
      ggplot2::aes(
        edge_width = jaccard,
        label = paste0("LE n=", overlap_n)
      ),
      color = "grey65",
      alpha = 0.75,
      label_colour = "grey30",
      label_size = 3.0,
      check_overlap = TRUE,
      angle_calc = "along",
      show.legend = FALSE
    ) +
    ggraph::scale_edge_width(range = c(0.5, 2.0)) +
    ggraph::geom_node_point(
      ggplot2::aes(
        size = node_size,
        fill = Direction
      ),
      shape = 21,
      color = "grey20",
      stroke = 0.35
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = display_label),
      repel = TRUE,
      size = 3.3,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "pCR_high" = group_cols[["pCR"]],
        "non_pCR_high" = group_cols[["non_pCR"]]
      ),
      labels = c(
        "pCR_high" = "pCR",
        "non_pCR_high" = "non-pCR"
      )
    ) +
    ggplot2::scale_size_continuous(
      name = "-log10(FDR)",
      range = c(4.0, 10.0),
      limits = c(0, 10)
    ) +
    ggplot2::labs(
      title = "Hypoxia/oxidative-stress community audit",
      subtitle = "Edges show retained leading-edge overlap; labels show the number of shared leading-edge genes",
      fill = "Enriched in"
    ) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        size = 8
      ),
      legend.position = "right"
    )

  print(p_hypoxia_oxidative_cluster_audit)

  ggplot2::ggsave(
    "figures/Fig5A_hypoxia_oxidative_cluster_audit.svg",
    p_hypoxia_oxidative_cluster_audit,
    width = 8.0,
    height = 5.5,
    device = "svg"
  )
}


fgsea_global_res %>%
  dplyr::filter(
    stringr::str_detect(
      display_label,
      stringr::regex(
        "G2.?M|hypoxia|oxidative|\\bROS\\b|adenoma up|epithelial.*mesenchymal|EMT",
        ignore_case = TRUE
      )
    )
  ) %>%
  dplyr::select(
    display_label,
    NES,
    padj,
    Direction,
    biological_axis
  ) %>%
  dplyr::arrange(padj, dplyr::desc(abs(NES))) %>%
  print(n = 100)


#-----------------------------------------------------------------#
# 13. Plot enrichment map
#-----------------------------------------------------------------#

vertices_tbl <- node_tbl %>%
  dplyr::mutate(name = node_id) %>%
  dplyr::select(
    name,
    dplyr::everything(),
    -node_id
  )

g <- igraph::graph_from_data_frame(
  d = edge_tbl,
  vertices = vertices_tbl,
  directed = FALSE
)

set.seed(20260703)

lay <- ggraph::create_layout(g, layout = "fr")
lay_df <- as.data.frame(lay)

needed_cols <- c(
  "network_cluster",
  "cluster_label",
  "cluster_direction",
  "cluster_direction_majority80",
  "Direction",
  "node_size",
  "biological_axis"
)

missing_cols <- setdiff(needed_cols, colnames(lay_df))

if (length(missing_cols) > 0) {
  stop(
    "Missing columns in layout: ",
    paste(missing_cols, collapse = ", ")
  )
}

cluster_summary <- lay_df %>%
  dplyr::filter(
    !is.na(network_cluster),
    !is.na(cluster_label),
    !is.na(cluster_direction)
  ) %>%
  dplyr::group_by(
    network_cluster,
    cluster_label,
    cluster_direction
  ) %>%
  dplyr::summarise(
    n_nodes = dplyr::n(),
    x = median(x, na.rm = TRUE),
    y = median(y, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_nodes))

cluster_summary_hull <- cluster_summary %>%
  dplyr::filter(n_nodes >= 2) %>%
  dplyr::select(network_cluster)

lay_hull <- lay_df %>%
  dplyr::semi_join(
    cluster_summary_hull,
    by = "network_cluster"
  )

cluster_summary_label <- cluster_summary %>%
  dplyr::filter(n_nodes >= 2)




p_enrichment_map <- ggraph::ggraph(lay) +
  ggforce::geom_mark_hull(
    data = lay_hull,
    ggplot2::aes(
      x = x,
      y = y,
      group = network_cluster,
      fill = cluster_direction
    ),
    inherit.aes = FALSE,
    alpha = 0.12,
    color = "grey60",
    linewidth = 0.35,
    concavity = 4,
    expand = grid::unit(2.5, "mm"),
    radius = grid::unit(2, "mm"),
    show.legend = FALSE
  ) +
  ggraph::geom_edge_link(
    ggplot2::aes(edge_width = jaccard),
    color = "grey70",
    alpha = 0.45,
    show.legend = FALSE
  ) +
  ggraph::scale_edge_width(range = c(0.2, 1.1)) +
  ggraph::geom_node_point(
    ggplot2::aes(
      size = node_size,
      fill = Direction
    ),
    shape = 21,
    color = "grey20",
    stroke = 0.25,
    alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    data = cluster_summary_label,
    ggplot2::aes(
      x = x,
      y = y,
      label = cluster_label
    ),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold",
    color = "black",
    segment.color = "grey45",
    segment.size = 0.25,
    box.padding = 0.35,
    point.padding = 0.2,
    seed = 20260703,
    min.segment.length = 0
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "pCR_high" = group_cols[["pCR"]],
      "non_pCR_high" = group_cols[["non_pCR"]],
      "mixed" = "grey70"
    ),
    labels = c(
      "pCR_high" = "pCR / pCR-dominant",
      "non_pCR_high" = "non-pCR / non-pCR-dominant",
      "mixed" = "Mixed community"
    ),
    breaks = c("pCR_high", "non_pCR_high", "mixed"),
    drop = FALSE
  ) +
  ggplot2::scale_size_continuous(
    name = "-log10(FDR)",
    range = c(2.0, 7.2),
    limits = c(0, 10)
  ) +
  ggplot2::labs(
    title = "Host transcriptional programs associated with TNT response",
    fill = "Direction"
  ) +
  ggplot2::theme_void(base_size = 10) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      hjust = 0.5,
      face = "bold",
      size = 12
    ),
    legend.position = "right"
  )

node_tbl %>%
  dplyr::filter(
    stringr::str_detect(
      display_label,
      stringr::regex("\\bMYC\\b", ignore_case = TRUE)
    )
  ) %>%
  dplyr::mutate(
    myc_term_type = dplyr::case_when(
      stringr::str_detect(
        display_label,
        stringr::regex("\\bMYC\\b.*\\bDN\\b|\\bDN\\b.*\\bMYC\\b|MYC.*DOWN|DOWN.*MYC", ignore_case = TRUE)
      ) ~ "MYC-down signature",
      
      stringr::str_detect(
        display_label,
        stringr::regex("\\bMYC\\b.*\\bUP\\b|\\bUP\\b.*\\bMYC\\b|MYC TARGET|HALLMARK MYC", ignore_case = TRUE)
      ) ~ "MYC-up / target signature",
      
      TRUE ~ "MYC-associated, unclear direction"
    )
  ) %>%
  dplyr::select(
    network_cluster,
    cluster_label,
    display_label,
    myc_term_type,
    NES,
    padj,
    Direction,
    biological_axis
  ) %>%
  dplyr::arrange(network_cluster, padj) %>%
  print(n = 100) %>% 
  as.data.frame()

print(p_enrichment_map)

ggplot2::ggsave(
  "figures/Fig5A_host_RNAseq_enrichment_map_mixed_direction.svg",
  p_enrichment_map,
  width = 9.0,
  height = 6.8,
  device = "svg"
)

#-----------------------------------------------------------------#
# 13B. Additional map: simple 80% majority rule for community hulls
#-----------------------------------------------------------------#

cluster_summary_majority80 <- lay_df %>%
  dplyr::filter(
    !is.na(network_cluster),
    !is.na(cluster_label),
    !is.na(cluster_direction_majority80)
  ) %>%
  dplyr::group_by(
    network_cluster,
    cluster_label,
    cluster_direction_majority80
  ) %>%
  dplyr::summarise(
    n_nodes = dplyr::n(),
    x = median(x, na.rm = TRUE),
    y = median(y, na.rm = TRUE),
    .groups = "drop"
  )

lay_hull_majority80 <- lay_df %>%
  dplyr::semi_join(
    cluster_summary_majority80 %>%
      dplyr::filter(n_nodes >= 2) %>%
      dplyr::select(network_cluster),
    by = "network_cluster"
  )

p_enrichment_map_majority80 <- ggraph::ggraph(lay) +
  ggforce::geom_mark_hull(
    data = lay_hull_majority80,
    ggplot2::aes(
      x = x,
      y = y,
      group = network_cluster,
      fill = cluster_direction_majority80
    ),
    inherit.aes = FALSE,
    alpha = 0.16,
    color = "grey55",
    linewidth = 0.40,
    concavity = 4,
    expand = grid::unit(2.5, "mm"),
    radius = grid::unit(2, "mm"),
    show.legend = TRUE
  ) +
  ggraph::geom_edge_link(
    ggplot2::aes(edge_width = jaccard),
    color = "grey70",
    alpha = 0.45,
    show.legend = FALSE
  ) +
  ggraph::scale_edge_width(range = c(0.2, 1.1)) +
  ggraph::geom_node_point(
    ggplot2::aes(
      size = node_size,
      fill = Direction
    ),
    shape = 21,
    color = "grey20",
    stroke = 0.25,
    alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    data = cluster_summary_majority80 %>%
      dplyr::filter(n_nodes >= 2),
    ggplot2::aes(
      x = x,
      y = y,
      label = cluster_label
    ),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold",
    color = "black",
    segment.color = "grey45",
    segment.size = 0.25,
    box.padding = 0.35,
    point.padding = 0.2,
    seed = 20260703,
    min.segment.length = 0
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "pCR_high" = group_cols[["pCR"]],
      "non_pCR_high" = group_cols[["non_pCR"]],
      "mixed" = "grey70"
    ),
    labels = c(
      "pCR_high" = "pCR node / pCR-dominant community",
      "non_pCR_high" = "non-pCR node / non-pCR-dominant community",
      "mixed" = "Mixed community (<80% agreement)"
    ),
    limits = c("pCR_high", "non_pCR_high", "mixed"),
    drop = FALSE
  ) +
  ggplot2::scale_size_continuous(
    name = "-log10(FDR)",
    range = c(2.0, 7.2),
    limits = c(0, 10)
  ) +
  ggplot2::labs(
    title = "Host transcriptional programs associated with TNT response",
    fill = "Node / community direction"
  ) +
  ggplot2::theme_void(base_size = 10) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      hjust = 0.5,
      face = "bold",
      size = 12
    ),
    legend.position = "right"
  )

print(p_enrichment_map_majority80)

ggplot2::ggsave(
  "figures/Fig5A_host_RNAseq_enrichment_map_majority80.svg",
  p_enrichment_map_majority80,
  width = 9.0,
  height = 6.8,
  device = "svg"
)

#-----------------------------------------------------------------#
# 14. Final direction sanity check
#-----------------------------------------------------------------#

node_tbl %>%
  dplyr::filter(
    stringr::str_detect(
      display_label,
      stringr::regex("adenoma|colorectal|colon|rectal|\\bCRC\\b", ignore_case = TRUE)
    )
  ) %>%
  dplyr::select(
    display_label,
    NES,
    padj,
    Direction,
    biological_axis,
    network_cluster,
    cluster_label
  ) %>%
  print(n = 50)


#-----------------------------------------------------------------#
# 15. Save downstream inputs for Fig. 5B and Fig. 5C
#-----------------------------------------------------------------#

dir.create(
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq",
  recursive = TRUE,
  showWarnings = FALSE
)

fig5_downstream_objects <- c(
  "gsea_plot_df",
  "gene_sets_use_long",
  "fgsea_res",
  "vst_int",
  "col_int",
  "group_cols"
)

if (exists("met_int")) {
  fig5_downstream_objects <- c(fig5_downstream_objects, "met_int")
} else {
  warning("met_int was not found. Fig. 5C metabolite-related analyses will require met_int later.")
}

missing_fig5_downstream_objects <- fig5_downstream_objects[
  !vapply(fig5_downstream_objects, exists, logical(1))
]

if (length(missing_fig5_downstream_objects) > 0) {
  stop(
    "Missing required downstream objects: ",
    paste(missing_fig5_downstream_objects, collapse = ", ")
  )
}

save(
  list = fig5_downstream_objects,
  file = "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5B_C_downstream_inputs.RData"
)

message(
  "Saved Fig. 5B/5C downstream inputs: ",
  "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq/Fig5B_C_downstream_inputs.RData"
)
