# Analyzed by Nayun Kim
# CNU Translational Microbiome Lab.

# ==================== 0. Setting ====================

# ---------- 0-1. Directory setting ----------

rm(list = ls())
setwd("D:/2-연구/2-CRC metagenomics/")
options(stringsAsFactors = FALSE)
set.seed(123)

library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggforce)
library(purrr)
library(broom)
library(ggrepel)

# PCoA
library(vegan)

# Heatmap
library(ComplexHeatmap)
library(circlize)

# ---------- 0-2. Rdata ----------
# Rdata
# load("260224 final script/260526 enzyme.RData")
load("input/metabolite_metadata.RData")



##### 1. global magnitude coherence

  # subject-level overall shift magnitude
  # modality-level overall shift coherence



  # microbiome shift magnitude
    # Relative abundance + Bray-Curtis
    # CLR + Aitchison

  # metabolome shift magnitude
    # absolute conc. + Euclidean
    # absolute conc. + Manhattan 
    # log10(x + pseudocount) + Euclidean
    # log10(x + pseudocount) + Manhattan


# 2. feature-level directional coherence = pairwise feature coherence
# 3. response-specific coherence = group-moditified coherence

colnames(m)
head(m)
# ---------- 0-3. Importing ----------
# # Before samples
# before_samples <- colnames(gb[, -1])
valid_samples <- metabolite_before_valid$SampleID

s[1:3, ]
t[1:3, ]
metabolite[1:3, ]

# Enzyme data

ko <- read_tsv("260224 final script/Input/merged_genefamilies_KO_named.tsv") %>% 
  dplyr::rename(GeneFamily_raw = `# Gene Family`) %>%
  
  tidyr::separate(
    GeneFamily_raw,
    into = c("KO", "Taxon"),
    sep = "\\|",
    fill = "right",
    extra = "merge",
    remove = FALSE
  ) %>%
  
  dplyr::mutate(
    Taxon = dplyr::na_if(Taxon, ""),
    
    ## 핵심: KO total row와 stratified row를 먼저 구분
    Feature_level = dplyr::case_when(
      is.na(Taxon) ~ "KO_total",
      TRUE ~ "KO_stratified"
    ),
    
    Genus_raw = stringr::str_extract(Taxon, "(?<=g__)[^\\.]+"),
    Species_raw = stringr::str_extract(Taxon, "(?<=s__)[^\\.]+"),
    
    ## KO total과 taxonomically unclassified를 명확히 분리
    Genus = dplyr::case_when(
      Feature_level == "KO_total" ~ "TOTAL",
      is.na(Genus_raw) ~ "UNCLASSIFIED_GENUS",
      TRUE ~ Genus_raw
    ),
    
    Species = dplyr::case_when(
      Feature_level == "KO_total" ~ "TOTAL",
      is.na(Species_raw) ~ "UNCLASSIFIED_SPECIES",
      TRUE ~ Species_raw
    ),
    
    Taxon_status = dplyr::case_when(
      Feature_level == "KO_total" ~ "unstratified_total",
      is.na(Genus_raw) & is.na(Species_raw) ~ "taxonomically_unclassified",
      is.na(Species_raw) ~ "species_unclassified",
      TRUE ~ "species_classified"
    )
  ) %>%
  
  dplyr::relocate(
    KO, Feature_level, Taxon_status, Genus, Species, Taxon,
    .after = GeneFamily_raw
  ) %>%
  
  dplyr::select(-Genus_raw, -Species_raw); head(ko)
  ko[10000:10003,]








# Enzyme data
ec <- read_tsv("260224 final script/Input/merged_genefamilies_EC_named.tsv") %>% 
  dplyr::rename(GeneFamily_raw = `# Gene Family`) %>%
  
  tidyr::separate(
    GeneFamily_raw,
    into = c("EC", "Taxon"),
    sep = "\\|",
    fill = "right",
    extra = "merge",
    remove = FALSE
  ) %>%
  
  dplyr::mutate(
    Taxon = dplyr::na_if(Taxon, ""),
    
    ## 핵심: KO total row와 stratified row를 먼저 구분
    Feature_level = dplyr::case_when(
      is.na(Taxon) ~ "EC_total",
      TRUE ~ "EC_stratified"
    ),
    
    Genus_raw = stringr::str_extract(Taxon, "(?<=g__)[^\\.]+"),
    Species_raw = stringr::str_extract(Taxon, "(?<=s__)[^\\.]+"),
    
    ## KO total과 taxonomically unclassified를 명확히 분리
    Genus = dplyr::case_when(
      Feature_level == "EC_total" ~ "TOTAL",
      is.na(Genus_raw) ~ "UNCLASSIFIED_GENUS",
      TRUE ~ Genus_raw
    ),
    
    Species = dplyr::case_when(
      Feature_level == "EC_total" ~ "TOTAL",
      is.na(Species_raw) ~ "UNCLASSIFIED_SPECIES",
      TRUE ~ Species_raw
    ),
    
    Taxon_status = dplyr::case_when(
      Feature_level == "EC_total" ~ "unstratified_total",
      is.na(Genus_raw) & is.na(Species_raw) ~ "taxonomically_unclassified",
      is.na(Species_raw) ~ "species_unclassified",
      TRUE ~ "species_classified"
    )
  ) %>%
  
  dplyr::relocate(
    EC, Feature_level, Taxon_status, Genus, Species, Taxon,
    .after = GeneFamily_raw
  ) %>%
  
  dplyr::select(-Genus_raw, -Species_raw)

  head(ec)
  ec[10000:10003,]


ec_valid_tmp <- ec %>% 
  filter(!EC %in% c("UNMAPPED", "UNGROUPED")) %>% 
  select(EC, Genus, Species, all_of(valid_samples)); ec_valid_tmp




# Sorting by prevalence
# 20 before-samples * 20% prevalence = 4
# 20 before-samples * 30% prevalence = 6
# 20 before-samples * 50% prevalence = 10
ec_prev <- ec_valid_tmp %>% 
  filter(Genus == "TOTAL") %>% 
  select(-Genus, -Species) %>% 
  group_by(EC) %>% 
  dplyr::slice_head(n = 1) %>% # find("slice") # "package:IRanges" "package:dplyr"
  ungroup() %>% 
  column_to_rownames("EC") %>% 
  t() %>% 
  as.data.frame(); ec_prev


prevalence <- colSums(ec_prev[, -1] > 0)
valid_ec_cols <- names(prevalence[prevalence >= 10])

ec_valid <- ec_valid_tmp %>% 
  filter(EC %in% valid_ec_cols); ec_valid


ec_valid_wide <- ec_valid %>% 
  filter(is.na(Genus)) %>% 
  select(-Genus, -Species) %>% 
  group_by(EC) %>% 
  dplyr::slice_head(n = 1) %>%
  column_to_rownames("EC") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleID") %>% 
  left_join(metabolite_before_valid %>%
              select(TRG_1, SampleID),
            by = "SampleID") %>% 
  relocate(c(TRG_1), .before = "SampleID"); ec_valid_wide

# Species level data frame
ec_species <- ec %>% 
  select(EC, Genus, Species, all_of(valid_samples)) %>% 
  filter(!is.na(Genus),
         !EC %in% c("UNGROUPED", "UNMAPPED"),
         !grepl("^\\d+(\\.\\d+)*:Deleted entry$", EC)) %>%
  mutate(Name = paste(EC, Species, sep = "|")) %>% 
  relocate(Name, .before = EC); ec_species

ec_s_top3 <- ec_species %>% 
  mutate(Total = rowSums(across(all_of(valid_samples)), na.rm = T)) %>% 
  group_by(EC) %>% 
  slice_max(order_by = Total,
            n = 3,
            with_ties = F) %>% 
  ungroup()

ec_s_valid <- ec_s_top3 %>% 
  filter(EC %in% valid_ec_cols, EC != "1.3.99.1: Deleted entry"); ec_s_valid

ec_s_valid_wide <- ec_s_valid %>% 
  select(-EC, -Genus, -Species) %>% 
  column_to_rownames("Name") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleID") %>% 
  left_join(metabolite_before_valid %>%
              select(TRG_1, SampleID),
            by = "SampleID") %>% 
  relocate(TRG_1, .before = "SampleID"); ec_s_valid_wide[1:6, 1:6]

valid_ec_s_cols <- unique(ec_s_valid$Name)

# Prevalence
s_prev <- s_trans %>% 
  rownames_to_column("SampleID") %>% 
  filter(SampleID %in% metabolite_before_valid$SampleID) %>% 
  pivot_longer(cols = -SampleID,
               names_to = "Species",
               values_to = "Abundance") %>% 
  group_by(Species) %>% 
  summarise(SampleCount = n_distinct(SampleID[Abundance > 0]),
            Prevalence = SampleCount / n_distinct(metabolite_before_valid$SampleID),
            .groups = "drop"); s_prev

# g_prev <- g_trans %>% 
#   rownames_to_column("SampleID") %>% 
#   filter(SampleID %in% metabolite_before_valid$SampleID) %>% 
#   pivot_longer(cols = -SampleID,
#                names_to = "Genus",
#                values_to = "Abundance") %>% 
#   group_by(Genus) %>% 
#   summarise(SampleCount = n_distinct(SampleID[Abundance > 0]),
#             Prevalence = SampleCount / n_distinct(metabolite_before_valid$SampleID),
#             .groups = "drop"); g_prev





# ---------- 0-4. Palette ----------
# 그룹별 색상
group_colors <- c("CR" = "#F7D9BC", "nonCR" = "#80461B")
stage_colors <- c("early" = "#1b9e77", "advanced" = "#d95f02")














# ==================== 1. Differential Enrichment ====================

# ---------- 1-1. Full screening ----------
# Data frame
ec_long <- ec_valid %>% 
  select(-Genus, -Species) %>% 
  group_by(EC) %>%
  slice(1) %>%
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1, Tstage), by = "SampleID"); ec_long

ec_stat <- ec_long %>% 
  group_by(EC) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  mutate(FDR = p.adjust(Wilcox, method = "BH")) %>% 
  relocate(FDR, .after = Wilcox) %>% 
  arrange(Wilcox); ec_stat

# write.csv(ec_stat, file = "Data/260526 Wilcoxon - enzyme.csv", row.names = FALSE)

ec_sig_names <- ec_stat %>% 
  arrange(Wilcox) %>% 
  head(20) %>% 
  pull(EC); ec_sig_names

ec_target_names <- c("2.7.2.1: Acetate kinase",
                     "6.3.4.21: Nicotinate phosphoribosyltransferase",
                     "2.3.1.8: Phosphate acetyltransferase",
                     "3.5.1.19: Nicotinamidase",
                     "4.2.1.54: Lactoyl-CoA dehydratase",
                     "2.8.3.1: Propionate CoA-transferase",
                     "4.1.1.74: Indolepyruvate decarboxylase",
                     "2.8.3.17: Cinnamoyl-CoA:phenyllactate CoA-transferase",
                     "2.7.2.7: Butyrate kinase",
                     "2.3.1.19: Phosphate butyryltransferase",
                     "2.8.3.8: Acetate CoA-transferase")

ec_subtarget_names <- c("6.3.1.5: NAD(+) synthase",
                        "2.4.2.19: Nicotinate-nucleotide diphosphorylase (carboxylating)",
                        "6.3.5.1: NAD(+) synthase (glutamine-hydrolyzing)",
                        "5.4.99.2: Methylmalonyl-CoA mutase",
                        "2.7.7.18: Nicotinate-nucleotide adenylyltransferase",
                        "1.1.1.157: 3-hydroxybutyryl-CoA dehydrogenase",
                        "2.5.1.72: Quinolinate synthase",
                        "1.2.1.4: Aldehyde dehydrogenase (NADP(+))",
                        "2.7.2.15: Propionate kinase",
                        "1.2.1.3: Aldehyde dehydrogenase (NAD(+))",
                        "4.2.1.55: 3-hydroxybutyryl-CoA dehydratase",
                        "3.5.1.4: Amidase",
                        "1.4.3.16: L-aspartate oxidase",
                        "1.2.1.5: Aldehyde dehydrogenase (NAD(P)(+))",
                        "2.6.1.57: Aromatic-amino-acid transaminase",
                        "2.3.1.9: Acetyl-CoA C-acetyltransferase",
                        "1.3.8.1: Short-chain acyl-CoA dehydrogenase",
                        "6.2.1.1: Acetate--CoA ligase",
                        "4.2.1.28: Propanediol dehydratase")

ec_order <- ec_stat %>% 
  filter(EC %in% ec_sig_names) %>% 
  arrange(Wilcox) %>% 
  pull(EC) %>% 
  str_wrap(width = 15); ec_order

ec_sig_long <- ec_valid %>% 
  filter(EC %in% ec_sig_names) %>% 
  group_by(EC) %>% 
  slice(1) %>% 
  select(-Genus, -Species) %>%
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before %>% select(SampleID, TRG_1), by = "SampleID") %>% 
  mutate(EC_wrap = str_wrap(EC, width = 15),
         EC_wrap = factor(EC_wrap, levels = ec_order)); ec_sig_long

p_ec_sig <- ec_sig_long %>% 
  ggplot(aes(x = TRG_1, y = Abundance, fill = TRG_1)) +
  facet_wrap(. ~ EC_wrap, scales = "free_y", nrow = 2, ncol = 10) + 
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CR", "nonCR"))) +
  theme_classic() + 
  theme(axis.title = element_text(size = rel(1.1)), 
        axis.text.x = element_text(size =rel(1.3), color = "black"), 
        legend.position = "bottom") +
  labs(x = "Disease", y = "Abundance", fill = "TRG",
       title = "Significant Enzymes"); p_ec_sig

ggsave("Figure/7-00. Significant enzymes.svg", device = "svg",
       plot = p_ec_sig, height = 12, width = 25)

ec_target_long <- ec_valid %>% 
  filter(EC %in% ec_target_names) %>% 
  group_by(EC) %>% 
  slice(1) %>% 
  select(-Genus, -Species) %>%
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before %>% select(SampleID, TRG_1), by = "SampleID"); ec_target_long

p_ec_target <- ec_target_long %>% 
  ggplot(aes(x = TRG_1, y = Abundance, fill = TRG_1)) +
  facet_wrap(. ~ str_wrap(EC, 15), scales = "free_y", nrow = 2, ncol = 6) + 
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CR", "nonCR"))) +
  theme_classic() + 
  theme(axis.title = element_text(size = rel(1.1)), 
        axis.text.x = element_text(size =rel(1.3), color = "black"), 
        legend.position = "bottom") +
  labs(x = "Disease", y = "Abundance", fill = "TRG",
       title = "Target Enzymes"); p_ec_target

ggsave("Figure/7-00. Target enzymes.svg", device = "svg",
       plot = p_ec_target, height = 12, width = 15)

ec_subtarget_long <- ec_valid %>% 
  filter(EC %in% ec_subtarget_names) %>% 
  group_by(EC) %>% 
  slice(1) %>% 
  select(-Genus, -Species) %>%
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before %>% select(SampleID, TRG_1), by = "SampleID"); ec_subtarget_long

p_ec_subtarget <- ec_subtarget_long %>% 
  ggplot(aes(x = TRG_1, y = Abundance, fill = TRG_1)) +
  facet_wrap(. ~ str_wrap(EC, 15), scales = "free_y", nrow = 2, ncol = 10) + 
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CR", "nonCR"))) +
  theme_classic() + 
  theme(axis.title = element_text(size = rel(1.1)), 
        axis.text.x = element_text(size =rel(1.3), color = "black"), 
        legend.position = "bottom") +
  labs(x = "Disease", y = "Abundance", fill = "TRG",
       title = "Sub-Target Enzymes"); p_ec_subtarget

ggsave("Figure/7-00. Subtarget enzymes.svg", device = "svg",
       plot = p_ec_subtarget, height = 12, width = 25)

# ---------- 1-2. Acetate kinase ----------
# Data frame
# Species: E. rectale
acetate_kinase_long <- ec_valid %>% 
  filter(EC == "2.7.2.1: Acetate kinase", !is.na(Species)) %>% 
  select(-Genus) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); acetate_kinase_long

acetate_kinase_stat <- acetate_kinase_long %>% 
  group_by(Species) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(s_prev %>% select(Species, Prevalence), by = "Species") %>% 
  arrange(Wilcox); acetate_kinase_stat

# Genus
acetate_kinase_long2 <- ec_valid %>% 
  filter(EC == "2.7.2.1: Acetate kinase", !is.na(Genus)) %>% 
  select(-Species) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); acetate_kinase_long2

acetate_kinase_stat2 <- acetate_kinase_long2 %>% 
  group_by(Genus) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(g_prev %>% select(Genus, Prevalence), by = "Genus") %>% 
  arrange(Wilcox); acetate_kinase_stat2

# write.csv(acetate_kinase_stat, file = "Data/260522 Wilcoxon - acetate kinase (species).csv", row.names = FALSE)

p_acetate_e.rectale <- acetate_kinase_long %>% 
  filter(Species == "Eubacterium_rectale") %>% 
  ggplot(aes(x = TRG_1, y = Abundance, fill = TRG_1)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CR", "nonCR"))) +
  theme_classic() + 
  theme(axis.title = element_text(size = rel(1.1)), 
        axis.text.x = element_text(size =rel(1.3), color = "black"), 
        legend.position = "bottom") +
  labs(x = "Disease", y = "Eubacterium rectale", fill = "TRG",
       title = "2.7.2.1: Acetate kinase"); p_acetate_e.rectale

p_acetate_r.inulinivorans <- acetate_kinase_long %>% 
  filter(Species == "Roseburia_inulinivorans") %>% 
  ggplot(aes(x = TRG_1, y = Abundance, fill = TRG_1)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CR", "nonCR"))) +
  theme_classic() + 
  theme(axis.title = element_text(size = rel(1.1)), 
        axis.text.x = element_text(size =rel(1.3), color = "black"), 
        legend.position = "bottom") +
  labs(x = "Disease", y = "Roseburia inulinivorans", fill = "TRG",
       title = "2.7.2.1: Acetate kinase"); p_acetate_r.inulinivorans

p_acetate_d.sp <- acetate_kinase_long %>% 
  filter(Species == "Dorea_sp_OM02_2LB") %>% 
  ggplot(aes(x = TRG_1, y = Abundance, fill = TRG_1)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CR", "nonCR"))) +
  theme_classic() + 
  theme(axis.title = element_text(size = rel(1.1)), 
        axis.text.x = element_text(size =rel(1.3), color = "black"), 
        legend.position = "bottom") +
  labs(x = "Disease", y = "Dorea sp. OM02-2LB", fill = "TRG",
       title = "2.7.2.1: Acetate kinase"); p_acetate_d.sp

p_acetate_f.parusnitzii <- acetate_kinase_long %>% 
  filter(Species == "Faecalibacterium_prausnitzii") %>% 
  ggplot(aes(x = TRG_1, y = Abundance, fill = TRG_1)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CR", "nonCR"))) +
  theme_classic() + 
  theme(axis.title = element_text(size = rel(1.1)), 
        axis.text.x = element_text(size =rel(1.3), color = "black"), 
        legend.position = "bottom") +
  labs(x = "Disease", y = "Faecalibacterium prausnitzii", fill = "TRG",
       title = "2.7.2.1: Acetate kinase"); p_acetate_f.parusnitzii

p_acetate_kinase <- ggarrange(p_acetate_e.rectale, p_acetate_r.inulinivorans, p_acetate_d.sp, p_acetate_f.parusnitzii,
                              nrow = 1, ncol = 4); p_acetate_kinase

ggsave("Figure/7-01. Acetate kinase (species).svg", device = "svg",
       plot = p_acetate_kinase, height = 6, width = 12)

# ---------- 1-3. Nicotinate phosphoribosyltransferase ----------
# Data frame
# Species: D. longicatena, R. bromii, E. retale, F. saccharivorans
nicotinate_transferase_long <- ec_valid %>% 
  filter(EC == "6.3.4.21: Nicotinate phosphoribosyltransferase", !is.na(Species)) %>% 
  select(-Genus) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); nicotinate_transferase_long

nicotinate_transferase_stat <- nicotinate_transferase_long %>% 
  group_by(Species) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(s_prev %>% select(Species, Prevalence), by = "Species") %>% 
  arrange(Wilcox); nicotinate_transferase_stat

# Genus
nicotinate_transferase_long2 <- ec_valid %>% 
  filter(EC == "6.3.4.21: Nicotinate phosphoribosyltransferase", !is.na(Genus)) %>% 
  select(-Species) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); nicotinate_transferase_long2

nicotinate_transferase_stat2 <- nicotinate_transferase_long2 %>% 
  group_by(Genus) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(g_prev %>% select(Genus, Prevalence), by = "Genus") %>% 
  arrange(Wilcox); nicotinate_transferase_stat2

# write.csv(nicotinate_transferase_stat, file = "Data/260522 Wilcoxon - nicotinate phosphoribosyltransferase (species).csv", row.names = FALSE)

p_nicotinate_e.rectale <- nicotinate_transferase_long %>% 
  filter(Species == "Eubacterium_rectale") %>% 
  ggplot(aes(x = TRG_1, y = Abundance, fill = TRG_1)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CR", "nonCR"))) +
  theme_classic() + 
  theme(axis.title = element_text(size = rel(1.1)), 
        axis.text.x = element_text(size =rel(1.3), color = "black"), 
        legend.position = "bottom") +
  labs(x = "Disease", y = "Eubacterium rectale", fill = "TRG",
       title = "6.3.4.21: Nicotinate phosphoribosyltransferase"); p_nicotinate_e.rectale

ggsave("Figure/7-02. Nicotinate phosphoribosyltransferase (species).svg", device = "svg",
       plot = p_nicotinate_e.rectale, height = 6, width = 3)

# ---------- 1-4. Phosphate acetyltransferase ----------
# Data frame
# Species: E. rectale
phosphate_transferase_long <- ec_valid %>% 
  filter(EC == "2.3.1.8: Phosphate acetyltransferase", !is.na(Species)) %>% 
  select(-Genus) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); phosphate_transferase_long

phosphate_transferase_stat <- phosphate_transferase_long %>% 
  group_by(Species) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(s_prev %>% select(Species, Prevalence), by = "Species") %>% 
  arrange(Wilcox); phosphate_transferase_stat

# Genus
phosphate_transferase_long2 <- ec_valid %>% 
  filter(EC == "2.3.1.8: Phosphate acetyltransferase", !is.na(Genus)) %>% 
  select(-Species) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); phosphate_transferase_long2

phosphate_transferase_stat2 <- phosphate_transferase_long2 %>% 
  group_by(Genus) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(g_prev %>% select(Genus, Prevalence), by = "Genus") %>% 
  arrange(Wilcox); phosphate_transferase_stat2

# write.csv(phosphate_transferase_stat, file = "Data/260522 Wilcoxon - phosphate acetyltransferase (species).csv", row.names = FALSE)

p_phosphate_e.rectale <- phosphate_transferase_long %>% 
  filter(Species == "Eubacterium_rectale") %>% 
  ggplot(aes(x = TRG_1, y = Abundance, fill = TRG_1)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CR", "nonCR"))) +
  theme_classic() + 
  theme(axis.title = element_text(size = rel(1.1)), 
        axis.text.x = element_text(size =rel(1.3), color = "black"), 
        legend.position = "bottom") +
  labs(x = "Disease", y = "Eubacterium rectale", fill = "TRG",
       title = "2.3.1.8: Phosphate acetyltransferase"); p_phosphate_e.rectale

p_phosphate_r.inulinivorans <- phosphate_transferase_long %>% 
  filter(Species == "Roseburia_inulinivorans") %>% 
  ggplot(aes(x = TRG_1, y = Abundance, fill = TRG_1)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CR", "nonCR"))) +
  theme_classic() + 
  theme(axis.title = element_text(size = rel(1.1)), 
        axis.text.x = element_text(size =rel(1.3), color = "black"), 
        legend.position = "bottom") +
  labs(x = "Disease", y = "Roseburia inulinivorans", fill = "TRG",
       title = "2.3.1.8: Phosphate acetyltransferase"); p_phosphate_r.inulinivorans

p_phosphate_b.plebeius <- phosphate_transferase_long %>% 
  filter(Species == "Bacteroides_plebeius") %>% 
  ggplot(aes(x = TRG_1, y = Abundance, fill = TRG_1)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CR", "nonCR"))) +
  theme_classic() + 
  theme(axis.title = element_text(size = rel(1.1)), 
        axis.text.x = element_text(size =rel(1.3), color = "black"), 
        legend.position = "bottom") +
  labs(x = "Disease", y = "Bacteroides plebeius", fill = "TRG",
       title = "2.3.1.8: Phosphate acetyltransferase"); p_phosphate_b.plebeius

p_phosphate_transferase <- ggarrange(p_phosphate_e.rectale, p_phosphate_r.inulinivorans, p_phosphate_b.plebeius,
                                     nrow = 1, ncol = 3); p_phosphate_transferase

ggsave("Figure/7-03. Phosphate acetyltransferase (species).svg", device = "svg",
       plot = p_phosphate_transferase, height = 6, width = 9)

# ---------- 1-5. NAD(+) synthase ----------
# Data frame
# Species: D. longicatena, R. bromii, E. rectale, F. saccharivorans
nad_synthase_long <- ec_valid %>% 
  filter(EC == "6.3.1.5: NAD(+) synthase", !is.na(Species)) %>% 
  select(-Genus) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); nad_synthase_long

nad_synthase_stat <- nad_synthase_long %>% 
  group_by(Species) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(s_prev %>% select(Species, Prevalence), by = "Species") %>% 
  arrange(Wilcox); nad_synthase_stat

# Genus
nad_synthase_long2 <- ec_valid %>% 
  filter(EC == "6.3.1.5: NAD(+) synthase", !is.na(Genus)) %>% 
  select(-Species) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); nad_synthase_long2

nad_synthase_stat2 <- nad_synthase_long2 %>% 
  group_by(Genus) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(g_prev %>% select(Genus, Prevalence), by = "Genus") %>% 
  arrange(Wilcox); nad_synthase_stat2

# write.csv(nad_synthase_stat, file = "Data/260522 Wilcoxon - nad synthase (species).csv", row.names = FALSE)

# ---------- 1-6. Nicotinate-nucleotide diphosphorylase ----------
# Data frame
# Species: D. longicatena, R. bromii, E. rectale, F. saccharivorans
nicotinate_nucleotide_long <- ec_valid %>% 
  filter(EC == "2.4.2.19: Nicotinate-nucleotide diphosphorylase (carboxylating)", !is.na(Species)) %>% 
  select(-Genus) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); nicotinate_nucleotide_long

nicotinate_nucleotide_stat <- nicotinate_nucleotide_long %>% 
  group_by(Species) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(s_prev %>% select(Species, Prevalence), by = "Species") %>% 
  arrange(Wilcox); nicotinate_nucleotide_stat

# Genus
nicotinate_nucleotide_long2 <- ec_valid %>% 
  filter(EC == "2.4.2.19: Nicotinate-nucleotide diphosphorylase (carboxylating)", !is.na(Genus)) %>% 
  select(-Species) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); nicotinate_nucleotide_long2

nicotinate_nucleotide_stat2 <- nicotinate_nucleotide_long2 %>% 
  group_by(Genus) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(g_prev %>% select(Genus, Prevalence), by = "Genus") %>% 
  arrange(Wilcox); nicotinate_nucleotide_stat2

# write.csv(nicotinate_nucleotide_stat, file = "Data/260522 Wilcoxon - nicotinate-nucleotide diphosphorylase (species).csv", row.names = FALSE)

p_nicotinate2_e.rectale <- nicotinate_nucleotide_long %>% 
  filter(Species == "Eubacterium_rectale") %>% 
  ggplot(aes(x = TRG_1, y = Abundance, fill = TRG_1)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("CR", "nonCR"))) +
  theme_classic() + 
  theme(axis.title = element_text(size = rel(1.1)), 
        axis.text.x = element_text(size =rel(1.3), color = "black"), 
        legend.position = "bottom") +
  labs(x = "Disease", y = "Eubacterium rectale", fill = "TRG",
       title = "2.4.2.19: Nicotinate-nucleotide diphosphorylase (carboxylating)"); p_nicotinate2_e.rectale

ggsave("Figure/7-04. Nicotinate-nucleotide diphosphorylase (species).svg", device = "svg",
       plot = p_nicotinate2_e.rectale, height = 6, width = 3)

# ---------- 1-7. Nicotinamidase ----------
# Data frame
# Species: D. longicatena, R. bromii, E. rectale, F. saccharivorans
nicotinamidase_long <- ec_valid %>% 
  filter(EC == "3.5.1.19: Nicotinamidase", !is.na(Species)) %>% 
  select(-Genus) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); nicotinamidase_long

nicotinamidase_stat <- nicotinamidase_long %>% 
  group_by(Species) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(s_prev %>% select(Species, Prevalence), by = "Species") %>% 
  arrange(Wilcox); nicotinamidase_stat

# Genus
nicotinamidase_long2 <- ec_valid %>% 
  filter(EC == "3.5.1.19: Nicotinamidase", !is.na(Genus)) %>% 
  select(-Species) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); nicotinamidase_long2

nicotinamidase_stat2 <- nicotinamidase_long2 %>% 
  group_by(Genus) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(g_prev %>% select(Genus, Prevalence), by = "Genus") %>% 
  arrange(Wilcox); nicotinamidase_stat2

# ---------- 1-8. Lactoyl-CoA dehydratase ----------
# Data frame
# Species: ?
lactoyl_dehydratase_long <- ec_valid %>% 
  filter(EC == "4.2.1.54: Lactoyl-CoA dehydratase", !is.na(Species)) %>% 
  select(-Genus) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); lactoyl_dehydratase_long

lactoyl_dehydratase_stat <- lactoyl_dehydratase_long %>% 
  group_by(Species) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(s_prev %>% select(Species, Prevalence), by = "Species") %>% 
  arrange(Wilcox); lactoyl_dehydratase_stat

# Genus
lactoyl_dehydratase_long2 <- ec_valid %>% 
  filter(EC == "4.2.1.54: Lactoyl-CoA dehydratase", !is.na(Genus)) %>% 
  select(-Species) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); lactoyl_dehydratase_long2

lactoyl_dehydratase_stat2 <- lactoyl_dehydratase_long2 %>% 
  group_by(Genus) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(g_prev %>% select(Genus, Prevalence), by = "Genus") %>% 
  arrange(Wilcox); lactoyl_dehydratase_stat2

# ---------- 1-9. Propionate CoA-transferase ----------
# Data frame
# Species: Dorea, E. rectale
propionate_transferase_long <- ec_valid %>% 
  filter(EC == "2.8.3.1: Propionate CoA-transferase", !is.na(Species)) %>% 
  select(-Genus) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); propionate_transferase_long

propionate_transferase_stat <- propionate_transferase_long %>% 
  group_by(Species) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(s_prev %>% select(Species, Prevalence), by = "Species") %>% 
  arrange(Wilcox); propionate_transferase_stat

# Genus
propionate_transferase_long2 <- ec_valid %>% 
  filter(EC == "2.8.3.1: Propionate CoA-transferase", !is.na(Genus)) %>% 
  select(-Species) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); propionate_transferase_long2

propionate_transferase_stat2 <- propionate_transferase_long2 %>% 
  group_by(Genus) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(g_prev %>% select(Genus, Prevalence), by = "Genus") %>% 
  arrange(Wilcox); propionate_transferase_stat2

# ---------- 1-10. Indolepyruvate decarboxylase ----------
# Data frame
# Species: Dorea, E. rectale
indolepyruvate_decarboxylase_long <- ec_valid %>% 
  filter(EC == "4.1.1.74: Indolepyruvate decarboxylase", !is.na(Species)) %>% 
  select(-Genus) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); indolepyruvate_decarboxylase_long

indolepyruvate_decarboxylase_stat <- indolepyruvate_decarboxylase_long %>% 
  group_by(Species) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(s_prev %>% select(Species, Prevalence), by = "Species") %>% 
  arrange(Wilcox); indolepyruvate_decarboxylase_stat

# Genus
indolepyruvate_decarboxylase_long2 <- ec_valid %>% 
  filter(EC == "4.1.1.74: Indolepyruvate decarboxylase", !is.na(Genus)) %>% 
  select(-Species) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); indolepyruvate_decarboxylase_long2

indolepyruvate_decarboxylase_stat2 <- indolepyruvate_decarboxylase_long2 %>% 
  group_by(Genus) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(g_prev %>% select(Genus, Prevalence), by = "Genus") %>% 
  arrange(Wilcox); indolepyruvate_decarboxylase_stat2

# ---------- 1-11. Cinnamoyl-CoA ----------
# Data frame
# Species: Faecalibacterium, E. rectale
cinnamoyl_coa_long <- ec_valid %>% 
  filter(EC == "2.8.3.17: Cinnamoyl-CoA:phenyllactate CoA-transferase", !is.na(Species)) %>% 
  select(-Genus) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); cinnamoyl_coa_long

cinnamoyl_coa_stat <- cinnamoyl_coa_long %>% 
  group_by(Species) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(s_prev %>% select(Species, Prevalence), by = "Species") %>% 
  arrange(Wilcox); cinnamoyl_coa_stat

# Genus
cinnamoyl_coa_long2 <- ec_valid %>% 
  filter(EC == "2.8.3.17: Cinnamoyl-CoA:phenyllactate CoA-transferase", !is.na(Genus)) %>% 
  select(-Species) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); cinnamoyl_coa_long2

cinnamoyl_coa_stat2 <- cinnamoyl_coa_long2 %>% 
  group_by(Genus) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(g_prev %>% select(Genus, Prevalence), by = "Genus") %>% 
  arrange(Wilcox); cinnamoyl_coa_stat2

# ---------- 1-12. Butyrate kinase ----------
# Data frame
# Species: C. comes, E. rectale
butyrate_kinase_long <- ec_valid %>% 
  filter(EC == "2.7.2.7: Butyrate kinase", !is.na(Species)) %>% 
  select(-Genus) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); butyrate_kinase_long

butyrate_kinase_stat <- butyrate_kinase_long %>% 
  group_by(Species) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(s_prev %>% select(Species, Prevalence), by = "Species") %>% 
  arrange(Wilcox); butyrate_kinase_stat

# Genus
butyrate_kinase_long2 <- ec_valid %>% 
  filter(EC == "2.7.2.7: Butyrate kinase", !is.na(Genus)) %>% 
  select(-Species) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); butyrate_kinase_long2

butyrate_kinase_stat2 <- butyrate_kinase_long2 %>% 
  group_by(Genus) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(g_prev %>% select(Genus, Prevalence), by = "Genus") %>% 
  arrange(Wilcox); butyrate_kinase_stat2

# write.csv(butyrate_kinase_stat, file = "Data/260522 Wilcoxon - butyrate kinase (species).csv", row.names = FALSE)

# ---------- 1-13. Phosphate butyryltransferase ----------
# Data frame
# Species: C. comes, E. rectale
phosphate_butyryltransferase_long <- ec_valid %>% 
  filter(EC == "2.3.1.19: Phosphate butyryltransferase", !is.na(Species)) %>% 
  select(-Genus) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); phosphate_butyryltransferase_long

phosphate_butyryltransferase_stat <- phosphate_butyryltransferase_long %>% 
  group_by(Species) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(s_prev %>% select(Species, Prevalence), by = "Species") %>% 
  arrange(Wilcox); phosphate_butyryltransferase_stat

# Genus
phosphate_butyryltransferase_long2 <- ec_valid %>% 
  filter(EC == "2.3.1.19: Phosphate butyryltransferase", !is.na(Genus)) %>% 
  select(-Species) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); phosphate_butyryltransferase_long2

phosphate_butyryltransferase_stat2 <- phosphate_butyryltransferase_long2 %>% 
  group_by(Genus) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(g_prev %>% select(Genus, Prevalence), by = "Genus") %>% 
  arrange(Wilcox); phosphate_butyryltransferase_stat2

# ---------- 1-14. Acetate CoA-transferase ----------
# Data frame
# Species: E. rectale, C. comes
acetate_transferase_long <- ec_valid %>% 
  filter(EC == "2.8.3.8: Acetate CoA-transferase", !is.na(Species)) %>% 
  select(-Genus) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); acetate_transferase_long

acetate_transferase_stat <- acetate_transferase_long %>% 
  group_by(Species) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(s_prev %>% select(Species, Prevalence), by = "Species") %>% 
  arrange(Wilcox); acetate_transferase_stat

# write.csv(acetate_transferase_stat, file = "Data/260522 Wilcoxon - acetate coa-transferase (species).csv", row.names = FALSE)

# Genus
acetate_transferase_long2 <- ec_valid %>% 
  filter(EC == "2.8.3.8: Acetate CoA-transferase", !is.na(Genus)) %>% 
  select(-Species) %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); acetate_transferase_long2

acetate_transferase_stat2 <- acetate_transferase_long2 %>% 
  group_by(Genus) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  left_join(g_prev %>% select(Genus, Prevalence), by = "Genus") %>% 
  arrange(Wilcox); acetate_transferase_stat2

# ---------- 1-15. Envfit PCoA ----------
# EC significant
ec_list <- c("2.7.2.1: Acetate kinase", "6.3.4.21: Nicotinate phosphoribosyltransferase", "2.3.1.8: Phosphate acetyltransferase", "2.4.2.19: Nicotinate-nucleotide diphosphorylase (carboxylating)")

ec_sig <- ec_valid %>% 
  filter(is.na(Genus) & is.na(Species) & EC %in% ec_list) %>% 
  select(-Genus, -Species) %>% 
  group_by(EC) %>% 
  slice(1) %>% 
  column_to_rownames("EC") %>% 
  t() %>% 
  as.data.frame(); ec_sig
  

# Envfit calculation
ec_meta_bac_df <- cbind(meta_sig, bac_sig2, ec_sig); ec_meta_bac_df

set.seed(42)
ef_ec_meta_bac <- envfit(pcoa_bc$points, 
                      ec_meta_bac_df,
                      permutations = 999); ef_ec_meta_bac
#                                                                     Dim1     Dim2     r2 Pr(>r)    
# Butyrate                                                        -0.86913 -0.49458 0.7786  0.001 ***
# Acetate                                                         -0.85939 -0.51131 0.9566  0.001 ***
# Propionate                                                      -0.80711 -0.59040 0.9120  0.001 ***
# Indole_acetic_acid                                               0.36299  0.93179 0.0848  0.423    
# Indole_lactic_acid                                              -0.85527  0.51818 0.1375  0.248    
# Indolepropionic_acid                                            -0.78075 -0.62484 0.2193  0.130    
# Nicotinic_acid                                                  -0.94349  0.33139 0.4610  0.006 ** 
# Coprococcus_comes                                               -0.98825 -0.15285 0.1593  0.174    
# Eubacterium_rectale                                             -0.98297 -0.18374 0.1502  0.256    
# Dorea_longicatena                                               -0.94143 -0.33721 0.1805  0.160    
# Phocaeicola_vulgatus                                             0.66164  0.74982 0.0297  0.683    
# Roseburia_inulinivorans                                         -0.39865  0.91710 0.0713  0.509    
# Eggerthella_lenta                                                0.67802 -0.73504 0.2912  0.085 .  
# 2.3.1.8: Phosphate acetyltransferase                            -0.88243 -0.47044 0.0459  0.684    
# 2.4.2.19: Nicotinate-nucleotide diphosphorylase (carboxylating) -0.21003  0.97769 0.0067  0.949    
# 2.7.2.1: Acetate kinase                                         -0.63966  0.76866 0.1406  0.275    
# 6.3.4.21: Nicotinate phosphoribosyltransferase                  -0.99794 -0.06413 0.1810  0.189

# Vector coordinates 추출
vec_ec_meta_bac <- scores(ef_ec_meta_bac, display = "vectors")

vec_ec_meta_bac <- as.data.frame(vec_ec_meta_bac) %>% 
  rename(PCoA1 = Dim1,
         PCoA2 = Dim2) %>% 
  rownames_to_column("Feature") %>% 
  mutate(pval = ef_ec_meta_bac$vectors$pvals,
         r2 = ef_ec_meta_bac$vectors$r,
         group = case_when(str_detect(Feature, "Indole") ~ "Tryptophan metabolites",
                           Feature == "Nicotinic_acid" ~ "Tryptophan metabolites",
                           Feature %in% c("Butyrate", "Acetate", "Propionate") ~ "SCFA",
                           Feature %in% ec_list ~ "EC",
                           TRUE ~ "Bacteria"),
         label = str_replace(Feature, "_", " ")) %>% 
  mutate(label = case_when(label == "Indolepropionic acid" ~ "Indole propionic acid",
                           label == "Indole lactic_acid" ~ "Indole lactic acid",
                           label == "Indole acetic_acid" ~ "Indole acetic acid",
                           label == "2.3.1.8: Phosphate acetyltransferase" ~ "Phosphate acetyltransferase",
                           label == "2.4.2.19: Nicotinate-nucleotide diphosphorylase (carboxylating)" ~ "Nicotinate-nucleotide diphosphorylase",
                           label == "2.7.2.1: Acetate kinase" ~ "Acetate kinase",
                           label == "6.3.4.21: Nicotinate phosphoribosyltransferase" ~ "Nicotinate phosphoribosyltransferase",
                           TRUE ~ label))

# Arrow scaling (optional)
arrow_mult <- 0.6

vec_ec_meta_bac <- vec_ec_meta_bac %>%
  mutate(PCoA1 = PCoA1 * arrow_mult,
         PCoA2 = PCoA2 * arrow_mult)

### PCoA plot
p_pcoa_ec_meta_bac1 <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, fill = Group)) +
  
  # envfit vectors
  geom_segment(data = vec_ec_meta_bac,
               aes(x = 0, y = 0, xend = PCoA1, yend = PCoA2),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 0.8,
               alpha = 0.4,
               color = "black") +
  
  geom_text_repel(data = vec_ec_meta_bac,
                  aes(x = PCoA1, y = PCoA2, label = label, color = group),
                  inherit.aes = FALSE,
                  size = 3.5,
                  fontface = "bold",
                  segment.color = NA,
                  show.legend = TRUE,
                  
                  box.padding = 0.5,      # label끼리 거리
                  point.padding = 0.5,    # 점/벡터 끝과 거리
                  force = 3,              # 밀어내는 힘
                  max.overlaps = Inf,     # 겹쳐도 삭제 안함
                  min.segment.length = 0, # 짧은 연결선도 허용
                  seed = 123) +           # 위치 고정
  
  # point coordinating
  geom_point(shape  = 21,
             size   = 3.5,
             color  = "white",
             stroke = 0.9) +
  
  # PERMANOVA 결과 annotation
  annotate("text", x = Inf, y = -0.32,
           hjust = 1.05, vjust = 1.3,
           label = sprintf("PERMANOVA: R² = %.4f, %s", r2_bc, fmt_pval(pv_bc)),
           size = 3.5) +
  
  scale_fill_manual(name = "TRG", values  = group_colors) +
  scale_color_manual(name = "",
                     values = c("SCFA" = "#6FA67C",
                                "Tryptophan metabolites" = "#6D7FA1",
                                "EC" = "#9A84B8",
                                "Bacteria" = "#C98C5A")) +
  
  # geom_hline(yintercept = 0, linetype = "dashed", color = "gray85", linewidth = 0.3) +
  # geom_vline(xintercept = 0, linetype = "dashed", color = "gray85", linewidth = 0.3) +
  
  labs(x = sprintf("PCoA1 (%.1f%%)", pcoa1_var),
       y = sprintf("PCoA2 (%.1f%%)", pcoa2_var),
       title = "Metabolite & Bacteria & Enzyme (Wilcoxon)") +
  
  theme_classic2() +
  theme(aspect.ratio = 1,
        legend.position  = "bottom",
        legend.title     = element_text(face = "bold"),
        axis.title       = element_text(face = "bold")) +
  
  guides(fill = guide_legend(override.aes = list(color = NA)),
         color = guide_legend(override.aes = list(linewidth = 1.2, alpha = 1))); p_pcoa_ec_meta_bac1

# ==================== 2. Correlation ====================

# ---------- 2-1. Full screening ----------
# Data frame
ec_wide <- ec_valid %>% 
  select(-Genus, -Species) %>% 
  group_by(EC) %>% 
  slice(1) %>% 
  column_to_rownames("EC") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleID"); ec_wide

meta_wide <- metabolite_before_valid %>% 
  select(-TRG_1, -Tstage); meta_wide

# Correlation data frame
ec_meta_df <- ec_wide %>% 
  left_join(meta_wide, by = "SampleID") %>%
  column_to_rownames("SampleID") %>% 
  as.matrix(); ec_meta_df

# Correlation
corr_ec_meta <- Hmisc::rcorr(ec_meta_df, type = "spearman")
corr_ec_meta2 <- Hmisc::rcorr(ec_meta_df, type = "pearson")

# Result table
ec_var <- colnames(ec_wide)[-1]; ec_var
meta_var <- colnames(meta_wide)[-1]; meta_var

corr_mat <- corr_ec_meta$r
p_mat <- corr_ec_meta$P

corr_mat2 <- corr_ec_meta2$r
p_mat2 <- corr_ec_meta2$P

corr_ec_meta_res <- expand.grid(EC = ec_var,
                                Metabolite = meta_var,
                                stringsAsFactors = FALSE) %>% 
  mutate(Corr_s = map2_dbl(EC, Metabolite, \(p, f) corr_mat[p, f]),
         Pval_s = map2_dbl(EC, Metabolite, \(p, f) p_mat[p, f]),
         Corr_p = map2_dbl(EC, Metabolite, \(p, f) corr_mat2[p, f]),
         Pval_p = map2_dbl(EC, Metabolite, \(p, f) p_mat2[p, f])) %>% 
  mutate(Pval_merged = (Pval_s + Pval_p) / 2) %>% 
  arrange(Pval_merged) %>% 
  drop_na(Pval_merged)
  
head(corr_ec_meta_res, 10)
#                                           EC                Metabolite       Corr         Pval       FDR
# 1  2.3.1.222: Phosphate propanoyltransferase     Glycodeoxycholic_acid  0.8476872 2.377738e-06 0.1195803
# 2             2.4.2.3: Uridine phosphorylase          Deoxycholic_acid  0.8451128 2.737644e-06 0.1195803
# 3                          3.5.4.43: NO_NAME Glycoursodeoxycholic_acid -0.8069498 1.713844e-05 0.3722181
# 4                      2.4.1.4: Amylosucrase     Glycodeoxycholic_acid  0.7965401 2.639494e-05 0.3722181
# 5    1.21.4.1: D-proline reductase (dithiol)     Glycodeoxycholic_acid  0.7957879 2.720569e-05 0.3722181
# 6             2.4.2.3: Uridine phosphorylase     Glycodeoxycholic_acid  0.7927793 3.066797e-05 0.3722181
# 7                          2.1.3.15: NO_NAME             Glutamic_acid  0.7909774 3.291873e-05 0.3722181
# 8                5.4.2.7: Phosphopentomutase     Taurodeoxycholic_acid  0.7900846 3.408590e-05 0.3722181
# 9  1.14.13.39: Nitric-oxide synthase (NADPH)               Aminophenol  0.7844584 4.229743e-05 0.4105671
# 10           3.4.13.19: Membrane dipeptidase      Indolepropionic_acid  0.7759398 5.796529e-05 0.5063848

#                                                            EC            Metabolite     Corr_s       Pval_s     Corr_p       Pval_p  Pval_merged
# 1                   2.3.1.222: Phosphate propanoyltransferase Glycodeoxycholic_acid  0.8476872 2.377738e-06  0.8355713 4.519200e-06 3.448469e-06
# 2                              2.4.2.3: Uridine phosphorylase      Deoxycholic_acid  0.8451128 2.737644e-06  0.8030180 2.023444e-05 1.148604e-05
# 3                                       2.4.1.4: Amylosucrase Glycodeoxycholic_acid  0.7965401 2.639494e-05  0.7673496 7.859470e-05 5.249482e-05
# 4                                           2.1.3.15: NO_NAME         Glutamic_acid  0.7909774 3.291873e-05  0.7598479 1.014945e-04 6.720662e-05
# 5                                       2.4.1.4: Amylosucrase      Deoxycholic_acid  0.7368421 2.108078e-04  0.8175669 1.073664e-05 1.107722e-04
# 6                     1.21.4.1: D-proline reductase (dithiol)      Deoxycholic_acid  0.7323308 2.412164e-04  0.7685421 7.539979e-05 1.583081e-04
# 7  2.1.2.11: 3-methyl-2-oxobutanoate hydroxymethyltransferase Glycodeoxycholic_acid -0.7544190 1.214468e-04 -0.7228952 3.171068e-04 2.192768e-04
# 8         2.1.1.107: Uroporphyrinogen-III C-methyltransferase            Tryptophan  0.7548872 1.196028e-04  0.7180910 3.629983e-04 2.413006e-04
# 9        2.1.1.217: tRNA (adenine(22)-N(1))-methyltransferase Glycodeoxycholic_acid  0.7115457 4.345361e-04  0.7453353 1.623752e-04 2.984557e-04
# 10                               3.4.24.78: GPR endopeptidase         Glutamic_acid  0.7218045 3.270659e-04  0.7248521 2.998860e-04 3.134760e-04

corr_ec_meta_file <- corr_ec_meta_res %>% 
  left_join(ec_long %>% 
              group_by(EC) %>% 
              summarise(Abundance = sum(Abundance), .groups = "drop") %>% 
              rename(ECAbund = Abundance), by = "EC") %>% 
  left_join(meta_long %>% 
              rename(MetAbund = Abundance), by = "Metabolite") %>% 
  relocate(ECAbund, .after = EC) %>% 
  relocate(MetAbund, .after = Metabolite); corr_ec_meta_file

# write.csv(corr_ec_meta_file, file = "Data/260526 Correlation - enzyme n metabolite.csv", row.names = FALSE)

ec_meta_plot <- as.data.frame(ec_meta_df) %>% 
  rownames_to_column("SampleID") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1), by = "SampleID"); ec_meta_plot

# Significant no.1
p_nad_glutamine_niacin <- ec_meta_plot %>% 
  ggplot(aes(x = log10(Nicotinic_acid + 1), y = log10(`6.3.5.1: NAD(+) synthase (glutamine-hydrolyzing)` + 1))) +
  geom_point(size = rel(3), aes(color = TRG_1)) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  stat_cor(aes(label = paste("Spearman:", after_stat(r.label), after_stat(p.label))),
           method = "spearman", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.95) +
  stat_cor(aes(label = paste("Pearson:", after_stat(r.label), after_stat(p.label))),
           method = "pearson", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.90) +
  geom_smooth(method = "lm", color = "gray70", se = FALSE, alpha = 0.2) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(x = "Log10(Nicotinic acid + 1)",
       y = "Log10(EC:6.3.5.1 + 1)",
       title = "NAD(+) synthase (glutamine-hydrolyzing)"); p_nad_glutamine_niacin

### Spearman correlation p-value
# Significant no.2
p_nad_synthase_niacin <- ec_meta_plot %>% 
  ggplot(aes(x = log10(Nicotinic_acid + 1), y = log10(`6.3.1.5: NAD(+) synthase` + 1))) +
  geom_point(size = rel(3), aes(color = TRG_1)) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  stat_cor(aes(label = paste("Spearman:", after_stat(r.label), after_stat(p.label))),
           method = "spearman", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.95) +
  stat_cor(aes(label = paste("Pearson:", after_stat(r.label), after_stat(p.label))),
           method = "pearson", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.90) +
  geom_smooth(method = "lm", color = "gray70", se = FALSE, alpha = 0.2) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(x = "Log10(Nicotinic acid + 1)",
       y = "Log10(EC:6.3.1.5 + 1)",
       title = "NAD(+) synthase"); p_nad_synthase_niacin

# Significant no.3
p_nicotinamidase_ipa <- ec_meta_plot %>% 
  ggplot(aes(x = log10(Indolepropionic_acid + 1), y = log10(`3.5.1.19: Nicotinamidase` + 1))) +
  geom_point(size = rel(3), aes(color = TRG_1)) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  stat_cor(aes(label = paste("Spearman:", after_stat(r.label), after_stat(p.label))),
           method = "spearman", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.95) +
  stat_cor(aes(label = paste("Pearson:", after_stat(r.label), after_stat(p.label))),
           method = "pearson", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.90) +
  geom_smooth(method = "lm", color = "gray70", se = FALSE, alpha = 0.2) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(x = "Log10(Indole propionic acid + 1)",
       y = "Log10(EC:3.5.1.19 + 1)",
       title = "Nicotinamidase"); p_nicotinamidase_ipa

# Significant no.4
p_acetate_coa_propionate <- ec_meta_plot %>% 
  ggplot(aes(x = log10(Propionate + 1), y = log10(`2.8.3.8: Acetate CoA-transferase` + 1))) +
  geom_point(size = rel(3), aes(color = TRG_1)) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  stat_cor(aes(label = paste("Spearman:", after_stat(r.label), after_stat(p.label))),
           method = "spearman", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.95) +
  stat_cor(aes(label = paste("Pearson:", after_stat(r.label), after_stat(p.label))),
           method = "pearson", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.90) +
  geom_smooth(method = "lm", color = "gray70", se = FALSE, alpha = 0.2) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(x = "Log10(Propionate + 1)",
       y = "Log10(EC:2.8.3.8 + 1)",
       title = "Acetate CoA-transferase"); p_acetate_coa_propionate

# Significant no.5
p_cinnamoyl_coa_niacin <- ec_meta_plot %>% 
  ggplot(aes(x = log10(Nicotinic_acid + 1), y = log10(`2.8.3.17: Cinnamoyl-CoA:phenyllactate CoA-transferase` + 1))) +
  geom_point(size = rel(3), aes(color = TRG_1)) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  stat_cor(aes(label = paste("Spearman:", after_stat(r.label), after_stat(p.label))),
           method = "spearman", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.95) +
  stat_cor(aes(label = paste("Pearson:", after_stat(r.label), after_stat(p.label))),
           method = "pearson", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.90) +
  geom_smooth(method = "lm", color = "gray70", se = FALSE, alpha = 0.2) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(x = "Log10(Nicotinic acid + 1)",
       y = "Log10(EC:2.8.3.17 + 1)",
       title = "Cinnamoyl-CoA: phenyllactate CoA-transferase"); p_cinnamoyl_coa_niacin

p_ec_metabolite <- ggarrange(p_nad_glutamine_niacin, p_nad_synthase_niacin, p_nicotinamidase_ipa,
                             p_acetate_coa_propionate, p_cinnamoyl_coa_niacin); p_ec_metabolite

ggsave("Figure/7-05. EC-metabolite scatter - screening.svg", device = "svg", 
       plot = p_ec_metabolite, width = 12, height = 8)

### Butyrate-related
# 1) EC 2.8.3.8: Acetate CoA-transferase
p_acetate_coa_butyrate <- ec_meta_plot %>% 
  ggplot(aes(x = log10(Butyrate + 1), y = log10(`2.8.3.8: Acetate CoA-transferase` + 1))) +
  geom_point(size = rel(3), aes(color = TRG_1)) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  stat_cor(aes(label = paste("Spearman:", after_stat(r.label), after_stat(p.label))),
           method = "spearman", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.95) +
  stat_cor(aes(label = paste("Pearson:", after_stat(r.label), after_stat(p.label))),
           method = "pearson", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.90) +
  geom_smooth(method = "lm", color = "gray70", se = FALSE, alpha = 0.2) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(x = "Log10(Butyrate + 1)",
       y = "Log10(EC:2.8.3.8 + 1)",
       title = "Acetate CoA-transferase"); p_acetate_coa_butyrate

# 2) EC 2.3.1.19: Phosphate butyryltransferase
p_phosphate_butyryl_butyrate <- ec_meta_plot %>% 
  ggplot(aes(x = log10(Butyrate + 1), y = log10(`2.3.1.19: Phosphate butyryltransferase` + 1))) +
  geom_point(size = rel(3), aes(color = TRG_1)) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  stat_cor(aes(label = paste("Spearman:", after_stat(r.label), after_stat(p.label))),
           method = "spearman", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.95) +
  stat_cor(aes(label = paste("Pearson:", after_stat(r.label), after_stat(p.label))),
           method = "pearson", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.90) +
  geom_smooth(method = "lm", color = "gray70", se = FALSE, alpha = 0.2) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(x = "Log10(Butyrate + 1)",
       y = "Log10(EC:2.3.1.19 + 1)",
       title = "Phosphate butyryltransferase"); p_phosphate_butyryl_butyrate

# 3) EC 2.7.2.7: Butyrate kinase
p_butyrate_kinase_butyrate <- ec_meta_plot %>% 
  ggplot(aes(x = log10(Butyrate + 1), y = log10(`2.7.2.7: Butyrate kinase` + 1))) +
  geom_point(size = rel(3), aes(color = TRG_1)) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  stat_cor(aes(label = paste("Spearman:", after_stat(r.label), after_stat(p.label))),
           method = "spearman", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.95) +
  stat_cor(aes(label = paste("Pearson:", after_stat(r.label), after_stat(p.label))),
           method = "pearson", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.90) +
  geom_smooth(method = "lm", color = "gray70", se = FALSE, alpha = 0.2) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(x = "Log10(Butyrate + 1)",
       y = "Log10(EC:2.7.2.7 + 1)",
       title = "Butyrate kinase"); p_butyrate_kinase_butyrate

p_ec_butyrate <- ggarrange(p_acetate_coa_butyrate, p_phosphate_butyryl_butyrate, p_butyrate_kinase_butyrate,
                           nrow = 1, ncol = 3); p_ec_butyrate

ggsave("Figure/7-06. EC-butyrate scatter - screening.svg", device = "svg", 
       plot = p_ec_butyrate, width = 12, height = 4)

### Acetate-related
# 1) EC 2.3.1.8: Phosphate acetyltransferase
p_phosphate_acetyl_acetate <- ec_meta_plot %>% 
  ggplot(aes(x = log10(Acetate + 1), y = log10(`2.3.1.8: Phosphate acetyltransferase` + 1))) +
  geom_point(size = rel(3), aes(color = TRG_1)) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  stat_cor(aes(label = paste("Spearman:", after_stat(r.label), after_stat(p.label))),
           method = "spearman", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.95) +
  stat_cor(aes(label = paste("Pearson:", after_stat(r.label), after_stat(p.label))),
           method = "pearson", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.90) +
  geom_smooth(method = "lm", color = "gray70", se = FALSE, alpha = 0.2) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(x = "Log10(Acetate + 1)",
       y = "Log10(EC:2.3.1.8 + 1)",
       title = "Phosphate acetyltransferase"); p_phosphate_acetyl_acetate

# 2) EC 2.7.2.1: Acetate kinase
p_acetate_kinase_acetate <- ec_meta_plot %>% 
  ggplot(aes(x = log10(Acetate + 1), y = log10(`2.7.2.1: Acetate kinase` + 1))) +
  geom_point(size = rel(3), aes(color = TRG_1)) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  stat_cor(aes(label = paste("Spearman:", after_stat(r.label), after_stat(p.label))),
           method = "spearman", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.95) +
  stat_cor(aes(label = paste("Pearson:", after_stat(r.label), after_stat(p.label))),
           method = "pearson", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.90) +
  geom_smooth(method = "lm", color = "gray70", se = FALSE, alpha = 0.2) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(x = "Log10(Acetate + 1)",
       y = "Log10(EC:2.7.2.1 + 1)",
       title = "Acetate kinase"); p_acetate_kinase_acetate

# 3) EC 2.8.3.8: Acetate CoA-transferase
p_acetate_coa_acetate <- ec_meta_plot %>% 
  ggplot(aes(x = log10(Acetate + 1), y = log10(`2.8.3.8: Acetate CoA-transferase` + 1))) +
  geom_point(size = rel(3), aes(color = TRG_1)) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  stat_cor(aes(label = paste("Spearman:", after_stat(r.label), after_stat(p.label))),
           method = "spearman", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.95) +
  stat_cor(aes(label = paste("Pearson:", after_stat(r.label), after_stat(p.label))),
           method = "pearson", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.90) +
  geom_smooth(method = "lm", color = "gray70", se = FALSE, alpha = 0.2) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(x = "Log10(Acetate + 1)",
       y = "Log10(EC:2.8.3.8 + 1)",
       title = "Acetate CoA-transferase"); p_acetate_coa_acetate

p_ec_acetate <- ggarrange(p_phosphate_acetyl_acetate, p_acetate_kinase_acetate, p_acetate_coa_acetate,
                          nrow = 1, ncol = 3); p_ec_acetate

ggsave("Figure/7-07. EC-acetate scatter - screening.svg", device = "svg", 
       plot = p_ec_acetate, width = 12, height = 4)

### Propionate-related
# 1) EC 2.8.3.1: Propionate CoA-transferase
p_propionate_coa_propionate <- ec_meta_plot %>% 
  ggplot(aes(x = log10(Propionate + 1), y = log10(`2.8.3.1: Propionate CoA-transferase` + 1))) +
  geom_point(size = rel(3), aes(color = TRG_1)) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  stat_cor(aes(label = paste("Spearman:", after_stat(r.label), after_stat(p.label))),
           method = "spearman", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.95) +
  stat_cor(aes(label = paste("Pearson:", after_stat(r.label), after_stat(p.label))),
           method = "pearson", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.90) +
  geom_smooth(method = "lm", color = "gray70", se = FALSE, alpha = 0.2) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(x = "Log10(Propionate + 1)",
       y = "Log10(EC:2.8.3.1 + 1)",
       title = "Propionate CoA-transferase"); p_propionate_coa_propionate

# 2) EC 4.2.1.54: Lactoyl-CoA dehydratase
p_lactoyl_coa_propionate <- ec_meta_plot %>% 
  ggplot(aes(x = log10(Propionate + 1), y = log10(`4.2.1.54: Lactoyl-CoA dehydratase` + 1))) +
  geom_point(size = rel(3), aes(color = TRG_1)) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  stat_cor(aes(label = paste("Spearman:", after_stat(r.label), after_stat(p.label))),
           method = "spearman", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.95) +
  stat_cor(aes(label = paste("Pearson:", after_stat(r.label), after_stat(p.label))),
           method = "pearson", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.90) +
  geom_smooth(method = "lm", color = "gray70", se = FALSE, alpha = 0.2) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(x = "Log10(Propionate + 1)",
       y = "Log10(EC:4.2.1.54 + 1)",
       title = "Lactoyl-CoA dehydratase"); p_lactoyl_coa_propionate

p_ec_propionate <- ggarrange(p_propionate_coa_propionate, p_lactoyl_coa_propionate,
                             nrow = 1, ncol = 2); p_ec_propionate

ggsave("Figure/7-08. EC-propionate scatter - screening.svg", device = "svg", 
       plot = p_ec_propionate, width = 8, height = 4)

### Indole-related
# 1) EC 2.8.3.17: Cinnamoyl-CoA:phenyllactate CoA-transferase (IPA)
p_cinnamoyl_coa_ipa <- ec_meta_plot %>% 
  ggplot(aes(x = log10(Indolepropionic_acid + 1), y = log10(`2.8.3.17: Cinnamoyl-CoA:phenyllactate CoA-transferase` + 1))) +
  geom_point(size = rel(3), aes(color = TRG_1)) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  stat_cor(aes(label = paste("Spearman:", after_stat(r.label), after_stat(p.label))),
           method = "spearman", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.95) +
  stat_cor(aes(label = paste("Pearson:", after_stat(r.label), after_stat(p.label))),
           method = "pearson", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.90) +
  geom_smooth(method = "lm", color = "gray70", se = FALSE, alpha = 0.2) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(x = "Log10(Indole propionic acid + 1)",
       y = "Log10(EC:2.8.3.17 + 1)",
       title = "Cinnamoyl-CoA:phenyllactate CoA-transferase"); p_cinnamoyl_coa_ipa

# 2) EC 2.8.3.17: Cinnamoyl-CoA:phenyllactate CoA-transferase (ILA)
p_cinnamoyl_coa_ila <- ec_meta_plot %>% 
  ggplot(aes(x = log10(Indole_lactic_acid + 1), y = log10(`2.8.3.17: Cinnamoyl-CoA:phenyllactate CoA-transferase` + 1))) +
  geom_point(size = rel(3), aes(color = TRG_1)) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  stat_cor(aes(label = paste("Spearman:", after_stat(r.label), after_stat(p.label))),
           method = "spearman", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.95) +
  stat_cor(aes(label = paste("Pearson:", after_stat(r.label), after_stat(p.label))),
           method = "pearson", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.90) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(x = "Log10(Indole lactic acid + 1)",
       y = "Log10(EC:2.8.3.17 + 1)",
       title = "Cinnamoyl-CoA:phenyllactate CoA-transferase"); p_cinnamoyl_coa_ila

# 3) EC 4.1.1.74: Indolepyruvate decarboxylase (IAA)
p_indolepyruvate_decarboxylase_iaa <- ec_meta_plot %>% 
  ggplot(aes(x = log10(Indole_acetic_acid + 1), y = log10(`4.1.1.74: Indolepyruvate decarboxylase` + 1))) +
  geom_point(size = rel(3), aes(color = TRG_1)) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  stat_cor(aes(label = paste("Spearman:", after_stat(r.label), after_stat(p.label))),
           method = "spearman", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.95) +
  stat_cor(aes(label = paste("Pearson:", after_stat(r.label), after_stat(p.label))),
           method = "pearson", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.90) +
  geom_smooth(method = "lm", color = "gray70", se = FALSE, alpha = 0.2) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(x = "Log10(Indole acetic acid + 1)",
       y = "Log10(EC:4.1.1.74 + 1)",
       title = "Indolepyruvate decarboxylase"); p_indolepyruvate_decarboxylase_iaa

p_ec_indole <- ggarrange(p_cinnamoyl_coa_ipa, p_cinnamoyl_coa_ila, p_indolepyruvate_decarboxylase_iaa,
                         nrow = 1, ncol = 3); p_ec_indole

ggsave("Figure/7-09. EC-indole scatter - screening.svg", device = "svg", 
       plot = p_ec_indole, width = 12, height = 4)

### Niacin-related
# 1) EC 3.5.1.19: Nicotinamidase
p_nicotinamidase_niacin <- ec_meta_plot %>% 
  ggplot(aes(x = log10(Nicotinic_acid + 1), y = log10(`3.5.1.19: Nicotinamidase` + 1))) +
  geom_point(size = rel(3), aes(color = TRG_1)) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  stat_cor(aes(label = paste("Spearman:", after_stat(r.label), after_stat(p.label))),
           method = "spearman", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.95) +
  stat_cor(aes(label = paste("Pearson:", after_stat(r.label), after_stat(p.label))),
           method = "pearson", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.90) +
  geom_smooth(method = "lm", color = "gray70", se = FALSE, alpha = 0.2) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(x = "Log10(Nicotinic acid + 1)",
       y = "Log10(EC:3.5.1.19 + 1)",
       title = "Nicotinamidase"); p_nicotinamidase_niacin

# 2) EC 6.3.4.21: Nicotinate phosphoribosyltransferase
p_nicotinate_phophoribosyl_niacin <- ec_meta_plot %>% 
  ggplot(aes(x = log10(Nicotinic_acid + 1), y = log10(`6.3.4.21: Nicotinate phosphoribosyltransferase` + 1))) +
  geom_point(size = rel(3), aes(color = TRG_1)) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  stat_cor(aes(label = paste("Spearman:", after_stat(r.label), after_stat(p.label))),
           method = "spearman", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.95) +
  stat_cor(aes(label = paste("Pearson:", after_stat(r.label), after_stat(p.label))),
           method = "pearson", output.type = "text", label.x.npc = 0.02, label.y.npc = 0.90) +
  geom_smooth(method = "lm", color = "gray70", se = FALSE, alpha = 0.2) +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "none") +
  labs(x = "Log10(Nicotinic acid + 1)",
       y = "Log10(EC:6.3.4.21 + 1)",
       title = "Nicotinate phosphoribosyltransferase"); p_nicotinate_phophoribosyl_niacin

p_ec_niacin <- ggarrange(p_nicotinamidase_niacin, p_nicotinate_phophoribosyl_niacin,
                         nrow = 1, ncol = 2); p_ec_niacin

ggsave("Figure/7-10. EC-niacin scatter - screening.svg", device = "svg", 
       plot = p_ec_niacin, width = 8, height = 4)

# ---------- 2-2. NAD(+) synthase (glutamine-hydrolyzing) ----------
nad_synthase_wide <- ec_valid %>% 
  filter(EC == "6.3.5.1: NAD(+) synthase (glutamine-hydrolyzing)" & !is.na(Species)) %>% 
  select(-EC, -Genus) %>%
  column_to_rownames("Species") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleID"); nad_synthase_wide

# Correlation data frame
nad_synthase_meta_df <- nad_synthase_wide %>% 
  left_join(meta_wide, by = "SampleID") %>%
  column_to_rownames("SampleID") %>% 
  as.matrix(); nad_synthase_meta_df

# Correlation
corr_nad_synthase_meta <- Hmisc::rcorr(nad_synthase_meta_df, type = "spearman")

# Result table
nad_synthase_var <- colnames(nad_synthase_wide)[-1]; nad_synthase_var
meta_var <- colnames(meta_wide)[-1]; meta_var

corr_mat <- corr_nad_synthase_meta$r
p_mat <- corr_nad_synthase_meta$P

corr_nad_synthase_meta_res <- expand.grid(Species = nad_synthase_var,
                                            Metabolite = meta_var,
                                            stringsAsFactors = FALSE) %>% 
  mutate(Corr = map2_dbl(Species, Metabolite, \(p, f) corr_mat[p, f]),
         Pval = map2_dbl(Species, Metabolite, \(p, f) p_mat[p, f])) %>% 
  arrange(Pval, desc(abs(Corr))) %>% 
  drop_na(Pval)

head(corr_nad_synthase_meta_res, 10)
#                      Species                                         Metabolite      Corr         Pval
# 1         Tyzzerella_nexilis Tauroursodeoxycholic_acid_Taurohyodeoxycholic_acid 0.8788717 3.407493e-07
# 2          Absiella_dolichum Tauroursodeoxycholic_acid_Taurohyodeoxycholic_acid 0.7846600 4.197598e-05
# 3      Ruminococcus_lactaris                                      Glutamic_acid 0.7510043 1.356489e-04
# 4    Clostridium_perfringens                                   Taurocholic_acid 0.7411126 1.851048e-04
# 5           Blautia_hansenii                               Ursodeoxycholic_acid 0.7254763 2.945642e-04
# 6   Collinsella_intestinalis                               Ursodeoxycholic_acid 0.7254763 2.945642e-04
# 7      Megasphaera_stantonii                               Ursodeoxycholic_acid 0.7254763 2.945642e-04
# 8  Veillonellaceae_bacterium                               Ursodeoxycholic_acid 0.7254763 2.945642e-04
# 9      Ruminococcus_callidus                               Indolepropionic_acid 0.7007965 5.779932e-04
# 10   Roseburia_inulinivorans                         Taurochenodeoxycholic_acid 0.6977642 6.250705e-04

corr_nad_synthase_meta_file <- corr_nad_synthase_meta_res %>% 
  left_join(nad_synthase_wide %>% 
              summarise(across(-SampleID, ~sum(.x, na.rm = TRUE))) %>% 
              pivot_longer(cols = everything(),
                           names_to = "Species",
                           values_to = "SpeAbund"), by = "Species") %>% 
  left_join(meta_long %>% 
              rename(MetAbund = Abundance), by = "Metabolite") %>% 
  relocate(SpeAbund, .after = Species) %>% 
  relocate(MetAbund, .after = Metabolite); corr_nad_synthase_meta_file

# write.csv(corr_nad_synthase_meta_file, file = "Data/260523 Spearman - nad synthase (glutamine-hydrolyzing).csv", row.names = FALSE)

# ---------- 2-3. NAD(+) synthase ----------
nad_synthase2_wide <- ec_valid %>% 
  filter(EC == "6.3.1.5: NAD(+) synthase" & !is.na(Species)) %>% 
  select(-EC, -Genus) %>%
  column_to_rownames("Species") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleID"); nad_synthase2_wide

# Correlation data frame
nad_synthase2_meta_df <- nad_synthase2_wide %>% 
  left_join(meta_wide, by = "SampleID") %>%
  column_to_rownames("SampleID") %>% 
  as.matrix(); nad_synthase2_meta_df

# Correlation
corr_nad_synthase2_meta <- Hmisc::rcorr(nad_synthase2_meta_df, type = "spearman")

# Result table
nad_synthase2_var <- colnames(nad_synthase2_wide)[-1]; nad_synthase2_var
meta_var <- colnames(meta_wide)[-1]; meta_var

corr_mat <- corr_nad_synthase2_meta$r
p_mat <- corr_nad_synthase2_meta$P

corr_nad_synthase2_meta_res <- expand.grid(Species = nad_synthase2_var,
                                          Metabolite = meta_var,
                                          stringsAsFactors = FALSE) %>% 
  mutate(Corr = map2_dbl(Species, Metabolite, \(p, f) corr_mat[p, f]),
         Pval = map2_dbl(Species, Metabolite, \(p, f) p_mat[p, f])) %>% 
  arrange(Pval, desc(abs(Corr))) %>% 
  drop_na(Pval)

head(corr_nad_synthase2_meta_res, 10)
#                     Species            Metabolite       Corr         Pval
# 1   Fusobacterium_nucleatum         Acetaminophen  0.8096692 1.524539e-05
# 2     Lactobacillus_ruminis              Dopamine  0.7534448 1.253622e-04
# 3  Fusobacterium_mortiferum  Ursodeoxycholic_acid  0.7254763 2.945642e-04
# 4     Veillonella_infantium      Glycocholic_acid  0.7011544 5.726415e-04
# 5   Streptococcus_cristatus Taurolithocholic_acid  0.6907253 7.470216e-04
# 6      Enterococcus_faecium Taurodeoxycholic_acid  0.6825835 9.126286e-04
# 7         Eggerthella_lenta              Butyrate -0.6690119 1.257406e-03
# 8     Mogibacterium_timidum         Acetaminophen  0.6544086 1.744752e-03
# 9     Staphylococcus_aureus  Ursodeoxycholic_acid  0.6491103 1.956794e-03
# 10      Leuconostoc_citreum Taurolithocholic_acid  0.6388723 2.427889e-03

corr_nad_synthase2_meta_file <- corr_nad_synthase2_meta_res %>% 
  left_join(nad_synthase2_wide %>% 
              summarise(across(-SampleID, ~sum(.x, na.rm = TRUE))) %>% 
              pivot_longer(cols = everything(),
                           names_to = "Species",
                           values_to = "SpeAbund"), by = "Species") %>% 
  left_join(meta_long %>% 
              rename(MetAbund = Abundance), by = "Metabolite") %>% 
  relocate(SpeAbund, .after = Species) %>% 
  relocate(MetAbund, .after = Metabolite); corr_nad_synthase2_meta_file

# write.csv(corr_nad_synthase2_meta_file, file = "Data/260523 Spearman - nad synthase.csv", row.names = FALSE)

# ---------- 2-4. Nicotinamidase ----------
nicotinamidase_wide <- ec_valid %>% 
  filter(EC == "3.5.1.19: Nicotinamidase" & !is.na(Species)) %>% 
  select(-EC, -Genus) %>%
  column_to_rownames("Species") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleID"); nicotinamidase_wide

# Correlation data frame
nicotinamidase_meta_df <- nicotinamidase_wide %>% 
  left_join(meta_wide, by = "SampleID") %>%
  column_to_rownames("SampleID") %>% 
  as.matrix(); nicotinamidase_meta_df

# Correlation
corr_nicotinamidase_meta <- Hmisc::rcorr(nicotinamidase_meta_df, type = "spearman")

# Result table
nicotinamidase_var <- colnames(nicotinamidase_wide)[-1]; nicotinamidase_var
meta_var <- colnames(meta_wide)[-1]; meta_var

corr_mat <- corr_nicotinamidase_meta$r
p_mat <- corr_nicotinamidase_meta$P

corr_nicotinamidase_meta_res <- expand.grid(Species = nicotinamidase_var,
                                           Metabolite = meta_var,
                                           stringsAsFactors = FALSE) %>% 
  mutate(Corr = map2_dbl(Species, Metabolite, \(p, f) corr_mat[p, f]),
         Pval = map2_dbl(Species, Metabolite, \(p, f) p_mat[p, f])) %>% 
  arrange(Pval, desc(abs(Corr))) %>% 
  drop_na(Pval)

head(corr_nicotinamidase_meta_res, 10)
#                       Species                                         Metabolite       Corr        Pval
# 1       Staphylococcus_aureus                               Ursodeoxycholic_acid  0.6491103 0.001956794
# 2         Lactobacillus_sakei                                   Deoxycholic_acid  0.6245997 0.003239482
# 3      Streptococcus_gordonii                                           Butyrate  0.6144073 0.003947528
# 4    Streptococcus_pneumoniae                              Taurolithocholic_acid  0.6050329 0.004707501
# 5    Streptococcus_pneumoniae                            gamma_Aminobutyric_acid  0.5988042 0.005276308
# 6      Streptococcus_gordonii                                         Propionate  0.5853801 0.006696018
# 7     Bifidobacterium_bifidum                                           Valerate -0.5770875 0.007719704
# 8        Streptococcus_oralis                                    Lithocholi_acid  0.5662416 0.009248169
# 9       Staphylococcus_aureus Tauroursodeoxycholic_acid_Taurohyodeoxycholic_acid  0.5439896 0.013157302
# 10 Streptococcus_vestibularis                                         Tryptophan  0.5439788 0.013159475

corr_nicotinamidase_meta_file <- corr_nicotinamidase_meta_res %>% 
  left_join(nicotinamidase_wide %>% 
              summarise(across(-SampleID, ~sum(.x, na.rm = TRUE))) %>% 
              pivot_longer(cols = everything(),
                           names_to = "Species",
                           values_to = "SpeAbund"), by = "Species") %>% 
  left_join(meta_long %>% 
              rename(MetAbund = Abundance), by = "Metabolite") %>% 
  relocate(SpeAbund, .after = Species) %>% 
  relocate(MetAbund, .after = Metabolite); corr_nicotinamidase_meta_file

# write.csv(corr_nicotinamidase_meta_file, file = "Data/260523 Spearman - nicotinamidase.csv", row.names = FALSE)

# ---------- 2-5. Acetate CoA-transferase ----------
acetate_transferase_wide <- ec_valid %>% 
  filter(EC == "2.8.3.8: Acetate CoA-transferase" & !is.na(Species)) %>% 
  select(-EC, -Genus) %>%
  column_to_rownames("Species") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleID"); acetate_transferase_wide

# Correlation data frame
acetate_transferase_meta_df <- acetate_transferase_wide %>% 
  left_join(meta_wide, by = "SampleID") %>%
  column_to_rownames("SampleID") %>% 
  as.matrix(); acetate_transferase_meta_df

# Correlation
corr_acetate_transferase_meta <- Hmisc::rcorr(acetate_transferase_meta_df, type = "spearman")

# Result table
acetate_transferase_var <- colnames(acetate_transferase_wide)[-1]; acetate_transferase_var
meta_var <- colnames(meta_wide)[-1]; meta_var

corr_mat <- corr_acetate_transferase_meta$r
p_mat <- corr_acetate_transferase_meta$P

corr_acetate_transferase_meta_res <- expand.grid(Species = acetate_transferase_var,
                                            Metabolite = meta_var,
                                            stringsAsFactors = FALSE) %>% 
  mutate(Corr = map2_dbl(Species, Metabolite, \(p, f) corr_mat[p, f]),
         Pval = map2_dbl(Species, Metabolite, \(p, f) p_mat[p, f])) %>% 
  arrange(Pval, desc(abs(Corr))) %>% 
  drop_na(Pval)

head(corr_acetate_transferase_meta_res, 10)
#                             Species                                         Metabolite      Corr         Pval
# 1             Ruminococcus_lactaris                                      Glutamic_acid 0.7549569 0.0001193303
# 2             Megasphaera_stantonii                               Ursodeoxycholic_acid 0.7254763 0.0002945642
# 3             Megasphaera_stantonii Tauroursodeoxycholic_acid_Taurohyodeoxycholic_acid 0.6079883 0.0044558872
# 4             Klebsiella_pneumoniae                                      Acetaminophen 0.5851135 0.0067271051
# 5             Ruminococcus_lactaris                                             Indole 0.5770875 0.0077197039
# 6            Escherichia_fergusonii Tauroursodeoxycholic_acid_Taurohyodeoxycholic_acid 0.5680934 0.0089710796
# 7              Streptococcus_oralis                                    Lithocholi_acid 0.5614020 0.0100053991
# 8                Klebsiella_oxytoca                              Chenodeoxycholic_acid 0.5569178 0.0107512446
# 9                Klebsiella_oxytoca                              Taurolithocholic_acid 0.5440613 0.0131428441
# 10 Streptococcus_dysgalactiae_group                                      Acetaminophen 0.5406549 0.0138435858

corr_acetate_transferase_meta_file <- corr_acetate_transferase_meta_res %>% 
  left_join(acetate_transferase_wide %>% 
              summarise(across(-SampleID, ~sum(.x, na.rm = TRUE))) %>% 
              pivot_longer(cols = everything(),
                           names_to = "Species",
                           values_to = "SpeAbund"), by = "Species") %>% 
  left_join(meta_long %>% 
              rename(MetAbund = Abundance), by = "Metabolite") %>% 
  relocate(SpeAbund, .after = Species) %>% 
  relocate(MetAbund, .after = Metabolite); corr_acetate_transferase_meta_file

# write.csv(corr_acetate_transferase_meta_file, file = "Data/260523 Spearman - acetate-coa transferase.csv", row.names = FALSE)

# ---------- 2-6. Cinnamoyl-CoA ----------
cinnamoyl_coa_wide <- ec_valid %>% 
  filter(EC == "2.8.3.17: Cinnamoyl-CoA:phenyllactate CoA-transferase" & !is.na(Species)) %>% 
  select(-EC, -Genus) %>%
  column_to_rownames("Species") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleID"); cinnamoyl_coa_wide

# Correlation data frame
cinnamoyl_coa_meta_df <- cinnamoyl_coa_wide %>% 
  left_join(meta_wide, by = "SampleID") %>%
  column_to_rownames("SampleID") %>% 
  as.matrix(); cinnamoyl_coa_meta_df

# Correlation
corr_cinnamoyl_coa_meta <- Hmisc::rcorr(cinnamoyl_coa_meta_df, type = "spearman")

# Result table
cinnamoyl_coa_var <- colnames(cinnamoyl_coa_wide)[-1]; cinnamoyl_coa_var
meta_var <- colnames(meta_wide)[-1]; meta_var

corr_mat <- corr_cinnamoyl_coa_meta$r
p_mat <- corr_cinnamoyl_coa_meta$P

corr_cinnamoyl_coa_meta_res <- expand.grid(Species = cinnamoyl_coa_var,
                                                 Metabolite = meta_var,
                                                 stringsAsFactors = FALSE) %>% 
  mutate(Corr = map2_dbl(Species, Metabolite, \(p, f) corr_mat[p, f]),
         Pval = map2_dbl(Species, Metabolite, \(p, f) p_mat[p, f])) %>% 
  arrange(Pval, desc(abs(Corr))) %>% 
  drop_na(Pval)

head(corr_cinnamoyl_coa_meta_res, 10)
#                        Species              Metabolite       Corr        Pval
# 1            Eggerthella_lenta                Butyrate -0.6637480 0.001417808
# 2        Mogibacterium_timidum           Acetaminophen  0.6544086 0.001744752
# 3       Mogibacterium_diversum           Glutamic_acid -0.5646321 0.009494633
# 4  Peptostreptococcus_stomatis gamma_Aminobutyric_acid  0.5589841 0.010402153
# 5            Eggerthella_lenta   Taurolithocholic_acid -0.5332614 0.015468229
# 6            Eggerthella_lenta      Indole_lactic_acid -0.5245093 0.017585662
# 7            Eggerthella_lenta              Propionate -0.5223942 0.018130616
# 8            Eggerthella_lenta                 Acetate -0.5177849 0.019365023
# 9        Clostridium_citroniae          Nicotinic_acid -0.4935676 0.026994322
# 10           Eggerthella_lenta        Xanthurenic_acid -0.4870558 0.029404011

corr_cinnamoyl_coa_meta_file <- corr_cinnamoyl_coa_meta_res %>% 
  left_join(cinnamoyl_coa_wide %>% 
              summarise(across(-SampleID, ~sum(.x, na.rm = TRUE))) %>% 
              pivot_longer(cols = everything(),
                           names_to = "Species",
                           values_to = "SpeAbund"), by = "Species") %>% 
  left_join(meta_long %>% 
              rename(MetAbund = Abundance), by = "Metabolite") %>% 
  relocate(SpeAbund, .after = Species) %>% 
  relocate(MetAbund, .after = Metabolite); corr_cinnamoyl_coa_meta_file

# write.csv(corr_cinnamoyl_coa_meta_file, file = "Data/260523 Spearman - cinnamoyl-coa.csv", row.names = FALSE)

# ---------- 2-7. Envfit PCoA ----------
# EC significant
ec_list2 <- c("6.3.5.1: NAD(+) synthase (glutamine-hydrolyzing)", "6.3.1.5: NAD(+) synthase",
              "3.5.1.19: Nicotinamidase", "2.8.3.8: Acetate CoA-transferase", "2.8.3.17: Cinnamoyl-CoA:phenyllactate CoA-transferase")

ec_sig2 <- ec_valid %>% 
  filter(is.na(Genus) & is.na(Species) & EC %in% ec_list2) %>% 
  select(-Genus, -Species) %>% 
  group_by(EC) %>% 
  slice(1) %>% 
  column_to_rownames("EC") %>% 
  t() %>% 
  as.data.frame(); ec_sig2


# Envfit calculation
ec_meta_bac_df2 <- cbind(meta_sig, bac_sig2, ec_sig2); ec_meta_bac_df2

set.seed(42)
ef_ec_meta_bac2 <- envfit(pcoa_bc$points, 
                          ec_meta_bac_df2,
                          permutations = 999); ef_ec_meta_bac2
#                                                           Dim1     Dim2     r2 Pr(>r)    
# Butyrate                                              -0.86913 -0.49458 0.7786  0.001 ***
# Acetate                                               -0.85939 -0.51131 0.9566  0.001 ***
# Propionate                                            -0.80711 -0.59040 0.9120  0.001 ***
# Indole_acetic_acid                                     0.36299  0.93179 0.0848  0.423    
# Indole_lactic_acid                                    -0.85527  0.51818 0.1375  0.248    
# Indolepropionic_acid                                  -0.78075 -0.62484 0.2193  0.130    
# Nicotinic_acid                                        -0.94349  0.33139 0.4610  0.006 ** 
# Coprococcus_comes                                     -0.98825 -0.15285 0.1593  0.174    
# Eubacterium_rectale                                   -0.98297 -0.18374 0.1502  0.256    
# Dorea_longicatena                                     -0.94143 -0.33721 0.1805  0.160    
# Phocaeicola_vulgatus                                   0.66164  0.74982 0.0297  0.683    
# Roseburia_inulinivorans                               -0.39865  0.91710 0.0713  0.509    
# Eggerthella_lenta                                      0.67802 -0.73504 0.2912  0.085 .  
# 2.8.3.17: Cinnamoyl-CoA:phenyllactate CoA-transferase  0.54530 -0.83824 0.1055  0.349    
# 2.8.3.8: Acetate CoA-transferase                       0.99891 -0.04658 0.1652  0.186    
# 3.5.1.19: Nicotinamidase                              -0.34015 -0.94037 0.4282  0.053 .  
# 6.3.1.5: NAD(+) synthase                              -0.19062 -0.98166 0.1515  0.219    
# 6.3.5.1: NAD(+) synthase (glutamine-hydrolyzing)      -0.53479  0.84499 0.1342  0.274

# Vector coordinates 추출
vec_ec_meta_bac2 <- scores(ef_ec_meta_bac2, display = "vectors")

vec_ec_meta_bac2 <- as.data.frame(vec_ec_meta_bac2) %>% 
  rename(PCoA1 = Dim1,
         PCoA2 = Dim2) %>% 
  rownames_to_column("Feature") %>% 
  mutate(pval = ef_ec_meta_bac2$vectors$pvals,
         r2 = ef_ec_meta_bac2$vectors$r,
         group = case_when(str_detect(Feature, "Indole") ~ "Tryptophan metabolites",
                           Feature == "Nicotinic_acid" ~ "Tryptophan metabolites",
                           Feature %in% c("Butyrate", "Acetate", "Propionate") ~ "SCFA",
                           Feature %in% ec_list2 ~ "EC",
                           TRUE ~ "Bacteria"),
         label = str_replace(Feature, "_", " ")) %>% 
  mutate(label = case_when(label == "Indolepropionic acid" ~ "Indole propionic acid",
                           label == "Indole lactic_acid" ~ "Indole lactic acid",
                           label == "Indole acetic_acid" ~ "Indole acetic acid",
                           label == "6.3.5.1: NAD(+) synthase (glutamine-hydrolyzing)" ~ "NAD(+) synthase (glutamine-hydrolyzing)",
                           label == "6.3.1.5: NAD(+) synthase" ~ "NAD(+) synthase",
                           label == "3.5.1.19: Nicotinamidase" ~ "Nicotinamidase",
                           label == "2.8.3.8: Acetate CoA-transferase" ~ "Acetate CoA-transferase",
                           label == "2.8.3.17: Cinnamoyl-CoA:phenyllactate CoA-transferase" ~ "Cinnamoyl-CoA:phenyllactate CoA-transferase",
                           TRUE ~ label))

# Arrow scaling (optional)
arrow_mult <- 0.6

vec_ec_meta_bac2 <- vec_ec_meta_bac2 %>%
  mutate(PCoA1 = PCoA1 * arrow_mult,
         PCoA2 = PCoA2 * arrow_mult)

### PCoA plot
p_pcoa_ec_meta_bac2 <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, fill = Group)) +
  
  # envfit vectors
  geom_segment(data = vec_ec_meta_bac2,
               aes(x = 0, y = 0, xend = PCoA1, yend = PCoA2),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 0.8,
               alpha = 0.4,
               color = "black") +
  
  geom_text_repel(data = vec_ec_meta_bac2,
                  aes(x = PCoA1, y = PCoA2, label = label, color = group),
                  inherit.aes = FALSE,
                  size = 3.5,
                  fontface = "bold",
                  segment.color = NA,
                  show.legend = TRUE,
                  
                  box.padding = 0.5,      # label끼리 거리
                  point.padding = 0.5,    # 점/벡터 끝과 거리
                  force = 3,              # 밀어내는 힘
                  max.overlaps = Inf,     # 겹쳐도 삭제 안함
                  min.segment.length = 0, # 짧은 연결선도 허용
                  seed = 123) +           # 위치 고정
  
  # point coordinating
  geom_point(shape  = 21,
             size   = 3.5,
             color  = "white",
             stroke = 0.9) +
  
  # PERMANOVA 결과 annotation
  annotate("text", x = Inf, y = -0.32,
           hjust = 1.05, vjust = 1.3,
           label = sprintf("PERMANOVA: R² = %.4f, %s", r2_bc, fmt_pval(pv_bc)),
           size = 3.5) +
  
  scale_fill_manual(name = "TRG", values  = group_colors) +
  scale_color_manual(name = "",
                     values = c("SCFA" = "#6FA67C",
                                "Tryptophan metabolites" = "#6D7FA1",
                                "EC" = "#9A84B8",
                                "Bacteria" = "#C98C5A")) +
  
  # geom_hline(yintercept = 0, linetype = "dashed", color = "gray85", linewidth = 0.3) +
  # geom_vline(xintercept = 0, linetype = "dashed", color = "gray85", linewidth = 0.3) +
  
  labs(x = sprintf("PCoA1 (%.1f%%)", pcoa1_var),
       y = sprintf("PCoA2 (%.1f%%)", pcoa2_var),
       title = "Metabolite & Bacteria & Enzyme (Spearman)") +
  
  theme_classic2() +
  theme(aspect.ratio = 1,
        legend.position  = "bottom",
        legend.title     = element_text(face = "bold"),
        axis.title       = element_text(face = "bold")) +
  
  guides(fill = guide_legend(override.aes = list(color = NA)),
         color = guide_legend(override.aes = list(linewidth = 1.2, alpha = 1))); p_pcoa_ec_meta_bac2

p_pcoa_ec_meta_bac <- ggarrange(p_pcoa_ec_meta_bac1, p_pcoa_ec_meta_bac2); p_pcoa_ec_meta_bac

# ==================== 3. Heatmap ====================

# ---------- 3-1. Significant enzymes ----------
### 1. 데이터 준비 및 전처리 (언더바 제거 포함)
ec_before_sig <- ec_stat %>% filter(Wilcox < 0.05)
wilcox_met <- ec_before_sig %>% select(EC, Wilcox) 

# 1) 원본 컬럼명으로 metabolite_before_valid에서 데이터 추출
target_ec_orig <- wilcox_met$EC
group_var1     <- factor(metabolite_before_valid[[1]],
                         levels = c("CR", "nonCR"))
group_var2     <- factor(metabolite_before_valid[[2]],
                         levels = c("early", "advanced"))

ec_df <- ec_wide %>% 
  select(SampleID,all_of(target_ec_orig)) %>%
  column_to_rownames("SampleID")
head(ec_df)

# 2) 시각화를 위해 모든 이름의 언더바("_")를 스페이스(" ")로 치환
colnames(ec_df) <- str_remove(colnames(ec_df), "^\\d+(\\.\\d+)*:\\s*")
wilcox_met$EC   <- str_remove(wilcox_met$EC, "^\\d+(\\.\\d+)*:\\s*")

# 4) 대사물질별(Column-wise) Z-score 스케일링 후 전치(Transpose)
ec_scaled <- scale(ec_df) 
heat_ec   <- t(ec_scaled) # 결과: 행(대사물질) x 열(샘플)

### 2. 샘플 그룹화 및 Clustering 순서 정렬
sample_info <- data.frame(TRG_1 = group_var1,
                          Tstage = group_var2,
                          row.names = colnames(heat_ec))

idx_CR    <- which(sample_info$TRG_1 == "CR")
idx_nonCR <- which(sample_info$TRG_1 == "nonCR")

heat_CR    <- heat_mat[, idx_CR,    drop = FALSE]
heat_nonCR <- heat_mat[, idx_nonCR, drop = FALSE]

# 그룹 내 계층적 군집화
col_order_CR    <- colnames(heat_CR)[hclust(dist(t(heat_CR)))$order]
col_order_nonCR <- colnames(heat_nonCR)[hclust(dist(t(heat_nonCR)))$order]

final_order <- c(col_order_CR, col_order_nonCR)

heat_final     <- heat_ec[, final_order, drop = FALSE]
annotation_col <- sample_info[final_order, , drop = FALSE]

### 3. Annotation 설정 (TRG_1 테두리 및 수평 레전드 추가)
ann_colors = list(
  TRG_1 = c(CR = "#F7D9BC", nonCR = "#80461B"),
  Tstage = c(early = "#1b9e77", advanced = "#d95f02")
)

# TRG_1 레전드를 직접 생성하여 완벽하게 왼쪽에 강제 배치하기 위해 
# 기본 Annotation에서는 레전드를 숨김(show_legend = FALSE) 처리합니다.
ha <- HeatmapAnnotation(df                   = annotation_col,
                        col                  = ann_colors,
                        annotation_name_side = "left",
                        show_legend          = c(TRG_1 = FALSE, Tstage = FALSE),     
                        gp                   = gpar(col = "black", lwd = 0.3)) # TRG_1 타일에 검은색 테두리 추가

### 4. 메인 히트맵 (Z-score) 생성
col_fun_main <- colorRamp2(c(-2, 0, 2),
                           c("#4575b4", "#F8F8F8", "#d73027"))

set.seed(123)
ht_main <- Heatmap(matrix                   = heat_final,
                   name                     = "z-score",
                   col                      = col_fun_main,
                   top_annotation           = ha,
                   cluster_columns          = FALSE,
                   cluster_rows             = TRUE,
                   clustering_distance_rows = "manhattan",
                   clustering_method_rows   = "ward.D2",
                   show_column_names        = FALSE,
                   show_row_names           = TRUE,
                   row_names_gp             = gpar(fontsize = 10),
                   row_names_max_width      = unit(12, "cm"),
                   border                   = TRUE,
                   rect_gp                  = gpar(col = "black", lwd = 0.3),
                   
                   # Z-score 레전드 수평 설정
                   heatmap_legend_param = list(
                     title            = "z-score",
                     title_position   = "topcenter",
                     legend_direction = "horizontal",
                     at               = c(-2, -1, 0, 1, 2),
                     labels           = c("-2", "-1", "0", "1", "2"))) # column_title 인자 삭제 완료
ht_main

### 5. p-value 타일 히트맵 생성 (오류 완벽 수정본)
# 핵심 수정: ht_main의 순서를 강제로 뽑지 않고, 
# 데이터 행렬(heat_final)의 '초기 행 순서'를 그대로 사용합니다.
base_ec <- rownames(heat_final)

# p-value 데이터 전처리
pval_data <- wilcox_met %>%
  mutate(log10_p  = -log10(Wilcox),
         sig_star = case_when(Wilcox < 0.005 ~ "**",
                              Wilcox < 0.01 ~ "*",
                              TRUE          ~ "")) %>%
  column_to_rownames("EC")

# base_mets (초기 순서)와 100% 동일한 순서로 p-value 행렬 생성
# 이제 ht_main의 행렬과 pval_mat의 행렬 순서가 완벽히 일치합니다.
pval_mat <- matrix(pval_data[base_ec, "log10_p"], ncol = 1, dimnames = list(base_ec, "-log10(p)"))
sig_text <- matrix(pval_data[base_ec, "sig_star"], ncol = 1, dimnames = list(base_ec, "-log10(p)"))

pval_range   <- range(pval_mat, na.rm = TRUE)
col_fun_pval <- colorRamp2(c(pval_range[1], mean(pval_range), pval_range[2]),
                           c("white", "#6BAED6", "#08306B"))

ht_pval <- Heatmap(matrix          = pval_mat,
                   name            = "-log10(p)",
                   col             = col_fun_pval,
                   cluster_rows    = FALSE, # 군집화는 ht_main이 알아서 통제하므로 FALSE
                   cluster_columns = FALSE,
                   width           = unit(0.5, "cm"),
                   show_row_names  = TRUE,  # 최종 대사물질 이름은 여기서 출력
                   row_names_gp    = gpar(fontsize = 10),
                   row_names_max_width      = unit(12, "cm"),
                   # row_order = met_order <--- 이부분을 삭제했습니다! (자동 정렬 유도)
                   
                   # logP 레전드 수평 설정
                   heatmap_legend_param = list(
                     title            = "-log10(p)",
                     title_position   = "topcenter",
                     legend_direction = "horizontal"),
                   
                   cell_fun = function(j, i, x, y, w, h, fill) {
                     if (!is.na(fill)) {
                       grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "grey60", lwd = 0.4))
                     } else {
                       grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "grey60", lwd = 0.4))
                     }
                     if (!is.na(sig_text[i, j]) && sig_text[i, j] != "") {
                       grid.text(sig_text[i, j], x = x, y = y, 
                                 gp = gpar(fontsize = 10, fontface = "bold", col = "black"))
                     }
                   }
)
ht_pval

### 6. 커스텀 레전드 배치 및 최종 출력
# TRG_1 레전드를 최좌측에 강제 배치하기 위해 수동으로 생성
lgd_trg <- Legend(
  title            = "TRG_1",
  title_position   = "topcenter",
  direction        = "horizontal",
  legend_gp        = gpar(fill = c("#F7D9BC", "#80461B")),
  labels           = c("CR", "nonCR"),
  border           = "black" # 레전드 내부 박스 테두리
)

lgd_stage <- Legend(
  title            = "Tstage",
  title_position   = "topcenter",
  direction        = "horizontal",
  legend_gp        = gpar(fill = c("#1b9e77", "#d95f02")),
  labels           = c("early", "advanced"),
  border           = "black" # 레전드 내부 박스 테두리
)

# 플롯 병합 및 그리기
# merge_legend = TRUE 를 사용하여 하단에 모든 레전드가 가로로 정렬되게 함
draw(
  ht_main + ht_pval,
  heatmap_legend_side  = "bottom",
  annotation_legend_side = "bottom",
  merge_legend         = TRUE,
  heatmap_legend_list  = list(lgd_trg, lgd_stage), # TRG_1 커스텀 레전드를 리스트 최우선으로 삽입 (최좌측 배치)
  padding              = unit(c(5, 10, 5, 5), "mm")
)

# 파일 저장
svg("Figure/7-11. Enzyme heatmap samples.svg", 
    width = 10, height = 8)
draw(ht_main + ht_pval,
     heatmap_legend_side = "bottom", 
     merge_legend = TRUE,
     heatmap_legend_list = list(lgd_trg),
     padding = unit(c(5, 10, 5, 5), "mm"))
dev.off()

# ---------- 3-2. Target enzymes ----------
### 1. 데이터 준비 및 전처리 (언더바 제거 포함)
ec_before_sig <- ec_stat %>% filter(EC %in% ec_target_names)
wilcox_met <- ec_before_sig %>% select(EC, Wilcox) 

# 1) 원본 컬럼명으로 metabolite_before_valid에서 데이터 추출
target_ec_orig <- wilcox_met$EC
group_var1     <- factor(metabolite_before_valid[[1]],
                         levels = c("CR", "nonCR"))
group_var2     <- factor(metabolite_before_valid[[2]],
                         levels = c("early", "advanced"))

ec_df <- ec_wide %>% 
  select(SampleID,all_of(target_ec_orig)) %>%
  column_to_rownames("SampleID")
head(ec_df)

# 2) 시각화를 위해 모든 이름의 언더바("_")를 스페이스(" ")로 치환
colnames(ec_df) <- str_remove(colnames(ec_df), "^\\d+(\\.\\d+)*:\\s*")
wilcox_met$EC   <- str_remove(wilcox_met$EC, "^\\d+(\\.\\d+)*:\\s*")

# 3) 효소별(Column-wise) Z-score 스케일링 후 전치(Transpose)
ec_scaled <- scale(ec_df) 
heat_ec   <- t(ec_scaled) # 결과: 행(대사물질) x 열(샘플)

### 2. 샘플 그룹화 및 Clustering 순서 정렬
sample_info <- data.frame(TRG_1 = group_var1,
                          Tstage = group_var2,
                          row.names = colnames(heat_ec))

idx_CR    <- which(sample_info$TRG_1 == "CR")
idx_nonCR <- which(sample_info$TRG_1 == "nonCR")

heat_CR    <- heat_mat[, idx_CR,    drop = FALSE]
heat_nonCR <- heat_mat[, idx_nonCR, drop = FALSE]

# 그룹 내 계층적 군집화
col_order_CR    <- colnames(heat_CR)[hclust(dist(t(heat_CR)))$order]
col_order_nonCR <- colnames(heat_nonCR)[hclust(dist(t(heat_nonCR)))$order]

final_order <- c(col_order_CR, col_order_nonCR)

heat_final     <- heat_ec[, final_order, drop = FALSE]
annotation_col <- sample_info[final_order, , drop = FALSE]

### 3. Annotation 설정 (TRG_1 테두리 및 수평 레전드 추가)
ann_colors = list(
  TRG_1 = c(CR = "#F7D9BC", nonCR = "#80461B"),
  Tstage = c(early = "#1b9e77", advanced = "#d95f02")
)

# TRG_1 레전드를 직접 생성하여 완벽하게 왼쪽에 강제 배치하기 위해 
# 기본 Annotation에서는 레전드를 숨김(show_legend = FALSE) 처리합니다.
ha <- HeatmapAnnotation(df                   = annotation_col,
                        col                  = ann_colors,
                        annotation_name_side = "left",
                        show_legend          = c(TRG_1 = FALSE, Tstage = FALSE),     
                        gp                   = gpar(col = "black", lwd = 0.3)) # TRG_1 타일에 검은색 테두리 추가

### 4. 메인 히트맵 (Z-score) 생성
col_fun_main <- colorRamp2(c(-2, 0, 2),
                           c("#4575b4", "#F8F8F8", "#d73027"))

set.seed(123)
ht_main <- Heatmap(matrix                   = heat_final,
                   name                     = "z-score",
                   col                      = col_fun_main,
                   top_annotation           = ha,
                   cluster_columns          = FALSE,
                   cluster_rows             = TRUE,
                   clustering_distance_rows = "manhattan",
                   clustering_method_rows   = "ward.D2",
                   show_column_names        = FALSE,
                   show_row_names           = TRUE,
                   row_names_gp             = gpar(fontsize = 10),
                   row_names_max_width      = unit(12, "cm"),
                   border                   = TRUE,
                   rect_gp                  = gpar(col = "black", lwd = 0.3),
                   
                   # Z-score 레전드 수평 설정
                   heatmap_legend_param = list(
                     title            = "z-score",
                     title_position   = "topcenter",
                     legend_direction = "horizontal",
                     at               = c(-2, -1, 0, 1, 2),
                     labels           = c("-2", "-1", "0", "1", "2"))) # column_title 인자 삭제 완료
ht_main

### 5. p-value 타일 히트맵 생성 (오류 완벽 수정본)
# 핵심 수정: ht_main의 순서를 강제로 뽑지 않고, 
# 데이터 행렬(heat_final)의 '초기 행 순서'를 그대로 사용합니다.
base_ec <- rownames(heat_final)

# p-value 데이터 전처리
pval_data <- wilcox_met %>%
  mutate(log10_p  = -log10(Wilcox),
         sig_star = case_when(Wilcox < 0.01 ~ "**",
                              Wilcox < 0.05 ~ "*",
                              TRUE          ~ "")) %>%
  column_to_rownames("EC")

# base_mets (초기 순서)와 100% 동일한 순서로 p-value 행렬 생성
# 이제 ht_main의 행렬과 pval_mat의 행렬 순서가 완벽히 일치합니다.
pval_mat <- matrix(pval_data[base_ec, "log10_p"], ncol = 1, dimnames = list(base_ec, "-log10(p)"))
sig_text <- matrix(pval_data[base_ec, "sig_star"], ncol = 1, dimnames = list(base_ec, "-log10(p)"))

pval_range   <- range(pval_mat, na.rm = TRUE)
col_fun_pval <- colorRamp2(c(pval_range[1], mean(pval_range), pval_range[2]),
                           c("white", "#6BAED6", "#08306B"))

ht_pval <- Heatmap(matrix          = pval_mat,
                   name            = "-log10(p)",
                   col             = col_fun_pval,
                   cluster_rows    = FALSE, # 군집화는 ht_main이 알아서 통제하므로 FALSE
                   cluster_columns = FALSE,
                   width           = unit(0.5, "cm"),
                   show_row_names  = TRUE,  # 최종 대사물질 이름은 여기서 출력
                   row_names_gp    = gpar(fontsize = 10),
                   row_names_max_width      = unit(12, "cm"),
                   # row_order = met_order <--- 이부분을 삭제했습니다! (자동 정렬 유도)
                   
                   # logP 레전드 수평 설정
                   heatmap_legend_param = list(
                     title            = "-log10(p)",
                     title_position   = "topcenter",
                     legend_direction = "horizontal"),
                   
                   cell_fun = function(j, i, x, y, w, h, fill) {
                     if (!is.na(fill)) {
                       grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "grey60", lwd = 0.4))
                     } else {
                       grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "grey60", lwd = 0.4))
                     }
                     if (!is.na(sig_text[i, j]) && sig_text[i, j] != "") {
                       grid.text(sig_text[i, j], x = x, y = y, 
                                 gp = gpar(fontsize = 10, fontface = "bold", col = "black"))
                     }
                   }
)
ht_pval

### 6. 커스텀 레전드 배치 및 최종 출력
# TRG_1 레전드를 최좌측에 강제 배치하기 위해 수동으로 생성
lgd_trg <- Legend(
  title            = "TRG_1",
  title_position   = "topcenter",
  direction        = "horizontal",
  legend_gp        = gpar(fill = c("#F7D9BC", "#80461B")),
  labels           = c("CR", "nonCR"),
  border           = "black" # 레전드 내부 박스 테두리
)

lgd_stage <- Legend(
  title            = "Tstage",
  title_position   = "topcenter",
  direction        = "horizontal",
  legend_gp        = gpar(fill = c("#1b9e77", "#d95f02")),
  labels           = c("early", "advanced"),
  border           = "black" # 레전드 내부 박스 테두리
)

# 플롯 병합 및 그리기
# merge_legend = TRUE 를 사용하여 하단에 모든 레전드가 가로로 정렬되게 함
draw(
  ht_main + ht_pval,
  heatmap_legend_side  = "bottom",
  annotation_legend_side = "bottom",
  merge_legend         = TRUE,
  heatmap_legend_list  = list(lgd_trg, lgd_stage), # TRG_1 커스텀 레전드를 리스트 최우선으로 삽입 (최좌측 배치)
  padding              = unit(c(5, 10, 5, 5), "mm")
)

# 파일 저장
svg("Figure/7-12. Target Enzyme heatmap samples.svg", 
    width = 8, height = 4.2)
draw(ht_main + ht_pval,
     heatmap_legend_side = "bottom", 
     merge_legend = TRUE,
     heatmap_legend_list = list(lgd_trg),
     padding = unit(c(5, 10, 5, 5), "mm"))
dev.off()

# ==================== 4. Basic ====================

# ---------- 4-1. Volcano plot ----------
### Total
# 분석 결과를 담을 데이터프레임 초기화
res_before <- data.frame(
  EC = character(),
  T_test_p_value = numeric(),
  Wilcoxon_p_value = numeric(),
  CR_mean = numeric(),
  nonCR_mean = numeric(),
  CR_count = integer(),
  nonCR_count = integer(),
  Log2FC = numeric(),
  stringsAsFactors = FALSE
)

ec_valid_total <- ec_valid %>% 
  filter(is.na(Genus)) %>% 
  select(-Genus, -Species) %>% 
  group_by(EC) %>% 
  slice(1) %>%
  ungroup() %>% 
  column_to_rownames("EC") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleID") %>% 
  left_join(metabolite_before_valid %>%
              select(TRG_1, Tstage, SampleID),
            by = "SampleID") %>% 
  relocate(c(TRG_1, Tstage), .before = "SampleID"); ec_valid_total

# t-test & Wilcoxon test 
for (ec in valid_ec_rows) {
  data_CR <- ec_valid_total[[ec]][ec_valid_total$TRG_1 == "CR"]
  data_nonCR <- ec_valid_total[[ec]][ec_valid_total$TRG_1 == "nonCR"]
  
  # 그룹별 평균값과 개수 계산
  CR_mean <- mean(data_CR, na.rm = TRUE)
  nonCR_mean <- mean(data_nonCR, na.rm = TRUE)
  
  CR_count <- sum(data_CR > 0)
  nonCR_count <- sum(data_nonCR > 0)
  
  # Log2FC 계산 
  Log2FC <- log2((CR_mean) / (nonCR_mean))
  
  # 분산 계산 및 검증
  var_CR <- var(data_CR, na.rm = TRUE)
  var_nonCR <- var(data_nonCR, na.rm = TRUE)
  
  # Wilcoxon test와 t-test 수행, 분산이 0이 아닌 경우
  if (var_CR > 0 & var_nonCR > 0) {
    t_test_result <- t.test(data_CR, data_nonCR, var.equal = TRUE)
    wilcox_test_result <- wilcox.test(data_CR, data_nonCR)
    
    # 결과 저장
    res_before <- rbind(res_before, data.frame(
      EC = ec,
      T_test_p_value = t_test_result$p.value,
      Wilcoxon_p_value = wilcox_test_result$p.value,
      CR_mean = CR_mean,
      nonCR_mean = nonCR_mean,
      CR_count = CR_count,
      nonCR_count = nonCR_count,
      Log2FC = Log2FC,
      stringsAsFactors = FALSE
    ))
  } else {
    # 분산이 0인 경우 NA로 처리
    res_before <- rbind(res_before, data.frame(
      EC = ec,
      T_test_p_value = NA,
      Wilcoxon_p_value = NA,
      CR_mean = CR_mean,
      nonCR_mean = nonCR_mean,
      CR_count = CR_count,
      nonCR_count = nonCR_count,
      Log2FC = Log2FC,
      stringsAsFactors = FALSE
    ))
  }
}

res_before %>% 
  mutate(merged_p = sqrt(T_test_p_value*Wilcoxon_p_value)) %>% 
  arrange(merged_p) %>% 
  filter(merged_p < 0.2) %>% 
  select(EC, merged_p, CR_mean, nonCR_mean, Log2FC)  

#                                                                              EC     merged_p      CR_mean   nonCR_mean     Log2FC
# 1                                    3.6.1.31: Phosphoribosyl-ATP diphosphatase 0.0008193257 160.06336641 124.59361754  0.3614130
# 2                                 2.3.1.28: Chloramphenicol O-acetyltransferase 0.0027269003  33.23798627  57.10455364 -0.7807728
# 3                                1.1.1.17: Mannitol-1-phosphate 5-dehydrogenase 0.0101585058  28.76194930  63.38759042 -1.1400390
# 4             3.2.1.96: Mannosyl-glycoprotein endo-beta-N-acetylglucosaminidase 0.0106615090  10.70943817  31.74765600 -1.5677673
# 5                  2.3.1.180: Beta-ketoacyl-[acyl-carrier-protein] synthase III 0.0110948710 176.70627884 154.35122411  0.1951364
# 6                                        1.1.1.103: L-threonine 3-dehydrogenase 0.0116136912  45.77314570  87.94107291 -0.9420357
# 7                    3.1.7.2: Guanosine-3',5'-bis(diphosphate) 3'-diphosphatase 0.0145512045  12.44662748  27.87089538 -1.1630045
# 8                                               2.7.1.85: Beta-glucoside kinase 0.0147232089  29.66908717  14.33692677  1.0492247
# 9                                             3.8.1.3: Haloacetate dehalogenase 0.0207781997   8.43043014   3.28746181  1.3586321
# 10                               4.1.1.83: 4-hydroxyphenylacetate decarboxylase 0.0221050567  44.18423927  74.79767945 -0.7594617

res_before %>% 
  arrange(Wilcoxon_p_value) %>% 
  filter(Wilcoxon_p_value < 0.1)

#                                                                            EC T_test_p_value Wilcoxon_p_value      CR_mean   nonCR_mean CR_count nonCR_count     Log2FC
# 1                                  3.6.1.31: Phosphoribosyl-ATP diphosphatase   0.0008414227      0.000797809 160.06336641 124.59361754       11           9  0.3614130
# 2                               2.3.1.28: Chloramphenicol O-acetyltransferase   0.0017541406      0.004239105  33.23798627  57.10455364       11           9 -0.7807728
# 3                2.3.1.180: Beta-ketoacyl-[acyl-carrier-protein] synthase III   0.0165401852      0.007442248 176.70627884 154.35122411       11           9  0.1951364
# 4                                           3.8.1.3: Haloacetate dehalogenase   0.0445417520      0.009692784   8.43043014   3.28746181       11           9  1.3586321
# 5           3.2.1.96: Mannosyl-glycoprotein endo-beta-N-acetylglucosaminidase   0.0091173063      0.012467254  10.70943817  31.74765600       11           9 -1.5677673
# 6                              1.1.1.17: Mannitol-1-phosphate 5-dehydrogenase   0.0064916377      0.015896642  28.76194930  63.38759042       11           9 -1.1400390
# 7                                             2.7.1.85: Beta-glucoside kinase   0.0136363944      0.015896642  29.66908717  14.33692677       11           9  1.0492247

-log10(0.03) # 1.522879
-log10(0.02) # 1.69897

# Volcano plot 
p_ec_volcano <- res_before %>% 
  filter(!is.na(Wilcoxon_p_value)) %>% 
  mutate(EC = str_remove(EC, "^\\d+(\\.\\d+)*:\\s*"),
         sig = ifelse(-log10(Wilcoxon_p_value) >= 1.5, 
                      "sig", "ns"),
         direction = ifelse(Log2FC > 0, 
                            "increase", "decrease"),
         sig_dir = paste0(sig, "_", direction), 
         sig_dir = ifelse(grepl("sig", sig_dir), sig_dir, "ns"),
         sig_dir = case_when(sig_dir == "sig_increase" ~ "pCR enriched",
                             sig_dir == "sig_decrease" ~ "non-pCR enriched",
                             TRUE ~ sig_dir),
         ID = ifelse(sig == "sig", EC, NA)) %>% 
  ggplot(aes(x = Log2FC, y = -log10(Wilcoxon_p_value))) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 1.5, linetype = "dotted") +
  geom_point(aes(color = sig_dir), size = rel(2.5)) +
  geom_text_repel(aes(label = ID), size = rel(4)) +
  theme_classic() +
  scale_x_continuous(limits = c(-9, 9), breaks = seq(-9, 9, by = 3)) +
  scale_color_manual(values = c("pCR enriched" = "#7BB9D0",
                                "non-pCR enriched" = "#E36E65",
                                ns = "gray70")) +
  theme(aspect.ratio = 1, 
        legend.position = "bottom") +
  labs(x = "Log2 Fold Change", y = "-Log10 P-value"); p_ec_volcano

ggsave("Figure/7-13. EC volcano.svg", device = "svg", 
       plot = p_ec_volcano, width = 6, height = 6)

### Species
# 분석 결과를 담을 데이터프레임 초기화
res_s_before <- data.frame()

# Wilcoxon test 
for (ec in valid_ec_s_cols) {
  
  data_CR <- ec_s_valid_wide[[ec]][ec_s_valid_wide$TRG_1 == "CR"]
  data_nonCR <- ec_s_valid_wide[[ec]][ec_s_valid_wide$TRG_1 == "nonCR"]
  
  CR_mean <- mean(data_CR, na.rm = TRUE)
  nonCR_mean <- mean(data_nonCR, na.rm = TRUE)
  
  CR_count <- sum(data_CR > 0, na.rm = T)
  nonCR_count <- sum(data_nonCR > 0, na.rm = T)
  
  Log2FC <- log2((CR_mean + 1e-6) /
                   (nonCR_mean + 1e-6))
  
  wilcox_test_result <- tryCatch(
    wilcox.test(data_CR, data_nonCR),
    error = function(e) NULL
  )
  
  res_s_before <- rbind(
    res_s_before,
    data.frame(
      EC = ec,
      EC_label = str_remove(ec, "^\\d+(\\.\\d+)*:\\s*"),
      Wilcoxon_p_value = ifelse(is.null(wilcox_test_result), NA, wilcox_test_result$p.value),
      CR_mean = CR_mean,
      nonCR_mean = nonCR_mean,
      CR_count = CR_count,
      nonCR_count = nonCR_count,
      Log2FC = Log2FC
    )
  )
}

res_s_before %>% 
  select(-EC_label) %>% 
  arrange(Wilcoxon_p_value) %>% 
  filter(Wilcoxon_p_value < 0.02)

#                                                                              EC Wilcoxon_p_value    CR_mean nonCR_mean CR_count nonCR_count    Log2FC
# 1                        3.8.1.3: Haloacetate dehalogenase|Ruminococcus_torques       0.00976397  3.3962940  1.1638733       11           7  1.545027
# 2                   3.5.1.44: Protein-glutamine glutaminase|Eubacterium_rectale       0.01213906 11.8320073  3.7087473       11           7  1.673691
# 3                  3.5.1.53: N-carbamoylputrescine amidase|Ruminococcus_torques       0.01213906  3.4395589  0.9064918       11           7  1.923857
# 4        4.1.2.14: 2-dehydro-3-deoxy-phosphogluconate aldolase|Blautia_wexlerae       0.01489907  4.7440188  0.9358701       11           6  2.341728
# 5  2.6.1.59: dTDP-4-amino-4,6-dideoxygalactose transaminase|Eubacterium_rectale       0.01634915  3.0161018  0.1499278        7           1  4.330336
# 6 1.7.2.2: Nitrite reductase (cytochrome; ammonia-forming)|Bacteroides_plebeius       0.01666428  0.6873693  1.4372376        4           8 -1.064140

-log10(0.02) # 1.69897

# Volcano plot 
p_ec_s_volcano <- res_s_before %>% 
  filter(!is.na(Wilcoxon_p_value)) %>% 
  mutate(sig = ifelse(-log10(Wilcoxon_p_value) >= 1.75, 
                      "sig", "ns"),
         direction = ifelse(Log2FC > 0, 
                            "increase", "decrease"),
         sig_dir = paste0(sig, "_", direction), 
         sig_dir = ifelse(grepl("sig", sig_dir), sig_dir, "ns"),
         sig_dir = case_when(sig_dir == "sig_increase" ~ "pCR enriched",
                             sig_dir == "sig_decrease" ~ "non-pCR enriched",
                             TRUE ~ sig_dir),
         ID = ifelse(sig == "sig", EC_label, NA)) %>% 
  ggplot(aes(x = Log2FC, y = -log10(Wilcoxon_p_value))) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 1.75, linetype = "dotted") +
  geom_point(aes(color = sig_dir), size = rel(2.5)) +
  geom_text_repel(aes(label = ID), size = rel(4)) +
  theme_classic() +
  scale_color_manual(values = c("pCR enriched" = "#7BB9D0",
                                "non-pCR enriched" = "#E36E65",
                                ns = "gray70")) +
  theme(aspect.ratio = 1, 
        legend.position = "bottom") +
  labs(x = "Log2 Fold Change", y = "-Log10 P-value"); p_ec_s_volcano

ggsave("Figure/7-13. EC volcano (species).svg", device = "svg", 
       plot = p_ec_s_volcano, width = 6, height = 6)

# ---------- 4-2. MA plot ----------
### Total
ma_df <- ec_valid_total %>% 
  pivot_longer(cols = -c(SampleID, TRG_1, Tstage),
               names_to = "EC",
               values_to = "Abundance") %>% 
  group_by(EC, TRG_1) %>% 
  summarise(Mean = mean(Abundance, na.rm = T), .groups = "drop") %>% 
  pivot_wider(names_from = TRG_1,
              values_from = Mean) %>% 
  mutate(Log10Mean = log10((CR + nonCR) / 2),
         Log2FC = log2(CR) - log2(nonCR)) %>% 
  left_join(res_s_before %>%
              select(EC, Wilcoxon_p_value) %>% 
              distinct(EC, .keep_all = T),
            by = "EC") %>% 
  # significance
  mutate(EC = str_remove(EC, "^\\d+(\\.\\d+)*:\\s*"),
         sig = ifelse(-log10(Wilcoxon_p_value) >= 1.0 & Log10Mean >= 0,
                      "sig", "ns"),
         direction = case_when(Log2FC > 0 ~ "CR enriched",
                               Log2FC < 0 ~ "non-CR enriched",
                               TRUE       ~ "ns"),
         sig_dir = ifelse(sig == "sig", direction, "ns"))

ma_label <- ma_df %>% 
  filter(sig == "sig") %>% 
  group_by(direction) %>% 
  slice_max(order_by = abs(Log2FC), n = 5, with_ties = F) %>% 
  ungroup()

p_ec_ma <- ma_df %>% 
  ggplot(aes(x = Log10Mean, y = Log2FC)) +
  geom_point(data = ~filter(.x, sig == "ns"),
             color = "gray70",
             alpha = 0.5,
             size = rel(2.3)) +
  geom_point(data = ~filter(.x, sig == "sig"),
             aes(color = sig_dir),
             size = rel(2.8)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_text_repel(data = ma_label, aes(label = EC), size = rel(4)) +
  scale_color_manual(values = c("CR enriched" = "#7BB9D0",
                                "non-CR enriched" = "#E36E65",
                                "ns" = "gray70")) +
  theme(aspect.ratio = 1.7, 
        legend.position = "none") +
  theme_classic() + 
  labs(x = "Log10 Mean Abundance", y = "Log2 Fold Change"); p_ec_ma

ggsave("Figure/7-14. EC ma plot.svg", device = "svg", 
       plot = p_ec_ma, width = 9, height = 6)

### Species level
ma_s_df <- ec_s_valid_wide %>% 
  pivot_longer(cols = -c(SampleID, TRG_1, Tstage),
               names_to = "EC",
               values_to = "Abundance") %>% 
  group_by(EC, TRG_1) %>% 
  summarise(Mean = mean(Abundance, na.rm = T), .groups = "drop") %>% 
  pivot_wider(names_from = TRG_1,
              values_from = Mean) %>% 
  mutate(Log10Mean = log10((CR + nonCR) / 2 + 1e-6),
         Log2FC = log2((CR + 1e-6) / (nonCR + 1e-6))) %>% 
  left_join(res_s_before %>%
              select(EC, Wilcoxon_p_value) %>% 
              distinct(EC, .keep_all = T),
            by = "EC") %>% 
  # significance
  mutate(EC = str_remove(EC, "^\\d+(\\.\\d+)*:\\s*"),
         sig = ifelse(-log10(Wilcoxon_p_value) >= 1 & Log10Mean >= 0,
                      "sig", "ns"),
         direction = case_when(Log2FC > 0 ~ "CR enriched",
                               Log2FC < 0 ~ "non-CR enriched",
                               TRUE       ~ "ns"),
         sig_dir = ifelse(sig == "sig", direction, "ns"))

ma_s_label <- ma_s_df %>% 
  filter(sig == "sig") %>% 
  group_by(direction) %>% 
  slice_max(order_by = abs(Log2FC), n = 3, with_ties = F) %>% 
  ungroup()

# Labeling 기준: top 3 for each
# Coloring 기준: Wilcoxon -log10 > 1
p_ec_s_ma <- ma_s_df %>% 
  drop_na() %>% 
  ggplot(aes(x = Log10Mean, y = Log2FC)) +
  geom_point(data = ~filter(.x, sig == "ns"),
             color = "gray70",
             alpha = 0.5,
             size = rel(2.3)) +
  geom_point(data = ~filter(.x, sig == "sig"),
             aes(color = sig_dir),
             size = rel(2.8)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_text_repel(data = ma_s_label, aes(label = EC), size = rel(4),
                  max.overlaps = Inf, force = 2, force_pull = 0.5) +
  scale_color_manual(values = c("CR enriched" = "#7BB9D0",
                                "non-CR enriched" = "#E36E65",
                                "ns" = "gray70")) +
  theme(aspect.ratio = 1.7, 
        legend.position = "none") +
  theme_classic() + 
  labs(x = "Log10 Mean Abundance", y = "Log2 Fold Change"); p_ec_s_ma

ggsave("Figure/7-14. EC ma plot (species).svg", device = "svg", 
       plot = p_ec_s_ma, width = 10, height = 7)

# ---------- 4-3. PCoA plot ----------
group_var1 <- factor(ec_valid_wide[[1]])              # TRG_1  (CR / nonCR)
group_var2 <- factor(ec_valid_wide[[2]])              # Tstage (early / advanced)
ec_df      <- ec_valid_wide[, -c(1:3), drop = FALSE]


# distance 계산
dist_ec <- vegdist(ec_df, method = "bray")

# PERMANOVA 계산
set.seed(42)

perm_ec <- adonis2(
  dist_ec ~ group_var,
  permutations = 9999,
  by           = "margin"
); perm_ec

# 주요 수치 추출
r2_ec  <- round(perm_ec["group_var", "R2"], 4)
pv_ec  <- perm_ec["group_var", "Pr(>F)"]

# 분산 설명 비율 (시각화용)
pcoa_ec <- cmdscale(dist_ec, k = 2, eig = TRUE)

eig_vals   <- pcoa_ec$eig
eig_vals   <- ifelse(eig_vals < 0, 0, eig_vals)   # 음수 고유값 보정
pcoa1_var  <- round(eig_vals[1] / sum(eig_vals) * 100, 1)
pcoa2_var  <- round(eig_vals[2] / sum(eig_vals) * 100, 1)

pcoa_df <- data.frame(
  PCoA1  = pcoa_ec$points[, 1],
  PCoA2  = pcoa_ec$points[, 2],
  Group  = group_var,
  SampleID = ec_valid_wide$SampleID
)

p_pcoa_ec <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, fill = Group)) +
  
  # stat_ellipse(geom      = "polygon",
  #              type      = "t",
  #              level     = 0.90,
  #              alpha     = 0.08,
  #              linewidth = 0.5,
  #              linetype  = "dotted") +
  
  geom_point(shape  = 21,
             size   = 3.5,
             color  = "white",
             stroke = 0.9) +
  
  geom_text_repel(aes(label = SampleID),
                  size          = 2.8,
                  color         = "gray30",
                  box.padding   = 0.4,
                  max.overlaps  = 20,
                  segment.color = "gray70",
                  segment.size  = 0.3,
                  show.legend   = FALSE) +
  
  # PERMANOVA 결과 annotation
  annotate("label", x = Inf, y = Inf,
           hjust = 1.05, vjust = 1.3,
           label = sprintf("PERMANOVA (Bray-Curtis)\nR2 = %.4f\n%s", r2_ec, fmt_pval(pv_ec)),
           size = 3.2,
           fill = "transparent",
           linewidth = 0.3,
           fontface = "plain") +
  
  scale_fill_manual(values  = group_colors) +
  
  # geom_hline(yintercept = 0, linetype = "dashed", color = "gray85", linewidth = 0.3) +
  # geom_vline(xintercept = 0, linetype = "dashed", color = "gray85", linewidth = 0.3) +
  
  labs(x = sprintf("PCoA1 (%.1f%%)", pcoa1_var),
       y = sprintf("PCoA2 (%.1f%%)", pcoa2_var)) +
  
  theme_classic2() +
  theme(aspect.ratio = 1, 
        legend.position  = "bottom",
        legend.title     = element_text(face = "bold"),
        axis.title       = element_text(face = "bold")); p_pcoa_ec

ggsave("Figure/7-15. EC pcoa.svg", device = "svg", 
       plot = p_pcoa_ec, width = 6, height = 6)

# ==================== 5. Correlation ====================

# ---------- 5-1. Wilcoxon test ----------
# Data frame
ec_s_long <- ec_s_valid %>% 
  pivot_longer(cols = all_of(metabolite_before_valid$SampleID),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1, Tstage), by = "SampleID"); ec_s_long

ec_s_stat <- ec_s_long %>% 
  group_by(Name) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  mutate(FDR = p.adjust(Wilcox, method = "BH")) %>% 
  relocate(FDR, .after = Wilcox) %>% 
  arrange(Wilcox) %>% 
  drop_na(); ec_s_stat

# write.csv(ec_s_stat, file = "Data/260527 Wilcoxon - enzyme (species).csv", row.names = FALSE)

# ---------- 5-2. Spearman test ----------
# Data frame
meta_cols <- c("TRG_1", "Tstage", "SampleID")

ec_mat <- ec_s_valid_wide %>% 
  select(-all_of(meta_cols))

# ==================== 6. Wilcoxon Test ====================

# ---------- 6-1. Total ----------
# 임시 데이터
ec_total_wide <- ec %>% 
  filter(!EC %in% c("UNMAPPED", "UNGROUPED")) %>% 
  select(EC, Genus, Species, all_of(valid_samples)) %>% 
  filter(is.na(Genus)) %>% 
  select(-Genus, -Species) %>% 
  group_by(EC) %>% 
  slice(1) %>%
  ungroup() %>%
  column_to_rownames("EC") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(
    metabolite_before_valid %>%
      select(
        SampleID,
        TRG_1,
        Tstage
      ),
    by = "SampleID"
  ) %>%
  
  relocate(
    TRG_1,
    Tstage,
    .before = SampleID
  ); ec_total_wide

prevalence <- colSums(
  ec_total_wide[, -(1:3)] > 0,
  na.rm = TRUE
)

valid_ec_cols <- names(
  prevalence[prevalence >= 10]
)

ec_total_long <- ec_total_wide %>%
  select(-TRG_1, -Tstage) %>% 
  column_to_rownames("SampleID") %>% 
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("EC") %>% 
  filter(EC %in% valid_ec_cols) %>% 
  pivot_longer(
    cols = all_of(valid_samples),
    names_to = "SampleID",
    values_to = "Abundance"
  ) %>% 
  left_join(
    metabolite_before_valid %>%
      select(
        SampleID,
        TRG_1,
        Tstage
      ),
    by = "SampleID"
  ) %>%
  
  relocate(
    TRG_1,
    Tstage,
    .before = SampleID
  ); ec_total_long

# Wilcoxon 계산
ec_total_stat <- ec_total_wide %>% 
  pivot_longer(cols = all_of(valid_ec_cols),
               names_to = "EC",
               values_to = "Abundance") %>% 
  group_by(EC) %>% 
  summarise(
    
    # Wilcoxon test
    Wilcoxon_p_value = wilcox.test(
      Abundance[TRG_1 == "CR"],
      Abundance[TRG_1 == "nonCR"]
    )$p.value,
    
    # Group abundance
    CR_abundance = mean(
      Abundance[TRG_1 == "CR"],
      na.rm = TRUE
    ),
    
    nonCR_abundance = mean(
      Abundance[TRG_1 == "nonCR"],
      na.rm = TRUE
    ),
    
    # Total mean abundance
    Mean_abundance = mean(
      Abundance,
      na.rm = TRUE
    )
  ) %>% 
  mutate(
    # Fold change
    Log2FC = log2(
      (CR_abundance + 1e-6) /
        (nonCR_abundance + 1e-6)
    )
  ) %>% 
  arrange(Wilcoxon_p_value); ec_total_stat

# write.csv(ec_total_stat, file = "Data/260527 Wilcoxon - total enzyme.csv", row.names = FALSE)

ec_total_stat_mean1.0 <- ec_total_stat %>% 
  filter(Mean_abundance >= 1.0); ec_total_stat_mean1.0

# write.csv(ec_total_stat_mean1.0, file = "Data/260527 Wilcoxon - total enzyme - mean abundance 1.csv", row.names = FALSE)

# ---------- 6-2. Stratified ----------
# 임시데이터
ec_strat_tmp <- ec %>% 
  
  select(
    EC,
    Genus,
    Species,
    all_of(valid_samples)
  ) %>% 
  
  filter(
    !is.na(Genus),
    !EC %in% c("UNGROUPED", "UNMAPPED"),
    !grepl("^\\d+(\\.\\d+)*: Deleted entry$", EC)
  ) %>%
  
  mutate(
    Name = paste(EC, Species, sep = "|")
  ) %>% 
  
  relocate(Name, .before = EC); ec_strat_tmp

# Prevalence 필터링
ec_strat_valid <- ec_strat_tmp %>% 
  
  filter(
    EC %in% valid_ec_cols
  )

# Wide matrix
ec_strat_wide <- ec_strat_valid %>% 
  
  select(-EC, -Genus, -Species) %>% 
  
  column_to_rownames("Name") %>% 
  
  t() %>% 
  
  as.data.frame() %>% 
  
  rownames_to_column("SampleID") %>% 
  
  left_join(
    metabolite_before_valid %>%
      select(
        SampleID,
        TRG_1,
        Tstage
      ),
    
    by = "SampleID"
  ) %>% 
  
  relocate(
    TRG_1,
    Tstage,
    .before = SampleID
  ); ec_strat_wide[1:6, 1:6]

# Long format 변환
ec_strat_long <- ec_strat_wide %>% 
  select(-TRG_1, -Tstage) %>% 
  column_to_rownames("SampleID") %>% 
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("EC") %>% 
  pivot_longer(
    cols = all_of(valid_samples),
    names_to = "SampleID",
    values_to = "Abundance"
  ) %>% 
  left_join(
    metabolite_before_valid %>%
      select(
        SampleID,
        TRG_1,
        Tstage
      ),
    by = "SampleID"
  ) %>%
  
  relocate(
    TRG_1,
    Tstage,
    .before = SampleID
  ); ec_strat_long

# Statistics
ec_strat_stat_tmp <- ec_strat_long %>% 
  
  group_by(EC) %>% 
  
  summarise(
    
    # Wilcoxon
    Wilcoxon_p_value = wilcox.test(
      Abundance[TRG_1 == "CR"],
      Abundance[TRG_1 == "nonCR"]
    )$p.value,
    
    # Group abundance
    CR_abundance = mean(
      Abundance[TRG_1 == "CR"],
      na.rm = TRUE
    ),
    
    nonCR_abundance = mean(
      Abundance[TRG_1 == "nonCR"],
      na.rm = TRUE
    ),
    
    # Overall abundance
    Mean_abundance = mean(
      Abundance,
      na.rm = TRUE
    )
    
  ) %>% 
  
  mutate(
    
    # Fold change
    Log2FC = log2(
      (CR_abundance + 1e-6) /
        (nonCR_abundance + 1e-6)
    )
    
  ) %>% 
  
  arrange(Wilcoxon_p_value); ec_strat_stat_tmp

ec_strat_stat <- ec_strat_stat %>% 
  
  separate(
    EC,
    into = c("EC", "Species"),
    sep = "\\|"
  ) %>% 
  drop_na(Wilcoxon_p_value); ec_strat_stat

# write.csv(ec_strat_stat2, file = "Data/260527 Wilcoxon - stratified enzyme.csv", row.names = FALSE)

ec_strat_stat_mean1.0 <- ec_strat_stat2 %>% 
  filter(Mean_abundance >= 1.0)

# write.csv(ec_strat_stat_mean1.0, file = "Data/260527 Wilcoxon - stratified enzyme total mean abund 1.csv", row.names = FALSE)

# ---------- 6-3. Volcano plot ----------
ec_total_stat_mean1.0 %>%  
  arrange(Wilcoxon_p_value) %>% 
  filter(Wilcoxon_p_value < 0.02)

#   EC                                                                Wilcoxon_p_value CR_abundance nonCR_abundance Mean_abundance Log2FC
# 1 3.6.1.31: Phosphoribosyl-ATP diphosphatase                                0.000798       160.            125.           144.    0.361
# 2 2.3.1.28: Chloramphenicol O-acetyltransferase                             0.00424         33.2            57.1           44.0  -0.781
# 3 2.3.1.180: Beta-ketoacyl-[acyl-carrier-protein] synthase III              0.00744        177.            154.           167.    0.195
# 4 3.8.1.3: Haloacetate dehalogenase                                         0.00969          8.43            3.29           6.12  1.36 
# 5 3.2.1.96: Mannosyl-glycoprotein endo-beta-N-acetylglucosaminidase         0.0125          10.7            31.7           20.2  -1.57 
# 6 1.1.1.17: Mannitol-1-phosphate 5-dehydrogenase                            0.0159          28.8            63.4           44.3  -1.14 
# 7 2.7.1.85: Beta-glucoside kinase                                           0.0159          29.7            14.3           22.8   1.05 
# 8 3.2.1.122: Maltose-6'-phosphate glucosidase                               0.0159          18.1            31.2           24.0  -0.789

-log10(0.02) # 1.69897

p_ec_total_volcano <- ec_total_stat_mean1.0 %>% 
  mutate(EC = str_remove(EC, "^\\d+(\\.\\d+)*:\\s*"),
         sig = ifelse(-log10(Wilcoxon_p_value) >= 1.5, "sig", "ns"),
         direction = ifelse(Log2FC > 0, "increase", "decrease"),
         sig_dir = paste0(sig, "_", direction), 
         sig_dir = ifelse(grepl("sig", sig_dir), sig_dir, "ns"),
         sig_dir = case_when(sig_dir == "sig_increase" ~ "pCR enriched",
                             sig_dir == "sig_decrease" ~ "non-pCR enriched",
                             TRUE ~ sig_dir),
         ID = ifelse(sig == "sig", EC, NA)) %>% 
  ggplot(aes(x = Log2FC, y = -log10(Wilcoxon_p_value))) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 1.5, linetype = "dotted") +
  geom_point(aes(color = sig_dir), size = rel(2.5)) +
  geom_text_repel(aes(label = ID), size = rel(4)) +
  theme_classic() +
  scale_x_continuous(limits = c(-9, 9), breaks = seq(-9, 9, by = 3)) +
  scale_color_manual(values = c("pCR enriched" = "#7BB9D0",
                                "non-pCR enriched" = "#E36E65",
                                ns = "gray70")) +
  theme(aspect.ratio = 1, 
        legend.position = "bottom") +
  labs(title = "Total Enzymes (Abundance ≥ 1)",
       x = "Log2 Fold Change", y = "-Log10 P-value"); p_ec_total_volcano

ggsave("Figure/7-13. EC volcano (total).jpg", device = "jpg", 
       plot = p_ec_total_volcano, width = 8, height = 8)

# ---------- 6-4. MA plot ----------
ec_strat_stat_mean1.0 %>%  
  arrange(Wilcoxon_p_value) %>% 
  filter(Wilcoxon_p_value < 0.004)

-log10(0.004) # 2.39794

p_ec_strat_ma <- ec_strat_stat_mean1.0 %>% 
  mutate(EC = str_remove(EC, "^\\d+(\\.\\d+)*:\\s*"),
         Name = paste0(EC, "| ", Species),
         sig = ifelse(-log10(Wilcoxon_p_value) >= 1.0, "sig", "ns"),
         sig2 = ifelse(-log10(Wilcoxon_p_value) >= 2.5, "sig", "ns"),
         direction = ifelse(Log2FC > 0, "increase", "decrease"),
         sig_dir = paste0(sig, "_", direction), 
         sig_dir = ifelse(grepl("sig", sig_dir), sig_dir, "ns"),
         sig_dir = case_when(sig_dir == "sig_increase" ~ "pCR enriched",
                             sig_dir == "sig_decrease" ~ "non-pCR enriched",
                             TRUE ~ sig_dir),
         ID = ifelse(sig2 == "sig", Name, NA),
         Log10Mean = log10(Mean_abundance)) %>% 
  ggplot(aes(x = Log10Mean, y = Log2FC)) +
  geom_point(data = ~filter(.x, sig == "ns"),
             color = "gray70",
             alpha = 0.5,
             size = rel(2.3)) +
  geom_point(data = ~filter(.x, sig == "sig"),
             aes(color = sig_dir),
             size = rel(2.8),
             alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_text_repel(aes(label = ID), size = rel(4), max.overlaps = Inf) +
  scale_color_manual(values = c("pCR enriched" = "#7BB9D0",
                                "non-pCR enriched" = "#E36E65",
                                "ns" = "gray70")) +
  theme(aspect.ratio = 1.7, 
        legend.position = "none") +
  theme_classic() + 
  labs(title = "Stratified Enzymes (Abundance ≥ 1)",
       x = "Log10 Mean Abundance", y = "Log2 Fold Change"); p_ec_strat_ma

ggsave("Figure/7-14. EC ma (stratified).jpg", device = "jpg", 
       plot = p_ec_strat_ma, width = 10, height = 6)

# ---------- 6-5. Significant heatmap ----------
### 1. 데이터 준비 및 전처리 (언더바 제거 포함)
ec_before_sig <- ec_total_stat_mean1.0 %>% filter(Wilcoxon_p_value < 0.05) %>% rename(Wilcox = Wilcoxon_p_value); ec_before_sig
wilcox_met <- ec_before_sig %>% select(EC, Wilcox) 

# 1) 원본 컬럼명으로 metabolite_before_valid에서 데이터 추출
target_ec_orig <- wilcox_met$EC
group_var0     <- factor(metabolite_before_valid[[3]],
                         levels = c("0", "1", "2"))
group_var1     <- factor(metabolite_before_valid[[1]],
                         levels = c("CR", "nonCR"))
group_var2     <- factor(metabolite_before_valid[[2]],
                         levels = c("early", "advanced"))

ec_df <- ec_wide %>% 
  select(SampleID,all_of(target_ec_orig)) %>%
  column_to_rownames("SampleID")
head(ec_df)

# 2) 시각화를 위해 모든 이름의 언더바("_")를 스페이스(" ")로 치환
colnames(ec_df) <- str_remove(colnames(ec_df), "^\\d+(\\.\\d+)*:\\s*")
wilcox_met$EC   <- str_remove(wilcox_met$EC, "^\\d+(\\.\\d+)*:\\s*")

# 4) 효소별(Column-wise) Z-score 스케일링 후 전치(Transpose)
ec_scaled <- scale(ec_df) 
heat_ec   <- t(ec_scaled) # 결과: 행(대사물질) x 열(샘플)

### 2. 샘플 그룹화 및 Clustering 순서 정렬
sample_info <- data.frame(TRG_score = group_var0,
                          TRG_1 = group_var1,
                          Tstage = group_var2,
                          row.names = colnames(heat_ec))

idx_CR    <- which(sample_info$TRG_1 == "CR")
idx_nonCR <- which(sample_info$TRG_1 == "nonCR")

heat_CR    <- heat_mat[, idx_CR,    drop = FALSE]
heat_nonCR <- heat_mat[, idx_nonCR, drop = FALSE]

# 그룹 내 계층적 군집화
col_order_CR    <- colnames(heat_CR)[hclust(dist(t(heat_CR)))$order]
col_order_nonCR <- colnames(heat_nonCR)[hclust(dist(t(heat_nonCR)))$order]

final_order <- c(col_order_CR, col_order_nonCR)

heat_final     <- heat_ec[, final_order, drop = FALSE]
annotation_col <- sample_info[final_order, , drop = FALSE]

### 3. Annotation 설정 (TRG_1 테두리 및 수평 레전드 추가)
ann_colors = list(
  TRG_1 = c(CR = "#F7D9BC", nonCR = "#80461B"),
  Tstage = c(early = "#1b9e77", advanced = "#d95f02"),
  TRG_score = c("0" = "white", "1" = "grey50", "2" = "grey20")
)

# TRG_1 레전드를 직접 생성하여 완벽하게 왼쪽에 강제 배치하기 위해 
# 기본 Annotation에서는 레전드를 숨김(show_legend = FALSE) 처리합니다.
ha <- HeatmapAnnotation(df                   = annotation_col,
                        col                  = ann_colors,
                        annotation_name_side = "left",
                        show_legend          = c(TRG_score = FALSE, TRG_1 = FALSE, Tstage = FALSE),     
                        gp                   = gpar(col = "black", lwd = 0.3)) # TRG_1 타일에 검은색 테두리 추가

### 4. 메인 히트맵 (Z-score) 생성
col_fun_main <- colorRamp2(c(-2, 0, 2),
                           c("#4575b4", "#F8F8F8", "#d73027"))

set.seed(123)
ht_main <- Heatmap(matrix                   = heat_final,
                   name                     = "z-score",
                   col                      = col_fun_main,
                   top_annotation           = ha,
                   cluster_columns          = FALSE,
                   cluster_rows             = TRUE,
                   clustering_distance_rows = "manhattan",
                   clustering_method_rows   = "ward.D2",
                   show_column_names        = FALSE,
                   show_row_names           = TRUE,
                   row_names_gp             = gpar(fontsize = 10),
                   row_names_max_width      = unit(12, "cm"),
                   border                   = TRUE,
                   rect_gp                  = gpar(col = "black", lwd = 0.3),
                   
                   # Z-score 레전드 수평 설정
                   heatmap_legend_param = list(
                     title            = "z-score",
                     title_position   = "topcenter",
                     legend_direction = "horizontal",
                     at               = c(-2, -1, 0, 1, 2),
                     labels           = c("-2", "-1", "0", "1", "2"))) # column_title 인자 삭제 완료
ht_main

### 5. p-value 타일 히트맵 생성 (오류 완벽 수정본)
# 핵심 수정: ht_main의 순서를 강제로 뽑지 않고, 
# 데이터 행렬(heat_final)의 '초기 행 순서'를 그대로 사용합니다.
base_ec <- rownames(heat_final)

# p-value 데이터 전처리
pval_data <- wilcox_met %>%
  mutate(log10_p  = -log10(Wilcox),
         sig_star = case_when(Wilcox < 0.005 ~ "**",
                              Wilcox < 0.01 ~ "*",
                              TRUE          ~ "")) %>%
  column_to_rownames("EC")

# base_mets (초기 순서)와 100% 동일한 순서로 p-value 행렬 생성
# 이제 ht_main의 행렬과 pval_mat의 행렬 순서가 완벽히 일치합니다.
pval_mat <- matrix(pval_data[base_ec, "log10_p"], ncol = 1, dimnames = list(base_ec, "-log10(p)"))
sig_text <- matrix(pval_data[base_ec, "sig_star"], ncol = 1, dimnames = list(base_ec, "-log10(p)"))

pval_range   <- range(pval_mat, na.rm = TRUE)
col_fun_pval <- colorRamp2(c(pval_range[1], mean(pval_range), pval_range[2]),
                           c("white", "#6BAED6", "#08306B"))

ht_pval <- Heatmap(matrix          = pval_mat,
                   name            = "-log10(p)",
                   col             = col_fun_pval,
                   cluster_rows    = FALSE, # 군집화는 ht_main이 알아서 통제하므로 FALSE
                   cluster_columns = FALSE,
                   width           = unit(0.5, "cm"),
                   show_row_names  = TRUE,  # 최종 대사물질 이름은 여기서 출력
                   row_names_gp    = gpar(fontsize = 10),
                   row_names_max_width      = unit(12, "cm"),
                   # row_order = met_order <--- 이부분을 삭제했습니다! (자동 정렬 유도)
                   
                   # logP 레전드 수평 설정
                   heatmap_legend_param = list(
                     title            = "-log10(p)",
                     title_position   = "topcenter",
                     legend_direction = "horizontal"),
                   
                   cell_fun = function(j, i, x, y, w, h, fill) {
                     if (!is.na(fill)) {
                       grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "grey60", lwd = 0.4))
                     } else {
                       grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "grey60", lwd = 0.4))
                     }
                     if (!is.na(sig_text[i, j]) && sig_text[i, j] != "") {
                       grid.text(sig_text[i, j], x = x, y = y, 
                                 gp = gpar(fontsize = 10, fontface = "bold", col = "black"))
                     }
                   }
)
ht_pval

### 6. 커스텀 레전드 배치 및 최종 출력
# TRG_1 레전드를 최좌측에 강제 배치하기 위해 수동으로 생성
lgd_trg <- Legend(
  title            = "TRG_1",
  title_position   = "topcenter",
  direction        = "horizontal",
  legend_gp        = gpar(fill = c("#F7D9BC", "#80461B")),
  labels           = c("CR", "nonCR"),
  border           = "black" # 레전드 내부 박스 테두리
)

lgd_stage <- Legend(
  title            = "Tstage",
  title_position   = "topcenter",
  direction        = "horizontal",
  legend_gp        = gpar(fill = c("#1b9e77", "#d95f02")),
  labels           = c("early", "advanced"),
  border           = "black" # 레전드 내부 박스 테두리
)

lgd_score <- Legend(
  title            = "TRG_score",
  title_position   = "topcenter",
  direction        = "horizontal",
  legend_gp        = gpar(fill = c("white", "grey50", "grey20")),
  labels           = c("0", "1", "2"),
  border           = "black" # 레전드 내부 박스 테두리
)

# 플롯 병합 및 그리기
# merge_legend = TRUE 를 사용하여 하단에 모든 레전드가 가로로 정렬되게 함
draw(
  ht_main + ht_pval,
  heatmap_legend_side  = "bottom",
  annotation_legend_side = "bottom",
  merge_legend         = TRUE,
  heatmap_legend_list  = list(lgd_score, lgd_trg, lgd_stage), # TRG_1 커스텀 레전드를 리스트 최우선으로 삽입 (최좌측 배치)
  padding              = unit(c(5, 10, 5, 5), "mm")
)

# 파일 저장
svg("Figure/7-11. Enzyme heatmap samples (total).svg", 
    width = 10, height = 8)
draw(ht_main + ht_pval,
     heatmap_legend_side = "bottom", 
     merge_legend = TRUE,
     heatmap_legend_list = list(lgd_trg),
     padding = unit(c(5, 10, 5, 5), "mm"))
dev.off()

# ---------- 6-6. Target enzymes ----------
### 1. 데이터 준비 및 전처리 (언더바 제거 포함)
ec_before_sig <- ec_total_stat_mean1.0 %>% filter(EC %in% ec_target_names) %>% rename(Wilcox = Wilcoxon_p_value); ec_before_sig
wilcox_met <- ec_before_sig %>% select(EC, Wilcox) 

# 1) 원본 컬럼명으로 metabolite_before_valid에서 데이터 추출
target_ec_orig <- wilcox_met$EC
group_var1     <- factor(metabolite_before_valid[[1]],
                         levels = c("CR", "nonCR"))
group_var2     <- factor(metabolite_before_valid[[2]],
                         levels = c("early", "advanced"))

ec_df <- ec_wide %>% 
  select(SampleID,all_of(target_ec_orig)) %>%
  column_to_rownames("SampleID")
head(ec_df)

# 2) 시각화를 위해 모든 이름의 언더바("_")를 스페이스(" ")로 치환
colnames(ec_df) <- str_remove(colnames(ec_df), "^\\d+(\\.\\d+)*:\\s*")
wilcox_met$EC   <- str_remove(wilcox_met$EC, "^\\d+(\\.\\d+)*:\\s*")

# 3) 효소별(Column-wise) Z-score 스케일링 후 전치(Transpose)
ec_scaled <- scale(ec_df) 
heat_ec   <- t(ec_scaled) # 결과: 행(대사물질) x 열(샘플)

### 2. 샘플 그룹화 및 Clustering 순서 정렬
sample_info <- data.frame(TRG_1 = group_var1,
                          Tstage = group_var2,
                          row.names = colnames(heat_ec))

idx_CR    <- which(sample_info$TRG_1 == "CR")
idx_nonCR <- which(sample_info$TRG_1 == "nonCR")

heat_CR    <- heat_mat[, idx_CR,    drop = FALSE]
heat_nonCR <- heat_mat[, idx_nonCR, drop = FALSE]

# 그룹 내 계층적 군집화
col_order_CR    <- colnames(heat_CR)[hclust(dist(t(heat_CR)))$order]
col_order_nonCR <- colnames(heat_nonCR)[hclust(dist(t(heat_nonCR)))$order]

final_order <- c(col_order_CR, col_order_nonCR)

heat_final     <- heat_ec[, final_order, drop = FALSE]
annotation_col <- sample_info[final_order, , drop = FALSE]

### 3. Annotation 설정 (TRG_1 테두리 및 수평 레전드 추가)
ann_colors = list(
  TRG_1 = c(CR = "#F7D9BC", nonCR = "#80461B"),
  Tstage = c(early = "#1b9e77", advanced = "#d95f02")
)

# TRG_1 레전드를 직접 생성하여 완벽하게 왼쪽에 강제 배치하기 위해 
# 기본 Annotation에서는 레전드를 숨김(show_legend = FALSE) 처리합니다.
ha <- HeatmapAnnotation(df                   = annotation_col,
                        col                  = ann_colors,
                        annotation_name_side = "left",
                        show_legend          = c(TRG_1 = FALSE, Tstage = FALSE),     
                        gp                   = gpar(col = "black", lwd = 0.3)) # TRG_1 타일에 검은색 테두리 추가

### 4. 메인 히트맵 (Z-score) 생성
col_fun_main <- colorRamp2(c(-2, 0, 2),
                           c("#4575b4", "#F8F8F8", "#d73027"))

set.seed(123)
ht_main <- Heatmap(matrix                   = heat_final,
                   name                     = "z-score",
                   col                      = col_fun_main,
                   top_annotation           = ha,
                   cluster_columns          = FALSE,
                   cluster_rows             = TRUE,
                   clustering_distance_rows = "manhattan",
                   clustering_method_rows   = "ward.D2",
                   show_column_names        = FALSE,
                   show_row_names           = TRUE,
                   row_names_gp             = gpar(fontsize = 10),
                   row_names_max_width      = unit(12, "cm"),
                   border                   = TRUE,
                   rect_gp                  = gpar(col = "black", lwd = 0.3),
                   
                   # Z-score 레전드 수평 설정
                   heatmap_legend_param = list(
                     title            = "z-score",
                     title_position   = "topcenter",
                     legend_direction = "horizontal",
                     at               = c(-2, -1, 0, 1, 2),
                     labels           = c("-2", "-1", "0", "1", "2"))) # column_title 인자 삭제 완료
ht_main

### 5. p-value 타일 히트맵 생성 (오류 완벽 수정본)
# 핵심 수정: ht_main의 순서를 강제로 뽑지 않고, 
# 데이터 행렬(heat_final)의 '초기 행 순서'를 그대로 사용합니다.
base_ec <- rownames(heat_final)

# p-value 데이터 전처리
pval_data <- wilcox_met %>%
  mutate(log10_p  = -log10(Wilcox),
         sig_star = case_when(Wilcox < 0.05 ~ "**",
                              Wilcox < 0.1 ~ "*",
                              TRUE          ~ "")) %>%
  column_to_rownames("EC")

# base_mets (초기 순서)와 100% 동일한 순서로 p-value 행렬 생성
# 이제 ht_main의 행렬과 pval_mat의 행렬 순서가 완벽히 일치합니다.
pval_mat <- matrix(pval_data[base_ec, "log10_p"], ncol = 1, dimnames = list(base_ec, "-log10(p)"))
sig_text <- matrix(pval_data[base_ec, "sig_star"], ncol = 1, dimnames = list(base_ec, "-log10(p)"))

pval_range   <- range(pval_mat, na.rm = TRUE)
col_fun_pval <- colorRamp2(c(pval_range[1], mean(pval_range), pval_range[2]),
                           c("white", "#6BAED6", "#08306B"))

ht_pval <- Heatmap(matrix          = pval_mat,
                   name            = "-log10(p)",
                   col             = col_fun_pval,
                   cluster_rows    = FALSE, # 군집화는 ht_main이 알아서 통제하므로 FALSE
                   cluster_columns = FALSE,
                   width           = unit(0.5, "cm"),
                   show_row_names  = TRUE,  # 최종 대사물질 이름은 여기서 출력
                   row_names_gp    = gpar(fontsize = 10),
                   row_names_max_width      = unit(12, "cm"),
                   # row_order = met_order <--- 이부분을 삭제했습니다! (자동 정렬 유도)
                   
                   # logP 레전드 수평 설정
                   heatmap_legend_param = list(
                     title            = "-log10(p)",
                     title_position   = "topcenter",
                     legend_direction = "horizontal"),
                   
                   cell_fun = function(j, i, x, y, w, h, fill) {
                     if (!is.na(fill)) {
                       grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "grey60", lwd = 0.4))
                     } else {
                       grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "grey60", lwd = 0.4))
                     }
                     if (!is.na(sig_text[i, j]) && sig_text[i, j] != "") {
                       grid.text(sig_text[i, j], x = x, y = y, 
                                 gp = gpar(fontsize = 10, fontface = "bold", col = "black"))
                     }
                   }
)
ht_pval

### 6. 커스텀 레전드 배치 및 최종 출력
# TRG_1 레전드를 최좌측에 강제 배치하기 위해 수동으로 생성
lgd_trg <- Legend(
  title            = "TRG_1",
  title_position   = "topcenter",
  direction        = "horizontal",
  legend_gp        = gpar(fill = c("#F7D9BC", "#80461B")),
  labels           = c("CR", "nonCR"),
  border           = "black" # 레전드 내부 박스 테두리
)

lgd_stage <- Legend(
  title            = "Tstage",
  title_position   = "topcenter",
  direction        = "horizontal",
  legend_gp        = gpar(fill = c("#1b9e77", "#d95f02")),
  labels           = c("early", "advanced"),
  border           = "black" # 레전드 내부 박스 테두리
)

# 플롯 병합 및 그리기
# merge_legend = TRUE 를 사용하여 하단에 모든 레전드가 가로로 정렬되게 함
draw(
  ht_main + ht_pval,
  heatmap_legend_side  = "bottom",
  annotation_legend_side = "bottom",
  merge_legend         = TRUE,
  heatmap_legend_list  = list(lgd_trg, lgd_stage), # TRG_1 커스텀 레전드를 리스트 최우선으로 삽입 (최좌측 배치)
  padding              = unit(c(5, 10, 5, 5), "mm")
)

# 파일 저장
svg("Figure/7-12. Target Enzyme heatmap samples (total).svg", 
    width = 8, height = 4.2)
draw(ht_main + ht_pval,
     heatmap_legend_side = "bottom", 
     merge_legend = TRUE,
     heatmap_legend_list = list(lgd_trg),
     padding = unit(c(5, 10, 5, 5), "mm"))
dev.off()

# ---------- 6-7. Regression scatter ----------
# Data frame
bug_names <- c("Eubacterium_rectale", "Roseburia_inulinivorans", "Faecalibacterium_prausnitzii", "Coprococcus_comes", "Coprococcus_catus")
ec_names  <- c("6.3.4.21: Nicotinate phosphoribosyltransferase",
               "2.3.1.8: Phosphate acetyltransferase",
               "2.7.2.1: Acetate kinase",
               "2.8.3.8: Acetate CoA-transferase",
               "2.7.2.7: Butyrate kinase")
met_names <- c("Nicotinic_acid", "Acetate", "Butyrate")

scatter_df <- metabolite_before_valid %>% 
  select(SampleID, all_of(met_names)) %>% 
  left_join(ec_total_wide %>% # Log-transformed 값 사용
              select(SampleID, all_of(ec_names)) %>% 
              mutate(`6.3.4.21: Nicotinate phosphoribosyltransferase` = log10(`6.3.4.21: Nicotinate phosphoribosyltransferase` + 1),
                     `2.3.1.8: Phosphate acetyltransferase` = log10(`2.3.1.8: Phosphate acetyltransferase` + 1),
                     `2.7.2.1: Acetate kinase` = log10(`2.7.2.1: Acetate kinase` + 1),
                     `2.8.3.8: Acetate CoA-transferase` = log10(`2.8.3.8: Acetate CoA-transferase` + 1),
                     `2.7.2.7: Butyrate kinase` = log10(`2.7.2.7: Butyrate kinase` + 1)),
            by = "SampleID") %>% 
  left_join(s_trans %>% # Raw 값 사용
              rownames_to_column("SampleID") %>% 
              select(SampleID, all_of(bug_names)),
            by = "SampleID") %>% 
  left_join(metabolite_before_valid %>% # Raw 값 사용
              select(SampleID, TRG_1),
            by = "SampleID"); scatter_df

# Nicotinate phosphoribosyltransferase - Eubacterium_rectale - Nicotinic_acid
scale_factor <- max(scatter_df$`6.3.4.21: Nicotinate phosphoribosyltransferase`, na.rm = T) /
  max(scatter_df$Nicotinic_acid, na.rm = T); scale_factor

scale_factor <- scale_factor * 1.5

regression_df1 <- scatter_df %>% 
  rename(meta_value = Nicotinic_acid,
         bac_value  = Eubacterium_rectale,
         enz_log    = `6.3.4.21: Nicotinate phosphoribosyltransferase`) %>% 
  mutate(meta_scaled = meta_value * scale_factor) %>% 
  select(SampleID, TRG_1, meta_value, meta_scaled, bac_value, enz_log); head(regression_df1)

meta_rng <- range(regression_df1$meta_value, na.rm = T)
enz_rng <- range(regression_df1$enz_log, na.rm = T)

p_regression1 <- regression_df1 %>% 
  ggplot(aes(x = bac_value)) +
  
  # Enzyme
  geom_point(
    aes(y = enz_log),
    size = rel(3),
    alpha = 0.8,
    color = "seagreen"
  ) +
  
  geom_smooth(
    aes(y = enz_log),
    method = "lm",
    color = "seagreen",
    se = T,
    alpha = 0.1
  ) +
  
  # Metabolite
  geom_point(
    aes(y = meta_scaled),
    size = rel(3),
    alpha = 0.7,
    color = "coral"
  ) +
  
  geom_smooth(
    aes(y = meta_scaled),
    method = "lm",
    color = "coral",
    linetype = "dashed",
    se = T,
    alpha = 0.1
  ) +
  
  # Correlation annotation
  stat_cor(
    aes(
      y = enz_log,
      label = paste(
        "Enzyme Spearman:",
        after_stat(r.label),
        after_stat(p.label)
      )
    ),
    method = "spearman",
    output.type = "text",
    label.x.npc = 0.02,
    label.y.npc = 0.95,
    size = 3.8
  ) +
  
  stat_cor(
    aes(
      y = meta_scaled,
      label = paste(
        "Metabolite Spearman:",
        after_stat(r.label),
        after_stat(p.label)
      )
    ),
    method = "spearman",
    output.type = "text",
    label.x.npc = 0.02,
    label.y.npc = 0.90,
    size = 3.8
  ) +
  
  # scale_color_manual(
  #   values = group_colors
  # ) +
  
  # X-axis
  scale_x_continuous(
    expand = c(0.005, 0.005),
    name = "Eubacterium rectale (%)"
  ) +
  
  # Dual y-axis
  scale_y_continuous(
    expand = c(0.005, 0.005),
    
    # Y-axis (left)
    name = "Log10(Nicotinate phosphoribosyltransferase + 1)",
    
    # Y-axis (right)
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Nicotinic acid"
    )
  ) +
  
  theme_classic() +
  
  theme(
    axis.title.y = element_text(color = "seagreen"),
    axis.title.y.right = element_text(color = "coral"),
    aspect.ratio = 1,
    legend.position = "bottom"
  ) +
  
  labs(title = "Nicotinate phosphoribosyltransferase - E. rectale"); p_regression1

# Phosphate acetyltransferase - Roseburia inulinivorans - Acetate
scale_factor <- max(scatter_df$`2.3.1.8: Phosphate acetyltransferase`, na.rm = T) /
  max(scatter_df$Acetate, na.rm = T); scale_factor

scale_factor <- scale_factor * 1.5

regression_df2 <- scatter_df %>% 
  rename(meta_value = Acetate,
         bac_value  = Roseburia_inulinivorans,
         enz_log    = `2.3.1.8: Phosphate acetyltransferase`) %>% 
  mutate(meta_scaled = meta_value * scale_factor) %>% 
  select(SampleID, TRG_1, meta_value, meta_scaled, bac_value, enz_log); head(regression_df2)

meta_rng <- range(regression_df2$meta_value, na.rm = T)
enz_rng <- range(regression_df2$enz_log, na.rm = T)

p_regression2 <- regression_df2 %>% 
  ggplot(aes(x = bac_value)) +
  
  # Enzyme
  geom_point(
    aes(y = enz_log),
    size = rel(3),
    alpha = 0.8,
    color = "seagreen"
  ) +
  
  geom_smooth(
    aes(y = enz_log),
    method = "lm",
    color = "seagreen",
    se = T,
    alpha = 0.1
  ) +
  
  # Metabolite
  geom_point(
    aes(y = meta_scaled),
    size = rel(3),
    alpha = 0.7,
    color = "coral"
  ) +
  
  geom_smooth(
    aes(y = meta_scaled),
    method = "lm",
    color = "coral",
    linetype = "dashed",
    se = T,
    alpha = 0.1
  ) +
  
  # Correlation annotation
  stat_cor(
    aes(
      y = enz_log,
      label = paste(
        "Enzyme Spearman:",
        after_stat(r.label),
        after_stat(p.label)
      )
    ),
    method = "spearman",
    output.type = "text",
    label.x.npc = 0.02,
    label.y.npc = 0.95,
    size = 3.8
  ) +
  
  stat_cor(
    aes(
      y = meta_scaled,
      label = paste(
        "Metabolite Spearman:",
        after_stat(r.label),
        after_stat(p.label)
      )
    ),
    method = "spearman",
    output.type = "text",
    label.x.npc = 0.02,
    label.y.npc = 0.90,
    size = 3.8
  ) +
  
  # scale_color_manual(
  #   values = group_colors
  # ) +
  
  # X-axis
  scale_x_continuous(
    expand = c(0.005, 0.005),
    name = "Roseburia inulinivorans (%)"
  ) +
  
  # Dual y-axis
  scale_y_continuous(
    expand = c(0.005, 0.005),
    
    # Y-axis (left)
    name = "Log10(Phosphate acetyltransferase + 1)",
    
    # Y-axis (right)
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Acetate"
    )
  ) +
  
  theme_classic() +
  
  theme(
    axis.title.y = element_text(color = "seagreen"),
    axis.title.y.right = element_text(color = "coral"),
    aspect.ratio = 1,
    legend.position = "bottom"
  ) +
  
  labs(title = "Phosphate acetyltransferase - R. inulinivorans"); p_regression2

# Acetate kinase - Roseburia inulinivorans - Acetate
scale_factor <- max(scatter_df$`2.7.2.1: Acetate kinase`, na.rm = T) /
  max(scatter_df$Acetate, na.rm = T); scale_factor

scale_factor <- scale_factor * 1.5

regression_df3_1 <- scatter_df %>% 
  rename(meta_value = Acetate,
         bac_value  = Roseburia_inulinivorans,
         enz_log    = `2.7.2.1: Acetate kinase`) %>% 
  mutate(meta_scaled = meta_value * scale_factor) %>% 
  select(SampleID, TRG_1, meta_value, meta_scaled, bac_value, enz_log); head(regression_df3_1)

meta_rng <- range(regression_df3_1$meta_value, na.rm = T)
enz_rng <- range(regression_df3_1$enz_log, na.rm = T)

p_regression3_1 <- regression_df3_1 %>% 
  ggplot(aes(x = bac_value)) +
  
  # Enzyme
  geom_point(
    aes(y = enz_log),
    size = rel(3),
    alpha = 0.8,
    color = "seagreen"
  ) +
  
  geom_smooth(
    aes(y = enz_log),
    method = "lm",
    color = "seagreen",
    se = T,
    alpha = 0.1
  ) +
  
  # Metabolite
  geom_point(
    aes(y = meta_scaled),
    size = rel(3),
    alpha = 0.7,
    color = "coral"
  ) +
  
  geom_smooth(
    aes(y = meta_scaled),
    method = "lm",
    color = "coral",
    linetype = "dashed",
    se = T,
    alpha = 0.1
  ) +
  
  # Correlation annotation
  stat_cor(
    aes(
      y = enz_log,
      label = paste(
        "Enzyme Spearman:",
        after_stat(r.label),
        after_stat(p.label)
      )
    ),
    method = "spearman",
    output.type = "text",
    label.x.npc = 0.02,
    label.y.npc = 0.95,
    size = 3.8
  ) +
  
  stat_cor(
    aes(
      y = meta_scaled,
      label = paste(
        "Metabolite Spearman:",
        after_stat(r.label),
        after_stat(p.label)
      )
    ),
    method = "spearman",
    output.type = "text",
    label.x.npc = 0.02,
    label.y.npc = 0.90,
    size = 3.8
  ) +
  
  # scale_color_manual(
  #   values = group_colors
  # ) +
  
  # X-axis
  scale_x_continuous(
    expand = c(0.005, 0.005),
    name = "Roseburia inulinivorans (%)"
  ) +
  
  # Dual y-axis
  scale_y_continuous(
    expand = c(0.005, 0.005),
    
    # Y-axis (left)
    name = "Log10(Acetate kinase + 1)",
    
    # Y-axis (right)
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Acetate"
    )
  ) +
  
  theme_classic() +
  
  theme(
    aspect.ratio = 1,
    legend.position = "bottom"
  ) +
  
  labs(title = "Acetate kinase - R. inulinivorans"); p_regression3_1

# Acetate kinase - Eubacterium rectale - Acetate
scale_factor <- max(scatter_df$`2.7.2.1: Acetate kinase`, na.rm = T) /
  max(scatter_df$Acetate, na.rm = T); scale_factor

scale_factor <- scale_factor * 1.5

regression_df3_2 <- scatter_df %>% 
  rename(meta_value = Acetate,
         bac_value  = Eubacterium_rectale,
         enz_log    = `2.7.2.1: Acetate kinase`) %>% 
  mutate(meta_scaled = meta_value * scale_factor) %>% 
  select(SampleID, TRG_1, meta_value, meta_scaled, bac_value, enz_log); head(regression_df3_2)

meta_rng <- range(regression_df3_2$meta_value, na.rm = T)
enz_rng <- range(regression_df3_2$enz_log, na.rm = T)

p_regression3_2 <- regression_df3_2 %>% 
  ggplot(aes(x = bac_value)) +
  
  # Enzyme
  geom_point(
    aes(y = enz_log),
    size = rel(3),
    alpha = 0.8,
    color = "seagreen"
  ) +
  
  geom_smooth(
    aes(y = enz_log),
    method = "lm",
    color = "seagreen",
    se = T,
    alpha = 0.1
  ) +
  
  # Metabolite
  geom_point(
    aes(y = meta_scaled),
    size = rel(3),
    alpha = 0.7,
    color = "coral"
  ) +
  
  geom_smooth(
    aes(y = meta_scaled),
    method = "lm",
    color = "coral",
    linetype = "dashed",
    se = T,
    alpha = 0.1
  ) +
  
  # Correlation annotation
  stat_cor(
    aes(
      y = enz_log,
      label = paste(
        "Enzyme Spearman:",
        after_stat(r.label),
        after_stat(p.label)
      )
    ),
    method = "spearman",
    output.type = "text",
    label.x.npc = 0.02,
    label.y.npc = 0.95,
    size = 3.8
  ) +
  
  stat_cor(
    aes(
      y = meta_scaled,
      label = paste(
        "Metabolite Spearman:",
        after_stat(r.label),
        after_stat(p.label)
      )
    ),
    method = "spearman",
    output.type = "text",
    label.x.npc = 0.02,
    label.y.npc = 0.90,
    size = 3.8
  ) +
  
  # scale_color_manual(
  #   values = group_colors
  # ) +
  
  # X-axis
  scale_x_continuous(
    expand = c(0.005, 0.005),
    name = "Eubacterium rectale (%)"
  ) +
  
  # Dual y-axis
  scale_y_continuous(
    expand = c(0.005, 0.005),
    
    # Y-axis (left)
    name = "Log10(Acetate kinase + 1)",
    
    # Y-axis (right)
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Acetate"
    )
  ) +
  
  theme_classic() +
  
  theme(
    aspect.ratio = 1,
    legend.position = "bottom"
  ) +
  
  labs(title = "Acetate kinase - E. rectale"); p_regression3_2

# Acetate kinase - Faecalibacterium prausnitzii - Acetate
scale_factor <- max(scatter_df$`2.7.2.1: Acetate kinase`, na.rm = T) /
  max(scatter_df$Acetate, na.rm = T); scale_factor

scale_factor <- scale_factor * 1.5

regression_df3_3 <- scatter_df %>% 
  rename(meta_value = Acetate,
         bac_value  = Faecalibacterium_prausnitzii,
         enz_log    = `2.7.2.1: Acetate kinase`) %>% 
  mutate(meta_scaled = meta_value * scale_factor) %>% 
  select(SampleID, TRG_1, meta_value, meta_scaled, bac_value, enz_log); head(regression_df3_3)

meta_rng <- range(regression_df3_3$meta_value, na.rm = T)
enz_rng <- range(regression_df3_3$enz_log, na.rm = T)

p_regression3_3 <- regression_df3_3 %>% 
  ggplot(aes(x = bac_value)) +
  
  # Enzyme
  geom_point(
    aes(y = enz_log),
    size = rel(3),
    alpha = 0.8,
    color = "seagreen"
  ) +
  
  geom_smooth(
    aes(y = enz_log),
    method = "lm",
    color = "seagreen",
    se = T,
    alpha = 0.1
  ) +
  
  # Metabolite
  geom_point(
    aes(y = meta_scaled),
    size = rel(3),
    alpha = 0.7,
    color = "coral"
  ) +
  
  geom_smooth(
    aes(y = meta_scaled),
    method = "lm",
    color = "coral",
    linetype = "dashed",
    se = T,
    alpha = 0.1
  ) +
  
  # Correlation annotation
  stat_cor(
    aes(
      y = enz_log,
      label = paste(
        "Enzyme Spearman:",
        after_stat(r.label),
        after_stat(p.label)
      )
    ),
    method = "spearman",
    output.type = "text",
    label.x.npc = 0.02,
    label.y.npc = 0.95,
    size = 3.8
  ) +
  
  stat_cor(
    aes(
      y = meta_scaled,
      label = paste(
        "Metabolite Spearman:",
        after_stat(r.label),
        after_stat(p.label)
      )
    ),
    method = "spearman",
    output.type = "text",
    label.x.npc = 0.02,
    label.y.npc = 0.90,
    size = 3.8
  ) +
  
  # scale_color_manual(
  #   values = group_colors
  # ) +
  
  # X-axis
  scale_x_continuous(
    expand = c(0.005, 0.005),
    name = "Faecalibacterium prausnitzii (%)"
  ) +
  
  # Dual y-axis
  scale_y_continuous(
    expand = c(0.005, 0.005),
    
    # Y-axis (left)
    name = "Log10(Acetate kinase + 1)",
    
    # Y-axis (right)
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Acetate"
    )
  ) +
  
  theme_classic() +
  
  theme(
    aspect.ratio = 1,
    legend.position = "bottom"
  ) +
  
  labs(title = "Acetate kinase - F. prausnitzii"); p_regression3_3

# Acetate CoA-transferase - Coprococcus catus - Butyrate
scale_factor <- max(scatter_df$`2.8.3.8: Acetate CoA-transferase`, na.rm = T) /
  max(scatter_df$Butyrate, na.rm = T); scale_factor

scale_factor <- scale_factor * 1.5

regression_df4 <- scatter_df %>% 
  rename(meta_value = Butyrate,
         bac_value  = Coprococcus_catus,
         enz_log    = `2.8.3.8: Acetate CoA-transferase`) %>% 
  mutate(meta_scaled = meta_value * scale_factor) %>% 
  select(SampleID, TRG_1, meta_value, meta_scaled, bac_value, enz_log); head(regression_df4)

meta_rng <- range(regression_df4$meta_value, na.rm = T)
enz_rng <- range(regression_df4$enz_log, na.rm = T)

p_regression4 <- regression_df4 %>% 
  ggplot(aes(x = bac_value)) +
  
  # Enzyme
  geom_point(
    aes(y = enz_log),
    size = rel(3),
    alpha = 0.8,
    color = "seagreen"
  ) +
  
  geom_smooth(
    aes(y = enz_log),
    method = "lm",
    color = "seagreen",
    se = T,
    alpha = 0.1
  ) +
  
  # Metabolite
  geom_point(
    aes(y = meta_scaled),
    size = rel(3),
    alpha = 0.7,
    color = "coral"
  ) +
  
  geom_smooth(
    aes(y = meta_scaled),
    method = "lm",
    color = "coral",
    linetype = "dashed",
    se = T,
    alpha = 0.1
  ) +
  
  # Correlation annotation
  stat_cor(
    aes(
      y = enz_log,
      label = paste(
        "Enzyme Spearman:",
        after_stat(r.label),
        after_stat(p.label)
      )
    ),
    method = "spearman",
    output.type = "text",
    label.x.npc = 0.02,
    label.y.npc = 0.95,
    size = 3.8
  ) +
  
  stat_cor(
    aes(
      y = meta_scaled,
      label = paste(
        "Metabolite Spearman:",
        after_stat(r.label),
        after_stat(p.label)
      )
    ),
    method = "spearman",
    output.type = "text",
    label.x.npc = 0.02,
    label.y.npc = 0.90,
    size = 3.8
  ) +
  
  # scale_color_manual(
  #   values = group_colors
  # ) +
  
  # X-axis
  scale_x_continuous(
    expand = c(0.001, 0.001),
    name = "Coprococcus catus (%)"
  ) +
  
  # Dual y-axis
  scale_y_continuous(
    expand = c(0.001, 0.001),
    
    # Y-axis (left)
    name = "Log10(Acetate CoA-transferase + 1)",
    
    # Y-axis (right)
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Butyrate"
    )
  ) +
  
  theme_classic() +
  
  theme(
    aspect.ratio = 1,
    legend.position = "bottom"
  ) +
  
  labs(title = "Acetate CoA-transferase - C. catus"); p_regression4

# Acetate CoA-transferase - Coprococcus catus - Butyrate
scale_factor <- max(scatter_df$`2.7.2.7: Butyrate kinase`, na.rm = T) /
  max(scatter_df$Butyrate, na.rm = T); scale_factor

scale_factor <- scale_factor * 1.5

regression_df5 <- scatter_df %>% 
  rename(meta_value = Butyrate,
         bac_value  = Coprococcus_comes,
         enz_log    = `2.7.2.7: Butyrate kinase`) %>% 
  mutate(meta_scaled = meta_value * scale_factor) %>% 
  select(SampleID, TRG_1, meta_value, meta_scaled, bac_value, enz_log); head(regression_df5)

meta_rng <- range(regression_df5$meta_value, na.rm = T)
enz_rng <- range(regression_df5$enz_log, na.rm = T)

p_regression5 <- regression_df5 %>% 
  ggplot(aes(x = bac_value)) +
  
  # Enzyme
  geom_point(
    aes(y = enz_log),
    size = rel(3),
    alpha = 0.8,
    color = "seagreen"
  ) +
  
  geom_smooth(
    aes(y = enz_log),
    method = "lm",
    color = "seagreen",
    se = T,
    alpha = 0.1
  ) +
  
  # Metabolite
  geom_point(
    aes(y = meta_scaled),
    size = rel(3),
    alpha = 0.7,
    color = "coral"
  ) +
  
  geom_smooth(
    aes(y = meta_scaled),
    method = "lm",
    color = "coral",
    linetype = "dashed",
    se = T,
    alpha = 0.1
  ) +
  
  # Correlation annotation
  stat_cor(
    aes(
      y = enz_log,
      label = paste(
        "Enzyme Spearman:",
        after_stat(r.label),
        after_stat(p.label)
      )
    ),
    method = "spearman",
    output.type = "text",
    label.x.npc = 0.02,
    label.y.npc = 0.95,
    size = 3.8
  ) +
  
  stat_cor(
    aes(
      y = meta_scaled,
      label = paste(
        "Metabolite Spearman:",
        after_stat(r.label),
        after_stat(p.label)
      )
    ),
    method = "spearman",
    output.type = "text",
    label.x.npc = 0.02,
    label.y.npc = 0.90,
    size = 3.8
  ) +
  
  # scale_color_manual(
  #   values = group_colors
  # ) +
  
  # X-axis
  scale_x_continuous(
    expand = c(0.001, 0.001),
    name = "Coprococcus comes (%)"
  ) +
  
  # Dual y-axis
  scale_y_continuous(
    expand = c(0.001, 0.001),
    
    # Y-axis (left)
    name = "Log10(Butyrate kinase + 1)",
    
    # Y-axis (right)
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Butyrate"
    )
  ) +
  
  theme_classic() +
  
  theme(
    aspect.ratio = 1,
    legend.position = "bottom"
  ) +
  
  labs(title = "Butyrate kinase - C. comes"); p_regression5

p_regression_a <- ggarrange(p_regression1, p_regression2); p_regression_a
p_regression_b <- ggarrange(p_regression3_1, p_regression3_2, p_regression3_3, nrow = 1, ncol = 3); p_regression_b
p_regression_c <- ggarrange(p_regression4, p_regression5); p_regression_c

ggsave("Figure/7-16. EC regression scatter (1).svg", device = "svg", 
       plot = p_regression_a, width = 10, height = 5)
ggsave("Figure/7-16. EC regression scatter (2).svg", device = "svg", 
       plot = p_regression_b, width = 15, height = 5)
ggsave("Figure/7-16. EC regression scatter (3).svg", device = "svg", 
       plot = p_regression_c, width = 10, height = 5)

# ==================== 7. Spearman Correlation ====================

# ---------- 7-1. Total ----------
# EC 풍부도 데이터
ec_total_matrix <- ec_total_wide %>% 
  select(-TRG_1, -Tstage, -TRG_score) %>% 
  column_to_rownames("SampleID"); ec_total_matrix[1:6, 1:6]

# Metabolite 풍부도 데이터
metabolite_valid_matrix <- metabolite_before_valid %>% 
  select(-TRG_1) %>% 
  column_to_rownames("SampleID"); metabolite_valid_matrix[1:6, 1:6]

# Correlation input
ec_metabolite_matrix <- bind_cols(ec_total_matrix, metabolite_matrix) %>% 
  as.matrix()

# Correlation 계산
corr_s <- Hmisc::rcorr(ec_metabolite_matrix, type = "spearman")
corr_p <- Hmisc::rcorr(ec_metabolite_matrix, type = "pearson")

# 결과 정리
ec_var <- colnames(ec_total_matrix)
meta_var <- colnames(metabolite_matrix)

corr_ec_total_meta <- expand.grid(EC = ec_var,
                                  Metabolite = meta_var,
                                  stringsAsFactors = FALSE) %>% 
  mutate(Corr_s = map2_dbl(EC, Metabolite, ~ corr_s$r[.x, .y]),
         Pval_s = map2_dbl(EC, Metabolite, ~ corr_s$P[.x, .y]),
         Corr_p = map2_dbl(EC, Metabolite, ~ corr_p$r[.x, .y]),
         Pval_p = map2_dbl(EC, Metabolite, ~ corr_p$P[.x, .y])) %>% 
  mutate(Pval_merged = (Pval_s + Pval_p)/2) %>% 
  arrange(Pval_merged); corr_ec_total_meta

# write.csv(corr_ec_total_meta, file = "Data/260601 Spearman - total enzyme n metabolite.csv", row.names = FALSE)

corr_ec_total_meta0.3 <- corr_ec_total_meta %>% 
  filter(abs(Corr_s) > 0.3 & abs(Corr_p) > 0.3); corr_ec_total_meta0.3

# write.csv(corr_ec_total_meta0.3, file = "Data/260601 Spearman - total enzyme n metabolite corr 0.3.csv", row.names = FALSE)

# ---------- 7-2. Stratified ----------
# EC 풍부도 데이터
ec_strat_matrix <- ec_strat_wide %>% 
  select(-TRG_1, -Tstage, -TRG_score) %>% 
  column_to_rownames("SampleID"); ec_strat_matrix[1:6, 1:6]

# Metabolite 풍부도 데이터
metabolite_valid_matrix <- metabolite_before_valid %>% 
  select(-TRG_1) %>% 
  column_to_rownames("SampleID"); metabolite_valid_matrix[1:6, 1:6]

# Correlation input
ec_strat_meta_matrix <- bind_cols(ec_strat_matrix, metabolite_matrix) %>% 
  as.matrix()

# Correlation 계산
corr_s <- Hmisc::rcorr(ec_strat_meta_matrix, type = "spearman")
corr_p <- Hmisc::rcorr(ec_strat_meta_matrix, type = "pearson")

# 결과 정리
ec_var <- colnames(ec_total_matrix)
meta_var <- colnames(metabolite_matrix)

corr_ec_total_meta <- expand.grid(EC = ec_var,
                                  Metabolite = meta_var,
                                  stringsAsFactors = FALSE) %>% 
  mutate(Corr_s = map2_dbl(EC, Metabolite, ~ corr_s$r[.x, .y]),
         Pval_s = map2_dbl(EC, Metabolite, ~ corr_s$P[.x, .y]),
         Corr_p = map2_dbl(EC, Metabolite, ~ corr_p$r[.x, .y]),
         Pval_p = map2_dbl(EC, Metabolite, ~ corr_p$P[.x, .y])) %>% 
  mutate(Pval_merged = (Pval_s + Pval_p)/2) %>% 
  arrange(Pval_merged); corr_ec_total_meta

# write.csv(corr_ec_total_meta, file = "Data/260601 Spearman - total enzyme n metabolite.csv", row.names = FALSE)

corr_ec_total_meta0.3 <- corr_ec_total_meta %>% 
  filter(abs(Corr_s) > 0.3 & abs(Corr_p) > 0.3); corr_ec_total_meta0.3

# write.csv(corr_ec_total_meta0.3, file = "Data/260601 Spearman - total enzyme n metabolite corr 0.3.csv", row.names = FALSE)

# ==================================================
save.image(file = "RData/260526 enzyme.RData")




###------------------ JW Huh ------------------###



# canonical histidine degradation vs. ImP brach

# histidine --> urocanate --> ImP
# histidine --> urocanate --> 4-imidazolone-5-propanoate

# 1. L-histidine --> trans-urocanate
#    - histidine ammonia-lyase; hutH; 4.3.1.3; K01745
#
# 2-1. trans-urocanate --> imidazole propionate (dihydrourocanate)
#    - urocanate reductase; urdA; 1.3.99.33
# 
# 2-2. trans-urocanate --> 4-imidazolone-5-propanoate
#    - urocanate hydratase; hutU; 4.2.1.49; K01712
#
# 3-2. 4-imidazolone-5-propanoate --> N-formimino-L-glutamate
#    - imidazolonepropionase; hutI; 3.5.2.7; K01468
#
# 4-2. FIGLU --> glutamate
#    - hutG and many others
ec <- read_tsv("260224 final script/Input/merged_genefamilies_EC_named.tsv",
               show_col_types = FALSE) %>%
  rename(Feature = `# Gene Family`) %>%
  mutate(
    has_taxon = str_detect(Feature, "\\|"),
    EC = str_replace(Feature, "\\|.*$", ""),
    Taxon = if_else(
      has_taxon,
      str_replace(Feature, "^[^|]+\\|", ""),
      NA_character_
    ),
    Taxon_type = case_when(
      !has_taxon ~ "total",
      Taxon == "unclassified" ~ "unclassified",
      TRUE ~ "classified"
    ),
    Genus = if_else(
      Taxon_type == "classified",
      str_extract(Taxon, "(?<=g__)[^\\.]+"),
      NA_character_
    ),
    Species = if_else(
      Taxon_type == "classified",
      str_extract(Taxon, "(?<=s__)[^\\.]+"),
      NA_character_
    )
  ) %>%
  relocate(EC, Taxon_type, has_taxon, Taxon, Genus, Species, .after = Feature)
  # head(ec); dim(ec) ; colnames(ec)


# Before only
m_before <- m %>%
  filter(TNT == "Before") %>%
  distinct(SampleID, .keep_all = TRUE); dim(m_before) # 26 samples * 29 vars

  m_before %>% count(TRG_1) # CR 11, nonCR 15

head(m_before)
# Remove 'unmapped' and 'ungrouped'
ec_valid_tmp <- ec %>%
  filter(!EC %in% c("UNMAPPED", "UNGROUPED")) %>%
  select(EC, Taxon_type, Taxon, 
         Genus, Species, Feature, m_before$SampleID)
  # dim(ec_valid_tmp) # 173608 * 32


ec_valid_0.2 <- ec_valid_tmp[rowSums(ec_valid_tmp[,-c(1:6)]>0) >= 26*0.2,]


dim(ec_valid_0.2) # 62271 EC

### Wilcoxon

# Long-form conversion

ec_valid_long <- ec_valid_0.2 %>% 
  pivot_longer(
    
    cols = m_before$SampleID, 
    names_to = "SampleID",
    values_to = "Abundance"
    
    ) %>% 
  inner_join(m_before, by = "SampleID")
  # head(ec_valid_long, 20)
  # dim(ec_valid_long) # 1619046 * 36



# library(dplyr)
# library(tidyr)
# library(stringr)

# ec_wilcox_all <- ec_valid_long %>%
#   filter(TRG_1 %in% c("CR", "nonCR")) %>%
#   group_by(Feature, EC, Taxon_type, Taxon, Genus, Species) %>%
#   summarise(
#     n_total = n_distinct(SampleID),
#     n_CR = n_distinct(SampleID[TRG_1 == "CR"]),
#     n_nonCR = n_distinct(SampleID[TRG_1 == "nonCR"]),
#     
#     n_pos_CR = sum(TRG_1 == "CR" & Abundance > 0, na.rm = TRUE),
#     n_pos_nonCR = sum(TRG_1 == "nonCR" & Abundance > 0, na.rm = TRUE),
#     
#     prev_CR = n_pos_CR / n_CR,
#     prev_nonCR = n_pos_nonCR / n_nonCR,
#     
#     median_CR = median(Abundance[TRG_1 == "CR"], na.rm = TRUE),
#     median_nonCR = median(Abundance[TRG_1 == "nonCR"], na.rm = TRUE),
#     
#     mean_CR = mean(Abundance[TRG_1 == "CR"], na.rm = TRUE),
#     mean_nonCR = mean(Abundance[TRG_1 == "nonCR"], na.rm = TRUE),
#     
#     median_diff = median_CR - median_nonCR,
#     mean_diff = mean_CR - mean_nonCR,
#     
#     p_value = {
#       x <- Abundance[TRG_1 == "CR"]
#       y <- Abundance[TRG_1 == "nonCR"]
#       
#       x <- x[!is.na(x)]
#       y <- y[!is.na(y)]
#       
#       if (length(x) >= 2 &&
#           length(y) >= 2 &&
#           length(unique(c(x, y))) > 1) {
#         tryCatch(
#           wilcox.test(x, y, exact = FALSE)$p.value,
#           error = function(e) NA_real_
#         )
#       } else {
#         NA_real_
#       }
#     },
#     .groups = "drop"
#   ) %>%
#   mutate(
#     q_value = p.adjust(p_value, method = "BH"),
#     direction = case_when(
#       median_CR > median_nonCR ~ "CR higher",
#       median_CR < median_nonCR ~ "nonCR higher",
#       TRUE ~ "No median difference"
#     )
#   ) %>%
#   arrange(q_value, p_value)



# Total only
ec_wilcox_total <- ec_wilcox_all %>%
  filter(Taxon_type == "total")

# species-stratified
ec_wilcox_species <- ec_wilcox_all %>%
  filter(Taxon_type == "classified")

ec_wilcox_species %>% 
  arrange(p_value) %>% 
  as.data.frame() %>% 
  head(30)


# Vars to test - 1) Imidazole propionate-related
selected_1_ImP <- c("4.3.1.3",
                    "1.3.99.33",
                    "4.2.1.49",
                    "3.5.2.7")

res_selected_1_ImP <- ec_wilcox_all %>%
  filter(
    Reduce(
      `|`,
      lapply(selected_1_ImP, function(k) {
        grepl(k, Feature, fixed = TRUE, ignore.case = TRUE) |
          grepl(k, EC, fixed = TRUE, ignore.case = TRUE)
      })
    )
  )

res_selected_1_ImP %>% 
  arrange(p_value) %>% 
  # filter(p_value < 0.1, 
  #        direction != "No median difference") %>% 
  as.data.frame() %>% 
  filter(Taxon_type == "total")
  head(30)



# Vars to test - 2) 
selected_2_gtx <- c("fragilysin", 
                    "colibactin", 
                    "cytolethal",
                    "3.4.24.74",
                    "3.5.2.6")

res_selected_2_gtx <- ec_wilcox_all %>%
  filter(
    Reduce(
      `|`,
      lapply(selected_2_gtx, function(k) {
        grepl(k, Feature, fixed = TRUE, ignore.case = TRUE) |
          grepl(k, EC, fixed = TRUE, ignore.case = TRUE) |
          grepl(k, Taxon, fixed = TRUE, ignore.case = TRUE)
      })
    )
  )

