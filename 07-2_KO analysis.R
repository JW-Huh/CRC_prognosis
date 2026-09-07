# Analyzed by Nayun Kim
# CNU Translational Microbiome Lab.

# ==================== 0. Setting ====================

# ---------- 0-1. Directory setting ----------
setwd("C:/Users/user/Desktop/NAYUN/CRC TNT/")
options(java.parameters = "-Xmx64g", stringsAsFactors = F)
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
load("RData/260522 enzyme.RData")

# ---------- 0-3. Importing ----------
# Before samples
before_samples <- colnames(gb[, -1])
valid_samples <- metabolite_before_valid$SampleID

# Raw data
ko <- read_tsv("Input/merged_genefamilies_KO_named.tsv") %>% 
  separate(`# Gene Family`,
           into = c("KO", "Taxon"),
           sep = "\\|",
           fill = "right",
           extra = "merge") %>% 
  mutate(Genus  = str_extract(Taxon, "(?<=g__)[^\\.]+"),
         Species = str_extract(Taxon, "(?<=s__)[^\\.]+")) %>%
  as.data.frame() %>% 
  relocate(c(Genus, Species), .after = KO) %>% 
  select(-Taxon); head(ko)

ko_valid_tmp <- ko %>% 
  filter(!KO %in% c("UNMAPPED", "UNGROUPED")) %>% 
  select(KO, Genus, Species, all_of(valid_samples)); ko_valid_tmp

# Sorting by prevalence
# 20 before-samples * 50% prevalence = 10
ko_prev <- ko_valid_tmp %>% 
  filter(is.na(Genus)) %>% 
  select(-Genus, -Species) %>% 
  group_by(KO) %>% 
  slice(1) %>% 
  column_to_rownames("KO") %>% 
  t() %>% 
  as.data.frame(); ko_prev

prevalence <- colSums(ko_prev[, -1] > 0)
valid_ko_cols <- names(prevalence[prevalence >= 10])

# Species-level data frame (cols = samples, rows = KO)
ko_valid_stratified <- ko_valid_tmp %>% 
  filter(KO %in% valid_ko_cols); ko_valid_stratified

# Summarized data frame (cols = KO, rows = samples)
ko_valid_wide <- ko_valid_tmp %>% 
  filter(KO %in% valid_ko_cols, is.na(Genus)) %>% 
  select(-Genus, -Species) %>% 
  group_by(KO) %>% 
  slice(1) %>% 
  column_to_rownames("KO") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleID") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1, Tstage), by = "SampleID") %>% 
  relocate(c(TRG_1, Tstage), .before = "SampleID"); ko_valid_wide

ko_valid_long <- ko_valid_stratified %>% 
  filter(KO %in% valid_ko_cols, is.na(Genus)) %>% 
  select(-Genus, -Species) %>% 
  group_by(KO) %>% 
  slice(1) %>% 
  pivot_longer(cols = all_of(valid_samples),
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  left_join(metabolite_before_valid %>% select(SampleID, TRG_1, Tstage), by = "SampleID") %>% 
  relocate(c(TRG_1, Tstage), .before = "SampleID"); ko_valid_long

# ==================== 1. Diversity ====================

# ---------- 4-3. PCoA plot ----------
group_var1 <- factor(ko_valid_wide[[1]])              # TRG_1  (CR / nonCR)
group_var2 <- factor(ko_valid_wide[[2]])              # Tstage (early / advanced)
ko_df      <- ko_valid_wide[, -c(1:3), drop = FALSE]

# distance 계산
dist_ko <- vegdist(ko_df, method = "bray")

# PERMANOVA 계산
set.seed(42)
perm_ko <- adonis2(
  dist_ko ~ group_var1,
  permutations = 9999,
  by           = "margin"
); perm_ko

# 주요 수치 추출
r2_ko  <- round(perm_ko["group_var", "R2"], 4)
pv_ko  <- perm_ko["group_var", "Pr(>F)"]

# 분산 설명 비율 (시각화용)
pcoa_ko <- cmdscale(dist_ko, k = 2, eig = TRUE)

eig_vals   <- pcoa_ko$eig
eig_vals   <- ifelse(eig_vals < 0, 0, eig_vals)   # 음수 고유값 보정
pcoa1_var  <- round(eig_vals[1] / sum(eig_vals) * 100, 1)
pcoa2_var  <- round(eig_vals[2] / sum(eig_vals) * 100, 1)

pcoa_df <- data.frame(
  PCoA1  = pcoa_ko$points[, 1],
  PCoA2  = pcoa_ko$points[, 2],
  Group  = group_var,
  SampleID = ko_valid_wide$SampleID
)

p_pcoa_ko <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, fill = Group)) +
  
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
           label = sprintf("PERMANOVA (Bray-Curtis)\nR2 = %.4f\n%s", r2_ko, fmt_pval(pv_ko)),
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
        axis.title       = element_text(face = "bold")); p_pcoa_ko

ggsave("Figure/8-01. KO pcoa.svg", device = "svg", 
       plot = p_pcoa_ko, width = 6, height = 6)

# ==================== 2. Differential Enrichment ====================

# ---------- 2-1. Volcano plot ----------
# 분석 결과를 담을 데이터프레임 초기화
volcano_df <- data.frame(
  KO = character(),
  T_test_p_value = numeric(),
  Wilcoxon_p_value = numeric(),
  CR_mean = numeric(),
  nonCR_mean = numeric(),
  CR_count = integer(),
  nonCR_count = integer(),
  Log2FC = numeric(),
  stringsAsFactors = FALSE
)

# t-test & Wilcoxon test 
for (ko in valid_ko_cols) {
  data_CR <- ko_valid_wide[[ko]][ko_valid_wide$TRG_1 == "CR"]
  data_nonCR <- ko_valid_wide[[ko]][ko_valid_wide$TRG_1 == "nonCR"]
  
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
    volcano_df <- rbind(volcano_df, data.frame(
      KO = ko,
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
    volcano_df <- rbind(volcano_df, data.frame(
      KO = ko,
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

volcano_df %>% 
  mutate(merged_p = sqrt(T_test_p_value*Wilcoxon_p_value)) %>% 
  arrange(merged_p) %>% 
  filter(merged_p < 0.2) %>% 
  select(KO, merged_p, CR_mean, nonCR_mean, Log2FC)  

#                                                                       KO    merged_p      CR_mean   nonCR_mean     Log2FC
# 1                                                        K06919: NO_NAME 0.002854147  61.21987038  29.78473636  1.0394268 - putative DNA primase/helicase
# 2      K13280: signal peptidase, endoplasmic reticulum-type [EC:3.4.-.-] 0.003679655  11.57114422   4.59424729  1.3326311
# 3                          K03892: ArsR family transcriptional regulator 0.003958610   5.83692137  20.00105940 -1.7767969
# 4                               K01195: beta-glucuronidase [EC:3.2.1.31] 0.011089583   9.37415413   4.50342157  1.0576670
# 5                                 K06406: stage V sporulation protein AD 0.012000035   9.25511117   3.80057793  1.2840315
# 6                        K03613: electron transport complex protein RnfE 0.014456884  13.38689429  26.57587504 -0.9892959
# 7                                                        K05942: NO_NAME 0.014640754  10.26052711   3.34400897  1.6174542 - citrate (Re)-synthase [EC:2.3.3.3]
# 8                                K02411: flagellar assembly protein FliH 0.017288171   6.77870834   2.02036792  1.7463924
# 9                                                    K01200: pullulanase 0.017945985  16.99036449   6.94843526  1.2899568
# 10                                                       K06926: NO_NAME 0.020222181   7.75100809   2.45810956  1.6568347 - uncharacterized protein * (GenBank) Predicted ATPases

volcano_df %>% 
  arrange(Wilcoxon_p_value) %>% 
  filter(Wilcoxon_p_value < 0.01)

#                                                                  KO T_test_p_value Wilcoxon_p_value   CR_mean nonCR_mean CR_count nonCR_count    Log2FC
# 1 K13280: signal peptidase, endoplasmic reticulum-type [EC:3.4.-.-]    0.005922278      0.002286259 11.571144   4.594247       11           9  1.332631 
# 2                                                   K06919: NO_NAME    0.002591342      0.003143606 61.219870  29.784736       11           9  1.039427 - putative DNA primase/helicase
# 3                     K03892: ArsR family transcriptional regulator    0.002770561      0.005656109  5.836921  20.001059       11           9 -1.776797
# 4                                               K01200: pullulanase    0.043274339      0.007442248 16.990364   6.948435       11           9  1.289957
# 5                            K06406: stage V sporulation protein AD    0.019349105      0.007442248  9.255111   3.800578       11           9  1.284032
# 6                                                   K06926: NO_NAME    0.042189800      0.009692784  7.751008   2.458110       11           8  1.656835 - uncharacterized protein
# 7                                                 K02406: flagellin    0.102919061      0.009763970  5.224712   1.765992       11           7  1.564872

-log10(0.01) # 2

# Volcano plot 
p_ko_volcano <- volcano_df %>% 
  filter(!is.na(Wilcoxon_p_value)) %>% 
  mutate(KO = str_remove(KO, "^K\\d+:\\s*"),
         KO =  str_remove(KO, "\\s*\\[.*?\\]$"),
         sig = ifelse(-log10(Wilcoxon_p_value) >= 1.5, 
                      "sig", "ns"),
         direction = ifelse(Log2FC > 0, 
                            "increase", "decrease"),
         sig_dir = paste0(sig, "_", direction), 
         sig_dir = ifelse(grepl("sig", sig_dir), sig_dir, "ns"),
         sig_dir = case_when(sig_dir == "sig_increase" ~ "pCR enriched",
                             sig_dir == "sig_decrease" ~ "non-pCR enriched",
                             TRUE ~ sig_dir),
         ID = ifelse(sig == "sig", KO, NA)) %>% 
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
  labs(x = "Log2 Fold Change", y = "-Log10 P-value"); p_ko_volcano

ggsave("Figure/8-02. KO volcano.svg", device = "svg", 
       plot = p_ko_volcano, width = 6, height = 6)

# ---------- 2-3. Full screening ----------
# Data frame
ko_stat <- ko_valid_long %>% 
  group_by(KO) %>% 
  summarise(Wilcox = wilcox.test(Abundance[TRG_1 == "CR"], Abundance[TRG_1 == "nonCR"])$p.value,
            Abun_CR = sum(Abundance[TRG_1 == "CR"], na.rm = T),
            Abun_nonCR = sum(Abundance[TRG_1 == "nonCR"], na.rm = T)) %>% 
  arrange(Wilcox); ko_stat

# write.csv(ko_stat, file = "Data/260527 Wilcoxon - ko.csv", row.names = FALSE)
