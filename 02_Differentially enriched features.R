
# Sample Filter: TNT - Before

rm(list = ls())
options(java.parameters = "-Xmx64g", stringsAsFactors = F)
setwd("C:/Users/user/Desktop/윤채빈/Data/CRC Metagenomics")

library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggforce)
library(gridExtra)


load("input/R_image/7-2. after-taxonomy-barplot.RData")
load("input/R_image/7-3. after-differential-enrichment-testing.RData")



###############  Differentially enriched bacteria (DEB) ###############

### Identifying DEB through Wilcoxon rank-sum test 


###############  1. Strain  ###############

# Prevalence filtering (>=10%)
tb_0.1 = tb %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Strain) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.1) %>% 
  pull(Strain) %>% 
  sort() ; tb_0.1


# Make an input file for differential enrichment testing 
# Long-form conversion
# merge with metadata 
tb_input = tb %>% 
  filter(Strain %in% tb_0.1) %>% 
  pivot_longer(cols = -Strain,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  merge(mb2 %>% 
          select(SampleID,
                 TRG_score, TRG_1, TRG_2, TRG_3,
                 Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                 Pre_Op_Nstage, Pre_Op_Nstage_bin,
                 BMI, BMI_bin, CEA, CEA_bin),
        by = "SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR")))

min_abund = tb_input %>% 
  filter(abundance > 0) %>% 
  pull(abundance) %>% min() # 2e-05


# TRG_1: CR vs. nonCR
# Wilcoxon test (CR vs. nonCR)
wilcox_t = tb_input %>% 
  group_by(Strain) %>% 
  summarise(p_value = wilcox.test(abundance ~ TRG_1, exact = F)$p.value,
            .groups = "drop") %>% 
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>% 
  arrange(p_value) ; wilcox_t

sig_strain = wilcox_t %>% 
  filter(p_value < 0.05) %>% 
  arrange(p_value) %>% 
  pull(Strain)
# 49 strains, P < 0.1
# 22 strains, P < 0.05

sig_strain = tb_input %>% 
  filter(Strain %in% sig_strain) %>% 
  group_by(Strain) %>% 
  summarise(avg_abund = mean(abundance)) %>% 
  arrange(-avg_abund) %>% 
  pull(Strain) ; sig_strain

p4_TRG_1_strain = tb_input %>% 
  filter(Strain %in% sig_strain) %>% 
  mutate(Strain = factor(Strain, levels = sig_strain),
         TRG_1 = factor(TRG_1, levels = c("nonCR", "CR"))) %>% 
  ggplot(aes(TRG_1, log10(abundance + min_abund/2))) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  facet_wrap(~Strain, ncol = 2) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0) +
  theme_classic() +
  coord_flip(ylim = c(-5, 1.5)) +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(0.95)),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = rel(1.05)),
        axis.text.y = element_text(size = rel(1.3)),
        strip.background = element_rect(color = NA),
        aspect.ratio = 0.2) +
  labs(x = "Tumor Regression Grade",
       y = "Log10 (Abundance + 1e-5)") ; p4_TRG_1_strain

# ggsave("figure/04-3-1_TRG_1-strain.svg",
#        plot = p4_TRG_1_strain, width = 7, height = 12)


### Bad bugs
# "Fusobacterium_nucleatum|SGB6011",
# "Intestinimonas_butyriciproducens|SGB15126"    
# "Hungatella_hathewayi|SGB4742"
# "GGB9719_SGB15272|SGB15272"
# "Faecalicatena_contorta|SGB4614" 
# "Gemella_sanguinis|SGB7298_group"
# "Leuconostoc_gelidum|SGB7127" 

# Bad indicators reported by my previous paper (Huh et al. 2022 Microbiome)
bad_t1 = c("Dialister_invisus|SGB5825_group", 
           "Alistipes_senegalensis|SGB2296",
           "Fusobacterium_nucleatum|SGB6011",
           "Fusobacterium_nucleatum|SGB6013",
           "Fusobacterium_nucleatum|SGB6014")

# CRC signature bacteria
bad_t2 = c("Fusobacterium_nucleatum|SGB6011",
           "Fusobacterium_nucleatum|SGB6013",
           "Fusobacterium_nucleatum|SGB6014",
           "Gemella_morbillorum|SGB7295",
           "Peptostreptococcus_stomatis|SGB748",
           "Parvimonas_micra|SGB6653",
           "Hungatella_hathewayi|SGB4741",
           "Hungatella_hathewayi|SGB4742")

# Colitis-inducing bacteria 
"Klebsiella_pneumoniae|SGB10115_group"

bad_t0 = c(bad_t2, bad_t1, "Klebsiella_pneumoniae|SGB10115_group") %>% unique()

# Significant bacteria 중 CR에 많은 것 (nonCR에도 3개 이상 존재)
sig_good_t1 = c("Eubacterium_rectale|SGB4933",
                "Bifidobacterium_longum|SGB17248", 
                "Roseburia_inulinivorans|SGB4940",
                "Ruminococcus_sp_AF13_28|SGB4834", 
                "Clostridiales_bacterium_KLE1615|SGB5090", 
                "Clostridium_SGB4750|SGB4750",
                "Bifidobacterium_dentium|SGB17234",
                "GGB2653_SGB3574|SGB3574",
                "Lancefieldella_parvula|SGB966",
                "Actinomyces_graevenitzii|SGB17130_group")

# Significant bacteria 중 nonCR에 많은 것
sig_bad_t1 = c("GGB9719_SGB15272|SGB15272", 
               "Faecalicatena_contorta|SGB4614", 
               "Gemella_sanguinis|SGB7298_group", 
               "Intestinimonas_butyriciproducens|SGB15126")

# TRG_1
p4_TRG_1_strain3 = tb_input %>% 
  filter(Strain %in% bad_t0) %>% 
  mutate(Strain = factor(Strain, levels = bad_t0)) %>% 
  ggplot(aes(TRG_1, log10(abundance + min_abund/2))) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  facet_wrap(~Strain, ncol = 2) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.05)),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = rel(1.05)),
        axis.text.y = element_text(size = rel(1.3)),
        strip.background = element_rect(color = NA),
        aspect.ratio = 0.15) +
  labs(x = "Tumor Regression Grade",
       y = "Log10 (Abundance + 1e-5)") ; p4_TRG_1_strain3

# ggsave("figure/04-3-2_TRG_1-bad_known_strains.svg",
#        plot = p4_TRG_1_strain3, width = 9, height = 6.5)

p4_TRG_1_strain2 = tb_input %>% 
  filter(Strain %in% sig_bad_t1) %>% 
  mutate(Strain = factor(Strain, levels = sig_bad_t1)) %>% 
  ggplot(aes(TRG_1, log10(abundance + min_abund/2))) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  facet_wrap(~Strain, ncol = 1) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.05)),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = rel(1.05)),
        axis.text.y = element_text(size = rel(1.3)),
        strip.background = element_rect(color = NA),
        aspect.ratio = 0.15) +
  labs(x = "Tumor Regression Grade",
       y = "Log10 (Abundance + 1e-5)") ; p4_TRG_1_strain2

# ggsave("figure/04-3-3_TRG_1-bad_sig.svg",
#        plot = p4_TRG_1_strain2, width = 6, height = 7.5)

p4_TRG_1_strain4 = tb_input %>% 
  filter(Strain %in% c(sig_good_t1, sig_bad_t1)) %>% 
  mutate(Strain = factor(Strain, levels = c(sig_good_t1, sig_bad_t1))) %>% 
  ggplot(aes(TRG_1, log10(abundance + min_abund/2))) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  facet_wrap(~Strain, ncol = 2, dir = "v") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 4.5, method = "wilcox",
                     hide.ns = T, tip.length = 0) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.3)),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.6)),
        strip.background = element_rect(color = NA),
        aspect.ratio = 0.15) +
  labs(y = "Log10 (abundance + 1e-5)") ; p4_TRG_1_strain4

# ggsave("figure/04-3-4_TRG_1-strain_sig_strain_all.svg",
#        plot = p4_TRG_1_strain3, width = 8.5, height = 6)


# "Top-10 Bad bugs"
# TRG_1
sig_bad_top10 =
  c("Faecalicatena_contorta|SGB4614", 
    "Hungatella_hathewayi|SGB4742", 
    "Eisenbergiella_tayi|SGB4988",
    "GGB9719_SGB15272|SGB15272",
    "Solobacterium_moorei|SGB6826",
    "Gemella_sanguinis|SGB7298_group", 
    "Leuconostoc_gelidum|SGB7127",
    "Intestinimonas_butyriciproducens|SGB15126", 
    "GGB45613_SGB63326|SGB63326",
    "Fusobacterium_nucleatum|SGB6011")

p4_TRG_1_strain_bad = tb_input %>% 
  filter(Strain %in% sig_bad_top10) %>% 
  mutate(T_bin = ifelse(Pre_Op_Tstage >= 3, "advanced", "early"),
         Stage_response = factor(paste0(T_bin, "_", TRG_1),
                                 levels = rev(c("early_CR", "early_nonCR",
                                                "advanced_CR", "advanced_nonCR")))) %>% 
  group_by(Stage_response, TRG_1, SampleID) %>% 
  summarise(abundance = sum(abundance)) %>% 
  ggplot(aes(Stage_response, log10(abundance + min_abund/2))) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("early_CR", "early_nonCR"), 
                                        c("advanced_CR", "advanced_nonCR"),
                                        c("early_nonCR", "advanced_nonCR")),
                     size = 4.5, method = "wilcox", 
                     label.y = c(-0.9, 0.2, 0.6),
                     tip.length = 0.01,
                     hide.ns = T) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.1)),
        axis.text.x = element_text(size = rel(1.3)),
        strip.background = element_rect(color = NA)) +
  labs(y = "Log10 (abundance + 1e-5)") ; p4_TRG_1_strain_bad

# ggsave("figure/04-3-5_TRG_1-strain-bad10.svg",
#        plot = p4_TRG_1_strain_bad, width = 6, height = 2.5)


# Correlation filtering 
cor_t = tb_input %>%
  group_by(Strain) %>%
  group_modify(~{
    df_tmp = .x %>%
      select(abundance, TRG_score) %>%
      mutate(
        abundance = as.numeric(as.character(abundance)),
        TRG_score = as.numeric(as.character(TRG_score))
      ) %>%
      na.omit()
    
    if (nrow(df_tmp) >= 3 &&
        length(unique(df_tmp$abundance)) > 1 &&
        length(unique(df_tmp$TRG_score)) > 1) {
      
      rc = Hmisc::rcorr(as.matrix(df_tmp), type = "spearman")
      tibble(
        n = nrow(df_tmp),
        spearman_r = rc$r["abundance", "TRG_score"],
        p_value = rc$P["abundance", "TRG_score"]
      )
    } else {
      tibble(n = nrow(df_tmp), spearman_r = NA_real_, p_value = NA_real_)
    }
  }) %>%
  ungroup(); head(cor_t) 

cor_t %>% 
  filter(p_value < 0.05) %>% 
  arrange(p_value)

# TRG_1
p4_TRG_score_strain = tb_input %>% 
  filter(Strain %in% sig_strain) %>% 
  mutate(Strain = factor(Strain, levels = sig_strain)) %>% 
  ggplot(aes(TRG_score, log10(abundance + 1.5e-05))) +
  geom_jitter(aes(fill = TRG_score), shape = 21, size = 2) +
  geom_smooth(method = "lm", se = 0.95, color = "#898989", alpha = 0.3) +
  scale_fill_gradient(low = "#50C878",
                      high = "#ff746c") +
  facet_wrap(~Strain, nrow = 4) +
  stat_cor(method = "spearman") +
  theme_pubr() +
  theme(legend.position = "none",
        strip.background = element_rect(color = NA),
        axis.title = element_text(size = rel(1.05)),
        aspect.ratio = 0.9) +
  labs(x = "Tumor Regression Grade",
       y = "Log10 (Abundance + pseudo)") ; p4_TRG_score_strain

# ggsave("figure/04-3-6_TRG_score-strain.svg",
#        plot = p4_TRG_score_strain, width = 12, height = 9)


## Eubacterium rectale visualization
library(ggbeeswarm)

E.rectale_tile = tb_input %>%
  filter(Strain == "Eubacterium_rectale|SGB4933") %>% 
  mutate(
    TRG_score = factor(TRG_score, levels = 0:3),
    Pre_Op_Tstage = factor(Pre_Op_Tstage, levels = c("2", "3", "4"))
  )

bubble_df = E.rectale_tile %>%
  group_by(TRG_score, Pre_Op_Tstage) %>%
  mutate(
    id_within_tile = row_number(),
    row = ((id_within_tile - 1) %/% 3),
    col = ((id_within_tile - 1) %% 3)
  ) %>%
  ungroup() %>%
  mutate(
    x_adj = as.numeric(TRG_score) - 0.3 + col * 0.15,
    y_adj = as.numeric(Pre_Op_Tstage) - 0.3 + row * 0.2
  )

# 빈 타일 격자 생성
tile_grid = expand.grid(
  TRG_score = factor(0:3),
  Pre_Op_Tstage = factor(c("2", "3", "4"), levels = c("2", "3", "4"))
)

p4_TRG_score_E.rectale = ggplot() +
  geom_tile(data = tile_grid, 
            aes(x = as.numeric(TRG_score), 
                y = as.numeric(Pre_Op_Tstage)),
            fill = "white", color = "grey80", width = 0.95, height = 0.95) +
  
  geom_point(data = bubble_df,
             aes(x = x_adj, y = y_adj, size = abundance, fill = abundance),
             shape = 21, alpha = 0.9) +
  
  scale_size_continuous(range = c(1.5, 7), name = "E. rectale") +
  scale_fill_viridis_c(name = "Abundance (%)") +
  
  scale_x_continuous(breaks = 1:4, labels = levels(E.rectale_tile$TRG_score),
                     name = "TRG Score") +
  scale_y_continuous(breaks = 1:3, labels = levels(E.rectale_tile$Pre_Op_Tstage),
                     name = "Preoperative T stage") +
  
  coord_fixed() +
  theme_minimal(base_size = 14) + 
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = rel(1.05), vjust = 1),
    axis.text.y = element_text(size = rel(1.05), hjust = 1),
    legend.position = "right"
  ) ; p4_TRG_score_E.rectale

# ggsave("figure/04-3-7_TRG_score-strain-E.rectale.svg",
#        plot = p4_TRG_score_E.rectale, width = 7.5, height = 4.5)



### In ongoing samples

# Prevalence filtering (>=10%)
to = t %>% select(c("Strain", mo$SampleID))

to_0.1 = to %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Strain) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.1) %>% 
  pull(Strain) %>% 
  sort() ; to_0.1


# Make an input file
to_input = to %>% 
  filter(Strain %in% to_0.1) %>% 
  pivot_longer(cols = -Strain,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  merge(mo %>% 
          select(SampleID,
                 TRG_score, TRG_1, TRG_2, TRG_3,
                 Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                 Pre_Op_Nstage, Pre_Op_Nstage_bin,
                 BMI, BMI_bin, CEA, CEA_bin),
        by = "SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR")))

min_abund = to_input %>% 
  filter(abundance > 0) %>% 
  pull(abundance) %>% min() # 4e-05


# TRG_1: CR vs. nonCR
# Wilcoxon test (CR vs. nonCR)
wilcox_to = to_input %>% 
  group_by(Strain) %>% 
  summarise(p_value = wilcox.test(abundance ~ TRG_1, exact = F)$p.value,
            .groups = "drop") %>% 
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>% 
  arrange(p_value) ; wilcox_to

sig_strain_o = wilcox_to %>% 
  filter(p_value < 0.05) %>% 
  arrange(p_value) %>% 
  pull(Strain)
# 19 strains, P < 0.1
# 9 strains, P < 0.05

sig_strain_o = to_input %>% 
  filter(Strain %in% sig_strain_o) %>% 
  group_by(Strain) %>% 
  summarise(avg_abund = mean(abundance)) %>% 
  arrange(-avg_abund) %>% 
  pull(Strain) ; sig_strain_o

p5_TRG_1_strain_o = to_input %>% 
  filter(Strain %in% sig_strain_o) %>% 
  mutate(Strain = factor(Strain, levels = sig_strain_o)) %>% 
  ggplot(aes(TRG_1, log10(abundance + min_abund/2))) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  facet_wrap(~Strain, ncol = 2) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0) +
  theme_classic() +
  coord_flip(ylim = c(-5, 1.5)) +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(0.95)),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = rel(1.05)),
        axis.text.y = element_text(size = rel(1.3)),
        strip.background = element_rect(color = NA),
        aspect.ratio = 0.2) +
  labs(x = "Tumor Regression Grade",
       y = "Log10 (Abundance + 2e-5)") ; p5_TRG_1_strain_o

# ggsave("figure/05-3-1_TRG_1-strain_ongoing.svg",
#        plot = p5_TRG_1_strain_o, width = 7, height = 6.5)



### In paired samples

tp = t %>% select(c("Strain", mp$SampleID))

# TRG_1
t %>% 
  filter(Strain %in% c(wilcox_to %>% filter(p_value < 0.1) %>% arrange(p_value) %>% pull(Strain), 
                       wilcox_t %>% filter(p_value < 0.1) %>% arrange(p_value) %>% pull(Strain))) %>% 
  pivot_longer(cols = -Strain,
               names_to = "SampleID",
               values_to = "Abundance") %>% 
  merge(., mp, by.x = "SampleID") %>% 
  group_by(Strain, TRG_1, TNT) %>% 
  summarise(Abundance = mean(Abundance, na.rm = T), .groups = "drop") %>% 
  pivot_wider(names_from = c(TNT, TRG_1), values_from = Abundance) %>% 
  mutate(
    # Before에서 CR / noncR
    Before_FC = log10((`Before_CR`+1e-05) / (`Before_nonCR`+1e-05)),
    
    # Ongoing에서 CR/ nonCR
    Ongoing_FC = log10((`Ongoing_CR`+1e-05) / (`Ongoing_nonCR`+1e-05)),
    
    # CR에서 Ongoing / Before
    CR_FC = log10((`Ongoing_CR`+1e-05) / (`Before_CR`+1e-05)),
    
    # nonCR에서 Ongoing / Before
    nonCR_FC = log10((`Ongoing_nonCR`+1e-05) / (`Before_nonCR`+1e-05))
  ) %>% 
  select(Strain, Before_CR, Before_nonCR, Before_FC, Ongoing_CR, Ongoing_nonCR, Ongoing_FC, CR_FC, nonCR_FC)

# Strains selected based on fold-change
FC_strains = c("Fournierella_massiliensis|SGB72786",
               "Veillonella_dispar|SGB6952",
               "Ruminococcus_gnavus|SGB4571",
               "GGB9597_SGB15022|SGB15022",
               "GGB9708_SGB15234|SGB15234",
               "GGB9635_SGB15106|SGB15106",
               "Roseburia_inulinivorans|SGB4940",
               "Lacrimispora_amygdalina|SGB4716",
               "GGB35068_SGB47850|SGB47850",
               "Collinsella_aerofaciens|SGB14483",
               "Allisonella_histaminiformans|SGB5843",
               "Agathobaculum_butyriciproducens|SGB14993",
               "Gemella_sanguinis|SGB7298_group",
               "Clostridiales_bacterium_KLE1615|SGB5090",
               "Coprococcus_comes|SGB4577",
               "Propionibacterium_acidifaciens|SGB15903",
               "Rothia_mucilaginosa|SGB16986",
               "Faecalibacterium_prausnitzii|SGB15316",
               "Gemella_morbillorum|SGB7295",
               "Candidatus_Nanosynsacchari_sp_TM7_ANC_38_39_G1_1|SGB19882",
               "TM7_phylum_sp_oral_taxon_348|SGB19880",
               "Intestinimonas_butyriciproducens|SGB15126")


# Make an input file
# Prevalence >= 0.1, P < 0.05
# FC 차이 나는 strain
# CRC or colitis inducing bacteria 
# Before sample & CR, nonCR의 abundance가 0 이상인 sample이 3개 이상
tp_input = tp %>%
  filter(Strain %in% c(sig_strain, FC_strains, bad_t0)) %>% 
  pivot_longer(cols = -Strain,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  merge(mp %>% select(SampleID, TRG_1, TNT), by = "SampleID") %>% 
  group_by(Strain) %>% 
  summarise(abundance_nonzero_CR = sum(abundance > 0 & TRG_1 == "CR" & TNT == "Before"),
            abundance_nonzero_nonCR = sum(abundance > 0 & TRG_1 == "nonCR" & TNT == "Before"),
            .groups = "drop") %>% 
  filter(abundance_nonzero_CR >= 3 & abundance_nonzero_nonCR >= 3) %>% 
  pull(Strain)

# tp_input = tp %>%
#   filter(Strain %in% tp_input) %>% 
#   pivot_longer(cols = -Strain,
#                names_to = "SampleID",
#                values_to = "abundance") %>% 
#   merge(mp %>% 
#           select(SampleID, SNU_ID, TNT,
#                  TRG_score, TRG_1, TRG_2, TRG_3,
#                  Pre_Op_Tstage, Pre_Op_Tstage_bin, 
#                  Pre_Op_Nstage, Pre_Op_Nstage_bin),
#         by = "SampleID") %>% 
#   mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR"))) %>% 
#   group_by(Strain, TRG_1, TNT) %>%
#   summarise(total_abundance = sum(abundance + 1e-5), .groups = "drop") %>%  # pseudo-count 추가
#   pivot_wider(names_from = TNT, values_from = total_abundance) %>%
#   mutate(
#     log10FC = log10((Ongoing + 1e-05) / (Before + 1e-05))
#   ) %>%
#   select(Strain, TRG_1, log10FC) %>%
#   pivot_wider(names_from = TRG_1, values_from = log10FC, names_prefix = "log10FC_") %>%
#   mutate(Strain = fct_reorder(Strain, log10FC_CR)) %>%
#   pivot_longer(cols = starts_with("log10FC_"), names_to = "Group", values_to = "log10FC") %>%
#   mutate(Group = recode(Group, log10FC_CR = "CR", log10FC_nonCR = "nonCR"))

# p5_Fold_change = tp_input %>% 
#   mutate(Group = factor(Group, levels = c("nonCR", "CR"))) %>% 
#   ggplot(aes(log10FC, Strain, fill = Group)) +
#   geom_bar(stat = "identity", color = "black",
#            position = position_dodge(width = 1)) +
#   geom_errorbar(aes(xmin = log10FC - ci95, xmax = log10FC + ci95),
#                 width = 0.2, position = position_dodge(width = 1)) +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   scale_fill_manual(values = c("CR" = "#F7D9BC",
#                                "nonCR" = "#80461B"),
#                     breaks = c("CR", "nonCR"),
#                     labels = c("CR", "nonCR")) +
#   labs(title = "Log10 Fold Change",
#        x = "Log10 Fold Change (Ongoing / Before)",
#        y = "Strain",
#        fill = "TRG_1") +
#   theme_classic() ; p5_Fold_change

# ggsave("figure/06-3-1_FC-strain_paired (total).svg",
#        plot = p5_Fold_change, width = 10, height = 9)

tp_input = tp %>%
  filter(Strain %in% tp_input) %>% 
  pivot_longer(cols = -Strain,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  left_join(mp %>% select(SNU_ID, SampleID, TRG_1, TNT), by = "SampleID") %>%
  select(-SampleID) %>% 
  pivot_wider(names_from = TNT, values_from = abundance) %>%
  mutate(log10FC = log10((Ongoing + 1e-05) / (Before + 1e-05))) %>%
  group_by(Strain, TRG_1)

tp_input_summary = tp_input %>% 
  summarise(
    mean_log10FC = mean(log10FC, na.rm = T),
    sd_log10FC = sd(log10FC, na.rm = T),
    n = n(),
    se_log10FC = sd_log10FC / sqrt(n),
    ci95 = 1.96 * se_log10FC,
    .groups = "drop") %>% 
  arrange(desc(TRG_1 == "CR"), desc(mean_log10FC))

tp_input_wilcox = tp_input %>% 
  group_by(Strain) %>% 
  summarise(
    wilcox_test = list(wilcox.test(log10FC ~ TRG_1, data = cur_data(), exact = FALSE)),
    p_value = wilcox_test[[1]]$p.value,
    statistic = wilcox_test[[1]]$statistic,
    .groups = "drop"
  ) %>%
  select(Strain, p_value, statistic) %>% 
  arrange(p_value)

tp_input_summary = tp_input_summary %>% 
  left_join(tp_input_wilcox %>% select(Strain, p_value), by = "Strain") ; rm(tp_input_wilcox)

p5_Fold_change = tp_input_summary %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("nonCR", "CR")),
         Strain = factor(Strain, levels = rev(unique(tp_input_summary$Strain)))) %>% 
  ggplot(aes(mean_log10FC, Strain, fill = TRG_1)) +
  geom_errorbar(aes(xmin = mean_log10FC - ci95, xmax = mean_log10FC + ci95),
                width = 0.2, position = position_dodge(width = 1), alpha = 0.5) +
  geom_point(aes(fill = TRG_1), shape = 21, size = 4,
             stat = "identity", position = position_dodge(width = 1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B"),
                    breaks = c("CR", "nonCR"),
                    labels = c("CR", "nonCR")) +
  labs(title = "Log10 Fold Change",
       x = "Log10 Fold Change (Ongoing / Before)",
       y = "Strain",
       fill = "TRG_1") +
  theme_classic() ; p5_Fold_change

# ggsave("figure/06-3-1_FC-strain_paired (mean).svg",
#        plot = p5_Fold_change, width = 10, height = 9)


FC_8strains = c("Lacrimispora_amygdalina|SGB4716",
                "Coprococcus_comes|SGB4577",
                "Faecalibacterium_prausnitzii|SGB15316",
                "Agathobaculum_butyriciproducens|SGB14993",
                "Clostridium_SGB4750|SGB4750",
                "Eubacterium_rectale|SGB4933",
                "Ruminococcus_sp_AF13_28|SGB4834",
                "Lancefieldella_parvula|SGB966")

p5_Fold_change_8 = tp_input_summary %>% 
  filter(Strain %in% FC_8strains) %>% 
  mutate(Strain = factor(Strain, levels = rev(FC_8strains))) %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("nonCR", "CR"))) %>% 
  ggplot(aes(mean_log10FC, Strain, fill = TRG_1)) +
  geom_bar(stat = "identity", color = "black",
           position = position_dodge(width = 1)) +
  geom_errorbar(aes(xmin = mean_log10FC - ci95, xmax = mean_log10FC + ci95),
                width = 0.2, position = position_dodge(width = 1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B"),
                    breaks = c("CR", "nonCR"),
                    labels = c("CR", "nonCR")) +
  labs(title = "Log10 Fold Change",
       x = "Log10 Fold Change (Ongoing / Before)",
       y = "Strain",
       fill = "TRG_1") +
  theme_classic() ; p5_Fold_change_8

# ggsave("figure/06-3-1_FC-8strain_paired.svg",
#        plot = p5_Fold_change_8, width = 10, height = 4)


# nonCR-enrichment 10 strains
sig_bad_top10

# sig_bad_top10_input = tp %>%
#   filter(Strain %in% sig_bad_top10) %>%
#   pivot_longer(cols = -Strain,
#                names_to = "SampleID",
#                values_to = "abundance") %>%
#   merge(mp %>%
#           select(SampleID, SNU_ID, TNT,
#                  TRG_score, TRG_1, TRG_2, TRG_3,
#                  Pre_Op_Tstage, Pre_Op_Tstage_bin,
#                  Pre_Op_Nstage, Pre_Op_Nstage_bin),
#         by = "SampleID") %>%
#   mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR"))) %>%
#   group_by(Strain, TRG_1, TNT) %>%
#   summarise(total_abundance = sum(abundance + 1e-5), .groups = "drop") %>%  # pseudo-count 추가
#   pivot_wider(names_from = TNT, values_from = total_abundance) %>%
#   mutate(
#     log10FC = log10((Ongoing + 1e-05) / (Before + 1e-05))
#   ) %>%
#   select(Strain, TRG_1, log10FC) %>%
#   pivot_wider(names_from = TRG_1, values_from = log10FC, names_prefix = "log10FC_") %>%
#   mutate(Strain = fct_reorder(Strain, log10FC_CR)) %>%
#   pivot_longer(cols = starts_with("log10FC_"), names_to = "Group", values_to = "log10FC") %>%
#   mutate(Group = recode(Group, log10FC_CR = "CR", log10FC_nonCR = "nonCR"))

# p5_Fold_change_nonCR = sig_bad_top10_input %>%
#   mutate(Group = factor(Group, levels = c("nonCR", "CR"))) %>%
#   ggplot(aes(log10FC, Strain, fill = Group)) +
#   geom_bar(stat = "identity", color = "black",
#            position = position_dodge(width = 1)) +
#   geom_errorbar(aes(xmin = log10FC - ci95, xmax = log10FC + ci95),
#                 width = 0.2, position = position_dodge(width = 1)) +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   scale_fill_manual(values = c("CR" = "#F7D9BC",
#                                "nonCR" = "#80461B"),
#                     breaks = c("CR", "nonCR"),
#                     labels = c("CR", "nonCR")) +
#   labs(title = "Log10 Fold Change",
#        x = "Log10 Fold Change (Ongoing / Before)",
#        y = "Strain",
#        fill = "TRG_1") +
#   theme_classic() ; p5_Fold_change_nonCR

# ggsave("figure/06-3-1_FC-strain_paired_nonCR (total).svg",
#        plot = p5_Fold_change_nonCR, width = 10, height = 9)

sig_bad_top10_input = tp %>%
  filter(Strain %in% sig_bad_top10) %>% 
  pivot_longer(cols = -Strain,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  left_join(mp %>% select(SNU_ID, SampleID, TRG_1, TNT), by = "SampleID") %>%
  select(-SampleID) %>% 
  pivot_wider(names_from = TNT, values_from = abundance) %>%
  mutate(log10FC = log10((Ongoing + 1e-05) / (Before + 1e-05))) %>%
  group_by(Strain, TRG_1)

sig_bad_top10_summary = sig_bad_top10_input %>% 
  summarise(
    mean_log10FC = mean(log10FC, na.rm = T),
    sd_log10FC = sd(log10FC, na.rm = T),
    n = n(),
    se_log10FC = sd_log10FC / sqrt(n),
    ci95 = 1.96 * se_log10FC,
    .groups = "drop") %>% 
  arrange(desc(TRG_1 == "CR"), desc(mean_log10FC))

sig_bad_top10_wilcox = sig_bad_top10_input %>% 
  group_by(Strain) %>% 
  summarise(
    wilcox_test = list(wilcox.test(log10FC ~ TRG_1, data = cur_data(), exact = FALSE)),
    p_value = wilcox_test[[1]]$p.value,
    statistic = wilcox_test[[1]]$statistic,
    .groups = "drop"
  ) %>%
  select(Strain, p_value, statistic) %>% 
  arrange(p_value)

sig_bad_top10_summary = sig_bad_top10_summary %>% 
  left_join(sig_bad_top10_wilcox %>% select(Strain, p_value), by = "Strain") ; rm(sig_bad_top10_wilcox)


p5_Fold_change_nonCR = sig_bad_top10_summary %>% 
  filter(mean_log10FC != 0) %>% 
  mutate(TRG_1 = factor(TRG_1, level = c("nonCR", "CR")),
         Strain = factor(Strain, levels = rev(unique(sig_bad_top10_summary$Strain)))) %>% 
  ggplot(aes(mean_log10FC, Strain, fill = TRG_1)) +
  geom_errorbar(aes(xmin = mean_log10FC - ci95, xmax = mean_log10FC + ci95),
                width = 0.2, position = position_dodge(width = 1)) +
  geom_point(aes(fill = TRG_1), shape = 21, size = 4,
             stat = "identity", position = position_dodge(width = 1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B"),
                    breaks = c("CR", "nonCR"),
                    labels = c("CR", "nonCR")) +
  labs(title = "Log10 Fold Change",
       x = "Log10 Fold Change (Ongoing / Before)",
       y = "Strain",
       fill = "TRG_1") +
  theme_classic() ; p5_Fold_change_nonCR

# ggsave("figure/06-3-1_FC-strain_paired_nonCR (mean).svg",
#        plot = p5_Fold_change_nonCR, width = 10, height = 9)


# Shared strains
# Prevalence ≥ 10%, Higher abundance ≥ Mean higher abundance * 0.1
tp %>%
  rowwise() %>%
  mutate(prevalence = sum(c_across(-Strain) > 0)) %>%
  ungroup() %>%
  filter(prevalence >= nrow(m) * 0.1) %>%
  pivot_longer(cols = -Strain,
               names_to = "SampleID",
               values_to = "abundance") %>%
  merge(mp %>% select(SampleID, TNT, TRG_1), by = "SampleID") %>%
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR"))) %>%
  group_by(Strain, TNT, TRG_1) %>%
  summarise(mean_abundance = mean(abundance, na.rm = T),
            prevalence = sum(abundance > 0, na.rm = T)) %>%
  pivot_wider(names_from = c("TNT", "TRG_1"),
              values_from = c(mean_abundance, prevalence),
              names_glue = "{TNT}_{TRG_1}_{.value}") %>%
  write_csv(., "strain_paired (25.11.17).csv")

# 98 strains
paired_strain = readxl::read_xlsx("strain_paired (25.11.17).xlsx", sheet = "Feature selection") %>% 
  select(Strain) %>% 
  unique() ; paired_strain


# Abundant, Prevalent (Abund. ≥ 0.1, Prev. ≥ 8)
abund_prev_strain = c("Agathobaculum_butyriciproducens|SGB14993",
                      "Faecalibacterium_prausnitzii|SGB15316",
                      "Roseburia_inulinivorans|SGB4940",
                      "Clostridiales_bacterium_KLE1615|SGB5090",
                      "Eubacterium_rectale|SGB4933",
                      "Bifidobacterium_longum|SGB17248",
                      "Collinsella_aerofaciens|SGB14483",
                      "Coprococcus_comes|SGB4577")

tp_abund_prev = tp %>% 
  filter(Strain %in% abund_prev_strain) %>% 
  pivot_longer(cols = -Strain,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  merge(mp %>% 
          select(SampleID, SNU_ID, TNT,
                 TRG_score, TRG_1, TRG_2, TRG_3,
                 Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                 Pre_Op_Nstage, Pre_Op_Nstage_bin,
                 BMI, BMI_bin, CEA, CEA_bin),
        by = "SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR")))

# Agathobaculum_butyriciproducens|SGB14993
p5_Agathobaculum_butyriciproducens = tp_abund_prev %>% 
  filter(Strain == "Agathobaculum_butyriciproducens|SGB14993") %>% 
  mutate(TRG = factor(TRG_1, levels = c("CR", "nonCR")),
         TNT = factor(TNT, levels = c("Before", "Ongoing")),
         Group = paste(TNT, TRG_1, sep = "_"),
         Group = factor(Group, levels = c("Before_CR", "Ongoing_CR", 
                                          "Before_nonCR", "Ongoing_nonCR"))) %>% 
  ggplot(aes(Group, abundance)) +
  geom_point(aes(color = TRG_1)) +
  geom_line(aes(group = SNU_ID), color = "gray60", alpha = 0.5,
            arrow = arrow(type = "open", length = unit(0.05, "inches"))) +
  geom_boxplot(aes(fill = TRG_1),
               alpha = 0.3, outlier.alpha = 0,
               width = 0.6) +
  stat_compare_means(comparisons = list(c("Before_CR", "Ongoing_CR"),
                                        c("Before_nonCR", "Ongoing_nonCR")),
                     tip.length = 0.02,
                     method = "wilcox", paired = T) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Sample", y = "Abundance", 
       title = "Agathobaculum butyriciproducens (SGB14993)") ; p5_Agathobaculum_butyriciproducens

# Faecalibacterium_prausnitzii|SGB15316
p5_Faecalibacterium_prausnitzii = tp_abund_prev %>% 
  filter(Strain == "Faecalibacterium_prausnitzii|SGB15316") %>% 
  mutate(TRG = factor(TRG_1, levels = c("CR", "nonCR")),
         TNT = factor(TNT, levels = c("Before", "Ongoing")),
         Group = paste(TNT, TRG_1, sep = "_"),
         Group = factor(Group, levels = c("Before_CR", "Ongoing_CR", 
                                          "Before_nonCR", "Ongoing_nonCR"))) %>% 
  ggplot(aes(Group, abundance)) +
  geom_point(aes(color = TRG_1)) +
  geom_line(aes(group = SNU_ID), color = "gray60", alpha = 0.5,
            arrow = arrow(type = "open", length = unit(0.05, "inches"))) +
  geom_boxplot(aes(fill = TRG_1),
               alpha = 0.3, outlier.alpha = 0,
               width = 0.6) +
  stat_compare_means(comparisons = list(c("Before_CR", "Ongoing_CR"),
                                        c("Before_nonCR", "Ongoing_nonCR")),
                     tip.length = 0.02,
                     method = "wilcox", paired = T) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Sample", y = "Abundance", 
       title = "Faecalibacterium prausnitzii (SGB15316)") ; p5_Faecalibacterium_prausnitzii

# Roseburia_inulinivorans|SGB4940
p5_Roseburia_inulinivorans = tp_abund_prev %>% 
  filter(Strain == "Roseburia_inulinivorans|SGB4940") %>% 
  mutate(TRG = factor(TRG_1, levels = c("CR", "nonCR")),
         TNT = factor(TNT, levels = c("Before", "Ongoing")),
         Group = paste(TNT, TRG_1, sep = "_"),
         Group = factor(Group, levels = c("Before_CR", "Ongoing_CR", 
                                          "Before_nonCR", "Ongoing_nonCR"))) %>% 
  ggplot(aes(Group, abundance)) +
  geom_point(aes(color = TRG_1)) +
  geom_line(aes(group = SNU_ID), color = "gray60", alpha = 0.5,
            arrow = arrow(type = "open", length = unit(0.05, "inches"))) +
  geom_boxplot(aes(fill = TRG_1),
               alpha = 0.3, outlier.alpha = 0,
               width = 0.6) +
  stat_compare_means(comparisons = list(c("Before_CR", "Ongoing_CR"),
                                        c("Before_nonCR", "Ongoing_nonCR")),
                     tip.length = 0.02,
                     method = "wilcox", paired = T) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Sample", y = "Abundance", 
       title = "Roseburia inulinivorans (SGB4940)") ; p5_Roseburia_inulinivorans

# Clostridiales_bacterium_KLE1615|SGB5090
p5_Clostridiales_bacterium_KLE1615 = tp_abund_prev %>% 
  filter(Strain == "Clostridiales_bacterium_KLE1615|SGB5090") %>% 
  mutate(TRG = factor(TRG_1, levels = c("CR", "nonCR")),
         TNT = factor(TNT, levels = c("Before", "Ongoing")),
         Group = paste(TNT, TRG_1, sep = "_"),
         Group = factor(Group, levels = c("Before_CR", "Ongoing_CR", 
                                          "Before_nonCR", "Ongoing_nonCR"))) %>% 
  ggplot(aes(Group, abundance)) +
  geom_point(aes(color = TRG_1)) +
  geom_line(aes(group = SNU_ID), color = "gray60", alpha = 0.5,
            arrow = arrow(type = "open", length = unit(0.05, "inches"))) +
  geom_boxplot(aes(fill = TRG_1),
               alpha = 0.3, outlier.alpha = 0,
               width = 0.6) +
  stat_compare_means(comparisons = list(c("Before_CR", "Ongoing_CR"),
                                        c("Before_nonCR", "Ongoing_nonCR")),
                     tip.length = 0.02,
                     method = "wilcox", paired = T) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Sample", y = "Abundance", 
       title = "Clostridiales bacterium KLE1615 (SGB5090)") ; p5_Clostridiales_bacterium_KLE1615

# Eubacterium_rectale|SGB4933
p5_Eubacterium_rectale = tp_abund_prev %>% 
  filter(Strain == "Eubacterium_rectale|SGB4933") %>% 
  mutate(TRG = factor(TRG_1, levels = c("CR", "nonCR")),
         TNT = factor(TNT, levels = c("Before", "Ongoing")),
         Group = paste(TNT, TRG_1, sep = "_"),
         Group = factor(Group, levels = c("Before_CR", "Ongoing_CR", 
                                          "Before_nonCR", "Ongoing_nonCR"))) %>% 
  ggplot(aes(Group, abundance)) +
  geom_point(aes(color = TRG_1)) +
  geom_line(aes(group = SNU_ID), color = "gray60", alpha = 0.5,
            arrow = arrow(type = "open", length = unit(0.05, "inches"))) +
  geom_boxplot(aes(fill = TRG_1),
               alpha = 0.3, outlier.alpha = 0,
               width = 0.6) +
  stat_compare_means(comparisons = list(c("Before_CR", "Ongoing_CR"),
                                        c("Before_nonCR", "Ongoing_nonCR")),
                     tip.length = 0.02,
                     method = "wilcox", paired = T) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Sample", y = "Abundance", 
       title = "Eubacterium rectale (SGB4933)") ; p5_Eubacterium_rectale

# Bifidobacterium_longum|SGB17248
p5_Bifidobacterium_longum = tp_abund_prev %>% 
  filter(Strain == "Bifidobacterium_longum|SGB17248") %>% 
  mutate(TRG = factor(TRG_1, levels = c("CR", "nonCR")),
         TNT = factor(TNT, levels = c("Before", "Ongoing")),
         Group = paste(TNT, TRG_1, sep = "_"),
         Group = factor(Group, levels = c("Before_CR", "Ongoing_CR", 
                                          "Before_nonCR", "Ongoing_nonCR"))) %>% 
  ggplot(aes(Group, abundance)) +
  geom_point(aes(color = TRG_1)) +
  geom_line(aes(group = SNU_ID), color = "gray60", alpha = 0.5,
            arrow = arrow(type = "open", length = unit(0.05, "inches"))) +
  geom_boxplot(aes(fill = TRG_1),
               alpha = 0.3, outlier.alpha = 0,
               width = 0.6) +
  stat_compare_means(comparisons = list(c("Before_CR", "Ongoing_CR"),
                                        c("Before_nonCR", "Ongoing_nonCR")),
                     tip.length = 0.02,
                     method = "wilcox", paired = T) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Sample", y = "Abundance", 
       title = "Bifidobacterium longum (SGB17248)") ; p5_Bifidobacterium_longum

# Collinsella_aerofaciens|SGB14483
p5_Collinsella_aerofaciens = tp_abund_prev %>% 
  filter(Strain == "Collinsella_aerofaciens|SGB14483") %>% 
  mutate(TRG = factor(TRG_1, levels = c("CR", "nonCR")),
         TNT = factor(TNT, levels = c("Before", "Ongoing")),
         Group = paste(TNT, TRG_1, sep = "_"),
         Group = factor(Group, levels = c("Before_CR", "Ongoing_CR", 
                                          "Before_nonCR", "Ongoing_nonCR"))) %>% 
  ggplot(aes(Group, abundance)) +
  geom_point(aes(color = TRG_1)) +
  geom_line(aes(group = SNU_ID), color = "gray60", alpha = 0.5,
            arrow = arrow(type = "open", length = unit(0.05, "inches"))) +
  geom_boxplot(aes(fill = TRG_1),
               alpha = 0.3, outlier.alpha = 0,
               width = 0.6) +
  stat_compare_means(comparisons = list(c("Before_CR", "Ongoing_CR"),
                                        c("Before_nonCR", "Ongoing_nonCR")),
                     tip.length = 0.02,
                     method = "wilcox", paired = T) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Sample", y = "Abundance", 
       title = "Collinsella aerofaciens (SGB14483)") ; p5_Collinsella_aerofaciens

# Coprococcus_comes|SGB4577
p5_Coprococcus_comes = tp_abund_prev %>% 
  filter(Strain == "Coprococcus_comes|SGB4577") %>% 
  mutate(TRG = factor(TRG_1, levels = c("CR", "nonCR")),
         TNT = factor(TNT, levels = c("Before", "Ongoing")),
         Group = paste(TNT, TRG_1, sep = "_"),
         Group = factor(Group, levels = c("Before_CR", "Ongoing_CR", 
                                          "Before_nonCR", "Ongoing_nonCR"))) %>% 
  ggplot(aes(Group, abundance)) +
  geom_point(aes(color = TRG_1)) +
  geom_line(aes(group = SNU_ID), color = "gray60", alpha = 0.5,
            arrow = arrow(type = "open", length = unit(0.05, "inches"))) +
  geom_boxplot(aes(fill = TRG_1),
               alpha = 0.3, outlier.alpha = 0,
               width = 0.6) +
  stat_compare_means(comparisons = list(c("Before_CR", "Ongoing_CR"),
                                        c("Before_nonCR", "Ongoing_nonCR")),
                     tip.length = 0.02,
                     method = "wilcox", paired = T) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Sample", y = "Abundance", 
       title = "Coprococcus comes (SGB4577)") ; p5_Coprococcus_comes

p5_paired_strains = grid.arrange(p5_Agathobaculum_butyriciproducens,
                                 p5_Faecalibacterium_prausnitzii,
                                 p5_Roseburia_inulinivorans,
                                 p5_Clostridiales_bacterium_KLE1615,
                                 p5_Eubacterium_rectale,
                                 p5_Bifidobacterium_longum,
                                 p5_Collinsella_aerofaciens,
                                 p5_Coprococcus_comes,
                                 ncol = 4, nrow = 2) ; p5_paired_strains

# ggsave("figure/06-3-1_tendency of strains.svg",
#        plot = p5_paired_strains, width = 16, height = 10)


### Bad bacteria
# 이미 대장암에서 나쁘다고 알려진 bacteria
bad_bac = c("Peptostreptococcus_stomatis|SGB748",
            "Gemella_morbillorum|SGB7295",
            "Parvimonas_micra|SGB6653",
            "Hungatella_hathewayi|SGB4742",
            "Fusobacterium_nucleatum|SGB6011",
            "Solobacterium_moorei|SGB6826")

bad_bac_input = tp %>% 
  filter(Strain %in% bad_bac) %>% 
  pivot_longer(cols = -Strain,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  left_join(mp %>% select(SNU_ID, SampleID, TRG_1, TNT), by = "SampleID") %>% 
  select(-SampleID) %>% 
  pivot_wider(names_from = TNT, values_from = abundance) %>% 
  mutate(log10FC = log10((Ongoing + 1e-05) / (Before + 1e-05))) %>% 
  group_by(Strain)

bad_bac_summary = bad_bac_input %>% 
  summarise(
    mean_log10FC = mean(log10FC, na.rm = T),
    sd_log10FC = sd(log10FC, na.rm = T),
    n = n(),
    se_log10FC = sd_log10FC / sqrt(n),
    ci95 = 1.96 * se_log10FC,
    .groups = "drop"
  ) %>% 
  arrange(desc(mean_log10FC))

bad_bac_wilcox = bad_bac_input %>% 
  group_by(Strain) %>% 
  summarise(
    wilcox_test = list(wilcox.test(log10FC ~ TRG_1, data = cur_data(), exact = F)),
    p_value = wilcox_test[[1]]$p.value,
    statistic = wilcox_test[[1]]$statistic,
    .groups = "drop"
  ) %>% 
  select(Strain, p_value, statistic) %>% 
  arrange(p_value)

bad_bac_summary = bad_bac_summary %>% 
  left_join(bad_bac_wilcox %>% select(Strain, p_value), by = "Strain") ; rm(bad_bac_wilcox)


p5_Fold_change_bad_bac = bad_bac_summary %>% 
  mutate(Strain = factor(Strain, levels = rev(unique(bad_bac_summary$Strain)))) %>% 
  ggplot(aes(mean_log10FC, Strain)) +
  geom_errorbar(aes(xmin = mean_log10FC - ci95, xmax = mean_log10FC + ci95),
                width = 0.2, position = position_dodge(width = 1)) +
  geom_point(shape = 21, size = 4, fill = "steelblue",
             stat = "identity", position = position_dodge(width = 1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(aes(x = max(mean_log10FC + ci95 + 0.1), label = sprintf("%.3f", p_value)),
            position = position_dodge(width = 1), hjust = 0, size = 3) +
  labs(title = "Log10 Fold Change",
       x = "Log10 Fold Change (Ongoing / Before)", y = "Strain") +
  theme_classic() ; p5_Fold_change_bad_bac

# ggsave("figure/06-3-1_FC-strain_paired_bad_bac (mean).svg",
#        plot = p5_Fold_change_bad_bac, width = 10, height = 6)



###############  2. Species  ###############

# Prevalence filtering
sb_0.2 = sb %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Species) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.2) %>% 
  pull(Species) %>% 
  sort() ; sb_0.2


# Make an input file for differential enrichment testing 
sb_input = sb %>% 
  filter(Species %in% sb_0.2) %>% 
  pivot_longer(cols = -Species,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  merge(mb2 %>% 
          select(SampleID,
                 TRG_score, TRG_1, TRG_2, TRG_3,
                 Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                 Pre_Op_Nstage, Pre_Op_Nstage_bin,
                 BMI, BMI_bin, CEA, CEA_bin),
        by = "SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR")))

min_abund = sb_input %>% 
  filter(abundance >0) %>% 
  pull(abundance) %>% min() # 2e-05


# TRG_1: CR vs. nonCR
# Wilcoxon test (CR vs. nonCR)
wilcox_s = sb_input %>% 
  group_by(Species) %>% 
  summarise(p_value = wilcox.test(abundance ~ TRG_1, exact = F)$p.value,
            .groups = "drop") %>% 
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>% 
  arrange(p_value) ; wilcox_s

sig_species = wilcox_s %>% 
  filter(p_value < 0.1) %>% 
  arrange(p_value) %>% 
  pull(Species) ; sig_species
# 25 species, P < 0.1
# 12 species, P < 0.05

sig_species = sb_input %>% 
  filter(Species %in% sig_species) %>% 
  group_by(Species) %>% 
  summarise(avg_abund = mean(abundance)) %>% 
  arrange(-avg_abund) %>% 
  pull(Species) ; sig_species

p4_TRG_1_species = sb_input %>% 
  filter(Species %in% sig_species) %>% 
  mutate(Species = factor(Species, levels = sig_species)) %>% 
  ggplot(aes(TRG_1, log10(abundance + min_abund/2))) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  facet_wrap(~Species, ncol = 3) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(0.95)),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = rel(1.05)),
        axis.text.y = element_text(size = rel(1.3)),
        strip.background = element_rect(color = NA),
        aspect.ratio = 0.2) +
  labs(x = "Tumor Regression Grade",
       y = "Log10 (Abundance + 1e-5)") ; p4_TRG_1_species

# ggsave("figure/04-3-8_TRG_1-species.svg",
#        plot = p4_TRG_1_species, width = 12, height = 9)



###############  3. Genus  ###############

gb_0.2 = gb %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Genus) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.2) %>% 
  pull(Genus) %>% 
  sort() ; gb_0.2


# Make an input file for differential enrichment testing
gb_input = gb %>% 
  filter(Genus %in% gb_0.2) %>% 
  pivot_longer(cols = -Genus,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  merge(mb2 %>% 
          select(SampleID,
                 TRG_score, TRG_1, TRG_2, TRG_3,
                 Pre_Op_Tstage, Pre_Op_Tstage_bin,
                 Pre_Op_Nstage, Pre_Op_Nstage_bin,
                 BMI, BMI_bin, CEA, CEA_bin),
        by = "SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR")))

min_abund = gb_input %>% 
  filter(abundance > 0) %>% 
  pull(abundance) %>% min() # 2e-05


# TRG_1: CR vs. nonCR
# Wilcoxon test
wilcox_g = gb_input %>% 
  group_by(Genus) %>% 
  summarise(p_value = wilcox.test(abundance ~ TRG_1, exact = F)$p.value,
            .groups = "drop") %>% 
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>% 
  arrange(p_value) ; wilcox_g

sig_genus = wilcox_g %>% 
  filter(p_value < 0.1) %>% 
  arrange(p_value) %>% 
  pull(Genus) ; sig_genus
# 14 genera, P < 0.1
# 4 genera, P < 0.05

sig_genus = gb_input %>% 
  filter(Genus %in% sig_genus) %>% 
  group_by(Genus) %>% 
  summarise(avg_abund = mean(abundance)) %>% 
  arrange(-avg_abund) %>% 
  pull(Genus) ; sig_genus

p4_TRG_1_genus = gb_input %>% 
  filter(Genus %in% sig_genus) %>% 
  mutate(Genus = factor(Genus, levels = sig_genus)) %>% 
  ggplot(aes(TRG_1, log10(abundance + min_abund/2))) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  facet_wrap(~Genus, ncol = 2) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(0.95)),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = rel(1.05)),
        axis.text.y = element_text(size = rel(1.3)),
        strip.background = element_rect(color = NA),
        aspect.ratio = 0.2) +
  labs(x = "Tumor Regression Grade",
       y = "Log10 (Abundance + pseudo)") ; p4_TRG_1_genus

# ggsave("figure/04-3-9_TRG_1-genus.svg",
#        plot = p4_TRG_1_genus, width = 8, height = 8)



###############  4. Family  ###############

# Prevalence filtering
fb_0.2 = fb %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Family) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.2) %>% 
  pull(Family) %>% 
  sort() ; fb_0.2


# Make an input file for differential enrichment testing 
fb_input = fb %>% 
  filter(Family %in% fb_0.2) %>% 
  pivot_longer(cols = -Family,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  merge(mb2 %>% 
          select(SampleID,
                 TRG_score, TRG_1, TRG_2, TRG_3,
                 Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                 Pre_Op_Nstage, Pre_Op_Nstage_bin,
                 BMI, BMI_bin, CEA, CEA_bin),
        by = "SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR")))

min_abund = fb_input %>% 
  filter(abundance > 0) %>% 
  pull(abundance) %>% min() # 9e-05


# TRG_1: CR vs. nonCR
# Wilcoxon test 
wilcox_f = fb_input %>% 
  group_by(Family) %>% 
  summarise(p_value = wilcox.test(abundance ~ TRG_1, exact = F)$p.value,
            .groups = "drop") %>% 
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>% 
  arrange(p_value) ; wilcox_f

sig_family = wilcox_f %>% 
  filter(p_value < 0.1) %>% 
  arrange(p_value) %>% 
  pull(Family) ; sig_family
# 5 families, P < 0.1
# 0 family, P < 0.05

sig_family = fb_input %>% 
  filter(Family %in% sig_family) %>% 
  group_by(Family) %>% 
  summarise(avg_abund = mean(abundance)) %>% 
  arrange(-avg_abund) %>% 
  pull(Family) ; sig_family

p4_TRG_1_family = fb_input %>% 
  filter(Family %in% sig_family) %>% 
  mutate(Family = factor(Family, levels = sig_family)) %>% 
  ggplot(aes(TRG_1, log10(abundance + min_abund/2))) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  facet_wrap(~Family, ncol = 2) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0) +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(0.95)),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = rel(1.05)),
        axis.text.y = element_text(size = rel(1.3)),
        strip.background = element_rect(color = NA),
        aspect.ratio = 0.2) +
  labs(x = "Tumor Regression Grade",
       y = "Log10 (Abundance + pseudo)") ; p4_TRG_1_family

# ggsave("figure/04-3-10_TRG_1-family.svg",
#        plot = p4_TRG_1_family, width = 8, height = 3.5)



###############  Differentially-enriched pathway (DEPs)  ###############

# Import
path_abund = read_tsv("input/merged_pathabundance.tsv")
# path_abund$Pathway %>% head(20)

# 각 pathway가 유래한 미생물 정보를 제거하여 단순화
path_abund_simple = path_abund %>% filter(!grepl("\\|", Pathway))
# head(path_abund_simple)

path_cover = read_tsv("input/merged_pathcoverage.tsv")
path_cover_simple = path_cover %>% filter(!grepl("\\|", Pathway))

# Pre-processing
# Merging pathway abundance & coverage data
path = merge(
  
  # Pathway abundance 
  path_abund_simple %>% 
    pivot_longer(cols = -Pathway,
                 names_to = "SampleID",
                 values_to = "Abundance"),
  
  # Pathway coverage
  path_cover_simple %>% 
    pivot_longer(cols = -Pathway,
                 names_to = "SampleID",
                 values_to = "Coverage"),
  
  by = c("Pathway", "SampleID")) %>% 
  merge(., mb2, by.x = "SampleID") %>% 
  mutate(
    # Multiplying
    Abund_Cov = Abundance * Coverage,
    
    # Coverage threshold
    Coverage_cat1 = ifelse(Coverage >= 0.9, 1, 0),
    Coverage_cat2 = ifelse(Coverage >= 0.5, 1, 0),
    Coverage_cat3 = ifelse(Coverage >= 0.2, 1, 0),
    Coverage_cat4 = ifelse(Coverage >= 0.1, 1, 0),
    
    Abundance_0.9 = Abundance * Coverage_cat1,
    Abundance_0.5 = Abundance * Coverage_cat2,
    Abundance_0.2 = Abundance * Coverage_cat3,
    Abundance_0.1 = Abundance * Coverage_cat4) ; head(path)


# Prevalence filtering
rowSums(path_abund_simple[, -1] > 0) %>% summary()  
# 1st Qu  17.75
# Median  39.00
# 3rd Qu  42.00
# Max     42.00


##### Prevalence filtering

# >= 50% ;  430 pathways
pathway_0.5 = path_abund_simple$Pathway[rowSums(path_abund_simple[, -1] > 0) >= 
                                          (nrow(mb2) * 0.5)]

# >= 20% ; 475 pathways
pathway_0.2 = path_abund_simple$Pathway[rowSums(path_abund_simple[,-1] > 0) >= 
                                          (nrow(mb2) * 0.2)]

##### Coverage filtering 
# 최소 샘플 1개에서는 100% coverage인 경우 = 1 / 26 ≒ 0.04
pathway_cover_0.04 =
  path_cover_simple[apply(path_cover_simple[, -1], 1, mean) >= 1/nrow(mb2), ] 

pathway_filtered = intersect(pathway_0.2, 
                             pathway_cover_0.04$Pathway[-c(1:2)]) # 145 pathways

path_sorted = path %>% 
  filter(Pathway %in% pathway_filtered) ; dim(path_sorted)  # 3770 rows

head(path_sorted)


# Abundance vs. Abundance_0.5
path_sorted %>% 
  ggplot(aes(Abundance, Abundance_0.1)) +
  geom_point() +
  geom_abline() +
  theme_classic() +
  theme(aspect.ratio = 1)

sum(path_sorted$Abundance > path_sorted$Abundance_0.9) # 1412
sum(path_sorted$Abundance > path_sorted$Abundance_0.5) # 1363
sum(path_sorted$Abundance > path_sorted$Abundance_0.2) # 1342
sum(path_sorted$Abundance > path_sorted$Abundance_0.1) # 1337
sum(path_sorted$Abundance == path_sorted$Abundance_0.1) # 2433

# 각 샘플마다 coverage가 특정 threshold를 넘지 못하는 경우, abundance를 0으로 설정
# coverage 0.5를 기준으로 할 경우 1363개의 pathway abundance가 0으로 설정
# 이를 0.1로 낮춰도 여전히 1337개가 0으로 설정됨 
# -> 엄격한 기준(0.5)과 관대한 기준(0.1)에서 큰 차이가 보이지 않음
# 이 경우 충분한 coverage를 보이는 샘플로 진행하는게 적절하다고 판단



###############  1. Identifying DEPs  ###############

# The half of Minimum abundance
epsilon = path_sorted %>% 
  filter(Abund_Cov > 0) %>% 
  pull(Abundance) %>% 
  sort(decreasing = F) %>% 
  .[[1]]/2 ; epsilon  # 3473.534/2 = 1736.767


# Performed Wilcoxon test for 145 filtered pathways 
# Abundance will be ignored if coverage in the sample is below 0.5 (using Abundance_0.5)
# Add epsilon value to calculate fold change (to avoid dividing by 0)

res_wilcox = merge(
  
  # TRG_1
  path_sorted %>% 
    group_by(Pathway) %>% 
    summarise(
      P_TRG_1 = wilcox.test(Abundance_0.5 ~ TRG_1,
                            exact = F)$p.value,
      FC_TRG_1 =
        (mean(Abundance_0.5[TRG_1 == "CR"], na.rm = T) + epsilon) /
        (mean(Abundance_0.5[TRG_1 == "nonCR"], na.rm = T) + epsilon),
      .groups = "drop"),
  
  # TRG_2
  path_sorted %>% 
    group_by(Pathway) %>% 
    summarise(P_TRG_2 = wilcox.test(Abundance_0.5 ~ TRG_2,
                                    exact = F)$p.value,
              FC_TRG_2 =
                (mean(Abundance_0.5[TRG_2 == "Regress"], na.rm = T) + epsilon) /
                (mean(Abundance_0.5[TRG_2 == "Bad"], na.rm = T) + epsilon),
              .groups = "drop"),
  by = "Pathway") %>% 
  
  merge(.,
        
        # TRG_3
        path_sorted %>% 
          group_by(Pathway) %>% 
          summarise(P_TRG_3 = wilcox.test(Abundance_0.5 ~ TRG_3,
                                          exact = F)$p.value,
                    FC_TRG_3 =
                      (mean(Abundance_0.5[TRG_3 == "NR"], na.rm = T) + epsilon) /
                      (mean(Abundance_0.5[TRG_3 == "R"], na.rm = T) + epsilon),
                    .groups = "drop"),
        by = "Pathway") %>% 
  
  mutate(P_overall = (P_TRG_1 * P_TRG_2 * P_TRG_3)^(1/3)) %>% 
  merge(.,
        
        # Add average Abundance & Coverage information
        path_sorted %>% 
          group_by(Pathway) %>% 
          summarise(Abundance = mean(Abundance),
                    Coverage = mean(Coverage),
                    Abundance_0.9 = mean(Abundance_0.9),
                    Abundance_0.5 = mean(Abundance_0.5)),
        by = "Pathway") ; head(res_wilcox)


# TRG_1
sig_TRG_1 = res_wilcox %>% 
  arrange(P_TRG_1) %>% 
  select(Pathway, 
         P_TRG_1, FC_TRG_1, 
         Abundance, Coverage) %>% 
  mutate(L2FC = log2(FC_TRG_1)) %>% 
  select(-FC_TRG_1) %>% 
  arrange(P_TRG_1) %>%
  filter(P_TRG_1 < 0.1); head(sig_TRG_1) ; dim(sig_TRG_1)

# With coverage
sig_TRG_1_pathways = sig_TRG_1 %>% 
  arrange(-Abundance) %>% 
  pull(Pathway) ; sig_TRG_1_pathways

p4_TRG_1_path_wCov = path_sorted %>% 
  filter(Pathway %in% sig_TRG_1_pathways) %>% 
  mutate(Pathway2 = factor(Pathway, levels = sig_TRG_1_pathways)) %>% 
  ggplot(aes(TRG_1, Abundance_0.5)) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  facet_wrap(~Pathway2, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     method = "wilcox.test",
                     size = rel(3.5), tip.length = 0) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = rel(0.7)),
        axis.text.y = element_text(size = rel(1.05)),
        strip.text = element_text(size = rel(1.05)),
        axis.title.y = element_blank(),
        strip.background = element_rect(color = NA)) +
  labs(y = "Functional abundance of Microbial pathways") ; p4_TRG_1_path_wCov

# ggsave("figure/04-3-11_significant-pathways_TRG_1-wCoverage.svg",
#        plot = p4_TRG_1_path_wCov, width = 14, height = 5)


unique(path[grepl("imidazole", path$Pathway), ] %>% pull(Pathway))
unique(path[grepl("propionate", path$Pathway), ] %>% pull(Pathway))
unique(path[grepl("propanoate", path$Pathway), ] %>% pull(Pathway))
unique(path[grepl("urocanate", path$Pathway), ] %>% pull(Pathway))



###############  2. Sucrose pathway  ###############

# Any pathway contain a term "sucrose"
path_sucrose = path_sorted %>% 
  filter(grepl("sucrose", Pathway, ignore.case = TRUE)) %>% 
  pull(Pathway) %>% 
  unique()

# Ordering by abundance
path_sucrose = path_sorted %>% 
  filter(Pathway %in% path_sucrose) %>% 
  group_by(Pathway) %>% 
  summarise(Abund_sum = sum(Abundance_0.5)) %>% 
  arrange(-Abund_sum) %>% 
  pull(Pathway) ; path_sucrose
# [1] "PWY-7238: sucrose biosynthesis II"                       
# [2] "PWY-5384: sucrose degradation IV (sucrose phosphorylase)"
# [3] "PWY-621: sucrose degradation III (sucrose invertase)" 


# With coverage
p4_sucrose_TRG_1_wCoverage = path_sorted %>% 
  filter(Pathway %in% path_sucrose) %>% 
  mutate(Pathway2 = factor(Pathway, levels = path_sucrose)) %>% 
  ggplot(aes(TRG_1, Abundance_0.5)) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.5, aes(fill = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  facet_wrap(~Pathway2, ncol = 1) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     method = "wilcox.test",
                     size = rel(3), tip.length = 0) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = rel(0.7)),
        axis.text.y = element_text(size = rel(1.05)),
        strip.text = element_text(size = rel(1.05)),
        strip.background = element_rect(color = NA)) +
  labs(y = "Functional abundance of microbial pathway",
       x = "Tumor Regression Grade") ; p4_sucrose_TRG_1_wCoverage

# ggsave("figure/04-3-12_sucrose-pathways_TRG1-wCoverage.svg",
#        plot = p4_sucrose_TRG_1_wCoverage, width = 6, height = 4)


path_sucrose_wide = path_sorted %>% 
  filter(Pathway %in% path_sucrose) %>% 
  pivot_wider(id_cols = c("SampleID", "TRG_1", "TRG_score"),
              names_from = Pathway,
              values_from = Abundance_0.5) %>% 
  mutate(degradation = `PWY-621: sucrose degradation III (sucrose invertase)`+
           `PWY-5384: sucrose degradation IV (sucrose phosphorylase)`)

# Spearman, p-value
path_sucrose_stat = cor.test(path_sucrose_wide$`PWY-7238: sucrose biosynthesis II`,
                             path_sucrose_wide$degradation, 
                             method = "spearman") ; path_sucrose_stat

path_sucrose_spearman = path_sucrose_stat$estimate ; path_sucrose_spearman # 0.2420154 
path_sucrose_p = path_sucrose_stat$p.value ; path_sucrose_p # 0.2335914


# 시각화: 이를 바탕으로 categorization 기준 설정
p4_sucrose_TRG_1_categorization = path_sucrose_wide %>% 
  ggplot(aes( `PWY-7238: sucrose biosynthesis II`,
              degradation)) +
  geom_point(aes(color = TRG_1)) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  annotate("text", x = min(path_sucrose_wide$`PWY-7238: sucrose biosynthesis II` + 12000), 
           y = max(path_sucrose_wide$degradation), 
           label = paste("Spearman r =", round(path_sucrose_spearman, 2), 
                         "\np-value =", round(path_sucrose_p, 4)), 
           hjust = 0, vjust = 1, size = 3, color = "black") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "PWY-621: sucrose degradation III (sucrose invertase) + \n
       PWY-5384: sucrose degradation IV (sucrose phosphorylase)"); p4_sucrose_TRG_1_categorization

# ggsave("figure/04-3-12_sucrose-TRG_1-scatter.svg",
#        plot = p4_sucrose_TRG_1_categorization, width = 5, height = 5)


path_sucrose_wide = path_sucrose_wide %>% 
  mutate(
    
    Sucrose_degradation = ifelse(
      degradation < 10000, 1, 0
    ),
    
    Sucrose_biosynthesis = ifelse(
      `PWY-7238: sucrose biosynthesis II` > 30000, 1, 0
    ),
    
    Sucrose_point = Sucrose_degradation + Sucrose_biosynthesis,
    
    Sucrose_bin = ifelse(Sucrose_point == 2, "Y", "N")
    
  )

table(path_sucrose_wide$Sucrose_bin,
      path_sucrose_wide$TRG_score)
#    0  1  2  3
# N  1  0  7  3
# Y 10  3  1  1

path_sucrose_wilcox = wilcox.test(TRG_score ~ Sucrose_bin, data = path_sucrose_wide)
path_sucrose_wilcox$p.value # 0.0009721879


library(ggbeeswarm)
library(waffle)
library(hrbrthemes)

# Make a sucrose dataframe
sucrose_df = as.data.frame(table(path_sucrose_wide$Sucrose_bin,
                                 path_sucrose_wide$TRG_score))
colnames(sucrose_df) = c("Sucrose_bin", "TRG_score", "Frequency")


# TRG_score를 숫자형으로 변환
is.numeric(sucrose_df$TRG_score) # FALSE 
sucrose_df$TRG_score = as.numeric(as.character(sucrose_df$TRG_score))
score_colors = c("0" = "#FFEDA0", 
                 "1" = "#FEB24C",
                 "2" = "#FC4E2A",
                 "3" = "#BD0026")


# waffle plot
p4_sucrose_TRG_1_waffle = sucrose_df %>% 
  group_by(Sucrose_bin) %>% 
  mutate(label = paste(TRG_score, Sucrose_bin, sep = "_")) %>% 
  ggplot(aes(fill = factor(TRG_score), values = Frequency)) +
  geom_waffle(n_rows = 5, size = 2, colour = "white", flip = T) +
  scale_fill_manual(values = score_colors,
                    name = "TRG Score") +
  facet_wrap(~Sucrose_bin, nrow = 1, scales = "fixed") +
  coord_equal() +
  annotate("text", x = 0.5, y = 3.5,
           label = paste("p-value = ", round(path_sucrose_wilcox$p.value, 4)),
           size = 3, color = "black", hjust = 0) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = rel(1.2)),
        strip.text= element_text(size = rel(2), face = "bold")) ; p4_sucrose_TRG_1_waffle

# ggsave("figure/04-3-12_sucrose-TRG_1-waffle.svg",
#        plot = p4_sucrose_TRG_1_waffle, width = 8, height = 4.5)


# Ongoing
# Any pathway contain a term "sucrose"
path_sucrose_o = path_sorted_o %>% 
  filter(grepl("sucrose", Pathway, ignore.case = TRUE)) %>% 
  pull(Pathway) %>% 
  unique()

# Ordering by abundance
path_sucrose_o = path_sorted_o %>% 
  filter(Pathway %in% path_sucrose_o) %>% 
  group_by(Pathway) %>% 
  summarise(Abund_sum = sum(Abundance_0.5)) %>% 
  arrange(-Abund_sum) %>% 
  pull(Pathway) ; path_sucrose_o
# [1] "PWY-7238: sucrose biosynthesis II"                       
# [2] "PWY-5384: sucrose degradation IV (sucrose phosphorylase)"
# [3] "PWY-621: sucrose degradation III (sucrose invertase)" 


path_sucrose_wide_o = path_sorted_o %>% 
  filter(Pathway %in% path_sucrose_o) %>% 
  pivot_wider(id_cols = c("SampleID", "TRG_1", "TRG_score"),
              names_from = Pathway,
              values_from = Abundance_0.5) %>% 
  mutate(degradation = `PWY-621: sucrose degradation III (sucrose invertase)`+
           `PWY-5384: sucrose degradation IV (sucrose phosphorylase)`)

# Spearman, p-value
path_sucrose_stat_o = cor.test(path_sucrose_wide_o$`PWY-7238: sucrose biosynthesis II`,
                               path_sucrose_wide_o$degradation, 
                               method = "spearman") ; path_sucrose_stat_o

path_sucrose_spearman_o = path_sucrose_stat_o$estimate ; path_sucrose_spearman_o # 0.2333533   
path_sucrose_p_o = path_sucrose_stat_o$p.value ; path_sucrose_p_o # 0.3844116


# 시각화: 이를 바탕으로 categorization 기준 설정
p5_sucrose_TRG_1_categorization = path_sucrose_wide_o %>% 
  ggplot(aes( `PWY-7238: sucrose biosynthesis II`,
              degradation)) +
  geom_point(aes(color = TRG_1)) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  annotate("text", x = min(path_sucrose_wide$`PWY-7238: sucrose biosynthesis II` + 12000), 
           y = max(path_sucrose_wide_o$degradation), 
           label = paste("Spearman r =", round(path_sucrose_spearman_o, 2), 
                         "\np-value =", round(path_sucrose_p_o, 4)), 
           hjust = 0, vjust = 1, size = 3, color = "black") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y = "PWY-621: sucrose degradation III (sucrose invertase) + \n
       PWY-5384: sucrose degradation IV (sucrose phosphorylase)"); p5_sucrose_TRG_1_categorization

# ggsave("figure/05-3-12_sucrose-TRG_1-scatter_ongoing.svg",
#        plot = p5_sucrose_TRG_1_categorization, width = 5, height = 5)


path_sucrose_wide_o = path_sucrose_wide_o %>% 
  mutate(
    
    Sucrose_degradation = ifelse(
      degradation < 10000, 1, 0
    ),
    
    Sucrose_biosynthesis = ifelse(
      `PWY-7238: sucrose biosynthesis II` > 40000, 1, 0
    ),
    
    Sucrose_point = Sucrose_degradation + Sucrose_biosynthesis,
    
    Sucrose_bin = ifelse(Sucrose_point == 2, "Y", "N")
    
  )

table(path_sucrose_wide_o$Sucrose_bin,
      path_sucrose_wide_o$TRG_score)
#   0 1 2 3
# N 5 1 4 1
# Y 2 1 1 1

path_sucrose_wilcox_o = wilcox.test(TRG_score ~ Sucrose_bin, data = path_sucrose_wide_o)
path_sucrose_wilcox_o$p.value # 0.9041631


library(ggbeeswarm)
library(waffle)
library(hrbrthemes)

# Make a sucrose dataframe
sucrose_df_o = as.data.frame(table(path_sucrose_wide_o$Sucrose_bin,
                                 path_sucrose_wide_o$TRG_score))
colnames(sucrose_df_o) = c("Sucrose_bin", "TRG_score", "Frequency")


# TRG_score를 숫자형으로 변환
is.numeric(sucrose_df_o$TRG_score) # FALSE 
sucrose_df_o$TRG_score = as.numeric(as.character(sucrose_df_o$TRG_score))
score_colors = c("0" = "#FFEDA0", 
                 "1" = "#FEB24C",
                 "2" = "#FC4E2A",
                 "3" = "#BD0026")


# waffle plot
p5_sucrose_TRG_1_waffle = sucrose_df_o %>% 
  group_by(Sucrose_bin) %>% 
  mutate(label = paste(TRG_score, Sucrose_bin, sep = "_")) %>% 
  ggplot(aes(fill = factor(TRG_score), values = Frequency)) +
  geom_waffle(n_rows = 5, size = 2, colour = "white", flip = T) +
  scale_fill_manual(values = score_colors,
                    name = "TRG Score") +
  facet_wrap(~Sucrose_bin, nrow = 1, scales = "fixed") +
  coord_equal() +
  annotate("text", x = 0.5, y = 3.5,
           label = paste("p-value = ", round(path_sucrose_wilcox_o$p.value, 4)),
           size = 3, color = "black", hjust = 0) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = rel(1.2)),
        strip.text= element_text(size = rel(2), face = "bold")) ; p5_sucrose_TRG_1_waffle

# ggsave("figure/05-3-12_sucrose-TRG_1-waffle_ongoing.svg",
#        plot = p5_sucrose_TRG_1_waffle, width = 8, height = 4.5)



###############  3. Histidine pathway  ###############

# Any pathway contain a term "histidine"
path_histidine = path_sorted %>% 
  filter(grepl("histidine", Pathway, ignore.case = TRUE)) %>% 
  pull(Pathway) %>% 
  unique()

# Ordering by abundance
path_histidine = path_sorted %>% 
  filter(Pathway %in% path_histidine) %>% 
  group_by(Pathway) %>% 
  summarise(Abund_sum = sum(Abundance_0.5)) %>% 
  arrange(-Abund_sum) %>% 
  pull(Pathway) ; path_histidine
# [1] "HISTSYN-PWY: L-histidine biosynthesis"
# [2] "HISDEG-PWY: L-histidine degradation I"


# With coverage
p4_histidine_TRG_1_wCoverage = path_sorted %>% 
  filter(Pathway %in% path_histidine) %>% 
  mutate(Pathway2 = factor(Pathway, levels = path_histidine)) %>% 
  ggplot(aes(TRG_1, Abundance_0.5)) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.5, aes(fill = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  facet_wrap(~Pathway2, ncol = 1) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     method = "wilcox.test",
                     size = rel(3), tip.length = 0) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = rel(0.7)),
        axis.text.y = element_text(size = rel(1.05)),
        strip.text = element_text(size = rel(1.05)),
        strip.background = element_rect(color = NA)) +
  labs(y = "Functional abundance of microbial pathway",
       x = "Tumor Regression Grade") ; p4_histidine_TRG_1_wCoverage

# ggsave("figure/04-3-13_histidine-pathways_TRG1-wCoverage.svg",
#        plot = p4_histidine_TRG_1_wCoverage, width = 6, height = 4)


path_histidine_wide = path_sorted %>% 
  filter(Pathway %in% path_histidine) %>% 
  pivot_wider(id_cols = c("SampleID", "TRG_1", "TRG_score"),
              names_from = Pathway,
              values_from = Abundance_0.5)

# Spearman, p-value
path_histidine_stat = cor.test(path_histidine_wide$`HISTSYN-PWY: L-histidine biosynthesis`,
                               path_histidine_wide$`HISDEG-PWY: L-histidine degradation I`, 
                               method = "spearman") ; path_histidine_stat

path_histidine_spearman = path_histidine_stat$estimate ; path_histidine_spearman # -0.4909452 
path_histidine_p = path_histidine_stat$p.value ; path_histidine_p # 0.01087472


# 시각화: 이를 바탕으로 categorization 기준 설정
p4_histidine_TRG_1_categorization = path_histidine_wide %>% 
  ggplot(aes( `HISTSYN-PWY: L-histidine biosynthesis`,
              `HISDEG-PWY: L-histidine degradation I`)) +
  geom_point(aes(color = TRG_1)) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  annotate("text", x = min(path_histidine_wide$`HISTSYN-PWY: L-histidine biosynthesis`), 
           y = max(path_histidine_wide$`HISDEG-PWY: L-histidine degradation I`), 
           label = paste("Spearman r =", round(path_histidine_spearman, 2), 
                         "\np-value =", round(path_histidine_p, 4)), 
           hjust = 0, vjust = 1, size = 3, color = "black") +
  theme_classic() +
  theme(legend.position = "none") ; p4_histidine_TRG_1_categorization

# ggsave("figure/04-3-13_histidine-TRG_1-scatter.svg",
#        plot = p4_histidine_TRG_1_categorization, width = 5, height = 5)


path_histidine_wide = path_histidine_wide %>% 
  mutate(
    
    Histidine_degradation = ifelse(
      `HISDEG-PWY: L-histidine degradation I` < 5000, 1, 0
    ),
    
    Histidine_biosynthesis = ifelse(
      `HISTSYN-PWY: L-histidine biosynthesis` > 15000, 1, 0
    ),
    
    Histidine_point = Histidine_degradation + Histidine_biosynthesis,
    
    Histidine_bin = ifelse(Histidine_point == 2, "Y", "N")
    
  )

table(path_histidine_wide$Histidine_bin,
      path_histidine_wide$TRG_score)
#    0  1  2  3
# N  0  1  3  1
# Y 11  2  5  3

path_histidine_wilcox = wilcox.test(TRG_score ~ Histidine_bin, data = path_histidine_wide)
path_histidine_wilcox$p.value # 0.08491633


library(ggbeeswarm)
library(waffle)
library(hrbrthemes)


# Make a histidine dataframe
histidine_df = as.data.frame(table(path_histidine_wide$Histidine_bin,
                                   path_histidine_wide$TRG_score))
colnames(histidine_df) = c("Histidine_bin", "TRG_score", "Frequency")


# TRG_score를 숫자형으로 변환
is.numeric(histidine_df$TRG_score) # FALSE 
histidine_df$TRG_score = as.numeric(as.character(histidine_df$TRG_score))
score_colors = c("0" = "#FFEDA0", 
                 "1" = "#FEB24C",
                 "2" = "#FC4E2A",
                 "3" = "#BD0026")


# waffle plot
p4_histidine_TRG_1_waffle = histidine_df %>% 
  group_by(Histidine_bin) %>% 
  mutate(label = paste(TRG_score, Histidine_bin, sep = "_")) %>% 
  ggplot(aes(fill = factor(TRG_score), values = Frequency)) +
  geom_waffle(n_rows = 5, size = 2, colour = "white", flip = T) +
  scale_fill_manual(values = score_colors,
                    name = "TRG Score") +
  facet_wrap(~Histidine_bin, nrow = 1, scales = "fixed") +
  coord_equal() +
  annotate("text", x = 0.5, y = 3.5,
           label = paste("p-value = ", round(path_histidine_wilcox$p.value, 4)),
           size = 3, color = "black", hjust = 0) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = rel(1.2)),
        strip.text= element_text(size = rel(2), face = "bold")) ; p4_histidine_TRG_1_waffle

# ggsave("figure/04-3-13_histidine-TRG_1-waffle.svg",
#        plot = p4_histidine_TRG_1_waffle, width = 8, height = 4.5)



###############  4. Thiamine pathway  ###############

# Any pathway contain a term "thiamine"
path_thiamine = path_sorted %>% 
  filter(grepl("thiamine", Pathway, ignore.case = TRUE)) %>% 
  pull(Pathway) %>% 
  unique()

# Ordering by abundance
path_thiamine = path_sorted %>% 
  filter(Pathway %in% path_thiamine) %>% 
  group_by(Pathway) %>% 
  summarise(Abund_sum = sum(Abundance_0.5)) %>% 
  arrange(-Abund_sum) %>% 
  pull(Pathway) ; path_thiamine
# [1] "PWY-7357: thiamine phosphate formation from pyrithiamine and oxythiamine (yeast)" 
# [2] "PWY-6892: thiazole component of thiamine diphosphate biosynthesis I"              
# [3] "PWY-6897: thiamine diphosphate salvage II"                                        
# [4] "THISYNARA-PWY: superpathway of thiamine diphosphate biosynthesis III (eukaryotes)"
# [5] "THISYN-PWY: superpathway of thiamine diphosphate biosynthesis I"


# With coverage
p4_thiamine_TRG_1_wCoverage = path_sorted %>% 
  filter(Pathway %in% path_thiamine) %>% 
  mutate(Pathway2 = factor(Pathway, levels = path_thiamine)) %>% 
  ggplot(aes(TRG_1, Abundance_0.5)) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.5, aes(fill = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  facet_wrap(~Pathway2, ncol = 1) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     method = "wilcox.test",
                     size = rel(3), tip.length = 0) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = rel(0.7)),
        axis.text.y = element_text(size = rel(1.05)),
        strip.text = element_text(size = rel(1.05)),
        strip.background = element_rect(color = NA)) +
  labs(y = "Functional abundance of microbial pathway",
       x = "Tumor Regression Grade") ; p4_thiamine_TRG_1_wCoverage

# ggsave("figure/04-3-14_thiamine-pathways_TRG1-wCoverage.svg",
#        plot = p4_thiamine_TRG_1_wCoverage, width = 6, height = 6)


path_thiamine_wide = path_sorted %>% 
  filter(Pathway %in% path_thiamine) %>% 
  pivot_wider(id_cols = c("SampleID", "TRG_1", "TRG_score"),
              names_from = Pathway,
              values_from = Abundance_0.5)

# Spearman, p-value
path_thiamine_stat = cor.test(path_thiamine_wide$`PWY-6897: thiamine diphosphate salvage II`,
                              path_thiamine_wide$`PWY-6892: thiazole component of thiamine diphosphate biosynthesis I`,
                              method = "spearman") ; path_thiamine_stat

path_thiamine_spearman = path_thiamine_stat$estimate ; path_thiamine_spearman # -0.3149162 
path_thiamine_p = path_thiamine_stat$p.value ; path_thiamine_p # 0.1171231


# 시각화: 이를 바탕으로 categorization 기준 설정
p4_thiamine_TRG_1_categorization = path_thiamine_wide %>% 
  ggplot(aes( `PWY-6897: thiamine diphosphate salvage II`,
              `PWY-6892: thiazole component of thiamine diphosphate biosynthesis I`)) +
  geom_point(aes(color = TRG_1)) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  annotate("text", x = max(path_thiamine_wide$`PWY-6897: thiamine diphosphate salvage II`-2500), 
           y = max(path_thiamine_wide$`PWY-6892: thiazole component of thiamine diphosphate biosynthesis I`), 
           label = paste("Spearman r =", round(path_thiamine_spearman, 2), 
                         "\np-value =", round(path_thiamine_p, 4)), 
           hjust = 0, vjust = 1, size = 3, color = "black") +
  theme_classic() +
  theme(legend.position = "none") ; p4_thiamine_TRG_1_categorization

# ggsave("figure/04-3-14_thiamine-TRG_1-scatter.svg",
#        plot = p4_thiamine_TRG_1_categorization, width = 5, height = 5)


path_thiamine_wide = path_thiamine_wide %>% 
  mutate(
    
    Thiamine_salvage = ifelse(
      `PWY-6897: thiamine diphosphate salvage II` > 10000, 1, 0
    ),
    
    Thiazole_biosynthesis = ifelse(
      `PWY-6892: thiazole component of thiamine diphosphate biosynthesis I` < 10000, 1, 0
    ),
    
    Thiamine_point = Thiamine_salvage + Thiazole_biosynthesis,
    
    Thiamine_bin = ifelse(Thiamine_point == 2, "Y", "N")
    
  )

table(path_thiamine_wide$Thiamine_bin,
      path_thiamine_wide$TRG_score)
#   0 1 2 3
# N 6 3 8 4
# Y 5 0 0 0

path_thiamine_wilcox = wilcox.test(TRG_score ~ Thiamine_bin, data = path_thiamine_wide)
path_thiamine_wilcox$p.value # 0.01077809


library(ggbeeswarm)
library(waffle)
library(hrbrthemes)


# Make a thiamine dataframe
thiamine_df = as.data.frame(table(path_thiamine_wide$Thiamine_bin,
                                  path_thiamine_wide$TRG_score))
colnames(thiamine_df) = c("Thiamine_bin", "TRG_score", "Frequency")


# TRG_score를 숫자형으로 변환
is.numeric(thiamine_df$TRG_score) # FALSE 
thiamine_df$TRG_score = as.numeric(as.character(thiamine_df$TRG_score))
score_colors = c("0" = "#FFEDA0", 
                 "1" = "#FEB24C",
                 "2" = "#FC4E2A",
                 "3" = "#BD0026")


# waffle plot
p4_thiamine_TRG_1_waffle = thiamine_df %>% 
  group_by(Thiamine_bin) %>% 
  mutate(label = paste(TRG_score, Thiamine_bin, sep = "_")) %>% 
  ggplot(aes(fill = factor(TRG_score), values = Frequency)) +
  geom_waffle(n_rows = 5, size = 2, colour = "white", flip = T) +
  scale_fill_manual(values = score_colors,
                    name = "TRG Score") +
  facet_wrap(~Thiamine_bin, nrow = 1, scales = "fixed") +
  coord_equal() +
  annotate("text", x = 0.5, y = 6.0,
           label = paste("p-value = ", round(path_thiamine_wilcox$p.value, 4)),
           size = 3, color = "black", hjust = 0) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = rel(1.2)),
        strip.text= element_text(size = rel(2), face = "bold")) ; p4_thiamine_TRG_1_waffle

# ggsave("figure/04-3-14_thiamine-TRG_1-waffle.svg",
#        plot = p4_thiamine_TRG_1_waffle, width = 8, height = 4.5)



# Ongoing

# Any pathway contain a term "thiamine"
path_thiamine_o = path_sorted_o %>% 
  filter(grepl("thiamine", Pathway, ignore.case = TRUE)) %>% 
  pull(Pathway) %>% 
  unique()

# Ordering by abundance
path_thiamine_o = path_sorted_o %>% 
  filter(Pathway %in% path_thiamine_o) %>% 
  group_by(Pathway) %>% 
  summarise(Abund_sum = sum(Abundance_0.5)) %>% 
  arrange(-Abund_sum) %>% 
  pull(Pathway) ; path_thiamine_o
# [1] "PWY-7357: thiamine phosphate formation from pyrithiamine and oxythiamine (yeast)" 
# [2] "PWY-6892: thiazole component of thiamine diphosphate biosynthesis I"              
# [3] "PWY-6897: thiamine diphosphate salvage II"                                        
# [4] "THISYNARA-PWY: superpathway of thiamine diphosphate biosynthesis III (eukaryotes)"


path_thiamine_wide_o = path_sorted_o %>% 
  filter(Pathway %in% path_thiamine_o) %>% 
  pivot_wider(id_cols = c("SampleID", "TRG_1", "TRG_score"),
              names_from = Pathway,
              values_from = Abundance_0.5)

# Spearman, p-value
path_thiamine_stat_o = cor.test(path_thiamine_wide_o$`PWY-6897: thiamine diphosphate salvage II`,
                                path_thiamine_wide_o$`PWY-6892: thiazole component of thiamine diphosphate biosynthesis I`,
                                method = "spearman") ; path_thiamine_stat_o

path_thiamine_spearman_o = path_thiamine_stat_o$estimate ; path_thiamine_spearman_o # -0.687189 
path_thiamine_p_o = path_thiamine_stat_o$p.value ; path_thiamine_p_o # 0.003269681


# 시각화: 이를 바탕으로 categorization 기준 설정
p5_thiamine_TRG_1_categorization = path_thiamine_wide_o %>% 
  ggplot(aes( `PWY-6897: thiamine diphosphate salvage II`,
              `PWY-6892: thiazole component of thiamine diphosphate biosynthesis I`)) +
  geom_point(aes(color = TRG_1)) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  annotate("text", x = max(path_thiamine_wide_o$`PWY-6897: thiamine diphosphate salvage II`-4000), 
           y = max(path_thiamine_wide_o$`PWY-6892: thiazole component of thiamine diphosphate biosynthesis I`), 
           label = paste("Spearman r =", round(path_thiamine_spearman_o, 2), 
                         "\np-value =", round(path_thiamine_p_o, 4)), 
           hjust = 0, vjust = 1, size = 3, color = "black") +
  theme_classic() +
  theme(legend.position = "none") ; p5_thiamine_TRG_1_categorization

# ggsave("figure/05-3-14_thiamine-TRG_1-scatter_ongoing.svg",
#        plot = p5_thiamine_TRG_1_categorization, width = 5, height = 5)


path_thiamine_wide_o = path_thiamine_wide_o %>% 
  mutate(
    
    Thiamine_salvage = ifelse(
      `PWY-6897: thiamine diphosphate salvage II` > 8000, 1, 0
    ),
    
    Thiazole_biosynthesis = ifelse(
      `PWY-6892: thiazole component of thiamine diphosphate biosynthesis I` < 15000, 1, 0
    ),
    
    Thiamine_point = Thiamine_salvage + Thiazole_biosynthesis,
    
    Thiamine_bin = ifelse(Thiamine_point == 2, "Y", "N")
    
  )

table(path_thiamine_wide_o$Thiamine_bin,
      path_thiamine_wide_o$TRG_score)
#   0 1 2 3
# N 3 2 4 1
# Y 4 0 1 1

path_thiamine_wilcox_o = wilcox.test(TRG_score ~ Thiamine_bin, data = path_thiamine_wide_o)
path_thiamine_wilcox_o$p.value # 0.4196973


library(ggbeeswarm)
library(waffle)
library(hrbrthemes)


# Make a thiamine dataframe
thiamine_df_o = as.data.frame(table(path_thiamine_wide_o$Thiamine_bin,
                                    path_thiamine_wide_o$TRG_score))
colnames(thiamine_df_o) = c("Thiamine_bin", "TRG_score", "Frequency")


# TRG_score를 숫자형으로 변환
is.numeric(thiamine_df_o$TRG_score) # FALSE 
thiamine_df_o$TRG_score = as.numeric(as.character(thiamine_df_o$TRG_score))
score_colors = c("0" = "#FFEDA0", 
                 "1" = "#FEB24C",
                 "2" = "#FC4E2A",
                 "3" = "#BD0026")


# waffle plot
p5_thiamine_TRG_1_waffle = thiamine_df_o %>% 
  group_by(Thiamine_bin) %>% 
  mutate(label = paste(TRG_score, Thiamine_bin, sep = "_")) %>% 
  ggplot(aes(fill = factor(TRG_score), values = Frequency)) +
  geom_waffle(n_rows = 5, size = 2, colour = "white", flip = T) +
  scale_fill_manual(values = score_colors,
                    name = "TRG Score") +
  facet_wrap(~Thiamine_bin, nrow = 1, scales = "fixed") +
  coord_equal() +
  annotate("text", x = 0.5, y = 3.0,
           label = paste("p-value = ", round(path_thiamine_wilcox_o$p.value, 4)),
           size = 3, color = "black", hjust = 0) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = rel(1.2)),
        strip.text= element_text(size = rel(2), face = "bold")) ; p5_thiamine_TRG_1_waffle

# ggsave("figure/05-3-14_thiamine-TRG_1-waffle_ongoing.svg",
#        plot = p5_thiamine_TRG_1_waffle, width = 8, height = 4.5)



###############  5. Rhamnose pathway  ###############

# Any pathway contain a term "rhamnose"
path_rhamnose = path_sorted %>% 
  filter(grepl("rhamnose", Pathway, ignore.case = TRUE)) %>% 
  pull(Pathway) %>% 
  unique()

# Ordering by abundance
path_rhamnose = path_sorted %>% 
  filter(Pathway %in% path_rhamnose) %>% 
  group_by(Pathway) %>% 
  summarise(Abund_sum = sum(Abundance_0.5)) %>% 
  arrange(-Abund_sum) %>% 
  pull(Pathway) ; path_rhamnose
# [1] "DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis"
# [2] "RHAMCAT-PWY: L-rhamnose degradation I"


# With coverage
p4_rhamnose_TRG_1_wCoverage = path_sorted %>% 
  filter(Pathway %in% path_rhamnose) %>% 
  mutate(Pathway2 = factor(Pathway, levels = path_rhamnose)) %>% 
  ggplot(aes(TRG_1, Abundance_0.5)) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.5, aes(fill = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  facet_wrap(~Pathway2, ncol = 1) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     method = "wilcox.test",
                     size = rel(3), tip.length = 0) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = rel(0.7)),
        axis.text.y = element_text(size = rel(1.05)),
        strip.text = element_text(size = rel(1.05)),
        strip.background = element_rect(color = NA)) +
  labs(y = "Functional abundance of microbial pathway",
       x = "Tumor Regression Grade") ; p4_rhamnose_TRG_1_wCoverage

# ggsave("figure/04-3-15_rhamnose-pathways_TRG1-wCoverage.svg",
#        plot = p4_rhamnose_TRG_1_wCoverage, width = 6, height = 3)


path_rhamnose_wide = path_sorted %>% 
  filter(Pathway %in% path_rhamnose) %>% 
  pivot_wider(id_cols = c("SampleID", "TRG_1", "TRG_score"),
              names_from = Pathway,
              values_from = Abundance_0.5)

# Spearman, p-value
path_rhamnose_stat = cor.test(path_rhamnose_wide$`DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis`,
                              path_rhamnose_wide$`RHAMCAT-PWY: L-rhamnose degradation I`,
                              method = "spearman") ; path_rhamnose_stat

path_rhamnose_spearman = path_rhamnose_stat$estimate ; path_rhamnose_spearman # -0.1236053 
path_rhamnose_p = path_rhamnose_stat$p.value ; path_rhamnose_p # 0.5474517


# 시각화: 이를 바탕으로 categorization 기준 설정
p4_rhamnose_TRG_1_categorization = path_rhamnose_wide %>% 
  ggplot(aes (`DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis`,
              `RHAMCAT-PWY: L-rhamnose degradation I`)) +
  geom_point(aes(color = TRG_1)) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  annotate("text", x = min(path_rhamnose_wide$`DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis`), 
           y = max(path_rhamnose_wide$`RHAMCAT-PWY: L-rhamnose degradation I`), 
           label = paste("Spearman r =", round(path_rhamnose_spearman, 2), 
                         "\np-value =", round(path_rhamnose_p, 4)), 
           hjust = 0, vjust = 1, size = 3, color = "black") +
  theme_classic() +
  theme(legend.position = "none") ; p4_rhamnose_TRG_1_categorization

# ggsave("figure/04-3-15_rhamnose-TRG_1-scatter.svg",
#        plot = p4_rhamnose_TRG_1_categorization, width = 5, height = 5)


path_rhamnose_wide = path_rhamnose_wide %>% 
  mutate(
    
    Rhamnose_biosynthesis = ifelse(
      `DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis` < 35000, 1, 0
    ),
    
    Rhamnose_degradation = ifelse(
      `RHAMCAT-PWY: L-rhamnose degradation I` > 5000, 1, 0
    ),
    
    Rhamnose_point = Rhamnose_biosynthesis + Rhamnose_degradation,
    
    Rhamnose_bin = ifelse(Rhamnose_point == 2, "Y", "N")
    
  )

table(path_rhamnose_wide$Rhamnose_bin,
      path_rhamnose_wide$TRG_score)
#   0 1 2 3
# N 5 1 8 4
# Y 6 2 0 0

path_rhamnose_wilcox = wilcox.test(TRG_score ~ Rhamnose_bin, data = path_rhamnose_wide)
path_rhamnose_wilcox$p.value # 0.007417135


library(ggbeeswarm)
library(waffle)
library(hrbrthemes)


# Make a rhamnose dataframe
rhamnose_df = as.data.frame(table(path_rhamnose_wide$Rhamnose_bin,
                                  path_rhamnose_wide$TRG_score))
colnames(rhamnose_df) = c("Rhamnose_bin", "TRG_score", "Frequency")


# TRG_score를 숫자형으로 변환
is.numeric(rhamnose_df$TRG_score) #FALSE
rhamnose_df$TRG_score = as.numeric(as.character(rhamnose_df$TRG_score))
score_colors = c("0" = "#FFEDA0", 
                 "1" = "#FEB24C",
                 "2" = "#FC4E2A",
                 "3" = "#BD0026")


# waffle plot
p4_rhamnose_TRG_1_waffle = rhamnose_df %>% 
  group_by(Rhamnose_bin) %>% 
  mutate(label = paste(TRG_score, Rhamnose_bin, sep = "_")) %>% 
  ggplot(aes(fill = factor(TRG_score), values = Frequency)) +
  geom_waffle(n_rows = 5, size = 2, colour = "white", flip = T) +
  scale_fill_manual(values = score_colors,
                    name = "TRG score") +
  facet_wrap(~Rhamnose_bin, nrow = 1, scales = "fixed") +
  coord_equal() +
  annotate("text", x = 0.5, y = 5.0,
           label = paste("p-value = ", round(path_rhamnose_wilcox$p.value, 4)),
           size = 3, color = "black", hjust = 0) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = rel(1.2)),
        strip.text= element_text(size = rel(2), face = "bold")) ; p4_rhamnose_TRG_1_waffle

# ggsave("figure/04-3-15_rhamnose-TRG_1-waffle.svg",
#        plot = p4_rhamnose_TRG_1_waffle, width = 8, height = 4.5)


# Ongoing
# Any pathway contain a term "rhamnose"
path_rhamnose_o = path_sorted_o %>% 
  filter(grepl("rhamnose", Pathway, ignore.case = TRUE)) %>% 
  pull(Pathway) %>% 
  unique()

# Ordering by abundance
path_rhamnose_o = path_sorted_o %>% 
  filter(Pathway %in% path_rhamnose_o) %>% 
  group_by(Pathway) %>% 
  summarise(Abund_sum = sum(Abundance_0.5)) %>% 
  arrange(-Abund_sum) %>% 
  pull(Pathway) ; path_rhamnose_o
# [1] "DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis"
# [2] "RHAMCAT-PWY: L-rhamnose degradation I"


# With coverage
p5_rhamnose_TRG_1_wCoverage = path_sorted_o %>% 
  filter(Pathway %in% path_rhamnose_o) %>% 
  mutate(Pathway2 = factor(Pathway, levels = path_rhamnose_o)) %>% 
  ggplot(aes(TRG_1, Abundance_0.5)) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.5, aes(fill = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  facet_wrap(~Pathway2, ncol = 1) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     method = "wilcox.test",
                     size = rel(3), tip.length = 0) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = rel(0.7)),
        axis.text.y = element_text(size = rel(1.05)),
        strip.text = element_text(size = rel(1.05)),
        strip.background = element_rect(color = NA)) +
  labs(y = "Functional abundance of microbial pathway",
       x = "Tumor Regression Grade") ; p5_rhamnose_TRG_1_wCoverage

# ggsave("figure/05-3-15_rhamnose-pathways_TRG1-wCoverage.svg",
#        plot = p5_rhamnose_TRG_1_wCoverage, width = 6, height = 3)


path_rhamnose_wide_o = path_sorted_o %>% 
  filter(Pathway %in% path_rhamnose_o) %>% 
  pivot_wider(id_cols = c("SampleID", "TRG_1", "TRG_score"),
              names_from = Pathway,
              values_from = Abundance_0.5)

# Spearman, p-value
path_rhamnose_stat_o = cor.test(path_rhamnose_wide_o$`DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis`,
                                path_rhamnose_wide_o$`RHAMCAT-PWY: L-rhamnose degradation I`,
                                method = "spearman") ; path_rhamnose_stat_o

path_rhamnose_spearman_o = path_rhamnose_stat_o$estimate ; path_rhamnose_spearman_o # -0.2762698  
path_rhamnose_p_o = path_rhamnose_stat_o$p.value ; path_rhamnose_p_o # 0.3003076


# 시각화: 이를 바탕으로 categorization 기준 설정
p5_rhamnose_TRG_1_categorization = path_rhamnose_wide_o %>% 
  ggplot(aes (`DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis`,
              `RHAMCAT-PWY: L-rhamnose degradation I`)) +
  geom_point(aes(color = TRG_1)) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  annotate("text", x = min(path_rhamnose_wide_o$`DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis`), 
           y = max(path_rhamnose_wide_o$`RHAMCAT-PWY: L-rhamnose degradation I`), 
           label = paste("Spearman r =", round(path_rhamnose_spearman_o, 2), 
                         "\np-value =", round(path_rhamnose_p_o, 4)), 
           hjust = 0, vjust = 1, size = 3, color = "black") +
  theme_classic() +
  theme(legend.position = "none") ; p5_rhamnose_TRG_1_categorization

# ggsave("figure/05-3-15_rhamnose-TRG_1-scatter.svg",
#        plot = p5_rhamnose_TRG_1_categorization, width = 5, height = 5)


path_rhamnose_wide_o = path_rhamnose_wide_o %>% 
  mutate(
    
    Rhamnose_biosynthesis = ifelse(
      `DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis` < 26000, 1, 0
    ),
    
    Rhamnose_degradation = ifelse(
      `RHAMCAT-PWY: L-rhamnose degradation I` > 6000, 1, 0
    ),
    
    Rhamnose_point = Rhamnose_biosynthesis + Rhamnose_degradation,
    
    Rhamnose_bin = ifelse(Rhamnose_point == 2, "Y", "N")
    
  )

table(path_rhamnose_wide_o$Rhamnose_bin,
      path_rhamnose_wide_o$TRG_score)
#   0 1 2 3
# N 6 1 4 2
# Y 1 1 1 0

path_rhamnose_wilcox_o = wilcox.test(TRG_score ~ Rhamnose_bin, data = path_rhamnose_wide_o)
path_rhamnose_wilcox_o$p.value # 0.9430058


library(ggbeeswarm)
library(waffle)
library(hrbrthemes)


# Make a rhamnose dataframe
rhamnose_df_o = as.data.frame(table(path_rhamnose_wide_o$Rhamnose_bin,
                                    path_rhamnose_wide_o$TRG_score))
colnames(rhamnose_df_o) = c("Rhamnose_bin", "TRG_score", "Frequency")


# TRG_score를 숫자형으로 변환
is.numeric(rhamnose_df_o$TRG_score) #FALSE
rhamnose_df_o$TRG_score = as.numeric(as.character(rhamnose_df_o$TRG_score))
score_colors = c("0" = "#FFEDA0", 
                 "1" = "#FEB24C",
                 "2" = "#FC4E2A",
                 "3" = "#BD0026")


# waffle plot
p5_rhamnose_TRG_1_waffle = rhamnose_df_o %>% 
  group_by(Rhamnose_bin) %>% 
  mutate(label = paste(TRG_score, Rhamnose_bin, sep = "_")) %>% 
  ggplot(aes(fill = factor(TRG_score), values = Frequency)) +
  geom_waffle(n_rows = 5, size = 2, colour = "white", flip = T) +
  scale_fill_manual(values = score_colors,
                    name = "TRG score") +
  facet_wrap(~Rhamnose_bin, nrow = 1, scales = "fixed") +
  coord_equal() +
  annotate("text", x = 0.5, y = 4.0,
           label = paste("p-value = ", round(path_rhamnose_wilcox_o$p.value, 4)),
           size = 3, color = "black", hjust = 0) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = rel(1.2)),
        strip.text= element_text(size = rel(2), face = "bold")) ; p5_rhamnose_TRG_1_waffle

# ggsave("figure/05-3-15_rhamnose-TRG_1-waffle.svg",
#        plot = p5_rhamnose_TRG_1_waffle, width = 8, height = 4.5)



###############  6. Isoprene pathway  ###############

# Any pathway contain a term "isoprene"
path_isoprene = path_sorted %>% 
  filter(grepl("isoprene", Pathway, ignore.case = TRUE)) %>% 
  pull(Pathway) %>% 
  unique()

# Ordering by abundance
path_isoprene = path_sorted %>% 
  filter(Pathway %in% path_isoprene) %>% 
  group_by(Pathway) %>% 
  summarise(Abund_sum = sum(Abundance_0.5)) %>% 
  arrange(-Abund_sum) %>% 
  pull(Pathway) ; path_isoprene
# [1] "PWY-6270: isoprene biosynthesis I"


# With coverage
p4_isoprene_TRG_1_wCoverage = path_sorted %>% 
  filter(Pathway %in% path_isoprene) %>% 
  mutate(Pathway2 = factor(Pathway, levels = path_isoprene)) %>% 
  ggplot(aes(TRG_1, Abundance_0.5)) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.5, aes(fill = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  facet_wrap(~Pathway2, ncol = 1) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     method = "wilcox.test",
                     size = rel(3), tip.length = 0) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = rel(0.7)),
        axis.text.y = element_text(size = rel(1.05)),
        strip.text = element_text(size = rel(1.05)),
        strip.background = element_rect(color = NA)) +
  labs(y = "Functional abundance of microbial pathway",
       x = "Tumor Regression Grade") ; p4_isoprene_TRG_1_wCoverage

# ggsave("figure/04-3-16_isoprene-pathways_TRG1-wCoverage.svg",
#        plot = p4_isoprene_TRG_1_wCoverage, width = 6, height = 3)


path_isoprene_wide = path_sorted %>% 
  filter(Pathway %in% path_isoprene) %>% 
  pivot_wider(id_cols = c("SampleID", "TRG_1", "TRG_score"),
              names_from = Pathway,
              values_from = Abundance_0.5)


# 시각화: 이를 바탕으로 categorization 기준 설정
p4_isoprene_TRG_1_categorization = path_isoprene_wide %>% 
  ggplot(aes (1:nrow(path_isoprene_wide),`PWY-6270: isoprene biosynthesis I`)) +
  geom_point(aes(color = TRG_1)) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  theme_classic() +
  labs(x = "Sample #") +
  theme(legend.position = "none") ; p4_isoprene_TRG_1_categorization

# ggsave("figure/04-3-16_isoprene-TRG_1-scatter.svg",
#        plot = p4_isoprene_TRG_1_categorization, width = 5, height = 5)


path_isoprene_wide = path_isoprene_wide %>% 
  mutate(
    
    Isoprene_biosynthesis = ifelse(
      `PWY-6270: isoprene biosynthesis I` > 15000, 1, 0
    ),
    
    Isoprene_point = Isoprene_biosynthesis,
    
    Isoprene_bin = ifelse(Isoprene_point == 1, "Y", "N")
    
  )

table(path_isoprene_wide$Isoprene_bin,
      path_isoprene_wide$TRG_score)
#   0 1 2 3
# N 4 0 1 1
# Y 7 3 7 3

path_isoprene_wilcox = wilcox.test(TRG_score ~ Isoprene_bin, data = path_isoprene_wide)
path_isoprene_wilcox$p.value # 0.3841635


library(ggbeeswarm)
library(waffle)
library(hrbrthemes)


# Make a rhamnose dataframe
isoprene_df = as.data.frame(table(path_isoprene_wide$Isoprene_bin,
                                  path_isoprene_wide$TRG_score))
colnames(isoprene_df) = c("Isoprene_bin", "TRG_score", "Frequency")


# TRG_score를 숫자형으로 변환
is.numeric(isoprene_df$TRG_score) #FALSE
isoprene_df$TRG_score = as.numeric(as.character(isoprene_df$TRG_score))
score_colors = c("0" = "#FFEDA0", 
                 "1" = "#FEB24C",
                 "2" = "#FC4E2A",
                 "3" = "#BD0026")


# waffle plot
p4_isoprene_TRG_1_waffle = isoprene_df %>% 
  group_by(Isoprene_bin) %>% 
  mutate(label = paste(TRG_score, Isoprene_bin, sep = "_")) %>% 
  ggplot(aes(fill = factor(TRG_score), values = Frequency)) +
  geom_waffle(n_rows = 5, size = 2, colour = "white", flip = T) +
  scale_fill_manual(values = score_colors,
                    name = "TRG score") +
  facet_wrap(~Isoprene_bin, nrow = 1, scales = "fixed") +
  coord_equal() +
  annotate("text", x = 0.5, y = 5.0,
           label = paste("p-value = ", round(path_isoprene_wilcox$p.value, 4)),
           size = 3, color = "black", hjust = 0) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = rel(1.2)),
        strip.text= element_text(size = rel(2), face = "bold")) ; p4_isoprene_TRG_1_waffle

# ggsave("figure/04-3-16_isoprene-TRG_1-waffle.svg",
#        plot = p4_isoprene_TRG_1_waffle, width = 8, height = 4.5)



###############  7. Methylerythritol phosphate (MEP) pathway  ###############

# Any pathway contain a term "methylerythritol"
path_methylerythritol = path_sorted %>% 
  filter(grepl("methylerythritol", Pathway, ignore.case = TRUE)) %>% 
  pull(Pathway) %>% 
  unique()

# Ordering by abundance
path_methylerythritol = path_sorted %>% 
  filter(Pathway %in% path_methylerythritol) %>% 
  group_by(Pathway) %>% 
  summarise(Abund_sum = sum(Abundance_0.5)) %>% 
  arrange(-Abund_sum) %>% 
  pull(Pathway) ; path_methylerythritol
# [1] "NONMEVIPP-PWY: methylerythritol phosphate pathway I"
# [2] "PWY-7560: methylerythritol phosphate pathway II" 


# With coverage
p4_methylerythritol_TRG_1_wCoverage = path_sorted %>% 
  filter(Pathway %in% path_methylerythritol) %>% 
  mutate(Pathway2 = factor(Pathway, levels = path_methylerythritol)) %>% 
  ggplot(aes(TRG_1, Abundance_0.5)) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.5, aes(fill = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  facet_wrap(~Pathway2, ncol = 1) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     method = "wilcox.test",
                     size = rel(3), tip.length = 0) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = rel(0.7)),
        axis.text.y = element_text(size = rel(1.05)),
        strip.text = element_text(size = rel(1.05)),
        strip.background = element_rect(color = NA)) +
  labs(y = "Functional abundance of microbial pathway",
       x = "Tumor Regression Grade") ; p4_methylerythritol_TRG_1_wCoverage

# ggsave("figure/04-3-17_methylerythritol-pathways_TRG1-wCoverage.svg",
#        plot = p4_methylerythritol_TRG_1_wCoverage, width = 6, height = 3)


path_methylerythritol_wide = path_sorted %>% 
  filter(Pathway %in% path_methylerythritol) %>% 
  pivot_wider(id_cols = c("SampleID", "TRG_1", "TRG_score"),
              names_from = Pathway,
              values_from = Abundance_0.5)

# Spearman, p-value
path_methylerythritol_stat = cor.test(path_methylerythritol_wide$`NONMEVIPP-PWY: methylerythritol phosphate pathway I`,
                                      path_methylerythritol_wide$`PWY-7560: methylerythritol phosphate pathway II`,
                                      method = "spearman") ; path_methylerythritol_stat

path_methylerythritol_spearman = path_methylerythritol_stat$estimate ; path_methylerythritol_spearman # 0.3996593
path_methylerythritol_p = path_methylerythritol_stat$p.value ; path_methylerythritol_p # 0.04308998


# 시각화: 이를 바탕으로 categorization 기준 설정
p4_methylerythritol_TRG_1_categorization = path_methylerythritol_wide %>% 
  ggplot(aes (`NONMEVIPP-PWY: methylerythritol phosphate pathway I`,
              `PWY-7560: methylerythritol phosphate pathway II`)) +
  geom_point(aes(color = TRG_1)) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  annotate("text", x = min(path_methylerythritol_wide$`NONMEVIPP-PWY: methylerythritol phosphate pathway I`), 
           y = max(path_methylerythritol_wide$`PWY-7560: methylerythritol phosphate pathway II`), 
           label = paste("Spearman r =", round(path_methylerythritol_spearman, 2), 
                         "\np-value =", round(path_methylerythritol_p, 4)), 
           hjust = 0, vjust = 1, size = 3, color = "black") +
  theme_classic() +
  theme(legend.position = "none") ; p4_methylerythritol_TRG_1_categorization

# ggsave("figure/04-3-17_methylerythritolTRG_1-scatter.svg",
#        plot = p4_methylerythritol_TRG_1_categorization, width = 5, height = 5)


path_methylerythritol_wide = path_methylerythritol_wide %>% 
  mutate(
    
    Methylerythritol_biosynthesis = ifelse(
      `NONMEVIPP-PWY: methylerythritol phosphate pathway I` > 20000, 1, 0
    ),
    
    Methylerythritol_degradation = ifelse(
      `PWY-7560: methylerythritol phosphate pathway II` > 13000, 1, 0
    ),
    
    Methylerythritol_point = Methylerythritol_biosynthesis + Methylerythritol_degradation,
    
    Methylerythritol_bin = ifelse(Methylerythritol_point == 2, "Y", "N")
    
  )

table(path_methylerythritol_wide$Methylerythritol_bin,
      path_methylerythritol_wide$TRG_score)
#   0 1 2 3
# N 5 0 1 0
# Y 6 3 7 4

path_methylerythritol_wilcox = wilcox.test(TRG_score ~ Methylerythritol_bin, data = path_methylerythritol_wide)
path_methylerythritol_wilcox$p.value # 0.0391301


library(ggbeeswarm)
library(waffle)
library(hrbrthemes)

# Make a sucrose dataframe
methylerythritol_df = as.data.frame(table(path_methylerythritol_wide$Methylerythritol_bin,
                                          path_methylerythritol_wide$TRG_score))
colnames(methylerythritol_df) = c("Methylerythritol_bin", "TRG_score", "Frequency")


# TRG_score를 숫자형으로 변환
is.numeric(methylerythritol_df$TRG_score) # FALSE 
methylerythritol_df$TRG_score = as.numeric(as.character(methylerythritol_df$TRG_score))
score_colors = c("0" = "#FFEDA0", 
                 "1" = "#FEB24C",
                 "2" = "#FC4E2A",
                 "3" = "#BD0026")


# waffle plot
p4_methylerythritol_TRG_1_waffle = methylerythritol_df %>% 
  group_by(Methylerythritol_bin) %>% 
  mutate(label = paste(TRG_score, Methylerythritol_bin, sep = "_")) %>% 
  ggplot(aes(fill = factor(TRG_score), values = Frequency)) +
  geom_waffle(n_rows = 5, size = 2, colour = "white", flip = T) +
  scale_fill_manual(values = score_colors,
                    name = "TRG Score") +
  facet_wrap(~Methylerythritol_bin, nrow = 1, scales = "fixed") +
  coord_equal() +
  annotate("text", x = 0.5, y = 5.0,
           label = paste("p-value = ", round(path_methylerythritol_wilcox$p.value, 4)),
           size = 3, color = "black", hjust = 0) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = rel(1.2)),
        strip.text= element_text(size = rel(2), face = "bold")) ; p4_methylerythritol_TRG_1_waffle

# ggsave("figure/04-3-17_methylerythritol-TRG_1-waffle.svg",
#        plot = p4_methylerythritol_TRG_1_waffle, width = 8, height = 4.5)



### In ongoing samples

# Pre-processing
path_o = merge(
  
  # Pathway abundance 
  path_abund_simple %>% 
    pivot_longer(cols = -Pathway,
                 names_to = "SampleID",
                 values_to = "Abundance"),
  
  # Pathway coverage
  path_cover_simple %>% 
    pivot_longer(cols = -Pathway,
                 names_to = "SampleID",
                 values_to = "Coverage"),
  
  by = c("Pathway", "SampleID")) %>% 
  merge(., mo, by.x = "SampleID") %>% 
  mutate(
    # Multiplying
    Abund_Cov = Abundance * Coverage,
    
    # Coverage threshold
    Coverage_cat1 = ifelse(Coverage >= 0.9, 1, 0),
    Coverage_cat2 = ifelse(Coverage >= 0.5, 1, 0),
    Coverage_cat3 = ifelse(Coverage >= 0.2, 1, 0),
    Coverage_cat4 = ifelse(Coverage >= 0.1, 1, 0),
    
    Abundance_0.9 = Abundance * Coverage_cat1,
    Abundance_0.5 = Abundance * Coverage_cat2,
    Abundance_0.2 = Abundance * Coverage_cat3,
    Abundance_0.1 = Abundance * Coverage_cat4) ; head(path_o)


# Prevalence filtering
rowSums(path_abund_simple[, -1] > 0) %>% summary()  
# 1st Qu  17.75
# Median  39.00
# 3rd Qu  42.00
# Max     42.00


##### Prevalence filtering

# >= 50% ;  455 pathways
pathway_0.5_o = path_abund_simple$Pathway[rowSums(path_abund_simple[, -1] > 0) >= 
                                            (nrow(mo) * 0.5)]

# >= 20% ; 486 pathways
pathway_0.2_o = path_abund_simple$Pathway[rowSums(path_abund_simple[,-1] > 0) >= 
                                            (nrow(mo) * 0.2)]

##### Coverage filtering 
# 최소 샘플 1개에서는 100% coverage인 경우
pathway_cover_0.04_o =
  path_cover_simple[apply(path_cover_simple[, -1], 1, mean) >= 1/nrow(mo), ] 

pathway_filtered_o = intersect(pathway_0.2_o, 
                               pathway_cover_0.04_o$Pathway[-c(1:2)]) # 136 pathways

path_sorted_o = path_o %>% 
  filter(Pathway %in% pathway_filtered_o) ; dim(path_sorted_o)  # 2176 rows

head(path_sorted_o)


# Abundance vs. Abundance_0.5
path_sorted_o %>% 
  ggplot(aes(Abundance, Abundance_0.1)) +
  geom_point() +
  geom_abline() +
  theme_classic() +
  theme(aspect.ratio = 1)

sum(path_sorted_o$Abundance > path_sorted_o$Abundance_0.9) # 814
sum(path_sorted_o$Abundance > path_sorted_o$Abundance_0.5) # 796
sum(path_sorted_o$Abundance > path_sorted_o$Abundance_0.2) # 792
sum(path_sorted_o$Abundance > path_sorted_o$Abundance_0.1) # 789
sum(path_sorted_o$Abundance == path_sorted_o$Abundance_0.1) # 1387


##### Identifying DEPs

# The half of minimum abundance
epsilon = path_sorted_o %>% 
  filter(Abund_Cov > 0) %>% 
  pull(Abundance) %>% 
  sort(decreasing = F) %>% 
  .[[1]]/2 ; epsilon  # 5256.234/2 = 2628.117


# Performed Wilcoxon test for 136 filtered pathways 
# Abundance will be ignored if coverage in the sample is below 0.5 (using Abundance_0.5)
# Add epsilon value to calculate fold change (to avoid dividing by 0)

res_wilcox_o = merge(
  
  # TRG_1
  path_sorted_o %>% 
    group_by(Pathway) %>% 
    summarise(
      P_TRG_1 = wilcox.test(Abundance_0.5 ~ TRG_1,
                            exact = F)$p.value,
      FC_TRG_1 =
        (mean(Abundance_0.5[TRG_1 == "CR"], na.rm = T) + epsilon) /
        (mean(Abundance_0.5[TRG_1 == "nonCR"], na.rm = T) + epsilon),
      .groups = "drop"),
  
  # TRG_2
  path_sorted_o %>% 
    group_by(Pathway) %>% 
    summarise(P_TRG_2 = wilcox.test(Abundance_0.5 ~ TRG_2,
                                    exact = F)$p.value,
              FC_TRG_2 =
                (mean(Abundance_0.5[TRG_2 == "Regress"], na.rm = T) + epsilon) /
                (mean(Abundance_0.5[TRG_2 == "Bad"], na.rm = T) + epsilon),
              .groups = "drop"),
  by = "Pathway") %>% 
  
  merge(.,
        
        # TRG_3
        path_sorted_o %>% 
          group_by(Pathway) %>% 
          summarise(P_TRG_3 = wilcox.test(Abundance_0.5 ~ TRG_3,
                                          exact = F)$p.value,
                    FC_TRG_3 =
                      (mean(Abundance_0.5[TRG_3 == "NR"], na.rm = T) + epsilon) /
                      (mean(Abundance_0.5[TRG_3 == "R"], na.rm = T) + epsilon),
                    .groups = "drop"),
        by = "Pathway") %>% 
  
  mutate(P_overall = (P_TRG_1 * P_TRG_2 * P_TRG_3)^(1/3)) %>% 
  merge(.,
        
        # Add average Abundance & Coverage information
        path_sorted_o %>% 
          group_by(Pathway) %>% 
          summarise(Abundance = mean(Abundance),
                    Coverage = mean(Coverage),
                    Abundance_0.9 = mean(Abundance_0.9),
                    Abundance_0.5 = mean(Abundance_0.5)),
        by = "Pathway") ; head(res_wilcox_o)


sig_TRG_1_o = res_wilcox_o %>% 
  arrange(P_TRG_1) %>% 
  select(Pathway, 
         P_TRG_1, FC_TRG_1, 
         Abundance, Coverage) %>% 
  mutate(L2FC = log2(FC_TRG_1)) %>% 
  select(-FC_TRG_1) %>% 
  arrange(P_TRG_1) %>%
  filter(P_TRG_1 < 0.1); head(sig_TRG_1_o) ; dim(sig_TRG_1_o)

#   Pathway                                                      P_TRG_1    Abundance Coverage  L2FC
# 1 PWY66-429: fatty acid biosynthesis initiation (mitochondria) 0.05673822 22326.205 0.9999978 -0.1220855
# 2         PWY-621: sucrose degradation III (sucrose invertase) 0.06047561  8353.675 0.2360563 -1.9043715
# 3       SALVADEHYPOX-PWY: adenosine nucleotides degradation II 0.06047561  8849.156 0.2500000 -2.0952317
# 4                                 PWY-5686: UMP biosynthesis I 0.07194424 27678.286 0.9375000  0.2043975
# 5                                PWY-7790: UMP biosynthesis II 0.07194424 27678.286 0.9375000  0.2043975
# 6                               PWY-7791: UMP biosynthesis III 0.07194424 27678.286 0.9375000  0.2043975
# 7            PWY-7851: coenzyme A biosynthesis II (eukaryotic) 0.09033759 22192.033 0.9375000  0.2304592


# With coverage
sig_TRG_1_pathways_o = sig_TRG_1_o %>% 
  arrange(-Abundance) %>% 
  pull(Pathway) ; sig_TRG_1_pathways_o
# [1] "PWY-5686: UMP biosynthesis I"                                 
# [2] "PWY-7790: UMP biosynthesis II"                               
# [3] "PWY-7791: UMP biosynthesis III"                              
# [4] "PWY66-429: fatty acid biosynthesis initiation (mitochondria)" 
# [5] "PWY-7851: coenzyme A biosynthesis II (eukaryotic)"            
# [6] "SALVADEHYPOX-PWY: adenosine nucleotides degradation II"      
# [7] "PWY-621: sucrose degradation III (sucrose invertase)"   


p5_TRG_1_path_wCov_o = path_sorted_o %>% 
  filter(Pathway %in% sig_TRG_1_pathways_o) %>% 
  mutate(Pathway2 = factor(Pathway, levels = sig_TRG_1_pathways_o)) %>% 
  ggplot(aes(TRG_1, Abundance_0.5)) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  facet_wrap(~Pathway2, nrow = 2, scales = "free") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     method = "wilcox.test",
                     size = rel(3.5), tip.length = 0) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = rel(0.7)),
        axis.text.y = element_text(size = rel(1.05)),
        strip.text = element_text(size = rel(1.05)),
        axis.title.y = element_blank(),
        strip.background = element_rect(color = NA)) +
  labs(y = "Functional abundance of Microbial pathways") ; p5_TRG_1_path_wCov_o

# ggsave("figure/05-3-2_significant-pathways_TRG_1-wCoverage_ongoing.svg",
#        plot = p5_TRG_1_path_wCov_o, width = 12, height = 8)


p5_TRG_1_path_coa_wCov_o = path_sorted_o %>% 
  filter(Pathway == "COA-PWY: coenzyme A biosynthesis I (prokaryotic)") %>% 
  ggplot(aes(TRG_1, Abundance_0.5)) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  facet_wrap(~Pathway, nrow = 2, scales = "free") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     method = "wilcox.test",
                     size = rel(3.5), tip.length = 0) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = rel(0.7)),
        axis.text.y = element_text(size = rel(1.05)),
        strip.text = element_text(size = rel(1.05)),
        axis.title.y = element_blank(),
        strip.background = element_rect(color = NA)) +
  labs(y = "Functional abundance of Microbial pathways") ; p5_TRG_1_path_coa_wCov_o

# ggsave("figure/05-3-2_CoA-PWY_TRG_1-wCoverage_ongoing.svg",
#        plot = p5_TRG_1_path_coa_wCov_o, width = 3, height = 4)



### In paired samples

path_p = merge(
  path_abund_simple %>% 
    pivot_longer(cols = -Pathway,
                 names_to = "SampleID",
                 values_to = "Abundance"),
  
  path_cover_simple %>% 
    pivot_longer(cols = -Pathway,
                 names_to = "SampleID",
                 values_to = "Coverage"),
  by = c("Pathway", "SampleID")) %>% 
  merge(., mp, by.x = "SampleID") %>% 
  mutate(
    # Multiplying
    Abund_Cov = Abundance * Coverage,
    
    # Coverage threshold
    Coverage_cat1 = ifelse(Coverage >= 0.9, 1, 0),
    Coverage_cat2 = ifelse(Coverage >= 0.5, 1, 0),
    Coverage_cat3 = ifelse(Coverage >= 0.2, 1, 0),
    Coverage_cat4 = ifelse(Coverage >= 0.1, 1, 0),
    
    Abundance_0.9 = Abundance * Coverage_cat1,
    Abundance_0.5 = Abundance * Coverage_cat2,
    Abundance_0.2 = Abundance * Coverage_cat3,
    Abundance_0.1 = Abundance * Coverage_cat4) ; head(path_p)


# TRG_1
# Prevalence ≥ 20%, Coverage ≥ 0.5, 최소 샘플 1개에서는 coverage 100%, P < 0.2
path_abund_simple %>% 
  filter(Pathway %in% c(res_wilcox %>% filter(P_TRG_1 < 0.2) %>% arrange(P_TRG_1) %>% pull(Pathway), 
                        res_wilcox_o %>% filter(P_TRG_1 < 0.2) %>% arrange(P_TRG_1) %>% pull(Pathway))) %>% 
  pivot_longer(cols = -Pathway,
               names_to = "SampleID",
               values_to = "Abund_Cov") %>% 
  merge(., mp, by.x = "SampleID") %>% 
  group_by(Pathway, TRG_1, TNT) %>% 
  summarise(Abund_Cov = mean(Abund_Cov, na.rm = T), .groups = "drop") %>% 
  pivot_wider(names_from = c(TNT, TRG_1), values_from = Abund_Cov) %>% 
  mutate(
    # Before에서 CR / noncR
    Before_FC = log10((`Before_CR`+1e-05) / (`Before_nonCR`+1e-05)),
    
    # Ongoing에서 CR/ nonCR
    Ongoing_FC = log10((`Ongoing_CR`+1e-05) / (`Ongoing_nonCR`+1e-05)),
    
    # CR에서 Ongoing / Before
    CR_FC = log10((`Ongoing_CR`+1e-05) / (`Before_CR`+1e-05)),
    
    # nonCR에서 Ongoing / Before
    nonCR_FC = log10((`Ongoing_nonCR`+1e-05) / (`Before_nonCR`+1e-05))
  ) %>% 
  select(Pathway, Before_CR, Before_nonCR, Before_FC, Ongoing_CR, Ongoing_nonCR, Ongoing_FC, CR_FC, nonCR_FC)


path_input = res_wilcox %>% 
  filter(P_TRG_1 < 0.2) %>% 
  pull(Pathway)

# path_input = path_input %>%
#   group_by(Pathway, TRG_1, TNT) %>%
#   summarise(total_abund_cov = sum(Abund_Cov + epsilon), .groups = "drop") %>%  # pseudo-count 추가
#   pivot_wider(names_from = TNT, values_from = total_abund_cov) %>%
#   mutate(
#     log10FC = log10((Ongoing + epsilon) / (Before + epsilon))
#   ) %>%
#   select(Pathway, TRG_1, log10FC) %>%
#   pivot_wider(names_from = TRG_1, values_from = log10FC, names_prefix = "log10FC_") %>%
#   mutate(Pathway = fct_reorder(Pathway, log10FC_CR)) %>%
#   pivot_longer(cols = starts_with("log10FC_"), names_to = "Group", values_to = "log10FC") %>%
#   mutate(Group = recode(Group, log10FC_CR = "CR", log10FC_nonCR = "nonCR"))

# p5_Fold_change_path = path_input %>% 
#   mutate(Group = factor(Group, levels = c("nonCR", "CR"))) %>% 
#   ggplot(aes(log10FC, Pathway, fill = Group)) +
#   geom_bar(stat = "identity", color = "black",
#            position = position_dodge(width = 1)) +
#   geom_errorbar(aes(xmin = log10FC - ci95, xmax = log10FC + ci95),
#                 width = 0.2, position = position_dodge(width = 1)) +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   scale_fill_manual(values = c("CR" = "#F7D9BC",
#                                "nonCR" = "#80461B"),
#                     breaks = c("CR", "nonCR"),
#                     labels = c("CR", "nonCR")) +
#   labs(title = "Log10 Fold Change",
#        x = "Log10 Fold Change (Ongoing / Before)",
#        y = "Pathway",
#        fill = "TRG_1") +
#   theme_classic() ; p5_Fold_change_path

# ggsave("figure/06-3-2_FC-pathway_paired (total).svg",
#        plot = p5_Fold_change_path, width = 10, height = 7)

path_input = path_p %>% 
  filter(Pathway %in% path_input) %>% 
  select(SNU_ID, TRG_1, TNT, Abund_Cov, Pathway) %>%
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR"))) %>% 
  pivot_wider(names_from = TNT,
              values_from = Abund_Cov) %>% 
  mutate(log10FC = log10((Ongoing + 1736.767) / (Before + 1736.767))) %>% 
  group_by(Pathway, TRG_1)

path_input_summary = path_input %>% 
  summarise(
    mean_log10FC = mean(log10FC, na.rm = T),
    sd_log10FC = sd(log10FC, na.rm = T),
    n = n(),
    se_log10FC = sd_log10FC / sqrt(n),
    ci95 = 1.96 * se_log10FC,
    .groups = "drop") %>% 
  arrange(desc(TRG_1 == "CR"), desc(mean_log10FC))

path_input_wilcox = path_input %>% 
  group_by(Pathway) %>% 
  summarise(
    wilcox_test = list(wilcox.test(log10FC ~ TRG_1, data = cur_data(), exact = FALSE)),
    p_value = wilcox_test[[1]]$p.value,
    statistic = wilcox_test[[1]]$statistic,
    .groups = "drop"
  ) %>% 
  select(Pathway, p_value, statistic) %>% 
  arrange(p_value)

path_input_summary = path_input_summary %>% 
  left_join(path_input_wilcox %>% select(Pathway, p_value), by = "Pathway") ; rm(path_input_wilcox)

p5_Fold_change_path = path_input_summary %>% 
  filter(mean_log10FC != 0) %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("nonCR", "CR")),
         Pathway = factor(Pathway, levels = rev(unique(path_input_summary$Pathway)))) %>% 
  ggplot(aes(mean_log10FC, Pathway, fill = TRG_1)) +
  geom_errorbar(aes(xmin = mean_log10FC - ci95, xmax = mean_log10FC + ci95),
                width = 0.2, position = position_dodge(width = 1), alpha = 0.5) +
  geom_point(aes(fill = TRG_1), shape = 21, size = 4,
             stat = "identity", position = position_dodge(width = 1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B"),
                    breaks = c("CR", "nonCR"),
                    labels = c("CR", "nonCR")) +
  labs(title = "Log10 Fold Change",
       x = "Log10 Fold Change (Ongoing / Before)",
       y = "Pathway",
       fill = "TRG_1") +
  theme_classic() ; p5_Fold_change_path

# ggsave("figure/06-3-2_FC-pathway_paired (mean).svg",
#        plot = p5_Fold_change_path, width = 12, height = 7)


# Shared pathways
# prev ≥ 0.1
path_p %>%
  rowwise() %>%
  group_by(Pathway) %>%
  summarise(prevalence = sum(Abund_Cov > 0), .groups = "drop") %>%
  filter(prevalence >= nrow(m) * 0.1) %>% 
  left_join(path_p, by = "Pathway") %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR"))) %>%
  group_by(Pathway, TNT, TRG_1) %>%
  summarise(mean_abund_cov = mean(Abund_Cov, na.rm = T),
            prevalence = sum(Abund_Cov > 0, na.rm = T)) %>%
  pivot_wider(names_from = c("TNT", "TRG_1"),
              values_from = c(mean_abund_cov, prevalence),
              names_glue = "{TNT}_{TRG_1}_{.value}") %>%
  write_csv(., "pathway_paired (25.11.19).csv")


# Prevalent pathways (Prev. ≥ 8)
prev_pathway = c("PWY-7560: methylerythritol phosphate pathway II",
                 "HISTSYN-PWY: L-histidine biosynthesis",
                 "PWY-6270: isoprene biosynthesis I",
                 "PWY-6892: thiazole component of thiamine diphosphate biosynthesis I")

path_prev = path_p %>% 
  filter(Pathway %in% prev_pathway)

# PWY-7560: methylerythritol phosphate pathway II
p5_PWY7560 = path_prev %>% 
  filter(Pathway == "PWY-7560: methylerythritol phosphate pathway II") %>% 
  mutate(TRG = factor(TRG_1, levels = c("CR", "nonCR")),
         TNT = factor(TNT, levels = c("Before", "Ongoing")),
         Group = paste(TNT, TRG_1, sep = "_"),
         Group = factor(Group, levels = c("Before_CR", "Ongoing_CR", 
                                          "Before_nonCR", "Ongoing_nonCR"))) %>% 
  ggplot(aes(Group, Abund_Cov)) +
  geom_point(aes(color = TRG_1)) +
  geom_line(aes(group = SNU_ID), color = "gray60", alpha = 0.5,
            arrow = arrow(type = "open", length = unit(0.05, "inches"))) +
  geom_boxplot(aes(fill = TRG_1),
               alpha = 0.3, outlier.alpha = 0,
               width = 0.6) +
  stat_compare_means(comparisons = list(c("Before_CR", "Ongoing_CR"),
                                        c("Before_nonCR", "Ongoing_nonCR")),
                     tip.length = 0.02,
                     method = "wilcox", paired = T) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Sample", y = "Abundance", 
       title = "PWY-7560: methylerythritol phosphate pathway II") ; p5_PWY7560

# HISTSYN-PWY: L-histidine biosynthesis
p5_HISTSYN_PWY = path_prev %>% 
  filter(Pathway == "HISTSYN-PWY: L-histidine biosynthesis") %>% 
  mutate(TRG = factor(TRG_1, levels = c("CR", "nonCR")),
         TNT = factor(TNT, levels = c("Before", "Ongoing")),
         Group = paste(TNT, TRG_1, sep = "_"),
         Group = factor(Group, levels = c("Before_CR", "Ongoing_CR", 
                                          "Before_nonCR", "Ongoing_nonCR"))) %>% 
  ggplot(aes(Group, Abund_Cov)) +
  geom_point(aes(color = TRG_1)) +
  geom_line(aes(group = SNU_ID), color = "gray60", alpha = 0.5,
            arrow = arrow(type = "open", length = unit(0.05, "inches"))) +
  geom_boxplot(aes(fill = TRG_1),
               alpha = 0.3, outlier.alpha = 0,
               width = 0.6) +
  stat_compare_means(comparisons = list(c("Before_CR", "Ongoing_CR"),
                                        c("Before_nonCR", "Ongoing_nonCR")),
                     tip.length = 0.02,
                     method = "wilcox", paired = T) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Sample", y = "Abundance", 
       title = "HISTSYN-PWY: L-histidine biosynthesis") ; p5_HISTSYN_PWY

# PWY-6270: isoprene biosynthesis I
p5_PWY6270 = path_prev %>% 
  filter(Pathway == "PWY-6270: isoprene biosynthesis I") %>% 
  mutate(TRG = factor(TRG_1, levels = c("CR", "nonCR")),
         TNT = factor(TNT, levels = c("Before", "Ongoing")),
         Group = paste(TNT, TRG_1, sep = "_"),
         Group = factor(Group, levels = c("Before_CR", "Ongoing_CR", 
                                          "Before_nonCR", "Ongoing_nonCR"))) %>% 
  ggplot(aes(Group, Abund_Cov)) +
  geom_point(aes(color = TRG_1)) +
  geom_line(aes(group = SNU_ID), color = "gray60", alpha = 0.5,
            arrow = arrow(type = "open", length = unit(0.05, "inches"))) +
  geom_boxplot(aes(fill = TRG_1),
               alpha = 0.3, outlier.alpha = 0,
               width = 0.6) +
  stat_compare_means(comparisons = list(c("Before_CR", "Ongoing_CR"),
                                        c("Before_nonCR", "Ongoing_nonCR")),
                     tip.length = 0.02,
                     method = "wilcox", paired = T) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Sample", y = "Abundance", 
       title = "PWY-6270: isoprene biosynthesis I") ; p5_PWY6270

# PWY-6892: thiazole component of thiamine diphosphate biosynthesis I
p5_PWY6892 = path_prev %>% 
  filter(Pathway == "PWY-6892: thiazole component of thiamine diphosphate biosynthesis I") %>% 
  mutate(TRG = factor(TRG_1, levels = c("CR", "nonCR")),
         TNT = factor(TNT, levels = c("Before", "Ongoing")),
         Group = paste(TNT, TRG_1, sep = "_"),
         Group = factor(Group, levels = c("Before_CR", "Ongoing_CR", 
                                          "Before_nonCR", "Ongoing_nonCR"))) %>% 
  ggplot(aes(Group, Abund_Cov)) +
  geom_point(aes(color = TRG_1)) +
  geom_line(aes(group = SNU_ID), color = "gray60", alpha = 0.5,
            arrow = arrow(type = "open", length = unit(0.05, "inches"))) +
  geom_boxplot(aes(fill = TRG_1),
               alpha = 0.3, outlier.alpha = 0,
               width = 0.6) +
  stat_compare_means(comparisons = list(c("Before_CR", "Ongoing_CR"),
                                        c("Before_nonCR", "Ongoing_nonCR")),
                     tip.length = 0.02,
                     method = "wilcox", paired = T) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Sample", y = "Abundance", 
       title = "PWY-6892: thiazole component of thiamine diphosphate biosynthesis I") ; p5_PWY6892

p5_paired_pathways = grid.arrange(p5_PWY7560,
                                  p5_HISTSYN_PWY,
                                  p5_PWY6270,
                                  p5_PWY6892,
                                  ncol = 4, nrow = 1) ; p5_paired_pathways

# ggsave("figure/06-3-2_tendency of pathways.svg",
#        plot = p5_paired_pathways, width = 16, height = 5)


### Selected pathways
path_select = c("PWY-5384: sucrose degradation IV (sucrose phosphorylase)",
                "HISTSYN-PWY: L-histidine biosynthesis",
                "RHAMCAT-PWY: L-rhamnose degradation I",
                "PWY-2941: L-lysine biosynthesis II",
                "GLCMANNANAUT-PWY: superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation")

path_select_input = path_p %>% 
  filter(Pathway %in% path_select) %>% 
  select(Pathway, SNU_ID, TRG_1, TNT, Abund_Cov) %>% 
  pivot_wider(names_from = TNT, values_from = Abund_Cov)


# PWY-5384: sucrose degradation IV (sucrose phosphorylase)
p5_PWY_5384 = path_select_input %>% 
  filter(Pathway == "PWY-5384: sucrose degradation IV (sucrose phosphorylase)") %>% 
  ggplot(aes(x = Before, y = Ongoing)) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.7,
              linetype = "dashed", size = 0.4) +
  geom_point(aes(fill = TRG_1), 
             size = 3, shape = 21, color = "gray40") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  theme_pubr() +
  labs(x = "Before", y = "Ongoing", color = "TRG_1",
       title = "PWY-5384: sucrose degradation IV (sucrose phosphorylase)") ; p5_PWY_5384

# HISTSYN-PWY: L-histidine biosynthesis
p5_HISTSYN_PWY = path_select_input %>% 
  filter(Pathway == "HISTSYN-PWY: L-histidine biosynthesis") %>% 
  ggplot(aes(x = Before, y = Ongoing)) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.7,
              linetype = "dashed", size = 0.4) +
  geom_point(aes(fill = TRG_1), 
             size = 3, shape = 21, color = "gray40") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  theme_pubr() +
  labs(x = "Before", y = "Ongoing", color = "TRG_1",
       title = "HISTSYN-PWY: L-histidine biosynthesis") ; p5_HISTSYN_PWY

# RHAMCAT-PWY: L-rhamnose degradation I
p5_RHAMCAT_PWY = path_select_input %>% 
  filter(Pathway == "RHAMCAT-PWY: L-rhamnose degradation I") %>% 
  ggplot(aes(x = Before, y = Ongoing)) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.7,
              linetype = "dashed", size = 0.4) +
  geom_point(aes(fill = TRG_1), 
             size = 3, shape = 21, color = "gray40") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  theme_pubr() +
  labs(x = "Before", y = "Ongoing", color = "TRG_1",
       title = "RHAMCAT-PWY: L-rhamnose degradation I") ; p5_RHAMCAT_PWY

# PWY-2941: L-lysine biosynthesis II
p5_PWY_2941 = path_select_input %>% 
  filter(Pathway == "PWY-2941: L-lysine biosynthesis II") %>% 
  ggplot(aes(x = Before, y = Ongoing)) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.7,
              linetype = "dashed", size = 0.4) +
  geom_point(aes(fill = TRG_1), 
             size = 3, shape = 21, color = "gray40") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  theme_pubr() +
  labs(x = "Before", y = "Ongoing", color = "TRG_1",
       title = "PWY-2941: L-lysine biosynthesis II") ; p5_PWY_2941

# GLCMANNANAUT-PWY: superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation
p5_GLCMANNANAUT_PWY = path_select_input %>% 
  filter(Pathway == "GLCMANNANAUT-PWY: superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation") %>% 
  ggplot(aes(x = Before, y = Ongoing)) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.7,
              linetype = "dashed", size = 0.4) +
  geom_point(aes(fill = TRG_1), 
             size = 3, shape = 21, color = "gray40") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  theme_pubr() +
  labs(x = "Before", y = "Ongoing", color = "TRG_1",
       title = "GLCMANNANAUT-PWY: superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation") ; p5_GLCMANNANAUT_PWY

p5_selected_path = ggarrange(p5_PWY_5384,
                             p5_HISTSYN_PWY,
                             p5_RHAMCAT_PWY,
                             p5_PWY_2941,
                             p5_GLCMANNANAUT_PWY,
                             nrow = 2, ncol = 3) ; p5_selected_path

# ggsave("figure/06-3-2_selected pathways_scatter.svg",
#        plot = p5_selected_path, width = 15, height = 10)


selected_shared_path = c("PWY-621: sucrose degradation III (sucrose invertase)",
                         "PWY-6892: thiazole component of thiamine diphosphate biosynthesis I",
                         "PANTO-PWY: phosphopantothenate biosynthesis I")

path_select_shared_input = path_p %>% 
  filter(Pathway %in% selected_shared_path) %>% 
  select(Pathway, SNU_ID, TRG_1, TNT, Abund_Cov) %>% 
  pivot_wider(names_from = TNT, values_from = Abund_Cov)


# PWY-621: sucrose degradation III (sucrose invertase)
p5_PWY621 = path_select_shared_input %>% 
  filter(Pathway == "PWY-621: sucrose degradation III (sucrose invertase)") %>% 
  ggplot(aes(x = Before, y = Ongoing)) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.7,
              linetype = "dashed", size = 0.4) +
  geom_point(aes(fill = TRG_1), 
             size = 3, shape = 21, color = "gray40") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  theme_pubr() +
  labs(x = "Before", y = "Ongoing", color = "TRG_1",
       title = "PWY-621: sucrose degradation III (sucrose invertase)") ; p5_PWY621

# PWY-6892: thiazole component of thiamine diphosphate biosynthesis I
p5_PWY6892 = path_select_shared_input %>% 
  filter(Pathway == "PWY-6892: thiazole component of thiamine diphosphate biosynthesis I") %>% 
  ggplot(aes(x = Before, y = Ongoing)) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.7,
              linetype = "dashed", size = 0.4) +
  geom_point(aes(fill = TRG_1), 
             size = 3, shape = 21, color = "gray40") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  theme_pubr() +
  labs(x = "Before", y = "Ongoing", color = "TRG_1",
       title = "PWY-6892: thiazole component of thiamine diphosphate biosynthesis I") ; p5_PWY6892

# PANTO-PWY: phosphopantothenate biosynthesis I
p5_PANTO_PWY = path_select_shared_input %>% 
  filter(Pathway == "PANTO-PWY: phosphopantothenate biosynthesis I") %>% 
  ggplot(aes(x = Before, y = Ongoing)) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.7,
              linetype = "dashed", size = 0.4) +
  geom_point(aes(fill = TRG_1), 
             size = 3, shape = 21, color = "gray40") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  theme_pubr() +
  labs(x = "Before", y = "Ongoing", color = "TRG_1",
       title = "PANTO-PWY: phosphopantothenate biosynthesis I") ; p5_PANTO_PWY

p5_selected_shared_path = ggarrange(p5_PWY621,
                                    p5_PWY6892,
                                    p5_PANTO_PWY,
                                    nrow = 1, ncol = 3) ; p5_selected_shared_path

# ggsave("figure/06-3-2_selected shared pathways_scatter.svg",
#        plot = p5_selected_shared_path, width = 15, height = 5)



###############  8. Coenzyme A  ###############

# Ongoing

# Any pathway contain a term "methylerythritol"
path_coenzymeA_o = path_sorted_o %>% 
  filter(grepl("coenzyme A", Pathway, ignore.case = TRUE)) %>% 
  pull(Pathway) %>% 
  unique()

# Ordering by abundance
path_coenzymeA_o = path_sorted_o %>% 
  filter(Pathway %in% path_coenzymeA_o) %>% 
  group_by(Pathway) %>% 
  summarise(Abund_sum = sum(Abundance_0.5)) %>% 
  arrange(-Abund_sum) %>% 
  pull(Pathway) ; path_coenzymeA_o
# [1] "COA-PWY-1: superpathway of coenzyme A biosynthesis III (mammals)"  
# [2] "COA-PWY: coenzyme A biosynthesis I (prokaryotic)"                  
# [3] "PWY-7851: coenzyme A biosynthesis II (eukaryotic)"                 
# [4] "PANTOSYN-PWY: superpathway of coenzyme A biosynthesis I (bacteria)"


path_coenzymeA_wide_o = path_sorted_o %>% 
  filter(Pathway %in% path_coenzymeA_o) %>% 
  pivot_wider(id_cols = c("SampleID", "TRG_1", "TRG_score"),
              names_from = Pathway,
              values_from = Abundance_0.5)


# Spearman, p-value
path_coenzymeA_stat_o = cor.test(path_coenzymeA_wide_o$`COA-PWY: coenzyme A biosynthesis I (prokaryotic)`,
                                 path_coenzymeA_wide_o$`PANTOSYN-PWY: superpathway of coenzyme A biosynthesis I (bacteria)`,
                                 method = "spearman") ; path_coenzymeA_stat_o

path_coenzymeA_spearman_o = path_coenzymeA_stat_o$estimate ; path_coenzymeA_spearman_o # 0.08388523 
path_coenzymeA_p_o = path_coenzymeA_stat_o$p.value ; path_coenzymeA_p_o # 0.7574225


# 시각화: 이를 바탕으로 categorization 기준 설정
p5_coenzymeA_TRG_1_categorization = path_coenzymeA_wide_o %>% 
  ggplot(aes (`COA-PWY: coenzyme A biosynthesis I (prokaryotic)`,
              `PANTOSYN-PWY: superpathway of coenzyme A biosynthesis I (bacteria)`)) +
  geom_point(aes(color = TRG_1)) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  annotate("text", x = min(path_coenzymeA_wide_o$`COA-PWY: coenzyme A biosynthesis I (prokaryotic)`), 
           y = max(path_coenzymeA_wide_o$`PANTOSYN-PWY: superpathway of coenzyme A biosynthesis I (bacteria)`), 
           label = paste("Spearman r =", round(path_coenzymeA_spearman_o, 2), 
                         "\np-value =", round(path_coenzymeA_p_o, 4)), 
           hjust = 0, vjust = 1, size = 3, color = "black") +
  theme_classic() +
  theme(legend.position = "none") ; p5_coenzymeA_TRG_1_categorization

# ggsave("figure/05-3-17_coenzymeA_TRG_1-scatter_ongoing.svg",
#        plot = p5_coenzymeA_TRG_1_categorization, width = 5, height = 5)


path_coenzymeA_wide_o = path_coenzymeA_wide_o %>% 
  mutate(
    
    coenzymeA_biosynthesis = ifelse(
      `COA-PWY: coenzyme A biosynthesis I (prokaryotic)` > 20000, 1, 0
    ),
    
    coenzymeA_degradation = ifelse(
      `PANTOSYN-PWY: superpathway of coenzyme A biosynthesis I (bacteria)` > 20000, 1, 0
    ),
    
    coenzymeA_point = coenzymeA_biosynthesis + coenzymeA_degradation,
    
    coenzymeA_bin = ifelse(coenzymeA_point == 2, "Y", "N")
    
  )

table(path_coenzymeA_wide_o$coenzymeA_bin,
      path_coenzymeA_wide_o$TRG_score)
#   0 1 2 3
# N 6 2 2 2
# Y 1 0 3 0

path_coenzymeA_wilcox_o = wilcox.test(TRG_score ~ coenzymeA_bin, data = path_coenzymeA_wide_o)
path_coenzymeA_wilcox_o$p.value # 0.4784073


library(ggbeeswarm)
library(waffle)
library(hrbrthemes)

# Make a sucrose dataframe
coenzymeA_df_o = as.data.frame(table(path_coenzymeA_wide_o$coenzymeA_bin,
                                     path_coenzymeA_wide_o$TRG_score))
colnames(coenzymeA_df_o) = c("coenzymeA_bin", "TRG_score", "Frequency")


# TRG_score를 숫자형으로 변환
is.numeric(coenzymeA_df_o$TRG_score) # FALSE 
coenzymeA_df_o$TRG_score = as.numeric(as.character(coenzymeA_df_o$TRG_score))
score_colors = c("0" = "#FFEDA0", 
                 "1" = "#FEB24C",
                 "2" = "#FC4E2A",
                 "3" = "#BD0026")


# waffle plot
p5_coenzymeA_TRG_1_waffle = coenzymeA_df_o %>% 
  group_by(coenzymeA_bin) %>% 
  mutate(label = paste(TRG_score, coenzymeA_bin, sep = "_")) %>% 
  ggplot(aes(fill = factor(TRG_score), values = Frequency)) +
  geom_waffle(n_rows = 5, size = 2, colour = "white", flip = T) +
  scale_fill_manual(values = score_colors,
                    name = "TRG Score") +
  facet_wrap(~coenzymeA_bin, nrow = 1, scales = "fixed") +
  coord_equal() +
  annotate("text", x = 0.5, y = 4.0,
           label = paste("p-value = ", round(path_coenzymeA_wilcox_o$p.value, 4)),
           size = 3, color = "black", hjust = 0) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = rel(1.2)),
        strip.text= element_text(size = rel(2), face = "bold")) ; p5_coenzymeA_TRG_1_waffle

# ggsave("figure/05-3-17_coenzyme A-TRG_1-waffle_ongoing.svg",
#        plot = p5_coenzymeA_TRG_1_waffle, width = 8, height = 4.5)



###############  9. Phosphopantothenate  ###############

# Ongoing

# Any pathway contain a term "panto"
path_panto_o = path_sorted_o %>% 
  filter(grepl("panto", Pathway, ignore.case = TRUE)) %>% 
  pull(Pathway) %>% 
  unique()

# Ordering by abundance
path_panto_o = path_sorted_o %>% 
  filter(Pathway %in% path_panto_o) %>% 
  group_by(Pathway) %>% 
  summarise(Abund_sum = sum(Abundance_0.5)) %>% 
  arrange(-Abund_sum) %>% 
  pull(Pathway) ; path_panto_o
# [1] "PANTO-PWY: phosphopantothenate biosynthesis I"                     
# [2] "PANTOSYN-PWY: superpathway of coenzyme A biosynthesis I (bacteria)"


path_panto_wide_o = path_sorted_o %>% 
  filter(Pathway %in% path_panto_o) %>% 
  pivot_wider(id_cols = c("SampleID", "TRG_1", "TRG_score"),
              names_from = Pathway,
              values_from = Abundance_0.5)


# Spearman, p-value
path_panto_stat_o = cor.test(path_panto_wide_o$`PANTO-PWY: phosphopantothenate biosynthesis I`,
                             path_panto_wide_o$`PANTOSYN-PWY: superpathway of coenzyme A biosynthesis I (bacteria)`,
                             method = "spearman") ; path_panto_stat_o

path_panto_spearman_o = path_panto_stat_o$estimate ; path_panto_spearman_o # 0.6313467  
path_panto_p_o = path_panto_stat_o$p.value ; path_panto_p_o # 0.008716314


# 시각화: 이를 바탕으로 categorization 기준 설정
p5_panto_TRG_1_categorization = path_panto_wide_o %>% 
  ggplot(aes (`PANTO-PWY: phosphopantothenate biosynthesis I`,
              `PANTOSYN-PWY: superpathway of coenzyme A biosynthesis I (bacteria)`)) +
  geom_point(aes(color = TRG_1)) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  annotate("text", x = min(path_panto_wide_o$`PANTO-PWY: phosphopantothenate biosynthesis I`), 
           y = max(path_panto_wide_o$`PANTOSYN-PWY: superpathway of coenzyme A biosynthesis I (bacteria)`), 
           label = paste("Spearman r =", round(path_panto_spearman_o, 2), 
                         "\np-value =", round(path_panto_p_o, 4)), 
           hjust = 0, vjust = 1, size = 3, color = "black") +
  theme_classic() +
  theme(legend.position = "none") ; p5_panto_TRG_1_categorization

# ggsave("figure/05-3-17_panto_TRG_1-scatter_ongoing.svg",
#        plot = p5_panto_TRG_1_categorization, width = 5, height = 5)


path_panto_wide_o = path_panto_wide_o %>% 
  mutate(
    
    panto_biosynthesis = ifelse(
      `PANTO-PWY: phosphopantothenate biosynthesis I` > 18000, 1, 0
    ),
    
    panto_degradation = ifelse(
      `PANTOSYN-PWY: superpathway of coenzyme A biosynthesis I (bacteria)` > 17000, 1, 0
    ),
    
    panto_point = panto_biosynthesis + panto_degradation,
    
    panto_bin = ifelse(panto_point == 2, "Y", "N")
    
  )

table(path_panto_wide_o$panto_bin,
      path_panto_wide_o$TRG_score)
#   0 1 2 3
# N 6 2 1 2
# Y 1 0 4 0

path_panto_wilcox_o = wilcox.test(TRG_score ~ panto_bin, data = path_panto_wide_o)
path_panto_wilcox_o$p.value # 0.2785261


library(ggbeeswarm)
library(waffle)
library(hrbrthemes)

# Make a sucrose dataframe
panto_df_o = as.data.frame(table(path_panto_wide_o$panto_bin,
                                 path_panto_wide_o$TRG_score))
colnames(panto_df_o) = c("panto_bin", "TRG_score", "Frequency")


# TRG_score를 숫자형으로 변환
is.numeric(panto_df_o$TRG_score) # FALSE 
panto_df_o$TRG_score = as.numeric(as.character(panto_df_o$TRG_score))
score_colors = c("0" = "#FFEDA0", 
                 "1" = "#FEB24C",
                 "2" = "#FC4E2A",
                 "3" = "#BD0026")


# waffle plot
p5_panto_TRG_1_waffle = panto_df_o %>% 
  group_by(panto_bin) %>% 
  mutate(label = paste(TRG_score, panto_bin, sep = "_")) %>% 
  ggplot(aes(fill = factor(TRG_score), values = Frequency)) +
  geom_waffle(n_rows = 5, size = 2, colour = "white", flip = T) +
  scale_fill_manual(values = score_colors,
                    name = "TRG Score") +
  facet_wrap(~panto_bin, nrow = 1, scales = "fixed") +
  coord_equal() +
  annotate("text", x = 0.5, y = 4.0,
           label = paste("p-value = ", round(path_panto_wilcox_o$p.value, 4)),
           size = 3, color = "black", hjust = 0) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = rel(1.2)),
        strip.text= element_text(size = rel(2), face = "bold")) ; p5_panto_TRG_1_waffle

# ggsave("figure/05-3-17_panto-TRG_1-waffle_ongoing.svg",
#        plot = p5_panto_TRG_1_waffle, width = 8, height = 4.5)



###############  Pathway by bacteria  ###############

# HISTSYN-PWY: L-histidine biosynthesis
# Among 7 significant pathways (P-value = 0.0968)
# highest mean abundance: 19303.062 
# highest mean coverage: 0.96 (vs. 0.69 of PWY-6892, 0.35 of RHAMCAT-PWY, 0.26 of PWY-5384, 
#                                  0.82 of PWY-6270, 0.15 of PWY-621, 0.82 of PWY-7560)
# Log2 FC: 0.260


# Filtering 
path_histidine_filtered = path_abund %>%
  filter(str_detect(Pathway, "HISTSYN-PWY: L-histidine biosynthesis")) %>% 
  filter(Pathway != "HISTSYN-PWY: L-histidine biosynthesis")

path_histidine_filtered_cover = path_cover %>%
  filter(str_detect(Pathway, "HISTSYN-PWY: L-histidine biosynthesis"))%>% 
  filter(Pathway != "HISTSYN-PWY: L-histidine biosynthesis")

path_histidine_filtered_0.5 = path_histidine_filtered

# coverage가 0.5 미만인 경우는 abundance를 0으로 간주
path_histidine_filtered_0.5[path_histidine_filtered_cover < 0.5] = 0


path_histidine_filtered_long = path_histidine_filtered_0.5 %>% 
  pivot_longer(
    
    cols = -Pathway, 
    names_to = "SampleID", 
    values_to = "Abundance_0.5"
    
  ) %>% 
  filter(Abundance_0.5 > 0) %>%
  merge(., mb2 %>% select(SampleID, 
                          Age, Sex,
                          TRG_score, TRG_1, TRG_2, TRG_3, 
                          Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                          Pre_Op_Nstage, Pre_Op_Nstage_bin),
        by = "SampleID") %>% 
  separate(Pathway, into = c("Pathway_ID", "Taxonomy"), sep = "\\|", fill = "right") %>%
  
  # 2. Taxonomy 분할: g__와 s__가 있는 경우에만 genus/species 추출
  mutate(
    Genus = ifelse(str_detect(Taxonomy, "g__"),
                   str_extract(Taxonomy, "g__([^\\.]+)") %>% str_remove("g__"),
                   "unclassified"),
    Species = ifelse(str_detect(Taxonomy, "s__"),
                     str_extract(Taxonomy, "s__.+") %>% str_remove("s__"),
                     "unclassified")
  ) ; head(path_histidine_filtered_long)


# By genus & TRG_1
path_histidine_genus = path_histidine_filtered_long %>% 
  filter(Genus != "unclassified") %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance_0.5), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice(1:6) %>%
  pull(Genus)

path_histidine_barplot = path_histidine_filtered_long %>% 
  mutate(Genus_grouped = case_when(
    Genus %in% path_histidine_genus ~ Genus,
    Genus == "unclassified" ~ "unclassified",
    TRUE ~ "Others"
  ))


path_histidine_barplot_grouped = path_histidine_barplot %>%
  group_by(Genus_grouped, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15))


# 4. stacked barplot 그리기
genus_colors = c(
  "Prevotella"              = "#78C38F",  # 밝은 민트-그린 (긍정적 강조)
  "Eubacterium"             = "#66B2E4",  # 부드러운 파스텔 블루
  "Subdoligranulum"         = "#E68A71",  # 페일 청록 (pastel turquoise)
  "Dialister"               = "#E06646",  # 따뜻한 오렌지-살구톤
  "Firmicutes_unclassified" = "#B34040",  # 진한 적갈색
  "Parabacteroides"         = "#CC5500",  # 고동색-주황 (burnt orange)
  "Others"                  = "#BFBFBF",   # 아주 옅은 회색
  "unclassified"            = "#DDDDDD"  # 연한 회색
)


p4_pathway_histidine_genus = path_histidine_barplot_grouped %>% 
  mutate(Genus_grouped = factor(Genus_grouped, 
                                levels = rev(c("Prevotella",
                                               "Eubacterium",
                                               "Subdoligranulum",
                                               "Dialister",
                                               "Firmicutes_unclassified",
                                               "Parabacteroides",
                                               "Others",
                                               "unclassified")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional abundance (HISTSYN-PWY)",
       fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors) ; p4_pathway_histidine_genus

# ggsave("figure/04-3-18_histidine-genus.svg",
#        plot = p4_pathway_histidine_genus, width = 5, height = 5.5)


genus_colors = c(
  "Prevotella"              = "#78C38F", 
  "Eubacterium"             = "#66B2E4",  
  "Subdoligranulum"         = "#E68A71",  
  "Dialister"               = "#E06646",  
  "Firmicutes_unclassified" = "#B34040",  
  "Parabacteroides"         = "#CC5500",  
  "Others"                  = "#BFBFBF"  
)


p4_pathway_histidine_genus_woUnclassified = path_histidine_barplot_grouped %>% 
  filter(Genus_grouped != "unclassified") %>% 
  mutate(Genus_grouped = factor(Genus_grouped, 
                                levels = rev(c("Prevotella",
                                               "Eubacterium",
                                               "Subdoligranulum",
                                               "Dialister",
                                               "Firmicutes_unclassified",
                                               "Parabacteroides",
                                               "Others")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional abundance (HISTSYN-PWY)",
       fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors) ; p4_pathway_histidine_genus_woUnclassified

# ggsave("figure/04-3-18_histidine-genus-woUnclassified.svg",
#        plot = p4_pathway_histidine_genus_woUnclassified, width = 5, height = 5.5)



# PWY-6892: thiazole component of thiamine diphosphate biosynthesis I
# Filtering 
path_thiamine_filtered = path_abund %>%
  filter(str_detect(Pathway, "PWY-6892: thiazole component of thiamine diphosphate biosynthesis I")) %>% 
  filter(Pathway != "PWY-6892: thiazole component of thiamine diphosphate biosynthesis I")

path_thiamine_filtered_cover = path_cover %>%
  filter(str_detect(Pathway, "PWY-6892: thiazole component of thiamine diphosphate biosynthesis I"))%>% 
  filter(Pathway != "PWY-6892: thiazole component of thiamine diphosphate biosynthesis I")

path_thiamine_filtered_0.5 = path_thiamine_filtered

# coverage가 0.5 미만인 경우는 abundance를 0으로 간주
path_thiamine_filtered_0.5[path_thiamine_filtered_cover < 0.5] = 0


path_thiamine_filtered_long = path_thiamine_filtered_0.5 %>% 
  pivot_longer(
    
    cols = -Pathway, 
    names_to = "SampleID", 
    values_to = "Abundance_0.5"
    
  ) %>% 
  filter(Abundance_0.5 > 0) %>%
  merge(., mb2 %>% select(SampleID, 
                          Age, Sex,
                          TRG_score, TRG_1, TRG_2, TRG_3, 
                          Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                          Pre_Op_Nstage, Pre_Op_Nstage_bin),
        by = "SampleID") %>% 
  separate(Pathway, into = c("Pathway_ID", "Taxonomy"), sep = "\\|", fill = "right") %>%
  
  # 2. Taxonomy 분할: g__와 s__가 있는 경우에만 genus/species 추출
  mutate(
    Genus = ifelse(str_detect(Taxonomy, "g__"),
                   str_extract(Taxonomy, "g__([^\\.]+)") %>% str_remove("g__"),
                   "unclassified"),
    Species = ifelse(str_detect(Taxonomy, "s__"),
                     str_extract(Taxonomy, "s__.+") %>% str_remove("s__"),
                     "unclassified")
  ) ; head(path_thiamine_filtered_long)


# By genus & TRG_1
path_thiamine_genus = path_thiamine_filtered_long %>% 
  filter(Genus != "unclassified") %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance_0.5), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice(1:6) %>%
  pull(Genus)

path_thiamine_barplot = path_thiamine_filtered_long %>% 
  mutate(Genus_grouped = case_when(
    Genus %in% path_thiamine_genus ~ Genus,
    Genus == "unclassified" ~ "unclassified",
    TRUE ~ "Others"
  ))

path_thiamine_barplot %>% 
  filter(Genus == "Anaerostipes") %>% 
  select(Species) # Anaerostipes_hadrus

path_thiamine_barplot %>% 
  filter(Genus == "Klebsiella") %>% 
  select(Species) # Klebsiella_pneumoniae, Klebsiella_oxytoca


path_thiamine_barplot_grouped = path_thiamine_barplot %>%
  group_by(Genus_grouped, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15))


# 4. stacked barplot 그리기
genus_colors = c(
  "Blautia"             = "#9A4E97",
  "Anaerostipes"        = "#A3C37B",  
  "Coprococcus"         = "#E3A15C",  
  "Klebsiella"          = "#D84B6A",  
  "Escherichia"         = "#4C91E0",  
  "Eubacterium"         = "#66B2E4",  
  "Others"              = "#BFBFBF", 
  "unclassified"        = "#DDDDDD" 
)


p4_pathway_thiamine_genus = path_thiamine_barplot_grouped %>% 
  mutate(Genus_grouped = factor(Genus_grouped, 
                                levels = rev(c("Blautia",
                                               "Anaerostipes",
                                               "Coprococcus",
                                               "Klebsiella",
                                               "Escherichia",
                                               "Eubacterium",
                                               "Others",
                                               "unclassified")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional abundance (PWY-6892)",
       fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors) ; p4_pathway_thiamine_genus

# ggsave("figure/04-3-19_thiamine-genus.svg",
#        plot = p4_pathway_thiamine_genus, width = 5, height = 5.5)


genus_colors = c(
  "Blautia"             = "#9A4E97",
  "Anaerostipes"        = "#A3C37B",  
  "Coprococcus"         = "#E3A15C",  
  "Klebsiella"          = "#D84B6A",  
  "Escherichia"         = "#4C91E0",  
  "Eubacterium"         = "#66B2E4", 
  "Others"              = "#BFBFBF"
)


p4_pathway_thiamine_genus_woUnclassified = path_thiamine_barplot_grouped %>% 
  filter(Genus_grouped != "unclassified") %>% 
  mutate(Genus_grouped = factor(Genus_grouped, 
                                levels = rev(c("Blautia",
                                               "Anaerostipes",
                                               "Coprococcus",
                                               "Klebsiella",
                                               "Escherichia",
                                               "Eubacterium",
                                               "Others")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional abundance (PWY-6892)",
       fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors) ; p4_pathway_thiamine_genus_woUnclassified

# ggsave("figure/04-3-19_thiamine-genus-woUnclassified.svg",
#        plot = p4_pathway_thiamine_genus_woUnclassified, width = 5, height = 5.5)



# PWY-621: sucrose degradation III (sucrose invertase)
# Filtering
path_sucrose_filtered = path_abund %>%
  filter(str_detect(Pathway, "PWY-621: sucrose degradation III \\(sucrose invertase\\)")) %>% 
  filter(Pathway != "PWY-621: sucrose degradation III (sucrose invertase)")

path_sucrose_filtered_cover = path_cover %>% 
  filter(str_detect(Pathway, "PWY-621: sucrose degradation III \\(sucrose invertase\\)")) %>% 
  filter(Pathway != "PWY-621: sucrose degradation III (sucrose invertase)")

path_sucrose_filtered_0.5 = path_sucrose_filtered

# coverage가 0.5 미만인 경우는 abudnance를 0으로 간주
path_sucrose_filtered_0.5[path_sucrose_filtered_cover < 0.5] = 0


path_sucrose_filtered_long = path_sucrose_filtered_0.5 %>% 
  pivot_longer(
    
    cols = -Pathway,
    names_to = "SampleID",
    values_to = "Abundance_0.5"
    
  ) %>% 
  filter(Abundance_0.5 > 0) %>% 
  merge(., mb2 %>% select(SampleID, 
                          Age, Sex,
                          TRG_score, TRG_1, TRG_2, TRG_3, 
                          Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                          Pre_Op_Nstage, Pre_Op_Nstage_bin),
        by = "SampleID") %>% 
  separate(Pathway, into = c("Pathway_ID", "Taxonomy"), sep = "\\|", fill = "right") %>% 
  
  # 2. Taxonomy 분할: g__와 s__가 있는 경우에만 genus/species 추출
  mutate(
    Genus = ifelse(str_detect(Taxonomy, "g__"),
                   str_extract(Taxonomy, "g__([^\\.]+)") %>% str_remove("g__"),
                   "unclassified"),
    Species = ifelse(str_detect(Taxonomy, "s__"),
                     str_extract(Taxonomy, "s__.+") %>% str_remove("s__"),
                     "unclassified")
  ) ; head(path_sucrose_filtered_long)


# By genus & TRG_1
path_sucrose_genus = path_sucrose_filtered_long %>% 
  filter(Genus != "unclassified") %>% 
  group_by(Genus) %>% 
  summarise(Total = sum(Abundance_0.5), .groups = "drop") %>% 
  arrange(desc(Total)) %>% 
  slice(1:6) %>% 
  pull(Genus)

path_sucrose_barplot = path_sucrose_filtered_long %>% 
  mutate(Genus_grouped = case_when(
    Genus %in% path_sucrose_genus ~ Genus,
    Genus == "unclassified" ~ "unclassified",
    TRUE ~ "Others"
  ))

path_sucrose_barplot %>% 
  filter(Genus == "Anaerostipes") %>% 
  select(Species) # Anaerostipes_hadrus

path_sucrose_barplot %>% 
  filter(Genus == "Klebsiella") %>% 
  select(Species)


path_sucrose_barplot_grouped = path_sucrose_barplot %>% 
  group_by(Genus_grouped, TRG_1) %>% 
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15))


# 4. stacked barplot 그리기
genus_colors = c(
  "Anaerostipes"        = "#A3C37B", 
  "Streptococcus"       = "#D96C6C",
  "Escherichia"         = "#4C91E0",  
  "Bacteroides"         = "#F1A7C1",
  "Klebsiella"          = "#D84B6A",  
  "Bifidobacterium"     = "#F0C16B",
  "Others"              = "#BFBFBF", 
  "unclassified"        = "#DDDDDD"
)


p4_pathway_sucrose_genus = path_sucrose_barplot_grouped %>% 
  mutate(Genus_grouped = factor(Genus_grouped, 
                                levels = rev(c("Anaerostipes",
                                               "Streptococcus",
                                               "Escherichia",
                                               "Bacteroides",
                                               "Klebsiella",
                                               "Bifidobacterium",
                                               "Others",
                                               "unclassified")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional abundance (PWY-621)",
       fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors) ; p4_pathway_sucrose_genus

# ggsave("figure/04-3-19_sucrose-genus.svg",
#        plot = p4_pathway_sucrose_genus, width = 5, height = 5.5)


genus_colors = c(
  "Anaerostipes"        = "#A3C37B", 
  "Streptococcus"       = "#D96C6C",
  "Escherichia"         = "#4C91E0",  
  "Bacteroides"         = "#F1A7C1",
  "Klebsiella"          = "#D84B6A",  
  "Bifidobacterium"     = "#F0C16B",
  "Others"              = "#BFBFBF"
)


p4_pathway_sucrose_genus_woUnclassified = path_sucrose_barplot_grouped %>% 
  filter(Genus_grouped != "unclassified") %>% 
  mutate(Genus_grouped = factor(Genus_grouped, 
                                levels = rev(c("Anaerostipes",
                                               "Streptococcus",
                                               "Escherichia",
                                               "Bacteroides",
                                               "Klebsiella",
                                               "Bifidobacterium",
                                               "Others")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional abundance (PWY-621)",
       fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors) ; p4_pathway_sucrose_genus_woUnclassified

# ggsave("figure/04-3-19_sucrose-genus-woUnclassified.svg",
#        plot = p4_pathway_sucrose_genus_woUnclassified, width = 5, height = 5.5)



# PWY-5384: sucrose degradation IV (sucrose phosphorylase)
# Filtering
path_sucrose4_filtered = path_abund %>%
  filter(str_detect(Pathway, "PWY-5384: sucrose degradation IV \\(sucrose phosphorylase\\)")) %>% 
  filter(Pathway != "PWY-5384: sucrose degradation IV (sucrose phosphorylase)")

path_sucrose4_filtered_cover = path_cover %>%
  filter(str_detect(Pathway, "PWY-5384: sucrose degradation IV \\(sucrose phosphorylase\\)"))%>% 
  filter(Pathway != "PWY-5384: sucrose degradation IV (sucrose phosphorylase)")

path_sucrose4_filtered_0.5 = path_sucrose4_filtered

# coverage가 0.5 미만인 경우는 abudnance를 0으로 간주
path_sucrose4_filtered_0.5[path_sucrose4_filtered_cover < 0.5] = 0


path_sucrose4_filtered_long = path_sucrose4_filtered_0.5 %>% 
  pivot_longer(
    
    cols = -Pathway,
    names_to = "SampleID",
    values_to = "Abundance_0.5"
    
  ) %>% 
  filter(Abundance_0.5 > 0) %>% 
  merge(., mb2 %>% select(SampleID, 
                          Age, Sex,
                          TRG_score, TRG_1, TRG_2, TRG_3, 
                          Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                          Pre_Op_Nstage, Pre_Op_Nstage_bin),
        by = "SampleID") %>% 
  separate(Pathway, into = c("Pathway_ID", "Taxonomy"), sep = "\\|", fill = "right") %>% 
  
  # 2. Taxonomy 분할: g__와 s__가 있는 경우에만 genus/species 추출
  mutate(
    Genus = ifelse(str_detect(Taxonomy, "g__"),
                   str_extract(Taxonomy, "g__([^\\.]+)") %>% str_remove("g__"),
                   "unclassified"),
    Species = ifelse(str_detect(Taxonomy, "s__"),
                     str_extract(Taxonomy, "s__.+") %>% str_remove("s__"),
                     "unclassified")
  ) ; head(path_sucrose4_filtered_long)


# By genus & TRG_1
path_sucrose4_genus = path_sucrose4_filtered_long %>% 
  filter(Genus != "unclassified") %>% 
  group_by(Genus) %>% 
  summarise(Total = sum(Abundance_0.5), .groups = "drop") %>% 
  arrange(desc(Total)) %>% 
  slice(1:6) %>% 
  pull(Genus)

path_sucrose4_barplot = path_sucrose4_filtered_long %>% 
  mutate(Genus_grouped = case_when(
    Genus %in% path_sucrose4_genus ~ Genus,
    Genus == "unclassified" ~ "unclassified",
    TRUE ~ "Others"
  ))

path_sucrose4_barplot %>% 
  filter(Genus == "Anaerostipes") %>% 
  select(Species)


path_sucrose4_barplot_grouped = path_sucrose4_barplot %>% 
  group_by(Genus_grouped, TRG_1) %>% 
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15)) 


# 4. stacked barplot 그리기
genus_colors = c(
  "Bifidobacterium"     = "#F0C16B",
  "Anaerostipes"        = "#A3C37B", 
  "Escherichia"         = "#4C91E0",  
  "Streptococcus"       = "#D96C6C",
  "Klebsiella"          = "#D84B6A",  
  "Roseburia"           = "#6F8A3F",
  "Others"              = "#BFBFBF", 
  "unclassified"        = "#DDDDDD"
)


p4_pathway_sucrose4_genus = path_sucrose4_barplot_grouped %>% 
  mutate(Genus_grouped = factor(Genus_grouped, 
                                levels = rev(c("Bifidobacterium",
                                               "Anaerostipes",
                                               "Escherichia",
                                               "Streptococcus",
                                               "Klebsiella",
                                               "Roseburia",
                                               "Others",
                                               "unclassified")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional abundance (PWY-5384)",
       fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors) ; p4_pathway_sucrose4_genus

# ggsave("figure/04-3-19_sucrose4-genus.svg",
#        plot = p4_pathway_sucrose4_genus, width = 5, height = 5.5)


genus_colors = c(
  "Bifidobacterium"     = "#F0C16B",
  "Anaerostipes"        = "#A3C37B", 
  "Escherichia"         = "#4C91E0",  
  "Streptococcus"       = "#D96C6C",
  "Klebsiella"          = "#D84B6A",  
  "Roseburia"           = "#6F8A3F",
  "Others"              = "#BFBFBF"
)


p4_pathway_sucrose4_genus_woUnclassified = path_sucrose4_barplot_grouped %>% 
  filter(Genus_grouped != "unclassified") %>% 
  mutate(Genus_grouped = factor(Genus_grouped, 
                                levels = rev(c("Bifidobacterium",
                                               "Anaerostipes",
                                               "Escherichia",
                                               "Streptococcus",
                                               "Klebsiella",
                                               "Roseburia",
                                               "Others")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional abundance (PWY-5384)",
       fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors) ; p4_pathway_sucrose4_genus_woUnclassified

# ggsave("figure/04-3-19_sucrose4-genus-woUnclassified.svg",
#        plot = p4_pathway_sucrose4_genus_woUnclassified, width = 5, height = 5.5)



# PWY-7238: sucrose biosynthesis II
# Filtering
path_sucrose_syn_filtered = path_abund %>% 
  filter(str_detect(Pathway, "PWY-7238: sucrose biosynthesis II")) %>% 
  filter(Pathway != "PWY-7238: sucrose biosynthesis II")

path_sucrose_syn_filtered_cover = path_cover %>%
  filter(str_detect(Pathway, "PWY-7238: sucrose biosynthesis II"))%>% 
  filter(Pathway != "PWY-7238: sucrose biosynthesis II")

path_sucrose_syn_filtered_0.5 = path_sucrose_syn_filtered

# coverage가 0.5 미만인 경우는 abudnance를 0으로 간주
path_sucrose_syn_filtered_0.5[path_sucrose_syn_filtered_cover < 0.5] = 0


path_sucrose_syn_filtered_long = path_sucrose_syn_filtered_0.5 %>% 
  pivot_longer(
    
    cols = -Pathway,
    names_to = "SampleID",
    values_to = "Abundance_0.5"
    
  ) %>% 
  filter(Abundance_0.5 > 0) %>% 
  merge(., mb2 %>% select(SampleID, 
                          Age, Sex,
                          TRG_score, TRG_1, TRG_2, TRG_3, 
                          Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                          Pre_Op_Nstage, Pre_Op_Nstage_bin),
        by = "SampleID") %>% 
  separate(Pathway, into = c("Pathway_ID", "Taxonomy"), sep = "\\|", fill = "right") %>% 
  
  # 2. Taxonomy 분할: g__와 s__가 있는 경우에만 genus/species 추출
  mutate(
    Genus = ifelse(str_detect(Taxonomy, "g__"),
                   str_extract(Taxonomy, "g__([^\\.]+)") %>% str_remove("g__"),
                   "unclassified"),
    Species = ifelse(str_detect(Taxonomy, "s__"),
                     str_extract(Taxonomy, "s__.+") %>% str_remove("s__"),
                     "unclassified")
  ) ; head(path_sucrose_syn_filtered_long)


# By genus & TRG_1
path_sucrose_syn_genus = path_sucrose_syn_filtered_long %>% 
  filter(Genus != "unclassified") %>% 
  group_by(Genus) %>% 
  summarise(Total = sum(Abundance_0.5), .groups = "drop") %>% 
  arrange(desc(Total)) %>% 
  slice(1:6) %>% 
  pull(Genus)

path_sucrose_syn_barplot = path_sucrose_syn_filtered_long %>% 
  mutate(Genus_grouped = case_when(
    Genus %in% path_sucrose_syn_genus ~ Genus,
    Genus == "unclassified" ~ "unclassified",
    TRUE ~ "Others"
  ))


path_sucrose_syn_barplot_grouped = path_sucrose_syn_barplot %>% 
  group_by(Genus_grouped, TRG_1) %>% 
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15)) 


# 4. stacked barplot 그리기
genus_colors = c(
  "Blautia"             = "#9A4E97",
  "Faecalibacterium"    = "#8BC34A",
  "Bifidobacterium"     = "#F0C16B",
  "Collinsella"         = "#1C85C4", 
  "Eubacterium"         = "#66B2E4",  
  "Dorea"               = "#FF9800",
  "Others"              = "#BFBFBF", 
  "unclassified"        = "#DDDDDD"
)


p4_pathway_sucrose_syn_genus = path_sucrose_syn_barplot_grouped %>% 
  mutate(Genus_grouped = factor(Genus_grouped, 
                                levels = rev(c("Blautia",
                                               "Faecalibacterium",
                                               "Bifidobacterium",
                                               "Collinsella",
                                               "Eubacterium",
                                               "Dorea",
                                               "Others",
                                               "unclassified")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional abundance (PWY-7238)",
       fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors) ; p4_pathway_sucrose_syn_genus

# ggsave("figure/04-3-19_sucrose_syn-genus.svg",
#        plot = p4_pathway_sucrose_syn_genus, width = 5, height = 5.5)


genus_colors = c(
  "Blautia"             = "#9A4E97",
  "Faecalibacterium"    = "#8BC34A",
  "Bifidobacterium"     = "#F0C16B",
  "Collinsella"         = "#1C85C4", 
  "Eubacterium"         = "#66B2E4",  
  "Dorea"               = "#FF9800",
  "Others"              = "#BFBFBF"
)


p4_pathway_sucrose_syn_genus_woUnclassified = path_sucrose_syn_barplot_grouped %>% 
  filter(Genus_grouped != "unclassified") %>% 
  mutate(Genus_grouped = factor(Genus_grouped, 
                                levels = rev(c("Blautia",
                                               "Faecalibacterium",
                                               "Bifidobacterium",
                                               "Collinsella",
                                               "Eubacterium",
                                               "Dorea",
                                               "Others")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional abundance (PWY-7238)",
       fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors) ; p4_pathway_sucrose_syn_genus_woUnclassified

# ggsave("figure/04-3-19_sucrose_syn-genus-woUnclassified.svg",
#        plot = p4_pathway_sucrose_syn_genus_woUnclassified, width = 5, height = 5.5)



# RHAMCAT-PWY: L-rhamnose degradation I
# Filtering 
path_rhamnose_filtered = path_abund %>%
  filter(str_detect(Pathway, "RHAMCAT-PWY: L-rhamnose degradation I")) %>% 
  filter(Pathway != "RHAMCAT-PWY: L-rhamnose degradation I")

path_rhamnose_filtered_cover = path_cover %>%
  filter(str_detect(Pathway, "RHAMCAT-PWY: L-rhamnose degradation I"))%>% 
  filter(Pathway != "RHAMCAT-PWY: L-rhamnose degradation I")

path_rhamnose_filtered_0.5 = path_rhamnose_filtered

# coverage가 0.5 미만인 경우는 abundance를 0으로 간주
path_rhamnose_filtered_0.5[path_rhamnose_filtered_cover < 0.5] = 0

path_rhamnose_filtered_long = path_rhamnose_filtered_0.5 %>% 
  pivot_longer(
    
    cols = -Pathway, 
    names_to = "SampleID", 
    values_to = "Abundance_0.5"
    
  ) %>% 
  filter(Abundance_0.5 > 0) %>%
  merge(., mb2 %>% select(SampleID, 
                          Age, Sex,
                          TRG_score, TRG_1, TRG_2, TRG_3, 
                          Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                          Pre_Op_Nstage, Pre_Op_Nstage_bin),
        by = "SampleID") %>% 
  separate(Pathway, into = c("Pathway_ID", "Taxonomy"), sep = "\\|", fill = "right") %>%
  
  # 2. Taxonomy 분할: g__와 s__가 있는 경우에만 genus/species 추출
  mutate(
    Genus = ifelse(str_detect(Taxonomy, "g__"),
                   str_extract(Taxonomy, "g__([^\\.]+)") %>% str_remove("g__"),
                   "unclassified"),
    Species = ifelse(str_detect(Taxonomy, "s__"),
                     str_extract(Taxonomy, "s__.+") %>% str_remove("s__"),
                     "unclassified")
  ) ; head(path_rhamnose_filtered_long)


# By genus & TRG_1
path_rhamnose_genus = path_rhamnose_filtered_long %>% 
  filter(Genus != "unclassified") %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance_0.5), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice(1:6) %>%
  pull(Genus)

path_rhamnose_barplot = path_rhamnose_filtered_long %>% 
  mutate(Genus_grouped = case_when(
    Genus %in% path_rhamnose_genus ~ Genus,
    Genus == "unclassified" ~ "unclassified",
    TRUE ~ "Others"
  ))

path_rhamnose_barplot %>% 
  filter(Genus == "Klebsiella") %>% 
  select(Species) %>% 
  unique()

path_rhamnose_barplot %>% 
  filter(Genus == "Blautia") %>% 
  select(Species) %>% 
  unique()

path_rhamnose_barplot %>% 
  filter(Genus == "Bacteroides") %>% 
  select(Species)


path_rhamnose_barplot_grouped = path_rhamnose_barplot %>%
  group_by(Genus_grouped, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15)) 


# 4. stacked barplot 그리기
genus_colors = c(
  "Bacteroides"         = "#F1A7C1",
  "Enterococcus"        = "#66A3D2",  
  "Blautia"             = "#9A4E97",
  "Klebsiella"          = "#D84B6A",  
  "Lactococcus"         = "#A0D468",  
  "Raoultella"          = "#E8C13B",  
  "Others"              = "#BFBFBF", 
  "unclassified"        = "#DDDDDD" 
)


p4_pathway_rhamnose_genus = path_rhamnose_barplot_grouped %>% 
  mutate(Genus_grouped = factor(Genus_grouped, 
                                levels = rev(c("Bacteroides",
                                               "Enterococcus",
                                               "Blautia",
                                               "Klebsiella",
                                               "Lactococcus",
                                               "Raoultella",
                                               "Others",
                                               "unclassified")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional abundance (RHAMCAT-PWY)",
       fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors) ; p4_pathway_rhamnose_genus

# ggsave("figure/04-3-19_rhamnose-genus.svg",
#        plot = p4_pathway_rhamnose_genus, width = 5, height = 5.5)


genus_colors = c(
  "Bacteroides"         = "#F1A7C1",
  "Enterococcus"        = "#66A3D2",  
  "Blautia"             = "#9A4E97",
  "Klebsiella"          = "#D84B6A",  
  "Lactococcus"         = "#A0D468",  
  "Raoultella"          = "#E8C13B",  
  "Others"              = "#BFBFBF"
)


p4_pathway_rhamnose_genus_woUnclassified = path_rhamnose_barplot_grouped %>% 
  filter(Genus_grouped != "unclassified") %>% 
  mutate(Genus_grouped = factor(Genus_grouped, 
                                levels = rev(c("Bacteroides",
                                               "Enterococcus",
                                               "Blautia",
                                               "Klebsiella",
                                               "Lactococcus",
                                               "Raoultella",
                                               "Others")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional abundance (RHAMCAT-PWY)",
       fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors) ; p4_pathway_rhamnose_genus_woUnclassified

# ggsave("figure/04-3-19_rhamnose-genus_woUnclassified.svg",
#        plot = p4_pathway_rhamnose_genus_woUnclassified, width = 5, height = 5.5)



# DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis
# Filtering
path_rhamnose_syn_filtered = path_abund %>%
  filter(str_detect(Pathway, "DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis")) %>% 
  filter(Pathway != "DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis")

path_rhamnose_syn_filtered_cover = path_cover %>%
  filter(str_detect(Pathway, "DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis"))%>% 
  filter(Pathway != "DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis")

path_rhamnose_syn_filtered_0.5 = path_rhamnose_syn_filtered

# coverage가 0.5 미만인 경우는 abudnance를 0으로 간주
path_rhamnose_syn_filtered_0.5[path_rhamnose_syn_filtered_cover < 0.5] = 0


path_rhamnose_syn_filtered_long = path_rhamnose_syn_filtered_0.5 %>% 
  pivot_longer(
    
    cols = -Pathway,
    names_to = "SampleID",
    values_to = "Abundance_0.5"
    
  ) %>% 
  filter(Abundance_0.5 > 0) %>% 
  merge(., mb2 %>% select(SampleID, 
                          Age, Sex,
                          TRG_score, TRG_1, TRG_2, TRG_3, 
                          Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                          Pre_Op_Nstage, Pre_Op_Nstage_bin),
        by = "SampleID") %>% 
  separate(Pathway, into = c("Pathway_ID", "Taxonomy"), sep = "\\|", fill = "right") %>% 
  
  # 2. Taxonomy 분할: g__와 s__가 있는 경우에만 genus/species 추출
  mutate(
    Genus = ifelse(str_detect(Taxonomy, "g__"),
                   str_extract(Taxonomy, "g__([^\\.]+)") %>% str_remove("g__"),
                   "unclassified"),
    Species = ifelse(str_detect(Taxonomy, "s__"),
                     str_extract(Taxonomy, "s__.+") %>% str_remove("s__"),
                     "unclassified")
  ) ; head(path_rhamnose_syn_filtered_long)


# By genus & TRG_1
path_rhamnose_syn_genus = path_rhamnose_syn_filtered_long %>% 
  filter(Genus != "unclassified") %>% 
  group_by(Genus) %>% 
  summarise(Total = sum(Abundance_0.5), .groups = "drop") %>% 
  arrange(desc(Total)) %>% 
  slice(1:6) %>% 
  pull(Genus)

path_rhamnose_syn_barplot = path_rhamnose_syn_filtered_long %>% 
  mutate(Genus_grouped = case_when(
    Genus %in% path_rhamnose_syn_genus ~ Genus,
    Genus == "unclassified" ~ "unclassified",
    TRUE ~ "Others"
  ))

path_rhamnose_syn_barplot %>% 
  filter(Genus == "Bacteroides") %>% 
  select(Species) %>% 
  unique()


path_rhamnose_syn_barplot_grouped = path_rhamnose_syn_barplot %>% 
  group_by(Genus_grouped, TRG_1) %>% 
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15)) 


# 4. stacked barplot 그리기
genus_colors = c(
  "Bacteroides"         = "#F1A7C1",
  "Blautia"             = "#9A4E97", 
  "Faecalibacterium"    = "#8BC34A",  
  "Collinsella"         = "#1C85C4",
  "Prevotella"          = "#78C38F",  
  "Lactobacillus"       = "#E17A6A",
  "Others"              = "#BFBFBF", 
  "unclassified"        = "#DDDDDD"
)


p4_pathway_rhamnose_syn_genus = path_rhamnose_syn_barplot_grouped %>% 
  mutate(Genus_grouped = factor(Genus_grouped, 
                                levels = rev(c("Bacteroides",
                                               "Blautia",
                                               "Faecalibacterium",
                                               "Collinsella",
                                               "Prevotella",
                                               "Lactobacillus",
                                               "Others",
                                               "unclassified")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional abundance (DTDPRHAMSYN-PWY)",
       fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors) ; p4_pathway_rhamnose_syn_genus

# ggsave("figure/04-3-19_rhamnose_syn-genus.svg",
#        plot = p4_pathway_rhamnose_syn_genus, width = 5, height = 5.5)


genus_colors = c(
  "Bacteroides"         = "#F1A7C1",
  "Blautia"             = "#9A4E97", 
  "Faecalibacterium"    = "#8BC34A",  
  "Collinsella"         = "#1C85C4",
  "Prevotella"          = "#78C38F",  
  "Lactobacillus"       = "#E17A6A",
  "Others"              = "#BFBFBF"
)


p4_pathway_rhamnose_syn_genus_woUnclassified = path_rhamnose_syn_barplot_grouped %>% 
  filter(Genus_grouped != "unclassified") %>% 
  mutate(Genus_grouped = factor(Genus_grouped, 
                                levels = rev(c("Bacteroides",
                                               "Blautia",
                                               "Faecalibacterium",
                                               "Collinsella",
                                               "Prevotella",
                                               "Lactobacillus",
                                               "Others")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional abundance (DTDPRHAMSYN-PWY)",
       fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors) ; p4_pathway_rhamnose_syn_genus_woUnclassified

# ggsave("figure/04-3-19_rhamnose_syn-genus-woUnclassified.svg",
#        plot = p4_pathway_rhamnose_syn_genus_woUnclassified, width = 5, height = 5.5)



# PWY-7560: methylerythritol phosphate pathway II
# Filtering
path_mep_filtered = path_abund %>%
  filter(str_detect(Pathway, "PWY-7560: methylerythritol phosphate pathway II")) %>% 
  filter(Pathway != "PWY-7560: methylerythritol phosphate pathway II")

path_mep_filtered_cover = path_cover %>%
  filter(str_detect(Pathway, "PWY-7560: methylerythritol phosphate pathway II"))%>% 
  filter(Pathway != "PWY-7560: methylerythritol phosphate pathway II")

path_mep_filtered_0.5 = path_mep_filtered

# coverage가 0.5 미만인 경우는 abundance를 0으로 간주
path_mep_filtered_0.5[path_mep_filtered_cover < 0.5] = 0


path_mep_filtered_long = path_mep_filtered_0.5 %>% 
  pivot_longer(
    
    cols = -Pathway, 
    names_to = "SampleID", 
    values_to = "Abundance_0.5"
    
  ) %>% 
  filter(Abundance_0.5 > 0) %>%
  merge(., mb2 %>% select(SampleID, 
                          Age, Sex,
                          TRG_score, TRG_1, TRG_2, TRG_3, 
                          Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                          Pre_Op_Nstage, Pre_Op_Nstage_bin),
        by = "SampleID") %>% 
  separate(Pathway, into = c("Pathway_ID", "Taxonomy"), sep = "\\|", fill = "right") %>%
  
  # 2. Taxonomy 분할: g__와 s__가 있는 경우에만 genus/species 추출
  mutate(
    Genus = ifelse(str_detect(Taxonomy, "g__"),
                   str_extract(Taxonomy, "g__([^\\.]+)") %>% str_remove("g__"),
                   "unclassified"),
    Species = ifelse(str_detect(Taxonomy, "s__"),
                     str_extract(Taxonomy, "s__.+") %>% str_remove("s__"),
                     "unclassified")
  ) ; head(path_mep_filtered_long)


# By genus & TRG_1
path_mep_genus = path_mep_filtered_long %>% 
  filter(Genus != "unclassified") %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance_0.5), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice(1:6) %>%
  pull(Genus)

path_mep_barplot = path_mep_filtered_long %>% 
  mutate(Genus_grouped = case_when(
    Genus %in% path_mep_genus ~ Genus,
    Genus == "unclassified" ~ "unclassified",
    TRUE ~ "Others"
  ))

path_mep_barplot %>% 
  filter(Genus == "Blautia") %>% 
  select(Species)


path_mep_barplot_grouped = path_mep_barplot %>%
  group_by(Genus_grouped, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15)) 


# 4. stacked barplot 그리기
genus_colors = c(
  "Blautia"             = "#9A4E97",
  "Bacteroides"         = "#F1A7C1",
  "Eubacterium"         = "#66B2E4",
  "Rothia"              = "#B5D9B0", 
  "Others"              = "#BFBFBF", 
  "unclassified"        = "#DDDDDD" 
)


p4_pathway_mep_genus = path_mep_barplot_grouped %>% 
  mutate(Genus_grouped = factor(Genus_grouped, 
                                levels = rev(c("Blautia",
                                               "Bacteroides",
                                               "Eubacterium",
                                               "Rothia",
                                               "Others",
                                               "unclassified")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional abundance (PWY-7560)",
       fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors) ; p4_pathway_mep_genus

# ggsave("figure/04-3-19_MEP-genus.svg",
#        plot = p4_pathway_mep_genus, width = 5, height = 5.5)


genus_colors = c(
  "Blautia"             = "#9A4E97",
  "Bacteroides"         = "#F1A7C1",
  "Eubacterium"         = "#66B2E4",
  "Rothia"              = "#B5D9B0", 
  "Others"              = "#BFBFBF"
)


p4_pathway_mep_genus_woUnclassified = path_mep_barplot_grouped %>% 
  filter(Genus_grouped != "unclassified") %>% 
  mutate(Genus_grouped = factor(Genus_grouped, 
                                levels = rev(c("Blautia",
                                               "Bacteroides",
                                               "Eubacterium",
                                               "Rothia",
                                               "Others")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional abundance (PWY-7560)",
       fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors) ; p4_pathway_mep_genus_woUnclassified

# ggsave("figure/04-3-19_MEP-genus-woUnclassified.svg",
#        plot = p4_pathway_mep_genus_woUnclassified, width = 5, height = 5.5)



# PWY-6270: isoprene biosynthesis I
# Filtering 
path_isoprene_filtered = path_abund %>%
  filter(str_detect(Pathway, "PWY-6270: isoprene biosynthesis I")) %>% 
  filter(Pathway != "PWY-6270: isoprene biosynthesis I")

path_isoprene_filtered_cover = path_cover %>%
  filter(str_detect(Pathway, "PWY-6270: isoprene biosynthesis I"))%>% 
  filter(Pathway != "PWY-6270: isoprene biosynthesis I")

path_isoprene_filtered_0.5 = path_isoprene_filtered

# coverage가 0.5 미만인 경우는 abundance를 0으로 간주
path_isoprene_filtered_0.5[path_isoprene_filtered_cover < 0.5] = 0


path_isoprene_filtered_long = path_isoprene_filtered_0.5 %>% 
  pivot_longer(
    
    cols = -Pathway, 
    names_to = "SampleID", 
    values_to = "Abundance_0.5"
    
  ) %>% 
  filter(Abundance_0.5 > 0) %>%
  merge(., mb2 %>% select(SampleID, 
                          Age, Sex,
                          TRG_score, TRG_1, TRG_2, TRG_3, 
                          Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                          Pre_Op_Nstage, Pre_Op_Nstage_bin),
        by = "SampleID") %>% 
  separate(Pathway, into = c("Pathway_ID", "Taxonomy"), sep = "\\|", fill = "right") %>%
  
  # 2. Taxonomy 분할: g__와 s__가 있는 경우에만 genus/species 추출
  mutate(
    Genus = ifelse(str_detect(Taxonomy, "g__"),
                   str_extract(Taxonomy, "g__([^\\.]+)") %>% str_remove("g__"),
                   "unclassified"),
    Species = ifelse(str_detect(Taxonomy, "s__"),
                     str_extract(Taxonomy, "s__.+") %>% str_remove("s__"),
                     "unclassified")
  ) ; head(path_isoprene_filtered_long)


# By genus & TRG_1
path_isoprene_genus = path_isoprene_filtered_long %>% 
  filter(Genus != "unclassified") %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance_0.5), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice(1:6) %>%
  pull(Genus)

path_isoprene_barplot = path_isoprene_filtered_long %>% 
  mutate(Genus_grouped = case_when(
    Genus %in% path_isoprene_genus ~ Genus,
    Genus == "unclassified" ~ "unclassified",
    TRUE ~ "Others"
  ))


path_isoprene_barplot_grouped = path_isoprene_barplot %>%
  group_by(Genus_grouped, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15)) 


# 4. stacked barplot 그리기
genus_colors = c(
  "Bacteroides"         = "#F1A7C1",
  "Blautia"             = "#9A4E97",
  "Eubacterium"         = "#66B2E4",
  "Rothia"              = "#B5D9B0", 
  "Others"              = "#BFBFBF", 
  "unclassified"        = "#DDDDDD" 
)


p4_pathway_isoprene_genus = path_isoprene_barplot_grouped %>% 
  mutate(Genus_grouped = factor(Genus_grouped, 
                                levels = rev(c("Bacteroides",
                                               "Blautia",
                                               "Eubacterium",
                                               "Rothia",
                                               "Others",
                                               "unclassified")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional abundance (PWY-6270)",
       fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors) ; p4_pathway_isoprene_genus

# ggsave("figure/04-3-19_isoprene-genus.svg",
#        plot = p4_pathway_isoprene_genus, width = 5, height = 5.5)


genus_colors = c(
  "Bacteroides"         = "#F1A7C1",
  "Blautia"             = "#9A4E97",
  "Eubacterium"         = "#66B2E4",
  "Rothia"              = "#B5D9B0", 
  "Others"              = "#BFBFBF"
)


p4_pathway_isoprene_genus_woUnclassified = path_isoprene_barplot_grouped %>% 
  filter(Genus_grouped != "unclassified") %>% 
  mutate(Genus_grouped = factor(Genus_grouped, 
                                levels = rev(c("Bacteroides",
                                               "Blautia",
                                               "Eubacterium",
                                               "Rothia",
                                               "Others")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Genus_grouped)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional abundance (PWY-7560)",
       fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = genus_colors) ; p4_pathway_isoprene_genus_woUnclassified

# ggsave("figure/04-3-19_isoprene-genus-woUnclassified.svg",
#        plot = p4_pathway_isoprene_genus_woUnclassified, width = 5, height = 5.5)



##### Species ###

# Anaerostipes
path_thiamine_species_A = path_thiamine_barplot %>%
  filter(Genus == "Anaerostipes") %>% 
  group_by(Species, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15))

path_sucrose_species_A = path_sucrose_barplot %>% 
  filter(Genus == "Anaerostipes") %>% 
  group_by(Species, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15))

path_sucrose4_species_A = path_sucrose4_barplot %>% 
  filter(Genus == "Anaerostipes") %>% 
  group_by(Species, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15))

species_colors = c(
  "Anaerostipes_hadrus"   = "#A3C37B"
)

p4_pathway_thiamine_species_A = path_thiamine_species_A %>% 
  mutate(Species = factor(Species)) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional Abundance (PWY-6892_Anaerostipes)",
       fill = "Species") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_thiamine_species_A

p4_pathway_sucrose_species_A = path_sucrose_species_A %>% 
  mutate(Species = factor(Species)) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional Abundance (PWY-621_Anaerostipes)",
       fill = "Species") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_sucrose_species_A

p4_pathway_sucrose4_species_A = path_sucrose4_species_A %>% 
  mutate(Species = factor(Species)) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional Abundance (PWY-5384_Anaerostipes)",
       fill = "Species") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_sucrose4_species_A


path_species_abund_A = sb[grep("^Anaerostipes_hadrus", sb$Species), ] %>% 
  pivot_longer(cols = -Species,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  merge(mb2 %>% 
          select(SampleID,
                 TRG_score, TRG_1, TRG_2, TRG_3,
                 Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                 Pre_Op_Nstage, Pre_Op_Nstage_bin,
                 BMI, BMI_bin, CEA, CEA_bin),
        by = "SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR"))) %>% 
  group_by(Species, TRG_1) %>% 
  summarise(total_abund = sum(abundance)) %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", total_abund / 11, total_abund / 15))
  

genus_colors = c(
  "Anaerostipes"   = "#A3C37B",
  "Blautia"        = "#9A4E97",
  "Bacteroides"    = "#F1A7C1",
  "Klebsiella"     = "#D84B6A",
  "Prevotella"     = "#78C38F",
  "Dialister"      = "#E06646"
)

p4_pathway_species_abund_A = path_species_abund_A %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(y = "Species Abundance", fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_species_abund_A

p4_pathway_species_Anaerostipes = ggarrange(p4_pathway_thiamine_species_A,
                                            p4_pathway_sucrose_species_A,
                                            p4_pathway_sucrose4_species_A,
                                            p4_pathway_species_abund_A,
                                            nrow = 1) ; p4_pathway_species_Anaerostipes

# ggsave("figure/04-3-20_pathway_species_Anaerostipes.svg",
#        plot = p4_pathway_species_Anaerostipes, width = 20, height = 5.5)


# Klebsiella
path_thiamine_species_K = path_thiamine_barplot %>%
  filter(Genus == "Klebsiella") %>% 
  group_by(Species, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15))

path_sucrose_species_K = path_sucrose_barplot %>%
  filter(Genus == "Klebsiella") %>% 
  group_by(Species, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15))

path_rhamnose_species_K = path_rhamnose_barplot %>%
  filter(Genus == "Klebsiella") %>% 
  group_by(Species, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15))

species_colors = c(
  "Klebsiella_pneumoniae"        = "#D84B6A",
  "Klebsiella_oxytoca"           = "#F3A8A0",
  "Klebsiella_quasipneumoniae"   = "#F1C2D4",
  "Klebsiella_aerogenes"         = "#F28D8C",
  "Others"                       = "#BFBFBF"
)

p4_pathway_thiamine_species_K = path_thiamine_species_K %>% 
  mutate(Species = factor(Species,
                          levels = rev(c("Klebsiella_pneumoniae",
                                         "Klebsiella_oxytoca")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional Abundance (PWY-6892_Klebsiella)",
       fill = "Species") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_thiamine_species_K

p4_pathway_sucrose_species_K = path_sucrose_species_K %>% 
  mutate(Species = factor(Species,
                          levels = rev(c("Klebsiella_pneumoniae",
                                         "Klebsiella_quasipneumoniae")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional Abundance (PWY-621_Klebsiella)",
       fill = "Species") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_sucrose_species_K

p4_pathway_rhamnose_species_K = path_rhamnose_species_K %>% 
  mutate(Species = factor(Species,
                          levels = rev(c("Klebsiella_pneumoniae",
                                         "Klebsiella_quasipneumoniae",
                                         "Klebsiella_aerogenes")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional Abundance (RHAMCAT-PWY_Klebsiella)",
       fill = "Species") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_rhamnose_species_K


path_species_abund_K = sb[grep("^Klebsiella_", sb$Species), ] %>% 
  pivot_longer(cols = -Species,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  merge(mb2 %>% 
          select(SampleID,
                 TRG_score, TRG_1, TRG_2, TRG_3,
                 Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                 Pre_Op_Nstage, Pre_Op_Nstage_bin,
                 BMI, BMI_bin, CEA, CEA_bin),
        by = "SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR"))) %>% 
  mutate(Species = case_when(
    Species == "Klebsiella_pneumoniae"      ~ "Klebsiella_pneumoniae",
    Species == "Klebsiella_oxytoca"         ~ "Klebsiella_oxytoca",
    Species == "Klebsiella_quasipneumoniae" ~ "Klebsiella_quasipneumoniae",
    Species == "Klebsiella_aerogenes"       ~ "Klebsiella_aerogenes",
    TRUE ~ "Others"  # 해당하지 않으면 "Others"
  )) %>% 
  group_by(Species, TRG_1) %>% 
  summarise(total_abund = sum(abundance)) %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", total_abund / 11, total_abund / 15))

genus_colors = c(
  "Anaerostipes"   = "#A3C37B",
  "Blautia"        = "#9A4E97",
  "Bacteroides"    = "#F1A7C1",
  "Klebsiella"     = "#D84B6A",
  "Prevotella"     = "#78C38F",
  "Dialister"      = "#E06646"
)

p4_pathway_species_abund_K = path_species_abund_K %>% 
  mutate(Species = factor(Species,
                          level = rev(c("Klebsiella_pneumoniae",
                                        "Klebsiella_oxytoca",
                                        "Klebsiella_quasipneumoniae",
                                        "Klebsiella_aerogenes",
                                        "Others")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(y = "Species Abundance", fill = "Genus") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_species_abund_K

p4_pathway_species_Klebsiella = ggarrange(p4_pathway_thiamine_species_K,
                                          p4_pathway_sucrose_species_K,
                                          p4_pathway_rhamnose_species_K,
                                          p4_pathway_species_abund_K,
                                          nrow = 1) ; p4_pathway_species_Klebsiella

# ggsave("figure/04-3-20_pathway_species_Klebsiella.svg",
#        plot = p4_pathway_species_Klebsiella, width = 20, height = 5.5)


# Blautia
path_rhamnose_species_B = path_rhamnose_barplot %>%
  filter(Genus == "Blautia") %>% 
  group_by(Species, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15))

path_mep_species_B = path_mep_barplot %>%
  filter(Genus == "Blautia") %>% 
  group_by(Species, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15))

path_thiamine_species_B = path_thiamine_barplot %>%
  filter(Genus == "Blautia") %>% 
  group_by(Species, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15))

species_colors = c(
  "Ruminococcus_gnavus"         = "#8E4B97",
  "Blautia_obeum"               = "#D16BB1",
  "Blautia_hydrogenotrophica"   = "#5A2E88",
  "Blautia_sp"                  = "#7D39A1",
  "Blautia_wexlerae"            = "#B257C4",
  "Others"                      = "#BFBFBF"
)

p4_pathway_rhamnose_species_B = path_rhamnose_species_B %>% 
  mutate(Species = factor(Species)) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional Abundance (RHAMCAT-PWY_Blautia)",
       fill = "Species") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_rhamnose_species_B

p4_pathway_mep_species_B = path_mep_species_B %>% 
  mutate(Species = factor(Species,
                          level = rev(c("Blautia_obeum",
                                        "Blautia_hydrogenotrophica")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional Abundance (PWY-7560_Blautia)",
       fill = "Species") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_mep_species_B

p4_pathway_thiamine_species_B = path_thiamine_species_B %>% 
  mutate(Species = factor(Species,
                          level = rev(c("Blautia_wexlerae",
                                        "Blautia_obeum",
                                        "Blautia_sp")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional Abundance (PWY-6892_Blautia)",
       fill = "Species") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_thiamine_species_B


path_species_abund_B = sb[grep("^Blautia|^Ruminococcus_gnavus", sb$Species), ] %>% 
  pivot_longer(cols = -Species,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  merge(mb2 %>% 
          select(SampleID,
                 TRG_score, TRG_1, TRG_2, TRG_3,
                 Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                 Pre_Op_Nstage, Pre_Op_Nstage_bin,
                 BMI, BMI_bin, CEA, CEA_bin),
        by = "SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR"))) %>% 
  mutate(Species = case_when(
    Species == "Blautia_wexlerae"          ~ "Blautia_wexlerae",
    Species == "Blautia_sp"                ~ "Blautia_sp",
    Species == "Ruminococcus_gnavus"       ~ "Ruminococcus_gnavus",
    Species == "Blautia_obeum"             ~ "Blautia_obeum",
    Species == "Blautia_hydrogenotrophica" ~ "Blautia_hydrogenotrophica",
    TRUE ~ "Others"  # 해당하지 않으면 "Others"
  )) %>% 
  group_by(Species, TRG_1) %>% 
  summarise(total_abund = sum(abundance)) %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", total_abund / 11, total_abund / 15))

genus_colors = c(
  "Anaerostipes"   = "#A3C37B",
  "Blautia"        = "#9A4E97",
  "Bacteroides"    = "#F1A7C1",
  "Klebsiella"     = "#D84B6A",
  "Prevotella"     = "#78C38F",
  "Dialister"      = "#E06646"
)

p4_pathway_species_abund_B = path_species_abund_B %>% 
  mutate(Species = factor(Species,
                          level = rev(c("Blautia_wexlerae", 
                                        "Blautia_obeum",
                                        "Blautia_sp",
                                        "Ruminococcus_gnavus",
                                        "Blautia_hydrogenotrophica",
                                        "Others")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(y = "Species Abundance", fill = "Species") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_species_abund_B

p4_pathway_species_Blautia = ggarrange(p4_pathway_thiamine_species_B,
                                       p4_pathway_rhamnose_species_B,
                                       p4_pathway_mep_species_B,
                                       p4_pathway_species_abund_B,
                                       nrow = 1) ; p4_pathway_species_Blautia

# ggsave("figure/04-3-20_pathway_species_Blautia.svg",
#        plot = p4_pathway_species_Blautia, width = 20, height = 5.5)


# Bacteroides
path_rhamnose_syn_species_Bac_list = path_rhamnose_syn_barplot %>%
  filter(Genus == "Bacteroides") %>% 
  group_by(Species) %>% 
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  filter(Total_Abund > 1500) %>% 
  pull(Species)

path_rhamnose_species_Bac = path_rhamnose_barplot %>%
  filter(Genus == "Bacteroides") %>% 
  group_by(Species, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15))

path_rhamnose_syn_species_Bac = path_rhamnose_syn_barplot %>%
  filter(Genus == "Bacteroides") %>% 
  group_by(Species, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15))


path_rhamnose_species_Bac_barplot = path_rhamnose_species_Bac %>% 
  mutate(Species_grouped = case_when(
    Species %in% path_rhamnose_syn_species_Bac_list ~ Species,
    TRUE ~ "Others"
  ))

path_rhamnose_syn_species_Bac_barplot = path_rhamnose_syn_species_Bac %>% 
  mutate(Species_grouped = case_when(
    Species %in% path_rhamnose_syn_species_Bac_list ~ Species,
    TRUE ~ "Others"
  ))

species_colors = c(
  "Bacteroides_vulgatus"         = "#F1A7C1", # [고정] 중간 핑크
  "Phocaeicola_vulgatus"         = "#F1A7C1",
  "Bacteroides_thetaiotaomicron" = "#FFD1DC", # 아주 연한 핑크
  "Bacteroides_xylanisolvens"    = "#D16A7D", # 톤다운 로즈 (중간-진함)
  "Bacteroides_ovatus"           = "#FAD0C4", # 연한 살구 코랄 (톤 변화)
  "Bacteroides_uniformis"        = "#EC8B99", # 선명한 코랄 핑크
  "Bacteroides_caccae"           = "#FBCFE8", # 밝은 라벤더 핑크 (톤 변화)
  "Bacteroides_faecis"           = "#A95C68", # 딥 로즈 (확실히 어두움)
  "Bacteroides_fragilis"         = "#E56B6F", # [고정] 코랄 로즈
  "Bacteroides_cellulosilyticus" = "#FFB3BA", # 파스텔 핑크
  "Bacteroides_dorei"            = "#C14450", # 루비 레드 (확실히 진함)
  "Phocaeicola_dorei"            = "#C14450",
  "Bacteroides_plebeius"         = "#C9A9A9", # 중간 톤 핑크
  "Phocaeicola_plebeius"         = "#C9A9A9",
  "Bacteroides_nordii"           = "#8E5B66", # 더스티 퍼플 로즈 (가장 어두움)
  "Others"                       = "#BFBFBF"  # 회색
)

p4_pathway_rhamnose_species_Bac = path_rhamnose_species_Bac_barplot %>% 
  mutate(Species = factor(Species_grouped,
                          levels = rev(c("Bacteroides_vulgatus",
                                         "Bacteroides_thetaiotaomicron",
                                         "Bacteroides_xylanisolvens",
                                         "Bacteroides_ovatus",
                                         "Bacteroides_uniformis",
                                         "Bacteroides_caccae",
                                         "Bacteroides_faecis",
                                         "Others")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional Abundance (RHAMCAT-PWY_Bacteroides)",
       fill = "Species") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_rhamnose_species_Bac

p4_pathway_rhamnose_syn_species_Bac = path_rhamnose_syn_species_Bac_barplot %>% 
  mutate(Species = factor(Species_grouped,
                          level = rev(c("Bacteroides_vulgatus",
                                        "Bacteroides_thetaiotaomicron",
                                        "Bacteroides_xylanisolvens",
                                        "Bacteroides_ovatus",
                                        "Bacteroides_uniformis",
                                        "Bacteroides_caccae",
                                        "Bacteroides_fragilis",
                                        "Bacteroides_cellulosilyticus", 
                                        "Bacteroides_dorei", 
                                        "Bacteroides_plebeius", 
                                        "Bacteroides_nordii",
                                        "Others")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional Abundance (DTDPRHAMSYN-PWY_Bacteroides)",
       fill = "Species") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_rhamnose_syn_species_Bac


path_species_abund_Bac = sb[grep("^Bacteroides|Phocaeicola_vulgatus|Phocaeicola_dorei|Phocaeicola_plebeius",
                                 sb$Species), ] %>% 
  pivot_longer(cols = -Species,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  merge(mb2 %>% 
          select(SampleID,
                 TRG_score, TRG_1, TRG_2, TRG_3,
                 Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                 Pre_Op_Nstage, Pre_Op_Nstage_bin,
                 BMI, BMI_bin, CEA, CEA_bin),
        by = "SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR"))) %>% 
  mutate(Species = case_when(
    Species == "Phocaeicola_vulgatus"      ~ "Phocaeicola_vulgatus",
    Species == "Bacteroides_thetaiotaomicron"  ~ "Bacteroides_thetaiotaomicron",
    Species == "Bacteroides_xylanisolvens" ~ "Bacteroides_xylanisolvens",
    Species == "Bacteroides_ovatus"       ~ "Bacteroides_ovatus",
    Species == "Bacteroides_uniformis" ~ "Bacteroides_uniformis",
    Species == "Bacteroides_caccae" ~ "Bacteroides_caccae",
    Species == "Bacteroides_fragilis" ~ "Bacteroides_fragilis",
    Species == "Bacteroides_cellulosilyticus" ~ "Bacteroides_cellulosilyticus",
    Species == "Phocaeicola_dorei" ~ "Phocaeicola_dorei",
    Species == "Phocaeicola_plebeius" ~ "Phocaeicola_plebeius",
    Species == "Bacteroides_nordii" ~ "Bacteroides_nordii",
    TRUE ~ "Others"  # 해당하지 않으면 "Others"
  )) %>% 
  group_by(Species, TRG_1) %>% 
  summarise(total_abund = sum(abundance)) %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", total_abund / 11, total_abund / 15))

genus_colors = c(
  "Anaerostipes"   = "#A3C37B",
  "Blautia"        = "#9A4E97",
  "Bacteroides"    = "#F1A7C1",
  "Klebsiella"     = "#D84B6A",
  "Prevotella"     = "#78C38F",
  "Dialister"      = "#E06646"
)

p4_pathway_species_abund_Bac = path_species_abund_Bac %>% 
  mutate(Species = factor(Species,
                          level = rev(c("Phocaeicola_vulgatus",
                                        "Bacteroides_thetaiotaomicron",
                                        "Bacteroides_xylanisolvens",
                                        "Bacteroides_ovatus",
                                        "Bacteroides_uniformis",
                                        "Bacteroides_caccae",
                                        "Bacteroides_fragilis",
                                        "Bacteroides_cellulosilyticus",
                                        "Phocaeicola_dorei",
                                        "Phocaeicola_plebeius",
                                        "Bacteroides_nordii",
                                        "Others")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(y = "Species Abundance", fill = "Species") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_species_abund_Bac

p4_pathway_species_Bacteroides = ggarrange(p4_pathway_rhamnose_species_Bac,
                                           p4_pathway_rhamnose_syn_species_Bac,
                                           p4_pathway_species_abund_Bac,
                                           nrow = 1) ; p4_pathway_species_Bacteroides

# ggsave("figure/04-3-20_pathway_species_Bacteroides.svg",
#        plot = p4_pathway_species_Bacteroides, width = 15, height = 5.5)


# Prevotella
path_histidine_species_P = path_histidine_barplot %>%
  filter(Genus == "Prevotella") %>% 
  group_by(Species, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15))

path_rhamnose_syn_species_P = path_rhamnose_syn_barplot %>%
  filter(Genus == "Prevotella") %>% 
  group_by(Species, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15))

species_colors = c(
  "Prevotella_sp_CAG_5226"     = "#78C38F",
  "Prevotella_copri"           = "#2A8E74",
  "Segatella_copri"            = "#2A8E74",
  "Prevotella_stercorea"       = "#8ACF6A",
  "Prevotella_sp_CAG_279"      = "#6E9F75",
  "Others"                     = "#BFBFBF"
)

p4_pathway_histidine_species_P = path_histidine_species_P %>% 
  mutate(Species = factor(Species)) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional Abundance (HISTSYN-PWY_Prevotella)",
       fill = "Species") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_histidine_species_P

p4_pathway_rhamnose_syn_species_P = path_rhamnose_syn_species_P %>% 
  mutate(Species = factor(Species,
                          level = rev(c("Prevotella_copri",
                                        "Prevotella_stercorea",
                                        "Prevotella_sp_CAG_279")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional Abundance (DTDPRHAMSYN-PWY_Prevotella)",
       fill = "Species") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_rhamnose_syn_species_P


path_species_abund_P = sb[grep("^Prevotella|Segatella_copri", sb$Species), ] %>% 
  pivot_longer(cols = -Species,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  merge(mb2 %>% 
          select(SampleID,
                 TRG_score, TRG_1, TRG_2, TRG_3,
                 Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                 Pre_Op_Nstage, Pre_Op_Nstage_bin,
                 BMI, BMI_bin, CEA, CEA_bin),
        by = "SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR"))) %>% 
  mutate(Species = case_when(
    Species == "Prevotella_sp_CAG_5226"      ~ "Others",
    Species == "Prevotella_copri" ~ "Prevotella_copri",
    Species == "Segatella_copri"         ~ "Segatella_copri",
    Species == "Prevotella_stercorea" ~ "Prevotella_stercorea",
    Species == "Prevotella_sp_CAG_279"       ~ "Others",
    TRUE ~ "Others"  # 해당하지 않으면 "Others"
  )) %>% 
  group_by(Species, TRG_1) %>% 
  summarise(total_abund = sum(abundance)) %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", total_abund / 11, total_abund / 15))

genus_colors = c(
  "Anaerostipes"   = "#A3C37B",
  "Blautia"        = "#9A4E97",
  "Bacteroides"    = "#F1A7C1",
  "Klebsiella"     = "#D84B6A",
  "Prevotella"     = "#78C38F",
  "Dialister"      = "#E06646"
)

p4_pathway_species_abund_P = path_species_abund_P %>% 
  mutate(Species = factor(Species,
                          level = rev(c("Segatella_copri",
                                        "Prevotella_stercorea",
                                        "Others")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(y = "Species Abundance", fill = "Species") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_species_abund_P

p4_pathway_species_Prevotella = ggarrange(p4_pathway_histidine_species_P,
                                          p4_pathway_rhamnose_syn_species_P,
                                          p4_pathway_species_abund_P,
                                          nrow = 1) ; p4_pathway_species_Prevotella

# ggsave("figure/04-3-20_pathway_species_Prevotella.svg",
#        plot = p4_pathway_species_Prevotella, width = 15, height = 5.5)


# Dialister
path_histidine_species_D = path_histidine_barplot %>%
  filter(Genus == "Dialister") %>% 
  group_by(Species, TRG_1) %>%
  summarise(Total_Abund = sum(Abundance_0.5), .groups = "drop") %>% 
  mutate(normal.abund = ifelse(TRG_1 == "CR", Total_Abund / 11, Total_Abund / 15))

species_colors = c(
  "Dialister_succinatiphilus"  = "#E06646",
  "Dialister_sp_CAG_357"       = "#f29a80",
  "Dialister_sp_CAG_486"       = "#fdccbe",
  "Others"                     = "#BFBFBF"
)

p4_pathway_histidine_species_D = path_histidine_species_D %>% 
  mutate(Species = factor(Species,
                          levels = rev(c("Dialister_succinatiphilus",
                                         "Dialister_sp_CAG_357",
                                         "Dialister_sp_CAG_486")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "Functional Abundance (HISTSYN-PWY_Dialister)",
       fill = "Species") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_histidine_species_D


path_species_abund_D = sb[grep("^Dialister", sb$Species), ] %>% 
  pivot_longer(cols = -Species,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  merge(mb2 %>% 
          select(SampleID,
                 TRG_score, TRG_1, TRG_2, TRG_3,
                 Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                 Pre_Op_Nstage, Pre_Op_Nstage_bin,
                 BMI, BMI_bin, CEA, CEA_bin),
        by = "SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR"))) %>% 
  mutate(Species = case_when(
    Species == "Dialister_succinatiphilus"      ~ "Dialister_succinatiphilus",
    Species == "Dialister_sp_CAG_357"         ~ "Dialister_sp_CAG_357",
    Species == "Dialister_sp_CAG_486" ~ "Dialister_sp_CAG_486",
    TRUE ~ "Others"  # 해당하지 않으면 "Others"
  )) %>% 
  group_by(Species, TRG_1) %>% 
  summarise(total_abund = sum(abundance)) %>% 
mutate(normal.abund = ifelse(TRG_1 == "CR", total_abund / 11, total_abund / 15))

genus_colors = c(
  "Anaerostipes"   = "#A3C37B",
  "Blautia"        = "#9A4E97",
  "Bacteroides"    = "#F1A7C1",
  "Klebsiella"     = "#D84B6A",
  "Prevotella"     = "#78C38F",
  "Dialister"      = "#E06646"
)

p4_pathway_species_abund_D = path_species_abund_D %>% 
  mutate(Species = factor(Species,
                          level = rev(c("Dialister_succinatiphilus",
                                        "Others")))) %>% 
  ggplot(aes(x = TRG_1, y = normal.abund, fill = Species)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(y = "Genus Abundance", fill = "Species") +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(size = rel(0.7)),
    legend.position = "right"
  ) +
  scale_fill_manual(values = species_colors) ; p4_pathway_species_abund_D

p4_pathway_species_Dialister = ggarrange(p4_pathway_histidine_species_D,
                                         p4_pathway_species_abund_D,
                                         nrow = 1) ; p4_pathway_species_Dialister

# ggsave("figure/04-3-20_pathway_species_Dialister.svg",
#        plot = p4_pathway_species_Dialister, width = 10, height = 5.5)




# Remove unnecessary objects
rm(list = ls(pattern = "^p4_")); rm(list = ls(pattern = "^p5_"))
save.image(file = "input/R_image/7-3. after-differential-enrichment-testing.RData")



