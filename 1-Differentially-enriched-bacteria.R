


### CRC Read-based metagenomics analysis ###


### 1-1. Setting

rm(list = ls()) 
options( java.parameters = "-Xmx64g" ) 



getwd() # 현재 R세션의 작업 디렉토리 확인
setwd( dir = "C:/Users/user/Desktop/충남대/4. 연구/008-CRC metagenomics" ) 


### 1-2. Load 'packages'

# Import one by one (you can use 'require()' function as well)

library(tidyverse) 
library(magrittr) # 데이터처리와 분석파이프라인 효율증대 (%>% 제공)
library(ggpubr)   # ggplot2 기반 출판물 수준의 그래픽 생성을 도움
library(ggforce)  # ggplot2 확장 패키지


# save.image(file = "Differentially-enriched-bacteria.RData")  
# load("Differentially-enriched-bacteria.RData")




#--------------- 1. Import ---------------# 

m <- read.csv("input/CRC_meta.csv") %>% 
  filter(!is.na(SampleID),
         !is.na(TRG_score),
         !SampleID %in% c(1, 6, 30)) %>% 
  mutate(TRG_2 = ifelse(TRG_score > 0, "nonCR", "CR"),
         SampleID_2 = paste0("Sample_", SampleID)); head(m)

f <- read.csv("input/family.csv") %>% 
  select(c("clade_name", m$SampleID_2)) %>% 
  mutate(Family = str_extract(clade_name, "(?<=f__).*")) %>% 
  select(-clade_name) ; head(g)

g <- read.csv("input/genus.csv") %>% 
  select(c("clade_name", m$SampleID_2)) %>% 
  mutate(Genus = str_extract(clade_name, "(?<=g__).*")) %>% 
  select(-clade_name) ; head(g)

s <- read.csv("input/species.csv") %>% 
  select(c("clade_name", m$SampleID_2)) %>% 
  mutate(Species = str_extract(clade_name, "(?<=s__).*")) %>% 
  select(-clade_name) ; head(s)

t <- read.csv("input/strain.csv") %>% 
  select(c("clade_name", m$SampleID_2))%>% 
  mutate(Strain = str_extract(clade_name, "(?<=s__).*")) %>% 
  select(-clade_name) ; head(t)

# p_abund <- read_tsv("input/pathway_abundance_matrix_simplified.tsv") %>% 
#   select(c("clade_name", m$SampleID_2))
# 
# p_cover <- read_tsv("input/pathway_coverage_matrix_simplified.tsv") %>% 
#   select(c("clade_name", m$SampleID_2))




#--------------- 2. Differentially-enriched bacteria ---------------# 


  
#----- 2-1. Family -----#
  
# Prevalence filtering 
f2 <- f %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Family) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence > nrow(m)*0.5)



# Merging dataframes
merged_f <- f %>%
  filter(Family %in% f2$Family) %>% 
  pivot_longer(cols = -Family, 
               names_to = "SampleID_2",
               values_to = "abundance") %>%
  merge(m %>% select(SampleID_2, TRG_score, TRG_2), by = "SampleID_2") %>% 
  mutate(TRG_2 = factor(TRG_2, levels = c("CR", "nonCR")))

# merged_f$abundance[merged_f$abundance > 0 ] %>% min() # 6e-05

# Wilcoxon test for TRG_2
wilcox_f <- merged_f %>%
  group_by(Family) %>%
  summarise(p_value = wilcox.test(abundance ~ TRG_2, exact = FALSE)$p.value,
            .groups = 'drop') %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>% 
  arrange(p_value)

sig_family <- wilcox_f %>% 
  filter(p_value < 0.05) %>% 
  pull(Family) 

sig_family <- merged_f %>% 
  filter(Family %in% sig_family) %>% 
  group_by(Family) %>% 
  summarise(avg_abund = mean(abundance)) %>%
  arrange(-avg_abund) %>% 
  pull(Family)


# TRG_2: CR vs. nonCR
p1_TRG_2_family <- merged_f %>% 
  filter(Family %in% sig_family) %>% 
  mutate(Family = factor(Family, levels = sig_family)) %>% 
  ggplot(aes(TRG_2, log10(abundance + 1.5e-05)))  +
  geom_jitter(aes(color = TRG_2)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0) +
  facet_wrap(~Family, nrow = 3)  +
  scale_color_manual(values = c("CR" = "#00BE67",
                                "nonCR" = "coral3")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")), 
                     size = 4,
                     method = "wilcox",
                     hide.ns = T, 
                     tip.length = 0) +
  theme_pubr() +
  coord_flip() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.1))) +
  labs(x = "Tumor Regression Grade",
       y = "Log10 (Abundance + pseudo)")

ggsave("figures/01-TRG_2-family.svg", device = "svg", plot = p1_TRG_2_family,
       width = 5, height = 5)


# Continuous scale using TRG_score
p1_TRG_score_family <- merged_f %>% 
  filter(Family %in% sig_family) %>% 
  mutate(Family = factor(Family, levels = sig_family)) %>% 
  ggplot(aes(TRG_score, log10(abundance + 1.5e-05)))  +
  geom_jitter(aes(fill = TRG_score), shape = 21, size = 2) +
  geom_smooth(method = "lm", se = 0.95, color = "gray50", alpha = 0.3) +
  scale_fill_gradient(low = "#00BE67",
                      high = "coral3") +
  facet_wrap(~Family, 
             nrow = 1) +
  stat_cor(method = "spearman") +
  theme_pubr() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.05)),
        axis.title = element_text(size = rel(1.2)),
        aspect.ratio = 0.9) +
  labs(x = "Tumor Regression Grade",
       y = "Log10 (Abundance + pseudo)")

ggsave("figures/01-TRG_score-family.svg", device = "svg", 
       plot = p1_TRG_score_family,
       width = 8, height = 3)




#----- 2-2. Genus -----#

# Prevalence filtering 
g2 <- g %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Genus) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence > nrow(m)*0.5)



# Merging dataframes
merged_g <- g %>%
  filter(Genus %in% g2$Genus) %>% 
  pivot_longer(cols = -Genus, 
               names_to = "SampleID_2",
               values_to = "abundance") %>%
  merge(m %>% select(SampleID_2, TRG_score, TRG_2), by = "SampleID_2") %>% 
  mutate(TRG_2 = factor(TRG_2, levels = c("CR", "nonCR")))

# merged_f$abundance[merged_f$abundance > 0 ] %>% min() # 6e-05

# Wilcoxon test for TRG_2
wilcox_g <- merged_g %>%
  group_by(Genus) %>%
  summarise(p_value = wilcox.test(abundance ~ TRG_2, exact = FALSE)$p.value,
            .groups = 'drop') %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>% 
  arrange(p_value)

sig_genus <- wilcox_g %>% 
  filter(p_value < 0.05) %>% 
  pull(Genus) 

sig_genus <- merged_g %>% 
  filter(Genus %in% sig_genus) %>% 
  group_by(Genus) %>% 
  summarise(avg_abund = mean(abundance)) %>%
  arrange(-avg_abund) %>% 
  pull(Genus)


# TRG_2: CR vs. nonCR
p1_TRG_2_genus <- merged_g %>% 
  filter(Genus %in% sig_genus) %>% 
  mutate(Genus = factor(Genus, levels = sig_genus)) %>% 
  ggplot(aes(TRG_2, log10(abundance + 1.5e-05)))  +
  geom_jitter(aes(color = TRG_2)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0) +
  facet_wrap(~Genus, nrow = 5)  +
  scale_color_manual(values = c("CR" = "#00BE67",
                                "nonCR" = "coral3")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")), 
                     size = 4,
                     method = "wilcox",
                     hide.ns = T, 
                     tip.length = 0) +
  theme_pubr() +
  coord_flip() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.1))) +
  labs(x = "Tumor Regression Grade",
       y = "Log10 (Abundance + pseudo)")

ggsave("figures/01-TRG_2-genus.svg", device = "svg", plot = p1_TRG_2_genus,
       width = 8, height = 6)


# Continuous scale using TRG_score
p1_TRG_score_genus <- merged_g %>% 
  filter(Genus %in% sig_genus) %>% 
  mutate(Genus = factor(Genus, levels = sig_genus)) %>% 
  ggplot(aes(TRG_score, log10(abundance + 1.5e-05)))  +
  geom_jitter(aes(fill = TRG_score), shape = 21, size = 2) +
  geom_smooth(method = "lm", se = 0.95, color = "gray50", alpha = 0.3) +
  scale_fill_gradient(low = "#00BE67",
                      high = "coral3") +
  facet_wrap(~Genus, 
             nrow = 2) +
  stat_cor(method = "spearman") +
  theme_pubr() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.05)),
        axis.title = element_text(size = rel(1.2)),
        aspect.ratio = 0.9) +
  labs(x = "Tumor Regression Grade",
       y = "Log10 (Abundance + pseudo)")

ggsave("figures/01-TRG_score-genus.svg", device = "svg", 
       plot = p1_TRG_score_genus,
       width = 11, height = 5)




#----- 2-3. Species -----# 

# Prevalence filtering 
s2 <- s %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Species) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence > nrow(m)*0.5)



# Merging dataframes
merged_s <- s %>%
  filter(Species %in% s2$Species) %>% 
  pivot_longer(cols = -Species, 
               names_to = "SampleID_2",
               values_to = "abundance") %>%
  merge(m %>% select(SampleID_2, TRG_score, TRG_2), by = "SampleID_2") %>% 
  mutate(TRG_2 = factor(TRG_2, levels = c("CR", "nonCR")))

# merged_s$abundance[merged_s$abundance > 0 ] %>% min() # 3e-05

# Wilcoxon test for TRG_2
wilcox_s <- merged_s %>%
  group_by(Species) %>%
  summarise(p_value = wilcox.test(abundance ~ TRG_2, exact = FALSE)$p.value,
            .groups = 'drop') %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>% 
  arrange(p_value)

sig_species <- wilcox_s %>% 
  filter(p_value < 0.05) %>% 
  pull(Species) 

sig_species <- merged_s %>% 
  filter(Species %in% sig_species) %>% 
  group_by(Species) %>% 
  summarise(avg_abund = mean(abundance)) %>%
  arrange(-avg_abund) %>% 
  pull(Species)


# TRG_2: CR vs. nonCR
p1_TRG_2_species <- merged_s %>% 
  filter(Species %in% sig_species) %>% 
  mutate(Species = factor(Species, levels = sig_species)) %>% 
  ggplot(aes(TRG_2, log10(abundance + 1.5e-05)))  +
  geom_jitter(aes(color = TRG_2)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0) +
  facet_wrap(~Species, nrow = 8)  +
  scale_color_manual(values = c("CR" = "#00BE67",
                                "nonCR" = "coral3")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")), 
                     size = 4,
                     method = "wilcox",
                     hide.ns = T, 
                     tip.length = 0) +
  theme_pubr() +
  coord_flip() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.1))) +
  labs(x = "Tumor Regression Grade",
       y = "Log10 (Abundance + pseudo)")

ggsave("figures/01-TRG_2-species.svg", device = "svg", plot = p1_TRG_2_species,
       width = 9, height = 9)



# Continuous scale using TRG_score
p1_TRG_score_species <- merged_s %>% 
  filter(Species %in% sig_species) %>% 
  mutate(Species = factor(Species, levels = sig_species)) %>% 
  ggplot(aes(TRG_score, log10(abundance + 1.5e-05)))  +
  geom_jitter(aes(fill = TRG_score), shape = 21, size = 2) +
  geom_smooth(method = "lm", se = 0.95, color = "gray50", alpha = 0.3) +
  scale_fill_gradient(low = "#00BE67",
                      high = "coral3") +
  facet_wrap(~Species, 
             nrow = 2) +
  stat_cor(method = "spearman") +
  theme_pubr() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.05)),
        axis.title = element_text(size = rel(1.2)),
        aspect.ratio = 0.9) +
  labs(x = "Tumor Regression Grade",
       y = "Log10 (Abundance + pseudo)")

ggsave("figures/01-TRG_score-species.svg", device = "svg", 
       plot = p1_TRG_score_species,
       width = 20, height = 6)



#----- 2-4. Strain -----#

# Prevalence filtering 
t2 <- t %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Strain) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence > nrow(m)*0.5)



# Merging dataframes
merged_t <- t %>%
  filter(Strain %in% t2$Strain) %>% 
  pivot_longer(cols = -Strain, 
               names_to = "SampleID_2",
               values_to = "abundance") %>%
  merge(m %>% select(SampleID_2, TRG_score, TRG_2), by = "SampleID_2") %>% 
  mutate(TRG_2 = factor(TRG_2, levels = c("CR", "nonCR")))

# merged_t$abundance[merged_t$abundance > 0 ] %>% min() # 3e-05

# Wilcoxon test for TRG_2
wilcox_t <- merged_t %>%
  group_by(Strain) %>%
  summarise(p_value = wilcox.test(abundance ~ TRG_2, exact = FALSE)$p.value,
            .groups = 'drop') %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>% 
  arrange(p_value)

sig_strain <- wilcox_t %>% 
  filter(p_value < 0.05) %>% 
  pull(Strain) 

sig_strain <- merged_t %>% 
  filter(Strain %in% sig_strain) %>% 
  group_by(Strain) %>% 
  summarise(avg_abund = mean(abundance)) %>%
  arrange(-avg_abund) %>% 
  pull(Strain)


# TRG_2: CR vs. nonCR
p1_TRG_2_strain <- merged_t %>% 
  filter(Strain %in% sig_strain) %>% 
  mutate(Strain = factor(Strain, levels = sig_strain)) %>% 
  ggplot(aes(TRG_2, log10(abundance + 1.5e-05)))  +
  geom_jitter(aes(color = TRG_2)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0) +
  facet_wrap(~Strain, nrow = 8)  +
  scale_color_manual(values = c("CR" = "#00BE67",
                                "nonCR" = "coral3")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")), 
                     size = 4,
                     method = "wilcox",
                     hide.ns = T, 
                     tip.length = 0) +
  theme_pubr() +
  coord_flip() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.1)),
        axis.title = element_text(size = rel(1.1))) +
  labs(x = "Tumor Regression Grade",
       y = "Log10 (Abundance + pseudo)")

ggsave("figures/01-TRG_2-strain.svg", device = "svg", plot = p1_TRG_2_strain,
       width = 12, height = 8)



# Continuous scale using TRG_score
p1_TRG_score_strain <- merged_t %>% 
  filter(Strain %in% sig_strain) %>% 
  mutate(Strain = factor(Strain, levels = sig_strain)) %>% 
  ggplot(aes(TRG_score, log10(abundance + 1.5e-05)))  +
  geom_jitter(aes(fill = TRG_score), shape = 21, size = 2) +
  geom_smooth(method = "lm", se = 0.95, color = "gray50", alpha = 0.3) +
  scale_fill_gradient(low = "#00BE67",
                      high = "coral3") +
  facet_wrap(~Strain, 
             nrow = 2) +
  stat_cor(method = "spearman") +
  theme_pubr() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.05)),
        axis.title = element_text(size = rel(1.2)),
        aspect.ratio = 0.9) +
  labs(x = "Tumor Regression Grade",
       y = "Log10 (Abundance + pseudo)")

ggsave("figures/01-TRG_score-strain.svg", device = "svg", 
       plot = p1_TRG_score_strain,
       width = 20, height = 5)


