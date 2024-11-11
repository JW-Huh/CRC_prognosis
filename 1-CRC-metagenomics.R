


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


# save.image(file = "Biobigdata-1.RData")  
# load("Biobigdata-1.RData")




#--------------- 1. Import ---------------# 

m <- read.csv("input/CRC_meta.csv") %>% 
  filter(!is.na(SampleID),
         !is.na(TRG_score),
         !SampleID %in% c(1, 6, 30)) %>% 
  mutate(TRG_2 = ifelse(TRG_score > 0, "nonCR", "CR"),
         TRG_3 = ifelse(TRG_score > 1, "bad", "nearCR"),
         SampleID_2 = paste0("Sample_", SampleID)); head(m)

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

head(m)

 # Prevalence filtering 

s2 <- s %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Species) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence > 21)


merged_s <- s %>%
  filter(Species %in% s2$Species) %>% 
  pivot_longer(cols = -Species, 
               names_to = "SampleID_2",
               values_to = "abundance") %>%
  merge(m %>% select(SampleID_2, TRG_score, TRG_2, TRG_3), by = "SampleID_2") %>% 
  mutate(TRG_2 = factor(TRG_2, levels = c("CR", "nonCR")),
         TRG_3 = factor(TRG_3, levels = c("nearCR", "bad")))



merged_s$abundance[merged_s$abundance > 0 ] %>% min() # 3e-05

min(merged_s$abundance)

# TRG_2: CR vs. nonCR
merged_s %>% 
  ggplot(aes(TRG_2, log10(abundance + 1.5e-05)))  +
  geom_jitter(aes(color = TRG_2)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0) +
  facet_wrap(~Species, 
             # scales = "free_y", 
             nrow = 16) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")), 
                     method = "wilcox",
                     hide.ns = T, 
                     tip.length = 0) +
  theme_minimal() +
  coord_flip() +
  theme(legend.position = "none")

# TRG_3: nearCR vs. bad
merged_s %>% 
  ggplot(aes(TRG_3, log10(abundance + 1.5e-05)))  +
  geom_jitter(aes(color = TRG_3)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0) +
  facet_wrap(~Species, 
             # scales = "free_y", 
             nrow = 16) +
  stat_compare_means(comparisons = list(c("nearCR", "bad")), 
                     method = "wilcox",
                     hide.ns = T, 
                     tip.length = 0) +
  theme_minimal() +
  coord_flip() +
  theme(legend.position = "none")


merged_s %>% 
  ggplot(aes(TRG_score, log10(abundance + 1.5e-05)))  +
  geom_jitter(aes(color = TRG_score)) +
  geom_smooth(method = "lm") +
  scale_color_gradient() +
  facet_wrap(~Species, 
             nrow = 5) +
  stat_cor(method = "spearman") +
  theme_minimal() +
  theme(legend.position = "none",
        aspect.ratio = 1.2)



# TRG_2 범주에 따른 Wilcoxon 검정
trg2_results <- merged_s %>%
  filter(TRG_2 %in% c("nonCR", "CR")) %>%
  group_by(Species) %>%
  summarise(p_value = wilcox.test(abundance ~ TRG_2, 
                                  data = .,
                                  exact = FALSE)$p.value,
            .groups = 'drop') %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>% 
  arrange(FDR)

trg2_results %>% 
  filter(FDR < 0.1)

# TRG_3 범주에 따른 Wilcoxon 검정
trg3_results <- merged_s %>%
  filter(TRG_3 %in% c("bad", "nearCR")) %>%
  group_by(clade_name) %>%
  summarise(p_value = wilcox.test(abundance ~ TRG_3, 
                                  data = ., 
                                  exact = FALSE)$p.value,
            .groups = 'drop') %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>% 
  arrange(FDR)



#-------------# 
t2 <- t %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-clade_name) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence > 8)


merged_t <- t %>%
  filter(clade_name %in% t2$clade_name) %>% 
  pivot_longer(cols = -clade_name, 
               names_to = "SampleID_2",
               values_to = "abundance") %>%
  left_join(m %>% select(SampleID_2, TRG_2, TRG_3), by = "SampleID_2")


# TRG_2 범주에 따른 Wilcoxon 검정
trg2_results <- merged_t %>%
  filter(TRG_2 %in% c("nonCR", "CR")) %>%
  group_by(clade_name) %>%
  summarise(p_value = wilcox.test(abundance ~ TRG_2, 
                                  data = .,
                                  exact = FALSE)$p.value,
            .groups = 'drop') %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>% 
  arrange(FDR)

s %>% 
  filter(clade_name == 
           "k__Bacteria|p__Actinobacteria|c__Actinomycetia|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_oris")

t %>% 
  filter(clade_name == 
           "k__Bacteria|p__Actinobacteria|c__Actinomycetia|o__Actinomycetales|f__Actinomycetaceae|g__Actinomyces|s__Actinomyces_oris|t__SGB15878")

# TRG_3 범주에 따른 Wilcoxon 검정
trg3_results <- merged_t %>%
  filter(TRG_3 %in% c("bad", "nearCR")) %>%
  group_by(clade_name) %>%
  summarise(p_value = wilcox.test(abundance ~ TRG_3, 
                                  data = ., 
                                  exact = FALSE)$p.value,
            .groups = 'drop') %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>% 
  arrange(FDR)



#--------------- 3. Differentially-enriched pathways ---------------# 
