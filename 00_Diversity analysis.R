
# Sample Filter: TNT - Before

rm(list = ls())
options(java.parameters = "-Xmx64g", stringsAsFactors = F)
setwd("D:/2-연구/2-CRC metagenomics")


library(vegan)
library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggforce)





###############  Import  ###############

# load("260224 final Input file/RData/7-1. after-diversity-analysis.RData")

m = readxl::read_xlsx("260224 final Input file/metadata_최종.xlsx") %>% 
  mutate(SampleID = paste0("Sample_", SampleID),
         TRG = ifelse(is.na(TRG), "nearCR", TRG),
         Pre_Op_Tstage_bin = ifelse(Pre_Op_Tstage >= 3, "1", "0"),
         Pre_Op_Nstage_bin = ifelse(Pre_Op_Nstage >= 1, "1", "0"),
         Age_bin = ifelse(Age >= 60, "1", "0"),
         CEA = ifelse(CEA == "<0.5", 0.5, CEA),
         CEA = as.numeric(CEA),
         CEA_bin = ifelse(CEA < 5.0, "0", "1"),
         BMI_bin = ifelse(BMI > 25.0, "1", "0")) %>% 
  filter(!is.na(TRG_score))

table(m$TRG_score)
table(m$cTstage)
table(m$cNstage)
table(m$CEA_bin)
# 42 Samples (CR 18, nearCR 5, PR 13, Poor 6)
# T-stage(Tumor stage): 종양의 크기, 퍼진 정도 (2: 6, 3: 25, 4: 11)
# N-stage(Lymph Node stage): 림프절에 암이 전이된 정도
# cTstage (2:cT2 6, 3:cT3 25, 4:cT4a 8, 5:cT4b 3)


# CR: "#9FD0E4"
# nonCR: "#F2A7A0"
# 
# CR: "#7FCDBB"
# nonCR: "#F29C99"
# 
# CR: "#5FB8A5"
# nonCR: "#E98580"
# 
# CR: "#4FAE9A"
# nonCR: "#DE7872"
# 
# CR: "#63B3A3"
# nonCR: "#D97C78"
# 
# CR: "#479C8C"
# nonCR: "#CF6B66"
# 
# CR: "#6CC3B0"
# nonCR: "#EB8D89"
# 
# CR: "#AFC6E9"
# nonCR: "#E8A2A8"



# Abundance table 
f <- read.csv("260224 final Input file/Abundance file/family.csv") %>% 
  select(c("clade_name", m$SampleID)) %>% 
  mutate(Family = str_extract(clade_name, "(?<=f__).*")) %>% 
  select(-clade_name) ; head(f)

g <- read.csv("260224 final Input file/Abundance file/genus.csv") %>% 
  select(c("clade_name", m$SampleID)) %>% 
  mutate(Genus = str_extract(clade_name, "(?<=g__).*")) %>% 
  select(-clade_name) ; head(g)

s <- read.csv("260224 final Input file/Abundance file/species.csv") %>% 
  select(c("clade_name", m$SampleID)) %>% 
  mutate(Species = str_extract(clade_name, "(?<=s__).*")) %>% 
  select(-clade_name) ; head(s)
  s[1:3, ]

t <- read.csv("260224 final Input file/Abundance file/strain.csv") %>% 
  select(c("clade_name", m$SampleID)) %>% 
  mutate(Strain = str_extract(clade_name, "(?<=s__).*")) %>%
  mutate(Strain = gsub("t__", "", Strain)) %>% 
  select(-clade_name) ; head(t)

### Calculate Alpha-diversity indice 
m = data.frame(
  ObservedStrain = colSums(t[, -ncol(t)] > 0),
  Shannon = diversity(t(t[, -ncol(t)]), index = "shannon"),
  InvSimpson = diversity(t(t[, -ncol(t)]), index = "invsimpson")) %>%
  rownames_to_column("SampleID") %>%
  merge(m, ., by = "SampleID")


# m = data.frame(
#   ObservedStrain = colSums(s[, -ncol(s)] > 0),
#   Shannon = diversity(t(s[, -ncol(s)]), index = "shannon"),
#   InvSimpson = diversity(t(s[, -ncol(s)]), index = "invsimpson")) %>% 
#   rownames_to_column("SampleID") %>% 
#   merge(m, ., by = "SampleID")


# Sorting TNT before samples  
mb = m %>% filter(TNT == "Before")
mo = m %>% filter(TNT == "Ongoing")

table(mb$Pre_Op_Tstage)
table(mo$Pre_Op_Tstage)

table(mb$TRG_score, mb$Pre_Op_Tstage)
table(mo$TRG_score, mo$Pre_Op_Tstage)

# Sorting "Before" samples  
fb = f %>% select(c("Family", mb$SampleID))
gb = g %>% select(c("Genus", mb$SampleID))
sb = s %>% select(c("Species", mb$SampleID))
tb = t %>% select(c("Strain", mb$SampleID))


# Clinical indicators p-value
# Continuous
wilcox.test(Age ~ TRG_1, mb)
wilcox.test(CEA ~ TRG_1, mb)
wilcox.test(CEA ~ TRG_1, mb)

# Categorical
fisher.test(table(mb$Age_bin, mb$TRG_1))
fisher.test(table(mb$Sex, mb$TRG_1))
fisher.test(table(mb$BMI_bin, mb$TRG_1))
fisher.test(table(mb$CEA_bin, mb$TRG_1))
fisher.test(table(mb$cNstage, mb$TRG_1))


chisq.test(table(mb$ASA, mb$TRG_1))
chisq.test(table(mb$cTstage, mb$TRG_1))



# Sorting TNT before samples

###############  Diversity analysis - Strain  ###############

# At the strain levels


###############  1. Before  ###############

### Beta-diversity: PCoA
t_dist = tb[, -1] %>% 
  t() %>% 
  as.matrix() %>% 
  vegdist()

t_ord = cmdscale(t_dist, k = 10, eig = T)

mb2 = t_ord$points %>% 
  data.frame() %>% 
  rownames_to_column("SampleID") %>% 
  merge(., mb, by = "SampleID") %>% 
  arrange(SubjectID, TNT) %>% 
  mutate(SampleID = factor(SampleID, levels = m$SampleID)) %>% 
  arrange(SampleID) ; head(mb2)

# Variance explained
var_exp1 = round(t_ord$eig/sum(t_ord$eig) * 100, 1)[1] ; var_exp1  # 14
var_exp2 = round(t_ord$eig/sum(t_ord$eig) * 100, 1)[2] ; var_exp2  # 10.5


# 1. 고정된 feature 변수: 거리 행렬 구성용
features = c("X1", "X2")

# 2. 그룹 변수 리스트 (모두 binary)
group_vars = c("Age_bin", "Sex", "Pre_Op_Tstage_bin", "Pre_Op_Nstage_bin", 
               "TRG_1", "TRG_2", "TRG_3", "CEA_bin", "BMI_bin")


# 거리 행렬 고정
temp_data = mb2 %>% 
  select(all_of(features)) %>%
  drop_na()

t_dist = dist(temp_data)  # 또는 vegdist() for Bray-Curtis

# 결과 저장용
res_group = data.frame(); set.seed(123) ; for (gvar in group_vars) {
  
  # 그룹 정보 추출
  group_data = mb2 %>%
    select(Group = !!sym(gvar)) %>%
    drop_na()
  
  # adonis2 수행
  adonis_res = adonis2(t_dist ~ Group, data = group_data, permutations = 999)
  adonis_p = adonis_res$`Pr(>F)`[1]
  
  # 결과 저장
  res_group = rbind(res_group,
                    data.frame(Group_Variable = gvar,
                               p_adonis = adonis_p))
} ; res_group

# Group_Variable         p_adonis
# 1           Age_bin    0.800
# 2               Sex    0.514
# 3 Pre_Op_Tstage_bin    0.920
# 4 Pre_Op_Nstage_bin    0.422
# 5             TRG_1    0.333
# 6             TRG_2    0.159*
# 7             TRG_3    0.989
# 8           CEA_bin    0.423
# 9           BMI_bin    0.099*


# Wilcoxon test & Welch's t-test 
res_t = data.frame(Variable = character(),
                   Group = character(),
                   p_wilcox = numeric(),
                   p_ttest = numeric(),
                   stringsAsFactors = F)

# 그룹 변수와 피처 변수에 대해 반복
set.seed(123); for (gvar in group_vars) {
  for (var in features) {
    
    # 테스트용 데이터 구성
    ptest = mb2 %>%
      select(Y = !!sym(var), Group = !!sym(gvar)) %>%
      drop_na()
    
    # Wilcoxon test
    wilcox_res = wilcox.test(Y ~ Group, data = ptest)
    wilcox_p = wilcox_res$p.value
    
    # Welch’s t-test
    ttest_res = t.test(Y ~ Group, data = ptest)
    ttest_p = ttest_res$p.value
    
    # 결과 저장
    res_t = rbind(res_t,
                  data.frame(Variable = var,
                             Group = gvar,
                             p_wilcox = wilcox_p,
                             p_ttest = ttest_p))
  }
};res_t

# Variable              Group   p_wilcox     p_ttest
# 1        X1           Age_bin 0.63401962   0.68614245
# 2        X2           Age_bin 0.63401962   0.62004377
# 3        X1               Sex 0.89712646   0.72219774
# 4        X2               Sex 0.21992588   0.21712243
# 5        X1 Pre_Op_Tstage_bin 0.91826087   0.74616453
# 6        X2 Pre_Op_Tstage_bin 0.81070234   0.74825234
# 7        X1 Pre_Op_Nstage_bin 0.09461538*  0.08406584*
# 8        X2 Pre_Op_Nstage_bin 0.88000000   0.88006281
# 9        X1             TRG_1 0.68321676   0.26445561
# 10       X2             TRG_1 0.18043892*  0.23999161
# 11       X1             TRG_2 0.70450666   0.26178111
# 12       X2             TRG_2 0.08492705*  0.10354959*
# 13       X1             TRG_3 1.00000000   0.89272087
# 14       X2             TRG_3 0.91826087   0.92756209
# 15       X1           CEA_bin 0.28307692   0.26547560
# 16       X2           CEA_bin 0.51451505   0.51180570
# 17       X1           BMI_bin 0.95500152   0.30856192
# 18       X2           BMI_bin 0.03516570** 0.05025731**



# ##### Faith's phylogenetic diversity
# # -->  No difference between CR vs. nonCR --> not going to use this
# 
# tree <- read.tree("260224 final Input file/mpa_vJan21_CHOCOPhlAnSGB_202103.nwk")
#   length(tree$tip.label)
#   head(tree$tip.label, 20)
#   tail(tree$tip.label, 20)
# 
# mpa <-  read.csv("260224 final Input file/Abundance file/strain.csv",
#                  check.names = F, comment.char = "#")
# 
# tax_col <- colnames(mpa)[1]
# 
# 
# # SGB row만 남기기
# mpa_sgb <- mpa %>%
#   dplyr::filter(grepl("t__SGB", .data[[tax_col]]))
# 
# # 첫 컬럼에서 숫자 ID만 추출
# mpa_sgb[[tax_col]] <- gsub(".*t__SGB", "", mpa_sgb[[tax_col]])
# 
# rownames(mpa_sgb) <- mpa_sgb[[tax_col]]
# mpa_sgb <- mpa_sgb[, -1, drop = FALSE]
# 
# # 공통 taxon 확인
# common_taxa <- intersect(rownames(mpa_sgb), tree$tip.label)
# 
# length(common_taxa)
# head(common_taxa)
# 
# # abundance와 tree 맞추기
# mpa_sgb2 <- mpa_sgb[common_taxa, , drop = FALSE]
# tree2 <- ape::keep.tip(tree, common_taxa)
# 
# # tree 순서에 맞춤
# mpa_sgb2 <- mpa_sgb2[tree2$tip.label, , drop = FALSE]
# 
# # presence/absence
# otu_pa <- t(mpa_sgb2 > 0)
# 
# # Faith's PD
# faith_df <- picante::pd(otu_pa, tree2, include.root = TRUE) %>%
#   rownames_to_column("SampleID") %>%
#   rename(FaithPD = PD, FaithSR = SR)
# 
# 
# m <- m %>%
#   left_join(faith_df, by = "SampleID")
# 
# mb <- m %>% filter(TNT == "Before")
# 
# 
# p_faith <- mb %>%
#   ggplot(aes(TRG_1, FaithPD)) +
#   geom_jitter(aes(color = TRG_1), width = 0.12, size = 2) +
#   geom_boxplot(aes(fill = TRG_1), alpha = 0.3, outlier.alpha = 0, width = 0.6) +
#   stat_compare_means(comparisons = list(c("CR", "nonCR")),
#                      method = "wilcox", tip.length = 0.02) +
#   scale_fill_manual(values = c("CR" = "#9FD0E4",
#                                "nonCR" = "#F2A7A0")) +
#   scale_color_manual(values = c("CR" = "#9FD0E4",
#                                 "nonCR" = "#F2A7A0")) +
#   theme_pubr() +
#   theme(legend.position = "none") +
#   labs(x = "Tumor Regression Grade",
#        y = "Faith's phylogenetic diversity")



# Alpha diversity: Shannon diversity by TRG score (complete, near, partial, poor)
p4_shannon_TRG = mb2 %>%
  mutate(TRG = factor(TRG, levels = c("CR", "nearCR", "PR", "Poor"))) %>%
  ggplot(aes(TRG, Shannon)) +
  geom_jitter(aes(color = TRG)) +
  geom_boxplot(aes(fill = TRG),
               alpha = 0.3, outlier.alpha = 0, width = 0.6) +
  stat_compare_means(comparisons = list(c("CR", "nearCR"),
                                        c("CR", "PR"),
                                        c("CR", "Poor")),
                     tip.length = 0.02, method = "wilcox") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nearCR" = "#FFB343",
                               "PR" = "#FF7518",
                               "Poor" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nearCR" = "#FFB343",
                                "PR" = "#FF7518",
                                "Poor" = "#80461B")) +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1.1)),
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Tumor Regression Grade",
       y = "Shannon diversity") +
  coord_cartesian(ylim = c(2, 5.5)) ; p4_shannon_TRG

# ggsave("figure/04-1-1_Shannon-strain.svg",
#        plot = p4_shannon_TRG, width = 3, height = 5)

# TRG_1
p4_shannon_TRG_cat = mb2 %>% 
  ggplot(aes(TRG_1, Shannon)) +
  geom_jitter(aes(color = TRG_1)) +
  geom_boxplot(aes(fill = TRG_1), 
               alpha = 0.3, outlier.alpha = 0,
               width = 0.6) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     tip.length = 0.02,
                     method = "wilcox") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Tumor Regression Grade",
       y = "Shannon diversity") +
  coord_cartesian(ylim = c(2, 5)) ; p4_shannon_TRG_cat

# ggsave("figure/04-1-1_Shannon-strain_cat.svg",
#        plot = p4_shannon_TRG_cat, width = 1.8, height = 5)


# Alpha diversity: Observed strain
p4_observed_TRG = mb2 %>% 
  mutate(TRG = factor(TRG, levels = c("CR", "nearCR", "PR", "Poor"))) %>% 
  ggplot(aes(TRG, ObservedStrain)) +
  geom_boxplot(aes(fill = TRG),
               alpha = 0.3, outlier.alpha = 0, width = 0.6) +
  geom_jitter(aes(color = TRG)) +
  stat_compare_means(comparisons = list(c("CR", "nearCR"),
                                        c("CR", "PR"),
                                        c("CR", "Poor")),
                     tip.length = 0.02, method = "wilcox") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nearCR" = "#FFB343",
                               "PR" = "#FF7518",
                               "Poor" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nearCR" = "#FFB343",
                                "PR" = "#FF7518",
                                "Poor" = "#80461B")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Tumor Regression Grade", y = "Observed Strain #") ; p4_observed_TRG

# ggsave("figure/04-1-2_Observed-strain.svg",
#        plot = p4_observed_TRG, width = 3, height = 5)

# TRG_1
p4_observed_TRG_cat = mb2 %>% 
  ggplot(aes(TRG_1, ObservedStrain)) +
  geom_boxplot(aes(fill = TRG_1),
               alpha = 0.3, outlier.alpha = 0, width = 0.6) +
  geom_jitter(aes(color = TRG_1)) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     method = "wilcox", tip.length = 0.02) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Tumor Regression Grade", y = "Observed Strain #") +
  coord_cartesian(ylim = c(50, 500)) ; p4_observed_TRG_cat

# ggsave("figure/04-1-2_Observed-strain_cat.svg",
#        plot = p4_observed_TRG_cat, width = 1.8, height = 5)


# Inverted Simpson
p4_invsimpson_TRG = mb2 %>% 
  mutate(TRG = factor(TRG, levels = c("CR", "nearCR", "PR", "Poor"))) %>% 
  ggplot(aes(TRG, InvSimpson)) +
  geom_boxplot(aes(fill = TRG),
               alpha = 0.3, outlier.alpha = 0, width = 0.6) +
  geom_jitter(aes(color = TRG)) +
  stat_compare_means(comparisons = list(c("CR", "nearCR"),
                                        c("CR", "PR"),
                                        c("CR", "Poor")),
                     tip.length = 0.02, method = "wilcox") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nearCR" = "#FFB343",
                               "PR" = "#FF7518",
                               "Poor" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nearCR" = "#FFB343",
                                "PR" = "#FF7518",
                                "Poor" = "#80461B")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Tumor Regression Grade", y = "Inverted Simpson") ; p4_invsimpson_TRG

# ggsave("figure/04-1-3_Inverted Simpson-strain.svg",
#        plot = p4_invsimpson_TRG, width = 3, height = 5)

# TRG_1
p4_invsimpson_TRG_cat = mb2 %>% 
  ggplot(aes(TRG_1, InvSimpson)) +
  geom_boxplot(aes(fill = TRG_1),
               alpha = 0.3, outlier.alpha = 0, width = 0.6) +
  geom_jitter(aes(color = TRG_1)) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     tip.length = 0.02, method = "wilcox") +
  # scale_fill_manual(values = c("CR" = "#F7D9BC",
  #                              "nonCR" = "#80461B")) +
  # scale_color_manual(values = c("CR" = "#F7D9BC",
  #                               "nonCR" = "#80461B")) +
  scale_fill_manual(values = c("CR" = "#4FAE9A",
                               "nonCR" = "#DE7872")) +
  scale_color_manual(values = c("CR" = "#4FAE9A",
                                "nonCR" = "#DE7872")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Tumor Regression Grade", y = "Inverted Simpson") ; p4_invsimpson_TRG_cat

# ggsave("figure/04-1-3_Inverted Simpson-strain_cat.svg",
#        plot = p4_invsimpson_TRG_cat, width = 2, height = 5)


# Pielou's evenness
even_t = t(tb[, -1])

H = diversity(even_t, index = "shannon")
S = specnumber(even_t)

evenness_t = H / log(S)
even_df_t = data.frame(
  SampleID = rownames(even_t),
  Evenness = evenness_t
)

mb2 = even_df_t %>% left_join(mb2, by = "SampleID"); rm(even_t, even_df_t, evenness_t, H, S)

p4_evenness_TRG = mb2 %>% 
  mutate(TRG = factor(TRG, levels = c("CR", "nearCR", "PR", "Poor"))) %>% 
  ggplot(aes(TRG, Evenness)) +
  geom_boxplot(aes(fill = TRG),
               alpha = 0.3, outlier.alpha = 0, width = 0.6) +
  geom_jitter(aes(color = TRG)) +
  stat_compare_means(comparisons = list(c("CR", "nearCR"),
                                        c("CR", "PR"),
                                        c("CR", "Poor")),
                     tip.length = 0.02, method = "wilcox") +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nearCR" = "#FFB343",
                               "PR" = "#FF7518",
                               "Poor" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nearCR" = "#FFB343",
                                "PR" = "#FF7518",
                                "Poor" = "#80461B")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Tumor Regression Grade", y = "Pielou's evenness") ; p4_evenness_TRG

# ggsave("figure/04-1-4_Evenness-strain.svg",
#        plot = p4_evenness_TRG, width = 3, height = 5)

# TRG_1
# CR: "#63B3A3"
# nonCR: "#D97C78"
# 
# CR: "#4FAE9A"
# nonCR: "#DE7872"

p4_evenness_TRG_cat = mb2 %>% 
  ggplot(aes(TRG_1, Evenness)) +
  geom_boxplot(aes(fill = TRG_1),
               alpha = 0.5, outlier.alpha = 0, width = 0.6) +
  geom_jitter(aes(color = TRG_1)) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     tip.length = 0.02, method = "wilcox") +
  # scale_fill_manual(values = c("CR" = "#F7D9BC",
  #                              "nonCR" = "#80461B")) +
  # scale_color_manual(values = c("CR" = "#F7D9BC",
  #                               "nonCR" = "#80461B")) +
  scale_fill_manual(values = c("CR" = "#4FAE9A",
                               "nonCR" = "#DE7872")) +
  scale_color_manual(values = c("CR" = "#4FAE9A",
                                "nonCR" = "#DE7872")) +
  theme_pubr() +
  theme(legend.position = "none", 
        axis.title = element_text(size = rel(1.1)), 
        axis.text = element_text(size = rel(0.95))) +
  labs(x = "Tumor Regression Grade", y = "Pielou's evenness") +
  scale_y_continuous(breaks = c(0.6, 0.7, 0.8)) ; p4_evenness_TRG_cat

# ggsave("figure/04-1-4_Evenness-strain_cat.svg",
#        plot = p4_evenness_TRG_cat, width = 2, height = 5)






##### Violin + boxplot
library(ggbeeswarm)

## 1) Evenness
med_CR <- median(plot_df$Evenness[plot_df$TRG_1 == "CR"], na.rm = TRUE)
med_nonCR <- median(plot_df$Evenness[plot_df$TRG_1 == "nonCR"], na.rm = TRUE)
delta_med <- round(med_nonCR - med_CR, 3)

plot_df <- mb2 %>%
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR"))) %>%
  select(TRG_1, Evenness) %>%
  na.omit()

y_max <- max(plot_df$Evenness, na.rm = TRUE)
y_bracket <- y_max + 0.015
y_text <- y_max + 0.022

p4_evenness_final <- ggplot(
  plot_df,
  aes(x = TRG_1, y = Evenness, fill = TRG_1, color = TRG_1)
) +
  
  # violin
  geom_violin(
    width = 0.8,
    trim = FALSE,
    alpha = 0.35,   # 이전보다 조금 더 진하게
    color = NA
  ) +
  
  # boxplot: fill 없음, whisker 제거
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    fill = NA,
    color = "black",
    linewidth = 0.5,
    coef = 0
  ) +

  # raw points
  ggbeeswarm::geom_quasirandom(
    width = 0.10,
    size = 2.3,
    alpha = 0.9
  ) +
  
  # p-value bracket
  annotate("segment", x = 1, xend = 2, y = y_bracket, yend = y_bracket, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1, y = y_bracket - 0.004, yend = y_bracket, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2, y = y_bracket - 0.004, yend = y_bracket, linewidth = 0.4) +
  annotate("text", x = 1.5, y = y_text, label = "P = 0.047", size = 3.8) +
  
  scale_fill_manual(values = c("CR" = "#4FAE9A",
                               "nonCR" = "#DE7872")) +
  scale_color_manual(values = c("CR" = "#4FAE9A",
                                "nonCR" = "#DE7872")) +
  
  scale_y_continuous(
    breaks = c(0.6, 0.7, 0.8),
    limits = c(min(plot_df$Evenness, na.rm = TRUE) - 0.02, y_text + 0.01)
  ) +
  
  labs(
    x = "Tumor Regression Grade",
    y = "Pielou's evenness"
  ) +
  
  theme_pubr() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    axis.text.y = element_text(size = rel(0.9))
    ); p4_evenness_final


ggsave("figures/04-1-4_Evenness-strain_cat_final.svg",
       plot = p4_evenness_final, width = 2.5, height = 4.5)



## 2) Inverted Simpson
med_CR <- median(plot_df2$InvSimpson[plot_df2$TRG_1 == "CR"], na.rm = TRUE)
med_nonCR <- median(plot_df2$InvSimpson[plot_df2$TRG_1 == "nonCR"], na.rm = TRUE)
delta_med <- round(med_nonCR - med_CR, 3)


plot_df2 <- mb2 %>%
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR"))) %>%
  select(TRG_1, InvSimpson) %>%
  na.omit()

y_max2 <- max(plot_df2$InvSimpson, na.rm = TRUE)
y_bracket2 <- y_max2 + 0.25
y_text2 <- y_max2 + 0.38

p4_invsimpson_final <- ggplot(
  plot_df2,
  aes(x = TRG_1, y = InvSimpson, fill = TRG_1, color = TRG_1)
) +
  
  # violin
  geom_violin(
    width = 0.8,
    trim = FALSE,
    alpha = 0.35,
    color = NA
  ) +
  
  # boxplot: fill 없음, whisker 제거
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    fill = NA,
    color = "black",
    linewidth = 0.5,
    coef = 0
  ) +

  # raw points
  ggbeeswarm::geom_quasirandom(
    width = 0.10,
    size = 2.3,
    alpha = 0.9
  ) +
  
  # p-value bracket
  annotate("segment", x = 1, xend = 2, y = y_bracket2, yend = y_bracket2, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1, y = y_bracket2 - 0.08, yend = y_bracket2, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2, y = y_bracket2 - 0.08, yend = y_bracket2, linewidth = 0.4) +
  annotate("text", x = 1.5, y = y_text2, label = "P = 0.087", size = 3.8) +
  
  scale_fill_manual(values = c("CR" = "#4FAE9A",
                               "nonCR" = "#DE7872")) +
  scale_color_manual(values = c("CR" = "#4FAE9A",
                                "nonCR" = "#DE7872")) +
  
  scale_y_continuous(
    limits = c(min(plot_df2$InvSimpson, na.rm = TRUE) - 0.3, y_text2 + 0.15)
  ) +
  
  labs(
    x = "Tumor Regression Grade",
    y = "Inverse Simpson index"
  ) +
  
  theme_pubr() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    axis.text.y = element_text(size = rel(0.9))
  )

p4_invsimpson_final

ggsave("figures/04-1-5_InvSimpson-strain_cat_final.svg",
       plot = p4_invsimpson_final, width = 2.5, height = 4.5)





# ######## Gardner–Altman plot: evaluate 95% CI and effect size
# # install.packages("dabestr")
# library(dabestr)
# 
# evenness_df <- mb2 %>%
#   mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR"))) %>%
#   select(TRG_1, Evenness) %>%
#   na.omit()
# 
# # 1) dabest object 생성
# dabest_obj <- dabestr::load(
#   data = evenness_df,
#   x = TRG_1,
#   y = Evenness,
#   idx = c("CR", "nonCR")
# )
# 
# # 2) effect size 계산
# dabest_mean <- dabestr::mean_diff(dabest_obj)
# 
# # 3) Gardner–Altman plot
# p_evenness_GA <- dabestr::dabest_plot(
#   dabest_mean,
#   float_contrast = TRUE
# )
# 
# p_evenness_GA
# 
# # 4) Median difference
# dabest_median <- dabestr::median_diff(dabest_obj)
# 
# p_evenness_GA_median <- dabestr::dabest_plot(
#   dabest_median,
#   float_contrast = TRUE
# )
# 
# p_evenness_GA_median






############# PCoA - Beta diversity #############
# ============================================================
# 0. Packages
# ============================================================

# ============================================================
# 0. Packages
# ============================================================
library(ggbeeswarm)
library(cowplot)
library(scales)

# ============================================================
# 1. Prepare Before-only species matrix
# ============================================================
prepare_species_before <- function(df,
                                   feature_col = "Species",
                                   meta_df,
                                   sample_col = "SampleID",
                                   group_col = "TRG_1",
                                   time_col = "Time",
                                   before_label = "Before") {
  
  if (!feature_col %in% names(df)) stop(paste0("feature_col '", feature_col, "' not found in species object."))
  if (!sample_col %in% names(meta_df)) stop(paste0("sample_col '", sample_col, "' not found in metadata."))
  if (!group_col %in% names(meta_df)) stop(paste0("group_col '", group_col, "' not found in metadata."))
  if (!time_col %in% names(meta_df)) stop(paste0("time_col '", time_col, "' not found in metadata."))
  
  meta_sub <- meta_df %>%
    filter(.data[[time_col]] == before_label,
           .data[[group_col]] %in% c("CR", "nonCR")) %>%
    distinct(.data[[sample_col]], .keep_all = TRUE)
  
  # feature x sample matrix
  mat <- df %>%
    mutate("{feature_col}" := make.unique(as.character(.data[[feature_col]]))) %>%
    tibble::column_to_rownames(feature_col) %>%
    as.matrix()
  
  mode(mat) <- "numeric"
  
  common_samples <- intersect(colnames(mat), meta_sub[[sample_col]])
  if (length(common_samples) == 0) {
    stop("No overlapping SampleIDs between species matrix and metadata.")
  }
  
  mat <- mat[, common_samples, drop = FALSE]
  
  # sample x feature
  samp_by_feat <- t(mat)
  
  # all-zero sample 제거
  keep_sample <- rowSums(samp_by_feat, na.rm = TRUE) > 0
  samp_by_feat <- samp_by_feat[keep_sample, , drop = FALSE]
  
  # all-zero feature 제거
  keep_feature <- colSums(samp_by_feat, na.rm = TRUE) > 0
  samp_by_feat <- samp_by_feat[, keep_feature, drop = FALSE]
  
  meta_sub <- meta_sub %>%
    filter(.data[[sample_col]] %in% rownames(samp_by_feat)) %>%
    arrange(match(.data[[sample_col]], rownames(samp_by_feat)))
  
  samp_by_feat <- samp_by_feat[meta_sub[[sample_col]], , drop = FALSE]
  
  list(
    samp_by_feat = samp_by_feat,
    meta = meta_sub
  )
}

# ============================================================
# 2. Bray-Curtis + PCoA + PERMANOVA + betadisper
# ============================================================
run_species_bray_pcoa <- function(samp_by_feat, meta_df,
                                  sample_col = "SampleID",
                                  group_col = "TRG_1",
                                  permutations = 999,
                                  seed = 1234) {
  
  set.seed(seed)
  
  dist_obj <- vegdist(samp_by_feat, method = "bray")
  pcoa_res <- ape::pcoa(dist_obj)
  
  coords <- as.data.frame(pcoa_res$vectors[, 1:2, drop = FALSE])
  colnames(coords) <- c("PCo1", "PCo2")
  coords[[sample_col]] <- rownames(coords)
  
  var_exp <- round(100 * pcoa_res$values$Relative_eig[1:2], 1)
  
  sample_ids <- attr(dist_obj, "Labels")
  meta_use <- meta_df %>%
    filter(.data[[sample_col]] %in% sample_ids) %>%
    distinct(.data[[sample_col]], .keep_all = TRUE) %>%
    arrange(match(.data[[sample_col]], sample_ids))
  
  adonis_res <- adonis2(
    as.formula(paste("dist_obj ~", group_col)),
    data = meta_use,
    permutations = permutations
  )
  
  bd <- betadisper(dist_obj, group = meta_use[[group_col]], type = "centroid")
  
  set.seed(seed)
  bd_perm <- permutest(bd, permutations = permutations)
  
  list(
    dist = dist_obj,
    pcoa = pcoa_res,
    coords = coords,
    var_exp = var_exp,
    adonis = adonis_res,
    betadisper = bd,
    betadisper_perm = bd_perm,
    meta_use = meta_use
  )
}

# ============================================================
# 3. Convex hull helper
# ============================================================
get_convex_hull_df <- function(plot_df,
                               group_col = "TRG_plot",
                               x_col = "PCo1",
                               y_col = "PCo2") {
  
  plot_df %>%
    group_by(.data[[group_col]]) %>%
    filter(n() >= 3) %>%
    slice(chull(.data[[x_col]], .data[[y_col]])) %>%
    ungroup()
}

# ============================================================
# 4. Make plot components
#    - main
#    - top strip
#    - right strip
#    - dispersion
# ============================================================
make_species_pcoa_components <- function(res, meta_df,
                                         sample_col = "SampleID",
                                         group_col = "TRG_1") {
  
  plot_df <- res$coords %>%
    left_join(meta_df, by = sample_col) %>%
    mutate(
      TRG_plot = case_when(
        .data[[group_col]] == "CR" ~ "pCR",
        .data[[group_col]] == "nonCR" ~ "non_pCR",
        TRUE ~ as.character(.data[[group_col]])
      ),
      TRG_plot = factor(TRG_plot, levels = c("pCR", "non_pCR")),
      top_y = ifelse(TRG_plot == "pCR", 0.492, 0.508),
      right_x = ifelse(TRG_plot == "pCR", 0.492, 0.508)
    )
  
  hull_df <- get_convex_hull_df(plot_df, group_col = "TRG_plot",
                                x_col = "PCo1", y_col = "PCo2")
  
  xlim_main <- range(plot_df$PCo1, na.rm = TRUE)
  ylim_main <- range(plot_df$PCo2, na.rm = TRUE)
  
  p_pco1 <- wilcox.test(PCo1 ~ TRG_plot, data = plot_df, exact = FALSE)$p.value
  p_pco2 <- wilcox.test(PCo2 ~ TRG_plot, data = plot_df, exact = FALSE)$p.value
  
  med1_df <- plot_df %>%
    group_by(TRG_plot) %>%
    summarise(
      med = median(PCo1, na.rm = TRUE),
      top_y = mean(top_y),
      .groups = "drop"
    )
  
  med2_df <- plot_df %>%
    group_by(TRG_plot) %>%
    summarise(
      med = median(PCo2, na.rm = TRUE),
      right_x = mean(right_x),
      .groups = "drop"
    )
  
  R2 <- res$adonis$R2[1]
  p_perm <- res$adonis$`Pr(>F)`[1]
  
  # -----------------------------
  # Top strip
  # -----------------------------
  p_top <- ggplot(plot_df, aes(x = PCo1, y = top_y, color = TRG_plot)) +
    geom_point(
      shape = 124,
      size = 6.2,
      alpha = 0.9,
      position = position_jitter(height = 0.0025, width = 0)
    ) +
    geom_segment(
      data = med1_df,
      aes(x = med, xend = med, y = top_y - 0.012, yend = top_y + 0.012),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.85
    ) +
    annotate(
      "text",
      x = mean(xlim_main),
      y = 0.545,
      label = paste0("PCo1: P = ", formatC(p_pco1, format = "f", digits = 3)),
      size = 3.6
    ) +
    scale_x_continuous(
      limits = xlim_main,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(limits = c(0.46, 0.55)) +
    scale_color_manual(values = c("pCR" = "#4FAE9A",
                                  "non_pCR" = "#DE7872")) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(5, 8, 0, 8)
    )
  
  # -----------------------------
  # Right strip
  # -----------------------------
  p_right <- ggplot(plot_df, aes(x = right_x, y = PCo2, color = TRG_plot)) +
    geom_point(
      shape = 95,
      size = 8.3,
      alpha = 0.9,
      position = position_jitter(width = 0.0025, height = 0)
    ) +
    geom_segment(
      data = med2_df,
      aes(x = right_x - 0.012, xend = right_x + 0.012, y = med, yend = med),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.85
    ) +
    annotate(
      "text",
      x = 0.50,
      y = max(ylim_main) + 0.08 * diff(ylim_main),
      label = paste0("PCo2: P = ", formatC(p_pco2, format = "f", digits = 3)),
      size = 3.4
    ) +
    scale_x_continuous(limits = c(0.46, 0.55)) +
    scale_y_continuous(
      limits = c(min(ylim_main), max(ylim_main) + 0.12 * diff(ylim_main)),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_color_manual(values = c("pCR" = "#4FAE9A",
                                  "non_pCR" = "#DE7872")) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(8, 5, 8, 0)
    )
  
  # -----------------------------
  # Main panel
  # -----------------------------
  p_main <- ggplot(plot_df, aes(x = PCo1, y = PCo2, color = TRG_plot)) +
    geom_polygon(
      data = hull_df,
      aes(x = PCo1, y = PCo2, fill = TRG_plot, color = TRG_plot, group = TRG_plot),
      inherit.aes = FALSE,
      alpha = 0.06,
      linewidth = 0.8
    ) +
    geom_point(size = 2.9, alpha = 0.95) +
    annotate(
      "text",
      x = xlim_main[1] + 0.03 * diff(xlim_main),
      y = ylim_main[2] - 0.05 * diff(ylim_main),
      hjust = 0, vjust = 1,
      label = paste0("PERMANOVA: R² = ", sprintf("%.3f", R2),
                     ", P = ", formatC(p_perm, format = "f", digits = 3)),
      size = 4.0
    ) +
    scale_color_manual(values = c("pCR" = "#4FAE9A",
                                  "non_pCR" = "#DE7872")) +
    scale_fill_manual(values = c("pCR" = "#4FAE9A",
                                 "non_pCR" = "#DE7872")) +
    labs(
      x = paste0("PCo1 (", res$var_exp[1], "%)"),
      y = paste0("PCo2 (", res$var_exp[2], "%)")
    ) +
    # coord_fixed(ratio = 1) +
    theme_classic(base_size = 12) +
    theme(
      aspect.ratio = 1,
      legend.title = element_blank(),
      legend.position = c(0.87, 0.90),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = alpha("white", 0.75), color = NA),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    )
  
  # -----------------------------
  # Dispersion
  # -----------------------------
  sample_ids <- attr(res$dist, "Labels")
  
  meta_use <- meta_df %>%
    filter(.data[[sample_col]] %in% sample_ids) %>%
    distinct(.data[[sample_col]], .keep_all = TRUE) %>%
    arrange(match(.data[[sample_col]], sample_ids))
  
  disp_df <- meta_use %>%
    mutate(distance_to_centroid = as.numeric(res$betadisper$distances)) %>%
    mutate(
      TRG_plot = case_when(
        .data[[group_col]] == "CR" ~ "pCR",
        .data[[group_col]] == "nonCR" ~ "non_pCR",
        TRUE ~ as.character(.data[[group_col]])
      ),
      TRG_plot = factor(TRG_plot, levels = c("pCR", "non_pCR"))
    )
  
  p_wilcox <- wilcox.test(distance_to_centroid ~ TRG_plot, data = disp_df, exact = FALSE)$p.value
  p_bd <- res$betadisper_perm$tab[1, "Pr(>F)"]
  y_max <- max(disp_df$distance_to_centroid, na.rm = TRUE)
  
  p_disp <- ggplot(disp_df, aes(x = TRG_plot, y = distance_to_centroid,
                                color = TRG_plot, fill = TRG_plot)) +
    geom_boxplot(width = 0.48, alpha = 0.14, outlier.shape = NA, color = "black") +
    ggbeeswarm::geom_quasirandom(width = 0.08, size = 2.4, alpha = 0.9) +
    annotate("segment", x = 1, xend = 2, y = y_max * 1.07, yend = y_max * 1.07, linewidth = 0.4) +
    annotate("segment", x = 1, xend = 1, y = y_max * 1.03, yend = y_max * 1.07, linewidth = 0.4) +
    annotate("segment", x = 2, xend = 2, y = y_max * 1.03, yend = y_max * 1.07, linewidth = 0.4) +
    annotate(
      "text",
      x = 1.5, y = y_max * 1.12,
      label = paste0("permutest P = ", formatC(p_bd, format = "f", digits = 3),
                     "\nWilcoxon P = ", formatC(p_wilcox, format = "f", digits = 3)),
      size = 3.7
    ) +
    scale_color_manual(values = c("pCR" = "#4FAE9A",
                                  "non_pCR" = "#DE7872")) +
    scale_fill_manual(values = c("pCR" = "#4FAE9A",
                                 "non_pCR" = "#DE7872")) +
    labs(x = NULL, y = "Distance to centroid") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )
  
  list(
    plot_df = plot_df,
    hull_df = hull_df,
    main = p_main,
    top = p_top,
    right = p_right,
    dispersion = p_disp
  )
}

# ============================================================
# 5. Combine with cowplot
# ============================================================
combine_pcoa_with_strips <- function(main_plot, top_plot, right_plot) {
  ggdraw() +
    draw_plot(main_plot,  x = 0.00, y = 0.00, width = 0.80, height = 0.80) +
    draw_plot(top_plot,   x = 0.00, y = 0.80, width = 0.80, height = 0.18) +
    draw_plot(right_plot, x = 0.80, y = 0.00, width = 0.18, height = 0.80)
}

# ============================================================
# 6. Run species-level analysis
# ============================================================

time_col_to_use <- if ("Time" %in% names(m)) "Time" else "TNT"

sp_before <- prepare_species_before(
  df = s,
  feature_col = "Species",
  meta_df = m,
  sample_col = "SampleID",
  group_col = "TRG_1",
  time_col = time_col_to_use,
  before_label = "Before"
)

sp_bray <- run_species_bray_pcoa(
  samp_by_feat = sp_before$samp_by_feat,
  meta_df = sp_before$meta,
  sample_col = "SampleID",
  group_col = "TRG_1",
  permutations = 999,
  seed = 1234
)

sp_comp <- make_species_pcoa_components(
  res = sp_bray,
  meta_df = sp_bray$meta_use,
  sample_col = "SampleID",
  group_col = "TRG_1"
)

# ============================================================
# 7. Outputs
# ============================================================

# 각각 따로 보기
sp_comp$main
sp_comp$top
sp_comp$right
sp_comp$dispersion

# R 안에서 임시 조합
sp_combined <- combine_pcoa_with_strips(
  main_plot = sp_comp$main,
  top_plot = sp_comp$top,
  right_plot = sp_comp$right
)

sp_combined

# ============================================================
# 8. Optional save
# ============================================================

# 메인, strip, divergence를 따로 저장하는 것이 가장 권장됩니다.
ggsave("../figures/species_bray_main.svg", device = "svg", 
       sp_comp$main, width = 4.8, height = 4.8)
ggsave("../figures/species_bray_top.svg", device = "svg", 
       sp_comp$top, width = 4.8, height = 1.0)
ggsave("../figures/species_bray_right.svg", device = "svg", 
       sp_comp$right, width = 1.0, height = 4.8)
ggsave("../figures/species_bray_dispersion.svg", device = "svg", 
       sp_comp$dispersion, width = 3.2, height = 3.6)

# 조합본 저장
ggsave("../figures/species_bray_combined.svg", device = "svg", 
       sp_combined, width = 6.0, height = 6.0)




t[1:3,]
t %>% filter(Strain == "Escherichia_coli|SGB10068")
t %>% filter(grepl("fragilis", Strain)) # Bacteroides_fragilis|SGB1855
t %>% filter(Strain == "Escherichia_coli|SGB10068")

#---------------------------#

p4_pcoa_TRG = mb2 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = TRG), size = 4, shape = 21) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nearCR" = "#FFB343",
                               "PR" = "#FF7518",
                               "Poor" = "#80461B")) +
  theme_pubr() +
  theme(aspect.ratio = 1, 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1))) +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)")) ; p4_pcoa_TRG

# ggsave("figure/04-1-5_PCoA-strain.svg",
#        plot = p4_pcoa_TRG, width = 8.5, height = 6)


# PERMANOVA result 
res_group

# TRG_1
p4_pcoa_TRG_cat = mb2 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = TRG_1), size = 4, shape = 21) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = min(mb2$X1), y = max(mb2$X2), 
           label = paste("PERMANOVA P-value =", 
                         res_group$p_adonis[5]), # TRG_1
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "Tumor Regression Grade") ; p4_pcoa_TRG_cat

# ggsave("figure/04-1-5_PCoA-strain-cat.svg",
#        plot = p4_pcoa_TRG_cat, width = 8.5, height = 6)

p4_hist_t_x1_TRG_1 = mb2 %>% 
  select(SampleID, X1, TRG_1) %>% 
  as.data.frame() %>% 
  ggplot(aes(X1, color = TRG_1)) +
  geom_histogram(aes(fill = TRG_1, color = TRG_1), 
                 position = "identity", linewidth = 1, alpha = 0.5) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  annotate("text", x = 0.45, y = 2.8, 
           label = paste("p-value =", round(res_t$p_wilcox[9], 3)), size = 4) +
  scale_x_continuous(breaks = NULL, expand = c(0,0)) +
  scale_y_continuous(breaks = NULL, expand = c(0,0)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) ; p4_hist_t_x1_TRG_1

# ggsave("figure/04-1-5_X1-histogram.svg",
#        plot = p4_hist_t_x1_TRG_1, width = 5, height = 3)

p4_hist_t_x2_TRG_1 = mb2 %>% 
  select(SampleID, X2, TRG_1) %>% 
  as.data.frame() %>% 
  ggplot(aes(X2, color = TRG_1)) +
  geom_histogram(aes(fill = TRG_1, color = TRG_1),
                 position  = "identity", linewidth = 1, alpha = 0.5) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  annotate("text", x = 0.25, y = 1.8, label = paste("p-value =", round(res_t$p_wilcox[2], 3))) +
  scale_x_continuous(breaks = NULL, expand = c(0,0)) +
  scale_y_continuous(breaks = NULL, expand = c(0,0)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  coord_flip() ; p4_hist_t_x2_TRG_1

# ggsave("figure/04-1-5_X2-histogram.svg",
#        plot = p4_hist_t_x2_TRG_1, width = 3, height = 5)



# PCoA plot - Bacteria
gb$Genus

genus_interest = c("Prevotella", 
                   "Bacteroides", 
                   "Phocaeicola", 
                   "Faecalibacterium", 
                   "Ruminococcus",
                   "Bifidobacterium",
                   "Akkermansia",
                   "Alistipes",
                   "Fusobacterium",
                   "Klebsiella",
                   "Campylobacter",
                   "Parvimonas",
                   "Porphyromonas",
                   "Gemella",
                   "Peptostreptococcus")

mb3 = gb %>% 
  filter(Genus %in% genus_interest) %>% 
  column_to_rownames("Genus") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleID") %>% 
  left_join(mb2, ., by = "SampleID")


# 새로운 binary 변수 추가
for (genus in genus_interest) {
  # 해당 genus 값의 median 계산 (NA 제외)
  med_val = median(mb3[[genus]], na.rm = TRUE)
  
  # 새로운 binary 변수 생성: median 이상이면 1, 미만이면 0
  mb3[[paste0(genus, "_bin")]] = ifelse(mb3[[genus]] > med_val, "High", "Low")
}


# 결과 저장용 데이터프레임 초기화
res_permanova = data.frame()

# 반복 수행: genus마다
set.seed(123); for (genus in genus_interest) {
  
  # binary 그룹 변수명
  bin_var = paste0(genus, "_bin")
  
  # 분석 대상 데이터 (NA 제거)
  temp_data = mb3 %>%
    select(X1, X2, Group = !!sym(bin_var)) %>%
    drop_na()
  
  # X1, X2에 대한 distance matrix (Euclidean)
  t_dist = dist(temp_data %>% select(X1, X2))
  
  # PERMANOVA 수행
  adonis_res = adonis2(t_dist ~ Group, data = temp_data, permutations = 999)
  adonis_p = adonis_res$`Pr(>F)`[1]
  
  # 결과 저장
  res_permanova = rbind(res_permanova,
                        data.frame(Genus = genus,
                                   GroupVar = bin_var,
                                   p_adonis = adonis_p))
}

# 결과 출력
res_permanova

# Genus               GroupVar p_adonis
# 1          Prevotella         Prevotella_bin    0.004**
# 2         Bacteroides        Bacteroides_bin    0.035**
# 3         Phocaeicola        Phocaeicola_bin    0.838
# 4    Faecalibacterium   Faecalibacterium_bin    0.142*
# 5        Ruminococcus       Ruminococcus_bin    0.595
# 6     Bifidobacterium    Bifidobacterium_bin    0.004**
# 7         Akkermansia        Akkermansia_bin    0.517
# 8           Alistipes          Alistipes_bin    0.879
# 9       Fusobacterium      Fusobacterium_bin    0.072*
# 10         Klebsiella         Klebsiella_bin    0.092*
# 11      Campylobacter      Campylobacter_bin    0.001**
# 12         Parvimonas         Parvimonas_bin    0.253
# 13      Porphyromonas      Porphyromonas_bin    0.333
# 14            Gemella            Gemella_bin    0.353
# 15 Peptostreptococcus Peptostreptococcus_bin    0.028**


# PCoA - Prevotella
p4_pcoa_prevotella_cat = mb3 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Prevotella_bin), size = 4, shape = 21) +
  scale_fill_manual(values = c("High" = "#00bb77",
                               "Low" = "gray90")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = 0.2, y = max(mb3$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[1]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "Prevotella") ; p4_pcoa_prevotella_cat

# ggsave("figure/04-1-6_PCoA-strain-Prevotella-cat.svg",
#        plot = p4_pcoa_prevotella_cat, width = 6.5, height = 5.5)

p4_pcoa_prevotella_con = mb3 %>%
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = log2(Prevotella+0.00001)), size = 4, shape = 21) +
  scale_fill_viridis_c() +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = 0.2, y = max(mb3$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[1]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "Prevotella") ; p4_pcoa_prevotella_con

# ggsave("figure/04-1-7_PCoA-strain-Prevotella-continuous.svg",
#        plot = p4_pcoa_prevotella_con, width = 6.5, height = 5.5)


# PCoA - Bacteroides
p4_pcoa_bacteroides_cat = mb3 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Bacteroides_bin), size = 4, shape = 21) +
  scale_fill_manual(values = c("High" = "#00bb77",
                               "Low" = "gray90")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = 0.2, y = max(mb3$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[2]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "Bacteroides") ; p4_pcoa_bacteroides_cat

# ggsave("figure/04-1-6_PCoA-strain-Bacteroides-cat.svg",
#        plot = p4_pcoa_bacteroides_cat, width = 6.5, height = 5.5)

p4_pcoa_bacteroides_con = mb3 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = log2(Bacteroides+0.00001)), size = 4, shape = 21) +
  scale_fill_viridis_c() +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = 0.2, y = max(mb3$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[2]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "Bacteroides") ; p4_pcoa_bacteroides_con

# ggsave("figure/04-1-7_PCoA-strain-Bacteroides-continuous.svg",
#        plot = p4_pcoa_bacteroides_con, width = 6.5, height = 5.5)


# PCoA plot - Phocaeicola
p4_pcoa_phocaeicola_cat = mb3 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Phocaeicola_bin), size = 4, shape = 21) +
  scale_fill_manual(values = c("High" = "#00bb77",
                               "Low" = "gray90")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = 0.2, y = max(mb3$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[3]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "Phocaeicola") ; p4_pcoa_phocaeicola_cat

# ggsave("figure/04-1-6_PCoA-strain-Phocaeicola-cat.svg",
#        plot = p4_pcoa_phocaeicola_cat, width = 6.5, height = 5.5)

p4_pcoa_phocaeicola_con = mb3 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = log2(Phocaeicola+0.00001)), size = 4, shape = 21) +
  scale_fill_viridis_c() +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = 0.2, y = max(mb3$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[3]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "Phocaeicola") ; p4_pcoa_phocaeicola_con

# ggsave("figure/04-1-7_PCoA-strain-Phocaeicola-continuous.svg",
#        plot = p4_pcoa_phocaeicola_con, width = 6.5, height = 5.5)


# PCoA plot - Faecalibacterium
p4_pcoa_faecalibacterium_cat = mb3 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Faecalibacterium_bin), size = 4, shape = 21) +
  scale_fill_manual(values = c("High" = "#00bb77",
                               "Low" = "gray90")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = 0.15, y = max(mb3$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[4]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "Faecalibacterium") ; p4_pcoa_faecalibacterium_cat

# ggsave("figure/04-1-6_PCoA-strain-Faecalibacterium-cat.svg",
#        plot = p4_pcoa_faecalibacterium_cat, width = 6.5, height = 5.5)


# PCoA plot - Bifidobacterium
p4_pcoa_bifidobacterium_cat = mb3 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Bifidobacterium_bin), size = 4, shape = 21) +
  scale_fill_manual(values = c("High" = "#00bb77",
                               "Low" = "gray90")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = 0.15, y = max(mb3$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[6]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "Bifidobacterium") ; p4_pcoa_bifidobacterium_cat

# ggsave("figure/04-1-6_PCoA-strain-Bifidobacterium-cat.svg",
#        plot = p4_pcoa_bifidobacterium_cat, width = 6.5, height = 5.5)


# PCoA plot - Fusobacterium
p4_pcoa_fusobacterium_cat = mb3 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Fusobacterium_bin), size = 4, shape = 21) +
  scale_fill_manual(values = c("High" = "#00bb77",
                               "Low" = "gray90")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = 0.15, y = max(mb3$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[9]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "Fusobacterium") ; p4_pcoa_fusobacterium_cat

# ggsave("figure/04-1-6_PCoA-strain-Fusobacterium-cat.svg",
#        plot = p4_pcoa_fusobacterium_cat, width = 6.5, height = 5.5)


# PCoA plot - Peptostreptococcus
p4_pcoa_peptostreptococcus_cat = mb3 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Peptostreptococcus_bin), size = 4, shape = 21) +
  scale_fill_manual(values = c("High" = "#00bb77",
                               "Low" = "gray90")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = 0.1, y = max(mb3$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[15]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "Peptostreptococcus") ; p4_pcoa_peptostreptococcus_cat

# ggsave("figure/04-1-6_PCoA-strain-Peptostreptococcus-cat.svg",
#        plot = p4_pcoa_peptostreptococcus_cat, width = 6.5, height = 5.5)


# PCoA plot - Campylobacter
p4_pcoa_campylobacter_cat = mb3 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Campylobacter_bin), size = 4, shape = 21) +
  scale_fill_manual(values = c("High" = "#00bb77",
                               "Low" = "gray90")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = 0.15, y = max(mb3$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[11]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "Campylobacter") ; p4_pcoa_campylobacter_cat

# ggsave("figure/04-1-6_PCoA-strain-Campylobacter-cat.svg",
#        plot = p4_pcoa_campylobacter_cat, width = 6.5, height = 5.5)


### PCoA plot - Strain (ratio)
strains_interst = c("Eubacterium_rectale|SGB4933",
                    "Roseburia_inulinivorans|SGB4940",
                    "Faecalibacterium_prausnitzii|SGB15316",
                    "Faecalibacterium_prausnitzii|SGB15317",
                    "Faecalibacterium_prausnitzii|SGB15318",
                    "Faecalibacterium_prausnitzii|SGB15322",
                    "Faecalibacterium_prausnitzii|SGB15323",
                    "Faecalibacterium_prausnitzii|SGB15332",
                    "Faecalibacterium_prausnitzii|SGB15342")

mb6 = tb %>% 
  filter(Strain %in% strains_interst) %>% 
  column_to_rownames("Strain") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleID") %>% 
  left_join(mb2, ., by = "SampleID")

mb6 = mb6 %>% 
  mutate(Faecalibacterium_prausnitzii_t = rowSums(select(., starts_with("Faecalibacterium_prausnitzii"))))

# 비율 계산
mb6 = mb6 %>% 
  mutate(Er_Ri_t = log10 ((`Eubacterium_rectale|SGB4933` + 1e-05) / (`Roseburia_inulinivorans|SGB4940` + 1e-05)),
         Fp_Ri_t = log10 ((`Faecalibacterium_prausnitzii_t` + 1e-05) / (`Roseburia_inulinivorans|SGB4940` + 1e-05)))

log_ratio = c("Er_Ri_t", "Fp_Ri_t")

# Median 기준 binary
for (log in log_ratio) {
  med_val = median(mb6[[log]], na.rm = T)
  cat(log, "median: ", med_val, "\n")
  mb6[[paste0(log, "_bin")]] = ifelse(mb6[[log]] > med_val, "High", "Low")
}

# PERMANOVA
res_permanova = data.frame()

set.seed(123) ; for (log in log_ratio) {
  bin_var = paste0(log, "_bin")
  
  temp_data = mb6 %>% 
    select(X1, X2, Group = !!sym(bin_var)) %>% 
    drop_na()
  
  t_dist = dist(temp_data %>% select(X1, X2))
  
  adonis_res = adonis2(t_dist ~ Group, data = temp_data, permutations = 999)
  adonis_p = adonis_res$`Pr(>F)`[1]
  
  res_permanova = rbind(res_permanova,
                        data.frame(Log_ratio = log,
                                   GroupVar = bin_var,
                                   p_adonis = adonis_p))
}

# Log_ratio    GroupVar p_adonis
# 1   Er_Ri_t Er_Ri_t_bin    0.035**
# 2   Fp_Ri_t Fp_Ri_t_bin    0.462


# PCoA - E. rectale / R. inulinivorans
p4_pcoa_Er_Ri = mb6 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Er_Ri_t_bin), size = 4, shape = 21) +
  scale_fill_manual(values = c("High" = "#FFB3A7",
                               "Low" = "gray90")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = min(mb6$X1), y = max(mb6$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[1]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "E. rectale-to-R. inulinivorans ratio") ; p4_pcoa_Er_Ri

# ggsave("figure/04-1-9_PCoA-strain-Er_Ri.svg",
#        plot = p4_pcoa_Er_Ri, width = 6.5, height = 5.5)


# PCoA - F. prausnitzii / R. inulinivorans
p4_pcoa_Fp_Ri = mb6 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Fp_Ri_t_bin), size = 4, shape = 21) +
  scale_fill_manual(values = c("High" = "#FFB3A7",
                               "Low" = "gray90")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = min(mb6$X1), y = max(mb6$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[2]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "F. prausnitzii-to-R. inulinivorans ratio") ; p4_pcoa_Fp_Ri

# ggsave("figure/04-1-9_PCoA-strain-Fp_Ri.svg",
#        plot = p4_pcoa_Fp_Ri, width = 7.0, height = 5.5)


# 각 변수에 대한 TRG_1의 비율
Er_Ri_t_TRG_1 = as.data.frame(table(mb6$TRG_1, mb6$Er_Ri_t_bin))
colnames(Er_Ri_t_TRG_1) = c("TRG_1", "Er_Ri_bin", "Count")

Er_Ri_t_p_value_CR = chisq.test(table(mb6$Er_Ri_t_bin[mb6$TRG_1 == "CR"]))$p.value
Er_Ri_t_p_value_nonCR = chisq.test(table(mb6$Er_Ri_t_bin[mb6$TRG_1 == "nonCR"]))$p.value

Er_Ri_t_TRG_1 = Er_Ri_t_TRG_1 %>% 
  group_by(TRG_1) %>% 
  mutate(p_value = ifelse(TRG_1 == "CR",
                          paste0("Chi-sq p = ", round(Er_Ri_t_p_value_CR, 3)),
                          paste0("Chi-sq p = ", round(Er_Ri_t_p_value_nonCR, 3))))

p4_Er_Ri_TRG_1 = Er_Ri_t_TRG_1 %>% 
  mutate(TRG_1 = factor(TRG_1, levels = rev(c("CR", "nonCR")))) %>% 
  mutate(Er_Ri_bin = factor(Er_Ri_bin, levels = rev(c("High", "Low")))) %>% 
  ggplot(aes(TRG_1, Count, fill = Er_Ri_bin)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(data = Er_Ri_t_TRG_1 %>% group_by(TRG_1) %>%
              summarise(total_count = sum(Count),
                        p_value = first(p_value)),
            aes(x = TRG_1, y = 3, label = p_value),
            inherit.aes = F, size = 5, color = "black") +
  scale_fill_manual(values = c("High" = "#FFB3A7",
                               "Low" = "gray90")) +
  theme_pubr() +
  labs(x = "TRG_1", y = "Count", 
       fill = "E. rectale-to-R. inulinivorans ratio") +
  coord_flip() ; p4_Er_Ri_TRG_1

# ggsave("figure/04-1-9_TRG_1-strain-Er_Ri (y = TRG).svg",
#        plot = p4_Er_Ri_TRG_1, width = 5, height = 3, limitsize = F)


Fp_Ri_t_TRG_1 = as.data.frame(table(mb6$TRG_1, mb6$Fp_Ri_t_bin))
colnames(Fp_Ri_t_TRG_1) = c("TRG_1", "Fp_Ri_bin", "Count")

Fp_Ri_t_p_value_CR = chisq.test(table(mb6$Fp_Ri_t_bin[mb6$TRG_1 == "CR"]))$p.value
Fp_Ri_t_p_value_nonCR = chisq.test(table(mb6$Fp_Ri_t_bin[mb6$TRG_1 == "nonCR"]))$p.value

Fp_Ri_t_TRG_1 = Fp_Ri_t_TRG_1 %>% 
  group_by(TRG_1) %>% 
  mutate(p_value = ifelse(TRG_1 == "CR",
                          paste0("Chi-sq p = ", round(Fp_Ri_t_p_value_CR, 3)),
                          paste0("Chi-sq p = ", round(Fp_Ri_t_p_value_nonCR, 3))))

p4_Fp_Ri_TRG_1 = Fp_Ri_t_TRG_1 %>% 
  mutate(TRG_1 = factor(TRG_1, levels = rev(c("CR", "nonCR")))) %>% 
  mutate(Fp_Ri_bin = factor(Fp_Ri_bin, levels = rev(c("High", "Low")))) %>% 
  ggplot(aes(TRG_1, Count, fill = Fp_Ri_bin)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(data = Fp_Ri_t_TRG_1 %>% group_by(TRG_1) %>%
              summarise(total_count = sum(Count),
                        p_value = first(p_value)),
            aes(x = TRG_1, y = 3, label = p_value),
            inherit.aes = F, size = 5, color = "black") +
  scale_fill_manual(values = c("High" = "#FFB3A7",
                               "Low" = "gray90")) +
  theme_pubr() +
  labs(x = "TRG_1", y = "Count",
       fill = "F. prausnitzii-to-R. inulinivorans ratio") +
  coord_flip() ; p4_Fp_Ri_TRG_1

# ggsave("figure/04-1-9_TRG_1-strain-Fp_Ri (y = TRG).svg",
#        plot = p4_Fp_Ri_TRG_1, width = 5, height = 3)


### PCoA plot - Species (ratio)
s_dist = sb[, -1] %>% 
  t() %>% 
  as.matrix() %>% 
  vegdist()

s_ord = cmdscale(s_dist, k = 10, eig = T)

mb4 = s_ord$points %>% 
  data.frame() %>% 
  rownames_to_column("SampleID") %>% 
  merge(., mb, by = "SampleID") %>% 
  arrange(SubjectID, TNT) %>% 
  mutate(SampleID = factor(SampleID, levels = m$SampleID)) %>% 
  arrange(SampleID) ; head(mb4)

# Variance explained
var_exp1 = round(s_ord$eig/sum(s_ord$eig) * 100, 1)[1] ; var_exp1  # 15
var_exp2 = round(s_ord$eig/sum(s_ord$eig) * 100, 1)[2] ; var_exp2  # 10.8

# 1. 고정된 feature 변수: 거리 행렬 구성용
features = c("X1", "X2")

# 2. 그룹 변수 리스트 (모두 binary)
group_vars = c("Age_bin", "Sex", "Pre_Op_Tstage_bin", "Pre_Op_Nstage_bin", 
               "TRG_1", "TRG_2", "TRG_3", "CEA_bin", "BMI_bin")


# 거리 행렬 고정
temp_data = mb4 %>% 
  select(all_of(features)) %>%
  drop_na()

s_dist = dist(temp_data)  # 또는 vegdist() for Bray-Curtis

# 결과 저장용
res_group = data.frame(); set.seed(123) ; for (svar in group_vars) {
  
  # 그룹 정보 추출
  group_data = mb4 %>%
    select(Group = !!sym(svar)) %>%
    drop_na()
  
  # adonis2 수행
  adonis_res = adonis2(s_dist ~ Group, data = group_data, permutations = 999)
  adonis_p = adonis_res$`Pr(>F)`[1]
  
  # 결과 저장
  res_group = rbind(res_group,
                    data.frame(Group_Variable = svar,
                               p_adonis = adonis_p))
} ; res_group


sb$Species

species_interest = c("Eubacterium_rectale",
                     "Roseburia_inulinivorans",
                     "Faecalibacterium_prausnitzii")

mb5 = sb %>% 
  filter(Species %in% species_interest) %>% 
  column_to_rownames("Species") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("SampleID") %>% 
  left_join(mb4, ., by = "SampleID")

# 비율 계산
mb5 = mb5 %>% 
  mutate(Er_Ri = log10 ((Eubacterium_rectale + 1e-05) / (Roseburia_inulinivorans + 1e-05)),
         Fp_Ri = log10 ((Faecalibacterium_prausnitzii + 1e-05) / (Roseburia_inulinivorans + 1e-05)),
         Er_Fp_Ri = log10 (((Eubacterium_rectale + 1e-05) + (Faecalibacterium_prausnitzii + 1e-05)) / (Roseburia_inulinivorans + 1e-05)))

log_ratio = c("Er_Ri", "Fp_Ri", "Er_Fp_Ri")


### Median 기준 binary
for (log in log_ratio) {
  med_val = median(mb5[[log]], na.rm = T)
  cat(log, "median:", med_val, "\n")
  mb5[[paste0(log, "_bin")]] = ifelse(mb5[[log]] > med_val, "High", "Low")
}

# PERMANOVA
res_permanova = data.frame()

set.seed(123); for (log in log_ratio) {
  
  # binary 그룹 변수명
  bin_var = paste0(log, "_bin")
  
  # 분석 대상 데이터 (NA 제거)
  temp_data = mb5 %>%
    select(X1, X2, Group = !!sym(bin_var)) %>%
    drop_na()
  
  # X1, X2에 대한 distance matrix (Euclidean)
  s_dist = dist(temp_data %>% select(X1, X2))
  
  # PERMANOVA 수행
  adonis_res = adonis2(s_dist ~ Group, data = temp_data, permutations = 999)
  adonis_p = adonis_res$`Pr(>F)`[1]
  
  # 결과 저장
  res_permanova = rbind(res_permanova,
                        data.frame(Log_ratio = log,
                                   GroupVar = bin_var,
                                   p_adonis = adonis_p))
}

# Log_ratio     GroupVar p_adonis
# 1     Er_Ri    Er_Ri_bin    0.039**
# 2     Fp_Ri    Fp_Ri_bin    0.383
# 3  Er_Fp_Ri Er_Fp_Ri_bin    0.631


# PCoA - E. rectale / R. inulinivorans
p4_pcoa_Er_Ri = mb5 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Er_Ri_bin), size = 4, shape = 21) +
  scale_fill_manual(values = c("High" = "#FFB3A7",
                               "Low" = "gray90")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = min(mb5$X1), y = max(mb5$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[1]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "E. rectale-to-R. inulinivorans ratio") ; p4_pcoa_Er_Ri

# ggsave("figure/04-1-9_PCoA-species-Er_Ri.svg",
#        plot = p4_pcoa_Er_Ri, width = 6.5, height = 5.5)


# PCoA - F. prausnitzii / R. inulinivorans
p4_pcoa_Fp_Ri = mb5 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Fp_Ri_bin), size = 4, shape = 21) +
  scale_fill_manual(values = c("High" = "#FFB3A7",
                               "Low" = "gray90")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = min(mb5$X1), y = max(mb5$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[2]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "F. prausnitzii-to-R. inulinivorans ratio") ; p4_pcoa_Fp_Ri

# ggsave("figure/04-1-9_PCoA-species-Fp_Ri.svg",
#        plot = p4_pcoa_Fp_Ri, width = 7.0, height = 5.5)


# PCoA - (E. rectale + F. prausnitzii) / R. inulinivorans
p4_pcoa_Er_Fp_Ri = mb5 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Er_Fp_Ri_bin), size = 4, shape = 21) +
  scale_fill_manual(values = c("High" = "#FFB3A7",
                               "Low" = "gray90")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = min(mb5$X1), y = max(mb5$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[3]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "E. rectale + F. prausnitzii-to-R. inulinivorans ratio") ; p4_pcoa_Er_Fp_Ri

# ggsave("figure/04-1-9_PCoA-species-Er_Fp_Ri.svg",
#        plot = p4_pcoa_Er_Fp_Ri, width = 8, height = 5.5)


# 각 변수에 대한 TRG_1의 비율
Er_Ri_TRG_1 = as.data.frame(table(mb5$TRG_1, mb5$Er_Ri_bin))
colnames(Er_Ri_TRG_1) = c("TRG_1", "Er_Ri_bin", "Count")

p4_Er_Ri_TRG_1 = Er_Ri_TRG_1 %>% 
  mutate(TRG_1 = factor(TRG_1, levels = rev(c("CR", "nonCR")))) %>% 
  mutate(Er_Ri_bin = factor(Er_Ri_bin, levels = rev(c("High", "Low")))) %>% 
  ggplot(aes(TRG_1, Count, fill = Er_Ri_bin)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("High" = "#FFB3A7",
                               "Low" = "gray90")) +
  theme_pubr() +
  labs(x = "TRG_1", y = "Count", fill = "E. rectale-to-R. inulinivorans ratio") +
  coord_flip() ; p4_Er_Ri_TRG_1

# ggsave("figure/04-1-9_TRG_1-species-Er_Ri (y = TRG).svg",
#        plot = p4_Er_Ri_TRG_1, width = 5, height = 3)


Fp_Ri_TRG_1 = as.data.frame(table(mb5$TRG_1, mb5$Fp_Ri_bin))
colnames(Fp_Ri_TRG_1) = c("TRG_1", "Fp_Ri_bin", "Count")

p4_Fp_Ri_TRG_1 = Fp_Ri_TRG_1 %>% 
  mutate(TRG_1 = factor(TRG_1, levels = rev(c("CR", "nonCR")))) %>% 
  mutate(Fp_Ri_bin = factor(Fp_Ri_bin, levels = rev(c("High", "Low")))) %>% 
  ggplot(aes(TRG_1, Count, fill = Fp_Ri_bin)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("High" = "#FFB3A7",
                               "Low" = "gray90")) +
  theme_pubr() +
  labs(x = "TRG_1", y = "Count",
       fill = "F. prausnitzii-to-R. inulinivorans ratio") +
  coord_flip() ; p4_Fp_Ri_TRG_1

# ggsave("figure/04-1-9_TRG_1-species-Fp_Ri (y = TRG).svg",
#        plot = p4_Fp_Ri_TRG_1, width = 5, height = 3)


Er_Fp_Ri_TRG_1 = as.data.frame(table(mb5$TRG_1, mb5$Er_Fp_Ri_bin))
colnames(Er_Fp_Ri_TRG_1) = c("TRG_1", "Er_Fp_Ri_bin", "Count")

p4_Er_Fp_Ri_TRG_1 = Er_Fp_Ri_TRG_1 %>% 
  mutate(TRG_1 = factor(TRG_1, levels = rev(c("CR", "nonCR")))) %>% 
  mutate(Er_Fp_Ri_bin = factor(Er_Fp_Ri_bin, levels = rev(c("High", "Low")))) %>% 
  ggplot(aes(TRG_1, Count, fill = Er_Fp_Ri_bin)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("High" = "#FFB3A7",
                               "Low" = "gray90")) +
  theme_pubr() +
  labs(x = "TRG_1", y = "Count",
       fill = "(E. rectale + F. prausnitzii)-to-R. inulinivorans Ratio") +
  coord_flip() ; p4_Er_Fp_Ri_TRG_1

# ggsave("figure/04-1-9_TRG_1-species-Er_Fp_Ri (y = TRG_1).svg",
#        plot = p4_Er_Fp_Ri_TRG_1, width = 5, height = 3)


### 분포 기준 binary
ggplot(mb5, aes(x = Er_Ri)) + 
  geom_histogram(fill = "blue", color = "black", alpha = 0.7) +
  labs(x = "Er_Ri", y = "Count") +
  theme_pubr() # criteria: -1

ggplot(mb5, aes(x = Fp_Ri)) + 
  geom_histogram(fill = "blue", color = "black", alpha = 0.7) +
  labs(x = "Fp_Ri", y = "Count") +
  theme_pubr() # criteria: 1.8

thresholds = list("Er_Ri" = -1, "Fp_Ri" = 1.8)
log_ratio = c("Er_Ri", "Fp_Ri")

for (log in log_ratio) {
  thresholds_value = thresholds[[log]]
  mb5[[paste0(log, "_bin2")]] = ifelse(mb5[[log]] > thresholds_value, "High", "Low")
}

# PERMANOVA
res_permanova = data.frame()

set.seed(123) ; for (log in log_ratio) {
  bin_var = paste0(log, "_bin2")
  
  # 분석 대상 데이터 (NA 제거)
  temp_data = mb5 %>%
    select(X1, X2, Group = !!sym(bin_var)) %>%
    drop_na()
  
  # X1, X2에 대한 distance matrix (Euclidean)
  s_dist = dist(temp_data %>% select(X1, X2))
  
  # PERMANOVA 수행
  adonis_res = adonis2(s_dist ~ Group, data = temp_data, permutations = 999)
  adonis_p = adonis_res$`Pr(>F)`[1]
  
  # 결과 저장
  res_permanova = rbind(res_permanova,
                        data.frame(Log_ratio = log,
                                   GroupVar = bin_var,
                                   p_adonis = adonis_p))
}

# Log_ratio   GroupVar p_adonis
# 1     Er_Ri Er_Ri_bin2    0.330
# 2     Fp_Ri Fp_Ri_bin2    0.314


# PCoA - E. rectale / R. inulinivorans
p4_pcoa_Er_Ri = mb5 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Er_Ri_bin2), size = 4, shape = 21) +
  scale_fill_manual(values = c("High" = "#FFB3A7",
                               "Low" = "gray90")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = min(mb5$X1), y = max(mb5$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[1]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "E. rectale-to-R. inulinivorans ratio") ; p4_pcoa_Er_Ri

# ggsave("figure/04-1-9_PCoA-species-Er_Ri (-1).svg",
#        plot = p4_pcoa_Er_Ri, width = 6.8, height = 5.5)


Er_Ri_TRG_1 = as.data.frame(table(mb5$TRG_1, mb5$Er_Ri_bin2))
colnames(Er_Ri_TRG_1) = c("TRG_1", "Er_Ri_bin", "Count")

p4_Er_Ri_TRG_1 = Er_Ri_TRG_1 %>% 
  mutate(TRG_1 = factor(TRG_1, levels = rev(c("CR", "nonCR")))) %>% 
  mutate(Er_Ri_bin = factor(Er_Ri_bin, levels = rev(c("High", "Low")))) %>% 
  ggplot(aes(TRG_1, Count, fill = Er_Ri_bin)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("High" = "#FFB3A7",
                               "Low" = "gray90")) +
  theme_pubr() +
  labs(x = "TRG_1", y = "Count", fill = "E. rectale-to-R. inulinivorans") +
  coord_flip() ; p4_Er_Ri_TRG_1

# ggsave("figure/04-1-9_TRG_1-species-Er_Ri (-1).svg",
#        plot = p4_Er_Ri_TRG_1, width = 5, height = 3)


# PCoA - F. prausnitzii / R. inulinivorans
p4_pcoa_Fp_Ri = mb5 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Fp_Ri_bin2), size = 4, shape = 21) +
  scale_fill_manual(values = c("High" = "#FFB3A7",
                               "Low" = "gray90")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = min(mb5$X1), y = max(mb5$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[2]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "F. prausnitzii-to-R. inulinivorans ratio") ; p4_pcoa_Fp_Ri

# ggsave("figure/04-1-9_PCoA-species-Fp_Ri (1.8).svg",
#        plot = p4_pcoa_Fp_Ri, width = 7.0, height = 5.5)


Fp_Ri_TRG_1 = as.data.frame(table(mb5$TRG_1, mb5$Fp_Ri_bin2))
colnames(Fp_Ri_TRG_1) = c("TRG_1", "Fp_Ri_bin", "Count")

p4_Fp_Ri_TRG_1 = Fp_Ri_TRG_1 %>% 
  mutate(TRG_1 = factor(TRG_1, levels = rev(c("CR", "nonCR")))) %>% 
  mutate(Fp_Ri_bin = factor(Fp_Ri_bin, levels = rev(c("High", "Low")))) %>% 
  ggplot(aes(TRG_1, Count, fill = Fp_Ri_bin)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("High" = "#FFB3A7",
                               "Low" = "gray90")) +
  theme_pubr() +
  labs(x = "TRG_1", y = "Count",
       fill = "F. prausnitzii / R. inulinivorans Ratio") +
  coord_flip() ; p4_Fp_Ri_TRG_1

# ggsave("figure/04-1-9_TRG_1-species-Fp_Ri (1.8).svg",
#        plot = p4_Fp_Ri_TRG_1, width = 5, height = 3)


# Continuous
res_permanova = data.frame()

set.seed(123) ; for (log in log_ratio) {
  
  # 분석 대상 데이터 (NA 제거)
  temp_data = mb5 %>%
    select(X1, X2, Group = !!sym(log)) %>%
    drop_na()
  
  # X1, X2에 대한 distance matrix (Euclidean)
  s_dist = dist(temp_data %>% select(X1, X2))
  
  # PERMANOVA 수행
  adonis_res = adonis2(s_dist ~ Group, data = temp_data, permutations = 999)
  adonis_p = adonis_res$`Pr(>F)`[1]
  
  # 결과 저장
  res_permanova = rbind(res_permanova,
                        data.frame(Log_ratio = log,
                                   p_adonis = adonis_p))
}

# Log_ratio p_adonis
# 1     Er_Ri    0.100*
# 2     Fp_Ri    0.357


p4_pcoa_Er_Ri = mb5 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Er_Ri), size = 4, shape = 21) +
  scale_fill_gradient2(low = "#00BFFF", mid = "gray90", high = "#FFEB3B") +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = min(mb5$X1), y = max(mb5$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[1]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "Er_Ri") ; p4_pcoa_Er_Ri

# ggsave("figure/04-1-9_PCoA-species-Er_Ri (con).svg",
#        plot = p4_pcoa_Er_Ri, width = 7.0, height = 5.5)


p4_pcoa_Fp_Ri = mb5 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Fp_Ri), size = 4, shape = 21) +
  scale_fill_gradient2(low = "#00BFFF", mid = "gray90", high = "#FFEB3B", midpoint = 1.614811) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = min(mb5$X1), y = max(mb5$X2), 
           label = paste0("PERMANOVA, P=", res_permanova$p_adonis[2]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "Fp_Ri") ; p4_pcoa_Fp_Ri

# ggsave("figure/04-1-9_PCoA-species-Fp_Ri (con).svg",
#        plot = p4_pcoa_Fp_Ri, width = 7.0, height = 5.5)



### PCoA plot - Clinical indicator
# Sex
# M 16 ; F 10

res_group
# Group_Variable p_adonis
# 1           Age_bin    0.800
# 2               Sex    0.514
# 3 Pre_Op_Tstage_bin    0.920
# 4 Pre_Op_Nstage_bin    0.422
# 5             TRG_1    0.333
# 6             TRG_2    0.159
# 7             TRG_3    0.989
# 8           CEA_bin    0.423
# 9           BMI_bin    0.099

p4_pcoa_sex_cat = mb3 %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Sex), size = 4, shape = 21) +
  scale_fill_manual(values = c("M" = "#5a7ef0",
                               "F" = "#f05a7e"),
                    labels = c("M" = "Male", "F" = "Female")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = 0.2, y = max(mb3$X2), 
           label = paste0("PERMANOVA P=", res_group$p_adonis[2]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "Sex") ; p4_pcoa_sex_cat

# ggsave("figure/04-1-8_PCoA-strain-Sex-cat.svg",
#        plot = p4_pcoa_sex_cat, width = 6.5, height = 5.5)


# Age
# Young(0) 9 ; Old(1) 17

p4_pcoa_age_cat = mb3 %>% 
  mutate(Age_bin = factor(Age_bin, levels = c("0", "1"), 
                          labels = c("<60", ">=60"))) %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Age_bin), size = 4, shape = 21) +
  scale_fill_manual(values = c(">=60" = "#A19AD3",
                               "<60" = "#FFF574")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = 0.2, y = max(mb3$X2), 
           label = paste0("PERMANOVA P=", res_group$p_adonis[1]),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "Age") ; p4_pcoa_age_cat

# ggsave("figure/04-1-8_PCoA-strain-Age-cat.svg",
#        plot = p4_pcoa_age_cat, width = 6.5, height = 5.5)


# Preoperative T stage
# Stage 2(0) 4 ; Stage 3-4(1) 22

p4_pcoa_Tstage_cat = mb3 %>% 
  mutate(Pre_Op_Tstage_bin = factor(Pre_Op_Tstage_bin,
                                    levels = c("0", "1"), 
                                    labels = c("early", "advanced"))) %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Pre_Op_Tstage_bin), size = 4, shape = 21) +
  scale_fill_manual(values = c("early" = "#1b9e77ff",
                               "advanced" = "#d95f02ff")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = 0.15, y = max(mb3$X2), 
           label = paste0("PERMANOVA P=", res_group$p_adonis[3]),
           hjust = 0, vjust = 1, size = 4, color = "black") + 
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "Preoperative T stage") ; p4_pcoa_Tstage_cat

# ggsave("figure/04-1-8_PCoA-strain-Tstage-cat.svg",
#        plot = p4_pcoa_Tstage_cat, width = 6.5, height = 5.5)


# Preoperative N stage
# Stage 0(0) 3 ; Stage 1(1) 23

p4_pcoa_Nstage_cat = mb3 %>% 
  mutate(Pre_Op_Nstage_bin = factor(Pre_Op_Nstage_bin,
                                    levels = c("0", "1"), 
                                    labels = c("early", "advanced"))) %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = Pre_Op_Nstage_bin), size = 4, shape = 21) +
  scale_fill_manual(values = c("early" = "#04c112ff",
                               "advanced" = "#ff7f50ff")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = 0.1, y = max(mb3$X2), 
           label = paste0("PERMANOVA P=", res_group$p_adonis[4]),
           hjust = 0, vjust = 1, size = 4, color = "black") + 
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "Preoperative N stage") ; p4_pcoa_Nstage_cat

# ggsave("figure/04-1-8_PCoA-strain-Nstage-cat.svg",
#        plot = p4_pcoa_Nstage_cat, width = 6.5, height = 5.5)


# CEA
# CEA < 5.0(0) 22 ; CEA >= 5.0(1) 4
p4_pcoa_CEA_cat = mb3 %>% 
  mutate(CEA_bin = factor(CEA_bin,
                          levels = c("0", "1"), 
                          labels = c("Low", "High"))) %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = CEA_bin), size = 4, shape = 21) +
  scale_fill_manual(values = c("Low" = "#81C784",
                               "High" = "#FF9800")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = 0.1, y = max(mb3$X2), 
           label = paste0("PERMANOVA P=", res_group$p_adonis[8]),
           hjust = 0, vjust = 1, size = 4, color = "black") + 
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "CEA") ; p4_pcoa_CEA_cat

# ggsave("figure/04-1-8_PCoA-strain-CEA-cat.svg",
#        plot = p4_pcoa_CEA_cat, width = 6.5, height = 5.5)


# BMI
# BMI <= 25.0(0) 19 ; BMI >25.0(1) 7
p4_pcoa_BMI_cat = mb3 %>% 
  mutate(BMI_bin = factor(BMI_bin,
                          levels = c("0", "1"), 
                          labels = c("≤ 25", "> 25"))) %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = BMI_bin), size = 4, shape = 21) +
  scale_fill_manual(values = c("≤ 25" = "#F7CDCD",
                               "> 25" = "#ED6B6B")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = 0.1, y = max(mb3$X2), 
           label = paste0("PERMANOVA P=", res_group$p_adonis[9]),
           hjust = 0, vjust = 1, size = 4, color = "black") + 
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "BMI") ; p4_pcoa_BMI_cat

# ggsave("figure/04-1-8_PCoA-strain-BMI-cat.svg",
#        plot = p4_pcoa_BMI_cat, width = 6.5, height = 5.5)


# TNT (Before, Ongoing)
# Before 26, Ongoing 16
t_dist_all = t[, -ncol(t)] %>% 
  t() %>% 
  as.matrix() %>% 
  vegdist()  # distance 계산

t_ord_all = cmdscale(t_dist_all, k = 10, eig = T)  # 고차원 -> 저차원

m2 = t_ord_all$points %>% 
  data.frame() %>% 
  rownames_to_column("SampleID") %>% 
  merge(., m, by = "SampleID") %>% 
  arrange(SubjectID, TNT) %>% 
  mutate(SampleID = factor(SampleID)) %>% 
  arrange(SampleID) ; head(m2)

# Variance explained
var_exp1 = round(t_ord_all$eig/sum(t_ord_all$eig) * 100, 1)[1] ; var_exp1  # 11.3
var_exp2 = round(t_ord_all$eig/sum(t_ord_all$eig) * 100, 1)[2] ; var_exp2  # 9.1

# Many different combinations
# PCo1 to PCo10 ~ TRG categories
vars = c("X1", "X2") ; groups = paste0("TRG_", 1:3)

# Make a result container
res_t = data.frame(Variable = character(),
                   Group = character(),
                   p_adonis = numeric(),
                   p_wilcox = numeric(),
                   p_ttest = numeric(),
                   stringsAsFactors = F)


# 각 변수와 그룹에 대해 반복 수행
set.seed(123) ; for(var in vars) {
  
  # Make temporary dataset
  temp_data = m2 %>% 
    select(X1, X2, Group = TNT) %>% 
    na.omit()
  
  
  # adonis2 검정 (PERMANOVA)
  
  adonis_res = adonis2(t_dist_all ~ Group,
                       data = temp_data,
                       permutations = 999)
  
  adonis_p = adonis_res$`Pr(>F)`[1]
  
  
  ptest = m2 %>% 
    select(Y = !!sym(var), Group = TRG_1) %>% 
    na.omit()
  
  
  # Wilcoxon Rank sum test
  wilcox_res = wilcox.test(ptest$Y ~ ptest$Group)
  wilcox_p = wilcox_res$p.value
  
  
  # t-test
  ttest_res = t.test(ptest$Y ~ ptest$Group)
  ttest_p = ttest_res$p.value
  
  
  # 결과 데이터프레임에 추가
  res_t = rbind(res_t,
                data.frame(Variable = var, 
                           Group = "TNT",
                           p_adonis = adonis_p,
                           p_wilcox = wilcox_p,
                           p_ttest = ttest_p))
} ; res_t

p4_pcoa_TNT_cat = m2 %>% 
  mutate(TNT = factor(TNT)) %>% 
  ggplot(aes(X1, X2)) +
  geom_point(aes(fill = TNT), size = 4, shape = 21) +
  scale_fill_manual(values = c("Before" = "#1f77b4",
                               "Ongoing" = "#FB2C36")) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        legend.position = "right") +
  annotate("text", x = min(m2$X1), y = max(m2$X2), 
           label = paste0("PERMANOVA P=", res_t$p_adonis[1]),
           hjust = 0, vjust = 1, size = 4, color = "black") + 
  labs(x = paste0("PCo1 (", var_exp1, " %)"),
       y = paste0("PCo2 (", var_exp2, " %)"),
       fill = "TNT") ; p4_pcoa_TNT_cat



###############  2. Paired  ###############

# Filtering samples (SNU_ID가 Before, Ongoing에 모두 존재)
mp = m %>% 
  group_by(SNU_ID) %>% 
  filter(n() > 1) %>% 
  ungroup

table(mp$TRG_score)
table(mp$cTstage)
# 32 Samples (CR 14. nearCR 4, PR 10, Poor 4)
# cTstage (cT2: 4, cT3: 20, cT4a: 6, cT4b: 2)


# Abundance table
tp = t %>% select(c("Strain", mp$SampleID))


## Beta-diversity: PCoA
t_dist = tp[, -1] %>% 
  t() %>% 
  as.matrix() %>% 
  vegdist()

t_ord = cmdscale(t_dist, k = 10, eig = T)

mp2 = t_ord$points %>% 
  data.frame() %>% 
  rownames_to_column("SampleID") %>% 
  merge(., mp, by = "SampleID") %>% 
  arrange(SubjectID, TNT) %>% 
  mutate(SampleID = factor(SampleID, levels = m$SampleID)) %>% 
  arrange(SampleID) ; head(mp2)

# Variance explained
var_exp1 = round(t_ord$eig/sum(t_ord$eig) * 100, 1)[1] ; var_exp1  # 12.6
var_exp2 = round(t_ord$eig/sum(t_ord$eig) * 100, 1)[2] ; var_exp2  # 10.3


# 1. 고정된 feature 변수: 거리 행렬 구성용
features = c("X1", "X2")

# 2. 그룹 변수 리스트 (모두 binary)
group_vars = c("Age_bin", "Sex", "Pre_Op_Tstage_bin", "Pre_Op_Nstage_bin", 
               "TRG_1", "TRG_2", "TRG_3", "CEA_bin", "BMI_bin")


# 거리 행렬 고정
temp_data = mp2 %>% 
  select(all_of(features)) %>%
  drop_na()

t_dist = dist(temp_data)  # 또는 vegdist() for Bray-Curtis

# 결과 저장용
res_group_p = data.frame(); set.seed(123) ; for (gvar in group_vars) {
  
  # 그룹 정보 추출
  group_data = mp2 %>%
    select(Group = !!sym(gvar)) %>%
    drop_na()
  
  # adonis2 수행
  adonis_res = adonis2(t_dist ~ Group, data = group_data, permutations = 999)
  adonis_p = adonis_res$`Pr(>F)`[1]
  
  # 결과 저장
  res_group_p = rbind(res_group_p,
                      data.frame(Group_Variable = gvar,
                                 p_adonis = adonis_p))
} ; res_group_p

# Group_Variable p_adonis
# 1           Age_bin    0.097*
# 2               Sex    0.335
# 3 Pre_Op_Tstage_bin    0.303
# 4 Pre_Op_Nstage_bin    0.742
# 5             TRG_1    0.331
# 6             TRG_2    0.038**
# 7             TRG_3    0.542
# 8           CEA_bin    0.004**
# 9           BMI_bin    0.528


# Wilcoxon test & Welch's t-test 
res_t_p = data.frame(Variable = character(),
                     Group = character(),
                     p_wilcox = numeric(),
                     p_ttest = numeric(),
                     stringsAsFactors = F)

# 그룹 변수와 피처 변수에 대해 반복
set.seed(123); for (gvar in group_vars) {
  for (var in features) {
    
    # 테스트용 데이터 구성
    ptest = mp2 %>%
      select(Y = !!sym(var), Group = !!sym(gvar)) %>%
      drop_na()
    
    # Wilcoxon test
    wilcox_res = wilcox.test(Y ~ Group, data = ptest)
    wilcox_p = wilcox_res$p.value
    
    # Welch’s t-test
    ttest_res = t.test(Y ~ Group, data = ptest)
    ttest_p = ttest_res$p.value
    
    # 결과 저장
    res_t_p = rbind(res_t_p,
                    data.frame(Variable = var,
                               Group = gvar,
                               p_wilcox = wilcox_p,
                               p_ttest = ttest_p))
  }
};res_t_p

# Variable             Group   p_wilcox    p_ttest
# 1        X1           Age_bin 0.09890503*    0.04219736**
# 2        X2           Age_bin 0.95522078     0.93630688
# 3        X1               Sex 0.98471254     0.71653897
# 4        X2               Sex 0.16951608*    0.13701203*
# 5        X1 Pre_Op_Tstage_bin 0.20850945     0.06483322*
# 6        X2 Pre_Op_Tstage_bin 0.42363737     0.45034507
# 7        X1 Pre_Op_Nstage_bin 0.98123797     0.23449756
# 8        X2 Pre_Op_Nstage_bin 0.86928819     0.82472747
# 9        X1             TRG_1 0.95522078     0.28448209
# 10       X2             TRG_1 0.35710143     0.29308158
# 11       X1             TRG_2 0.44153916     0.13887238*
# 12       X2             TRG_2 0.13496013*    0.05006114*
# 13       X1             TRG_3 1.00000000     0.60151207
# 14       X2             TRG_3 0.84760845     0.76789939
# 15       X1           CEA_bin 1.544926e-05** 7.936319e-07**
# 16       X2           CEA_bin 0.4356008**    0.1433232*
# 17       X1           BMI_bin 0.654656       0.1515174*
# 18       X2           BMI_bin 0.5878909      0.5842439


# Alpha diversity: Shannon diversity
# TRG_1
p5_shannon_TRG_cat_p = mp2 %>%
  mutate(TRG = factor(TRG_1, levels = c("CR", "nonCR")),
         TNT = factor(TNT, levels = c("Before", "Ongoing")),
         Group = paste(TNT, TRG_1, sep = "_"),
         Group = factor(Group, levels = c("Before_CR", "Ongoing_CR", "Before_nonCR", "Ongoing_nonCR"))) %>% 
  ggplot(aes(Group, Shannon)) +
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
  labs(x = "Sample", y = "Shannon diversity") +
  coord_cartesian(ylim = c(2, 5)) ; p5_shannon_TRG_cat_p

# ggsave("figure/05-1-1_Shannon-strain_cat_paired.svg",
#        plot = p5_shannon_TRG_cat_p, width = 3, height = 5)


# Alpha diversity: Observed strain
# TRG_1
p5_observed_TRG_cat_p = mp2 %>%
  mutate(TRG = factor(TRG_1, levels = c("CR", "nonCR")),
         TNT = factor(TNT, levels = c("Before", "Ongoing")),
         Group = paste(TNT, TRG_1, sep = "_"),
         Group = factor(Group, levels = c("Before_CR", "Ongoing_CR", "Before_nonCR", "Ongoing_nonCR"))) %>% 
  ggplot(aes(Group, ObservedStrain)) +
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
  labs(x = "Sample", y = "Observed strains") +
  coord_cartesian(ylim = c(100, 500)) ; p5_observed_TRG_cat_p

# ggsave("figure/05-1-2_Observed-strain_cat_paired.svg",
#        plot = p5_observed_TRG_cat_p, width = 3, height = 5)


# Alpha diversity: Inverted Simpson
# TRG_1
p5_invsimpson_TRG_cat_p = mp2 %>%
  mutate(TRG = factor(TRG_1, levels = c("CR", "nonCR")),
         TNT = factor(TNT, levels = c("Before", "Ongoing")),
         Group = paste(TNT, TRG_1, sep = "_"),
         Group = factor(Group, levels = c("Before_CR", "Ongoing_CR", "Before_nonCR", "Ongoing_nonCR"))) %>% 
  ggplot(aes(Group, InvSimpson)) +
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
  labs(x = "Sample", y = "Inverted Simpson") ; p5_invsimpson_TRG_cat_p

# ggsave("figure/05-1-3_Inverted Simpson-strain_cat_paired.svg",
#        plot = p5_invsimpson_TRG_cat_p, width = 3, height = 5)


# Alpha diversity: Pielou's evenness
even_t = t(tp[, -1])

H = diversity(even_t, index = "shannon")
S = specnumber(even_t)

evenness_t = H / log(S)
even_df_t = data.frame(
  SampleID = rownames(even_t),
  Evenness = evenness_t
)

mp2 = even_df_t %>% left_join(mp2, by = "SampleID") ; rm(even_t, even_df_t, evenness_t, H, S)

#TRG_1
p5_evenness_TRG_cat_p = mp2 %>%
  mutate(TRG = factor(TRG_1, levels = c("CR", "nonCR")),
         TNT = factor(TNT, levels = c("Before", "Ongoing")),
         Group = paste(TNT, TRG_1, sep = "_"),
         Group = factor(Group, levels = c("Before_CR", "Ongoing_CR", "Before_nonCR", "Ongoing_nonCR"))) %>% 
  ggplot(aes(Group, Evenness)) +
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
  labs(x = "Sample", y = "Evenness") ; p5_evenness_TRG_cat_p

# ggsave("figure/05-1-4_Evenness-strain_cat_paired.svg",
#        plot = p5_evenness_TRG_cat_p, width = 3, height = 5)




# Remove unnecessary objects
rm(list = ls(pattern = "^p4")) ; rm(list = ls(pattern = "^p5_"))
rm(adonis_p, adonis_res, ttest_p, ttest_res, wilcox_p, wilcox_res, temp_data, t_ord, t_ord_all,
   var, gvar, med_val, t_dist, t_dist_all, features, genus, ptest, group_data, bin_var)

save.image(file = "input/R_image/7-1. after-diversity-analysis.RData")

