
# Sample Filter: TNT - Before

rm(list = ls())
options(java.parameters = "-Xmx64g", stringsAsFactors = F)
setwd("C:/Users/user/Desktop/윤채빈/Data/CRC Metagenomics")

library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggforce)


load("input/R_image/7-1. after-diversity-analysis.RData")
load("input/R_image/7-2. after-taxonomy-barplot.RData")



###############  Taxonoimc barplot  ###############

# At the Family & Genus levels


###############  1. Family  ###############

fb_abund = fb %>% 
  column_to_rownames("Family") %>% 
  rowSums() %>% 
  sort(decreasing = T) %>% 
  data.frame() %>% 
  set_colnames("Total.abund") %>% 
  mutate(avg.abund = round(Total.abund / nrow(mb), 3)) %>% 
  rownames_to_column("Family")

fb_top10 = fb_abund %>% head(10) %>% pull(Family) ; fb_top10
# [1] "Lachnospiraceae"            "Oscillospiraceae"           "Bifidobacteriaceae"        
# [4] "Bacteroidaceae"             "Erysipelotrichaceae"        "Clostridiaceae"            
# [7] "Prevotellaceae"             "Streptococcaceae"           "Eubacteriales_unclassified"
# [10] "Coriobacteriaceae" 

fb_input = fb %>% 
  filter(Family %in% fb_top10) %>% 
  column_to_rownames("Family")

fb_input["Others", ] = apply(fb_input, 2, 
                             function(column) {100 - sum(column)})

fb_input = fb_input %>% 
  rownames_to_column("Family") %>% 
  pivot_longer(cols = -Family, names_to = "SampleID", values_to = "Abund") %>% 
  left_join(., mb %>% 
              select(SampleID, TRG, TRG_1, TRG_2, TRG_3, 
                     Pre_Op_Tstage, Pre_Op_Tstage_bin))


col_fb_top10 = c("Others" = "gray",
                 "Lachnospiraceae" = "#FF6B6B",
                 "Oscillospiraceae" = "#f7d9bc",
                 "Bifidobacteriaceae" = "#E0B0FF",
                 "Bacteroidaceae" = "#F88379",
                 "Erysipelotrichaceae" = "#FBEC5D",
                 "Clostridiaceae" = "#B5C7EB",
                 "Prevotellaceae" = "#00bb77",
                 "Streptococcaceae" = "#E89EB8",
                 "Eubacteriales_unclassified" = "#6395ee",
                 "Coriobacteriaceae" = "#FFFF99")


# Taxonomy barplot by TRG score
p4_taxonomy_TRG_f = fb_input %>% 
  group_by(Family, TRG) %>% 
  summarise(Abund = sum(Abund)) %>% 
  mutate(Family = factor(Family, levels = rev(c(fb_top10, "Others")))) %>% 
  ggplot(aes(TRG, Abund)) +
  geom_bar(aes(fill = Family), stat = "identity", position = "fill") +
  scale_fill_manual(values = col_fb_top10) +
  theme_pubr() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = rel(0.8)),
        axis.title.x = element_blank()) +
  labs(y = "Relative abundance") ; p4_taxonomy_TRG_f

# ggsave("figure/04-2-1_barplot-family.svg",
#        plot = p4_taxonomy_TRG_f, width = 5, height = 5.5)

# Taxonomy barplot by TRG_1
p4_taxonomy_TRG_cat_f = fb_input %>% 
  group_by(Family, TRG_1) %>% 
  summarise(Abund = sum(Abund)) %>% 
  mutate(Family = factor(Family, levels = rev(c(fb_top10, "Others")))) %>% 
  ggplot(aes(TRG_1, Abund)) +
  geom_bar(aes(fill = Family), stat = "identity", position = "fill") +
  scale_fill_manual(values = col_fb_top10) +
  theme_pubr() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = rel(0.8)),
        axis.title.x = element_blank()) +
  labs(y = "Relative abundance") ; p4_taxonomy_TRG_cat_f

# ggsave("figure/04-2-2_barplot-family_cat.svg",
#        plot = p4_taxonomy_TRG_cat_f, width = 4.5, height = 5.5)

fb_input %>% 
  group_by(Family, TRG_1) %>% 
  summarise(Abund = sum(Abund)) %>% 
  mutate(Family = factor(Family, levels = rev(c(fb_top10, "Others")))) %>% 
  pivot_wider(names_from = TRG_1, values_from = Abund) %>% 
  arrange(desc(Family))

fb_input %>% 
  group_by(Family) %>% 
  summarise(p_value = wilcox.test(Abund[TRG_1 == "CR"] / sum(Abund),
                                  Abund[TRG_1 == "nonCR"] / sum(Abund))$p.value) %>% 
  mutate(p_adj = p.adjust(p_value, method = "fdr"),
         Family = factor(Family, levels = rev(c(fb_top10, "Others")))) %>% 
  arrange(desc(Family))
  

# By sample - TRG_1
p4_taxonomy_TRG_cat_f_by_sample_CR = fb_input %>% 
  filter(TRG_1 == "CR") %>% 
  mutate(SampleID = factor(SampleID, levels = fb_input %>%
                             filter(Family == "Lachnospiraceae") %>%
                             arrange(desc(Abund)) %>%
                             pull(SampleID))) %>%
  mutate(Family = factor(Family, levels = rev(c(fb_top10, "Others")))) %>% 
  ggplot(aes(SampleID, Abund)) +
  geom_bar(aes(fill = Family), stat = "identity", position = "fill") +
  scale_fill_manual(values = col_fb_top10) +
  theme_pubr() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = rel(0.8)),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "Relative abundance") ; p4_taxonomy_TRG_cat_f_by_sample_CR

# ggsave("figure/04-2-2_barplot-family_sample_CR.svg",
#        plot = p4_taxonomy_TRG_cat_f_by_sample_CR, width = 7.5, height = 5.5)

p4_taxonomy_TRG_cat_f_by_sample_nonCR = fb_input %>% 
  filter(TRG_1 == "nonCR") %>% 
  mutate(SampleID = factor(SampleID, levels = fb_input %>%
                             filter(Family == "Lachnospiraceae") %>%
                             arrange(desc(Abund)) %>%
                             pull(SampleID))) %>%
  mutate(Family = factor(Family, levels = rev(c(fb_top10, "Others")))) %>% 
  ggplot(aes(SampleID, Abund)) +
  geom_bar(aes(fill = Family), stat = "identity", position = "fill") +
  scale_fill_manual(values = col_fb_top10) +
  theme_pubr() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = rel(0.8)),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "Relative abundance") ; p4_taxonomy_TRG_cat_f_by_sample_nonCR

# ggsave("figure/04-2-2_barplot-family_sample_nonCR.svg",
#        plot = p4_taxonomy_TRG_cat_f_by_sample_nonCR, width = 7.5, height = 5.5)



### In ongoing samples

# Sorting "Ongoing" samples
fo = f %>% select(c("Family", mo$SampleID))

fo_abund = fo %>% 
  column_to_rownames("Family") %>% 
  rowSums() %>% 
  sort(decreasing = T) %>% 
  data.frame() %>% 
  set_colnames("Total.abund") %>% 
  mutate(avg.abund = round(Total.abund / nrow(mo), 3)) %>% 
  rownames_to_column("Family")

fo_top10 = fo_abund %>% head(10) %>% pull(Family) ; fo_top10
# [1] "Lachnospiraceae"     "Oscillospiraceae"    "Bifidobacteriaceae"  
# [4] "Erysipelotrichaceae" "Enterobacteriaceae"  "Clostridiaceae"      
# [7] "Lactobacillaceae"    "Coriobacteriaceae"   "Bacteroidaceae"      
# [10] "Selenomonadaceae" 

fo_input = fo %>% 
  filter(Family %in% fo_top10) %>% 
  column_to_rownames("Family")

fo_input["Others", ] = apply(fo_input, 2,
                             function(column) {100 - sum(column)})

fo_input = fo_input %>% 
  rownames_to_column("Family") %>% 
  pivot_longer(cols = -Family, names_to = "SampleID", values_to = "Abund") %>% 
  left_join(., mo %>% 
              select(SampleID, TRG, TRG_1, TRG_2, TRG_3,
                     Pre_Op_Tstage, Pre_Op_Tstage_bin))


col_fo_top10 = c("Others" = "gray",
                 "Lachnospiraceae" = "#FF6B6B",
                 "Oscillospiraceae" = "#f7d9bc",
                 "Bifidobacteriaceae" = "#E0B0FF",
                 "Bacteroidaceae" = "#F88379",
                 "Erysipelotrichaceae" = "#FBEC5D",
                 "Clostridiaceae" = "#B5C7EB",
                 "Coriobacteriaceae" = "#FFFF99",
                 "Enterobacteriaceae" = "#C3D9A4",
                 "Lactobacillaceae" = "#80CBC4",
                 "Selenomonadaceae" = "#8C99D6")

fo_top10 = factor(names(col_fo_top10)[2:11], levels = names(col_fo_top10)[2:11]) %>% as.character()


# Taxonomy barplot by TRG score
p5_taxonomy_TRG_f_o = fo_input %>% 
  group_by(Family, TRG) %>% 
  summarise(Abund = sum(Abund)) %>% 
  mutate(Family = factor(Family, levels = rev(c(fo_top10, "Others")))) %>% 
  ggplot(aes(TRG, Abund)) +
  geom_bar(aes(fill = Family), stat = "identity", position = "fill") +
  scale_fill_manual(values = col_fo_top10) +
  theme_pubr() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = rel(0.8)),
        axis.title.x = element_blank()) +
  labs(y = "Relative abundance") ; p5_taxonomy_TRG_f_o

# ggsave("figure/05-2-1_barplot-family_ongoing.svg",
#        plot = p4_taxonomy_TRG_f_o, width = 5, height = 5.5)

# Taxonomy barplot by TRG_1
p5_taxonomy_TRG_cat_f_o = fo_input %>% 
  group_by(Family, TRG_1) %>% 
  summarise(Abund = sum(Abund)) %>% 
  mutate(Family = factor(Family, levels = rev(c(fo_top10, "Others")))) %>% 
  ggplot(aes(TRG_1, Abund)) +
  geom_bar(aes(fill = Family), stat = "identity", position = "fill") +
  scale_fill_manual(values = col_fo_top10) +
  theme_pubr() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = rel(0.8)),
        axis.title.x = element_blank()) +
  labs(y = "Relative abundance") ; p5_taxonomy_TRG_cat_f_o

# ggsave("figure/05-2-2_barplot-family_cat_ongoing.svg",
#        plot = p5_taxonomy_TRG_cat_f_o, width = 4.5, height = 5.5)

fo_input %>% 
  group_by(Family, TRG_1) %>% 
  summarise(Abund = sum(Abund)) %>% 
  mutate(Family = factor(Family, levels = c(fo_top10, "Others"))) %>% 
  pivot_wider(names_from = TRG_1, values_from = Abund) %>% 
  arrange(Family)

fo_input %>% 
  group_by(Family) %>% 
  summarise(p_value = wilcox.test(Abund[TRG_1 == "CR"] / sum(Abund),
                                  Abund[TRG_1 == "nonCR"] / sum(Abund))$p.value) %>% 
  mutate(p_adj = p.adjust(p_value, method = "fdr"),
         Family = factor(Family, levels = c(fo_top10, "Others"))) %>% 
  arrange(Family)

# By sample - TRG_1
p5_taxonomy_TRG_cat_f_by_sample_CR_o = fo_input %>% 
  filter(TRG_1 == "CR") %>% 
  mutate(SampleID = factor(SampleID, levels = fo_input %>%
                             filter(Family == "Lachnospiraceae") %>%
                             arrange(desc(Abund)) %>%
                             pull(SampleID))) %>%
  mutate(Family = factor(Family, levels = rev(c(fo_top10, "Others")))) %>% 
  ggplot(aes(SampleID, Abund)) +
  geom_bar(aes(fill = Family), stat = "identity", position = "fill") +
  scale_fill_manual(values = col_fo_top10) +
  theme_pubr() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = rel(0.8)),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "Relative abundance") ; p5_taxonomy_TRG_cat_f_by_sample_CR_o

# ggsave("figure/05-2-2_barplot-family_sample_CR.svg",
#        plot = p5_taxonomy_TRG_cat_f_by_sample_CR_o, width = 6, height = 5.5)

p5_taxonomy_TRG_cat_f_by_sample_nonCR_o = fo_input %>% 
  filter(TRG_1 == "nonCR") %>% 
  mutate(SampleID = factor(SampleID, levels = fo_input %>%
                             filter(Family == "Lachnospiraceae") %>%
                             arrange(desc(Abund)) %>%
                             pull(SampleID))) %>%
  mutate(Family = factor(Family, levels = rev(c(fo_top10, "Others")))) %>% 
  ggplot(aes(SampleID, Abund)) +
  geom_bar(aes(fill = Family), stat = "identity", position = "fill") +
  scale_fill_manual(values = col_fo_top10) +
  theme_pubr() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = rel(0.8)),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "Relative abundance") ; p5_taxonomy_TRG_cat_f_by_sample_nonCR_o

# ggsave("figure/05-2-2_barplot-family_sample_nonCR.svg",
#        plot = p5_taxonomy_TRG_cat_f_by_sample_nonCR_o, width = 6, height = 5.5)



###############  2. Genus  ###############

gb_abund = gb %>% 
  column_to_rownames("Genus") %>% 
  rowSums() %>% 
  sort(decreasing = T) %>% 
  data.frame() %>% 
  set_colnames("Total.abund") %>% 
  mutate(avg.abund = round(Total.abund / nrow(mb), 3)) %>% 
  rownames_to_column("Genus")

gb_top10 = gb_abund %>% head(10) %>% pull(Genus) ; gb_top10
# [1] "Blautia"              "Faecalibacterium"      "Bifidobacterium"             
# [4] "Anaerostipes"         "Bacteroides"           "Mediterraneibacter"          
# [7] "Ruminococcus"         "Dorea"                 "Lachnospiraceae_unclassified"
# [10] "Fusicatenibacter"

gb_input = gb %>% 
  filter(Genus %in% gb_top10) %>% 
  column_to_rownames("Genus")

gb_input["Others", ] = apply(gb_input, 2, 
                             function(column) { 100 - sum(column)})

gb_input = gb_input %>% 
  rownames_to_column("Genus") %>% 
  pivot_longer(cols = -Genus, names_to = "SampleID", values_to = "Abund") %>% 
  left_join(., mb %>% 
              select(SampleID, TRG, TRG_1, TRG_2, TRG_3, 
                     Pre_Op_Tstage, Pre_Op_Tstage_bin))


col_gb_top10 = c("Others" = "gray",
                 "Blautia" = "#FF6B6B",
                 "Faecalibacterium" = "#f7d9bc",
                 "Bifidobacterium" = "#E0B0FF",
                 "Anaerostipes" = "#FBEC5D",
                 "Bacteroides" = "#F88379",
                 "Mediterraneibacter" = "#B5C7EB",
                 "Ruminococcus" = "#00bb77",
                 "Dorea" = "#E89EB8",
                 "Lachnospiraceae_unclassified" = "#FFFF99",
                 "Fusicatenibacter" = "#6395ee")


# Taxonomy barplot by TRG score
p4_taxonomy_TRG_g = gb_input %>% 
  group_by(Genus, TRG) %>% 
  summarise(Abund = sum(Abund)) %>% 
  mutate(Genus = factor(Genus, levels = rev(c(gb_top10, "Others")))) %>% 
  ggplot(aes(TRG, Abund)) +
  geom_bar(aes(fill = Genus), stat = "identity", position = "fill") +
  scale_fill_manual(values = col_gb_top10) +
  theme_pubr() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = rel(0.8)),
        axis.title.x = element_blank()) +
  labs(y = "Relative abundance") ; p4_taxonomy_TRG_g

# ggsave("figure/04-2-3_barplot-genus.svg",
#        plot = p4_taxonomy_TRG_g, width = 5, height = 5.5)

# Taxonomy barplot by TRG_1
p4_taxonomy_TRG_cat_g = gb_input %>% 
  group_by(Genus, TRG_1) %>% 
  summarise(Abund = sum(Abund)) %>% 
  mutate(Genus = factor(Genus, levels = rev(c(gb_top10, "Others")))) %>% 
  ggplot(aes(TRG_1, Abund)) +
  geom_bar(aes(fill = Genus), stat = "identity", position = "fill") +
  scale_fill_manual(values = col_gb_top10) +
  theme_pubr() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = rel(0.8)),
        axis.title.x = element_blank()) +
  labs(y = "Relative abundance") ; p4_taxonomy_TRG_cat_g

# ggsave("figure/04-2-4_barplot-genus_cat.svg",
#        plot = p4_taxonomy_TRG_cat_g, width = 4.5, height = 5.5)



### In ongoing samples

# Sorting "Ongoing" samples
go = g %>% select(c("Genus", mo$SampleID))

go_abund = go %>% 
  column_to_rownames("Genus") %>% 
  rowSums() %>% 
  sort(decreasing = T) %>% 
  data.frame() %>% 
  set_colnames("Total.abund") %>% 
  mutate(avg.abund = round(Total.abund / nrow(mo), 3)) %>% 
  rownames_to_column("Genus")

go_top10 = go_abund %>% head(10) %>% pull(Genus) ; go_top10
# [1] "Blautia"           "Bifidobacterium"         "Faecalibacterium"            
# [4] "Anaerostipes"      "Mediterraneibacter"      "Lachnospiraceae_unclassified"
# [7] "Ruminococcus"      "Anaerobutyricum"         "Dorea"                       
# [10] "Escherichia"

go_input = go %>% 
  filter(Genus %in% go_top10) %>% 
  column_to_rownames("Genus")

go_input["Others", ] = apply(go_input, 2,
                             function(column) {100 - sum(column)})

go_input = go_input %>% 
  rownames_to_column("Genus") %>% 
  pivot_longer(cols = -Genus, names_to = "SampleID", values_to = "Abund") %>% 
  left_join(., mo %>% 
              select(SampleID, TRG, TRG_1, TRG_2, TRG_3,
                     Pre_Op_Tstage, Pre_Op_Tstage_bin))


col_go_top10 = c("Others" = "gray",
                 "Blautia" = "#FF6B6B",
                 "Faecalibacterium" = "#f7d9bc",
                 "Bifidobacterium" = "#E0B0FF",
                 "Anaerostipes" = "#FBEC5D",
                 "Mediterraneibacter" = "#B5C7EB",
                 "Ruminococcus" = "#00bb77",
                 "Dorea" = "#E89EB8",
                 "Anaerobutyricum" = "#FFFF99",
                 "Lachnospiraceae_unclassified" = "#FF6B6B",
                 "Escherichia" = "#8C99D6")

go_top10 = factor(names(col_go_top10)[2:11], levels = names(col_go_top10)[2:11]) %>% as.character()


# Taxonomy barplot by TRG score
p5_taxonomy_TRG_g_o = go_input %>% 
  group_by(Genus, TRG) %>% 
  summarise(Abund = sum(Abund)) %>% 
  mutate(Genus = factor(Genus, levels = rev(c(go_top10, "Others")))) %>% 
  ggplot(aes(TRG, Abund)) +
  geom_bar(aes(fill = Genus), stat = "identity", position = "fill") +
  scale_fill_manual(values = col_go_top10) +
  theme_pubr() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = rel(0.8)),
        axis.title.x = element_blank()) +
  labs(y = "Relative abundance") ; p5_taxonomy_TRG_g_o

# ggsave("figure/05-2-3_barplot-genus_ongoing.svg",
#        plot = p5_taxonomy_TRG_g_o, width = 5, height = 5.5)

# Taxonomy barplot by TRG_1
p5_taxonomy_TRG_cat_g_o = go_input %>% 
  group_by(Genus, TRG_1) %>% 
  summarise(Abund = sum(Abund)) %>% 
  mutate(Genus = factor(Genus, levels = rev(c(go_top10, "Others")))) %>% 
  ggplot(aes(TRG_1, Abund)) +
  geom_bar(aes(fill = Genus), stat = "identity", position = "fill") +
  scale_fill_manual(values = col_go_top10) +
  theme_pubr() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = rel(0.8)),
        axis.title.x = element_blank()) +
  labs(y = "Relative abundance") ; p5_taxonomy_TRG_cat_g_o

# ggsave("figure/05-2-4_barplot-genus_cat_ongoing.svg",
#        plot = p5_taxonomy_TRG_cat_g_o, width = 4.5, height = 5.5)



# Remove unnecessary objects
rm(list = ls(pattern = "^p4_")) ; rm(list = ls(pattern = "^p5_"))
save.image(file = "input/R_image/7-2. after-taxonomy-barplot.RData")



