
# Sample Filter: TNT - Before

rm(list = ls())
options(java.parameters = "-Xmx64g", stringsAsFactors = F)
setwd("C:/Users/user/Desktop/윤채빈/Data/CRC Metagenomics")

library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggforce)
library(gridExtra)
library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)
library(circlize)
library(grid)


load("input/R_image/7-3. after-differential-enrichment-testing.RData")
load("input/R_image/7-4. after-heatmap.RData")



#################### Heatmap ####################

##### Strain #####

### At the Strain level by TRG group

head(tb_input)

# 1. Sorting DEB (P<0.1)
sig_strain_0.1 = wilcox_t %>% 
  filter(p_value < 0.1) %>% 
  arrange(p_value) %>% 
  pull(Strain)


# 2. Calculate mean abundance for the DEB
heat_tb = tb_input %>% 
  filter(Strain %in% sig_strain_0.1) %>% 
  group_by(Strain, TRG_1) %>% 
  summarise(Mean_Abundance = mean(abundance, na.rm = T), .groups = "drop") %>% 
  arrange(Strain, TRG_1)


# 3. Conversion to wide format
heat_tb_wide = heat_tb %>% 
  pivot_wider(names_from = TRG_1, values_from = Mean_Abundance) %>% 
  column_to_rownames("Strain") %>% 
  as.matrix()


# 4. Half of min. abundance
min_abund = tb_input$abundance[tb_input$abundance > 0] %>% min() # 2e-05


# 5. Log10 with pseudo-abundance
heat_tb_wide_log = log10(heat_tb_wide + min_abund/2)


# 6. Bacteria-wise scaling
heat_tb_scaled = t(scale(t(heat_tb_wide_log))) # z-score per row


# 7. 별표 표시용 annotation matrix 생성
get_asterisks = function(p) {
  if (is.na(p)) return("")
  if(p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("")
}

annots = sapply(rownames(heat_tb_wide_log), function(strain) {
  p = wilcox_t$p_value[wilcox_t$Strain == strain]
  get_asterisks(p)
})

annotation_matrix = matrix("",
                           nrow = nrow(heat_tb_wide_log),
                           ncol = ncol(heat_tb_wide_log))

rownames(annotation_matrix) = rownames(heat_tb_wide_log)
colnames(annotation_matrix) = colnames(heat_tb_wide_log)

annotation_matrix[, "CR"] = annots # CR 열에만 표시


# 8. Draw
color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(300)

p4_heat_t2 = pheatmap(
  mat = heat_tb_wide_log,
  cluster_rows = T,
  cluster_cols = F,
  clustering_distance_rows = "manhattan",  # 또는 "euclidean", "spearman", "pearson" 등
  clustering_method_rows = "ward.D2",      # 또는 "complete", "average", "single" 등
  color = color,
  scale = "none",
  angle_col = 0,
  display_numbers = annotation_matrix,
  number_color = "black",
  fontsize_number = 14,
  border_color = "gray70"
) ; p4_heat_t2

# ggsave("figure/04-4-1_heatmap_strain_pval0.1_TRG1.svg",
#        p4_heat_t2, width = 5.5, height = 10)



### At the Strain level by sample

# pseudo-abundance = 1e-05
min_abund/2

# 1. Wide-format 변환 (Strain x Sample matrix)
heat_tb_sample = tb_input %>% 
  filter(Strain %in% sig_strain_0.1) %>% 
  group_by(Strain, SampleID) %>% 
  summarise(abundance = sum(as.numeric(abundance), na.rm = T), .groups = "drop") %>% 
  pivot_wider(names_from = SampleID,
              values_from = abundance,
              values_fill = 0)


# 2. 행 이름 설정 + matrix 변환
heat_tb_mat = heat_tb_sample %>% 
  column_to_rownames("Strain") %>% 
  as.matrix() ; dim(heat_tb_mat) # 49 strains in 26 samples


# 3. Log10 변환 + pseudo-abundance 추가
heat_tb_log = log10(heat_tb_mat + min_abund/2)


# 4. Strain별 z-score scaling
heat_tb_scaled = t(scale(t(heat_tb_log)))
  # summary(t(heat_tb_scaled)) # Mean 0


# 5. 샘플 그룹 및 T-stage 정보 추출
  # mb2$Pre_Op_Tstage_bin %>% table()
  # 0  1 
  # 4 22 → 4 early & 22 advanced

sample_info = tb_input %>% 
  distinct(SampleID, TRG_1, Pre_Op_Tstage_bin) %>% 
  filter(SampleID %in% colnames(heat_tb_scaled)) %>% 
  mutate(Tstage = factor(Pre_Op_Tstage_bin,
                         levels = c(0, 1),
                         label = c("early", "advanced"))) %>% 
  select(SampleID, TRG_1, Tstage) %>% 
  column_to_rownames("SampleID")


# 6. CR / nonCR 나눈 후 각 그룹 안에서 clustering
heat_tb_CR = heat_tb_scaled[, sample_info$TRG_1 == "CR"]
heat_tb_nonCR = heat_tb_scaled[, sample_info$TRG_1 == "nonCR"]

col_order_tb_CR = colnames(heat_tb_CR)[hclust(dist(t(heat_tb_CR)))$order]
col_order_tb_nonCR = colnames(heat_tb_nonCR)[hclust(dist(t(heat_tb_nonCR)))$order]

final_tb_order = c(col_order_tb_CR, col_order_tb_nonCR)


# 7. 열 순서 및 주석 재정렬
heat_tb_final = heat_tb_scaled[, final_tb_order]
annotation_col = sample_info[final_tb_order, , drop = F]

col_fun_main = circlize::colorRamp2(
  c(-2, 0, 2),
  # Blue - Light Grey (almost white) - Red
  c("#4575b4", "#F8F8F8", "#d73027")
)


# 8. 색상 설정
ann_colors = list(
  TRG_1 = c(CR = "#F7D9BC", nonCR = "#80461B"),
  Tstage = c(early = "#1b9e77", advanced = "#d95f02")
)


# 기본 ht_main
set.seed(123) ; ht_main = Heatmap(
  matrix = heat_tb_final,
  name = "z-score",
  top_annotation = HeatmapAnnotation(df = annotation_col,
                                     col = ann_colors,
                                     annotation_name_side = "left"),
  cluster_columns = F,
  cluster_rows = T,
  clustering_distance_rows = "manhattan", # 또는 "euclidean", "spearman", "pearson" 등
  clustering_method_rows = "ward.D2",     # 또는 "complete", "average", "single" 등
  
  show_column_names = F,
  show_row_names = T,
  row_names_gp = gpar(fontsize = 8),
  border = T,
  heatmap_legend_param = list(
    title = "z-score",
    title_position = "topcenter",
    legend_direction = "horizontal",
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2")
  ),
  col = col_fun_main,
  column_title = "Strain Abundance (log10 + z-score)\nGrouped CR vs nonCr, hierarchical clusteirng",
  rect_gp = gpar(col = "black", lwd = 0.3) # 테두리 추가 (선 색 + 두께)
) ; ht_main

# svg("figure/04-4-2_heatmap_samples_strain.svg", width = 12, height = 15) ; draw(
#   ht_main,
#   heatmap_legend_side = "bottom",       # 또는 "left", "top"으로 조정 가능
#   annotation_legend_side = "bottom",    # T_stage, TRG_1 legend는 그대로
#   padding = unit(c(5, 10, 5, 5), "mm")  # top, right, bottom, left padding


### Add a heatmap to display p-value

# 1. Strain 순서 추출 (위에서 아래로 보이는 순서 기준)
strain_order = rownames(heat_tb_final)[order.dendrogram(row_dend(ht_main))]


# 2. p-value matrix
pval_data = wilcox_t %>% 
  filter(Strain %in% strain_order) %>% 
  mutate(
    log10_p = -log10(p_value),
    sig_star = case_when(
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      T ~ "")) %>% 
  column_to_rownames("Strain")

pval_mat = matrix(
  pval_data[strain_order, "log10_p", drop = T],
  ncol = 1,
  dimnames = list(strain_order, "-log10(p)"))

summary(pval_mat)
range(pval_mat, na.rm = T) # 1.012480 ~ 2.516706

sig_text = matrix(
  pval_data[strain_order, "sig_star", drop = T],
  ncol = 1,
  dimnames = list(strain_order, "-log10(p)"))


# 3. 색상 스케일 (정상 순서로: 낮음 → 진청색)
col_fun_pval = colorRamp2(
  c(1.0, 1.8, 2.6),  # Adjusted to your actual value range
  c("white", "#6BAED6", "#08306B")  # white → soft blue → deep blue
)


# 4. Heatmap 객체
ht_pval = Heatmap(
  matrix = pval_mat,
  name = "-log10(p)",
  col = col_fun_pval,
  cluster_rows = F,
  cluster_columns = F,
  width = unit(0.5, "cm"),
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  row_order = strain_order,
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (!is.na(fill)) {
      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "grey60", lwd = 0.4))
    } else {
      grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "grey60", lwd = 0.4))
    }
    grid.text(sig_text[i, j], x = x, y = y, gp = gpar(fontsize = 10))
  }
)

# svg("figure/04-4-2_heatmap_samples_strain_pval.svg",
#     width = 4, height = 13); draw(ht_pval, heatmap_legend_side = "right"); dev.off()



### In ongoing samples

head(to_input)

# 1. Sorting DEB (P < 0.1)
sig_strain_0.1_o = wilcox_to %>% 
  filter(p_value < 0.1) %>% 
  arrange(p_value) %>% 
  pull(Strain)


# 2. Calculate mean abundance for the DEB
heat_to = to_input %>% 
  filter(Strain %in% sig_strain_0.1_o) %>% 
  group_by(Strain, TRG_1) %>% 
  summarise(Mean_Abundance = mean(abundance, na.rm = T), .groups = "drop") %>% 
  arrange(Strain, TRG_1)


# 3. Conversion to wide format
heat_to_wide = heat_to %>% 
  pivot_wider(names_from = TRG_1, values_from = Mean_Abundance) %>% 
  column_to_rownames("Strain") %>% 
  as.matrix()

  
# 4. Half of min. abundance
min_abund = to_input$abundance[to_input$abundance > 0] %>% min() # 4e-05


# 5. Log10 with pseudo-abundance
heat_to_wide_log = log10(heat_to_wide + min_abund/2)


# 6. Bacteria-wise scaling
heat_to_scaled = t(scale(t(heat_to_wide_log))) # z-score per row


# 7. 별표 표시용 annotation matrix 생성
get_asterisks = function(p) {
  if(is.na(p)) return("")
  if(p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("")
}

annots = sapply(rownames(heat_to_wide_log), function(strain) {
  p = wilcox_to$p_value[wilcox_to$Strain == strain]
  get_asterisks(p)
})

annotation_matrix = matrix("",
                           nrow = nrow(heat_to_wide_log),
                           ncol = ncol(heat_to_wide_log))

rownames(annotation_matrix) = rownames(heat_to_wide_log)
colnames(annotation_matrix) = colnames(heat_to_wide_log)

annotation_matrix[, "CR"] = annots # CR 열에만 표시


# 8. Draw
color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(300)

p5_heat_t2_on = pheatmap(
  mat = heat_to_wide_log,
  cluster_rows = T,
  cluster_cols = F,
  clustering_distance_rows = "manhattan", # 또는 "euclidean", "spearman", "pearson" 등
  clustering_method_rows = "ward.D2",     # 또는 "complete", "average", "single" 등
  color = color,
  scale = "none",
  angle_col = 0,
  display_numbers = annotation_matrix,
  number_color = "black",
  fontsize_number = 14,
  border_color = "grey70"
) ; p5_heat_t2_on

# ggsave("figure/05-4-1_heatmap_strain_pval0.1_TRG1_ongoing.svg",
#        p5_heat_t2_on, width = 6.5, height = 5)


# pseudo-abundance = 4e-05
min_abund/2

# 1. Wide-format 변환 (Strain x Sample matrix)
heat_to_sample = to_input %>% 
  filter(Strain %in% sig_strain_0.1_o) %>% 
  group_by(Strain, SampleID) %>% 
  summarise(abundance = sum(as.numeric(abundance), na.rm = T), .groups = "drop") %>% 
  pivot_wider(names_from = SampleID,
              values_from = abundance,
              values_fill = 0)


# 2. 행 이름 설정 + matrix 변환
heat_to_mat = heat_to_sample %>% 
  column_to_rownames("Strain") %>% 
  as.matrix() ; dim(heat_to_mat) # 19 strains in 16 samples


# 3. Log10 변환 + pseudo-abundance 추가
heat_to_log = log10(heat_to_mat + min_abund/2)


# 4. Strain별 z-score scaling
heat_to_scaled = t(scale(t(heat_to_log)))
# summary(t(heat_to_scaled)) # Mean 0


# 5. 샘플 그룹 및 T_stage 정보 추출
  # mo$Pre_Op_Tstage_bin %>% table()
  # 0  1 
  # 2 14 → 2 early & 14 advanced

sample_info_o = to_input %>% 
  distinct(SampleID, TRG_1, Pre_Op_Tstage_bin) %>% 
  filter(SampleID %in% colnames(heat_to_scaled)) %>% 
  mutate(Tstage = factor(Pre_Op_Tstage_bin,
                         levels = c(0, 1),
                         label = c("early", "advanced"))) %>% 
  select(SampleID, TRG_1, Tstage) %>% 
  column_to_rownames("SampleID")


# 6. CR / nonCR 나눈 후 각 그룹 안에서 clustering
heat_to_CR = heat_to_scaled[, sample_info_o$TRG_1 == "CR"]
heat_to_nonCR = heat_to_scaled[, sample_info_o$TRG_1 == "nonCR"]

col_order_to_CR = colnames(heat_to_CR)[hclust(dist(t(heat_to_CR)))$order]
col_order_to_nonCR = colnames(heat_to_nonCR)[hclust(dist(t(heat_to_nonCR)))$order]

final_to_order = c(col_order_to_CR, col_order_to_nonCR)


# 7. 열 순서 및 주석 재정렬
heat_to_final = heat_to_scaled[, final_to_order]
annotation_col = sample_info_o[final_to_order, , drop = F]

col_fun_main = circlize::colorRamp2(
  c(-2, 0, 2),
  # Blue - Light Grey (almost white) - Red
  c("#4575b4", "#F8F8F8", "#d73027"))


# 8. 색상 설정
ann_colors = list(
  TRG_1 = c(CR = "#F7D9BC", nonCR = "#80461B"),
  Tstage = c(early = "#1b9e77", advanced = "#d95f02"))


# 기본 ht_main
set.seed(123) ; ht_main = Heatmap(
  matrix = heat_to_final,
  name = "z-score",
  top_annotation = HeatmapAnnotation(df = annotation_col, 
                                     col = ann_colors,
                                     annotation_name_side = "left"),
  cluster_columns = F,
  cluster_rows = T,
  clustering_distance_rows = "manhattan", # 또는 "euclidean", "spearman", "pearson" 등
  clustering_method_rows = "ward.D2",     # 또는 "complete", "average", "single" 등
  
  show_column_names = F,
  show_row_names = T,
  row_names_gp = gpar(fontsize = 8),
  border = T, 
  heatmap_legend_param = list(
    title = "z-score",
    title_position = "topcenter",
    legend_direction = "horizontal",
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2")
  ),
  col = col_fun_main,
  column_title = "Strain Abundance (log10 + z-score)\nGrouped CR vs nonCR, hierarchical clustering",
  rect_gp = gpar(col = "black", lwd = 0.3)  # 테두리 추가 (선 색 + 두께)
) ; ht_main

# svg("figure/05-4-2_heatmap_samples_strain_ongoing.svg", width = 9, height = 7) ; draw(
#   ht_main,
#   heatmap_legend_side = "bottom",       # 또는 "left", "top"으로 조정 가능
#   annotation_legend_side = "bottom",    # T_stage, TRG_1 legend는 그대로
#   padding = unit(c(5, 10, 5, 5), "mm")  # top, right, bottom, left padding
# ) ; dev.off()


### Add a heatmap to display p-value

# 1. Strain 순서 추출 (위에서 아래로 보이는 순서 기준)
strain_order_o = rownames(heat_to_final)[order.dendrogram(row_dend(ht_main))]


# 2. p-value matrix
pval_data_o = wilcox_to %>% 
  filter(Strain %in% strain_order_o) %>% 
  mutate(
    log10_p = -log10(p_value),
    sig_star = case_when(
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "")) %>% 
  column_to_rownames("Strain")

pval_mat_o = matrix(
  pval_data_o[strain_order_o, "log10_p", drop = T],
  ncol = 1,
  dimnames = list(strain_order_o, "-log10(p)"))

summary(pval_mat_o)
range(pval_mat_o, na.rm = T)

sig_text_o = matrix(
  pval_data_o[strain_order_o, "sig_star", drop = T],
  ncol = 1,
  dimnames = list(strain_order_o, "-log10(p)"))


# 3. 색상 스케일 (정상 순서로: 낮음 → 진청색)
col_fun_pval = colorRamp2(
  c(1.0, 1.8, 2.6),
  c("white", "#6BAED6", "#08306B")  # white → soft blue → deep blue
)


# 4. Heatmap 객체
ht_pval = Heatmap(
  matrix = pval_mat_o,
  name = "-log10(p)",
  col = col_fun_pval,
  cluster_rows = F,
  cluster_columns = F,
  width = unit(0.5, "cm"),
  show_row_names = T,
  row_names_gp = gpar(fontsize = 10),
  row_order = strain_order_o,
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (!is.na(fill)) {
      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "grey60", lwd = 0.4))
    } else {
      grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "grey60", lwd = 0.4))
    }
    grid.text(sig_text_o[i, j], x = x, y = y, gp = gpar(fontsize = 10))
  }
) ; ht_pval

# svg("figure/05-4-2_heatmap_samples_strain_pval_ongoing.svg",
#     width = 4, height = 5); draw(ht_pval, heatmap_legend_side = "right"); dev.off()



##### Species #####

### At the Species level by TRG group

head(sb_input)

# 1. Sorting DEG (P<0.1)
sig_species_0.1 = wilcox_s %>% 
  filter(p_value < 0.1) %>% 
  arrange(p_value) %>% 
  pull(Species)


# 2. Calculate mean abundance for the DEB
heat_sb = sb_input %>% 
  filter(Species %in% sig_species_0.1) %>% 
  group_by(Species, TRG_1) %>% 
  summarise(Mean_Abundance = mean(abundance, na.rm = T), .groups = "drop") %>% 
  arrange(Species, TRG_1)


# 3. Conversion to wide format
heat_sb_wide = heat_sb %>% 
  pivot_wider(names_from = TRG_1, values_from = Mean_Abundance) %>% 
  column_to_rownames("Species") %>% 
  as.matrix()


# 4. Half of min. abundance
min_abund = sb_input$abundance[sb_input$abundance > 0] %>% min # 2e-05


# 5. Log10 with pseudo-abundance
heat_sb_wide_log = log10(heat_sb_wide + min_abund/2)


# 6. Bacteriqa-wise Scaling
heat_sb_scaled = t(scale(t(heat_sb_wide_log))) # z-score per row


# 7. 별표 표시용 annotation matrix 생성
get_asterisks = function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("")
}

annots = sapply(rownames(heat_sb_wide_log), function(species) {
  p = wilcox_s$p_value[wilcox_s$Species == species]
  get_asterisks(p)
})

annotation_matrix = matrix("",
                           nrow = nrow(heat_sb_wide_log),
                           ncol = ncol(heat_sb_wide_log))

rownames(annotation_matrix) = rownames(heat_sb_wide_log)
colnames(annotation_matrix) = colnames(heat_sb_wide_log)

annotation_matrix[, "CR"] = annots # CR 열에만 표시


# 8. Draw
color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(300)

p4_heat_s2 = pheatmap(
  mat = heat_sb_wide_log,
  cluster_rows = T,
  cluster_cols = F,
  
  clustering_distance_rows = "manhattan", # 또는 "euclidean", "correlatoin" 등
  clustering_method_rows = "average",     # 또는 "complete", "ward.D2", "single" 등
  
  color = color,
  scale = "none",
  angel_col = 0,
  display_numbers = annotation_matrix,
  number_color = "black",
  fontsize_number = 14,
  border_color = "grey70"
) ; p4_heat_s2

# ggsave("figure/04-4-3_heatmap_species_pval0.1_TRG1.svg",
#        p4_heat_s2, width = 5.5, height = 10)



### At the Species level by sample

# pseudo-abundance = 1e-05
min_abund/2

# 1. Wide-format 변환 (Species x Sample matrix)
heat_sb_sample = sb_input %>% 
  filter(Species %in% sig_species_0.1) %>% 
  group_by(Species, SampleID) %>% 
  summarise(abundance = sum(as.numeric(abundance), na.rm = T), .groups = "drop") %>% 
  pivot_wider(names_from = SampleID,
              values_from = abundance,
              values_fill = 0)


# 2. 행 이름 설정 + matrix 변환
heat_sb_mat = heat_sb_sample %>% 
  column_to_rownames("Species") %>% 
  as.matrix() ; dim(heat_sb_mat) # 25 species in 26 samples


# 3. Log10 변환 + pseudo-abundance 추가
heat_sb_log = log10(heat_sb_mat + min_abund/2)


# 4. Species별 z-score scaling
heat_sb_scaled = t(scale(t(heat_sb_log)))
# summary(t(heat_sb_scaled)) # Mean 0


# 5. 샘플 그룹 및 T_stage 정보 추출
sample_info = sb_input %>% 
  distinct(SampleID, TRG_1, Pre_Op_Tstage_bin) %>% 
  filter(SampleID %in% colnames(heat_sb_scaled)) %>% 
  mutate(Tstage = factor(Pre_Op_Tstage_bin,
                         levels = c(0, 1),
                         label = c("early", "advanced"))) %>% 
  select(SampleID, TRG_1, Tstage) %>% 
  column_to_rownames("SampleID")


# 6. CR / nonCR 나눈 후 각 그룹 안에서 clustering
heat_sb_CR = heat_sb_scaled[, sample_info$TRG_1 == "CR"]
heat_sb_nonCR = heat_sb_scaled[, sample_info$TRG_1 == "nonCR"]

col_order_sb_CR = colnames(heat_sb_CR)[hclust(dist(t(heat_sb_CR)))$order]
col_order_sb_nonCR = colnames(heat_sb_nonCR)[hclust(dist(t(heat_sb_nonCR)))$order]

final_sb_order = c(col_order_sb_CR, col_order_sb_nonCR)


# 7. 열 순서 및 주석 재정렬
heat_sb_final = heat_sb_scaled[, final_sb_order]
annotation_col = sample_info[final_sb_order, , drop = F]

col_fun_main = circlize::colorRamp2(
  c(-2, 0 ,2),
  # Blue - Light Grey (almost white) - Red
  c("#4575b4", "#F8F8F8", "#d73027")
)


# 8. 색상 설정
ann_colors = list(
  TRG_1 = c(CR = "#F7D9BC", nonCR = "#80461B"),
  Tstage = c(early = "#1b9e77", advanced = "#d95f02")
)


# 기본 ht_main
set.seed(123) ; ht_main = Heatmap(
  matrix = heat_sb_final,
  name = "z-score",
  top_annotation = HeatmapAnnotation(df = annotation_col,
                                     col = ann_colors,
                                     annotation_name_side = "left"),
  cluster_columns = F,
  cluster_rows = T,
  
  clustering_distance_rows = "manhattan", # 또는 "euclidean", "correlatoin" 등
  clustering_method_rows = "average",     # 또는 "complete", "ward.D2", "single" 등
  
  show_column_names = F,
  show_row_names = T,
  row_names_gp = gpar(fontsize = 8),
  border = T,
  heatmap_legend_param = list(
    title = "z-score",
    title_position = "topcenter",
    legend_direction = "horizontal",
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2")
  ),
  col = col_fun_main,
  column_title = "Species Abundance (log10 + z-score)\nGrouped CR vs nonCR, hierarchical clustering",
  rect_gp = gpar(col = "black", lwd = 0.3)  # 테두리 추가 (선 색 + 두께)
) ; ht_main

# svg("figure/04-4-4_heatmap_samples_species.svg", width = 12, height = 10) ; draw(
#   ht_main,
#   heatmap_legend_side = "bottom",       # 또는 "left", "top"으로 조정 가능
#   annotation_legend_side = "bottom",    # T_stage, TRG_1 legend는 그대로
#   padding = unit(c(5, 10, 5, 5), "mm")  # top, right, bottom, left padding
# ) ; dev.off()


### Add p-value

# 1. Species 순서 추출 (위에서 아래로 보이는 순서 기준)
species_order = rownames(heat_sb_final)[order.dendrogram(row_dend(ht_main))]


# 2. p-value matrix
pval_data = wilcox_s %>% 
  filter(Species %in% species_order) %>% 
  mutate(
    log10_p = -log10(p_value),
    sig_star = case_when(
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>% 
  column_to_rownames("Species")

pval_mat = matrix(
  pval_data[species_order, "log10_p", drop = T],
  ncol = 1,
  dimnames = list(species_order, "-log10(p)")
)

summary(pval_mat)
range(pval_mat, na.rm = T)

sig_text = matrix(
  pval_data[species_order, "sig_star", drop = T],
  ncol = 1,
  dimnames = list(species_order, "-log10(p)")
)


# 3. 색상 스케일 (정상 순서로: 낮음 → 진청색)
col_fun_pval = colorRamp2(
  c(1.0, 1.7, 2.5),
  c("white", "#6BAED6", "#08306B") # white → soft blue → deep blue
)


# 4. Heatmap 객체
ht_pval = Heatmap(
  matrix = pval_mat,
  name = "-log10(p)",
  col = col_fun_pval,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  width = unit(0.5, "cm"),
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  row_order = species_order,
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (!is.na(fill)) {
      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "grey60", lwd = 0.4))
    } else {
      grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "grey60", lwd = 0.4))
    }
    grid.text(sig_text[i, j], x = x, y = y, gp = gpar(fontsize = 10))
  }
)

# svg("figure/04-4-4_heatmap_samples_species_pval.svg",width = 4, height = 7); draw(
#   ht_pval, heatmap_legend_side = "right"); dev.off()



### In ongoing samples

# Prevalence filtering (>=10%)
so = s %>% select(c("Species", mo$SampleID))

so_0.2 = so %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Species) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.2) %>% 
  pull(Species) %>% 
  sort() ; so_0.2


# Make an input file
so_input = so %>% 
  filter(Species %in% so_0.2) %>% 
  pivot_longer(cols = -Species,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  merge(mo %>% 
          select(SampleID,
                 TRG_score, TRG_1, TRG_2, TRG_3,
                 Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                 Pre_Op_Nstage, Pre_Op_Nstage_bin),
        by = "SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR")))

head(so_input)


# Wilcoxon test (CR vs. nonCR)
wilcox_so = so_input %>% 
  group_by(Species) %>% 
  summarise(p_value = wilcox.test(abundance ~ TRG_1, exact = F)$p.value,
            .groups = "drop") %>% 
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>% 
  arrange(p_value) ; wilcox_so


# 1. Sorting DEB (P < 0.1)
sig_species_0.1_o = wilcox_so %>% 
  filter(p_value < 0.1) %>% 
  arrange(p_value) %>% 
  pull(Species)


# 2. Calculate mean abundance for the DEB
heat_so = so_input %>% 
  filter(Species %in% sig_species_0.1_o) %>% 
  group_by(Species, TRG_1) %>% 
  summarise(Mean_abundance = mean(abundance, na.rm = T), .groups = "drop") %>% 
  arrange(Species, TRG_1)


# 3. Conversion to wide format
heat_so_wide = heat_so %>% 
  pivot_wider(names_from = TRG_1, values_from = Mean_abundance) %>% 
  column_to_rownames("Species") %>% 
  as.matrix()


# 4. Half of min. abundance
min_abund = so_input$abundance[so_input$abundance > 0] %>% min() # 6e-05


# 5. Log10 with pseudo-abundance
heat_so_wide_log = log10(heat_so_wide + min_abund/2)


# 6. Bacteria-wise scaling
heat_so_scaled = t(scale(t(heat_so_wide_log))) # z-score per row


# 7. 별표 표시용 annotation matrix 생성
get_asterisks = function(p) {
  if(is.na(p)) return("")
  if(p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("")
}

annots = sapply(rownames(heat_so_wide_log), function(species) {
  p = wilcox_so$p_value[wilcox_so$Species == species]
  get_asterisks(p)
})

annotation_matrix = matrix("",
                           nrow = nrow(heat_so_wide_log),
                           ncol = ncol(heat_so_wide_log))

rownames(annotation_matrix) = rownames(heat_so_wide_log)
colnames(annotation_matrix) = colnames(heat_so_wide_log)

annotation_matrix[, "CR"] = annots # CR 열에만 표시


# 8. Draw
color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(300)

p5_heat_so_on = pheatmap(
  mat = heat_so_wide_log,
  cluster_rows = T,
  cluster_cols = F,
  clustering_distance_rows = "manhattan", # 또는 "euclidean", "spearman", "pearson" 등
  clustering_method_rows = "ward.D2",     # 또는 "complete", "average", "single" 등
  color = color,
  scale = "none",
  angle_col = 0,
  display_numbers = annotation_matrix,
  number_color = "black",
  fontsize_number = 14,
  border_color = "grey70"
) ; p5_heat_so_on

# ggsave("figure/05-4-3_heatmap_species_pval0.1_TRG1_ongoing.svg",
#        p5_heat_so_on, width = 6.5, height = 4)


# pseudo-abundance = 4e-05
min_abund/2

# 1. Wide-format 변환 (Species x Sample matrix)
heat_so_sample = so_input %>% 
  filter(Species %in% sig_species_0.1_o) %>% 
  group_by(Species, SampleID) %>% 
  summarise(abundance = sum(as.numeric(abundance), na.rm = T), .groups = "drop") %>% 
  pivot_wider(names_from = SampleID,
              values_from = abundance,
              values_fill = 0)


# 2. 행 이름 설정 + matrix 변환
heat_so_mat = heat_so_sample %>% 
  column_to_rownames("Species") %>% 
  as.matrix() ; dim(heat_so_mat) # 17 species in 16 samples


# 3. Log10 변환 + pseudo-abundance 추가
heat_so_log = log10(heat_so_mat + min_abund/2)


# 4. Species별 z-score scaling
heat_so_scaled = t(scale(t(heat_so_log)))
# summary(t(heat_so_scaled)) # Mean 0


# 5. 샘플 그룹 및 T_stage 정보 추출
sample_info_o = so_input %>% 
  distinct(SampleID, TRG_1, Pre_Op_Tstage_bin) %>% 
  filter(SampleID %in% colnames(heat_so_scaled)) %>% 
  mutate(Tstage = factor(Pre_Op_Tstage_bin,
                         levels = c(0, 1),
                         label = c("early", "advanced"))) %>% 
  select(SampleID, TRG_1, Tstage) %>% 
  column_to_rownames("SampleID")


# 6. CR / nonCR 나눈 후 각 그룹 안에서 clustering
heat_so_CR = heat_so_scaled[, sample_info_o$TRG_1 == "CR"]
heat_so_nonCR = heat_so_scaled[, sample_info_o$TRG_1 == "nonCR"]

col_order_so_CR = colnames(heat_so_CR)[hclust(dist(t(heat_so_CR)))$order]
col_order_so_nonCR = colnames(heat_so_nonCR)[hclust(dist(t(heat_so_nonCR)))$order]

final_so_order = c(col_order_so_CR, col_order_so_nonCR)


# 7. 열 순서 및 주석 재정렬
heat_so_final = heat_so_scaled[, final_so_order]
annotation_col = sample_info_o[final_so_order, , drop = F]

col_fun_main = circlize::colorRamp2(
  c(-2, 0, 2),
  # Blue - Light Grey (almost white) - Red
  c("#4575b4", "#F8F8F8", "#d73027")
)


# 8. 색상 설정
ann_colors = list(
  TRG_1 = c(CR = "#F7D9BC", nonCR = "#80461B"),
  Tstage = c(early = "#1b9e77", advanced = "#d95f02")
)


# 기본 ht_main
set.seed(123) ; ht_main = Heatmap(
  matrix = heat_so_final,
  name = "z-score",
  top_annotation = HeatmapAnnotation(df = annotation_col, 
                                     col = ann_colors,
                                     annotation_name_side = "left"),
  cluster_columns = F,
  cluster_rows = T,
  clustering_distance_rows = "manhattan", # 또는 "euclidean", "spearman", "pearson" 등
  clustering_method_rows = "ward.D2",     # 또는 "complete", "average", "single" 등
  
  show_column_names = F,
  show_row_names = T,
  row_names_gp = gpar(fontsize = 8),
  border = T, 
  heatmap_legend_param = list(
    title = "z-score",
    title_position = "topcenter",
    legend_direction = "horizontal",
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2")
  ),
  col = col_fun_main,
  column_title = "Species Abundance (log10 + z-score)\nGrouped CR vs nonCR, hierarchical clustering",
  rect_gp = gpar(col = "black", lwd = 0.3)  # 테두리 추가 (선 색 + 두께)
) ; ht_main

# svg("figure/05-4-4_heatmap_samples_species_ongoing.svg", width = 9, height = 7) ; draw(
#   ht_main,
#   heatmap_legend_side = "bottom",       # 또는 "left", "top"으로 조정 가능
#   annotation_legend_side = "bottom",    # T_stage, TRG_1 legend는 그대로
#   padding = unit(c(5, 10, 5, 5), "mm")  # top, right, bottom, left padding
# ) ; dev.off()


### Add a heatmap to display p-value

# 1. Species 순서 추출 (위에서 아래로 보이는 순서 기준)
species_order_o = rownames(heat_so_final)[order.dendrogram(row_dend(ht_main))]


# 2. p-value matrix
pval_data_o = wilcox_so %>% 
  filter(Species %in% species_order_o) %>% 
  mutate(
    log10_p = -log10(p_value),
    sig_star = case_when(
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>% 
  column_to_rownames("Species")

pval_mat_o = matrix(
  pval_data_o[species_order_o, "log10_p", drop = T],
  ncol = 1,
  dimnames = list(species_order_o, "-log10(p)")
)

summary(pval_mat_o)
range(pval_mat_o, na.rm = T)

sig_text_o = matrix(
  pval_data_o[species_order_o, "sig_star", drop = T],
  ncol = 1,
  dimnames = list(species_order_o, "-log10(p)")
)


# 3. 색상 스케일 (정상 순서로: 낮음 → 진청색)
col_fun_pval = colorRamp2(
  c(1.0, 1.7, 2.5),
  c("white", "#6BAED6", "#08306B")  # white → soft blue → deep blue
)


# 4. Heatmap 객체
ht_pval = Heatmap(
  matrix = pval_mat_o,
  name = "-log10(p)",
  col = col_fun_pval,
  cluster_rows = F,
  cluster_columns = F,
  width = unit(0.5, "cm"),
  show_row_names = T,
  row_names_gp = gpar(fontsize = 10),
  row_order = species_order_o,
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (!is.na(fill)) {
      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "grey60", lwd = 0.4))
    } else {
      grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "grey60", lwd = 0.4))
    }
    grid.text(sig_text_o[i, j], x = x, y = y, gp = gpar(fontsize = 10))
  }
) ; ht_pval

# svg("figure/05-4-4_heatmap_samples_species_pval_ongoing.svg",
#     width = 4, height = 5); draw(ht_pval, heatmap_legend_side = "right"); dev.off()



##### Microbial Pathway #####

# Pathway coverage
  # >= 0.5 분야 표준처럼 자주 인용됨
  # >= 0.8
  # >= 0.9 매우 보수적, false positive 최소화

# Pathways filtering
pathway_0.2       # prevalence >= 20%; 475
pathway_0.5       # prevalence >= 50%; 430
pathway_filtered  # prevalence >= 20% & mean coverage 0.04 ; 145 pathways
  # 여기서 mean coverage 0.04는 총 26개의 샘플 중 최소 한 개가 100% coverage일 경우
  # 4%의 평균 coverage가 나오는 것을 기준으로 설정

# TRG_1에 대한 wilcoxon rank-sum test 결과 저장
sig_TRG_1_pathways # 7 pathways, P<0.1
sig_TRG_1_pathways_0.2 = res_wilcox %>% 
  filter(P_TRG_1 < 0.2) %>% 
  pull(Pathway) # 11 pathways, P < 0.2

# Conversion into wider matrix: Pathway x Sample
heat_path = path %>%
  filter(Pathway %in% sig_TRG_1_pathways_0.2) %>%
  select(Pathway, SampleID, Abundance_0.5) %>%
  group_by(Pathway, SampleID) %>%
  summarise(Abundance_0.5 = mean(Abundance_0.5, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = SampleID, values_from = Abundance_0.5, values_fill = 0)

heat_path_mat = heat_path %>%
  column_to_rownames("Pathway") %>%
  as.matrix()

colSums(heat_path_mat)
rowSums(heat_path_mat) %>% as.vector()


# Log10 + z-score
min_abund_path = path$Abundance_0.5[path$Abundance_0.5 > 0] %>% min() # 4553.008
heat_path_mat_log = log10(heat_path_mat + min_abund_path/2) # 4553.008/2 = 2276.504
heat_path_mat_scaled = t(scale(t(heat_path_mat_log)))

# Sample metadata
sample_info = path %>%
  distinct(SampleID, TRG_1, Pre_Op_Tstage_bin) %>%
  filter(SampleID %in% colnames(heat_path_mat_scaled)) %>%
  mutate(Tstage = factor(Pre_Op_Tstage_bin, 
                         levels = c(0, 1), 
                         label = c("early", "advanced"))) %>%
  select(SampleID, TRG_1, Tstage) %>%
  column_to_rownames("SampleID")


# Column order by TRG_1 clustering
heat_CR = heat_path_mat_scaled[, sample_info$TRG_1 == "CR"]
heat_nonCR = heat_path_mat_scaled[, sample_info$TRG_1 == "nonCR"]

col_order_CR = colnames(heat_CR)[hclust(dist(t(heat_CR)))$order]
col_order_nonCR = colnames(heat_nonCR)[hclust(dist(t(heat_nonCR)))$order]

final_order = c(col_order_CR, col_order_nonCR)

heat_path_final_mat = heat_path_mat_scaled[, final_order]

annotation_col = sample_info[final_order, , drop = F]

# Heatmap
col_fun_main = colorRamp2(c(-2, 0, 2), c("#4575b4", "#F8F8F8", "#d73027"))
ann_colors = list(
  TRG_1 = c(CR = "#F7D9BC", nonCR = "#80461B"),
  Tstage = c(early = "#1b9e77", advanced = "#d95f02"))


set.seed(123); ht_sig = Heatmap(
  matrix = heat_path_final_mat,
  name = "z-score",
  top_annotation = HeatmapAnnotation(df = annotation_col, 
                                     col = ann_colors,
                                     annotation_name_side = "left"),
  cluster_columns = F,
  cluster_rows = T,
  
  clustering_distance_rows = "manhattan", # 또는 "euclidean", "correlatoin" 등
  clustering_method_rows = "average",   # 또는 "complete", "ward.D2", "single" 등
  
  show_column_names = F,
  show_row_names = T,
  row_names_gp = gpar(fontsize = 8),
  border = T, 
  heatmap_legend_param = list(
    title = "z-score",
    title_position = "topcenter",
    legend_direction = "horizontal",
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2")),
  col = col_fun_main,
  column_title = "Pathway Abundance (log10 + z-score)\nGrouped CR vs nonCR, hierarchical clustering",
  rect_gp = gpar(col = "black", lwd = 0.3)  # 테두리 추가 (선 색 + 두께)
) ; ht_sig

# svg("figure/04-4-5_heatmap_samples_pathway.svg", width = 12, height = 5); draw(
#   ht_sig,
#   heatmap_legend_side = "bottom",
#   annotation_legend_side = "bottom"
# ); dev.off()


### p-value tile for pathway

# 1. draw()로부터 row 순서 추출
row_index = row_order(draw(ht_sig))
pathway_order = rownames(heat_path_final_mat)[row_index]


# 2. p-value 정리
pval_df = res_wilcox %>%
  mutate(
    log10p = -log10(P_TRG_1),
    sig_star = case_when(
      P_TRG_1 < 0.05 ~ "**",
      P_TRG_1 < 0.1 ~ "*",
      TRUE ~ "")) %>%
  column_to_rownames("Pathway")


# 3. Matrix 구성
pval_mat = matrix(pval_df[pathway_order, "log10p"], ncol = 1,
                  dimnames = list(pathway_order, "-log10(p)"))

sig_text = matrix(pval_df[pathway_order, "sig_star"], ncol = 1,
                  dimnames = list(pathway_order, "-log10(p)"))


# 4. 색상 정의
col_fun_pval = colorRamp2(
  c(0.7, 1.0, 1.7),
  c("white", "#6BAED6", "#08306B")
)


# 5. 히트맵 생성
ht_pval = Heatmap(
  matrix = pval_mat,
  name = "-log10(p)",
  col = col_fun_pval,
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = T,
  row_names_gp = gpar(fontsize = 8),
  row_names_max_width = unit(15, "cm"),  
  width = unit(0.5, "cm"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    # 배경 타일
    grid.rect(x = x, y = y, width = w, height = h,
              gp = gpar(fill = fill, col = "black", lwd = 0.6))
    
    # 별표 표시
    star = sig_text[i, j]
    if (!is.na(star) && star != "") {
      grid.text(label = star, x = x, y = y,
                gp = gpar(fontsize = 10, fontface = "bold", col = "black"))
    }
  }
) ; ht_pval

# svg("figure/04-4-5_heatmap_samples_pathway_pval.svg", width = 8, height = 5) ; draw(
#   ht_pval, heatmap_legend_side = "bottom"); dev.off()


### Coverage tile for pathway

# 1. Coverage 계산
coverage_df = path %>%
  filter(Pathway %in% pathway_order) %>%
  group_by(Pathway) %>%
  summarise(mean_coverage = mean(Coverage, na.rm = TRUE), .groups = "drop") %>%
  column_to_rownames("Pathway")


# 2. matrix로 변환
coverage_mat = matrix(
  coverage_df[pathway_order, "mean_coverage"],
  ncol = 1,
  dimnames = list(pathway_order, "Coverage"))


# 3. 색상 정의 (하양 → 노랑 → 빨강)
col_fun_cov = colorRamp2(
  c(0.0, 0.5, 1.0),
  c("#FFF7BC", "#FEE08B", "#D73027")
)
  

# 4. Heatmap 생성
ht_cov = Heatmap(
  matrix = coverage_mat,
  name = "Coverage",
  col = col_fun_cov,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  width = unit(0.5, "cm"),
  row_names_max_width = unit(18, "cm"),  
  cell_fun = function(j, i, x, y, w, h, fill) {
    # 배경 타일
    grid.rect(x = x, y = y, width = w, height = h,
              gp = gpar(fill = fill, col = "black", lwd = 0.5))
  }
) ; ht_cov

# svg("figure/04-4-5_heatmap_samples_pathway_coverage.svg", width = 8, height = 5); draw(
#   ht_cov, heatmap_legend_side = "bottom"); dev.off()


### Log-fold tile for pathway

# 1. Log Fold change 계산
fc_df_CR = path %>% 
  filter(Pathway %in% pathway_order) %>% 
  group_by(Pathway, TRG_1) %>% 
  filter(TRG_1 == "CR") %>% 
  summarise(mean_abund = mean(Abundance_0.5, na.rm = T))

fc_df_nonCR = path %>% 
  filter(Pathway %in% pathway_order) %>% 
  group_by(Pathway, TRG_1) %>% 
  filter(TRG_1 == "nonCR") %>% 
  summarise(mean_abund = mean(Abundance_0.5, na.rm = T))

fc_df = fc_df_CR %>%
  left_join(fc_df_nonCR, by = "Pathway", suffix = c("_CR", "_nonCR")) %>% 
  select(-c(TRG_1_CR, TRG_1_nonCR)) %>% 
  mutate(Fold_change = (mean_abund_CR + 2276.504) / (mean_abund_nonCR + 2276.504),
         Log_Fold_change = log10(Fold_change)) %>% 
  column_to_rownames("Pathway"); rm(fc_df_CR); rm(fc_df_nonCR)


# 2. matrix로 변환
fc_mat = matrix(
  fc_df[pathway_order, "Log_Fold_change"],
  ncol = 1,
  dimnames = list(pathway_order, "Fold_change")
)


# 3. 색상 정의
col_fun_cov = colorRamp2(
  c(-0.4, 0, 0.3),
  c("#AF7AC5", "#FFFFFF", "#48C9B0")
)

# 4. Heatmap 생성
ht_fc = Heatmap(
  matrix = fc_mat,
  name = "Fold_change",
  col = col_fun_cov,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  width = unit(0.5, "cm"),
  row_names_max_width = unit(18, "cm"),  
  cell_fun = function(j, i, x, y, w, h, fill) {
    # 배경 타일
    grid.rect(x = x, y = y, width = w, height = h,
              gp = gpar(fill = fill, col = "black", lwd = 0.5))
  }
) ; ht_fc

# svg("figure/04-4-5_heatmap_samples_pathway_fold-change.svg", width = 8, height = 5); draw(
#   ht_fc, heatmap_legend_side = "bottom"); dev.off()



# Boxplot
pathway_df_0.1 = heat_path %>% 
  filter(Pathway %in% sig_TRG_1_pathways) %>% 
  pivot_longer(cols = -Pathway) %>% 
  mutate(SampleID = name) %>% 
  left_join(mb2 %>% 
              select(SampleID, TRG_1, TRG_score, Pre_Op_Tstage_bin),
            by = "SampleID") %>% 
  mutate(abundance = value)


p4_pathways_1phosphate = pathway_df_0.1 %>% 
  filter(Pathway == "PWY-7560: methylerythritol phosphate pathway II") %>% 
  ggplot(aes(TRG_1, abundance)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  geom_jitter(aes(color = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0.02) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.2)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(size = rel(1.3))) +
  labs(y = "PWY-7560: methylerythritol phosphate pathway II") ; p4_pathways_1phosphate


p4_pathways_2histidine = pathway_df_0.1 %>% 
  filter(Pathway == "HISTSYN-PWY: L-histidine biosynthesis") %>% 
  ggplot(aes(TRG_1, abundance)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  geom_jitter(aes(color = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0.02) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.2)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(size = rel(1.3))) +
  labs(y = "HISTSYN-PWY: L-histidine biosynthesis") ; p4_pathways_2histidine


p4_pathways_3isoprene = pathway_df_0.1 %>% 
  filter(Pathway == "PWY-6270: isoprene biosynthesis I") %>% 
  ggplot(aes(TRG_1, abundance)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  geom_jitter(aes(color = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0.02) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.2)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(size = rel(1.3))) +
  labs(y = "PWY-6270: isoprene biosynthesis I") ; p4_pathways_3isoprene


p4_pathways_4thiamine = pathway_df_0.1 %>% 
  filter(Pathway == "PWY-6892: thiazole component of thiamine diphosphate biosynthesis I") %>% 
  ggplot(aes(TRG_1, abundance)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  geom_jitter(aes(color = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0.02) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.2)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(size = rel(1.3))) +
  labs(y = "PWY-6892: thiazole component of thiamine diphosphate biosynthesis I") ; p4_pathways_4thiamine


p4_pathways_5sucrose = pathway_df_0.1 %>% 
  filter(Pathway == "PWY-5384: sucrose degradation IV (sucrose phosphorylase)") %>% 
  ggplot(aes(TRG_1, abundance)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  geom_jitter(aes(color = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0.02) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.2)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(size = rel(1.3))) +
  labs(y = "PWY-5384: sucrose degradation IV (sucrose phosphorylase)") ; p4_pathways_5sucrose


p4_pathways_6sucrose = pathway_df_0.1 %>% 
  filter(Pathway == "PWY-621: sucrose degradation III (sucrose invertase)") %>% 
  ggplot(aes(TRG_1, abundance)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  geom_jitter(aes(color = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0.02) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.2)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(size = rel(1.3))) +
  labs(y = "PWY-621: sucrose degradation III (sucrose invertase)") ; p4_pathways_6sucrose


p4_pathways_7rhamnose = pathway_df_0.1 %>% 
  filter(Pathway == "RHAMCAT-PWY: L-rhamnose degradation I") %>% 
  ggplot(aes(TRG_1, abundance)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  geom_jitter(aes(color = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0.02) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.2)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(size = rel(1.3))) +
  labs(y = "RHAMCAT-PWY: L-rhamnose degradation I") ; p4_pathways_7rhamnose

# ggsave("figure/04-4-6_TRG_1-7pathways.svg",
#        plot = ggarrange(p4_pathways_1phosphate,
#                         p4_pathways_2histidine,
#                         p4_pathways_3isoprene,
#                         p4_pathways_4thiamine,
#                         p4_pathways_5sucrose,
#                         p4_pathways_6sucrose,
#                         p4_pathways_7rhamnose,
#                         ncol = 7, nrow = 1),
#        width = 14, height = 4)



### In ongoing samples

# Pathways filtering
pathway_0.2_o      # Prevalence >= 20%: 486
pathway_0.5_o      # Prevalence >= 50%: 455
pathway_filtered_o # Prevalence >= 20% & mean coverage 0.04: 136 pathways

# TRG_1에 대한 wilcoxon rank-sum test 결과 저장
sig_TRG_1_pathways_o
sig_TRG_1_pathways_0.2_o = res_wilcox_o %>% 
  filter(P_TRG_1 < 0.2) %>% 
  pull(Pathway) # 15 pathways, P < 0.2

# Conversion into wider matrix: Pathway x Sample
heat_path_o = path_o %>% 
  filter(Pathway %in% sig_TRG_1_pathways_0.2_o) %>% 
  select(Pathway, SampleID, Abundance_0.5) %>% 
  group_by(Pathway, SampleID) %>% 
  summarise(Abundance_0.5 = sum(Abundance_0.5, na.rm = T), .groups = "drop") %>% 
  pivot_wider(names_from = SampleID, values_from = Abundance_0.5, values_fill = 0)

heat_path_mat_o = heat_path_o %>% 
  column_to_rownames("Pathway") %>% 
  as.matrix()

colSums(heat_path_mat_o)
rowSums(heat_path_mat_o) %>% as.vector()

# Log10 + z-score
min_abund_path = path_o$Abundance_0.5[path_o$Abundance_0.5 > 0] %>% min() # 5256.234
heat_path_mat_log_o = log10(heat_path_mat_o + min_abund_path/2) # 5256.234/2 = 2628.117
heat_path_mat_scaled_o = t(scale(t(heat_path_mat_log_o)))

# Sample metadata
sample_info_o = path_o %>% 
  distinct(SampleID, TRG_1, Pre_Op_Tstage_bin) %>% 
  filter(SampleID %in% colnames(heat_path_mat_scaled_o)) %>% 
  mutate(Tstage = factor(Pre_Op_Tstage_bin,
                         levels = c(0, 1),
                         label = c("early", "advanced"))) %>% 
  select(SampleID, TRG_1, Tstage) %>% 
  column_to_rownames("SampleID")

# Column order by TRG_1 clustering
heat_CR = heat_path_mat_scaled_o[, sample_info_o$TRG_1 == "CR"]
heat_nonCR = heat_path_mat_scaled_o[, sample_info_o$TRG_1 == "nonCR"]

col_order_CR_o = colnames(heat_CR)[hclust(dist(t(heat_CR)))$order]
col_order_nonCR_o = colnames(heat_nonCR)[hclust(dist(t(heat_nonCR)))$order]

final_order_o = c(col_order_CR_o, col_order_nonCR_o)

heat_path_final_mat_o = heat_path_mat_scaled_o[, final_order_o]

annotation_col_o = sample_info_o[final_order_o, , drop = F]

# Heatmap
col_fun_main = colorRamp2(c(-2, 0, 2), c("#4575b4", "#F8F8F8", "#d73027"))
ann_colors = list(
  TRG_1 = c(CR = "#F7D9BC", nonCR = "#80461B"),
  Tstage = c(early = "#1b9e77", advanced = "#d95f02")
)


set.seed(123) ; ht_sig = Heatmap(
  matrix = heat_path_final_mat_o,
  name = "z-score",
  top_annotation = HeatmapAnnotation(df = annotation_col_o, 
                                     col = ann_colors,
                                     annotation_name_side = "left"),
  cluster_columns = F,
  cluster_rows = T,
  
  clustering_distance_rows = "manhattan", # 또는 "euclidean", "correlatoin" 등
  clustering_method_rows = "average",   # 또는 "complete", "ward.D2", "single" 등
  
  show_column_names = F,
  show_row_names = T,
  row_names_gp = gpar(fontsize = 8),
  border = T, 
  heatmap_legend_param = list(
    title = "z-score",
    title_position = "topcenter",
    legend_direction = "horizontal",
    at = c(-2, -1, 0, 1, 2),
    labels = c("-2", "-1", "0", "1", "2")),
  col = col_fun_main,
  column_title = "Pathway Abundance (log10 + z-score)\nGrouped CR vs nonCR, hierarchical clustering",
  rect_gp = gpar(col = "black", lwd = 0.3)  # 테두리 추가 (선 색 + 두께)
) ; ht_sig

# svg("figure/05-4-5_heatmap_samples_pathway_ongoing.svg", width = 9, height = 5); draw(
#   ht_sig,
#   heatmap_legend_side = "bottom",
#   annotation_legend_side = "bottom"
# ); dev.off()


### p-value tile for pathway

# 1. draw()로부터 row 순서 추출
row_index_o = row_order(draw(ht_sig))
pathway_order_o = rownames(heat_path_final_mat_o)[row_index_o]


# 2. p-value 정리
pval_df_o = res_wilcox_o %>% 
  mutate(
    log10p = -log10(P_TRG_1),
    sig_star = case_when(
      P_TRG_1 < 0.05 ~ "**",
      P_TRG_1 < 0.1 ~ "*",
      TRUE ~ ""
    )
  ) %>% 
  column_to_rownames("Pathway")


# 3. Matrix 구성
pval_mat_o = matrix(pval_df_o[pathway_order_o, "log10p"], ncol = 1,
                    dimnames = list(pathway_order_o, "-log10(p)"))

sig_text_o = matrix(pval_df_o[pathway_order_o, "sig_star"], ncol = 1,
                    dimnames = list(pathway_order_o, "-log10(p)"))


# 4. 색상 정의
col_fun_pval = colorRamp2(
  c(0.7, 1.0, 1.7),
  c("white", "#6BAED6", "#08306B")
)


# 5. 히트맵 생성
ht_pval = Heatmap(
  matrix = pval_mat_o,
  name = "-log10(p)",
  col = col_fun_pval,
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = T,
  row_names_gp = gpar(fontsize = 8),
  row_names_max_width = unit(15, "cm"),  
  width = unit(0.5, "cm"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    # 배경 타일
    grid.rect(x = x, y = y, width = w, height = h,
              gp = gpar(fill = fill, col = "black", lwd = 0.6))
    
    # 별표 표시
    star = sig_text_o[i, j]
    if (!is.na(star) && star != "") {
      grid.text(label = star, x = x, y = y,
                gp = gpar(fontsize = 10, fontface = "bold", col = "black"))
    }
  }
) ; ht_pval

# svg("figure/05-4-5_heatmap_samples_pathway_pval_ongoing.svg", width = 6.8, height = 5) ; draw(
#   ht_pval, heatmap_legend_side = "bottom"); dev.off()


### Coverage tile for pathway

# 1. Coverage 계산
coverage_df_o = path_o %>% 
  filter(Pathway %in% pathway_order_o) %>% 
  group_by(Pathway) %>% 
  summarise(Mean_coverage = mean(Coverage, na.rm = T), .groups = "drop") %>% 
  column_to_rownames("Pathway")


# 2. Matrix로 변환
coverage_mat_o = matrix(
  coverage_df_o[pathway_order_o, "Mean_coverage"],
  ncol = 1,
  dimnames = list(pathway_order_o, "Coverage")
)


# 3. 색상 정의 (하양 → 노랑 → 빨강)
col_fun_cov = colorRamp2(
  c(0.0, 0.5, 1.0),
  c("#FFF7BC", "#FEE08B", "#D73027")
)


# 4. Heatmap 생성
ht_cov = Heatmap(
  matrix = coverage_mat_o,
  name = "Coverage",
  col = col_fun_cov,
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = T,
  row_names_gp = gpar(fontsize = 9),
  width = unit(0.5, "cm"),
  row_names_max_width = unit(18, "cm"),  
  cell_fun = function(j, i, x, y, w, h, fill) {
    # 배경 타일
    grid.rect(x = x, y = y, width = w, height = h,
              gp = gpar(fill = fill, col = "black", lwd = 0.5))
  }
) ; ht_cov

# svg("figure/05-4-5_heatmap_samples_pathway_coverage_ongoing.svg",
#     width = 8, height = 5); draw(ht_cov, heatmap_legend_side = "bottom"); dev.off()


### Log-fold tile for pathway

# 1. Log Fold change 계산
fc_df_CR_o = path_o %>% 
  filter(Pathway %in% pathway_order_o) %>% 
  group_by(Pathway, TRG_1) %>% 
  filter(TRG_1 == "CR") %>% 
  summarise(mean_abund = mean(Abundance_0.5, na.rm = T))

fc_df_nonCR_o = path_o %>% 
  filter(Pathway %in% pathway_order_o) %>% 
  group_by(Pathway, TRG_1) %>% 
  filter(TRG_1 == "nonCR") %>% 
  summarise(mean_abund = mean(Abundance_0.5, na.rm = T))

fc_df_o = fc_df_CR_o %>%
  left_join(fc_df_nonCR_o, by = "Pathway", suffix = c("_CR", "_nonCR")) %>% 
  select(-c(TRG_1_CR, TRG_1_nonCR)) %>% 
  mutate(Fold_change = (mean_abund_CR + 2628.117) / (mean_abund_nonCR + 2628.117),
         Log_Fold_change = log10(Fold_change)) %>% 
  column_to_rownames("Pathway"); rm(fc_df_CR_o); rm(fc_df_nonCR_o)


# 2. matrix로 변환
fc_mat_o = matrix(
  fc_df_o[pathway_order_o, "Log_Fold_change"],
  ncol = 1,
  dimnames = list(pathway_order_o, "Fold_change")
)


# 3. 색상 정의
col_fun_cov = colorRamp2(
  c(-0.7, 0, 0.1),
  c("#AF7AC5", "#FFFFFF", "#48C9B0")
)

# 4. Heatmap 생성
ht_fc_o = Heatmap(
  matrix = fc_mat_o,
  name = "Fold_change",
  col = col_fun_cov,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  width = unit(0.5, "cm"),
  row_names_max_width = unit(18, "cm"),  
  cell_fun = function(j, i, x, y, w, h, fill) {
    # 배경 타일
    grid.rect(x = x, y = y, width = w, height = h,
              gp = gpar(fill = fill, col = "black", lwd = 0.5))
  }
) ; ht_fc_o

# svg("figure/05-4-5_heatmap_samples_pathway_fold-change_ongoing.svg", width = 8, height = 5); draw(
#   ht_fc_o, heatmap_legend_side = "bottom"); dev.off()



# Boxplot
pathway_df_0.1_o = heat_path_o %>% 
  filter(Pathway %in% sig_TRG_1_pathways_o) %>% 
  pivot_longer(cols = -Pathway) %>% 
  mutate(SampleID = name) %>% 
  left_join(mo %>% 
              select(SampleID, TRG_1, TRG_score, Pre_Op_Tstage_bin),
            by = "SampleID") %>% 
  mutate(abundance = value)


p5_pathways_1UMP = pathway_df_0.1_o %>% 
  filter(Pathway == "PWY-5686: UMP biosynthesis I") %>% 
  ggplot(aes(TRG_1, abundance)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  geom_jitter(aes(color = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0.02) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.2)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(size = rel(1.3))) +
  labs(y = "PWY-5686: UMP biosynthesis I") ; p5_pathways_1UMP


p5_pathways_2UMP = pathway_df_0.1_o %>% 
  filter(Pathway == "PWY-7790: UMP biosynthesis II") %>% 
  ggplot(aes(TRG_1, abundance)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  geom_jitter(aes(color = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0.02) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.2)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(size = rel(1.3))) +
  labs(y = "PWY-7790: UMP biosynthesis II") ; p5_pathways_2UMP


p5_pathways_3UMP = pathway_df_0.1_o %>% 
  filter(Pathway == "PWY-7791: UMP biosynthesis III") %>% 
  ggplot(aes(TRG_1, abundance)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  geom_jitter(aes(color = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0.02) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.2)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(size = rel(1.3))) +
  labs(y = "PWY-7791: UMP biosynthesis III") ; p5_pathways_3UMP


p5_pathways_4fattyacid = pathway_df_0.1_o %>% 
  filter(Pathway == "PWY66-429: fatty acid biosynthesis initiation (mitochondria)") %>% 
  ggplot(aes(TRG_1, abundance)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  geom_jitter(aes(color = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0.02) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.2)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(size = rel(1.3))) +
  labs(y = "PWY66-429: fatty acid biosynthesis initiation (mitochondria)") ; p5_pathways_4fattyacid


p5_pathways_5coenzymeA = pathway_df_0.1_o %>% 
  filter(Pathway == "PWY-7851: coenzyme A biosynthesis II (eukaryotic)") %>% 
  ggplot(aes(TRG_1, abundance)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  geom_jitter(aes(color = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0.02) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.2)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(size = rel(1.3))) +
  labs(y = "PWY-7851: coenzyme A biosynthesis II (eukaryotic)") ; p5_pathways_5coenzymeA


p5_pathways_6adenosine = pathway_df_0.1_o %>% 
  filter(Pathway == "SALVADEHYPOX-PWY: adenosine nucleotides degradation II") %>% 
  ggplot(aes(TRG_1, abundance)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  geom_jitter(aes(color = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0.02) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.2)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(size = rel(1.3))) +
  labs(y = "SALVADEHYPOX-PWY: adenosine nucleotides degradation II") ; p5_pathways_6adenosine


p5_pathways_7sucrose = pathway_df_0.1_o %>% 
  filter(Pathway == "PWY-621: sucrose degradation III (sucrose invertase)") %>% 
  ggplot(aes(TRG_1, abundance)) +
  geom_boxplot(alpha = 0.5, outlier.alpha = 0, aes(fill = TRG_1)) +
  geom_jitter(aes(color = TRG_1)) +
  scale_fill_manual(values = c("CR" = "#F7D9BC",
                               "nonCR" = "#80461B")) +
  scale_color_manual(values = c("CR" = "#F7D9BC",
                                "nonCR" = "#80461B")) +
  stat_compare_means(comparisons = list(c("CR", "nonCR")),
                     size = 3.5, method = "wilcox",
                     hide.ns = T, tip.length = 0.02) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = rel(1.2)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(size = rel(1.3))) +
  labs(y = "PWY-621: sucrose degradation III (sucrose invertase)") ; p5_pathways_7sucrose

# ggsave("figure/05-4-6_TRG_1-7pathways_ongoing.svg",
#        plot = ggarrange(p5_pathways_1UMP,
#                         p5_pathways_2UMP,
#                         p5_pathways_3UMP,
#                         p5_pathways_4fattyacid,
#                         p5_pathways_5coenzymeA,
#                         p5_pathways_6adenosine,
#                         p5_pathways_7sucrose,
#                         ncol = 7, nrow = 1),
#        width = 14, height = 4)




# Remove unnecessary objects
rm(list = ls(pattern = "^p4_")); rm(list = ls(pattern = "^p5_"))
rm(list = ls(pattern = "^heat_")); rm(list = ls(pattern = "^ht_"))
rm(list = ls(pattern = "^annotation"))
save.image(file = "input/R_image/7-4. after-heatmap.RData")





