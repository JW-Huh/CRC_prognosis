
# Sample Filter: TNT - Before

rm(list = ls())
options(java.parameters = "-Xmx64g", stringsAsFactors = F)
setwd("C:/Users/user/Desktop/윤채빈/Data/CRC Metagenomics")

library(circlize)
library(vegan)
library(tidyverse)
library(magrittr)
library(ggpubr)
library(ggforce)
library(pROC)
library(randomForest)
library(caret)
library(ggplot2)
library(gridExtra)
library(purrr)


load("input/R_image/7-4. after-heatmap.RData")
load("input/R_image/7-5. after-random-forest-final_Stratified K-Fold.RData")



#################### 1. Pathway ####################

# Prev. >= 10%, 최소 샘플 1개에서는 100% coverage
# Wilcoxon p-value top 10%

# Before samples

pathway_0.1 = path_abund_simple$Pathway[rowSums(path_abund_simple[,-1] > 0) >= 
                                          (nrow(m) * 0.1)]

pathway_filtered_0.1 = intersect(pathway_0.1, 
                                 pathway_cover_0.04$Pathway[-c(1:2)])

path_sorted_0.1 = path %>% 
  filter(Pathway %in% pathway_filtered_0.1)

RF_path = path_sorted_0.1 %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Pathway,
              values_from = Abundance_0.5) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_path) # 146 columns (145 pathways)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(fastshap)

# 1. 데이터 준비
RF_work = RF_path
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# Top-K 변수 선택 개수
n_top_features = trunc(ncol(RF_work %>% select(-TRG_1)) * 0.1) # 14

all_p_values = sapply(colnames(RF_work)[colnames(RF_work) != "TRG_1"], function(f) {
  wilcox.test(RF_work[[f]] ~ RF_work$TRG_1)$p.value
})

top_features_overall_path_before = names(sort(all_p_values))[1:n_top_features]

cat("\n\n============================================",
    "\n [Key Biomarkers identified by Wilcoxon]",
    "\n These pathways showed the most robust difference.",
    "\n======================================================\n",
    paste(top_features_overall_path_before, collapse = "\n"))

# 시각화
par(mfrow=c(3,5), mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
for(f in top_features_overall_path_before[1:n_top_features]){
  boxplot(RF_work[[f]] ~ RF_work$TRG_1, 
          main = f, 
          col = c("#80461B", "#F7D9BC"),
          ylab = "Abundance")
}
dev.off()

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1
shap_matrix = matrix(NA, nrow = N, ncol = n_top_features)
colnames(shap_matrix) = top_features_overall_path_before

# SHAP 예측 함수 정의
shap_fun = function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, "CR"]
}

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_path_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  # (4) SHAP 계산 (필요 시 주석 해제)
  # shap_value = explain(rf_model,
  #                      X = train_data[, top_features_overall_path_before],
  #                      newdata = test_data[, top_features_overall_path_before],
  #                      pred_wrapper = shap_fun,
  #                      nsim = 100)
  # 
  # shap_matrix[test_idx, ] = as.matrix(shap_value)
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj = roc(response = true_labels, 
              predictor = k_fold_pred, 
              levels = c("nonCR", "CR"), 
              direction = "<",
              quiet = TRUE)

auc_val = auc(roc_obj) ; auc_val # 0.7939  

imp_df_before_path = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_before_path

shap_before_path = shap_matrix

# 시각화
p4_g_roc_impro = ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0("Top ", n_top_features, " Pathways Selected before CV Loop"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/04-5-1_ROC curve-pathway_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

pathway_0.1_o = path_abund_simple$Pathway[rowSums(path_abund_simple[,-1] > 0) >= 
                                            (nrow(m) * 0.1)]

pathway_filtered_0.1_o = intersect(pathway_0.1_o, 
                                   pathway_cover_0.04_o$Pathway[-c(1:2)])

path_sorted_0.1_o = path_o %>% 
  filter(Pathway %in% pathway_filtered_0.1_o)

RF_path_o = path_sorted_0.1_o %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Pathway,
              values_from = Abundance_0.5) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_path_o) # 137 columns (136 pathways)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(fastshap)

# 1. 데이터 준비
RF_work = RF_path_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# Top-K 변수 선택 개수
n_top_features = trunc(ncol(RF_work %>% select(-TRG_1)) * 0.1) # 13

all_p_values = sapply(colnames(RF_work)[colnames(RF_work) != "TRG_1"], function(f) {
  wilcox.test(RF_work[[f]] ~ RF_work$TRG_1)$p.value
})

top_features_overall_path_ongoing = names(sort(all_p_values))[1:n_top_features]

cat("\n\n============================================",
    "\n [Key Biomarkers identified by Wilcoxon]",
    "\n These pathways showed the most robust difference.",
    "\n======================================================\n",
    paste(top_features_overall_path_ongoing, collapse = "\n"))

# 시각화
par(mfrow=c(3,5), mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
for(f in top_features_overall_path_ongoing[1:n_top_features]){
  boxplot(RF_work[[f]] ~ RF_work$TRG_1, 
          main = f, 
          col = c("#80461B", "#F7D9BC"),
          ylab = "Abundance")
}
dev.off()

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1
shap_matrix_o = matrix(NA, nrow = N, ncol = n_top_features)
colnames(shap_matrix_o) = top_features_overall_path_ongoing

# SHAP 예측 함수 정의
shap_fun = function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, "CR"]
}

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_path_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  # (4) SHAP 계산 (필요 시 주석 해제)
  # shap_value = explain(rf_model,
  #                      X = train_data[, top_features_overall_path_ongoing],
  #                      newdata = test_data[, top_features_overall_path_ongoing],
  #                      pred_wrapper = shap_fun,
  #                      nsim = 100)
  # 
  # shap_matrix_o[test_idx, ] = as.matrix(shap_value)
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj = roc(response = true_labels, 
              predictor = k_fold_pred, 
              levels = c("nonCR", "CR"), 
              direction = "<",
              quiet = TRUE)

auc_val = auc(roc_obj) ; auc_val # 0.6587  

imp_df_ongoing_path = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_ongoing_path

shap_ongoing_path = shap_matrix_o

# 시각화
p5_g_roc_impro = ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0("Top ", n_top_features, " Pathways Selected before CV Loop"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/04-5-1_ROC curve-pathway_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)



#################### 2. Strain ####################

# Prev. >= 10%,
# Wilcoxon p-value top 10%

# Before samples

tb_0.1 = tb %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Strain) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.1) %>% 
  pull(Strain) %>% 
  sort()

tb_input_0.1 = tb %>% 
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

RF_strain = tb_input_0.1 %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Strain,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_strain) # 498 columns (497 strains)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(fastshap)

# 1. 데이터 준비
RF_work = RF_strain
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# Top-K 변수 선택 개수
n_top_features = trunc(ncol(RF_work %>% select(-TRG_1)) * 0.1) # 49

all_p_values = sapply(colnames(RF_work)[colnames(RF_work) != "TRG_1"], function(f) {
  wilcox.test(RF_work[[f]] ~ RF_work$TRG_1)$p.value
})

top_features_overall_strain_before = names(sort(all_p_values))[1:n_top_features]

cat("\n\n============================================",
    "\n [Key Biomarkers identified by Wilcoxon]",
    "\n These strains showed the most robust difference.",
    "\n======================================================\n",
    paste(top_features_overall_strain_before, collapse = "\n"))

# 시각화
par(mfrow=c(5,10), mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
for(f in top_features_overall_strain_before[1:n_top_features]){
  boxplot(RF_work[[f]] ~ RF_work$TRG_1, 
          main = f, 
          col = c("#80461B", "#F7D9BC"),
          ylab = "Abundance")
}
dev.off()

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1
shap_matrix = matrix(NA, nrow = N, ncol = n_top_features)
colnames(shap_matrix) = top_features_overall_strain_before

# SHAP 예측 함수 정의
shap_fun = function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, "CR"]
}

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_strain_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  # (4) SHAP 계산 (필요 시 주석 해제)
  # shap_value = explain(rf_model,
  #                      X = train_data[, top_features_overall_strain_before],
  #                      newdata = test_data[, top_features_overall_strain_before],
  #                      pred_wrapper = shap_fun,
  #                      nsim = 100)
  # 
  # shap_matrix[test_idx, ] = as.matrix(shap_value)
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj = roc(response = true_labels, 
              predictor = k_fold_pred, 
              levels = c("nonCR", "CR"), 
              direction = "<",
              quiet = TRUE)

auc_val = auc(roc_obj) ; auc_val # 0.8788 (SHAP 계산 전 구한 값, 다른 taxa level과 동일한 조건)

imp_df_before_strain = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_before_strain

shap_before_strain = shap_matrix

# 시각화
p4_g_roc_impro = ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0("Top ", n_top_features, " Strains Selected before CV Loop"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/04-5-2_ROC curve-strain_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

to_0.1 = to %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Strain) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.1) %>% 
  pull(Strain) %>% 
  sort()

to_input_0.1 = to %>% 
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

RF_strain_o = to_input_0.1 %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Strain,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_strain_o) # 306 columns (305 strains)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(fastshap)

# 1. 데이터 준비
RF_work = RF_strain_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# Top-K 변수 선택 개수
n_top_features = trunc(ncol(RF_work %>% select(-TRG_1)) * 0.1) # 30

all_p_values = sapply(colnames(RF_work)[colnames(RF_work) != "TRG_1"], function(f) {
  wilcox.test(RF_work[[f]] ~ RF_work$TRG_1)$p.value
})

top_features_overall_strain_ongoing = names(sort(all_p_values))[1:n_top_features]

cat("\n\n============================================",
    "\n [Key Biomarkers identified by Wilcoxon]",
    "\n These strains showed the most robust difference.",
    "\n======================================================\n",
    paste(top_features_overall_strain_ongoing, collapse = "\n"))

# 시각화
par(mfrow=c(5,6), mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
for(f in top_features_overall_strain_ongoing[1:n_top_features]){
  boxplot(RF_work[[f]] ~ RF_work$TRG_1, 
          main = f, 
          col = c("#80461B", "#F7D9BC"),
          ylab = "Abundance")
}
dev.off()

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1
shap_matrix_o = matrix(NA, nrow = N, ncol = n_top_features)
colnames(shap_matrix_o) = top_features_overall_strain_ongoing

# SHAP 예측 함수 정의
shap_fun = function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, "CR"]
}

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_strain_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  # (4) SHAP 계산 (필요 시 주석 해제)
  # shap_value = explain(rf_model,
  #                      X = train_data[, top_features_overall_strain_ongoing],
  #                      newdata = test_data[, top_features_overall_strain_ongoing],
  #                      pred_wrapper = shap_fun,
  #                      nsim = 100)
  # 
  # shap_matrix_o[test_idx, ] = as.matrix(shap_value)
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj = roc(response = true_labels, 
              predictor = k_fold_pred, 
              levels = c("nonCR", "CR"), 
              direction = "<",
              quiet = TRUE)

auc_val = auc(roc_obj) ; auc_val # 1 (SHAP 계산 전 구한 값, 다른 taxa level과 동일한 조건)

imp_df_ongoing_strain = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_ongoing_strain

shap_ongoing_strain = shap_matrix_o

# 시각화
p5_g_roc_impro = ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0("Top ", n_top_features, " Strains Selected before CV Loop"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/04-5-2_ROC curve-strain_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)



#################### 3. Species ####################

# Prev. >= 10%,

# Wilcoxon p-value top 10%

# Before samples

sb_0.1 = sb %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Species) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.1) %>% 
  pull(Species) %>% 
  sort()

sb_input_0.1 = sb %>% 
  filter(Species %in% sb_0.1) %>% 
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

RF_species = sb_input_0.1 %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Species,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_species) # 460 columns (459 species)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(fastshap)

# 1. 데이터 준비
RF_work = RF_species
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# Top-K 변수 선택 개수
n_top_features = trunc(ncol(RF_work %>% select(-TRG_1)) * 0.1) # 45

all_p_values = sapply(colnames(RF_work)[colnames(RF_work) != "TRG_1"], function(f) {
  wilcox.test(RF_work[[f]] ~ RF_work$TRG_1)$p.value
})

top_features_overall_species_before = names(sort(all_p_values))[1:n_top_features]

cat("\n\n============================================",
    "\n [Key Biomarkers identified by Wilcoxon]",
    "\n These species showed the most robust difference.",
    "\n======================================================\n",
    paste(top_features_overall_species_before, collapse = "\n"))

# 시각화
par(mfrow=c(5,10), mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
for(f in top_features_overall_species_before[1:n_top_features]){
  boxplot(RF_work[[f]] ~ RF_work$TRG_1, 
          main = f, 
          col = c("#80461B", "#F7D9BC"),
          ylab = "Abundance")
}
dev.off()

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1
shap_matrix = matrix(NA, nrow = N, ncol = n_top_features)
colnames(shap_matrix) = top_features_overall_species_before

# SHAP 예측 함수 정의
shap_fun = function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, "CR"]
}

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_species_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  # (4) SHAP 계산 (필요 시 주석 해제)
  # shap_value = explain(rf_model,
  #                      X = train_data[, top_features_overall_species_before],
  #                      newdata = test_data[, top_features_overall_species_before],
  #                      pred_wrapper = shap_fun,
  #                      nsim = 100)
  # 
  # shap_matrix[test_idx, ] = as.matrix(shap_value)
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj = roc(response = true_labels, 
              predictor = k_fold_pred, 
              levels = c("nonCR", "CR"), 
              direction = "<",
              quiet = TRUE)

auc_val = auc(roc_obj) ; auc_val # 0.897  

imp_df_before_species = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_before_species

shap_before_species = shap_matrix

# 시각화
p4_g_roc_impro = ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0("Top ", n_top_features, " Species Selected before CV Loop"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/04-5-3_ROC curve-species_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

so_0.1 = so %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Species) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.1) %>% 
  pull(Species) %>% 
  sort()

so_input_0.1 = so %>% 
  filter(Species %in% so_0.1) %>% 
  pivot_longer(cols = -Species,
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

RF_species_o = so_input_0.1 %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Species,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_species_o) # 289 columns (288 species)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(fastshap)

# 1. 데이터 준비
RF_work = RF_species_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# Top-K 변수 선택 개수
n_top_features = trunc(ncol(RF_work %>% select(-TRG_1)) * 0.1) # 28

all_p_values = sapply(colnames(RF_work)[colnames(RF_work) != "TRG_1"], function(f) {
  wilcox.test(RF_work[[f]] ~ RF_work$TRG_1)$p.value
})

top_features_overall_species_ongoing = names(sort(all_p_values))[1:n_top_features]

cat("\n\n============================================",
    "\n [Key Biomarkers identified by Wilcoxon]",
    "\n These species showed the most robust difference.",
    "\n======================================================\n",
    paste(top_features_overall_species_ongoing, collapse = "\n"))

# 시각화
par(mfrow=c(5,6), mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
for(f in top_features_overall_species_ongoing[1:n_top_features]){
  boxplot(RF_work[[f]] ~ RF_work$TRG_1, 
          main = f, 
          col = c("#80461B", "#F7D9BC"),
          ylab = "Abundance")
}
dev.off()

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1
shap_matrix_o = matrix(NA, nrow = N, ncol = n_top_features)
colnames(shap_matrix_o) = top_features_overall_species_ongoing

# SHAP 예측 함수 정의
shap_fun = function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, "CR"]
}

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_species_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  # (4) SHAP 계산 (필요 시 주석 해제)
  # shap_value = explain(rf_model,
  #                      X = train_data[, top_features_overall_species_ongoing],
  #                      newdata = test_data[, top_features_overall_species_ongoing],
  #                      pred_wrapper = shap_fun,
  #                      nsim = 100)
  # 
  # shap_matrix_o[test_idx, ] = as.matrix(shap_value)
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj = roc(response = true_labels, 
              predictor = k_fold_pred, 
              levels = c("nonCR", "CR"), 
              direction = "<",
              quiet = TRUE)

auc_val = auc(roc_obj) ; auc_val # 1  

imp_df_ongoing_species = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_ongoing_species

shap_ongoing_species = shap_matrix_o

# 시각화
p5_g_roc_impro = ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0("Top ", n_top_features, " Species Selected before CV Loop"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/04-5-3_ROC curve-species_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)



#################### 4. Genus ####################

# Prev. >= 10%,
# Wilcoxon p-value top 10%

# Before samples

gb_0.1 = gb %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Genus) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.1) %>% 
  pull(Genus) %>% 
  sort()

gb_input_0.1 = gb %>% 
  filter(Genus %in% gb_0.1) %>% 
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

RF_genus = gb_input_0.1 %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Genus,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_genus) # 264 columns (263 genera)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_genus
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# Top-K 변수 선택 개수
n_top_features = trunc(ncol(RF_work %>% select(-TRG_1)) * 0.1) # 26

all_p_values = sapply(colnames(RF_work)[colnames(RF_work) != "TRG_1"], function(f) {
  wilcox.test(RF_work[[f]] ~ RF_work$TRG_1)$p.value
})

top_features_overall_genus_before = names(sort(all_p_values))[1:n_top_features]

cat("\n\n============================================",
    "\n [Key Biomarkers identified by Wilcoxon]",
    "\n These genera showed the most robust difference.",
    "\n======================================================\n",
    paste(top_features_overall_genus_before, collapse = "\n"))

# 시각화
par(mfrow=c(6,5), mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
for(f in top_features_overall_genus_before[1:n_top_features]){
  boxplot(RF_work[[f]] ~ RF_work$TRG_1, 
          main = f, 
          col = c("#80461B", "#F7D9BC"),
          ylab = "Abundance")
}
dev.off()

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_genus_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj = roc(response = true_labels, 
              predictor = k_fold_pred, 
              levels = c("nonCR", "CR"), 
              direction = "<",
              quiet = TRUE)

auc_val = auc(roc_obj) ; auc_val # 0.8242  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p4_g_roc_impro = ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0("Top ", n_top_features, " Genera Selected before CV Loop"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/04-5-4_ROC curve-genus_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

go_0.1 = go %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Genus) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.1) %>% 
  pull(Genus) %>% 
  sort()

go_input_0.1 = go %>% 
  filter(Genus %in% go_0.1) %>% 
  pivot_longer(cols = -Genus,
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

RF_genus_o = go_input_0.1 %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Genus,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_genus_o) # 175 columns (174 genera)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_genus_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# Top-K 변수 선택 개수
n_top_features = trunc(ncol(RF_work %>% select(-TRG_1)) * 0.1) # 17

all_p_values = sapply(colnames(RF_work)[colnames(RF_work) != "TRG_1"], function(f) {
  wilcox.test(RF_work[[f]] ~ RF_work$TRG_1)$p.value
})

top_features_overall_genus_ongoing = names(sort(all_p_values))[1:n_top_features]

cat("\n\n============================================",
    "\n [Key Biomarkers identified by Wilcoxon]",
    "\n These genera showed the most robust difference.",
    "\n======================================================\n",
    paste(top_features_overall_genus_ongoing, collapse = "\n"))

# 시각화
par(mfrow=c(4,5), mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
for(f in top_features_overall_genus_ongoing[1:n_top_features]){
  boxplot(RF_work[[f]] ~ RF_work$TRG_1, 
          main = f, 
          col = c("#80461B", "#F7D9BC"),
          ylab = "Abundance")
}
dev.off()

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_genus_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj = roc(response = true_labels, 
              predictor = k_fold_pred, 
              levels = c("nonCR", "CR"), 
              direction = "<",
              quiet = TRUE)

auc_val = auc(roc_obj) ; auc_val # 0.873  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p5_g_roc_impro = ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0("Top ", n_top_features, " Genera Selected before CV Loop"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/04-5-4_ROC curve-genus_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)



#################### 5. Family ####################

# Prev. >= 10%,
# Wilcoxon p-value top 10%

# Before samples

fb_0.1 = fb %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Family) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.1) %>% 
  pull(Family) %>% 
  sort()

fb_input_0.1 = fb %>% 
  filter(Family %in% fb_0.1) %>% 
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

RF_family = fb_input_0.1 %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Family,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_family) # 103 columns (102 families)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_family
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# Top-K 변수 선택 개수
n_top_features = trunc(ncol(RF_work %>% select(-TRG_1)) * 0.1) # 10

all_p_values = sapply(colnames(RF_work)[colnames(RF_work) != "TRG_1"], function(f) {
  wilcox.test(RF_work[[f]] ~ RF_work$TRG_1)$p.value
})

top_features_overall_family_before = names(sort(all_p_values))[1:n_top_features]

cat("\n\n============================================",
    "\n [Key Biomarkers identified by Wilcoxon]",
    "\n These families showed the most robust difference.",
    "\n======================================================\n",
    paste(top_features_overall_family_before, collapse = "\n"))

# 시각화
par(mfrow=c(2,5), mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
for(f in top_features_overall_family_before[1:n_top_features]){
  boxplot(RF_work[[f]] ~ RF_work$TRG_1, 
          main = f, 
          col = c("#80461B", "#F7D9BC"),
          ylab = "Abundance")
}
dev.off()

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_family_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj = roc(response = true_labels, 
              predictor = k_fold_pred, 
              levels = c("nonCR", "CR"), 
              direction = "<",
              quiet = TRUE)

auc_val = auc(roc_obj) ; auc_val # 0.5879  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p4_g_roc_impro = ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0("Top ", n_top_features, " Families Selected before CV Loop"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/04-5-5_ROC curve-family_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

fo_0.1 = fo %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Family) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.1) %>% 
  pull(Family) %>% 
  sort()

fo_input_0.1 = fo %>% 
  filter(Family %in% fo_0.1) %>% 
  pivot_longer(cols = -Family,
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

RF_family_o = fo_input_0.1 %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Family,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_family_o) # 74 columns (73 families)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_family_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# Top-K 변수 선택 개수
n_top_features = trunc(ncol(RF_work %>% select(-TRG_1)) * 0.1) # 7

all_p_values = sapply(colnames(RF_work)[colnames(RF_work) != "TRG_1"], function(f) {
  wilcox.test(RF_work[[f]] ~ RF_work$TRG_1)$p.value
})

top_features_overall_family_ongoing = names(sort(all_p_values))[1:n_top_features]

cat("\n\n============================================",
    "\n [Key Biomarkers identified by Wilcoxon]",
    "\n These families showed the most robust difference.",
    "\n======================================================\n",
    paste(top_features_overall_family_ongoing, collapse = "\n"))

# 시각화
par(mfrow=c(2,5), mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
for(f in top_features_overall_family_ongoing[1:n_top_features]){
  boxplot(RF_work[[f]] ~ RF_work$TRG_1, 
          main = f, 
          col = c("#80461B", "#F7D9BC"),
          ylab = "Abundance")
}
dev.off()

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_family_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj = roc(response = true_labels, 
              predictor = k_fold_pred, 
              levels = c("nonCR", "CR"), 
              direction = "<",
              quiet = TRUE)

auc_val = auc(roc_obj) ; auc_val # 0.746  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p5_g_roc_impro = ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0("Top ", n_top_features, " Families Selected before CV Loop"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/04-5-5_ROC curve-family_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)



#################### 6. Order ####################

# Prev. >= 10%,
# Wilcoxon p-value top 10%

# Before samples
o = read.csv("input/order.csv") %>% 
  select(c("Order", m$SampleID)) ; head(o)

ob = o %>% select(c("Order", mb$SampleID))

ob_0.1 = ob %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Order) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.1) %>% 
  pull(Order) %>% 
  sort()

ob_input_0.1 = ob %>% 
  filter(Order %in% ob_0.1) %>% 
  pivot_longer(cols = -Order,
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

RF_order = ob_input_0.1 %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Order,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_order) # 80 columns (79 orders)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_order
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# Top-K 변수 선택 개수
n_top_features = trunc(ncol(RF_work %>% select(-TRG_1)) * 0.1) # 7

all_p_values = sapply(colnames(RF_work)[colnames(RF_work) != "TRG_1"], function(f) {
  wilcox.test(RF_work[[f]] ~ RF_work$TRG_1)$p.value
})

top_features_overall_order_before = names(sort(all_p_values))[1:n_top_features]

cat("\n\n============================================",
    "\n [Key Biomarkers identified by Wilcoxon]",
    "\n These orders showed the most robust difference.",
    "\n======================================================\n",
    paste(top_features_overall_order_before, collapse = "\n"))

# 시각화
par(mfrow=c(2,5), mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
for(f in top_features_overall_order_before[1:n_top_features]){
  boxplot(RF_work[[f]] ~ RF_work$TRG_1, 
          main = f, 
          col = c("#80461B", "#F7D9BC"),
          ylab = "Abundance")
}
dev.off()

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_order_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj = roc(response = true_labels, 
              predictor = k_fold_pred, 
              levels = c("nonCR", "CR"), 
              direction = "<",
              quiet = TRUE)

auc_val = auc(roc_obj) ; auc_val # 0.6727  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p4_g_roc_impro = ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0("Top ", n_top_features, " Orders Selected before CV Loop"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/04-5-6_ROC curve-order_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

oo = o %>% select(c("Order", mo$SampleID))

oo_0.1 = oo %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Order) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.1) %>% 
  pull(Order) %>% 
  sort()

oo_input_0.1 = oo %>% 
  filter(Order %in% oo_0.1) %>% 
  pivot_longer(cols = -Order,
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

RF_order_o = oo_input_0.1 %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Order,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  select(-SampleID) ; dim(RF_order_o) # 54 columns (53 orders)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_order_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# Top-K 변수 선택 개수
n_top_features = trunc(ncol(RF_work %>% select(-TRG_1)) * 0.1) # 5

all_p_values = sapply(colnames(RF_work)[colnames(RF_work) != "TRG_1"], function(f) {
  wilcox.test(RF_work[[f]] ~ RF_work$TRG_1)$p.value
})

top_features_overall_order_ongoing = names(sort(all_p_values))[1:n_top_features]

cat("\n\n============================================",
    "\n [Key Biomarkers identified by Wilcoxon]",
    "\n These orders showed the most robust difference.",
    "\n======================================================\n",
    paste(top_features_overall_order_ongoing, collapse = "\n"))

# 시각화
par(mfrow=c(1,5), mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
for(f in top_features_overall_order_ongoing[1:n_top_features]){
  boxplot(RF_work[[f]] ~ RF_work$TRG_1, 
          main = f, 
          col = c("#80461B", "#F7D9BC"),
          ylab = "Abundance")
}
dev.off()

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_order_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj = roc(response = true_labels, 
              predictor = k_fold_pred, 
              levels = c("nonCR", "CR"), 
              direction = "<",
              quiet = TRUE)

auc_val = auc(roc_obj) ; auc_val # 0.746  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p5_g_roc_impro = ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0("Top ", n_top_features, " Orders Selected before CV Loop"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/04-5-6_ROC curve-order_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)



#################### 7. Class ####################

# Prev. >= 10%,
# Wilcoxon p-value top 10%

# Before samples

c = read.csv("input/class.csv") %>% 
  select(c("Class", m$SampleID)) ; head(c)

cb = c %>% select(c("Class", mb$SampleID))

cb_0.1 = cb %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Class) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.1) %>% 
  pull(Class) %>% 
  sort() ; cb_0.1

cb_input_0.1 = cb %>% 
  filter(Class %in% cb_0.1) %>% 
  pivot_longer(cols = -Class,
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

RF_class = cb_input_0.1 %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Class,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_class) # 70 columns (69 classes)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_class
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# Top-K 변수 선택 개수
n_top_features = trunc(ncol(RF_work %>% select(-TRG_1)) * 0.1) # 6

all_p_values = sapply(colnames(RF_work)[colnames(RF_work) != "TRG_1"], function(f) {
  wilcox.test(RF_work[[f]] ~ RF_work$TRG_1)$p.value
})

top_features_overall_class_before = names(sort(all_p_values))[1:n_top_features]

cat("\n\n============================================",
    "\n [Key Biomarkers identified by Wilcoxon]",
    "\n These classes showed the most robust difference.",
    "\n======================================================\n",
    paste(top_features_overall_class_before, collapse = "\n"))

# 시각화
par(mfrow=c(2,3), mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
for(f in top_features_overall_class_before[1:n_top_features]){
  boxplot(RF_work[[f]] ~ RF_work$TRG_1, 
          main = f, 
          col = c("#80461B", "#F7D9BC"),
          ylab = "Abundance")
}
dev.off()

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_class_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj = roc(response = true_labels, 
              predictor = k_fold_pred, 
              levels = c("nonCR", "CR"), 
              direction = "<",
              quiet = TRUE)

auc_val = auc(roc_obj) ; auc_val # 0.7212  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p4_g_roc_impro = ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0("Top ", n_top_features, " Classes Selected before CV Loop"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/04-5-7_ROC curve-class_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

co = c %>% select(c("Class", mo$SampleID))

co_0.1 = co %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Class) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.1) %>% 
  pull(Class) %>% 
  sort()

co_input_0.1 = co %>% 
  filter(Class %in% co_0.1) %>% 
  pivot_longer(cols = -Class,
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

RF_class_o = co_input_0.1 %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Class,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  select(-SampleID) ; dim(RF_class_o) # 45 columns (44 classes)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_class_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# Top-K 변수 선택 개수
n_top_features = trunc(ncol(RF_work %>% select(-TRG_1)) * 0.1) # 4

all_p_values = sapply(colnames(RF_work)[colnames(RF_work) != "TRG_1"], function(f) {
  wilcox.test(RF_work[[f]] ~ RF_work$TRG_1)$p.value
})

top_features_overall_class_ongoing = names(sort(all_p_values))[1:n_top_features]

cat("\n\n============================================",
    "\n [Key Biomarkers identified by Wilcoxon]",
    "\n These classes showed the most robust difference.",
    "\n======================================================\n",
    paste(top_features_overall_class_ongoing, collapse = "\n"))

# 시각화
par(mfrow=c(1,4), mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
for(f in top_features_overall_class_ongoing[1:n_top_features]){
  boxplot(RF_work[[f]] ~ RF_work$TRG_1, 
          main = f, 
          col = c("#80461B", "#F7D9BC"),
          ylab = "Abundance")
}
dev.off()

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_class_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj = roc(response = true_labels, 
              predictor = k_fold_pred, 
              levels = c("nonCR", "CR"), 
              direction = "<",
              quiet = TRUE)

auc_val = auc(roc_obj) ; auc_val # 0.6825  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p5_g_roc_impro = ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0("Top ", n_top_features, " Classes Selected before CV Loop"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/04-5-7_ROC curve-class_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)



#################### 8. Phylum ####################

# Prev. >= 10%,
# Wilcoxon p-value top 10%

# Before samples

p = read.csv("input/phylum.csv") %>% 
  select(c("clade_name", m$SampleID)) %>% 
  mutate(Phylum = str_extract(clade_name, "(?<=p__).*")) %>% 
  select(-clade_name) ; head(p)

pb = p %>% select(c("Phylum", mb$SampleID))

pb_0.1 = pb %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Phylum) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.1) %>% 
  pull(Phylum) %>% 
  sort() ; pb_0.1

pb_input_0.1 = pb %>% 
  filter(Phylum %in% pb_0.1) %>% 
  pivot_longer(cols = -Phylum,
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

RF_phylum = pb_input_0.1 %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Phylum,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_phylum) # 11 columns (10 phyla)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_phylum
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# Top-K 변수 선택 개수
n_top_features = trunc(ncol(RF_work %>% select(-TRG_1)) * 0.1) # 1

all_p_values = sapply(colnames(RF_work)[colnames(RF_work) != "TRG_1"], function(f) {
  wilcox.test(RF_work[[f]] ~ RF_work$TRG_1)$p.value
})

top_features_overall_phylum_before = names(sort(all_p_values))[1:n_top_features]

cat("\n\n============================================",
    "\n [Key Biomarkers identified by Wilcoxon]",
    "\n These phyla showed the most robust difference.",
    "\n======================================================\n",
    paste(top_features_overall_phylum_before, collapse = "\n"))

# 시각화
par(mfrow=c(1,1), mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
for(f in top_features_overall_phylum_before[1:n_top_features]){
  boxplot(RF_work[[f]] ~ RF_work$TRG_1, 
          main = f, 
          col = c("#80461B", "#F7D9BC"),
          ylab = "Abundance")
}
dev.off()

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_phylum_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj = roc(response = true_labels, 
              predictor = k_fold_pred, 
              levels = c("nonCR", "CR"), 
              direction = "<",
              quiet = TRUE)

auc_val = auc(roc_obj) ; auc_val # 0.6848  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p4_g_roc_impro = ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0("Top ", n_top_features, " Phyla Selected before CV Loop"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/04-5-8_ROC curve-phylum_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

po = p %>% select(c("Phylum", mo$SampleID))

po_0.1 = po %>% 
  rowwise() %>% 
  mutate(prevalence = sum(c_across(-Phylum) > 0)) %>% 
  ungroup() %>% 
  filter(prevalence >= nrow(m) * 0.1) %>% 
  pull(Phylum) %>% 
  sort()

po_input_0.1 = po %>% 
  filter(Phylum %in% po_0.1) %>% 
  pivot_longer(cols = -Phylum,
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

RF_phylum_o = po_input_0.1 %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Phylum,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  select(-SampleID) ; dim(RF_phylum_o) # 10 columns (9 phyla)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_phylum_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# Top-K 변수 선택 개수
n_top_features = round(ncol(RF_work %>% select(-TRG_1)) * 0.1) # 1 (0.9개라 반올림하여 1개 선택)

all_p_values = sapply(colnames(RF_work)[colnames(RF_work) != "TRG_1"], function(f) {
  wilcox.test(RF_work[[f]] ~ RF_work$TRG_1)$p.value
})

top_features_overall_phylum_ongoing = names(sort(all_p_values))[1:n_top_features]

cat("\n\n============================================",
    "\n [Key Biomarkers identified by Wilcoxon]",
    "\n These phyla showed the most robust difference.",
    "\n======================================================\n",
    paste(top_features_overall_phylum_ongoing, collapse = "\n"))

# 시각화
par(mfrow=c(1,1), mar=c(2, 2, 2, 1), oma=c(1, 1, 1, 1))
for(f in top_features_overall_phylum_ongoing[1:n_top_features]){
  boxplot(RF_work[[f]] ~ RF_work$TRG_1, 
          main = f, 
          col = c("#80461B", "#F7D9BC"),
          ylab = "Abundance")
}
dev.off()

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_phylum_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj = roc(response = true_labels, 
              predictor = k_fold_pred, 
              levels = c("nonCR", "CR"), 
              direction = "<",
              quiet = TRUE)

auc_val = auc(roc_obj) ; auc_val # 0.6508  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p5_g_roc_impro = ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0("Top ", n_top_features, " Phyla Selected before CV Loop"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/04-5-8_ROC curve-phylum_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)



#################### 9. Clinical indicator ####################

# Before samples

RF_meta = mb2 %>% 
  select(SampleID, TRG_1, Pre_Op_Tstage, Pre_Op_Nstage, BMI, CEA) %>% 
  mutate(TRG_1 = factor(TRG_1),
         Pre_Op_Tstage = factor(Pre_Op_Tstage),
         Pre_Op_Nstage = factor(Pre_Op_Nstage),
         CEA_bin = ifelse(CEA < 5.0, "<5", "5+")) %>% 
  select(-c(BMI, CEA)) %>% 
  column_to_rownames("SampleID")


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_meta
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
top_features_overall_meta_before = colnames(RF_meta)[colnames(RF_meta) != "TRG_1"]

cat("\n\n============================================",
    "\n [Key Biomarkers identified by Wilcoxon]",
    "\n These clniical indicators showed the most robust difference.",
    "\n======================================================\n",
    paste(top_features_overall_meta_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_meta_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj = roc(response = true_labels, 
              predictor = k_fold_pred, 
              levels = c("nonCR", "CR"), 
              direction = "<",
              quiet = TRUE)

auc_val = auc(roc_obj) ; auc_val # 0.7121  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p4_g_roc_impro = ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/04-5-9_ROC curve-clinical indicator_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

RF_meta_o = mo %>% 
  select(SampleID, TRG_1, Pre_Op_Tstage, Pre_Op_Nstage, BMI, CEA) %>% 
  mutate(TRG_1 = factor(TRG_1),
         Pre_Op_Tstage = factor(Pre_Op_Tstage),
         Pre_Op_Nstage = factor(Pre_Op_Nstage),
         CEA_bin = ifelse(CEA < 5.0, "<5", "5+")) %>% 
  select(-c(BMI, CEA)) %>% 
  column_to_rownames("SampleID")


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_meta_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
top_features_overall_meta_ongoing = colnames(RF_meta)[colnames(RF_meta) != "TRG_1"]

cat("\n\n============================================",
    "\n [Key Biomarkers identified by Wilcoxon]",
    "\n These clniical indicators showed the most robust difference.",
    "\n======================================================\n",
    paste(top_features_overall_meta_ongoing, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_meta_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj = roc(response = true_labels, 
              predictor = k_fold_pred, 
              levels = c("nonCR", "CR"), 
              direction = "<",
              quiet = TRUE)

auc_val = auc(roc_obj) ; auc_val # 0.881  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p5_g_roc_impro = ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/04-5-9_ROC curve-clinical indicator_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)



### Performance comparison of random forest models

# Before

RF_model_AUC = readxl::read_xlsx("RF_results.xlsx", sheet = "Stratified K-Fold", skip = 2)

p4_RF_all_level_before = RF_model_AUC[c(2:8, 1, 9), 1:2] %>% 
  mutate(Variable = factor(Variable, levels = rev(Variable))) %>% 
  ggplot(aes(x = Variable, y = `Baseline - AUC`, fill = "#E76F51")) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  geom_text(aes(label = `Baseline - AUC`),
            position = position_stack(vjust = 1),
            color = "black", size = 4) +
  labs(title = "AUC - Baseline", x = NULL, y = "Mean AUC") +
  theme_classic(base_size = 13) +
  theme(axis.text.y = element_text(size = rel(1.1), vjust = 0.5, 
                                   color= "black", hjust = 1),
        legend.position = "none") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  coord_flip(ylim = c(0.3, NA)) ; p4_RF_all_level_before

# ggsave("figure/04-5-0_AUC barplot_before.svg",
#        p4_RF_all_level_before, width = 4, height = 5)

# Ongoing

p5_RF_all_level_ongoing = RF_model_AUC[c(2:8, 1, 9), c(1, 3)] %>% 
  mutate(Variable = factor(Variable, levels = rev(Variable))) %>% 
  ggplot(aes(x = Variable, y = `After RT - AUC`, fill = "#E76F51")) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  geom_text(aes(label = `After RT - AUC`),
            position = position_stack(vjust = 1),
            color = "black", size = 4) +
  labs(title = "AUC - After RT", x = NULL, y = "Mean AUC") +
  theme_classic(base_size = 13) +
  theme(axis.text.y = element_text(size = rel(1.1), vjust = 0.5, 
                                   color= "black", hjust = 1),
        legend.position = "none") + 
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  coord_flip(ylim = c(0.3, NA)) ; p5_RF_all_level_ongoing

# ggsave("figure/04-5-0_AUC barplot_ongoing.svg",
#        p5_RF_all_level_ongoing, width = 4, height = 5)



#################### 10. Combined Feature ####################

##### Strain Before + Ongoing

# Before samples

all_strain = union(top_features_overall_strain_before,
                   top_features_overall_strain_ongoing) %>% 
  gsub("\\.", "|", .)

RF_strain_before_ongoing = RF_strain %>% 
  select(TRG_1, any_of(all_strain)) %>% 
  mutate(TRG_1 = factor(TRG_1))


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_strain_before_ongoing
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
top_features_overall_strain_all = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [All features selected from either 'Before' or 'Ongoing' samples]",
    "\n======================================================\n",
    paste(top_features_overall_strain_all, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_strain_all, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.8606  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(top_features_overall_strain_all), 
                         " All Strains Selected from Either 'Before' or 'Ongoing' Samples"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/04-5-10_ROC curve-strain_before + ongoing.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)



##### Strain + Clinical Indicator

# Before samples

mb2_selected = mb2 %>% 
  select(SampleID, Pre_Op_Tstage, Pre_Op_Nstage, CEA_bin) %>% 
  column_to_rownames("SampleID")

RF_strain_meta = RF_strain %>% 
  select(TRG_1, any_of(gsub("\\.", "|", top_features_overall_strain_before))) %>% 
  rownames_to_column("SampleID") %>% 
  inner_join(mb2_selected %>% 
               rownames_to_column("SampleID"),
             by = "SampleID") %>% 
  column_to_rownames("SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1),
         Pre_Op_Tstage = factor(Pre_Op_Tstage),
         Pre_Op_Nstage = factor(Pre_Op_Nstage),
         CEA_bin = factor(CEA_bin))


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_strain_meta
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
top_features_overall_strain_meta_before = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Selected features + cTstage + cNstage + CEA]",
    "\n======================================================\n",
    paste(top_features_overall_strain_meta_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_strain_meta_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.8848  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(top_features_overall_strain_before),
                         " Strains + cTstage + cNstage + CEA"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/04-5-10_ROC curve-strain + clinical indicator_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)



##### Strain + Pathway + Clinical Indicator

# Before samples

top_features_overall_path_before_mods = top_features_overall_path_before
top_features_overall_path_before_mods = top_features_overall_path_before_mods %>% 
  str_replace("\\.", "-") %>%
  str_replace("\\.\\.", ": ") %>%
  str_replace("\\.\\.", " (") %>%
  str_replace_all("\\.", " ") %>%
  str_replace("L.", "L-") %>% 
  str_replace("\\ $", ")") %>% 
  recode(.,
         "GL-MANNANAUT-PWY: superpathway of N acetylglucosamine (N acetylmannosamine and N acetylneuraminate degradation" = 
           "GLCMANNANAUT-PWY: superpathway of N-acetylglucosamine, N-acetylmannosamine and N-acetylneuraminate degradation",
         "PWY-5676: acetyl CoA fermentation to butanoate II" = 
           "PWY-5676: acetyl-CoA fermentation to butanoate II")

mb2_selected = mb2 %>% 
  select(SampleID, Pre_Op_Tstage, Pre_Op_Nstage, CEA_bin) %>% 
  column_to_rownames("SampleID")

RF_strain_path_meta = RF_strain_meta %>% 
  rownames_to_column("SampleID") %>% 
  inner_join(RF_path %>% 
               select(-TRG_1) %>% 
               select(any_of(top_features_overall_path_before_mods)) %>% 
               rownames_to_column("SampleID"),
             by = "SampleID") %>% 
  column_to_rownames("SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1),
         Pre_Op_Tstage = factor(Pre_Op_Tstage),
         Pre_Op_Nstage = factor(Pre_Op_Nstage),
         CEA_bin = factor(CEA_bin))


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_strain_path_meta
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
top_features_overall_strain_path_meta_before = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Selected features + cTstage + cNstage + CEA]",
    "\n======================================================\n",
    paste(top_features_overall_strain_path_meta_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_strain_path_meta_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.8909  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(top_features_overall_strain_before), " Strains + ",
                         length(top_features_overall_path_before), " Pathways + cTstage + cNstage + CEA"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/04-5-10_ROC curve-strain + pathway + clinical indicator_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)



##### Species Before + Ongoing

# Before samples

all_species = union(top_features_overall_species_before,
                    top_features_overall_species_ongoing) %>% 
  gsub("\\.", "|", .)

RF_species_before_ongoing = RF_species %>% 
  select(TRG_1, any_of(all_species)) %>% 
  mutate(TRG_1 = factor(TRG_1))


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_species_before_ongoing
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
top_features_overall_species_all = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [All features selected from either 'Before' or 'Ongoing' samples]",
    "\n======================================================\n",
    paste(top_features_overall_species_all, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_species_all, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.8909  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(top_features_overall_species_all), 
                         " All Species Selected from Either 'Before' or 'Ongoing' Samples"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/04-5-11_ROC curve-species_before + ongoing.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)



##### Species + Clinical Indicator

# Before samples

mb2_selected = mb2 %>% 
  select(SampleID, Pre_Op_Tstage, Pre_Op_Nstage, CEA_bin) %>% 
  column_to_rownames("SampleID")

RF_species_meta = RF_species %>% 
  select(TRG_1, any_of(gsub("\\.", "|", top_features_overall_species_before))) %>% 
  rownames_to_column("SampleID") %>% 
  inner_join(mb2_selected %>% 
               rownames_to_column("SampleID"),
             by = "SampleID") %>% 
  column_to_rownames("SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1),
         Pre_Op_Tstage = factor(Pre_Op_Tstage),
         Pre_Op_Nstage = factor(Pre_Op_Nstage),
         CEA_bin = factor(CEA_bin))


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_species_meta
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
top_features_overall_species_meta_before = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Selected features + cTstage + cNstage + CEA]",
    "\n======================================================\n",
    paste(top_features_overall_species_meta_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_species_meta_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.8909  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(top_features_overall_species_before),
                         " Species + cTstage + cNstage + CEA"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/04-5-11_ROC curve-species + clinical indicator_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)



##### Species + Pathway + Clinical Indicator

# Before samples

mb2_selected = mb2 %>% 
  select(SampleID, Pre_Op_Tstage, Pre_Op_Nstage, CEA_bin) %>% 
  column_to_rownames("SampleID")

RF_species_path_meta = RF_species_meta %>% 
  rownames_to_column("SampleID") %>% 
  inner_join(RF_path %>% 
               select(-TRG_1) %>% 
               select(any_of(top_features_overall_path_before_mods)) %>% 
               rownames_to_column("SampleID"),
             by = "SampleID") %>% 
  column_to_rownames("SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1),
         Pre_Op_Tstage = factor(Pre_Op_Tstage),
         Pre_Op_Nstage = factor(Pre_Op_Nstage),
         CEA_bin = factor(CEA_bin))


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_species_path_meta
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
top_features_overall_species_path_meta_before = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Selected features + cTstage + cNstage + CEA]",
    "\n======================================================\n",
    paste(top_features_overall_species_path_meta_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_species_path_meta_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.8909  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(top_features_overall_species_before), " Species + ",
                         length(top_features_overall_path_before), " Pathways + cTstage + cNstage + CEA"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/04-5-11_ROC curve-species + pathway + clinical indicator_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)



#===========================================================================
### Shared Features --------------------------------------------------------

#################### 1. Shared Pathway ####################

shared_path = intersect(top_features_overall_path_before, 
                        top_features_overall_path_ongoing) # 2 pathways

shared_path = recode(shared_path,
                     "PWY.621..sucrose.degradation.III..sucrose.invertase." =
                       "PWY-621: sucrose degradation III (sucrose invertase)",
                     "PANTO.PWY..phosphopantothenate.biosynthesis.I" =
                       "PANTO-PWY: phosphopantothenate biosynthesis I")


# Before samples

RF_shared_path = path %>% 
  filter(Pathway %in% shared_path) %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Pathway, 
              values_from = Abundance_0.5) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_shared_path) # 3 columns (2 pathways)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_shared_path
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
shared_pathway_before = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Shared features]",
    "\n======================================================\n",
    paste(shared_pathway_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(shared_pathway_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.8  

imp_df_before = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_before

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_path), " Shared Pathways"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/05-5-1_ROC curve-shared pathway_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

RF_shared_path_o = path_o %>% 
  filter(Pathway %in% shared_path) %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Pathway, 
              values_from = Abundance_0.5) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_shared_path_o) # 3 columns (2 pathways)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_shared_path_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
shared_pathway_ongoing = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Shared features]",
    "\n======================================================\n",
    paste(shared_pathway_ongoing, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(shared_pathway_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_ongoing = roc(response = true_labels, 
                      predictor = k_fold_pred, 
                      levels = c("nonCR", "CR"), 
                      direction = "<",
                      quiet = TRUE)

auc_val_ongoing = auc(roc_obj_ongoing) ; auc_val_ongoing # 0.5  

imp_df_ongoing = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_ongoing

# 시각화
p5_g_roc_impro = ggroc(roc_obj_ongoing, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_path), " Shared Features"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/05-5-1_ROC curve-shared pathway_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)


# Before, Ongoing Plot
roc_list_all = list(Before = roc_obj_before,
                    Ongoing = roc_obj_ongoing)

p5_g_roc_impro = ggroc(roc_list_all, legacy.axes = TRUE, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Before"  = "#E76F51",
                                "Ongoing" = "#2A9D8F")) +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  annotate("text", x = 0.70, y = 0.15, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#2A9D8F") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_path), " Shared Pathways"),
       x = "1 - Specificity", y = "Sensitivity", color = "TNT") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/05-5-1_ROC curve-shared pathway_all.svg",
#        plot = p5_g_roc_impro, width = 6, height = 5)


### Feature importance, Abundance, Prevalence

# Before
path_shap = colMeans(abs(shap_before_path)) %>% 
  enframe(name = "Feature", value = "Avg_SHAP")
path_importance = imp_df_before_path
path_imp_list = path_importance$Feature

p4_imp_path_mda = path_importance %>% 
  filter(Feature %in% gsub("[^A-Za-z0-9]", ".", shared_path)) %>% 
  select(-Avg_MeanDecreaseGini) %>% 
  mutate(Feature = factor(Feature, levels = rev(path_imp_list))) %>% 
  ggplot(aes(Feature, Avg_MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "#4682B4", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 7.5)) +
  labs(title = "Mean Decrease Accuracy") +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) ; p4_imp_path_mda

p4_imp_path_mdg = path_importance %>% 
  filter(Feature %in% gsub("[^A-Za-z0-9]", ".", shared_path)) %>% 
  select(-Avg_MeanDecreaseAccuracy) %>% 
  mutate(Feature = factor(Feature, levels = rev(path_imp_list))) %>% 
  ggplot(aes(Feature, Avg_MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "#6A994E", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1.5)) +
  labs(title = "Mean Decrease Gini") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) ; p4_imp_path_mdg

p4_imp_path_shap = path_shap %>% 
  filter(Feature %in% gsub("[^A-Za-z0-9]", ".", shared_path)) %>%
  mutate(Feature = factor(Feature, levels = rev(path_imp_list))) %>% 
  ggplot(aes(Feature, Avg_SHAP)) +
  geom_bar(stat = "identity", fill = "#B22222", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 0.05)) +
  labs(title = "SHAP") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) ; p4_imp_path_shap

p4_imp_path_abund = path %>% 
  filter(Pathway %in% shared_path) %>% 
  group_by(Pathway, TRG_1) %>% 
  summarise(Mean_abundance = mean(Abundance_0.5, na.rm = T),
            `Log10(Mean_abundance)` = log10(Mean_abundance + 1e-05),
            .groups = "drop") %>% 
  pivot_wider(names_from = TRG_1,
              values_from = c(Mean_abundance, `Log10(Mean_abundance)`),
              names_glue = "{.value}_{TRG_1}") %>% 
  pivot_longer(cols = c(`Log10(Mean_abundance)_CR`, `Log10(Mean_abundance)_nonCR`),
               names_to = "Group", values_to = "Log10(Mean_abundance)") %>%
  mutate(Group = factor(Group, levels = c("Log10(Mean_abundance)_nonCR", "Log10(Mean_abundance)_CR")),
         Pathway = factor(Pathway, levels = shared_path),
         `Log10(Mean_abundance)` = `Log10(Mean_abundance)`+ 5.1) %>% 
  ggplot(aes(x = Pathway, y = `Log10(Mean_abundance)`, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  scale_fill_manual(values = c("Log10(Mean_abundance)_nonCR" = "#80461B", "Log10(Mean_abundance)_CR" = "#F7D9BC"),
                    labels = c("nonCR", "CR")) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 9.5)) +
  labs(title = "Log10(Abund. + P)") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") ; p4_imp_path_abund

p4_imp_path_prev = path %>% 
  filter(Pathway %in% shared_path) %>% 
  group_by(Pathway, TRG_1) %>% 
  summarise(Prevalence = sum(Abundance_0.5 > 0),
            .groups = "drop") %>% 
  pivot_wider(names_from = TRG_1,
              values_from = Prevalence,
              names_glue = "{.value}_{TRG_1}") %>% 
  pivot_longer(cols = c(Prevalence_CR, Prevalence_nonCR),
               names_to = "Group", values_to = "Prevalence") %>%
  mutate(Group = factor(Group, levels = c("Prevalence_nonCR", "Prevalence_CR")),
         Pathway = factor(Pathway, levels = shared_path)) %>% 
  ggplot(aes(x = Pathway, y = Prevalence, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  scale_fill_manual(values = c("Prevalence_nonCR" = "#80461B", "Prevalence_CR" = "#F7D9BC"),
                    labels = c("nonCR", "CR")) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 15)) +
  labs(title = "Prevalence") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") ; p4_imp_path_prev


# Ongoing
path_shap_o = colMeans(abs(shap_ongoing_path)) %>% 
  enframe(name = "Feature", value = "Avg_SHAP")
path_importance_o = imp_df_ongoing_path

p5_imp_path_mda = path_importance_o %>% 
  filter(Feature %in% gsub("[^A-Za-z0-9]", ".", shared_path)) %>% 
  select(-Avg_MeanDecreaseGini) %>% 
  mutate(Feature = factor(Feature, levels = rev(path_imp_list))) %>% 
  ggplot(aes(Feature, Avg_MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "#4682B4", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 7.5)) +
  labs(title = "Mean Decrease Accuracy") +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) ; p5_imp_path_mda

p5_imp_path_mdg = path_importance_o %>% 
  filter(Feature %in% gsub("[^A-Za-z0-9]", ".", shared_path)) %>% 
  select(-Avg_MeanDecreaseAccuracy) %>% 
  mutate(Feature = factor(Feature, levels = rev(path_imp_list))) %>% 
  ggplot(aes(Feature, Avg_MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "#6A994E", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1.5)) +
  labs(title = "Mean Decrease Gini") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) ; p5_imp_path_mdg

p5_imp_path_shap = path_shap_o %>% 
  filter(Feature %in% gsub("[^A-Za-z0-9]", ".", shared_path)) %>%
  mutate(Feature = factor(Feature, levels = rev(path_imp_list))) %>% 
  ggplot(aes(Feature, Avg_SHAP)) +
  geom_bar(stat = "identity", fill = "#B22222", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 0.05)) +
  labs(title = "SHAP") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) ; p5_imp_path_shap

p5_imp_path_abund = path_o %>% 
  filter(Pathway %in% shared_path) %>% 
  group_by(Pathway, TRG_1) %>% 
  summarise(Mean_abundance = mean(Abundance_0.5, na.rm = T),
            `Log10(Mean_abundance)` = log10(Mean_abundance + 1736.767),
            .groups = "drop") %>% 
  pivot_wider(names_from = TRG_1,
              values_from = c(Mean_abundance, `Log10(Mean_abundance)`),
              names_glue = "{.value}_{TRG_1}") %>% 
  pivot_longer(cols = c(`Log10(Mean_abundance)_CR`, `Log10(Mean_abundance)_nonCR`),
               names_to = "Group", values_to = "Log10(Mean_abundance)") %>%
  mutate(Group = factor(Group, levels = c("Log10(Mean_abundance)_nonCR", "Log10(Mean_abundance)_CR")),
         Pathway = factor(Pathway, levels = shared_path),
         `Log10(Mean_abundance)` = `Log10(Mean_abundance)` + 5.1) %>% 
  ggplot(aes(x = Pathway, y = `Log10(Mean_abundance)`, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  scale_fill_manual(values = c("Log10(Mean_abundance)_nonCR" = "#80461B", "Log10(Mean_abundance)_CR" = "#F7D9BC"),
                    labels = c("nonCR", "CR")) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 9.5)) +
  labs(title = "Log10(Abund. + P)") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") ; p5_imp_path_abund

p5_imp_path_prev = path_o %>% 
  filter(Pathway %in% shared_path) %>% 
  group_by(Pathway, TRG_1) %>% 
  summarise(Prevalence = sum(Abundance_0.5 > 0),
            .groups = "drop") %>% 
  pivot_wider(names_from = TRG_1,
              values_from = Prevalence,
              names_glue = "{.value}_{TRG_1}") %>% 
  pivot_longer(cols = c(Prevalence_CR, Prevalence_nonCR),
               names_to = "Group", values_to = "Prevalence") %>%
  mutate(Group = factor(Group, levels = c("Prevalence_nonCR", "Prevalence_CR")),
         Pathway = factor(Pathway, levels = shared_path)) %>% 
  ggplot(aes(x = Pathway, y = Prevalence, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  scale_fill_manual(values = c("Prevalence_nonCR" = "#80461B", "Prevalence_CR" = "#F7D9BC"),
                    labels = c("nonCR", "CR")) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 15)) +
  labs(title = "Prevalence") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") ; p5_imp_path_prev

p5_imp_path = grid.arrange(p4_imp_path_mda,
                           p4_imp_path_mdg, 
                           p4_imp_path_shap,
                           p4_imp_path_abund,
                           p4_imp_path_prev,
                           p5_imp_path_mda,
                           p5_imp_path_mdg, 
                           p5_imp_path_shap,
                           p5_imp_path_abund,
                           p5_imp_path_prev,
                           nrow = 2, widths = c(2, 1, 1, 1, 1))

# ggsave("figure/05-5-1_Importance barplot-shared pathway_all.svg",
#        plot = p5_imp_path, width = 20, height = 10)



#################### 2. Shared Strain ####################

shared_strain = intersect(top_features_overall_strain_before, 
                          top_features_overall_strain_ongoing) # 12 strains

shared_strain = gsub("\\.", "|", shared_strain)


# Before samples

RF_shared_strain = tb %>% 
  filter(Strain %in% shared_strain) %>% 
  pivot_longer(cols = -Strain,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  merge(mb %>% 
          select(SampleID,
                 TRG_score, TRG_1, TRG_2, TRG_3,
                 Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                 Pre_Op_Nstage, Pre_Op_Nstage_bin,
                 BMI, BMI_bin, CEA, CEA_bin),
        by = "SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR"))) %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Strain,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_shared_strain) # 13 columns (12 strains)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_shared_strain
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
shared_strain_before = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Shared features]",
    "\n======================================================\n",
    paste(shared_strain_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(shared_strain_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.8909  

imp_df_before = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_before

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_strain), " Shared Strains"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/05-5-2_ROC curve-shared strain_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

RF_shared_strain_o = to %>% 
  filter(Strain %in% shared_strain) %>% 
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
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR"))) %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Strain,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_shared_strain_o) # 13 columns (12 strains)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_shared_strain_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
shared_strain_ongoing = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Shared features]",
    "\n======================================================\n",
    paste(shared_strain_ongoing, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(shared_strain_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_ongoing = roc(response = true_labels, 
                      predictor = k_fold_pred, 
                      levels = c("nonCR", "CR"), 
                      direction = "<",
                      quiet = TRUE)

auc_val_ongoing = auc(roc_obj_ongoing) ; auc_val_ongoing # 0.9524  

imp_df_ongoing = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_ongoing

# 시각화
p5_g_roc_impro = ggroc(roc_obj_ongoing, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_strain), " Shared Strains"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/05-5-2_ROC curve-shared strain_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)


# Before, Ongoing Plot
roc_list_all = list(Before = roc_obj_before,
                    Ongoing = roc_obj_ongoing)

p5_g_roc_impro = ggroc(roc_list_all, legacy.axes = TRUE, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Before"  = "#E76F51",
                                "Ongoing" = "#2A9D8F")) +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  annotate("text", x = 0.70, y = 0.15, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#2A9D8F") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_strain), " Shared Strains"),
       x = "1 - Specificity", y = "Sensitivity", color = "TNT") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/05-5-2_ROC curve-shared strain_all.svg",
#        plot = p5_g_roc_impro, width = 6, height = 5)


### Feature importance, Abundance, Prevalence

# Before
strain_shap = colMeans(abs(shap_before_strain)) %>% 
  enframe(name = "Feature", value = "Avg_SHAP")
strain_importance = imp_df_before_strain
strain_imp_list = strain_importance$Feature %>% 
  gsub("\\.", "|", .)

p4_imp_strain_mda = strain_importance %>% 
  mutate(Feature = gsub("\\.", "|", Feature)) %>% 
  filter(Feature %in% shared_strain) %>% 
  select(-Avg_MeanDecreaseGini) %>% 
  mutate(Feature = factor(Feature, levels = rev(strain_imp_list))) %>% 
  ggplot(aes(Feature, Avg_MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "#4682B4", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(-0.35, 7.0)) +
  labs(title = "Mean Decrease Accuracy") +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) ; p4_imp_strain_mda

p4_imp_strain_mdg = strain_importance %>% 
  mutate(Feature = gsub("\\.", "|", Feature)) %>% 
  filter(Feature %in% shared_strain) %>% 
  select(-Avg_MeanDecreaseAccuracy) %>% 
  mutate(Feature = factor(Feature, levels = rev(strain_imp_list))) %>% 
  ggplot(aes(Feature, Avg_MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "#6A994E", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 0.55)) +
  labs(title = "Mean Decrease Gini") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) ; p4_imp_strain_mdg

p4_imp_strain_shap = strain_shap %>% 
  mutate(Feature = gsub("\\.", "|", Feature)) %>% 
  filter(Feature %in% shared_strain) %>% 
  mutate(Feature = factor(Feature, levels = rev(strain_imp_list))) %>% 
  ggplot(aes(Feature, Avg_SHAP)) +
  geom_bar(stat = "identity", fill = "#B22222", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 0.03)) +
  labs(title = "SHAP") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) ; p4_imp_strain_shap

p4_imp_strain_abund = tb_input_0.1 %>% 
  filter(Strain %in% shared_strain) %>% 
  group_by(Strain, TRG_1) %>% 
  summarise(Mean_abundance = mean(abundance, na.rm = T),
            `Log10(Mean_abundance)` = log10(Mean_abundance + 1e-05),
            .groups = "drop") %>% 
  pivot_wider(names_from = TRG_1,
              values_from = c(Mean_abundance, `Log10(Mean_abundance)`),
              names_glue = "{.value}_{TRG_1}") %>% 
  pivot_longer(cols = c(`Log10(Mean_abundance)_CR`, `Log10(Mean_abundance)_nonCR`),
               names_to = "Group", values_to = "Log10(Mean_abundance)") %>%
  mutate(Group = factor(Group, levels = c("Log10(Mean_abundance)_nonCR", "Log10(Mean_abundance)_CR")),
         Strain = factor(Strain, levels = rev(strain_imp_list)),
         `Log10(Mean_abundance)` = `Log10(Mean_abundance)` + 5.1) %>% 
  ggplot(aes(x = Strain, y = `Log10(Mean_abundance)`, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  scale_fill_manual(values = c("Log10(Mean_abundance)_nonCR" = "#80461B", "Log10(Mean_abundance)_CR" = "#F7D9BC"),
                    labels = c("nonCR", "CR")) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 5.5)) +
  labs(title = "Log10(Abund. + P)") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") ; p4_imp_strain_abund

p4_imp_strain_prev = tb_input_0.1 %>% 
  filter(Strain %in% shared_strain) %>% 
  group_by(Strain, TRG_1) %>% 
  summarise(Prevalence = sum(abundance > 0),
            .groups = "drop") %>% 
  pivot_wider(names_from = TRG_1,
              values_from = Prevalence,
              names_glue = "{.value}_{TRG_1}") %>% 
  pivot_longer(cols = c(Prevalence_CR, Prevalence_nonCR),
               names_to = "Group", values_to = "Prevalence") %>%
  mutate(Group = factor(Group, levels = c("Prevalence_nonCR", "Prevalence_CR")),
         Strain = factor(Strain, levels = rev(strain_imp_list))) %>% 
  ggplot(aes(x = Strain, y = Prevalence, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  scale_fill_manual(values = c("Prevalence_nonCR" = "#80461B", "Prevalence_CR" = "#F7D9BC"),
                    labels = c("nonCR", "CR")) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 14)) +
  labs(title = "Prevalence") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") ; p4_imp_strain_prev


# Ongoing
strain_shap_o = colMeans(abs(shap_ongoing_strain)) %>% 
  enframe(name = "Feature", value = "Avg_SHAP")
strain_importance_o = imp_df_ongoing_strain

p5_imp_strain_mda = strain_importance_o %>% 
  mutate(Feature = gsub("\\.", "|", Feature)) %>% 
  filter(Feature %in% shared_strain) %>% 
  select(-Avg_MeanDecreaseGini) %>% 
  mutate(Feature = factor(Feature, levels = rev(strain_imp_list))) %>% 
  ggplot(aes(Feature, Avg_MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "#4682B4", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(-0.35, 7.0)) +
  labs(title = "Mean Decrease Accuracy") +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) ; p5_imp_strain_mda

p5_imp_strain_mdg = strain_importance_o %>% 
  mutate(Feature = gsub("\\.", "|", Feature)) %>% 
  filter(Feature %in% shared_strain) %>% 
  select(-Avg_MeanDecreaseAccuracy) %>% 
  mutate(Feature = factor(Feature, levels = rev(strain_imp_list))) %>% 
  ggplot(aes(Feature, Avg_MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "#6A994E", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 0.55)) +
  labs(title = "Mean Decrease Gini") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) ; p5_imp_strain_mdg

p5_imp_strain_shap = strain_shap_o %>% 
  mutate(Feature = gsub("\\.", "|", Feature)) %>% 
  filter(Feature %in% shared_strain) %>% 
  mutate(Feature = factor(Feature, levels = rev(strain_imp_list))) %>% 
  ggplot(aes(Feature, Avg_SHAP)) +
  geom_bar(stat = "identity", fill = "#B22222", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 0.03)) +
  labs(title = "SHAP") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) ; p5_imp_strain_shap

p5_imp_strain_abund = to_input_0.1 %>% 
  filter(Strain %in% shared_strain) %>% 
  group_by(Strain, TRG_1) %>% 
  summarise(Mean_abundance = mean(abundance, na.rm = T),
            `Log10(Mean_abundance)` = log10(Mean_abundance + 1e-05),
            .groups = "drop") %>% 
  pivot_wider(names_from = TRG_1,
              values_from = c(Mean_abundance, `Log10(Mean_abundance)`),
              names_glue = "{.value}_{TRG_1}") %>% 
  pivot_longer(cols = c(`Log10(Mean_abundance)_CR`, `Log10(Mean_abundance)_nonCR`),
               names_to = "Group", values_to = "Log10(Mean_abundance)") %>%
  mutate(Group = factor(Group, levels = c("Log10(Mean_abundance)_nonCR", "Log10(Mean_abundance)_CR")),
         Strain = factor(Strain, levels = rev(strain_imp_list)),
         `Log10(Mean_abundance)` = `Log10(Mean_abundance)` + 5.1) %>% 
  ggplot(aes(x = Strain, y = `Log10(Mean_abundance)`, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  scale_fill_manual(values = c("Log10(Mean_abundance)_nonCR" = "#80461B", "Log10(Mean_abundance)_CR" = "#F7D9BC"),
                    labels = c("nonCR", "CR")) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 5.5)) +
  labs(title = "Log10(Abund. + P)") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") ; p5_imp_strain_abund

p5_imp_strain_prev = to_input_0.1 %>% 
  filter(Strain %in% shared_strain) %>% 
  group_by(Strain, TRG_1) %>% 
  summarise(Prevalence = sum(abundance > 0),
            .groups = "drop") %>% 
  pivot_wider(names_from = TRG_1,
              values_from = Prevalence,
              names_glue = "{.value}_{TRG_1}") %>% 
  pivot_longer(cols = c(Prevalence_CR, Prevalence_nonCR),
               names_to = "Group", values_to = "Prevalence") %>%
  mutate(Group = factor(Group, levels = c("Prevalence_nonCR", "Prevalence_CR")),
         Strain = factor(Strain, levels = rev(strain_imp_list))) %>% 
  ggplot(aes(x = Strain, y = Prevalence, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  scale_fill_manual(values = c("Prevalence_nonCR" = "#80461B", "Prevalence_CR" = "#F7D9BC"),
                    labels = c("nonCR", "CR")) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 14)) +
  labs(title = "Prevalence") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") ; p5_imp_strain_prev

p5_imp_strain = grid.arrange(p4_imp_strain_mda,
                             p4_imp_strain_mdg, 
                             p4_imp_strain_shap,
                             p4_imp_strain_abund,
                             p4_imp_strain_prev,
                             p5_imp_strain_mda,
                             p5_imp_strain_mdg, 
                             p5_imp_strain_shap,
                             p5_imp_strain_abund,
                             p5_imp_strain_prev,
                             nrow = 2, widths = c(2, 1, 1, 1, 1))

# ggsave("figure/05-5-2_Importance barplot-shared strain_all.svg",
#        plot = p5_imp_strain, width = 20, height = 10)



#################### 3. Shared Species ####################

shared_species = intersect(top_features_overall_species_before, 
                           top_features_overall_species_ongoing) # 10 species


# Before samples

RF_shared_species = sb %>% 
  filter(Species %in% shared_species) %>% 
  pivot_longer(cols = -Species,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  merge(mb %>% 
          select(SampleID,
                 TRG_score, TRG_1, TRG_2, TRG_3,
                 Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                 Pre_Op_Nstage, Pre_Op_Nstage_bin,
                 BMI, BMI_bin, CEA, CEA_bin),
        by = "SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR"))) %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Species,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_shared_species) # 11 columns (10 species)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_shared_species
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
shared_species_before = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Shared features]",
    "\n======================================================\n",
    paste(shared_species_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(shared_species_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.9333  

imp_df_before = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_before

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_species), " Shared Species"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/05-5-3_ROC curve-shared species_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

RF_shared_species_o = so %>% 
  filter(Species %in% shared_species) %>% 
  pivot_longer(cols = -Species,
               names_to = "SampleID",
               values_to = "abundance") %>% 
  merge(mo %>% 
          select(SampleID,
                 TRG_score, TRG_1, TRG_2, TRG_3,
                 Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                 Pre_Op_Nstage, Pre_Op_Nstage_bin,
                 BMI, BMI_bin, CEA, CEA_bin),
        by = "SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1, levels = c("CR", "nonCR"))) %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Species,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_shared_species_o) # 11 columns (10 species)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_shared_species_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
shared_species_ongoing = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Shared features]",
    "\n======================================================\n",
    paste(shared_species_ongoing, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(shared_species_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_ongoing = roc(response = true_labels, 
                      predictor = k_fold_pred, 
                      levels = c("nonCR", "CR"), 
                      direction = "<",
                      quiet = TRUE)

auc_val_ongoing = auc(roc_obj_ongoing) ; auc_val_ongoing # 0.9683  

imp_df_ongoing = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_ongoing

# 시각화
p5_g_roc_impro = ggroc(roc_obj_ongoing, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_species), " Shared Species"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/05-5-3_ROC curve-shared species_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)


# Before, Ongoing Plot
roc_list_all = list(Before = roc_obj_before,
                    Ongoing = roc_obj_ongoing)

p5_g_roc_impro = ggroc(roc_list_all, legacy.axes = TRUE, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Before"  = "#E76F51",
                                "Ongoing" = "#2A9D8F")) +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  annotate("text", x = 0.70, y = 0.15, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#2A9D8F") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_species), " Shared Species"),
       x = "1 - Specificity", y = "Sensitivity", color = "TNT") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/05-5-3_ROC curve-shared species_all.svg",
#        plot = p5_g_roc_impro, width = 6, height = 5)


### Feature importance, Abundance, Prevalence

# Before
species_shap = colMeans(abs(shap_before_species)) %>% 
  enframe(name = "Feature", value = "Avg_SHAP")
species_importance = imp_df_before_species
species_imp_list = species_importance$Feature %>% 
  gsub("\\.", "|", .)

p4_imp_species_mda = species_importance %>% 
  mutate(Feature = gsub("\\.", "|", Feature)) %>% 
  filter(Feature %in% shared_species) %>% 
  select(-Avg_MeanDecreaseGini) %>% 
  mutate(Feature = factor(Feature, levels = rev(species_imp_list))) %>% 
  ggplot(aes(Feature, Avg_MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "#4682B4", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 7.5)) +
  labs(title = "Mean Decrease Accuracy") +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) ; p4_imp_species_mda

p4_imp_species_mdg = species_importance %>% 
  mutate(Feature = gsub("\\.", "|", Feature)) %>% 
  filter(Feature %in% shared_species) %>% 
  select(-Avg_MeanDecreaseAccuracy) %>% 
  mutate(Feature = factor(Feature, levels = rev(species_imp_list))) %>% 
  ggplot(aes(Feature, Avg_MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "#6A994E", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1.5)) +
  labs(title = "Mean Decrease Gini") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) ; p4_imp_species_mdg

p4_imp_species_shap = species_shap %>% 
  mutate(Feature = gsub("\\.", "|", Feature)) %>% 
  filter(Feature %in% shared_species) %>% 
  mutate(Feature = factor(Feature, levels = rev(species_imp_list))) %>% 
  ggplot(aes(Feature, Avg_SHAP)) +
  geom_bar(stat = "identity", fill = "#B22222", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 0.05)) +
  labs(title = "SHAP") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) ; p4_imp_species_shap

p4_imp_species_abund = sb_input_0.1 %>% 
  filter(Species %in% shared_species) %>% 
  group_by(Species, TRG_1) %>% 
  summarise(Mean_abundance = mean(abundance, na.rm = T),
            `Log10(Mean_abundance)` = log10(Mean_abundance + 1e-05),
            .groups = "drop") %>% 
  pivot_wider(names_from = TRG_1,
              values_from = c(Mean_abundance, `Log10(Mean_abundance)`),
              names_glue = "{.value}_{TRG_1}") %>% 
  pivot_longer(cols = c(`Log10(Mean_abundance)_CR`, `Log10(Mean_abundance)_nonCR`),
               names_to = "Group", values_to = "Log10(Mean_abundance)") %>%
  mutate(Group = factor(Group, levels = c("Log10(Mean_abundance)_nonCR", "Log10(Mean_abundance)_CR")),
         Species = factor(Species, levels = rev(species_imp_list)),
         `Log10(Mean_abundance)` = `Log10(Mean_abundance)` + 5.1) %>% 
  ggplot(aes(x = Species, y = `Log10(Mean_abundance)`, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  scale_fill_manual(values = c("Log10(Mean_abundance)_nonCR" = "#80461B", "Log10(Mean_abundance)_CR" = "#F7D9BC"),
                    labels = c("nonCR", "CR")) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 9.5)) +
  labs(title = "Log10(Abund. + P)") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") ; p4_imp_species_abund

p4_imp_species_prev = sb_input_0.1 %>% 
  filter(Species %in% shared_species) %>% 
  group_by(Species, TRG_1) %>% 
  summarise(Prevalence = sum(abundance > 0),
            .groups = "drop") %>% 
  pivot_wider(names_from = TRG_1,
              values_from = Prevalence,
              names_glue = "{.value}_{TRG_1}") %>% 
  pivot_longer(cols = c(Prevalence_CR, Prevalence_nonCR),
               names_to = "Group", values_to = "Prevalence") %>%
  mutate(Group = factor(Group, levels = c("Prevalence_nonCR", "Prevalence_CR")),
         Species = factor(Species, levels = rev(species_imp_list))) %>% 
  ggplot(aes(x = Species, y = Prevalence, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  scale_fill_manual(values = c("Prevalence_nonCR" = "#80461B", "Prevalence_CR" = "#F7D9BC"),
                    labels = c("nonCR", "CR")) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 15)) +
  labs(title = "Prevalence") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") ; p4_imp_species_prev


# Ongoing
species_shap_o = colMeans(abs(shap_ongoing_species)) %>% 
  enframe(name = "Feature", value = "Avg_SHAP")
species_importance_o = imp_df_ongoing_species

p5_imp_species_mda = species_importance_o %>% 
  mutate(Feature = gsub("\\.", "|", Feature)) %>% 
  filter(Feature %in% shared_species) %>% 
  select(-Avg_MeanDecreaseGini) %>% 
  mutate(Feature = factor(Feature, levels = rev(species_imp_list))) %>% 
  ggplot(aes(Feature, Avg_MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "#4682B4", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 7.5)) +
  labs(title = "Mean Decrease Accuracy") +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) ; p5_imp_species_mda

p5_imp_species_mdg = species_importance_o %>% 
  mutate(Feature = gsub("\\.", "|", Feature)) %>% 
  filter(Feature %in% shared_species) %>% 
  select(-Avg_MeanDecreaseAccuracy) %>% 
  mutate(Feature = factor(Feature, levels = rev(species_imp_list))) %>% 
  ggplot(aes(Feature, Avg_MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "#6A994E", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1.5)) +
  labs(title = "Mean Decrease Gini") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) ; p5_imp_species_mdg

p5_imp_species_shap = species_shap_o %>% 
  mutate(Feature = gsub("\\.", "|", Feature)) %>% 
  filter(Feature %in% shared_species) %>% 
  mutate(Feature = factor(Feature, levels = rev(species_imp_list))) %>% 
  ggplot(aes(Feature, Avg_SHAP)) +
  geom_bar(stat = "identity", fill = "#B22222", color = "black") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 0.05)) +
  labs(title = "SHAP") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) ; p5_imp_species_shap

p5_imp_species_abund = so_input_0.1 %>% 
  filter(Species %in% shared_species) %>% 
  group_by(Species, TRG_1) %>% 
  summarise(Mean_abundance = mean(abundance, na.rm = T),
            `Log10(Mean_abundance)` = log10(Mean_abundance + 1e-05),
            .groups = "drop") %>% 
  pivot_wider(names_from = TRG_1,
              values_from = c(Mean_abundance, `Log10(Mean_abundance)`),
              names_glue = "{.value}_{TRG_1}") %>% 
  pivot_longer(cols = c(`Log10(Mean_abundance)_CR`, `Log10(Mean_abundance)_nonCR`),
               names_to = "Group", values_to = "Log10(Mean_abundance)") %>%
  mutate(Group = factor(Group, levels = c("Log10(Mean_abundance)_nonCR", "Log10(Mean_abundance)_CR")),
         Species = factor(Species, levels = rev(species_imp_list)),
         `Log10(Mean_abundance)` = `Log10(Mean_abundance)` + 5.1) %>% 
  ggplot(aes(x = Species, y = `Log10(Mean_abundance)`, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  scale_fill_manual(values = c("Log10(Mean_abundance)_nonCR" = "#80461B", "Log10(Mean_abundance)_CR" = "#F7D9BC"),
                    labels = c("nonCR", "CR")) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 9.5)) +
  labs(title = "Log10(Abund. + P)") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") ; p5_imp_species_abund

p5_imp_species_prev = so_input_0.1 %>% 
  filter(Species %in% shared_species) %>% 
  group_by(Species, TRG_1) %>% 
  summarise(Prevalence = sum(abundance > 0),
            .groups = "drop") %>% 
  pivot_wider(names_from = TRG_1,
              values_from = Prevalence,
              names_glue = "{.value}_{TRG_1}") %>% 
  pivot_longer(cols = c(Prevalence_CR, Prevalence_nonCR),
               names_to = "Group", values_to = "Prevalence") %>%
  mutate(Group = factor(Group, levels = c("Prevalence_nonCR", "Prevalence_CR")),
         Species = factor(Species, levels = rev(species_imp_list))) %>% 
  ggplot(aes(x = Species, y = Prevalence, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  scale_fill_manual(values = c("Prevalence_nonCR" = "#80461B", "Prevalence_CR" = "#F7D9BC"),
                    labels = c("nonCR", "CR")) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 15)) +
  labs(title = "Prevalence") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") ; p5_imp_species_prev

p5_imp_species = grid.arrange(p4_imp_species_mda,
                              p4_imp_species_mdg, 
                              p4_imp_species_shap,
                              p4_imp_species_abund,
                              p4_imp_species_prev,
                              p5_imp_species_mda,
                              p5_imp_species_mdg, 
                              p5_imp_species_shap,
                              p5_imp_species_abund,
                              p5_imp_species_prev,
                              nrow = 2, widths = c(2, 1, 1, 1, 1))

# ggsave("figure/05-5-3_Importance barplot-shared species_all.svg",
#        plot = p5_imp_species, width = 20, height = 10)



#################### 4. Clinical Indicator ####################

# Before samples

RF_meta = mb2 %>% 
  select(SampleID, TRG_1, Pre_Op_Tstage, Pre_Op_Nstage, BMI, CEA) %>% 
  mutate(TRG_1 = factor(TRG_1),
         Pre_Op_Tstage = factor(Pre_Op_Tstage),
         Pre_Op_Nstage = factor(Pre_Op_Nstage),
         CEA_bin = ifelse(CEA < 5.0, "<5", "5+")) %>% 
  select(-c(BMI, CEA)) %>% 
  column_to_rownames("SampleID")


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_meta
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
top_features_overall_meta_before = colnames(RF_meta)[colnames(RF_meta) != "TRG_1"]

cat("\n\n============================================",
    "\n [Clinical indicators]",
    "\n======================================================\n",
    paste(top_features_overall_meta_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_meta_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx]  = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.7121  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/05-5-4_ROC curve-clinical indicator_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

RF_meta_o = mo %>% 
  select(SampleID, TRG_1, Pre_Op_Tstage, Pre_Op_Nstage, BMI, CEA) %>% 
  mutate(TRG_1 = factor(TRG_1),
         Pre_Op_Tstage = factor(Pre_Op_Tstage),
         Pre_Op_Nstage = factor(Pre_Op_Nstage),
         CEA_bin = ifelse(CEA < 5.0, "<5", "5+")) %>% 
  select(-c(BMI, CEA)) %>% 
  column_to_rownames("SampleID")


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_meta_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
top_features_overall_meta_ongoing = colnames(RF_meta)[colnames(RF_meta) != "TRG_1"]

cat("\n\n============================================",
    "\n [Clinical indicators]",
    "\n======================================================\n",
    paste(top_features_overall_meta_ongoing, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(top_features_overall_meta_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_ongoing = roc(response = true_labels, 
                      predictor = k_fold_pred, 
                      levels = c("nonCR", "CR"), 
                      direction = "<",
                      quiet = TRUE)

auc_val_ongoing = auc(roc_obj_ongoing) ; auc_val_ongoing # 0.881  

imp_df = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df

# 시각화
p5_g_roc_impro = ggroc(roc_obj_ongoing, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/05-5-4_ROC curve-clinical indicator_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)


# Before, Ongoing Plot
roc_list_all = list(Before = roc_obj_before,
                    Ongoing = roc_obj_ongoing)

p5_g_roc_impro = ggroc(roc_list_all, legacy.axes = TRUE, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Before"  = "#E76F51",
                                "Ongoing" = "#2A9D8F")) +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  annotate("text", x = 0.70, y = 0.15, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#2A9D8F") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = "cTstage + cNstage + CEA",
       x = "1 - Specificity", y = "Sensitivity", color = "TNT") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/05-5-4_ROC curve-clinical indicator_all.svg",
#        plot = p5_g_roc_impro, width = 6, height = 5)



#################### 5. Combined Feature ####################

##### Shared Strain + Clinical Indicator

# Before samples

mb2_selected = mb2 %>% 
  select(SampleID, Pre_Op_Tstage, Pre_Op_Nstage, CEA_bin) %>% 
  column_to_rownames("SampleID")

RF_shared_strain_meta = RF_shared_strain %>% 
  rownames_to_column("SampleID") %>% 
  inner_join(mb2_selected %>% 
               rownames_to_column("SampleID"),
             by = "SampleID") %>% 
  column_to_rownames("SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1),
         Pre_Op_Tstage = factor(Pre_Op_Tstage),
         Pre_Op_Nstage = factor(Pre_Op_Nstage),
         CEA_bin = factor(CEA_bin))


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_shared_strain_meta
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
shared_strain_meta_before = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Shared features + cTstage + cNstage + CEA]",
    "\n======================================================\n",
    paste(shared_strain_meta_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(shared_strain_meta_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.8848  

imp_df_before = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_before

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_strain), " Shared Strains + cTstage + cNstage + CEA"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/05-5-4_ROC curve-shared strain + clinical indicator_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

mo2_selected = mo %>% 
  select(SampleID, Pre_Op_Tstage, Pre_Op_Nstage, CEA_bin) %>% 
  column_to_rownames("SampleID")

RF_shared_strain_meta_o = RF_shared_strain_o %>% 
  rownames_to_column("SampleID") %>% 
  inner_join(mo2_selected %>% 
               rownames_to_column("SampleID"),
             by = "SampleID") %>% 
  column_to_rownames("SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1),
         Pre_Op_Tstage = factor(Pre_Op_Tstage),
         Pre_Op_Nstage = factor(Pre_Op_Nstage),
         CEA_bin = factor(CEA_bin))


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_shared_strain_meta_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
shared_strain_meta_ongoing = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Shared features + cTstage + cNstage + CEA]",
    "\n======================================================\n",
    paste(shared_strain_meta_ongoing, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(shared_strain_meta_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_ongoing = roc(response = true_labels, 
                      predictor = k_fold_pred, 
                      levels = c("nonCR", "CR"), 
                      direction = "<",
                      quiet = TRUE)

auc_val_ongoing = auc(roc_obj_ongoing) ; auc_val_ongoing # 0.9841  

imp_df_ongoing = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_ongoing

# 시각화
p5_g_roc_impro = ggroc(roc_obj_ongoing, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_strain), " Shared Strain + cTstage + cNstage + CEA"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/05-5-4_ROC curve-shared strain + clinical indicator_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)


# Before, Ongoing Plot
roc_list_all = list(Before = roc_obj_before,
                    Ongoing = roc_obj_ongoing)

p5_g_roc_impro = ggroc(roc_list_all, legacy.axes = TRUE, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Before"  = "#E76F51",
                                "Ongoing" = "#2A9D8F")) +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  annotate("text", x = 0.70, y = 0.15, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#2A9D8F") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_strain), " Shared Strains + cTstage + cNstage + CEA"),
       x = "1 - Specificity", y = "Sensitivity", color = "TNT") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/05-5-4_ROC curve-shared strain + clinical indicator_all.svg",
#        plot = p5_g_roc_impro, width = 6, height = 5)



##### Shared Strain + Shared Pathway + Clinical Indicator

# Before samples

mb2_selected = mb2 %>% 
  select(SampleID, Pre_Op_Tstage, Pre_Op_Nstage, CEA_bin) %>% 
  column_to_rownames("SampleID")

RF_shared_strain_path_meta = RF_shared_strain %>% 
  rownames_to_column("SampleID") %>% 
  inner_join(RF_shared_path %>% 
               select(-TRG_1) %>% 
               rownames_to_column("SampleID"),
             by = "SampleID") %>% 
  inner_join(mb2_selected %>% 
               rownames_to_column("SampleID"),
             by = "SampleID") %>% 
  column_to_rownames("SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1),
         Pre_Op_Tstage = factor(Pre_Op_Tstage),
         Pre_Op_Nstage = factor(Pre_Op_Nstage),
         CEA_bin = factor(CEA_bin))


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_shared_strain_path_meta
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
shared_strain_path_meta_before = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Shared features + cTstage + cNstage + CEA]",
    "\n======================================================\n",
    paste(shared_strain_path_meta_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(123)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(shared_strain_path_meta_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.9091  

imp_df_before = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_before

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_strain), " Shared Strains + ",
                         length(shared_path), " Shared Pathways + cTstage + cNstage + CEA"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/05-5-5_ROC curve-shared strain + shared pathway + clinical indicator_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

mo2_selected = mo %>% 
  select(SampleID, Pre_Op_Tstage, Pre_Op_Nstage, CEA_bin) %>% 
  column_to_rownames("SampleID")

RF_shared_strain_path_meta_o = RF_shared_strain_o %>% 
  rownames_to_column("SampleID") %>% 
  inner_join(RF_shared_path_o %>% 
               select(-TRG_1) %>% 
               rownames_to_column("SampleID"),
             by = "SampleID") %>% 
  inner_join(mo2_selected %>% 
               rownames_to_column("SampleID"),
             by = "SampleID") %>% 
  column_to_rownames("SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1),
         Pre_Op_Tstage = factor(Pre_Op_Tstage),
         Pre_Op_Nstage = factor(Pre_Op_Nstage),
         CEA_bin = factor(CEA_bin))


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_shared_strain_path_meta_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
shared_strain_path_meta_ongoing = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Shared features + cTstage + cNstage + CEA]",
    "\n======================================================\n",
    paste(shared_strain_path_meta_ongoing, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(shared_strain_path_meta_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_ongoing = roc(response = true_labels, 
                      predictor = k_fold_pred, 
                      levels = c("nonCR", "CR"), 
                      direction = "<",
                      quiet = TRUE)

auc_val_ongoing = auc(roc_obj_ongoing) ; auc_val_ongoing # 0.9524  

imp_df_ongoing = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_ongoing

# 시각화
p5_g_roc_impro = ggroc(roc_obj_ongoing, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_strain), " Shared Strains + ",
                         length(shared_path), " Shared Pathways + cTstage + cNstage + CEA"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/05-5-5_ROC curve-shared strain + shared pathway + clinical indicator_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)


# Before, Ongoing Plot
roc_list_all = list(Before = roc_obj_before,
                    Ongoing = roc_obj_ongoing)

p5_g_roc_impro = ggroc(roc_list_all, legacy.axes = TRUE, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Before"  = "#E76F51",
                                "Ongoing" = "#2A9D8F")) +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  annotate("text", x = 0.70, y = 0.15, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#2A9D8F") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_strain), " Shared Strains + ",
                         length(shared_path), " Shared Pathways + cTstage + cNstage + CEA"),
       x = "1 - Specificity", y = "Sensitivity", color = "TNT") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/05-5-5_ROC curve-shared strain + shared pathway + clinical indicator_all.svg",
#        plot = p5_g_roc_impro, width = 6, height = 5)



##### Shared Species + Clinical Indicator

# Before samples

mb2_selected = mb2 %>% 
  select(SampleID, Pre_Op_Tstage, Pre_Op_Nstage, CEA_bin) %>% 
  column_to_rownames("SampleID")

RF_shared_species_meta = RF_shared_species %>% 
  rownames_to_column("SampleID") %>% 
  inner_join(mb2_selected %>% 
               rownames_to_column("SampleID"),
             by = "SampleID") %>% 
  column_to_rownames("SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1),
         Pre_Op_Tstage = factor(Pre_Op_Tstage),
         Pre_Op_Nstage = factor(Pre_Op_Nstage),
         CEA_bin = factor(CEA_bin))


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_shared_species_meta
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
shared_species_meta_before = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Shared features + cTstage + cNstage + CEA]",
    "\n======================================================\n",
    paste(shared_species_meta_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(shared_species_meta_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.9333  

imp_df_before = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_before

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_species), " Shared Species + cTstage + cNstage + CEA"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/05-5-6_ROC curve-shared species + clinical indicator_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

mo2_selected = mo %>% 
  select(SampleID, Pre_Op_Tstage, Pre_Op_Nstage, CEA_bin) %>% 
  column_to_rownames("SampleID")

RF_shared_species_meta_o = RF_shared_species_o %>% 
  rownames_to_column("SampleID") %>% 
  inner_join(mo2_selected %>% 
               rownames_to_column("SampleID"),
             by = "SampleID") %>% 
  column_to_rownames("SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1),
         Pre_Op_Tstage = factor(Pre_Op_Tstage),
         Pre_Op_Nstage = factor(Pre_Op_Nstage),
         CEA_bin = factor(CEA_bin))


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_shared_species_meta_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
shared_species_meta_ongoing = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Shared features + cTstage + cNstage + CEA]",
    "\n======================================================\n",
    paste(shared_species_meta_ongoing, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(shared_species_meta_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_ongoing = roc(response = true_labels, 
                      predictor = k_fold_pred, 
                      levels = c("nonCR", "CR"), 
                      direction = "<",
                      quiet = TRUE)

auc_val_ongoing = auc(roc_obj_ongoing) ; auc_val_ongoing # 0.9841  

imp_df_ongoing = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_ongoing

# 시각화
p5_g_roc_impro = ggroc(roc_obj_ongoing, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_species), " Shared Species + cTstage + cNstage + CEA"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/05-5-6_ROC curve-shared species + clinical indicator_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)


# Before, Ongoing Plot
roc_list_all = list(Before = roc_obj_before,
                    Ongoing = roc_obj_ongoing)

p5_g_roc_impro = ggroc(roc_list_all, legacy.axes = TRUE, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Before"  = "#E76F51",
                                "Ongoing" = "#2A9D8F")) +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  annotate("text", x = 0.70, y = 0.15, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#2A9D8F") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_species), " Shared Species + cTstage + cNstage + CEA"),
       x = "1 - Specificity", y = "Sensitivity", color = "TNT") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/05-5-6_ROC curve-shared species + clinical indicator_all.svg",
#        plot = p5_g_roc_impro, width = 6, height = 5)



##### Shared Species + Shared Pathway + Clinical Indicator

# Before samples

mb2_selected = mb2 %>% 
  select(SampleID, Pre_Op_Tstage, Pre_Op_Nstage, CEA_bin) %>% 
  column_to_rownames("SampleID")

RF_shared_species_path_meta = RF_shared_species %>% 
  rownames_to_column("SampleID") %>% 
  inner_join(RF_shared_path %>% 
               select(-TRG_1) %>% 
               rownames_to_column("SampleID"),
             by = "SampleID") %>% 
  inner_join(mb2_selected %>% 
               rownames_to_column("SampleID"),
             by = "SampleID") %>% 
  column_to_rownames("SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1),
         Pre_Op_Tstage = factor(Pre_Op_Tstage),
         Pre_Op_Nstage = factor(Pre_Op_Nstage),
         CEA_bin = factor(CEA_bin))


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_shared_species_path_meta
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
shared_species_path_meta_before = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Shared features + cTstage + cNstage + CEA]",
    "\n======================================================\n",
    paste(shared_species_path_meta_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(shared_species_path_meta_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.9394  

imp_df_before = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_before

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_species), " Shared Species + ",
                         length(shared_path), " Shared Pathways + cTstage + cNstage + CEA"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/05-5-7_ROC curve-shared species + shared pathway + clinical indicator_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

mo2_selected = mo %>% 
  select(SampleID, Pre_Op_Tstage, Pre_Op_Nstage, CEA_bin) %>% 
  column_to_rownames("SampleID")

RF_shared_species_path_meta_o = RF_shared_species_o %>% 
  rownames_to_column("SampleID") %>% 
  inner_join(RF_shared_path_o %>% 
               select(-TRG_1) %>% 
               rownames_to_column("SampleID"),
             by = "SampleID") %>% 
  inner_join(mo2_selected %>% 
               rownames_to_column("SampleID"),
             by = "SampleID") %>% 
  column_to_rownames("SampleID") %>% 
  mutate(TRG_1 = factor(TRG_1),
         Pre_Op_Tstage = factor(Pre_Op_Tstage),
         Pre_Op_Nstage = factor(Pre_Op_Nstage),
         CEA_bin = factor(CEA_bin))


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_shared_species_path_meta_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
shared_species_path_meta_ongoing = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Shared features + cTstage + cNstage + CEA]",
    "\n======================================================\n",
    paste(shared_species_path_meta_ongoing, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(shared_species_path_meta_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_ongoing = roc(response = true_labels, 
                      predictor = k_fold_pred, 
                      levels = c("nonCR", "CR"), 
                      direction = "<",
                      quiet = TRUE)

auc_val_ongoing = auc(roc_obj_ongoing) ; auc_val_ongoing # 0.9683  

imp_df_ongoing = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_ongoing

# 시각화
p5_g_roc_impro = ggroc(roc_obj_ongoing, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_species), " Shared Species + ",
                         length(shared_path), " Shared Pathways + cTstage + cNstage + CEA"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/05-5-7_ROC curve-shared species + shared pathway + clinical indicator_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)


# Before, Ongoing Plot
roc_list_all = list(Before = roc_obj_before,
                    Ongoing = roc_obj_ongoing)

p5_g_roc_impro = ggroc(roc_list_all, legacy.axes = TRUE, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Before"  = "#E76F51",
                                "Ongoing" = "#2A9D8F")) +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  annotate("text", x = 0.70, y = 0.15, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#2A9D8F") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(shared_species), " Shared Species + ",
                         length(shared_path), " Shared Pathways + cTstage + cNstage + CEA"),
       x = "1 - Specificity", y = "Sensitivity", color = "TNT") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/05-5-7_ROC curve-shared species + shared pathway + clinical indicator_all.svg",
#        plot = p5_g_roc_impro, width = 6, height = 5)



#===========================================================================
### Pathway-associated species ---------------------------------------------

#################### 1. Pathway-associated species ####################

# Pathway 관련 species (36)
path_s = c(unique(c(path_thiamine_species_A$Species,
                    path_sucrose_species_A$Species,
                    path_sucrose4_species_A$Species)),
           
           unique(c(path_thiamine_species_K$Species,
                    path_sucrose_species_K$Species,
                    path_rhamnose_species_K$Species)),
           
           unique(c(path_rhamnose_species_B$Species,
                    path_mep_species_B$Species)),
           
           unique(c(path_rhamnose_species_Bac$Species,
                    path_rhamnose_syn_species_Bac$Species)),
           
           unique(c(path_rhamnose_syn_species_P$Species,
                    path_histidine_species_P$Species)),
           
           unique(c(path_histidine_species_D$Species)),
           
           # Pathway data에서의 이름과 species data에서의 이름이 다름
           "Phocaeicola_vulgatus",
           "Phocaeicola_dorei",
           "Phocaeicola_plebeius",
           "Phocaeicola_massiliensis",
           "Phocaeicola_coprocola",
           "Phocaeicola_coprophilus",
           "Segatella_copri")

# Pathway data에서는 찾았지만 species data에서는 발견 X
# Bacteroides_sp_CAG_144, Bacteroides_sp_OM08_11
# Prevotella_sp_CAG_279, Prevotella_sp_CAG_5226
# Dialister_sp_CAG_357, Dialister_sp_CAG_486


# Before samples

path_s_input = sb %>% 
  filter(Species %in% path_s) %>% 
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

RF_path_species = path_s_input %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Species,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_path_species) # 31 columns (30 species)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_path_species
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
pathway_s_before = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Pathway-associated features]",
    "\n======================================================\n",
    paste(pathway_s_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(pathway_s_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.1879  

imp_df_before = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_before

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(pathway_s_before), " Pathway-associated Species"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/06-5-1_ROC curve-pathway-associated-species_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

path_s_input_o = so %>% 
  filter(Species %in% path_s) %>% 
  pivot_longer(cols = -Species,
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

RF_path_species_o = path_s_input_o %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Species,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_path_species_o) # 31 columns (30 species)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_path_species_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
pathway_s_ongoing = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Pathway-associated features]",
    "\n======================================================\n",
    paste(pathway_s_ongoing, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(pathway_s_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_ongoing = roc(response = true_labels, 
                      predictor = k_fold_pred, 
                      levels = c("nonCR", "CR"), 
                      direction = "<",
                      quiet = TRUE)

auc_val_ongoing = auc(roc_obj_ongoing) ; auc_val_ongoing # 0.3651  

imp_df_ongoing = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_ongoing

# 시각화
p5_g_roc_impro = ggroc(roc_obj_ongoing, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(pathway_s_ongoing), " Pathway-associated Species"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/06-5-1_ROC curve-pathway-associated-species_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)


# Before, Ongoing Plot
roc_list_all = list(Before = roc_obj_before,
                    Ongoing = roc_obj_ongoing)

p5_g_roc_impro = ggroc(roc_list_all, legacy.axes = TRUE, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Before"  = "#E76F51",
                                "Ongoing" = "#2A9D8F")) +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  annotate("text", x = 0.70, y = 0.15, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#2A9D8F") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(pathway_s_ongoing), " Pathway-associated Species"),
       x = "1 - Specificity", y = "Sensitivity", color = "TNT") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/06-5-1_ROC curve-pathway-associated-species_all.svg",
#        plot = p5_g_roc_impro, width = 6, height = 5)



#################### 2. Pathway associated Species x Pathway pair ####################

# Species_Pathway
barplot_path = c("PWY-6892: thiazole component of thiamine diphosphate biosynthesis I",
                 "PWY-621: sucrose degradation III (sucrose invertase)",
                 "PWY-5384: sucrose degradation IV (sucrose phosphorylase)",
                 "RHAMCAT-PWY: L-rhamnose degradation I",
                 "PWY-7560: methylerythritol phosphate pathway II",
                 "DTDPRHAMSYN-PWY: dTDP-&beta;-L-rhamnose biosynthesis",
                 "HISTSYN-PWY: L-histidine biosynthesis")

barplot_species = c(unique(c(path_thiamine_species_A$Species,
                             path_sucrose_species_A$Species,
                             path_sucrose4_species_A$Species)),
                    
                    unique(c(path_thiamine_species_K$Species,
                             path_sucrose_species_K$Species,
                             path_rhamnose_species_K$Species)),
                    
                    unique(c(path_rhamnose_species_B$Species,
                             path_mep_species_B$Species)),
                    
                    unique(c(path_rhamnose_species_Bac$Species,
                             path_rhamnose_syn_species_Bac$Species)),
                    
                    unique(c(path_rhamnose_syn_species_P$Species,
                             path_histidine_species_P$Species)),
                    
                    unique(c(path_histidine_species_D$Species)))


# Before samples (47)

data_config = list(
  path_thiamine_filtered_long     = c(path_thiamine_species_A$Species, path_thiamine_species_K$Species),
  path_sucrose_filtered_long      = c(path_sucrose_species_A$Species, path_sucrose_species_K$Species),
  path_sucrose4_filtered_long     = path_sucrose4_species_A$Species,
  path_rhamnose_filtered_long     = c(path_rhamnose_species_K$Species, path_rhamnose_species_B$Species, path_rhamnose_species_Bac$Species),
  path_mep_filtered_long          = path_mep_species_B$Species,
  path_rhamnose_syn_filtered_long = c(path_rhamnose_syn_species_Bac$Species, path_rhamnose_syn_species_P$Species),
  path_histidine_filtered_long    = c(path_histidine_species_P$Species, path_histidine_species_D$Species)
)

path_species_bind = imap_dfr(data_config, function(species, df_name) {
  get(df_name) %>% 
    filter(Pathway_ID %in% barplot_path) %>% 
    filter(Species %in% species)
}) %>% 
  mutate(species_path = paste0(Species, "_", Pathway_ID))

RF_species_and_path = path_species_bind %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = species_path,
              values_from = Abundance_0.5,
              values_fill = 0) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_species_and_path) # 48 columns (47 species_pathway)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_species_and_path
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
# 공통 feature
species_pathway_pair_before = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Pathway-associated features X Pathway]",
    "\n======================================================\n",
    paste(species_pathway_pair_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(species_pathway_pair_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.3727  

imp_df_before = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_before

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(species_pathway_pair_before), " Species X Pathway Pairs"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/06-5-2_ROC curve-pathway-associated-species x pathway_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples (39)

process_path = function(df, species, mo_data, path) {
  df %>%
    pivot_longer(cols = -Pathway, names_to = "SampleID", values_to = "Abundance_0.5") %>%
    filter(Abundance_0.5 > 0) %>%
    
    left_join(mo_data %>% select(SampleID, Age, Sex, TRG_score, TRG_1, TRG_2, TRG_3, 
                                 Pre_Op_Tstage, Pre_Op_Tstage_bin, 
                                 Pre_Op_Nstage, Pre_Op_Nstage_bin), 
              by = "SampleID") %>%
    
    separate(Pathway, into = c("Pathway_ID", "Taxonomy"), sep = "\\|", fill = "right") %>%
    mutate(Genus = if_else(str_detect(Taxonomy, "g__"),
                           str_remove(str_extract(Taxonomy, "g__[^.]+"), "g__"),
                           "unclassified"),
           Species = if_else(str_detect(Taxonomy, "s__"),
                             str_remove(str_extract(Taxonomy, "s__.+"), "s__"),
                             "unclassified")) %>%
    
    filter(Pathway_ID %in% path) %>%
    filter(Species %in% species) %>%
    mutate(species_path = paste0(Species, "_", Pathway_ID))
}

data_config_o = list(
  list(df = path_thiamine_filtered_0.5,     sp = c(path_thiamine_species_A$Species, path_thiamine_species_K$Species)),
  list(df = path_sucrose_filtered_0.5,      sp = c(path_sucrose_species_A$Species, path_sucrose_species_K$Species)),
  list(df = path_sucrose4_filtered_0.5,     sp = path_sucrose4_species_A$Species),
  list(df = path_rhamnose_filtered_0.5,     sp = c(path_rhamnose_species_K$Species, path_rhamnose_species_B$Species, path_rhamnose_species_Bac$Species)),
  list(df = path_rhamnose_syn_filtered_0.5, sp = c(path_rhamnose_syn_species_Bac$Species, path_rhamnose_syn_species_P$Species)),
  list(df = path_mep_filtered_0.5,          sp = path_mep_species_B$Species),
  list(df = path_histidine_filtered_0.5,    sp = c(path_histidine_species_P$Species, path_histidine_species_D$Species))
)

path_species_bind_on = map_dfr(data_config_o, function(comb) {
  process_path(comb$df, comb$sp, mo, barplot_path)
}) %>% 
  drop_na() %>% 
  relocate(SampleID)

RF_species_and_path_o = path_species_bind_on %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = species_path,
              values_from = Abundance_0.5,
              values_fill = 0) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  select(-SampleID) ; dim(RF_species_and_path_o) # 40 columns (39 species_pathway)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_species_and_path_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
species_pathway_pair_ongoing = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Pathway-associated features X Pathway]",
    "\n======================================================\n",
    paste(species_pathway_pair_ongoing, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(species_pathway_pair_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_ongoing = roc(response = true_labels, 
                      predictor = k_fold_pred, 
                      levels = c("nonCR", "CR"), 
                      direction = "<",
                      quiet = TRUE)

auc_val_ongoing = auc(roc_obj_ongoing) ; auc_val_ongoing # 0.5635  

imp_df_ongoing = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_ongoing

# 시각화
p5_g_roc_impro = ggroc(roc_obj_ongoing, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(species_pathway_pair_ongoing), " Species X Pathway Pairs"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/06-5-2_ROC curve-pathway-associated-species x pathway_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)


# Before, Ongoing Plot
roc_list_all = list(Before = roc_obj_before,
                    Ongoing = roc_obj_ongoing)

p5_g_roc_impro = ggroc(roc_list_all, legacy.axes = TRUE, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Before"  = "#E76F51",
                                "Ongoing" = "#2A9D8F")) +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  annotate("text", x = 0.70, y = 0.15, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#2A9D8F") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = paste0(length(species_pathway_pair_before), " Species X Pathway Pairs (Before)", "\n",
                         length(species_pathway_pair_ongoing), " Species X Pathway Pairs (Ongoing)"),
       x = "1 - Specificity", y = "Sensitivity", color = "TNT") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/06-5-2_ROC curve-pathway-associated-species x pathway_all.svg",
#        plot = p5_g_roc_impro, width = 6, height = 5)



#################### 3. Specific Species  ####################

### Anaerostipes hadrus

# Before samples

RF_path_species_A_hadrus = path_s_input %>% 
  filter(Species == "Anaerostipes_hadrus") %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Species,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_path_species_A_hadrus)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_path_species_A_hadrus
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
A_hadrus_before = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Selected features]",
    "\n======================================================\n",
    paste(A_hadrus_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(A_hadrus_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.6212  

imp_df_before = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_before

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = "Anaerostipes hadrus",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/06-5-3_ROC curve-Anaerostipes hadrus_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

RF_path_species_A_hadrus_o = path_s_input_o %>% 
  filter(Species == "Anaerostipes_hadrus") %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Species,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_path_species_A_hadrus_o)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_path_species_A_hadrus_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
A_hadrus_ongoing = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Selected features]",
    "\n======================================================\n",
    paste(A_hadrus_ongoing, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(A_hadrus_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_ongoing = roc(response = true_labels, 
                      predictor = k_fold_pred, 
                      levels = c("nonCR", "CR"), 
                      direction = "<",
                      quiet = TRUE)

auc_val_ongoing = auc(roc_obj_ongoing) ; auc_val_ongoing # 0.4603  

imp_df_ongoing = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_ongoing

# 시각화
p5_g_roc_impro = ggroc(roc_obj_ongoing, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = "Anaerostipes hadrus",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/06-5-3_ROC curve-Anaerostipes hadrus_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)


# Before, Ongoing Plot
roc_list_all = list(Before = roc_obj_before,
                    Ongoing = roc_obj_ongoing)

p5_g_roc_impro = ggroc(roc_list_all, legacy.axes = TRUE, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Before"  = "#E76F51",
                                "Ongoing" = "#2A9D8F")) +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  annotate("text", x = 0.70, y = 0.15, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#2A9D8F") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = "Anaerostipes hadrus",
       x = "1 - Specificity", y = "Sensitivity", color = "TNT") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/06-5-3_ROC curve-Anaerostipes hadrus_all.svg",
#        plot = p5_g_roc_impro, width = 6, height = 5)



### Bacteroides vulgatus

# Before samples

RF_path_species_B_vulgatus = path_s_input %>% 
  filter(Species == "Phocaeicola_vulgatus") %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Species,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_path_species_B_vulgatus)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_path_species_B_vulgatus
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
B_vulgatus_before = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Selected features]",
    "\n======================================================\n",
    paste(B_vulgatus_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(B_vulgatus_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.6091  

imp_df_before = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_before

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = "Bacteroides vulgatus (Phocaeicola vulgatus)",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/06-5-3_ROC curve-Bacteroides vulgatus_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

RF_path_species_B_vulgatus_o = path_s_input_o %>% 
  filter(Species == "Phocaeicola_vulgatus") %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Species,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_path_species_B_vulgatus_o)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_path_species_B_vulgatus_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
B_vulgatus_ongoing = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Selected features]",
    "\n======================================================\n",
    paste(B_vulgatus_ongoing, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(B_vulgatus_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_ongoing = roc(response = true_labels, 
                      predictor = k_fold_pred, 
                      levels = c("nonCR", "CR"), 
                      direction = "<",
                      quiet = TRUE)

auc_val_ongoing = auc(roc_obj_ongoing) ; auc_val_ongoing # 0.4444  

imp_df_ongoing = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_ongoing

# 시각화
p5_g_roc_impro = ggroc(roc_obj_ongoing, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = "Bacteroides vulgatus (Phocaeicola vulgatus)",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/06-5-3_ROC curve-Bacteroides vulgatus_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)


# Before, Ongoing Plot
roc_list_all = list(Before = roc_obj_before,
                    Ongoing = roc_obj_ongoing)

p5_g_roc_impro = ggroc(roc_list_all, legacy.axes = TRUE, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Before"  = "#E76F51",
                                "Ongoing" = "#2A9D8F")) +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  annotate("text", x = 0.70, y = 0.15, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#2A9D8F") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = "Bacteroides vulgatus (Phocaeicola vulgatus)",
       x = "1 - Specificity", y = "Sensitivity", color = "TNT") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/06-5-3_ROC curve-Bacteroides vulgatus_all.svg",
#        plot = p5_g_roc_impro, width = 6, height = 5)



### Prevotella copri

# Before samples

RF_path_species_P_copri = path_s_input %>% 
  filter(Species == "Segatella_copri") %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Species,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_path_species_P_copri)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_path_species_P_copri
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
P_copri_before = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Selected features]",
    "\n======================================================\n",
    paste(P_copri_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(P_copri_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.3303  

imp_df_before = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_before

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = "Prevotella copri (Segatella copri)",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/06-5-3_ROC curve-Prevotella copri_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

RF_path_species_P_copri_o = path_s_input_o %>% 
  filter(Species == "Segatella_copri") %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Species,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_path_species_P_copri_o)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_path_species_P_copri_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
P_copri_ongoing = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Selected features]",
    "\n======================================================\n",
    paste(P_copri_ongoing, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(P_copri_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_ongoing = roc(response = true_labels, 
                      predictor = k_fold_pred, 
                      levels = c("nonCR", "CR"), 
                      direction = "<",
                      quiet = TRUE)

auc_val_ongoing = auc(roc_obj_ongoing) ; auc_val_ongoing # 0.6587  

imp_df_ongoing = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_ongoing

# 시각화
p5_g_roc_impro = ggroc(roc_obj_ongoing, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = "Prevotella copri (Segatella copri)",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/06-5-3_ROC curve-Prevotella copri_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)


# Before, Ongoing Plot
roc_list_all = list(Before = roc_obj_before,
                    Ongoing = roc_obj_ongoing)

p5_g_roc_impro = ggroc(roc_list_all, legacy.axes = TRUE, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Before"  = "#E76F51",
                                "Ongoing" = "#2A9D8F")) +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  annotate("text", x = 0.70, y = 0.15, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#2A9D8F") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = "Prevotella copri (Segatella copri)",
       x = "1 - Specificity", y = "Sensitivity", color = "TNT") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/06-5-3_ROC curve-Prevotella copri_all.svg",
#        plot = p5_g_roc_impro, width = 6, height = 5)



### Anaerostipes hadrus + Bacteroides vulgatus

# Before samples

RF_path_species_A_hadrus_B_vulgatus = path_s_input %>% 
  filter(Species %in% c("Anaerostipes_hadrus", "Phocaeicola_vulgatus")) %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Species,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_path_species_A_hadrus_B_vulgatus)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_path_species_A_hadrus_B_vulgatus
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
A_hadrus_B_vulgatus_before = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Selected features]",
    "\n======================================================\n",
    paste(A_hadrus_B_vulgatus_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(A_hadrus_B_vulgatus_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.6364  

imp_df_before = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_before

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = "Anaerostipes hadrus + Bacteroides vulgatus (Phocaeicola vulgatus)",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/06-5-4_ROC curve-Anaerostipes hadrus + Bacteroides vulgatus_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

RF_path_species_A_hadrus_B_vulgatus_o = path_s_input_o %>% 
  filter(Species %in% c("Anaerostipes_hadrus", "Phocaeicola_vulgatus")) %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Species,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_path_species_A_hadrus_B_vulgatus_o)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_path_species_A_hadrus_B_vulgatus_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
A_hadrus_B_vulgatus_ongoing = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Selected features]",
    "\n======================================================\n",
    paste(A_hadrus_B_vulgatus_ongoing, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(A_hadrus_B_vulgatus_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_ongoing = roc(response = true_labels, 
                      predictor = k_fold_pred, 
                      levels = c("nonCR", "CR"), 
                      direction = "<",
                      quiet = TRUE)

auc_val_ongoing = auc(roc_obj_ongoing) ; auc_val_ongoing # 0.3175  

imp_df_ongoing = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_ongoing

# 시각화
p5_g_roc_impro = ggroc(roc_obj_ongoing, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = "Anaerostipes hadrus + Bacteroides vulgatus (Phocaeicola vulgatus)",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/06-5-4_ROC curve-Anaerostipes hadrus + Bacteroides vulgatus_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)


# Before, Ongoing Plot
roc_list_all = list(Before = roc_obj_before,
                    Ongoing = roc_obj_ongoing)

p5_g_roc_impro = ggroc(roc_list_all, legacy.axes = TRUE, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Before"  = "#E76F51",
                                "Ongoing" = "#2A9D8F")) +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  annotate("text", x = 0.70, y = 0.15, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#2A9D8F") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = "Anaerostipes hadrus + Bacteroides vulgatus (Phocaeicola vulgatus)",
       x = "1 - Specificity", y = "Sensitivity", color = "TNT") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/06-5-4_ROC curve-Anaerostipes hadrus + Bacteroides vulgatus_all.svg",
#        plot = p5_g_roc_impro, width = 6, height = 5)



### Anaerostipes hadrus + Bacteroides vulgatus + Prevotella copri

# Before samples

RF_path_species_A_hadrus_B_vulgatus_P_copri = path_s_input %>% 
  filter(Species %in% c("Anaerostipes_hadrus", "Phocaeicola_vulgatus", "Segatella_copri")) %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Species,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_path_species_A_hadrus_B_vulgatus_P_copri)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_path_species_A_hadrus_B_vulgatus_P_copri
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
A_hadrus_B_vulgatus_P_copri_before = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Selected features]",
    "\n======================================================\n",
    paste(A_hadrus_B_vulgatus_P_copri_before, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 5
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(A_hadrus_B_vulgatus_P_copri_before, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_before = roc(response = true_labels, 
                     predictor = k_fold_pred, 
                     levels = c("nonCR", "CR"), 
                     direction = "<",
                     quiet = TRUE)

auc_val_before = auc(roc_obj_before) ; auc_val_before # 0.4273  

imp_df_before = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_before

# 시각화
p4_g_roc_impro = ggroc(roc_obj_before, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = "Anaerostipes hadrus + Bacteroides vulgatus (Phocaeicola vulgatus) + Prevotella copri (Segatella copri)",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p4_g_roc_impro

# ggsave("figure/06-5-4_ROC curve-Anaerostipes hadrus + Bacteroides vulgatus + Prevotella copri_before.svg",
#        plot = p4_g_roc_impro, width = 5, height = 5)


# Ongoing samples

RF_path_species_A_hadrus_B_vulgatus_P_copri_o = path_s_input_o %>% 
  filter(Species %in% c("Anaerostipes_hadrus", "Phocaeicola_vulgatus", "Segatella_copri")) %>% 
  pivot_wider(id_cols = c(SampleID, TRG_1),
              names_from = Species,
              values_from = abundance) %>% 
  mutate(TRG_1 = factor(TRG_1)) %>% 
  column_to_rownames(var = "SampleID") ; dim(RF_path_species_A_hadrus_B_vulgatus_P_copri_o)


### Stratified K-Fold Cross Validation ###

# Wilcoxon 기반 feature selection
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)

# 1. 데이터 준비
RF_work = RF_path_species_A_hadrus_B_vulgatus_P_copri_o
colnames(RF_work) = make.names(colnames(RF_work))
RF_work$TRG_1 = factor(RF_work$TRG_1, levels = c("nonCR", "CR"))

# 2. Feature selection
A_hadrus_B_vulgatus_P_copri_ongoing = names(RF_work)[names(RF_work) != "TRG_1"]

cat("\n\n============================================",
    "\n [Selected features]",
    "\n======================================================\n",
    paste(A_hadrus_B_vulgatus_P_copri_ongoing, collapse = "\n"))

# 3. Stratified K-Fold CV
N = nrow(RF_work)
k_folds = 3
importance_list = list()
k_fold_pred = numeric(N)
true_labels = RF_work$TRG_1

# Stratified Fold 생성
set.seed(2026)

folds = createFolds(RF_work$TRG_1, k = k_folds, list = TRUE, returnTrain = FALSE)

# 진행바
pb = txtProgressBar(min = 0, max = k_folds, style = 3)

cat(paste0("\n>> Stratified ", k_folds, "-Fold CV 수행 중...\n"))

for (i in 1:k_folds) {
  
  # (1) 데이터 분할
  test_idx = folds[[i]]
  train_data = RF_work[-test_idx, ]
  test_data  = RF_work[test_idx, ]
  
  # (2) 모델 학습
  # formula 생성: TRG_1 ~ top_feature1 + top_feature2 ...
  rf_formula = as.formula(paste("TRG_1 ~", paste(A_hadrus_B_vulgatus_P_copri_ongoing, collapse = "+")))
  
  # Class Imbalance 해결을 위한 sampsize 조정 (옵션)
  # min_class_n = min(table(train_data$TRG_1))
  
  rf_model = randomForest(rf_formula, 
                          data = train_data, 
                          ntree = 1000,
                          importance = TRUE)
  
  importance_list[[i]] = as.data.frame(importance(rf_model))
  # strata = train_data$TRG_1, # 불균형 심하면 주석 해제
  # sampsize = rep(min_class_n, 2)) # 불균형 심하면 주석 해제
  
  # (3) 예측
  k_fold_pred[test_idx] = predict(rf_model, test_data, type = "prob")[, "CR"]
  
  setTxtProgressBar(pb, i)
}
close(pb)

# 4. 결과 확인
roc_obj_ongoing = roc(response = true_labels, 
                      predictor = k_fold_pred, 
                      levels = c("nonCR", "CR"), 
                      direction = "<",
                      quiet = TRUE)

auc_val_ongoing = auc(roc_obj_ongoing) ; auc_val_ongoing # 0.4762  

imp_df_ongoing = do.call(rbind, importance_list) %>% 
  as.data.frame() %>% 
  mutate(Feature = rep(rownames(importance_list[[1]]), length(importance_list))) %>% 
  group_by(Feature) %>% 
  summarise(Avg_MeanDecreaseAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE),
            Avg_MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Avg_MeanDecreaseAccuracy)) ; imp_df_ongoing

# 시각화
p5_g_roc_impro = ggroc(roc_obj_ongoing, legacy.axes = TRUE, size = 1.2, colour = "#E76F51") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = "Anaerostipes hadrus + Bacteroides vulgatus (Phocaeicola vulgatus) + Prevotella copri (Segatella copri)",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/06-5-4_ROC curve-Anaerostipes hadrus + Bacteroides vulgatus + Prevotella copri_ongoing.svg",
#        plot = p5_g_roc_impro, width = 5, height = 5)


# Before, Ongoing Plot
roc_list_all = list(Before = roc_obj_before,
                    Ongoing = roc_obj_ongoing)

p5_g_roc_impro = ggroc(roc_list_all, legacy.axes = TRUE, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Before"  = "#E76F51",
                                "Ongoing" = "#2A9D8F")) +
  annotate("text", x = 0.70, y = 0.20, 
           label = paste0("AUC = ", round(auc_val_before, 3)), 
           size = 5, fontface = "bold", color = "#E76F51") +
  annotate("text", x = 0.70, y = 0.15, 
           label = paste0("AUC = ", round(auc_val_ongoing, 3)), 
           size = 5, fontface = "bold", color = "#2A9D8F") +
  labs(title = "Improved Prediction (Feature Selection)",
       subtitle = "Anaerostipes hadrus + Bacteroides vulgatus (Phocaeicola vulgatus) + Prevotella copri (Segatella copri)",
       x = "1 - Specificity", y = "Sensitivity", color = "TNT") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust=0.5)) ; p5_g_roc_impro

# ggsave("figure/06-5-4_ROC curve-Anaerostipes hadrus + Bacteroides vulgatus + Prevotella copri_all.svg",
#        plot = p5_g_roc_impro, width = 6, height = 5)




# Remove unnecessary objects
rm(list = ls(pattern = "^p4_")); rm(list = ls(pattern = "^p5_"))
save.image(file = "input/R_image/7-5. after-random-forest-final_Stratified K-Fold.RData")



