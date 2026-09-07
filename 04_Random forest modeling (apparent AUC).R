# ==================== 0. Setting ====================

# ---------- 0-1. Directory setting ----------
setwd("C:/Users/sissi/Desktop/NAYUN/CRC TNT/")
options(java.parameters = "-Xmx64g", stringsAsFactors = F)

library(dplyr)
library(tidyr)
library(tibble)
library(pROC)
library(randomForest)
library(rsample)

# ---------- 0-2. Rdata ----------
# Rdata
# load("RData/05_RF input 260722.RData")
load("RData/05_RF modeling 260804.RData")

# ---------- 0-3. Function ----------
# fastshap API function
shap_fun <- function(object, newdata) predict(object, newdata = newdata, type = "prob")[, "CR"]

# ---------- 0-4. Palette ----------
cols_input <- c("Strain"     = "#1F78B4",
                "Species"    = "#2CB1C9",
                "Genus"      = "#5FC8B6",
                "Family"     = "#98DDD9",
                "Phylum"     = "#C7E6F2",
                "Pathway"    = "#F2B233",
                "KO"         = "#D48806",
                "Metabolite" = "#5FBF66",
                "CI"         = "#AEAEAE")

# ==================== 1. Setting ====================

# ---------- 1-1. Common setting ----------
SEED      <- 123  # 재현성 시드
N_TREE    <- 1000 # Random forest 트리 개수
V_FOLD    <- 5    # Stratified CV fold 수 (group_vfold_cv 사용)
NSIM_SHAP <- 100  # fastshap Monte Carlo 반복 횟수

# Grid feature selection cutoff
grid_cutoff <- expand.grid(prev_cutoff   = c(0.10, 0.20, 0.30, 0.40, 0.50), # Prevalence: Prev_CR 또는 Prev_nonCR 중 하나 이상에서 기준 만족
                           d_cutoff      = c(0.10, 0.20, 0.30),             # Effect size: |D| = |2*AUC-1|
                           cons_cutoff   = c(0.70, 0.80, 0.90),             # Consistency: direction consistency, Before/Ongoing 방향 일치 비율
                           pval_cutoff   = 0.01,                            # Binomial p-value: consistency ratio에 대한 p-value
                           wilcox_cutoff = c(0.01, 0.05, 0.10, 0.20),       # Wilcoxon p-value: feature abundance에 대한 p-value
                           n_top         = c(10, 20, 30, 40, 50),           # Features: random forest 모델링에 사용할 feature 개수
                           stringsAsFactors = FALSE)

# ---------- 1-2. Input setting ----------
# m, mb_pair, mo_pair 는 이미 load된 상태 (SNU_ID, SampleID, TRG_1 포함)

b_mat    <- met_b_mat   # Before 시점, rownames = SNU_ID
o_mat    <- met_o_mat   # Ongoing 시점, rownames = SNU_ID
samp_mat <- met_mat     # 샘플 단위 전체 매트릭스, rownames = SampleID

# 세 매트릭스에 공통으로 존재하는 feature(컬럼)만 사용
common_cols <- Reduce(intersect, list(colnames(samp_mat), colnames(b_mat), colnames(o_mat)))
b_mat       <- b_mat[,    common_cols, drop = FALSE]
o_mat       <- o_mat[,    common_cols, drop = FALSE]
samp_mat    <- samp_mat[, common_cols, drop = FALSE]
feat_cands  <- common_cols

patient_meta <- mm %>%
  select(SNU_ID, TRG_1) %>% distinct() %>%
  mutate(SNU_ID = as.character(SNU_ID))

sample_meta <- mm %>%
  select(SampleID, SNU_ID, TRG_1) %>%
  mutate(SNU_ID = as.character(SNU_ID))

# ==================== 2. Modeling ====================

# ---------- 2-1. Feature scoring ----------
# 처음 1회만 계산하는 global feature scoring

cr_all  <- patient_meta %>% filter(TRG_1 == "CR")    %>% pull(SNU_ID)
ncr_all <- patient_meta %>% filter(TRG_1 == "nonCR") %>% pull(SNU_ID)

cr_b  <- intersect(cr_all,  rownames(b_mat))
ncr_b <- intersect(ncr_all, rownames(b_mat))
cr_o  <- intersect(cr_all,  rownames(o_mat))
ncr_o <- intersect(ncr_all, rownames(o_mat))

score_list <- vector("list", length(feat_cands))

for (i in seq_along(feat_cands)) {
  
  feat <- feat_cands[i]
  
  # ── Before / Ongoing 각각에서 CR - nonCR pairwise 차이 ──────────────────────
  b_d <- as.vector(outer(b_mat[cr_b, feat], b_mat[ncr_b, feat], "-"))
  o_d <- as.vector(outer(o_mat[cr_o, feat], o_mat[ncr_o, feat], "-"))
  
  b_pos <- sum(b_d > 0); b_neg <- sum(b_d < 0); b_n <- b_pos + b_neg
  o_pos <- sum(o_d > 0); o_neg <- sum(o_d < 0); o_n <- o_pos + o_neg
  
  all_cr  <- c(b_mat[cr_b,  feat], o_mat[cr_o,  feat])
  all_ncr <- c(b_mat[ncr_b, feat], o_mat[ncr_o, feat])
  
  # ── Prevalence (양수 비율) ───────────────────────────────────────────────
  prev_cr  <- mean(all_cr  > 0, na.rm = TRUE)
  prev_ncr <- mean(all_ncr > 0, na.rm = TRUE)
  
  # ── AUC -> D (effect size, = 2*AUC - 1) ─────────────────────────────────
  mat_gt <- outer(all_cr, all_ncr, ">")
  mat_eq <- outer(all_cr, all_ncr, "==")
  mat_00 <- outer(all_cr == 0, all_ncr == 0, "&")
  valid  <- !mat_00
  total  <- sum(valid)
  auc <- if (total > 0) (sum(mat_gt & valid) + 0.5 * sum(mat_eq & valid & !mat_00)) / total else NA_real_
  D   <- if (!is.na(auc)) 2 * auc - 1 else NA_real_
  
  # ── Wilcoxon p-value ──────────────────────────────────────────────────
  wilcox_p <- tryCatch(
    suppressWarnings(wilcox.test(all_cr, all_ncr, exact = FALSE)$p.value),
    error = function(e) NA_real_
  )
  
  # ── Consistency / binomial p-value ───────────────────────────────────
  # Before, Ongoing 두 시점의 방향(부호)이 일치할 때만 유효한 feature로 인정
  if (b_n == 0 || o_n == 0) {
    cons <- NA_real_; pval <- NA_real_
  } else {
    b_dir <- ifelse(b_pos >= b_neg, 1L, -1L)
    o_dir <- ifelse(o_pos >= o_neg, 1L, -1L)
    if (b_dir != o_dir) {
      cons <- NA_real_; pval <- NA_real_
    } else {
      b_pv <- binom.test(max(b_pos, b_neg), b_n, 0.5, "greater")$p.value
      o_pv <- binom.test(max(o_pos, o_neg), o_n, 0.5, "greater")$p.value
      pval <- sqrt(b_pv * o_pv)                                     # 두 시점 p-value 기하평균
      cons <- (max(b_pos, b_neg) / b_n + max(o_pos, o_neg) / o_n) / 2
    }
  }
  
  score_list[[i]] <- data.frame(
    Feature = feat, Prev_CR = prev_cr, Prev_nonCR = prev_ncr,
    D = D, Consistency = cons, Pvalue = pval, Wilcoxon_p = wilcox_p
  )
}

scores_global <- bind_rows(score_list)
# columns: Feature, Prev_CR, Prev_nonCR, D, Consistency, Pvalue, Wilcoxon_p

# ---------- 2-2. Stratified 5-fold CV ----------
# 처음 1회만 계산 후 모든 grid 조합/최종 실행에 사용

RF_data <- samp_mat %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(sample_meta %>% select(SampleID, SNU_ID, TRG_1), by = "SampleID") %>%
  column_to_rownames("SampleID")

set.seed(SEED)
folds <- rsample::group_vfold_cv(RF_data, group = SNU_ID, v = V_FOLD, strata = TRG_1)

N <- nrow(RF_data)
true_labels <- RF_data$TRG_1

# ---------- 2-3. Grid search ----------
# 조합별 AUC 저장
grid_auc <- rep(NA_real_, nrow(grid_cutoff))

for (gi in seq_len(nrow(grid_cutoff))) {
  
  g <- grid_cutoff[gi, ]
  
  # grid search에서의 cutoff 조합으로 global feature selection
  sel_g <- scores_global %>%
    filter(!is.na(D), !is.na(Consistency), !is.na(Pvalue)) %>%
    filter(Prev_CR >= g$prev_cutoff | Prev_nonCR >= g$prev_cutoff) %>%
    filter(abs(D) >= g$d_cutoff, Consistency >= g$cons_cutoff, Pvalue < g$pval_cutoff) %>%
    filter(!is.na(Wilcoxon_p), Wilcoxon_p < g$wilcox_cutoff) %>%
    arrange(Pvalue)
  
  if (nrow(sel_g) > g$n_top) sel_g <- sel_g %>% slice_head(n = g$n_top)
  
  fc_g <- sel_g$Feature
  if (length(fc_g) == 0) next
  
  # Cutoff 조합이 정해진 후 seed 고정 (feature 조합에 따른 차이만 확인하기 위함)
  set.seed(SEED)
  
  # 정해진 feature set으로 5-fold CV
  k_pred_g <- rep(NA_real_, N)
  names(k_pred_g) <- rownames(RF_data)
  
  for (i in seq_len(nrow(folds))) {
    tr_data <- rsample::analysis(folds$splits[[i]])
    te_data <- rsample::assessment(folds$splits[[i]])
    
    mod <- randomForest::randomForest(x = tr_data[, fc_g, drop = FALSE],
                                      y = factor(tr_data$TRG_1, levels = c("nonCR", "CR")),
                                      ntree = N_TREE)
    pp  <- predict(mod, newdata = te_data[, fc_g, drop = FALSE], type = "prob")[, "CR"]
    idx <- match(rownames(te_data), rownames(RF_data))
    k_pred_g[idx] <- pp
  }
  
  valid_g <- !is.na(k_pred_g)
  if (length(unique(true_labels[valid_g])) < 2) next
  
  roc_g <- tryCatch(pROC::roc(true_labels[valid_g], k_pred_g[valid_g],
                              levels = c("nonCR", "CR"), direction = "<", quiet = TRUE),
                    error = function(e) NULL)
  if (!is.null(roc_g)) grid_auc[gi] <- as.numeric(pROC::auc(roc_g))
  
  cat(sprintf("  [%d/%d] n_feat=%d AUC=%s\n",
              gi, nrow(grid_cutoff), length(fc_g),
              ifelse(is.na(grid_auc[gi]), "NA", sprintf("%.3f", grid_auc[gi]))))
}

grid_result <- cbind(grid_cutoff, AUC = grid_auc)

# ---------- 2-4. Best cutoff ----------
best_g <- grid_result %>%
  filter(!is.na(AUC)) %>%
  arrange(desc(AUC), n_top) %>%
  slice(1)

cat("\n>> Best cutoff 조합:\n"); print(best_g)

# ---------- 2-5. Final model ----------
# 최종적으로 정해진 feature set
sel_final <- scores_global %>%
  filter(!is.na(D), !is.na(Consistency), !is.na(Pvalue)) %>%
  filter(Prev_CR >= best_g$prev_cutoff | Prev_nonCR >= best_g$prev_cutoff) %>%
  filter(abs(D) >= best_g$d_cutoff, Consistency >= best_g$cons_cutoff,
         Pvalue < best_g$pval_cutoff) %>%
  filter(!is.na(Wilcoxon_p), Wilcoxon_p < best_g$wilcox_cutoff) %>%
  arrange(Pvalue)

if (nrow(sel_final) > best_g$n_top) sel_final <- sel_final %>% slice_head(n = best_g$n_top)

# 모델 학습에 사용할 feature 이름만 추출
fc <- sel_final$Feature
cat(sprintf("최종 선택된 feature 수: %d\n", length(fc)))

# Best cutoff가 정해진 후 seed 고정 (grid search와 최대한 유사한 결과를 내기 위함)
# 여기서는 importance와 SHAP 계산이 추가되기 때문에 난수를 추가로 소비
# 완전히 똑같은 결과는 아닐 수 있으나 "동일한 규칙으로 고정한 seed"
set.seed(SEED)

k_pred <- rep(NA_real_, N)
names(k_pred) <- rownames(RF_data)

importance_log <- vector("list", nrow(folds))   # fold별 MeanDecreaseAccuracy/Gini
shap_log       <- vector("list", nrow(folds))   # fold별 SHAP (long format)

for (i in seq_len(nrow(folds))) {
  tr_data <- rsample::analysis(folds$splits[[i]])
  te_data <- rsample::assessment(folds$splits[[i]])
  
  mod <- randomForest::randomForest(x = tr_data[, fc, drop = FALSE],
                                    y = factor(tr_data$TRG_1, levels = c("nonCR", "CR")),
                                    ntree = N_TREE,
                                    importance = TRUE)
  pp  <- predict(mod, newdata = te_data[, fc, drop = FALSE], type = "prob")[, "CR"]
  idx <- match(rownames(te_data), rownames(RF_data))
  k_pred[idx] <- pp
  
  # Feature importance (fold별)
  importance_log[[i]] <- as.data.frame(randomForest::importance(mod)) %>%
    rownames_to_column("Feature") %>%
    transmute(Feature, MDA = MeanDecreaseAccuracy, MDG = MeanDecreaseGini, Fold = i)
  
  # SHAP (fold별 test set 대상)
  shap_val <- fastshap::explain(object       = mod,
                                X            = tr_data[, fc, drop = FALSE],
                                newdata      = te_data[, fc, drop = FALSE],
                                pred_wrapper = shap_fun,
                                nsim         = NSIM_SHAP)
  
  shap_log[[i]] <- as.data.frame(shap_val) %>%
    mutate(SampleID = rownames(te_data)) %>%
    pivot_longer(-SampleID, names_to = "Feature", values_to = "SHAP") %>%
    left_join(te_data[, fc, drop = FALSE] %>% as.data.frame() %>%
                mutate(SampleID = rownames(te_data)) %>%
                pivot_longer(-SampleID, names_to = "Feature", values_to = "Abundance"),
              by = c("SampleID", "Feature")) %>%
    mutate(Fold = i)
  
  cat(sprintf("  Fold %d/%d 완료 (importance+SHAP 계산 포함)\n", i, nrow(folds)))
}

# Fold별 결과 취합
importance_all <- bind_rows(importance_log)
shap_all       <- bind_rows(shap_log)

# Feature별 mean |SHAP| (fold/sample 모두 pooled)
shap_summary <- shap_all %>%
  group_by(Feature) %>%
  summarise(MeanAbsSHAP = mean(abs(SHAP), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MeanAbsSHAP))

# Fold마다 스케일이 다를 수 있어 fold 내 percent_rank로 정규화 후 평균
importance_scaled <- importance_all %>%
  group_by(Fold) %>%
  mutate(MDA_scaled = percent_rank(MDA), MDG_scaled = percent_rank(MDG)) %>%
  ungroup()

importance_summary <- importance_scaled %>%
  group_by(Feature) %>%
  summarise(MDA = mean(MDA_scaled), MDG = mean(MDG_scaled), .groups = "drop")

feature_summary <- shap_summary %>%
  left_join(importance_summary, by = "Feature") %>%
  left_join(sel_final %>% select(Feature, Consistency, Pvalue), by = "Feature") %>%
  arrange(desc(MeanAbsSHAP))

valid_idx <- !is.na(k_pred)

roc_obj <- pROC::roc(true_labels[valid_idx], k_pred[valid_idx],
                     levels = c("nonCR", "CR"), direction = "<", quiet = TRUE)

apparent_auc <- as.numeric(pROC::auc(roc_obj))

cat(sprintf("\n[RF_Global | Direction] Apparent AUC = %.3f\n", apparent_auc))
cat(sprintf("Best cutoff: prev=%.2f d=%.2f cons=%.2f pval=%.2f wilcox=%.2f n_top=%d\n",
            best_g$prev_cutoff, best_g$d_cutoff, best_g$cons_cutoff,
            best_g$pval_cutoff, best_g$wilcox_cutoff, best_g$n_top))

# ==================== 3. Results ====================

# ---------- 3-1. Saving ----------
# Input 종류 바꿔서 저장
input_name <- "Metabolite"

pred_df <- data.frame(SampleID = rownames(RF_data),
                      SNU_ID   = RF_data$SNU_ID,
                      TRG_1    = true_labels,
                      Pred_CR  = k_pred)

# ── 이번 실행(=이 input) 결과를 하나의 list로 묶음 ────────────────────────────
result_this_input <- list(Input             = input_name,
                          Model             = "RF_Global",
                          FS                = "Direction",
                          Grid              = grid_result,        # 3600개 cutoff 조합별 AUC 전체
                          Best_criteria     = best_g,             # grid 중 채택된 best cutoff
                          Scores_global     = scores_global,      # 전체 환자 대상 Direction scoring 원본 (전체 feature)
                          Selected_features = fc,                 # best cutoff로 최종 선택된 feature 이름 벡터
                          Selected_scores   = sel_final,          # 위 feature들의 Consistency/Pvalue
                          AUC               = apparent_auc,       # apparent AUC (5-fold pooled)
                          ROC               = roc_obj,
                          Pred              = pred_df,            # 샘플별 예측 확률
                          Importance        = importance_all,     # fold별 importance (raw)
                          Importance_summary= importance_summary, # feature별 평균(percent_rank) importance
                          SHAP              = shap_all,           # fold별 SHAP (long format)
                          SHAP_summary      = shap_summary,       # feature별 mean |SHAP|
                          Feature_summary   = feature_summary)    # SHAP+importance 합쳐 정렬한 요약표

# 여러 input(omics)을 돌아가며 실행할 때 all_results에 누적
# (스크립트를 input만 바꿔 다시 돌리면 all_results[[input_name]]에 계속 쌓임)
if (!exists("all_results")) all_results <- list()
all_results[[input_name]] <- result_this_input

### CI AUC 계산 결과 추가
# Clinical stage와 CEA 값으로 baseline 확인 목적
ci_features <- c("Pre_Op_Tstage", "Pre_Op_Nstage", "CEA_bin")

RF_ci <- m %>%
  select(SNU_ID, TRG_1, all_of(ci_features)) %>%
  mutate(SNU_ID = as.character(SNU_ID),
         TRG_1  = factor(TRG_1, levels = c("nonCR", "CR"))) %>%
  distinct(SNU_ID, .keep_all = TRUE)

set.seed(SEED)
folds_ci <- rsample::group_vfold_cv(RF_ci, group = SNU_ID, v = V_FOLD, strata = TRG_1)

k_pred_ci <- rep(NA_real_, nrow(RF_ci))
names(k_pred_ci) <- RF_ci$SNU_ID

# Best cutoff 없이 한 번만 seed 고정
set.seed(SEED)
for (i in seq_len(nrow(folds_ci))) {
  
  tr_data <- rsample::analysis(folds_ci$splits[[i]])
  te_data <- rsample::assessment(folds_ci$splits[[i]])
  
  mod_ci <- randomForest::randomForest(x = tr_data[, ci_features, drop = FALSE],
                                       y = tr_data$TRG_1,
                                       ntree = N_TREE)
  
  pp  <- predict(mod_ci, newdata = te_data[, ci_features, drop = FALSE], type = "prob")[, "CR"]
  idx <- match(te_data$SNU_ID, RF_ci$SNU_ID)
  k_pred_ci[idx] <- pp
  
  cat(sprintf("  [CI] Fold %d/%d 완료\n", i, nrow(folds_ci)))
}

valid_ci <- !is.na(k_pred_ci)

roc_ci <- pROC::roc(RF_ci$TRG_1[valid_ci], k_pred_ci[valid_ci],
                    levels = c("nonCR", "CR"), direction = "<", quiet = TRUE)
auc_ci <- as.numeric(pROC::auc(roc_ci))

cat(sprintf("\n[CI] Apparent AUC = %.3f (n = %d)\n", auc_ci, sum(valid_ci)))

pred_df_ci <- data.frame(SNU_ID  = RF_ci$SNU_ID,
                         TRG_1   = RF_ci$TRG_1,
                         Pred_CR = k_pred_ci)

# 같은 형태로 맞춰서 all_results에 추가
result_ci <- list(Input = "CI",
                  Model = "RF_CI",
                  FS    = NA,
                  AUC   = auc_ci,
                  ROC   = roc_ci,
                  Pred  = pred_df_ci)

if (!exists("all_results")) all_results <- list()
all_results[["CI"]] <- result_ci

# 최종 RDS 파일 저장
saveRDS(all_results, file = "RData/RDS/260907 RF_Global_Direction_all_results.rds")

# ---------- 3-2. ROC curve ----------
roc_list <- lapply(all_results, function(res) res$ROC)

# Legend에 input 이름, AUC, n 표시
names(roc_list) <- sapply(all_results, function(res) {
  sprintf("%s\n(AUC = %.3f; n = %d)", res$Input, res$AUC, length(res$Selected_features))
}) 

p_roc_all <- ggroc(roc_list, legacy.axes = TRUE, linewidth = 1.2, alpha = 0.75) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  scale_color_manual(values = unname(cols_input[names(all_results)]), name = "Input data") +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_pubr() +
  theme(aspect.ratio = 1,
        panel.border = element_rect(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.95, 0.05),
        legend.justification.inside = c(1, 0)); p_roc_all

ggsave("Figure/05 Prediction/01-5-ROC_curve-global-final.svg",
       plot = p_roc_all, width = 9, height = 6)

# ---------- 3-3. SHAP barplot ----------
TOP_N <- 10   # input당 표시할 top feature 개수

tile_vars   <- c("MDA", "MDG", "Consistency", "neg_log10_p")
tile_titles <- c("MDA (%)", "MDG (%)", "Consistency", "-log10(p)")

panel_list <- list()

for (nm in setdiff(names(all_results), "CI")) {
  
  res <- all_results[[nm]]
  
  anno <- res$Feature_summary %>%
    arrange(desc(MeanAbsSHAP)) %>%
    slice_head(n = TOP_N) %>%
    mutate(Feature = factor(Feature, levels = rev(Feature)),
           neg_log10_p = -log10(Pvalue))
  
  # SHAP bar
  p_bar <- ggplot(anno, aes(x = MeanAbsSHAP, y = Feature)) +
    geom_col(width = 0.7, fill = cols_input[[nm]]) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = "Mean |SHAP|", y = NULL, title = nm) +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
  
  # Importance tile (MDA / MDG / Consistency / Binomial-p)
  tile_list <- list()
  
  for (j in seq_along(tile_vars)) {
    v <- tile_vars[j]
    
    p_tile <- ggplot(anno, aes(x = 1, y = Feature, fill = .data[[v]])) +
      geom_tile() +
      labs(x = tile_titles[j], y = NULL) +
      theme_bw(base_size = 10) +
      theme(axis.text.y = element_blank(), axis.ticks = element_blank(),
            axis.text.x = element_blank(), panel.grid = element_blank(),
            legend.position = "bottom",
            plot.margin = ggplot2::margin(5.5, 0, 5.5, 0))
    
    if (v == "MDA")         p_tile <- p_tile + scale_fill_gradient(low = "white", high = "#2C7FB8", limits = c(0, 1))
    if (v == "MDG")         p_tile <- p_tile + scale_fill_gradient(low = "white", high = "#2F9EAA", limits = c(0, 1))
    if (v == "Consistency") p_tile <- p_tile + scale_fill_gradient(low = "white", high = "#4C956C",  limits = c(0.7, 1))
    if (v == "neg_log10_p") p_tile <- p_tile + scale_fill_gradient(low = "white", high = "#D6406B", limits = c(0, 15))
    
    tile_list[[j]] <- p_tile
  }
  
  # bar + tile들을 한 줄로 병합
  panel_list[[nm]] <- patchwork::wrap_plots(c(list(p_bar), tile_list), nrow = 1,
                                            widths = c(6, rep(0.6, length(tile_list))))
}

# 세로로 panel 쌓아서 출력
p_shap_facet_all <- patchwork::wrap_plots(panel_list, ncol = 1) &
  theme(legend.position = "none"); p_shap_facet_all

ggsave("Figure/05 Prediction/02-2-SHAP-best_model-top10.svg",
       plot = p_shap_facet_all, width = 7, height = 3.2 * length(panel_list))

ggsave("Figure/05 Prediction/02-3-SHAP-best_model-for_legend.svg",
       plot = panel_list$Strain, width = 30, height = 2)

# ========================================
save.image(file = "RData/05_RF modeling 260907.RData")
