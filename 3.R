rm(list = ls())
library(tidyverse)
library(caret)
library(DALEX)
library(pROC)
library(ggvenn)

setwd("D:\\009.Spesis\\GEO\\GSE65682")

dat_expr <- read.table("GEOdata\\exp.txt", header=TRUE, sep="\t", row.names=1, stringsAsFactors=FALSE, check.names=FALSE)
dat_group <- read.table("GEOdata\\group.txt", header=TRUE, sep="\t", row.names=1, stringsAsFactors=FALSE)
dat_group$sample <- rownames(dat_group)

ch_genes <- read.csv("D:\\009.Spesis\\GEO\\GSE65682\\GEOdata/common_genes.csv")
dat_expr <- dat_expr[, dat_group$sample]

dat_expr_ch <- dat_expr[ch_genes$Gene.Symbol, ] %>% t() %>% as.data.frame()
dat_expr_ch$group <- ifelse(dat_group$group == 'normal', '0', '1')

set.seed(2025)
idx <- createDataPartition(dat_expr_ch$group, p = 0.7, list = FALSE) %>% as.vector()
dat_train <- dat_expr_ch[idx, ]
dat_test <- dat_expr_ch[-idx, ]

tr_control <- trainControl(method = 'repeatedcv', repeats = 5)

model_svm <- train(group ~ ., data = dat_train, method = 'svmRadial', prob.model = TRUE, tuneLength = 1, trControl = tr_control)
model_rf  <- train(group ~ ., data = dat_train, method = 'rf', ntree = 100, prob.model = TRUE, tuneLength = 1, trControl = tr_control)
model_xgb <- train(group ~ ., data = dat_train, method = 'xgbTree', tuneLength = 1, trControl = tr_control)
model_glm <- train(group ~ ., data = dat_train, method = 'glm', family = 'binomial', tuneLength = 1, trControl = tr_control)

predicted_fun <- function(object, newdata) {
  predict(object, newdata, type = 'prob')[,2]
}

explainer_svm <- DALEX::explain(model_svm, label = 'SVM', data = dat_test, y = as.numeric(dat_test$group), predict_function = predicted_fun)
explainer_rf  <- DALEX::explain(model_rf,  label = 'RF',  data = dat_test, y = as.numeric(dat_test$group), predict_function = predicted_fun)
explainer_xgb <- DALEX::explain(model_xgb, label = 'XGB', data = dat_test, y = as.numeric(dat_test$group), predict_function = predicted_fun)
explainer_glm <- DALEX::explain(model_glm, label = 'GLM', data = dat_test, y = as.numeric(dat_test$group), predict_function = predicted_fun)

mp_svm <- model_performance(explainer_svm)
mp_rf  <- model_performance(explainer_rf)
mp_xgb <- model_performance(explainer_xgb)
mp_glm <- model_performance(explainer_glm)

pdf(file = '8b-residualtest.pdf', width = 6, height = 5)
plot(mp_svm, mp_rf, mp_xgb, mp_glm)
dev.off()

pdf(file = '8c-residual_boxplottest.pdf', width = 6, height = 5)
plot(mp_svm, mp_rf, mp_xgb, mp_glm, geom = 'boxplot')
dev.off()

dat_roc <- data.frame(
  group = dat_test$group,
  SVM = predict(model_svm, dat_test, type = 'prob')$`1`,
  RF  = predict(model_rf, dat_test, type = 'prob')$`1`,
  XGB = predict(model_xgb, dat_test, type = 'prob')$`1`,
  GLM = predict(model_glm, dat_test, type = 'prob')$`1`
)

roc_res_model <- roc(group ~ SVM + RF + XGB + GLM, data = dat_roc, auc = TRUE, ci = TRUE)

dat_roc_plot_model <- rbind(
  data.frame(specificity = roc_res_model$SVM$specificities, sensitivity = roc_res_model$SVM$sensitivities, model = 'SVM') %>% arrange(desc(specificity)),
  data.frame(specificity = roc_res_model$RF$specificities, sensitivity = roc_res_model$RF$sensitivities, model = 'RF') %>% arrange(desc(specificity)),
  data.frame(specificity = roc_res_model$XGB$specificities, sensitivity = roc_res_model$XGB$sensitivities, model = 'XGB') %>% arrange(desc(specificity)),
  data.frame(specificity = roc_res_model$GLM$specificities, sensitivity = roc_res_model$GLM$sensitivities, model = 'GLM') %>% arrange(desc(specificity))
)

dat_roc_plot_model$model <- factor(dat_roc_plot_model$model, levels = c('SVM', 'RF', 'XGB', 'GLM'))

p <- ggplot(dat_roc_plot_model, aes(1 - specificity, sensitivity, color = model)) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, color = 'grey', lty = 'dashed') +
  annotate('text', x = 0.75, y = 0.2, label = sprintf('AUC(SVM) = %.3f (%.3f-%.3f)', roc_res_model$SVM$auc, roc_res_model$SVM$ci[1], roc_res_model$SVM$ci[3]), color = '#4DBBD5') +
  annotate('text', x = 0.75, y = 0.15, label = sprintf('AUC(RF) = %.3f (%.3f-%.3f)', roc_res_model$RF$auc, roc_res_model$RF$ci[1], roc_res_model$RF$ci[3]), color = '#E64B35') +
  annotate('text', x = 0.75, y = 0.1, label = sprintf('AUC(XGB) = %.3f (%.3f-%.3f)', roc_res_model$XGB$auc, roc_res_model$XGB$ci[1], roc_res_model$XGB$ci[3]), color = '#3C5488') +
  annotate('text', x = 0.75, y = 0.05, label = sprintf('AUC(GLM) = %.3f (%.3f-%.3f)', roc_res_model$GLM$auc, roc_res_model$GLM$ci[1], roc_res_model$GLM$ci[3]), color = '#F39B7F') +
  theme_bw() +
  coord_fixed() +
  labs(x = '1-Specificity (FPR)', y = 'Sensitivity (TPR)') +
  scale_color_manual(values = c('#4DBBD5', '#E64B35', '#3C5488', '#F39B7F'))

ggsave(file = '8d-ROC_modeltest.pdf', p, height = 6, width = 7)

vi_svm <- variable_importance(explainer_svm, loss_function = loss_root_mean_square)
vi_rf  <- variable_importance(explainer_rf,  loss_function = loss_root_mean_square)
vi_xgb <- variable_importance(explainer_xgb, loss_function = loss_root_mean_square)
vi_glm <- variable_importance(explainer_glm, loss_function = loss_root_mean_square)

pdf(file = '8a-gene 6_importancetest1.pdf', width = 6, height = 12)
plot(vi_svm, vi_rf, vi_xgb, vi_glm, max_vars = 15)
dev.off()

top_gene_svm <- vi_svm %>% filter(!variable %in% c('_baseline_', '_full_model_', 'group')) %>% group_by(variable) %>% summarise(m = median(dropout_loss)) %>% arrange(desc(m)) %>% pull(variable, n = 15)
top_gene_rf  <- vi_rf  %>% filter(!variable %in% c('_baseline_', '_full_model_', 'group')) %>% group_by(variable) %>% summarise(m = median(dropout_loss)) %>% arrange(desc(m)) %>% pull(variable, n = 15)
top_gene_xgb <- vi_xgb %>% filter(!variable %in% c('_baseline_', '_full_model_', 'group')) %>% group_by(variable) %>% summarise(m = median(dropout_loss)) %>% arrange(desc(m)) %>% pull(variable, n = 15)
top_gene_glm <- vi_glm %>% filter(!variable %in% c('_baseline_', '_full_model_', 'group')) %>% group_by(variable) %>% summarise(m = median(dropout_loss)) %>% arrange(desc(m)) %>% pull(variable, n = 15)

p <- ggvenn(list(SVM = top_gene_svm, RF = top_gene_rf, XGB = top_gene_xgb, GLM = top_gene_glm), c('SVM', 'RF', 'XGB', 'GLM'),
            show_percentage = FALSE,
            fill_alpha = 0.5,
            stroke_color = NA,
            fill_color = c('#4DBBD5', '#E64B35', '#3C5488', '#F39B7F'))

ggsave(file = '8e-model_gene_venntest.pdf', p, width = 5, height = 5)

gene_model <- Reduce(intersect, list(top_gene_svm, top_gene_rf, top_gene_xgb, top_gene_glm))
write.csv(gene_model, file = "hub_gene.csv", row.names = FALSE)

dat_test1 <- dat_test[, c('PAG1', 'group')]
roc_res_PAG1 <- roc(group ~ PAG1, data = dat_test1, auc = TRUE, ci = TRUE)
dat_roc_plot_PAG1 <- data.frame(specificity = roc_res_PAG1$specificities, sensitivity = roc_res_PAG1$sensitivities) %>%
  arrange(desc(specificity), sensitivity)

p <- ggplot(dat_roc_plot_PAG1, aes(1 - specificity, sensitivity)) +
  geom_line(color = '#4DBBD5') +
  geom_abline(slope = 1, intercept = 0, color = 'grey', lty = 'dashed') +
  annotate('text', x = 0.8, y = 0.1,
           label = sprintf('AUC = %.3f (%.3f-%.3f)', roc_res_PAG1$auc, roc_res_PAG1$ci[1], roc_res_PAG1$ci[3]),
           color = '#4DBBD5') +
  theme_bw() +
  coord_fixed() +
  labs(x = '1-Specificity (FPR)', y = 'Sensitivity (TPR)')

ggsave(file = '8f-ROC_PAG1test.pdf', p, height = 6, width = 6)
