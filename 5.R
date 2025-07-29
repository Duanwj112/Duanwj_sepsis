rm(list = ls())

library(tidyverse)
library(rms)
library(forestplot)
library(ggDCA)
library(rmda)
library(pROC)
library(nomogramFormula)

setwd("D:\\009.Spesis\\GEO\\GSE65682")

dat_expr_train <- read.table("GEOdata\\exp.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
dat_group_train <- read.table("GEOdata\\group.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
dat_group_train$sample <- rownames(dat_group_train)

key_genes <- read.csv("D:\\009.Spesis\\GEO\\hub_gene.csv")$x

dat_expr_train_key <- dat_expr_train[key_genes, ] %>%
  na.omit() %>% t() %>% as.data.frame()
dat_expr_train_key$group <- ifelse(dat_group_train$group == 'normal', 0, 1)

ddist <- datadist(dat_expr_train_key)
options(datadist = 'ddist')

fit_log_lrm <- lrm(as.formula(paste0('group ~ ', paste(key_genes, collapse = ' + '))), data = dat_expr_train_key, x = TRUE, y = TRUE)

nomo <- nomogram(
  fit_log_lrm,
  fun = plogis,
  fun.at = c(0.01, 0.1, 0.5, 0.9, 0.99),
  lp = TRUE,
  funlabel = 'Risk'
)

pdf(file = 'Fig 5/4A_diag_nomogram.pdf', width = 21, height = 10)
plot(nomo, col.grid = c('grey50', 'lightgrey'))
dev.off()

set.seed(2024)
dat_cal <- calibrate(fit_log_lrm, method = 'boot', B = 500)
dat_cal <- dat_cal[, 3:1] %>% as.data.frame()
colnames(dat_cal) <- c('bias_actual', 'apparent_actual', 'pre')
dat_cal <- data.frame(
  actual = c(dat_cal$bias_actual, dat_cal$apparent_actual),
  pre = rep(dat_cal$pre, 2),
  group = c(rep('Bias Corrected', nrow(dat_cal)), rep('Apparent', nrow(dat_cal)))
)

p <- ggplot(dat_cal, aes(pre, actual, color = group)) +
  geom_abline(slope = 1, intercept = 0, color = 'grey', lty = 2) +
  geom_line() +
  theme_bw() +
  theme(legend.title = element_blank()) +
  labs(x = 'Predicted Probability', y = 'Actual Probability') +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  coord_fixed() +
  scale_color_manual(values = c('#4DBBD5', '#E64B35'))

ggsave(file = 'Fig 5/4B_diag_calibration.pdf', p, height = 8, width = 8)

p <- dca(fit_log_lrm) %>%
  ggplot(aes(thresholds, NB, color = model)) +
  theme_bw() +
  scale_color_manual(values = c('#4DBBD5', '#E64B35', '#3C5488'), labels = c('Logistic Model', 'All', 'None')) +
  scale_linetype_manual(values = rep(1, 3)) +
  guides(lty = 'none') +
  labs(color = element_blank())

ggsave(file = 'Fig 5/4C_diag_DCA.pdf', p, height = 6, width = 7)

dca <- decision_curve(as.formula(paste0('group ~ ', paste(key_genes, collapse = ' + '))), dat_expr_train_key)

pdf(file = 'Fig 5/4D_CIC.pdf', height = 6, width = 7)
plot_clinical_impact(dca, col = c('#E64B35', '#4DBBD5'))
dev.off()

nomo_formula <- formula_rd(nomo)
nomo_point_train <- points_cal(nomo_formula$formula, as.matrix(dat_expr_train_key))

dat_roc_nomo_train <- data.frame(
  row.names = rownames(dat_group_train), 
  value = nomo_point_train, 
  group = factor(dat_group_train$group, levels = c('normal', 'spesis'))
)

roc_res_nomo_train <- roc(group ~ value, data = dat_roc_nomo_train, auc = TRUE, ci = TRUE)

dat_roc_plot_nomo_train <- data.frame(specificity = roc_res_nomo_train$specificities, 
                                      sensitivity = roc_res_nomo_train$sensitivities) %>%
  dplyr::arrange(desc(specificity), sensitivity)

p <- ggplot(dat_roc_plot_nomo_train, aes(1 - specificity, sensitivity)) + 
  geom_line(color = '#4DBBD5') +
  geom_abline(slope = 1, intercept = 0, color = 'grey', lty = 'dashed') +
  geom_area(fill = '#4DBBD5', alpha = 0.2) +
  annotate('text', label = paste0(
    'AUC = ',
    round(roc_res_nomo_train$auc, 3),
    ' (',
    round(roc_res_nomo_train$ci[1], 3),
    '-',
    round(roc_res_nomo_train$ci[3], 3),
    ')'
  ),
  x = 0.8, y = 0.1, color = '#4DBBD5') +
  theme_bw() +
  coord_fixed() +
  labs(x = '1-Specificity (FPR)', y = 'Sensitivity (TPR)')

ggsave(file = 'Fig 5/4E_ROC_nomo.pdf', p, height = 6, width = 6)

setwd("D:\\009.Spesis\\GEO\\GSE26378")

dat_expr_GSE26378 <- read.table("GEOdata\\exp.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
dat_group_GSE26378 <- read.table("GEOdata\\group.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
dat_group_GSE26378$sample <- rownames(dat_group_GSE26378)

nomo_point_GSE26378 <- points_cal(nomo_formula$formula, t(dat_expr_GSE26378))

dat_roc_nomo_GSE26378 <- data.frame(
  row.names = rownames(dat_group_GSE26378), 
  value = nomo_point_GSE26378, 
  group = factor(dat_group_GSE26378$group, levels = c('normal', 'spesis'))
)

roc_res_nomo_GSE26378 <- roc(group ~ value, data = dat_roc_nomo_GSE26378, auc = TRUE, ci = TRUE)

dat_roc_plot_nomo_GSE26378 <- data.frame(specificity = roc_res_nomo_GSE26378$specificities, 
                                         sensitivity = roc_res_nomo_GSE26378$sensitivities) %>%
  dplyr::arrange(desc(specificity), sensitivity)

p <- ggplot(dat_roc_plot_nomo_GSE26378, aes(1 - specificity, sensitivity)) + 
  geom_line(color = '#4DBBD5') + 
  geom_abline(slope = 1, intercept = 0, color = 'grey', lty = 'dashed') +
  geom_area(fill = '#4DBBD5', alpha = 0.2) + 
  annotate('text', label = paste0(
    'AUC = ',
    round(roc_res_nomo_GSE26378$auc, 3),
    ' (',
    round(roc_res_nomo_GSE26378$ci[1], 3),
    '-',
    round(roc_res_nomo_GSE26378$ci[3], 3),
    ')'
  ),
  x = 0.8, y = 0.1, color = '#4DBBD5') +
  theme_bw() + 
  coord_fixed() + 
  labs(x = '1-Specificity (FPR)', y = 'Sensitivity (TPR)')

ggsave(file = 'D:\\009.Spesis\\GEO\\GSE65682\\Fig 5/4F_ROC_GSE26378.pdf', p, height = 6, width = 6)
