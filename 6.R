rm(list = ls())

library(CIBERSORT)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggcorrplot)
library(psych)
library(data.table)

setwd("D:\\009.Spesis\\GEO\\GSE65682")

dat_risk <- readRDS('GEOdata/dat_risk.Rds')
rownames(dat_risk) <- dat_risk[[1]]
dat_risk <- dat_risk[, -1]  

dat_expr <- fread(input = "GEOdata/exp.txt", header = TRUE, data.table = FALSE)
rownames(dat_expr) <- dat_expr[, 1]
dat_expr <- dat_expr[, -1]

dat_res_cibersort <- cibersort(
  LM22,
  as.matrix(dat_expr),
  perm = 100,
  QN = TRUE
) %>% as.data.frame()

rownames(dat_res_cibersort) <- make.unique(substr(rownames(dat_res_cibersort), 1, 12))
dat_res_cibersort <- dat_res_cibersort[, 1:22]
dat_res_cibersort <- dat_res_cibersort[, colSums(dat_res_cibersort) > 0]

rownames(dat_res_cibersort) <- make.unique(substr(rownames(dat_res_cibersort), 1, 12))
dat_res_cibersort <- dat_res_cibersort[rownames(dat_risk), ]

long_dat_cibersort <- dat_res_cibersort
long_dat_cibersort$group <- dat_risk$riskgroup
long_dat_cibersort <- reshape2::melt(long_dat_cibersort)
colnames(long_dat_cibersort) <- c('group', 'cell', 'score')
long_dat_cibersort$group <- factor(long_dat_cibersort$group, levels = c('Low', 'High'))

p <- ggplot(long_dat_cibersort, aes(cell, score, fill = group)) +
  geom_boxplot(outlier.color = NA) +
  stat_compare_means(
    aes(group = group),
    label = 'p.signif',
    method = 'wilcox.test',
    paired = FALSE,
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c('***', '**', '*', 'ns')
    )
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = element_blank(), y = 'Immune Infiltration', fill = 'Group') +
  scale_fill_manual(values = c('#4DBBD5', '#E64B35'))

ggsave(file = 'CIBERSORT_boxplot_all1.pdf', p, width = 20, height = 5)

diagnostic_genes <- read.csv("D:\\009.Spesis\\GEO\\hub_gene.csv")
dat_expr_filtered <- dat_expr[rownames(dat_expr) %in% diagnostic_genes$x, ]

colnames(dat_expr) <- substr(colnames(dat_expr), 1, 12)
dat_expr <- dat_expr[, rownames(dat_risk)]
dat_expr <- t(dat_expr) %>% data.frame()

dat_expr_cor <- data.frame(
  row.names = rownames(dat_expr),
  GZMA = dat_expr$GZMA,
  PAG1 = dat_expr$PAG1,
  ITM2A = dat_expr$ITM2A,
  G0S2 = dat_expr$G0S2,
  RiskScore = dat_risk$riskscore
)
saveRDS(dat_expr_cor, file = 'dat_expr_cor.Rds')

cor_res_cibersort_cell_gene <- psych::corr.test(
  dat_res_cibersort,
  dat_expr_cor,
  method = 'spearman',
  adjust = 'BH'
)

mat_cor_cibersort_cell_gene <- cor_res_cibersort_cell_gene$r
mat_p_cibersort_cell_gene <- cor_res_cibersort_cell_gene$p

p <- ggcorrplot(
  mat_cor_cibersort_cell_gene,
  method = 'square',
  type = 'full',
  colors = c('#4DBBD5', 'white', '#E64B35'),
  p.mat = mat_p_cibersort_cell_gene,
  pch.col = NA,
  ggtheme = theme_void(),
  legend.title = 'Correlation',
  outline.color = NA
) +
  theme(axis.text.y = element_text(hjust = 1)) +
  geom_text(aes(label = ifelse(
    pvalue < 0.001, '***',
    ifelse(pvalue < 0.01, '**',
           ifelse(pvalue < 0.05, '*', NA))
  )))

ggsave(file = 'corr_heatmap_CIBERSORT_cell_gene.pdf', p, width = 10, height = 6)

# ssGSEA analysis
rm(list = ls())
library(GSVA)

dat_expr <- fread(input = "GEOdata/exp.txt", header = TRUE, data.table = FALSE)
rownames(dat_expr) <- dat_expr[, 1]
dat_expr <- dat_expr[, -1]

dat_risk <- readRDS('GEOdata/dat_risk.Rds')
rownames(dat_risk) <- dat_risk[[1]]
dat_risk <- dat_risk[, -1]  

dat_cell_marker <- read.csv('GEOdata/cellMarker.csv')
dat_cell_marker <- lapply(split(dat_cell_marker, dat_cell_marker$Celltype), function(x)
  unique(x$Metagene))

dat_res_ssgsea <- ssgseaParam(as.matrix(dat_expr), dat_cell_marker) %>% gsva() %>% as.data.frame()

colnames(dat_res_ssgsea) <- substr(colnames(dat_res_ssgsea), 1, 12)
dat_res_ssgsea <- dat_res_ssgsea[, rownames(dat_risk)]

long_dat_ssgsea <- dat_res_ssgsea %>% t() %>% as.data.frame()
long_dat_ssgsea$group <- dat_risk$riskgroup
long_dat_ssgsea <- reshape2::melt(long_dat_ssgsea)
colnames(long_dat_ssgsea) <- c('group', 'cell', 'score')
long_dat_ssgsea$group <- factor(long_dat_ssgsea$group, levels = c('Low', 'High'))

p <- ggplot(long_dat_ssgsea, aes(cell, score, fill = group)) +
  geom_boxplot(outlier.color = NA) +
  stat_compare_means(
    aes(group = group),
    label = 'p.signif',
    method = 'wilcox.test',
    paired = FALSE,
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c('***', '**', '*', 'ns')
    )
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = element_blank(), y = 'Immune Infiltration', fill = 'Group') +
  scale_fill_manual(values = c('#4DBBD5', '#E64B35'))

ggsave(file = 'ssGSEA_boxplot.pdf', p, width = 18, height = 5)


# ssGSEA correlation heatmap
dat_expr_cor <- readRDS('dat_expr_cor.Rds')

colnames(dat_res_ssgsea) <- substr(colnames(dat_res_ssgsea), 1, 12)
dat_res_ssgsea <- dat_res_ssgsea[, rownames(dat_risk)]

cor_res_ssgsea_cell_gene <- psych::corr.test(
  t(dat_res_ssgsea),
  dat_expr_cor,
  method = 'spearman',
  adjust = 'BH'
)

mat_cor_ssgsea_cell_gene <- cor_res_ssgsea_cell_gene$r
mat_p_ssgsea_cell_gene <- cor_res_ssgsea_cell_gene$p

p <- ggcorrplot(
  mat_cor_ssgsea_cell_gene,
  method = 'square',
  type = 'full',
  colors = c('#4DBBD5', 'white', '#E64B35'),
  p.mat = mat_p_ssgsea_cell_gene,
  pch.col = NA,
  ggtheme = theme_void(),
  legend.title = 'Correlation',
  outline.color = NA
) +
  theme(axis.text.y = element_text(hjust = 1)) +
  geom_text(aes(label = ifelse(
    pvalue < 0.001, '***',
    ifelse(pvalue < 0.01, '**',
           ifelse(pvalue < 0.05, '*', NA))
  )))

ggsave(file = 'corr_heatmap_ssGSEA_cell_gene.pdf', p, width = 10, height = 6)
