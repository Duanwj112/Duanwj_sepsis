rm(list = ls())
library(tidyverse)
library(limma)
library(ggplot2)

setwd("E:\\009.Spesis\\GEO\\GSE65682")

dat_expr <- read.table("GEOdata/exp.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
dat_group <- read.table("GEOdata/group.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
dat_group$sample <- rownames(dat_group)
dat_group <- dat_group[, c("sample", setdiff(colnames(dat_group), "sample"))]

design <- model.matrix( ~ 0 + factor(dat_group$group))
colnames(design) <- levels(factor(dat_group$group))
rownames(design) <- dat_group$sample

fit <- lmFit(dat_expr, design)
contrasts <- c('spesis-normal')
cont_matrix <- makeContrasts(contrasts = contrasts, levels = design)
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2)

dat_res_diff <- topTable(fit2, coef = contrasts, n = Inf) %>% na.omit()
dat_res_diff$P.Value[dat_res_diff$P.Value == 0] <- min(dat_res_diff$P.Value[dat_res_diff$P.Value != 0])
dat_res_diff$adj.P.Val[dat_res_diff$adj.P.Val == 0] <- min(dat_res_diff$adj.P.Val[dat_res_diff$adj.P.Val != 0])

logFC_cutoff <- 1
dat_res_diff$change <- ifelse(
  dplyr::between(dat_res_diff$logFC, -logFC_cutoff, logFC_cutoff) |
    dat_res_diff$adj.P.Val >= 0.05,
  'Not',
  ifelse(dat_res_diff$logFC >= logFC_cutoff, 'Up', 'Down')
)
gene_DEGs <- dat_res_diff %>%
  dplyr::filter(change != 'Not') %>%
  dplyr::arrange(desc(logFC)) %>%
  rownames()

dat_vocano <- dat_res_diff[, c('logFC', 'adj.P.Val', 'change')]
dat_vocano$logP <- -log10(dat_vocano$adj.P.Val)
dat_vocano$change <- factor(dat_vocano$change, levels = c('Down', 'Not', 'Up'))

p <- ggplot(dat_vocano, aes(logFC, logP, color = change)) +
  geom_point(alpha = 0.6) +
  theme_bw() +
  labs(x = 'LogFC', y = '-Log10(adj. p)', color = 'Change') +
  geom_hline(yintercept = -log10(0.05), lty = 2) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = 2) +
  scale_x_continuous(limits = c(-11, 11)) +
  scale_color_manual(values = c('#4DBBD5', 'grey', '#E64B35'))

ggsave(file = 'GSE65682_limma.pdf', p, width = 7, height = 6)
write.csv(dat_res_diff, file = "GSE65682_limma_result.csv", row.names = TRUE)
write.csv(data.frame(gene = gene_DEGs), file = "GSE65682_limma_filter.csv", row.names = FALSE)

logFC_cutoff <- 0.58
dat_res_diff$change <- ifelse(
  dplyr::between(dat_res_diff$logFC, -logFC_cutoff, logFC_cutoff) |
    dat_res_diff$adj.P.Val >= 0.05,
  'Not',
  ifelse(dat_res_diff$logFC >= logFC_cutoff, 'Up', 'Down')
)
gene_DEGs <- dat_res_diff %>%
  dplyr::filter(change != 'Not') %>%
  dplyr::arrange(desc(logFC)) %>%
  rownames()

dat_vocano <- dat_res_diff[, c('logFC', 'adj.P.Val', 'change')]
dat_vocano$logP <- -log10(dat_vocano$adj.P.Val)
dat_vocano$change <- factor(dat_vocano$change, levels = c('Down', 'Not', 'Up'))

p <- ggplot(dat_vocano, aes(logFC, logP, color = change)) +
  geom_point(alpha = 0.6) +
  theme_bw() +
  labs(x = 'LogFC', y = '-Log10(adj. p)', color = 'Change') +
  geom_hline(yintercept = -log10(0.05), lty = 2) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = 2) +
  scale_x_continuous(limits = c(-11, 11)) +
  scale_color_manual(values = c('#4DBBD5', 'grey', '#E64B35'))

ggsave(file = 'GSE65682_limma_log2FC 0.58.pdf', p, width = 7, height = 6)
write.csv(dat_res_diff, file = "GSE65682_limma_result_log2FC 0.58.csv", row.names = TRUE)
write.csv(data.frame(gene = gene_DEGs), file = "GSE65682_limma_filter_log2FC 0.58.csv", row.names = FALSE)
