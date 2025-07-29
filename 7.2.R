rm(list = ls())

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(enrichplot)
library(ggplot2)
library(ggridges)
library(limma)
library(edgeR)
library(RColorBrewer)
library(data.table)

setwd("D:/009.Spesis/GEO/GSE65682")

risk <- read.csv("4gene_model_risk_group.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(risk) <- risk$sample
risk$sample <- NULL

dat_expr <- fread("GEOdata/exp.txt", header = TRUE, data.table = FALSE)
rownames(dat_expr) <- dat_expr[, 1]
dat_expr <- dat_expr[, -1]

dat_expr <- dat_expr[, !duplicated(substr(colnames(dat_expr), 1, 12))]
colnames(dat_expr) <- substr(colnames(dat_expr), 1, 12)
dat_expr <- dat_expr[, rownames(risk)]

expr_log <- log2(dat_expr + 1)

design <- model.matrix(~ 0 + factor(risk$riskgroup, levels = c("Low", "High")))
colnames(design) <- c("Low", "High")

fit <- lmFit(expr_log, design)

contrast.matrix <- makeContrasts(High - Low, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

res_limma <- topTable(fit2, number = Inf, sort.by = "none")

dat_res_diff <- res_limma %>%
  rownames_to_column("gene") %>%
  as.data.frame() %>%
  na.omit()

write.csv(dat_res_diff, file = "Fig6_limma_diffexpr.csv", row.names = FALSE, quote = FALSE)

logFC_all_gene <- set_names(dat_res_diff$logFC, dat_res_diff$gene)
logFC_all_gene <- sort(logFC_all_gene, decreasing = TRUE)

geneset_c2 <- read.gmt("GEOdata/c2.all.v2024.1.Hs.symbols.gmt.txt")

set.seed(2025)
res_gsea <- GSEA(
  geneList      = logFC_all_gene,
  TERM2GENE     = geneset_c2,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  seed          = 2024
)

dat_res_gsea <- as.data.frame(res_gsea)
write.csv(dat_res_gsea, file = "GSEA_result.csv", row.names = FALSE, quote = FALSE)

p_up <- gseaplot2(
  res_gsea,
  geneSetID    = c(
    "ZHOU_INFLAMMATORY_RESPONSE_LPS_UP",
    "HAMAI_APOPTOSIS_VIA_TRAIL_UP",
    "REACTOME_NEUTROPHIL_DEGRANULATION",
    "REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING"
  ),
  pvalue_table = FALSE,
  color        = c("#4DBBD5", "#E64B35", "#3C5488", "#692F7C")
)
ggsave(file = "GSEA_plot_UP.pdf", p_up, width = 8, height = 8)

p_down <- gseaplot2(
  res_gsea,
  geneSetID    = c(
    "WP_TCELL_RECEPTOR_SIGNALING",
    "BIOCARTA_TCR_PATHWAY",
    "PID_CD8_TCR_PATHWAY"
  ),
  pvalue_table = FALSE,
  color        = c("#4DBBD5", "#E64B35", "#3C5488")
)
ggsave(file = "GSEA_plot_Down.pdf", p_down, width = 8, height = 8)
