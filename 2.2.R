rm(list = ls())
library(tidyverse)
library(broom)
library(forestplot)
library(broom)

setwd("D:\\009.Spesis\\GEO\\GSE65682")

filter_res <- read.csv("common_genes_40.csv",
                       header = TRUE,
                       stringsAsFactors = FALSE)

exp <- read.table("GEOdata/exp.txt",
                  header = TRUE,
                  sep = "\t",
                  row.names = 1,
                  stringsAsFactors = FALSE)

group_df <- read.table("GEOdata/group.txt",
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE)

genes <- filter_res$x
dat <- as.data.frame(t(exp[genes, ])) %>%
  rownames_to_column("sample") %>%
  left_join(group_df, by="sample") %>%
  mutate(group=factor(group,levels=c("normal","spesis")))

genes <- filter_res$x

uni_res <- map_dfr(genes, function(gene) {
  fml <- as.formula(paste0("group ~ ", gene))
  fit <- glm(fml, data = dat, family = binomial)
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == gene) %>%
    transmute(
      gene    = gene,
      OR      = estimate,
      CI_low  = conf.low,
      CI_high = conf.high,
      p.value = p.value
    )
})

tabletext_uni <- cbind(
  c("Gene", uni_res$gene),
  c("OR (95% CI)",
    paste0(
      round(uni_res$OR, 2),
      " (",
      round(uni_res$CI_low, 2),
      "-", 
      round(uni_res$CI_high, 2),
      ")"
    )
  ),
  c("p-value",
    formatC(uni_res$p.value, format = "f", digits = 3)
  )
)

write.csv(
  uni_res,
  file     = "uni_res_GSE65682.csv",
  row.names = FALSE,
  quote     = TRUE
)

library(forestplot)
library(grid)

setwd("D:/009.Spesis/GEO/GSE65682")
dat_cox <- read.csv("uni_res_GSE65682_40.csv", 
                    header = TRUE, 
                    stringsAsFactors = FALSE)
colnames(dat_cox)[1] <- "gene"
dat_cox$OR      <- as.numeric(dat_cox$OR)
dat_cox$OR1     <- as.numeric(dat_cox$OR1)
dat_cox$OR2     <- as.numeric(dat_cox$OR2)
dat_cox$p.value <- as.numeric(dat_cox$p.value)

labeltext <- cbind(
  Gene  = dat_cox$gene,
  OR_CI = sprintf("%.2f (%.2f–%.2f)", 
                  dat_cox$OR, dat_cox$OR1, dat_cox$OR2),
  Pval  = formatC(dat_cox$p.value, 
                  format = "e", digits = 2)
)
labeltext <- rbind(
  c("Gene", "OR (95% CI)", "p-value"),
  labeltext
)

or_vals  <- dat_cox$OR
ci_lower <- dat_cox$OR1
ci_upper <- dat_cox$OR2

lims <- quantile(c(ci_lower, ci_upper), 
                 probs = c(0.05, 0.95), 
                 na.rm = TRUE)
xlim <- lims
clip <- lims

log_ticks  <- pretty(log10(c(ci_lower, ci_upper)), n = 5)
xticks_all <- 10^log_ticks
xticks <- xticks_all[xticks_all >= clip[1] & xticks_all <= clip[2]]
if (length(xticks) < 2) xticks <- pretty(xlim, n = 5)

txtgp <- fpTxtGp(
  label = gpar(fontsize = 10),
  ticks = gpar(fontsize =  9),
  xlab  = gpar(fontsize = 11)
)

n_row <- nrow(labeltext)
hrzl <- c(
  list(
    "1" = gpar(lwd = 1.5, col = "#000000"),
    "2" = gpar(lwd = 1,   col = "#CCCCCC")
  ),
  setNames(
    list(gpar(lwd = 1.5, col = "#000000")),
    as.character(n_row + 1)
  )
)

pdf("plot.pdf", width = 16, height = 8, paper = "special", onefile = FALSE)
forestplot(
  labeltext   = labeltext,
  mean        = c(NA, or_vals),
  lower       = c(NA, ci_lower),
  upper       = c(NA, ci_upper),
  clip        = clip,
  ci.vertices = TRUE,
  fn.ci_norm  = fpDrawNormalCI,
  xlog        = TRUE,
  xlab        = "Odds Ratio",
  xticks      = xticks,
  zero        = 1,
  colwidths   = c(0.15, 0.15, 0.15, 0.55),
  graphwidth  = unit(0.55, "npc"),
  colgap      = unit(1, "mm"),
  boxsize     = 0.25,
  txt_gp      = txtgp,
  hrzl_lines  = hrzl,
  lineheight  = unit(7, "mm"),
  graph.pos   = 4,
  align       = c("l", "c", "c")
)
dev.off()
