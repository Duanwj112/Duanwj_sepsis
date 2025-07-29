rm(list = ls())

library(tidyverse)
library(broom)
library(forestplot)

setwd("D:\\009.Spesis\\GEO\\GSE65682")

filter_res <- read.csv("D:\\009.Spesis\\GEO\\hub_gene.csv", header = TRUE, stringsAsFactors = FALSE)

exp <- read.table("GEOdata/exp.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

group_df <- read.table("GEOdata/group.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

genes <- filter_res$x
dat <- as.data.frame(t(exp[genes, ])) %>%
  rownames_to_column("sample") %>%
  left_join(group_df, by = "sample") %>%
  mutate(group = factor(group, levels = c("normal", "spesis")))

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

multi_fit <- glm(
  as.formula(paste0("group ~ ", paste(genes, collapse = " + "))),
  data   = dat,
  family = binomial
)

multi_res <- tidy(multi_fit, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term %in% genes) %>%
  transmute(
    gene    = term,
    OR      = estimate,
    CI_low  = conf.low,
    CI_high = conf.high,
    p.value = p.value
  )

tabletext_multi <- cbind(
  c("Gene", multi_res$gene),
  c("OR (95% CI)",
    paste0(
      round(multi_res$OR, 2),
      " (",
      round(multi_res$CI_low, 2),
      "-",
      round(multi_res$CI_high, 2),
      ")"
    )
  ),
  c("p-value",
    formatC(multi_res$p.value, format = "f", digits = 3)
  )
)

write.csv(uni_res, file = "uni_res_GSE65682.csv", row.names = FALSE, quote = TRUE)
write.csv(multi_res, file = "multi_res_GSE65682.csv", row.names = FALSE, quote = TRUE)

multi_coefs <- broom::tidy(multi_fit)
coef_df <- multi_coefs %>%
  select(term, estimate) %>%
  rename(gene = term, coef = estimate)

saveRDS(coef_df, file = "coef_logistic.Rds")

coef_vec   <- setNames(coef_df$coef, coef_df$gene)
expr_mat   <- as.matrix(dat[genes])
intercept  <- coef_vec["(Intercept)"]
linear_score <- as.numeric(intercept + expr_mat %*% coef_vec[genes])
names(linear_score) <- rownames(dat)

pred_prob <- 1 / (1 + exp(-linear_score))
names(pred_prob) <- rownames(dat)
dat$riskscore <- pred_prob

cutoff_med <- median(pred_prob, na.rm = TRUE)
dat$riskgroup <- ifelse(dat$riskscore >= cutoff_med, "High", "Low")
dat$riskgroup <- factor(dat$riskgroup, levels = c("Low", "High"))

output_df <- dat %>%
  select(sample, riskscore, riskgroup)

saveRDS(output_df, file = "dat_risk.Rds")
