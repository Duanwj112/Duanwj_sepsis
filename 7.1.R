library(Nebulosa)
library(Seurat)
library(patchwork)

sc_data <- readRDS('sc_Rds/sc_data.Rds')

options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
options("repos" = c(CRAN = "https://mirrors.westlake.edu.cn/CRAN/"))

sc_data <- RunTSNE(sc_data, dims = 1:30)

p4 <- plot_density(
  sc_data,
  c("G0S2", "GZMA", "ITM2A", "PAG1"),
  reduction = "tsne",
  joint     = TRUE,
  combine   = FALSE
)

combined_plot <- wrap_plots(p4)

ggsave(
  filename = "Fig 1/Fig 4E2_combined.pdf",
  plot     = combined_plot,
  device   = "pdf",
  width    = 12,
  height   = 8,
  dpi      = 300
)
