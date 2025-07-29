rm(list = ls())
library(tidyverse)
library(Seurat)
library(scDblFinder)
library(harmony)
library(SingleR)
library(HGNChelper)
library(AUCell)
library(monocle)
library(igraph)
library(CellChat)

# Read 10X data
sc_integer <- lapply(dir('GSE167363_RAW/', full.names = TRUE), function(x) Read10X(x, gene.column = 2))
sample_id <- list.files('GSE167363_RAW/') %>% strsplit('_') %>% sapply(`[`, 1) %>% unique()
for (i in seq_along(sample_id)) {
  sc_integer[[i]]@Dimnames[[2]] <- paste0(sc_integer[[i]]@Dimnames[[2]], '_', sample_id[i])
}

sc_integer <- lapply(sc_integer, function(x) {
  CreateSeuratObject(counts = x, min.cells = 3, min.features = 200)
})

for (i in seq_along(sample_id)) {
  sc_integer[[i]]$sample <- sample_id[i]
  sc_integer[[i]]$group <- ifelse(
    sample_id[i] %in% c('GSM5102900','GSM5102901'), 'Normal',
    ifelse(
      sample_id[i] %in% c('GSM5511351','GSM5102902','GSM5511352','GSM5102903'), 'Dead_sepsis',
      ifelse(
        sample_id[i] %in% c('GSM5102904','GSM5511353','GSM5511355','GSM5102905','GSM5511354','GSM5511356'), 'Alive_sepsis',
        NA
      )
    )
  )
}

sc_integer <- Reduce(merge, sc_integer) %>% JoinLayers()
saveRDS(sc_integer, file = 'sc_Rds/sc_integer.Rds')

sc_filtered <- sc_integer
sc_filtered$log10GenesPerUMI <- log10(sc_filtered$nFeature_RNA) / log10(sc_filtered$nCount_RNA)
sc_filtered$mitoRatio <- PercentageFeatureSet(sc_filtered, pattern = '^MT-') / 100
sc_filtered <- subset(sc_filtered,
                      nCount_RNA > 200 & nCount_RNA < 5000 &
                        nFeature_RNA > 200 & log10GenesPerUMI > 0.8 &
                        mitoRatio < 0.2)

mat_count <- LayerData(sc_filtered, assay = 'RNA', layer = 'counts')
sc_filtered <- CreateSeuratObject(mat_count[rowSums(mat_count > 0) >= 100, ], meta.data = sc_filtered@meta.data)

sc_filtered <- NormalizeData(sc_filtered)
sc_filtered <- ScaleData(sc_filtered, features = rownames(sc_filtered))

sce <- as.SingleCellExperiment(sc_filtered)
sce <- scDblFinder(sce, samples = 'sample')
sc_filtered <- as.Seurat(sce)
sc_filtered <- subset(sc_filtered, subset = scDblFinder.class == 'singlet')
saveRDS(sc_filtered, file = 'sc_Rds/sc_filtered.Rds')

VlnPlot(sc_filtered,
        features = c('nCount_RNA','nFeature_RNA','log10GenesPerUMI','mitoRatio'),
        pt.size = 0.01,
        cols = c('#4DBBD5','#E64B35','#00A087'),
        ncol = 4,
        group.by = 'group')

rm(list = ls())
sc_cc <- readRDS('sc_Rds/sc_filtered.Rds')
s_genes <- CaseMatch(cc.genes$s.genes, rownames(sc_cc))
g2m_genes <- CaseMatch(cc.genes$g2m.genes, rownames(sc_cc))

sc_cc <- CellCycleScoring(sc_cc,
                          s.features = s_genes,
                          g2m.features = g2m_genes,
                          set.ident = TRUE)
sc_cc$CC.Difference <- sc_cc$S.Score - sc_cc$G2M.Score

options(future.globals.maxSize = 16 * 1024^3)
sc_cc <- SCTransform(sc_cc,
                     vars.to.regress = c('S.Score','G2M.Score','CC.Difference','mitoRatio'),
                     verbose = TRUE)
sc_cc <- RunPCA(sc_cc, seed.use = 2024, features = rownames(sc_cc))
saveRDS(sc_cc, file = 'sc_Rds/sc_cc.Rds')
DimPlot(sc_cc, reduction = 'pca', group.by = 'Phase')

rm(list = ls())
sc_data <- readRDS('sc_Rds/sc_cc.Rds')
LabelPoints(VariableFeaturePlot(sc_data),
            points = head(VariableFeatures(sc_data),10),
            repel = TRUE)

Idents(sc_data) <- 'group'
FeatureScatter(sc_data,
               feature1 = 'nCount_RNA',
               feature2 = 'nFeature_RNA',
               pt.size = 0.1)

sc_data <- RunHarmony(sc_data,
                      group.by.vars = 'sample',
                      plot_convergence = TRUE)
sc_data <- RunPCA(sc_data, seed.use = 2024, features = VariableFeatures(sc_data))
ElbowPlot(sc_data, ndims = 50)

sc_data <- FindNeighbors(sc_data,
                         reduction = 'pca',
                         dims = 1:30,
                         verbose = TRUE)
sc_data <- FindClusters(sc_data,
                        verbose = TRUE,
                        resolution = c(0.6,0.8,1.0,1.2,1.4))
sc_data <- RunUMAP(sc_data,
                   seed.use = 2024,
                   reduction = 'harmony',
                   dims = 1:30)

DimPlot(sc_data, reduction = 'umap', group.by = 'SCT_snn_res.1', label = TRUE, repel = TRUE)
DimPlot(sc_data, reduction = 'umap', group.by = 'sample')
DimPlot(sc_data, reduction = 'umap', group.by = 'group')

hpca <- celldex::HumanPrimaryCellAtlasData()
testdata <- GetAssayData(sc_data, layer = 'data')

library(HGNChelper)
source('input/gene_sets_prepare.R')
source('input/sctype_score_.R')
db_ <- 'input/ScTypeDB_full.xlsx'
tissue <- 'Immune system'
gs_list <- gene_sets_prepare(db_, tissue)

es_max <- sctype_score(sc_data[['SCT']]@scale.data,
                       scaled = TRUE,
                       gs = gs_list$gs_positive,
                       gs2 = gs_list$gs_negative)

cL_res <- do.call(rbind, lapply(unique(sc_data$SCT_snn_res.0.8), function(cl) {
  es_max_cl <- sort(rowSums(es_max[, rownames(sc_data@meta.data)[sc_data$SCT_snn_res.0.8==cl]]),
                    decreasing = TRUE)
  head(data.frame(cluster = cl,
                  type    = names(es_max_cl),
                  scores  = es_max_cl,
                  ncells  = sum(sc_data$SCT_snn_res.0.8==cl)),
       10)
}))

sctype_score <- cL_res %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = scores)
sctype_score$type[as.numeric(as.character(sctype_score$scores)) < sctype_score$ncells/4] <- 'Unknown'

sc_data$cell_anno <- left_join(sc_data@meta.data,
                               sctype_score[,c('cluster','type')],
                               by = c('SCT_snn_res.1'='cluster'))$type

DimPlot(sc_data, reduction = 'umap', group.by = 'cell_anno', label = TRUE, repel = TRUE)

memory.limit(size = 132000)
ggsave('Fig_1A1.png', plot = last_plot(), width = 8, height = 6, dpi = 300)

ggplot(sc_data@meta.data, aes(cell_anno, fill = cell_anno)) +
  geom_bar(stat = 'count') +
  geom_text(data = data.frame(cell_anno = names(table(sc_data$cell_anno)),
                              n        = as.numeric(table(sc_data$cell_anno))),
            aes(cell_anno, n, label = n)) +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = NULL) +
  guides(fill = FALSE)

ggplot(sc_data@meta.data, aes(sample, fill = cell_anno)) +
  geom_bar(position = 'fill') +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = NULL)

ggplot(sc_data@meta.data, aes(group, fill = cell_anno)) +
  geom_bar(position = 'fill') +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = NULL)

saveRDS(sc_data, file = 'sc_Rds/sc_data.Rds')

rm(list = ls())
sc_data <- readRDS('sc_Rds/sc_data.Rds')
Idents(sc_data) <- 'cell_anno'
DefaultAssay(sc_data) <- 'RNA'
cells_to_keep <- rownames(sc_data@meta.data[sc_data@meta.data$group=='LUAD',])
sc_data <- subset(sc_data, cells = cells_to_keep)

sc_marker <- FindAllMarkers(sc_data, min.pct = 0.1, logfc.threshold = 0.5)
sc_marker <- sc_marker %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5)

marker1 <- sc_marker %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
marker5 <- sc_marker %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

FeaturePlot(sc_data, features = marker1$gene)
ggsave('Fig_1B1.png', plot = last_plot(), width = 8, height = 6, dpi = 300)

VlnPlot(sc_data, features = marker1$gene, pt.size = 0)
ggsave('Fig_1B2.png', plot = last_plot(), width = 8, height = 6, dpi = 300)

DotPlot(sc_data, features = unique(marker5$gene)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

saveRDS(sc_marker, file = 'sc_Rds/sc_marker.Rds')

rm(list = ls())
sc_data_AUCell <- readRDS('sc_Rds/sc_data.Rds')
gene_anoikis <- read.csv('E:/009.Phenotype_summary/Ferroptosis.csv') %>%
  filter(Category=='Protein Coding') %>%
  pull(Gene.Symbol) %>%
  unique()

Idents(sc_data_AUCell) <- 'cell_anno'
gene_anoikis <- list(Anoikis = gene_anoikis)

mat_expr <- t(FetchData(sc_data_AUCell, vars = rownames(sc_data_AUCell)))
cell_ranking <- AUCell_buildRankings(mat_expr)
cell_AUC <- AUCell_calcAUC(gene_anoikis, cell_ranking)
auc_score <- t(getAUC(cell_AUC)) %>% as.data.frame()
sc_data_AUCell$anoikis_auc <- auc_score$Anoikis

FeaturePlot(sc_data_AUCell, features = 'anoikis_auc')
VlnPlot(sc_data_AUCell, features = 'anoikis_auc')

sc_data_AUCell$anoikis_group <- ifelse(sc_data_AUCell$anoikis_auc >= median(sc_data_AUCell$anoikis_auc),
                                       'high_anoikis','low_anoikis')
DimPlot(sc_data_AUCell, reduction = 'umap', group.by = 'anoikis_group')
saveRDS(sc_data_AUCell, file = 'sc_Rds/sc_anoikis_score.Rds')

rm(list = ls())
sc_data_pseudo <- readRDS('sc_Rds/sc_data.Rds')
mat_expr <- as(t(FetchData(sc_data_pseudo, vars = rownames(sc_data_pseudo))), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = sc_data_pseudo@meta.data)
fd <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = rownames(sc_data_pseudo),
                                                  row.names = rownames(sc_data_pseudo)))
cds <- newCellDataSet(mat_expr, phenoData = pd, featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size()) %>%
  estimateSizeFactors() %>%
  estimateDispersions() %>%
  detectGenes(min_expr = 0.1)

disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table,
                                 mean_expression >= 0.1 &
                                   dispersion_empirical >= 0.5 * dispersion_fit)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)

diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = '~cell_anno')
ordering_genes <- rownames(subset(diff_test_res, qval < 0.01))
saveRDS(ordering_genes, file = 'sc_Rds/ordering_genes.Rds')

ordering_genes <- readRDS('sc_Rds/ordering_genes.Rds')
cds <- cds %>%
  setOrderingFilter(ordering_genes) %>%
  reduceDimension(max_components = 2, reduction_method = 'DDRTree') %>%
  orderCells()
saveRDS(cds, file = 'sc_Rds/cds.Rds')

cds <- readRDS('sc_Rds/cds.Rds')
plot_cell_trajectory(cds, color_by = 'Pseudotime')
plot_cell_trajectory(cds, color_by = 'State')
plot_cell_trajectory(cds, color_by = 'cell_anno')
plot_cell_trajectory(cds, color_by = 'cell_anno') + facet_wrap(~cell_anno)

gene_anoikis <- read.csv('input/GeneCards-TEX.csv') %>%
  filter(Category=='Protein Coding') %>%
  pull(Gene.Symbol) %>%
  intersect(rownames(cds))

plot_pseudotime_heatmap(cds[gene_anoikis[1:10],],
                        num_clusters = 3,
                        show_rownames = TRUE,
                        return_heatmap = TRUE,
                        cores = 1)

rm(list = ls())
sc_data_cellchat <- readRDS('sc_Rds/sc_data.Rds')
mat_input <- LayerData(sc_data_cellchat, assay = 'RNA', layer = 'counts')
meta <- data.frame(row.names = rownames(sc_data_cellchat@meta.data),
                   cell_anno = sc_data_cellchat$cell_anno)
cellchat <- createCellChat(object = mat_input, meta = meta, group.by = 'cell_anno') %>%
  addMeta(meta = meta) %>%
  setIdent(ident.use = 'cell_anno')
group_size <- as.numeric(table(cellchat@idents))
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
future::plan('default')
cellchat <- cellchat %>%
  identifyOverExpressedGenes() %>%
  identifyOverExpressedInteractions() %>%
  projectData(PPI.human) %>%
  computeCommunProb(raw.use = TRUE) %>%
  filterCommunication(min.cells = 10) %>%
  computeCommunProbPathway() %>%
  aggregateNet() %>%
  netAnalysis_computeCentrality(slot.name = 'netP')
saveRDS(cellchat, file = 'sc_Rds/cellchat.Rds')

netVisual_circle(cellchat@net$count, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE)
pdf(file = 'communication_interaction_plot1.pdf', height = 10, width = 12)
netVisual_circle(cellchat@net$count, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE)
dev.off()

mat0 <- matrix(0, nrow = nrow(cellchat@net$weight), ncol = ncol(cellchat@net$weight),
               dimnames = dimnames(cellchat@net$weight))
mat0[5,] <- cellchat@net$weight[5,]
pdf(file = 'single_cell_communication_plot.pdf', height = 10, width = 12)
netVisual_circle(mat0, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE)
dev.off()

pdf(file = 'cell_communication_heatmap.pdf', height = 10, width = 12)
netVisual_heatmap(cellchat)
dev.off()

pdf(file = 'cellular_centrality_heatmap.pdf', height = 10, width = 12)
netAnalysis_signalingRole_network(cellchat)
dev.off()

pdf(file = 'outgoing_signaling_heatmap.pdf', height = 10, width = 12)
netAnalysis_signalingRole_heatmap(cellchat, pattern = 'outgoing')
dev.off()

pdf(file = 'incoming_signaling_heatmap.pdf', height = 10, width = 12)
netAnalysis_signalingRole_heatmap(cellchat, pattern = 'incoming')
dev.off()

pdf(file = 'communication_bubble_plot.pdf', height = 20, width = 20)
netVisual_bubble(cellchat, remove.isolate = FALSE)
dev.off()
