# run libraries
library(Seurat) 
library(ggplot2)
library(patchwork)
library(clustree)
library(harmony)

# load data
merged_singlets <- readRDS("merged_singlets.rds")

# ========================== seurat workflow ===================================
merged_singlets <- NormalizeData(merged_singlets,
                                 normalization.method = "LogNormalize", scale.factor = 10000)
merged_singlets <- FindVariableFeatures(merged_singlets)
merged_singlets <- ScaleData(merged_singlets)
merged_singlets <- RunPCA(merged_singlets)
ElbowPlot(merged_singlets)
merged_singlets <- RunUMAP(merged_singlets, dims = 1:20, reduction = "pca")

# DimPlots before integration
Dataset <- DimPlot(merged_singlets, reduction = "umap", group.by = "Dataset", raster = FALSE) +
  ggtitle("Dataset") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0, 'cm'),
        legend.key.height = unit(0, 'cm'),
        legend.key.width= unit(0, 'cm'),
        legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3)))

Sample <- DimPlot(merged_singlets, reduction = "umap", group.by = "Sample", raster = FALSE) +
  ggtitle("Sample") + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0, 'cm'),
        legend.key.height = unit(0, 'cm'),
        legend.key.width= unit(0, 'cm'),
        legend.text = element_text(size = 5)) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 3)), )

plots <- list(cnd, Dataset, Sample)
wrap_plots(plots, ncol = 2, nrow = 2) &
  plot_annotation(title = "HL before Harmony")

# ========================== harmony integration ===============================
harmony.all.sample <- merged_singlets %>% 
  RunHarmony(group.by.vars = "Sample",
             max.iter.harmony = 20,
             plot_convergence = FALSE) 

harmony.all.sample <- harmony.all.sample %>% 
  RunUMAP(reduction = "harmony", dims = 1:20, reduction.name = "Harmony_UMAP") %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20)

# Clustering with louvain (algorithm 1). each time you run clustering, the data is stored in meta.data columns
for (res in c(0, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 1)) {
  print(paste0("resolution: ", res))
  harmony.all.sample <- FindClusters(harmony.all.sample, resolution = res, algorithm = 1)
}

# clustree package helps to visualize how cells are distributed between clusters depending on resolution
# to demonstrate stability: node_colour = "sc3_stability"
clustree::clustree(harmony.all.sample@meta.data, prefix = "RNA_snn_res.") +
  theme(legend.key.size = unit(0.3, 'cm'))

# change data in seurat_cluster to RNA_snn_res.0.25 (in default seurat_cluster resolution is equal to 0.8)
harmony.all.sample@meta.data$seurat_clusters <- harmony.all.sample@meta.data$RNA_snn_res.0.25
harmony.all.sample@meta.data <- harmony.all.sample@meta.data[, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 11)]
harmony.all.sample@meta.data$seurat_clusters <- as.numeric(harmony.all.sample@meta.data$seurat_clusters)

saveRDS(harmony.all.sample, "harmony.all.sample.0.25.rds")
harmony.all.sample.0.25 <- readRDS("harmony.all.sample.0.25.rds")

# DimPlot after harmony integration
Dataset <- DimPlot(harmony.all.sample.0.25, reduction = "Harmony_UMAP", group.by = "Dataset", raster = FALSE) +
  ggtitle("Dataset") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0, 'cm'),
        legend.key.height = unit(0, 'cm'),
        legend.key.width= unit(0, 'cm'),
        legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3)))

Sample <- DimPlot(harmony.all.sample.0.25, reduction = "Harmony_UMAP", group.by = "Sample", raster = FALSE) +  
  ggtitle("Sample") + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(0, 'cm'),
        legend.key.height = unit(0, 'cm'),
        legend.key.width= unit(0, 'cm'),
        legend.text = element_text(size = 5)) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 3)), )

plots <- list(Dataset, seurat_clusters, Sample, cellType)
wrap_plots(plots, ncol = 2, nrow = 2) &
  plot_annotation(title = "HL after harmony.Sample res: 0.25")

# ========================== find markers ======================================
# FindConservedMarkers (total Endo cells vs. all other cells) 
harmony.all.sample.0.25 <- readRDS("rds/harmony.all.sample.0.25.rds") 
Idents(harmony.all.sample.0.25) <- harmony.all.sample.0.25@meta.data$seurat_clusters

harmony.all.sample.0.25 <- RenameIdents(object = harmony.all.sample.0.25, "16" = "endo cells", "7" = "endo cells")

HL.FindConservedMarkers.AllEndoCells <- FindConservedMarkers(harmony.all.sample.0.25,
                                                             ident.1 = "endo cells",
                                                             grouping.var = "Condition",
                                                             only.pos = TRUE)

saveRDS(HL.FindConservedMarkers.AllEndoCells, "rds/HL.FindConservedMarkers.AllEndoCells.rds")

# FeaturePlot LRRC32 expression in lung cells
# all lung cells
FeaturePlot(harmony.all.sample.0.25,
            features = "LRRC32",
            reduction = "Harmony_UMAP",
            # split.by = "Condition",
            keep.scale = "feature",
            min.cutoff = "q10",
            raster = FALSE,
            cols = c("grey", "blue")
)

# VlnPlot LRRC32 in all lung clusters 
plot <- VlnPlot(harmony.all.sample.0.25,
                features = c("LRRC32"), 
                split.by = "Condition",
                split.plot = TRUE,
                group.by = "seurat_clusters", 
                # pt.size = 0,
                raster = FALSE, 
                assay = "RNA", 
                slot = "data") &
  ggtitle("LRRC32 Expression in Lung Clusters") 
wrap_plots(plots = plot, ncol = 1) 

# FindMarkers (Endo.cells_IPF and Endo.cells_CTRL) --- differential gene expression
harmony.all.sample.0.25 <- readRDS("rds/harmony.all.sample.0.25.rds") 

harmony.all.sample.0.25$clust.cnd <- paste0(
  harmony.all.sample.0.25@meta.data$seurat_clusters, "_", harmony.all.sample.0.25@meta.data$Condition)
Idents(harmony.all.sample.0.25) <- harmony.all.sample.0.25$clust.cnd

harmony.all.sample.0.25 <- RenameIdents(object = harmony.all.sample.0.25, "16_CTRL" = "endo_CTRL", 
                                        "16_IPF" = "endo_IPF", "7_CTRL" = "endo_CTRL", "7_IPF" = "endo_IPF")

Markers <- FindMarkers(harmony.all.sample.0.25,
                       logfc.threshold = 0.25,
                       ident.1 = "endo_IPF",
                       ident.2 = "endo_CTRL",
                       test.use = "MAST",
                       slot = "data")

saveRDS(Markers, "test/HL.FindMarkers_AllEndoIPF.vs.AllEndoCTRL.rds")

# featurePlot LRRC32 expression in total endo IPF vs. CTRL
# CTRL vs. IPF cells
plot <- FeaturePlot(harmony.all.sample.0.25,
                    features = "LRRC32",
                    reduction = "Harmony_UMAP",
                    split.by = "Condition",
                    keep.scale = "feature",
                    min.cutoff = "q10",
                    raster = FALSE,
                    cols = c("grey", "blue")
) 
print(wrap_plots(plot, nrow = 2))
