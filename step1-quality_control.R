# run libraries 
library(Seurat)
library(celldex)
library(SingleR)
library(DoubletFinder)
library(patchwork)

# read merged and clean data
merged_seurat <- readRDS("merged_seurat.rds")

# calculate mitPercent 
merged_seurat@meta.data["mitoPercent"] <- 
  PercentageFeatureSet(merged_seurat, pattern = "^MT-")

# Quality control (charts)
# Density plot
ggplot(merged_seurat@meta.data, aes(x = mitoPercent)) + geom_density() +
  geom_vline(xintercept = 15, color = "red", linetype = 2)

# evaluating the quality of data
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)

plot1 <- FeatureScatter(merged_seurat,
                        feature1 = "nCount_RNA", feature2 = "mitoPercent", 
                        raster=FALSE) + NoLegend()
plot2 <- FeatureScatter(merged_seurat,
                        feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        raster=FALSE) + NoLegend()
plot1 + plot2

# Filtering data based on the QC charts
merged_seurat_filtered <- subset(merged_seurat, subset =
                                   nCount_RNA > 100 & nCount_RNA < 40000 &
                                   nFeature_RNA > 200 & nFeature_RNA < 6000 &
                                   mitoPercent < 15)

# annotate cells automatically based on the ref
singleR_labels <- SingleR(test = as.SingleCellExperiment(merged_seurat_filtered),
                          ref = ref, labels = ref$label.main) 
merged_seurat_filtered@meta.data["singleR_labels"] <- singleR_labels$labels

saveRDS(merged_seurat_filtered, "merged_seurat_filtered")
