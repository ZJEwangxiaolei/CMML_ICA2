# This code shows the performance of DecontX(default) and scCDC(preclustered by DecontX)  
setwd("/Users/wangxiaolei/Downloads")
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(celda)
library(dplyr)
library(scCDC)

# Load data
pbmc <- Read10X("./PBMC4K/filtered_gene_bc_matrices/GRCh38")

# Seurat workflow
pbmc <- CreateSeuratObject(pbmc)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nCount_RNA < 20000 & nFeature_RNA < 3000 & percent.mt < 10)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc)
pbmc <- RunUMAP(pbmc, dims = 1:15)

# Run DecontX
set.seed(12345)
decontx_results <- decontX(x = counts_matrix, z = umap_coords)
# Add contamination estimates to Seurat metadata
pbmc$DecontX_Contamination <- decontx_results$contamination
head(pbmc@meta.data)

# Visualize contamination levels on UMAP
FeaturePlot(pbmc_DBX, features = "DecontX_Contamination", reduction
            = "umap") +
  scale_colour_gradientn(colours =
                           rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu"))) +
  labs(title = "Estimated Ambient RNA Contamination (DecontX)")
# Plot distribution of contamination scores
ggplot(pbmc@meta.data, aes(x = DecontX_Contamination)) +
  geom_histogram(bins = 50) +
  theme_classic() +
  labs(title = "Distribution of DecontX Contamination Scores", x
       = "Contamination Fraction")

# Example: Access decontaminated counts
decontaminated_counts <- decontx_results$decontXcounts
dim(decontaminated_counts)

# Run scCDC
library(scCDC)
set.seed(12345)

# Step 1: Detect global contamination-causing genes (GCGs) 
GCGs <- ContaminationDetection(pbmc)

# Step 2: Quantify contamination ratios
contamination_ratio = ContaminationQuantification(pbmc,rownames(GCGs))
str(contamination_ratio)

# Step 3: Correct contamination
seurat_obj_corrected = ContaminationCorrection(pbmc,rownames(GCGs))

# Set corrected assay as default
DefaultAssay(seuratobj_corrected) = "Corrected"

# Extract corrected count matrix
corrected_count_matrix = data.frame(seurat_obj_corrected@assays[["Corrected"]]@layers$counts)

dim(corrected_count_matrix)

# Data normalization
seurat_obj_corrected <- NormalizeData(seurat_obj_corrected, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj_corrected <- FindVariableFeatures(seurat_obj_corrected, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_obj_corrected)
seurat_obj_corrected <- ScaleData(seurat_obj_corrected, features = all.genes)
seurat_obj_corrected <- RunPCA(seurat_obj_corrected, features = VariableFeatures(object = seurat_obj_corrected))
ElbowPlot(seurat_obj_corrected)
seurat_obj_corrected <- FindNeighbors(seurat_obj_corrected, dims = 1:15)
seurat_obj_corrected <- FindClusters(seurat_obj_corrected, resolution = 0.2, verbose = FALSE)
seurat_obj_corrected <- RunUMAP(seurat_obj_corrected, dims = 1:15)

gcg_features <- head(rownames(GCGs), 2)
gcg_features

