# This code shows the performance of DecontX and scCDC when preclustered
# Some additional exploration of doublet cells was done by DoubletFinder
setwd("/Users/wangxiaolei/Downloads")
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(celda)
library(dplyr)
pbmc <- Read10X("./PBMC4K/filtered_gene_bc_matrices/GRCh38")
##Seurat workflow
pbmc <- CreateSeuratObject(pbmc)
pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern =
                                              "^MT")
VlnPlot(pbmc, features = c("nFeature_RNA",
"nCount_RNA","percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nCount_RNA<20000 & nFeature_RNA
               <3000 &percent.mt < 10)
pbmc <- NormalizeData(pbmc,normalization.method =
                        "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst",
                             nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object =
                                                   pbmc))
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.7)
pbmc <- RunUMAP(pbmc, dims = 1:15)
# Visualize initial UMAP
DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend()

# Sweep parameters to find optimal pK
sweep.res.list_pbmc <- paramSweep(pbmc, PCs = 1:15, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT =
                                     FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
# Visualize pK selection plot
ggplot(bcmvn_pbmc, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() +
  theme_classic() +
  labs(title = "pK Optimization Plot", x = "pK Value", y =
         "BCmvn Metric")
# Select the pK value with the highest BCmetric
optimal_pk <-
  as.numeric(as.character(bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)]))
print(paste("Optimal pK selected:", optimal_pk))

# Get the number of cells
n_cells <- ncol(pbmc)
print(paste("Number of cells:", n_cells))
# Estimate expected doublet rate (e.g., 4% for 5k cells)
doublet_rate_estimate <- 0.04
nExp_poi <- round(doublet_rate_estimate * n_cells)
print(paste("Estimated number of doublets:", nExp_poi))

# Run DoubletFinder with the optimal pK and estimated doublet rate
pbmc <- doubletFinder(pbmc, PCs = 1:15, pN = 0.25, pK =
                        optimal_pk, nExp = nExp_poi, sct = FALSE)
# Check the metadata names
head(pbmc@meta.data)
df_classification_col <-
  colnames(pbmc@meta.data)[grepl("DF.classifications",
                                 colnames(pbmc@meta.data))]
print(paste("DoubletFinder classification column:",
            df_classification_col))

# Visualize predicted doublets on UMAP
DimPlot(pbmc, reduction = "umap", group.by =
          df_classification_col) +
  labs(title = paste("DoubletFinder Results (pK=", optimal_pk,
                     ", Rate=", doublet_rate_estimate, ")"))
# See the proportion of identified doublets
table(pbmc@meta.data[[df_classification_col]])

# Try a LOWER doublet rate
nExp_poi_low <- round(doublet_rate_estimate * 0.5 * n_cells)
pbmc <- doubletFinder(pbmc, PCs = 1:15, pN = 0.25, pK =
                        optimal_pk, nExp = nExp_poi_low, sct = FALSE)
df_classification_col_low <-
  colnames(pbmc@meta.data)[grepl("DF.classifications",
                                 colnames(pbmc@meta.data))]
# Try a HIGHER doublet rate
nExp_poi_high <- round(doublet_rate_estimate * 2.0 * n_cells)
nExp_poi_high <- min(nExp_poi_high, round(0.25 * n_cells))
pbmc <- doubletFinder(pbmc, PCs = 1:15, pN = 0.25, pK =
                        optimal_pk, nExp = nExp_poi_high, sct = FALSE)
df_classification_col_high <-
  colnames(pbmc@meta.data)[grepl("DF.classifications",
                                 colnames(pbmc@meta.data))]

#Part2

# Get raw counts and UMAP coordinates
counts_matrix <- GetAssayData(pbmc,slot = "counts")
umap_coords <- as.vector(pbmc@active.ident)
# Check dimensions
dim(counts_matrix)
dim(umap_coords)

# Run DecontX
set.seed(12345)
decontx_results <- decontX(x = counts_matrix, z = umap_coords)
# Add contamination estimates to Seurat metadata
pbmc_DBX <- pbmc
pbmc_DBX$DecontX_Contamination <- decontx_results$contamination
head(pbmc_DBX@meta.data)

# Visualize contamination levels on UMAP
FeaturePlot(pbmc_DBX, features = "DecontX_Contamination", reduction
            = "umap") +
  scale_colour_gradientn(colours =
                           rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu"))) +
  labs(title = "Estimated Ambient RNA Contamination (DecontX)")
# Plot distribution of contamination scores
ggplot(pbmc_DBX@meta.data, aes(x = DecontX_Contamination)) +
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


FeaturePlot(pbmc, features = c("MALAT1","SERF2"))
FeaturePlot(seurat_obj_corrected, features = c("MALAT1","SERF2"))

GCGs_DBX <- ContaminationDetection(pbmc_DBX)
GCGs_corrected <- ContaminationDetection(seurat_obj_corrected)

dbx_top100 <- head(GCGs_DBX$mean_distance, 100)
cor_top100 <- head(GCGs_corrected$mean_distance, 100)
pbmc_top100 <- head(GCGs$mean_distance, 100)

df_boxplot <- data.frame(
  mean_distance = c(pbmc_top100,dbx_top100, cor_top100),
  group = rep(c("Before Correction","DecontX pre-clustered", "scCDC"), each=100)
)


ggplot(df_boxplot, aes(x = group, y = mean_distance, fill = group)) +
  geom_boxplot() +
  labs(
    title = "Top 100 mean_distance Boxplot",
    x = "Group",             
    y = "mean_distance"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",          
    axis.text.x = element_text(size = 14) 
  )

