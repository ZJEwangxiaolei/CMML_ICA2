# The input file was from >cellbender remove-background --input "C:\Users\Admin\Desktop\GRCh38" --output "C:\Users\Admin\Desktop\out_GRCh38.h5"

data.file <- 'out_GRCh38.h5'
data.data <- Read10X_h5(filename = data.file, use.names = TRUE)

# create Seurat object
obj <- CreateSeuratObject(counts = data.data)
obj

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 1000)
obj <- ScaleData(obj, features = VariableFeatures(object = obj))
obj <- RunPCA(obj, npcs = 30)  
obj <- RunUMAP(obj, dims = 1:15)  

DimPlot(obj, reduction = "umap", label = TRUE)