
library("Seurat")
library("arrow")
options(Seurat.object.assay.version = "v5")

matrix <- readRDS("E:/Mouse_M1/Mouse_M1_10xV3_Matrix.RDS")
meta <- arrow::read_feather("E:/Mouse_M1/Mouse_M1_10xV3_Metadata.feather")
obj <- CreateSeuratObject(counts = matrix, meta.data = meta)

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

# VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
