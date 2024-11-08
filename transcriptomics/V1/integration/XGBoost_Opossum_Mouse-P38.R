# XGBoost classification

library(Seurat)
library(reticulate)
library(scrubletR)
library(ggplot2)
library(cowplot)
library(dplyr)
library(SeuratDisk)

# Source scripts with functions
source("C:/Ryan/GitHub/trachtenberg-lab/transcriptomics/tools/seurat_functions.R")
source("C:/Ryan/GitHub/trachtenberg-lab/transcriptomics/tools/seurat_integration_functions.R")
source("C:/Ryan/GitHub/trachtenberg-lab/transcriptomics/xgboost/xgboost_train.R")
source("C:/Ryan/GitHub/trachtenberg-lab/transcriptomics/xgboost/plottingFxns.R")

obj.mouse <- readRDS("E:/Transcriptomics_V1/Mouse/seurat/mouse_v1_P38_glutamatergic_processed.rds")
Idents(obj.mouse) <- "subclass"
obj.mouse$species <- "Mouse"
obj.opossum <- readRDS("E:/Transcriptomics_V1/Opossum/seurat/opossum_v1_glutamatergic_processed.rds")
Idents(obj.opossum) <- "subclass"
obj.opossum$species <- "Opossum"

common.features <- intersect(rownames(obj.mouse), rownames(obj.opossum))
obj.combined <- merge(obj.mouse[common.features,], y = obj.opossum[common.features,])
obj.list <- SplitObject(obj.combined, split.by = "species")

obj.list[["Mouse"]] <- SCTransform(obj.list[["Mouse"]], vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 1, algorithm = 4, method = "igraph")

obj.list[["Opossum"]] <- SCTransform(obj.list[["Opossum"]], vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 1, algorithm = 4, method = "igraph")

integration.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 5000)
# obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = integration.features)
# anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", anchor.features = integration.features)
# obj.combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# obj.sct.list <- SplitObject(obj.combined.sct, split.by = "dataset")

# Define a wrapper function to train model on common highly variable genes
TrainXGBoost = function(train, test, train.clusters = "", test.clusters = "", features = c()){
  
  Idents(train) = train.clusters
  Idents(test) = test.clusters
  
  AC_model <- TrainModel(train, training_genes = features)
  
  return(AC_model)
}

obj.list.sub <- list()
obj.list.sub[["Mouse"]] <- SubsampleObject(obj.list[["Mouse"]], "subclass", 200)
obj.list.sub[["Opossum"]] <- SubsampleObject(obj.list[["Opossum"]], "subclass", 200)

# Train model
MO_model <- TrainXGBoost(obj.list.sub[["Opossum"]], 
                         obj.list.sub[["Mouse"]], 
                         train.clusters = "subclass", 
                         test.clusters = "subclass", 
                         features = integration.features)

## Save model
# saveRDS(MO_model, "E:/xgboost/data/train_Mouse_test_Opossum.rds")

## Load trained model
# MO_model = readRDS("E:/xgboost/data/train_Mouse_test_Opossum.rds")

# Apply model
Idents(obj.list.sub[["Mouse"]]) <- "subclass"
Idents(obj.list.sub[["Opossum"]]) <- "subclass"
confusion_matrix <- BuildConfusionMatrix(obj.list.sub[["Mouse"]], obj.list.sub[["Opossum"]], model = MO_model)

# Visualize mappings
plotConfusionMatrix(confusion_matrix, plot.return = TRUE, stagger.threshold = 15)

# ggsave("G:/Shared drives/Opossum transcriptomics/figures/Confusion_Matrix_XGBoost.png", conf_plot, width=7.2, height=4.5, dpi=500)
