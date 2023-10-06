# XGBoost classification

library(Seurat)
library(reticulate)
library(scrubletR)
library(ggplot2)
library(cowplot)
library(dplyr)
library(SeuratDisk)

# Source scripts with functions
source("E:/xgboost/xgboost_train.R")
source("E:/xgboost/plottingFxns.R")

obj.mouse <- LoadH5Seurat("E:/Mouse_M1/seurat/mouse_m1_gabaergic.h5seurat")
Idents(obj.mouse) <- "type"
obj.treeshrew <- LoadH5Seurat("E:/Tree_Shrew_M1/seurat/treeshrew_m1_gabaergic_sct.h5seurat")
Idents(obj.treeshrew) <- "ts_integrated_snn_res.1"
obj.treeshrew$dataset <- "Tree_Shrew_M1"

common.features <- intersect(rownames(obj.mouse), rownames(obj.treeshrew))
obj.combined <- merge(obj.mouse[common.features,], y = obj.treeshrew[common.features,])
obj.list <- SplitObject(obj.combined, split.by = "dataset")

obj.list[["Mouse_M1"]] <- SCTransform(obj.list[["Mouse_M1"]], vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 1, algorithm = 4, method = "igraph")

obj.list[["Tree_Shrew_M1"]] <- SCTransform(obj.list[["Tree_Shrew_M1"]], vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 1, algorithm = 4, method = "igraph")

integration.features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = integration.features)
anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", anchor.features = integration.features)
obj.combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

obj.sct.list <- SplitObject(obj.combined.sct, split.by = "dataset")

# Define a wrapper function to train model on common highly variable genes
TrainXGBoost = function(train, test, train.clusters = "", test.clusters = "", features = c()){
  
  Idents(train) = train.clusters
  Idents(test) = test.clusters
  
  AC_model <- TrainModel(train, training_genes = features)
  
  return(AC_model)
}

# Train model
obj.sct.list[["Mouse_M1"]] <- obj.sct.list[["Mouse_M1"]][, obj.sct.list[["Mouse_M1"]]$type != "NA"]

MO_model <- TrainXGBoost(obj.sct.list[["Mouse_M1"]], 
                         obj.sct.list[["Tree_Shrew_M1"]], 
                         train.clusters = "type", 
                         test.clusters = "ts_integrated_snn_res.1", 
                         features = integration.features) # 

## Save model
# saveRDS(MO_model, "E:/xgboost/data/train_Mouse_test_Opossum.rds")

## Load trained model
# MO_model = readRDS("E:/xgboost/data/train_Mouse_test_Opossum.rds")

# Apply model
Idents(obj.sct.list[["Mouse_M1"]]) <- "type"
Idents(obj.sct.list[["Tree_Shrew_M1"]]) <- "ts_integrated_snn_res.1"
confusion_matrix <- BuildConfusionMatrix(obj.sct.list[["Tree_Shrew_M1"]], obj.sct.list[["Mouse_M1"]], model = MO_model)

# Visualize mappings
plotConfusionMatrix(confusion_matrix, plot.return = TRUE, stagger.threshold = 15)

# ggsave("G:/Shared drives/Opossum transcriptomics/figures/Confusion_Matrix_XGBoost.png", conf_plot, width=7.2, height=4.5, dpi=500)
