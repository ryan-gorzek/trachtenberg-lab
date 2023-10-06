# Example XGBoost classification

# Source scripts with functions
source("xgboost_train.R")
source("plottingFxns.R")

# Define a wrapper function to train model on common highly variable genes
TrainXGBoost = function(train, test, train.clusters = "seurat_clusters", test.clusters = "seurat_clusters"){
  
  Idents(train) = "seurat_clusters"
  Idents(test) = "seurat_clusters"
  
  DefaultAssay(train) = "RNA"
  DefaultAssay(test) = "RNA"
  
  train = FindVariableFeatures(train, selection.method = "vst", nfeatures = 2000)
  test = FindVariableFeatures(test, selection.method = "vst", nfeatures = 2000)
  common_HVGs <- intersect(VariableFeatures(train), VariableFeatures(test))
  message("Using ", length(common_HVGs), " common highly variable genes...\n")
  
  AC_model <- TrainModel(train, training_genes = common_HVGs)
  
  return(AC_model)
}

# Train model
AC_model <- TrainXGBoost(MouseACref, PengAC)

# Save model
saveRDS(AC_model, "../../data/train_MouseACref_test_Peng.rds")

# Load trained model
AC_model = readRDS("../../data/train_MouseACref_test_Peng.rds")

# Apply model
train_MouseACref_test_Peng <- BuildConfusionMatrix(PengAC, MouseACref, model = AC_model)

# Visualize mappings
plotConfusionMatrix(train_MouseACref_test_Peng, plot.return = TRUE, stagger.threshold = 15)