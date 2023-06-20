
library("Seurat")
library("anndata")
options(Seurat.object.assay.version = "v5")

data <- read_h5ad("E:/Mouse_V1/P28NR/Mouse_V1_P28NR_Glut.h5ad")
obj <- CreateSeuratObject(counts = t(data$layers[["counts"]]),
                          assay = "RNA",
                          meta.data = data$obs)
obj <- SetAssayData(object = obj, assay = "RNA", slot = "data",
                    new.data = t(data$layers[["logCPT"]]))
obj <- SetAssayData(object = obj, assay = "RNA", slot = "scale.data",
                    new.data = t(data$layers[["z-score"]]))
obj[["RNA"]] <- AddMetaData(object = obj[["RNA"]], metadata = data$var)

Idents(obj) <- obj[["type"]]

hvg <- obj[["RNA"]][["highly_variable"]]
VariableFeatures(obj) <- rownames(obj)[hvg[, "highly_variable"]]

embedding <- matrix(unlist(data$obsm["X_umap"]), ncol = 2, byrow = FALSE)
rownames(embedding) <- data$obs_names
colnames(embedding) <- c("umap_1", "umap_2")
obj[["umap"]] <- CreateDimReducObject(embedding, assay = "RNA", key = "umap_")

embedding <- matrix(unlist(data$obsm["X_pca"]), ncol = 50, byrow = FALSE)
loading <- matrix(unlist(data$varm["PCs"]), ncol = 50, byrow = FALSE)
rownames(embedding) <- data$obs_names
colnames(embedding) <- sprintf("pc_%s", seq(1:50))
rownames(loading) <- data$var_names
colnames(loading) <- sprintf("pc_%s", seq(1:50))
obj[["pca"]] <- CreateDimReducObject(embedding, loading,
                                     assay = "RNA", key = "pca_")

embedding <- matrix(unlist(data$obsm["X_harmony"]), ncol = 50, byrow = FALSE)
rownames(embedding) <- data$obs_names
colnames(embedding) <- sprintf("pc_%s", seq(1:50))
obj[["hmny"]] <- CreateDimReducObject(embedding, assay = "RNA", key = "hmny_")
