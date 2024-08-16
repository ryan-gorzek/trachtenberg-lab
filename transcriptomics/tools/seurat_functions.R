
PreprocessData <- function(sample_IDs, data_path, project_name, mapping_path) {
  
  # Load the data.
  print("Loading 10x data...")
  objs <- c()
  genes <- list()
  
  for (sample in sample_IDs) {
    
    temp.obj.path <- paste(data_path, sample, "/outs/filtered_feature_bc_matrix/", sep = "")
    temp.obj.data <- Read10X(temp.obj.path)
    temp.obj <- CreateSeuratObject(counts = temp.obj.data, project = project_name)
    temp.obj$sample <- sample
    temp.obj <- scrublet_R(seurat_obj = temp.obj)
    objs <- append(objs, temp.obj)
    
  }
  
  obj <- merge(objs[[1]], y = objs[2:length(objs)], add.cell.ids = sample_IDs, project = project_name)
  genes$pre.map <- rownames(obj)
  
  # Map gene names/IDs if specified.
  print("Mapping genes...")
  if (!is.na(mapping_path)) {
    genes.mapping <- read.csv(mapping_path)
    para.idx <- genes.mapping$Gene.stable.ID %in% unique(genes.mapping$Gene.stable.ID[duplicated(genes.mapping$Gene.stable.ID)])
    genes.mapping <- genes.mapping[!para.idx,]
    genes.mapping.other <- as.list(genes.mapping[, 3])
    genes.mapping.self <- as.list(genes.mapping[, 1])
    ids.mapping.self <- as.list(genes.mapping[, 2])
    genes.self <- rownames(obj)
    for (gene in genes.mapping.other) {
      
      idx.other <- which(genes.mapping.other %in% gene)
      
      if (length(idx.other) == 1) {
        
        gene.self <- genes.mapping.self[idx.other]
        id.self <- ids.mapping.self[idx.other]
        
        if (gene.self == "") {
          
          idx.self <- which(genes.self %in% id.self)
          genes.self[idx.self] <- gene
          
        } else {
          
          idx.self <- which(genes.self %in% gene.self)
          genes.self[idx.self] <- gene
          
        }
      }
    }
  }
  else { genes.self <- rownames(obj) }
  
  # Rebuild Seurat object.
  print("Rebuilding object...")
  obj.df <- as.data.frame(as.matrix(obj[["RNA"]]@counts))
  rownames(obj.df) <- genes.self
  obj.temp <- CreateSeuratObject(counts = obj.df, meta.data = obj[[]])
  obj <- obj.temp
  genes$post.map <- rownames(obj)
  
  # Filter cells and genes.
  print("Filtering cells and genes...")
  obj.prefilt <- obj
  cell_mask <- Reduce(intersect,list(WhichCells(obj, expression = nFeature_RNA > 700),
                                     WhichCells(obj, expression = nFeature_RNA < 6500),
                                     WhichCells(obj, expression = nCount_RNA < 40000)))
  
  gene_mask <- rownames(obj)[Matrix::rowSums(obj[["RNA"]]@counts > 0) > 8]
  
  obj <- subset(obj, features = gene_mask, cells = cell_mask)
  
  genes$post.filt <- rownames(obj)
  
  return(list(obj = obj, obj.prefilt = obj.prefilt, genes = genes))
  
}

PlotQC <- function(data) {
  
  # Create a custom theme for the violin plots
  custom_theme_vln <- theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.position = "none", # Remove legends
    axis.title.x = element_blank(), # Remove x-axis labels
    axis.text.x = element_text(size = 8) # Remove x-axis tick labels
  )
  
  # Create a custom theme for the scatter plots
  custom_theme_scatter <- theme(
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 6),
    legend.position = "none", # Remove legends
    axis.title.x = element_text(size = 8), # Add x-axis labels back
    axis.text.x = element_text(size = 6) # Add x-axis tick labels back
  )
  
  # Create individual plots for pre-filtered data with the custom theme
  vln_plot1_pre <- VlnPlot(data$obj.prefilt, features = c("nFeature_RNA"), group.by = "sample", raster = FALSE) + ylim(0, 10000) + custom_theme_vln
  vln_plot2_pre <- VlnPlot(data$obj.prefilt, features = c("nCount_RNA"), group.by = "sample", raster = FALSE) + ylim(0, 40000) + custom_theme_vln
  scatter_plot_pre <- FeatureScatter(data$obj.prefilt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = FALSE) + xlim(0, 40000) + ylim(0, 10000) + custom_theme_scatter + coord_fixed(ratio = 4)
  
  # Create individual plots for post-filtered data with the custom theme
  vln_plot1_post <- VlnPlot(data$obj, features = c("nFeature_RNA"), group.by = "sample", raster = FALSE) + ylim(0, 10000) + custom_theme_vln
  vln_plot2_post <- VlnPlot(data$obj, features = c("nCount_RNA"), group.by = "sample", raster = FALSE) + ylim(0, 40000) + custom_theme_vln
  scatter_plot_post <- FeatureScatter(data$obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = FALSE) + xlim(0, 40000) + ylim(0, 10000) + custom_theme_scatter + coord_fixed(ratio = 4)
  
  final_layout_feature <- vln_plot1_pre | vln_plot1_post
  final_layout_count <- vln_plot2_pre | vln_plot2_post
  final_layout_scatter <- scatter_plot_pre | scatter_plot_post
  
  # Display the final layout
  print(final_layout_feature)
  print(final_layout_count)
  print(final_layout_scatter)
  
}

ClusterSCT <- function(obj, resolutions) {
  
  obj <- SCTransform(obj, vst.flavor = "v2", return.only.var.genes = FALSE, verbose = FALSE) %>%
         RunPCA(npcs = 30, verbose = FALSE) %>%
         FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
         RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)
  
  for (res in resolutions) {
    obj <- FindClusters(obj, resolution = res, algorithm = 4, method = "igraph")
  }
  
  return(obj)
  
}

NormalizePCA <- function(obj, nfeatures = 3000, npcs = 30) {
  
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
  all.genes <- rownames(obj)
  obj <- ScaleData(obj, features = all.genes)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = npcs)
  return(obj)
  
}

PlotClusters <- function(obj, group.id) {
  
  if (!missing(group.id)) {
    Idents(obj) <- group.id
  }
  obj$active.ident <- obj@active.ident
  dimplot1 <- DimPlot(obj, reduction = "umap", label = TRUE, raster = FALSE) + NoLegend() + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
  dimplot2 <- DimPlot(obj, reduction = "umap", group.by = "sample", raster = FALSE, shuffle = TRUE) + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
  cluster.sample <- table(obj$sample, obj$active.ident) %>%
    as.data.frame.matrix() %>%
    rownames_to_column(var = "sample")
  cluster.sample[-1] <- lapply(cluster.sample[-1], function(x) x/sum(x))
  cluster.sample <- cluster.sample %>%
    pivot_longer(
      cols = -c("sample"),
      names_to = "cluster",
      values_to = "count"
    )
  cluster.sample$cluster <- factor(cluster.sample$cluster, levels = unique(cluster.sample$cluster))
  barplot1 <- ggplot(cluster.sample, aes(x=cluster, y=count, fill=sample)) +
    geom_bar(stat="identity") +
    theme_minimal()
  dimplot3 <- DimPlot(obj, reduction = "umap", group.by = "predicted_doublets", raster = FALSE) + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
  # Summarize doublets by cluster
  df <- obj[[]]
  df$active.ident <- obj@active.ident
  df_summary <- df %>%
    dplyr::group_by(active.ident, predicted_doublets) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::mutate(fraction = count / sum(count))
  # Create the stacked bar plot
  barplot2 <- ggplot(df_summary, aes(x = active.ident, y = fraction, fill = predicted_doublets)) +
    geom_bar(stat = "identity") +
    labs(x = "Clusters", y = "Doublet Fraction", fill = "Value") +
    theme_minimal()
  featplot1 <- FeaturePlot(obj, "nFeature_RNA", raster = FALSE) + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
  vlnplot1 <- VlnPlot(obj, "nFeature_RNA")
  featplot2 <- FeaturePlot(obj, "nCount_RNA", raster = FALSE) + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
  vlnplot2 <- VlnPlot(obj, "nCount_RNA")
  
  print(dimplot1)
  print(dimplot2)
  print(barplot1)
  print(dimplot3)
  print(barplot2)
  print(featplot1)
  print(vlnplot1)
  print(featplot2)
  print(vlnplot2)
  
}

PlotFeatures <- function(obj, features) {
  
  if (class(features) == "list") { features <- unname(unlist(features)) }
  
  for (feat in features) {
    if (feat %in% rownames(obj)) {
      featplot <- FeaturePlot(obj, reduction = "umap", features = feat, raster = FALSE) + NoLegend() + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
      featplot <- suppressWarnings(featplot)
      print(featplot)
    }
  }
  
}

RemoveByUMAP <- function(seurat_obj, umap1_range, umap2_range) {
  # Check if UMAP reduction exists
  if (!"umap" %in% names(seurat_obj@reductions)) {
    stop("UMAP reduction not found in the Seurat object. Please ensure UMAP has been run.")
  }
  
  # Extract UMAP coordinates
  umap_coords <<- Embeddings(seurat_obj, reduction = "umap")
  
  # Find cells within the specified UMAP1 and UMAP2 range
  cells_to_remove <<- which(
    umap_coords[, 1] >= umap1_range[1] & umap_coords[, 1] <= umap1_range[2] &
      umap_coords[, 2] >= umap2_range[1] & umap_coords[, 2] <= umap2_range[2]
  )
  
  # Get cell names to remove
  cells_to_remove <<- rownames(umap_coords)[cells_to_remove]
  
  # Remove cells from Seurat object
  cell_names <- colnames(seurat_obj)
  seurat_obj <- subset(seurat_obj, cells = cell_names[!cell_names %in% cells_to_remove])
  
  return(seurat_obj)
}

SubclassByIdent <- function(obj, subclass.idx) {
  
  for (id in names(subclass.idx)) {
    
    Idents(obj) <- id
    
    id.num <- strsplit(id, ".", fixed = TRUE)[[1]]
    if (length(id.num) == 2) { sbcl.col <- paste0("subclass.", id.num[2]) }
    else if (length(id.num) == 3) { sbcl.col <- paste0("subclass.", id.num[2], ".", id.num[3]) }
    else { sbcl.col <- paste0("subclass.", id.num[1]) }
    
    if ((sbcl.col %in% colnames(obj[[]]) == FALSE)) { obj[[sbcl.col]] <- NA }
    
    for (sbcl in names(subclass.idx[[id]])) {
      
      if (!is.null(subclass.idx[[id]][[sbcl]])) {
        
        cell.names <- WhichCells(obj, ident = subclass.idx[[id]][[sbcl]])
        obj[[sbcl.col]][cell.names,] <- sbcl
        
      }
    }
  }
  
  return(obj)
  
}

IdentMarkerDict <- function(obj, subclass.labels, ident.labels, save.path) {
  
  marker.dict <- list()
  
  for (sbcl in subclass.labels) {
    
    marker.dict[[sbcl]] <- list()
    
    for (id in ident.labels) {
      
      id.num <- strsplit(id, ".", fixed = TRUE)[[1]]
      if (length(id.num) == 2) { sbcl.col <- paste0("subclass.", id.num[2]) }
      else if (length(id.num) == 3) { sbcl.col <- paste0("subclass.", id.num[2], ".", id.num[3]) }
      else { sbcl.col <- paste0("subclass.", id.num[1]) }
      
      Idents(obj) <- sbcl.col
      if (sbcl %in% levels(obj)) {
        
        obj.sbcl.id <- subset(obj, idents = sbcl)
        
        DefaultAssay(obj.sbcl.id) <- "SCT"
        Idents(obj.sbcl.id) <- id
        if (length(levels(obj.sbcl.id)) > 1) {
          
          marker.dict[[sbcl]][[id]] <- FindAllMarkers(obj.sbcl.id, only.pos = TRUE, logfc.threshold = 0.1)
          
        }
      }
    }
  }
  
  saveRDS(marker.dict, save.path)
  return(marker.dict)
  
}

SubclassMarkerDict <- function(obj, subclass.labels, save.path) {
  
  marker.dict <- list()
  
  for (sbcl in subclass.labels) {
    
    marker.dict[[sbcl]] <- list()
    
    DefaultAssay(obj) <- "SCT"
    Idents(obj) <- sbcl
      
    marker.dict[[sbcl]] <- FindAllMarkers(obj, only.pos = TRUE, logfc.threshold = 0.1)

  }
  
  saveRDS(marker.dict, save.path)
  return(marker.dict)
  
}

SaveDotPlots <- function(obj, markers, subclass.labels, ident.labels, savepath, ens.id) {
  
  for (sbcl in subclass.labels) {
    
    folder_path <- paste0(savepath, gsub("/", "", sbcl), "/")
    make_folder(folder_path)
    
    for (id in ident.labels) {
      
      id.num <- strsplit(id, ".", fixed = TRUE)[[1]]
      if (length(id.num) == 2) { sbcl.col <- paste0("subclass.", id.num[2]) }
      else if (length(id.num) == 3) { sbcl.col <- paste0("subclass.", id.num[2], ".", id.num[3]) }
      else { sbcl.col <- paste0("subclass.", id.num[1]) }
      
      Idents(obj) <- sbcl.col
      if (sbcl %in% levels(obj)) {
        
        obj.sbcl.id <- subset(obj, idents = sbcl)
        
        DefaultAssay(obj.sbcl.id) <- "SCT"
        Idents(obj.sbcl.id) <- id
        if (!all(is.na(as.numeric(levels(obj.sbcl.id))))) {
          levels(obj.sbcl.id) <- factor(rev(as.character(sort(as.numeric(levels(obj.sbcl.id))))))
        }
        else { levels(obj.sbcl.id) <- factor(rev(levels(obj.sbcl.id))) }
        if (length(levels(obj.sbcl.id)) > 1) {

          id.path <- paste0(folder_path, gsub("/", "", id), "/")
          make_folder(id.path)
          
          marker_sets <- list(c(1:20))
          for (set in marker_sets) {
            
            all.markers.FC <- list()
            all.markers.PD <- list()
            
            for (type in levels(obj.sbcl.id)) {
              
              all.markers <- markers[[sbcl]][[id]]
              gene.counts <- table(all.markers$gene)
              unique.markers <- all.markers[all.markers$gene %in% names(gene.counts)[gene.counts == 1],]
              type.markers <- unique.markers[unique.markers$cluster == type,]
              type.markers$pct.diff <- type.markers$pct.1 - type.markers$pct.2
              all.markers.FC[[type]] <- top_genes_desc(type.markers, "avg_log2FC", set)
              all.markers.PD[[type]] <- top_genes_desc(type.markers, "pct.diff", set)
              
            }
            
            # make plots
            plot.FC <- DotPlot(obj.sbcl.id, features = rev(all.markers.FC), cols = c("lightgrey", "red"), scale = FALSE) + 
              theme(axis.text.x = element_text(angle = 90)) + 
              theme(panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA))
            plot.FC$data$features.plot <- factor(plot.FC$data$features.plot, levels = levels(plot.FC$data$features.plot), 
                                                 labels = unname(strip_if_contains(unlist(rev(all.markers.FC)), ens.id, paste0(ens.id, "000000"))))
            ggsave(paste0(id.path, sbcl.col, "_DotPlot_FC.png"), plot = plot.FC, width = 24, dpi = 300)
            
            plot.PD <- DotPlot(obj.sbcl.id, features = rev(all.markers.PD), cols = c("lightgrey", "red"), scale = FALSE) +
              theme(axis.text.x = element_text(angle = 90)) +
              theme(panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA))
            plot.PD$data$features.plot <- factor(plot.PD$data$features.plot, levels = levels(plot.PD$data$features.plot), 
                                                 labels = unname(strip_if_contains(unlist(rev(all.markers.PD)), ens.id, paste0(ens.id, "000000"))))
            ggsave(paste0(id.path, sbcl.col, "_DotPlot_PD.png"), plot = plot.PD, width = 24, dpi = 300)
            
          }
        }
      }
    }
  }
  
}

SaveFeaturePlots <- function(obj, markers, subclass.labels, ident.labels, savepath) {
  
  get_plot_limits <- function(plot) {
    gg_build <- ggplot_build(plot)
    xlims <- gg_build$layout$panel_scales_x[[1]]$range$range
    ylims <- gg_build$layout$panel_scales_y[[1]]$range$range
    return(list(xlims = xlims, ylims = ylims))
  }
  
  for (sbcl in subclass.labels) {
    
    folder_path <- paste0(savepath, gsub("/", "", sbcl), "/")
    make_folder(folder_path)
    
    for (id in ident.labels) {
      
      id.num <- strsplit(id, ".", fixed = TRUE)[[1]]
      if (length(id.num) == 2) { sbcl.col <- paste0("subclass.", id.num[2]) }
      else if (length(id.num) == 3) { sbcl.col <- paste0("subclass.", id.num[2], ".", id.num[3]) }
      else { sbcl.col <- paste0("subclass.", id.num[1]) }
      
      Idents(obj) <- sbcl.col
      if (sbcl %in% levels(obj)) {
        
        obj.sbcl.id <- subset(obj, idents = sbcl)
        
        DefaultAssay(obj.sbcl.id) <- "SCT"
        Idents(obj.sbcl.id) <- id
        if (!all(is.na(as.numeric(levels(obj.sbcl.id))))) {
          levels(obj.sbcl.id) <- factor(rev(as.character(sort(as.numeric(levels(obj.sbcl.id))))))
        }
        else { levels(obj.sbcl.id) <- factor(rev(levels(obj.sbcl.id))) }
        if (length(levels(obj.sbcl.id)) > 1) {
          
          id.path <- paste0(folder_path, gsub("/", "", id), "/")
          make_folder(id.path)
          
          marker_sets <- list(c(1:20))
          for (set in marker_sets) {
            
            all.markers.FC <- list()
            all.markers.PD <- list()
            
            for (type in levels(obj.sbcl.id)) {
              
              all.markers <- markers[[sbcl]][[id]]
              gene.counts <- table(all.markers$gene)
              unique.markers <- all.markers[all.markers$gene %in% names(gene.counts)[gene.counts == 1],]
              type.markers <- unique.markers[unique.markers$cluster == type,]
              type.markers$pct.diff <- type.markers$pct.1 - type.markers$pct.2
              all.markers.FC[[type]] <- top_genes_desc(type.markers, "avg_log2FC", set)
              all.markers.PD[[type]] <- top_genes_desc(type.markers, "pct.diff", set)
              
              # make feature plots
              features.FC <- unlist(all.markers.FC[[type]])
              features.PD <- unlist(all.markers.PD[[type]])
              
              for (feature_set in list(features.FC, features.PD)) {
                feature_subset <- feature_set[1:min(20, length(feature_set))]
                feature_subset <- feature_subset[!is.na(feature_subset)]
                if (length(feature_subset) > 0) {
                  plots <- FeaturePlot(obj.sbcl.id, features = feature_subset, cols = c("lightgrey", "red"), ncol = 5)
                  
                  # Extract plot limits
                  plot_limits <- get_plot_limits(plots[[1]])
                  xlims <- plot_limits$xlims
                  ylims <- plot_limits$ylims
                  
                  x.range <- diff(xlims)
                  y.range <- diff(ylims)
                  
                  max.range <- max(x.range, y.range)
                  
                  if (x.range < max.range) {
                    xlims <- mean(xlims) + c(-1, 1) * (max.range / 2)
                  }
                  
                  if (y.range < max.range) {
                    ylims <- mean(ylims) + c(-1, 1) * (max.range / 2)
                  }
                  
                  # Adjust the plots with new limits and coord_equal
                  for (i in 1:length(plots$patches$plots)) {
                    plots[[i]] <- plots[[i]] + coord_equal(xlim = xlims, ylim = ylims)
                  }
                  
                  file_prefix <- ifelse(identical(feature_set, features.FC), "FeaturePlot_FC", "FeaturePlot_PD")
                  ggsave(paste0(id.path, sbcl.col, "_", file_prefix, "_id-", gsub("/", "", type), ".png"), plot = plots, width = 24, height = 16, dpi = 300)
                }
              }
            }
          }
        }
      }
    }
  }
}

top_genes_two <- function(df, sort_column_asc, sort_column_desc, idx) {
  
  # Sort the dataframe by the specified columns: one ascending and one descending
  sorted_df <- df[order(df[[sort_column_asc]], -df[[sort_column_desc]]), ]
  
  # Get the top n rownames
  top_n <- sorted_df$gene[idx]
  
  return(top_n[!is.na(top_n)])
}

top_genes_asc <- function(df, sort_column_asc, idx) {
  
  # Sort the dataframe by the specified columns: one ascending and one descending
  sorted_df <- df[order(df[[sort_column_asc]]), ]
  
  # Get the top n rownames
  top_n <- sorted_df$gene[idx]
  
  return(top_n[!is.na(top_n)])
}

top_genes_desc <- function(df, sort_column_desc, idx) {
  
  # Sort the dataframe by the specified columns: one ascending and one descending
  sorted_df <- df[order(-df[[sort_column_desc]]), ]
  
  # Get the top n rownames
  top_n <- sorted_df$gene[idx]
  
  return(top_n[!is.na(top_n)])
}

make_folder <- function(folder_path) {
  if (!dir.exists(folder_path)) {
    # Create the folder
    dir.create(folder_path)
  }
}

strip_if_contains <- function(vector, search, strip) {
  # Initialize an empty vector to store results
  result_vector <- vector
  
  # Loop through each element of the vector
  for (i in seq_along(vector)) {
    if (grepl(search, vector[i])) {
      result_vector[i] <- paste0("U#", sub(paste0("^", strip), "", vector[i]))
    }
  }
  
  return(result_vector)
}

sort_idents <- function(vec) {
  suffixes <<- sub(".*_", "", vec)
  
  # Determine if suffixes are numeric or alphabetic
  if (all(grepl("^[0-9]+$", suffixes))) {
    sorted_vec <<- vec[order(as.numeric(suffixes))]
  } else {
    sorted_vec <<- vec[order(suffixes)]
  }
  
  return(factor(rev(sorted_vec), levels = rev(sorted_vec)))
}

sort_by_reference <- function(vec, ref) {
  # Find the intersection of vec and ref
  common_elements <- intersect(ref, vec)
  
  # Order vec according to the position of common elements in ref
  sorted_vec <- vec[order(match(vec, common_elements, nomatch = Inf))]
  
  # Remove elements that are not in common_elements
  sorted_vec <- sorted_vec[sorted_vec %in% common_elements]

  return(factor(rev(sorted_vec), levels = rev(sorted_vec)))
}

SaveIdentConfusionMatrices <- function(obj, subclass.labels, ident.labels, savepath) {
  
  DefaultAssay(obj) <- "SCT"
  
  for (sbcl in subclass.labels) {
    
    folder_path <- paste0(savepath, gsub("/", "", sbcl), "/")
    make_folder(folder_path)
    
    for (id in ident.labels) {
      
      id.num <- strsplit(id, ".", fixed = TRUE)[[1]]
      if (length(id.num) == 2) { sbcl.col <- paste0("subclass.", id.num[2]) }
      else if (length(id.num) == 3) { sbcl.col <- paste0("subclass.", id.num[2], ".", id.num[3]) }
      else { sbcl.col <- paste0("subclass.", id.num[1]) }
      
      Idents(obj) <- sbcl.col
      if (sbcl %in% levels(obj)) {
        
        obj.sbcl.id <- subset(obj, idents = sbcl)
        
        DefaultAssay(obj.sbcl.id) <- "SCT"
        Idents(obj.sbcl.id) <- id
        levels(obj.sbcl.id) <- sort_idents(levels(obj.sbcl.id))
        if (length(levels(obj.sbcl.id)) > 1) {
          
          id.path <- paste0(folder_path, gsub("/", "", id), "/")
          make_folder(id.path)

          mdl <- TrainModel(obj.sbcl.id, training_genes = VariableFeatures(obj))
          
          if (!is.null(mdl$confusion)) {
            confusion_plot <- mdl$confusion +
              coord_equal() + 
              theme(axis.text.x = element_text(angle = 90),
                    panel.background = element_rect(fill = "white", color = NA), 
                    plot.background = element_rect(fill = "white", color = NA))
            
            ggsave(paste0(id.path, sbcl.col, "_XGBoost.png"), plot = confusion_plot, width = 5, height = 5, dpi = 300)
          }
        }
      }
    }
  }
}

SaveSubclassConfusionMatrices <- function(obj, subclass.cols, subclass.order, savepath) {
  
  DefaultAssay(obj) <- "SCT"
  
  for (sbcl in subclass.cols) {
    
    folder_path <- paste0(savepath)
    make_folder(folder_path)
      
    DefaultAssay(obj) <- "SCT"
    Idents(obj) <- sbcl
    levels(obj) <- sort_by_reference(levels(obj), subclass.order)
    
    mdl <- TrainModel(obj, training_genes = VariableFeatures(obj))
          
    if (!is.null(mdl$confusion)) {
      confusion_plot <- mdl$confusion +
                        coord_equal() + 
                        theme(axis.text.x = element_text(angle = 90),
                              panel.background = element_rect(fill = "white", color = NA), 
                              plot.background = element_rect(fill = "white", color = NA))
            
      ggsave(paste0(folder_path, sbcl, "_XGBoost.png"), plot = confusion_plot, width = 5, height = 5, dpi = 300)

    }
  }
}

LabelCells <- function(obj, subclass_resolution) {

  # Iterate over each subclass and resolution
  for (subclass in names(subclass_resolution)) {
    resolution <- subclass_resolution[[subclass]]
    
    # Extract the subclass and cluster columns based on the resolution
    subclass_col <- paste0("subclass.", resolution)
    cluster_col <- paste0("SCT_snn_res.", resolution)
    
    # Get the cells belonging to the current subclass
    cells <- obj@meta.data[[subclass_col]] == subclass
    
    # Assign the subclass label to the 'subclass' column
    obj@meta.data$subclass[cells] <- subclass
    
    # Extract the cluster labels for these cells
    clusters <- obj@meta.data[[cluster_col]][cells]
    
    # Get the cluster sizes in decreasing order
    cluster_sizes <- sort(table(clusters), decreasing = TRUE)
    
    # Create a mapping from original cluster labels to new type labels
    cluster_to_type <- setNames(paste0(subclass, "_", seq_along(cluster_sizes)), names(cluster_sizes))
    
    # Assign the new type labels based on the cluster sizes
    if (sum(cluster_sizes > 0) > 1) {
      obj@meta.data$type[cells] <- cluster_to_type[as.character(clusters)]
    } else { obj@meta.data$type[cells] <- subclass }
    obj@meta.data$subclass.type[cells] <- subclass
  }
  
  if (any(is.na(obj@meta.data$type))) { warning("NAs present in type column...") }
  return(obj)

}

IdentBySample <- function(obj, y_limits = c(0, 0.50)) {
  
  # Assuming your dataframe is named df with columns 'subclass' and 'sample'
  df <- obj[[]]
  df$active.ident <- obj@active.ident
  df$smpl <- df$sample
  specified_order <- levels(obj)

  # Check for NAs and warn the user
  na_samples <- df %>% filter(is.na(active.ident)) %>% count(smpl)
  if(nrow(na_samples) > 0) {
    warning("The following samples contained NAs and were dropped: ",
            paste(na_samples$sample, " (", na_samples$n, ")", collapse = "\n"))
  }
  
  # Drop NAs
  df <- df %>% drop_na(active.ident)
  
  # Calculate the relative proportions of each active.ident by sample
  relative_proportions <- df %>%
    group_by(sample, active.ident) %>%
    dplyr::summarize(count = n(), .groups = 'drop') %>%
    ungroup() %>%
    group_by(sample) %>%
    mutate(Proportion = count / sum(count))
  
  # Calculate the median proportions across samples
  median_proportions <- relative_proportions %>%
    group_by(active.ident) %>%
    summarize(MedianProportion = median(Proportion))
  
  # Rename the columns for clarity
  colnames(median_proportions) <- c("active.ident", "MedianProportion")
  
  # Convert active.ident to a factor and specify the order
  median_proportions$active.ident <- factor(median_proportions$active.ident, levels = specified_order)
  relative_proportions$active.ident <- factor(relative_proportions$active.ident, levels = specified_order)
  
  # # Set y-axis limits
  # y_limits <- c(0, 0.60) # Modify these values as needed
  
  # Calculate the maximum proportion for each active.ident
  max_proportions <- relative_proportions %>%
    group_by(active.ident) %>%
    summarize(MaxProportion = max(Proportion))
  
  # Merge the max proportions with the median proportions
  plot_data <- merge(median_proportions, max_proportions, by = "active.ident")
  
  # Create the barplot
  p <- ggplot(plot_data, aes(x = active.ident, y = MedianProportion, fill = active.ident)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = scales::percent(MedianProportion, accuracy = 0.1), 
                  y = MaxProportion + 0.01), # Place text labels slightly above the highest scatter point
              vjust = -0.5) + # Adjust the position of the text
    geom_point(data = relative_proportions, aes(x = active.ident, y = Proportion), 
               color = "black") + # Scatter the different samples
    theme_minimal() +
    theme(legend.position = "none", 
          aspect.ratio = 1) + # Remove the legend and make the plot square
    xlab("") +
    ylab("Relative Proportion") +
    ggtitle("") +
    scale_y_continuous(labels = scales::percent, limits = y_limits) +
    scale_x_discrete(limits = specified_order) + 
    theme(axis.text.x = element_text(size=12), 
          axis.text.y = element_text(size=12), 
          axis.title.x = element_text(size=12), 
          axis.title.y = element_text(size=12))
  
  print(p)
  
}

PlotIdentGeneCounts <- function(nested_list, subclass, clustering_res) {

  # Extract dataframe
  df <- nested_list[[subclass]][[clustering_res]]
  if (!is.null(df)) {
    # Calculate the maximum avg_log2FC value
    max_log2FC <- max(df$avg_log2FC, na.rm = TRUE)
    
    # Create a sequence of avg_log2FC values from 0.2 to the maximum or 2 (whichever is smaller)
    log2FC_grid <- seq(0.2, min(max_log2FC, 2), length.out = 100)
    
    # Calculate the number of genes for each cluster as avg_log2FC varies
    gene_counts <- df %>%
      group_by(cluster) %>%
      do(data.frame(avg_log2FC = log2FC_grid,
                    count = sapply(log2FC_grid, function(x) sum(.$avg_log2FC > x, na.rm = TRUE)))) %>%
      ungroup()
    levels(gene_counts$cluster) <- rev(sort_idents(levels(gene_counts$cluster)))
    # Plot the number of genes
    ggplot(gene_counts, aes(x = avg_log2FC, y = count, color = cluster, group = cluster)) +
      geom_line(size = 1) + # Thicker lines
      labs(title = paste("Number of DE Genes for", subclass, "by", clustering_res),
           x = "avg_log2FC",
           y = "Number of Genes",
           color = "Cluster") +
      theme_minimal() +
      scale_x_continuous(limits = c(0.2, 2)) + # Set x-axis limits
      geom_vline(xintercept = 0.2, linetype = "dotted", color = "black") +
      annotate("text", x = 0.2, y = max(gene_counts$count) * 0.9, label = "0.2", hjust = -0.2, vjust = -0.5) +
      geom_vline(xintercept = 0.5, linetype = "dotted", color = "black") +
      annotate("text", x = 0.5, y = max(gene_counts$count) * 0.9, label = "0.5", hjust = -0.2, vjust = -0.5) +
      theme(
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
        axis.text.y = element_text(size = 12)   # Increase y-axis tick label size
      )
  }
}

PlotSubclassGeneCounts <- function(nested_list, subclass.col, subclass.order) {
  
  # Extract dataframe
  df <- nested_list[[subclass.col]]
  if (!is.null(df)) {
    # Calculate the maximum avg_log2FC value
    max_log2FC <- max(df$avg_log2FC, na.rm = TRUE)
    
    # Create a sequence of avg_log2FC values from 0.2 to the maximum or 2 (whichever is smaller)
    log2FC_grid <- seq(0.2, min(max_log2FC, 2), length.out = 100)
    
    # Calculate the number of genes for each cluster as avg_log2FC varies
    gene_counts <- df %>%
      group_by(cluster) %>%
      do(data.frame(avg_log2FC = log2FC_grid,
                    count = sapply(log2FC_grid, function(x) sum(.$avg_log2FC > x, na.rm = TRUE)))) %>%
      ungroup()
    levels(gene_counts$cluster) <- rev(sort_by_reference(levels(gene_counts$cluster), subclass.order))
    # Plot the number of genes
    ggplot(gene_counts, aes(x = avg_log2FC, y = count, color = cluster, group = cluster)) +
      geom_line(size = 1) + # Thicker lines
      labs(title = paste("Number of DE Genes for Subclasses"),
           x = "avg_log2FC",
           y = "Number of Genes",
           color = "Subclass") +
      theme_minimal() +
      scale_x_continuous(limits = c(0.2, 2)) + # Set x-axis limits
      geom_vline(xintercept = 0.2, linetype = "dotted", color = "black") +
      annotate("text", x = 0.2, y = max(gene_counts$count) * 0.9, label = "0.2", hjust = -0.2, vjust = -0.5) +
      geom_vline(xintercept = 0.5, linetype = "dotted", color = "black") +
      annotate("text", x = 0.5, y = max(gene_counts$count) * 0.9, label = "0.5", hjust = -0.2, vjust = -0.5) +
      theme(
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
        axis.text.y = element_text(size = 12)   # Increase y-axis tick label size
      )
  }
}

plotGeneFractions <- function(df, gene_list) {
  # Calculate the fraction of genes in each cluster that belong to the gene_list
  fraction_data <- df %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarize(
      total_genes = n(),
      matching_genes = sum(gene %in% gene_list),
      fraction = matching_genes / total_genes
    )
  
  # Plot the results
  ggplot(fraction_data, aes(x = cluster, y = fraction)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "",
         x = "Cluster",
         y = "Fraction of Genes") +
    ylim(0, 1) +
    theme_minimal()
}

PlotSubclassGeneCountCDF <- function(de_df, all_genes, sample.name, subclass.order, subclass_colors) {
  # Extract dataframes
  df <- de_df$subclass
  
  # Initialize an empty data frame to store cumulative counts
  cdf_data <- data.frame()
  
  for (i in 1:length(subclass.order)) {
    
    sbcl <- subclass.order[i]
    
    # Filter dataframes for the specific cluster pair
    df_filtered <- df %>% filter(cluster == sbcl)
    
    # Calculate the maximum avg_log2FC value
    max_log2FC <- max(df_filtered$avg_log2FC, na.rm = TRUE)
    
    # Create a sequence of avg_log2FC values from 0.2 to the maximum or 2 (whichever is smaller)
    log2FC_grid <- seq(0.2, min(max_log2FC, 2), length.out = 50)
    
    count <- c()
    for (l in log2FC_grid) {
      de_genes <- df_filtered$gene[df_filtered$avg_log2FC > l]
      count <- c(count, length(de_genes)) }
    
    # Calculate the cumulative number of shared genes as avg_log2FC varies
    shared_gene_counts <- data.frame(
      avg_log2FC = log2FC_grid,
      count = count,
      Cluster1 = sbcl,
      Color = subclass_colors[i]
    )

    # Combine with existing data
    cdf_data <- rbind(cdf_data, shared_gene_counts)
    
  }
  
  cdf_data$Cluster1 <- factor(cdf_data$Cluster1, levels = subclass.order[subclass.order %in% cdf_data$Cluster1])
  # Plot cumulative distribution function
  ggplot(cdf_data, aes(x = avg_log2FC, y = count, color = Cluster1, group = Cluster1)) +
    geom_line(size = 1) + # Use line size
    scale_color_manual(values = setNames(subclass_colors, subclass.order), guide = guide_legend(reverse = TRUE)) + # Set the colors using subclass_colors
    labs(title = paste0(""),
         x = "avg_log2FC",
         y = "Number of DE Genes",
         color = "Cluster Pair") +
    theme_minimal() +
    scale_x_continuous(limits = c(0.2, 2)) + # Set x-axis limits
    geom_vline(xintercept = 0.2, linetype = "dotted", color = "black") +
    annotate("text", x = 0.2, y = max(cdf_data$count) * 0.9, label = "0.2", hjust = -0.2, vjust = -0.5) +
    geom_vline(xintercept = 0.5, linetype = "dotted", color = "black") +
    annotate("text", x = 0.5, y = max(cdf_data$count) * 0.9, label = "0.5", hjust = -0.2, vjust = -0.5) +
    geom_vline(xintercept = 1, linetype = "dotted", color = "black") +
    annotate("text", x = 1, y = max(cdf_data$count) * 0.9, label = "1.0", hjust = -0.2, vjust = -0.5) +
    theme(
      axis.title.x = element_text(size = 14), # Increase x-axis label size
      axis.title.y = element_text(size = 14), # Increase y-axis label size
      axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
      axis.text.y = element_text(size = 12)   # Increase y-axis tick label size
    )
}

PlotSubclassGeneCountCDF_AtPoints <- function(de_df, all_genes, sample.name, subclass.order, subclass_colors, log2FC_points) {
  # Extract dataframes
  df <- de_df$subclass
  
  # Initialize an empty data frame to store cumulative counts
  cdf_data <- data.frame()
  
  for (i in 1:length(subclass.order)) {
    
    sbcl <- subclass.order[i]
    
    # Filter dataframes for the specific cluster pair
    df_filtered <- df %>% filter(cluster == sbcl)
    
    count <- c()
    for (l in log2FC_points) {
      de_genes <- df_filtered$gene[df_filtered$avg_log2FC > l]
      count <- c(count, length(de_genes))
    }
    
    # Calculate the cumulative number of shared genes as avg_log2FC varies
    shared_gene_counts <<- data.frame(
      avg_log2FC = log2FC_points,
      count = count,
      Cluster1 = sbcl,
      Color = subclass_colors[i]
    )
    
    # Combine with existing data
    cdf_data <- rbind(cdf_data, shared_gene_counts)
    
  }
  
  # Plot cumulative distribution function
  ggplot(cdf_data, aes(x = avg_log2FC, y = count, color = Cluster1, group = Cluster1)) +
    geom_line(size = 1) + # Use line size
    geom_point(size = 2) + # Add points at specified log2FC values
    scale_color_manual(values = setNames(subclass_colors, subclass.order)) + # Set the colors using subclass_colors
    labs(title = paste0(""),
         x = "avg_log2FC",
         y = "Cumulative Number of Shared Genes",
         color = "Cluster Pair") +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14), # Increase x-axis label size
      axis.title.y = element_text(size = 14), # Increase y-axis label size
      axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
      axis.text.y = element_text(size = 12)   # Increase y-axis tick label size
    )
}

ShuffleExpression <- function(seurat_obj, metadata_column, assay = "RNA", ignore_group = FALSE) {
  # Check if the metadata_column exists in the Seurat object metadata
  if (!(metadata_column %in% colnames(seurat_obj@meta.data))) {
    stop(paste("Column", metadata_column, "not found in Seurat object metadata"))
  }
  
  # Create a copy of the expression data from the counts slot
  expr_data <- as.matrix(GetAssayData(seurat_obj, assay = assay, slot = "counts"))
  
  DefaultAssay(seurat_obj) <- assay
  Idents(seurat_obj) <- metadata_column
  
  if (ignore_group) {
    # Shuffle expression values across all cells
    shuffled_expr <- t(apply(expr_data, 1, function(x) x[sample(length(x))]))
    colnames(shuffled_expr) <- colnames(expr_data)
    
    # Create a new Seurat object to ensure the original remains unchanged
    seurat_obj_copy <- seurat_obj
    seurat_obj_copy <- SetAssayData(seurat_obj_copy, assay = assay, slot = "counts", new.data = as(shuffled_expr, "dgCMatrix"))
    
    return(seurat_obj_copy)
  } else {
    # Get unique categories
    categories <- unique(seurat_obj@meta.data[[metadata_column]])
    
    # Initialize a list to store shuffled expression matrices
    shuffled_expr_list <- list()
    cell_names_list <- list()
    
    # Loop through each category and shuffle expression values within that category
    for (category in categories) {
      print(category)
      
      # Get the cells belonging to the current category
      category_cells <- WhichCells(seurat_obj, idents = category)
      
      # Get the expression matrix for the current category cells
      category_expr <- expr_data[, category_cells]
      
      # Shuffle expression values for each gene within the current category
      shuffled_expr <- t(apply(category_expr, 1, function(x) x[sample(length(x))]))
      
      # Add the shuffled expression matrix to the list
      shuffled_expr_list[[category]] <- shuffled_expr
      cell_names_list[[category]] <- colnames(category_expr)
    }
    
    # Concatenate the shuffled expression matrices
    shuffled_expr_concat <- do.call(cbind, shuffled_expr_list)
    
    # Combine the cell names to maintain order
    all_cell_names <- unlist(cell_names_list)
    
    # Reorder the columns to match the original cell names
    cell_order <- colnames(expr_data)
    shuffled_expr_final <- shuffled_expr_concat[, match(cell_order, all_cell_names)]
    colnames(shuffled_expr_final) <- cell_order
    
    # Create a new Seurat object to ensure the original remains unchanged
    seurat_obj_copy <- seurat_obj
    seurat_obj_copy <- SetAssayData(seurat_obj_copy, assay = assay, slot = "counts", new.data = as(shuffled_expr_final, "dgCMatrix"))
    
    return(seurat_obj_copy)
  }
}

SubsampleObject <- function(seurat_obj, metadata_col, cells_per_category) {
  # Extract metadata
  metadata <- seurat_obj@meta.data
  
  # Get unique values in the specified metadata column
  unique_values <- unique(metadata[[metadata_col]])
  
  # Initialize a list to store subsampled cell names
  subsampled_cells <- list()
  
  # Loop through each unique value in the metadata column
  for (value in unique_values) {
    # Get cells that belong to the current metadata category
    cells_in_category <- rownames(metadata[metadata[[metadata_col]] == value, ])
    
    # Determine the number of cells to sample
    n_cells_to_sample <- min(length(cells_in_category), cells_per_category)
    
    # Sample cells
    sampled_cells <- sample(cells_in_category, n_cells_to_sample)
    
    # Add sampled cells to the list
    subsampled_cells <- c(subsampled_cells, sampled_cells)
  }
  
  # Subset the Seurat object to include only the subsampled cells
  subsampled_seurat_obj <- subset(seurat_obj, cells = as.character(subsampled_cells))
  
  return(subsampled_seurat_obj)
}

SubsampleObjectMultipleIterations <- function(seurat_obj, metadata_col, cells_per_category, iterations) {
  # Extract metadata
  metadata <- seurat_obj@meta.data
  
  # Get unique values in the specified metadata column
  unique_values <- unique(metadata[[metadata_col]])
  
  # Initialize a list to store sampled cells for each iteration
  iterations_list <- vector("list", iterations)
  
  # Initialize a list to keep track of already sampled cells
  sampled_cells_overall <- vector("list", iterations)
  
  for (iter in 1:iterations) {
    # Initialize a vector to store sampled cells for this iteration
    subsampled_cells <- c()
    
    # Initialize a list to track sampled cells for this iteration
    sampled_cells_iter <- list()
    
    for (value in unique_values) {
      # Get cells that belong to the current metadata category
      cells_in_category <- rownames(metadata[metadata[[metadata_col]] == value, ])
      
      # Exclude cells that have already been sampled in previous iterations
      previously_sampled <- unlist(lapply(sampled_cells_overall, function(x) if (!is.null(x[[value]])) x[[value]] else c()))
      available_cells <- setdiff(cells_in_category, subsampled_cells) # Ensure no duplicates within iteration
      
      # Check if there are available cells to sample
      if (length(available_cells) > 0) {
        if (length(available_cells) >= cells_per_category) {
          # Sample without replacement
          sampled_cells <- sample(available_cells, cells_per_category, replace = FALSE)
        } else {
          # Sample all available cells
          sampled_cells <- available_cells
          
          # If needed, sample additional cells from the entire pool
          remaining_needed <- cells_per_category - length(available_cells)
          
          # Sample from the previously used cells plus the remaining available cells
          additional_pool <- setdiff(cells_in_category, sampled_cells)
          
          # Check if additional_pool is non-empty before sampling
          if (length(additional_pool) > 0) {
            additional_cells <- sample(additional_pool, remaining_needed, replace = FALSE)
            # Combine the available and additional sampled cells
            sampled_cells <- c(sampled_cells, additional_cells)
          }
        }
      } else {
        # If no available cells, sample with replacement from the entire category
        sampled_cells <- sample(cells_in_category, cells_per_category, replace = FALSE)
      }
      
      # Add sampled cells to the subsampled_cells vector
      subsampled_cells <- c(subsampled_cells, sampled_cells)
      
      # Store sampled cells for this category in the tracker for this iteration
      sampled_cells_iter[[value]] <- sampled_cells
    }
    
    # Store the subsampled cells for this iteration
    iterations_list[[iter]] <- subsampled_cells
    
    # Store the sampled cells for this iteration in the overall tracker
    sampled_cells_overall[[iter]] <- sampled_cells_iter
  }
  
  return(iterations_list)
}

SplitObjectHalf <- function(seurat_object, seed = 123) {
  # Get the number of cells
  num_cells <- ncol(seurat_object)
  
  # Generate a random split
  set.seed(seed)  # Setting seed for reproducibility
  random_split <- sample(1:num_cells, num_cells, replace = FALSE)
  
  # Split into two groups
  half1 <- random_split[1:(num_cells %/% 2)]
  half2 <- random_split[((num_cells %/% 2) + 1):num_cells]
  
  # Subset the Seurat object into two halves
  seurat_half1 <- subset(seurat_object, cells = colnames(seurat_object)[half1])
  seurat_half2 <- subset(seurat_object, cells = colnames(seurat_object)[half2])
  
  # Return a list containing the two halves
  return(list(obj1 = seurat_half1, obj2 = seurat_half2))
}

MinDistance <- function(seurat_obj, reference_df, pc1_colname = "PC_1", pc2_colname = "PC_2") {
  # Extract PC1 and PC2 coordinates from the Seurat object
  cell_coords <- as.data.frame(Embeddings(seurat_obj, "pca")[, c(pc1_colname, pc2_colname)])
  
  # Function to calculate the Euclidean distance between two points
  euclidean_dist <- function(x1, y1, x2, y2) {
    sqrt((x1 - x2)^2 + (y1 - y2)^2)
  }
  
  # Initialize a vector to store the minimum distances
  min_distances <- numeric(nrow(cell_coords))
  
  # Calculate the minimum distance for each cell
  for (i in 1:nrow(cell_coords)) {
    distances <- apply(reference_df, 1, function(ref_point) {
      euclidean_dist(cell_coords[i, "PC_1"], cell_coords[i, "PC_2"], ref_point["X..X"], ref_point["Y"])
    })
    min_distances[i] <- min(distances)
  }
  
  # Add the distances to the Seurat object's metadata
  seurat_obj <- AddMetaData(seurat_obj, metadata = min_distances, col.name = "min_distance_to_reference")
  
  return(seurat_obj)
}
