
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
  
  obj <- merge(objs[[1]], y = objs[[2]], add.cell.ids = sample_IDs, project = project_name)
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

PlotClusters <- function(obj, group.id) {
  
  if (!missing(group.id)) {
    Idents(obj) <- group.id
  }
  dimplot1 <- DimPlot(obj, reduction = "umap", label = TRUE, raster = FALSE) + NoLegend() + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
  dimplot2 <- DimPlot(obj, reduction = "umap", group.by = "sample", raster = FALSE, shuffle = TRUE) + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
  dimplot3 <- DimPlot(obj, reduction = "umap", group.by = "predicted_doublets", raster = FALSE) + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
  # Summarize doublets by cluster
  df <- obj[[]]
  df$active.ident <- obj@active.ident
  df_summary <- df %>%
    group_by(active.ident, predicted_doublets) %>%
    summarise(count = n()) %>%
    mutate(fraction = count / sum(count))
  # Create the stacked bar plot
  barplot <-ggplot(df_summary, aes(x = active.ident, y = fraction, fill = predicted_doublets)) +
    geom_bar(stat = "identity") +
    labs(x = "Clusters", y = "Doublet Fraction", fill = "Value") +
    theme_minimal()
  featplot1 <- FeaturePlot(obj, "nFeature_RNA", raster = FALSE) + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
  vlnplot1 <- VlnPlot(obj, "nFeature_RNA")
  featplot2 <- FeaturePlot(obj, "nCount_RNA", raster = FALSE) + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
  vlnplot2 <- VlnPlot(obj, "nCount_RNA")
  
  print(dimplot1)
  print(dimplot2)
  print(dimplot3)
  print(barplot)
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

# LabelCells <- function(obj, column, lst) {
#   
#   
#   
#   
# }

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

MarkerDict <- function(obj, subclass.labels, ident.labels) {
  
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
  
  return(marker.dict)
  
}

SaveDotPlots <- function(obj, markers, subclass.labels, ident.labels, savepath, ens.id) {
  
  for (sbcl in subclass.labels) {
    
    folder_path <- paste0(savepath, sbcl, "/")
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
          
          id.path <- paste0(folder_path, id, "/")
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

SaveFeaturePlots <- function(obj, markers, subclass.labels, ident.labels, savepath, ens.id) {
  
  for (sbcl in subclass.labels) {
    
    folder_path <- paste0(savepath, sbcl, "/")
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
          
          id.path <- paste0(folder_path, id, "/")
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
            
            # make feature plots
            features.FC <- unlist(rev(all.markers.FC))
            features.PD <- unlist(rev(all.markers.PD))
            
            for (feature_set in list(features.FC, features.PD)) {
              feature_subset <- feature_set[1:min(20, length(feature_set))]
              feature_labels <- unname(strip_if_contains(feature_subset, ens.id, paste0(ens.id, "000000")))
              
              plots <- FeaturePlot(obj.sbcl.id, features = feature_subset, cols = c("lightgrey", "red"), ncol = 5)
              
              for (j in 1:length(plots)) {
                plots[[j]] <- plots[[j]] + 
                  ggtitle(feature_labels[j]) +
                  coord_equal()
              }
              
              file_prefix <- ifelse(identical(feature_set, features.FC), "FeaturePlot_FC", "FeaturePlot_PD")
              ggsave(paste0(id.path, sbcl.col, "_", file_prefix, "_id-", type, ".png"), plot = plots, width = 24, height = 16, dpi = 300)
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

# Function to generate the plot
plot_gene_counts <- function(nested_list, subclass, clustering_res) {
  # Extract dataframe
  df <- nested_list[[subclass]][[clustering_res]]
  
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
  
  # Plot the number of genes
  ggplot(gene_counts, aes(x = avg_log2FC, y = count, color = cluster, group = cluster)) +
    geom_line(size = 1) + # Thicker lines
    labs(title = paste("Number of DE Genes for", subclass, "at", clustering_res),
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

# Function to generate heatmap
plot_intersection_heatmap <- function(list1, list2, subclass1, subclass2, ident1, ident2) {
  # Extract dataframes
  df1 <- list1[[subclass1]][[ident1]]
  df2 <- list2[[subclass2]][[ident2]]
  
  # Find intersecting genes
  intersecting_genes <<- intersect(df1$gene, df2$gene)
  
  # Filter dataframes to only include intersecting genes
  df1_filtered <<- df1 %>% filter(gene %in% intersecting_genes)
  df2_filtered <<- df2 %>% filter(gene %in% intersecting_genes)
  
  # Create grid
  cluster_combinations <- expand.grid(
    Cluster1 = unique(df1$cluster), 
    Cluster2 = unique(df2$cluster)
  )
  
  # Count intersections for each combo
  intersection_counts <- cluster_combinations %>%
    rowwise() %>%
    mutate(Count = sum(df1_filtered$cluster == Cluster1 & df2_filtered$cluster == Cluster2)) %>%
    ungroup()
  
  # Plot heatmap
  ggplot(intersection_counts, aes(x = Cluster1, y = Cluster2, fill = Count)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = Count), color = "black") +
    scale_fill_gradient(low = "white", high = "red") +
    scale_y_discrete(limits = rev(levels(intersection_counts$Cluster2))) + # Reverse y-axis order
    labs(title = "", # paste("", subclass, "at", clustering_res)
         x = "",
         y = "",
         fill = "# DE Genes") +
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      # axis.title.x = element_text(size = 14), # Increase x-axis label size
      # axis.title.y = element_text(size = 14), # Increase y-axis label size
      axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
      axis.text.y = element_text(size = 12)   # Increase y-axis tick label size
    ) # Make cells square
}
