
IntegrateObjects <- function(seurat_obj1, seurat_obj2, resolutions = 1, nfeatures = 3000, subsample = FALSE) {
  # Ensure the two objects have the same set of features (genes)
  common_features <- intersect(rownames(seurat_obj1), rownames(seurat_obj2))
  seurat_obj1 <- subset(seurat_obj1, features = common_features)
  seurat_obj2 <- subset(seurat_obj2, features = common_features)
  
  # If subsample is TRUE, randomly subsample the larger object to match the smaller one
  if (!is.logical(subsample)) {
    seurat_obj1 <- seurat_obj1[, sample(colnames(seurat_obj1), subsample)]
    seurat_obj2 <- seurat_obj2[, sample(colnames(seurat_obj2), subsample)]
  } else if (subsample) {
    if (ncol(seurat_obj1) > ncol(seurat_obj2)) {
      seurat_obj1 <- seurat_obj1[, sample(colnames(seurat_obj1), ncol(seurat_obj2))]
    } else if (ncol(seurat_obj2) > ncol(seurat_obj1)) {
      seurat_obj2 <- seurat_obj2[, sample(colnames(seurat_obj2), ncol(seurat_obj1))]
    }
  }
  
  # List of objects to integrate
  seurat_list <- list(seurat_obj1, seurat_obj2)
  
  # Perform SCTransform v2 on each object
  seurat_list <- lapply(seurat_list, function(x) {
    x <- SCTransform(x, vst.flavor = "v2", return.only.var.genes = FALSE, verbose = FALSE) %>%
         RunPCA(npcs = 30, verbose = FALSE) %>%
         RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)
    return(x)
  })
  
  # Select integration features
  integration_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = nfeatures)
  
  # Prepare objects for integration
  seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = integration_features)
  
  # Find integration anchors
  anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = integration_features)
  
  # Integrate data
  integrated_seurat <- IntegrateData(anchorset = anchors, normalization.method = "SCT", features.to.integrate = common_features)
  
  # Run PCA and UMAP on the integrated object
  integrated_seurat <- RunPCA(integrated_seurat, verbose = FALSE) %>%
                       RunUMAP(dims = 1:30, verbose = FALSE) %>%
                       FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE)
  
  for (res in resolutions) {
    integrated_seurat <- FindClusters(integrated_seurat, resolution = res, algorithm = 4, method = "igraph")
  }
  
  return(integrated_seurat)
}

PlotIntegration <- function(obj, integvar, integclust, subclass.order) {
  
  dimplot1 <- DimPlot(obj, reduction = "umap", group.by = integvar, label = FALSE, raster = FALSE, shuffle = TRUE) + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
  print(dimplot1)
  
  for (icl in integclust) {
    
    dimplot1 <- DimPlot(obj, reduction = "umap", split.by = integvar, group.by = icl, label = TRUE, raster = FALSE) + NoLegend() + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
    dimplot2 <- DimPlot(obj, reduction = "umap", split.by = integvar, group.by = "subclass", label = TRUE, raster = FALSE) + NoLegend() + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
    heatmap1 <- IntegratedClusterOverlapHeatmap(obj, integvar, "subclass", icl, subclass.order)
    dimplot3 <- DimPlot(obj, reduction = "umap", split.by = integvar, group.by = "type", label = TRUE, raster = FALSE) + NoLegend() + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
    heatmap2 <- IntegratedClusterOverlapHeatmap(obj, integvar, "type", icl, subclass.order)
    heatmaps3 <- IntegratedClusterMakeupHeatmap(obj, integvar, "type", icl, subclass.order)
    heatmaps4 <- IdentToIntegratedClusterHeatmap(obj, integvar, "type", icl, subclass.order)
    
    print(dimplot1)
    print(dimplot2)
    print(heatmap1)
    print(dimplot3)
    print(heatmap2)
    grid.arrange(grobs = heatmaps3, ncol = 2)
    grid.arrange(grobs = heatmaps4, ncol = 2)
  
  }
}

MapObjects <- function(seurat_obj1, seurat_obj2, idents, assay = "SCT") {
  
  objs <- list(seurat_obj1, seurat_obj2)
  # Perform SCTransform v2 on each object
  if (assay == "integrated") {
    objs <- lapply(objs, function(x) {
      # SCTransform(x, vst.flavor = "v2", return.only.var.genes = FALSE, verbose = FALSE) %>%
      x <- RunPCA(x, npcs = 30, assay = "integrated", verbose = FALSE) %>%
           RunUMAP(reduction = "pca", dims = 1:30, assay = "integrated", return.model = TRUE, verbose = FALSE)
      return(x)
    })
  } else if (assay == "SCT") {
    objs <- lapply(objs, function(x) {
      x <- SCTransform(x, vst.flavor = "v2", return.only.var.genes = FALSE, verbose = FALSE) %>%
           RunPCA(npcs = 30, verbose = FALSE) %>%
           RunUMAP(reduction = "pca", dims = 1:30, return.model = TRUE, verbose = FALSE)
      return(x)
    })
  }
  objs.idx <- list(qry = c(1, 2), ref = c(2, 1))
  objs.mapped <- list()
  
  refdata <- list()
  for (id in idents) { refdata[[id]] <- id }
  
  for (idx in c(1, 2)) {

    # Transfer data from reference to query
    reference <- objs[[objs.idx$ref[idx]]]
    query <- objs[[objs.idx$qry[idx]]]
    anchors.query <- FindTransferAnchors(reference = reference, query = query, reference.reduction = "pca", dims = 1:30)
    # if (nrow(anchors.query@anchors) < 50) { k.weight = 10 } # floor(nrow(anchors.query@anchors) * 0.25)
    # else { k.weight = 50 }
    objs.mapped[[idx]] <- MapQuery(anchorset = anchors.query, 
                                   reference = reference, query = query, 
                                   refdata = refdata, 
                                   reference.reduction = "pca", 
                                   reduction.model = "umap")
  
  }
  return(objs.mapped)
}

PlotMapping <- function(objs, idents = c("subclass", "type"), ident.order = NULL, title.key = "species") {
  
  for (obj.idx in c(1, 2)) {
    obj <- objs[[obj.idx]]
    for (id in idents) {
      id.levels <- as.character(unlist(unique(objs[[setdiff(c(1, 2), obj.idx)]][[id]])))
      for (red in c("umap", "ref.umap")) {
        for (cl in c("%s", "predicted.%s", "predicted.%s.score")) {
          if (cl == "predicted.%s.score") {
            plot <- FeaturePlot(obj, sprintf(cl, id), reduction = red) + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
          }
          else {
            plot <- DimPlot(obj, reduction = red, group.by = sprintf(cl, id), label = TRUE, raster = FALSE) + NoLegend() + xlim(-18, 18) + ylim(-18, 18) + coord_equal()
          }
          print(plot + ggtitle(paste(obj[[title.key]][1,], sprintf(cl, id), "on", red)))
        }
      }
      plot <- PlotMappedLabelsHeatmap(obj, id, id.levels, normalize = "row", ident.order = ident.order)
      print(plot)
      plot <- PlotMappingQualityHeatmap(obj, id, id.levels, sprintf("predicted.%s.score", id), ident.order = ident.order)
      print(plot)
    }
  }
  
}

IntegratedClusterOverlapHeatmap <- function(integrated.obj, integvar.col, ident.col, cluster.col, primary_order,
                                     col.low = "white", col.high = "red", x.lab.rot = TRUE) {
  
  # Extract the relevant columns
  metadata <- integrated.obj@meta.data %>% select(all_of(c(integvar.col, cluster.col, ident.col)))
  
  # Rename columns for convenience
  colnames(metadata) <- c("integvar", "cluster", "ident")
  
  # Get unique integvar levels
  integvar_levels <- unique(metadata$integvar)
  
  # Initialize an empty list to store overlap matrices for each pair of integvar levels
  overlap_matrices <- list()
  
  # Sorting function
  sort_ident <- function(ident, primary_order) {
    primary <- sapply(ident, function(x) str_extract(x, paste(primary_order, collapse = "|")))
    suffix <- sapply(ident, function(x) str_extract(x, "(?<=_)[A-Za-z0-9]+$"))
    suffix_numeric <- suppressWarnings(as.numeric(suffix))
    suffix[is.na(suffix_numeric)] <- paste0("Z", suffix[is.na(suffix_numeric)])  # Add "Z" prefix to non-numeric suffixes to sort them correctly
    suffix_numeric[is.na(suffix_numeric)] <- Inf
    df <- data.frame(ident = ident, primary = primary, suffix = suffix, suffix_numeric = suffix_numeric)
    df <- df %>% arrange(match(primary, primary_order), suffix_numeric, suffix)
    return(df$ident)
  }
  
  # Loop through each pair of integvar levels
  for (i in 1:(length(integvar_levels) - 1)) {
    for (j in (i + 1):length(integvar_levels)) {
      integvar1 <- integvar_levels[i]
      integvar2 <- integvar_levels[j]
      
      # Filter the metadata for the two integvar levels
      data_integvar1 <- metadata %>% filter(integvar == integvar1)
      data_integvar2 <- metadata %>% filter(integvar == integvar2)
      
      # Get unique ident levels for each integvar level
      ident_levels1 <- unique(data_integvar1$ident)
      ident_levels2 <- unique(data_integvar2$ident)
      
      # Sort ident levels
      sorted_ident_levels1 <- sort_ident(ident_levels1, primary_order)
      sorted_ident_levels2 <- sort_ident(ident_levels2, primary_order)
      
      # Get unique clusters
      clusters <- unique(metadata$cluster)
      
      # Initialize the overlap matrix
      overlap_matrix <- matrix(0, nrow = length(sorted_ident_levels1), ncol = length(sorted_ident_levels2), 
                               dimnames = list(rev(sorted_ident_levels1), sorted_ident_levels2))
      
      # Calculate overlap fractions
      for (ident1 in sorted_ident_levels1) {
        for (ident2 in sorted_ident_levels2) {
          for (cluster in clusters) {
            fraction_ident1 <- sum(data_integvar1$ident == ident1 & data_integvar1$cluster == cluster) / sum(data_integvar1$ident == ident1)
            fraction_ident2 <- sum(data_integvar2$ident == ident2 & data_integvar2$cluster == cluster) / sum(data_integvar2$ident == ident2)
            overlap_matrix[ident1, ident2] <- overlap_matrix[ident1, ident2] + min(fraction_ident1, fraction_ident2)
          }
        }
      }
      
      # Convert overlap fractions to percentages
      overlap_matrix <- overlap_matrix * 100
      
      # Store the overlap matrix in the list
      overlap_matrices[[paste(integvar1, integvar2, sep = "_vs_")]] <- overlap_matrix
    }
  }
  
  # Plot the heatmap for each pair of integvar levels
  for (name in names(overlap_matrices)) {
    overlap_matrix <- overlap_matrices[[name]]
    
    # Melt the matrix for ggplot
    melted <- melt(overlap_matrix)
    colnames(melted) <- c("row", "col", "Percentage")
    
    # Create the heatmap plot
    p <- ggplot(melted, aes(y = row, x = col)) + 
      geom_tile(aes(fill = Percentage)) + 
      scale_fill_gradient(low = col.low, high = col.high, limits = c(0, 100)) + 
      geom_text(aes(label = sprintf("%.0f", Percentage)), size = 5) +
      theme_bw() + 
      xlab(integvar_levels[2]) + 
      ylab(integvar_levels[1]) + 
      theme(axis.text.x = element_text(size=16, face="italic", hjust=1, angle = ifelse(x.lab.rot, 90, 0)),
            axis.text.y = element_text(size=16, face="italic"),
            axis.title.x = element_text(size=16),
            axis.title.y = element_text(size=16)) +
      coord_fixed()
    
    # Ensure correct rotation
    if (x.lab.rot) {
      p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    } else {
      p <- p + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
    }
    
    # print(p + ggtitle(paste("Overlap Heatmap:", name, ident.col, "at", cluster.col)))
    p <- p + ggtitle(paste(name, ident.col, "at", cluster.col))
  }
  return(p)
}

IntegratedClusterMakeupHeatmap <- function(integrated.obj, integvar.col, ident.col, cluster.col, primary_order, 
                                           col.low = "white", col.high = "red", x.lab.rot = TRUE) {
  # Extract the relevant columns
  metadata <- integrated.obj@meta.data %>% select(all_of(c(integvar.col, cluster.col, ident.col)))
  
  # Rename columns for convenience
  colnames(metadata) <- c("integvar", "cluster", "ident")
  
  # Get unique integvar levels
  integvar_levels <- unique(metadata$integvar)
  
  # Sorting function
  sort_ident <- function(ident, primary_order) {
    primary <- sapply(ident, function(x) str_extract(x, paste(primary_order, collapse = "|")))
    suffix <- sapply(ident, function(x) str_extract(x, "(?<=_)[A-Za-z0-9]+$"))
    suffix_numeric <- suppressWarnings(as.numeric(suffix))
    suffix[is.na(suffix_numeric)] <- paste0("Z", suffix[is.na(suffix_numeric)])  # Add "Z" prefix to non-numeric suffixes to sort them correctly
    suffix_numeric[is.na(suffix_numeric)] <- Inf
    df <- data.frame(ident = ident, primary = primary, suffix = suffix, suffix_numeric = suffix_numeric)
    df <- df %>% arrange(match(primary, primary_order), suffix_numeric, suffix)
    return(df$ident)
  }
  
  # Initialize a list to store plots
  plot_list <- list()
  cluster_levels <- unique(metadata$cluster)
  cluster_levels <- sort(cluster_levels)
  
  # Plot the fraction of each cluster.col that comes from each ident.col, split by integvar
  for (iv in integvar_levels) {
    data_integvar <- metadata %>% filter(integvar == iv)
    ident_levels <- unique(data_integvar$ident)
    sorted_ident_levels <- sort_ident(ident_levels, primary_order)

    fraction_matrix <- matrix(0, nrow = length(cluster_levels), ncol = length(sorted_ident_levels), 
                              dimnames = list(cluster_levels, sorted_ident_levels))
    
    for (cluster in cluster_levels) {
      for (ident in sorted_ident_levels) {
        if (any(data_integvar$cluster == cluster)) {
        fraction_matrix[cluster, ident] <- sum(data_integvar$cluster == cluster & data_integvar$ident == ident) / sum(metadata$cluster == cluster)
        # fraction_matrix[cluster, ident] <- sum(data_integvar$cluster == cluster & data_integvar$ident == ident) / sum(data_integvar$ident == ident)
        }
        else { fraction_matrix[cluster, ident] <- 0 }
      }
    }
    
    # Melt the matrix for ggplot
    fraction_matrix <- fraction_matrix * 100
    melted_fraction <- melt(fraction_matrix)
    colnames(melted_fraction) <- c("row", "col", "Fraction")
    
    # Create the heatmap plot
    p_fraction <- ggplot(melted_fraction, aes(x = col, y = factor(row, levels = rev(cluster_levels)))) + 
      geom_tile(aes(fill = Fraction)) + 
      scale_fill_gradient(low = col.low, high = col.high, limits = c(0, 100)) + 
      geom_text(aes(label = sprintf("%.0f", Fraction)), size = 3) +
      theme_bw() + 
      xlab(iv) + 
      ylab("") + 
      theme(axis.text.x = element_text(size=12, face="italic", hjust=1, angle = ifelse(x.lab.rot, 90, 0)), 
            axis.text.y = element_text(size=12, face="italic"), 
            axis.title.x = element_text(size=12), 
            axis.title.y = element_text(size=12)) +
      coord_fixed(ratio = length(ident_levels) / length(cluster_levels))
    
    # Ensure correct rotation
    if (x.lab.rot) {
      p_fraction <- p_fraction + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    } else {
      p_fraction <- p_fraction + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
    }
    
    plot_list[[iv]] <- p_fraction + theme(legend.position = "none")
    if (iv == integvar_levels[1]) { plot_list[[iv]] <- plot_list[[iv]] + ggtitle("% Makeup of Integrated Clusters") }
    else { plot_list[[iv]] <- plot_list[[iv]] + ggtitle("") }
  }
  return(plot_list)
}

IdentToIntegratedClusterHeatmap <- function(integrated.obj, integvar.col, ident.col, cluster.col, primary_order, 
                                           col.low = "white", col.high = "red", x.lab.rot = TRUE) {
  # Extract the relevant columns
  metadata <- integrated.obj@meta.data %>% select(all_of(c(integvar.col, cluster.col, ident.col)))
  
  # Rename columns for convenience
  colnames(metadata) <- c("integvar", "cluster", "ident")
  
  # Get unique integvar levels
  integvar_levels <- unique(metadata$integvar)
  
  # Sorting function
  sort_ident <- function(ident, primary_order) {
    primary <- sapply(ident, function(x) str_extract(x, paste(primary_order, collapse = "|")))
    suffix <- sapply(ident, function(x) str_extract(x, "(?<=_)[A-Za-z0-9]+$"))
    suffix_numeric <- suppressWarnings(as.numeric(suffix))
    suffix[is.na(suffix_numeric)] <- paste0("Z", suffix[is.na(suffix_numeric)])  # Add "Z" prefix to non-numeric suffixes to sort them correctly
    suffix_numeric[is.na(suffix_numeric)] <- Inf
    df <- data.frame(ident = ident, primary = primary, suffix = suffix, suffix_numeric = suffix_numeric)
    df <- df %>% arrange(match(primary, primary_order), suffix_numeric, suffix)
    return(df$ident)
  }
  
  # Initialize a list to store plots
  plot_list <- list()
  cluster_levels <- unique(metadata$cluster)
  cluster_levels <- sort(cluster_levels)
  
  # Plot the fraction of each cluster.col that comes from each ident.col, split by integvar
  for (iv in integvar_levels) {
    data_integvar <- metadata %>% filter(integvar == iv)
    ident_levels <- unique(data_integvar$ident)
    sorted_ident_levels <- sort_ident(ident_levels, primary_order)
    
    fraction_matrix <- matrix(0, nrow = length(cluster_levels), ncol = length(sorted_ident_levels), 
                              dimnames = list(cluster_levels, sorted_ident_levels))
    
    for (cluster in cluster_levels) {
      for (ident in sorted_ident_levels) {
        if (any(data_integvar$cluster == cluster)) {
          fraction_matrix[cluster, ident] <- sum(data_integvar$cluster == cluster & data_integvar$ident == ident) / sum(data_integvar$ident == ident)
        }
        else { fraction_matrix[cluster, ident] <- 0 }
      }
    }
    
    # Melt the matrix for ggplot
    fraction_matrix <- fraction_matrix * 100
    melted_fraction <- melt(fraction_matrix)
    colnames(melted_fraction) <- c("row", "col", "Fraction")
    
    # Create the heatmap plot
    p_fraction <- ggplot(melted_fraction, aes(x = col, y = factor(row, levels = rev(cluster_levels)))) + 
      geom_tile(aes(fill = Fraction)) + 
      scale_fill_gradient(low = col.low, high = col.high, limits = c(0, 100)) + 
      geom_text(aes(label = sprintf("%.0f", Fraction)), size = 3) +
      theme_bw() + 
      xlab(iv) + 
      ylab("") + 
      theme(axis.text.x = element_text(size=12, face="italic", hjust=1, angle = ifelse(x.lab.rot, 90, 0)), 
            axis.text.y = element_text(size=12, face="italic"), 
            axis.title.x = element_text(size=12), 
            axis.title.y = element_text(size=12)) +
      coord_fixed(ratio = length(ident_levels) / length(cluster_levels))
    
    # Ensure correct rotation
    if (x.lab.rot) {
      p_fraction <- p_fraction + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    } else {
      p_fraction <- p_fraction + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
    }
    
    plot_list[[iv]] <- p_fraction + theme(legend.position = "none")
    if (iv == integvar_levels[1]) { plot_list[[iv]] <- plot_list[[iv]] + ggtitle(paste("% of", ident.col, "in Integrated Clusters")) }
    else { plot_list[[iv]] <- plot_list[[iv]] + ggtitle("") }
  }
  return(plot_list)
}

PlotIdentDEIntersection <- function(list1, list2, all_genes1, all_genes2, sample.name1, sample.name2, subclass1, subclass2, ident1, ident2, log2FC_threshold, percentage = FALSE) {
  # Extract dataframes
  df1 <- list1[[subclass1]][[ident1]]
  df2 <- list2[[subclass2]][[ident2]]
  
  # Filter by log2FC_threshold
  df1 <- df1 %>% filter(avg_log2FC > log2FC_threshold)
  df2 <- df2 %>% filter(avg_log2FC > log2FC_threshold)
  
  # Find intersecting genes
  intersecting_genes <- intersect(all_genes1, all_genes2)
  
  # Create grid
  cluster_combinations <- expand.grid(
    Cluster1 = unique(df1$cluster), 
    Cluster2 = unique(df2$cluster)
  )
  
  # Initialize an empty data frame to store intersection counts or percentages
  intersection_counts <- data.frame()
  
  # Calculate intersections for each cluster pair
  for (i in 1:nrow(cluster_combinations)) {
    cluster1 <- cluster_combinations$Cluster1[i]
    cluster2 <- cluster_combinations$Cluster2[i]
    
    # Filter dataframes for the specific cluster pair
    df1_filtered <- df1 %>% filter(cluster == cluster1)
    df2_filtered <- df2 %>% filter(cluster == cluster2)
    
    # Find intersecting DE genes for the specific cluster pair
    intersecting_de_genes <- intersect(df1_filtered$gene, df2_filtered$gene)
    
    # Count the number of intersecting DE genes
    count <- length(intersecting_de_genes)
    
    if (percentage) {
      # Calculate total intersecting DE genes for percentage calculation
      # total_intersecting_genes <<- length(intersect(df1_filtered$gene, intersecting_genes)) + length(intersect(df2_filtered$gene, intersecting_genes))
      total_intersecting_genes <- length(union(intersect(df1_filtered$gene, intersecting_genes), intersect(df2_filtered$gene, intersecting_genes)))
      # Calculate percentage
      percent <- (count / total_intersecting_genes) * 100
      intersection_counts <- rbind(intersection_counts, data.frame(Cluster1 = cluster1, Cluster2 = cluster2, Value = percent))
    } else {
      intersection_counts <- rbind(intersection_counts, data.frame(Cluster1 = cluster1, Cluster2 = cluster2, Value = count))
    }
  }
  
  fill_label <- if (percentage) "Percentage of DE Genes" else "# DE Genes"
  
  intersection_counts$Cluster1 <- factor(intersection_counts$Cluster1, levels = rev(sort_idents(levels(intersection_counts$Cluster1))))
  intersection_counts$Cluster2 <- factor(intersection_counts$Cluster2, levels = rev(sort_idents(levels(intersection_counts$Cluster2))))
  
  # Plot heatmap
  ggplot(intersection_counts, aes(x = Cluster1, y = Cluster2, fill = Value)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = round(Value, 2)), color = "black") +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, max(intersection_counts$Value))) +
    scale_y_discrete(limits = rev(levels(factor(intersection_counts$Cluster2)))) + # Reverse y-axis order
    labs(
      title = paste0("Shared DE Genes (log2FC > ", log2FC_threshold, ")"), 
      x = sample.name1,
      y = sample.name2,
      fill = fill_label
    ) +
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
      axis.text.y = element_text(size = 12)   # Increase y-axis tick label size
    ) # Make cells square
}

PlotSubclassDEIntersectionHeatmap <- function(list1, list2, all_genes1, all_genes2, sample.name1, sample.name2, subclass.order, log2FC_threshold, output_path, percentage = FALSE) {
  # Extract dataframes
  df1 <- list1$subclass
  df2 <- list2$subclass
  
  # Filter by log2FC_threshold
  df1 <- df1 %>% filter(avg_log2FC > log2FC_threshold & cluster %in% subclass.order)
  subclasses <- subclass.order[subclass.order %in% df1$cluster]
  df1$cluster <- factor(df1$cluster, levels = subclasses)
  df2 <- df2 %>% filter(avg_log2FC > log2FC_threshold & cluster %in% subclass.order)
  subclasses <- subclass.order[subclass.order %in% df2$cluster]
  df2$cluster <- factor(df2$cluster, levels = subclasses)
  
  # Find intersecting genes
  intersecting_genes <- intersect(all_genes1, all_genes2)
  
  # Create grid
  cluster_combinations <<- expand.grid(
    Cluster1 = unique(df1$cluster), 
    Cluster2 = unique(df2$cluster)
  )
  
  # Initialize an empty data frame to store intersection counts or percentages
  intersection_counts <- data.frame()
  
  # Calculate intersections for each cluster pair
  for (i in 1:nrow(cluster_combinations)) {
    cluster1 <- cluster_combinations$Cluster1[i]
    cluster2 <- cluster_combinations$Cluster2[i]
    
    # Filter dataframes for the specific cluster pair
    df1_filtered <- df1 %>% filter(cluster == cluster1)
    df2_filtered <- df2 %>% filter(cluster == cluster2)
    
    # Find intersecting DE genes for the specific cluster pair
    intersecting_de_genes <- intersect(df1_filtered$gene, df2_filtered$gene)
    all_de_genes <- union(df1_filtered$gene[df1_filtered$gene %in% intersecting_genes], df2_filtered$gene[df2_filtered$gene %in% intersecting_genes])
    
    # Count the number of intersecting DE genes
    count <- length(intersecting_de_genes)
    
    if (percentage) {
      # Calculate total intersecting DE genes for percentage calculation
      total_intersecting_genes <- length(union(intersect(df1_filtered$gene, intersecting_genes), intersect(df2_filtered$gene, intersecting_genes)))
      # Calculate percentage
      percent <- (count / total_intersecting_genes) * 100
      intersection_counts <- rbind(intersection_counts, data.frame(Cluster1 = cluster1, Cluster2 = cluster2, Value = percent))
    } else {
      intersection_counts <- rbind(intersection_counts, data.frame(Cluster1 = cluster1, Cluster2 = cluster2, Value = count))
    }
    
    # Save intersecting genes to a text file
    filename <- paste0(output_path, "/", gsub("/", "", cluster1), "_vs_", gsub("/", "", cluster2), "_intersecting_genes.txt")
    write.table(intersecting_de_genes, file = filename, quote = FALSE, row.names = FALSE, col.names = FALSE)
    filename <- paste0(output_path, "/", gsub("/", "", cluster1), "_vs_", gsub("/", "", cluster2), "_all_genes.txt")
    write.table(all_de_genes, file = filename, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  fill_label <- if (percentage) "Percentage of DE Genes" else "# DE Genes"
  
  intersection_counts$Cluster1 <- factor(intersection_counts$Cluster1, levels = rev(sort_by_reference(levels(intersection_counts$Cluster1), subclass.order)))
  intersection_counts$Cluster2 <- factor(intersection_counts$Cluster2, levels = rev(sort_by_reference(levels(intersection_counts$Cluster2), subclass.order)))
  
  # Plot heatmap
  ggplot(intersection_counts, aes(x = Cluster1, y = Cluster2, fill = Value)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = round(Value, 2)), color = "black") +
    scale_fill_gradient(low = "white", high = "red", limits = c(0, max(intersection_counts$Value))) +
    scale_y_discrete(limits = rev(levels(factor(intersection_counts$Cluster2)))) + # Reverse y-axis order
    labs(
      title = paste0("Shared DE Genes (log2FC > ", log2FC_threshold, ")"), 
      x = sample.name1,
      y = sample.name2,
      fill = fill_label
    ) +
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),  # Rotate x-axis tick labels
      axis.text.y = element_text(size = 12)   # Increase y-axis tick label size
    ) # Make cells square
}

PlotSubclassDEIntersectionTotalCDF <- function(list, ortho_genes, subclass.order, subclass_colors) {
  # Extract dataframes
  df <- list$subclass
  
  # Initialize an empty data frame to store cumulative counts
  cdf_data <- data.frame()
  
  # Calculate cumulative distribution for each specified cluster pair
  for (i in seq_along(subclass.order)) {
    sbcl <- subclass.order[[i]]
    
    # Filter dataframes for the specific cluster pair
    df_filtered <- df %>% filter(cluster == sbcl)
    
    # Calculate the maximum avg_log2FC value
    max_log2FC <- max(df_filtered$avg_log2FC, na.rm = TRUE)
    
    # Create a sequence of avg_log2FC values from 0.2 to the maximum or 2 (whichever is smaller)
    log2FC_grid <- seq(0.2, min(max_log2FC, 2), length.out = 50)
    
    count <- c()
    for (l in log2FC_grid) {
      de_genes <- df_filtered$gene[df_filtered$avg_log2FC > l]
      de_genes_ortho <- de_genes[de_genes %in% ortho_genes]
      count <- c(count, length(de_genes_ortho) / length(de_genes))
    }
    
    # Calculate the cumulative number of shared genes as avg_log2FC varies
    shared_gene_counts <- data.frame(
      avg_log2FC = log2FC_grid,
      count = count,
      Cluster = sbcl,
      Color = subclass_colors[i]
    )
    
    # Combine with existing data
    cdf_data <- rbind(cdf_data, shared_gene_counts)
    
  }
  
  cdf_data$Cluster <- factor(cdf_data$Cluster, levels = subclass.order)
  # Plot cumulative distribution function
  ggplot(cdf_data, aes(x = avg_log2FC, y = count, color = Cluster, group = Cluster)) +
    geom_line(size = 1) + # Use line size
    scale_color_manual(values = setNames(subclass_colors, subclass.order), guide = guide_legend(reverse = TRUE)) + # Set the colors using subclass_colors
    labs(title = paste0(""),
         x = "avg_log2FC",
         y = "Fraction of DE Genes with a 1:1 Ortholog",
         color = "Cluster Pair") +
    theme_minimal() +
    scale_x_continuous(limits = c(0.2, 2)) + # Set x-axis limits
    scale_y_continuous(limits = c(0, 1)) + # Set x-axis limits
    geom_vline(xintercept = 0.2, linetype = "dotted", color = "black") +
    annotate("text", x = 0.2, y = 0.95, label = "0.2", hjust = -0.2, vjust = -0.5) +
    geom_vline(xintercept = 0.5, linetype = "dotted", color = "black") +
    annotate("text", x = 0.5, y = 0.95, label = "0.5", hjust = -0.2, vjust = -0.5) +
    geom_vline(xintercept = 1, linetype = "dotted", color = "black") +
    annotate("text", x = 1, y = 0.95, label = "1", hjust = -0.2, vjust = -0.5) +
    theme(
      axis.title.x = element_text(size = 14), # Increase x-axis label size
      axis.title.y = element_text(size = 14), # Increase y-axis label size
      axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
      axis.text.y = element_text(size = 12)   # Increase y-axis tick label size
    )
}

PlotSubclassDEIntersectionCDF <- function(list1, list2, all_genes1, all_genes2, sample.name1, sample.name2, subclass.order, output_path, cluster_pairs, pair_colors, normalize.within = TRUE) {
  # Extract dataframes
  df1 <- list1$subclass
  df2 <- list2$subclass
  
  # Find intersecting genes
  intersecting_genes <- intersect(all_genes1, all_genes2)
  
  # Initialize an empty data frame to store cumulative counts
  cdf_data <- data.frame()
  
  # Calculate cumulative distribution for each specified cluster pair
  for (i in seq_along(cluster_pairs)) {
    pair <- cluster_pairs[[i]]
    cluster1 <- pair[1]
    cluster2 <- pair[2]
    
    # Filter dataframes for the specific cluster pair
    df1_filtered <- df1 %>% filter(cluster == cluster1)
    df2_filtered <- df2 %>% filter(cluster == cluster2)
    
    # Calculate the maximum avg_log2FC value
    max_log2FC <- max(c(df1_filtered$avg_log2FC, df2_filtered$avg_log2FC), na.rm = TRUE)
    
    # Create a sequence of avg_log2FC values from 0.2 to the maximum or 2 (whichever is smaller)
    log2FC_grid <- seq(0.2, min(max_log2FC, 2), length.out = 50)
    
    total_intersecting_genes_init <- length(union(intersect(df1_filtered$gene, intersecting_genes), intersect(df2_filtered$gene, intersecting_genes)))
    
    count <- c()
    for (l in log2FC_grid) {
      de_genes_1 <- df1_filtered$gene[df1_filtered$avg_log2FC > l]
      de_genes_2 <- df2_filtered$gene[df2_filtered$avg_log2FC > l]
      intersecting_de_genes <- intersect(de_genes_1, de_genes_2)
      total_intersecting_genes <- length(union(intersect(de_genes_1, intersecting_genes), intersect(de_genes_2, intersecting_genes)))
      if (normalize.within == TRUE) { 
        count <- c(count, length(intersecting_de_genes) / total_intersecting_genes) }
      else { count <- c(count, length(intersecting_de_genes) / total_intersecting_genes_init) }
      
    }
    
    # Calculate the cumulative number of shared genes as avg_log2FC varies
    shared_gene_counts <<- data.frame(
      avg_log2FC = log2FC_grid,
      count = count,
      Cluster1 = cluster1,
      Cluster2 = cluster2,
      Color = pair_colors[i]
    )
    
    # Combine with existing data
    cdf_data <- rbind(cdf_data, shared_gene_counts)
    
    # # Save intersecting genes to a text file
    # filename <- paste0(output_path, "/", gsub("/", "", cluster1), "_vs_", gsub("/", "", cluster2), "_intersecting_genes.txt")
    # write.table(intersecting_de_genes, file = filename, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  cdf_data$Cluster1 <- factor(cdf_data$Cluster1, levels = subclass.order[subclass.order %in% cdf_data$Cluster1])
  cdf_data$Cluster2 <- factor(cdf_data$Cluster2, levels = subclass.order[subclass.order %in% cdf_data$Cluster2])
  # Plot cumulative distribution function
  ggplot(cdf_data, aes(x = avg_log2FC, y = count, color = interaction(Cluster1, Cluster2), group = interaction(Cluster1, Cluster2))) +
    geom_line(size = 1) + # Use line size
    scale_color_manual(values = setNames(pair_colors, unique(interaction(cdf_data$Cluster1, cdf_data$Cluster2))), guide = guide_legend(reverse = TRUE)) + # Set the colors using pair_colors
    labs(title = paste0(""),
         x = "avg_log2FC",
         y = "Fraction of 1:1 DE Genes Shared Across Species",
         color = "Cluster Pair") +
    theme_minimal() +
    scale_x_continuous(limits = c(0.2, 2)) + # Set x-axis limits
    geom_vline(xintercept = 0.2, linetype = "dotted", color = "black") +
    annotate("text", x = 0.2, y = 0.5, label = "0.2", hjust = -0.2, vjust = -0.5) +
    geom_vline(xintercept = 0.5, linetype = "dotted", color = "black") +
    annotate("text", x = 0.5, y = 0.5, label = "0.5", hjust = -0.2, vjust = -0.5) +
    geom_vline(xintercept = 1, linetype = "dotted", color = "black") +
    annotate("text", x = 1, y = 0.5, label = "1", hjust = -0.2, vjust = -0.5) +
    theme(
      axis.title.x = element_text(size = 14), # Increase x-axis label size
      axis.title.y = element_text(size = 14), # Increase y-axis label size
      axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
      axis.text.y = element_text(size = 12)   # Increase y-axis tick label size
    )
}

PlotSubclassDEIntersectionCDF_TopX <- function(list1, list2, all_genes1, all_genes2, sample.name1, sample.name2, subclass.order, output_path, cluster_pairs, pair_colors, top_genes_seq) {
  # Extract dataframes
  df1 <- list1$subclass
  df2 <- list2$subclass
  
  # Find intersecting genes
  intersecting_genes <- intersect(all_genes1, all_genes2)
  
  # Initialize an empty data frame to store cumulative counts
  cdf_data <- data.frame()
  
  # Calculate cumulative distribution for each specified cluster pair
  for (i in seq_along(cluster_pairs)) {
    pair <- cluster_pairs[[i]]
    cluster1 <- pair[1]
    cluster2 <- pair[2]
    
    # Filter dataframes for the specific cluster pair
    df1_filtered <- df1 %>% filter(cluster == cluster1) %>% filter(gene %in% intersecting_genes)
    df2_filtered <- df2 %>% filter(cluster == cluster2) %>% filter(gene %in% intersecting_genes)
    
    # Rank genes by avg_log2FC
    df1_filtered <- df1_filtered %>% arrange(desc(avg_log2FC))
    df2_filtered <- df2_filtered %>% arrange(desc(avg_log2FC))
    
    # Calculate the cumulative fraction of shared genes as top X DE genes vary
    shared_gene_counts <- data.frame(
      top_genes = top_genes_seq,
      count = sapply(top_genes_seq, function(x) {
        top_genes1 <- head(df1_filtered$gene, x)
        top_genes2 <- head(df2_filtered$gene, x)
        length(intersect(top_genes1, top_genes2)) / x
      }),
      Cluster1 = cluster1,
      Cluster2 = cluster2,
      Color = pair_colors[i]
    )
    
    # Combine with existing data
    cdf_data <- rbind(cdf_data, shared_gene_counts)
    
  }
  
  cdf_data$Cluster1 <- factor(cdf_data$Cluster1, levels = subclass.order[subclass.order %in% cdf_data$Cluster1])
  cdf_data$Cluster2 <- factor(cdf_data$Cluster2, levels = subclass.order[subclass.order %in% cdf_data$Cluster2])
  cdf_data$ClusterPair <- interaction(cdf_data$Cluster1, cdf_data$Cluster2)
  # Plot cumulative distribution function
  ggplot(cdf_data, aes(x = top_genes, y = count, color = ClusterPair, group = ClusterPair)) +
    geom_line(size = 1) + # Use line size
    scale_color_manual(values = setNames(pair_colors, unique(cdf_data$ClusterPair)), guide = guide_legend(reverse = TRUE)) + # Set the colors using pair_colors
    labs(title = paste0(""),
         x = "Top X DE Genes",
         y = "Fraction of 1:1 DE Genes Shared Across Species",
         color = "Cluster Pair") +
    scale_x_reverse() + # Reverse the x-axis
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14), # Increase x-axis label size
      axis.title.y = element_text(size = 14), # Increase y-axis label size
      axis.text.x = element_text(size = 12),  # Increase x-axis tick label size
      axis.text.y = element_text(size = 12)   # Increase y-axis tick label size
    )
}

PlotIdentCrossConfusionMatrices <- function(obj1, obj2, sample.name1, sample.name2, assay = "SCT", subclass.labels, ident.labels, n_iters = 10, all.genes = FALSE, ident.genes = FALSE, upsample = FALSE, downsample = FALSE) {
  library(gridExtra)
  
  DefaultAssay(obj1) <- assay
  DefaultAssay(obj2) <- assay
  
  objs.list <- list(obj1, obj2)
  if (all.genes) {
    variable.features = intersect(rownames(obj1), rownames(obj2))
  }
  else {
    variable.features <<- SelectIntegrationFeatures(objs.list, nfeatures = 3000)
  }
  
  all_plots <- list()
  
  for (sbcl in subclass.labels) {
    
    all_plots[[sbcl]] <- list()
    
    for (id in ident.labels) {
      
      all_plots[[sbcl]][[id]] <- list(avg1 = NA, avg2 = NA, grid1 = NA, grid2 = NA)
      
      id.num <- strsplit(id, ".", fixed = TRUE)[[1]]
      if (length(id.num) == 2) { sbcl.col <- paste0("subclass.", id.num[2]) }
      else if (length(id.num) == 3) { sbcl.col <- paste0("subclass.", id.num[2], ".", id.num[3]) }
      else { sbcl.col <- paste0("subclass.", id.num[1]) }
      
      Idents(obj1) <- sbcl.col
      Idents(obj2) <- sbcl.col
      
      if (sbcl %in% levels(obj1) & sbcl %in% levels(obj2)) {
        
        obj1.sbcl.id <- subset(obj1, idents = sbcl)
        obj2.sbcl.id <- subset(obj2, idents = sbcl)
        
        if (ident.genes == TRUE && all.genes == FALSE && assay == "SCT") {
          print("Using IntegrationFeatures from the Ident level...")
          obj1.sbcl.id.temp <- SCTransform(obj1.sbcl.id, verbose = FALSE)
          obj2.sbcl.id.temp <- SCTransform(obj2.sbcl.id, verbose = FALSE)
          objs.list <- list(obj1.sbcl.id.temp, obj2.sbcl.id.temp)
          variable.features <- SelectIntegrationFeatures(objs.list, nfeatures = 3000)
        }
        else if (ident.genes == TRUE && all.genes == FALSE && assay == "integrated") {
          print("Using IntegrationFeatures from the Ident level...")
          obj1.sbcl.id.temp <- FindVariableFeatures(obj1.sbcl.id, assay = "integrated", nfeatures = 3000)
          obj2.sbcl.id.temp <- FindVariableFeatures(obj2.sbcl.id, assay = "integrated", nfeatures = 3000)
          objs.list <- list(obj1.sbcl.id.temp, obj2.sbcl.id.temp)
          variable.features <- SelectIntegrationFeatures(objs.list, nfeatures = 3000)
        }
        
        DefaultAssay(obj1.sbcl.id) <- assay
        DefaultAssay(obj2.sbcl.id) <- assay
        
        Idents(obj1.sbcl.id) <- id
        Idents(obj2.sbcl.id) <- id
        
        levels(obj1.sbcl.id) <- sort_idents(levels(obj1.sbcl.id))
        levels(obj2.sbcl.id) <- sort_idents(levels(obj2.sbcl.id))
        
        if (length(levels(obj1.sbcl.id)) > 1 & length(levels(obj2.sbcl.id)) > 1) {
          
          confusion_matrices1 <- list()
          confusion_matrices2 <- list()

          names1 <- unique(levels(Idents(obj1.sbcl.id)))
          names2 <- unique(levels(Idents(obj2.sbcl.id)))
          
          # Initialize matrices to store summed percentages
          confusion_matrix1_sum <- matrix(0, nrow = length(names2), ncol = length(names1),
                                          dimnames = list(names2, names1))
          confusion_matrix2_sum <- matrix(0, nrow = length(names1), ncol = length(names2),
                                          dimnames = list(names1, names2))
          
          valid_matrices1 <- 0
          valid_matrices2 <- 0
          
          for (i in 1:n_iters) {
            
            # Helper function to compute confusion matrix
            compute_confusion_matrix <- function(train_obj, test_obj) {
              mdl <- TrainModel(train_obj, assay = assay, training_genes = variable.features, upsample = upsample, downsample = downsample)
              confusion_matrix <- BuildConfusionMatrix(test_obj, train_obj, model = mdl, assay = assay)
              as.table(confusion_matrix)
            }
            
            confusion_matrix1 <- compute_confusion_matrix(obj1.sbcl.id, obj2.sbcl.id)
            confusion_matrix2 <- compute_confusion_matrix(obj2.sbcl.id, obj1.sbcl.id)
            
            if (!is.null(confusion_matrix1)) {
              cm1 <- as.matrix(confusion_matrix1)
              if (nrow(cm1) > 0 && ncol(cm1) > 0) {
                cm1_aligned <- align_matrix(cm1, names2, names1)
                confusion_matrix1_sum <- confusion_matrix1_sum + cm1_aligned
                confusion_matrices1[[valid_matrices1 + 1]] <- cm1_aligned
                valid_matrices1 <- valid_matrices1 + 1
              }
            }
            
            if (!is.null(confusion_matrix2)) {
              cm2 <- as.matrix(confusion_matrix2)
              if (nrow(cm2) > 0 && ncol(cm2) > 0) {
                cm2_aligned <- align_matrix(cm2, names1, names2)
                confusion_matrix2_sum <- confusion_matrix2_sum + cm2_aligned
                confusion_matrices2[[valid_matrices2 + 1]] <- cm2_aligned
                valid_matrices2 <- valid_matrices2 + 1
              }
            }
          }
          
          # Average the confusion matrices if there are valid matrices
          if (valid_matrices1 > 0) {
            confusion_matrix1_avg <- confusion_matrix1_sum / valid_matrices1
            avg_plot1 <- create_ident_confusion_matrix_plot(confusion_matrix1_avg, sample.name1, sample.name2)
            avg_plot1 <- avg_plot1 + ggtitle(paste0(length(variable.features), " Features"))
            grid_plot1 <- plot_ident_individual_confusion_matrices(confusion_matrices1, valid_matrices1, sample.name1, sample.name2)
            all_plots[[sbcl]][[id]][["avg1"]] <- avg_plot1
            all_plots[[sbcl]][[id]][["grid1"]] <- grid_plot1
          }
          
          if (valid_matrices2 > 0) {
            confusion_matrix2_avg <- confusion_matrix2_sum / valid_matrices2
            avg_plot2 <- create_ident_confusion_matrix_plot(confusion_matrix2_avg, sample.name2, sample.name1)
            avg_plot2 <- avg_plot2 + ggtitle(paste0(length(variable.features), " Features"))
            grid_plot2 <- plot_ident_individual_confusion_matrices(confusion_matrices2, valid_matrices2, sample.name2, sample.name1)
            all_plots[[sbcl]][[id]][["avg2"]] <- avg_plot2
            all_plots[[sbcl]][[id]][["grid2"]] <- grid_plot2
          }
        }
      }
    }
  }
  
  return(all_plots)
}

align_matrix <- function(mat, row_names, col_names) {
  aligned_mat <- matrix(0, nrow = length(row_names), ncol = length(col_names),
                        dimnames = list(row_names, col_names))
  mat_rownames <- rownames(mat)
  mat_colnames <- colnames(mat)
  for (r in mat_rownames) {
    for (c in mat_colnames) {
      aligned_mat[r, c] <- mat[r, c]
    }
  }
  return(aligned_mat)
}

sort_idents <- function(vec) {
  suffixes <- sub(".*_", "", vec)
  
  # Determine if suffixes are numeric or alphabetic
  if (all(grepl("^[0-9]+$", suffixes))) {
    sorted_vec <- vec[order(as.numeric(suffixes))]
  } else {
    sorted_vec <- vec[order(suffixes)]
  }
  
  return(factor(rev(sorted_vec), levels = rev(sorted_vec)))
}

create_ident_confusion_matrix_plot <- function(confusion_matrix, sample.name1, sample.name2) {
  row.levels <- sort_idents(rownames(confusion_matrix))
  col.levels <- sort_idents(colnames(confusion_matrix))
  
  if (!is.null(confusion_matrix)) {
    p <- plotConfusionMatrix(confusion_matrix, plot.return = TRUE, stagger.threshold = 15, 
                             row.levels = row.levels, col.levels = col.levels)
    confusion_plot <- p +
      coord_equal() +
      labs(
        x = paste0("Predicted (Train = ", sample.name1, ")"),
        y = paste0("True (Test = ", sample.name2, ")")
      ) +
      theme(axis.text.x = element_text(angle = 90),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA))
    return(confusion_plot)
  }
}

plot_ident_individual_confusion_matrices <- function(confusion_matrices, n_valid_iters, sample.name1, sample.name2) {
  plot_list <- list()
  
  for (i in 1:n_valid_iters) {
    confusion_matrix <- confusion_matrices[[i]]
    row.levels <- sort_idents(rownames(confusion_matrix))
    col.levels <- sort_idents(colnames(confusion_matrix))
    
    if (!is.null(confusion_matrix)) {
      p <- plotConfusionMatrix(confusion_matrix, plot.return = TRUE, stagger.threshold = 15, 
                               row.levels = row.levels, col.levels = col.levels)
      confusion_plot <- p +
        coord_equal() +
        labs(
          title = paste0("Iteration ", i)
        ) +
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_rect(fill = "white", color = NA),
              plot.background = element_rect(fill = "white", color = NA),
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
      plot_list[[i]] <- confusion_plot
    }
  }
  
  grid_plot <- grid.arrange(grobs = plot_list, ncol = 5)
  return(grid_plot)
}

PlotSubclassCrossConfusionMatrices <- function(obj1, obj2, sample.name1, sample.name2, assay = "SCT", subclass.order, n_iters = 10, all.genes = FALSE, upsample = FALSE, downsample = FALSE) {
  library(gridExtra)
  
  DefaultAssay(obj1) <- assay
  DefaultAssay(obj2) <- assay
  Idents(obj1) <- "subclass"
  Idents(obj2) <- "subclass"
  levels(obj1) <- sort_by_reference(levels(obj1), subclass.order)
  levels(obj2) <- sort_by_reference(levels(obj2), subclass.order)
  objs.list <- list(obj1, obj2)
  
  if (all.genes) {
    variable.features = intersect(rownames(obj1), rownames(obj2))
  }
  else {
    variable.features <<- SelectIntegrationFeatures(objs.list, nfeatures = 3000)
  }
  
  # Helper function to compute confusion matrix
  compute_confusion_matrix <- function(train_obj, test_obj) {
    mdl <- TrainModel(train_obj, assay = assay, training_genes = variable.features, upsample = upsample, downsample = downsample)
    confusion_matrix <- BuildConfusionMatrix(test_obj, train_obj, model = mdl, assay = assay)
    as.table(confusion_matrix)
  }
  
  # Initialize lists to store all row and column names
  all_row_names <- unique(c(levels(Idents(obj1)), levels(Idents(obj2))))
  all_col_names <- unique(c(levels(Idents(obj1)), levels(Idents(obj2))))
  
  # Function to align matrices to have the same row and column names
  align_matrix <- function(mat, row_names, col_names) {
    aligned_mat <- matrix(0, nrow = length(row_names), ncol = length(col_names),
                          dimnames = list(row_names, col_names))
    mat_rownames <- rownames(mat)
    mat_colnames <- colnames(mat)
    for (r in mat_rownames) {
      for (c in mat_colnames) {
        aligned_mat[r, c] <- mat[r, c]
      }
    }
    return(aligned_mat)
  }
  
  # Function to filter rows and columns with non-zero sums
  filter_non_zero <- function(mat) {
    row_sums <- rowSums(mat)
    col_sums <- colSums(mat)
    non_zero_rows <- rownames(mat)[row_sums != 0]
    non_zero_cols <- colnames(mat)[col_sums != 0]
    mat[non_zero_rows, non_zero_cols]
  }
  
  # Initialize matrices to store summed percentages and lists to store individual matrices
  confusion_matrix1_sum <- matrix(0, nrow = length(all_row_names), ncol = length(all_col_names),
                                  dimnames = list(all_row_names, all_col_names))
  confusion_matrix2_sum <- matrix(0, nrow = length(all_row_names), ncol = length(all_col_names),
                                  dimnames = list(all_row_names, all_col_names))
  
  confusion_matrices1 <- list()
  confusion_matrices2 <- list()
  
  # Compute confusion matrices for both directions and store them
  for (i in 1:n_iters) {
    cm1 <- compute_confusion_matrix(objs.list[[1]], objs.list[[2]])
    cm2 <- compute_confusion_matrix(objs.list[[2]], objs.list[[1]])
    cm1_aligned <- align_matrix(cm1, all_row_names, all_col_names)
    cm2_aligned <- align_matrix(cm2, all_row_names, all_col_names)
    confusion_matrix1_sum <- confusion_matrix1_sum + cm1_aligned
    confusion_matrix2_sum <- confusion_matrix2_sum + cm2_aligned
    confusion_matrices1[[i]] <- cm1_aligned
    confusion_matrices2[[i]] <- cm2_aligned
  }
  
  # Average the confusion matrices
  confusion_matrix1_avg <- confusion_matrix1_sum / n_iters
  confusion_matrix2_avg <- confusion_matrix2_sum / n_iters
  
  # Filter rows and columns with non-zero sums
  confusion_matrix1_avg <- filter_non_zero(confusion_matrix1_avg)
  confusion_matrix2_avg <- filter_non_zero(confusion_matrix2_avg)
  
  # Create the averaged confusion matrix plots
  plot1 <- create_confusion_matrix_plot(confusion_matrix1_avg, sample.name1, sample.name2, subclass.order)
  plot1 <- plot1 + ggtitle(paste0(length(variable.features), " Features"))
  plot2 <- create_confusion_matrix_plot(confusion_matrix2_avg, sample.name2, sample.name1, subclass.order)
  plot2 <- plot2 + ggtitle(paste0(length(variable.features), " Features"))
  
  # Create the individual confusion matrix plots
  grid_plot1 <- plot_individual_confusion_matrices(confusion_matrices1, n_iters, sample.name1, sample.name2, subclass.order, filter_non_zero)
  grid_plot2 <- plot_individual_confusion_matrices(confusion_matrices2, n_iters, sample.name2, sample.name1, subclass.order, filter_non_zero)
  
  # Return the plots in a list
  return(list(avg_confusion_plot1 = plot1, avg_confusion_plot2 = plot2, individual_grid_plot1 = grid_plot1, individual_grid_plot2 = grid_plot2))
}

create_confusion_matrix_plot <- function(confusion_matrix, sample.name1, sample.name2, subclass.order) {
  row.levels <- sort_by_reference(rownames(confusion_matrix), subclass.order)
  col.levels <- sort_by_reference(colnames(confusion_matrix), subclass.order)
  
  if (!is.null(confusion_matrix)) {
    p <- plotConfusionMatrix(confusion_matrix, plot.return = TRUE, stagger.threshold = 15, 
                             row.levels = row.levels, col.levels = col.levels)
    confusion_plot <- p +
      coord_equal() +
      labs(
        x = paste0("Predicted (Train = ", sample.name1, ")"),
        y = paste0("True (Test = ", sample.name2, ")")
      ) +
      theme(axis.text.x = element_text(angle = 90),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA))
    return(confusion_plot)
  }
}

plot_individual_confusion_matrices <- function(confusion_matrices, n_iters, sample.name1, sample.name2, subclass.order, filter_func) {
  plot_list <- list()
  
  for (i in 1:n_iters) {
    confusion_matrix <- filter_func(confusion_matrices[[i]])
    row.levels <- sort_by_reference(rownames(confusion_matrix), subclass.order)
    col.levels <- sort_by_reference(colnames(confusion_matrix), subclass.order)
    
    if (!is.null(confusion_matrix)) {
      p <- plotConfusionMatrix(confusion_matrix, plot.return = TRUE, stagger.threshold = 15, 
                               row.levels = row.levels, col.levels = col.levels)
      confusion_plot <- p +
        coord_equal() +
        labs(
          title = paste0("Iteration ", i)
        ) +
        theme(axis.text.x = element_text(angle = 90),
              panel.background = element_rect(fill = "white", color = NA),
              plot.background = element_rect(fill = "white", color = NA),
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
      plot_list[[i]] <- confusion_plot
    }
  }
  
  grid_plot <- grid.arrange(grobs = plot_list, ncol = 5)
  return(grid_plot)
}

PlotMappedLabelsHeatmap <- function(data, column_name, column_levels, normalize = NULL, ident.order = NULL, col.low = "white", col.high = "red", x.lab.rot = TRUE) {
  
  # Create confusion matrix
  confusion_matrix <- table(as.character(unlist(data[[column_name]])), as.character(unlist(data[[paste0("predicted.", column_name)]])))
  confusion_matrix <- as.matrix(confusion_matrix)
  
  add_zeros_to_table <- function(tbl, new_row_names, new_col_names) {
    # Add new rows of zeros
    for (row_name in new_row_names) {
      if (!any(row_name %in% rownames(tbl))) {
        tbl <- rbind(tbl, setNames(t(rep(0, ncol(tbl))), row_name))
        orig_names <- rownames(tbl)[rownames(tbl) != ""]
        rownames(tbl) <- c(orig_names, row_name)
      }
    }
    
    # Add new columns of zeros
    for (col_name in new_col_names) {
      if (!any(col_name %in% colnames(tbl))) {
        tbl <- cbind(tbl, setNames(rep(0, nrow(tbl)), col_name))
        orig_names <- colnames(tbl)[colnames(tbl) != ""]
        colnames(tbl) <- c(orig_names, col_name)
      }
    }
    
    return(tbl)
  }

  # Ensure all possible levels are present in the confusion matrix
  row_levels <- as.character(unlist(unique(data[[column_name]])))
  col_levels <- column_levels
  rows_to_add <- row_levels[row_levels %in% rownames(confusion_matrix) == FALSE]
  cols_to_add <- col_levels[col_levels %in% colnames(confusion_matrix) == FALSE]
  confusion_matrix <- add_zeros_to_table(confusion_matrix, rows_to_add, cols_to_add)
  
  confusion_df <- as.data.frame(confusion_matrix)
  
  # Melt the dataframe for ggplot2
  melted <- melt(confusion_matrix)
  colnames(melted) <- c("row", "col", "Count")
  melted$Count <- as.numeric(melted$Count)
  
  # Normalize if needed
  if (!is.null(normalize)) {
    if (normalize == "row") {
      melted <- ddply(melted, .(row), transform, Percentage = Count / sum(Count) * 100)
    } else if (normalize == "col") {
      melted <- ddply(melted, .(col), transform, Percentage = Count / sum(Count) * 100)
    }
  } else {
    melted$Percentage <- melted$Count
  }
  
  # Sorting function
  sort_ident <- function(ident, primary_order) {
    primary <- sapply(ident, function(x) str_extract(x, paste(primary_order, collapse = "|")))
    suffix <- sapply(ident, function(x) str_extract(x, "(?<=_)[A-Za-z0-9]+$"))
    suffix_numeric <- suppressWarnings(as.numeric(suffix))
    suffix[is.na(suffix_numeric)] <- paste0("Z", suffix[is.na(suffix_numeric)])  # Add "Z" prefix to non-numeric suffixes to sort them correctly
    suffix_numeric[is.na(suffix_numeric)] <- Inf
    df <- data.frame(ident = ident, primary = primary, suffix = suffix, suffix_numeric = suffix_numeric)
    df <- df %>% arrange(match(primary, primary_order), suffix_numeric, suffix)
    return(unlist(df$ident))
  }
  
  # Get unique levels
  row_levels <- unique(melted$row)
  col_levels <- unique(melted$col)
  
  # Sort and factorize identifiers
  if (is.null(ident.order)) {
    row_levels <- sort_ident(row_levels, unique(melted$row))
    col_levels <- sort_ident(col_levels, unique(melted$col))
  } else {
    row_levels <- sort_ident(row_levels, ident.order)
    col_levels <- sort_ident(col_levels, ident.order)
  }
  
  melted$row <- factor(melted$row, levels = row_levels)
  melted$col <- factor(melted$col, levels = col_levels)
  
  if (nrow(confusion_matrix) < 10) { fontsize = 5 }
  else { fontsize = 4 }
  
  # Plot heatmap
  p <- ggplot(melted, aes(y = factor(row, levels = rev(row_levels)), x = factor(col, levels = col_levels))) + 
    geom_tile(aes(fill = Percentage)) + 
    scale_fill_gradient(low = col.low, high = col.high, limits = c(0, max(melted$Percentage))) + 
    geom_text(aes(label = sprintf("%.0f", Percentage)), size = fontsize) +
    theme_bw() + 
    ylab(column_name) + 
    xlab(paste0("predicted_", column_name)) + 
    theme(
      axis.text.x = element_text(size = 16, face = "italic", hjust = 0, angle = ifelse(x.lab.rot, 90, 0)),
      axis.text.y = element_text(size = 16, face = "italic"),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16)
    ) +
    coord_fixed()

  # Ensure correct rotation
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
  }
  
  return(p)
}

PlotMappingQualityHeatmap <- function(data, column_name, column_levels, value_column, ident.order = NULL, col.low = "white", col.high = "red", x.lab.rot = TRUE) {

  # Create average matrix
  average_matrix <<- with(data, tapply(as.numeric(unlist(data[[value_column]])), 
                                       list(as.character(unlist(data[[column_name]])), 
                                       as.character(unlist(data[[paste0("predicted.", column_name)]]))), mean, na.rm = TRUE))
  average_matrix[is.na(average_matrix)] <- 0
  average_matrix <- average_matrix * 100

  add_zeros_to_table <- function(tbl, new_row_names, new_col_names) {
    # Add new rows of zeros
    for (row_name in new_row_names) {
      if (!any(row_name %in% rownames(tbl))) {
        tbl <- rbind(tbl, setNames(t(rep(0, ncol(tbl))), row_name))
        orig_names <- rownames(tbl)[rownames(tbl) != ""]
        rownames(tbl) <- c(orig_names, row_name)
      }
    }

    # Add new columns of zeros
    for (col_name in new_col_names) {
      if (!any(col_name %in% colnames(tbl))) {
        tbl <- cbind(tbl, setNames(rep(0, nrow(tbl)), col_name))
        orig_names <- colnames(tbl)[colnames(tbl) != ""]
        colnames(tbl) <- c(orig_names, col_name)
      }
    }

    return(tbl)
  }

  # Ensure all possible levels are present in the average matrix
  row_levels <<- as.character(unlist(unique(data[[column_name]])))
  col_levels <<- column_levels
  rows_to_add <<- row_levels[row_levels %in% rownames(average_matrix) == FALSE]
  cols_to_add <<- col_levels[col_levels %in% colnames(average_matrix) == FALSE]
  average_matrix <- add_zeros_to_table(average_matrix, rows_to_add, cols_to_add)

  # average_df <- as.data.frame(average_matrix)

  # Melt the dataframe for ggplot2
  melted <- melt(average_matrix)
  colnames(melted) <- c("row", "col", "Value")
  melted$Value <- as.numeric(melted$Value)

  # Sorting function
  sort_ident <- function(ident, primary_order) {
    primary <- sapply(ident, function(x) str_extract(x, paste(primary_order, collapse = "|")))
    suffix <- sapply(ident, function(x) str_extract(x, "(?<=_)[A-Za-z0-9]+$"))
    suffix_numeric <- suppressWarnings(as.numeric(suffix))
    suffix[is.na(suffix_numeric)] <- paste0("Z", suffix[is.na(suffix_numeric)])  # Add "Z" prefix to non-numeric suffixes to sort them correctly
    suffix_numeric[is.na(suffix_numeric)] <- Inf
    df <- data.frame(ident = ident, primary = primary, suffix = suffix, suffix_numeric = suffix_numeric)
    df <- df %>% arrange(match(primary, primary_order), suffix_numeric, suffix)
    return(unlist(df$ident))
  }

  # Get unique levels
  row_levels <- unique(melted$row)
  col_levels <- unique(melted$col)

  # Sort and factorize identifiers
  if (is.null(ident.order)) {
    row_levels <- sort_ident(row_levels, unique(melted$row))
    col_levels <- sort_ident(col_levels, unique(melted$col))
  } else {
    row_levels <- sort_ident(row_levels, ident.order)
    col_levels <- sort_ident(col_levels, ident.order)
  }

  melted$row <- factor(melted$row, levels = row_levels)
  melted$col <- factor(melted$col, levels = col_levels)

  if (nrow(average_matrix) < 10) { fontsize = 5 }
  else { fontsize = 4 }

  # Plot heatmap
  p <- ggplot(melted, aes(y = factor(row, levels = rev(row_levels)), x = factor(col, levels = col_levels))) +
    geom_tile(aes(fill = Value)) +
    scale_fill_gradient(low = col.low, high = col.high, limits = c(0, max(melted$Value, na.rm = TRUE))) +
    geom_text(aes(label = ifelse(Value == 0, "NA", sprintf("%.0f", Value))), size = fontsize) +
    theme_bw() +
    ylab(column_name) +
    xlab(paste0("predicted_", column_name)) +
    theme(
      axis.text.x = element_text(size = 16, face = "italic", hjust = 0, angle = ifelse(x.lab.rot, 90, 0)),
      axis.text.y = element_text(size = 16, face = "italic"),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16)
    ) +
    coord_fixed()

  # Ensure correct rotation
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
  }

  return(p)
}
