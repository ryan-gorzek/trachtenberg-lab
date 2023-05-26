
library(scAlign)
library(Seurat)
library(SingleCellExperiment)
library(gridExtra)
library(cowplot)
library(tensorflow)

## User paths
working.dir = "." #where our data file, kowalcyzk_gene_counts.rda is located
results.dir = "." #where the output should be stored

## Load in data
load(url('https://github.com/quon-titative-biology/examples/raw/master/scAlign_paired_alignment/kowalcyzk_gene_counts.rda'))

## Extract age and cell type labels
cell_age = unlist(lapply(strsplit(colnames(C57BL6_mouse_data), "_"), "[[", 1))
cell_type = gsub('HSC', '', unlist(lapply(strsplit(colnames(C57BL6_mouse_data), "_"), "[[", 2)))

## Separate young and old data
young_data = C57BL6_mouse_data[unique(row.names(C57BL6_mouse_data)),which(cell_age == "young")]
old_data   = C57BL6_mouse_data[unique(row.names(C57BL6_mouse_data)),which(cell_age == "old")]

## Set up young mouse Seurat object
youngMouseSeuratObj <- CreateSeuratObject(counts = young_data, project = "MOUSE_AGE", min.cells = 0)
youngMouseSeuratObj <- NormalizeData(youngMouseSeuratObj)
youngMouseSeuratObj <- ScaleData(youngMouseSeuratObj, do.scale=T, do.center=T, display.progress = T)

## Set up old mouse Seurat object
oldMouseSeuratObj <- CreateSeuratObject(counts = old_data, project = "MOUSE_AGE", min.cells = 0)
oldMouseSeuratObj <- NormalizeData(oldMouseSeuratObj)
oldMouseSeuratObj <- ScaleData(oldMouseSeuratObj, do.scale=T, do.center=T, display.progress = T)

## Gene selection
youngMouseSeuratObj <- FindVariableFeatures(youngMouseSeuratObj, do.plot = F, nFeature=3000)
oldMouseSeuratObj <- FindVariableFeatures(oldMouseSeuratObj, do.plot = F, nFeature=3000,)
genes.use = Reduce(intersect, list(VariableFeatures(youngMouseSeuratObj),
                                   VariableFeatures(oldMouseSeuratObj),
                                   rownames(youngMouseSeuratObj),
                                   rownames(oldMouseSeuratObj)))

## Combine our Seurat objects
hsc.combined <- merge(youngMouseSeuratObj, oldMouseSeuratObj, add.cell.ids = c("YOUNG", "OLD"), project = "HSC")
hsc.combined <- ScaleData(hsc.combined, do.scale=T, do.center=T, display.progress = T)
VariableFeatures(hsc.combined) = genes.use

## Run PCA and UMAP
hsc.combined = RunPCA(hsc.combined, do.print=FALSE)
hsc.combined = RunUMAP(hsc.combined, dims = 1:30)

## Plot UMAP results
plot.me <- data.frame(x=hsc.combined@reductions$umap@cell.embeddings[,1],
                      y=hsc.combined@reductions$umap@cell.embeddings[,2],
                      labels=Idents(hsc.combined),
                      stringsAsFactors=FALSE)
unaligned.plot <- ggplot(plot.me, aes(x=x, y=y, colour = labels)) +
  geom_point(size=2) +
  scale_colour_manual(values=c("blue", "red")) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.line = element_line(colour = 'black', linewidth=1))
plot(unaligned.plot)

## Create paired dataset SCE objects to pass into scAlignCreateObject
youngMouseSCE <- SingleCellExperiment(
  assays = list(counts = youngMouseSeuratObj@assays$RNA$counts[genes.use,],
                logcounts  = youngMouseSeuratObj@assays$RNA$data[genes.use,],
                scale.data = youngMouseSeuratObj@assays$RNA$scale.data[genes.use,])
)

oldMouseSCE <- SingleCellExperiment(
  assays = list(counts = oldMouseSeuratObj@assays$RNA$counts[genes.use,],
                logcounts  = oldMouseSeuratObj@assays$RNA$data[genes.use,],
                scale.data = oldMouseSeuratObj@assays$RNA$scale.data[genes.use,])
)

scAlignHSC = scAlignCreateObject(sce.objects = list("YOUNG"=youngMouseSCE, "OLD"=oldMouseSCE),
                                 labels = list(cell_type[which(cell_age == "young")], cell_type[which(cell_age == "old")]),
                                 data.use="scale.data",
                                 pca.reduce = TRUE,
                                 pcs.compute = 50,
                                 cca.reduce = TRUE,
                                 ccs.compute = 15,
                                 project.name = "scAlign_Kowalcyzk_HSC")

scAlignHSC = scAlign(scAlignHSC,
                     options=scAlignOptions(steps=5000, log.every=5000, norm=TRUE, batch.norm.layer=TRUE, early.stop=FALSE, architecture="small"),
                     encoder.data="scale.data",
                     decoder.data="logcounts",
                     supervised='none',
                     run.encoder=TRUE,
                     run.decoder=TRUE,
                     log.dir=file.path(results.dir, 'models','gene_input'),
                     device="CPU")

## Additional run of scAlign with PCA, the early.stopping heuristic terminates the training procedure too early with PCs as input so it is disabled.
scAlignHSC = scAlign(scAlignHSC,
                     options=scAlignOptions(steps=15000, log.every=1000, norm=TRUE, batch.norm.layer=TRUE, early.stop=FALSE),
                     encoder.data="PCA",
                     supervised='none',
                     run.encoder=TRUE,
                     run.decoder=FALSE,
                     log.dir=file.path(results.dir, 'models','pca_input'),
                     device="CPU")

## Additional run of scAlign with CCA
scAlignHSC = scAlign(scAlignHSC,
                     options=scAlignOptions(steps=5000, log.every=1000, norm=TRUE, batch.norm.layer=TRUE, early.stop=TRUE),
                     encoder.data="CCA",
                     supervised='none',
                     run.encoder=TRUE,
                     run.decoder=FALSE,
                     log.dir=file.path(results.dir, 'models','cca_input'),
                     device="CPU")

## Plot aligned data in tSNE space, when the data was processed in three different ways: 1) either using the original gene inputs, 2) after PCA dimensionality reduction for preprocessing, or 3) after CCA dimensionality reduction for preprocessing. Cells here are colored by input labels
set.seed(5678)
gene_plot = PlotTSNE(scAlignHSC, "ALIGNED-GENE", title="scAlign-Gene", perplexity=30)
pca_plot = PlotTSNE(scAlignHSC, "ALIGNED-PCA", title="scAlign-PCA", perplexity=30)
cca_plot = PlotTSNE(scAlignHSC, "ALIGNED-CCA", title="scAlign-CCA", perplexity=30)
legend = get_legend(PlotTSNE(scAlignHSC, "ALIGNED-GENE", title="scAlign-Gene", legend="right", max_iter=1))
combined_plot = grid.arrange(gene_plot, pca_plot, cca_plot, legend, nrow = 1, layout_matrix=cbind(1,1,1,2,2,2,3,3,3,4))

## Plot aligned data in tSNE space, when the data was processed in three different ways: 1) either using the original gene inputs, 2) after PCA dimensionality reduction for preprocessing, or 3) after CCA dimensionality reduction for preprocessing. Cells here are colored by dataset.
set.seed(5678)
gene_plot = PlotTSNE(scAlignHSC, "ALIGNED-GENE", title="scAlign-Gene", cols=c("red","blue"), labels.use="group.by", perplexity=30)
pca_plot = PlotTSNE(scAlignHSC, "ALIGNED-PCA", title="scAlign-PCA", cols=c("red","blue"), labels.use="group.by", perplexity=30)
cca_plot = PlotTSNE(scAlignHSC, "ALIGNED-CCA", title="scAlign-CCA", cols=c("red","blue"), labels.use="group.by", perplexity=30)
legend = get_legend(PlotTSNE(scAlignHSC, "ALIGNED-GENE", title="scAlign-Gene", cols=c("red","blue"), labels.use="group.by", legend="right", max_iter=1))
combined_plot = grid.arrange(gene_plot, pca_plot, cca_plot, legend, nrow = 1, layout_matrix=cbind(1,1,1,2,2,2,3,3,3,4))
