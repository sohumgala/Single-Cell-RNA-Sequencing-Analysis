# use this script as a template of how to perform analysis on a dataset that includes control and stimulated cells
# script is made using a workflow specified by Satija Lab (https://satijalab.org/seurat/)
# install necessary packages and load libraries
install.packages("Seurat")
install.packages("cowplot")
library(Seurat)
library(cowplot)

# INPUT: file path for expression matrices as either .txt or .tsv
# take care to load the control dataset and stimulated dataset in the correct places
ctrl.data <- read.table(file = "C://Users//Sohum//Documents//URECAStuff//immune_alignment_expression_matrices//immune_control_expression_matrix.txt", sep = "\t")
stim.data <- read.table(file = "C://Users//Sohum//Documents//URECAStuff//immune_alignment_expression_matrices//immune_stimulated_expression_matrix.txt", sep = "\t")

# Set up both of the Seurat objects in a way that specifies which is is the control and which is the stimulated
# INPUT: can change qc parameters or number of variable genes if desired
# Set up control object
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "CELLS_CTRL", min.cells = 5)
ctrl$stim <- "CTRL"
ctrl <- subset(ctrl, subset = nFeature_RNA > 500)
ctrl <- NormalizeData(ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)

# Set up stimulated object
stim <- CreateSeuratObject(counts = stim.data, project = "CELLS_STIM", min.cells = 5)
stim$stim <- "STIM"
stim <- subset(stim, subset = nFeature_RNA > 500)
stim <- NormalizeData(stim, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)

# set up single Seurat object to perform integrated analysis
# INPUT: can change the number of PCs to be used, but 1:20 is standard
cells.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
cells.combined <- IntegrateData(anchorset = cells.anchors, dims = 1:20)

DefaultAssay(cells.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
cells.combined <- ScaleData(cells.combined, verbose = FALSE)
cells.combined <- RunPCA(cells.combined, npcs = 30, verbose = FALSE)
cells.combined <- RunUMAP(cells.combined, reduction = "pca", dims = 1:20)
cells.combined <- FindNeighbors(cells.combined, reduction = "pca", dims = 1:20)
cells.combined <- FindClusters(cells.combined, resolution = 0.5)

# Create visualizations to see distribution of data
p1 <- DimPlot(cells.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(cells.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(cells.combined, reduction = "umap", split.by = "stim")

# use conserved cell type markers for each cluster in order to identify cell types for each cluster
# below is an example for cluster 6 (NK Cells)
# Using this example, must do this for all clusters to identify cell type and annotate
# INPUT: adjust parameters to cover all clusters
DefaultAssay(cells.combined) <- "RNA"
six.markers <- FindConservedMarkers(cells.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
FeaturePlot(cells.combined, features = "GNLY")
# Using the table and feature plot, assign cell type labels to clusters
# INPUT: cell type names for clusters
cells.combined <- RenameIdents(cells.combined, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", 
                                `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "NK", `7` = "T activated", `8` = "DC", `9` = "B Activated", 
                                `10` = "Mk", `11` = "pDC", `12` = "Eryth")

DimPlot(cells.combined, label = TRUE)
# can export this plot from plots window

# There are many plots that can be used to visualize differences in gene expression, they can be found on https://satijalab.org/seurat/
# this table is very useful as a list of genes for the entire dataset
# The following code is a template, and can also be used for looking at differential gene expression within a specific cluster
cells.combined$status <- cells.combined@meta.data$orig.ident
Idents(cells.combined) <- "status"
differential.markers <- FindMarkers(cells.combined, ident.1 = "CELLS_CTRL", ident.2 = "CELLS_STIM", verbose = FALSE)
# in this case, pct.1 will be the control and pct.2 will be the stimulated
# one can look for genes that are highly differentially expressed
# write the table to a file
# INPUT: output file name and directory
write.csv(differential.markers, "C://Users//Sohum//Documents//URECAStuff//Tables//differentialmarkers.csv", row.names = FALSE)


