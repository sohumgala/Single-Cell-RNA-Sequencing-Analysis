# install necessary packages and load libraries
# Seurat is a package that has tools for sc-RNAseq data analysis (Butler et al.)
# script is made using a workflow specified by Satija Lab (https://satijalab.org/seurat/)
# this script is a template to be used, values must be changed where there is "INPUT" commented
install.packages("Seurat")
install.packages("dplyr")
install.packages("patchwork")
library(Seurat)
library(dplyr)
library(patchwork)

# Load data and set up Seurat object
# INPUT: data.dir must be set to the file path of the folder that contains the 3 raw data matrices
# NOTE: file path formatting is different on different operating systems (could be double slash, single slash, or backslash)
cells.data <- Read10X(data.dir = "C://Users//Sohum//Documents//URECAStuff//pbmc3k_filtered_gene_bc_matrices//filtered_gene_bc_matrices//hg19")

# create Seurat object using data loaded above
# cells is an object of type Seurat where the data is in a format that is easily analyzed
# INPUT: project can be set to a specific project name
cells <- CreateSeuratObject(counts = cells.data, project = "cells", min.cells = 3, min.features = 200)

# Filter out cells that are of low quality/create noise
# High mitochondrial gene percentages can indicate low cell quality, so add column in metadata for percent.mt using [[]] operator
# INPUT: Depending on the species, the prefix for mitochondrial genes may be different than "^MT-" (as is the case for mice and rats)
cells[["percent.mt"]] <- PercentageFeatureSet(cells, pattern = "^MT-")

# other metrics are automatically stored in metadata, so we can visualize them as shown below in order to set cutoffs for filtered data
# violin plots
VlnPlot(cells, features = "nFeature_RNA")
VlnPlot(cells, features = "nCount_RNA")
VlnPlot(cells, features = "percent.mt")

# density plots
plot(density(cells@meta.data$nFeature_RNA))
plot(density(cells@meta.data$nCount_RNA))
plot(density(cells@meta.data$percent.mt))

# subset data based on above Quality Control metrics
# cutoffs are based on user's discretion, user must be mindful of how many cells are lost
# to check how many cells are lost, run "cells" and then "cells.filt" to see the difference in the number of samples
# INPUT: subset parameter must be adjusted based on above QC metrics
cells.filt <- subset(cells, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# normalize filtered data, parameters can be changed, but defaults usually work fine
cells.filt <- NormalizeData(cells.filt, normalization.method = "LogNormalize", scale.factor = 10000)

# identify the most variable genes across the dataset
# again, nfeatures parameter can be set, but default also works fine
cells.filt <- FindVariableFeatures(cells.filt, selection.method = "vst", nfeatures = 2000)

# scale data so that highly expressed genes do not dominate, another sort of normalization step
all.genes <- rownames(cells.filt)
cells.filt <- ScaleData(cells.filt, features = all.genes)

# linear dimension reduction (PCA)
cells.filt <- RunPCA(cells.filt, features = VariableFeatures(object = cells.filt))

# PCA groups data into principal components(PCs) that contain decreasingly significant data
# heatmaps can be used to visualize, below shows first 15
DimHeatmap(cells.filt, dims = 1:15, cells = 500, balanced = TRUE)

# the following functions can aid in deciding the number of PCs to use
# JackStrawPlot and ElbowPlot are useful for this purpose
# JackStraw is computationally expensive, so keep in mind
cells.filt <- JackStraw(cells.filt, num.replicate = 100)
cells.filt <- ScoreJackStraw(cells.filt, dims = 1:20)
JackStrawPlot(cells.filt, dims = 1:20)

ElbowPlot(cells.filt)

# Cluster the cells based on average gene expression
# INPUT: specify dims parameter as the number of PCs
# INPUT: specify resolution of the plot (lower resolution means fewer clusters with more cells), this may take experimentation
cells.filt <- FindNeighbors(cells.filt, dims = 1:10)
cells.filt <- FindClusters(cells.filt, resolution = 0.5)

# UMAP (Uniform Manifold Approximation and Projection) is used to visualize clusters
cells.filt <- RunUMAP(cells.filt, dims = 1:10)
DimPlot(cells.filt, reduction = "umap")

# There are 2 methods for finding cluster labels
# 1. Manually finding differentially expressed genes
cell.markers <- FindAllMarkers(cells.filt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cell.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#    INPUT: Write to csv in specified output directory
write.csv(cell.markers, "C://Users//Sohum//Documents//URECAStuff//Tables//cellmarkers.csv")
#    Now able to look at markers in a supervised way and use literature to assign labels
#    Other functions and visualizations that are useful for looking at canonical markers:
#    INPUT: Gene names
FeaturePlot(cells.filt, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
VlnPlot(cells.filt, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
#    This is useful when automated methods don't have enough input data or return ambiguous results

# 2. Use algorithms (GSVA and METANEIGHBOR) to partially automate this process
#    In order to do this, need to have a table of average gene expression of all genes per cluster
#    INPUT: Write to tsv in specified output directory
cluster.markers <- AverageExpression(cells.filt)
write.table(cluster.markers, "C://Users//Sohum//Documents//URECAStuff//Tables//clustermarkers.tsv", quote = F, sep = '\t', col.names = NA)
#    Next 2 Scripts will contain code for using this table and others for semi-automated clustering annotation

# A combination of these two methods may be used in order to label cell types most accurately


