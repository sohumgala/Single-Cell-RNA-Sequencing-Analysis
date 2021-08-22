# After cell types have been determined, it is easy to label the UMAP with appropriate cell types
# INPUT: cell ids based on previously specified annotation techniques
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(cells.filt)
cells.filt <- RenameIdents(cells.filt, new.cluster.ids)
DimPlot(cells.filt, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# Further customization of the plot is possible 

# save the plot by exporting from plots window
