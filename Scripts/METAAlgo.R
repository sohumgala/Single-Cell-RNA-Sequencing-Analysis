# Script (Diaz-Mejia et al) for GSVA annotation
# install necessary packages
install.packages("optparse")
install.packages("SummarizedExperiment")
install.packages("data.table")
install.packages("GSA")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MetaNeighbor")

suppressPackageStartupMessages(library(optparse))             
suppressPackageStartupMessages(library(SummarizedExperiment)) 
suppressPackageStartupMessages(library(data.table))           
suppressPackageStartupMessages(library(GSA))                  
suppressPackageStartupMessages(library(MetaNeighbor))         

oldw <- getOption("warn")
options( warn = -1 )

option_list <- list(
  # INPUT: replace default with tsv of average gene expression per cluster (from previous script)
  make_option(c("-i", "--infile_mat"), default = "C://Users//Sohum//Documents//URECAStuff//Tables//sampleclustermarkers.tsv"),
  # INPUT: replace default with file path to a .gmt file of canonical markers for cell types
  # able to create such a file in excel with a .gmt extension
  # first column is cell type, second column is blank
  # enter genes to be considered as canonical markers across the rows from the third column onwards
  # Example: https://github.com/jdime/scRNAseq_cell_cluster_labeling/tree/master/examples/inputs
  make_option(c("-c", "--infile_signature"), default = "C://Users//Sohum//Documents//URECAStuff//Tables//celltypemarkers.gmt"),
  # leave as gmt since file above is .gmt format
  make_option(c("-t", "--type_signature"), default = "gmt"),
  # if some clusters are known, can specify file path to .tsv with positive gold standards
  # clust1, clust2, clust3 names go down first column, and cell type go in the cooresponding second column
  # IMPORTANT: cell type names must be the same as specified in the .gmt file
  make_option(c("-g", "--infile_cell_types"), default = "NA"),
  # INPUT: replace default with output directory 
  make_option(c("-o", "--outdir"), default = "C://Users//Sohum//Documents//URECAStuff//METANEIGHBOR"),
  # INPUT: replace with a prefix on the output files if desired
  make_option(c("-p", "--prefix_outfiles"), default = "cells")
  
)

opt <- parse_args(OptionParser(option_list=option_list))

InfileMat       <- opt$infile_mat
InfileSignature <- opt$infile_signature
TypeSignature   <- opt$type_signature
InfileCellTypes <- opt$infile_cell_types
Outdir          <- opt$outdir
PrefixOutfiles  <- opt$prefix_outfiles

# INPUT: change temporary directory to a folder with a short file path, nothing will be saved here
Tempdir         <- "C://Users//Sohum//Documents//URECAStuff//Tempdir" 

StartTimeOverall<-Sys.time()

writeLines("\n*** Check that mandatory parameters are not 'NA' (default) ***\n")

ListMandatory<-list("infile_mat", "infile_signature", "type_signature", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
  }
}




dir.create(file.path(Outdir, "METANEIGHBOR"), showWarnings = F, recursive = T)
dir.create(file.path(Tempdir),        showWarnings = F, recursive = T)

writeLines("\n*** Load the matrix of genes (rows) vs. cell clusters (columns) ***\n")

### Load the matrix of genes (rows) vs. cell clusters (columns)
InfileMat.df<-data.frame(fread(InfileMat, sep="\t", na.strings=c("NA")), row.names=1)

writeLines("\n*** Load the cell type signatures ***\n")

if (length(grep('^gmt$', TypeSignature, perl = T, ignore.case = T))) {

  gmt1<-GSA.read.gmt(InfileSignature)
  gmt2<-list()
  for (i in 1:length(gmt1[[1]])){
    tmp<-unlist(gmt1[[1]][i])
    gmt2[[gmt1[[3]][i]]]<-tmp[which(tmp!="")]
  }
  AllGenesInSignature<-unique(unlist(x = gmt2, recursive = T, use.names = F))
  
  SigMat.df<-data.frame()
  for (gene_set_id in  names(gmt2)) {
    for (gene in AllGenesInSignature) {
      if (gene %in% gmt2[[gene_set_id]]) {
        SigMat.df[gene,gene_set_id] <- 1
      }else{
        SigMat.df[gene,gene_set_id] <- 0
      }
    }
  }
  
} else if (length(grep('^mat$', TypeSignature, perl = T, ignore.case = T))) {
  
  SigMat.df<-data.frame(fread(InfileSignature, sep="\t", na.strings=c("NA")), row.names=1)
  AllGenesInSignature<-rownames(SigMat.df)
  
} else {
  stop(paste("Parameter -c ", TypeSignature, " must be either  'gmt' or 'mat'. Use option -h for help.", sep = "", collapse = ""))
}
FullMat.df<-merge(SigMat.df, InfileMat.df, by = "row.names", all = T)
FullMat.df[is.na(FullMat.df)] <- 0
rownames(FullMat.df)<-FullMat.df[,"Row.names"]
FullMat.df$Row.names<-NULL

StudyIDs<-c(replicate(ncol(SigMat.df), "Signature"), replicate(ncol(InfileMat.df), "DataMatrix"))

writeLines("\n*** Load the cell type annotations ***\n")

InfileCellTypes.df<-data.frame(fread(InfileCellTypes, sep="\t", header = F, select = c(1:2)), row.names = 1)
CellTypes<-c(colnames(SigMat.df), c(colnames(InfileMat.df)))

writeLines("\n*** Create SummarizedExperiment object ***\n")

FullMat.df$chr<-1
FullMat.df$start<-1
FullMat.df$end<-1
FullMat.df$strand<-"*"

FullMat.seo <- makeSummarizedExperimentFromDataFrame(FullMat.df)
FullMat.seo

writeLines("\n*** Run MetaNeighbor US ***\n")

StartTimeMetaneighborUS<-Sys.time()
AUROC_scores_US = MetaNeighborUS(var_genes = AllGenesInSignature,
                                 dat = FullMat.seo, 
                                 study_id = StudyIDs,
                                 cell_type = CellTypes)

EndTimeMetaneighborUS<-Sys.time()

RowsDataMatrix<-grep(pattern = "^DataMatrix", rownames(AUROC_scores_US))
ColsDataMatrix<-grep(pattern = "^Signature",  colnames(AUROC_scores_US))
AUROC_scores_US_Datamatrix<-AUROC_scores_US[RowsDataMatrix,ColsDataMatrix]
colnames(AUROC_scores_US_Datamatrix)<-gsub(pattern = "^Signature\\|",  replacement = "", x = colnames(AUROC_scores_US_Datamatrix))
rownames(AUROC_scores_US_Datamatrix)<-gsub(pattern = "^DataMatrix\\|", replacement = "", x = rownames(AUROC_scores_US_Datamatrix))

OutfileEnrichmentScoresUS<-paste(Tempdir,"//",PrefixOutfiles,".MetaNeighborUS_AUROC.tsv", sep="")
Headers<-paste("MetaNeighborUS", paste(colnames(AUROC_scores_US_Datamatrix), sep="", collapse="\t"), sep="\t", collapse = "\t")
write.table(Headers,file = OutfileEnrichmentScoresUS, row.names = F, col.names = F, sep="\t", quote = F)
write.table(AUROC_scores_US_Datamatrix, file = OutfileEnrichmentScoresUS, row.names = T, col.names = F, sep="\t", quote = F, append = T)

writeLines("\n*** Report used options ***\n")

OutfileOptionsUsed<-paste(Tempdir,"//",PrefixOutfiles,".MetaNeighborUS_UsedOptions.txt", sep="")
TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=c(TimeOfRun,"\n","Options used:"))

for (optionInput in option_list) {
  write(file = OutfileOptionsUsed, x=(paste(optionInput@short_flag, optionInput@dest, opt[optionInput@dest], sep="\t", collapse="\t")),append = T)
}

writeLines("\n*** Report time used ***\n")

EndTimeOverall<-Sys.time()

TookTimeGsva    <-format(difftime(EndTimeMetaneighborUS, StartTimeMetaneighborUS, units = "secs"))
TookTimeOverall <-format(difftime(EndTimeOverall,        StartTimeOverall,        units = "secs"))

OutfileCPUusage<-paste(Tempdir,"/",PrefixOutfiles,".MetaNeighborUS_CPUusage.txt", sep="")
ReportTime<-c(
  paste("MetaNeighborUS",  TookTimeGsva,    collapse = "\t"),
  paste("overall",          TookTimeOverall, collapse = "\t")
)

write(file = OutfileCPUusage, x=c(ReportTime))

writeLines("\n*** Moving outfiles into outdir ***\n")
writeLines(paste(Outdir,"//METANEIGHBOR",sep="",collapse = ""))

outfiles_to_move <- list.files(Tempdir,pattern = paste(PrefixOutfiles, ".MetaNeighborUS_", sep=""), full.names = F)
sapply(outfiles_to_move,FUN=function(eachFile){
  
  file.copy(from=paste(Tempdir,"//",eachFile,sep=""),to=paste(Outdir,"//METANEIGHBOR",eachFile,sep=""),overwrite=T)
  file.remove(paste(Tempdir,"//",eachFile,sep=""))
})

options(warn = oldw)

print("END - All done!!! Took time:")
print(ReportTime)