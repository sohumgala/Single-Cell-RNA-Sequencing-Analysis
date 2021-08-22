# Script (Diaz-Mejia et al) for GSVA annotation
# install necessary packages
install.packages("optparse")
install.packages("parallel")
install.packages("data.table")
install.packages("GSA")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
BiocManager::install("GSVA")
BiocManager::install("qvalue")
install.packages("cluster")
install.packages("Seurat")

suppressPackageStartupMessages(library(optparse))     
suppressPackageStartupMessages(library(parallel))     
suppressPackageStartupMessages(library(data.table))   
suppressPackageStartupMessages(library(GSA))          
suppressPackageStartupMessages(library(GSVA))         
suppressPackageStartupMessages(library(qvalue))       
suppressPackageStartupMessages(library(cluster))      
suppressPackageStartupMessages(library(Seurat))       

oldw <- getOption("warn")
options( warn = -1 )



option_list <- list(
  # INPUT: replace default with tsv of average gene expression per cluster (from previous script)
  make_option(c("-i", "--infile_mat"), default = "C://Users//Sohum//Documents//URECAStuff//Tables//sampleclustermarkers.tsv"),
  # leave as "DGE" if file type is .tsv
  make_option(c("-t", "--infile_mat_type"), default = "DGE"),
  # INPUT: replace default with file path to a .gmt file of canonical markers for cell types
  # able to create such a file in excel or sheets, and download with a .gmt extension
  # first column is cell type, second column is blank
  # enter genes to be considered as canonical markers across the rows from the third column onwards
  # Example: https://github.com/jdime/scRNAseq_cell_cluster_labeling/tree/master/examples/inputs
  make_option(c("-c", "--infile_gmt"), default = "C://Users//Sohum//Documents//URECAStuff//Tables//celltypemarkers.gmt"),
  # INPUT: replace default with output directory 
  make_option(c("-o", "--outdir"), default = "C://Users//Sohum//Documents//URECAStuff//GSVA"),
  # INPUT: replace with a prefix on the output files if desired
  make_option(c("-p", "--prefix_outfiles"), default = "cells"),
  # can change p value and fdr cutoffs if desired, but defaults are generally accepted
  make_option(c("-e", "--pvalue_cutoff"), default = "0.05"),
  
  make_option(c("-f", "--fdr_cutoff"), default = "0.1")
  
)

opt <- parse_args(OptionParser(option_list=option_list))

InfileMat       <- opt$infile_mat
InfileMatType   <- opt$infile_mat_type
InfileGmt       <- opt$infile_gmt
Outdir          <- opt$outdir
PrefixOutfiles  <- opt$prefix_outfiles
PvalueCutoff    <- as.numeric(opt$pvalue_cutoff)
FdrCutoff       <- as.numeric(opt$fdr_cutoff)

# INPUT: change temporary directory to a folder with a short file path, nothing will be saved here 
Tempdir         <- "C://Users//Sohum//Documents//URECAStuff//Tempdir" 

StartTimeOverall<-Sys.time()



DefaultParameters <- list(
  DigitsForRound = 5 
)


writeLines("\n*** Check that mandatory parameters are not 'NA' (default) ***\n")

ListMandatory<-list("infile_mat", "infile_gmt", "outdir", "prefix_outfiles")
for (param in ListMandatory) {
  if (length(grep('^NA$',opt[[param]], perl = T))) {
    stop(paste("Parameter -", param, " can't be 'NA' (default). Use option -h for help.", sep = "", collapse = ""))
  }
}

dir.create(file.path(Outdir, "GSVA"), showWarnings = F, recursive = T)
dir.create(file.path(Tempdir),        showWarnings = F, recursive = T)

OutfileEnrichmentScores<-paste(Tempdir,"//",PrefixOutfiles,".GSVA_enrichment_scores.tsv", sep="")
OutfileEnrichScorsClust<-paste(Tempdir,"//",PrefixOutfiles,".GSVA_enrichment_scores_sorted.tsv", sep="")
OutfilePvalues         <-paste(Tempdir,"//",PrefixOutfiles,".GSVA_pvalues.tsv", sep="")
OutfileFdrvalues       <-paste(Tempdir,"//",PrefixOutfiles,".GSVA_fdr_values.tsv", sep="")
OutfileAllScores       <-paste(Tempdir,"//",PrefixOutfiles,".GSVA_all_scores_table.tsv", sep="")
OutfileFilteredES      <-paste(Tempdir,"//",PrefixOutfiles,".GSVA_filtered.tsv", sep="")
OutfileFinalLabel      <-paste(Tempdir,"//",PrefixOutfiles,".GSVA_final_label.tsv", sep="")

writeLines("\n*** Load data ***\n")



if (regexpr("^MTX$", InfileMatType, ignore.case = T)[1] == 1) {
  print("Loading MTX infiles")
  fullmat <- as.matrix(Read10X(data.dir = InfileMat))
}else if (regexpr("^DGE$", InfileMatType, ignore.case = T)[1] == 1) {
  print("Loading Digital Gene Expression matrix")
  fullmat <- as.matrix(data.frame(fread(InfileMat, check.names = F), row.names=1, check.names = F))
}else{
  stop(paste("Unexpected type of input: ", InfileMatType, "\n\nFor help type:\n\nRscript obtains_GSVA_for_MatrixColumns.R -h\n\n", sep=""))
}
rownames(fullmat) <- toupper(rownames(fullmat))

gmt1<-GSA.read.gmt(InfileGmt)
gmt2<-list()
for (i in 1:length(gmt1[[1]])){
  tmp<-unlist(gmt1[[1]][i])
  tmp <- toupper(tmp)
  gmt2[[gmt1[[3]][i]]]<-tmp[which(tmp!="")]
}

StartTimeGsva<-Sys.time()
EnrichmentScores<-gsva(expr=fullmat, gset.idx.list=gmt2, min.sz=1, max.sz=Inf, mx.diff=TRUE, verbose=T, parallel.sz=0)
SortedRowNames<-rownames(EnrichmentScores)
SortedColNames<-colnames(EnrichmentScores)
#
EnrichmentScores<-EnrichmentScores[SortedRowNames,SortedColNames]
write.table(data.frame("ENRICHMENT"=colnames(EnrichmentScores), t(round(x=EnrichmentScores, digits = DefaultParameters$DigitsForRound))), OutfileEnrichmentScores, row.names = F,sep="\t",quote = F)


Classes.list<-NULL
for (i in 1:nrow(EnrichmentScores)){
  tmp1<-unlist(strsplit(gmt1[[2]][match(rownames(EnrichmentScores)[i],gmt1[[3]])],"%"))[3]
  tmp2<-paste(unlist(gmt2[[match(rownames(EnrichmentScores)[i],names(gmt2))]]),collapse = ",")
  Classes.list<-rbind(Classes.list,c(tmp1,tmp2))
}
EndTimeGsva<-Sys.time()

writeLines("\n*** Get and print Enrichment, P-values and Q-values (FDR) ***\n")

HeaderCutoff<-paste(c("PassCutoff", "_p", PvalueCutoff, "_fdr", FdrCutoff) , sep="",collapse = "")
HeaderCutoff
HeadersForPandQvalues<-paste("CLASS","ColumnHeader","EnrichmentScore","p.Val","FDR", HeaderCutoff ,sep="\t",collapse = "")
write(x=HeadersForPandQvalues,file=OutfileAllScores)

for (columnNumber in 1:ncol(EnrichmentScores)){
  pvalues<-pnorm(-abs(scale(EnrichmentScores[,columnNumber])[,1]))
  qvalues<-qvalue(pvalues,pi0=1)$lfdr
  PassCutoffs<-ifelse((pvalues<=PvalueCutoff & qvalues<=FdrCutoff)==TRUE,1,0)
  concatenatedResults<-cbind(colnames(fullmat)[columnNumber], round(x=EnrichmentScores[,columnNumber], digits = DefaultParameters$DigitsForRound) , pvalues,qvalues,PassCutoffs)
  write.table(x=data.frame("CLASS"=rownames(concatenatedResults),concatenatedResults),file=OutfileAllScores, row.names = F,sep="\t",quote = F,col.names = F,append = T)
}

dataEPQ <- read.table(file = OutfileAllScores,row.names = NULL ,header = T)
ForPvaluesMat<- data.frame(x=dataEPQ[,"CLASS"], y=dataEPQ[,"ColumnHeader"], z=dataEPQ[,"p.Val"])
PvaluesMat<-xtabs(z~x+y, data=ForPvaluesMat)
PvaluesMat<-PvaluesMat[SortedRowNames,SortedColNames]
HeadersForPvalues<-paste(c("PVALUES", colnames(PvaluesMat)), sep="\t", collapse = "\t")
write(x=HeadersForPvalues,file=OutfilePvalues)
write.table(x=PvaluesMat, file=OutfilePvalues, row.names = T, sep="\t", quote = F, col.names = F, append = T)
ForFDRvaluesMat <- data.frame(x=dataEPQ[,"CLASS"], y=dataEPQ[,"ColumnHeader"], z=dataEPQ[,"FDR"])
FdrvaluesMat<-xtabs(z~x+y, data=ForFDRvaluesMat)
FdrvaluesMat<-FdrvaluesMat[SortedRowNames,SortedColNames]
HeadersForFdrvalues<-paste(c("FDR_VALUES", colnames(FdrvaluesMat)), sep="\t", collapse = "\t")
write(x=HeadersForFdrvalues,file=OutfileFdrvalues)
write.table(x=FdrvaluesMat, file=OutfileFdrvalues, row.names = T, sep="\t", quote = F, col.names = F, append = T)
FilteredESMatLogical<-(PvaluesMat<=PvalueCutoff & FdrvaluesMat<=FdrCutoff)
FilteredESMatLogical<-FilteredESMatLogical[SortedRowNames,SortedColNames]
FilteredESMatValues<-ifelse(FilteredESMatLogical==TRUE,EnrichmentScores,NA)
write.table(data.frame("ENRICHMENT_FILTERED"=rownames(EnrichmentScores), round(x=FilteredESMatValues, digits = DefaultParameters$DigitsForRound)),OutfileFilteredES, row.names = F,sep="\t",quote = F)


writeLines("\n*** Sort *.GSVA_enrichment_scores.tsv by similarity of prediction profiles ***\n")

predictions.mat<-as.matrix(data.frame(fread(OutfileEnrichmentScores, sep="\t", na.strings=c("NA")), row.names=1))
predictions.clusters.clust <- agnes(x = predictions.mat, metric = "manhattan")
predictions.clusters.order <- rownames(predictions.mat)[predictions.clusters.clust$order]

predictions.mat.t<-t(predictions.mat)
predictions.classes.clust <- agnes(x = predictions.mat.t, metric = "manhattan")
predictions.classes.order <- rownames(predictions.mat.t)[predictions.classes.clust$order]
predictions.mat.ordered   <- predictions.mat[predictions.clusters.order,predictions.classes.order]

write.table(data.frame("ENRICHMENT"=predictions.clusters.order, predictions.mat.ordered), OutfileEnrichScorsClust, row.names = F,sep="\t",quote = F)


write(paste(row.names(predictions.mat.ordered), colnames(predictions.mat.ordered)[max.col(predictions.mat.ordered, ties.method="first")], sep = "\t", collapse = "\n"),
      OutfileFinalLabel)


writeLines("\n*** Report used options ***\n")

OutfileOptionsUsed<-paste(Tempdir,"//",PrefixOutfiles,".GSVA_UsedOptions.txt", sep="")
TimeOfRun<-format(Sys.time(), "%a %b %d %Y %X")
write(file = OutfileOptionsUsed, x=c(TimeOfRun,"\n","Options used:"))

for (optionInput in option_list) {
  write(file = OutfileOptionsUsed, x=(paste(optionInput@short_flag, optionInput@dest, opt[optionInput@dest], sep="\t", collapse="\t")),append = T)
}

writeLines("\n*** Report time used ***\n")

EndTimeOverall<-Sys.time()

TookTimeGsva    <-format(difftime(EndTimeGsva,    StartTimeGsva,    units = "secs"))
TookTimeOverall <-format(difftime(EndTimeOverall, StartTimeOverall, units = "secs"))

OutfileCPUusage<-paste(Tempdir,"/",PrefixOutfiles,".GSVA_CPUusage.txt", sep="")
ReportTime<-c(
  paste("gsva",   TookTimeGsva,   collapse = "\t"),
  paste("overall",TookTimeOverall,collapse = "\t")
)

write(file = OutfileCPUusage, x=c(ReportTime))


writeLines("\n*** Moving outfiles into outdir ***\n")
writeLines(paste(Outdir,"//GSVA",sep="",collapse = ""))

outfiles_to_move <- list.files(Tempdir,pattern = paste(PrefixOutfiles, ".GSVA_", sep=""), full.names = F)
sapply(outfiles_to_move,FUN=function(eachFile){
  file.copy(from=paste(Tempdir,"//",eachFile,sep=""),to=paste(Outdir,"//GSVA",eachFile,sep=""),overwrite=T)
  file.remove(paste(Tempdir,"//",eachFile,sep=""))
})

options(warn = oldw)


print("END - All done!!! Took time:")
print(ReportTime)