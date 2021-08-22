# A Semi-Automated Single Cell RNA-Seq Analysis Pipeline

## Description

This folder contains scripts to guide scRNA-seq data analysis processes including clustering, cluster annotation, and differential gene expression analysis.  These scripts are intended for those with little to no experience performing analyses of this type, and are adapted from existing workflows (Butler et al., Nature Biotechnology 2018) (Diaz-Mejia et al., 2019). 

## Scripts

|        Script Name         |                           Task(s)                            |
| :------------------------: | :----------------------------------------------------------: |
| **`IntroQCClustering.R`**  |     Load dataset, perform quality control, cluster cells     |
|      **`GSVAAlgo.R`**      |   Run the GSVA Algorithm for automating cluster annotation   |
|      **`METAAlgo.R`**      | Run the METANEIGHBOR Algorithm for automating cluster annotation |
|      **`Annotate.R`**      |            Annotate clusters with cell type name             |
| **`IntegratedAnalysis.R`** | Perform analysis on dataset containing stimulated and control cells together |

## Introductory Notes

Downloading R and using the RStudio environment is the easiest way to use these scripts.  The first four scripts are to be used in order on a single dataset that contains only one group of cells that were all treated the same way.  The last script is to be used when an analysis of a dataset containing cells of two groups is required, such as control vs. stimulated.  Each of the scripts have comments that guide the user to understand what is happening, as well as input fields and values specific to their data.  Additionally, file handling varies from operating system (Windows vs. Linux), so each of the automated annotation algorithms contain two scripts.  Use the one whose name matches the operating system being used.

## Brief Notes for Each Script

### **`IntroQCClustering.R`**

In this script, statements should be run one at a time, since certain commands and parameters are dependent on results from statements that come before.  Wherever "INPUT" is commented, that signifies that a value in the following line(s) of code must be replaced.  Most commonly these are numerical parameters for mathematical functions or file paths for loading/writing data. 

### **`GSVAAlgo.R`** and **`METAAlgo.R`**

These are the scripts that can be used to automate the annotation of cell clusters with cell type names.  The  outputs of this algorithm should be used along with manual methods in order to most accurately annotate clusters.  There are several inputs that are required for each algorithm to work.  These are denoted by comments in the code.  The fields in option_list have to be set manually to file/folder paths or numeric parameters.  Sample file inputs can be found on [GitHub](https://github.com/jdime/scRNAseq_cell_cluster_labeling/tree/master/examples/inputs) and [Zenodo](https://zenodo.org/record/3369934#.XtZWfzpKi70) (Diaz-Mejia et al., 2019).  Additionally a temporary output directory with short file path must be specified.  

Unlike the previous script, once all the input fields have been set, the entire script can be run at once.  Though this is the case, it is also possible to run lines or blocks individually to ensure everything is working as intended.

### **`Annotate.R`**

This is a very short script, and contains a few commands that are necessary in order to create an annotated plot of clusters.  There are many directions an investigation can go from here, such as sub clustering or differential gene expression analysis.  Sub clustering involves sub setting the original Seurat object and performing a similar analysis as is done in the first script.

### **`IntegratedAnalysis.R`**

The purpose of this script is to perform an introductory analysis of a dataset with two groups of the same cells in order to understand some of the differences in gene expression.  Similar to the first script, statements should be run one at a time since parameters might need to be adjusted based on results obtained by previous commands.

## Original Source Scripts

[Seurat](https://satijalab.org/seurat/)

[Diaz-Mejia](https://github.com/jdime/scRNAseq_cell_cluster_labeling)