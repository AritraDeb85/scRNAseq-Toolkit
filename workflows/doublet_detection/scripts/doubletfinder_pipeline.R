############################################################
# doubletfinder_pipeline.R
#
# Author: Aritra Deb
# Version: 1.0.0
# Last updated: July 2026
#
# Developed and tested with:
#   Seurat v5
#   DoubletFinder (compatible with Seurat v5)
# Description:
# Automated DoubletFinder workflow for multi-sample
# single-cell RNA-seq datasets processed with Seurat.
#
# Inputs:
#   - 10x Genomics .h5 files
#   - Sample metadata
#
# Outputs:
#   - QC plots
#   - Doublet-filtered Seurat objects
#   - UMAP visualizations
#
# Tested on:
#   GSE305558 (human bone marrow scRNA-seq)
#
# This workflow can be adapted to other scRNA-seq datasets
# by modifying the sample information table and QC thresholds
# where appropriate.
############################################################
##############################################################
# Workflow Summary
#
# Input
#   ├── 10x Genomics .h5 files
#   └── Sample information table
#
# Processing
#   ├── Create Seurat object
#   ├── Quality control
#   ├── PCA
#   ├── Clustering
#   ├── Parameter optimization
#   ├── DoubletFinder
#   └── Remove predicted doublets
#
# Output
#   ├── QC figures
#   ├── Doublet UMAP
#   ├── Filtered Seurat object
#   └── Singlet UMAP
##############################################################
library(Seurat)
library(ggplot2)
library(SeuratDisk)
library(tidyverse)
library(hdf5r)
library(DoubletFinder)



##############################################################
# USER CONFIGURATION
#
# Before running this workflow on a new dataset:
#
# 1. Update the sample_info table below.
# 2. Modify QC thresholds if required.
# 3. Adjust QC thresholds and expected doublet-rate estimation
#    if required for your dataset.
# 4. Ensure input .h5 files are available in the working directory
#    or provide full file paths.
##############################################################
##Step 1. Create the sample information table
##############################################################
# Define the sample metadata for all samples to be processed.
# One row corresponds to one input sample.
sample_info <- data.frame(
  
  file = c(
    "GSM9179364_Young1_filtered_feature_bc_matrix.h5",
    "GSM9179365_Young2_filtered_feature_bc_matrix.h5",
    "GSM9179366_Young3_filtered_feature_bc_matrix.h5",
    "GSM9179367_Young4_filtered_feature_bc_matrix.h5",
    "GSM9179368_Young5_filtered_feature_bc_matrix.h5",
        
    "GSM9179369_Elderly1_filtered_feature_bc_matrix.h5",
    "GSM9179370_Elderly2_filtered_feature_bc_matrix.h5",
    "GSM9179371_Elderly3_filtered_feature_bc_matrix.h5",
    "GSM9179372_Elderly4_filtered_feature_bc_matrix.h5",
    "GSM9179373_Elderly5_filtered_feature_bc_matrix.h5",
    "GSM9179374_Elderly6_filtered_feature_bc_matrix.h5",
    "GSM9179375_Elderly7_filtered_feature_bc_matrix.h5"
      ),
  
  sample_title = c(
    "Young1","Young2","Young3","Young4",
    "Young5",
    
    "Elderly1","Elderly2","Elderly3","Elderly4",
    "Elderly5","Elderly6","Elderly7"
  ),
  
  sample_group = c(
    rep("Young",5),
    rep("Elderly",7)
  ),
  
   age = c(
    31, 28, 48, 34, 48,
    58, 66, 79, 63, 69, 71, 74
  ),

  sex = c(
    "Male","Female","Male","Male","Male",
    "Male","Female","Male","Male","Female","Male","Female"
  ),
  
  stringsAsFactors = FALSE
)
##############################################################
##Step 2. Write the function
##############################################################
run_doublet_detection <- function(
    h5_file,
    sample_group,
    sample_title
){
  
##############################################################
  ##Step 3. Read file automatically
##############################################################
  if(!file.exists(h5_file)){
    stop(
      paste(
        "Cannot find file:",
        h5_file
      )
    )
  }
  
  counts <- Read10X_h5(h5_file)
  
# Read10X_h5() may return either:
# - a sparse matrix, or
# - a list containing multiple assay types.
# Extract the Gene Expression matrix when multiple assays are present.
  if(is.list(counts)){
    
    cts <- counts$`Gene Expression`
    
  } else {
    
    cts <- counts
    
  }
  
  cat(
    "Genes:",
    nrow(cts),
    " Cells:",
    ncol(cts),
    "\n"
  )
##############################################################
  ##Step 4. Create Seurat object
##############################################################
  seu <- CreateSeuratObject(
    counts = cts,
    project = sample_title,
    min.cells = 3,
    min.features = 200
  )
  ##Take total cells loaded from unfiltered Seurat object
  raw_cell_number <- ncol(seu)
##############################################################
  ##Step 5. QC
##############################################################
  seu$mitoPercent <- PercentageFeatureSet(
    seu,
    pattern = "^MT-"
  )
  
  ##############################################################
  ## RAW QC PLOTS
  ##############################################################
  
  png(
    paste0(
      sample_title,
      "_QC_raw_violin.png"
    ),
    width = 1800,
    height = 600,
    res = 150
  )
  
  print(
    VlnPlot(
      seu,
      features = c(
        "nFeature_RNA",
        "nCount_RNA",
        "mitoPercent"
      ),
      ncol = 3,
      pt.size = 0.1
    )
  )
  
  dev.off()
  
  
  png(
    paste0(
      sample_title,
      "_QC_raw_scatter.png"
    ),
    width = 1200,
    height = 600,
    res = 150
  )
  
  p1 <- FeatureScatter(
    seu,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  )
  
  p2 <- FeatureScatter(
    seu,
    feature1 = "nCount_RNA",
    feature2 = "mitoPercent"
  )
  
  print(
    p1 + p2
  )
  
  dev.off()
  
  # Apply basic QC filtering.
# Thresholds should be adjusted according to dataset characteristics.
  seu <- subset(
    seu,
    subset =
      nCount_RNA > 800 &
      nFeature_RNA > 500 
#	  & mitoPercent < 10
  )
  
  cat(
    sample_title,
    ":",
    raw_cell_number,
    "raw cells ->",
    ncol(seu),
    "cells after QC\n"
  ) 
  
  ##############################################################
  ## POST-QC PLOTS
  ##############################################################
  
  png(
    paste0(
      sample_title,
      "_QC_filtered_violin.png"
    ),
    width = 1800,
    height = 600,
    res = 150
  )
  
  print(
    VlnPlot(
      seu,
      features = c(
        "nFeature_RNA",
        "nCount_RNA",
        "mitoPercent"
      ),
      ncol = 3,
      pt.size = 0.1
    )
  )
  
  dev.off()
  
  
  png(
    paste0(
      sample_title,
      "_QC_filtered_scatter.png"
    ),
    width = 1200,
    height = 600,
    res = 150
  )
  
  p1 <- FeatureScatter(
    seu,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  )
  
  p2 <- FeatureScatter(
    seu,
    feature1 = "nCount_RNA",
    feature2 = "mitoPercent"
  )
  
  print(
    p1 + p2
  )
  
  dev.off()
  
##############################################################
  ##Step 6. Standard preprocessing
##############################################################
  seu <- NormalizeData(seu)
  
  seu <- FindVariableFeatures(seu)
  
  seu <- ScaleData(seu)
  
  seu <- RunPCA(seu)
  
  seu <- FindNeighbors(
    seu,
    dims = 1:20
  )
  
  seu <- FindClusters(
    seu,
    resolution = 0.5
  )
  
  seu <- RunUMAP(
    seu,
    dims = 1:20
  )
  
##############################################################
  ##Step 7. pK calculation
##############################################################
  sweep.res <- paramSweep(
    seu,
    PCs = 1:20,
    sct = FALSE
  )
  
  sweep.stats <- summarizeSweep(
    sweep.res,
    GT = FALSE
  )
  
  bcmvn <- find.pK(
    sweep.stats
  )

## Select the pK value corresponding to the maximum BCmetric. 
  best_row <- which.max(bcmvn$BCmetric)
  
  pK <- as.numeric(
    as.character(
      bcmvn$pK[best_row]
    )
  )
  
  cat("Selected pK =", pK, "\n")
#  
##############################################################
  ##Step 8. Doublet estimation
##############################################################
  annotations <- as.character(
    Idents(seu)
  )
   
#  homotypic.prop <- modelHomotypic(
#    annotations
#  )

  cat("Running modelHomotypic...\n")
  
  homotypic.prop <- modelHomotypic(
    annotations
  )
  
  cat("homotypic.prop =", homotypic.prop, "\n")
  
## Estimate expected doublet rate
## based on the approximate number of loaded cells.
  if(raw_cell_number < 3000){
    expected_rate <- 0.02
  } else if(raw_cell_number < 6000){
    expected_rate <- 0.04
  } else if(raw_cell_number < 10000){
    expected_rate <- 0.05
  } else {
    expected_rate <- 0.06
  }
  
  nExp_poi <- round(
    expected_rate * raw_cell_number
  )
  
  nExp_poi.adj <- round(
    nExp_poi *
      (1 - homotypic.prop)
  )

  cat(
    "Raw cells:",
    raw_cell_number,
    "\n"
  )
  
  cat(
    "Expected doublets:",
    nExp_poi,
    "\n"
  )
  
  cat(
    "Adjusted doublets:",
    nExp_poi.adj,
    "\n"
  )
  
##############################################################
  ##Step 9. Run DoubletFinder
##############################################################
#  seu <- JoinLayers(seu)
  cat("Starting doubletFinder...\n")
  seu <- doubletFinder(
    seu,
    PCs = 1:20,
    pN = 0.25,
    pK = pK,
    nExp = nExp_poi.adj,
    reuse.pANN = NULL,
    sct = FALSE
  )
  
  saveRDS(
    seu,
    paste0(
      sample_title,
      "_with_doublet_calls.rds"
    )
  )
  
  ##############################################################
  ## DOUBLETFINDER UMAP
  ##############################################################
  
  df_col <- grep(
    "DF.classifications",
    colnames(seu@meta.data),
    value = TRUE
  )
  
  png(
    paste0(
      sample_title,
      "_doublet_umap.png"
    ),
    width = 1000,
    height = 800,
    res = 150
  )
  
  print(
    DimPlot(
      seu,
      group.by = df_col
    )
  )
  
  dev.off()
##############################################################
  ## Step 10. Automatically find DF column
##############################################################
  df_col <- grep(
    "DF.classifications",
    colnames(seu@meta.data),
    value = TRUE
  )
  
  print(df_col)
  
  print(
    colnames(seu@meta.data)
  )
  
  print(
    table(
      seu@meta.data[[df_col]]
    )
  )
  
##############################################################
  ## Step 11. Keep singlets
##############################################################
  singlet_cells <- rownames(
    seu@meta.data[
      seu@meta.data[[df_col]] == "Singlet",
      ,
      drop = FALSE
    ]
  )
  
  cat(
    "Keeping",
    length(singlet_cells),
    "singlets\n"
  )
  
  seu <- subset(
    seu,
    cells = singlet_cells
  )
  
##############################################################
  ##Step 12. Add metadata
##############################################################
  seu$sample_group <- sample_group
  
  seu$sample_title <- sample_title
  
##############################################################
  ##Step 13. Save automatically
##############################################################
  saveRDS(
    seu,
    paste0(
      sample_title,
      "_doublet_filtered.rds"
    )
  )
  
  ##############################################################
  ## FINAL SINGLET UMAP
  ##############################################################
  
  png(
    paste0(
      sample_title,
      "_final_singlets_umap.png"
    ),
    width = 1000,
    height = 800,
    res = 150
  )
  
  print(
    DimPlot(seu)
  )
  
  dev.off()
  
##############################################################
  ##Step 14. Close the function
##############################################################
  cat(
    sample_title,
    "finished successfully\n"
  )
  return(seu)
  
}

##############################################################
## OPTIONAL: Run the workflow on a single sample
##############################################################
##
##test_obj <- run_doublet_detection(
##  h5_file = sample_info$file[2],
##  sample_group = sample_info$sample_group[2],
##  sample_title = sample_info$sample_title[2]
##)
##

####################################################
## RUN ALL SAMPLES
####################################################

for(i in seq_len(nrow(sample_info))){
  
  cat(
    "\n====================================\n"
  )
  
  cat(
    "Processing:",
    sample_info$sample_title[i],
    "\n"
  )
  
  cat(
    "====================================\n"
  )
  
  tryCatch(
    
    {
      run_doublet_detection(
        h5_file      = sample_info$file[i],
        sample_group = sample_info$sample_group[i],
        sample_title = sample_info$sample_title[i]
      )
    },
    
    error = function(e){
      
      cat(
        "\nERROR in",
        sample_info$sample_title[i],
        "\n"
      )
      
      cat(
        conditionMessage(e),
        "\n"
      )
    }
    
  )
  
  gc()
  
}
cat("\nWorkflow completed successfully.\n")