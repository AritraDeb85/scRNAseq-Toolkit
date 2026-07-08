# scRNAseq-Toolkit
A collection of modular R workflows for single-cell RNA-seq (scRNA-seq) analysis, including quality control, doublet detection, batch correction, clustering, cell annotation, visualization, and other reproducible bioinformatics pipelines.

## Overview
This repository provides modular and reusable R workflows for common single-cell RNA-seq (scRNA-seq) preprocessing and downstream analysis tasks. The current focus is on automating computationally repetitive steps such as doublet detection with DoubletFinder and dataset integration with Harmony for large multi-sample scRNA-seq datasets processed in Seurat.

## Features
- Automated doublet detection workflow for multi-sample scRNA-seq datasets using Seurat and DoubletFinder
- Planned support for Harmony-based batch correction and dataset integration
- Designed for reusable, modular analysis of large scRNA-seq sample collections
- Generates QC plots, intermediate Seurat objects, and doublet-filtered outputs

## Repository Structure
scRNAseq-Toolkit/
├── README.md
├── workflows/
│   ├── doublet_detection/
│   ├── integration/
│   ├── quality_control/
│   ├── clustering/
│   └── annotation/
├── docs/
└── example_data/

## Dependencies
The workflows in this repository currently rely on the following R packages:

- Seurat
- SeuratDisk
- DoubletFinder
- hdf5r
- tidyverse
- harmony

## Workflows

### Quality Control
Planned workflow module for basic scRNA-seq quality control and filtering steps.

### Doublet Detection
This module provides an automated DoubletFinder workflow for multi-sample scRNA-seq datasets processed in Seurat. It is designed to perform quality control, estimate DoubletFinder parameters, identify predicted doublets, and generate singlet-filtered Seurat objects together with diagnostic QC and UMAP visualizations for downstream analysis.

### Batch Correction
Planned workflow module for Harmony-based batch correction and dataset integration.

### Clustering
Planned workflow module for clustering and downstream cell population analysis.

### Cell Annotation
Planned workflow module for marker-based or reference-based cell type annotation.

## Installation

## Example Dataset

## Outputs

## Citation

## Author
**Aritra Deb**  
Research Fellow in Bioinformatics / Computational Biology
