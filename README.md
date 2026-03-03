# Overview 
This pipeline processes 10x multiome output and performs:
      RNA preprocessing and normalization (SCTransform or standard workflow)
      ATAC chromatin assay construction
      RNA and ATAC quality control using MAD-based thresholds
      Dimensionality reduction (PCA, LSI)
      Multimodal interation using Weighted Nearest Neighbors (WNN)
      Clustering and UMAP visualization
      Optional peak calling with MACS2
      Clean export of processed objects and figures

# Design Principles
Reproducible: All dataset-specific parameters are controlled through a config file.
Modular: Each step (loading, QC, RNA processing, ATAC processing, integration) lives in its own script.
Dataset-agnostic: No hard-coded sample names, cluster IDs, or biological assumptions.
Extensible: Easily adapted for human or mouse data, MULTI-seq demultiplexing, peak-level re-quantification, downstream trajectory or motif analysis

# Repository Structure
.
├── config/
│   └── config.yaml
├── R/
│   ├── 00_packages.R
│   ├── 01_config.R
│   ├── 02_io_load.R
│   ├── 03_multiseq_demux.R
│   ├── 04_atac_assay.R
│   ├── 05_qc.R
│   ├── 06_rna_processing.R
│   ├── 07_multimodal_wnn.R
│   ├── 08_peak_calling.R
│   └── 99_utils.R
├── scripts/
│   └── run_pipeline.R
└── README.md

# Requirements
R ≥ 4.2 recommended
Core packages:
      Seurat
      Signac
      SoupX (optional)
      GenomeInfoDb
      EnsDb.Hsapiens.v86 (or mouse equivalent)
      ggplot2
      yaml
      MACS2 (optional, for peak calling)

# Input Requirements
Expected 10x multiome output:
filtered_feature_bc_matrix.h5
atac_fragments.tsv.gz
per_barcode_metrics.csv

Optional:
MULTI-seq FASTQs
Barcode reference CSV
All file paths are set in config.yaml.

# How to Run
1. Edit config/config.yaml with:
    Input directories
    Genome (hg38 or mm10)
    QC thresholds
    Clustering resolution
    Peak calling options
2. Run
    source("scripts/run_pipeline.R")
3. Output will be written to the specified output directory

# Example Workflow
Load 10x multiome output
Construct RNA and ATAC assays
Perform MAD-based filtering on:
      nFeature_RNA
      nCount_RNA
      percent.mt
      TSS enrichment
      nucleosome signal
Run SCTransform
Compute PCA (RNA) and LSI (ATAC)
Integrate modalities using WNN
Cluster and generate UMAP embeddings
Optionally call peaks and build a peak-level assay

# Extensions of this 
Snakemake or Nextflow integration
Docker containerization
Automated report generation
Batch integration across multiple samples
CI-based testing with toy datasets

