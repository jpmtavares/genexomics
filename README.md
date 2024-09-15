# GeneXOmics: Perturb-Seq Pipeline with Snakemake and Cell Ranger

This repository contains a Snakemake workflow for running a small **Perturb-Seq** pipeline using **10x Genomics Cell Ranger**. The pipeline performs both primary and secondary analyses, including steps to:

1. Copy FASTQ files into the project folder.
2. Run Cell Ranger `count` for primary and secondary gene expression analysis.
3. Run a custom Python script that generates PCA and UMAP plots.

## Table of Contents

1. [Project Files and Folders](#project-files-and-folders)
   - condaenv.yml
   - Dockerfile
   - Snakefile
   - snakeconfig.yml
   - scripts/
2. [System Requirements](#system-requirements)
   - [Docker Setup](#docker-setup)
   - [Cell Ranger Setup](#cell-ranger-setup)
     - Hardware Requirements
     - Software Requirements
3. [Data](#data)
   - Reference transcriptome
   - FASTQ files
4. [Inputs](#inputs)
   - [Cell Ranger](#cell-ranger)
   - [Matrix Plot Generation Script](#matrix-plot-generation-script)
5. [Outputs](#outputs)
   - [Cell Ranger Outputs Overview](#cell-ranger-outputs-overview)
   - [Cell Ranger Primary Analysis](#cell-ranger-primary-analysis)
   - [Cell Ranger Secondary Analysis](#cell-ranger-secondary-analysis)
     - Analysis Folder Tree Structure
   - [Matrix Plot Generation](#matrix-plot-generation)
6. [How to Run](#how-to-run)

## Project Files and Folders

- **`condaenv.yml`**: 
  - Defines a Conda environment that installs **Snakemake** for managing the workflow.
  
- **`Dockerfile`**: 
  - Builds a Docker image based on `litd/docker-cellranger` (which runs on CentOS), freely available on [Docker Hub](https://hub.docker.com/r/litd/docker-cellranger). 
  - Installs necessary dependencies (e.g., `gcc`, `make`, `git`, `openssl-devel`, etc.), 
  - Installs Miniconda to manage `conda` environments and `mamba` for package management, creating the environment specified in `condaenv.yml` for running the pipeline.

- **`Snakefile`**: 
  - Implements the pipeline workflow for the Perturb-Seq analysis:
    1. Copy FASTQ files into the project directory.
    2. Run **Cell Ranger** for primary and secondary analysis using `cellranger count`.
      - Note that since `cellranger count` does not use the `--nosecondary` option, the results of the secondary analysis will be saved in the `analysis/` folder.

- **`snakeconfig.yml`**: 
  - Contains Snakemake [configuration](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html) parameters for copying FASTQ files and running the `cellranger count` step.
  - Accounts for sample name and sample number using the [naming convention](https://support.illumina.com/help/BaseSpace_Sequence_Hub_OLH_009008_2/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm) from Illumina FASTQ files (required for cellranger):
      ```txt
      SampleName_SampleNumber_Lane_Read_001.fastq.gz
      ```
  - Gets raw FASTQ and reference transcriptome from **10x Genomics** test data for **Cell Ranger** test runs:
      ```yml
      SOURCE_FASTQ_DIR: "/opt/cellranger-8.0.1/external/cellranger_tiny_fastq/"
      TRANSCRIPTOME_DIR: "/opt/cellranger-8.0.1/external/cellranger_tiny_ref"
      ```
    - You can specify custom paths in **Snakemake** using the `--config` parameter. For example:
      ```bash
      snakemake --config SOURCE_FASTQ_DIR=/path/to/fastqs TRANSCRIPTOME_DIR=/path/to/ref/transcriptome
      ```

- **`scripts/`**:
  - Contains custom scripts used in various stages of the pipeline. So far, the only script is:
    - `matrix_plotting.py`: This script is responsible for reading the filtered gene expression matrices produced by Cell Ranger and generating quality control metrics, normalization, clustering, and visualizations (PCA and UMAP plots).

## System Requirements

### Docker Setup
The pipeline requires Docker for containerization. To build and run the Docker image, you'll need Docker installed on your system. Follow the official [Docker installation guide](https://docs.docker.com/engine/install/).

### Cell Ranger Setup
#### Hardware Requirements
According to 10x Genomics Cell Ranger [system requirements](https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-system-requirements), your system should meet the following minimum specifications:

- **Processor**: 8-core Intel or AMD (16 cores recommended).
- **Memory**: 64 GB RAM (128 GB recommended).
- **Disk Space**: 1 TB of free disk space.
- **Operating System**: 64-bit CentOS/RedHat 7.0 or Ubuntu 14.04.

#### Software Requirements
The majority of software dependencies are bundled within the `litd/docker-cellranger` Docker image. If your data needs to be demultiplexed and converted from base call (BCL) files into FASTQ files, `cellranger mkfastq` command requires **Illumina bcl2fastq** (v2.20 or later), which is also installed in the Docker image.

## Data

This pipeline was created using small-scale test datasets provided by **10x Genomics**, specifically designed for Cell Ranger `testrun`. These datasets include both reference transcriptome files and sample FASTQ files required to simulate the analysis.

### Reference transcriptome
```bash
/opt/cellranger-8.0.1/external/cellranger_tiny_ref/
```

### FASTQ files
```bash
/opt/cellranger-8.0.1/external/cellranger_tiny_fastq/
```
- `tinygex_S1_L001_I1_001.fastq.gz`
- `tinygex_S1_L001_R1_001.fastq.gz`
- `tinygex_S1_L001_R2_001.fastq.gz`
- `tinygex_S1_L002_I1_001.fastq.gz`
- `tinygex_S1_L002_R1_001.fastq.gz`
- `tinygex_S1_L002_R2_001.fastq.gz`

## Inputs

### Cell Ranger
#### Required
- **FASTQ files**: Essential for running the `cellranger count` pipeline. These files are either provided by your sequencing core or can be generated from BCL files using `cellranger mkfastq`.
- **Reference transcriptome**: A reference transcriptome is required, which can be [downloaded](https://www.10xgenomics.com/support/software/cell-ranger/downloads/eula?closeUrl=%2Fsupport%2Fsoftware%2Fcell-ranger&lastTouchOfferName=Cell%20Ranger&lastTouchOfferType=Software%20Download&product=chromium&redirectUrl=%2Fsupport%2Fsoftware%2Fcell-ranger%2Fdownloads%23reference-downloads) from **10x Genomics** (after acknowledging to their terms of use and privacy policies) or generated using `cellranger mkref` for custom references with a reference genome sequence (FASTA file) and gene annotations (GTF file).

#### Optional
- **Libraries CSV**: Needed for analyzing Feature Barcode libraries. The Libraries CSV file declares the input FASTQ data for the libraries that make up a Feature Barcode experiment.
- **Feature Reference CSV**: Required for processing Feature Barcode data. It indicates the unique Feature Barcode sequence and the location of each feature present in your experiment.

Please refer to the **10x Genomics** Cell Ranger Inputs [page](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-inputs-overview) to know more.

### Matrix Plot Generation Script
#### Required
- **MTX and TSV files**: The directory containing the processed matrix (MTX) and accompanying TSV files from Cell Ranger output (`filtered_feature_bc_matrix/`) is required. These files represent the filtered feature-barcode matrix for gene expression data.

## Outputs

### Cell Ranger Outputs Overview
The `cellranger count` pipeline generates the following files:

| **Output**          | **Description**                             | **File Path**                                                          |
|---------------------|------------------------------|-------------------------------------------------------------------------|
| **Run summary HTML**| Overview of results          | `/home/loka/tinygex_S1/cellranger/count/outs/web_summary.html`          |
| **Run summary CSV** | Summary statistics for the run | `/home/loka/tinygex_S1/cellranger/count/outs/metrics_summary.csv`        |
| **BAM**             | Aligned sequence data                     | `/home/loka/tinygex_S1/cellranger/count/outs/possorted_genome_bam.bam`   |
| **BAM BAI index**   | Index for aligned sequence data | `/home/loka/tinygex_S1/cellranger/count/outs/possorted_genome_bam.bam.bai` |
| BAM CSI index       |                              | null                                                                    |
| **Filtered feature-barcode matrices (MEX)**   |    | `/home/loka/tinygex_S1/cellranger/count/outs/filtered_feature_bc_matrix` |
| **Filtered feature-barcode matrices (HDF5)**  |    | `/home/loka/tinygex_S1/cellranger/count/outs/filtered_feature_bc_matrix.h5` |
| **Unfiltered feature-barcode matrices (MEX)** |    | `/home/loka/tinygex_S1/cellranger/count/outs/raw_feature_bc_matrix`      |
| **Unfiltered feature-barcode matrices (HDF5)**|    | `/home/loka/tinygex_S1/cellranger/count/outs/raw_feature_bc_matrix.h5`   |
| **Secondary analysis folder** | Outputs in CSV     | `/home/loka/tinygex_S1/cellranger/count/outs/analysis`                   |
| **Per-molecule read information** |                | `/home/loka/tinygex_S1/cellranger/count/outs/molecule_info.h5`           |
| CRISPR-specific analysis          |                | null                                                                    |
| Antibody aggregate barcodes       |                | null                                                                    |
| **Loupe Browser file** | CLOUPE for visualization in Loupe Browser | `/home/loka/tinygex_S1/cellranger/count/outs/cloupe.cloupe`              |
| Feature Reference      |                           | null                                                                    |
| Target Panel File      |                           | null                                                                    |
| Probe Set File         |                           | null                                                                    |

### Cell Ranger Primary Analysis
This is the initial stage of data processing performed by Cell Ranger. It includes:
  1. **Alignment**: Reads from single-cell RNA sequencing are aligned to a reference transcriptome.
  - generated files:
    - `possorted_genome_bam.bam`
    - `possorted_genome_bam.bam.bai`
  2. **Barcode Processing**: Identifies and assigns the reads to individual cells based on unique cell barcodes.
  - output folder: `filtered_feature_bc_matrix/`
  - generated files: 
    - `barcodes.tsv`: lists the cell barcodes used in the analysis.
  3. **Feature Extraction**: Extracts features like gene expression levels and other metrics for each cell.
  - output folder: `filtered_feature_bc_matrix/`
  - generated files: 
    - `features.tsv`: lists the features (e.g., genes) detected.
  4. **Gene Counting**: Quantifies the number of transcripts for each gene within each cell.
  - output folder: `filtered_feature_bc_matrix/`
  - generated files: 
      - `matrix.mtx`: a sparse matrix where rows represent genes, columns represent cells, and values represent gene counts.

The output of primary analysis is typically a matrix of gene expression counts per cell, `matrix.mtx`, which forms the basis for subsequent analyses.

### Cell Ranger Secondary Analysis
The secondary analysis builds upon the results of the primary analysis and involves:
  1. **Quality Control**: Assesses and filters out low-quality cells or genes based on various metrics.
    - generated files:
      - `metrics_summary.csv`
  2. **Dimensionality Reduction**: Reduces the complexity of the data to make it easier to visualize and analyze.
  - output folder: `analysis/pca/`
  - generated files:
    - `projection.csv`
    - `components.csv`
    - `variance.csv`
  3. **Visualization in 2-D space**: After running PCA, t-distributed Stochastic Neighbor Embedding (t-SNE) and Uniform Manifold Approximation and Projection (UMAP) are run to visualize cells in a 2-D space.
  - output folder:
    - `analysis/tsne/`
    - `analysis/umap/`
  - generated files:
    - `projection.csv`
  4. **Clustering**: Groups cells into clusters based on their gene expression profiles, helping to identify distinct cell types or states.
  - output folder: `analysis/clustering/`
  - generated files:
    - `clusters.csv`: lists the cluster assignments for each cell.
  5. **Differential Expression**: Identifies genes that are differentially expressed between different clusters or conditions.
  - output folder: `analysis/diffexp/`
  - generated files:
    - `differential_expression.csv`: contains tables of differentially expressed genes between clusters or conditions.

For more details on secondary analysis outputs, refer to the **10x Genomics** Cell Ranger Secondary Analysis Outputs [documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-secondary-analysis).

#### Analysis Folder Tree Structure
```
analysis/
├── clustering/
│   ├── gene_expression_graphclust/
│   │   └── clusters.csv
│   ├── gene_expression_kmeans/
│   │   ├── 2_clusters/
│   │   │   └── clusters.csv
│   │   ├── 3_clusters/
│   │   │   └── clusters.csv
│   │   ├── 4_clusters/
│   │   │   └── clusters.csv
│   │   ├── 5_clusters/
│   │   │   └── clusters.csv
│   │   ├── 6_clusters/
│   │   │   └── clusters.csv
│   │   ├── 7_clusters/
│   │   │   └── clusters.csv
│   │   ├── 8_clusters/
│   │   │   └── clusters.csv
│   │   ├── 9_clusters/
│   │   │   └── clusters.csv
│   │   └── 10_clusters/
│   │       └── clusters.csv
│
├── diffexp/
│   ├── gene_expression_graphclust/
│   │   └── differential_expression.csv
│   ├── gene_expression_kmeans/
│   │   ├── 2_clusters/
│   │   │   └── differential_expression.csv
│   │   ├── 3_clusters/
│   │   │   └── differential_expression.csv
│   │   ├── 4_clusters/
│   │   │   └── differential_expression.csv
│   │   ├── 5_clusters/
│   │   │   └── differential_expression.csv
│   │   ├── 6_clusters/
│   │   │   └── differential_expression.csv
│   │   ├── 7_clusters/
│   │   │   └── differential_expression.csv
│   │   ├── 8_clusters/
│   │   │   └── differential_expression.csv
│   │   ├── 9_clusters/
│   │   │   └── differential_expression.csv
│   │   └── 10_clusters/
│   │       └── differential_expression.csv
│
├── pca/
│   └── gene_expression_10_components/
│       ├── components.csv
│       ├── dispersion.csv
│       ├── features_selected.csv
│       ├── projection.csv
│       └── variance.csv
│
├── tsne/
│   └── gene_expression_2_components/
│       └── projection.csv
│
└── umap/
    └── gene_expression_2_components/
        └── projection.csv
```

### Matrix Plot Generation

In addition to the cellranger outputs, the pipeline includes a custom Python script `scripts/matrix_plotting.py` that generates and saves visualizations of the PCA and UMAP projections. These plots help with the interpretation of the secondary analysis results.

## How to Run

1. **Build Docker Image**:
   ```bash
   docker build -t loka-perturb-seq .
   ```
2. **Run the Snakemake Workflow**:
   ```bash
   docker run -v $(pwd):/home/loka loka-perturb-seq snakemake -j 4
   ```