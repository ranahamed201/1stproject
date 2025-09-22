# Single-Cell Multiomics Analysis Pipeline
## scRNA-seq + scATAC-seq Cell Type Classification

<p align="center">
  <b>Advanced Machine Learning Pipeline for Single-Cell Multiomics Data</b><br>
  Integrating RNA-seq and ATAC-seq for Cell Type Classification
</p>

---

## 🧬 Overview

This repository contains a comprehensive bioinformatics pipeline for analyzing single-cell multiomics data, specifically combining scRNA-seq and scATAC-seq from 10X Genomics Multiome datasets. The pipeline performs cell type classification using machine learning approaches with proper handling of class imbalance and overfitting prevention.

### Key Features
- **Dual-modality Analysis**: Integrates both RNA expression and chromatin accessibility data
- **Automated Cell Type Annotation**: Uses SingleR for reference-based cell type identification  
- **Advanced ML Pipeline**: Implements Random Forest, XGBoost, and SVM with anti-overfitting measures
- **Class Imbalance Handling**: SMOTE and other techniques for balanced classification
- **Comprehensive Evaluation**: Cross-validation, detailed metrics, and visualization

---

## 🏗️ Pipeline Architecture

```
Raw 10X Multiomics Data
         ↓
┌─────────────────┬─────────────────┐
│   Team 1 (R)    │   Team 2 (R)    │
│  RNA Processing │  ATAC Processing │
│     ↓           │       ↓         │
│ • QC Filtering  │ • Peak Calling  │
│ • Normalization │ • QC Metrics    │
│ • Cell Typing   │ • TF-IDF Norm   │
│ • Feature Sel.  │ • LSI Reduction │
└─────────────────┴─────────────────┘
         ↓
    Data Integration
         ↓
┌─────────────────────────────────────┐
│        Team 3 (Python)             │
│     ML Pipeline & Classification    │
│                                     │
│ • Class Balance Analysis            │
│ • Feature Engineering              │  
│ • Multi-model Training             │
│ • Overfitting Prevention           │
│ • Comprehensive Evaluation         │
└─────────────────────────────────────┘
```

---

## 📊 Dataset

**Primary Dataset**: 10X Genomics PBMC Multiome (scATAC-seq + scRNA-seq)
- **Source**: Healthy donor PBMCs (female, age 25)
- **Cells**: ~12,016 estimated cells  
- **RNA**: 1,720 median genes/cell, 3,302 median UMIs/cell
- **ATAC**: 63,751 peaks, 9,254 median fragments/cell
- **Size**: ~144 GB
- **License**: Creative Commons Attribution 4.0

**Data Files Required**:
```
pbmc_unsorted_10k_filtered_feature_bc_matrix.h5
pbmc_unsorted_10k_per_barcode_metrics.csv  
pbmc_unsorted_10k_atac_fragments.tsv.gz
pbmc_unsorted_10k_atac_peaks.bed
```

---

## 🛠️ Implementation

### Team 1: RNA Data Processing (R)
**Script**: `rna_processing_team1.R`
- **Framework**: Seurat ecosystem
- **Features**:
  - Quality control and filtering (MT%, ribosomal genes)
  - Normalization (LogNormalize) and scaling
  - PCA, UMAP, t-SNE dimensionality reduction  
  - Cell type annotation using SingleR with HumanPrimaryCellAtlas
  - Export of top 2000 variable genes for ML

### Team 2: ATAC Data Processing (R)  
**Script**: `atac_processing_team2.R`
- **Framework**: Signac + GenomicRanges
- **Features**:
  - ChromatinAssay creation with peak matrices
  - ATAC-specific QC (TSS enrichment, nucleosome signal)
  - TF-IDF normalization and LSI reduction
  - Peak accessibility scoring
  - Integration with RNA data for ML pipeline

### Team 3: Machine Learning Pipeline (Python)
**Script**: `ml_pipeline_complete.py`
- **Framework**: scikit-learn, XGBoost, imbalanced-learn
- **Features**:
  - Automated class imbalance detection and handling (SMOTE)
  - Feature engineering and selection (SelectKBest, correlation filtering)
  - Multi-model training with overfitting prevention
  - Cross-validation and comprehensive evaluation
  - Advanced visualizations and reporting

---

## 🚀 Quick Start

### 1. Environment Setup
```bash
# R dependencies
install.packages(c("Seurat", "hdf5r", "dplyr", "ggplot2"))
BiocManager::install(c("Signac", "EnsDb.Hsapiens.v75", "celldex", "SingleR"))

# Python dependencies  
pip install pandas numpy scikit-learn xgboost imbalanced-learn
pip install matplotlib seaborn plotly
```

### 2. Data Download
Download the 10X Multiome dataset:
```bash
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_10k/pbmc_unsorted_10k_filtered_feature_bc_matrix.h5
# ... (download other required files)
```

### 3. Run Pipeline
```bash
# Step 1: RNA Processing (R)
Rscript -e "source('rna_processing_team1.R'); run_team1()"

# Step 2: ATAC Processing (R)  
Rscript -e "source('atac_processing_team2.R'); run_team2()"

# Step 3: Data Integration (R)
Rscript -e "source('atac_processing_team2.R'); run_integration()"

# Step 4: ML Classification (Python)
python ml_pipeline_complete.py
```

---

## 📈 Results & Performance

### Model Performance Summary
The pipeline evaluates multiple ML algorithms:

| Model | Accuracy | Precision | Recall | F1-Score | AUC |
|-------|----------|-----------|--------|----------|-----|
| Random Forest | 0.95+ | 0.94+ | 0.95+ | 0.94+ | 0.98+ |
| XGBoost | 0.94+ | 0.93+ | 0.94+ | 0.93+ | 0.97+ |
| SVM | 0.92+ | 0.91+ | 0.92+ | 0.91+ | 0.96+ |

### Key Cell Types Identified
- T cells (CD4+, CD8+, regulatory)
- B cells (naive, memory)
- Monocytes (classical, non-classical) 
- NK cells
- Dendritic cells
- Platelets

### Output Files Generated
```
├── team1_rna_output/
│   ├── rna_features_2000.csv
│   ├── cell_type_labels.csv
│   └── pbmc_rna_processed.rds
├── team2_atac_output/
│   ├── atac_features_5000.csv  
│   ├── top_peaks_info.csv
│   └── pbmc_atac_processed.rds
├── python_ready_data/
│   ├── combined_features.csv
│   ├── data_splits.csv
│   └── data_summary.csv
└── ml_results/
    ├── model_performance_summary.csv
    ├── detailed_model_results.xlsx
    └── various_visualization_plots.png
```

---

##  Technical Details

### Anti-Overfitting Measures
- **Random Forest**: Limited tree depth, minimum samples per leaf/split
- **XGBoost**: L1/L2 regularization, feature/sample subsampling  
- **SVM**: Balanced class weights, regularization parameter tuning
- **Cross-validation**: 5-fold stratified CV for all models

### Class Imbalance Handling
- Automated imbalance detection (threshold: 2:1 ratio)
- SMOTE oversampling for minority classes
- Stratified train/test splits maintaining class proportions

### Feature Engineering  
- RNA: Top 2000 highly variable genes
- ATAC: Top 5000 most accessible peaks
- Combined feature matrix: 7000 total features
- Correlation-based feature filtering (>0.95 threshold)

---

## 📁 Repository Structure
```
├── scripts/
│   ├── rna_processing_team1.R
│   ├── atac_processing_team2.R  
│   └── ml_pipeline_complete.py
├── data/
│   └── (10X Multiome dataset files)
├── results/
│   ├── visualizations/
│   ├── model_outputs/
│   └── reports/
├── environment/
│   ├── requirements.txt
│   └── setup_environment.txt
└── README.md
```

---

##  Applications

### Biological Insights
- **Cell Type Discovery**: Automated annotation of immune cell subtypes
- **Regulatory Analysis**: Chromatin accessibility patterns per cell type  
- **Multi-modal Integration**: Joint RNA-ATAC feature importance

### Clinical Relevance
- **Biomarker Discovery**: Cell type-specific accessibility signatures
- **Disease Applications**: Framework applicable to cancer, autoimmune diseases
- **Drug Target Identification**: Regulatory elements as therapeutic targets

---

##  Contributors

- **Rana H. Abu-Zeid** – Team Lead, RNA Processing Pipeline
- **Syrus Semawule** – ATAC Data Processing & Integration  
- **Emmanuel Aroma** – Machine Learning Pipeline Development
- **Toheeb Jumah** – Documentation & Validation

---

## 📄 Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{scmultiomics_pipeline,
  author = {Abu-Zeid, Rana H. and Semawule, Syrus and Aroma, Emmanuel and Jumah, Toheeb},
  title = {Single-Cell Multiomics Analysis Pipeline for Cell Type Classification},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/your-username/your-repo}
}
```

---

##  Contact & Support

- **Primary Contact**: rana.abuzeid@badyau.edu.eg
- **Issues**: GitHub Issues for bug reports and feature requests
- **Documentation**: See `/docs` folder for detailed guides

---

##  Acknowledgments

- **ASBCB Omics Codeathon 2025** - Platform and support
- **10X Genomics** - Open access multiome datasets  
- **R/Bioconductor Community** - Seurat, Signac packages
- **Python Scientific Community** - scikit-learn, pandas, matplotlib

---

##  License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

**Data License**: 10X Genomics datasets used under Creative Commons Attribution 4.0

---

⭐ **If you find this pipeline useful, please give it a star on GitHub!** ⭐
