
# Methodology: scATAC-TF Cell Type Classification Pipeline

## Overview

This methodology describes a comprehensive machine learning approach for cell type classification using single-cell chromatin accessibility (scATAC-seq) and gene expression (scRNA-seq) data. The pipeline emphasizes automated cell annotation, class imbalance handling, and robust model evaluation.

---

## Phase 1: Data Preprocessing & Quality Control (R Implementation)

### Step 1: Data Acquisition
- **Input**: 10X Genomics multiome dataset (scATAC + scRNA)
- **Files**: `filtered_feature_bc_matrix.h5`, `per_barcode_summary_metrics.csv`, `atac_fragments.tsv.gz`
- **Validation**: Verify file integrity and expected dimensions

### Step 2: Quality Control Framework

#### RNA Quality Control (Team 1)
1. **Cell Filtering Criteria**:
   - Minimum genes per cell: 500
   - Maximum genes per cell: 6,000 (doublet threshold)
   - Minimum UMI counts: 800
   - Maximum mitochondrial gene percentage: 15%

2. **Gene Filtering**:
   - Complete removal of mitochondrial genes (^MT-)
   - Retain genes expressed in ≥10 cells
   - Remove ribosomal genes >50% of total expression

3. **Doublet Detection**:
   - Apply DoubletFinder with optimized parameters
   - Estimate doublet rate based on cell count (4-8%)
   - Remove predicted doublets maintaining cell type balance

#### ATAC Quality Control (Team 2)
1. **Cell Filtering Criteria**:
   - Minimum fragments per cell: 1,000
   - Maximum fragments per cell: 100,000
   - TSS enrichment score: >2
   - Nucleosome signal: <4
   - Fraction of reads in peaks: >15%

2. **Peak Filtering**:
   - Remove ENCODE blacklisted regions
   - Filter peaks accessible in <1% of cells
   - Apply TF-IDF normalization

### Step 3: Automated Cell Type Annotation

#### SingleR Implementation
1. **Reference Databases**:
   - Human Primary Cell Atlas (HPCA)
   - Monaco Immune Database (immune cell specific)

2. **Consensus Annotation**:
   - Primary: Monaco immune annotations
   - Fallback: HPCA annotations for missing cells
   - Standardization: Map to major PBMC types (T_cells, B_cells, Monocytes, NK_cells)

3. **Quality Assurance**:
   - Remove cell types with <10 cells
   - Validate against known marker genes
   - Generate confidence scores

### Step 4: Feature Engineering

#### RNA Features
- Select top 2,000 variable genes using variance stabilizing transformation
- Apply log-normalization with scale factor 10,000
- Generate normalized expression matrix

#### ATAC Features  
- Select top 5,000 most accessible peaks
- Apply TF-IDF normalization followed by LSI
- Calculate motif enrichment scores using JASPAR database

#### Integration
- Merge RNA and ATAC features for common cells
- Create combined feature matrix (2,000 genes + 5,000 peaks = 7,000 features)
- Ensure cell ID consistency between modalities

---

## Phase 2: Machine Learning Pipeline (Python Implementation)

### Step 5: Data Preparation

#### Feature Matrix Construction
- **Input**: Combined feature matrix [cells × 7,000 features]
- **Labels**: Automated cell type annotations
- **Validation**: Verify feature-label alignment

#### Train-Test Splitting
- **Method**: Stratified random split maintaining class proportions
- **Ratio**: 70% training, 30% testing
- **Seed**: Fixed random seed (42) for reproducibility

### Step 6: Class Imbalance Analysis & Correction

#### Imbalance Detection
1. **Metrics**:
   - Calculate class frequencies and proportions
   - Compute imbalance ratios (max_class_size / min_class_size)
   - Apply threshold: ratio >2.0 indicates significant imbalance

2. **Visualization**:
   - Bar plots of class distributions
   - Pie charts showing proportions
   - Imbalance ratio analysis

#### SMOTE Application
- **Condition**: Apply only if imbalance ratio >2.0
- **Parameters**: k_neighbors=3, random_state=42
- **Validation**: Compare class distributions before/after balancing
- **Visualization**: Side-by-side comparison plots

### Step 7: Feature Engineering & Selection

#### Preprocessing Steps
1. **Zero Variance Removal**: Drop features with no variation
2. **Correlation Filtering**: Remove features with correlation >0.95
3. **Univariate Selection**: SelectKBest with f_classif, k=1,000 features
4. **Standardization**: StandardScaler for SVM (fit on training only)

### Step 8: Model Development & Training

#### Model Configuration (Anti-Overfitting Focus)

**Random Forest**:
```python
RandomForestClassifier(
    n_estimators=100,
    max_depth=10,              # Prevent deep trees
    min_samples_split=20,      # Require more samples to split
    min_samples_leaf=10,       # Larger leaf nodes
    max_features='sqrt',       # Feature subsampling
    class_weight='balanced'    # Handle remaining imbalance
)
```

**XGBoost**:
```python
XGBClassifier(
    n_estimators=100,
    max_depth=6,               # Moderate depth
    learning_rate=0.1,         # Conservative learning rate
    subsample=0.8,             # Row subsampling
    colsample_bytree=0.8,      # Feature subsampling
    reg_alpha=0.1,             # L1 regularization
    reg_lambda=1.0             # L2 regularization
)
```

**Support Vector Machine**:
```python
SVC(
    C=1.0,                     # Regularization parameter
    kernel='rbf',
    gamma='scale',             # Automatic scaling
    class_weight='balanced',   # Handle imbalance
    probability=True           # Enable probability estimates
)
```

### Step 9: Model Evaluation Framework

#### Primary Metrics
- **Accuracy**: Overall classification correctness
- **Precision**: True positives / (True positives + False positives)
- **Recall**: True positives / (True positives + False negatives)
- **F1-Score**: Harmonic mean of precision and recall
- **AUC**: Area under ROC curve (multiclass: one-vs-rest)

#### Cross-Validation Analysis
- **Method**: 5-fold stratified cross-validation
- **Purpose**: Detect overfitting (train vs validation gap)
- **Threshold**: Gap >0.1 indicates significant overfitting
- **Reporting**: Mean and standard deviation across folds

#### Per-Class Analysis
- **Classification Reports**: Detailed metrics per cell type
- **Confusion Matrices**: Misclassification patterns
- **Class-Specific Performance**: Identify challenging cell types

### Step 10: Overfitting Detection & Prevention

#### Analysis Methods
1. **Learning Curves**: Plot training vs validation scores
2. **Cross-Validation Gaps**: Monitor train-validation differences  
3. **Feature Importance Stability**: Check consistency across folds

#### Prevention Strategies
1. **Regularization**: Built into model parameters
2. **Early Stopping**: Monitor validation performance
3. **Feature Selection**: Reduce dimensionality
4. **Ensemble Methods**: Reduce variance through averaging

---

## Phase 3: Results & Interpretation

### Step 11: Performance Visualization

#### Comparative Analysis
- **Model Comparison**: Bar plots of all metrics across models
- **Confusion Matrices**: Heatmaps for each model
- **ROC Curves**: Performance visualization (if binary)
- **Feature Importance**: Top contributing features per model

#### Class Imbalance Impact
- **Before/After SMOTE**: Distribution comparisons
- **Performance Impact**: Metrics with and without balancing
- **Class-Specific Improvements**: Per-class performance changes

### Step 12: Biological Validation

#### Marker Gene Analysis
- **Known Markers**: Validate against established cell type markers
- **Feature Importance**: Identify novel biomarkers
- **Biological Coherence**: Ensure predictions align with biology

#### External Validation
- **Literature Comparison**: Cross-reference with published studies
- **Database Validation**: Compare with ENCODE, ChIP-seq data
- **Pathway Analysis**: Assess biological pathway enrichment

### Step 13: Output Generation

#### Automated Reports
- **Performance Summary**: CSV with all model metrics
- **Detailed Analysis**: Excel workbook with per-class results
- **Visualization Suite**: High-quality plots for publication
- **JSON Report**: Machine-readable complete analysis

#### Reproducibility Assets
- **Parameter Logs**: All hyperparameters and random seeds
- **Environment Specs**: Package versions and system info
- **Data Splits**: Exact train/test cell assignments
- **Model Artifacts**: Saved trained models for reuse

---

## Quality Assurance Framework

### Validation Checkpoints
1. **Data Loading**: Verify dimensions and data types
2. **Preprocessing**: Check cell/feature retention rates
3. **Annotation**: Validate cell type assignments
4. **Feature Engineering**: Ensure proper normalization
5. **Model Training**: Monitor convergence and stability
6. **Evaluation**: Cross-validate all metrics

### Expected Performance Benchmarks
- **Classification Accuracy**: 85-95%
- **Cell Retention**: 70-80% after quality control
- **Feature Retention**: 15-25% after selection
- **Overfitting Gap**: <0.05 for production models
- **Cross-Validation Stability**: CV standard deviation <0.02

### Error Handling
- **Graceful Degradation**: Continue with warnings for non-critical errors
- **Automatic Fallbacks**: Alternative methods when primary fails
- **Comprehensive Logging**: Track all decisions and parameter choices
- **Validation Alerts**: Warn when benchmarks not met

This methodology ensures reproducible, robust, and biologically meaningful cell type classification while maintaining high computational standards and preventing common pitfalls in machine learning workflows.
