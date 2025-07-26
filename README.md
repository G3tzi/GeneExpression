# GeneClusterAnalysis

A comprehensive R package for gene expression data clustering and pathway enrichment analysis. This package provides a complete workflow for analyzing gene expression patterns, from data preprocessing to biological interpretation through pathway enrichment.

## Features

- **Data Loading**: Support for multiple input formats (CSV, TSV, ExpressionSet)
- **Data Preprocessing**: Multiple normalization methods (log2, quantile, z-score, DESeq2) with optional batch correction
- **Exploratory Analysis**: Outlier detection and data quality assessment tools
- **Multiple Clustering Methods**: Hierarchical clustering, k-means, and DBSCAN algorithms
- **Cluster Validation**: Within-cluster sum of squares, silhouette analysis, and gap statistics
- **Rich Visualizations**: Heatmaps, PCA plots, dendrograms, and cluster profile plots
- **Pathway Enrichment**: Over-representation analysis (ORA) with support for online gene set databases
- **Parallel Processing**: Support for parallel execution of computationally intensive tasks
- **Complete Pipeline**: Automated end-to-end workflow with customizable parameters

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biobase", "preprocessCore", "DESeq2", "sva"))

install.packages(c("cluster", "dbscan", "ComplexHeatmap", "ggplot2", 
                   "reshape2", "matrixStats", "proxy", "msigdbr",
                   "BiocParallel", "stringr"))
```

## Quick Start

### Option 1: Using the Complete Pipeline (Recommended)

```r
library(GeneClusterAnalysis)

results <- run_pipeline(
  data_path = "expression_data.csv",
  data_type = "csv",
  pheno_data = "sample_metadata.csv",
  clustering_methods = c("hierarchical", "kmeans"),
  k_clusters = 5,
  perform_enrichment = TRUE
)

clusters <- results$clustering_results
plots <- results$plots
enrichment <- results$enrichment_results
summary <- results$summary
```

### Option 2: Step-by-Step Analysis

```r
library(GeneClusterAnalysis)

eset <- read_expression("expression_data.csv", type = "csv")

normalized_eset <- normalize_data(eset, method = "log2")
cluster_result <- cluster_expression(normalized_eset, method = "kmeans", k = 5)

plot_heatmap(normalized_eset, clusters = cluster_result$clusters)
```

## Data Loading and Input Formats

The package supports multiple input formats through the `read_expression()` function:

### Supported File Formats

#### 1. CSV Files
```r
eset <- read_expression("expression_data.csv", type = "csv")

eset <- read_expression(
  path = "expression_data.csv", 
  type = "csv",
  pheno_data = "sample_metadata.csv"
)
```

#### 2. TSV Files
```r
eset <- read_expression("expression_data.tsv", type = "tsv")
```

#### 3. ExpressionSet from .RData
```r
eset <- read_expression("expression_data.RData", type = "eset")
```

### Expected Data Structure

#### Expression Data File Format
```
Gene_ID,Sample1,Sample2,Sample3,Sample4
GENE1,10.5,12.3,8.7,11.2
GENE2,5.2,6.8,4.9,7.1
GENE3,15.8,14.2,16.1,13.9
...
```

- **First column**: Gene identifiers (used as row names)
- **Remaining columns**: Sample expression values
- **Header**: Required (first row should contain sample names)

#### Sample Metadata File Format (Optional)
```
Sample_ID,Condition,Batch,Treatment
Sample1,Control,Batch1,Untreated
Sample2,Treated,Batch1,Drug_A
Sample3,Control,Batch2,Untreated
Sample4,Treated,Batch2,Drug_A
```

- **First column**: Sample identifiers (must match expression data columns)
- **Additional columns**: Sample attributes for visualization and batch correction

### Creating ExpressionSet Manually

```r
library(Biobase)

expr_matrix <- your_expression_matrix  # genes × samples
pheno_data <- your_sample_info         # sample metadata data.frame

eset <- ExpressionSet(
  assayData = expr_matrix,
  phenoData = AnnotatedDataFrame(pheno_data)
)
```

## Complete Pipeline Usage

The `run_pipeline()` function provides a comprehensive analysis workflow:

### Basic Pipeline

```r
results <- run_pipeline(
  data_path = "expression_data.csv",
  data_type = "csv"
)
```

### Advanced Pipeline Configuration

```r
results <- run_pipeline(
  # Data input
  data_path = "expression_data.csv",
  data_type = "csv",
  pheno_data = "sample_metadata.csv",
  output_dir = "my_analysis_results",
  
  # Normalization
  normalize = TRUE,
  norm_method = "quantile",
  batch_correction = "batch_column",
  filter_low_expr = TRUE,
  filter_threshold = 1,
  
  # Outlier detection
  detect_outliers = TRUE,
  outlier_method = "mad",
  outlier_threshold = 3,
  
  # Clustering
  clustering_methods = c("hierarchical", "kmeans", "dbscan"),
  target = "gene",
  k_clusters = 6,
  dbscan_eps = 0.5,
  
  # Validation
  perform_validation = TRUE,
  gap_statistic = TRUE,
  
  # Visualization
  create_plots = TRUE,
  heatmap_scale = "row",
  pca_color_by = "condition",
  
  # Enrichment analysis
  perform_enrichment = TRUE,
  gene_sets_source = "online",
  gene_sets_collection = "H",  # Hallmark pathways
  gene_sets_species = "Homo sapiens",
  
  # General
  verbose = TRUE,
  seed = 123
)
```

### Pipeline Output Structure

The pipeline returns a comprehensive results object:

```r
str(results, max.level = 2)

# $original_data          - Original ExpressionSet
# $processed_data         - Processed/normalized ExpressionSet
# $outlier_info          - Outlier detection results
# $clustering_results    - Results from all clustering methods
#   $hierarchical        - Hierarchical clustering results
#   $kmeans             - K-means clustering results
#   $dbscan             - DBSCAN clustering results (if used)
# $validation_results    - Clustering validation metrics
# $plots                - All generated visualizations
# $enrichment_results   - Pathway enrichment analysis
# $summary              - Analysis summary statistics
```

### Accessing Results

```r
# Get clustering results
hclust_clusters <- results$clustering_results$hierarchical$clusters
kmeans_clusters <- results$clustering_results$kmeans$clusters

# View plots
print(results$plots$pca)
print(results$plots$hierarchical_heatmap)

# Get enrichment results
enrichment_summary <- results$enrichment_results$kmeans$summary
detailed_enrichment <- results$enrichment_results$kmeans$detailed_results

# View analysis summary
print(results$summary)
```

### Pipeline Parameters

#### Data Input Parameters
- `data_path`: Path to expression data file
- `expression_set`: Alternative to file path - provide ExpressionSet directly
- `data_type`: File format ("csv", "tsv", "eset")
- `pheno_data`: Sample metadata (file path or data.frame)
- `output_dir`: Directory for saving results

#### Normalization Parameters
- `normalize`: Whether to perform normalization
- `norm_method`: Normalization method ("log2", "quantile", "zscore", "deseq")
- `batch_correction`: Column name for batch correction
- `filter_low_expr`: Filter low-expression genes
- `filter_threshold`: Minimum expression threshold
- `filter_min_samples`: Minimum samples meeting threshold

#### Clustering Parameters
- `clustering_methods`: Vector of methods ("hierarchical", "kmeans", "dbscan")
- `target`: Cluster "gene" or "sample"
- `k_clusters`: Number of clusters for hierarchical/k-means
- `hclust_method`: Linkage method for hierarchical clustering
- `distance`: Distance metric
- `dbscan_eps`: DBSCAN epsilon parameter
- `dbscan_minPts`: DBSCAN minimum points parameter

#### Enrichment Parameters
- `perform_enrichment`: Whether to perform pathway analysis
- `gene_sets_source`: Source for gene sets ("online", "gmt", "example", "custom")
- `gene_sets_collection`: MSigDB collection for online source
- `gene_sets_species`: Species for gene set databases
- `enrichment_p_threshold`: P-value threshold for significance

## Detailed Usage

### 1. Data Preprocessing and Normalization

```r
eset <- read_expression("data.csv", type = "csv")

normalized_eset <- normalize_data(
  eset = eset,
  method = "log2",              # Options: "log2", "quantile", "zscore", "deseq"
  filter = TRUE,                # Remove low-expression genes
  filter_threshold = 10,
  filter_min_samples = 3
)

# With batch correction
normalized_eset <- normalize_data(
  eset = eset,
  method = "quantile",
  batch = "batch_column_name"   # Column in pData(eset)
)
```

### 2. Exploratory Analysis

```r
# Detect and remove outlier genes
filtered_eset <- detect_outliers(
  eset = normalized_eset,
  method = "mad",               # Options: "mad", "iqr"
  threshold = 3,
  plot = TRUE,                  # Show distribution plot
  remove = TRUE                 # Remove outliers from output
)
```

### 3. Gene Clustering

#### Hierarchical Clustering
```r
hc_result <- cluster_expression(
  eset = filtered_eset,
  method = "hierarchical",
  k = 5,
  target = "gene",              # Options: "gene", "sample"
  hclust_method = "ward.D2",
  distance = "euclidean"
)

clusters <- hc_result$clusters
model <- hc_result$model
```

#### K-means Clustering
```r
km_result <- cluster_expression(
  eset = filtered_eset,
  method = "kmeans",
  k = 5,
  target = "gene",
  nstart = 25
)

# Find optimal k using elbow method
elbow_plot <- plot_elbow(filtered_eset, max_k = 15, target = "gene")
print(elbow_plot)
```

#### DBSCAN Clustering
```r
# Density-based clustering
dbscan_result <- cluster_expression(
  eset = filtered_eset,
  method = "dbscan",
  eps = 0.5,                    # Distance parameter
  minPts = 5,                   # Minimum points per cluster
  target = "gene",
  distance = "euclidean"
)
```

### 4. Cluster Validation

```r
expr_matrix <- exprs(filtered_eset)

# Within-cluster sum of squares
wss <- wss_score(t(expr_matrix), clusters)

# Silhouette analysis
sil <- silhouette_score(t(expr_matrix), clusters)
plot(sil)

# Gap statistic for optimal k
gap_stats <- gap_statistic(t(expr_matrix), max_k = 15, B = 50)
print(gap_stats)

# Parallel gap statistic (faster)
gap_stats_parallel <- gap_statistic_parallel(
  t(expr_matrix), 
  max_k = 15, 
  B = 50,
  BPPARAM = BiocParallel::MulticoreParam(workers = 4)
)
```

### 5. Visualization

```r
# Expression heatmap with cluster annotation
heatmap_plot <- plot_heatmap(
  eset = filtered_eset,
  clusters = clusters,
  scale = "row",                # Options: "none", "row", "column"
  show_row_names = FALSE
)

# PCA plot of samples
pca_plot <- plot_pca(
  eset = filtered_eset,
  color_by = "condition",       # Column in pData(eset)
  ellipse = TRUE
)

# Hierarchical clustering dendrogram
dendrogram <- plot_dendrogram(hc_result$model, main = "Gene Clustering")

# Cluster expression profiles
profile_plot <- plot_cluster_profiles(
  eset = filtered_eset,
  clusters = clusters,
  scale = TRUE
)

# DBSCAN results visualization
dbscan_plot <- plot_dbscan_results(
  eset = filtered_eset,
  clusters = dbscan_result$clusters,
  target = "gene"
)
```

### 6. Gene Set Enrichment Analysis

#### Load Gene Sets
```r
# Load gene sets from online databases (MSigDB)
hallmark_sets <- load_online_gene_sets(
  species = "Homo sapiens",
  collection = "H"              # Hallmark pathways
)

# GO Biological Process
go_bp_sets <- load_online_gene_sets(
  species = "Homo sapiens",
  collection = "C5",
  subcategory = "GO:BP"
)

# Load from GMT file
custom_sets <- load_gmt_gene_sets("path/to/your/file.gmt")
```

#### Perform Enrichment Analysis
```r
# Over-representation analysis
ora_results <- perform_ora(
  clusters = clusters,
  gene_sets = hallmark_sets,
  p_adjust_method = "BH",
  min_set_size = 10,
  max_set_size = 500
)

# Parallel ORA (faster for large datasets)
ora_parallel <- perform_ora_parallel(
  clusters = clusters,
  gene_sets = hallmark_sets,
  BPPARAM = BiocParallel::MulticoreParam(workers = 4)
)

# Summarize top enriched pathways
enrichment_summary <- summarize_ora_results(
  ora_results,
  top_n = 10,
  p_threshold = 0.05,
  min_enrichment_ratio = 1.5
)

# Get detailed results for specific cluster
cluster_1_results <- get_cluster_results(
  ora_results,
  cluster_id = 1,
  p_threshold = 0.05
)
```

## Parallel Processing

The package supports parallel processing for computationally intensive tasks:

```r
library(BiocParallel)

# Configure parallel backend
BPPARAM <- MulticoreParam(workers = 4)

kmeans_results <- kmeans_parallel_k(
  eset = your_eset,
  k_range = 2:10,
  BPPARAM = BPPARAM
)

gap_parallel <- gap_statistic_parallel(
  data = t(exprs(your_eset)),
  max_k = 15,
  BPPARAM = BPPARAM
)

ora_parallel <- perform_ora_parallel(
  clusters = clusters,
  gene_sets = gene_sets,
  BPPARAM = BPPARAM
)
```

## Data Requirements

### Input Data Format
The package expects data in Bioconductor **ExpressionSet** format:

```r
library(Biobase)

# Create ExpressionSet from matrix
expr_matrix <- your_expression_matrix  # genes x samples
pheno_data <- your_sample_info         # sample metadata

eset <- ExpressionSet(
  assayData = expr_matrix,
  phenoData = AnnotatedDataFrame(pheno_data)
)
```

### Gene Expression Matrix
- **Rows**: Genes/Features
- **Columns**: Samples
- **Values**: Expression levels (counts, FPKM, TPM, etc.)
- **Row names**: Gene identifiers (symbols, Ensembl IDs, etc.)
- **Column names**: Sample identifiers

## Output Structure

The package functions return structured results:

```r
# Clustering results
cluster_result <- list(
  model = clustering_model,      # Original clustering object
  clusters = named_vector        # Gene cluster assignments
)

# Enrichment results
enrichment_result <- list(
  cluster_1 = data.frame(...),  # ORA results for cluster 1
  cluster_2 = data.frame(...),  # ORA results for cluster 2
  # ... more clusters
)

# Pipeline results
pipeline_result <- list(
  original_data = eset,
  processed_data = processed_eset,
  outlier_info = outlier_summary,
  clustering_results = all_clustering_results,
  validation_results = validation_metrics,
  plots = list_of_plots,
  enrichment_results = enrichment_results,
  summary = analysis_summary
)
```

## Example Workflows

### Workflow 1: Basic Gene Clustering

```r
results <- run_pipeline(
  data_path = "expression_data.csv",
  data_type = "csv",
  clustering_methods = "kmeans",
  k_clusters = 5,
  perform_enrichment = FALSE,
  create_plots = TRUE
)

head(results$clustering_results$kmeans$clusters)

print(results$plots$pca)
```

### Workflow 2: Comprehensive Analysis with Batch Correction

```r
results <- run_pipeline(
  data_path = "expression_data.csv",
  data_type = "csv",
  pheno_data = "sample_metadata.csv",
  
  # Normalization with batch correction
  norm_method = "quantile",
  batch_correction = "batch",
  
  # Multiple clustering methods
  clustering_methods = c("hierarchical", "kmeans", "dbscan"),
  k_clusters = 6,
  dbscan_eps = 0.5,
  
  # Comprehensive validation
  perform_validation = TRUE,
  gap_statistic = TRUE,
  
  # Pathway enrichment
  perform_enrichment = TRUE,
  gene_sets_source = "online",
  gene_sets_collection = "C2",  # Curated pathways
  
  # Visualization
  pca_color_by = "treatment",
  heatmap_scale = "row",
  
  verbose = TRUE
)

# Compare clustering methods
sapply(results$clustering_results, function(x) length(unique(x$clusters)))

# View enrichment results
head(results$enrichment_results$kmeans$summary)
```

### Workflow 3: Custom Gene Sets Analysis

```r
# Analysis with custom gene sets
custom_pathways <- list(
  "Pathway1" = c("GENE1", "GENE2", "GENE3"),
  "Pathway2" = c("GENE4", "GENE5", "GENE6"),
  "Pathway3" = c("GENE7", "GENE8", "GENE9")
)

results <- run_pipeline(
  data_path = "expression_data.csv",
  data_type = "csv",
  clustering_methods = "hierarchical",
  k_clusters = 4,
  perform_enrichment = TRUE,
  gene_sets_source = "custom",
  custom_gene_sets = custom_pathways
)
```

## Gene Set Databases

The package supports multiple gene set sources:

- **MSigDB Collections**: Hallmark (H), Curated (C2), GO (C5), etc.
- **Custom GMT files**: User-provided gene sets
- **Online databases**: Automatic fetching via `msigdbr`

Supported species include human, mouse, and other model organisms supported by MSigDB.

## Tips and Best Practices

1. **Data Quality**: Always inspect your data and remove outliers before clustering
2. **Normalization**: Choose appropriate normalization based on your data type (counts vs. already normalized)
3. **Clustering Method**: 
   - Use hierarchical clustering for small datasets or when you need a dendrogram
   - Use k-means for medium to large datasets with spherical clusters
   - Use DBSCAN for datasets with irregular cluster shapes or when you expect noise
4. **Cluster Validation**: Use multiple metrics (elbow, gap statistic, silhouette) to determine optimal cluster number
5. **Gene Sets**: Use relevant gene sets for your organism and research context
6. **Parallel Processing**: Enable for large datasets to improve performance
7. **Pipeline Usage**: Use `run_pipeline()` for comprehensive analyses; use individual functions for custom workflows

## Dependencies

### Required Packages
- **Biobase**: ExpressionSet handling
- **cluster**: Clustering algorithms and validation
- **ggplot2**: Plotting
- **ComplexHeatmap**: Advanced heatmaps
- **stringr**: String manipulation

### Optional Packages
- **preprocessCore**: Quantile normalization
- **DESeq2**: DESeq2 normalization
- **sva**: Batch correction
- **dbscan**: DBSCAN clustering
- **msigdbr**: Online gene set databases
- **BiocParallel**: Parallel processing

## License

MIT + file LICENSE

## Support

For questions, issues, or feature requests:
- **Email**: [ghezzigabriele080@gmail.com]