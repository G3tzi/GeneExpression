pkgname <- "GeneExprClust"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('GeneExprClust')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("cluster_expression")
### * cluster_expression

flush(stderr()); flush(stdout())

### Name: cluster_expression
### Title: Unified clustering function with multiple methods
### Aliases: cluster_expression

### ** Examples

# Hierarchical clustering with automatic cluster assignment
data(eset)
hc_result <- cluster_expression(eset, method = "hierarchical", k = 4)
table(hc_result$clusters)

# K-means clustering
km_result <- cluster_expression(eset, method = "kmeans", k = 3, target = "sample")
print(km_result$clusters)

# DBSCAN clustering
db_result <- cluster_expression(eset, method = "dbscan", eps = 2.0, minPts = 5)
summary(as.factor(db_result$clusters))



cleanEx()
nameEx("create_example_gene_sets")
### * create_example_gene_sets

flush(stderr()); flush(stdout())

### Name: create_example_gene_sets
### Title: Create Simple Example Gene Sets
### Aliases: create_example_gene_sets

### ** Examples

data(eset)
# Create example gene sets
example_sets <- create_example_gene_sets(rownames(eset), n_pathways = 20)



cleanEx()
nameEx("dbscan_cluster")
### * dbscan_cluster

flush(stderr()); flush(stdout())

### Name: dbscan_cluster
### Title: Cluster genes using DBSCAN
### Aliases: dbscan_cluster

### ** Examples

# DBSCAN clustering of genes
data(eset)
db_genes <- dbscan_cluster(eset, eps = 2.5, minPts = 3, target = "gene")
table(db_genes$cluster)  # 0 indicates noise points

# DBSCAN clustering of samples with correlation distance
db_samples <- dbscan_cluster(eset, eps = 0.3, minPts = 2, 
                            distance = "correlation", target = "sample")
print(db_samples$cluster)

# Check for noise points (cluster = 0)
noise_genes <- names(db_genes$cluster)[db_genes$cluster == 0]
message(paste("Number of noise genes: ", length(noise_genes)))



cleanEx()
nameEx("detect_outliers")
### * detect_outliers

flush(stderr()); flush(stdout())

### Name: detect_outliers
### Title: Detect and optionally remove outlier genes based on variability
### Aliases: detect_outliers

### ** Examples

data(eset)
# Detect and remove outlier genes using MAD method
eset_filtered <- detect_outliers(eset, plot = TRUE)



cleanEx()
nameEx("eset")
### * eset

flush(stderr()); flush(stdout())

### Name: eset
### Title: Example Expression Set
### Aliases: eset
### Keywords: datasets

### ** Examples

data(eset)
str(eset)



cleanEx()
nameEx("extract_gene_names")
### * extract_gene_names

flush(stderr()); flush(stdout())

### Name: extract_gene_names
### Title: Extract Gene Names from ExpressionSet
### Aliases: extract_gene_names

### ** Examples

# Create a sample ExpressionSet
library(Biobase)

# Create expression matrix
expr_matrix <- matrix(rnorm(100), nrow = 10, ncol = 10)
rownames(expr_matrix) <- paste0("probe_", 1:10)
colnames(expr_matrix) <- paste0("sample_", 1:10)

# Create feature data with gene symbols
feature_data <- data.frame(
  SYMBOL = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", 
             "GENE6", "GENE7", "GENE8", "GENE9", "GENE10"),
  ENTREZID = paste0("ID_", 1:10),
  stringsAsFactors = FALSE
)
rownames(feature_data) <- rownames(expr_matrix)

# Create phenotype data
pheno_data <- data.frame(
  condition = rep(c("control", "treatment"), each = 5),
  stringsAsFactors = FALSE
)
rownames(pheno_data) <- colnames(expr_matrix)

# Create ExpressionSet
eset <- ExpressionSet(
  assayData = expr_matrix,
  phenoData = AnnotatedDataFrame(pheno_data),
  featureData = AnnotatedDataFrame(feature_data)
)

# Extract gene names using different methods
probe_ids <- extract_gene_names(eset)  # Default: featureNames
gene_symbols <- extract_gene_names(eset, id_type = "SYMBOL")
entrez_ids <- extract_gene_names(eset, id_type = "ENTREZID")



cleanEx()
nameEx("gap_statistic")
### * gap_statistic

flush(stderr()); flush(stdout())

### Name: gap_statistic
### Title: Compute Gap Statistic using k-means
### Aliases: gap_statistic

### ** Examples

# Load required package
library(cluster)

# Create sample data with natural clustering structure
set.seed(789)
# Three natural clusters
cluster1 <- matrix(rnorm(60, mean = 0, sd = 0.5), nrow = 15, ncol = 4)
cluster2 <- matrix(rnorm(60, mean = 2, sd = 0.5), nrow = 15, ncol = 4)
cluster3 <- matrix(rnorm(60, mean = -1.5, sd = 0.5), nrow = 15, ncol = 4)

data <- rbind(cluster1, cluster2, cluster3)
rownames(data) <- paste0("Gene", 1:45)
colnames(data) <- paste0("Sample", 1:4)

# Compute gap statistic
gap_result <- gap_statistic(data, max_k = 8, B = 10)

# View results
print(gap_result)

# Find optimal k (first local maximum)
optimal_k <- which.max(gap_result$gap)
message(paste("Suggested optimal k: ", optimal_k))

# Plot gap statistic
plot(1:nrow(gap_result), gap_result$gap, type = "b",
     xlab = "Number of clusters (k)", ylab = "Gap statistic",
     main = "Gap Statistic for Optimal k")
abline(v = optimal_k, col = "red", lty = 2)



cleanEx()
nameEx("gap_statistic_parallel")
### * gap_statistic_parallel

flush(stderr()); flush(stdout())

### Name: gap_statistic_parallel
### Title: Parallel Gap Statistic Computation
### Aliases: gap_statistic_parallel

### ** Examples

# Load required packages
library(cluster)
library(BiocParallel)

# Create larger sample dataset for parallel processing demonstration
set.seed(999)
# Create 4 natural clusters with more samples
n_genes <- 200
n_samples <- 20

cluster1 <- matrix(rnorm(n_genes/4 * n_samples, mean = 0), 
                  nrow = n_genes/4, ncol = n_samples)
cluster2 <- matrix(rnorm(n_genes/4 * n_samples, mean = 2), 
                  nrow = n_genes/4, ncol = n_samples)
cluster3 <- matrix(rnorm(n_genes/4 * n_samples, mean = -2), 
                  nrow = n_genes/4, ncol = n_samples)
cluster4 <- matrix(rnorm(n_genes/4 * n_samples, mean = 1.5), 
                  nrow = n_genes/4, ncol = n_samples)

large_data <- rbind(cluster1, cluster2, cluster3, cluster4)
rownames(large_data) <- paste0("Gene", 1:n_genes)
colnames(large_data) <- paste0("Sample", 1:n_samples)

# Set up parallel backend (adjust workers based on your system)
param <- MulticoreParam(workers = 2)  # Use 2 cores

# Compute gap statistic in parallel
system.time({
  gap_parallel_result <- gap_statistic_parallel(
    large_data, 
    max_k = 8, 
    B = 50,  # More bootstrap samples
    BPPARAM = param
  )
})

# View results
print(gap_parallel_result)

# Compare with sequential version timing
system.time({
  gap_sequential_result <- gap_statistic(large_data, max_k = 8, B = 50)
})

# Plot comparison
plot(1:nrow(gap_parallel_result), gap_parallel_result$gap, type = "b",
     xlab = "Number of clusters (k)", ylab = "Gap statistic",
     main = "Parallel Gap Statistic Results", col = "blue")



cleanEx()
nameEx("get_cluster_results")
### * get_cluster_results

flush(stderr()); flush(stdout())

### Name: get_cluster_results
### Title: Get Detailed Results for Specific Cluster
### Aliases: get_cluster_results

### ** Examples

## Not run: 
##D # Get results for specific cluster
##D cluster1_results <- get_cluster_results(ora_results, cluster_id = 1)
## End(Not run)



cleanEx()
nameEx("hierarchical_cluster")
### * hierarchical_cluster

flush(stderr()); flush(stdout())

### Name: hierarchical_cluster
### Title: Hierarchical clustering of genes or samples
### Aliases: hierarchical_cluster

### ** Examples

# Load required libraries
library(Biobase)

# Create sample ExpressionSet
set.seed(123)
expr_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
rownames(expr_matrix) <- paste0("Gene_", 1:100)
colnames(expr_matrix) <- paste0("Sample_", 1:10)
eset <- ExpressionSet(assayData = expr_matrix)

# Hierarchical clustering of genes
hc_genes <- hierarchical_cluster(eset, target = "gene")
plot(hc_genes, main = "Gene Clustering")

# Hierarchical clustering of samples
hc_samples <- hierarchical_cluster(eset, target = "sample", 
                                  hclust_method = "complete")
plot(hc_samples, main = "Sample Clustering")



cleanEx()
nameEx("kmeans_cluster")
### * kmeans_cluster

flush(stderr()); flush(stdout())

### Name: kmeans_cluster
### Title: K-means clustering of genes or samples
### Aliases: kmeans_cluster

### ** Examples

# K-means clustering of genes into 5 clusters
data(eset)
km_genes <- kmeans_cluster(eset, k = 5, target = "gene")
table(km_genes$cluster)



cleanEx()
nameEx("kmeans_parallel_k")
### * kmeans_parallel_k

flush(stderr()); flush(stdout())

### Name: kmeans_parallel_k
### Title: Parallel K-means with Multiple K Values
### Aliases: kmeans_parallel_k

### ** Examples

data(eset)
# Parallel k-means for multiple k values on genes
library(BiocParallel)
results <- kmeans_parallel_k(eset, k_range = 2:8, target = "gene")

# Extract results for k=5
k5_result <- results$k_5
table(k5_result$cluster)

# Compare total within-cluster sum of squares across k values
sapply(results, function(x) x$tot.withinss)



cleanEx()
nameEx("load_gmt_gene_sets")
### * load_gmt_gene_sets

flush(stderr()); flush(stdout())

### Name: load_gmt_gene_sets
### Title: Load Gene Sets from GMT File
### Aliases: load_gmt_gene_sets

### ** Examples

## Not run: 
##D # Load gene sets from GMT file
##D gene_sets <- load_gmt_gene_sets("pathways.gmt")
## End(Not run)



cleanEx()
nameEx("load_online_gene_sets")
### * load_online_gene_sets

flush(stderr()); flush(stdout())

### Name: load_online_gene_sets
### Title: Load Gene Sets from Online Databases
### Aliases: load_online_gene_sets

### ** Examples

# Get Hallmark pathways for human
hallmark_sets <- load_online_gene_sets("Homo sapiens", "H")
# Load Hallmark pathways for human
hallmark_sets <- load_online_gene_sets("Homo sapiens", "H")



cleanEx()
nameEx("normalize_data")
### * normalize_data

flush(stderr()); flush(stdout())

### Name: normalize_data
### Title: Normalize expression data in an ExpressionSet with optional
###   batch correction and filtering
### Aliases: normalize_data

### ** Examples

data(eset)
# Log2 normalization
eset_log2 <- normalize_data(eset, method = "log2")



cleanEx()
nameEx("ora_from_expression_set")
### * ora_from_expression_set

flush(stderr()); flush(stdout())

### Name: ora_from_expression_set
### Title: Perform ORA Analysis from ExpressionSet and Clusters
### Aliases: ora_from_expression_set

### ** Examples

## Not run: 
##D # ORA analysis from ExpressionSet
##D ora_results <- ora_from_expression_set(eset, clusters, gene_sets)
## End(Not run)



cleanEx()
nameEx("perform_ora")
### * perform_ora

flush(stderr()); flush(stdout())

### Name: perform_ora
### Title: Perform Over-Representation Analysis (ORA) for Gene Clusters
### Aliases: perform_ora

### ** Examples

# Basic ORA analysis
clusters <- c(1, 1, 2, 2, 3, 3, 1, 2)
names(clusters) <- c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", 
                     "GENE6", "GENE7", "GENE8")
gene_sets <- list(
  Pathway1 = c("GENE1", "GENE2", "GENE7", "GENE9", "GENE10"),
  Pathway2 = c("GENE3", "GENE4", "GENE8", "GENE11", "GENE12")
)
ora_results <- perform_ora(clusters, gene_sets)



cleanEx()
nameEx("perform_ora_parallel")
### * perform_ora_parallel

flush(stderr()); flush(stdout())

### Name: perform_ora_parallel
### Title: Parallel Over-Representation Analysis (ORA) for Gene Clusters
### Aliases: perform_ora_parallel

### ** Examples

## Not run: 
##D # Parallel ORA analysis
##D ora_results <- perform_ora_parallel(clusters, gene_sets)
## End(Not run)



cleanEx()
nameEx("plot_cluster_profiles")
### * plot_cluster_profiles

flush(stderr()); flush(stdout())

### Name: plot_cluster_profiles
### Title: Plot average gene expression profiles per cluster
### Aliases: plot_cluster_profiles

### ** Examples

expr_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
colnames(expr_matrix) <- paste0("Sample", 1:10)
rownames(expr_matrix) <- paste0("Gene", 1:100)
eset <- Biobase::ExpressionSet(assayData = expr_matrix)
clusters <- kmeans(expr_matrix, centers = 3)$cluster
plot_cluster_profiles(eset, clusters, scale = TRUE)



cleanEx()
nameEx("plot_dbscan_results")
### * plot_dbscan_results

flush(stderr()); flush(stdout())

### Name: plot_dbscan_results
### Title: Plot DBSCAN clustering results using PCA
### Aliases: plot_dbscan_results

### ** Examples

data(eset)
# Plot DBSCAN clustering results for genes
cluster_results <- c(1, 1, 2, 0, 2, 3)  # Example cluster assignments
names(cluster_results) <- rownames(eset)[1:6]
plot_dbscan_results(eset, cluster_results, target = "gene")



cleanEx()
nameEx("plot_dendrogram")
### * plot_dendrogram

flush(stderr()); flush(stdout())

### Name: plot_dendrogram
### Title: Plot dendrogram from hierarchical clustering model
### Aliases: plot_dendrogram

### ** Examples

expr_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
colnames(expr_matrix) <- paste0("Sample", 1:10)
rownames(expr_matrix) <- paste0("Gene", 1:100)
eset <- Biobase::ExpressionSet(assayData = expr_matrix)
hc <- hclust(dist(Biobase::exprs(eset)))
p <- plot_dendrogram(hc)
message(p)



cleanEx()
nameEx("plot_elbow")
### * plot_elbow

flush(stderr()); flush(stdout())

### Name: plot_elbow
### Title: Elbow plot for selecting number of clusters in k-means
### Aliases: plot_elbow

### ** Examples

data(eset)
# Generate elbow plot for gene clustering
elbow_plot <- plot_elbow(eset, max_k = 9, target = "gene")
message(elbow_plot)

# Generate elbow plot for sample clustering
elbow_plot_samples <- plot_elbow(eset, max_k = 8, target = "sample")
message(elbow_plot_samples)

# The "elbow" in the plot suggests optimal number of clusters



cleanEx()
nameEx("plot_heatmap")
### * plot_heatmap

flush(stderr()); flush(stdout())

### Name: plot_heatmap
### Title: Plot heatmap of expression data with optional clustering
###   annotation and scaling
### Aliases: plot_heatmap

### ** Examples

expr_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
colnames(expr_matrix) <- paste0("Sample", 1:10)
rownames(expr_matrix) <- paste0("Gene", 1:100)
gene_clusters <- kmeans(expr_matrix, centers = 3)$cluster
eset <- Biobase::ExpressionSet(
  assayData = expr_matrix
)
plot_heatmap(eset, clusters = gene_clusters, scale = "row")



cleanEx()
nameEx("plot_pca")
### * plot_pca

flush(stderr()); flush(stdout())

### Name: plot_pca
### Title: PCA plot of samples
### Aliases: plot_pca

### ** Examples

# Create example ExpressionSet
library(Biobase)

# Expression data (genes x samples)
expr_data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
rownames(expr_data) <- paste0("Gene", 1:100)
colnames(expr_data) <- paste0("Sample", 1:10)

# Phenotype data
pheno_data <- data.frame(
  sample_id = paste0("Sample", 1:10),
  treatment = rep(c("Control", "Treated"), each = 5),
  batch = rep(c("A", "B"), times = 5)
)
rownames(pheno_data) <- colnames(expr_data)

# Create ExpressionSet
eset <- ExpressionSet(
  assayData = expr_data,
  phenoData = AnnotatedDataFrame(pheno_data)
)

# Basic PCA plot
plot_pca(eset)

# PCA plot colored by treatment group with ellipses
plot_pca(eset, color_by = "treatment", ellipse = TRUE)

# PCA plot colored by batch
plot_pca(eset, color_by = "batch", main = "PCA by Batch")



cleanEx()
nameEx("read_expression")
### * read_expression

flush(stderr()); flush(stdout())

### Name: read_expression
### Title: Read gene expression data from file
### Aliases: read_expression

### ** Examples

# Load existing ExpressionSet from RData file
file_path <- system.file("data", "eset.rda", package = "GeneExprClust")
eset <- read_expression(file_path, type = "eset")

# Check the resulting ExpressionSet
print(eset)
dim(Biobase::exprs(eset))  # Get expression matrix dimensions
Biobase::pData(eset)       # Get phenotype data



cleanEx()
nameEx("run_pipeline")
### * run_pipeline

flush(stderr()); flush(stdout())

### Name: run_pipeline
### Title: Complete Gene Expression Analysis Pipeline
### Aliases: run_pipeline

### ** Examples

data(eset)
results <- run_pipeline(
  expression_set = eset,
  clustering_methods = c("hierarchical", "kmeans"),
  k_clusters = 4
)



cleanEx()
nameEx("silhouette_score")
### * silhouette_score

flush(stderr()); flush(stdout())

### Name: silhouette_score
### Title: Compute silhouette width for clustering
### Aliases: silhouette_score

### ** Examples

# Load required packages
library(cluster)
library(proxy)

# Create sample data with clear clusters
set.seed(456)
data <- rbind(
  matrix(rnorm(50, mean = 0), nrow = 10, ncol = 5),
  matrix(rnorm(50, mean = 3), nrow = 10, ncol = 5),
  matrix(rnorm(50, mean = -2), nrow = 10, ncol = 5)
)
rownames(data) <- paste0("Gene", 1:30)

# Create cluster assignments
clusters <- c(rep(1, 10), rep(2, 10), rep(3, 10))
names(clusters) <- rownames(data)

# Compute silhouette scores
sil_result <- silhouette_score(data, clusters)

# View summary
summary(sil_result)

# Plot silhouette
plot(sil_result)



cleanEx()
nameEx("summarize_ora_results")
### * summarize_ora_results

flush(stderr()); flush(stdout())

### Name: summarize_ora_results
### Title: Summarize ORA Results
### Aliases: summarize_ora_results

### ** Examples

## Not run: 
##D # Summarize ORA results
##D summary_table <- summarize_ora_results(ora_results, top_n = 5)
## End(Not run)



cleanEx()
nameEx("wss_score")
### * wss_score

flush(stderr()); flush(stdout())

### Name: wss_score
### Title: Compute Within-Cluster Sum of Squares (WSS)
### Aliases: wss_score

### ** Examples

# Create sample gene expression data
set.seed(123)
data <- matrix(rnorm(200), nrow = 20, ncol = 10)
rownames(data) <- paste0("Gene", 1:20)
colnames(data) <- paste0("Sample", 1:10)

# Create cluster assignments
clusters <- c(rep(1, 8), rep(2, 7), rep(3, 5))
names(clusters) <- rownames(data)

# Compute WSS
wss_result <- wss_score(data, clusters)
message(paste("Total WSS: ", round(wss_result, 2)))



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
