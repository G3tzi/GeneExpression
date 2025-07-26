## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6,
  warning = FALSE,
  message = FALSE
)

## ----load_packages------------------------------------------------------------
# Load required packages
library(GeneExprClust)
library(Biobase)
library(ggplot2)
library(ComplexHeatmap)

# Set random seed for reproducibility
set.seed(123)

## ----load_data----------------------------------------------------------------
# Example 1: Loading from CSV file
# eset <- read_expression("path/to/your/data.csv", type = "csv", 
#                        pheno_data = "path/to/sample_metadata.csv")

# Example 2: Creating sample data for demonstration
# Generate example expression data (1000 genes, 20 samples)
n_genes <- 1000
n_samples <- 20

# Create expression matrix with some structure
expr_matrix <- matrix(rnorm(n_genes * n_samples, mean = 8, sd = 2), 
                      nrow = n_genes, ncol = n_samples)

# Add some structured patterns for clustering
# Cluster 1: High expression in first 10 samples
expr_matrix[1:200, 1:10] <- expr_matrix[1:200, 1:10] + 3

# Cluster 2: High expression in last 10 samples  
expr_matrix[201:400, 11:20] <- expr_matrix[201:400, 11:20] + 2.5

# Cluster 3: Oscillating pattern
for(i in 401:600) {
  expr_matrix[i, ] <- expr_matrix[i, ] + 2 * sin(seq(0, 2*pi, length.out = n_samples))
}

# Add gene and sample names
rownames(expr_matrix) <- paste0("Gene_", sprintf("%04d", 1:n_genes))
colnames(expr_matrix) <- paste0("Sample_", sprintf("%02d", 1:n_samples))

# Create sample metadata
sample_data <- data.frame(
  condition = factor(rep(c("Control", "Treatment"), each = 10)),
  time_point = factor(rep(c("T1", "T2", "T3", "T4"), times = 5)),
  batch = factor(rep(c("Batch1", "Batch2"), times = 10)),
  row.names = colnames(expr_matrix)
)

# Create ExpressionSet
eset <- ExpressionSet(
  assayData = expr_matrix,
  phenoData = AnnotatedDataFrame(sample_data)
)

print(eset)

## ----data_overview------------------------------------------------------------
# Basic information about the dataset
message("Number of genes: ", nrow(eset))
message("Number of samples: ", ncol(eset))
message("Sample conditions: ", table(pData(eset)$condition))

# Expression data summary
summary(exprs(eset)[1:5, 1:5])

## ----normalization------------------------------------------------------------
# Apply log2 normalization with optional filtering
eset_norm <- normalize_data(
  eset,
  method = "zscore",
  filter = TRUE,
  filter_threshold = 1,
  filter_min_samples = 3
)

message("Genes after filtering: ", nrow(eset_norm))

# Compare before and after normalization
par(mfrow = c(1, 2))
boxplot(exprs(eset)[, 1:8], main = "Before Normalization",
        las = 2, cex.axis = 0.7)
boxplot(exprs(eset_norm)[, 1:8], main = "After zscore Normalization",
        las = 2, cex.axis = 0.7)

## ----batch_correction, eval=FALSE---------------------------------------------
# # If batch effects are present, correct them
# eset_batch_corrected <- normalize_data(
#   eset_norm,
#   method = "zscore",
#   batch = "batch"  # Column name in phenotype data
# )

## ----pca_analysis-------------------------------------------------------------
# PCA plot to visualize sample relationships
pca_plot <- plot_pca(eset_norm, color_by = "condition",
                     main = "PCA: Samples Colored by Condition")
message(pca_plot)

## ----outlier_detection--------------------------------------------------------
# Detect outlier genes based on variability
eset_clean <- detect_outliers(
  eset_norm,
  method = "mad",
  threshold = 3,
  plot = TRUE,
  remove = TRUE
)

message("Outlier genes detected and removed:", 
    length(attr(eset_clean, "outliers")))

## ----optimal_clusters---------------------------------------------------------
# Elbow plot for k-means
elbow_plot <- plot_elbow(eset_clean, max_k = 15, target = "gene")
message(elbow_plot)

# Gap statistic analysis
gap_results <- gap_statistic(t(exprs(eset_clean)), max_k = 12, B = 20)
print(head(gap_results))

# Plot gap statistic
plot(1:nrow(gap_results), gap_results$gap, type = "b", 
     xlab = "Number of Clusters", ylab = "Gap Statistic",
     main = "Gap Statistic for Optimal k")

## ----hierarchical_clustering--------------------------------------------------
# Perform hierarchical clustering
hc_result <- cluster_expression(
  eset_clean,
  method = "hierarchical",
  k = 5,
  target = "gene",
  hclust_method = "ward.D2",
  distance = "euclidean"
)

# Plot dendrogram
dendrogram_plot <- plot_dendrogram(hc_result$model, main = "Gene Hierarchical Clustering")
message(dendrogram_plot)

# Show cluster sizes
table(hc_result$clusters)

## ----kmeans_clustering--------------------------------------------------------
# Perform k-means clustering
kmeans_result <- cluster_expression(
  eset_clean,
  method = "kmeans",
  k = 5,
  target = "gene",
  nstart = 25
)

# Show cluster sizes
table(kmeans_result$clusters)

# Compare clustering methods
comparison_table <- table(
  Hierarchical = hc_result$clusters,
  KMeans = kmeans_result$clusters
)
print(comparison_table)

## ----cluster_validation-------------------------------------------------------
# Calculate Within-Cluster Sum of Squares (WSS)
wss_hc <- wss_score(exprs(eset_clean), hc_result$clusters)
wss_kmeans <- wss_score(exprs(eset_clean), kmeans_result$clusters)

message("WSS - Hierarchical: ", round(wss_hc, 2))
message("WSS - K-means: ", round(wss_kmeans, 2))

# Calculate Silhouette scores
sil_hc <- silhouette_score(exprs(eset_clean), hc_result$clusters)
sil_kmeans <- silhouette_score(exprs(eset_clean), kmeans_result$clusters)

message("Average Silhouette - Hierarchical: ", round(mean(sil_hc[, 3]), 3))
message("Average Silhouette - K-means: ", round(mean(sil_kmeans[, 3]), 3))

## ----heatmap_visualization----------------------------------------------------
# Create heatmap with cluster annotations
# Use hierarchical clustering results for demonstration
heatmap_plot <- plot_heatmap(
  eset_clean,
  clusters = hc_result$clusters,
  show_row_names = FALSE,
  scale = "row",
  main = "Gene Expression Heatmap with Clusters"
)

print(heatmap_plot)

## ----cluster_profiles---------------------------------------------------------
# Plot average expression profiles per cluster
profile_plot <- plot_cluster_profiles(
  eset_clean,
  clusters = hc_result$clusters,
  scale = TRUE
)
print(profile_plot)

## ----individual_cluster_heatmaps----------------------------------------------
# Focus on specific clusters
for(cluster_id in 1:3) {
  cluster_genes <- names(hc_result$clusters)[hc_result$clusters == cluster_id]
  
  if(length(cluster_genes) > 50) {
    cluster_genes <- sample(cluster_genes, 50)  # Sample for visualization
  }
  
  cluster_eset <- eset_clean[cluster_genes, ]
  
  heatmap_cluster <- plot_heatmap(
    cluster_eset,
    show_row_names = FALSE,
    scale = "row",
    main = paste("Cluster", cluster_id, "- Expression Heatmap")
  )
  
  print(heatmap_cluster)
}

## ----load_gene_sets-----------------------------------------------------------
# Option 1: Load from online databases (MSigDB)
# hallmark_sets <- load_online_gene_sets("Homo sapiens", "H")

# Option 2: Load from GMT file
# custom_sets <- load_gmt_gene_sets("path/to/gene_sets.gmt")

# Option 3: Create example gene sets for demonstration
gene_universe <- rownames(eset_clean)
example_gene_sets <- create_example_gene_sets(
  gene_universe = gene_universe,
  n_pathways = 100,
  pathway_size = 25
)

message("Number of gene sets loaded: ", length(example_gene_sets))
message("Example pathway sizes: ", sapply(example_gene_sets[1:5], length))

## ----perform_ora--------------------------------------------------------------
# Perform ORA analysis for each cluster
ora_results <- perform_ora(
  clusters = hc_result$clusters,
  gene_sets = example_gene_sets,
  universe = gene_universe,
  min_set_size = 5,
  max_set_size = 200,
  p_threshold = 0.05
)

# Summary of results
message("Clusters analyzed: ", length(ora_results))
for(cluster_name in names(ora_results)) {
  n_pathways <- nrow(ora_results[[cluster_name]])
  n_significant <- sum(ora_results[[cluster_name]]$p_adjusted < 0.05)
  message(cluster_name, "- Pathways tested:", n_pathways, 
      ", Significant:", n_significant)
}

## ----summarize_ora------------------------------------------------------------
# Create summary table of top enrichments
enrichment_summary <- summarize_ora_results(
  ora_results,
  top_n = 5,
  p_threshold = 0.05,
  min_enrichment_ratio = 1.5
)

# Display top enrichments
if(nrow(enrichment_summary) > 0) {
  print(head(enrichment_summary, 15))
} else {
  message("No significant enrichments found with current thresholds.\n")
}

## ----detailed_cluster_results-------------------------------------------------
# Get detailed results for cluster 1
if("cluster_1" %in% names(ora_results)) {
  cluster1_results <- get_cluster_results(
    ora_results, 
    cluster_id = 1, 
    p_threshold = 0.1,
    include_genes = FALSE
  )
  
  if(nrow(cluster1_results) > 0) {
    message("Top enriched pathways in Cluster 1:\n")
    print(head(cluster1_results[, c("pathway", "p_adjusted", "enrichment_ratio", 
                                   "genes_overlap", "genes_in_pathway")], 10))
  }
}

## ----cluster_stability--------------------------------------------------------
# Perform clustering with different parameters to assess stability
stability_results <- list()

for(k in 3:7) {
  kmeans_k <- cluster_expression(eset_clean, method = "kmeans", k = k, target = "gene")
  stability_results[[paste0("k", k)]] <- kmeans_k$clusters
}

# Compare cluster assignments across different k values
message("Cluster stability across different k values:\n")
for(i in 1:(length(stability_results)-1)) {
  for(j in (i+1):length(stability_results)) {
    # Calculate adjusted rand index or similar metric
    message("Comparing ", names(stability_results)[i], " vs ", names(stability_results)[j])
  }
}

## ----functional_annotation----------------------------------------------------
# Extract representative genes from each cluster
cluster_representatives <- list()

for(cluster_id in unique(hc_result$clusters)) {
  cluster_genes <- names(hc_result$clusters)[hc_result$clusters == cluster_id]
  
  # Calculate average expression for genes in this cluster
  cluster_expr <- exprs(eset_clean)[cluster_genes, , drop = FALSE]
  gene_means <- rowMeans(cluster_expr)
  
  # Get top expressed genes in this cluster
  top_genes <- names(sort(gene_means, decreasing = TRUE))[1:min(10, length(gene_means))]
  cluster_representatives[[paste0("Cluster_", cluster_id)]] <- top_genes
}

# Display representative genes
for(cluster_name in names(cluster_representatives)) {
  message(cluster_name, " representative genes:\n")
  message(paste(cluster_representatives[[cluster_name]], collapse = ", "), "\n\n")
}

## ----export_results, eval=FALSE-----------------------------------------------
# # Export cluster assignments
# write.csv(data.frame(
#   Gene = names(hc_result$clusters),
#   Cluster = hc_result$clusters
# ), "gene_cluster_assignments.csv", row.names = FALSE)
# 
# # Export enrichment results
# if(nrow(enrichment_summary) > 0) {
#   write.csv(enrichment_summary, "pathway_enrichment_results.csv", row.names = FALSE)
# }
# 
# # Export expression data for clusters
# for(cluster_id in unique(hc_result$clusters)) {
#   cluster_genes <- names(hc_result$clusters)[hc_result$clusters == cluster_id]
#   cluster_data <- exprs(eset_clean)[cluster_genes, ]
# 
#   write.csv(cluster_data, paste0("cluster_", cluster_id, "_expression.csv"))
# }

## ----session_info-------------------------------------------------------------
sessionInfo()

