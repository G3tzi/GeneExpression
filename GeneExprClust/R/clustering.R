#' Hierarchical clustering of genes or samples
#' @param eset ExpressionSet
#' @param hclust_method Linkage method (default: "ward.D2")
#' @param distance Distance metric (default: "euclidean")
#' @param target Cluster by "gene" (default) or "sample"
#' @return hclust object
#' @importFrom stats hclust dist
#' @importFrom Biobase exprs
#' @details If `target = "sample"`, samples (columns) are clustered by transposing the expression matrix. 
#' If `target = "gene"`, gene (rows) are clustered directly.
#' @examples
#' # Load required libraries
#' library(Biobase)
#' 
#' # Create sample ExpressionSet
#' set.seed(123)
#' expr_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' rownames(expr_matrix) <- paste0("Gene_", 1:100)
#' colnames(expr_matrix) <- paste0("Sample_", 1:10)
#' eset <- ExpressionSet(assayData = expr_matrix)
#' 
#' # Hierarchical clustering of genes
#' hc_genes <- hierarchical_cluster(eset, target = "gene")
#' plot(hc_genes, main = "Gene Clustering")
#' 
#' # Hierarchical clustering of samples
#' hc_samples <- hierarchical_cluster(eset, target = "sample", 
#'                                   hclust_method = "complete")
#' plot(hc_samples, main = "Sample Clustering")
#' @export
hierarchical_cluster <- function(eset, hclust_method = "ward.D2", distance = "euclidean", target = "gene") {
  expr <- exprs(eset)
  data <- if (target == "gene") expr else t(expr)
  dist_matrix <- dist(data, method = distance)
  hclust(dist_matrix, method = hclust_method)
}

#' K-means clustering of genes or samples
#' @param eset ExpressionSet
#' @param k Number of clusters
#' @param nstart Number of random starts (default: 25)
#' @param target Cluster by "gene" (default) or "sample"
#' @param ... Other params used by other methods
#' @return kmeans object
#' @importFrom stats kmeans
#' @importFrom Biobase exprs
#' @details If `target = "sample"`, samples (columns) are clustered by transposing the expression matrix. 
#' If `target = "gene"`, gene (rows) are clustered directly.
#' @examples
#' # K-means clustering of genes into 5 clusters
#' data(eset)
#' km_genes <- kmeans_cluster(eset, k = 5, target = "gene")
#' table(km_genes$cluster)
#' @export
kmeans_cluster <- function(eset, k = 5, nstart = 25, target = "gene", ...) {
  expr <- exprs(eset)
  data <- if (target == "gene") expr else t(expr)
  kmeans(data, centers = k, nstart = nstart)
}

#' Parallel K-means with Multiple K Values
#' @param eset ExpressionSet
#' @param k_range Vector of k values to try
#' @param target Cluster by "gene" (default) or "sample"
#' @param BPPARAM BiocParallel backend
#' @return List of kmeans results for each k
#' @examples
#' data(eset)
#' # Parallel k-means for multiple k values on genes
#' library(BiocParallel)
#' results <- kmeans_parallel_k(eset, k_range = 2:8, target = "gene")
#' 
#' # Extract results for k=5
#' k5_result <- results$k_5
#' table(k5_result$cluster)
#' 
#' # Compare total within-cluster sum of squares across k values
#' sapply(results, function(x) x$tot.withinss)
#' @export
kmeans_parallel_k <- function(eset, k_range = 2:10, target = "gene", 
                             BPPARAM = BiocParallel::bpparam()) {
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("BiocParallel package is required")
  }
  
  expr <- exprs(eset)
  data <- if (target == "gene") expr else t(expr)
  
  # Parallel execution across different k values
  results <- BiocParallel::bplapply(k_range, function(k) {
    kmeans(data, centers = k, nstart = 25)
  }, BPPARAM = BPPARAM)
  
  names(results) <- paste0("k_", k_range)
  return(results)
}

#' Cluster genes using DBSCAN
#'
#' @param eset ExpressionSet object
#' @param eps Numeric, DBSCAN epsilon parameter (required)
#' @param minPts Integer, minimum number of points per cluster (default 5)
#' @param distance Character, distance metric ("euclidean" or "correlation")
#' @param scale Logical, whether to z-score normalize expression values (default TRUE)
#' @param target Cluster "gene" (rows) or "sample" (columns)
#' @return List containing:
#' \itemize{
#'   \item cluster - Named vector of cluster assignments (0 = noise)
#'   \item eps - The eps parameter used
#'   \item minPts - The minPts parameter used
#' }
#' @importFrom Biobase exprs
#' @importFrom dbscan dbscan
#' @importFrom stats cor complete.cases dist as.dist
#' @examples
#' # DBSCAN clustering of genes
#' data(eset)
#' db_genes <- dbscan_cluster(eset, eps = 2.5, minPts = 3, target = "gene")
#' table(db_genes$cluster)  # 0 indicates noise points
#' 
#' # DBSCAN clustering of samples with correlation distance
#' db_samples <- dbscan_cluster(eset, eps = 0.3, minPts = 2, 
#'                             distance = "correlation", target = "sample")
#' print(db_samples$cluster)
#' 
#' # Check for noise points (cluster = 0)
#' noise_genes <- names(db_genes$cluster)[db_genes$cluster == 0]
#' message(paste("Number of noise genes: ", length(noise_genes)))
#' @export
dbscan_cluster <- function(
    eset, 
    eps, 
    minPts = 5, 
    distance = "euclidean", 
    scale = TRUE,
    target = "gene"
) {
  if (!inherits(eset, "ExpressionSet")) stop("Input must be an ExpressionSet.")
  if (!is.numeric(eps) || eps <= 0) stop("'eps' must be a positive numeric value.")
  if (!is.numeric(minPts) || minPts < 1) stop("'minPts' must be a positive integer.")
  if (!target %in% c("gene", "sample")) stop('target must be "gene" or "sample".')
  if (!distance %in% c("euclidean", "correlation")) {
    stop("Distance must be either 'euclidean' or 'correlation'")
  }
  
  mat <- Biobase::exprs(eset)
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    stop("Expression matrix must have at least 2 genes and 2 samples.")
  }

  if (target == "sample") mat <- t(mat)

  if (scale) {
    mat <- scale(mat)
    if (any(is.na(mat))) warning("Scaling introduced NA values - these will be removed.")
    mat <- mat[complete.cases(mat), ]
    if (nrow(mat) < 2) stop("Not enough complete cases after scaling.")
  }

  if (distance == "euclidean") {
    dist_matrix <- dist(mat)
  } else {
    dist_matrix <- as.dist(1 - cor(t(mat)))
  }

  clustering <- dbscan::dbscan(dist_matrix, eps = eps, minPts = minPts)

  labels <- clustering$cluster
  names(labels) <- if (target == "gene") rownames(mat) else colnames(eset)
  
  if (all(labels == 0)) warning("DBSCAN found no clusters - all points are noise.")
  
  return(list(
    cluster = labels,
    eps = eps,
    minPts = minPts,
    distance = distance,
    target = target
  ))
}

#' Unified clustering function with multiple methods
#' @param eset ExpressionSet
#' @param method One of "hierarchical", "kmeans", "dbscan"
#' @param k Number of clusters (used for hierarchical and kmeans)
#' @param target Cluster by "gene" (default) or "sample"
#' @param ... Additional parameters passed to method
#' @return A list with model object and cluster assignments
#' @importFrom stats cutree
#' @details If `target = "sample"`, samples (columns) are clustered by transposing the expression matrix. 
#' If `target = "gene"`, gene (rows) are clustered directly.
#' The `cluster` vector in the returned object is named by gene or sample names accordingly.
#' @examples
#' # Hierarchical clustering with automatic cluster assignment
#' data(eset)
#' hc_result <- cluster_expression(eset, method = "hierarchical", k = 4)
#' table(hc_result$clusters)
#' 
#' # K-means clustering
#' km_result <- cluster_expression(eset, method = "kmeans", k = 3, target = "sample")
#' print(km_result$clusters)
#' 
#' # DBSCAN clustering
#' db_result <- cluster_expression(eset, method = "dbscan", eps = 2.0, minPts = 5)
#' summary(as.factor(db_result$clusters))
#' @export
cluster_expression <- function(eset, method = c("hierarchical", "kmeans", "dbscan"), 
                               k = 3, target = "gene", ...) {
  method <- match.arg(method)
  if (method == "hierarchical") {
    model <- hierarchical_cluster(eset, target = target, ...)
    clusters <- cutree(model, k = k)
  } else if (method == "kmeans") {
    model <- kmeans_cluster(eset, k = k, target = target, ...)
    clusters <- model$cluster
  } else if (method == "dbscan") {
    model <- dbscan_cluster(eset, target = target, ...)
    clusters <- model$cluster
  }

  expr <- exprs(eset)
  if (target == "gene") {
    names(clusters) <- rownames(expr)
  } else {
    names(clusters) <- colnames(expr)
  }

  return(list(model = model, clusters = clusters))
}


#' Elbow plot for selecting number of clusters in k-means
#' @param eset ExpressionSet
#' @param max_k Maximum number of clusters to try
#' @param target Cluster by "gene" (default) or "sample"
#' @return A ggplot object
#' @importFrom Biobase exprs
#' @import ggplot2
#' @details Computes total within-cluster sum of squares (WSS) for k = 1 to max_k.
#' Plots WSS to help visually select the optimal number of clusters using the elbow method.
#' If `target = "gene"`, genes are clustered by transposing the expression matrix; 
#' otherwise samples are clustered directly.
#' @examples
#' data(eset)
#' # Generate elbow plot for gene clustering
#' elbow_plot <- plot_elbow(eset, max_k = 9, target = "gene")
#' message(elbow_plot)
#' 
#' # Generate elbow plot for sample clustering
#' elbow_plot_samples <- plot_elbow(eset, max_k = 8, target = "sample")
#' message(elbow_plot_samples)
#' 
#' # The "elbow" in the plot suggests optimal number of clusters
#' @export
plot_elbow <- function(eset, max_k = 10, target = "gene") {
  expr <- exprs(eset)
  data <- if (target == "gene") t(expr) else expr

  wss <- sapply(1:max_k, function(k) {
    kmeans(data, centers = k, nstart = 25)$tot.withinss
  })

  df <- data.frame(k = 1:max_k, WSS = wss)

  ggplot(df, aes(x = k, y = WSS)) +
    geom_line() +
    geom_point() +
    theme_minimal() +
    labs(title = "Elbow Method for Optimal k", x = "Number of Clusters", y = "Within-Cluster SS")
}
