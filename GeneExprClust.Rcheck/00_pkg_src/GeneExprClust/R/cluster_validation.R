#' Compute Within-Cluster Sum of Squares (WSS)
#'
#' @param data Numeric matrix or data frame of expression data (genes x samples)
#' @param clusters Named vector of cluster assignments
#' @return Total within-cluster sum of squares (numeric)
#' @examples
#' # Create sample gene expression data
#' set.seed(123)
#' data <- matrix(rnorm(200), nrow = 20, ncol = 10)
#' rownames(data) <- paste0("Gene", 1:20)
#' colnames(data) <- paste0("Sample", 1:10)
#' 
#' # Create cluster assignments
#' clusters <- c(rep(1, 8), rep(2, 7), rep(3, 5))
#' names(clusters) <- rownames(data)
#' 
#' # Compute WSS
#' wss_result <- wss_score(data, clusters)
#' message(paste("Total WSS: ", round(wss_result, 2)))
#' @export
wss_score <- function(data, clusters) {
  if (!is.matrix(data) && !is.data.frame(data)) stop("Data must be a matrix or data frame.")
  if (length(clusters) != nrow(data)) stop("Length of cluster vector must match number of rows in data.")
  if (any(is.na(clusters))) warning("Missing cluster assignments detected.")

  total_wss <- 0
  for (cl in unique(clusters)) {
    members <- data[clusters == cl, , drop = FALSE]
    if (nrow(members) < 2) next
    center <- colMeans(members)
    dists <- rowSums((t(t(members) - center))^2)
    total_wss <- total_wss + sum(dists)
  }
  return(total_wss)
}

#' Compute silhouette width for clustering
#'
#' @param data Numeric matrix or data frame (genes x samples)
#' @param clusters Named vector of cluster assignments
#' @param distance Character, distance metric (default = "euclidean")
#' @return Silhouette object
#' @importFrom cluster silhouette
#' @examples
#' # Load required packages
#' library(cluster)
#' library(proxy)
#' 
#' # Create sample data with clear clusters
#' set.seed(456)
#' data <- rbind(
#'   matrix(rnorm(50, mean = 0), nrow = 10, ncol = 5),
#'   matrix(rnorm(50, mean = 3), nrow = 10, ncol = 5),
#'   matrix(rnorm(50, mean = -2), nrow = 10, ncol = 5)
#' )
#' rownames(data) <- paste0("Gene", 1:30)
#' 
#' # Create cluster assignments
#' clusters <- c(rep(1, 10), rep(2, 10), rep(3, 10))
#' names(clusters) <- rownames(data)
#' 
#' # Compute silhouette scores
#' sil_result <- silhouette_score(data, clusters)
#' 
#' # View summary
#' summary(sil_result)
#' 
#' # Plot silhouette
#' plot(sil_result)
#' @export
silhouette_score <- function(data, clusters, distance = "euclidean") {
  if (!is.matrix(data) && !is.data.frame(data)) stop("Data must be a matrix or data frame.")
  if (length(clusters) != nrow(data)) stop("Length of cluster vector must match number of rows in data.")
  if (length(unique(clusters)) < 2) stop("At least 2 clusters are required for silhouette analysis.")

  d <- proxy::dist(data, method = distance)
  sil <- cluster::silhouette(clusters, d)
  return(sil)
}

#' Compute Gap Statistic using k-means
#'
#' @param data Numeric matrix or data frame (genes x samples)
#' @param max_k Maximum number of clusters to try (default 10)
#' @param B Number of bootstrap samples (default 20)
#' @return Data frame with gap statistic per K
#' @importFrom cluster clusGap
#' @examples
#' # Load required package
#' library(cluster)
#' 
#' # Create sample data with natural clustering structure
#' set.seed(789)
#' # Three natural clusters
#' cluster1 <- matrix(rnorm(60, mean = 0, sd = 0.5), nrow = 15, ncol = 4)
#' cluster2 <- matrix(rnorm(60, mean = 2, sd = 0.5), nrow = 15, ncol = 4)
#' cluster3 <- matrix(rnorm(60, mean = -1.5, sd = 0.5), nrow = 15, ncol = 4)
#' 
#' data <- rbind(cluster1, cluster2, cluster3)
#' rownames(data) <- paste0("Gene", 1:45)
#' colnames(data) <- paste0("Sample", 1:4)
#' 
#' # Compute gap statistic
#' gap_result <- gap_statistic(data, max_k = 8, B = 10)
#' 
#' # View results
#' print(gap_result)
#' 
#' # Find optimal k (first local maximum)
#' optimal_k <- which.max(gap_result$gap)
#' message(paste("Suggested optimal k: ", optimal_k))
#' 
#' # Plot gap statistic
#' plot(1:nrow(gap_result), gap_result$gap, type = "b",
#'      xlab = "Number of clusters (k)", ylab = "Gap statistic",
#'      main = "Gap Statistic for Optimal k")
#' abline(v = optimal_k, col = "red", lty = 2)
#' @export
gap_statistic <- function(data, max_k = 10, B = 20) {
  if (!is.matrix(data) && !is.data.frame(data)) stop("Data must be a matrix or data frame.")
  if (max_k < 2) stop("max_k must be at least 2.")
  if (B < 1) stop("Number of bootstrap iterations B must be >= 1.")

  gap <- cluster::clusGap(data, FUN = kmeans, K.max = max_k, B = B)
  return(as.data.frame(gap$Tab))
}

#' Parallel Gap Statistic Computation
#' @param data Numeric matrix or data frame (genes x samples)
#' @param max_k Maximum number of clusters to try (default 10)
#' @param B Number of bootstrap samples (default 20)
#' @param BPPARAM BiocParallel backend (default: bpparam())
#' @return Data frame with gap statistic per K
#' @importFrom BiocParallel bpparam bplapply
#' @importFrom cluster clusGap
#' @importFrom stats runif
#' @examples
#' # Load required packages
#' library(cluster)
#' library(BiocParallel)
#' 
#' # Create larger sample dataset for parallel processing demonstration
#' set.seed(999)
#' # Create 4 natural clusters with more samples
#' n_genes <- 200
#' n_samples <- 20
#' 
#' cluster1 <- matrix(rnorm(n_genes/4 * n_samples, mean = 0), 
#'                   nrow = n_genes/4, ncol = n_samples)
#' cluster2 <- matrix(rnorm(n_genes/4 * n_samples, mean = 2), 
#'                   nrow = n_genes/4, ncol = n_samples)
#' cluster3 <- matrix(rnorm(n_genes/4 * n_samples, mean = -2), 
#'                   nrow = n_genes/4, ncol = n_samples)
#' cluster4 <- matrix(rnorm(n_genes/4 * n_samples, mean = 1.5), 
#'                   nrow = n_genes/4, ncol = n_samples)
#' 
#' large_data <- rbind(cluster1, cluster2, cluster3, cluster4)
#' rownames(large_data) <- paste0("Gene", 1:n_genes)
#' colnames(large_data) <- paste0("Sample", 1:n_samples)
#' 
#' # Set up parallel backend (adjust workers based on your system)
#' param <- MulticoreParam(workers = 2)  # Use 2 cores
#' 
#' # Compute gap statistic in parallel
#' system.time({
#'   gap_parallel_result <- gap_statistic_parallel(
#'     large_data, 
#'     max_k = 8, 
#'     B = 50,  # More bootstrap samples
#'     BPPARAM = param
#'   )
#' })
#' 
#' # View results
#' print(gap_parallel_result)
#' 
#' # Compare with sequential version timing
#' system.time({
#'   gap_sequential_result <- gap_statistic(large_data, max_k = 8, B = 50)
#' })
#' 
#' # Plot comparison
#' plot(1:nrow(gap_parallel_result), gap_parallel_result$gap, type = "b",
#'      xlab = "Number of clusters (k)", ylab = "Gap statistic",
#'      main = "Parallel Gap Statistic Results", col = "blue")
#' @export
gap_statistic_parallel <- function(data, max_k = 10, B = 20, BPPARAM = BiocParallel::bpparam()) {
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("BiocParallel package is required for parallel processing")
  }

  if (!is.matrix(data)) data <- as.matrix(data)
  if (max_k < 2) stop("max_k must be at least 2.")
  if (B < 1) stop("Number of bootstrap iterations B must be >= 1.")
  
  n <- nrow(data)
  dist_obs <- sapply(1:max_k, function(k) {
    cl <- kmeans(data, centers = k, nstart = 10)
    sum(cl$withinss)
  })
  
  ref_data <- replicate(B, apply(data, 2, function(x) runif(n, min(x), max(x))), simplify = FALSE)
  
  dist_ref <- BiocParallel::bplapply(ref_data, function(rd) {
    sapply(1:max_k, function(k) {
      cl <- kmeans(rd, centers = k, nstart = 10)
      sum(cl$withinss)
    })
  }, BPPARAM = BPPARAM)

  dist_ref <- do.call(rbind, dist_ref)
  logW <- log(dist_obs)
  logW_ref <- apply(log(dist_ref), 2, mean)
  gap <- logW_ref - logW
  sk <- sqrt(1 + 1 / B) * apply(log(dist_ref), 2, sd)

  result <- data.frame(k = 1:max_k, gap = gap, sk = sk, logW = logW, logW_ref = logW_ref)
  return(result)
}

