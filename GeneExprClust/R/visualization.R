#' Plot heatmap of expression data with optional clustering annotation and scaling
#'
#' @param eset ExpressionSet object
#' @param clusters Named vector of cluster assignments for genes (optional)
#' @param main The graph's title
#' @param show_row_names Logical, whether to display gene names (default FALSE)
#' @param scale Character, one of "none", "row", or "column" for scaling expression data (default "none")
#' @return Heatmap object
#' @importFrom ComplexHeatmap Heatmap rowAnnotation
#' @importFrom Biobase exprs
#' @examples
#' expr_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' colnames(expr_matrix) <- paste0("Sample", 1:10)
#' rownames(expr_matrix) <- paste0("Gene", 1:100)
#' gene_clusters <- kmeans(expr_matrix, centers = 3)$cluster
#' eset <- Biobase::ExpressionSet(
#'   assayData = expr_matrix
#' )
#' plot_heatmap(eset, clusters = gene_clusters, scale = "row")
#' @export
plot_heatmap <- function(eset, clusters = NULL, main = NULL, show_row_names = FALSE, scale = c("none", "row", "column")) {
  if (!inherits(eset, "ExpressionSet")) {
    stop("Input must be an ExpressionSet object.")
  }

  mat <- Biobase::exprs(eset)
  scale <- match.arg(scale)

  if (scale == "row") mat <- t(scale(t(mat)))
  if (scale == "column") mat <- scale(mat)

  row_anno <- NULL
  if (!is.null(clusters)) {
    if (!is.vector(clusters) || is.null(names(clusters))) {
      warning("Clusters must be a named vector. Ignoring 'clusters'.")
    } else if (!all(rownames(mat) %in% names(clusters))) {
      warning("Some gene names in the expression matrix do not match cluster names. Unmatched genes will be ignored.")
      clusters <- clusters[intersect(rownames(mat), names(clusters))]
    }
    row_anno <- ComplexHeatmap::rowAnnotation(Cluster = factor(clusters))
  }

  ComplexHeatmap::Heatmap(mat,
                          name = "Expression",
                          show_row_names = show_row_names,
                          right_annotation = row_anno,
                          column_title = main)
}

#' PCA plot of samples
#'
#' @param eset ExpressionSet object
#' @param color_by Character, column name in pData(eset) to color samples by
#' @param main Character, main title for the plot (optional)
#' @param ellipse Logical, whether to draw group confidence ellipses (default TRUE)
#' @return ggplot object
#' @importFrom Biobase pData exprs
#' @importFrom stats sd prcomp
#' @import ggplot2
#' @examples
#' # Create example ExpressionSet
#' library(Biobase)
#' 
#' # Expression data (genes x samples)
#' expr_data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' rownames(expr_data) <- paste0("Gene", 1:100)
#' colnames(expr_data) <- paste0("Sample", 1:10)
#' 
#' # Phenotype data
#' pheno_data <- data.frame(
#'   sample_id = paste0("Sample", 1:10),
#'   treatment = rep(c("Control", "Treated"), each = 5),
#'   batch = rep(c("A", "B"), times = 5)
#' )
#' rownames(pheno_data) <- colnames(expr_data)
#' 
#' # Create ExpressionSet
#' eset <- ExpressionSet(
#'   assayData = expr_data,
#'   phenoData = AnnotatedDataFrame(pheno_data)
#' )
#' 
#' # Basic PCA plot
#' plot_pca(eset)
#' 
#' # PCA plot colored by treatment group with ellipses
#' plot_pca(eset, color_by = "treatment", ellipse = TRUE)
#' 
#' # PCA plot colored by batch
#' plot_pca(eset, color_by = "batch", main = "PCA by Batch")
#' @export
plot_pca <- function(eset, color_by = NULL, main = "PCA Plot", ellipse = FALSE) {
  if (!inherits(eset, "ExpressionSet")) {
    stop("Input must be an ExpressionSet.")
  }

  expr <- Biobase::exprs(eset)

  if (any(dim(expr) == 0)) stop("Empty expression matrix")
  if (any(!is.finite(expr))) stop("Non-finite values in expression matrix")
  if (sd(expr) == 0) stop("Zero variance in expression data")

  pca <- tryCatch({
    prcomp(t(expr), scale. = TRUE)
  }, error = function(e) {
    stop("PCA failed: ", e$message)
  })

  percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])

  if (!is.null(color_by)) {
    pd <- Biobase::pData(eset)
    if (color_by %in% colnames(pd)) {
      df$Group <- pd[[color_by]]
    } else {
      stop(sprintf("Column '%s' not found in phenotype data.", color_by))
    }
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)"),
      title = main
    ) +
    ggplot2::theme_minimal()

  if (!is.null(color_by)) {
    p <- p + ggplot2::aes(color = Group)
    if (ellipse) {
      p <- p + ggplot2::stat_ellipse(ggplot2::aes(group = Group), geom = "polygon", alpha = 0.2)
    }
  }

  return(p)
}

#' Plot dendrogram from hierarchical clustering model
#'
#' @param hclust_model Object returned by hclust()
#' @param main Title for the plot (optional)
#' @return ggplot object
#' @importFrom ggdendro ggdendrogram
#' @import ggplot2
#' @examples
#' expr_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' colnames(expr_matrix) <- paste0("Sample", 1:10)
#' rownames(expr_matrix) <- paste0("Gene", 1:100)
#' eset <- Biobase::ExpressionSet(assayData = expr_matrix)
#' hc <- hclust(dist(Biobase::exprs(eset)))
#' p <- plot_dendrogram(hc)
#' message(p)
#' @export
plot_dendrogram <- function(hclust_model, main = "Hierarchical Clustering Dendrogram") {
  if (!requireNamespace("ggdendro", quietly = TRUE)) {
    stop("ggdendro package required for plotting dendrograms")
  }
  
  p <- ggdendro::ggdendrogram(hclust_model, rotate = FALSE, size = 0.5) +
    ggplot2::labs(title = main) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
  
  return(p)
}

#' Plot average gene expression profiles per cluster
#'
#' @param eset ExpressionSet object
#' @param clusters Named vector of cluster assignments for genes
#' @param scale Logical, whether to z-score normalize gene expression before averaging (default FALSE)
#' @return ggplot object showing average expression profile per cluster across samples
#' @importFrom stats aggregate
#' @importFrom Biobase exprs
#' @import ggplot2 reshape2
#' @examples
#' expr_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' colnames(expr_matrix) <- paste0("Sample", 1:10)
#' rownames(expr_matrix) <- paste0("Gene", 1:100)
#' eset <- Biobase::ExpressionSet(assayData = expr_matrix)
#' clusters <- kmeans(expr_matrix, centers = 3)$cluster
#' plot_cluster_profiles(eset, clusters, scale = TRUE)
#' @export
plot_cluster_profiles <- function(eset, clusters, scale = FALSE) {
  if (!inherits(eset, "ExpressionSet")) stop("Input 'eset' must be an ExpressionSet.")
  if (is.null(names(clusters))) stop("Clusters vector must be named (gene names).")

  expr <- Biobase::exprs(eset)
  
  if (!all(rownames(expr) %in% names(clusters))) {
    stop("Gene names in expression set don't match cluster assignments.")
  }
  if (!all(names(clusters) %in% rownames(expr))) {
    stop("Gene names in cluster assignments don't match expression set.")
  }

  if (scale) expr <- t(scale(t(expr)))

  df <- data.frame(expr)
  df$Cluster <- clusters[rownames(df)]

  df_long <- reshape2::melt(df, id.vars = "Cluster", variable.name = "Sample", value.name = "Expression")
  if (nrow(df_long) == 0) stop("No data to plot.")

  agg <- aggregate(Expression ~ Cluster + Sample, df_long, function(x) mean(x, na.rm = TRUE))
  agg$Expression[is.nan(agg$Expression)] <- NA

  ggplot2::ggplot(agg, aes(x = Sample, y = Expression, group = Cluster, color = factor(Cluster))) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(color = "Cluster", x = "Sample", y = "Average Expression")
}

#' Plot DBSCAN clustering results using PCA
#'
#' Visualizes DBSCAN clustering results on genes or samples by performing PCA
#' for dimensionality reduction and plotting the first two principal components.
#' Noise points (cluster 0) are colored differently to distinguish from clusters.
#'
#' @param eset ExpressionSet object containing expression data
#' @param clusters Named vector of cluster assignments from DBSCAN (0 = noise)
#' @param target Character, either "gene" (default) or "sample" to indicate
#'   whether clustering was done on genes or samples
#'
#' @return A ggplot2 object showing the PCA plot colored by cluster assignment
#' @importFrom Biobase exprs
#' @importFrom scales hue_pal
#' @import ggplot2
#' @examples
#' data(eset)
#' # Plot DBSCAN clustering results for genes
#' cluster_results <- c(1, 1, 2, 0, 2, 3)  # Example cluster assignments
#' names(cluster_results) <- rownames(eset)[1:6]
#' plot_dbscan_results(eset, cluster_results, target = "gene")
#' @export
plot_dbscan_results <- function(eset, clusters, target = "gene") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting DBSCAN results")
  }
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase package required")
  }

  mat <- Biobase::exprs(eset)

  if (target == "gene") {
    data_mat <- mat
  } else {
    data_mat <- t(mat)
  }

  # Perform PCA for dimensionality reduction
  pca_res <- prcomp(data_mat, center = TRUE, scale. = TRUE)
  pcs <- as.data.frame(pca_res$x[, 1:2])

  # Use rownames of PCA matrix as the names to subset clusters
  pcs$name <- rownames(pcs)

  # Only keep clusters for the samples present in pcs (alignment)
  pcs$cluster <- factor(clusters[pcs$name], levels = c(0, sort(unique(clusters[clusters != 0]))))

  pcs$cluster_color <- ifelse(pcs$cluster == 0, "Noise", paste("Cluster", pcs$cluster))

  p <- ggplot2::ggplot(pcs, aes(x = PC1, y = PC2, color = cluster_color)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(values = c("Noise" = "gray50",
                                  scales::hue_pal()(length(unique(clusters[clusters != 0]))))) +
    theme_minimal() +
    labs(title = paste("DBSCAN Clustering on", target),
         color = "Cluster") +
    theme(legend.position = "right")

  return(p)
}
