create_test_data <- function() {
  expr_matrix <- matrix(c(
    rnorm(30, mean=10),
    rnorm(30, mean=5),
    rnorm(40, mean=1)
  ), nrow=10, byrow=TRUE)
  colnames(expr_matrix) <- paste0("Sample", 1:10)
  rownames(expr_matrix) <- paste0("Gene", 1:10)
  eset <- Biobase::ExpressionSet(assayData=expr_matrix)
  clusters <- rep(1:3, c(3, 3, 4))
  names(clusters) <- rownames(expr_matrix)
  list(eset=eset, clusters=clusters)
}

testthat::test_that("plot_cluster_profiles validates input correctly", {
  test_data <- create_test_data()
  
  testthat::expect_error(plot_cluster_profiles(matrix(1:10, nrow=2), "exprs"))
               
  bad_clusters <- test_data$clusters
  names(bad_clusters)[1] <- "WrongGene"
  testthat::expect_error(plot_cluster_profiles(test_data$eset, bad_clusters), "Gene names")
})
  
testthat::test_that("plot_cluster_profiles returns a ggplot object", {
  test_data <- create_test_data()
  p <- plot_cluster_profiles(test_data$eset, test_data$clusters)
  
  testthat::expect_s3_class(p, "ggplot")
  testthat::expect_s3_class(p$layers[[1]]$geom, "GeomLine")
  testthat::expect_s3_class(p$layers[[2]]$geom, "GeomPoint")
})

testthat::test_that("scaling works correctly", {
  test_data <- create_test_data()
  
  p_unscaled <- plot_cluster_profiles(test_data$eset, test_data$clusters, scale=FALSE)
  unscaled_data <- ggplot2::ggplot_build(p_unscaled)$plot$data
  
  p_scaled <- plot_cluster_profiles(test_data$eset, test_data$clusters, scale=TRUE)
  scaled_data <- ggplot2::ggplot_build(p_scaled)$plot$data
  
  testthat::expect_false(identical(unscaled_data$Expression, scaled_data$Expression))
  
  expr <- Biobase::exprs(test_data$eset)
  scaled_expr <- t(scale(t(expr)))
  testthat::expect_equal(mean(scaled_expr[1,]), 0, tolerance=1e-10)
  testthat::expect_equal(sd(scaled_expr[1,]), 1, tolerance=1e-10)
})

testthat::test_that("averaging works correctly with known input", {
  # Define known expression matrix
  expr <- matrix(c(
    1, 2, 3,
    4, 5, 6,
    7, 8, 9
  ), nrow = 3, byrow = TRUE)
  colnames(expr) <- c("S1", "S2", "S3")
  rownames(expr) <- c("G1", "G2", "G3")

  eset <- Biobase::ExpressionSet(assayData = expr)

  # Define clusters: G1 and G2 in cluster 1, G3 in cluster 2
  clusters <- c(1, 1, 2)
  names(clusters) <- rownames(expr)

  # Manually compute expected means
  # Cluster 1: mean of G1 and G2
  # Cluster 2: G3 values
  expected <- data.frame(
    Cluster = factor(c(1, 1, 1, 2, 2, 2)),
    Sample = factor(rep(c("S1", "S2", "S3"), each = 2)),
    Expression = c(
      mean(c(1, 4)), mean(c(2, 5)), mean(c(3, 6)),
      7, 8, 9
    )
  )

  # Get actual data
  p <- plot_cluster_profiles(eset, clusters)
  actual <- ggplot2::ggplot_build(p)$plot$data

  # Order both for comparison
  expected <- expected[order(expected$Cluster, expected$Sample), ]
  actual <- actual[order(actual$Cluster, actual$Sample), ]

  testthat::expect_equal(actual$Expression, expected$Expression)
})

testthat::test_that("plot contains correct elements", {
  test_data <- create_test_data()
  p <- plot_cluster_profiles(test_data$eset, test_data$clusters)
  pb <- ggplot2::ggplot_build(p)
  
  testthat::expect_equal(pb$plot$labels$x, "Sample")
  testthat::expect_equal(pb$plot$labels$y, "Average Expression")
  testthat::expect_equal(pb$plot$labels$colour, "Cluster")
  
  testthat::expect_equal(pb$plot$theme$axis.text.x$angle, 45)
  testthat::expect_equal(pb$plot$theme$axis.text.x$hjust, 1)
})

testthat::test_that("function handles edge cases", {
  # Single column
  expr_matrix <- matrix(rnorm(10), nrow=10, ncol=1)
  rownames(expr_matrix) <- paste0("G", 1:10)
  eset <- Biobase::ExpressionSet(assayData=expr_matrix)
  clusters <- rep(1:2, each=5)
  names(clusters) <- rownames(expr_matrix)
  p <- plot_cluster_profiles(eset, clusters)
  testthat::expect_s3_class(p, "ggplot")
})