test_that("cluster_expression returns expected structure for hierarchical", {
  mat <- matrix(rnorm(100), nrow = 10)
  eset <- Biobase::ExpressionSet(assayData = mat)
  
  res <- cluster_expression(eset, method = "hierarchical", k = 3, target = "gene")
  
  testthat::expect_type(res, "list")
  testthat::expect_s3_class(res$model, "hclust")
  testthat::expect_named(res, c("model", "clusters"))
  testthat::expect_equal(length(res$clusters), ncol(mat))
})

test_that("cluster_expression returns expected structure for kmeans", {
  mat <- matrix(rnorm(100), nrow = 10)
  eset <- Biobase::ExpressionSet(assayData = mat)
  
  res <- cluster_expression(eset, method = "kmeans", k = 4, target = "sample")
  
  testthat::expect_type(res, "list")
  testthat::expect_s3_class(res$model, "kmeans")
  testthat::expect_named(res, c("model", "clusters"))
  testthat::expect_equal(length(res$clusters), nrow(mat))
  testthat::expect_equal(length(unique(res$clusters)), 4)
})

test_that("cluster_expression works with dbscan", {
  # Create better test data with clear clusters
  set.seed(123)
  cluster1 <- matrix(rnorm(20, mean = 0, sd = 0.5), nrow = 5)
  cluster2 <- matrix(rnorm(20, mean = 5, sd = 0.5), nrow = 5)
  cluster3 <- matrix(rnorm(20, mean = 10, sd = 0.5), nrow = 5)
  mat <- rbind(cluster1, cluster2, cluster3)
  rownames(mat) <- paste0("Gene", 1:15)
  colnames(mat) <- paste0("Sample", 1:4)
  eset <- Biobase::ExpressionSet(assayData = mat)
  
  # Test with appropriate parameters
  res <- cluster_expression(
    eset, 
    method = "dbscan", 
    eps = 1.5,    # Increased from 2 to better match the data scale
    minPts = 3,   # Reasonable minimum cluster size
    target = "gene",
    distance = "euclidean"
  )
  
  testthat::expect_type(res, "list")
  testthat::expect_named(res, c("model", "clusters"))
  testthat::expect_equal(length(res$clusters), nrow(mat))
  testthat::expect_true(max(res$clusters) > 0, 
             info = "DBSCAN should find at least one cluster with this data")
})

test_that("cluster_expression errors with invalid method", {
  mat <- matrix(rnorm(100), nrow = 10)
  eset <- Biobase::ExpressionSet(assayData = mat)
  
  testthat::expect_error(cluster_expression(eset, method = "invalid"), regexp = "arg")
})

test_that("cluster_expression passes additional arguments", {
  mat <- matrix(rnorm(100), nrow = 10)
  eset <- Biobase::ExpressionSet(assayData = mat)
  
  res <- cluster_expression(eset, method = "kmeans", k = 3, nstart = 10)
  testthat::expect_s3_class(res$model, "kmeans")
})
