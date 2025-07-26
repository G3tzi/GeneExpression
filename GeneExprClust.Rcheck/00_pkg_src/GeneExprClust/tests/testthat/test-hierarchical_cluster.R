test_that("hierarchical_cluster returns hclust object for gene clustering", {
  mat <- matrix(rnorm(100), nrow = 10)
  eset <- Biobase::ExpressionSet(assayData = mat)
  
  hc <- hierarchical_cluster(eset, target = "gene")
  
  testthat::expect_s3_class(hc, "hclust")
  testthat::expect_equal(length(hc$order), ncol(mat))
})

test_that("hierarchical_cluster returns hclust object for sample clustering", {
  mat <- matrix(rnorm(100), nrow = 10)
  eset <- Biobase::ExpressionSet(assayData = mat)
  
  hc <- hierarchical_cluster(eset, target = "sample")
  
  testthat::expect_s3_class(hc, "hclust")
  testthat::expect_equal(length(hc$order), nrow(mat))
})

test_that("hierarchical_cluster handles different distance and linkage methods", {
  mat <- matrix(rnorm(100), nrow = 10)
  eset <- Biobase::ExpressionSet(assayData = mat)
  
  hc <- hierarchical_cluster(eset, hclust_method = "complete", distance = "manhattan", target = "sample")
  
  testthat::expect_s3_class(hc, "hclust")
})
