test_that("kmeans_cluster returns kmeans object for gene clustering", {
  mat <- matrix(rnorm(100), nrow = 10)
  eset <- Biobase::ExpressionSet(assayData = mat)
  
  km <- kmeans_cluster(eset, k = 3, target = "gene")
  
  testthat::expect_s3_class(km, "kmeans")
  testthat::expect_equal(length(km$cluster), ncol(mat))
  testthat::expect_equal(km$centers %>% nrow(), 3)
})

test_that("kmeans_cluster returns kmeans object for sample clustering", {
  mat <- matrix(rnorm(100), nrow = 10)
  eset <- Biobase::ExpressionSet(assayData = mat)
  
  km <- kmeans_cluster(eset, k = 4, target = "sample")
  
  testthat::expect_s3_class(km, "kmeans")
  testthat::expect_equal(length(km$cluster), nrow(mat))
  testthat::expect_equal(nrow(km$centers), 4)
})

test_that("kmeans_cluster respects nstart parameter", {
  mat <- matrix(rnorm(100), nrow = 10)
  eset <- Biobase::ExpressionSet(assayData = mat)
  
  # Just run to ensure no error; difficult to test nstart effect directly
  testthat::expect_s3_class(kmeans_cluster(eset, k = 3, nstart = 50), "kmeans")
})
