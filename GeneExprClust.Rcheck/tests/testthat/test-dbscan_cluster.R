test_that("dbscan_cluster returns dbscan object for gene clustering", {
  mat <- matrix(rnorm(100), nrow = 10)
  rownames(mat) <- paste0("Gene", 1:10)
  colnames(mat) <- paste0("Sample", 1:10)

  eset <- Biobase::ExpressionSet(assayData = mat)
  
  res <- dbscan_cluster(eset, eps = 0.7, minPts = 3, target = "gene")
  
  testthat::expect_type(res, "list")
  testthat::expect_named(res, c("cluster", "eps", "minPts", "distance", "target"))
  testthat::expect_equal(length(res$cluster), nrow(mat))
  testthat::expect_equal(names(res$cluster), rownames(mat))
})

test_that("dbscan_cluster returns dbscan object for sample clustering", {
  mat <- matrix(rnorm(100), nrow = 10)
  rownames(mat) <- paste0("Gene", 1:10)
  colnames(mat) <- paste0("Sample", 1:10)
  eset <- Biobase::ExpressionSet(assayData = mat)
  
  res <- dbscan_cluster(eset, eps = 0.5, minPts = 4, target = "sample")
  
  testthat::expect_type(res, "list")
  testthat::expect_named(res, c("cluster", "eps", "minPts", "distance", "target"))
  testthat::expect_equal(length(res$cluster), ncol(mat))               # number of samples
  testthat::expect_equal(names(res$cluster), colnames(mat))
})
